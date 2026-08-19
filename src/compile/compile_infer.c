/* Mathilda — Compile[]: the type-inference engine.
 *
 * A bottom-up walk over the Expr that computes each node's CompileType WITHOUT
 * emitting code — the emitter needs a result type before it lowers (e.g. to type
 * an If's result register before either branch is lowered, or to route an array
 * head).  It mirrors emit's result-type rules exactly and returns false for
 * anything outside the compilable subset.  Also here: the integer-closed head
 * table lookup (int_closed_head / int_only_head / all_args_int) that decides
 * when integer arguments keep a head's result an Integer, and the fixed-point
 * type helpers for the functional heads (Nest / Fold / FixedPoint / ...). */
#include "compile_emit.h"        /* Ctx / Val / FnSpec / LoopSpec / IntClosed + emit-core */
#include "../arithmetic.h"       /* is_rational */
#include "../symtab.h"           /* symtab_lookup, SymbolDef, ndarray_unary_kernel */
#include "compiled_function.h"   /* compiled_function_body — a CompiledFunction callee's body */
#include "../sym_names.h"        /* SYM_SameTest / SYM_List ... */
#include <string.h>
#include <stdlib.h>

/* Defined below infer_apply, which they all use; declared here because the
 * inference branches for Nest/Fold/FixedPoint/NestWhile need them. */
bool infer_apply(Ctx* c, const FnSpec* s, const CompileType* argt, int n,
                        CompileType* out);
int nest_fixed_type(Ctx* c, const FnSpec* s, CompileType t0);
int accum_fixed_type(Ctx* c, const FnSpec* s, CompileType t0,
                             const CompileType* rest, int nrest);
CompileType vec_elem_type(Ctx* c, const Expr* e);

/* Inference twin of emit_arr_unary: the array type a registered unary kernel
 * will produce over an operand of array type `ta`.  Both go through
 * nd_unary_elem, so the two cannot disagree about the element type -- which
 * they did, and which is what made Total[Floor[v]] answer with reinterpreted
 * integer bits. */
static bool infer_arr_unary(Ctx* c, const char* head, CompileType ta,
                            CompileType* out) {
    (void)c;
    SymbolDef* d = symtab_lookup(head);
    const NDUnaryKernel* k = d ? (const NDUnaryKernel*)d->ndarray_unary_kernel : NULL;
    CompileType er = nd_unary_elem(k, CT_ELEM(ta));
    if ((int)er < 0) return false;
    *out = CT_ARRAY(er, CT_RANK(ta));
    return true;
}
bool fp_opts(Expr* const* A, size_t na, size_t start,
                    const Expr** max_out, const Expr** same_out);

/* Pure type inference (no code emission) — needed to type an If's result
 * register before both branches are lowered.  Mirrors emit's result-type rules;
 * returns false for anything not compilable. */
/* Lift a scalar result element type back over an array operand.
 *
 * Several inference branches compute "the element type this head produces" and
 * then returned it directly, which silently turned an array-valued node into a
 * scalar one.  That was invisible while every array operation was delegated to
 * the ND layer — which decides its own result dtype — and became two distinct
 * bugs once fusion started trusting these types: the wrong output buffer, and
 * whole bodies (anything rooted at Sin/Exp/Log/ArcTan over an array) never
 * being recognised as fusable at all. */
static CompileType arr_like(CompileType operand, CompileType elem_result) {
    return CT_IS_ARRAY(operand) ? CT_ARRAY(elem_result, CT_RANK(operand)) : elem_result;
}
static CompileType elem_of(CompileType t) { return CT_IS_ARRAY(t) ? CT_ELEM(t) : t; }

/* Heads that are integer-CLOSED: given integer arguments the interpreter returns
 * an Integer, so a compiled body must too.  Every one of them also has a
 * registered real kernel — which is the trap, because reaching the kernel first
 * silently answers `35.` where the interpreter says `35`, and the numeric parity
 * tests cannot see the difference (see docs: result-HEAD parity).
 *
 * ONE table, consulted by both infer_type and emit, so the type a body is
 * declared to have and the opcode it actually runs cannot drift apart.
 *
 * Heads that look like they belong here and do NOT:
 *   LegendreP   `LegendreP[2, 2]` is 11/2 — Rational for some integer arguments.
 *   Beta        Rational.       HarmonicNumber  Rational.
 *   Divide      Rational; `Divide[7,3]` is 7/3 and no machine type holds it. */
static const IntClosed INT_CLOSED[] = {
    { "Factorial",  1, OP_FACT_I  },
    { "Gamma",      1, OP_GAMMA_I },
    { "Fibonacci",  1, OP_FIB_I   },
    { "LucasL",     1, OP_LUCAS_I },
    { "Binomial",   2, OP_BINOM_I },
    { "Pochhammer", 2, OP_POCH_I  },
};
const IntClosed* int_closed_head(const char* h, size_t na) {
    for (size_t i = 0; i < sizeof INT_CLOSED / sizeof INT_CLOSED[0]; i++)
        if (INT_CLOSED[i].arity == na && strcmp(INT_CLOSED[i].name, h) == 0)
            return &INT_CLOSED[i];
    return NULL;
}
/* Integer-ONLY heads: ones with no real counterpart at all.  Unlike INT_CLOSED
 * above there is no kernel behind them to fall through to — non-integer
 * arguments simply bail, exactly as they did before.
 *
 * `*pred` distinguishes the three that answer True/False (CT_BOOL) from the ones
 * that answer an Integer.  GCD and LCM are Flat and n-ary in the interpreter, so
 * any arity from 2 up is accepted and folded left. */
bool int_only_head(const char* h, size_t na, bool* pred) {
    *pred = false;
    if (na >= 2 && (strcmp(h, "GCD") == 0 || strcmp(h, "LCM") == 0)) return true;
    if ((na == 1 || na == 2)
        && (strcmp(h, "IntegerLength") == 0 || strcmp(h, "IntegerExponent") == 0)) return true;
    if (na == 3 && strcmp(h, "PowerMod") == 0) return true;
    if ((na == 1 && (strcmp(h, "EvenQ") == 0 || strcmp(h, "OddQ") == 0))
        || (na == 2 && strcmp(h, "Divisible") == 0)) { *pred = true; return true; }
    return false;
}

/* Peek at the argument types without emitting anything, so a head that is
 * integer-closed only on integer arguments can decline to the ordinary real
 * lowering for every other case. */
bool all_args_int(Ctx* c, Expr* const* A, size_t na) {
    for (size_t i = 0; i < na; i++) {
        CompileType t;
        if (!infer_type(c, A[i], &t) || t != CT_INT) return false;
    }
    return true;
}

bool infer_type(Ctx* c, const Expr* e, CompileType* out) {
    if (!e) return false;
    Slot imm; CompileType lt;
    if (literal(e, &imm, &lt)) { *out = lt; return true; }
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType st; if (scope_find(c, nm, &st, NULL) >= 0) { *out = st; return true; }
        int k = arg_find(c, nm);
        if (k >= 0) { *out = c->arg_types[k]; return true; }
        if (strcmp(nm, "True") == 0 || strcmp(nm, "False") == 0) { *out = CT_BOOL; return true; }
        if (strcmp(nm, "I") == 0) { *out = CT_COMPLEX; return true; }
        double cv; if (named_const(nm, &cv)) { *out = CT_REAL; return true; }
        Slot gi; CompileType gt;
        if (global_const(c, nm, &gi, &gt)) { *out = gt; return true; }
        return false;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    Expr** A = e->data.function.args; size_t na = e->data.function.arg_count;
    /* Slot[k] inside an inlined Function[body].  An index past the bound frame
     * falls through and fails, which is right: the interpreter would leave that
     * Slot unsubstituted, so its answer is not a machine number either. */
    if (h == SYM_Slot) { int k = fn_slot_index(c, A, na); if (k >= 0) { *out = c->slot[k].type; return true; } }
    /* Norm[{e1,...}, p] over a List literal -> scalar arithmetic: type the same
     * expansion emit_node lowers, so the two passes agree.  A declared-array
     * operand yields NULL and drops through to the array/reduction branches. */
    if (h == SYM_Norm) {
        Expr* ex = norm_try_expand(e);
        if (ex) { bool ok = infer_type(c, ex, out); expr_free(ex); return ok; }
    }
    /* Association read ops (B1): typed by the shared resolver.  Returns false for
     * a non-association Length/Values so the array branch below still handles it. */
    { CompileType at; if (try_infer_assoc(c, h, A, na, &at)) { *out = at; return true; } }
    CompileType ta, tb;
    #define IT(idx, dst) do { if (!infer_type(c, A[idx], &dst)) return false; } while (0)
    if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
        if (na == 0) { *out = CT_INT; return true; }
        IT(0, ta); for (size_t i = 1; i < na; i++) { IT(i, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; }
        *out = ta; return true;
    }
    if (strcmp(h, "Subtract") == 0 && na == 2) { IT(0, ta); IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; *out = ta; return true; }
    if (strcmp(h, "Minus") == 0 && na == 1)    { IT(0, ta); if (ta == CT_BOOL) return false; *out = ta; return true; }
    if (strcmp(h, "Divide") == 0 && na == 2)   { IT(0, ta); IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; if (ta < CT_REAL) ta = CT_REAL; *out = ta; return true; }
    if ((strcmp(h, "Mod") == 0 || strcmp(h, "Quotient") == 0) && na == 2) {
        IT(0, ta); IT(1, tb);
        /* An array operand is the registered binary kernel's business.  Both
         * heads have exact int64 arms and no complex one, so this branch used
         * to be the end of the road for them: Mod[x, 3] compiled and
         * Mod[v, 3] did not, and being outside the subset costs the whole
         * body, not just the Mod. */
        if (CT_IS_ARRAY(ta) || CT_IS_ARRAY(tb)) {
            if (CT_IS_ARRAY(ta) && CT_IS_ARRAY(tb)) return false;
            SymbolDef* dk = symtab_lookup(h);
            const NDBinaryKernel* kb = dk ? (const NDBinaryKernel*)dk->ndarray_binary_kernel : NULL;
            CompileType ea = CT_IS_ARRAY(ta) ? CT_ELEM(ta) : CT_ELEM(tb);
            CompileType es = CT_IS_ARRAY(ta) ? tb : ta;
            CompileType er = nd_binary_elem(kb, ea, es);
            if ((int)er < 0) return false;
            *out = CT_ARRAY(er, CT_RANK(CT_IS_ARRAY(ta) ? ta : tb));
            return true;
        }
        if (ta == CT_INT && tb == CT_INT) { *out = CT_INT; return true; }
        if (ta <= CT_REAL && tb <= CT_REAL && ta != CT_BOOL && tb != CT_BOOL) {
            /* Real Mod/Quotient: the interpreter evaluates both (Mod[2.5,1.2] is
             * 0.1, Quotient[7.5,2.] is 3), so declining here would drop a whole
             * body to the interpreter over an operation we can do. */
            *out = (h[0] == 'M') ? CT_REAL : CT_INT;
            return true;
        }
        return false;
    }
    if (strcmp(h, "Power") == 0 && na == 2) {
        IT(0, ta);
        if (A[1]->type == EXPR_INTEGER) {
            CompileType el = elem_of(ta);
            CompileType res = (el == CT_INT && A[1]->data.integer >= 0)
                            ? CT_INT : (el == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
            *out = arr_like(ta, res); return true;
        }
        int64_t rn, rd;
        if (is_rational(A[1], &rn, &rd) && rd == 2 && (rn == 1 || rn == -1)) {
            *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
            return true;
        }
        IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false;
        /* Integer base and integer exponent stay INTEGER, even though the
         * exponent's sign is only known per call: `2^n` is the Integer 1024 at
         * n = 10, and answering 1024. would be the wrong head.  OP_POW_II
         * abandons the call when the exponent turns out negative, because 2^-3
         * is the Rational 1/8 and no machine type holds it. */
        if (ta == CT_INT) { *out = CT_INT; return true; }
        if (ta < CT_REAL) ta = CT_REAL;
        *out = ta;
        return true;
    }
    uint16_t or_, oc_;
    if (na == 1 && (unary_math(h, &or_, &oc_) || strcmp(h, "Tanh") == 0)) {
        IT(0, ta);
        *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
        return true;
    }
    if (strcmp(h, "Log") == 0 && na == 2) {
        IT(0, ta); IT(1, tb); ta = num_common(ta, tb);
        if ((int)ta < 0) return false;
        *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
        return true;
    }
    /* Projections: the result is real even for a complex operand, and over an
     * array that means a REAL-element array, not an array of the operand's own
     * element type.  Getting this wrong is invisible while each operation is
     * delegated to the ND layer (which decides the result dtype itself), and
     * becomes a wrong answer the moment a fused loop allocates the output buffer
     * from this type instead. */
    if (strcmp(h, "Abs") == 0 && na == 1) {
        IT(0, ta);
        if (CT_IS_ARRAY(ta)) {
            CompileType el = CT_ELEM(ta) == CT_COMPLEX ? CT_REAL : CT_ELEM(ta);
            *out = CT_ARRAY(el, CT_RANK(ta)); return true;
        }
        *out = ta == CT_COMPLEX ? CT_REAL : ta; return true;
    }
    /* Integer-closed heads on integer arguments (see INT_CLOSED).  Must be tested
     * BEFORE the generic kernel dispatch below: each of these HAS a registered
     * real kernel, and reaching it first is exactly how `Binomial[7,3]` came back
     * as 35. instead of 35. */
    if (int_closed_head(h, na) && all_args_int(c, A, na)) { *out = CT_INT; return true; }
    /* Integer-ONLY heads: no real counterpart, so they exist over CT_INT or not
     * at all.  The three predicates answer True/False, hence CT_BOOL. */
    {
        bool pred = false;
        if (int_only_head(h, na, &pred) && all_args_int(c, A, na)) {
            *out = pred ? CT_BOOL : CT_INT; return true;
        }
    }
    /* Sign of a real is the INTEGER -1, 0 or 1 — `Sign[-2.5]` is `-1`, not
     * `-1.` — so the result type is CT_INT, not the argument's type.  Sign of a
     * complex is z/|z|, which is genuinely complex. */
    if (strcmp(h, "Sign") == 0 && na == 1) {
        IT(0, ta);
        if (CT_IS_ARRAY(ta)) return infer_arr_unary(c, h, ta, out);
        *out = (ta == CT_COMPLEX) ? CT_COMPLEX : CT_INT;
        return true;
    }
    if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0
         || strcmp(h, "IntegerPart") == 0) && na == 1) {
        IT(0, ta);
        /* The ARRAY case used to `return false` here, so a narrowing head could
         * only ever be a whole compiled body and never a subexpression:
         * Total[Floor[v]] had no inferable type.  It now reports what the
         * kernel will really produce, which for these four is an int64 buffer
         * — the same answer emit_arr_unary reaches. */
        if (CT_IS_ARRAY(ta)) return infer_arr_unary(c, h, ta, out);
        *out = CT_INT; return true; }
    if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 || strcmp(h, "Arg") == 0) && na == 1) {
        IT(0, ta);
        if (ta == CT_BOOL) return false;
        /* Over the integers all three stay integral: Re is the value, Im is 0,
         * and Arg is 0 (or Pi, which ARG_I hands back to the interpreter). */
        if (ta == CT_INT) { *out = CT_INT; return true; }
        *out = CT_IS_ARRAY(ta) ? CT_ARRAY(CT_REAL, CT_RANK(ta)) : CT_REAL;
        return true;
    }
    /* FractionalPart of an integer is the Integer 0, not 0. */
    if (strcmp(h, "FractionalPart") == 0 && na == 1) {
        IT(0, ta);
        if (ta == CT_INT) { *out = CT_INT; return true; }
        if (CT_IS_ARRAY(ta) || ta == CT_BOOL) return false;
        *out = ta; return true;
    }
    if (strcmp(h, "Conjugate") == 0 && na == 1) { IT(0, ta); *out = ta; return true; }
    if (strcmp(h, "HypergeometricPFQ") == 0 && na == 3) {
        const Expr* la = A[0];
        const Expr* lb = A[1];
        if (la->type != EXPR_FUNCTION || la->data.function.head->type != EXPR_SYMBOL
            || strcmp(la->data.function.head->data.symbol.name, "List") != 0
            || lb->type != EXPR_FUNCTION || lb->data.function.head->type != EXPR_SYMBOL
            || strcmp(lb->data.function.head->data.symbol.name, "List") != 0) return false;
        size_t np = la->data.function.arg_count, nq = lb->data.function.arg_count;
        if (np + nq + 2 > 8) return false;
        for (size_t i = 0; i < np; i++) {
            CompileType t;
            if (!infer_type(c, la->data.function.args[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
        }
        for (size_t j = 0; j < nq; j++) {
            CompileType t;
            if (!infer_type(c, lb->data.function.args[j], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
        }
        IT(2, ta); if (CT_IS_ARRAY(ta) || ta == CT_BOOL) return false;
        *out = CT_REAL; return true;
    }
    if (strcmp(h, "UnitStep") == 0 && na >= 1) {
        /* One array argument is the kernel's business, not the product-of-
         * comparisons lowering below: UnitStep's kernel is narrowing-only (no
         * real and no complex arm at all), so before this the array form had
         * no lowering and took the whole body down with it. */
        if (na == 1) { IT(0, ta); if (CT_IS_ARRAY(ta)) return infer_arr_unary(c, h, ta, out); }
        for (size_t i = 0; i < na; i++) { IT(i, ta); if (CT_IS_ARRAY(ta) || ta == CT_BOOL) return false; }
        *out = CT_INT; return true;                    /* UnitStep[0.5] is 1, not 1. */
    }
    if (strcmp(h, "UnitBox") == 0 && na == 1) {
        /* Single-argument box, matching the interpreter (multi-arg UnitBox is
         * left unevaluated there). An array argument routes to the narrowing
         * kernel; a scalar lowers to the pair of comparisons in the emitter and
         * answers with an exact Integer 0/1, so CT_INT like UnitStep. */
        IT(0, ta);
        if (CT_IS_ARRAY(ta)) return infer_arr_unary(c, h, ta, out);
        if (ta == CT_BOOL || ta == CT_COMPLEX) return false;
        *out = CT_INT; return true;                    /* UnitBox[0.5] is 1, not 1. */
    }
    if (strcmp(h, "Chop") == 0 && (na == 1 || na == 2)) {
        IT(0, ta);
        if (CT_IS_ARRAY(ta) || ta == CT_BOOL || ta == CT_COMPLEX) return false;
        if (na == 2 && A[1]->type != EXPR_REAL && A[1]->type != EXPR_INTEGER) return false;
        *out = ta;   /* CT_INT -> identity; CT_REAL -> Real */
        return true;
    }
    if (strcmp(h, "Clip") == 0 && na == 1) {          /* default bounds {-1, 1} */
        IT(0, ta);
        if (CT_IS_ARRAY(ta) || ta == CT_BOOL || ta == CT_COMPLEX) return false;
        *out = CT_REAL; return true;
    }
    if ((strcmp(h, "Clip") == 0 || strcmp(h, "Rescale") == 0) && na == 2) {
        const Expr* bnd = A[1];
        if (bnd->type != EXPR_FUNCTION || bnd->data.function.head->type != EXPR_SYMBOL
            || strcmp(bnd->data.function.head->data.symbol.name, "List") != 0
            || bnd->data.function.arg_count != 2) return false;
        IT(0, ta);
        CompileType tl, tu;
        if (!infer_type(c, bnd->data.function.args[0], &tl)
            || !infer_type(c, bnd->data.function.args[1], &tu)) return false;
        if (tl != CT_REAL || tu != CT_REAL || CT_IS_ARRAY(ta) || ta == CT_BOOL) return false;
        /* Rescale carries a complex argument through; Clip cannot (Min/Max need
         * an order, and the interpreter leaves complex Clip unevaluated). */
        if (ta == CT_COMPLEX) {
            if (h[0] == 'C') return false;
            *out = CT_COMPLEX; return true;
        }
        *out = CT_REAL; return true;
    }
    /* Max / Min over SCALARS fold pairwise.  An array operand is a different
     * operation entirely -- Max[v] is the largest ELEMENT, and Max[v, 3] is the
     * largest of the flattened whole -- so every array spelling is refused
     * here and left to the delegated reduction below (na == 1) or to the
     * interpreter (na >= 2).
     *
     * Refusing it is a CORRECTNESS fix, not a coverage one.  The pairwise fold
     * used to accept an array unchecked, so `Max[v]` lowered to the identity
     * (an empty fold returns its accumulator) and
     *     Compile[{{v, _Real, 1}}, Max[v] + 1.][{3., 1., 7., 2.}]
     * answered {4., 2., 8., 3.} where the interpreter answers 8.  Alone,
     * `Max[v]` happened to be rejected downstream, which is why nothing had
     * noticed: the wrong answer needed the head to appear inside a larger
     * expression.  Found by tools/compile_coverage.py on 2026-08-02. */
    if ((strcmp(h, "Max") == 0 || strcmp(h, "Min") == 0) && na >= 2) {
        IT(0, ta);
        if (CT_IS_ARRAY(ta)) return false;
        for (size_t i = 1; i < na; i++) {
            IT(i, tb);
            if (CT_IS_ARRAY(tb)) return false;
            ta = num_common(ta, tb);
            if ((int)ta < 0 || ta == CT_COMPLEX) return false;
        }
        *out = ta; return true;
    }
    if ((strcmp(h, "Max") == 0 || strcmp(h, "Min") == 0) && na == 1) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta)) { *out = ta; return true; }   /* Max[x] is x */
        /* falls through to the delegated-reduction block below */
    }
    if (strcmp(h, "ArcTan") == 0) {
        if (na == 1) {
            IT(0, ta);
            *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
            return true;
        }
        if (na == 2) {
            IT(0, ta); IT(1, tb); ta = num_common(ta, tb);
            if ((int)ta < 0) return false;
            *out = arr_like(ta, CT_REAL);
            return true;
        }
        return false;
    }
    if (na == 2 && (!strcmp(h, "Less") || !strcmp(h, "LessEqual") || !strcmp(h, "Greater") || !strcmp(h, "GreaterEqual") || !strcmp(h, "Equal") || !strcmp(h, "Unequal"))) { *out = CT_BOOL; return true; }
    /* Order[a,b] is the canonical comparison {1,0,-1}: always the INTEGER, like
     * Sign.  Only ordered scalars (Integer/Real, unified) compile; complex,
     * boolean and array args fall back to the interpreter. */
    if (na == 2 && !strcmp(h, "Order")) { IT(0, ta); IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0 || ta == CT_COMPLEX || CT_IS_ARRAY(ta)) return false; *out = CT_INT; return true; }
    if (!strcmp(h, "And") || !strcmp(h, "Or") || !strcmp(h, "Xor") || !strcmp(h, "Not")) { *out = CT_BOOL; return true; }
    if (strcmp(h, "If") == 0 && na == 3) {
        CompileType tt, te; if (!infer_type(c, A[1], &tt) || !infer_type(c, A[2], &te)) return false;
        if (tt == te) { *out = tt; return true; }
        *out = num_common(tt, te); return (int)*out >= 0;
    }
    /* Which / Switch / Piecewise: the ladder's result type, computed the same way
     * emit_ctrl does (shared collect_ladder), so the two cannot disagree. */
    if (strcmp(h, "Which") == 0 || strcmp(h, "Switch") == 0 || strcmp(h, "Piecewise") == 0) {
        int nr; bool hd; CompileType rt;
        if (!collect_ladder(c, h, A, na, NULL, &nr, &hd, &rt) || CT_IS_ARRAY(rt)) return false;
        *out = rt; return true;
    }
    if ((strcmp(h, "Sum") == 0 || strcmp(h, "Product") == 0) && na == 2) {
        LoopSpec s;
        /* s.var != NULL: the interpreter's Sum/Product reject a bare count
         * ({n}), so the compiled path must reject it too — only Do accepts it. */
        if (!loop_spec_parse(A[1], &s) || !s.var || !loop_spec_int_bounds(c, &s)
            || c->nscope >= CTX_MAX_SCOPE) return false;
        c->scope[c->nscope].name = s.var;
        c->scope[c->nscope].reg = 0; c->scope[c->nscope].type = CT_INT; c->nscope++;
        CompileType T; bool okT = infer_type(c, A[0], &T);
        c->nscope--;
        if (!okT || T == CT_BOOL) return false;
        *out = T; return true;
    }
    if (strcmp(h, "CompoundExpression") == 0 && na >= 1) return infer_type(c, A[na - 1], out);
    if ((strcmp(h, "With") == 0 || strcmp(h, "Module") == 0) && na == 2) {
        const Expr* L = A[0];
        if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
            || strcmp(L->data.function.head->data.symbol.name, "List") != 0
            || L->data.function.arg_count == 0
            || (int)(c->nscope + L->data.function.arg_count) > CTX_MAX_SCOPE) return false;
        int pushed = 0;
        for (size_t i = 0; i < L->data.function.arg_count; i++) {
            const Expr* spec = L->data.function.args[i];
            const char* vname = NULL; CompileType vt = CT_REAL;
            if (spec->type == EXPR_SYMBOL) vname = spec->data.symbol.name;
            else if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
                     && strcmp(spec->data.function.head->data.symbol.name, "Set") == 0
                     && spec->data.function.arg_count == 2 && spec->data.function.args[0]->type == EXPR_SYMBOL) {
                vname = spec->data.function.args[0]->data.symbol.name;
                if (!infer_type(c, spec->data.function.args[1], &vt)) { c->nscope -= pushed; return false; }
            } else { c->nscope -= pushed; return false; }
            c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = 0; c->scope[c->nscope].type = vt;
            c->nscope++; pushed++;
        }
        bool okb = infer_type(c, A[1], out);
        c->nscope -= pushed; return okb;
    }
    if ((!strcmp(h, "Set") || !strcmp(h, "AddTo") || !strcmp(h, "SubtractFrom")
         || !strcmp(h, "TimesBy") || !strcmp(h, "DivideBy")) && na == 2
        && A[0]->type == EXPR_SYMBOL) {
        CompileType vt; if (scope_find(c, A[0]->data.symbol.name, &vt, NULL) < 0) return false; *out = vt; return true;
    }
    /* Threaded list assignment {a, b, ...} = scalar: its value is the RHS (which
     * emit_list_thread broadcasts to each target). A List RHS (destructuring) or
     * an array RHS is not lowered -- report non-compilable so the whole body
     * bails to the interpreter, matching emit_ctrl. */
    if (!strcmp(h, "Set") && na == 2 && A[0]->type == EXPR_FUNCTION
        && A[0]->data.function.head->type == EXPR_SYMBOL
        && !strcmp(A[0]->data.function.head->data.symbol.name, "List")) {
        if (A[1]->type == EXPR_FUNCTION && A[1]->data.function.head->type == EXPR_SYMBOL
            && !strcmp(A[1]->data.function.head->data.symbol.name, "List")) return false;
        CompileType vt;
        if (!infer_type(c, A[1], &vt) || CT_IS_ARRAY(vt)) return false;
        *out = vt; return true;
    }
    if ((!strcmp(h, "Increment") || !strcmp(h, "Decrement")) && na == 1 && A[0]->type == EXPR_SYMBOL) {
        CompileType vt; if (scope_find(c, A[0]->data.symbol.name, &vt, NULL) < 0) return false; *out = vt; return true;
    }
    if ((!strcmp(h, "Do") && na >= 2) || (!strcmp(h, "While") && na == 2) || (!strcmp(h, "For") && na == 4)
        || (!strcmp(h, "If") && na == 2)) { *out = CT_INT; return true; }
    if (!strcmp(h, "Nest") && na == 3) {
        FnSpec s; CompileType tn, tx;
        if (!fn_resolve(A[0], 1, &s) || !infer_type(c, A[2], &tn) || tn != CT_INT
            || !infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0) return false;
        *out = (CompileType)tfp; return true;
    }
    /* Delegated structural heads — the rank rule lives with the table. */
    {
        const NdFnSpec* nf = nd_fn_lookup(h, na);
        if (nf) {
            IT(0, ta);
            for (int i = 0; i < nf->nextra; i++) {
                CompileType te;
                if (!infer_type(c, A[1 + i], &te) || te != CT_INT) return false;
            }
            CompileType rt = nd_fn_result(nf, ta);
            if ((int)rt < 0) return false;
            *out = rt; return true;
        }
    }
    /* Delegated two-array heads (COMPILE_MISSING.md §3).  Required, not optional:
     * compile_expr emits without pre-inferring, so a STANDALONE Dot[m,v] would
     * compile from the emit branch alone — but a COMPOSED body (Total[Dot[m,v]],
     * a.b + 1.0) infers its subterms, and without this it would bail wholesale,
     * which is exactly the cliff this lowering exists to remove. */
    {
        const NdFn2Spec* nf2 = nd_fn2_lookup(h, na);
        if (nf2) {
            CompileType ta, tb;
            if (!infer_type(c, A[0], &ta) || !infer_type(c, A[1], &tb)) return false;
            CompileType rt = nd_fn2_result(nf2, ta, tb);
            if ((int)rt < 0) return false;
            *out = rt; return true;
        }
    }
    /* Select / TakeWhile / LengthWhile / All-, Any-, NoneTrue — one predicate
     * loop; see the emit-side block. */
    if (na == 2 && (!strcmp(h, "Select") || !strcmp(h, "TakeWhile")
                    || !strcmp(h, "LengthWhile") || !strcmp(h, "AllTrue")
                    || !strcmp(h, "AnyTrue") || !strcmp(h, "NoneTrue"))) {
        FnSpec s; if (!fn_resolve(A[1], 1, &s)) return false;
        CompileType el = vec_elem_type(c, A[0]);
        if ((int)el < 0) return false;
        CompileType tb;
        if (!infer_apply(c, &s, &el, 1, &tb) || tb != CT_BOOL) return false;
        if (!strcmp(h, "Select") || !strcmp(h, "TakeWhile")) *out = CT_ARRAY(el, 1);
        else if (!strcmp(h, "LengthWhile"))                  *out = CT_INT;
        else                                                 *out = CT_BOOL;
        return true;
    }
    /* First[v] / Last[v] are v[[1]] and v[[-1]]. */
    if ((!strcmp(h, "First") || !strcmp(h, "Last")) && na == 1) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta)) return false;
        *out = CT_RANK(ta) == 1 ? CT_ELEM(ta) : CT_ARRAY(CT_ELEM(ta), CT_RANK(ta) - 1);
        return true;
    }
    /* Map / Scan over a rank-1 array — see the emit-side block for why the
     * result element type must equal the source's. */
    if ((!strcmp(h, "Map") || !strcmp(h, "Scan")) && na == 2) {
        bool scan = h[0] == 'S';
        FnSpec s; if (!fn_resolve(A[0], 1, &s)) return false;
        CompileType el = vec_elem_type(c, A[1]);
        if ((int)el < 0) return false;
        CompileType rt;
        if (!infer_apply(c, &s, &el, 1, &rt) || CT_IS_ARRAY(rt)) return false;
        if (scan) { *out = CT_INT; return true; }
        if (rt != el) return false;
        *out = CT_ARRAY(rt, 1); return true;
    }
    /* Fold[f, x0, v] / Fold[f, v]: the accumulator's widening fixed point with
     * the element type held fixed. */
    if (!strcmp(h, "Fold") && (na == 2 || na == 3)) {
        FnSpec s; if (!fn_resolve(A[0], 2, &s)) return false;
        CompileType el = vec_elem_type(c, A[na - 1]);
        if ((int)el < 0) return false;
        CompileType t0 = el;                       /* Fold[f, v] seeds from v[[1]] */
        if (na == 3 && !infer_type(c, A[1], &t0)) return false;
        int tfp = accum_fixed_type(c, &s, t0, &el, 1);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) return false;
        *out = (CompileType)tfp; return true;
    }
    /* NestList / FoldList — see the emit-side block for the element-type rule. */
    if ((!strcmp(h, "NestList") && na == 3)
        || (!strcmp(h, "FoldList") && (na == 2 || na == 3))) {
        bool nest = h[0] == 'N';
        FnSpec s; if (!fn_resolve(A[0], nest ? 1 : 2, &s)) return false;
        CompileType T;
        if (nest) {
            CompileType tn, tx;
            if (!infer_type(c, A[2], &tn) || tn != CT_INT
                || !infer_type(c, A[1], &tx)) return false;
            int tfp = nest_fixed_type(c, &s, tx);
            if (tfp < 0) return false;
            T = (CompileType)tfp;
        } else {
            CompileType el = vec_elem_type(c, A[na - 1]);
            if ((int)el < 0) return false;
            CompileType t0 = el;
            if (na == 3 && !infer_type(c, A[1], &t0)) return false;
            int tfp = accum_fixed_type(c, &s, t0, &el, 1);
            if (tfp < 0) return false;
            T = (CompileType)tfp;
        }
        if (!ct_is_elem(T)) return false;
        *out = CT_ARRAY(T, 1); return true;
    }
    /* FixedPointList / NestWhileList — a BUILT history, so Real/Complex only. */
    if ((!strcmp(h, "FixedPointList") && na >= 2 && na <= 4)
        || (!strcmp(h, "NestWhileList") && na == 3)) {
        bool fp = h[0] == 'F';
        FnSpec s, ss, ts; const Expr *mx = NULL, *st = NULL;
        if (!fn_resolve(A[0], 1, &s)) return false;
        if (fp) { if (!fp_opts(A, na, 2, &mx, &st)) return false; }
        else if (!fn_resolve(A[2], 1, &ts)) return false;
        CompileType tn;
        if (mx && (!infer_type(c, mx, &tn) || tn != CT_INT)) return false;
        CompileType tx; if (!infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0) return false;
        CompileType T = (CompileType)tfp;
        if (!ct_is_elem(T)) return false;
        CompileType tb;
        if (fp && st) {
            CompileType at[2] = { T, T };
            if (!fn_resolve(st, 2, &ss) || !infer_apply(c, &ss, at, 2, &tb) || tb != CT_BOOL)
                return false;
        }
        if (!fp && (!infer_apply(c, &ts, &T, 1, &tb) || tb != CT_BOOL)) return false;
        *out = CT_ARRAY(T, 1); return true;
    }
    if (!strcmp(h, "FixedPoint") && na >= 2 && na <= 4) {
        FnSpec s; const Expr *mx, *st;
        if (!fn_resolve(A[0], 1, &s) || !fp_opts(A, na, 2, &mx, &st)) return false;
        CompileType tn;
        if (mx && (!infer_type(c, mx, &tn) || tn != CT_INT)) return false;
        CompileType tx; if (!infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) return false;
        if (st) {                                  /* SameTest must yield a Bool */
            FnSpec ss; CompileType tb, at[2] = { (CompileType)tfp, (CompileType)tfp };
            if (!fn_resolve(st, 2, &ss) || !infer_apply(c, &ss, at, 2, &tb) || tb != CT_BOOL)
                return false;
        }
        *out = (CompileType)tfp; return true;
    }
    if (!strcmp(h, "NestWhile") && na == 3) {
        FnSpec s, ts; CompileType tx, tb;
        if (!fn_resolve(A[0], 1, &s) || !fn_resolve(A[2], 1, &ts)
            || !infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) return false;
        CompileType at = (CompileType)tfp;
        if (!infer_apply(c, &ts, &at, 1, &tb) || tb != CT_BOOL) return false;
        *out = at; return true;
    }
    /* Part (M3c).  A full run of scalar subscripts drops every axis and lands
     * back in the scalar lattice; anything else keeps at least one axis and
     * stays an array, with the rank counted the way build_axis_selector counts
     * it: a spec keeps its axis, a scalar subscript drops it, and axes past the
     * last subscript are implicit Alls and survive. */
    if (strcmp(h, "Part") == 0 && na >= 2) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta) || (int)(na - 1) > CT_RANK(ta)) return false;
        if (part_is_scalar_indexed(c, ta, (const Expr* const*)A + 1, na - 1)) {
            *out = CT_ELEM(ta); return true;
        }
        int keep = CT_RANK(ta) - (int)(na - 1);
        for (size_t i = 1; i < na; i++) if (subscript_is_literal_spec(A[i])) keep++;
        if (keep < 1 || keep > CT_MAX_RANK) return false;
        *out = CT_ARRAY(CT_ELEM(ta), keep);
        return true;
    }
    if (strcmp(h, "ConstantArray") == 0 && na == 2) {
        int rank; CompileType elem;
        if (!const_array_shape(c, (const Expr* const*)A, &rank, &elem)) return false;
        *out = CT_ARRAY(elem, rank); return true;
    }
    /* Range[hi] / Range[lo, hi] / Range[lo, hi, di] — see the emit-side block. */
    if (strcmp(h, "Range") == 0 && na >= 1 && na <= 3) {
        LoopSpec s;
        if (!range_spec(c, A, na, &s)) return false;
        *out = CT_ARRAY(CT_INT, 1); return true;
    }
    /* Table[body, spec...] — see the emit-side block for why the iterators must
     * be integer-bounded and the element type Real or Complex. */
    if (strcmp(h, "Table") == 0 && na >= 2 && (int)na - 1 <= CT_MAX_RANK) {
        int rank = (int)na - 1;
        LoopSpec sp[CT_MAX_RANK];
        for (int j = 0; j < rank; j++)
            if (!loop_spec_parse(A[j + 1], &sp[j]) || !loop_spec_int_bounds(c, &sp[j])) return false;
        if (c->nscope + rank > CTX_MAX_SCOPE) return false;
        int pushed = 0;
        for (int j = 0; j < rank; j++)
            if (sp[j].var) {
                c->scope[c->nscope].name = sp[j].var; c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = CT_INT; c->scope[c->nscope].built = false;
                c->nscope++; pushed++;
            }
        CompileType et; bool okT = infer_type(c, A[0], &et);
        c->nscope -= pushed;
        if (!okT || !ct_is_elem(et)) return false;
        *out = CT_ARRAY(et, rank); return true;
    }

    /* array -> scalar reductions (M3a): the only way an array type re-enters
     * the scalar lattice, so these must be inferable inside If/Sum/... */
    if ((strcmp(h, "Total") == 0 || strcmp(h, "Length") == 0) && na == 1) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta) || (h[0] == 'T' && CT_RANK(ta) != 1)) return false;
        *out = (h[0] == 'L') ? CT_INT : CT_ELEM(ta);
        return true;
    }
    /* Norm[array, p] two-arg form (V_NORM): valid literal p over a real array
     * operand -> Real.  The bare Norm[array] is the reduction below. */
    if (strcmp(h, "Norm") == 0 && na == 2) {
        CompileType ta; IT(0, ta);
        double immr;
        if (!nd_norm2_encode(A[1], ta, &immr)) return false;
        *out = CT_REAL; return true;
    }
    /* Mean / Median / Variance / StandardDeviation / RootMeanSquare / Max /
     * Min of a rank-1 array — the delegated reductions, table above. */
    {
        const NdRedSpec* nr = nd_red_lookup(h, na);
        if (nr) {
            IT(0, ta);
            if (nr->nextra == 1) {   /* RankedMin/Max: trailing rank must be an Integer */
                CompileType tk; IT(1, tk);
                if (tk != CT_INT) return false;
            }
            CompileType rt = nd_red_result(nr, ta);
            if ((int)rt < 0) return false;
            *out = rt; return true;
        }
    }
    if (na == 1) { SymbolDef* d = symtab_lookup(h); if (d && d->ndarray_unary_kernel) { const NDUnaryKernel* k = d->ndarray_unary_kernel; if (k->cplx || k->real) { IT(0, ta); if (CT_IS_ARRAY(ta)) { *out = CT_ARRAY(k->to_real ? CT_REAL : CT_ELEM(ta), CT_RANK(ta)); return true; } if (k->to_real) { *out = CT_REAL; return true; } if (ta == CT_COMPLEX) { if (!k->cplx) return false; *out = CT_COMPLEX; return true; } *out = (k->real_closed || k->real) ? CT_REAL : CT_COMPLEX; return true; }
        /* Narrowing-only kernel over an array: NDT_INT64 out.  Mirrors the
         * branch in emit_arr_unary, and must agree with it or a body that
         * infers here and emits there disagrees about the element type. */
        if (k->to_int_r) { IT(0, ta); if (CT_IS_ARRAY(ta) && CT_ELEM(ta) != CT_COMPLEX) { *out = CT_ARRAY(CT_INT, CT_RANK(ta)); return true; } } } }
    if (na >= 3 && na <= 8) {
        SymbolDef* d = symtab_lookup(h);
        if (d && d->ndarray_nary_kernel) {
            const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
            if (k->cplx && k->nargs == na) {
                for (size_t i = 0; i < na; i++) {
                    CompileType t;
                    if (!infer_type(c, A[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
                }
                *out = CT_REAL; return true;
            }
        }
    }
    if (na == 2) {
        SymbolDef* d = symtab_lookup(h);
        if (d && d->ndarray_binary_kernel) {
            const NDBinaryKernel* k = d->ndarray_binary_kernel;
            if (k->cplx || k->to_int) {
                IT(0, ta); IT(1, tb);
                CompileType t = num_common(ta, tb);
                if ((int)t < 0) return false;
                /* Over an array the result is an ARRAY, mirroring the unary
                 * branch above and the emit-side lowering.  Reporting a scalar
                 * here was harmless while every array op was delegated (the ND
                 * layer picks its own result dtype), and becomes a wrong output
                 * buffer as soon as a fused loop is sized from this type. */
                if (CT_IS_ARRAY(t)) {
                    if (CT_IS_ARRAY(ta) && CT_IS_ARRAY(tb)) return false;
                    CompileType ea = CT_IS_ARRAY(ta) ? CT_ELEM(ta) : CT_ELEM(tb);
                    CompileType es = CT_IS_ARRAY(ta) ? tb : ta;
                    CompileType er = nd_binary_elem(k, ea, es);
                    if ((int)er < 0) return false;
                    *out = CT_ARRAY(er, CT_RANK(t));
                    return true;
                }
                if (!k->cplx) return false;    /* scalar: no arm to call */
                *out = (t <= CT_REAL && k->real_closed) ? CT_REAL : CT_COMPLEX;
                return true;
            }
        }
    }
    #undef IT
    /* Inlined CompiledFunction call — must type here too, or a call inside an
     * If branch / Sum body could not be lowered.  Mirrors emit's binding. */
    {
        const CompiledFunction* cf = compiled_callee(c, h);
        if (cf && na >= 1 && compiled_function_num_args(cf) == na
            && na + (size_t)c->nscope <= CTX_MAX_SCOPE && c->inlining < 8) {
            const char* const* pn = compiled_function_arg_names(cf);
            const CompileType* pt = compiled_function_arg_types(cf);
            /* Type every argument in the CALLER's environment first — binding a
             * parameter before the later arguments are typed would let it
             * shadow a caller variable of the same name. */
            for (size_t i = 0; i < na; i++) {
                CompileType ta;
                if (!infer_type(c, A[i], &ta) || CT_IS_ARRAY(ta) || CT_IS_ARRAY(pt[i]))
                    return false;
            }
            int saved_scope = c->nscope;
            for (size_t i = 0; i < na; i++) {
                c->scope[c->nscope].name = pn[i];
                c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = pt[i];
                c->nscope++;
            }
            c->inlining++;
            bool okb = infer_type(c, compiled_function_body(cf), out);
            c->inlining--;
            c->nscope = saved_scope;
            return okb && !CT_IS_ARRAY(*out);
        }
    }
    return false;
}

/* Type `f[a1, ..., an]` for an already-resolved function value.  Mirrors
 * emit_apply's binding exactly — the two must agree or a construct types as one
 * thing and lowers as another. */
bool infer_apply(Ctx* c, const FnSpec* s, const CompileType* argt, int n,
                        CompileType* out) {
    switch (s->kind) {
        case FN_IDENTITY:
            if (n != 1) return false;
            *out = argt[0]; return true;

        case FN_LAMBDA: {
            if (n != s->nparams || c->nscope + n > CTX_MAX_SCOPE) return false;
            int saved_scope = c->nscope, saved_nslot = c->nslot;
            for (int i = 0; i < n; i++) {
                c->scope[c->nscope].name = s->pname[i];
                c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = argt[i];
                c->nscope++;
            }
            c->nslot = 0;                    /* see emit_apply for why */
            bool ok = infer_type(c, s->body, out);
            c->nslot = saved_nslot; c->nscope = saved_scope;
            return ok;
        }

        case FN_SLOTS: {
            if (n > FN_MAX_PARAMS) return false;
            int saved_nslot = c->nslot;
            struct { int reg; CompileType type; } saved[FN_MAX_PARAMS];
            memcpy(saved, c->slot, sizeof saved);
            for (int i = 0; i < n; i++) { c->slot[i].reg = 0; c->slot[i].type = argt[i]; }
            c->nslot = n;
            bool ok = infer_type(c, s->body, out);
            c->nslot = saved_nslot;
            memcpy(c->slot, saved, sizeof saved);
            return ok;
        }

        case FN_COMPOSE: {
            /* Composition[f1,...,fk][a...] = f1[f2[...fk[a...]]] — the INNERMOST
             * takes every argument, each outer takes one (src/eval.c:1449). */
            Expr* const* F = s->fexpr->data.function.args;
            size_t nf = s->fexpr->data.function.arg_count;
            CompileType t; FnSpec inner;
            if (!fn_resolve(F[nf - 1], n, &inner) || !infer_apply(c, &inner, argt, n, &t))
                return false;
            for (size_t k = nf - 1; k > 0; k--) {
                FnSpec g;
                if (!fn_resolve(F[k - 1], 1, &g) || !infer_apply(c, &g, &t, 1, &t)) return false;
            }
            *out = t; return true;
        }

        case FN_HEAD: {
            if (n > FN_MAX_PARAMS || c->nscope + n > CTX_MAX_SCOPE) return false;
            Expr* call = fn_head_call(s->head, n);
            if (!call) return false;
            int saved_scope = c->nscope;
            for (int i = 0; i < n; i++) {
                c->scope[c->nscope].name = fn_placeholder(i);
                c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = argt[i];
                c->nscope++;
            }
            bool ok = infer_type(c, call, out);
            c->nscope = saved_scope;
            expr_free(call);
            return ok;
        }
    }
    return false;
}

/* Fixed-point accumulator type for an iteration that feeds its own output back
 * in (Nest, Fold, FixedPoint, NestWhile): the register must hold a type wide
 * enough to absorb every iteration.  Starting from t0, widen until the output
 * type stops growing — the lattice is bounded, so this converges in a few
 * passes or not at all.
 *
 * Argument 0 is the accumulator; `rest` carries the types of any further
 * arguments, which do not vary (Fold's list element). */
int accum_fixed_type(Ctx* c, const FnSpec* s, CompileType t0,
                            const CompileType* rest, int nrest) {
    if (nrest < 0 || nrest + 1 > FN_MAX_PARAMS) return -1;
    CompileType t = t0;
    for (int iter = 0; iter < 4; iter++) {
        CompileType at[FN_MAX_PARAMS], tb;
        at[0] = t;
        for (int i = 0; i < nrest; i++) at[i + 1] = rest[i];
        if (!infer_apply(c, s, at, nrest + 1, &tb) || (int)tb < 0) return -1;
        if (tb == CT_BOOL && t != CT_BOOL) return -1;    /* can't fold a bool into a number */
        if (tb <= t) return (int)t;                       /* output coerces down into the accumulator */
        t = tb;                                           /* output widened it; grow and re-check */
    }
    return -1;
}
int nest_fixed_type(Ctx* c, const FnSpec* s, CompileType t0) {
    return accum_fixed_type(c, s, t0, NULL, 0);
}

/* The element type of a rank-1 array-valued expression, or CT_ERR.  Fold, Map
 * and the other list-consuming heads all need exactly this test. */
CompileType vec_elem_type(Ctx* c, const Expr* e) {
    CompileType t;
    if (!infer_type(c, e, &t) || !CT_IS_ARRAY(t) || CT_RANK(t) != 1) return CT_ERR;
    return CT_ELEM(t);
}

/* FixedPoint's trailing arguments: an application bound and/or SameTest -> s,
 * in either order, exactly as parse_fp_opts accepts them (src/funcprog.c:2511).
 * Returns false for a spelling the interpreter itself refuses. */
bool fp_opts(Expr* const* A, size_t na, size_t start,
                    const Expr** max_out, const Expr** same_out) {
    *max_out = NULL; *same_out = NULL;
    for (size_t i = start; i < na; i++) {
        const Expr* a = A[i];
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL
            && a->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            if (*same_out) return false;
            *same_out = a->data.function.args[1];
        } else if (!*max_out) {
            *max_out = a;        /* the bound; must infer as CT_INT to compile */
        } else return false;
    }
    return true;
}
