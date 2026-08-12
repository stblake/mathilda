/* Mathilda — Compile[]: elementary/special-function + projection heads (Log, Sin..Tanh, Abs/Sign/Floor/Ceiling/Round/Re/Im/Arg/Conjugate/IntegerPart/FractionalPart/UnitStep, Clip/Rescale, HypergeometricPFQ, Max/Min, ArcTan).
 *
 * One of emit_node's per-category dispatch functions: 1 = handled (ok iff
 * c->ok), -1 = handled but not lowerable, 0 = not one of these heads. */
#include "compile_emit.h"
#include "../arithmetic.h"
#include "../symtab.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

int emit_mathfn(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e; Slot imm;
    /* elementary unary math */
    uint16_t op_r, op_c;
    if (na == 1 && unary_math(h, &op_r, &op_c)) return emit_unary_math(c, h, A[0], op_r, op_c, out);
    if (strcmp(h, "Log") == 0 && na == 2) {   /* Log[b,x] = Log[x]/Log[b] */
        /* over an array, the ND layer's registered two-arg Log kernel maps the
         * whole buffer in one pass — the scalar lowering below would need an
         * array divide the ND layer does not have. */
        Val kv2;
        if (any_array_arg(c, A, na)) {
            if (try_kernel(c, h, A, na, &kv2)) { *out = kv2; return c->ok ? 1 : -1; }
            c->ok = false; return -1;
        }
        Val x; if (!emit_unary_math(c, h, A[1], OP_LOG_R, OP_LOG_C, &x)) return -1;
        Val b; if (!emit_unary_math(c, h, A[0], OP_LOG_R, OP_LOG_C, &b)) return -1;
        CompileType t = num_common(x.type, b.type); coerce(c, &x, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_DIV_C : OP_DIV_R, x, b, t);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Tanh") == 0 && na == 1) return emit_unary_math(c, h, A[0], OP_TANH_R, OP_TANH_C, out);

    /* Projections and rounding over an array.  Each of these heads is handled
     * below by a dedicated scalar branch that would never reach try_kernel, yet
     * every one of them has a registered ND unary kernel — so intercept the
     * array case here and map the whole buffer in one pass. */
    if (na == 1 && (!strcmp(h, "Abs") || !strcmp(h, "Sign") || !strcmp(h, "Floor")
                    || !strcmp(h, "Ceiling") || !strcmp(h, "Round") || !strcmp(h, "Re")
                    || !strcmp(h, "Im") || !strcmp(h, "Arg") || !strcmp(h, "Conjugate")
                    /* IntegerPart and UnitStep were missing from this list --
                     * the two narrowing heads whose scalar branch is dedicated
                     * AND whose kernel has no real or complex arm, so neither
                     * this interception nor try_kernel would take them and the
                     * array form had no lowering at all. */
                    || !strcmp(h, "IntegerPart") || !strcmp(h, "UnitStep"))
        && any_array_arg(c, A, na)) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (!CT_IS_ARRAY(a.type)) { c->ok = false; return -1; }
        return emit_arr_unary(c, h, a, out);
    }

    if (strcmp(h, "Abs") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (a.type == CT_COMPLEX) { *out = unop(c, OP_ABS_C, a, CT_REAL); return c->ok ? 1 : -1; }
        if (a.type == CT_INT)     { *out = unop(c, OP_ABS_I, a, CT_INT); return c->ok ? 1 : -1; }
        coerce(c, &a, CT_REAL); *out = unop(c, OP_ABS_R, a, CT_REAL); return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Sign") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (a.type == CT_INT)  { *out = unop(c, OP_SIGN_I, a, CT_INT); return c->ok ? 1 : -1; }
        /* CT_INT, not CT_REAL: `Sign[-2.5]` is the Integer -1 in the
         * interpreter, and a compiled path answering -1. would differ in HEAD
         * from the one it is supposed to be interchangeable with.  Same reason
         * UnitStep is lowered by hand. */
        if (a.type == CT_REAL) { *out = unop(c, OP_SIGN_R, a, CT_INT); return c->ok ? 1 : -1; }
        if (a.type == CT_COMPLEX) {
            /* z/|z|, and 0 at the origin — already in the shared kernel
             * registry, which this branch was shadowing by bailing first. */
            SymbolDef* d = symtab_lookup("Sign");
            const NDUnaryKernel* k = d ? (const NDUnaryKernel*)d->ndarray_unary_kernel : NULL;
            if (!k || !k->cplx) { c->ok = false; return -1; }
            *out = kern_unop(c, OP_KERN_CC, a, CT_COMPLEX, (const void*)k->cplx);
            return c->ok ? 1 : -1;
        }
        c->ok = false; return -1;
    }
    if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0
         || strcmp(h, "IntegerPart") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (a.type == CT_INT) { *out = a; return c->ok ? 1 : -1; }
        if (a.type != CT_REAL) coerce(c, &a, CT_REAL);
        /* IntegerPart is here rather than on its registered kernel because the
         * kernel returns a double and the interpreter returns an Integer —
         * `IntegerPart[2.5]` is `2`, not `2.`.  The kernel stays for the ARRAY
         * path, where a packed real buffer is the right answer. */
        uint16_t op = h[0] == 'F' ? OP_FLOOR_R : h[0] == 'C' ? OP_CEIL_R
                    : h[0] == 'R' ? OP_ROUND_R : OP_TRUNC_R;
        *out = unop(c, op, a, CT_INT); return c->ok ? 1 : -1;
    }
    /* FractionalPart of an integer is the Integer 0.  Its registered kernel
     * returns a double, so without this the head comes back as 0. — the same
     * trap IntegerPart is lowered by hand for. */
    if (strcmp(h, "FractionalPart") == 0 && na == 1) {
        CompileType at;
        if (infer_type(c, A[0], &at) && at == CT_INT) {
            Val a; if (!emit(c, A[0], &a)) return -1;
            free_if_tmp(c, a); memset(&imm, 0, sizeof imm); imm.i = 0;
            *out = emit_const(c, imm, CT_INT);
            return c->ok ? 1 : -1;
        }
    }
    if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 || strcmp(h, "Arg") == 0 || strcmp(h, "Conjugate") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (a.type == CT_INT && strcmp(h, "Conjugate") != 0) {
            /* Re[n] = n, Im[n] = the Integer 0, Arg[n] = the Integer 0 or Pi. */
            if (strcmp(h, "Re") == 0) { *out = a; return c->ok ? 1 : -1; }
            if (strcmp(h, "Im") == 0) {
                free_if_tmp(c, a); memset(&imm, 0, sizeof imm); imm.i = 0;
                *out = emit_const(c, imm, CT_INT); return c->ok ? 1 : -1;
            }
            *out = unop(c, OP_ARG_I, a, CT_INT); return c->ok ? 1 : -1;
        }
        if (a.type != CT_COMPLEX) {
            if (strcmp(h, "Im") == 0) { free_if_tmp(c, a); imm.r = 0.0; *out = emit_const(c, imm, CT_REAL); return c->ok ? 1 : -1; }
            /* Arg of a real is 0 or Pi, which the interpreter evaluates, so
             * widen and use the same opcode rather than bailing the whole body:
             * carg(x + 0i) gives exactly 0 for x > 0 and Pi for x < 0. */
            if (strcmp(h, "Arg") == 0) {
                coerce(c, &a, CT_COMPLEX);
                if (!c->ok) return -1;
                *out = unop(c, OP_ARG_C, a, CT_REAL);
                return c->ok ? 1 : -1;
            }
            *out = a; return c->ok ? 1 : -1;   /* Re/Conjugate of real = itself */
        }
        if (strcmp(h, "Re") == 0)   { *out = unop(c, OP_RE_C, a, CT_REAL); return c->ok ? 1 : -1; }
        if (strcmp(h, "Im") == 0)   { *out = unop(c, OP_IM_C, a, CT_REAL); return c->ok ? 1 : -1; }
        if (strcmp(h, "Arg") == 0)  { *out = unop(c, OP_ARG_C, a, CT_REAL); return c->ok ? 1 : -1; }
        *out = unop(c, OP_CONJ_C, a, CT_COMPLEX); return c->ok ? 1 : -1;
    }
    /* UnitStep / Clip / Rescale: lowered by hand rather than registered as
     * kernels, because for these the RESULT TYPE is the whole difficulty.
     *
     *   UnitStep returns an INTEGER (UnitStep[0.5] is 1, not 1.) — a real-valued
     *   kernel would answer with a different head from the interpreter.  A
     *   comparison already leaves 0/1 in the integer half of the slot, so the
     *   lowering is just the comparison, typed CT_INT.
     *
     *   Clip and Rescale take a LIST of bounds, which no kernel signature can
     *   express; a literal two-element List is destructured here instead.  Both
     *   are gated on the bounds being REAL: with exact bounds the interpreter
     *   returns an exact value where it clips (Clip[5, {1, 3}] is the Integer 3)
     *   and a Real where it does not, which no single compiled type can match. */
    /* HypergeometricPFQ[{a1..ap}, {b1..bq}, z] — the head that actually matters,
     * because the evaluator canonicalises Hypergeometric0F1, 1F1 and 2F1 into it
     * before anything downstream sees them.  Its LIST arguments are what no
     * kernel signature can express, so the two literal Lists are destructured
     * here and flattened into the consecutive block the n-ary opcode wants,
     * with p passed as the leading element so the kernel can find the split. */
    if (strcmp(h, "HypergeometricPFQ") == 0 && na == 3) {
        const Expr* la = A[0];
        const Expr* lb = A[1];
        if (la->type != EXPR_FUNCTION || la->data.function.head->type != EXPR_SYMBOL
            || strcmp(la->data.function.head->data.symbol.name, "List") != 0
            || lb->type != EXPR_FUNCTION || lb->data.function.head->type != EXPR_SYMBOL
            || strcmp(lb->data.function.head->data.symbol.name, "List") != 0)
            { c->ok = false; return -1; }
        size_t np = la->data.function.arg_count, nq = lb->data.function.arg_count;
        if (np + nq + 2 > 8) { c->ok = false; return -1; }   /* KERNN operand cap */

        const Expr* parts[8];
        size_t total = 0;
        for (size_t i = 0; i < np; i++) parts[total++] = la->data.function.args[i];
        for (size_t j = 0; j < nq; j++) parts[total++] = lb->data.function.args[j];
        parts[total++] = A[2];
        for (size_t i = 0; i < total; i++) {
            CompileType t;
            if (!infer_type(c, parts[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL)
                { c->ok = false; return -1; }
        }

        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) { c->ok = false; return -1; }
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx) { c->ok = false; return -1; }

        Slot z; memset(&z, 0, sizeof z);
        int base = alloc_temp(c);                     /* slot 0 holds p */
        for (size_t i = 0; i < total; i++) (void)alloc_temp(c);
        int after = c->temp_top;
        Slot pk; memset(&pk, 0, sizeof pk); pk.r = (double)np;
        ins(c, OP_CONST, (uint32_t)base, 0, 0, pk);
        for (size_t i = 0; i < total && c->ok; i++) {
            Val v;
            if (!emit(c, parts[i], &v)) return -1;
            coerce(c, &v, CT_REAL);
            if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)(base + 1 + (int)i), (uint32_t)v.reg, 0, z);
            c->temp_top = after;
        }
        if (!c->ok) return -1;
        c->temp_top = (base - c->nlocals) + 1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        ins_f(c, OP_KERNN, (uint16_t)(total + 1), (uint32_t)base, (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = CT_REAL;
        return c->ok ? 1 : -1;
    }

    if (strcmp(h, "UnitStep") == 0 && na >= 1) {
        Val acc; bool have = false;
        for (size_t i = 0; i < na; i++) {
            CompileType at;
            if (!infer_type(c, A[i], &at) || CT_IS_ARRAY(at) || at == CT_BOOL) { c->ok = false; return -1; }
            Val v; if (!emit(c, A[i], &v)) return -1;
            Val z;
            if (at == CT_INT) { Slot k; memset(&k, 0, sizeof k); z = emit_const(c, k, CT_INT); }
            else { coerce(c, &v, CT_REAL); Slot k; memset(&k, 0, sizeof k); k.r = 0.0; z = emit_const(c, k, CT_REAL); }
            if (!c->ok) return -1;
            Val ge = binop(c, (at == CT_INT) ? OP_GE_I : OP_GE_R, v, z, CT_INT);
            acc = have ? binop(c, OP_MUL_I, acc, ge, CT_INT) : ge;
            have = true;
        }
        if (!have) { c->ok = false; return -1; }
        *out = acc;
        return c->ok ? 1 : -1;
    }
    if ((strcmp(h, "Clip") == 0 || strcmp(h, "Rescale") == 0) && na == 2) {
        const Expr* bnd = A[1];
        if (bnd->type != EXPR_FUNCTION || bnd->data.function.head->type != EXPR_SYMBOL
            || strcmp(bnd->data.function.head->data.symbol.name, "List") != 0
            || bnd->data.function.arg_count != 2) { c->ok = false; return -1; }
        CompileType tx, tl, tu;
        if (!infer_type(c, A[0], &tx) || !infer_type(c, bnd->data.function.args[0], &tl)
            || !infer_type(c, bnd->data.function.args[1], &tu)) { c->ok = false; return -1; }
        /* Real bounds only — see the note above. */
        if (tl != CT_REAL || tu != CT_REAL || CT_IS_ARRAY(tx) || tx == CT_BOOL) { c->ok = false; return -1; }
        /* Rescale is just (x - lo)/(hi - lo), which is defined for a complex x
         * and which the interpreter evaluates (`Rescale[1. + I, {0., 2.}]` is
         * `0.5 + 0.5 I`).  Clip is NOT: the interpreter leaves `Clip[1. + I,
         * {0., 2.}]` unevaluated, because Min/Max need an order. */
        bool cx = (tx == CT_COMPLEX);
        if (cx && h[0] == 'C') { c->ok = false; return -1; }
        /* Temporaries are a STACK: binop pops its two operands and allocates the
         * destination, so both operands must be the top of it at that moment.
         * Everything below is therefore emitted in exactly the order it is
         * consumed — and `lo`, which Rescale needs twice, is emitted twice
         * rather than held across an intervening allocation.  (It is a literal
         * bound; the optimiser's value numbering folds the duplicate away.) */
        const Expr* elo = bnd->data.function.args[0];
        const Expr* ehi = bnd->data.function.args[1];
        Val x, lo, hi;
        if (!emit(c, A[0], &x)) return -1;
        coerce(c, &x, cx ? CT_COMPLEX : CT_REAL);
        if (!c->ok) return -1;

        if (h[0] == 'C') {                       /* Clip[x, {lo, hi}] = Min[Max[x, lo], hi] */
            if (!emit(c, elo, &lo)) return -1;
            Val mx = binop(c, OP_MAX_R, x, lo, CT_REAL);
            if (!emit(c, ehi, &hi)) return -1;
            *out = binop(c, OP_MIN_R, mx, hi, CT_REAL);
        } else {                                 /* Rescale[x, {lo, hi}] = (x - lo)/(hi - lo) */
            CompileType t = cx ? CT_COMPLEX : CT_REAL;
            if (!emit(c, elo, &lo)) return -1;
            if (cx) coerce(c, &lo, CT_COMPLEX);
            Val num = binop(c, cx ? OP_SUB_C : OP_SUB_R, x, lo, t);
            Val lo2, hi2;
            if (!emit(c, ehi, &hi2)) return -1;
            if (cx) coerce(c, &hi2, CT_COMPLEX);
            if (!emit(c, elo, &lo2)) return -1;
            if (cx) coerce(c, &lo2, CT_COMPLEX);
            Val den = binop(c, cx ? OP_SUB_C : OP_SUB_R, hi2, lo2, t);
            *out = binop(c, cx ? OP_DIV_C : OP_DIV_R, num, den, t);
        }
        return c->ok ? 1 : -1;
    }

    /* The scalar pairwise fold.  Array operands are refused for the reason
     * given beside the matching inference block: an array reaching the fold
     * makes Max[v] lower to the identity, which is a wrong answer rather than
     * a slow one.  na == 1 over an array falls through to the delegated
     * reduction further down. */
    if ((strcmp(h, "Max") == 0 || strcmp(h, "Min") == 0) && na >= 1
        && !any_array_arg(c, A, na)) {
        bool mx = h[1] == 'a';
        Val acc; if (!emit(c, A[0], &acc)) return -1;
        for (size_t i = 1; i < na; i++) {
            Val b; if (!emit(c, A[i], &b)) return -1;
            CompileType t = num_common(acc.type, b.type);
            if ((int)t < 0 || t == CT_COMPLEX) { c->ok = false; return -1; }
            coerce(c, &acc, t); coerce(c, &b, t);
            acc = binop(c, mx ? (t == CT_INT ? OP_MAX_I : OP_MAX_R) : (t == CT_INT ? OP_MIN_I : OP_MIN_R), acc, b, t);
        }
        *out = acc; return c->ok ? 1 : -1;
    }
    if (strcmp(h, "ArcTan") == 0) {
        if (na == 1) return emit_unary_math(c, h, A[0], OP_ATAN_R, OP_ATAN_C, out);
        if (na == 2) {
            Val kv2;
            if (any_array_arg(c, A, na)) {            /* ND two-arg ArcTan kernel */
                if (try_kernel(c, h, A, na, &kv2)) { *out = kv2; return c->ok ? 1 : -1; }
                c->ok = false; return -1;
            }
            Val x, y; if (!emit(c, A[0], &x) || !emit(c, A[1], &y)) return -1;
            coerce(c, &x, CT_REAL); coerce(c, &y, CT_REAL);
            *out = binop(c, OP_ATAN2_R, x, y, CT_REAL); return c->ok ? 1 : -1;   /* atan2(y,x): b=y top */
        }
        c->ok = false; return -1;
    }

    /* ConstantArray[v, dims] (M3c): a zeroed buffer, then a fill loop only when
     * the value is not already zero.  Emitting the fill as ordinary bytecode
     * (rather than as a memset variant per element type) means a non-constant
     * fill costs nothing extra and the optimiser hoists it like any other loop. */
    return 0;
}
