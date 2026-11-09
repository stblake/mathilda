/* Mathilda — Compile[]: elementary/special-function + projection heads (Log, Sin..Tanh, Abs/Sign/Floor/Ceiling/Round/Re/Im/Arg/Conjugate/IntegerPart/FractionalPart/UnitStep, Clip/Rescale, HypergeometricPFQ, Max/Min, ArcTan).
 *
 * One of emit_node's per-category dispatch functions: 1 = handled (ok iff
 * c->ok), -1 = handled but not lowerable, 0 = not one of these heads. */
#include "compile_emit.h"
#include "../arithmetic.h"
#include "../symtab.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../special_functions/sf_machine_meijerg.h" /* SF_MG_MAXP, sf_mg_shape_code */
#include "../special_functions/dirichlet.h"        /* dch_labels_real */
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Rank-1 array form of a hand-lowered elementwise head (Chop/Clip): emit the
 * synthesised Map[Function[H[#, rest...]], arr].  The Map lowering applies the
 * scalar arm to each element, so the array answer is the scalar one elementwise,
 * bit for bit, and every argument form rides the scalar arm.  Returns 1/-1. */
static int emit_map_over_array(Ctx* c, const char* h, const Expr* e,
                               Expr* const* A, size_t na, Val* out) {
    Expr* m = build_map_scalar_over_array(h, A, na);
    if (!m) { c->ok = false; return -1; }
    bool r = emit(c, m, out);
    if (!r && c->bail_node && expr_subtree_of(m, c->bail_node)) c->bail_node = (Expr*)e;
    expr_free(m);
    return r ? 1 : -1;
}

int emit_mathfn(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    Slot imm;
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
                     * array form had no lowering at all.  UnitBox is a third of
                     * exactly the same shape (narrowing-only kernel, dedicated
                     * scalar branch below). */
                    || !strcmp(h, "IntegerPart") || !strcmp(h, "UnitStep")
                    || !strcmp(h, "UnitBox") || !strcmp(h, "Unitize"))
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
            || la->data.function.head->data.symbol.name != SYM_List
            || lb->type != EXPR_FUNCTION || lb->data.function.head->type != EXPR_SYMBOL
            || lb->data.function.head->data.symbol.name != SYM_List)
            { c->ok = false; return -1; }
        size_t np = la->data.function.arg_count, nq = lb->data.function.arg_count;
        if (np + nq + 2 > KERNN_MAX_OPERANDS) { c->ok = false; return -1; }   /* KERNN operand cap */

        const Expr* parts[KERNN_MAX_OPERANDS];
        size_t total = 0;
        for (size_t i = 0; i < np; i++) parts[total++] = la->data.function.args[i];
        for (size_t j = 0; j < nq; j++) parts[total++] = lb->data.function.args[j];
        parts[total++] = A[2];
        bool cx = false;    /* any complex operand: the complex arm, OP_KERNN_C */
        for (size_t i = 0; i < total; i++) {
            CompileType t;
            if (!infer_type(c, parts[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL)
                { c->ok = false; return -1; }
            cx |= (t == CT_COMPLEX);
        }
        CompileType et = cx ? CT_COMPLEX : CT_REAL;

        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) { c->ok = false; return -1; }
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx) { c->ok = false; return -1; }

        Slot z; memset(&z, 0, sizeof z);
        int base = alloc_temp(c);                     /* slot 0 holds p */
        for (size_t i = 0; i < total; i++) (void)alloc_temp(c);
        int after = c->temp_top;
        Slot pk; memset(&pk, 0, sizeof pk);
        if (cx) pk.z = (double)np; else pk.r = (double)np;
        ins(c, OP_CONST, (uint32_t)base, 0, 0, pk);
        for (size_t i = 0; i < total && c->ok; i++) {
            Val v;
            if (!emit(c, parts[i], &v)) return -1;
            coerce(c, &v, et);
            if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)(base + 1 + (int)i), (uint32_t)v.reg, 0, z);
            c->temp_top = after;
        }
        if (!c->ok) return -1;
        c->temp_top = (base - c->nlocals) + 1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        ins_f(c, cx ? OP_KERNN_C : OP_KERNN, (uint16_t)(total + 1), (uint32_t)base,
              (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = et;
        return c->ok ? 1 : -1;
    }

    /* DirichletL[k, j, s] -- the run [k, j, s] for the n-ary kernel
     * ndk_DirichletL_n (sf_dirichlet_l_c: the scalar builtin's helper). Inside
     * a strip-mined loop s is a TILE: the run stops before s and VKERNN_* reads
     * s per element, the labels broadcast. Typed as compile_infer.c decides:
     * Real only for literal labels of a real character at real s. */
    /* WeierstrassP / PPrime / Zeta / Sigma [u, {g2, g3}] -- the literal List of
     * invariants destructured into the run [g2, g3, u] the n-ary kernel
     * (ndk_Weierstrass*_n, sf_machine_weierstrass: the scalar builtins' helper)
     * reads. Inside a strip-mined loop u is a TILE: the run stops before u and
     * VKERNN_* reads it per element, the invariants broadcast (the kernel's
     * plan memo then builds the lattice once). */
    if ((strcmp(h, "WeierstrassP") == 0 || strcmp(h, "WeierstrassPPrime") == 0 ||
         strcmp(h, "WeierstrassZeta") == 0 || strcmp(h, "WeierstrassSigma") == 0) && na == 2) {
        const Expr* L = A[1];
        if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL ||
            L->data.function.head->data.symbol.name != SYM_List || L->data.function.arg_count != 2)
            { c->ok = false; return -1; }
        CompileType t2, t3, tu;
        if (!infer_type(c, L->data.function.args[0], &t2) ||
            !infer_type(c, L->data.function.args[1], &t3) || !infer_type(c, A[0], &tu))
            { c->ok = false; return -1; }
        CompileType eu = CT_IS_ARRAY(tu) ? CT_ELEM(tu) : tu;
        if (CT_IS_ARRAY(t2) || CT_IS_ARRAY(t3) || t2 == CT_BOOL || t3 == CT_BOOL || eu == CT_BOOL)
            { c->ok = false; return -1; }
        bool cx = eu == CT_COMPLEX || t2 == CT_COMPLEX || t3 == CT_COMPLEX;
        CompileType et = cx ? CT_COMPLEX : CT_REAL;
        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) { c->ok = false; return -1; }
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx) { c->ok = false; return -1; }
        Slot z0; memset(&z0, 0, sizeof z0);
        int base = alloc_temp(c);
        (void)alloc_temp(c); (void)alloc_temp(c);
        int after = c->temp_top;
        for (int i = 0; i < 2 && c->ok; i++) {
            Val v;
            if (!emit(c, L->data.function.args[i], &v)) return -1;
            if (val_is_tile(v)) { c->ok = false; return -1; }   /* invariants broadcast only */
            coerce(c, &v, et);
            if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)(base + i), (uint32_t)v.reg, 0, z0);
            c->temp_top = after;
        }
        if (!c->ok) return -1;
        Val vu;
        if (!emit(c, A[0], &vu)) return -1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        if (val_is_tile(vu)) {
            coerce(c, &vu, et);
            if (!c->ok) return -1;
            int dst = alloc_tile(c);
            ins_f(c, cx ? OP_VKERNN_C : OP_VKERNN_R, 2, (uint32_t)dst,
                  (uint32_t)base, (uint32_t)vu.reg, kp);
            c->temp_top = (base - c->nlocals);
            out->reg = dst; out->tmp = true; out->type = et;
            return c->ok ? 1 : -1;
        }
        coerce(c, &vu, et);
        if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)(base + 2), (uint32_t)vu.reg, 0, z0);
        c->temp_top = (base - c->nlocals) + 1;
        ins_f(c, cx ? OP_KERNN_C : OP_KERNN, 3, (uint32_t)base, (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = et;
        return c->ok ? 1 : -1;
    }

    if (strcmp(h, "DirichletL") == 0 && na == 3) {
        CompileType tk, tj, ts;
        if (!infer_type(c, A[0], &tk) || !infer_type(c, A[1], &tj) || !infer_type(c, A[2], &ts)
            || CT_IS_ARRAY(tk) || CT_IS_ARRAY(tj) || CT_IS_ARRAY(ts)
            || tk == CT_BOOL || tj == CT_BOOL || ts == CT_BOOL
            || tk == CT_COMPLEX || tj == CT_COMPLEX) { c->ok = false; return -1; }
        bool cx = (ts == CT_COMPLEX) || !dch_labels_real(A[0], A[1]);
        CompileType et = cx ? CT_COMPLEX : CT_REAL;
        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) { c->ok = false; return -1; }
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx) { c->ok = false; return -1; }
        Slot z0; memset(&z0, 0, sizeof z0);
        int base = alloc_temp(c);
        (void)alloc_temp(c); (void)alloc_temp(c);
        int after = c->temp_top;
        for (int i = 0; i < 2 && c->ok; i++) {
            Val v;
            if (!emit(c, A[i], &v)) return -1;
            if (val_is_tile(v)) { c->ok = false; return -1; }   /* labels broadcast only */
            coerce(c, &v, et);
            if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)(base + i), (uint32_t)v.reg, 0, z0);
            c->temp_top = after;
        }
        if (!c->ok) return -1;
        Val vs;
        if (!emit(c, A[2], &vs)) return -1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        if (val_is_tile(vs)) {
            coerce(c, &vs, et);
            if (!c->ok) return -1;
            int dst = alloc_tile(c);
            ins_f(c, cx ? OP_VKERNN_C : OP_VKERNN_R, 2, (uint32_t)dst,
                  (uint32_t)base, (uint32_t)vs.reg, kp);
            c->temp_top = (base - c->nlocals);
            out->reg = dst; out->tmp = true; out->type = et;
            return c->ok ? 1 : -1;
        }
        coerce(c, &vs, et);
        if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)(base + 2), (uint32_t)vs.reg, 0, z0);
        c->temp_top = (base - c->nlocals) + 1;
        ins_f(c, cx ? OP_KERNN_C : OP_KERNN, 3, (uint32_t)base, (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = et;
        return c->ok ? 1 : -1;
    }

    /* MeijerG[{{a..},{a..}}, {{b..},{b..}}, z(, r)] -- the four parameter
     * Lists are destructured and flattened, behind a shape code and r, into the
     * run the n-ary kernel reads (sf_mg_eval_flat's layout):
     *     [n + 64 p + 4096 m, r, a_1..a_p, b_1..b_q, z].
     * Inside a strip-mined loop z is a TILE (MeijerG threads over z only, so it
     * is the one operand that may be array-valued): the run then stops before z
     * and VKERNN_* reads z from the tile, the parameters broadcast. */
    if (strcmp(h, "MeijerG") == 0 && (na == 3 || na == 4)) {
        const Expr* g[4];
        for (int s = 0; s < 2; s++) {
            const Expr* L = A[s];
            if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
                || L->data.function.head->data.symbol.name != SYM_List
                || L->data.function.arg_count != 2) { c->ok = false; return -1; }
            for (int t = 0; t < 2; t++) {
                const Expr* I2 = L->data.function.args[t];
                if (I2->type != EXPR_FUNCTION || I2->data.function.head->type != EXPR_SYMBOL
                    || I2->data.function.head->data.symbol.name != SYM_List)
                    { c->ok = false; return -1; }
                g[2 * s + t] = I2;
            }
        }
        size_t n = g[0]->data.function.arg_count, p = n + g[1]->data.function.arg_count;
        size_t m = g[2]->data.function.arg_count, q = m + g[3]->data.function.arg_count;
        if (p > SF_MG_MAXP || q > SF_MG_MAXP || p + q + 3 > KERNN_MAX_OPERANDS)
            { c->ok = false; return -1; }
        const Expr* parts[KERNN_MAX_OPERANDS];
        size_t np = 0;
        for (int s = 0; s < 4; s++)
            for (size_t i = 0; i < g[s]->data.function.arg_count; i++)
                parts[np++] = g[s]->data.function.args[i];
        bool cx = false;
        for (size_t i = 0; i < np; i++) {
            CompileType t;
            if (!infer_type(c, parts[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL)
                { c->ok = false; return -1; }
            cx |= (t == CT_COMPLEX);
        }
        if (na == 4) {
            CompileType tr;
            if (!infer_type(c, A[3], &tr) || CT_IS_ARRAY(tr) || tr == CT_BOOL || tr == CT_COMPLEX)
                { c->ok = false; return -1; }
        }
        CompileType tz;
        if (!infer_type(c, A[2], &tz) || CT_IS_ARRAY(tz) || tz == CT_BOOL) { c->ok = false; return -1; }
        cx |= (tz == CT_COMPLEX);
        CompileType et = cx ? CT_COMPLEX : CT_REAL;

        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) { c->ok = false; return -1; }
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx) { c->ok = false; return -1; }

        /* The scalar run: shape code, r, then the parameters; z appended below
         * when it is a scalar. */
        size_t nrun = np + 2;
        Slot z0; memset(&z0, 0, sizeof z0);
        int base = alloc_temp(c);
        for (size_t i = 1; i < nrun + 1; i++) (void)alloc_temp(c);   /* + z's slot */
        int after = c->temp_top;
        Slot ks; memset(&ks, 0, sizeof ks);
        double code = sf_mg_shape_code(n, p, m);
        if (cx) ks.z = code; else ks.r = code;
        ins(c, OP_CONST, (uint32_t)base, 0, 0, ks);
        if (na == 4) {
            Val v;
            if (!emit(c, A[3], &v)) return -1;
            if (val_is_tile(v)) { c->ok = false; return -1; }
            coerce(c, &v, et);
            if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)(base + 1), (uint32_t)v.reg, 0, z0);
            c->temp_top = after;
        } else {
            Slot one; memset(&one, 0, sizeof one);
            if (cx) one.z = 1.0; else one.r = 1.0;
            ins(c, OP_CONST, (uint32_t)(base + 1), 0, 0, one);
        }
        for (size_t i = 0; i < np && c->ok; i++) {
            Val v;
            if (!emit(c, parts[i], &v)) return -1;
            if (val_is_tile(v)) { c->ok = false; return -1; }   /* parameters broadcast only */
            coerce(c, &v, et);
            if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)(base + 2 + (int)i), (uint32_t)v.reg, 0, z0);
            c->temp_top = after;
        }
        if (!c->ok) return -1;
        Val vz;
        if (!emit(c, A[2], &vz)) return -1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        if (val_is_tile(vz)) {
            coerce(c, &vz, et);
            if (!c->ok) return -1;
            int dst = alloc_tile(c);
            ins_f(c, cx ? OP_VKERNN_C : OP_VKERNN_R, (uint16_t)nrun, (uint32_t)dst,
                  (uint32_t)base, (uint32_t)vz.reg, kp);
            c->temp_top = (base - c->nlocals);
            out->reg = dst; out->tmp = true; out->type = et;
            return c->ok ? 1 : -1;
        }
        coerce(c, &vz, et);
        if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)(base + (int)nrun), (uint32_t)vz.reg, 0, z0);
        c->temp_top = (base - c->nlocals) + 1;
        ins_f(c, cx ? OP_KERNN_C : OP_KERNN, (uint16_t)(nrun + 1), (uint32_t)base,
              (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = et;
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
    if (strcmp(h, "Unitize") == 0 && (na == 1 || na == 2)) {
        /* Unitize[x] = (x != 0) and Unitize[x, dx] = (|x| >= dx), each an exact
         * Integer 0/1 (CT_INT) -- the comparison result IS the answer, the same
         * arithmetic as the shared machine helpers in piecewise.c that the
         * interpreter's scalar and buffer paths run. A one-argument array was
         * intercepted above (the narrowing kernel); Unitize[v, dx] runs the
         * registered binary narrowing kernel over the buffer, which declines
         * (handing the call back) on a threshold that is not positive. */
        CompileType tx;
        if (!infer_type(c, A[0], &tx) || tx == CT_BOOL) { c->ok = false; return -1; }
        if (CT_IS_ARRAY(tx)) {
            Val kv;
            if (na == 2 && try_kernel(c, h, A, na, &kv)) { *out = kv; return c->ok ? 1 : -1; }
            c->ok = false; return -1;
        }
        Slot kz; memset(&kz, 0, sizeof kz);          /* 0 / 0. (all-zero bits) */
        if (na == 1) {
            Val v; if (!emit(c, A[0], &v)) return -1;
            if (v.type == CT_COMPLEX) v = unop(c, OP_ABS_C, v, CT_REAL);
            if (v.type == CT_INT) {
                *out = binop(c, OP_NE_I, v, emit_const(c, kz, CT_INT), CT_INT);
            } else {
                coerce(c, &v, CT_REAL);
                *out = binop(c, OP_NE_R, v, emit_const(c, kz, CT_REAL), CT_INT);
            }
            return c->ok ? 1 : -1;
        }
        CompileType td;
        if (!infer_type(c, A[1], &td) || CT_IS_ARRAY(td) || td == CT_BOOL
            || td == CT_COMPLEX) { c->ok = false; return -1; }
        bool both_int = (tx == CT_INT && td == CT_INT);
        /* Run-time guard: a threshold that is not positive has no compiled
         * answer (the interpreter issues Unitize::post and leaves the call
         * unevaluated), so OP_FAIL hands the call back. */
        {
            Val g; if (!emit(c, A[1], &g)) return -1;
            if (!both_int) coerce(c, &g, CT_REAL);
            if (!c->ok) return -1;
            Slot z; memset(&z, 0, sizeof z);
            int rz = alloc_temp(c), rc = alloc_temp(c);
            ins(c, OP_CONST, (uint32_t)rz, 0, 0, kz);
            ins(c, both_int ? OP_LE_I : OP_LE_R, (uint32_t)rc, (uint32_t)g.reg, (uint32_t)rz, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);   /* dx > 0 -> skip */
            ins(c, OP_FAIL, 0, 0, 0, z);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
            c->temp_top -= 2;
            free_if_tmp(c, g);
        }
        if (both_int) {
            /* (x >= dx) + (x <= -dx): no Abs, so INT64_MIN cannot overflow;
             * dx > 0 makes the two disjoint. */
            Val x1; if (!emit(c, A[0], &x1)) return -1;
            Val d1; if (!emit(c, A[1], &d1)) return -1;
            Val hi = binop(c, OP_GE_I, x1, d1, CT_INT);
            Val x2; if (!emit(c, A[0], &x2)) return -1;
            Val d2; if (!emit(c, A[1], &d2)) return -1;
            Val nd = unop(c, OP_NEG_I, d2, CT_INT);
            Val lo = binop(c, OP_LE_I, x2, nd, CT_INT);
            *out = binop(c, OP_ADD_I, hi, lo, CT_INT);
            return c->ok ? 1 : -1;
        }
        Val x; if (!emit(c, A[0], &x)) return -1;
        Val ax;
        if (x.type == CT_COMPLEX) ax = unop(c, OP_ABS_C, x, CT_REAL);
        else { coerce(c, &x, CT_REAL); ax = unop(c, OP_ABS_R, x, CT_REAL); }
        Val d; if (!emit(c, A[1], &d)) return -1;
        coerce(c, &d, CT_REAL);
        if (!c->ok) return -1;
        *out = binop(c, OP_GE_R, ax, d, CT_INT);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "UnitBox") == 0 && na == 1) {
        /* UnitBox[x] = (x >= -1/2) * (x <= 1/2), the closed box, answered as an
         * exact Integer 0/1 (CT_INT) exactly like the UnitStep idiom above. The
         * argument is emitted twice -- once per comparison -- because binop
         * consumes its operands; the optimiser folds the duplicate arg load. An
         * integer x is in the box iff x == 0, i.e. (x >= 0) * (x <= 0). */
        CompileType at;
        if (!infer_type(c, A[0], &at) || CT_IS_ARRAY(at) || at == CT_BOOL
            || at == CT_COMPLEX) { c->ok = false; return -1; }
        Val lo, hi;
        if (at == CT_INT) {
            Val va; if (!emit(c, A[0], &va)) return -1;
            Slot kz; memset(&kz, 0, sizeof kz);
            lo = binop(c, OP_GE_I, va, emit_const(c, kz, CT_INT), CT_INT);
            Val vb; if (!emit(c, A[0], &vb)) return -1;
            hi = binop(c, OP_LE_I, vb, emit_const(c, kz, CT_INT), CT_INT);
        } else {
            Val va; if (!emit(c, A[0], &va)) return -1;
            coerce(c, &va, CT_REAL);
            Slot klo; memset(&klo, 0, sizeof klo); klo.r = -0.5;
            lo = binop(c, OP_GE_R, va, emit_const(c, klo, CT_REAL), CT_INT);
            Val vb; if (!emit(c, A[0], &vb)) return -1;
            coerce(c, &vb, CT_REAL);
            Slot khi; memset(&khi, 0, sizeof khi); khi.r = 0.5;
            hi = binop(c, OP_LE_R, vb, emit_const(c, khi, CT_REAL), CT_INT);
        }
        if (!c->ok) return -1;
        *out = binop(c, OP_MUL_I, lo, hi, CT_INT);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Boole") == 0 && na == 1) {
        CompileType tb;
        if (!infer_type(c, A[0], &tb)) { c->ok = false; return -1; }
        /* A rank-1 bool array applies the scalar arm below to each element,
         * into an int64 buffer -- the interpreter's Boole buffer path. */
        if (CT_IS_ARRAY(tb)) return emit_elementwise_over_array(c, h, e, A[0], out);
        /* Boole of a non-Boolean is unevaluated in the interpreter: bail. */
        if (tb != CT_BOOL) { c->ok = false; return -1; }
        /* A truth value already lives in the integer slot as 0/1 (every
         * comparison and logic op writes RD.i), so Boole is a retyping MOVE --
         * the exact Integer the interpreter answers, with no arithmetic. */
        Val b; if (!emit(c, A[0], &b)) return -1;
        *out = unop(c, OP_MOVE, b, CT_INT);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Chop") == 0 && (na == 1 || na == 2)) {
        CompileType tx;
        if (!infer_type(c, A[0], &tx)) { c->ok = false; return -1; }
        /* A rank-1 array maps the scalar lowering over the buffer. */
        if (CT_IS_ARRAY(tx)) {
            if (CT_RANK(tx) != 1) { c->ok = false; return -1; }
            return emit_map_over_array(c, h, e, A, na, out);
        }
        if (tx == CT_BOOL || tx == CT_COMPLEX) { c->ok = false; return -1; }
        /* The tolerance (default 1e-10) is any real-valued scalar operand.  It
         * used to be accepted only as a machine Real/Integer literal, so an
         * exact one -- `Chop[x, 10^-6]`, `Chop[x, 1/1000]` -- sent the whole
         * body to the interpreter even though the general operand path already
         * lowers exact numeric constants into a real slot (literal() rounds a
         * Rational correctly; the Power lowering takes 10^-6) and the
         * optimiser folds a constant one.  Typed here, emitted below. */
        CompileType td = CT_REAL;
        if (na == 2) {
            if (!infer_type(c, A[1], &td)) { c->ok = false; return -1; }
            if (CT_IS_ARRAY(td) || td == CT_BOOL || td == CT_COMPLEX) { c->ok = false; return -1; }
        }
        /* A machine integer is exact and never chops -> identity. */
        if (tx == CT_INT) return emit(c, A[0], out) ? (c->ok ? 1 : -1) : -1;
        /* Chop[x] = x * (Abs[x] >= delta): keep x when |x|>=delta, else 0. -- a
         * branchless mask-multiply (the UnitStep idiom above). The chopped value
         * is a machine 0. (Real), the compiled counterpart of the interpreter's
         * exact Integer 0; a CompiledFunction returns machine types. `x` is
         * emitted twice (for Abs and for the product), which the optimiser folds
         * to one register for an argument and value-numbers for a computed x. */
        Val xa; if (!emit(c, A[0], &xa)) return -1;
        Val ax = unop(c, OP_ABS_R, xa, CT_REAL);
        Val dv;
        if (na == 2) {
            /* |delta|: the interpreter chops against the tolerance's magnitude. */
            Val d; if (!emit(c, A[1], &d)) return -1;
            coerce(c, &d, CT_REAL);
            if (!c->ok) return -1;
            dv = unop(c, OP_ABS_R, d, CT_REAL);
        } else {
            Slot dz; memset(&dz, 0, sizeof dz); dz.r = 1.0e-10;
            dv = emit_const(c, dz, CT_REAL);
        }
        Val mask = binop(c, OP_GE_R, ax, dv, CT_INT);
        coerce(c, &mask, CT_REAL);
        Val xb; if (!emit(c, A[0], &xb)) return -1;
        Val prod = binop(c, OP_MUL_R, mask, xb, CT_REAL);
        /* + 0.0: a chopped NEGATIVE x gives 0. * x == -0., where the
         * interpreter's chop (and its packed / NDArray buffer paths) answer
         * +0.  IEEE round-to-nearest makes -0. + 0. == +0., and x + 0. == x
         * for every unchopped x (only x == -0. itself changes, and it is
         * always chopped).  The optimiser has no x + 0 -> x identity fold, and
         * a contracted fma(mask, x, 0.) gives the same +0. */
        Slot pz; memset(&pz, 0, sizeof pz); pz.r = 0.0;
        Val zv = emit_const(c, pz, CT_REAL);
        *out = binop(c, OP_ADD_R, prod, zv, CT_REAL);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Clip") == 0 && na == 1) {
        /* Clip[x] default bounds {-1, 1} = Min[Max[x, -1.], 1.] -- the na==2 arm
         * below with the bounds synthesised. Real result (the machine 1. where
         * the interpreter's exact-bound scalar gives the Integer 1). */
        CompileType tx;
        if (!infer_type(c, A[0], &tx)) { c->ok = false; return -1; }
        if (CT_IS_ARRAY(tx)) {
            if (CT_RANK(tx) != 1) { c->ok = false; return -1; }
            return emit_map_over_array(c, h, e, A, na, out);
        }
        if (tx == CT_BOOL || tx == CT_COMPLEX) { c->ok = false; return -1; }
        Val x; if (!emit(c, A[0], &x)) return -1;
        coerce(c, &x, CT_REAL);
        if (!c->ok) return -1;
        Slot lo; memset(&lo, 0, sizeof lo); lo.r = -1.0;
        Val vlo = emit_const(c, lo, CT_REAL);
        Val mx = binop(c, OP_MAX_R, x, vlo, CT_REAL);
        Slot hi; memset(&hi, 0, sizeof hi); hi.r = 1.0;
        Val vhi = emit_const(c, hi, CT_REAL);
        *out = binop(c, OP_MIN_R, mx, vhi, CT_REAL);
        return c->ok ? 1 : -1;
    }
    if ((strcmp(h, "Clip") == 0 || strcmp(h, "Rescale") == 0) && na == 2) {
        /* A rank-1 array first arg maps the scalar lowering over the buffer. */
        { CompileType ta0;
          if (infer_type(c, A[0], &ta0) && CT_IS_ARRAY(ta0)) {
              if (CT_RANK(ta0) != 1) { c->ok = false; return -1; }
              return emit_map_over_array(c, h, e, A, na, out);
          } }
        const Expr* bnd = A[1];
        if (bnd->type != EXPR_FUNCTION || bnd->data.function.head->type != EXPR_SYMBOL
            || bnd->data.function.head->data.symbol.name != SYM_List
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
