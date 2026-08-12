/* Mathilda — Compile[]: array-construction + reduction heads (ConstantArray, Table, Range, Part, Total/Length, Mean/Median/...).
 *
 * One of emit_node's per-category dispatch functions: returns 1 if `h` is
 * one of its heads (result in *out, valid iff c->ok), -1 if it is one but
 * could not be lowered, 0 if `h` belongs to another category. */
#include "compile_emit.h"
#include "../arithmetic.h"
#include "../symtab.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include <string.h>
#include <stdlib.h>

int emit_array(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e;
    if (strcmp(h, "ConstantArray") == 0 && na == 2) {
        int rank; CompileType elem;
        if (!const_array_shape(c, (const Expr* const*)A, &rank, &elem)) { c->ok = false; return -1; }
        const Expr* d = A[1];
        const Expr* const* dexpr = &d;                 /* the bare-dimension case */
        if (rank > 1 || (d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL
                         && d->data.function.head->data.symbol.name == SYM_List))
            dexpr = (const Expr* const*)d->data.function.args;

        int base = -1;
        for (int i = 0; i < rank; i++) {               /* dims into consecutive regs */
            int r = alloc_temp(c);
            if (base < 0) base = r;
            Val dv;
            if (!emit(c, dexpr[i], &dv)) return -1;
            if (dv.type != CT_INT) { c->ok = false; return -1; }
            Slot z; memset(&z, 0, sizeof z);
            ins(c, OP_MOVE, (uint32_t)r, (uint32_t)dv.reg, 0, z);
            free_if_tmp(c, dv);
        }
        Slot el; memset(&el, 0, sizeof el); el.i = (long long)elem;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, (uint16_t)rank, (uint32_t)rout, (uint32_t)base, (uint32_t)base, el);

        Slot fz; memset(&fz, 0, sizeof fz);
        bool zero_fill = (A[0]->type == EXPR_REAL && A[0]->data.real == 0.0
                          && !signbit(A[0]->data.real));
        if (!zero_fill) {
            /* for k = 0 .. size-1: out[k] = v   (v hoisted by LICM if invariant) */
            int rn = alloc_temp(c), rk = alloc_temp(c);
            ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)rout, 0, fz);
            Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
            ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);
            int rc = alloc_temp(c);
            ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rk, (uint32_t)rn, fz);
            size_t jz = c->n;
            ins(c, OP_JZ, 0, (uint32_t)rc, 0, fz);
            c->temp_top--;                              /* rc dead after the guard */
            size_t body = c->n;
            Val fv;
            if (!emit(c, A[0], &fv)) return -1;
            coerce(c, &fv, elem);
            ins(c,a_store_op(elem),
                (uint32_t)rout, (uint32_t)rk, (uint32_t)fv.reg, fz);
            free_if_tmp(c, fv);
            Slot one; memset(&one, 0, sizeof one); one.i = 1;
            ins(c, OP_LOOP, (uint32_t)rk, (uint32_t)rn, (uint32_t)body, one);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
            c->temp_top -= 2;                           /* rn, rk */
        }
        c->temp_top -= rank;                            /* the dimension registers */
        /* A construction site: the interpreter's ConstantArray returns a List
         * whatever the arguments were, so the result kind must not follow them. */
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(elem, rank); out->built = true;
        return c->ok ? 1 : -1;
    }

    /* Table[body, spec1, ..., speck]: build a rank-k machine array.
     *
     * INTEGER iterators only, and that is not a convenience restriction.  The
     * interpreter walks a real iterator by repeated addition against a 1e-14
     * termination slack (src/list/table.c:90), so a closed-form `lo + k di`
     * differs from it in the last bits and, near the endpoint, in the element
     * COUNT.  Integer bounds reproduce that walk exactly.
     *
     * An INTEGER body is compiled, not refused: `Table[i, {i, 1, n}]` holds
     * exact Integers in the interpreter and a packed NDT_INT64 buffer holds
     * exactly those.  (Before the integer dtype existed this had to bail, since
     * a float64 buffer would have answered differently rather than faster.)
     *
     * One flat store index across k nested loops, innermost varying fastest —
     * which is row-major, and which is how the interpreter's nested rewrite
     * orders the elements too (see the multi-iterator Do). */
    if (strcmp(h, "Table") == 0 && na >= 2 && (int)na - 1 <= CT_MAX_RANK) {
        int rank = (int)na - 1;
        LoopSpec sp[CT_MAX_RANK];
        for (int j = 0; j < rank; j++)
            if (!loop_spec_parse(A[j + 1], &sp[j]) || !loop_spec_int_bounds(c, &sp[j]))
                { c->ok = false; return -1; }
        if (c->nscope + rank > CTX_MAX_SCOPE) { c->ok = false; return -1; }

        int pushed = 0;
        for (int j = 0; j < rank; j++)
            if (sp[j].var) {
                c->scope[c->nscope].name = sp[j].var; c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = CT_INT; c->scope[c->nscope].built = false;
                c->nscope++; pushed++;
            }
        CompileType et; bool okT = infer_type(c, A[0], &et);
        c->nscope -= pushed;
        if (!okT || !ct_is_elem(et)) { c->ok = false; return -1; }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;

        /* The dimension registers must be ONE contiguous run: A_NEW reads rank
         * of them starting at `base`. */
        int rn[CT_MAX_RANK], rlo[CT_MAX_RANK], base = -1;
        for (int j = 0; j < rank; j++) { int r = alloc_temp(c); if (base < 0) base = r; rn[j] = r; }
        for (int j = 0; j < rank; j++) rlo[j] = alloc_temp(c);
        for (int j = 0; j < rank; j++) {
            if (!emit_iter_count(c, &sp[j], rlo[j], rn[j])) return -1;
            emit_max_guard(c, rn[j], 1000000);
        }

        Slot el; memset(&el, 0, sizeof el); el.i = (long long)et;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, (uint16_t)rank, (uint32_t)rout, (uint32_t)base, (uint32_t)base, el);

        int rflat = alloc_temp(c);
        ins(c, OP_CONST, (uint32_t)rflat, 0, 0, k0);

        int riv[CT_MAX_RANK], rk[CT_MAX_RANK];
        size_t Lp[CT_MAX_RANK], jz[CT_MAX_RANK];
        int scope_entry = c->nscope;
        for (int j = 0; j < rank; j++) {
            riv[j] = alloc_temp(c); rk[j] = alloc_temp(c);
            ins(c, OP_MOVE, (uint32_t)riv[j], (uint32_t)rlo[j], 0, z);
            ins(c, OP_CONST, (uint32_t)rk[j], 0, 0, k0);
            if (sp[j].var) {
                c->scope[c->nscope].name = sp[j].var; c->scope[c->nscope].reg = riv[j];
                c->scope[c->nscope].type = CT_INT; c->scope[c->nscope].built = false;
                c->nscope++;
            }
            Lp[j] = c->n;
            int rc = alloc_temp(c);
            ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rk[j], (uint32_t)rn[j], z);
            jz[j] = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
            c->temp_top--;
        }

        Val bv;
        if (!emit(c, A[0], &bv)) { c->nscope = scope_entry; return -1; }
        if (CT_IS_ARRAY(bv.type)) { c->nscope = scope_entry; c->ok = false; return -1; }
        coerce(c, &bv, et);
        if (!c->ok) { c->nscope = scope_entry; return -1; }
        ins(c,a_store_op(et),
            (uint32_t)rout, (uint32_t)rflat, (uint32_t)bv.reg, z);
        free_if_tmp(c, bv);
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)rflat, 0, 0, one);

        for (int j = rank - 1; j >= 0; j--) {
            Slot sd; memset(&sd, 0, sizeof sd); sd.i = sp[j].di;
            ins(c, OP_INC_I, (uint32_t)riv[j], 0, 0, sd);
            ins(c, OP_INC_I, (uint32_t)rk[j], 0, 0, one);
            ins(c, OP_JMP, 0, 0, (uint32_t)Lp[j], z);
            if (c->ok) c->code[jz[j]].b = (uint32_t)c->n;
        }
        c->nscope = scope_entry;
        c->temp_top = base - c->nlocals;      /* every scalar temp is dead; rout is an arr */
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(et, rank); out->built = true;
        return c->ok ? 1 : -1;
    }

    /* Range[hi] / Range[lo, hi] / Range[lo, hi, di] — a rank-1 integer array.
     *
     * This is `Table[i, {i, lo, hi, di}]` with the body being the iterator
     * itself, so it shares the iteration machinery (emit_iter_count reproduces
     * the interpreter's element COUNT exactly) and simply stores the loop
     * variable instead of a body.
     *
     * INTEGER bounds only, for the reason Table gives: the interpreter walks a
     * real iterator by repeated addition against a termination slack, so a
     * closed-form step differs from it in the last bits and, at the endpoint, in
     * the number of elements. */
    if (strcmp(h, "Range") == 0 && na >= 1 && na <= 3) {
        LoopSpec s;
        if (!range_spec(c, A, na, &s)) { c->ok = false; return -1; }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        Slot one; memset(&one, 0, sizeof one); one.i = 1;

        int rn = alloc_temp(c), rlo = alloc_temp(c);
        if (!emit_iter_count(c, &s, rlo, rn)) return -1;
        emit_max_guard(c, rn, VM_ITER_SAFETY_CAP);

        Slot el; memset(&el, 0, sizeof el); el.i = (long long)CT_INT;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rn, (uint32_t)rn, el);

        int riv = alloc_temp(c), rk = alloc_temp(c);
        ins(c, OP_MOVE,  (uint32_t)riv, (uint32_t)rlo, 0, z);
        ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);

        size_t top = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rk, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        ins(c, OP_A_STORE_I, (uint32_t)rout, (uint32_t)rk, (uint32_t)riv, z);
        Slot sd; memset(&sd, 0, sizeof sd); sd.i = s.di;
        ins(c, OP_INC_I, (uint32_t)riv, 0, 0, sd);
        ins(c, OP_INC_I, (uint32_t)rk, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)top, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        c->temp_top = rn - c->nlocals;         /* every scalar temp is dead */
        out->reg = rout; out->tmp = true;
        out->type = CT_ARRAY(CT_INT, 1); out->built = true;
        return c->ok ? 1 : -1;
    }

    /* Part[a, subscripts...] (M3c).  See the block comment above emit_flat_index
     * for why this splits in two. */
    if (strcmp(h, "Part") == 0 && na >= 2) {
        Val a;
        if (!emit(c, A[0], &a)) return -1;
        if (!CT_IS_ARRAY(a.type) || (int)(na - 1) > CT_RANK(a.type)) { c->ok = false; return -1; }
        CompileType elem = CT_ELEM(a.type);
        Slot z; memset(&z, 0, sizeof z);

        if (part_is_scalar_indexed(c, a.type, (const Expr* const*)A + 1, na - 1)) {
            int ridx;
            if (!emit_flat_index(c, a, (const Expr* const*)A + 1, na - 1, &ridx)) return -1;
            int rd = alloc_temp(c);
            ins(c,a_load_op(elem),
                (uint32_t)rd, (uint32_t)a.reg, (uint32_t)ridx, z);
            /* rd sits above ridx, so drop ridx by relocating rd onto it. */
            ins(c, OP_MOVE, (uint32_t)ridx, (uint32_t)rd, 0, z);
            c->temp_top--;                              /* rd */
            free_if_tmp(c, a);
            out->reg = ridx; out->tmp = true; out->type = elem;
            return c->ok ? 1 : -1;
        }

        int keep = CT_RANK(a.type) - (int)(na - 1);
        for (size_t i = 1; i < na; i++) if (subscript_is_literal_spec(A[i])) keep++;
        if (keep < 1 || keep > CT_MAX_RANK) { c->ok = false; return -1; }

        int base = -1;
        PartSpec* ps = emit_partspec(c, (const Expr* const*)A + 1, na - 1, &base);
        if (!ps) return -1;
        if (!ctx_own_partspec(c, ps)) return -1;
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = ps;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_PART, (uint16_t)(na - 1), (uint32_t)rout, (uint32_t)base,
              (uint32_t)a.reg, ip);
        c->temp_top -= (int)(na - 1);                   /* the subscript registers */
        if (a.tmp) {
            /* The source was itself a temporary, and it sits BELOW the result in
             * the array stack, so it cannot simply be popped: the next alloc_arr
             * would hand out the register the result is living in.  Free it and
             * slide the result down into its slot, restoring LIFO. */
            ins(c, OP_ARR_FREE, (uint32_t)a.reg, 0, 0, z);
            ins(c, OP_A_XFER, (uint32_t)a.reg, (uint32_t)rout, 0, z);
            c->arr_top--;
            rout = a.reg;
        }
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(elem, keep); out->built = a.built;
        return c->ok ? 1 : -1;
    }

    /* Total[v] / Length[v]: the array -> scalar reductions.  Total delegates to
     * the NDArray reduction so its summation order — and therefore its
     * rounding — is identical to the interpreter's Total[]. */
    if ((strcmp(h, "Total") == 0 || strcmp(h, "Length") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        /* Length is dims[0] at any rank — the number of ROWS of a matrix, as in
         * the interpreter.  Total stays rank-1: ndred_total_all collapses every
         * axis, but Total[] reduces only the leading one, so at rank 2 the two
         * would disagree. */
        if (!CT_IS_ARRAY(a.type)
            || (h[0] == 'T' && CT_RANK(a.type) != 1)) { c->ok = false; return -1; }
        Slot z = { 0 };
        *out = (h[0] == 'L') ? arr_op(c, OP_V_LEN, a, arr_noop_val(), CT_INT, z)
                             : arr_op(c, OP_V_TOTAL, a, arr_noop_val(), CT_ELEM(a.type), z);
        return c->ok ? 1 : -1;
    }

    /* Norm[array, p] with a compile-time literal p (Norm[v,1] / Norm[v,Infinity]
     * / Norm[m,"Frobenius"] / ...): delegate to ndla_norm via V_NORM, the p baked
     * into imm.r.  A List-literal operand was already expanded to scalar ops in
     * emit_node, so the operand here is a genuine array; the bare Norm[array]
     * stays on the reduction table below (nd_red_lookup). */
    if (strcmp(h, "Norm") == 0 && na == 2) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        double immr;
        if (!nd_norm2_encode(A[1], a.type, &immr)) { c->ok = false; return -1; }
        Slot ip; memset(&ip, 0, sizeof ip); ip.r = immr;
        *out = arr_op(c, OP_V_NORM, a, arr_noop_val(), CT_REAL, ip);
        return c->ok ? 1 : -1;
    }

    /* Mean / Median / Variance / StandardDeviation / RootMeanSquare / Max /
     * Min of a rank-1 array, delegated to the interpreter's own reduction so
     * the compiled answer is bit-identical to the interpreted one.  RankedMin /
     * RankedMax (nextra == 1) carry a trailing integer rank in a second
     * register, emitted exactly as the A_NDFN trailing-int path does. */
    {
        const NdRedSpec* nr = nd_red_lookup(h, na);
        if (nr) {
            Val a; if (!emit(c, A[0], &a)) return -1;
            CompileType rt = nd_red_result(nr, a.type);
            if ((int)rt < 0) { c->ok = false; return -1; }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = nr;
            if (nr->nextra == 1) {
                Val nv; if (!emit(c, A[1], &nv)) return -1;
                if (nv.type != CT_INT) { c->ok = false; return -1; }
                *out = arr_op(c, OP_V_NDREDN, a, nv, rt, ip);   /* c->a=array, c->b=n */
            } else {
                *out = arr_op(c, OP_V_NDRED, a, arr_noop_val(), rt, ip);
            }
            return c->ok ? 1 : -1;
        }
    }

    return 0;
}
