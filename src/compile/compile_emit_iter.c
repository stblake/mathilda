/* Mathilda — Compile[]: iteration heads (Nest/Fold/NestList/FoldList/FixedPoint/NestWhile, Select/Map/Scan, First/Last).
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

int emit_iter(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e;
    if (strcmp(h, "Nest") == 0 && na == 3) {
        FnSpec fs;
        if (!fn_resolve(A[0], 1, &fs)) { c->ok = false; return -1; }
        CompileType tn; if (!infer_type(c, A[2], &tn) || tn != CT_INT) { c->ok = false; return -1; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return -1; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return -1; }
        CompileType t = (CompileType)tfp;
        Slot z = { 0 };
        /* Persistent registers FIRST, so freeing the init/count temps (which sit
         * above them) cannot clobber them. */
        int racc = alloc_temp(c), rn = alloc_temp(c), rcnt = alloc_temp(c);
        Val vx; if (!emit(c, A[1], &vx)) return -1; coerce(c, &vx, t); if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        Val vn; if (!emit(c, A[2], &vn)) return -1;
        ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
        /* Nest[f, x, -1] is UNEVALUATED (src/funcprog.c:2153); the counted loop
         * below would silently run zero times and hand back the seed. */
        emit_nonneg_guard(c, rn);
        Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, s0);
        size_t Lp = c->n;
        int rc = alloc_temp(c); ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z); c->temp_top--;
        Val acc = { racc, false, t, false }, vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        coerce(c, &vb, t); if (!c->ok) return -1;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        free_if_tmp(c, vb);
        Slot one; one.i = 1; ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;
        c->temp_top = (racc - c->nlocals) + 1;   /* keep racc; free rn, rcnt */
        out->reg = racc; out->tmp = true; out->type = t; return c->ok ? 1 : -1;
    }

    /* Reverse / Sort / Accumulate / Flatten / Transpose / Take / Drop over a
     * machine array, delegated to the interpreter's own entry point. */
    {
        const NdFnSpec* nf = nd_fn_lookup(h, na);
        if (nf) {
            Val a;
            if (!emit(c, A[0], &a)) return -1;
            CompileType rt = nd_fn_result(nf, a.type);
            if ((int)rt < 0) { c->ok = false; return -1; }
            Slot z = { 0 };
            int base = 0;
            for (int i = 0; i < nf->nextra; i++) {
                int r = alloc_temp(c);
                if (i == 0) base = r;
                Val v;
                if (!emit(c, A[1 + i], &v)) return -1;
                if (v.type != CT_INT) { c->ok = false; return -1; }
                ins(c, OP_MOVE, (uint32_t)r, (uint32_t)v.reg, 0, z);
                free_if_tmp(c, v);
            }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = nf;
            int rout = alloc_arr(c);
            ins_f(c, OP_A_NDFN, (uint16_t)nf->nextra, (uint32_t)rout,
                  (uint32_t)base, (uint32_t)a.reg, ip);
            c->temp_top -= nf->nextra;
            if (a.tmp) {   /* restore LIFO — the same slide the Part lowering does */
                ins(c, OP_ARR_FREE, (uint32_t)a.reg, 0, 0, z);
                ins(c, OP_A_XFER, (uint32_t)a.reg, (uint32_t)rout, 0, z);
                c->arr_top--;
                rout = a.reg;
            }
            out->reg = rout; out->tmp = true; out->type = rt; out->built = a.built;
            return c->ok ? 1 : -1;
        }
    }

    /* Dot / LinearSolve / Cross / LeastSquares / ListConvolve / ListCorrelate /
     * Join over TWO machine arrays, delegated (COMPILE_MISSING.md §3).  arr_op
     * reads both array registers, frees each per its tmp flag and allocates the
     * result — an array bank register (A_NDFN2) or a scalar temp (V_NDFN2, the
     * vector.vector Dot) depending on whether nd_fn2_result gave an array type. */
    {
        const NdFn2Spec* nf2 = nd_fn2_lookup(h, na);
        if (nf2) {
            Val a, b;
            if (!emit(c, A[0], &a)) return -1;
            if (!emit(c, A[1], &b)) return -1;
            CompileType rt = nd_fn2_result(nf2, a.type, b.type);
            if ((int)rt < 0) { c->ok = false; return -1; }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = nf2;
            *out = arr_op(c, CT_IS_ARRAY(rt) ? OP_A_NDFN2 : OP_V_NDFN2, a, b, rt, ip);
            return c->ok ? 1 : -1;
        }
    }

    /* Select / TakeWhile / LengthWhile / AllTrue / AnyTrue / NoneTrue over a
     * rank-1 array: one predicate loop, four things done with the answer.
     *
     * Select and TakeWhile produce a buffer whose length is only known once the
     * loop has run, so they allocate the upper bound the source gives them and
     * cut it to size (A_TRUNC).  An EMPTY result declines: the interpreter
     * cannot pack `{}` either, so it answers with a List, and a length-0 array
     * would be a different value.
     *
     * (These only became compilable once the interpreter itself grew NDArray
     * paths for them — every one used to return the call UNEVALUATED on a
     * packed argument, so there was nothing to be parity with.) */
    {
        int sel = -1;
        if (na == 2) {
            if      (strcmp(h, "Select") == 0)      sel = 0;
            else if (strcmp(h, "TakeWhile") == 0)   sel = 1;
            else if (strcmp(h, "LengthWhile") == 0) sel = 2;
            else if (strcmp(h, "AllTrue") == 0)     sel = 3;
            else if (strcmp(h, "AnyTrue") == 0)     sel = 4;
            else if (strcmp(h, "NoneTrue") == 0)    sel = 5;
        }
        if (sel >= 0) {
            bool keeps = (sel <= 1);            /* builds a buffer */
            bool pred_run = (sel <= 2);         /* scans a prefix / filters */
            FnSpec ps;
            if (!fn_resolve(A[1], 1, &ps)) { c->ok = false; return -1; }
            CompileType el = vec_elem_type(c, A[0]);
            if ((int)el < 0) { c->ok = false; return -1; }
            CompileType tb;
            if (!infer_apply(c, &ps, &el, 1, &tb) || tb != CT_BOOL) { c->ok = false; return -1; }

            Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
            Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
            int rn = alloc_temp(c), ri = alloc_temp(c), rk = alloc_temp(c);
            Val va;
            if (!emit(c, A[0], &va)) return -1;
            if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return -1; }
            ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)va.reg, 0, z);

            int rout = -1;
            if (keeps) {
                Slot el_i; memset(&el_i, 0, sizeof el_i); el_i.i = (long long)el;
                rout = alloc_arr(c);
                ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rn, (uint32_t)rn, el_i);
            }
            /* Boolean accumulator: All and None start True, Any starts False. */
            if (!pred_run) {
                Slot init; memset(&init, 0, sizeof init); init.i = (sel == 4) ? 0 : 1;
                ins(c, OP_CONST, (uint32_t)rk, 0, 0, init);
            } else ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);

            size_t Lp = c->n;
            int rc = alloc_temp(c);
            ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
            c->temp_top--;

            int relem = alloc_temp(c);
            int body_top = c->temp_top;
            ins(c,a_load_op(el),
                (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ri, z);
            Val ev = { relem, false, el, false }, vt;
            if (!emit_apply(c, &ps, &ev, 1, &vt)) return -1;
            if (vt.type != CT_BOOL) { c->ok = false; return -1; }

            size_t jstop = 0; bool has_stop = false;
            if (sel == 0) {                       /* Select: keep when true */
                size_t jskip = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
                ins(c,a_store_op(el),
                    (uint32_t)rout, (uint32_t)rk, (uint32_t)relem, z);
                ins(c, OP_INC_I, (uint32_t)rk, 0, 0, k1);
                if (c->ok) c->code[jskip].b = (uint32_t)c->n;
            } else if (sel == 1 || sel == 2) {
                /* A prefix, so the write index and the final length are both
                 * just `ri` — no separate counter to keep in step. */
                jstop = c->n; has_stop = true;
                ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
                if (sel == 1)
                    ins(c,a_store_op(el),
                        (uint32_t)rout, (uint32_t)ri, (uint32_t)relem, z);
            } else {
                /* Short-circuit.  AllTrue fires on a FALSE test, AnyTrue and
                 * NoneTrue on a TRUE one; the value they then answer with is
                 * True only for AnyTrue. */
                Slot outv; memset(&outv, 0, sizeof outv); outv.i = (sel == 4) ? 1 : 0;
                size_t jfalse = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
                size_t jcont = 0;
                if (sel == 3) {                   /* AllTrue: true -> keep going */
                    jcont = c->n; ins(c, OP_JMP, 0, 0, 0, z);
                    if (c->ok) c->code[jfalse].b = (uint32_t)c->n;   /* false -> fire */
                }
                ins(c, OP_CONST, (uint32_t)rk, 0, 0, outv);
                jstop = c->n; has_stop = true; ins(c, OP_JMP, 0, 0, 0, z);
                if (c->ok) {
                    if (sel == 3) c->code[jcont].b = (uint32_t)c->n;
                    else          c->code[jfalse].b = (uint32_t)c->n; /* false -> keep going */
                }
            }
            c->temp_top = body_top - 1;
            ins(c, OP_INC_I, (uint32_t)ri, 0, 0, k1);
            ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
            if (c->ok) {
                c->code[jz].b = (uint32_t)c->n;
                if (has_stop) c->code[jstop].b = (uint32_t)c->n;
            }

            if (keeps) {
                int rlen = (sel == 0) ? rk : ri;
                /* An empty result has no packed form; the interpreter answers
                 * with a List, so a length-0 array would not be the same value. */
                emit_nonzero_guard(c, rlen);
                ins_f(c, OP_A_TRUNC, 0, (uint32_t)rout, (uint32_t)rlen, 0, z);
                if (va.tmp) {          /* restore LIFO, as the Part lowering does */
                    ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
                    ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)rout, 0, z);
                    c->arr_top--;
                    rout = va.reg;
                }
                c->temp_top = (rn - c->nlocals);
                out->reg = rout; out->tmp = true; out->type = CT_ARRAY(el, 1);
                out->built = va.built;     /* filtered, not constructed */
                return c->ok ? 1 : -1;
            }
            free_if_tmp(c, va);
            /* LengthWhile's answer is the prefix length, which IS `ri`; the
             * predicates' is the boolean in `rk`.  Relocate it down onto rn so
             * the scratch above is reclaimed. */
            ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)(sel == 2 ? ri : rk), 0, z);
            c->temp_top = (rn - c->nlocals) + 1;
            out->reg = rn; out->tmp = true; out->built = false;
            out->type = (sel == 2) ? CT_INT : CT_BOOL;
            return c->ok ? 1 : -1;
        }
    }

    /* First[v] / Last[v]: exactly v[[1]] and v[[-1]], so lower them as that
     * rather than duplicating the indexing — the range check that makes
     * First[{}] decline comes along for free. */
    if ((strcmp(h, "First") == 0 || strcmp(h, "Last") == 0) && na == 1) {
        CompileType ta;
        if (!infer_type(c, A[0], &ta) || !CT_IS_ARRAY(ta)) { c->ok = false; return -1; }
        Expr* idx = expr_new_integer(h[0] == 'F' ? 1 : -1);
        Expr* args[2] = { expr_copy(A[0]), idx };
        Expr* part = expr_new_function(expr_new_symbol("Part"), args, 2);
        if (!part) { expr_free(idx); c->ok = false; return -1; }
        bool r = emit(c, part, out);
        if (!r && expr_subtree_of(part, c->bail_node)) c->bail_node = e;
        expr_free(part);
        return r;
    }

    /* Map[f, v] / Scan[f, v] over a rank-1 array.
     *
     * One element loop serves both: Map stores each result into a fresh buffer,
     * Scan drops it and answers like Do.  Rank 1 only — at rank >= 2 the
     * interpreter's Map applies f to each ROW (map_ndarray_axis,
     * src/funcprog.c:272), which is a different operation from elementwise.
     *
     * The result element type must equal the SOURCE's, because Map over a packed
     * array repacks with the source dtype (src/funcprog.c:289): Map[Abs, cv] on
     * a complex vector comes back complex-typed with real values.  Anything else
     * would answer with a different element type, so it declines. */
    if ((strcmp(h, "Map") == 0 || strcmp(h, "Scan") == 0) && na == 2) {
        bool scan = h[0] == 'S';
        FnSpec fs;
        if (!fn_resolve(A[0], 1, &fs)) { c->ok = false; return -1; }
        CompileType el = vec_elem_type(c, A[1]);
        if ((int)el < 0) { c->ok = false; return -1; }
        CompileType rt;
        if (!infer_apply(c, &fs, &el, 1, &rt) || CT_IS_ARRAY(rt)) { c->ok = false; return -1; }
        if (!scan && rt != el) { c->ok = false; return -1; }
        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;

        int rn = alloc_temp(c), ri = alloc_temp(c);
        Val va;
        if (!emit(c, A[1], &va)) return -1;
        if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return -1; }
        ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)va.reg, 0, z);
        /* An empty source makes the interpreter's Map return a plain List — the
         * result of ndarray_from_nested_list declining on List[] — where a
         * compiled buffer would be an empty NDArray. */
        if (!scan) emit_nonzero_guard(c, rn);

        /* Fast route.  When the body threads elementwise, Map[Function[u, body], v]
         * IS `body` with u bound to the whole array — and that is what the
         * existing elementwise fusion strip-mines into one tiled, optionally
         * threaded pass.  Measured 6.2x over the per-element loop below on
         * `#^2 + 1.` at 200k elements, for ten lines and no new machinery.
         *
         * Three conditions make the rewrite legal, and each corresponds to a way
         * the interpreter's Map differs from threading:
         *   - fuse_listable == 1: the interpreter threads a list exactly when the
         *     head is Listable, so `Map[Function[u, Total[u]], v]` (Total of each
         *     SCALAR element) must not become Total[v);
         *   - exactly one array leaf, and it is the parameter: a second array
         *     would make Map produce a list OF arrays, not an elementwise result;
         *   - the result is rank 1 with the promised element type.
         * Anything else rolls back and takes the general loop, which is correct
         * for every body — this is a cost split, not a subset restriction. */
        if (!scan && !(c->flags & COMPILE_NO_FUSE)
            && fs.kind == FN_LAMBDA && fs.nparams == 1 && c->nscope < CTX_MAX_SCOPE) {
            EmitMark mk = emit_mark(c);
            c->scope[c->nscope].name = fs.pname[0]; c->scope[c->nscope].reg = va.reg;
            c->scope[c->nscope].type = va.type;     c->scope[c->nscope].built = va.built;
            c->nscope++;
            FuseLeaves L; L.n = 0;
            bool legal = fuse_listable(c, fs.body) == 1
                      && fuse_collect(c, fs.body, &L) && L.n == 1
                      && L.name[0] == fs.pname[0];
            Val bv;
            bool got = legal && emit(c, fs.body, &bv) && CT_IS_ARRAY(bv.type)
                    && CT_RANK(bv.type) == 1 && CT_ELEM(bv.type) == rt;
            c->nscope--;
            if (got) {
                int res = bv.reg;
                if (va.tmp && bv.tmp && res != va.reg) {  /* restore LIFO, as Part does */
                    ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
                    ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)res, 0, z);
                    c->arr_top--;
                    res = va.reg;
                }
                c->temp_top = (rn - c->nlocals);          /* rn, ri unused on this route */
                out->reg = res; out->tmp = bv.tmp; out->type = bv.type;
                out->built = va.built;
                return c->ok ? 1 : -1;
            }
            emit_rollback(c, mk);
        }

        int rout = -1;
        if (!scan) {
            rout = alloc_arr(c);
            ins_f(c, OP_A_NEWLIKE, (uint16_t)(((unsigned)rt & 3u) << AF_R_SHIFT),
                  (uint32_t)rout, (uint32_t)va.reg, 0, z);
        }
        ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        int relem = alloc_temp(c);
        int body_top = c->temp_top;
        ins(c,a_load_op(el),
            (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ri, z);
        Val ev = { relem, false, el, false }, vb;
        if (!emit_apply(c, &fs, &ev, 1, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        if (!scan) {
            coerce(c, &vb, rt); if (!c->ok) return -1;
            ins(c,a_store_op(rt),
                (uint32_t)rout, (uint32_t)ri, (uint32_t)vb.reg, z);
        }
        c->temp_top = body_top - 1;                     /* body temps and relem */
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)ri, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        if (scan) {
            free_if_tmp(c, va);
            c->temp_top = (rn - c->nlocals);
            int r0 = alloc_temp(c); ins(c, OP_CONST, (uint32_t)r0, 0, 0, k0);
            out->reg = r0; out->tmp = true; out->type = CT_INT; out->built = false;
            return c->ok ? 1 : -1;
        }
        /* The source may itself be a temporary sitting BELOW the result in the
         * array stack, so it cannot just be popped: free it and slide the result
         * down into its slot, restoring LIFO.  Same dance as the Part lowering. */
        if (va.tmp) {
            ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
            ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)rout, 0, z);
            c->arr_top--;
            rout = va.reg;
        }
        c->temp_top = (rn - c->nlocals);                /* rn, ri dead; rout is an arr */
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(rt, 1);
        out->built = va.built;                          /* Map follows its source's kind */
        return c->ok ? 1 : -1;
    }

    /* Fold[f, x0, v] / Fold[f, v]: reduce a rank-1 array into a scalar
     * accumulator.  Nest's loop with an index and a per-iteration element. */
    if (strcmp(h, "Fold") == 0 && (na == 2 || na == 3)) {
        FnSpec fs;
        if (!fn_resolve(A[0], 2, &fs)) { c->ok = false; return -1; }
        CompileType el = vec_elem_type(c, A[na - 1]);
        if ((int)el < 0) { c->ok = false; return -1; }
        CompileType t0 = el;
        if (na == 3 && !infer_type(c, A[1], &t0)) { c->ok = false; return -1; }
        int tfp = accum_fixed_type(c, &fs, t0, &el, 1);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return -1; }
        CompileType T = (CompileType)tfp;
        Slot z = { 0 };

        /* Persistent registers first, so the operand lowerings below allocate
         * ABOVE them and freeing those temps cannot reach them. */
        int racc = alloc_temp(c), rn = alloc_temp(c), ri = alloc_temp(c);
        Val va;
        if (!emit(c, A[na - 1], &va)) return -1;
        if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return -1; }
        ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)va.reg, 0, z);

        uint16_t ldop = a_load_op(el);
        Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        if (na == 3) {
            Val vx; if (!emit(c, A[1], &vx)) return -1;
            coerce(c, &vx, T); if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);
        } else {
            /* Fold[f, {}] stays UNEVALUATED (src/funcprog.c:2282), so an empty
             * vector must fail the call rather than answer with a seed. */
            emit_nonzero_guard(c, rn);
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);
            int rs = alloc_temp(c);
            ins(c, ldop, (uint32_t)rs, (uint32_t)va.reg, (uint32_t)ri, z);
            Val sv = { rs, true, el, false };
            coerce(c, &sv, T); if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)sv.reg, 0, z);
            free_if_tmp(c, sv);
            Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k1);   /* the seed is consumed */
        }

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;                                  /* guard temp is dead */

        int relem = alloc_temp(c);
        int body_top = c->temp_top;
        ins(c, ldop, (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ri, z);
        Val argv[2] = { { racc, false, T, false }, { relem, false, el, false } }, vb;
        if (!emit_apply(c, &fs, argv, 2, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        coerce(c, &vb, T); if (!c->ok) return -1;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        c->temp_top = body_top - 1;                     /* drop the body's temps and relem */
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)ri, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        free_if_tmp(c, va);                             /* release the source buffer */
        c->temp_top = (racc - c->nlocals) + 1;          /* keep racc */
        out->reg = racc; out->tmp = true; out->type = T; out->built = false;
        return c->ok ? 1 : -1;
    }

    /* NestList[f, x, n] / FoldList[f, x0, v] / FoldList[f, v]: the same
     * accumulator loops as Nest and Fold, writing every iterate into a buffer
     * whose length is known before the loop starts (n + 1, or the source length
     * plus one for the seed).
     *
     * The element type must be Real or Complex — these BUILD their result, so
     * the ConstantArray rule applies: `NestList[2 # &, 1, 5]` holds exact
     * Integers in the interpreter and a packed buffer has no integer dtype. */
    if ((strcmp(h, "NestList") == 0 && na == 3)
        || (strcmp(h, "FoldList") == 0 && (na == 2 || na == 3))) {
        bool nest = h[0] == 'N';
        FnSpec fs;
        if (!fn_resolve(A[0], nest ? 1 : 2, &fs)) { c->ok = false; return -1; }

        CompileType el = CT_ERR, T;
        if (nest) {
            CompileType tn;
            if (!infer_type(c, A[2], &tn) || tn != CT_INT) { c->ok = false; return -1; }
            CompileType tx;
            if (!infer_type(c, A[1], &tx)) { c->ok = false; return -1; }
            int tfp = nest_fixed_type(c, &fs, tx);
            if (tfp < 0) { c->ok = false; return -1; }
            T = (CompileType)tfp;
        } else {
            el = vec_elem_type(c, A[na - 1]);
            if ((int)el < 0) { c->ok = false; return -1; }
            CompileType t0 = el;
            if (na == 3 && !infer_type(c, A[1], &t0)) { c->ok = false; return -1; }
            int tfp = accum_fixed_type(c, &fs, t0, &el, 1);
            if (tfp < 0) { c->ok = false; return -1; }
            T = (CompileType)tfp;
        }
        if (!ct_is_elem(T)) { c->ok = false; return -1; }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
        uint16_t stop = a_store_op(T);

        /* Persistent registers first; `rlen` must be the sole dimension operand
         * of A_NEW, which reads a contiguous run starting there. */
        int rlen = alloc_temp(c), racc = alloc_temp(c);
        int rn = alloc_temp(c), rcnt = alloc_temp(c), rone = alloc_temp(c);
        ins(c, OP_CONST, (uint32_t)rone, 0, 0, k1);

        Val va = { 0, false, CT_REAL, false };
        if (nest) {
            Val vn; if (!emit(c, A[2], &vn)) return -1;
            if (vn.type != CT_INT) { c->ok = false; return -1; }
            ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
            emit_nonneg_guard(c, rn);                    /* n < 0 is unevaluated */
            ins(c, OP_ADD_I, (uint32_t)rlen, (uint32_t)rn, (uint32_t)rone, z);
        } else {
            if (!emit(c, A[na - 1], &va)) return -1;
            if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return -1; }
            ins(c, OP_A_SIZE, (uint32_t)rlen, (uint32_t)va.reg, 0, z);
            if (na == 3) {
                /* seeded: the history is the seed plus one entry per element */
                ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)rlen, 0, z);
                ins(c, OP_ADD_I, (uint32_t)rlen, (uint32_t)rlen, (uint32_t)rone, z);
            } else {
                /* seedless: the first element IS the seed, so one fewer step.
                 * FoldList[f, {}] is `{}` in the interpreter, and an empty
                 * packed result is not the same value, so decline it. */
                emit_nonzero_guard(c, rlen);
                ins(c, OP_SUB_I, (uint32_t)rn, (uint32_t)rlen, (uint32_t)rone, z);
            }
        }
        emit_max_guard(c, rlen, 1 << 26);   /* bound the allocation, not the answer */

        Slot el_i; memset(&el_i, 0, sizeof el_i); el_i.i = (long long)T;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rlen, (uint32_t)rlen, el_i);

        /* seed into racc, and into slot 0 of the buffer */
        if (nest) {
            Val vx; if (!emit(c, A[1], &vx)) return -1;
            coerce(c, &vx, T); if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        } else if (na == 3) {
            Val vx; if (!emit(c, A[1], &vx)) return -1;
            coerce(c, &vx, T); if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        } else {
            int rs = alloc_temp(c);
            ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);
            ins(c,a_load_op(el),
                (uint32_t)rs, (uint32_t)va.reg, (uint32_t)rcnt, z);
            Val sv = { rs, true, el, false };
            coerce(c, &sv, T); if (!c->ok) return -1;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)sv.reg, 0, z);
            free_if_tmp(c, sv);
        }
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);
        ins(c, stop, (uint32_t)rout, (uint32_t)rcnt, (uint32_t)racc, z);

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        int body_top = c->temp_top, relem = -1;
        Val argv[2] = { { racc, false, T, false }, { 0, false, CT_REAL, false } }, vb;
        int nargs = 1;
        if (!nest) {
            /* Element k of the source pairs with history slot k+1: seeded folds
             * consume v[[k+1]] at step k, seedless ones v[[k+2]] (the first
             * element was the seed). */
            relem = alloc_temp(c);
            body_top = c->temp_top;
            int ridx = alloc_temp(c);
            if (na == 3) ins(c, OP_MOVE, (uint32_t)ridx, (uint32_t)rcnt, 0, z);
            else         ins(c, OP_ADD_I, (uint32_t)ridx, (uint32_t)rcnt, (uint32_t)rone, z);
            ins(c,a_load_op(el),
                (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ridx, z);
            c->temp_top--;                                /* ridx dead after the load */
            argv[1].reg = relem; argv[1].type = el;
            nargs = 2;
        }
        if (!emit_apply(c, &fs, argv, nargs, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        coerce(c, &vb, T); if (!c->ok) return -1;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        c->temp_top = nest ? body_top : body_top - 1;
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, k1);
        ins(c, stop, (uint32_t)rout, (uint32_t)rcnt, (uint32_t)racc, z);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        if (!nest && va.tmp) {           /* restore LIFO, as the Part lowering does */
            ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
            ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)rout, 0, z);
            c->arr_top--;
            rout = va.reg;
        }
        c->temp_top = (rlen - c->nlocals);
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(T, 1);
        /* NestList always constructs; FoldList's history is packed by the
         * interpreter with the SOURCE dtype, so it follows the source's kind. */
        out->built = nest ? true : va.built;
        return c->ok ? 1 : -1;
    }

    /* FixedPointList[f, x, ...] / NestWhileList[f, x, test]: the same loops as
     * their scalar twins, keeping every iterate.
     *
     * The length is not known until the loop has run, so the buffer GROWS
     * (A_PUSH) and is cut to size at the end (A_TRUNC).  Allocating the safety
     * cap up front instead would be 8 MB per call, and running the body twice
     * to count first would double every side effect a `Set` in it performs.
     * The capacity a particular run reached is never observable. */
    if ((strcmp(h, "FixedPointList") == 0 && na >= 2 && na <= 4)
        || (strcmp(h, "NestWhileList") == 0 && na == 3)) {
        bool fp = h[0] == 'F';
        FnSpec fs, ss, ts;
        const Expr *mx = NULL, *st = NULL;
        if (!fn_resolve(A[0], 1, &fs)) { c->ok = false; return -1; }
        if (fp) { if (!fp_opts(A, na, 2, &mx, &st)) { c->ok = false; return -1; } }
        else if (!fn_resolve(A[2], 1, &ts)) { c->ok = false; return -1; }

        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return -1; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0) { c->ok = false; return -1; }
        CompileType T = (CompileType)tfp;
        /* A built history, so the ConstantArray element-type rule applies. */
        if (!ct_is_elem(T)) { c->ok = false; return -1; }
        if (fp && st && !fn_resolve(st, 2, &ss)) { c->ok = false; return -1; }
        if (!fp) {   /* the while-test must yield a Bool on the accumulator */
            CompileType tb;
            if (!infer_apply(c, &ts, &T, 1, &tb) || tb != CT_BOOL) { c->ok = false; return -1; }
        }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
        uint16_t pflags = (uint16_t)(((unsigned)T & 3u) << AF_R_SHIFT);

        int rcap = alloc_temp(c), racc = alloc_temp(c);
        int rk = alloc_temp(c), rcnt = alloc_temp(c), rlim = alloc_temp(c);
        Slot cap0; memset(&cap0, 0, sizeof cap0); cap0.i = 16;   /* grown as needed */
        ins(c, OP_CONST, (uint32_t)rcap, 0, 0, cap0);

        Val vx; if (!emit(c, A[1], &vx)) return -1;
        coerce(c, &vx, T); if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        if (mx) {
            Val vn; if (!emit(c, mx, &vn)) return -1;
            if (vn.type != CT_INT) { c->ok = false; return -1; }
            ins(c, OP_MOVE, (uint32_t)rlim, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
            emit_nonneg_guard(c, rlim);
        } else {
            Slot cp; memset(&cp, 0, sizeof cp); cp.i = VM_ITER_SAFETY_CAP;
            ins(c, OP_CONST, (uint32_t)rlim, 0, 0, cp);
        }

        Slot el_i; memset(&el_i, 0, sizeof el_i); el_i.i = (long long)T;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rcap, (uint32_t)rcap, el_i);

        ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);
        ins_f(c, OP_A_PUSH, pflags, (uint32_t)rout, (uint32_t)rk, (uint32_t)racc, z);
        ins(c, OP_INC_I, (uint32_t)rk, 0, 0, k1);

        size_t Lp = c->n, jend = 0, jcap = 0;
        int body_top = c->temp_top;
        if (!fp) {                       /* NestWhileList tests BEFORE applying */
            Val acc0 = { racc, false, T, false }, vt;
            if (!emit_apply(c, &ts, &acc0, 1, &vt)) return -1;
            if (vt.type != CT_BOOL) { c->ok = false; return -1; }
            jend = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
            c->temp_top = body_top;
        }
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rlim, z);
        jcap = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        Val acc = { racc, false, T, false }, vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        coerce(c, &vb, T); if (!c->ok) return -1;
        int rsame = -1;
        if (fp) {                        /* compare before the accumulator moves */
            if (st) {
                Val sargv[2] = { { racc, false, T, false }, vb }, sv;
                sargv[1].tmp = false;
                if (!emit_apply(c, &ss, sargv, 2, &sv)) return -1;
                if (sv.type != CT_BOOL) { c->ok = false; return -1; }
                rsame = sv.reg;
            } else {
                rsame = alloc_temp(c);
                ins(c, emit_sameq_op(T), (uint32_t)rsame, (uint32_t)vb.reg, (uint32_t)racc, z);
            }
        }
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, k1);
        ins_f(c, OP_A_PUSH, pflags, (uint32_t)rout, (uint32_t)rk, (uint32_t)racc, z);
        ins(c, OP_INC_I, (uint32_t)rk, 0, 0, k1);
        if (fp) ins(c, OP_JZ, 0, (uint32_t)rsame, (uint32_t)Lp, z);  /* not same -> iterate */
        else    ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        c->temp_top = body_top;
        size_t jdone = c->n; ins(c, OP_JMP, 0, 0, 0, z);
        if (c->ok) c->code[jcap].b = (uint32_t)c->n;
        /* Reaching the cap of an UNBOUNDED run is where the interpreter gives
         * up too; with a user bound, falling through keeps what was collected. */
        if (!mx) ins(c, OP_FAIL, 0, 0, 0, z);
        if (c->ok) {
            c->code[jdone].b = (uint32_t)c->n;
            if (!fp) c->code[jend].b = (uint32_t)c->n;
        }
        ins_f(c, OP_A_TRUNC, 0, (uint32_t)rout, (uint32_t)rk, 0, z);

        c->temp_top = (rcap - c->nlocals);
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(T, 1); out->built = true;
        return c->ok ? 1 : -1;
    }

    /* FixedPoint[f, x] / FixedPoint[f, x, n] / SameTest -> s: iterate until two
     * successive values are SameQ.  See emit_sameq for why that is not Equal. */
    if (strcmp(h, "FixedPoint") == 0 && na >= 2 && na <= 4) {
        FnSpec fs, ss; const Expr *mx, *st;
        if (!fn_resolve(A[0], 1, &fs) || !fp_opts(A, na, 2, &mx, &st)) { c->ok = false; return -1; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return -1; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return -1; }
        CompileType T = (CompileType)tfp;
        if (st && !fn_resolve(st, 2, &ss)) { c->ok = false; return -1; }
        if (!st && T == CT_BOOL) { c->ok = false; return -1; }   /* no SameQ opcode for Bool */
        Slot z = { 0 };

        int racc = alloc_temp(c), rcnt = alloc_temp(c), rlim = alloc_temp(c);
        Val vx; if (!emit(c, A[1], &vx)) return -1;
        coerce(c, &vx, T); if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        if (mx) {
            Val vn; if (!emit(c, mx, &vn)) return -1;
            if (vn.type != CT_INT) { c->ok = false; return -1; }
            ins(c, OP_MOVE, (uint32_t)rlim, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
            /* A negative bound leaves the whole call unevaluated. */
            emit_nonneg_guard(c, rlim);
        } else {
            Slot cap; memset(&cap, 0, sizeof cap); cap.i = VM_ITER_SAFETY_CAP;
            ins(c, OP_CONST, (uint32_t)rlim, 0, 0, cap);
        }
        Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rlim, z);
        size_t jcap = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        int body_top = c->temp_top;
        Val acc = { racc, false, T, false }, vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        coerce(c, &vb, T); if (!c->ok) return -1;
        /* Compare BEFORE the accumulator is overwritten. */
        int rsame;
        if (st) {
            Val sargv[2] = { { racc, false, T, false }, vb }, sv;
            sargv[1].tmp = false;
            if (!emit_apply(c, &ss, sargv, 2, &sv)) return -1;
            if (sv.type != CT_BOOL) { c->ok = false; return -1; }
            rsame = sv.reg;
        } else {
            rsame = alloc_temp(c);
            ins(c, emit_sameq_op(T), (uint32_t)rsame, (uint32_t)vb.reg, (uint32_t)racc, z);
        }
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JZ, 0, (uint32_t)rsame, (uint32_t)Lp, z);   /* not same -> iterate */
        c->temp_top = body_top;
        size_t jend = c->n; ins(c, OP_JMP, 0, 0, 0, z);
        if (c->ok) c->code[jcap].b = (uint32_t)c->n;
        /* Reaching the cap of an UNBOUNDED run is where the interpreter gives up
         * too, so fail the call rather than answer with a non-fixed point.  With
         * a user bound, falling through returns f^n(x), which is what
         * FixedPoint[f, x, n] means. */
        if (!mx) ins(c, OP_FAIL, 0, 0, 0, z);
        if (c->ok) c->code[jend].b = (uint32_t)c->n;

        c->temp_top = (racc - c->nlocals) + 1;
        out->reg = racc; out->tmp = true; out->type = T; out->built = false;
        return c->ok ? 1 : -1;
    }

    /* NestWhile[f, x, test]: apply f while test holds on the CURRENT value, so
     * the test runs before the first application (src/funcprog.c:2323). */
    if (strcmp(h, "NestWhile") == 0 && na == 3) {
        FnSpec fs, ts;
        if (!fn_resolve(A[0], 1, &fs) || !fn_resolve(A[2], 1, &ts)) { c->ok = false; return -1; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return -1; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return -1; }
        CompileType T = (CompileType)tfp;
        Slot z = { 0 };

        int racc = alloc_temp(c), rcnt = alloc_temp(c), rlim = alloc_temp(c);
        Val vx; if (!emit(c, A[1], &vx)) return -1;
        coerce(c, &vx, T); if (!c->ok) return -1;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        Slot cap; memset(&cap, 0, sizeof cap); cap.i = VM_ITER_SAFETY_CAP;
        ins(c, OP_CONST, (uint32_t)rlim, 0, 0, cap);
        Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);

        size_t Lp = c->n;
        int body_top = c->temp_top;
        Val acc = { racc, false, T, false }, vt;
        if (!emit_apply(c, &ts, &acc, 1, &vt)) return -1;
        if (vt.type != CT_BOOL) { c->ok = false; return -1; }
        size_t jend = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
        c->temp_top = body_top;

        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rlim, z);
        size_t jcap = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        Val vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return -1;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return -1; }
        coerce(c, &vb, T); if (!c->ok) return -1;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        c->temp_top = body_top;
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) { c->code[jcap].b = (uint32_t)c->n; }
        ins(c, OP_FAIL, 0, 0, 0, z);                    /* cap: interpreter gives up too */
        if (c->ok) c->code[jend].b = (uint32_t)c->n;

        c->temp_top = (racc - c->nlocals) + 1;
        out->reg = racc; out->tmp = true; out->type = T; out->built = false;
        return c->ok ? 1 : -1;
    }
    return 0;
}
