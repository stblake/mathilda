/* Mathilda — Compile[]: control-flow heads (If, Sum/Product, CompoundExpression, With/Module, Set/AddTo/..., Do, While/For).
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

/* Thread an already-emitted scalar `val` into every target of a (possibly
 * nested) List LHS -- Wolfram Set semantics for a non-List RHS: {a, b} = c
 * binds a = c and b = c, and {a, {b, c}} = c recurses. Each leaf target must be
 * a scoped numeric local. `val` is BORROWED (left intact) so the caller can
 * return it as the Set's value and free it exactly once; the per-target widening
 * copy is coerced into its own temp and freed here (LIFO). Returns false with
 * c->ok = false on any target that is not a compilable local, so the whole body
 * bails to the interpreter (which threads correctly). */
static bool emit_list_thread(Ctx* c, const Expr* lhs, Val val) {
    Slot z; memset(&z, 0, sizeof z);
    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        const Expr* t = lhs->data.function.args[i];
        if (t->type == EXPR_SYMBOL) {
            CompileType vt;
            int vreg = scope_find(c, t->data.symbol.name, &vt, NULL);
            if (vreg < 0 || vt == CT_BOOL) { c->ok = false; return false; }
            Val cp = { val.reg, false, val.type, false };
            coerce(c, &cp, vt);                 /* widen a borrowed copy to the target */
            if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)vreg, (uint32_t)cp.reg, 0, z);
            if (cp.reg != val.reg) c->temp_top--;   /* free the coercion temp */
        } else if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL
                   && strcmp(t->data.function.head->data.symbol.name, "List") == 0) {
            if (!emit_list_thread(c, t, val)) return false;
        } else {
            c->ok = false; return false;         /* not a compilable target */
        }
    }
    return true;
}

int emit_ctrl(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e;
    /*
     * If[cond, then] -- no else branch. A STATEMENT, not a value: the
     * interpreter answers Null when cond is False, so this lowers exactly like
     * While (guard, jump over the body, answer the integer 0) and joins
     * stmt_valued_head so the 0 can never be the program's result.
     *
     * Omitting this was not a small gap. The compilable subset is a cliff, not a
     * slope -- one head outside it costs the WHOLE body -- and `If[test, var =
     * val]` is how anyone writes a running maximum or a conditional store. Both
     * a Sieve of Eratosthenes and a Collatz longest-chain search failed to
     * compile on this single form, and so ran interpreted: Collatz to 10^6 took
     * 240 s instead of 0.35 s. docs/design/compile.md section 11 specified it
     * ("If[c,t] -> f = Null only valid where the value is unused"); it was simply
     * never emitted.
     */
    if (strcmp(h, "If") == 0 && na == 2) {
        Slot z = { 0 };
        Val cond; if (!emit(c, A[0], &cond)) return -1;
        if (cond.type != CT_BOOL) { c->ok = false; return -1; }
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)cond.reg, 0, z);
        free_if_tmp(c, cond);
        Val th; if (!emit(c, A[1], &th)) return -1;
        free_if_tmp(c, th);                    /* value discarded, as While's body is */
        if (c->ok) c->code[jz].b = (uint32_t)c->n;   /* end label */
        int r0 = alloc_temp(c); Slot s0; s0.i = 0;
        ins(c, OP_CONST, (uint32_t)r0, 0, 0, s0);
        out->reg = r0; out->tmp = true; out->type = CT_INT;
        return c->ok ? 1 : -1;
    }

    /* If[cond, then, else]: branch control flow.  Result lands in one register
     * that whichever branch runs writes; the other branch is jumped over. */
    if (strcmp(h, "If") == 0 && na == 3) {
        CompileType tt, te;
        if (!infer_type(c, A[1], &tt) || !infer_type(c, A[2], &te)) { c->ok = false; return -1; }
        CompileType rt = (tt == te) ? tt : num_common(tt, te);
        /* An array result would have to be copied into the branch-join register;
         * a MOVE only duplicates the handle, so array-valued branches are not in
         * the M3a subset (see docs/design/compile_state.md). */
        if ((int)rt < 0 || CT_IS_ARRAY(rt)) { c->ok = false; return -1; }
        int rr = alloc_temp(c);                     /* persistent result reg */
        Slot z = { 0 };
        Val cond; if (!emit(c, A[0], &cond)) return -1;
        if (cond.type != CT_BOOL) { c->ok = false; return -1; }
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)cond.reg, 0, z);
        free_if_tmp(c, cond);
        Val th; if (!emit(c, A[1], &th)) return -1;
        if (CT_IS_ARRAY(th.type)) { c->ok = false; return -1; }
        coerce(c, &th, rt);
        ins(c, OP_MOVE, (uint32_t)rr, (uint32_t)th.reg, 0, z);
        free_if_tmp(c, th);
        size_t jmp = c->n; ins(c, OP_JMP, 0, 0, 0, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;   /* else label */
        Val el; if (!emit(c, A[2], &el)) return -1;
        if (CT_IS_ARRAY(el.type)) { c->ok = false; return -1; }
        coerce(c, &el, rt);
        ins(c, OP_MOVE, (uint32_t)rr, (uint32_t)el.reg, 0, z);
        free_if_tmp(c, el);
        if (c->ok) c->code[jmp].b = (uint32_t)c->n;   /* end label */
        out->reg = rr; out->tmp = true; out->type = rt;
        return c->ok ? 1 : -1;
    }

    /* Sum[body, {i, lo, hi}] / Product[...]: integer-counted accumulation loop.
     * The loop variable lives in a scoped register; the accumulator survives the
     * body's temporaries and is the result. */
    if ((strcmp(h, "Sum") == 0 || strcmp(h, "Product") == 0) && na == 2) {
        bool prod = h[0] == 'P';
        LoopSpec s;
        /* Named iterator required — see the matching note in infer_type. */
        if (!loop_spec_parse(A[1], &s) || !s.var || !loop_spec_int_bounds(c, &s)
            || c->nscope >= CTX_MAX_SCOPE) { c->ok = false; return -1; }  /* integer iteration only */
        /* body type with i bound as INT */
        c->scope[c->nscope].name = s.var; c->scope[c->nscope].reg = 0;
        c->scope[c->nscope].type = CT_INT; c->nscope++;
        CompileType T; bool okT = infer_type(c, A[0], &T);
        c->nscope--;
        /* An array accumulator would need a per-iteration copy, not a MOVE. */
        if (!okT || T == CT_BOOL || CT_IS_ARRAY(T)) { c->ok = false; return -1; }
        Slot z = { 0 };
        int racc = alloc_temp(c), rhi = alloc_temp(c), ri = alloc_temp(c);
        Slot iz; iz.i = 0; if (T == CT_INT) iz.i = prod ? 1 : 0; else if (T == CT_REAL) iz.r = prod ? 1.0 : 0.0; else iz.z = prod ? 1.0 : 0.0;
        ins(c, OP_CONST, (uint32_t)racc, 0, 0, iz);
        if (!emit_loop_bound(c, s.lo, ri)) return -1;
        if (!emit_loop_bound(c, s.hi, rhi)) return -1;
        c->scope[c->nscope].name = s.var; c->scope[c->nscope].reg = ri;
        c->scope[c->nscope].type = CT_INT; c->nscope++;
        size_t L = c->n;
        int rc = alloc_temp(c);
        ins(c, s.di > 0 ? OP_LE_I : OP_GE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;                                          /* free the guard temp */
        Val rb; if (!emit(c, A[0], &rb)) { c->nscope--; return -1; }
        if (CT_IS_ARRAY(rb.type)) { c->nscope--; c->ok = false; return -1; }
        coerce(c, &rb, T);
        uint16_t acc = prod ? (T == CT_INT ? OP_MUL_I : T == CT_REAL ? OP_MUL_R : OP_MUL_C)
                            : (T == CT_INT ? OP_ADD_I : T == CT_REAL ? OP_ADD_R : OP_ADD_C);
        ins(c, acc, (uint32_t)racc, (uint32_t)racc, (uint32_t)rb.reg, z);
        free_if_tmp(c, rb);
        /* Force-check the counter step (Sum/Product), same reason as Do: a
         * wrapped counter at the int64 edge never terminates. IF_FORCECHK. */
        Slot step; step.i = s.di; ins_f(c, OP_INC_I, IF_FORCECHK, (uint32_t)ri, 0, 0, step);
        ins(c, OP_JMP, 0, 0, (uint32_t)L, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;               /* loop-exit label */
        c->nscope--;
        c->temp_top -= 2;                                       /* free ri, rhi; racc is result */
        out->reg = racc; out->tmp = true; out->type = T;
        return c->ok ? 1 : -1;
    }

    /* ---- procedural constructs: local variables, mutation, loops ---- */

    /* CompoundExpression[e1,...,en]: run each for side effects, return the last. */
    if (strcmp(h, "CompoundExpression") == 0 && na >= 1) {
        for (size_t i = 0; i + 1 < na; i++) { Val v; if (!emit(c, A[i], &v)) return -1; free_if_tmp(c, v); }
        return emit(c, A[na - 1], out);
    }

    /* With[{v=init,...}, body] / Module[{v or v=init,...}, body]: numeric locals
     * as mutable registers.  The body result is relocated onto the first local's
     * register so the locals free cleanly (LIFO). */
    if ((strcmp(h, "With") == 0 || strcmp(h, "Module") == 0) && na == 2) {
        const Expr* L = A[0];
        if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
            || strcmp(L->data.function.head->data.symbol.name, "List") != 0
            || L->data.function.arg_count == 0
            || (int)(c->nscope + L->data.function.arg_count) > CTX_MAX_SCOPE) { c->ok = false; return -1; }
        size_t nl = L->data.function.arg_count;
        Slot z = { 0 };
        int base_reg = -1, pushed = 0;
        /* Array locals (M3c) get a register in the ARRAY bank on top of their
         * scalar slot, which stays allocated but unused: keeping every local's
         * scalar register in one run is what makes the single temp_top reset at
         * the end correct whatever mix of kinds the locals are. */
        int arr_entry = c->arr_top, arr_regs[16], narr = 0;
        for (size_t i = 0; i < nl; i++) {
            const Expr* spec = L->data.function.args[i];
            const char* vname = NULL; const Expr* init = NULL;
            if (spec->type == EXPR_SYMBOL) vname = spec->data.symbol.name;
            else if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
                     && strcmp(spec->data.function.head->data.symbol.name, "Set") == 0
                     && spec->data.function.arg_count == 2
                     && spec->data.function.args[0]->type == EXPR_SYMBOL) {
                vname = spec->data.function.args[0]->data.symbol.name; init = spec->data.function.args[1];
            } else { c->nscope -= pushed; c->ok = false; return -1; }
            int reg = alloc_temp(c);
            if (base_reg < 0) base_reg = reg;
            CompileType vt = CT_REAL;
            if (init) {
                Val iv; if (!emit(c, init, &iv)) { c->nscope -= pushed; return -1; }
                if (CT_IS_ARRAY(iv.type)) {
                    /* A local may be written through (u[[i]] = ...), so it has to
                     * OWN its array.  An initialiser that is already a temporary
                     * is adopted by leaving its slot allocated for the life of
                     * the scope; anything else — an argument, an enclosing local
                     * — is borrowed and is copied, which is exactly the value
                     * semantics the interpreter gives the same code. */
                    if (narr >= (int)(sizeof arr_regs / sizeof arr_regs[0])) {
                        c->nscope -= pushed; c->ok = false; return -1;
                    }
                    int areg;
                    if (iv.tmp && reg_is_owned_arr(iv.reg)) areg = iv.reg;
                    else {
                        areg = alloc_arr(c);
                        uint16_t f = (uint16_t)(((unsigned)CT_ELEM(iv.type) & 3u) << AF_R_SHIFT);
                        ins_f(c, OP_A_COPY, f, (uint32_t)areg, (uint32_t)iv.reg, 0, z);
                    }
                    arr_regs[narr++] = areg;
                    c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = areg;
                    c->scope[c->nscope].type = iv.type; c->scope[c->nscope].built = iv.built;
                    c->nscope++; pushed++;
                    continue;
                }
                vt = iv.type; ins(c, OP_MOVE, (uint32_t)reg, (uint32_t)iv.reg, 0, z); free_if_tmp(c, iv);
            } else { Slot s; s.r = 0.0; ins(c, OP_CONST, (uint32_t)reg, 0, 0, s); }
            c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = reg; c->scope[c->nscope].type = vt;
            c->nscope++; pushed++;
        }
        Val body; if (!emit(c, A[1], &body)) { c->nscope -= pushed; return -1; }
        c->nscope -= pushed;

        if (CT_IS_ARRAY(body.type)) {
            /* The result has to outlive the frees below, so every array local
             * EXCEPT the one carrying it is released, the bank is wound back to
             * where the scope started, and the value is moved down into the
             * first slot — a handle transfer, not a copy. */
            for (int i = narr - 1; i >= 0; i--)
                if (arr_regs[i] != body.reg) ins(c, OP_ARR_FREE, (uint32_t)arr_regs[i], 0, 0, z);
            c->arr_top = arr_entry;
            int rres = alloc_arr(c);
            if (!reg_is_owned_arr(body.reg)) {
                /* borrowed (an argument array): copy, or the caller's node would
                 * be freed twice — once by them and once as our result */
                uint16_t f = (uint16_t)(((unsigned)CT_ELEM(body.type) & 3u) << AF_R_SHIFT);
                ins_f(c, OP_A_COPY, f, (uint32_t)rres, (uint32_t)body.reg, 0, z);
            } else if (rres != body.reg) {
                ins(c, OP_A_XFER, (uint32_t)rres, (uint32_t)body.reg, 0, z);
            }
            c->temp_top = (base_reg - c->nlocals);   /* no scalar result to keep */
            out->reg = rres; out->tmp = true; out->type = body.type; out->built = body.built;
            return c->ok ? 1 : -1;
        }

        for (int i = narr - 1; i >= 0; i--) ins(c, OP_ARR_FREE, (uint32_t)arr_regs[i], 0, 0, z);
        c->arr_top = arr_entry;
        if (body.reg != base_reg) ins(c, OP_MOVE, (uint32_t)base_reg, (uint32_t)body.reg, 0, z);
        c->temp_top = (base_reg - c->nlocals) + 1;   /* free above base_reg; keep result */
        out->reg = base_reg; out->tmp = true; out->type = body.type;
        return c->ok ? 1 : -1;
    }

    /* Set / AddTo / SubtractFrom / TimesBy on a scoped local (return new value);
     * Increment / Decrement (v++ / v--) return the OLD value. */
    {
        int kind = -1;
        if (strcmp(h, "Set") == 0) kind = 0; else if (strcmp(h, "AddTo") == 0) kind = 1;
        else if (strcmp(h, "SubtractFrom") == 0) kind = 2; else if (strcmp(h, "TimesBy") == 0) kind = 3;
        else if (strcmp(h, "Increment") == 0) kind = 4; else if (strcmp(h, "Decrement") == 0) kind = 5;
        else if (strcmp(h, "DivideBy") == 0) kind = 6;
        /* 4 and 5 are the unary forms (v++ / v--); everything else takes a value. */
        #define IC_UNARY(k) ((k) == 4 || (k) == 5)
        /* Part assignment (M3c): u[[i, j]] = v, and the compound forms.
         *
         * The target must be an array the program OWNS, because the write goes
         * into the buffer in place.  An argument array is borrowed — the caller
         * still holds the node it passed in — so writing through it would mutate
         * a value the caller never offered up, and for a List argument (packed
         * into a temporary at the boundary) the write would silently vanish.
         * Copy it into a local first; that is what the interpreter's value
         * semantics do anyway. */
        if (kind >= 0 && !IC_UNARY(kind) && na == 2 && A[0]->type == EXPR_FUNCTION
            && A[0]->data.function.head->type == EXPR_SYMBOL
            && strcmp(A[0]->data.function.head->data.symbol.name, "Part") == 0
            && A[0]->data.function.arg_count >= 2) {
            const Expr* pt = A[0];
            const Expr* const* S = (const Expr* const*)pt->data.function.args;
            size_t ns = pt->data.function.arg_count - 1;
            Val arr;
            if (!emit(c, S[0], &arr)) return -1;
            if (!CT_IS_ARRAY(arr.type) || !reg_is_owned_arr(arr.reg) || arr.tmp
                || (int)ns > CT_RANK(arr.type)) { c->ok = false; return -1; }
            CompileType elem = CT_ELEM(arr.type);
            Slot z; memset(&z, 0, sizeof z);

            if (part_is_scalar_indexed(c, arr.type, S + 1, ns)) {
                int ridx;
                if (!emit_flat_index(c, arr, S + 1, ns, &ridx)) return -1;
                Val val;
                if (kind == 0) {                         /* plain Set */
                    if (!emit(c, A[1], &val)) return -1;
                } else {                                 /* AddTo / SubtractFrom / TimesBy */
                    int rold = alloc_temp(c);
                    ins(c,a_load_op(elem),
                        (uint32_t)rold, (uint32_t)arr.reg, (uint32_t)ridx, z);
                    Val cur = { rold, true, elem, false }, rhs;
                    if (!emit(c, A[1], &rhs)) return -1;
                    coerce(c, &rhs, elem);
                    uint16_t op = kind == 1 ? (elem == CT_COMPLEX ? OP_ADD_C : OP_ADD_R)
                                : kind == 2 ? (elem == CT_COMPLEX ? OP_SUB_C : OP_SUB_R)
                                : kind == 3 ? (elem == CT_COMPLEX ? OP_MUL_C : OP_MUL_R)
                                            : (elem == CT_COMPLEX ? OP_DIV_C : OP_DIV_R);
                    val = binop(c, op, cur, rhs, elem);
                }
                coerce(c, &val, elem);
                if (!c->ok) return -1;
                ins(c,a_store_op(elem),
                    (uint32_t)arr.reg, (uint32_t)ridx, (uint32_t)val.reg, z);
                /* Set returns the stored value; relocate it onto the index
                 * register so the whole subscript computation is reclaimed. */
                ins(c, OP_MOVE, (uint32_t)ridx, (uint32_t)val.reg, 0, z);
                pop_tmp(c, val);
                c->temp_top = (ridx - c->nlocals) + 1;
                out->reg = ridx; out->tmp = true; out->type = elem;
                return c->ok ? 1 : -1;
            }

            /* General spec: only a plain Set, because the compound forms would
             * have to read a whole slice, combine it and write it back — which
             * is Part-as-an-lvalue on a sub-array, not this. */
            if (kind != 0) { c->ok = false; return -1; }
            int base = -1;
            PartSpec* ps = emit_partspec(c, S + 1, ns, &base);
            if (!ps) return -1;
            if (!ctx_own_partspec(c, ps)) return -1;
            Val val;
            if (!emit(c, A[1], &val)) return -1;
            if (CT_IS_ARRAY(val.type)) ps->rhs_kind = AK_ARR;
            else {
                coerce(c, &val, elem);
                if (!c->ok) return -1;
                /* Same exactness rule as arr_prep's: writing 5 into an int64
                 * buffer must box as the Integer 5, or ndarray_part_set's
                 * mixed-exactness check sees an inexact right-hand side. */
                ps->rhs_kind = (elem == CT_COMPLEX) ? AK_COMPLEX
                             : (elem == CT_INT)     ? AK_INT
                                                    : AK_REAL;
            }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = ps;
            ins_f(c, OP_A_PARTSET, (uint16_t)ns, (uint32_t)arr.reg, (uint32_t)base,
                  (uint32_t)val.reg, ip);
            free_if_tmp(c, val);
            c->temp_top -= (int)ns;
            /* Set's value is the right-hand side; a general one is an array, and
             * an array result would need an owner, so this form is a statement:
             * it reports 0 and is only useful inside a CompoundExpression. */
            Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
            *out = emit_const(c, k0, CT_INT);
            return c->ok ? 1 : -1;
        }

        /* List LHS {a, b, ...} = rhs.  A non-List RHS THREADS (Wolfram Set
         * semantics: {a, b} = c binds a = c, b = c); this is lowered here. A
         * List RHS is DESTRUCTURING and an array RHS is the same shape -- neither
         * is lowered, so bail and let the interpreter (which handles both) take
         * the whole body. Only plain `Set` threads; AddTo/... on a list LHS is
         * not a Wolfram form. */
        if (kind == 0 && na == 2 && A[0]->type == EXPR_FUNCTION
            && A[0]->data.function.head->type == EXPR_SYMBOL
            && strcmp(A[0]->data.function.head->data.symbol.name, "List") == 0) {
            if (A[1]->type == EXPR_FUNCTION && A[1]->data.function.head->type == EXPR_SYMBOL
                && strcmp(A[1]->data.function.head->data.symbol.name, "List") == 0) {
                c->ok = false; return -1;            /* destructuring -> interpreter */
            }
            Val val; if (!emit(c, A[1], &val)) return -1;
            if (CT_IS_ARRAY(val.type)) { free_if_tmp(c, val); c->ok = false; return -1; }
            if (!emit_list_thread(c, A[0], val)) { free_if_tmp(c, val); return -1; }
            *out = val;                              /* Set returns the RHS value */
            return c->ok ? 1 : -1;
        }

        if (kind >= 0) {
            size_t want = IC_UNARY(kind) ? 1 : 2;
            if (na != want || A[0]->type != EXPR_SYMBOL) { c->ok = false; return -1; }
            CompileType vt; int vreg = scope_find(c, A[0]->data.symbol.name, &vt, NULL);
            if (vreg < 0 || vt == CT_BOOL) { c->ok = false; return -1; }   /* mutable numeric locals only */
            /* DivideBy always produces a Real (or Complex), never an Int. */
            if (kind == 6 && vt == CT_INT) { c->ok = false; return -1; }
            Slot z = { 0 };
            if (IC_UNARY(kind)) {
                int old = alloc_temp(c); ins(c, OP_MOVE, (uint32_t)old, (uint32_t)vreg, 0, z);
                if (vt == CT_INT) { Slot s; s.i = (kind == 4) ? 1 : -1; ins(c, OP_INC_I, (uint32_t)vreg, 0, 0, s); }
                else { int one = alloc_temp(c); Slot s; s.r = (kind == 4) ? 1.0 : -1.0; ins(c, OP_CONST, (uint32_t)one, 0, 0, s);
                       ins(c, vt == CT_COMPLEX ? OP_ADD_C : OP_ADD_R, (uint32_t)vreg, (uint32_t)vreg, (uint32_t)one, z); c->temp_top--; }
                out->reg = old; out->tmp = true; out->type = vt; return c->ok ? 1 : -1;
            }
            Val val; if (!emit(c, A[1], &val)) return -1;
            if (CT_IS_ARRAY(val.type)) { c->ok = false; return -1; }
            coerce(c, &val, vt); if (!c->ok) return -1;
            if (kind == 0) ins(c, OP_MOVE, (uint32_t)vreg, (uint32_t)val.reg, 0, z);
            else { uint16_t op = kind == 1 ? (vt == CT_INT ? OP_ADD_I : vt == CT_REAL ? OP_ADD_R : OP_ADD_C)
                              : kind == 2 ? (vt == CT_INT ? OP_SUB_I : vt == CT_REAL ? OP_SUB_R : OP_SUB_C)
                              : kind == 3 ? (vt == CT_INT ? OP_MUL_I : vt == CT_REAL ? OP_MUL_R : OP_MUL_C)
                                          : (vt == CT_COMPLEX ? OP_DIV_C : OP_DIV_R);
                   ins(c, op, (uint32_t)vreg, (uint32_t)vreg, (uint32_t)val.reg, z); }
            free_if_tmp(c, val);
            out->reg = vreg; out->tmp = false; out->type = vt; return c->ok ? 1 : -1;
        }
        #undef IC_UNARY
    }

    /* Do[body, {i, lo, hi}]: counted loop for side effects; returns a dummy 0. */
    /* Do[body, spec1, spec2, ...]: the interpreter nests the iterators with the
     * LAST varying fastest, so rewrite to the nested form and lower that.  Doing
     * it here rather than teaching the loop lowering about several specs keeps
     * one implementation of the loop and makes the multi-iterator form correct
     * by the same argument as the single one. */
    if (strcmp(h, "Do") == 0 && na > 2) {
        Expr** ia = malloc((na - 1) * sizeof(Expr*));
        if (!ia) { c->ok = false; return -1; }
        ia[0] = expr_copy((Expr*)A[0]);
        for (size_t i = 2; i < na; i++) ia[i - 1] = expr_copy((Expr*)A[i]);
        Expr* inner = expr_new_function(expr_new_symbol("Do"), ia, na - 1);
        free(ia);
        if (!inner) { c->ok = false; return -1; }
        Expr* oa[2] = { inner, expr_copy((Expr*)A[1]) };
        Expr* outer = expr_new_function(expr_new_symbol("Do"), oa, 2);
        if (!outer) { expr_free(inner); c->ok = false; return -1; }
        bool r = emit(c, outer, out);
        /* The rewrite is scaffolding, not the user's own tree: a bail inside it
         * would leave the diagnostic pointing into a node freed on the next
         * line, so blame the Do the user actually wrote. */
        if (!r && expr_subtree_of(outer, c->bail_node)) c->bail_node = e;
        expr_free(outer);
        return r;
    }
    if (strcmp(h, "Do") == 0 && na == 2) {
        LoopSpec s;
        if (!loop_spec_parse(A[1], &s) || !loop_spec_int_bounds(c, &s)
            || c->nscope >= CTX_MAX_SCOPE) { c->ok = false; return -1; }
        Slot z = { 0 };
        int rhi = alloc_temp(c), ri = alloc_temp(c);
        if (!emit_loop_bound(c, s.lo, ri)) return -1;
        if (!emit_loop_bound(c, s.hi, rhi)) return -1;
        int pushed = 0;
        if (s.var) {
            c->scope[c->nscope].name = s.var; c->scope[c->nscope].reg = ri;
            c->scope[c->nscope].type = CT_INT; c->nscope++; pushed = 1;
        }
        /* Unit step gets the fused loop instruction: OP_LOOP increments, tests
         * and branches in ONE, where the general shape below spends four
         * instructions per iteration on control alone (test, branch, increment,
         * back-edge).  On `Do[s = s + 1. i, {i, 1, n}]` that is 8 inner
         * instructions down to 5.  OP_LOOP compares `++i < a`, so the bound
         * register holds hi + 1; the entry guard still runs once, because a
         * loop whose range is empty must execute the body zero times and
         * OP_LOOP tests only at the bottom.
         *
         * A non-unit or negative step keeps the general form — OP_LOOP steps by
         * one, and the direction of the test is baked into it. */
        if (s.di == 1) {
            Slot one; memset(&one, 0, sizeof one); one.i = 1;
            int rc = alloc_temp(c), rend = alloc_temp(c);
            /* Entry guard, ONCE: OP_LOOP tests at the bottom, so an empty range
             * would otherwise run the body a first time. */
            ins(c, OP_LE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
            /* OP_LOOP compares `++i < a`, so the bound register holds hi + 1.
             * Force-check this add: if hi is INT64_MAX, hi + 1 overflows, and a
             * wrapped bound (INT64_MIN) would make the loop under-run. Bail to
             * the interpreter instead (see IF_FORCECHK). */
            ins(c, OP_CONST, (uint32_t)rc, 0, 0, one);
            ins_f(c, OP_ADD_I, IF_FORCECHK, (uint32_t)rend, (uint32_t)rhi, (uint32_t)rc, z);
            size_t Lb = c->n;
            Val bod; if (!emit(c, A[0], &bod)) { c->nscope -= pushed; return -1; }
            free_if_tmp(c, bod);
            ins(c, OP_LOOP, (uint32_t)ri, (uint32_t)rend, (uint32_t)Lb, one);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
            c->temp_top -= 2;                                /* rc, rend */
        } else {
            size_t Lp = c->n;
            int rc = alloc_temp(c);
            ins(c, s.di > 0 ? OP_LE_I : OP_GE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z); c->temp_top--;
            Val bod; if (!emit(c, A[0], &bod)) { c->nscope -= pushed; return -1; }
            free_if_tmp(c, bod);
            /* Force-check the counter step: in wrap mode i += di at the int64
             * edge would wrap and the `i <= hi` test never fails -> non-
             * terminating loop. Bail to the interpreter instead (IF_FORCECHK). */
            Slot step; step.i = s.di; ins_f(c, OP_INC_I, IF_FORCECHK, (uint32_t)ri, 0, 0, step);
            ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
        }
        c->nscope -= pushed; c->temp_top -= 2;
        int r0 = alloc_temp(c); Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)r0, 0, 0, s0);
        out->reg = r0; out->tmp = true; out->type = CT_INT; return c->ok ? 1 : -1;
    }

    /* While[test, body] / For[start, test, incr, body]: data-dependent loops. */
    if ((strcmp(h, "While") == 0 && na == 2) || (strcmp(h, "For") == 0 && na == 4)) {
        bool isfor = h[0] == 'F';
        Slot z = { 0 };
        if (isfor) { Val s; if (!emit(c, A[0], &s)) return -1; free_if_tmp(c, s); }
        size_t Lp = c->n;
        Val t; if (!emit(c, isfor ? A[1] : A[0], &t)) return -1;
        if (t.type != CT_BOOL) { c->ok = false; return -1; }
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)t.reg, 0, z); free_if_tmp(c, t);
        Val bod; if (!emit(c, isfor ? A[3] : A[1], &bod)) return -1; free_if_tmp(c, bod);
        if (isfor) { Val ic; if (!emit(c, A[2], &ic)) return -1; free_if_tmp(c, ic); }
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;
        int r0 = alloc_temp(c); Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)r0, 0, 0, s0);
        out->reg = r0; out->tmp = true; out->type = CT_INT; return c->ok ? 1 : -1;
    }

    /* Nest[f, x, n]: apply f n times, feeding each result back in.  The
     * accumulator lives in one persistent register (racc) typed to the
     * fixed-point type; the counted loop mirrors Do.  `f` is any function value
     * fn_resolve accepts — Function[u,body], #-slots, a bare head, Composition,
     * a CompiledFunction — so this one lowering covers all of them. */
    return 0;
}
