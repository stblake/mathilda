/* Mathilda — Compile[]: arithmetic heads (Plus, Times, Subtract, Minus, Divide, Mod, Quotient, Power).
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

int emit_arith(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e; Slot imm;
    if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
        bool mul = h[0] == 'T';
        if (na == 0) { imm.i = mul ? 1 : 0; *out = emit_const(c, imm, CT_INT); return c->ok ? 1 : -1; }
        Val acc; if (!emit(c, A[0], &acc)) return -1;
        for (size_t i = 1; i < na; i++) {
            Val b; if (!emit(c, A[i], &b)) return -1;
            CompileType t = num_common(acc.type, b.type);
            if ((int)t < 0) { c->ok = false; return -1; }
            if (CT_IS_ARRAY(t)) { acc = arr_ew(c, acc, b, t, !mul); if (!c->ok) return -1; continue; }
            coerce(c, &acc, t); coerce(c, &b, t);
            uint16_t op = mul ? (t == CT_INT ? OP_MUL_I : t == CT_REAL ? OP_MUL_R : OP_MUL_C)
                              : (t == CT_INT ? OP_ADD_I : t == CT_REAL ? OP_ADD_R : OP_ADD_C);
            acc = binop(c, op, acc, b, t);
        }
        *out = acc; return c->ok ? 1 : -1;
    }
    if ((strcmp(h, "Subtract") == 0) && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return -1;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return -1; }
        if (CT_IS_ARRAY(t)) {                       /* a - b = a + (-1) b */
            if (CT_IS_ARRAY(b.type)) {
                Val m1 = arr_real_const(c, -1.0);
                b = arr_ew(c, b, m1, b.type, false);
            } else {                                /* scalar subtrahend: negate in a register */
                arr_prep(c, &b, CT_ELEM(t));
                b = unop(c, b.type == CT_COMPLEX ? OP_NEG_C : OP_NEG_R, b, b.type);
            }
            if (!c->ok) return -1;
            *out = arr_ew(c, a, b, t, true); return c->ok ? 1 : -1;
        }
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_INT ? OP_SUB_I : t == CT_REAL ? OP_SUB_R : OP_SUB_C, a, b, t);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Minus") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (a.type == CT_BOOL) { c->ok = false; return -1; }
        if (CT_IS_ARRAY(a.type)) {                  /* -v = (-1) v */
            Val m1 = arr_real_const(c, -1.0);
            *out = arr_ew(c, a, m1, a.type, false); return c->ok ? 1 : -1;
        }
        *out = unop(c, a.type == CT_INT ? OP_NEG_I : a.type == CT_REAL ? OP_NEG_R : OP_NEG_C, a, a.type);
        return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Divide") == 0 && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return -1;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return -1; }
        if (CT_IS_ARRAY(t)) {
            /* a / b = a * b^-1: the ND layer has elementwise power and times,
             * but no divide. */
            CompileType et = CT_ELEM(t);
            if (CT_IS_ARRAY(b.type)) {
                Val em = arr_real_const(c, -1.0);
                Slot z = { 0 };
                b = arr_op(c, OP_V_POW, b, em, b.type, z);
            } else {
                arr_prep(c, &b, et);
                b = unop(c, b.type == CT_COMPLEX ? OP_INV_C : OP_INV_R, b, b.type);
            }
            if (!c->ok) return -1;
            *out = arr_ew(c, a, b, t, false); return c->ok ? 1 : -1;
        }
        if (t < CT_REAL) t = CT_REAL;
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_DIV_C : OP_DIV_R, a, b, t);
        return c->ok ? 1 : -1;
    }
    if ((strcmp(h, "Mod") == 0 || strcmp(h, "Quotient") == 0) && na == 2) {
        CompileType mt, nt;
        if (!infer_type(c, A[0], &mt) || !infer_type(c, A[1], &nt)) { c->ok = false; return -1; }
        /* An array operand goes to the registered binary kernel — see the
         * matching note in the inference. */
        if (CT_IS_ARRAY(mt) || CT_IS_ARRAY(nt)) {
            Val kv;
            if (try_kernel(c, h, A, na, &kv)) { *out = kv; return c->ok ? 1 : -1; }
            c->ok = false; return -1;
        }
        if (mt == CT_INT && nt == CT_INT) {
            Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return -1;
            *out = binop(c, h[0] == 'M' ? OP_MOD_I : OP_QUOT_I, a, b, CT_INT);
            return c->ok ? 1 : -1;
        }
        if (mt <= CT_REAL && nt <= CT_REAL && mt != CT_BOOL && nt != CT_BOOL) {
            if (h[0] == 'M') {          /* the registered machine kernel for Mod */
                Val kv;
                if (try_kernel(c, h, A, na, &kv)) { *out = kv; return c->ok ? 1 : -1; }
                c->ok = false; return -1;
            }
            /* Quotient[a,b] == Floor[a/b], and an INTEGER like the interpreter's. */
            Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return -1;
            coerce(c, &a, CT_REAL); coerce(c, &b, CT_REAL);
            Val q = binop(c, OP_DIV_R, a, b, CT_REAL);
            *out = unop(c, OP_FLOOR_R, q, CT_INT);
            return c->ok ? 1 : -1;
        }
        c->ok = false; return -1;
    }
    if (strcmp(h, "Power") == 0 && na == 2) {
        const Expr* base = A[0]; const Expr* ex = A[1];
        if (any_array_arg(c, A, na)) {
            /* base^exp over arrays: three shapes, all served by the NDArray
             * power fast paths (which promote to a complex dtype where a real
             * base leaves the real axis — caught by the element-type check). */
            Val a, b; if (!emit(c, base, &a) || !emit(c, ex, &b)) return -1;
            CompileType t = num_common(a.type, b.type);
            if ((int)t < 0 || !CT_IS_ARRAY(t)) { c->ok = false; return -1; }
            CompileType et = CT_ELEM(t);
            arr_prep(c, &a, et); arr_prep(c, &b, et);
            Slot z = { 0 };
            *out = arr_op(c, OP_V_POW, a, b, t, z);
            return c->ok ? 1 : -1;
        }
        if (ex->type == EXPR_INTEGER) {
            long long nexp = ex->data.integer;
            Val a; if (!emit(c, base, &a)) return -1;
            /* Integer exponents emit POWI directly rather than through unop(),
             * because the exponent rides in the immediate.  That means this site
             * has to make the tile decision itself — routing it through vec_unop
             * is what keeps a strip-mined `v^3` from reading a tile POINTER as a
             * double. */
            Slot s; memset(&s, 0, sizeof s); s.i = nexp;
            if (c->vector_mode && val_is_tile(a)) {
                CompileType rt = (a.type == CT_COMPLEX) ? CT_COMPLEX : CT_REAL;
                if (rt == CT_REAL) coerce(c, &a, CT_REAL);
                if (!c->ok) return -1;
                *out = vec_unop(c, rt == CT_COMPLEX ? OP_POWI_C : OP_POWI_R, a, rt, s);
                return c->ok ? 1 : -1;
            }
            if (a.type == CT_INT && nexp >= 0) { free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_I, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_INT; return c->ok ? 1 : -1; }
            if (a.type == CT_COMPLEX) { free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_C, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_COMPLEX; return c->ok ? 1 : -1; }
            coerce(c, &a, CT_REAL);
            free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_R, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_REAL; return c->ok ? 1 : -1;
        }
        int64_t rn, rd;
        if (is_rational(ex, &rn, &rd) && rd == 2 && (rn == 1 || rn == -1)) {
            Val a; if (!emit(c, base, &a)) return -1;
            if (a.type == CT_COMPLEX) { a = unop(c, OP_SQRT_C, a, CT_COMPLEX); }
            else { coerce(c, &a, CT_REAL); a = unop(c, OP_SQRT_R, a, CT_REAL); }
            if (rn == -1) a = unop(c, a.type == CT_COMPLEX ? OP_INV_C : OP_INV_R, a, a.type);
            *out = a; return c->ok ? 1 : -1;
        }
        Val a, b; if (!emit(c, base, &a) || !emit(c, ex, &b)) return -1;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return -1; }
        /* Integer base AND integer exponent, the exponent's value known only per
         * call.  `2^n` is an Integer for n >= 0 and a Rational below it, so the
         * only faithful lowering is one that computes exactly and abandons the
         * call on a negative exponent (and on 0^0, which is Indeterminate). */
        if (t == CT_INT && !val_is_tile(a) && !val_is_tile(b)) {
            *out = binop(c, OP_POW_II, a, b, CT_INT);
            return c->ok ? 1 : -1;
        }
        if (t < CT_REAL) t = CT_REAL;
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_POW_C : OP_POW_R, a, b, t);
        return c->ok ? 1 : -1;
    }

    return 0;
}
