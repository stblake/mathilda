/* Mathilda — Compile[]: comparison and boolean-logic head lowering.
 *
 * The comparison heads (Less..Unequal), Order, the n-ary And/Or/Xor and Not.
 * One of emit_node's per-category dispatch functions: returns 1 if `h` is one
 * of these heads (result in *out, valid iff c->ok), -1 if it is one but could
 * not be lowered, and 0 if `h` belongs to another category. */
#include "compile_emit.h"
#include "../arithmetic.h"
#include "../symtab.h"
#include <string.h>

int emit_logic(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e;
    /* comparisons (2-arg) -> Bool */
    struct { const char* name; uint16_t oi, orr, oc; } CMP[] = {
        { "Less", OP_LT_I, OP_LT_R, 0 }, { "LessEqual", OP_LE_I, OP_LE_R, 0 },
        { "Greater", OP_GT_I, OP_GT_R, 0 }, { "GreaterEqual", OP_GE_I, OP_GE_R, 0 },
        { "Equal", OP_EQ_I, OP_EQ_R, OP_EQ_C }, { "Unequal", OP_NE_I, OP_NE_R, OP_NE_C },
    };
    for (size_t k = 0; k < sizeof(CMP) / sizeof(CMP[0]); k++) {
        if (strcmp(h, CMP[k].name) == 0 && na == 2) {
            Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return -1;
            CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return -1; }
            coerce(c, &a, t); coerce(c, &b, t);
            uint16_t op = t == CT_COMPLEX ? CMP[k].oc : t == CT_INT ? CMP[k].oi : CMP[k].orr;
            if (!op) { c->ok = false; return -1; }   /* ordering of complex */
            *out = binop(c, op, a, b, CT_BOOL); return c->ok ? 1 : -1;
        }
    }
    /* Order[a,b] -> Sign[b - a]: +1 when a sorts before b (a < b), -1 when after,
     * 0 when equal.  Built from the existing SUB + SIGN opcodes; SIGN_* writes
     * the integer slot, so the result is the Integer {-1,0,1} that the
     * interpreter's Order returns (result-HEAD parity).  Mixed int/real args
     * unify to real; complex/array/bool args are rejected above in infer_type
     * and here, sending the call to the interpreter. */
    if (strcmp(h, "Order") == 0 && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return -1;
        CompileType t = num_common(a.type, b.type);
        if ((int)t < 0 || t == CT_COMPLEX || CT_IS_ARRAY(t)) { c->ok = false; return -1; }
        coerce(c, &a, t); coerce(c, &b, t);
        Val d = binop(c, t == CT_INT ? OP_SUB_I : OP_SUB_R, b, a, t);   /* b - a */
        *out = unop(c, t == CT_INT ? OP_SIGN_I : OP_SIGN_R, d, CT_INT);
        return c->ok ? 1 : -1;
    }
    /* boolean logic */
    if ((strcmp(h, "And") == 0 || strcmp(h, "Or") == 0 || strcmp(h, "Xor") == 0) && na >= 1) {
        uint16_t op = h[0] == 'A' ? OP_AND : h[0] == 'O' ? OP_OR : OP_XOR;
        Val acc; if (!emit(c, A[0], &acc)) return -1;
        if (acc.type != CT_BOOL) { c->ok = false; return -1; }
        for (size_t i = 1; i < na; i++) {
            Val b; if (!emit(c, A[i], &b)) return -1;
            if (b.type != CT_BOOL) { c->ok = false; return -1; }
            acc = binop(c, op, acc, b, CT_BOOL);
        }
        *out = acc; return c->ok ? 1 : -1;
    }
    if (strcmp(h, "Not") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return -1;
        if (a.type != CT_BOOL) { c->ok = false; return -1; }
        *out = unop(c, OP_NOT, a, CT_BOOL); return c->ok ? 1 : -1;
    }
    return 0;
}
