/* Mathilda — Compile[]: integer-closed and integer-only heads (Factorial/Gamma/Binomial/Pochhammer/Fibonacci/LucasL; GCD/LCM/IntegerLength/IntegerExponent/PowerMod/EvenQ/OddQ/Divisible).
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

int emit_int(Ctx* c, const char* h, const Expr* e, Expr** A, size_t na, Val* out) {
    (void)e; Slot imm;
    /* Integer-closed heads, before the generic kernel dispatch — see INT_CLOSED.
     * Declining here (rather than bailing) leaves every non-integer call to the
     * ordinary real kernel it already used. */
    {
        const IntClosed* ic = int_closed_head(h, na);
        if (ic && all_args_int(c, A, na)) {
            Val a; if (!emit(c, A[0], &a)) return -1;
            if (na == 1) { *out = unop(c, ic->op, a, CT_INT); return c->ok ? 1 : -1; }
            Val b; if (!emit(c, A[1], &b)) return -1;
            *out = binop(c, ic->op, a, b, CT_INT);
            return c->ok ? 1 : -1;
        }
    }
    /* Integer-only heads — see int_only_head. */
    {
        bool pred = false;
        if (int_only_head(h, na, &pred) && all_args_int(c, A, na)) {
            /* GCD / LCM: Flat and n-ary, folded left the way the interpreter's
             * own Flat attribute would associate them. */
            if (strcmp(h, "GCD") == 0 || strcmp(h, "LCM") == 0) {
                uint16_t op = (h[0] == 'G') ? OP_GCD_I : OP_LCM_I;
                Val acc; if (!emit(c, A[0], &acc)) return -1;
                for (size_t i = 1; i < na; i++) {
                    Val v; if (!emit(c, A[i], &v)) return -1;
                    acc = binop(c, op, acc, v, CT_INT);
                    if (!c->ok) return -1;
                }
                *out = acc; return c->ok ? 1 : -1;
            }
            if (strcmp(h, "IntegerLength") == 0 || strcmp(h, "IntegerExponent") == 0) {
                uint16_t op = (h[7] == 'L') ? OP_ILEN_I : OP_IEXP_I;
                Val a; if (!emit(c, A[0], &a)) return -1;
                Val b;
                if (na == 2) { if (!emit(c, A[1], &b)) return -1; }
                else {
                    /* The default base is 10 for both. */
                    memset(&imm, 0, sizeof imm); imm.i = 10;
                    b = emit_const(c, imm, CT_INT);
                    if (!c->ok) return -1;
                }
                *out = binop(c, op, a, b, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "PowerMod") == 0) {
                /* Ternary: three CONSECUTIVE registers, the K_NARY shape. Same
                 * allocate-then-MOVE dance as the n-ary kernel path, and for the
                 * same reason — each argument's own lowering needs temps above
                 * the run, so the run has to be reserved first. */
                int base = alloc_temp(c);
                (void)alloc_temp(c); (void)alloc_temp(c);
                int after = c->temp_top;
                Slot z; memset(&z, 0, sizeof z);
                for (size_t i = 0; i < 3 && c->ok; i++) {
                    Val v;
                    if (!emit(c, A[i], &v)) return -1;
                    ins(c, OP_MOVE, (uint32_t)(base + (int)i), (uint32_t)v.reg, 0, z);
                    c->temp_top = after;
                }
                if (!c->ok) return -1;
                c->temp_top = (base - c->nlocals) + 1;
                ins_f(c, OP_POWMOD_I, 3, (uint32_t)base, (uint32_t)base, 0, z);
                out->reg = base; out->tmp = true; out->type = CT_INT;
                return c->ok ? 1 : -1;
            }
            /* EvenQ / OddQ / Divisible are Mod-and-compare, so they need no
             * opcode of their own — and they inherit MOD_I's guard, which hands
             * `Divisible[n, 0]` to the interpreter rather than dividing by it. */
            {
                Val a; if (!emit(c, A[0], &a)) return -1;
                Val b;
                if (na == 2) { if (!emit(c, A[1], &b)) return -1; }
                else {
                    memset(&imm, 0, sizeof imm); imm.i = 2;
                    b = emit_const(c, imm, CT_INT);
                    if (!c->ok) return -1;
                }
                Val m = binop(c, OP_MOD_I, a, b, CT_INT);
                if (!c->ok) return -1;
                memset(&imm, 0, sizeof imm); imm.i = (strcmp(h, "OddQ") == 0) ? 1 : 0;
                Val k = emit_const(c, imm, CT_INT);
                if (!c->ok) return -1;
                *out = binop(c, OP_EQ_I, m, k, CT_BOOL);
                return c->ok ? 1 : -1;
            }
        }
    }
    return 0;
}
