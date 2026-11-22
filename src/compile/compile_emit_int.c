/* Mathilda — Compile[]: integer-closed and integer-only heads (Factorial/Gamma/Binomial/Pochhammer/Fibonacci/LucasL; GCD/LCM/IntegerLength/IntegerExponent/PowerMod/EvenQ/OddQ/Divisible/PrimePowerQ).
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
            /* GCD / LCM / BitXor / BitAnd / BitOr: Flat and n-ary, folded left the
             * way the interpreter's own Flat attribute would associate them. */
            uint16_t fold = strcmp(h, "GCD") == 0    ? OP_GCD_I
                          : strcmp(h, "LCM") == 0    ? OP_LCM_I
                          : strcmp(h, "BitXor") == 0 ? OP_BXOR_I
                          : strcmp(h, "BitAnd") == 0 ? OP_BAND_I
                          : strcmp(h, "BitOr") == 0  ? OP_BOR_I
                          : OP__COUNT;                   /* not a fold head */
            if (fold != OP__COUNT) {
                Val acc; if (!emit(c, A[0], &acc)) return -1;
                for (size_t i = 1; i < na; i++) {
                    Val v; if (!emit(c, A[i], &v)) return -1;
                    acc = binop(c, fold, acc, v, CT_INT);
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
            if (strcmp(h, "BitLength") == 0) {
                /* Unary: no base operand (BitLength is base 2 only). */
                Val a; if (!emit(c, A[0], &a)) return -1;
                *out = unop(c, OP_BLEN_I, a, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "BitNot") == 0) {
                /* Unary two's-complement NOT, total on int64. */
                Val a; if (!emit(c, A[0], &a)) return -1;
                *out = unop(c, OP_BNOT_I, a, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "BitSet") == 0) {
                /* Binary n, k -> n with bit k set; OP_BSET_I bails to the
                 * interpreter where there is no int64 answer (k < 0, or a
                 * bignum for n >= 0, k >= 63). */
                Val a, b;
                if (!emit(c, A[0], &a)) return -1;
                if (!emit(c, A[1], &b)) return -1;
                *out = binop(c, OP_BSET_I, a, b, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "BitClear") == 0) {
                /* Binary n, k -> n with bit k cleared; OP_BCLR_I bails to the
                 * interpreter where there is no int64 answer (k < 0, or a
                 * bignum for n < 0, k >= 63). */
                Val a, b;
                if (!emit(c, A[0], &a)) return -1;
                if (!emit(c, A[1], &b)) return -1;
                *out = binop(c, OP_BCLR_I, a, b, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "BitGet") == 0) {
                /* Binary n, k -> bit k of n as 0/1; OP_BGET_I bails to the
                 * interpreter only for k < 0 (BitGet stays unevaluated). */
                Val a, b;
                if (!emit(c, A[0], &a)) return -1;
                if (!emit(c, A[1], &b)) return -1;
                *out = binop(c, OP_BGET_I, a, b, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "BitShiftLeft") == 0 || strcmp(h, "BitShiftRight") == 0) {
                /* n, k -> n * 2^k (Left) or IntegerPart[n / 2^k] (Right); the
                 * one-argument form shifts by one place. Bails to the
                 * interpreter only where a left shift leaves int64. */
                uint16_t op = (h[8] == 'L') ? OP_BSHL_I : OP_BSHR_I;
                Val a, b;
                if (!emit(c, A[0], &a)) return -1;
                if (na == 2) { if (!emit(c, A[1], &b)) return -1; }
                else {
                    memset(&imm, 0, sizeof imm); imm.i = 1;
                    b = emit_const(c, imm, CT_INT);
                    if (!c->ok) return -1;
                }
                *out = binop(c, op, a, b, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "CoprimeQ") == 0) {
                /* CoprimeQ[a, b] is GCD[a, b] == 1; GCD_I's own guard hands
                 * an INT64_MIN operand back to the interpreter. */
                Val a, b;
                if (!emit(c, A[0], &a)) return -1;
                if (!emit(c, A[1], &b)) return -1;
                Val g = binop(c, OP_GCD_I, a, b, CT_INT);
                if (!c->ok) return -1;
                memset(&imm, 0, sizeof imm); imm.i = 1;
                Val one = emit_const(c, imm, CT_INT);
                if (!c->ok) return -1;
                *out = binop(c, OP_EQ_I, g, one, CT_BOOL);
                return c->ok ? 1 : -1;
            }
            if (strcmp(h, "PrimeQ") == 0) {
                /* Unary integer -> boolean; OP_PRIMEQ_I never bails. */
                Val a; if (!emit(c, A[0], &a)) return -1;
                *out = unop(c, OP_PRIMEQ_I, a, CT_BOOL); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "PrimePowerQ") == 0) {
                /* Unary integer -> boolean; OP_PPQ_I never bails. */
                Val a; if (!emit(c, A[0], &a)) return -1;
                *out = unop(c, OP_PPQ_I, a, CT_BOOL); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "DivisorSigma") == 0 || strcmp(h, "JacobiSymbol") == 0) {
                /* Binary integer -> integer; OP_DSIG_I bails to the
                 * interpreter wherever there is no int64 answer, OP_JACOBI_I
                 * never needs to. */
                uint16_t op = (h[0] == 'D') ? OP_DSIG_I : OP_JACOBI_I;
                Val a, b;
                if (!emit(c, A[0], &a)) return -1;
                if (!emit(c, A[1], &b)) return -1;
                *out = binop(c, op, a, b, CT_INT); return c->ok ? 1 : -1;
            }
            if (strcmp(h, "EulerPhi") == 0 || strcmp(h, "MoebiusMu") == 0) {
                /* Unary integer -> integer; bails past the int64 helper. */
                Val a; if (!emit(c, A[0], &a)) return -1;
                *out = unop(c, h[0] == 'E' ? OP_PHI_I : OP_MU_I, a, CT_INT);
                return c->ok ? 1 : -1;
            }
            if (strcmp(h, "PartitionsP") == 0 || strcmp(h, "PartitionsQ") == 0) {
                /* Unary integer -> integer; bails past the int64 table. */
                Val a; if (!emit(c, A[0], &a)) return -1;
                *out = unop(c, h[10] == 'P' ? OP_PARTP_I : OP_PARTQ_I, a, CT_INT);
                return c->ok ? 1 : -1;
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
