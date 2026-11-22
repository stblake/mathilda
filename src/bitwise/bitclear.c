/*
 * bitclear.c
 *
 * BitClear[n, k] -- the integer n with the bit k places from the right (the
 * coefficient of 2^k) set to 0. Integers are taken in two's complement with
 * an unlimited number of digits, so a non-negative n already has every high
 * bit off (BitClear[5, 100] is 5), and clearing a high bit of a negative n
 * subtracts 2^k, growing it to whatever size the bit needs (BitClear[-1, 0]
 * is -2, BitClear[-1, 64] is -1 - 2^64). GMP's mpz_clrbit has exactly that
 * semantics for either sign.
 *
 * Attributes Listable, Protected. The evaluator threads plain Lists before the
 * builtin runs, so it only ever sees one (n, k) pair.
 *
 * Everything else is left unevaluated SILENTLY, with no message, including
 * the wrong argument count: k < 0, a bignum k, a non-integer n or k (2.5, x),
 * BitClear[5], BitClear[5, 1, 2] -- the same contract as BitSet and BitGet.
 * There are no symbolic rules.
 *
 * Surfaces:
 *   - packed / NDArray: the int64 binary kernel ndk_BitClear_ii and its
 *     whole-buffer twins in src/ndinteger.c answer before this builtin runs.
 *   - Compile[]: OP_BCLR_I in src/compile/compile_emit_int.c.
 * All three, and the machine path below, call ndint_bitclear_i64, so every
 * surface agrees bit for bit; where it declines (a bignum answer) each
 * surface falls back to this builtin's GMP path.
 */

#include "bitwise.h"
#include "ndinteger.h"   /* ndint_bitclear_i64: the shared machine helper */

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

/* The largest k mpz_clrbit can honour on a negative operand, where clearing
 * the bit extends the number to k + 1 bits: an mpz holds at most INT_MAX
 * limbs, and mp_bitcnt_t is an unsigned long, which is only 32 bits on LLP64
 * Windows. A bit beyond either has no representation at all, so the call is
 * left unevaluated rather than letting GMP abort the process. */
static bool bitclear_k_representable(int64_t k) {
    if ((uint64_t)k > (uint64_t)ULONG_MAX) return false;
    return k < (int64_t)INT_MAX * GMP_NUMB_BITS;
}

Expr* builtin_bitclear(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* n = res->data.function.args[0];
    Expr* kx = res->data.function.args[1];
    if (kx->type != EXPR_INTEGER || kx->data.integer < 0) return NULL;
    int64_t k = kx->data.integer;

    if (n->type == EXPR_INTEGER) {
        int64_t r;
        if (ndint_bitclear_i64(n->data.integer, k, &r)) return expr_new_integer(r);
    } else if (n->type != EXPR_BIGINT) {
        return NULL;
    }

    /* A bignum n, or an int64 n < 0 whose bit k >= 63 makes a bignum. A
     * non-negative bignum with no bit at k is its own answer, whatever the
     * size of k, so it needs no representable-k check. */
    if (n->type == EXPR_BIGINT && mpz_sgn(n->data.bigint) >= 0 &&
        (uint64_t)k >= (uint64_t)mpz_sizeinbase(n->data.bigint, 2))
        return expr_copy(n);
    if (!bitclear_k_representable(k)) return NULL;

    /* Clear the bit in a fresh node; the operand may be shared and is never
     * written. Clearing a bit of a positive bignum can land it back in int64
     * range (BitClear[2^63, 63] is 0), which expr_bigint_normalize demotes. */
    Expr* r;
    if (n->type == EXPR_INTEGER) {
        r = expr_new_bigint_from_int64(n->data.integer);
    } else {
        r = expr_new_bigint_from_int64(0);
        mpz_set(r->data.bigint, n->data.bigint);
    }
    mpz_clrbit(r->data.bigint, (mp_bitcnt_t)k);
    return expr_bigint_normalize(r);
}
