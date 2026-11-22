#ifndef MATHILDA_NDINTEGER_H
#define MATHILDA_NDINTEGER_H

/* ---------------------------------------------------------------------------
 * ndinteger.c — exact-integer NDArray kernels for the number-theoretic heads.
 *
 * ndkernels.c owns the libc-expressible elementary and special functions: every
 * kernel there computes in `double` and its exactness question is only ever
 * "does this real answer round-trip". These are the opposite kind. GCD, LCM,
 * EulerPhi, MoebiusMu, DivisorSigma, IntegerLength, PowerMod, Prime and
 * IntegerDigits are defined on Z, their answers are exact Integers, and a
 * `double` intermediate is not an approximation of the answer but a different
 * answer. So they are written in int64 end to end, through the ci_*_i64
 * overflow-checked helpers, and ABANDON the whole array the moment a result
 * would not fit — the List path then re-runs the call and GMP answers exactly.
 * Never wrap, never round: the same contract as the narrowing kernels in
 * ndkernels.c and Total's int64 accumulate.
 *
 * WHY THIS EXISTS. `Mod[Range[10^6] 7919, 1000]` is how an integer array is
 * made, and until it kept its buffer every consumer downstream paid. The eighth
 * round fixed Mod and Quotient and left the rest of the integer domain with no
 * kernel at all: `GCD[cv, 1234]` and its seven siblings answered a visible
 * NDArray with the unevaluated call and a packed List with one Expr per
 * element. The natural next batch, and the same reasoning — an integer pipeline
 * is only as packed as its least-aware link.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

/* Register the unary/binary integer kernels (GCD, LCM, DivisorSigma, EulerPhi,
 * MoebiusMu, IntegerLength) on their symbols. Called from core_init. */
void ndinteger_init(void);

/* ---- head-level buffer paths ------------------------------------------- *
 * Three heads cannot be expressed as an element kernel and get their own entry
 * point, called from the head's builtin before it does anything else. Each
 * returns NULL for any shape it does not handle, leaving the builtin's ordinary
 * path untouched.
 *
 *   Prime      per-element Prime[n] is a Meissel-style prime count each time.
 *              Over an array the answer is a SIEVE: one pass to the largest
 *              requested index's upper bound, then a gather. That is a
 *              different algorithm, not a different loop, which is why it
 *              cannot be a kernel.
 *   PowerMod   ternary — PowerMod[b, e, m] — and NDUnary/NDBinaryKernel cannot
 *              express three arguments.
 *   IntegerDigits  the result is RAGGED (a list per element, of differing
 *              lengths), so no buffer holds it. The win is on the INPUT side:
 *              read the buffer directly instead of materialising 10^6 Integers
 *              only to read each one back out again.
 */
/* ---- DivisorSigma's machine helper -------------------------------------- *
 * sigma_k(n) in int64, shared by every machine surface of DivisorSigma (the
 * NDArray kernel and its sieve, Compile's OP_DSIG_I, the scalar builtin), so
 * all four agree by construction. Both return false -- DECLINE, never wrap --
 * for k < 0, for any overflow, and (ndint_divisor_sigma_i64) for n == 0 or
 * |n| past the trial-division ceiling; the caller then answers through GMP.
 * ndint_sigma_from_factors takes a factorisation the caller already has
 * (distinct primes p[0..np) with exponents e[0..np)). */
bool ndint_divisor_sigma_i64(int64_t k, int64_t n, int64_t* out);
bool ndint_sigma_from_factors(int64_t k, const int64_t* p, const int* e, int np,
                              int64_t* out);

/* EulerPhi / MoebiusMu / JacobiSymbol machine helpers, shared the same way
 * by the element kernels, Compile (OP_PHI_I / OP_MU_I / OP_JACOBI_I) and the
 * scalar builtins. The first two decline (false) past the trial-division
 * ceiling or for INT64_MIN; EulerPhi[0] and MoebiusMu[0] are 0. The Kronecker
 * symbol is total on int64. */
bool ndint_euler_phi_i64(int64_t n, int64_t* out);
/* Factor |n| (n != 0, |n| <= 10^12) by int64 trial division into distinct
 * primes p[0..*np) with exponents e[]; arrays need NDINT_MAX_FACTORS slots.
 * false for 0, INT64_MIN or past the ceiling. */
#define NDINT_MAX_FACTORS 64
bool ndint_factor_i64(int64_t n, int64_t* p, int* e, int* np);
bool ndint_moebius_mu_i64(int64_t n, int64_t* out);
int64_t ndint_kronecker_i64(int64_t a, int64_t n);

/* BitSet[n, k] in int64: n with bit k set, two's complement. Shared by the
 * scalar builtin (src/bitwise/bitset.c), the NDArray kernels and Compile's
 * OP_BSET_I, so every machine surface agrees bit for bit. Declines (false)
 * for k < 0 (BitSet stays unevaluated) and for n >= 0 with k >= 63 (the
 * answer is a bignum). A negative n already has every bit from 63 up set, so
 * BitSet[n, k] is n there. Inline: it sits in the VM's dispatch loop. */
static inline bool ndint_bitset_i64(int64_t n, int64_t k, int64_t* out) {
    if (k < 0) return false;
    if (k <= 62) { *out = (int64_t)((uint64_t)n | (UINT64_C(1) << k)); return true; }
    if (n < 0) { *out = n; return true; }
    return false;
}

/* BitClear[n, k] in int64: n with bit k cleared, two's complement. Shared by
 * the scalar builtin (src/bitwise/bitclear.c), the NDArray kernels and
 * Compile's OP_BCLR_I, so every machine surface agrees bit for bit. Declines
 * (false) for k < 0 (BitClear stays unevaluated) and for n < 0 with k >= 63
 * (clearing one of a negative n's high bits subtracts 2^k, a bignum answer).
 * A non-negative n has every bit from 63 up clear already, so BitClear[n, k]
 * is n there. Inline: it sits in the VM's dispatch loop. */
static inline bool ndint_bitclear_i64(int64_t n, int64_t k, int64_t* out) {
    if (k < 0) return false;
    if (k <= 62) { *out = (int64_t)((uint64_t)n & ~(UINT64_C(1) << k)); return true; }
    if (n >= 0) { *out = n; return true; }
    return false;
}

/* BitGet[n, k] in int64: bit k of n as 0 or 1, two's complement. Shared by
 * the scalar builtin (src/bitwise/bitget.c), the NDArray kernels and
 * Compile's OP_BGET_I. Declines (false) only for k < 0 (BitGet stays
 * unevaluated). Clamping k to 63 is exact: bit 63 is the sign bit, and in
 * unlimited two's complement every higher bit equals it. */
static inline bool ndint_bitget_i64(int64_t n, int64_t k, int64_t* out) {
    if (k < 0) return false;
    *out = (int64_t)(((uint64_t)n >> (k < 63 ? k : 63)) & 1);
    return true;
}

/* BitShiftLeft[n, k] in int64: n * 2^k, where a negative k shifts RIGHT and
 * truncates toward zero (IntegerPart[n / 2^-k], so BitShiftLeft[-5, -1] is
 * -2, not the arithmetic shift's -3). Shared by the scalar builtins
 * (src/bitwise/bitshiftleft.c), the NDArray kernels and Compile's
 * OP_BSHL_I / OP_BSHR_I, so every machine surface agrees bit for bit.
 * Declines (false) only when a left shift has no int64 answer (the result
 * is a bignum); every right shift is an int64. Shifting 0 is always 0.
 *
 * The right shift adds 2^s - 1 to a negative n before the arithmetic shift,
 * which rounds toward zero instead of toward -infinity. It cannot overflow:
 * the bias is added only to a negative n, and is at most 2^63 - 1. */
static inline bool ndint_bitshiftleft_i64(int64_t n, int64_t k, int64_t* out) {
    if (k >= 0) {
        if (k > 62) {
            /* 2^63 is no int64; only 0, and -1 with k == 63 (-2^63), fit. */
            if (n == 0 || (n == -1 && k == 63)) { *out = n == 0 ? 0 : INT64_MIN; return true; }
            return false;
        }
        int64_t r = (int64_t)((uint64_t)n << k);
        if ((r >> k) != n) return false;
        *out = r;
        return true;
    }
    uint64_t s = (uint64_t)0 - (uint64_t)k;   /* -k, exact even for INT64_MIN */
    if (s > 63) { *out = 0; return true; }
    int64_t bias = (n >> 63) & (int64_t)((UINT64_C(1) << s) - 1);
    *out = (n + bias) >> s;
    return true;
}

/* BitShiftRight[n, k] = BitShiftLeft[n, -k]. -INT64_MIN is not an int64, but
 * a left shift by 2^63 has an int64 answer only for n == 0. */
static inline bool ndint_bitshiftright_i64(int64_t n, int64_t k, int64_t* out) {
    if (k == INT64_MIN) {
        if (n != 0) return false;
        *out = 0;
        return true;
    }
    return ndint_bitshiftleft_i64(n, -k, out);
}

Expr* ndint_prime(Expr* res);          /* Prime[arr] */
Expr* ndint_powermod(Expr* res);       /* PowerMod[arr, e, m] */
Expr* ndint_integerdigits(Expr* res);  /* IntegerDigits[arr] */

/* ---- elementwise sign predicates ---------------------------------------- *
 * Positive / Negative / NonNegative / NonPositive. Listable, answering with a
 * List of True/False. Both sides now stay on the buffer: the INPUT is read
 * straight off it and the OUTPUT is a packed one-byte-per-element NDT_BOOL array
 * (performance.md §13 gap C.1, now closed), inheriting the input's presentation.
 * No per-element Expr and no evaluator round-trip on either side.
 *
 * `which` selects the comparison; see NDSignPred. Returns NULL unless the sole
 * argument is a rank-1 real or int64 NDArray (a complex or Indeterminate element
 * has no order, so those degrade to the List path element for element). */
typedef enum {
    NDSP_POSITIVE, NDSP_NEGATIVE, NDSP_NONNEGATIVE, NDSP_NONPOSITIVE
} NDSignPred;
Expr* ndint_sign_predicate(Expr* res, NDSignPred which);

/* PrimePowerQ over a rank-1 int64 buffer: one mth_prime_power_base_u64 per
 * element, answered as a packed NDT_BOOL List (the sign predicates' output
 * contract). NULL for a non-int64 dtype, rank > 1, an empty buffer, or an
 * option argument -- the caller degrades through ndarray_delist_and_reeval. */
Expr* ndint_prime_power_predicate(Expr* res);

/* PrimeQ[int64 buffer]: a packed bool List (a sieve for dense small values,
 * deterministic Miller-Rabin otherwise). NULL for any other shape. */
Expr* ndint_primeq(Expr* res);

/* CoprimeQ[a, b] / Divisible[n, m] over int64 buffers (one side may be an
 * exact Integer), answering with a packed NDT_BOOL List, as the sign
 * predicates and PrimePowerQ do. NULL for any shape it does not handle -- including a zero
 * divisor, whose ::divz message the element path must emit. */
typedef enum { NDBP_COPRIME, NDBP_DIVISIBLE } NDBinPred;
Expr* ndint_binary_predicate(Expr* res, NDBinPred which);

#endif /* MATHILDA_NDINTEGER_H */
