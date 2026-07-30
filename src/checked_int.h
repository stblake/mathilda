/* Mathilda — checked machine-integer arithmetic.
 *
 * The interpreter promotes an int64 that overflows to a GMP bigint, so any fast
 * path that wrapped would answer DIFFERENTLY from the interpreter on the same
 * input rather than merely faster — the one thing a fast path may never do.
 * Every operation that can overflow therefore tests for it and abandons, and the
 * caller re-runs the work through the interpreter, which gives the exact answer.
 * (This is also what the Wolfram Language does by default, under
 * RuntimeOptions -> "CatchMachineIntegerOverflow".)
 *
 * Each helper returns TRUE ON OVERFLOW, matching the __builtin_*_overflow
 * convention, and writes the wrapped result either way — callers must test
 * before using it.
 *
 * The builtins compile to the machine's own overflow flag plus a branch the
 * predictor never takes; the strict-C99 fallback below costs a compare or two
 * more.  Both are exercised: `-DCOMPILE_NO_OVERFLOW_BUILTIN` selects the
 * fallback, and tests/test_compile_overflow.c asserts the two agree at the
 * boundaries.
 *
 * Two families, same semantics:
 *   ci_add / ci_sub / ci_mul / ci_neg / ci_abs / ci_powi   `long long`
 *   ci_*_i64                                               `int64_t`
 * The compiler (src/compile/) uses the first because its Slot register is a
 * `long long`; the NDArray layer uses the second because a buffer is `int64_t*`.
 * They are the same operation and, on every target Mathilda builds for, the same
 * width — but C99 does not guarantee that, and passing an `int64_t*` where a
 * `long long*` is expected is a diagnosable error, so both spellings exist.
 *
 * This header was originally private to src/compile/compile_internal.h; it moved
 * out when the NDArray layer gained exact int64 buffers and needed the identical
 * contract (packed arrays, Phase 3).  compile_internal.h now includes it, so the
 * compiler's ~60 call sites are untouched.
 */
#ifndef MATHILDA_CHECKED_INT_H
#define MATHILDA_CHECKED_INT_H

#include <stdbool.h>
#include <stdint.h>
#include <limits.h>

#if defined(__GNUC__) && !defined(COMPILE_NO_OVERFLOW_BUILTIN)
#define ci_add(a, b, out)  __builtin_add_overflow((a), (b), (out))
#define ci_sub(a, b, out)  __builtin_sub_overflow((a), (b), (out))
#define ci_mul(a, b, out)  __builtin_mul_overflow((a), (b), (out))
#else
/* Strict C99: no GNU builtins, and signed overflow is undefined, so every test
 * has to be made on operands that have NOT yet overflowed. */
static inline bool ci_add_fn(long long a, long long b, long long* out) {
    bool ovf = (b > 0 && a > LLONG_MAX - b) || (b < 0 && a < LLONG_MIN - b);
    *out = ovf ? (long long)((unsigned long long)a + (unsigned long long)b) : a + b;
    return ovf;
}
static inline bool ci_sub_fn(long long a, long long b, long long* out) {
    bool ovf = (b < 0 && a > LLONG_MAX + b) || (b > 0 && a < LLONG_MIN + b);
    *out = ovf ? (long long)((unsigned long long)a - (unsigned long long)b) : a - b;
    return ovf;
}
static inline bool ci_mul_fn(long long a, long long b, long long* out) {
    bool ovf = false;
    if (a != 0 && b != 0) {
        if (a > 0) ovf = (b > 0) ? (a > LLONG_MAX / b) : (b < LLONG_MIN / a);
        else       ovf = (b > 0) ? (a < LLONG_MIN / b)
                                 : (a < LLONG_MAX / b);   /* both negative */
    }
    *out = ovf ? (long long)((unsigned long long)a * (unsigned long long)b) : a * b;
    return ovf;
}
#define ci_add(a, b, out)  ci_add_fn((a), (b), (out))
#define ci_sub(a, b, out)  ci_sub_fn((a), (b), (out))
#define ci_mul(a, b, out)  ci_mul_fn((a), (b), (out))
#endif

/* Negation and absolute value overflow at exactly one point: -LLONG_MIN is not
 * representable.  Easy to forget precisely because it is a single input. */
static inline bool ci_neg(long long a, long long* out) {
    *out = (a == LLONG_MIN) ? a : -a;
    return a == LLONG_MIN;
}
static inline bool ci_abs(long long a, long long* out) {
    *out = (a < 0) ? ((a == LLONG_MIN) ? a : -a) : a;
    return a == LLONG_MIN;
}

/* Binary exponentiation, overflow-checked at every step.  Mirrors the shape of
 * the unchecked `ipow_i` it replaces, so a folded constant and a computed value
 * still agree bit for bit.  `n` must be >= 0: a negative exponent on integers is
 * a Rational in the interpreter and is rejected before reaching here. */
static inline bool ci_powi(long long b, long long n, long long* out) {
    long long r = 1;
    while (n > 0) {
        if (n & 1) { if (ci_mul(r, b, &r)) return true; }
        n >>= 1;
        /* Squaring on the final round is dead work, and squaring a large base
         * is exactly where a spurious overflow would be reported for a result
         * that fits. */
        if (n > 0 && ci_mul(b, b, &b)) return true;
    }
    *out = r;
    return false;
}

/* ------------------------------------------------------------------ *
 *  int64_t spellings                                                  *
 * ------------------------------------------------------------------ *
 * Identical semantics on int64_t operands, for callers whose data is a
 * fixed-width buffer (NDT_INT64 NDArrays) rather than a `long long` register.
 * Written out rather than aliased so they stay correct on a target where
 * `long long` is wider than 64 bits. */
static inline bool ci_add_i64(int64_t a, int64_t b, int64_t* out) {
#if defined(__GNUC__) && !defined(COMPILE_NO_OVERFLOW_BUILTIN)
    return __builtin_add_overflow(a, b, out);
#else
    bool ovf = (b > 0 && a > INT64_MAX - b) || (b < 0 && a < INT64_MIN - b);
    *out = ovf ? (int64_t)((uint64_t)a + (uint64_t)b) : a + b;
    return ovf;
#endif
}
static inline bool ci_sub_i64(int64_t a, int64_t b, int64_t* out) {
#if defined(__GNUC__) && !defined(COMPILE_NO_OVERFLOW_BUILTIN)
    return __builtin_sub_overflow(a, b, out);
#else
    bool ovf = (b < 0 && a > INT64_MAX + b) || (b > 0 && a < INT64_MIN + b);
    *out = ovf ? (int64_t)((uint64_t)a - (uint64_t)b) : a - b;
    return ovf;
#endif
}
static inline bool ci_mul_i64(int64_t a, int64_t b, int64_t* out) {
#if defined(__GNUC__) && !defined(COMPILE_NO_OVERFLOW_BUILTIN)
    return __builtin_mul_overflow(a, b, out);
#else
    bool ovf = false;
    if (a != 0 && b != 0) {
        if (a > 0) ovf = (b > 0) ? (a > INT64_MAX / b) : (b < INT64_MIN / a);
        else       ovf = (b > 0) ? (a < INT64_MIN / b)
                                 : (a < INT64_MAX / b);   /* both negative */
    }
    *out = ovf ? (int64_t)((uint64_t)a * (uint64_t)b) : a * b;
    return ovf;
#endif
}
static inline bool ci_neg_i64(int64_t a, int64_t* out) {
    *out = (a == INT64_MIN) ? a : -a;
    return a == INT64_MIN;
}
static inline bool ci_abs_i64(int64_t a, int64_t* out) {
    *out = (a < 0) ? ((a == INT64_MIN) ? a : -a) : a;
    return a == INT64_MIN;
}

/* True when `v` is exactly representable as a double, i.e. |v| <= 2^53.  The
 * generic ndt_get/ndt_set pair routes through `double`, so an NDT_INT64 buffer
 * may only be read or written through it when this holds; otherwise the exact
 * accessors (ndt_get_i/ndt_set_i) must be used, or the operation must abandon.
 * See src/expr.h's NDType comment. */
static inline bool ci_fits_double(int64_t v) {
    return v <= ((int64_t)1 << 53) && v >= -((int64_t)1 << 53);
}

#endif /* MATHILDA_CHECKED_INT_H */
