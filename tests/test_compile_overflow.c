/* test_compile_overflow.c -- src/checked_int.h against exact GMP arithmetic.
 *
 * checked_int.h has two implementations of every overflow-checked helper: the
 * compiler builtins (__builtin_*_overflow) and a strict-C99 fallback selected by
 * -DCOMPILE_NO_OVERFLOW_BUILTIN. This file is built twice (compile_overflow_tests
 * and compile_overflow_fallback_tests, see CMakeLists.txt), so each build checks
 * one implementation against the same exact reference. Both pass, so they agree.
 *
 * The contract checked for every helper, on every operand pair drawn from a
 * boundary set (0, +-1, +-2, the int32 and int64 edges, the square-root edge
 * 3037000499/3037000500, powers of two):
 *   - the return value is true iff the exact result lies outside int64;
 *   - *out holds the exact result reduced mod 2^64 (the wrapped value), whether
 *     or not it overflowed. Compile[]'s wrap mode relies on that.
 * Both families are covered: ci_* on long long and ci_*_i64 on int64_t.
 *
 * No Mathilda runtime is needed; the test links only GMP. */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>
#include <gmp.h>
#include "checked_int.h"

static int g_fail = 0;
static long g_checks = 0;

/* int64 -> mpz without mpz_set_si, whose `long` is 32 bits on LLP64. */
static void z_set_i64(mpz_t z, int64_t v) {
    uint64_t mag = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_import(z, 1, 1, sizeof mag, 0, 0, &mag);
    if (v < 0) mpz_neg(z, z);
}

/* The exact value z as (fits int64?, z mod 2^64 as int64). */
static bool z_wrap(const mpz_t z, int64_t* wrapped) {
    mpz_t lo, hi, m;
    mpz_init(lo); mpz_init(hi); mpz_init(m);
    z_set_i64(lo, INT64_MIN);
    z_set_i64(hi, INT64_MAX);
    bool fits = mpz_cmp(z, lo) >= 0 && mpz_cmp(z, hi) <= 0;
    mpz_fdiv_r_2exp(m, z, 64);                 /* 0 <= m < 2^64 */
    uint64_t u = 0;
    size_t cnt = 0;
    mpz_export(&u, &cnt, 1, sizeof u, 0, 0, m);
    if (cnt == 0) u = 0;
    *wrapped = (int64_t)u;                     /* two's complement on every target */
    mpz_clear(lo); mpz_clear(hi); mpz_clear(m);
    return fits;
}

static void expect(const char* op, int64_t a, int64_t b, const mpz_t exact,
                   bool ovf, int64_t out) {
    int64_t want;
    bool fits = z_wrap(exact, &want);
    g_checks++;
    if (ovf != !fits || out != want) {
        g_fail++;
        if (g_fail <= 20)
            printf("FAIL %s(%lld, %lld): ovf=%d out=%lld; expected ovf=%d out=%lld\n",
                   op, (long long)a, (long long)b, (int)ovf, (long long)out,
                   (int)!fits, (long long)want);
    }
}

static const int64_t kVals[] = {
    0, 1, -1, 2, -2, 3, -3, 7, -7,
    INT32_MAX, INT32_MIN, (int64_t)INT32_MAX + 1, (int64_t)INT32_MIN - 1,
    INT64_C(4294967296), -INT64_C(4294967296),
    INT64_C(3037000499), INT64_C(3037000500), -INT64_C(3037000499), -INT64_C(3037000500),
    INT64_C(4611686018427387904), -INT64_C(4611686018427387904),   /* 2^62 */
    INT64_MAX, INT64_MAX - 1, INT64_MIN, INT64_MIN + 1,
    INT64_MAX / 2, INT64_MIN / 2, INT64_MAX / 3, INT64_MIN / 3,
};
#define NVALS ((int)(sizeof kVals / sizeof kVals[0]))

static void test_binary(void) {
    mpz_t za, zb, r;
    mpz_init(za); mpz_init(zb); mpz_init(r);
    for (int i = 0; i < NVALS; i++) {
        for (int j = 0; j < NVALS; j++) {
            int64_t a = kVals[i], b = kVals[j];
            z_set_i64(za, a); z_set_i64(zb, b);
            int64_t o; long long ol; bool f;

            mpz_add(r, za, zb);
            f = ci_add_i64(a, b, &o);                     expect("ci_add_i64", a, b, r, f, o);
            f = ci_add((long long)a, (long long)b, &ol);  expect("ci_add", a, b, r, f, (int64_t)ol);

            mpz_sub(r, za, zb);
            f = ci_sub_i64(a, b, &o);                     expect("ci_sub_i64", a, b, r, f, o);
            f = ci_sub((long long)a, (long long)b, &ol);  expect("ci_sub", a, b, r, f, (int64_t)ol);

            mpz_mul(r, za, zb);
            f = ci_mul_i64(a, b, &o);                     expect("ci_mul_i64", a, b, r, f, o);
            f = ci_mul((long long)a, (long long)b, &ol);  expect("ci_mul", a, b, r, f, (int64_t)ol);
        }
    }
    mpz_clear(za); mpz_clear(zb); mpz_clear(r);
}

static void test_unary(void) {
    mpz_t za, r;
    mpz_init(za); mpz_init(r);
    for (int i = 0; i < NVALS; i++) {
        int64_t a = kVals[i];
        z_set_i64(za, a);
        int64_t o; long long ol; bool f;

        mpz_neg(r, za);
        f = ci_neg_i64(a, &o);              expect("ci_neg_i64", a, 0, r, f, o);
        f = ci_neg((long long)a, &ol);      expect("ci_neg", a, 0, r, f, (int64_t)ol);

        mpz_abs(r, za);
        f = ci_abs_i64(a, &o);              expect("ci_abs_i64", a, 0, r, f, o);
        f = ci_abs((long long)a, &ol);      expect("ci_abs", a, 0, r, f, (int64_t)ol);
    }
    mpz_clear(za); mpz_clear(r);
}

/* ci_powi: n >= 0 only (a negative exponent never reaches it). The overflow
 * flag must be exact -- no spurious overflow from squaring past the last
 * round (e.g. 3037000500^1, (-2)^63). */
static void test_powi(void) {
    mpz_t za, r;
    mpz_init(za); mpz_init(r);
    for (int i = 0; i < NVALS; i++) {
        for (int64_t n = 0; n <= 65; n++) {
            int64_t b = kVals[i];
            z_set_i64(za, b);
            mpz_pow_ui(r, za, (unsigned long)n);
            int64_t o; long long ol; bool f;
            f = ci_powi_i64(b, n, &o);                       expect("ci_powi_i64", b, n, r, f, o);
            f = ci_powi((long long)b, (long long)n, &ol);    expect("ci_powi", b, n, r, f, (int64_t)ol);
        }
    }
    mpz_clear(za); mpz_clear(r);
}

int main(void) {
#if defined(__GNUC__) && !defined(COMPILE_NO_OVERFLOW_BUILTIN)
    const char* impl = "compiler builtins";
#else
    const char* impl = "strict-C99 fallback";
#endif
    test_binary();
    test_unary();
    test_powi();
    if (g_fail) {
        printf("checked_int.h (%s): %d of %ld checks FAILED\n", impl, g_fail, g_checks);
        return 1;
    }
    printf("checked_int.h (%s): all %ld checks passed\n", impl, g_checks);
    return 0;
}
