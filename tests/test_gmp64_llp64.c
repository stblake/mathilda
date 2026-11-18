/* Mathilda -- the int64/uint64 GMP/MPFR helpers (src/expr_gmp64.c) on their
 * LLP64 branch, checked against exact GMP on an LP64 host.
 *
 * On native Windows `long` is 32 bits, so GMP's `_si`/`_ui` calls truncate an
 * int64 above 2^31 (numeric-tower review F2).  The helpers take the `_si`/`_ui`
 * call when the value fits a `long` and an exact mpz otherwise.  That second
 * branch never runs where `long` is 64 bits, so this suite includes the helper
 * source with MTH_TEST_LLP64_PATH, which selects it with a 32-bit `long` range,
 * and checks every helper over the values on both sides of +-2^31, 2^32 and
 * +-2^63 against a reference built from the decimal string.
 *
 * Standalone: links GMP and MPFR only, not the kernel (the helpers are the
 * only definitions; expr.h's other declarations are never called).
 */
#define MTH_TEST_LLP64_PATH 1
#include "../src/expr_gmp64.c"

#include <inttypes.h>
#include <stdio.h>
#include <string.h>

static int g_fail = 0, g_checks = 0;

#define CHECK(cond, ...) do { g_checks++; if (!(cond)) { g_fail++; \
    printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

/* Exact references: through the decimal string, never through `long`. */
static void ref_i64(mpz_t z, int64_t v) {
    char buf[32]; snprintf(buf, sizeof buf, "%" PRId64, v); mpz_set_str(z, buf, 10);
}
static void ref_u64(mpz_t z, uint64_t v) {
    char buf[32]; snprintf(buf, sizeof buf, "%" PRIu64, v); mpz_set_str(z, buf, 10);
}

static const int64_t I64[] = {
    0, 1, -1, 7, -7, 2147483647, -2147483647, INT64_C(-2147483648), INT64_C(2147483648),
    INT64_C(-2147483649), INT64_C(4294967295), INT64_C(4294967296), INT64_C(-4294967296),
    INT64_C(4294967297), INT64_C(1099511627776), INT64_C(-1099511627777),
    INT64_C(9007199254740993), INT64_MAX, INT64_MAX - 1, INT64_MIN, INT64_MIN + 1,
};
static const uint64_t U64[] = {
    0, 1, 7, 2147483647u, 2147483648u, 4294967295u, UINT64_C(4294967296),
    UINT64_C(4294967297), UINT64_C(1099511627791), UINT64_C(18446744073709551557),
    UINT64_MAX, UINT64_MAX - 1, (uint64_t)INT64_MAX + 1,
};
#define N_I (sizeof I64 / sizeof I64[0])
#define N_U (sizeof U64 / sizeof U64[0])

static void test_mpz(void) {
    mpz_t a, r, x, y, t;
    mpz_inits(a, r, x, y, t, NULL);
    /* A multi-limb operand and a small one, both signs. */
    mpz_set_str(a, "-123456789012345678901234567890123", 10);
    for (int s = 0; s < 4; s++) {
        if (s == 1) mpz_neg(a, a);
        if (s == 2) mpz_set_si(a, 3);
        if (s == 3) mpz_set_si(a, -3);
        for (size_t i = 0; i < N_I; i++) {
            int64_t v = I64[i];
            ref_i64(y, v);
            mth_mpz_set_i64(x, v);
            CHECK(mpz_cmp(x, y) == 0, "set_i64(%" PRId64 ")", v);
            mpz_clear(x); mth_mpz_init_set_i64(x, v);
            CHECK(mpz_cmp(x, y) == 0, "init_set_i64(%" PRId64 ")", v);
            mth_mpz_mul_i64(r, a, v); mpz_mul(t, a, y);
            CHECK(mpz_cmp(r, t) == 0, "mul_i64(%" PRId64 ")", v);
            int c1 = mth_mpz_cmp_i64(a, v), c2 = mpz_cmp(a, y);
            CHECK((c1 > 0) == (c2 > 0) && (c1 < 0) == (c2 < 0), "cmp_i64(%" PRId64 ")", v);
            CHECK(mth_mpz_cmp_i64(y, v) == 0, "cmp_i64 equal (%" PRId64 ")", v);
            for (size_t j = 0; j < N_U; j++) {
                uint64_t d = U64[j];
                if (d == 0) continue;
                mpq_t q, qr; mpq_init(q); mpq_init(qr);
                mth_mpq_set_i64(q, v, d);
                mpz_set(mpq_numref(qr), y); ref_u64(mpq_denref(qr), v ? d : 1);  /* 0 is 0/1, as GMP */
                CHECK(mpz_cmp(mpq_numref(q), mpq_numref(qr)) == 0
                      && mpz_cmp(mpq_denref(q), mpq_denref(qr)) == 0,
                      "mpq_set_i64(%" PRId64 ", %" PRIu64 ")", v, d);
                mpq_clear(q); mpq_clear(qr);
            }
        }
        for (size_t i = 0; i < N_U; i++) {
            uint64_t v = U64[i];
            ref_u64(y, v);
            mth_mpz_set_u64(x, v);
            CHECK(mpz_cmp(x, y) == 0, "set_u64(%" PRIu64 ")", v);
            mpz_clear(x); mth_mpz_init_set_u64(x, v);
            CHECK(mpz_cmp(x, y) == 0, "init_set_u64(%" PRIu64 ")", v);
            CHECK(mth_mpz_get_u64(y) == v, "get_u64(%" PRIu64 ")", v);
            mpz_neg(t, y);
            CHECK(mth_mpz_get_u64(t) == v, "get_u64(-%" PRIu64 ") is |z|", v);
            mth_mpz_mul_u64(r, a, v); mpz_mul(t, a, y);
            CHECK(mpz_cmp(r, t) == 0, "mul_u64(%" PRIu64 ")", v);
            mth_mpz_add_u64(r, a, v); mpz_add(t, a, y);
            CHECK(mpz_cmp(r, t) == 0, "add_u64(%" PRIu64 ")", v);
            mth_mpz_sub_u64(r, a, v); mpz_sub(t, a, y);
            CHECK(mpz_cmp(r, t) == 0, "sub_u64(%" PRIu64 ")", v);
            int c1 = mth_mpz_cmp_u64(a, v), c2 = mpz_cmp(a, y);
            CHECK((c1 > 0) == (c2 > 0) && (c1 < 0) == (c2 < 0), "cmp_u64(%" PRIu64 ")", v);
            CHECK(mth_u64_fits_ulong(v) == (v <= UINT32_MAX), "fits_ulong(%" PRIu64 ")", v);
            if (v == 0) continue;
            /* fdiv: a mod v in [0, v). */
            mpz_fdiv_r(t, a, y);
            mpz_t m; mpz_init(m); ref_u64(m, mth_mpz_fdiv_u64(a, v));
            CHECK(mpz_cmp(m, t) == 0, "fdiv_u64(%" PRIu64 ")", v);
            /* divisible / divexact on a known multiple. */
            mpz_mul(t, a, y);
            CHECK(mth_mpz_divisible_u64_p(t, v), "divisible_u64(%" PRIu64 ")", v);
            mth_mpz_divexact_u64(m, t, v);
            CHECK(mpz_cmp(m, a) == 0, "divexact_u64(%" PRIu64 ")", v);
            if (v > 1) {
                mpz_add_ui(t, t, 1);
                CHECK(!mth_mpz_divisible_u64_p(t, v), "not divisible_u64(%" PRIu64 ")", v);
            }
            mpz_clear(m);
            mpq_t q, qr; mpq_init(q); mpq_init(qr);
            mth_mpq_set_u64(q, v, v == UINT64_MAX ? 3 : v + 1);
            ref_u64(mpq_numref(qr), v); ref_u64(mpq_denref(qr), v == UINT64_MAX ? 3 : v + 1);
            CHECK(mpq_equal(q, qr), "mpq_set_u64(%" PRIu64 ")", v);
            mpq_clear(q); mpq_clear(qr);
        }
    }
    mpz_clears(a, r, x, y, t, NULL);
}

#ifdef USE_MPFR
static void test_mpfr(void) {
    mpfr_t x, r, ref, a;
    mpfr_inits2(200, x, r, ref, a, NULL);      /* 200 bits: every int64 is exact */
    mpz_t y; mpz_init(y);
    mpfr_set_str(a, "1.25", 10, MPFR_RNDN);
    for (size_t i = 0; i < N_I; i++) {
        int64_t v = I64[i];
        ref_i64(y, v);
        mth_mpfr_set_i64(x, v, MPFR_RNDN); mpfr_set_z(ref, y, MPFR_RNDN);
        CHECK(mpfr_equal_p(x, ref), "mpfr_set_i64(%" PRId64 ")", v);
        mth_mpfr_mul_i64(r, a, v, MPFR_RNDN); mpfr_mul_z(ref, a, y, MPFR_RNDN);
        CHECK(mpfr_equal_p(r, ref), "mpfr_mul_i64(%" PRId64 ")", v);
        mth_mpfr_add_i64(r, a, v, MPFR_RNDN); mpfr_add_z(ref, a, y, MPFR_RNDN);
        CHECK(mpfr_equal_p(r, ref), "mpfr_add_i64(%" PRId64 ")", v);
        if (v != 0) {
            mth_mpfr_div_i64(r, a, v, MPFR_RNDN); mpfr_div_z(ref, a, y, MPFR_RNDN);
            CHECK(mpfr_equal_p(r, ref), "mpfr_div_i64(%" PRId64 ")", v);
        }
    }
    for (size_t i = 0; i < N_U; i++) {
        uint64_t v = U64[i];
        ref_u64(y, v);
        mth_mpfr_set_u64(x, v, MPFR_RNDN); mpfr_set_z(ref, y, MPFR_RNDN);
        CHECK(mpfr_equal_p(x, ref), "mpfr_set_u64(%" PRIu64 ")", v);
    }
    /* pow: 1.25^e against the exact mpq power for small |e|; for 2^32 + 1 the
     * truncated exponent would be 1, so 1.25^(2^32+1) must not equal 1.25. */
    static const int64_t E[] = { 0, 1, -1, 5, -5, 31, -31 };
    for (size_t i = 0; i < sizeof E / sizeof E[0]; i++) {
        mth_mpfr_pow_i64(r, a, E[i], MPFR_RNDN); mpfr_pow_si(ref, a, (long)E[i], MPFR_RNDN);
        CHECK(mpfr_equal_p(r, ref), "mpfr_pow_i64(%" PRId64 ")", E[i]);
    }
    mth_mpfr_pow_i64(r, a, INT64_C(4294967297), MPFR_RNDN);
    CHECK(mpfr_inf_p(r) || mpfr_cmp_d(r, 1e300) > 0, "mpfr_pow_i64(2^32+1) is huge, not 1.25");
    mpz_clear(y);
    mpfr_clears(x, r, ref, a, NULL);
}
#endif

int main(void) {
    test_mpz();
#ifdef USE_MPFR
    test_mpfr();
#endif
    printf("gmp64_llp64: %d checks, %d failures\n", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
