/* Mathilda -- int64/uint64 siblings of GMP's and MPFR's `_si`/`_ui` calls.
 *
 * GMP's mpz_set_si / mpz_mul_si / mpz_cmp_si / mpq_set_si / mpz_fdiv_ui / ...
 * and MPFR's mpfr_set_si / mpfr_mul_si / ... take a `long` or `unsigned long`.
 * That is 64 bits on LP64 (macOS, Linux) but only 32 on LLP64 (native
 * Windows), where passing an int64_t above 2^31 silently truncates.  Every
 * helper here is exact for the whole int64/uint64 range on every data model;
 * where `long` is 64 bits each is just the `_si`/`_ui` call.  Declared in
 * expr.h.  See docs/design/numeric_tower_review.md, finding F2.
 *
 * tests/test_gmp64_llp64.c includes this file with MTH_TEST_LLP64_PATH
 * defined, which forces the LLP64 branch (with a 32-bit `long` range) so the
 * Windows code is exercised on an LP64 host.
 */
#include "expr.h"
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>

/* MTH_LONG_IS_64: `long` holds every int64 (LP64), so each helper below is the
 * plain `_si`/`_ui` call.  Otherwise (LLP64 Windows, where `long` is 32 bits) a
 * value that fits `long` still takes the `_si`/`_ui` call and a wider one goes
 * through an exact mpz.  MTH_TEST_LLP64_PATH forces the second branch, with a
 * 32-bit `long` range, so the Windows code can be exercised on an LP64 host
 * (tests/test_gmp64_llp64.c builds a copy of these with it set). */
#if LONG_MAX >= INT64_MAX && !defined(MTH_TEST_LLP64_PATH)
#  define MTH_LONG_IS_64 1
#else
#  ifdef MTH_TEST_LLP64_PATH
#    define MTH_LONG_LO  INT32_MIN
#    define MTH_LONG_HI  INT32_MAX
#    define MTH_ULONG_HI UINT32_MAX
#  else
#    define MTH_LONG_LO  LONG_MIN
#    define MTH_LONG_HI  LONG_MAX
#    define MTH_ULONG_HI ULONG_MAX
#  endif
#  define MTH_FITS_LONG(v)  ((v) >= MTH_LONG_LO && (v) <= MTH_LONG_HI)
#  define MTH_FITS_ULONG(v) ((v) <= MTH_ULONG_HI)
#endif

bool mth_u64_fits_ulong(uint64_t v) {
#ifdef MTH_LONG_IS_64
    (void)v;
    return true;
#else
    return MTH_FITS_ULONG(v);
#endif
}

void mth_mpz_set_u64(mpz_t z, uint64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_set_ui(z, (unsigned long)v);
#else
    if (MTH_FITS_ULONG(v)) { mpz_set_ui(z, (unsigned long)v); return; }
    mpz_import(z, 1, -1, sizeof(uint64_t), 0, 0, &v);
#endif
}

void mth_mpz_set_i64(mpz_t z, int64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_set_si(z, (long)v);
#else
    if (MTH_FITS_LONG(v)) { mpz_set_si(z, (long)v); return; }
    /* Import the magnitude as a uint64 word, then apply the sign.  The
     * magnitude is computed in uint64 so INT64_MIN does not overflow. */
    uint64_t mag = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_import(z, 1, -1, sizeof(uint64_t), 0, 0, &mag);
    if (v < 0) mpz_neg(z, z);
#endif
}

void mth_mpz_init_set_i64(mpz_t z, int64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_init_set_si(z, (long)v);
#else
    mpz_init(z);
    mth_mpz_set_i64(z, v);
#endif
}

void mth_mpz_init_set_u64(mpz_t z, uint64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_init_set_ui(z, (unsigned long)v);
#else
    mpz_init(z);
    mth_mpz_set_u64(z, v);
#endif
}

void mth_mpz_mul_i64(mpz_t r, const mpz_t a, int64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_mul_si(r, a, (long)v);
#else
    if (MTH_FITS_LONG(v)) { mpz_mul_si(r, a, (long)v); return; }
    mpz_t t; mpz_init(t); mth_mpz_set_i64(t, v);
    mpz_mul(r, a, t);
    mpz_clear(t);
#endif
}

void mth_mpz_mul_u64(mpz_t r, const mpz_t a, uint64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_mul_ui(r, a, (unsigned long)v);
#else
    if (MTH_FITS_ULONG(v)) { mpz_mul_ui(r, a, (unsigned long)v); return; }
    mpz_t t; mpz_init(t); mth_mpz_set_u64(t, v);
    mpz_mul(r, a, t);
    mpz_clear(t);
#endif
}

void mth_mpz_add_u64(mpz_t r, const mpz_t a, uint64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_add_ui(r, a, (unsigned long)v);
#else
    if (MTH_FITS_ULONG(v)) { mpz_add_ui(r, a, (unsigned long)v); return; }
    mpz_t t; mpz_init(t); mth_mpz_set_u64(t, v);
    mpz_add(r, a, t);
    mpz_clear(t);
#endif
}

void mth_mpz_sub_u64(mpz_t r, const mpz_t a, uint64_t v) {
#ifdef MTH_LONG_IS_64
    mpz_sub_ui(r, a, (unsigned long)v);
#else
    if (MTH_FITS_ULONG(v)) { mpz_sub_ui(r, a, (unsigned long)v); return; }
    mpz_t t; mpz_init(t); mth_mpz_set_u64(t, v);
    mpz_sub(r, a, t);
    mpz_clear(t);
#endif
}

int mth_mpz_cmp_i64(const mpz_t a, int64_t v) {
#ifdef MTH_LONG_IS_64
    return mpz_cmp_si(a, (long)v);
#else
    if (MTH_FITS_LONG(v)) return mpz_cmp_si(a, (long)v);
    mpz_t t; mpz_init(t); mth_mpz_set_i64(t, v);
    int c = mpz_cmp(a, t);
    mpz_clear(t);
    return c;
#endif
}

int mth_mpz_cmp_u64(const mpz_t a, uint64_t v) {
#ifdef MTH_LONG_IS_64
    return mpz_cmp_ui(a, (unsigned long)v);
#else
    if (MTH_FITS_ULONG(v)) return mpz_cmp_ui(a, (unsigned long)v);
    mpz_t t; mpz_init(t); mth_mpz_set_u64(t, v);
    int c = mpz_cmp(a, t);
    mpz_clear(t);
    return c;
#endif
}

int mth_mpz_divisible_u64_p(const mpz_t a, uint64_t d) {
#ifdef MTH_LONG_IS_64
    return mpz_divisible_ui_p(a, (unsigned long)d);
#else
    if (MTH_FITS_ULONG(d)) return mpz_divisible_ui_p(a, (unsigned long)d);
    mpz_t t; mpz_init(t); mth_mpz_set_u64(t, d);
    int r = mpz_divisible_p(a, t);
    mpz_clear(t);
    return r;
#endif
}

void mth_mpz_divexact_u64(mpz_t r, const mpz_t a, uint64_t d) {
#ifdef MTH_LONG_IS_64
    mpz_divexact_ui(r, a, (unsigned long)d);
#else
    if (MTH_FITS_ULONG(d)) { mpz_divexact_ui(r, a, (unsigned long)d); return; }
    mpz_t t; mpz_init(t); mth_mpz_set_u64(t, d);
    mpz_divexact(r, a, t);
    mpz_clear(t);
#endif
}

uint64_t mth_mpz_fdiv_u64(const mpz_t a, uint64_t d) {
#ifdef MTH_LONG_IS_64
    return (uint64_t)mpz_fdiv_ui(a, (unsigned long)d);
#else
    if (MTH_FITS_ULONG(d)) return (uint64_t)mpz_fdiv_ui(a, (unsigned long)d);
    mpz_t t, m; mpz_init(t); mpz_init(m); mth_mpz_set_u64(t, d);
    mpz_fdiv_r(m, a, t);                        /* 0 <= m < d < 2^64 */
    uint64_t out = 0; size_t cnt = 0;
    mpz_export(&out, &cnt, -1, sizeof(uint64_t), 0, 0, m);
    mpz_clear(t); mpz_clear(m);
    return cnt ? out : 0;
#endif
}

uint64_t mth_mpz_get_u64(const mpz_t z) {
#ifdef MTH_LONG_IS_64
    return (uint64_t)mpz_get_ui(z);
#else
    /* The low 64 bits of |z|, as mpz_get_ui gives the low `long` bits. */
    if (mpz_sizeinbase(z, 2) <= 32) return (uint64_t)mpz_get_ui(z);
    mpz_t t; mpz_init(t);
    mpz_abs(t, z);
    mpz_tdiv_r_2exp(t, t, 64);
    uint64_t out = 0; size_t cnt = 0;
    mpz_export(&out, &cnt, -1, sizeof(uint64_t), 0, 0, t);
    mpz_clear(t);
    return out;
#endif
}

void mth_mpq_set_u64(mpq_t q, uint64_t num, uint64_t den) {
#ifdef MTH_LONG_IS_64
    mpq_set_ui(q, (unsigned long)num, (unsigned long)den);
#else
    if (MTH_FITS_ULONG(num) && MTH_FITS_ULONG(den)) {
        mpq_set_ui(q, (unsigned long)num, (unsigned long)den); return;
    }
    mth_mpz_set_u64(mpq_numref(q), num);        /* like mpq_set_ui: not   */
    mth_mpz_set_u64(mpq_denref(q), num ? den : 1); /* canonicalised; 0 is 0/1 */
#endif
}

void mth_mpq_set_i64(mpq_t q, int64_t num, uint64_t den) {
#ifdef MTH_LONG_IS_64
    mpq_set_si(q, (long)num, (unsigned long)den);
#else
    if (MTH_FITS_LONG(num) && MTH_FITS_ULONG(den)) {
        mpq_set_si(q, (long)num, (unsigned long)den); return;
    }
    mth_mpz_set_i64(mpq_numref(q), num);        /* like mpq_set_si: not   */
    mth_mpz_set_u64(mpq_denref(q), num ? den : 1); /* canonicalised; 0 is 0/1 */
#endif
}

#ifdef USE_MPFR
/* MPFR's intmax_t entry points are exact for every int64 on every data model
 * (expr.h includes <stdint.h> before <mpfr.h>, which enables them). */
int mth_mpfr_set_i64(mpfr_t x, int64_t v, mpfr_rnd_t rnd) {
#ifdef MTH_LONG_IS_64
    return mpfr_set_si(x, (long)v, rnd);
#else
    return mpfr_set_sj(x, (intmax_t)v, rnd);
#endif
}

int mth_mpfr_set_u64(mpfr_t x, uint64_t v, mpfr_rnd_t rnd) {
#ifdef MTH_LONG_IS_64
    return mpfr_set_ui(x, (unsigned long)v, rnd);
#else
    return mpfr_set_uj(x, (uintmax_t)v, rnd);
#endif
}

int mth_mpfr_mul_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd) {
#ifdef MTH_LONG_IS_64
    return mpfr_mul_si(r, a, (long)v, rnd);
#else
    if (MTH_FITS_LONG(v)) return mpfr_mul_si(r, a, (long)v, rnd);
    mpz_t t; mpz_init(t); mth_mpz_set_i64(t, v);
    int ter = mpfr_mul_z(r, a, t, rnd);         /* one rounding, as _si */
    mpz_clear(t);
    return ter;
#endif
}

int mth_mpfr_div_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd) {
#ifdef MTH_LONG_IS_64
    return mpfr_div_si(r, a, (long)v, rnd);
#else
    if (MTH_FITS_LONG(v)) return mpfr_div_si(r, a, (long)v, rnd);
    mpz_t t; mpz_init(t); mth_mpz_set_i64(t, v);
    int ter = mpfr_div_z(r, a, t, rnd);
    mpz_clear(t);
    return ter;
#endif
}

int mth_mpfr_add_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd) {
#ifdef MTH_LONG_IS_64
    return mpfr_add_si(r, a, (long)v, rnd);
#else
    if (MTH_FITS_LONG(v)) return mpfr_add_si(r, a, (long)v, rnd);
    mpz_t t; mpz_init(t); mth_mpz_set_i64(t, v);
    int ter = mpfr_add_z(r, a, t, rnd);
    mpz_clear(t);
    return ter;
#endif
}

int mth_mpfr_pow_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd) {
#ifdef MTH_LONG_IS_64
    return mpfr_pow_si(r, a, (long)v, rnd);
#else
    return mpfr_pow_sj(r, a, (intmax_t)v, rnd);
#endif
}
#endif /* USE_MPFR */
