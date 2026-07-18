/* Mathilda — Fourier / InverseFourier: the discrete Fourier transform.
 *
 * See fourier.h for the high-level description. The transform is defined,
 * with FourierParameters -> {a, b} (default {0, 1}), as
 *
 *   Fourier[list]_s        = N^{-(1-a)/2} Sum_r u_r Exp[+2 Pi I b (r-1)(s-1)/n]
 *   InverseFourier[list]_r = N^{-(1+a)/2} Sum_s v_s Exp[-2 Pi I b (r-1)(s-1)/n]
 *
 * where n is a dimension length and N the total element count. We factor this
 * as a *standard* unnormalised transform F[u]_k = Sum_j u_j Exp[base_sign 2 Pi
 * I jk/n] (base_sign = +1 for Fourier, -1 for InverseFourier) evaluated at the
 * gathered index (b*k) mod n, then scaled by the single factor N^{-(1∓a)/2}.
 * The b-gather (identity when b = 1) folds the FourierParameters `b` — including
 * the non-invertible |b| > 1 cases — onto a single transform. For a
 * multidimensional (rectangular) array the transform and the gather are applied
 * independently per axis (they are separable).
 *
 * Regimes:
 *   - machine  : FFTW nD plan (USE_FFTW) or a naive O(n^2) separable fallback;
 *   - arbitrary: an MPFR-complex FFT (radix-2 for powers of two, Bluestein
 *                otherwise) applied per axis;
 *   - symbolic : the exact transform built from Exp[..] roots of unity, handed
 *                to the evaluator to simplify.
 *
 * `b` must be a (machine) integer for the numeric regimes; the symbolic regime
 * accepts any exact `b`. Non-integer `b` on numeric input stays unevaluated.
 */

#include "fourier.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "arithmetic.h"
#include "numeric.h"
#include "common.h"
#include "ndarray.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifdef USE_MPFR
#include "numeric_complex.h"
#include <mpfr.h>
#endif

#ifdef USE_FFTW
#include <fftw3.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FOURIER_MAX_RANK NDARRAY_MAX_RANK

/* ==================================================================== *
 *  Small numeric-leaf helpers
 * ==================================================================== */

/* Extract a numeric leaf as a machine complex (re, im). Handles Integer,
 * BigInt, Real, MPFR, Rational[n,d] and Complex[..]; false for anything else. */
static bool leaf_to_cdouble(const Expr* e, double* re, double* im) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *re = (double)e->data.integer; *im = 0.0; return true;
        case EXPR_REAL:    *re = e->data.real;            *im = 0.0; return true;
        case EXPR_BIGINT:  *re = mpz_get_d(e->data.bigint); *im = 0.0; return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); *im = 0.0; return true;
#endif
        default: break;
    }
    int64_t rn, rd;
    if (is_rational(e, &rn, &rd)) { *re = (double)rn / (double)rd; *im = 0.0; return true; }
    Expr *cr, *ci;
    if (is_complex((Expr*)e, &cr, &ci)) {
        double r1, i1, r2, i2;
        if (!leaf_to_cdouble(cr, &r1, &i1)) return false;
        if (!leaf_to_cdouble(ci, &r2, &i2)) return false;
        /* (r1 + i1 I) + (r2 + i2 I) I = (r1 - i2) + (i1 + r2) I */
        *re = r1 - i2;
        *im = i1 + r2;
        return true;
    }
    return false;
}

/* A real numeric leaf as a double (for the FourierParameters `a`). */
static bool expr_to_double_num(const Expr* e, double* out) {
    double re, im;
    if (!leaf_to_cdouble(e, &re, &im)) return false;
    *out = re;
    return true;
}

/* An exact-integer leaf as a long (for the FourierParameters `b`). */
static bool expr_to_long_int(const Expr* e, long* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* ==================================================================== *
 *  Nested-array shape: rectangular validation + flatten
 * ==================================================================== */

static bool nd_verify(const Expr* e, const int64_t* dims, int level, int rank) {
    if (level == rank) return !head_is(e, SYM_List);      /* leaf must be scalar */
    if (!head_is(e, SYM_List)) return false;
    if ((int64_t)e->data.function.arg_count != dims[level]) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!nd_verify(e->data.function.args[i], dims, level + 1, rank)) return false;
    return true;
}

/* Fill dims[0..*rank-1] for a perfectly rectangular nested List of scalar
 * leaves. Returns false if `e` isn't a List, is ragged, or nests too deep. */
static bool nested_dims(const Expr* e, int64_t* dims, int* rank_out) {
    int rank = 0;
    const Expr* cur = e;
    while (head_is(cur, SYM_List)) {
        if (rank >= FOURIER_MAX_RANK) return false;
        dims[rank++] = (int64_t)cur->data.function.arg_count;
        if (cur->data.function.arg_count == 0) { *rank_out = rank; return true; }
        cur = cur->data.function.args[0];
    }
    if (rank == 0) return false;
    *rank_out = rank;
    return nd_verify(e, dims, 0, rank);
}

/* Collect leaves in row-major order into `out` (borrowed pointers). */
static void nd_flatten(const Expr* e, int level, int rank, Expr** out, size_t* idx) {
    if (level == rank) { out[(*idx)++] = (Expr*)e; return; }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nd_flatten(e->data.function.args[i], level + 1, rank, out, idx);
}

/* Rebuild a nested List from row-major leaves; ownership of `leaves` transfers
 * into the returned tree. */
static Expr* build_nested(Expr** leaves, const int64_t* dims, int rank,
                          int level, size_t* idx) {
    if (level == rank) return leaves[(*idx)++];
    int64_t n = dims[level];
    Expr** a = malloc((size_t)n * sizeof(Expr*));
    for (int64_t i = 0; i < n; i++)
        a[i] = build_nested(leaves, dims, rank, level + 1, idx);
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), a, (size_t)n);
    free(a);
    return l;
}

static void row_major_strides(const int64_t* dims, int rank, int64_t* strides) {
    strides[rank - 1] = 1;
    for (int i = rank - 2; i >= 0; i--) strides[i] = strides[i + 1] * dims[i + 1];
}

/* ==================================================================== *
 *  Regime classification
 * ==================================================================== */

typedef enum { REG_SYMBOLIC, REG_MACHINE, REG_ARB } Regime;

static Regime classify(Expr** leaves, size_t n, long* arb_bits_out) {
    bool has_mpfr = false;
    long min_bits = 0;
    for (size_t i = 0; i < n; i++) {
        if (!expr_is_numeric_like(leaves[i])) return REG_SYMBOLIC;
#ifdef USE_MPFR
        if (numeric_expr_is_mpfr(leaves[i])) {
            has_mpfr = true;
            long b = numeric_min_inexact_bits(leaves[i]);
            if (b > 0 && (min_bits == 0 || b < min_bits)) min_bits = b;
        }
#endif
    }
    if (has_mpfr) {
        if (min_bits < 53) min_bits = 53;
        *arb_bits_out = min_bits;
        return REG_ARB;
    }
    return REG_MACHINE;
}

/* Gathered source index along one axis: ((b*k) mod n + n) mod n. */
static int64_t gather_index(long b, int64_t k, int64_t n) {
    long bm = b % n;                       /* keep the product small */
    int64_t g = ((int64_t)bm * k) % n;
    if (g < 0) g += n;
    return g;
}

/* ==================================================================== *
 *  Machine-precision path
 * ==================================================================== */

#ifndef USE_FFTW
/* Naive separable standard transform (base_sign, b = 1) on an interleaved
 * (re, im) buffer — the fallback when FFTW is unavailable. */
static void naive_std_transform(double* buf, const int64_t* dims, int rank,
                                int base_sign) {
    int64_t strides[FOURIER_MAX_RANK];
    row_major_strides(dims, rank, strides);
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];

    for (int ax = 0; ax < rank; ax++) {
        int64_t n = dims[ax], s = strides[ax];
        double* cw = malloc((size_t)n * sizeof(double));
        double* sw = malloc((size_t)n * sizeof(double));
        double* tre = malloc((size_t)n * sizeof(double));
        double* tim = malloc((size_t)n * sizeof(double));
        for (int64_t k = 0; k < n; k++) {
            double ang = (double)base_sign * 2.0 * M_PI * (double)k / (double)n;
            cw[k] = cos(ang);
            sw[k] = sin(ang);
        }
        for (int64_t o = 0; o < N; o++) {
            if ((o / s) % n != 0) continue;          /* line base only */
            for (int64_t r = 0; r < n; r++) {
                int64_t idx = o + r * s;
                tre[r] = buf[2 * idx];
                tim[r] = buf[2 * idx + 1];
            }
            for (int64_t sidx = 0; sidx < n; sidx++) {
                double ar = 0.0, ai = 0.0;
                for (int64_t r = 0; r < n; r++) {
                    int64_t w = (r * sidx) % n;
                    double c = cw[w], sn = sw[w];
                    ar += tre[r] * c - tim[r] * sn;
                    ai += tre[r] * sn + tim[r] * c;
                }
                int64_t idx = o + sidx * s;
                buf[2 * idx] = ar;
                buf[2 * idx + 1] = ai;
            }
        }
        free(cw); free(sw); free(tre); free(tim);
    }
}
#endif

/* Real-collapse tolerance test + NDArray construction from an interleaved
 * (re, im) buffer `out` (length 2N, malloc'd — ownership transfers into the
 * returned EXPR_NDARRAY, NDT_COMPLEX64 or NDT_FLOAT64 when roundoff-real). */
static Expr* machine_build_ndarray(double* out, int64_t N, const int64_t* dims, int rank) {
    double max_mag = 0.0, max_im = 0.0;
    for (int64_t o = 0; o < N; o++) {
        double re = out[2 * o], im = out[2 * o + 1];
        double mag = hypot(re, im);
        if (mag > max_mag) max_mag = mag;
        double a = fabs(im);
        if (a > max_im) max_im = a;
    }
    double tol = 16.0 * (double)N * DBL_EPSILON;
    bool real_only = (max_mag == 0.0) || (max_im <= tol * max_mag);
    if (real_only) {
        double* rout = malloc((size_t)N * sizeof(double));
        for (int64_t o = 0; o < N; o++) rout[o] = out[2 * o];
        free(out);
        return expr_new_ndarray(rank, dims, rout, NDT_FLOAT64);
    }
    return expr_new_ndarray(rank, dims, out, NDT_COMPLEX64);
}

/* FFT (FFTW or naive) + b-gather + normalisation. Consumes `buf` (an
 * interleaved (re, im) double buffer of length 2N) and returns a fresh
 * malloc'd interleaved output buffer. */
static double* machine_transform_buf(double* buf, const int64_t* dims, int rank,
                                     double a, long b, int base_sign) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    double norm = pow((double)N, base_sign > 0 ? -(1.0 - a) / 2.0 : -(1.0 + a) / 2.0);
    double* std_buf;

#ifdef USE_FFTW
    fftw_complex* fb = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
    memcpy(fb, buf, (size_t)N * 2 * sizeof(double));
    free(buf);
    int idims[FOURIER_MAX_RANK];
    for (int i = 0; i < rank; i++) idims[i] = (int)dims[i];
    /* FFTW_BACKWARD (+1) matches our +2πi convention; FFTW_FORWARD (-1) the -. */
    int sign = (base_sign > 0) ? FFTW_BACKWARD : FFTW_FORWARD;
    fftw_plan p = fftw_plan_dft(rank, idims, fb, fb, sign, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    std_buf = (double*)fb;
#else
    naive_std_transform(buf, dims, rank, base_sign);
    std_buf = buf;
#endif

    int64_t strides[FOURIER_MAX_RANK];
    row_major_strides(dims, rank, strides);
    double* out = malloc((size_t)N * 2 * sizeof(double));
    for (int64_t o = 0; o < N; o++) {
        int64_t src = 0;
        for (int i = 0; i < rank; i++) {
            int64_t k = (o / strides[i]) % dims[i];
            src += gather_index(b, k, dims[i]) * strides[i];
        }
        out[2 * o]     = norm * std_buf[2 * src];
        out[2 * o + 1] = norm * std_buf[2 * src + 1];
    }

#ifdef USE_FFTW
    fftw_free(fb);
#else
    free(buf);
#endif
    return out;
}

/* List input: numericalise, pack, transform, then delist the packed result so
 * the surface value is a plain List (WMA prints `{...}`, Head is List). */
static Expr* machine_path(const Expr* nested, const int64_t* dims, int rank,
                          double a, long b, int base_sign) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    if (N <= 0) return NULL;

    Expr* nnum = numericalize(nested, numeric_machine_spec());
    if (!nnum) return NULL;
    int64_t d2[FOURIER_MAX_RANK];
    int r2 = 0;
    if (!nested_dims(nnum, d2, &r2) || r2 != rank) { expr_free(nnum); return NULL; }
    for (int i = 0; i < rank; i++) if (d2[i] != dims[i]) { expr_free(nnum); return NULL; }

    Expr** leaves = malloc((size_t)N * sizeof(Expr*));
    size_t idx = 0;
    nd_flatten(nnum, 0, rank, leaves, &idx);
    double* buf = malloc((size_t)N * 2 * sizeof(double));
    for (int64_t o = 0; o < N; o++) {
        double re, im;
        if (!leaf_to_cdouble(leaves[o], &re, &im)) {
            free(leaves); free(buf); expr_free(nnum); return NULL;
        }
        buf[2 * o] = re;
        buf[2 * o + 1] = im;
    }
    free(leaves);
    expr_free(nnum);

    double* out = machine_transform_buf(buf, dims, rank, a, b, base_sign);
    Expr* nd = machine_build_ndarray(out, N, dims, rank);
    Expr* lst = ndarray_to_nested_list(nd);
    expr_free(nd);
    return lst;
}

/* NDArray fast path: read the packed buffer directly (no List round trip) and
 * return an NDArray. An NDArray is always machine-precision numeric. */
static Expr* machine_path_ndarray(const Expr* nd, double a, long b, int base_sign) {
    int rank = nd->data.ndarray.rank;
    const int64_t* dims = nd->data.ndarray.dims;
    NDType dt = nd->data.ndarray.dtype;
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    if (N <= 0) return NULL;

    /* Pack every element to an interleaved (re, im) double via the dtype-
     * agnostic ndt_get choke point (real dtypes give im = 0). */
    double* buf = malloc((size_t)N * 2 * sizeof(double));
    for (int64_t o = 0; o < N; o++)
        ndt_get(nd->data.ndarray.data, (size_t)o, dt, &buf[2 * o], &buf[2 * o + 1]);

    double* out = machine_transform_buf(buf, dims, rank, a, b, base_sign);
    return machine_build_ndarray(out, N, dims, rank);
}

/* ==================================================================== *
 *  Arbitrary-precision path (MPFR-complex FFT: radix-2 + Bluestein)
 * ==================================================================== */

#ifdef USE_MPFR

static bool is_pow2(size_t n) { return n && ((n & (n - 1)) == 0); }

static size_t next_pow2(size_t n) {
    size_t m = 1;
    while (m < n) m <<= 1;
    return m;
}

static size_t bit_reverse(size_t x, int bits) {
    size_t r = 0;
    for (int i = 0; i < bits; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}

/* In-place radix-2 Cooley–Tukey FFT of a power-of-two length `m`.
 * Computes a_k <- Sum_j a_j Exp[sign 2 Pi I jk/m] (unnormalised). */
static void fft_pow2(ncpx* a, size_t m, int sign, mpfr_prec_t wp) {
    if (m <= 1) return;
    int bits = 0;
    while (((size_t)1 << bits) < m) bits++;
    for (size_t i = 0; i < m; i++) {
        size_t j = bit_reverse(i, bits);
        if (i < j) {
            ncpx tmp; ncpx_init(&tmp, wp);
            ncpx_set(&tmp, &a[i]); ncpx_set(&a[i], &a[j]); ncpx_set(&a[j], &tmp);
            ncpx_clear(&tmp);
        }
    }
    /* Full twiddle table W[k] = Exp[sign 2 Pi I k/m], k = 0 .. m/2-1. */
    ncpx* W = malloc((m / 2) * sizeof(ncpx));
    mpfr_t two_pi, ang;
    mpfr_init2(two_pi, wp); mpfr_init2(ang, wp);
    mpfr_const_pi(two_pi, MPFR_RNDN);
    mpfr_mul_2ui(two_pi, two_pi, 1, MPFR_RNDN);       /* 2*pi */
    for (size_t k = 0; k < m / 2; k++) {
        ncpx_init(&W[k], wp);
        mpfr_mul_si(ang, two_pi, sign, MPFR_RNDN);
        mpfr_mul_ui(ang, ang, (unsigned long)k, MPFR_RNDN);
        mpfr_div_ui(ang, ang, (unsigned long)m, MPFR_RNDN);
        mpfr_sin_cos(W[k].im, W[k].re, ang, MPFR_RNDN);
    }
    mpfr_clear(two_pi); mpfr_clear(ang);

    ncpx t, u;
    ncpx_init(&t, wp); ncpx_init(&u, wp);
    for (size_t len = 2; len <= m; len <<= 1) {
        size_t half = len / 2, step = m / len;
        for (size_t i = 0; i < m; i += len) {
            for (size_t j = 0; j < half; j++) {
                ncpx_mul(&t, &W[j * step], &a[i + j + half], wp);
                ncpx_set(&u, &a[i + j]);
                ncpx_add(&a[i + j], &u, &t);
                ncpx_sub(&a[i + j + half], &u, &t);
            }
        }
    }
    ncpx_clear(&t); ncpx_clear(&u);
    for (size_t k = 0; k < m / 2; k++) ncpx_clear(&W[k]);
    free(W);
}

/* Bluestein (chirp-z) transform for arbitrary length `n`. */
static void fft_bluestein(ncpx* a, size_t n, int sign, mpfr_prec_t wp) {
    size_t m = next_pow2(2 * n - 1);

    /* Chirp w_j = Exp[sign Pi I (j^2 mod 2n)/n]. */
    ncpx* w = malloc(n * sizeof(ncpx));
    mpfr_t pi, ang;
    mpfr_init2(pi, wp); mpfr_init2(ang, wp);
    mpfr_const_pi(pi, MPFR_RNDN);
    for (size_t j = 0; j < n; j++) {
        ncpx_init(&w[j], wp);
        unsigned long e = (unsigned long)(((unsigned long long)j * j) % (2ULL * n));
        mpfr_mul_si(ang, pi, sign, MPFR_RNDN);
        mpfr_mul_ui(ang, ang, e, MPFR_RNDN);
        mpfr_div_ui(ang, ang, (unsigned long)n, MPFR_RNDN);
        mpfr_sin_cos(w[j].im, w[j].re, ang, MPFR_RNDN);
    }
    mpfr_clear(pi); mpfr_clear(ang);

    ncpx* A = malloc(m * sizeof(ncpx));
    ncpx* B = malloc(m * sizeof(ncpx));
    for (size_t i = 0; i < m; i++) { ncpx_init(&A[i], wp); ncpx_init(&B[i], wp); ncpx_set_ui(&A[i], 0); ncpx_set_ui(&B[i], 0); }

    /* A[j] = a[j]*w[j]; B[0]=1, B[j]=B[m-j]=conj-chirp w[j]^{-1} = Exp[-..]. */
    for (size_t j = 0; j < n; j++) ncpx_mul(&A[j], &a[j], &w[j], wp);
    /* conj of w[j] is Exp[-sign Pi I e/n]; since |w|=1, w^{-1} = conj(w). */
    for (size_t j = 0; j < n; j++) {
        ncpx cj; ncpx_init(&cj, wp);
        mpfr_set(cj.re, w[j].re, MPFR_RNDN);
        mpfr_neg(cj.im, w[j].im, MPFR_RNDN);
        ncpx_set(&B[j], &cj);
        if (j != 0) ncpx_set(&B[m - j], &cj);
        ncpx_clear(&cj);
    }

    /* Circular convolution via a consistent forward/inverse FFT pair. */
    fft_pow2(A, m, +1, wp);
    fft_pow2(B, m, +1, wp);
    for (size_t i = 0; i < m; i++) ncpx_mul(&A[i], &A[i], &B[i], wp);
    fft_pow2(A, m, -1, wp);
    mpfr_t inv_m; mpfr_init2(inv_m, wp);
    mpfr_set_ui(inv_m, (unsigned long)m, MPFR_RNDN);
    mpfr_ui_div(inv_m, 1, inv_m, MPFR_RNDN);
    for (size_t i = 0; i < m; i++) ncpx_scale(&A[i], &A[i], inv_m);
    mpfr_clear(inv_m);

    /* a[k] = w[k] * conv[k]. */
    for (size_t k = 0; k < n; k++) ncpx_mul(&a[k], &w[k], &A[k], wp);

    for (size_t i = 0; i < m; i++) { ncpx_clear(&A[i]); ncpx_clear(&B[i]); }
    for (size_t j = 0; j < n; j++) ncpx_clear(&w[j]);
    free(A); free(B); free(w);
}

static void fft_mpfr(ncpx* a, size_t n, int sign, mpfr_prec_t wp) {
    if (n <= 1) return;
    if (is_pow2(n)) fft_pow2(a, n, sign, wp);
    else            fft_bluestein(a, n, sign, wp);
}

static Expr* arb_path(Expr** leaves, const int64_t* dims, int rank,
                      double a, long b, int base_sign, long target_bits) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    if (N <= 0) return NULL;

    long guard = 32;
    { int64_t t = N; while (t > 1) { guard += 2; t >>= 1; } }
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);

    /* Load leaves into an ncpx buffer at working precision. */
    ncpx* buf = malloc((size_t)N * sizeof(ncpx));
    for (int64_t o = 0; o < N; o++) {
        ncpx_init(&buf[o], wp);
        mpfr_t re, im; mpfr_init2(re, wp); mpfr_init2(im, wp);
        bool inexact = false;
        if (!get_approx_mpfr(leaves[o], re, im, &inexact)) {
            mpfr_clear(re); mpfr_clear(im);
            for (int64_t j = 0; j <= o; j++) ncpx_clear(&buf[j]);
            free(buf);
            return NULL;
        }
        mpfr_set(buf[o].re, re, MPFR_RNDN);
        mpfr_set(buf[o].im, im, MPFR_RNDN);
        mpfr_clear(re); mpfr_clear(im);
    }

    /* Standard transform per axis (strided lines). */
    int64_t strides[FOURIER_MAX_RANK];
    row_major_strides(dims, rank, strides);
    for (int ax = 0; ax < rank; ax++) {
        int64_t n = dims[ax], s = strides[ax];
        if (n <= 1) continue;
        ncpx* line = malloc((size_t)n * sizeof(ncpx));
        for (int64_t k = 0; k < n; k++) ncpx_init(&line[k], wp);
        for (int64_t o = 0; o < N; o++) {
            if ((o / s) % n != 0) continue;
            for (int64_t r = 0; r < n; r++) ncpx_set(&line[r], &buf[o + r * s]);
            fft_mpfr(line, (size_t)n, base_sign, wp);
            for (int64_t r = 0; r < n; r++) ncpx_set(&buf[o + r * s], &line[r]);
        }
        for (int64_t k = 0; k < n; k++) ncpx_clear(&line[k]);
        free(line);
    }

    /* Normalisation factor N^{-(1∓a)/2} at working precision. */
    mpfr_t norm, Nm, expo;
    mpfr_init2(norm, wp); mpfr_init2(Nm, wp); mpfr_init2(expo, wp);
    mpfr_set_si(Nm, (long)N, MPFR_RNDN);
    mpfr_set_d(expo, base_sign > 0 ? -(1.0 - a) / 2.0 : -(1.0 + a) / 2.0, MPFR_RNDN);
    mpfr_pow(norm, Nm, expo, MPFR_RNDN);
    mpfr_clear(Nm); mpfr_clear(expo);

    /* b-gather into fresh leaves at the target precision. */
    Expr** out = malloc((size_t)N * sizeof(Expr*));
    for (int64_t o = 0; o < N; o++) {
        int64_t src = 0;
        for (int i = 0; i < rank; i++) {
            int64_t k = (o / strides[i]) % dims[i];
            src += gather_index(b, k, dims[i]) * strides[i];
        }
        mpfr_t rre, rim; mpfr_init2(rre, target_bits); mpfr_init2(rim, target_bits);
        mpfr_mul(rre, buf[src].re, norm, MPFR_RNDN);
        mpfr_mul(rim, buf[src].im, norm, MPFR_RNDN);
        /* Chop imaginary roundoff below the achievable precision so a
         * numerically-real coefficient prints real (as WMA's precision
         * tracking does): |im| <= |value| * 2^-target_bits. */
        if (!mpfr_zero_p(rim)) {
            mpfr_t mag, thr;
            mpfr_init2(mag, target_bits); mpfr_init2(thr, target_bits);
            mpfr_hypot(mag, rre, rim, MPFR_RNDN);
            mpfr_mul_2si(thr, mag, -(long)target_bits, MPFR_RNDN);
            if (mpfr_cmpabs(rim, thr) <= 0) mpfr_set_zero(rim, 1);
            mpfr_clear(mag); mpfr_clear(thr);
        }
        out[o] = numeric_mpfr_make_complex(rre, rim);
        mpfr_clear(rre); mpfr_clear(rim);
    }
    mpfr_clear(norm);
    for (int64_t o = 0; o < N; o++) ncpx_clear(&buf[o]);
    free(buf);

    size_t idx = 0;
    Expr* tree = build_nested(out, dims, rank, 0, &idx);
    free(out);
    return tree;
}

#endif /* USE_MPFR */

/* ==================================================================== *
 *  Symbolic path (exact roots of unity)
 * ==================================================================== */

static Expr* mk_fn(const char* head_sym, Expr** args, size_t n) {
    Expr* r = expr_new_function(expr_new_symbol(head_sym), args, n);
    free(args);
    return r;
}

/* Build Exp[2 Pi I b (coeff/n)], with coeff already carrying the sign. */
static Expr* build_root_of_unity(int64_t coeff, int64_t n, const Expr* b_expr) {
    Expr** ta = malloc(5 * sizeof(Expr*));
    ta[0] = expr_new_integer(2);
    ta[1] = expr_new_symbol(SYM_Pi);
    ta[2] = make_complex(expr_new_integer(0), expr_new_integer(1));   /* I */
    ta[3] = expr_copy((Expr*)b_expr);
    ta[4] = make_rational(coeff, n);
    Expr* prod = mk_fn(SYM_Times, ta, 5);
    Expr** ea = malloc(sizeof(Expr*));
    ea[0] = prod;
    return mk_fn(SYM_Exp, ea, 1);
}

/* Apply the symbolic 1-D transform along one axis, evaluating as we go. */
static void symbolic_axis(Expr** work, const int64_t* dims, int rank,
                          int ax, const int64_t* strides, int64_t N,
                          const Expr* b_expr, int base_sign) {
    int64_t n = dims[ax], s = strides[ax];
    (void)rank;
    if (n <= 1) return;
    for (int64_t o = 0; o < N; o++) {
        if ((o / s) % n != 0) continue;
        Expr** newline = malloc((size_t)n * sizeof(Expr*));
        for (int64_t sidx = 0; sidx < n; sidx++) {
            Expr** terms = malloc((size_t)n * sizeof(Expr*));
            for (int64_t r = 0; r < n; r++) {
                int64_t coeff = (int64_t)base_sign * r * sidx;
                Expr* root = build_root_of_unity(coeff, n, b_expr);
                Expr** ma = malloc(2 * sizeof(Expr*));
                ma[0] = expr_copy(work[o + r * s]);
                ma[1] = root;
                terms[r] = mk_fn(SYM_Times, ma, 2);
            }
            Expr* sum = mk_fn(SYM_Plus, terms, (size_t)n);
            newline[sidx] = eval_and_free(sum);
        }
        for (int64_t r = 0; r < n; r++) {
            expr_free(work[o + r * s]);
            work[o + r * s] = newline[r];
        }
        free(newline);
    }
}

static Expr* symbolic_path(Expr** leaves, const int64_t* dims, int rank,
                           const Expr* a_expr, const Expr* b_expr, int base_sign,
                           int64_t N) {
    int64_t strides[FOURIER_MAX_RANK];
    row_major_strides(dims, rank, strides);

    Expr** work = malloc((size_t)N * sizeof(Expr*));
    for (int64_t o = 0; o < N; o++) work[o] = expr_copy(leaves[o]);

    for (int ax = 0; ax < rank; ax++)
        symbolic_axis(work, dims, rank, ax, strides, N, b_expr, base_sign);

    /* Normalisation exponent: forward (a-1)/2, inverse (-a-1)/2. */
    Expr** pa = malloc(2 * sizeof(Expr*));
    if (base_sign > 0) {
        pa[0] = expr_copy((Expr*)a_expr);
        pa[1] = expr_new_integer(-1);
    } else {
        Expr** na = malloc(2 * sizeof(Expr*));
        na[0] = expr_new_integer(-1);
        na[1] = expr_copy((Expr*)a_expr);
        pa[0] = mk_fn(SYM_Times, na, 2);
        pa[1] = expr_new_integer(-1);
    }
    Expr* plus = mk_fn(SYM_Plus, pa, 2);
    Expr** ea = malloc(2 * sizeof(Expr*));
    ea[0] = make_rational(1, 2);
    ea[1] = plus;
    Expr* expo = mk_fn(SYM_Times, ea, 2);
    Expr** pw = malloc(2 * sizeof(Expr*));
    pw[0] = expr_new_integer((int64_t)N);
    pw[1] = expo;
    Expr* norm = eval_and_free(mk_fn(SYM_Power, pw, 2));

    for (int64_t o = 0; o < N; o++) {
        Expr** ma = malloc(2 * sizeof(Expr*));
        ma[0] = expr_copy(norm);
        ma[1] = work[o];
        work[o] = eval_and_free(mk_fn(SYM_Times, ma, 2));
    }
    expr_free(norm);

    size_t idx = 0;
    Expr* tree = build_nested(work, dims, rank, 0, &idx);
    free(work);
    return tree;
}

/* ==================================================================== *
 *  Option / argument parsing + dispatch
 * ==================================================================== */

/* Is `e` a FourierParameters -> {a, b} option rule? Fills a_out/b_out. */
static bool parse_fourier_option(Expr* e, Expr** a_out, Expr** b_out) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed)
        return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_FourierParameters)
        return false;
    Expr* rhs = e->data.function.args[1];
    if (!head_is(rhs, SYM_List) || rhs->data.function.arg_count != 2) return false;
    *a_out = rhs->data.function.args[0];
    *b_out = rhs->data.function.args[1];
    return true;
}

/* Compute the transform of `data` (a validated nested list / NDArray). */
static Expr* fourier_compute(Expr* data, Expr* a_expr, Expr* b_expr, int inverse) {
    int base_sign_nd = inverse ? -1 : +1;

    /* NDArray fast path: an NDArray is always machine-precision numeric, so it
     * bypasses the List round trip entirely and returns an NDArray directly. */
    if (data->type == EXPR_NDARRAY) {
        if (data->data.ndarray.rank > FOURIER_MAX_RANK) return NULL;
        double a;
        long b;
        if (!expr_to_double_num(a_expr, &a) || !expr_to_long_int(b_expr, &b))
            return NULL;   /* non-integer b: stay unevaluated */
        return machine_path_ndarray(data, a, b, base_sign_nd);
    }

    int64_t dims[FOURIER_MAX_RANK];
    int rank = 0;
    if (!nested_dims(data, dims, &rank)) return NULL;
    int64_t N = 1;
    for (int i = 0; i < rank; i++) {
        if (dims[i] <= 0) return NULL;   /* reject empty */
        N *= dims[i];
    }

    Expr** leaves = malloc((size_t)N * sizeof(Expr*));
    size_t idx = 0;
    nd_flatten(data, 0, rank, leaves, &idx);

    long arb_bits = 0;
    Regime reg = classify(leaves, (size_t)N, &arb_bits);
    int base_sign = base_sign_nd;

    Expr* result = NULL;
    if (reg == REG_SYMBOLIC) {
        result = symbolic_path(leaves, dims, rank, a_expr, b_expr, base_sign, N);
    } else {
        double a;
        long b;
        if (!expr_to_double_num(a_expr, &a) || !expr_to_long_int(b_expr, &b)) {
            /* Non-integer b (or unreadable a) on numeric data: stay unevaluated. */
            free(leaves); return NULL;
        }
        if (reg == REG_MACHINE) {
            result = machine_path(data, dims, rank, a, b, base_sign);
        } else { /* REG_ARB */
#ifdef USE_MPFR
            result = arb_path(leaves, dims, rank, a, b, base_sign, arb_bits);
#else
            result = NULL;
#endif
        }
    }

    free(leaves);
    return result;
}

static Expr* fourier_core(Expr* res, int inverse) {
    const char* name = inverse ? "InverseFourier" : "Fourier";
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return builtin_arg_error(name, argc, 1, 2);

    Expr* data = res->data.function.args[0];

    /* Default {a, b} = {0, 1}; scan trailing args for the option and an
     * optional position list. */
    Expr* a_default = expr_new_integer(0);
    Expr* b_default = expr_new_integer(1);
    Expr* a_expr = a_default;
    Expr* b_expr = b_default;
    Expr* poslist = NULL;

    for (size_t i = 1; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        Expr *ao, *bo;
        if (parse_fourier_option(arg, &ao, &bo)) {
            a_expr = ao;
            b_expr = bo;
        } else if (head_is(arg, SYM_List) && poslist == NULL) {
            poslist = arg;                         /* position-selection form */
        } else {
            expr_free(a_default); expr_free(b_default);
            return NULL;                           /* unrecognised argument */
        }
    }

    Expr* full = fourier_compute(data, a_expr, b_expr, inverse);
    expr_free(a_default);
    expr_free(b_default);
    if (!full) return NULL;

    if (poslist == NULL) return full;

    /* Position form: Extract[full, poslist]. Delist a packed NDArray first so
     * Extract's list-of-positions form applies (position selection yields a
     * small list, so packing offers nothing here). */
    if (full->type == EXPR_NDARRAY) {
        Expr* l = ndarray_to_nested_list(full);
        expr_free(full);
        full = l;
    }
    Expr** xa = malloc(2 * sizeof(Expr*));
    xa[0] = full;
    xa[1] = expr_copy(poslist);
    Expr* xcall = expr_new_function(expr_new_symbol("Extract"), xa, 2);
    free(xa);
    return eval_and_free(xcall);
}

Expr* builtin_fourier(Expr* res)         { return fourier_core(res, 0); }
Expr* builtin_inverse_fourier(Expr* res) { return fourier_core(res, 1); }

/* ==================================================================== *
 *  FourierDCT / FourierDST — real discrete cosine / sine transforms
 * ====================================================================
 *
 * Four WMA types each (I..IV), selectable by integer 1..4 or the strings
 * "I".."IV"; the default is type II. Each is the real orthonormal (unitary)
 * transform v = M u with the real matrix (0-indexed i = s-1, j = r-1; length n):
 *
 *   DCT-I  (n>=2): sqrt(2/(n-1)) c_j cos(pi i j/(n-1)),  c_0=c_{n-1}=1/2 else 1
 *   DCT-II       : (1/sqrt(n)) cos(pi (j+1/2) i / n)
 *   DCT-III      : (1/sqrt(n)) d_j cos(pi j (i+1/2)/n),  d_0=1 else 2
 *   DCT-IV       : sqrt(2/n) cos(pi (j+1/2)(i+1/2)/n)
 *   DST-I        : sqrt(2/(n+1)) sin(pi (j+1)(i+1)/(n+1))
 *   DST-II       : (1/sqrt(n)) sin(pi (j+1/2)(i+1)/n)
 *   DST-III      : (1/sqrt(n)) ( 2 sin(pi (j+1)(i+1/2)/n) for j<n-1; (-1)^i for j=n-1 )
 *   DST-IV       : sqrt(2/n) sin(pi (j+1/2)(i+1/2)/n)
 *
 * The matrix is real and separable, so the transform is applied independently
 * per axis and real input yields real output. Complex input is handled by
 * applying M to the real and imaginary parts independently.
 *
 * Regimes mirror Fourier: exact input is numericalised with N and takes the
 * machine path (FFTW real-to-real REDFT/RODFT plan — O(n log n) — or a direct
 * O(n^2) matrix fallback); MPFR input takes a direct O(n^2) MPFR matrix path;
 * genuinely symbolic input stays unevaluated (WMA behaviour).
 */

/* Parse the type argument: Integer 1..4 or the strings "I".."IV". */
static bool parse_dct_type(const Expr* e, int* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        int64_t v = e->data.integer;
        if (v >= 1 && v <= 4) { *out = (int)v; return true; }
        return false;
    }
    if (e->type == EXPR_STRING) {
        const char* s = e->data.string;
        if (!s) return false;
        if (!strcmp(s, "I"))   { *out = 1; return true; }
        if (!strcmp(s, "II"))  { *out = 2; return true; }
        if (!strcmp(s, "III")) { *out = 3; return true; }
        if (!strcmp(s, "IV"))  { *out = 4; return true; }
    }
    return false;
}

/* Axis-length validity: DCT-I needs n >= 2 (it divides by n-1); all else n >= 1. */
static bool dct_axis_ok(int type, bool sine, int64_t n) {
    if (n < 1) return false;
    if (!sine && type == 1 && n < 2) return false;
    return true;
}

#ifdef USE_FFTW
/* Map (type, sine) to the FFTW real-to-real transform kind. */
static fftw_r2r_kind dct_fftw_kind(int type, bool sine) {
    if (!sine) {
        switch (type) { case 1: return FFTW_REDFT00; case 2: return FFTW_REDFT10;
                        case 3: return FFTW_REDFT01; default: return FFTW_REDFT11; }
    }
    switch (type) { case 1: return FFTW_RODFT00; case 2: return FFTW_RODFT10;
                    case 3: return FFTW_RODFT01; default: return FFTW_RODFT11; }
}

/* Per-axis scale mapping FFTW's unnormalised REDFT/RODFT output onto the WMA
 * orthonormal convention (the total scale is the product over axes). */
static double dct_axis_scale(int type, bool sine, int64_t n) {
    double dn = (double)n;
    switch (type) {
        case 1: return sine ? 1.0 / sqrt(2.0 * (dn + 1.0)) : 1.0 / sqrt(2.0 * (dn - 1.0));
        case 2: return 1.0 / (2.0 * sqrt(dn));
        case 3: return 1.0 / sqrt(dn);
        default: return 1.0 / sqrt(2.0 * dn);
    }
}
#else
/* Single matrix entry M[i][j] (machine double) per the formulas above. */
static double dct_entry(int i, int j, int n, int type, bool sine) {
    double dn = (double)n, di = (double)i, dj = (double)j;
    if (!sine) {
        switch (type) {
            case 1: { double c = (j == 0 || j == n - 1) ? 0.5 : 1.0;
                      return sqrt(2.0 / (dn - 1.0)) * c * cos(M_PI * di * dj / (dn - 1.0)); }
            case 2: return (1.0 / sqrt(dn)) * cos(M_PI * (dj + 0.5) * di / dn);
            case 3: { double d = (j == 0) ? 1.0 : 2.0;
                      return (1.0 / sqrt(dn)) * d * cos(M_PI * dj * (di + 0.5) / dn); }
            case 4: return sqrt(2.0 / dn) * cos(M_PI * (dj + 0.5) * (di + 0.5) / dn);
        }
    } else {
        switch (type) {
            case 1: return sqrt(2.0 / (dn + 1.0)) * sin(M_PI * (dj + 1.0) * (di + 1.0) / (dn + 1.0));
            case 2: return (1.0 / sqrt(dn)) * sin(M_PI * (dj + 0.5) * (di + 1.0) / dn);
            case 3: if (j == n - 1) return (1.0 / sqrt(dn)) * ((i % 2 == 0) ? 1.0 : -1.0);
                    return (1.0 / sqrt(dn)) * 2.0 * sin(M_PI * (dj + 1.0) * (di + 0.5) / dn);
            case 4: return sqrt(2.0 / dn) * sin(M_PI * (dj + 0.5) * (di + 0.5) / dn);
        }
    }
    return 0.0;
}

/* Direct O(n^2) matrix-vector transform along one axis (strided lines). The
 * matrix already carries the full WMA normalisation. */
static void dct_apply_axis_machine(double* data, const int64_t* dims, int rank,
                                   int ax, int type, bool sine) {
    int64_t strides[FOURIER_MAX_RANK];
    row_major_strides(dims, rank, strides);
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    int64_t n = dims[ax], s = strides[ax];
    double* M = malloc((size_t)n * (size_t)n * sizeof(double));
    for (int64_t i = 0; i < n; i++)
        for (int64_t j = 0; j < n; j++)
            M[i * n + j] = dct_entry((int)i, (int)j, (int)n, type, sine);
    double* tmp  = malloc((size_t)n * sizeof(double));
    double* outl = malloc((size_t)n * sizeof(double));
    for (int64_t o = 0; o < N; o++) {
        if ((o / s) % n != 0) continue;
        for (int64_t r = 0; r < n; r++) tmp[r] = data[o + r * s];
        for (int64_t i = 0; i < n; i++) {
            double acc = 0.0;
            for (int64_t j = 0; j < n; j++) acc += M[i * n + j] * tmp[j];
            outl[i] = acc;
        }
        for (int64_t i = 0; i < n; i++) data[o + i * s] = outl[i];
    }
    free(M); free(tmp); free(outl);
}
#endif

/* Consume the interleaved (re, im) buffer `buf` (length 2N) and return a fresh
 * interleaved output buffer, or NULL if any axis length is invalid for `type`. */
static double* dct_transform_buf(double* buf, const int64_t* dims, int rank,
                                 int type, bool sine) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) {
        if (!dct_axis_ok(type, sine, dims[i])) { free(buf); return NULL; }
        N *= dims[i];
    }
    double* re = malloc((size_t)N * sizeof(double));
    double* im = malloc((size_t)N * sizeof(double));
    bool has_im = false;
    for (int64_t o = 0; o < N; o++) {
        re[o] = buf[2 * o];
        im[o] = buf[2 * o + 1];
        if (im[o] != 0.0) has_im = true;
    }
    free(buf);

    double scale = 1.0;
#ifdef USE_FFTW
    int idims[FOURIER_MAX_RANK];
    fftw_r2r_kind kind[FOURIER_MAX_RANK];
    fftw_r2r_kind k = dct_fftw_kind(type, sine);
    for (int i = 0; i < rank; i++) { idims[i] = (int)dims[i]; kind[i] = k; }
    for (int i = 0; i < rank; i++) scale *= dct_axis_scale(type, sine, dims[i]);
    fftw_plan p = fftw_plan_r2r(rank, idims, re, re, kind, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    if (has_im) {
        fftw_plan p2 = fftw_plan_r2r(rank, idims, im, im, kind, FFTW_ESTIMATE);
        fftw_execute(p2);
        fftw_destroy_plan(p2);
    }
#else
    for (int ax = 0; ax < rank; ax++) {
        dct_apply_axis_machine(re, dims, rank, ax, type, sine);
        if (has_im) dct_apply_axis_machine(im, dims, rank, ax, type, sine);
    }
#endif

    double* out = malloc((size_t)N * 2 * sizeof(double));
    for (int64_t o = 0; o < N; o++) {
        out[2 * o]     = scale * re[o];
        out[2 * o + 1] = has_im ? scale * im[o] : 0.0;
    }
    free(re); free(im);
    return out;
}

/* List path: numericalise, pack, transform, delist to a plain nested List. */
static Expr* dct_machine_path(const Expr* nested, const int64_t* dims, int rank,
                              int type, bool sine) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    if (N <= 0) return NULL;

    Expr* nnum = numericalize(nested, numeric_machine_spec());
    if (!nnum) return NULL;
    int64_t d2[FOURIER_MAX_RANK];
    int r2 = 0;
    if (!nested_dims(nnum, d2, &r2) || r2 != rank) { expr_free(nnum); return NULL; }
    for (int i = 0; i < rank; i++) if (d2[i] != dims[i]) { expr_free(nnum); return NULL; }

    Expr** leaves = malloc((size_t)N * sizeof(Expr*));
    size_t idx = 0;
    nd_flatten(nnum, 0, rank, leaves, &idx);
    double* buf = malloc((size_t)N * 2 * sizeof(double));
    for (int64_t o = 0; o < N; o++) {
        double re, im;
        if (!leaf_to_cdouble(leaves[o], &re, &im)) {
            free(leaves); free(buf); expr_free(nnum); return NULL;
        }
        buf[2 * o] = re;
        buf[2 * o + 1] = im;
    }
    free(leaves);
    expr_free(nnum);

    double* out = dct_transform_buf(buf, dims, rank, type, sine);
    if (!out) return NULL;
    Expr* nd = machine_build_ndarray(out, N, dims, rank);
    Expr* lst = ndarray_to_nested_list(nd);
    expr_free(nd);
    return lst;
}

/* NDArray fast path: read the packed buffer directly and return an NDArray. */
static Expr* dct_machine_path_ndarray(const Expr* nd, int type, bool sine) {
    int rank = nd->data.ndarray.rank;
    const int64_t* dims = nd->data.ndarray.dims;
    NDType dt = nd->data.ndarray.dtype;
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    if (N <= 0) return NULL;

    double* buf = malloc((size_t)N * 2 * sizeof(double));
    for (int64_t o = 0; o < N; o++)
        ndt_get(nd->data.ndarray.data, (size_t)o, dt, &buf[2 * o], &buf[2 * o + 1]);

    double* out = dct_transform_buf(buf, dims, rank, type, sine);
    if (!out) return NULL;
    return machine_build_ndarray(out, N, dims, rank);
}

#ifdef USE_MPFR
/* Build the n*n MPFR transform matrix at working precision `wp`. Each entry is
 * mpfr_init2'd; the caller clears and frees the array. */
static mpfr_t* dct_matrix_build_mpfr(int n, int type, bool sine, mpfr_prec_t wp) {
    mpfr_t* M = malloc((size_t)n * (size_t)n * sizeof(mpfr_t));
    mpfr_t base, pi, ang, val;
    mpfr_init2(base, wp); mpfr_init2(pi, wp); mpfr_init2(ang, wp); mpfr_init2(val, wp);
    mpfr_const_pi(pi, MPFR_RNDN);

    /* Base coefficient (common to every entry of this type). */
    if (!sine) {
        switch (type) {
            case 1: mpfr_set_ui(base, 2, MPFR_RNDN); mpfr_div_ui(base, base, (unsigned long)(n - 1), MPFR_RNDN); mpfr_sqrt(base, base, MPFR_RNDN); break;
            case 4: mpfr_set_ui(base, 2, MPFR_RNDN); mpfr_div_ui(base, base, (unsigned long)n, MPFR_RNDN); mpfr_sqrt(base, base, MPFR_RNDN); break;
            default: mpfr_set_ui(base, 1, MPFR_RNDN); { mpfr_t t; mpfr_init2(t, wp); mpfr_sqrt_ui(t, (unsigned long)n, MPFR_RNDN); mpfr_div(base, base, t, MPFR_RNDN); mpfr_clear(t); } break; /* 1/sqrt(n) */
        }
    } else {
        switch (type) {
            case 1: mpfr_set_ui(base, 2, MPFR_RNDN); mpfr_div_ui(base, base, (unsigned long)(n + 1), MPFR_RNDN); mpfr_sqrt(base, base, MPFR_RNDN); break;
            case 4: mpfr_set_ui(base, 2, MPFR_RNDN); mpfr_div_ui(base, base, (unsigned long)n, MPFR_RNDN); mpfr_sqrt(base, base, MPFR_RNDN); break;
            default: mpfr_set_ui(base, 1, MPFR_RNDN); { mpfr_t t; mpfr_init2(t, wp); mpfr_sqrt_ui(t, (unsigned long)n, MPFR_RNDN); mpfr_div(base, base, t, MPFR_RNDN); mpfr_clear(t); } break;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mpfr_init2(M[i * n + j], wp);
            unsigned long num = 0, den = 1;
            int use_cos = 1;
            double factor = 1.0;
            int special = 0;   /* DST-III last column: value = base * (-1)^i */
            if (!sine) {
                switch (type) {
                    case 1: num = (unsigned long)i * (unsigned long)j; den = (unsigned long)(n - 1);
                            factor = (j == 0 || j == n - 1) ? 0.5 : 1.0; break;
                    case 2: num = (unsigned long)(2 * j + 1) * (unsigned long)i; den = (unsigned long)(2 * n); break;
                    case 3: num = (unsigned long)j * (unsigned long)(2 * i + 1); den = (unsigned long)(2 * n);
                            factor = (j == 0) ? 1.0 : 2.0; break;
                    case 4: num = (unsigned long)(2 * j + 1) * (unsigned long)(2 * i + 1); den = (unsigned long)(4 * n); break;
                }
            } else {
                use_cos = 0;
                switch (type) {
                    case 1: num = (unsigned long)(j + 1) * (unsigned long)(i + 1); den = (unsigned long)(n + 1); break;
                    case 2: num = (unsigned long)(2 * j + 1) * (unsigned long)(i + 1); den = (unsigned long)(2 * n); break;
                    case 3: if (j == n - 1) { special = 1; }
                            else { num = (unsigned long)(j + 1) * (unsigned long)(2 * i + 1); den = (unsigned long)(2 * n); factor = 2.0; }
                            break;
                    case 4: num = (unsigned long)(2 * j + 1) * (unsigned long)(2 * i + 1); den = (unsigned long)(4 * n); break;
                }
            }
            if (special) {
                mpfr_set(M[i * n + j], base, MPFR_RNDN);
                if (i % 2 != 0) mpfr_neg(M[i * n + j], M[i * n + j], MPFR_RNDN);
                continue;
            }
            mpfr_mul_ui(ang, pi, num, MPFR_RNDN);
            mpfr_div_ui(ang, ang, den, MPFR_RNDN);
            if (use_cos) mpfr_cos(val, ang, MPFR_RNDN);
            else         mpfr_sin(val, ang, MPFR_RNDN);
            mpfr_mul(M[i * n + j], base, val, MPFR_RNDN);
            if (factor != 1.0) mpfr_mul_d(M[i * n + j], M[i * n + j], factor, MPFR_RNDN);
        }
    }
    mpfr_clear(base); mpfr_clear(pi); mpfr_clear(ang); mpfr_clear(val);
    return M;
}

static Expr* dct_arb_path(Expr** leaves, const int64_t* dims, int rank,
                          int type, bool sine, long target_bits) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) {
        if (!dct_axis_ok(type, sine, dims[i])) return NULL;
        N *= dims[i];
    }
    if (N <= 0) return NULL;

    long guard = 32;
    { int64_t t = N; while (t > 1) { guard += 2; t >>= 1; } }
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);

    mpfr_t* re = malloc((size_t)N * sizeof(mpfr_t));
    mpfr_t* im = malloc((size_t)N * sizeof(mpfr_t));
    for (int64_t o = 0; o < N; o++) {
        mpfr_init2(re[o], wp); mpfr_init2(im[o], wp);
        mpfr_t rr, ii; mpfr_init2(rr, wp); mpfr_init2(ii, wp);
        bool inexact = false;
        if (!get_approx_mpfr(leaves[o], rr, ii, &inexact)) {
            mpfr_clear(rr); mpfr_clear(ii);
            for (int64_t j = 0; j <= o; j++) { mpfr_clear(re[j]); mpfr_clear(im[j]); }
            free(re); free(im);
            return NULL;
        }
        mpfr_set(re[o], rr, MPFR_RNDN); mpfr_set(im[o], ii, MPFR_RNDN);
        mpfr_clear(rr); mpfr_clear(ii);
    }

    int64_t strides[FOURIER_MAX_RANK];
    row_major_strides(dims, rank, strides);
    for (int ax = 0; ax < rank; ax++) {
        int64_t n = dims[ax], s = strides[ax];
        mpfr_t* M = dct_matrix_build_mpfr((int)n, type, sine, wp);
        mpfr_t *lre = malloc((size_t)n * sizeof(mpfr_t)), *lim = malloc((size_t)n * sizeof(mpfr_t));
        mpfr_t *ore = malloc((size_t)n * sizeof(mpfr_t)), *oim = malloc((size_t)n * sizeof(mpfr_t));
        mpfr_t prod; mpfr_init2(prod, wp);
        for (int64_t r = 0; r < n; r++) { mpfr_init2(lre[r], wp); mpfr_init2(lim[r], wp); mpfr_init2(ore[r], wp); mpfr_init2(oim[r], wp); }
        for (int64_t o = 0; o < N; o++) {
            if ((o / s) % n != 0) continue;
            for (int64_t r = 0; r < n; r++) { mpfr_set(lre[r], re[o + r * s], MPFR_RNDN); mpfr_set(lim[r], im[o + r * s], MPFR_RNDN); }
            for (int64_t i = 0; i < n; i++) {
                mpfr_set_zero(ore[i], 1); mpfr_set_zero(oim[i], 1);
                for (int64_t j = 0; j < n; j++) {
                    mpfr_mul(prod, M[i * n + j], lre[j], MPFR_RNDN); mpfr_add(ore[i], ore[i], prod, MPFR_RNDN);
                    mpfr_mul(prod, M[i * n + j], lim[j], MPFR_RNDN); mpfr_add(oim[i], oim[i], prod, MPFR_RNDN);
                }
            }
            for (int64_t i = 0; i < n; i++) { mpfr_set(re[o + i * s], ore[i], MPFR_RNDN); mpfr_set(im[o + i * s], oim[i], MPFR_RNDN); }
        }
        for (int64_t r = 0; r < n; r++) { mpfr_clear(lre[r]); mpfr_clear(lim[r]); mpfr_clear(ore[r]); mpfr_clear(oim[r]); }
        free(lre); free(lim); free(ore); free(oim);
        mpfr_clear(prod);
        for (int64_t t = 0; t < n * n; t++) mpfr_clear(M[t]);
        free(M);
    }

    Expr** out = malloc((size_t)N * sizeof(Expr*));
    for (int64_t o = 0; o < N; o++) {
        mpfr_t rre, rim; mpfr_init2(rre, target_bits); mpfr_init2(rim, target_bits);
        mpfr_set(rre, re[o], MPFR_RNDN); mpfr_set(rim, im[o], MPFR_RNDN);
        out[o] = numeric_mpfr_make_complex(rre, rim);
        mpfr_clear(rre); mpfr_clear(rim);
    }
    for (int64_t o = 0; o < N; o++) { mpfr_clear(re[o]); mpfr_clear(im[o]); }
    free(re); free(im);

    size_t idx = 0;
    Expr* tree = build_nested(out, dims, rank, 0, &idx);
    free(out);
    return tree;
}
#endif /* USE_MPFR */

/* Dispatch a validated data argument to the appropriate regime. */
static Expr* dct_compute(Expr* data, int type, bool sine) {
    if (data->type == EXPR_NDARRAY) {
        if (data->data.ndarray.rank > FOURIER_MAX_RANK) return NULL;
        return dct_machine_path_ndarray(data, type, sine);
    }

    int64_t dims[FOURIER_MAX_RANK];
    int rank = 0;
    if (!nested_dims(data, dims, &rank)) return NULL;
    int64_t N = 1;
    for (int i = 0; i < rank; i++) {
        if (dims[i] <= 0) return NULL;
        N *= dims[i];
    }

    Expr** leaves = malloc((size_t)N * sizeof(Expr*));
    size_t idx = 0;
    nd_flatten(data, 0, rank, leaves, &idx);

    long arb_bits = 0;
    Regime reg = classify(leaves, (size_t)N, &arb_bits);

    Expr* result = NULL;
    if (reg == REG_SYMBOLIC) {
        result = NULL;                       /* symbolic input stays unevaluated */
    } else if (reg == REG_MACHINE) {
        result = dct_machine_path(data, dims, rank, type, sine);
    } else { /* REG_ARB */
#ifdef USE_MPFR
        result = dct_arb_path(leaves, dims, rank, type, sine, arb_bits);
#else
        result = NULL;
#endif
    }
    free(leaves);
    return result;
}

static Expr* dct_core(Expr* res, bool sine) {
    const char* name = sine ? "FourierDST" : "FourierDCT";
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error(name, argc, 1, 2);

    Expr* data = res->data.function.args[0];
    int type = 2;                            /* FourierDCT[list] == type II */
    if (argc == 2) {
        if (!parse_dct_type(res->data.function.args[1], &type)) return NULL;
    }
    return dct_compute(data, type, sine);
}

Expr* builtin_fourier_dct(Expr* res) { return dct_core(res, false); }
Expr* builtin_fourier_dst(Expr* res) { return dct_core(res, true);  }

/* ==================================================================== *
 *  Registration
 * ==================================================================== */

void fourier_init(void) {
    symtab_add_builtin("Fourier", builtin_fourier);
    symtab_add_builtin("InverseFourier", builtin_inverse_fourier);
    symtab_add_builtin("FourierDCT", builtin_fourier_dct);
    symtab_add_builtin("FourierDST", builtin_fourier_dst);
    symtab_get_def("Fourier")->attributes |= ATTR_PROTECTED;
    symtab_get_def("InverseFourier")->attributes |= ATTR_PROTECTED;
    symtab_get_def("FourierDCT")->attributes |= ATTR_PROTECTED;
    symtab_get_def("FourierDST")->attributes |= ATTR_PROTECTED;
    /* FourierParameters is an option symbol. */
    symtab_get_def("FourierParameters")->attributes |= ATTR_PROTECTED;
}
