/* ndinteger.c — exact-integer NDArray kernels for the number-theoretic heads.
 * See ndinteger.h for what separates these from ndkernels.c's `double` set. */

#include "ndinteger.h"
#include "ndarray.h"
#include "expr.h"
#include "pack.h"
#include "symtab.h"
#include "sym_names.h"
#include "checked_int.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 *  int64 factorisation, by trial division.
 *
 *  The scalar builtins factor through GMP, which is the right choice for one
 *  number of unbounded size and the wrong one for 200000 of them: the per-call
 *  mpz_init/clear alone dominates. Stripping the factors of 2 and then dividing
 *  by the odd numbers needs no allocation at all and settles a 6-digit element
 *  in well under a microsecond.
 *
 *  It is also why there is a size ceiling. Trial division is O(sqrt(n)), so a
 *  single 10^18 element would cost 10^9 iterations and be far slower than the
 *  GMP path it replaced. NDI_FACTOR_MAX caps sqrt(n) at 10^6 iterations; past
 *  it the kernel declines, which abandons the whole array and lets the List
 *  path answer through GMP exactly as before. A ceiling that turns a fast path
 *  off is a performance decision; the ANSWER is identical either way.
 * -------------------------------------------------------------------------- */
#define NDI_FACTOR_MAX ((int64_t)1000000000000LL)   /* 10^12; sqrt = 10^6 */
#define NDI_MAX_FACTORS 64                          /* 2^64 > any int64 */

/* Factor |n| (n != 0) into distinct primes and exponents. Returns false when n
 * is 0 or beyond NDI_FACTOR_MAX. n == 1 gives the empty factorisation. */
static bool ndi_factor(int64_t n, int64_t* p, int* e, int* np) {
    if (n == 0) return false;
    int64_t m;
    if (ci_abs_i64(n, &m)) return false;            /* INT64_MIN has no |n| */
    if (m > NDI_FACTOR_MAX) return false;
    int k = 0;
    if ((m & 1) == 0) {
        p[k] = 2; e[k] = 0;
        while ((m & 1) == 0) { e[k]++; m >>= 1; }
        k++;
    }
    for (int64_t d = 3; d * d <= m; d += 2) {
        if (m % d) continue;
        if (k >= NDI_MAX_FACTORS) return false;
        p[k] = d; e[k] = 0; k++;
        while (m % d == 0) { e[k-1]++; m /= d; }
    }
    if (m > 1) {
        if (k >= NDI_MAX_FACTORS) return false;
        p[k] = m; e[k] = 1; k++;
    }
    *np = k;
    return true;
}

/* ---- unary integer kernels ---------------------------------------------- */

/* IntegerLength[n]: decimal digits of |n|. IntegerLength[0] is 0, matching
 * builtin_integerlength (and Mathematica). */
static bool ndk_IntegerLength_ii(int64_t n, int64_t* out) {
    int64_t m;
    if (ci_abs_i64(n, &m)) return false;
    int64_t d = 0;
    while (m > 0) { d++; m /= 10; }
    *out = d;
    return true;
}

/* EulerPhi[n] = |n| * prod (1 - 1/p). EulerPhi[0] is 0 and the sign of n is
 * ignored, both matching the scalar builtin. Divide before multiplying so the
 * running value never exceeds |n| and cannot overflow. */
static bool ndk_EulerPhi_ii(int64_t n, int64_t* out) {
    if (n == 0) { *out = 0; return true; }
    int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS]; int np;
    if (!ndi_factor(n, p, e, &np)) return false;
    int64_t m;
    if (ci_abs_i64(n, &m)) return false;
    int64_t phi = m;
    for (int i = 0; i < np; i++) phi = phi / p[i] * (p[i] - 1);
    *out = phi;
    return true;
}

/* MoebiusMu[n]: 0 if any exponent >= 2, else (-1)^(distinct primes). Declines
 * n == 0, where the scalar builtin leaves the call unevaluated -- so the whole
 * array degrades and the List path reproduces that exactly. Sign ignored. */
static bool ndk_MoebiusMu_ii(int64_t n, int64_t* out) {
    if (n == 0) return false;
    int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS]; int np;
    if (!ndi_factor(n, p, e, &np)) return false;
    for (int i = 0; i < np; i++) if (e[i] >= 2) { *out = 0; return true; }
    *out = (np % 2 == 0) ? 1 : -1;
    return true;
}

/* ---- binary integer kernels --------------------------------------------- */

static int64_t ndi_gcd_u(int64_t a, int64_t b) {
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* GCD[a, b] over the non-negative representatives: GCD[-12, 18] is 6 and
 * GCD[0, 5] is 5, both matching the scalar builtin. */
static bool ndk_GCD_ii(int64_t a, int64_t b, int64_t* out) {
    int64_t x, y;
    if (ci_abs_i64(a, &x) || ci_abs_i64(b, &y)) return false;
    *out = ndi_gcd_u(x, y);
    return true;
}

/* LCM[a, b] = |a|/gcd * |b|, in that order so the product is formed from the
 * already-reduced factor. LCM[0, n] is 0. Overflow declines. */
static bool ndk_LCM_ii(int64_t a, int64_t b, int64_t* out) {
    int64_t x, y;
    if (ci_abs_i64(a, &x) || ci_abs_i64(b, &y)) return false;
    if (x == 0 || y == 0) { *out = 0; return true; }
    int64_t g = ndi_gcd_u(x, y);
    return !ci_mul_i64(x / g, y, out);
}

/* DivisorSigma[k, n] = prod_i (1 + p^k + ... + p^(e_i k)).
 *
 * Argument order is the CALL order, which ndarray_map_binary composes for
 * either operand position -- so `k` is always the first argument and `n` the
 * second, whichever of the two was the array. Declines n == 0 (unevaluated in
 * the scalar builtin) and k < 0 (sigma_-1[6] is 2, a Rational, which no int64
 * buffer holds). The sign of n is ignored, matching the builtin. */
static bool ndk_DivisorSigma_ii(int64_t k, int64_t n, int64_t* out) {
    if (n == 0 || k < 0) return false;
    int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS]; int np;
    if (!ndi_factor(n, p, e, &np)) return false;
    int64_t total = 1;
    for (int i = 0; i < np; i++) {
        int64_t term = 1, pk = 1;
        for (int j = 0; j < e[i]; j++) {
            int64_t step;
            if (ci_powi_i64(p[i], k, &step)) return false;
            if (ci_mul_i64(pk, step, &pk)) return false;
            if (ci_add_i64(term, pk, &term)) return false;
        }
        if (ci_mul_i64(total, term, &total)) return false;
    }
    *out = total;
    return true;
}

/* Integer-only descriptors: no `cplx`, no `real`, `real_closed` false. Every
 * non-integer branch of ndarray_map_unary / _binary tests `cplx` and declines,
 * so a Real buffer degrades to the List path — which is what these heads want
 * anyway (EulerPhi[3.] is not EulerPhi[3]). */
static const NDUnaryKernel NDKU_IntegerLength =
    { NULL, NULL, false, false, NULL, ndk_IntegerLength_ii, true };
static const NDUnaryKernel NDKU_EulerPhi =
    { NULL, NULL, false, false, NULL, ndk_EulerPhi_ii, true };
static const NDUnaryKernel NDKU_MoebiusMu =
    { NULL, NULL, false, false, NULL, ndk_MoebiusMu_ii, true };
static const NDBinaryKernel NDKB_GCD =
    { NULL, false, NULL, ndk_GCD_ii, true };
static const NDBinaryKernel NDKB_LCM =
    { NULL, false, NULL, ndk_LCM_ii, true };
static const NDBinaryKernel NDKB_DivisorSigma =
    { NULL, false, NULL, ndk_DivisorSigma_ii, true };

/* ---------------------------------------------------------------------------
 *  Head-level buffer paths.
 * -------------------------------------------------------------------------- */

/* The sole argument, when it is a rank-1-or-higher int64 NDArray. */
static const Expr* ndi_int_arg(const Expr* res, size_t which, size_t argc) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != argc) return NULL;
    const Expr* a = res->data.function.args[which];
    if (!is_ndarray(a) || a->data.ndarray.dtype != NDT_INT64) return NULL;
    return a;
}

/* An exact non-negative int64 from an Integer argument. */
static bool ndi_exact_int(const Expr* e, int64_t* out) {
    if (!e || e->type != EXPR_INTEGER) return false;
    *out = e->data.integer;
    return true;
}

/* ---- Prime[arr] --------------------------------------------------------- *
 * One sieve for the whole array, not one prime count per element.
 *
 * Prime[n] alone is a Meissel/Lucy count, which is the right algorithm for a
 * single large index and the wrong one repeated 200000 times. Given the largest
 * index in the array, a sieve of Eratosthenes to Rosser's bound
 *     p_n < n (ln n + ln ln n)      (n >= 6)
 * produces every requested prime in one linear pass, and the gather is an
 * indexed read. The array as a whole then costs about what its largest single
 * element used to.
 *
 * Declines an index below 1 (Prime[0] is unevaluated, with a message, in the
 * scalar builtin) and an upper bound past what a byte sieve should allocate,
 * where the builtin's own path is the better one anyway. */
/* THE SIEVE IS NOT UNCONDITIONALLY BETTER, and the condition is the array's
 * LENGTH, not its values. One sieve to `limit` costs O(limit log log limit)
 * whatever the array holds, where the scalar path costs one Meissel count per
 * element -- so for a two-element array asking for a large index the builtin's
 * own path wins, and taking the sieve anyway would be a regression paid in
 * hundreds of megabytes.
 *
 *   under NDI_SIEVE_FREE   a few milliseconds and a few MB. Cheaper than a
 *                          single Meissel count at that size, so no amount of
 *                          elements is needed to justify it.
 *   up to NDI_SIEVE_MAX    only once there are enough elements to amortise it.
 *   past NDI_SIEVE_MAX     never; decline, and each element takes the Meissel
 *                          count exactly as before.
 */
#define NDI_SIEVE_FREE  ((int64_t)4000000)    /* ~4 MB, always worth it */
#define NDI_SIEVE_MAX   ((int64_t)50000000)   /* ~50 MB, the hard ceiling */
#define NDI_SIEVE_AMORT ((size_t)64)          /* elements needed past _FREE */

Expr* ndint_prime(Expr* res) {
    const Expr* a = ndi_int_arg(res, 0, 1);
    if (!a) return NULL;
    size_t sz = ndarray_size(a);
    if (sz == 0) return NULL;
    const void* in = a->data.ndarray.data;

    int64_t nmax = 0;
    for (size_t i = 0; i < sz; i++) {
        int64_t v = ndt_get_i(in, i, NDT_INT64);
        if (v < 1) return NULL;                    /* Prime[0] stays symbolic */
        if (v > nmax) nmax = v;
    }

    /* Rosser & Schoenfeld's upper bound on the n-th prime, with the small-n
     * cases (where ln ln n is <= 0) handled by a flat floor. */
    double dn = (double)nmax;
    int64_t limit;
    if (nmax < 6) limit = 13;
    else {
        double ln = log(dn), lnln = log(ln);
        limit = (int64_t)(dn * (ln + lnln)) + 16;
    }
    if (limit > NDI_SIEVE_MAX) return NULL;
    if (limit > NDI_SIEVE_FREE && sz < NDI_SIEVE_AMORT) return NULL;

    unsigned char* comp = calloc((size_t)limit + 1, 1);
    if (!comp) return NULL;
    for (int64_t i = 2; i * i <= limit; i++)
        if (!comp[i]) for (int64_t j = i * i; j <= limit; j += i) comp[j] = 1;

    /* Index the primes in order; primes[k] is the (k+1)-th prime. */
    int64_t* primes = malloc(sizeof(int64_t) * (size_t)nmax);
    if (!primes) { free(comp); return NULL; }
    int64_t found = 0;
    for (int64_t i = 2; i <= limit && found < nmax; i++)
        if (!comp[i]) primes[found++] = i;
    free(comp);
    if (found < nmax) { free(primes); return NULL; }   /* bound was too tight */

    int64_t* out = malloc(sizeof(int64_t) * sz);
    if (!out) { free(primes); return NULL; }
    for (size_t i = 0; i < sz; i++)
        out[i] = primes[ndt_get_i(in, i, NDT_INT64) - 1];
    free(primes);
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                 out, NDT_INT64);
}

/* ---- PowerMod[arr, e, m] ------------------------------------------------ *
 * Ternary, so no NDBinaryKernel can carry it. Restricted to a NON-NEGATIVE
 * exponent: PowerMod[b, -1, m] is a modular inverse, which exists only when
 * gcd(b, m) is 1, and an array where some elements have an inverse and some do
 * not has no uniform answer — the List path handles that case element by
 * element, and answers unevaluated where it must.
 *
 * The modulus is bounded so that the squaring step cannot overflow int64: with
 * |m| <= 3037000499 (< 2^31.5), every intermediate product is under 2^63. */
#define NDI_POWERMOD_MAX ((int64_t)3037000499LL)

Expr* ndint_powermod(Expr* res) {
    const Expr* a = ndi_int_arg(res, 0, 3);
    if (!a) return NULL;
    int64_t e, m;
    if (!ndi_exact_int(res->data.function.args[1], &e) || e < 0) return NULL;
    if (!ndi_exact_int(res->data.function.args[2], &m)) return NULL;
    if (m == 0) return NULL;
    int64_t am;
    if (ci_abs_i64(m, &am) || am > NDI_POWERMOD_MAX) return NULL;

    size_t sz = ndarray_size(a);
    const void* in = a->data.ndarray.data;
    int64_t* out = malloc(sizeof(int64_t) * (sz ? sz : 1));
    if (!out) return NULL;
    for (size_t i = 0; i < sz; i++) {
        int64_t b = ndt_get_i(in, i, NDT_INT64) % am;
        if (b < 0) b += am;
        int64_t r = 1 % am, k = e;
        while (k > 0) {
            if (k & 1) r = (r * b) % am;
            k >>= 1;
            if (k > 0) b = (b * b) % am;
        }
        /* builtin_powermod takes mpz_abs of the modulus before it starts, so a
         * negative modulus gives the same non-negative residue: PowerMod[3, 2,
         * -5] is 4, not -1. Verified against the scalar path rather than
         * assumed, since Mod's own sign convention is the opposite one. */
        out[i] = r;
    }
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                 out, NDT_INT64);
}

/* ---- IntegerDigits[arr] ------------------------------------------------- *
 * The result is ragged, so the win is entirely on the input side: read each
 * element straight out of the buffer rather than materialising the whole array
 * into Integers first and then reading each one back. Rank 1 only — a
 * higher-rank argument's result nests one level deeper per axis, which the
 * List path already builds correctly and is not the shape anyone measures. */
Expr* ndint_integerdigits(Expr* res) {
    const Expr* a = ndi_int_arg(res, 0, 1);
    if (!a || a->data.ndarray.rank != 1) return NULL;
    size_t sz = ndarray_size(a);
    const void* in = a->data.ndarray.data;

    Expr** rows = malloc(sizeof(Expr*) * (sz ? sz : 1));
    if (!rows) return NULL;
    for (size_t i = 0; i < sz; i++) {
        int64_t v = ndt_get_i(in, i, NDT_INT64), m;
        if (ci_abs_i64(v, &m)) {                    /* INT64_MIN: let GMP have it */
            for (size_t j = 0; j < i; j++) expr_free(rows[j]);
            free(rows);
            return NULL;
        }
        int64_t d[20]; int nd = 0;
        if (m == 0) d[nd++] = 0;
        while (m > 0) { d[nd++] = m % 10; m /= 10; }
        Expr** digits = malloc(sizeof(Expr*) * (size_t)nd);
        if (!digits) {
            for (size_t j = 0; j < i; j++) expr_free(rows[j]);
            free(rows);
            return NULL;
        }
        for (int j = 0; j < nd; j++) digits[j] = expr_new_integer(d[nd - 1 - j]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), digits, (size_t)nd);
        free(digits);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, sz);
    free(rows);
    return out;
}

/* ---- Positive / Negative / NonNegative / NonPositive -------------------- */

Expr* ndint_sign_predicate(Expr* res, NDSignPred which) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* a = res->data.function.args[0];
    if (!is_ndarray(a)) return NULL;
    NDType dt = a->data.ndarray.dtype;
    /* A complex buffer has no order, and the scalar builtins leave Positive[I]
     * unevaluated rather than answering False -- so only the ordered dtypes. */
    if (ndt_is_complex(dt)) return NULL;
    /* Rank > 1 threads to a nested List of True/False, which the List path
     * already builds; the flat vector is the shape that is measured. */
    if (a->data.ndarray.rank != 1) return NULL;

    size_t sz = ndarray_size(a);
    if (sz == 0) return NULL;           /* empty: let the List path build {} */
    const void* in = a->data.ndarray.data;

    /* The answer is now a one-byte-per-element bool BUFFER, not 10^6 shared
     * True/False Expr nodes (performance.md §13 gap C.1). It always presents as a
     * List -- a PACKED bool List, observably identical to the plain List of
     * True/False the predicate returned before this change, whichever surface the
     * input came from. (Unlike Sin, which inherits its input's presentation, a
     * sign predicate has always answered with a List; keeping that is what lets
     * the visible surface agree with the packed one and preserves the contract
     * the tests pin.) */
    void* buf = NULL;
    int64_t dims[1] = { (int64_t)sz };
    Expr* out = ndbuild_open_like(a, 1, dims, NDT_BOOL, &buf);
    if (!out) return NULL;
    out->data.ndarray.present_as = NDA_HEAD_LIST;   /* always a List, never NDArray[...] */
    pack_g_any_created = true;                       /* a packed list now exists */
    uint8_t* ob = (uint8_t*)buf;
    for (size_t i = 0; i < sz; i++) {
        bool val;
        if (dt == NDT_INT64) {
            int64_t v = ndt_get_i(in, i, dt);
            val = (which == NDSP_POSITIVE)    ? (v > 0)
                : (which == NDSP_NEGATIVE)    ? (v < 0)
                : (which == NDSP_NONNEGATIVE) ? (v >= 0)
                                              : (v <= 0);
        } else {
            double re, im;
            ndt_get(in, i, dt, &re, &im);
            /* An indeterminate element is unordered: Positive[Indeterminate]
             * stays symbolic rather than answering False, so abandon the array
             * and let the List path reproduce that element for element. */
            if (re != re) { expr_free(out); return NULL; }
            val = (which == NDSP_POSITIVE)    ? (re > 0.0)
                : (which == NDSP_NEGATIVE)    ? (re < 0.0)
                : (which == NDSP_NONNEGATIVE) ? (re >= 0.0)
                                              : (re <= 0.0);
        }
        ob[i] = val ? 1 : 0;
    }
    return out;
}

/* ---- registration ------------------------------------------------------- */

void ndinteger_init(void) {
    symtab_set_ndarray_unary_kernel("IntegerLength", &NDKU_IntegerLength);
    symtab_set_ndarray_unary_kernel("EulerPhi",      &NDKU_EulerPhi);
    symtab_set_ndarray_unary_kernel("MoebiusMu",     &NDKU_MoebiusMu);
    symtab_set_ndarray_binary_kernel("GCD",          &NDKB_GCD);
    symtab_set_ndarray_binary_kernel("LCM",          &NDKB_LCM);
    symtab_set_ndarray_binary_kernel("DivisorSigma", &NDKB_DivisorSigma);
}
