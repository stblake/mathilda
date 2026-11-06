/* Mathilda -- MathieuCharacteristicA[r, q] and MathieuCharacteristicB[r, q],
 * the characteristic values a_r(q) and b_r(q) of Mathieu's equation
 *
 *     y'' + (a - 2 q cos 2z) y = 0.
 *
 * For an integer r, a_r is the value with an even periodic solution ce_r
 * (period pi for even r, 2 pi for odd r) and b_r (r >= 1) the value with an odd
 * one se_r (DLMF 28.2(vi)); a_{-r} = a_r, b_{-r} = b_r, and there is no b_0
 * (MathieuCharacteristicB::zord). For a non-integer r = nu, a_nu(q) = b_nu(q)
 * is the value whose Floquet exponent is nu (DLMF 28.12): continuous in q from
 * a_nu(0) = nu^2, even in nu and in q.
 *
 * METHOD. Each class is the eigenvalue problem of a symmetric tridiagonal
 * matrix in the Fourier coefficients of the solution (DLMF 28.4):
 *
 *   a_{2m}     cos 2kz, k >= 0:      diag (2k)^2,        off-diagonal q,
 *                                    sqrt(2) q in the first row (symmetrised)
 *   a_{2m+1}   cos (2k+1)z, k >= 0:  diag (2k+1)^2, (0,0) entry 1 + q
 *   b_{2m+1}   sin (2k+1)z, k >= 0:  diag (2k+1)^2, (0,0) entry 1 - q
 *   b_{2m+2}   sin (2k+2)z, k >= 0:  diag (2k+2)^2
 *   a_nu       e^{i (nu+2k) z}, k in Z: diag (nu+2k)^2, nu = |r| mod 2
 *
 * and the wanted value is the m-th smallest eigenvalue of its class (for a_nu,
 * the (2m + [nu > 1])-th, m = floor(|r|/2): the rank of nu^2 among the
 * diagonal at q = 0). For real q the off-diagonals are nonzero and the
 * eigenvalues of an unreduced symmetric tridiagonal matrix are simple, so the
 * rank is preserved from q = 0 for every real q -- the ordering argument is
 * exact, not a tracking heuristic. The matrix is truncated at K where the
 * eigenvalue's truncation error, ~ q c_K c_{K+1} / |c|^2 (the dropped
 * coupling, to first order), is below 2^-(bits + 30): beyond the point where
 * the diagonal exceeds the eigenvalue, |c_{k+1}/c_k| <= q/((2k)^2 - a - q), so
 * the product P of those ratios needs q P^2 < 2^-(bits + 30).
 *
 *   machine    Sturm-sequence bisection in doubles (sf_mathieu_char_r): the
 *              ONE machine helper behind the scalar REPL path, the NDArray /
 *              packed map and Compile[], so they agree bit for bit.
 *   arbitrary  the machine value refined by Newton on the characteristic
 *              polynomial det(T - x) (its logarithmic derivative is a sum
 *              over the LDL^T pivots) in MPFR with the precision doubling,
 *              then CERTIFIED by two Sturm counts that put exactly the right
 *              number of eigenvalues below x - d and x + d; bisection when the
 *              certificate fails.
 *
 * Complex q (whose eigenvalue curves meet at branch points, so the rank
 * argument above no longer holds) and complex r stay unevaluated.
 *
 * Attributes: Listable, NHoldFirst, NumericFunction, Protected. NHoldFirst
 * keeps N from turning an integer order into a Real: N[MathieuCharacteristicA
 * [0, 25], 30] numericalizes q alone.
 *
 * Memory contract: builtins take ownership of res but never free it; they
 * return NULL (unevaluated) or an owned result.
 */
#include "mathieu_char.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "attr.h"
#include "expr.h"
#include "message.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Largest truncation; beyond it (|q| or |r| ~ 10^7) the call declines. */
#define MC_MAX_N 20000

/* Matrix classes. */
enum { MC_A_EVEN, MC_A_ODD, MC_B_ODD, MC_B_EVEN, MC_NONINT };

/* The class and wanted 0-based eigenvalue index for an integer order
 * (|r| = s), or -1 for b_0. */
static int mc_int_class(int odd, int64_t s, int64_t* idx) {
    if (!odd) {
        *idx = s / 2;
        return (s % 2 == 0) ? MC_A_EVEN : MC_A_ODD;
    }
    if (s == 0) return -1;
    if (s % 2 == 1) { *idx = (s - 1) / 2; return MC_B_ODD; }
    *idx = (s - 2) / 2;
    return MC_B_EVEN;
}

/* Half-width of the truncation: the coefficient index past which the product
 * P of the decay ratios q/((2k)^2 - lam - q) has q P^2 < 2^-(bits + 30)
 * (see the header). Returns -1 when it would exceed MC_MAX_N. */
static int64_t mc_size(double s, double aq, long bits) {
    double lam = s * s + 2.0 * aq + 2.0;
    double kc = ceil(0.5 * sqrt(lam + 2.0 * aq)) + 2.0;
    if (!(kc < MC_MAX_N)) return -1;
    int64_t k = (int64_t)kc;
    double lg = 0.0, target = -0.5 * ((double)bits + 30.0 + log2(1.0 + aq));
    while (lg > target) {
        double d = 4.0 * (double)k * (double)k - lam - aq;
        if (d < 1.0) d = 1.0;
        lg += log2(aq / d);
        if (++k > MC_MAX_N) return -1;
    }
    return k + 1;
}

/* ---- machine helper -------------------------------------------------------- */

/* Number of eigenvalues < x of the symmetric tridiagonal (d, e2 = off^2). */
static int64_t mc_count_d(const double* d, const double* e2, int64_t n, double x,
                          double pivmin) {
    int64_t c = 0;
    double t = d[0] - x;
    if (fabs(t) < pivmin) t = -pivmin;
    if (t < 0.0) c++;
    for (int64_t i = 1; i < n; i++) {
        t = d[i] - x - e2[i - 1] / t;
        if (fabs(t) < pivmin) t = -pivmin;
        if (t < 0.0) c++;
    }
    return c;
}

/* Fill the class matrix (diagonal d, squared off-diagonal e2) in doubles.
 * Returns the order, or -1. For MC_NONINT `nu` is |r| mod 2. */
static int64_t mc_build_d(int cls, double nu, double q, int64_t K, double** dp, double** ep) {
    int64_t n = cls == MC_NONINT ? 2 * K + 2 : K + 1;
    double* d = (double*)malloc((size_t)n * sizeof(double));
    double* e2 = (double*)malloc((size_t)n * sizeof(double));
    if (!d || !e2) { free(d); free(e2); return -1; }
    for (int64_t i = 0; i < n; i++) {
        double v;
        switch (cls) {
            case MC_A_EVEN: v = 2.0 * (double)i; break;
            case MC_A_ODD: case MC_B_ODD: v = 2.0 * (double)i + 1.0; break;
            case MC_B_EVEN: v = 2.0 * (double)i + 2.0; break;
            default: v = nu + 2.0 * (double)(i - K - 1); break;
        }
        d[i] = v * v;
        e2[i] = q * q;
    }
    if (cls == MC_A_EVEN && n > 1) e2[0] = 2.0 * q * q;
    if (cls == MC_A_ODD) d[0] += q;
    if (cls == MC_B_ODD) d[0] -= q;
    *dp = d; *ep = e2;
    return n;
}

/* The idx-th smallest eigenvalue by bisection on the Sturm count. */
static double mc_bisect_d(const double* d, const double* e2, int64_t n, int64_t idx) {
    double emax = 0.0, lo = HUGE_VAL, hi = -HUGE_VAL;
    for (int64_t i = 0; i < n; i++) {
        double r = 0.0;
        if (i > 0) r += sqrt(e2[i - 1]);
        if (i + 1 < n) r += sqrt(e2[i]);
        if (d[i] - r < lo) lo = d[i] - r;
        if (d[i] + r > hi) hi = d[i] + r;
        if (i + 1 < n && e2[i] > emax) emax = e2[i];
    }
    double pivmin = DBL_MIN * (emax > 1.0 ? emax : 1.0);
    lo -= 1.0; hi += 1.0;
    for (int it = 0; it < 4000; it++) {
        double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) break;
        if (mc_count_d(d, e2, n, mid, pivmin) <= idx) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* Shared classification of (odd, r) for both precisions: the class, the
 * index, nu (non-integer only) and s = |r|. False for b_0 or an order too
 * large. */
static bool mc_classify_d(int odd, double r, int* cls, int64_t* idx, double* nu) {
    double s = fabs(r);
    if (!(s < 1e7)) return false;
    if (s == floor(s)) {
        int c = mc_int_class(odd, (int64_t)s, idx);
        if (c < 0) return false;
        *cls = c; *nu = 0.0;
        return true;
    }
    double m = floor(0.5 * s);
    *nu = s - 2.0 * m;
    *cls = MC_NONINT;
    *idx = 2 * (int64_t)m + (*nu > 1.0 ? 1 : 0);
    return true;
}

bool sf_mathieu_char_r(int odd, double r, double q, double* out) {
    if (!isfinite(r) || !isfinite(q)) return false;
    int cls; int64_t idx; double nu;
    if (!mc_classify_d(odd, r, &cls, &idx, &nu)) return false;
    if (q == 0.0) { *out = r * r; return true; }
    int64_t K = mc_size(fabs(r), fabs(q), 53);
    if (K < 0) return false;
    if (cls != MC_NONINT && K < idx + 2) K = idx + 2;
    double *d, *e2;
    int64_t n = mc_build_d(cls, nu, q, K, &d, &e2);
    if (n < 0) return false;
    /* For MC_NONINT the rank idx = 2m + [nu > 1] counts the diagonal entries
     * (nu + 2k)^2 < r^2, k in [-(m + 1), m - 1] -- all inside the truncation
     * [-K - 1, K], since K > m + 1 -- so it is the rank in the matrix too. */
    double v = mc_bisect_d(d, e2, n, idx);
    free(d); free(e2);
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool sf_mathieu_char_a_r(double r, double q, double* out) { return sf_mathieu_char_r(0, r, q, out); }
bool sf_mathieu_char_b_r(double r, double q, double* out) { return sf_mathieu_char_r(1, r, q, out); }

/* ---- arbitrary precision ----------------------------------------------------- */

#ifdef USE_MPFR
/* Eigenvalues < x of (d, e2) in MPFR at the working precision of t. */
static int64_t mc_count_m(mpfr_t* d, mpfr_t* e2, int64_t n, const mpfr_t x, mpfr_t t,
                          mpfr_t u, const mpfr_t pivmin) {
    int64_t c = 0;
    mpfr_sub(t, d[0], x, MPFR_RNDN);
    for (int64_t i = 0;; ) {
        if (mpfr_cmpabs(t, pivmin) < 0) mpfr_neg(t, pivmin, MPFR_RNDN);
        if (mpfr_sgn(t) < 0) c++;
        if (++i >= n) break;
        mpfr_div(u, e2[i - 1], t, MPFR_RNDN);
        mpfr_sub(t, d[i], x, MPFR_RNDN);
        mpfr_sub(t, t, u, MPFR_RNDN);
    }
    return c;
}

/* One Newton step on det(T - x): x -= 1 / sum_i p_i'/p_i, the pivots
 * p_0 = d_0 - x, p_i = d_i - x - e2_{i-1}/p_{i-1} and their x-derivatives
 * p_i' = -1 + e2_{i-1} p_{i-1}'/p_{i-1}^2. False on a zero pivot. */
static bool mc_newton_m(mpfr_t* d, mpfr_t* e2, int64_t n, mpfr_t x, mpfr_prec_t wp,
                        double* step) {
    mpfr_t p, dp, s, t, u;
    mpfr_inits2(wp, p, dp, s, t, u, (mpfr_ptr)0);
    bool ok = true;
    mpfr_sub(p, d[0], x, MPFR_RNDN);
    mpfr_set_si(dp, -1, MPFR_RNDN);
    mpfr_set_zero(s, 1);
    for (int64_t i = 0;; ) {
        if (mpfr_zero_p(p)) { ok = false; break; }
        mpfr_div(t, dp, p, MPFR_RNDN);
        mpfr_add(s, s, t, MPFR_RNDN);
        if (++i >= n) break;
        /* dp_i = -1 + e2 * dp / p^2 = -1 + (e2 / p) * (dp / p) */
        mpfr_div(u, e2[i - 1], p, MPFR_RNDN);      /* u = e2/p        */
        mpfr_mul(dp, u, t, MPFR_RNDN);             /* e2 dp / p^2     */
        mpfr_sub_si(dp, dp, 1, MPFR_RNDN);
        mpfr_sub(p, d[i], x, MPFR_RNDN);
        mpfr_sub(p, p, u, MPFR_RNDN);
    }
    if (ok && !mpfr_zero_p(s) && mpfr_number_p(s)) {
        mpfr_ui_div(s, 1, s, MPFR_RNDN);
        mpfr_prec_round(x, wp, MPFR_RNDN);
        mpfr_sub(x, x, s, MPFR_RNDN);
        /* the step relative to max(|x|, 1) */
        double ax = fabs(mpfr_get_d(x, MPFR_RNDN));
        *step = fabs(mpfr_get_d(s, MPFR_RNDN)) / (ax > 1.0 ? ax : 1.0);
    } else {
        ok = false;
    }
    mpfr_clears(p, dp, s, t, u, (mpfr_ptr)0);
    return ok;
}

/* a_r / b_r at `bits` bits. `r` an exact integer order (r_int set) or an MPFR
 * real order; `q` an MPFR real. `seed` the machine value. */
static Expr* mc_char_mpfr(int odd, bool r_int, int64_t rint, const mpfr_t r,
                          const mpfr_t q, long bits) {
    double rd = r_int ? (double)rint : mpfr_get_d(r, MPFR_RNDN);
    double qd = mpfr_get_d(q, MPFR_RNDN);
    if (!isfinite(qd) || !isfinite(rd)) return NULL;
    int cls; int64_t idx; double nud;
    if (!mc_classify_d(odd, rd, &cls, &idx, &nud)) return NULL;
    /* A non-integer MPFR order whose double rounds to an integer (r within
     * 2^-53 of an integer): the classification above must use the MPFR
     * value. Decline rather than mis-classify. */
    if (!r_int && mpfr_integer_p(r)) {
        if (!(fabs(rd) < 1e7)) return NULL;
        r_int = true; rint = (int64_t)rd;
        int c = mc_int_class(odd, rint < 0 ? -rint : rint, &idx);
        if (c < 0) return NULL;
        cls = c;
    } else if (!r_int && cls != MC_NONINT) {
        return NULL;
    }
    double seed;
    if (!sf_mathieu_char_r(odd, rd, qd, &seed)) return NULL;
    int64_t K = mc_size(fabs(rd), fabs(qd), bits);
    if (K < 0) return NULL;
    if (cls != MC_NONINT && K < idx + 2) K = idx + 2;
    int64_t n = cls == MC_NONINT ? 2 * K + 2 : K + 1;
    /* Guard bits: the Sturm recurrence loses ~log2(||T||) bits. */
    double scale = 4.0 * (double)(K + 2) * (double)(K + 2) + 4.0 * fabs(qd);
    mpfr_prec_t wp = (mpfr_prec_t)bits + 32 + (mpfr_prec_t)ceil(log2(scale + 2.0));

    mpfr_t* d = (mpfr_t*)malloc((size_t)n * sizeof(mpfr_t));
    mpfr_t* e2 = (mpfr_t*)malloc((size_t)n * sizeof(mpfr_t));
    if (!d || !e2) { free(d); free(e2); return NULL; }
    mpfr_t nu, v, q2, x, t, u, pivmin, delta, xa;
    mpfr_inits2(wp, nu, v, q2, x, t, u, pivmin, delta, xa, (mpfr_ptr)0);
    if (cls == MC_NONINT) {
        /* nu = |r| - 2 floor(|r|/2) */
        mpfr_abs(nu, r, MPFR_RNDN);
        mpfr_div_2ui(v, nu, 1, MPFR_RNDN);
        mpfr_floor(v, v);
        mpfr_mul_2ui(v, v, 1, MPFR_RNDN);
        mpfr_sub(nu, nu, v, MPFR_RNDN);
    }
    mpfr_sqr(q2, q, MPFR_RNDN);
    for (int64_t i = 0; i < n; i++) {
        mpfr_init2(d[i], wp);
        mpfr_init2(e2[i], wp);
        switch (cls) {
            case MC_A_EVEN: mpfr_set_si(v, 2 * i, MPFR_RNDN); break;
            case MC_A_ODD: case MC_B_ODD: mpfr_set_si(v, 2 * i + 1, MPFR_RNDN); break;
            case MC_B_EVEN: mpfr_set_si(v, 2 * i + 2, MPFR_RNDN); break;
            default:
                mpfr_set_si(v, 2 * (i - K - 1), MPFR_RNDN);
                mpfr_add(v, v, nu, MPFR_RNDN);
                break;
        }
        mpfr_sqr(d[i], v, MPFR_RNDN);
        mpfr_set(e2[i], q2, MPFR_RNDN);
    }
    if (cls == MC_A_EVEN && n > 1) mpfr_mul_2ui(e2[0], q2, 1, MPFR_RNDN);
    if (cls == MC_A_ODD) mpfr_add(d[0], d[0], q, MPFR_RNDN);
    if (cls == MC_B_ODD) mpfr_sub(d[0], d[0], q, MPFR_RNDN);
    /* The index within the matrix. */
    int64_t j = idx;

    /* Newton from the machine seed with the precision doubling: the seed is
     * good to ~50 bits and each step doubles that, so one step per doubled
     * precision -- and more at wp while the step is still above the target
     * (a seed near 0 carries fewer correct bits). */
    mpfr_set_prec(x, 53);
    mpfr_set_d(x, seed, MPFR_RNDN);
    bool ok = true;
    double step = 1.0;
    for (mpfr_prec_t p = 96;; p *= 2) {
        if (p > wp) p = wp;
        ok = mc_newton_m(d, e2, n, x, p, &step);
        if (!ok || p == wp) break;
    }
    for (int it = 0; ok && it < 6 && !(step < ldexp(1.0, -(int)bits - 4)); it++)
        ok = mc_newton_m(d, e2, n, x, wp, &step);

    /* Certificate: exactly j eigenvalues below x - delta, at least j + 1
     * below x + delta, delta = 2^-(bits + 2) max(|x|, 1). */
    mpfr_set_ui(pivmin, 1, MPFR_RNDN);
    mpfr_mul_2si(pivmin, pivmin, -(long)wp * 4, MPFR_RNDN);
    mpfr_prec_round(x, wp, MPFR_RNDN);
    if (ok) {
        mpfr_abs(delta, x, MPFR_RNDN);
        if (mpfr_cmp_ui(delta, 1) < 0) mpfr_set_ui(delta, 1, MPFR_RNDN);
        mpfr_mul_2si(delta, delta, -(long)(bits + 2), MPFR_RNDN);
        mpfr_sub(xa, x, delta, MPFR_RNDN);
        ok = mc_count_m(d, e2, n, xa, t, u, pivmin) == j;
        mpfr_add(xa, x, delta, MPFR_RNDN);
        if (ok) ok = mc_count_m(d, e2, n, xa, t, u, pivmin) >= j + 1;
    }
    if (!ok) {
        /* Bisection from the Gershgorin bracket. */
        mpfr_t lo, hi;
        mpfr_inits2(wp, lo, hi, (mpfr_ptr)0);
        double glo = HUGE_VAL, ghi = -HUGE_VAL;
        for (int64_t i = 0; i < n; i++) {
            double dd = mpfr_get_d(d[i], MPFR_RNDN), rr = 3.0 * fabs(qd);
            if (dd - rr < glo) glo = dd - rr;
            if (dd + rr > ghi) ghi = dd + rr;
        }
        mpfr_set_d(lo, glo - 1.0, MPFR_RNDD);
        mpfr_set_d(hi, ghi + 1.0, MPFR_RNDU);
        for (long it = 0; it < (long)wp + 80; it++) {
            mpfr_add(xa, lo, hi, MPFR_RNDN);
            mpfr_div_2ui(xa, xa, 1, MPFR_RNDN);
            if (mpfr_equal_p(xa, lo) || mpfr_equal_p(xa, hi)) break;
            if (mc_count_m(d, e2, n, xa, t, u, pivmin) <= j) mpfr_set(lo, xa, MPFR_RNDN);
            else mpfr_set(hi, xa, MPFR_RNDN);
        }
        mpfr_add(x, lo, hi, MPFR_RNDN);
        mpfr_div_2ui(x, x, 1, MPFR_RNDN);
        mpfr_clears(lo, hi, (mpfr_ptr)0);
    }
    Expr* out = NULL;
    if (mpfr_number_p(x)) {
        mpfr_t o;
        mpfr_init2(o, (mpfr_prec_t)bits);
        mpfr_set(o, x, MPFR_RNDN);
        out = expr_new_mpfr_move(o);
    }
    for (int64_t i = 0; i < n; i++) { mpfr_clear(d[i]); mpfr_clear(e2[i]); }
    free(d); free(e2);
    mpfr_clears(nu, v, q2, x, t, u, pivmin, delta, xa, (mpfr_ptr)0);
    return out;
}
#endif

/* ---- builtins ------------------------------------------------------------------ */

/* A real number leaf (exact Integer / BigInt / Rational, or inexact). */
static bool mc_is_real_number(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION:
            return e->data.function.head->type == EXPR_SYMBOL &&
                   e->data.function.head->data.symbol.name == SYM_Rational;
        default:
            return false;
    }
}

static bool mc_is_zero(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer == 0) ||
           (e->type == EXPR_REAL && e->data.real == 0.0)
#ifdef USE_MPFR
           || (e->type == EXPR_MPFR && mpfr_zero_p(e->data.mpfr))
#endif
        ;
}

static bool mc_is_int_zero(const Expr* r) {
    if (r->type == EXPR_INTEGER) return r->data.integer == 0;
    if (r->type == EXPR_REAL) return r->data.real == 0.0;
#ifdef USE_MPFR
    if (r->type == EXPR_MPFR) return mpfr_zero_p(r->data.mpfr);
#endif
    return false;
}

static Expr* mc_builtin(Expr* res, int odd) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* r = res->data.function.args[0];
    const Expr* q = res->data.function.args[1];
    if (!mc_is_real_number(r) || !mc_is_real_number(q)) return NULL;
    if (odd && mc_is_int_zero(r)) {
        mth_message("MathieuCharacteristicB", "zord",
                    "There is no zero-order MathieuCharacteristicB.");
        return NULL;
    }
    bool inexact = arg_is_inexact(r) || arg_is_inexact(q);
    if (!inexact) {
        /* Exact: only q = 0 has a closed form, a_r(0) = b_r(0) = r^2. */
        if (mc_is_zero(q)) {
            Expr* args[2] = { expr_copy((Expr*)r), expr_new_integer(2) };
            return expr_new_function(expr_new_symbol(SYM_Power), args, 2);
        }
        return NULL;
    }
    long bits = 0;
    for (int i = 0; i < 2; i++) {
        long b = numeric_min_inexact_bits(res->data.function.args[i]);
        if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    if (r->type == EXPR_REAL || q->type == EXPR_REAL) bits = 53;
    if (bits <= 53) {
        /* The common shapes (an Integer or Real order, a Real q) read
         * directly; anything else through numericalize. */
        if ((r->type == EXPR_INTEGER || r->type == EXPR_REAL) && q->type == EXPR_REAL) {
            double rd = r->type == EXPR_INTEGER ? (double)r->data.integer : r->data.real, v;
            return sf_mathieu_char_r(odd, rd, q->data.real, &v) ? expr_new_real(v) : NULL;
        }
        NumericSpec ms = numeric_machine_spec();
        Expr* rn = numericalize(r, ms);
        Expr* qn = numericalize(q, ms);
        Expr* out = NULL;
        double v;
        if (rn && qn && rn->type == EXPR_REAL && qn->type == EXPR_REAL &&
            sf_mathieu_char_r(odd, rn->data.real, qn->data.real, &v))
            out = expr_new_real(v);
        if (rn) expr_free(rn);
        if (qn) expr_free(qn);
        return out;
    }
#ifdef USE_MPFR
    NumericSpec ps;
    ps.mode = NUMERIC_MODE_MPFR;
    ps.bits = bits;
    Expr* qn = numericalize(q, ps);
    Expr* rn = NULL;
    Expr* out = NULL;
    bool r_int = r->type == EXPR_INTEGER;
    int64_t rint = r_int ? r->data.integer : 0;
    if (!r_int) rn = numericalize(r, ps);
    if (qn && qn->type == EXPR_MPFR && (r_int || (rn && rn->type == EXPR_MPFR))) {
        mpfr_t rz;
        mpfr_init2(rz, 64);
        mpfr_set_si(rz, (long)rint, MPFR_RNDN);
        out = mc_char_mpfr(odd, r_int, rint, r_int ? rz : rn->data.mpfr, qn->data.mpfr, bits);
        mpfr_clear(rz);
    }
    if (qn) expr_free(qn);
    if (rn) expr_free(rn);
    return out;
#else
    return NULL;
#endif
}

Expr* builtin_mathieucharacteristica(Expr* res) { return mc_builtin(res, 0); }
Expr* builtin_mathieucharacteristicb(Expr* res) { return mc_builtin(res, 1); }

void mathieu_char_init(void) {
    const unsigned attrs = ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_NUMERICFUNCTION |
                           ATTR_PROTECTED;
    symtab_add_builtin("MathieuCharacteristicA", builtin_mathieucharacteristica);
    symtab_get_def("MathieuCharacteristicA")->attributes |= attrs;
    symtab_set_docstring("MathieuCharacteristicA",
        "MathieuCharacteristicA[r, q]\n"
        "\tThe characteristic value a_r(q) of Mathieu's equation\n"
        "\ty'' + (a - 2 q cos 2z) y = 0: for integer r the value with an even\n"
        "\tperiodic solution, for non-integer r the value with Floquet exponent r.\n"
        "\tReal r and q; machine or arbitrary precision.");
    symtab_add_builtin("MathieuCharacteristicB", builtin_mathieucharacteristicb);
    symtab_get_def("MathieuCharacteristicB")->attributes |= attrs;
    symtab_set_docstring("MathieuCharacteristicB",
        "MathieuCharacteristicB[r, q]\n"
        "\tThe characteristic value b_r(q) of Mathieu's equation: for integer\n"
        "\tr >= 1 the value with an odd periodic solution; for non-integer r it\n"
        "\tequals MathieuCharacteristicA[r, q]. Real r and q.");
}
