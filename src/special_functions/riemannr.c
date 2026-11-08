/* Mathilda -- RiemannR[x], Riemann's prime-counting function.
 *
 *   R(x) = Sum_{k>=1} MoebiusMu[k]/k LogIntegral[x^(1/k)]
 *        = 1 + Sum_{k>=1} (Log x)^k / (k k! Zeta[k + 1])        (Gram series)
 *
 * Every evaluation uses the Gram series in L = Log x: it is entire in L, so it
 * converges for every x != 0 (complex x through the principal Log), where the
 * Moebius sum converges only conditionally and very slowly.  Its derivative
 * comes out of the same pass,
 *
 *   R'(x) = 1/(x L) Sum_{k>=1} L^k / (k! Zeta[k + 1]),
 *
 * and is what Derivative[1][RiemannR] evaluates to at an inexact point, so
 * FindRoot can take Newton steps on R instead of finite differences.
 *
 * The difficulty is cancellation.  For L > 0 every term is positive and the
 * sum is as well-conditioned as its largest term.  For L < 0 (0 < x < 1) the
 * terms alternate, grow to about e^|L| / |L|^(3/2) near k = |L|, and cancel
 * down to R(x), which is O(1) or smaller (R(e^-u) ~ 21.9/u^3): at the largest
 * zero of R, u = 34142, the terms reach 10^14820 and the sum is 10^-13.  The
 * arbitrary-precision path therefore works with
 *
 *   W = p + g + log2(max term) + log2(terms) + 16     bits,
 *
 * and is ADAPTIVE in g: an answer smaller than its error budget (the sum
 * cancelled further than the guard bits allowed for) is recomputed with
 * enough extra bits to give p correct bits relative to |R|.
 *
 * 1/Zeta[k+1] costs nothing like a full-precision zeta per term.  Term k only
 * needs s_k = t_k / Zeta[k+1] to the ABSOLUTE accuracy of the sum, i.e. to
 * q_k = log2|t_k| - (target exponent) bits, and
 *
 *   1/Zeta[s] = Sum_{n>=1} MoebiusMu[n] / n^s
 *
 * reaches that after n^(k+1) > 2^q_k: for the large k that dominate the cost
 * one or two corrections (n = 2 is an exact shift).  Only small k, where that
 * Dirichlet series would be long, call mpfr_zeta_ui -- at q_k bits, not W.
 *
 * Evaluation by argument class:
 *   RiemannR[1] = 1 exactly; any other exact argument stays unevaluated (as in
 *     the reference: RiemannR[2], RiemannR[0], RiemannR[x]).
 *   machine real x >= 1/e   ->  rr_machine_real: the Gram series in double-
 *     double (t_k and the partial sum), shared with the NDArray / packed map
 *     and Compile[]; 1/e is where the alternating series starts to cancel.
 *   machine complex x, |Log x| <= 3/2  ->  rr_machine_complex, the same.
 *   any other machine argument  ->  the MPFR path at 53 bits, rounded.
 *   MPFR real / complex     ->  the MPFR path at the argument's precision.
 *   0. -> 0. (the limit x -> 0+).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "riemannr.h"

#include <math.h>
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdbool.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"   /* is_complex, make_complex */
#include "attr.h"
#include "common.h"       /* builtin_arg_error */
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "numeric.h"
#include "sf_scalar.h"
#include "sym_names.h"
#include "symtab.h"
#include "sym_intern.h"     /* intern_symbol */
#include "deriv.h"            /* deriv_register_numeric */
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

/* ------------------------------------------------------------------ */
/* Machine kernels                                                    */
/* ------------------------------------------------------------------ */

/* Zeta[n] for n = 2..64 rounded to double; Zeta[n] = 1 + 2^-n + ... rounds
 * to 1 + 2^-n (to within half an ulp) beyond the table. */
static const double rr_zeta_tab[65] = {
    0.0, 0.0, 1.64493406684822643647, 1.20205690315959428540,
    1.08232323371113819152, 1.03692775514336992633, 1.01734306198444913971, 1.00834927738192282684,
    1.00407735619794433938, 1.00200839282608221442, 1.00099457512781808534, 1.00049418860411946456,
    1.00024608655330804830, 1.00012271334757848915, 1.00006124813505870483, 1.00003058823630702049,
    1.00001528225940865187, 1.00000763719763789976, 1.00000381729326499984, 1.00000190821271655394,
    1.00000095396203387280, 1.00000047693298678781, 1.00000023845050272773, 1.00000011921992596531,
    1.00000005960818905126, 1.00000002980350351465, 1.00000001490155482837, 1.00000000745071178984,
    1.00000000372533402479, 1.00000000186265972351, 1.00000000093132743242, 1.00000000046566290650,
    1.00000000023283118337, 1.00000000011641550173, 1.00000000005820772088, 1.00000000002910385044,
    1.00000000001455192189, 1.00000000000727595984, 1.00000000000363797955, 1.00000000000181898965,
    1.00000000000090949478, 1.00000000000045474738, 1.00000000000022737368, 1.00000000000011368684,
    1.00000000000005684342, 1.00000000000002842171, 1.00000000000001421085, 1.00000000000000710543,
    1.00000000000000355271, 1.00000000000000177636, 1.00000000000000088818, 1.00000000000000044409,
    1.00000000000000022204, 1.00000000000000011102, 1.00000000000000005551, 1.00000000000000002776,
    1.00000000000000001388, 1.00000000000000000694, 1.00000000000000000347, 1.00000000000000000173,
    1.00000000000000000087, 1.00000000000000000043, 1.00000000000000000022, 1.00000000000000000011,
    1.00000000000000000005,
};

static double rr_zeta_d(long n) {
    if (n <= 64) return rr_zeta_tab[n];
    return 1.0 + ldexp(1.0, (int)-n);
}

/* Double-double helpers (Dekker; no fma, so the result does not depend on
 * whether the compiler contracts -- the makefile builds -ffp-contract=off). */
typedef struct { double hi, lo; } rr_dd;

static void rr_two_sum(double a, double b, double* s, double* e) {
    double t = a + b, bb = t - a;
    *s = t;
    *e = (a - (t - bb)) + (b - bb);
}
static void rr_split(double a, double* h, double* l) {
    double c = 134217729.0 * a;          /* 2^27 + 1 */
    double big = c - a;
    *h = c - big;
    *l = a - *h;
}
static void rr_two_prod(double a, double b, double* p, double* e) {
    double ah, al, bh, bl;
    *p = a * b;
    rr_split(a, &ah, &al);
    rr_split(b, &bh, &bl);
    *e = ((ah * bh - *p) + ah * bl + al * bh) + al * bl;
}
static rr_dd rr_dd_mul_d(rr_dd a, double b) {
    double p, e;
    rr_two_prod(a.hi, b, &p, &e);
    e += a.lo * b;
    rr_dd r; rr_two_sum(p, e, &r.hi, &r.lo);
    return r;
}
static rr_dd rr_dd_div_d(rr_dd a, double b) {
    double q1 = a.hi / b;
    rr_dd p = { 0.0, 0.0 };
    rr_two_prod(q1, b, &p.hi, &p.lo);
    double s, e;
    rr_two_sum(a.hi, -p.hi, &s, &e);
    e -= p.lo;
    e += a.lo;
    double q2 = (s + e) / b;
    rr_dd r; rr_two_sum(q1, q2, &r.hi, &r.lo);
    return r;
}
static rr_dd rr_dd_add(rr_dd a, rr_dd b) {
    double s, e;
    rr_two_sum(a.hi, b.hi, &s, &e);
    e += a.lo + b.lo;
    rr_dd r; rr_two_sum(s, e, &r.hi, &r.lo);
    return r;
}

static rr_dd rr_dd_mul_dd(rr_dd a, rr_dd b) {
    double p, e;
    rr_two_prod(a.hi, b.hi, &p, &e);
    e += a.hi * b.lo + a.lo * b.hi;
    rr_dd r; rr_two_sum(p, e, &r.hi, &r.lo);
    return r;
}

/* log x in double-double, absolute error ~3e-17 for every positive normal or
 * subnormal x: x = m 2^e with m in [sqrt(1/2), sqrt(2)), log x = e ln 2 +
 * log1p(m - 1) (m - 1 is exact), ln 2 split hi + lo.  A plain log(x) is off
 * by half an ulp of |log x| -- 6e-14 at x = 10^300 -- and R(x) inherits that
 * error relative (dR/R ~ dL for large L), 300 ulps. */
static rr_dd rr_log_dd(double x) {
    static const double ln2_hi = 0x1.62e42fefa39efp-1;
    static const double ln2_lo = 0x1.abc9e3b39803fp-56;
    int e;
    double m = frexp(x, &e);               /* x = m 2^e, m in [1/2, 1) */
    if (m < 0.70710678118654752440) { m *= 2.0; e -= 1; }
    double lm = log1p(m - 1.0);
    rr_dd L;
    rr_two_prod((double)e, ln2_hi, &L.hi, &L.lo);
    L.lo += (double)e * ln2_lo;
    rr_dd t = { lm, 0.0 };
    return rr_dd_add(L, t);
}

/* The real Gram series for L = log x >= -1 in double-double.  The recurrence
 * t_k = t_{k-1} L (k-1)/k^2 is carried to ~2^-100 relative per step, so the
 * ~e L steps before the peak at large x cost no accuracy; the only rounding
 * that reaches the answer is Zeta[k+1]'s own (half an ulp, on terms of one
 * sign when L > 0) and the final one.  Writes R and, when dout != NULL,
 * Sum L^k/(k! Zeta[k+1]) (the derivative's numerator). */
static bool rr_gram_dd(rr_dd L, double* rout, double* dout) {
    if (!(L.hi >= -1.0) || !(L.hi <= 709.79)) return false;
    rr_dd t = L;                        /* t_1 = L */
    rr_dd S = { 0.0, 0.0 }, D = { 0.0, 0.0 };
    double aL = fabs(L.hi);
    for (long k = 1; k < 4000; k++) {
        if (k > 1) {
            t = rr_dd_mul_dd(t, L);
            t = rr_dd_mul_d(t, (double)(k - 1));
            t = rr_dd_div_d(t, (double)k * (double)k);
        }
        rr_dd s = rr_dd_div_d(t, rr_zeta_d(k + 1));
        S = rr_dd_add(S, s);
        if (dout) D = rr_dd_add(D, rr_dd_mul_d(s, (double)k));
        if ((double)k > 2.0 * aL + 2.0) {
            double scale = fabs(S.hi) > 1.0 ? fabs(S.hi) : 1.0;
            double dscale = fabs(D.hi) > 1.0 ? fabs(D.hi) : 1.0;
            if (fabs(t.hi) < 0x1p-60 * scale &&
                (!dout || fabs(t.hi) * (double)k < 0x1p-60 * dscale))
                break;
        }
    }
    rr_dd one = { 1.0, 0.0 };
    rr_dd R = rr_dd_add(one, S);
    double r = R.hi + R.lo;
    if (!isfinite(r)) return false;
    *rout = r;
    if (dout) { *dout = D.hi + D.lo; if (!isfinite(*dout)) return false; }
    return true;
}

bool sf_riemannr_r(double x, double* out) {
    if (!(x > 0.0) || !isfinite(x)) return false;
    if (x == 1.0) { *out = 1.0; return true; }
    return rr_gram_dd(rr_log_dd(x), out, NULL);
}

/* The complex Gram series in plain double for |L| <= 3/2 (cancellation then
 * costs at most a couple of bits), L = Log z principal.  Declines a positive
 * real z, which the real kernel or the MPFR path owns, so a real argument
 * never comes back as a Complex. */
bool sf_riemannr_c(double zr, double zi, double* ore, double* oim) {
    if (!isfinite(zr) || !isfinite(zi)) return false;
    if (zi == 0.0 && zr >= 0.0) return false;
    double Lr = log(hypot(zr, zi)), Li = atan2(zi, zr);
    if (!(hypot(Lr, Li) <= 1.5)) return false;
    double tr = Lr, ti = Li, Sr = 0.0, Si = 0.0;
    for (long k = 1; k < 200; k++) {
        if (k > 1) {
            double f = (double)(k - 1) / ((double)k * (double)k);
            double nr = (tr * Lr - ti * Li) * f, ni = (tr * Li + ti * Lr) * f;
            tr = nr; ti = ni;
        }
        double z = rr_zeta_d(k + 1);
        Sr += tr / z; Si += ti / z;
        if (k > 6 && hypot(tr, ti) < 0x1p-60) break;
    }
    *ore = 1.0 + Sr; *oim = Si;
    return isfinite(*ore) && isfinite(*oim);
}

/* NDArray / Compile complex ABI. */
bool riemannr_machine_complex(double are, double aim, double* ore, double* oim) {
    if (aim == 0.0 && are > 0.0) {
        if (!sf_riemannr_r(are, ore)) return false;
        *oim = 0.0;
        return true;
    }
    return sf_riemannr_c(are, aim, ore, oim);
}

/* ------------------------------------------------------------------ */
/* Arbitrary precision                                                */
/* ------------------------------------------------------------------ */
#ifdef USE_MPFR

/* MoebiusMu[n] for n <= RR_MU_MAX, sieved once (riemannr_init). */
#define RR_MU_MAX 1024
static signed char rr_mu[RR_MU_MAX + 1];

static void rr_mu_sieve(void) {
    for (int n = 0; n <= RR_MU_MAX; n++) rr_mu[n] = 1;
    rr_mu[0] = 0;
    for (int p = 2; p <= RR_MU_MAX; p++) {
        bool prime = true;
        for (int d = 2; d * d <= p; d++) if (p % d == 0) { prime = false; break; }
        if (!prime) continue;
        for (int m = p; m <= RR_MU_MAX; m += p) rr_mu[m] = (signed char)-rr_mu[m];
        for (long m = (long)p * p; m <= RR_MU_MAX; m += (long)p * p) rr_mu[m] = 0;
    }
}

/* log2|a + b i| as a double estimate (exponent-level accuracy suffices). */
static double rr_log2abs(const mpfr_t a, const mpfr_t b, bool cplx) {
    double ea = mpfr_zero_p(a) ? -1e300 : (double)mpfr_get_exp(a);
    double eb = (!cplx || mpfr_zero_p(b)) ? -1e300 : (double)mpfr_get_exp(b);
    return ea > eb ? ea : eb;
}

/* log2 |t_k| = log2 (|L|^k / (k k!)) by lgamma: the planning estimate. */
static double rr_log2_term(double aL, long k) {
    return ((double)k * log(aL) - lgamma((double)k + 1.0) - log((double)k)) / M_LN2;
}

/* Complex-or-real in-place multiply u <- u y (y = yr + i yi), tmp scratch. */
static void rr_cmul(mpfr_t ur, mpfr_t ui, const mpfr_t yr, const mpfr_t yi, bool cplx,
                    mpfr_t tmp1, mpfr_t tmp2) {
    if (!cplx) { mpfr_mul(ur, ur, yr, MPFR_RNDN); return; }
    mpfr_mul(tmp1, ur, yr, MPFR_RNDN);
    mpfr_mul(tmp2, ui, yi, MPFR_RNDN);
    mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);          /* re */
    mpfr_mul(tmp2, ur, yi, MPFR_RNDN);
    mpfr_mul(ui, ui, yr, MPFR_RNDN);
    mpfr_add(ui, ui, tmp2, MPFR_RNDN);              /* im */
    mpfr_set(ur, tmp1, MPFR_RNDN);
}

/* u <- u (k-1) / (n k^2), both parts.  One mpfr_div_ui by n k^2 when that
 * fits an unsigned long (the division is the inner loop's dominant cost;
 * unsigned long is 32 bits on LLP64, hence the test, not an assumption). */
static void rr_step_scalar(mpfr_t ur, mpfr_t ui, bool cplx, long k, long n) {
    unsigned long uk = (unsigned long)k, un = (unsigned long)n;
    bool one_div = uk <= ULONG_MAX / uk && uk * uk <= ULONG_MAX / un;
    for (int part = 0; part < (cplx ? 2 : 1); part++) {
        mpfr_ptr u = part ? ui : ur;
        mpfr_mul_ui(u, u, uk - 1, MPFR_RNDN);
        if (one_div) {
            mpfr_div_ui(u, u, uk * uk * un, MPFR_RNDN);
        } else {
            mpfr_div_ui(u, u, uk, MPFR_RNDN);
            mpfr_div_ui(u, u, uk, MPFR_RNDN);
            if (un > 1) mpfr_div_ui(u, u, un, MPFR_RNDN);
        }
    }
}

/* Past its peak a term needs fewer bits every step (it is only ever added to
 * the sum); drop the excess so the long tail runs at its own size. */
static void rr_trim(mpfr_t ur, mpfr_t ui, bool cplx, double need) {
    mpfr_prec_t want = (mpfr_prec_t)(need > 64.0 ? need : 64.0);
    if (want + 256 < mpfr_get_prec(ur)) {
        mpfr_prec_round(ur, want, MPFR_RNDN);
        if (cplx) mpfr_prec_round(ui, want, MPFR_RNDN);
    }
}

/* (Sr, Si) += sgn (ur, ui); (Dr, Di) += sgn k (ur, ui) when want_d. */
static void rr_accum(mpfr_t Sr, mpfr_t Si, mpfr_t Dr, mpfr_t Di, const mpfr_t ur,
                     const mpfr_t ui, bool cplx, int sgn, long k, bool want_d, mpfr_t tmp) {
    for (int part = 0; part < (cplx ? 2 : 1); part++) {
        mpfr_srcptr u = part ? ui : ur;
        if (sgn > 0) mpfr_add(part ? Si : Sr, part ? Si : Sr, u, MPFR_RNDN);
        else         mpfr_sub(part ? Si : Sr, part ? Si : Sr, u, MPFR_RNDN);
        if (want_d) {
            mpfr_set_prec(tmp, mpfr_get_prec(u));
            mpfr_mul_ui(tmp, u, (unsigned long)k, MPFR_RNDN);
            if (sgn > 0) mpfr_add(part ? Di : Dr, part ? Di : Dr, tmp, MPFR_RNDN);
            else         mpfr_sub(part ? Di : Dr, part ? Di : Dr, tmp, MPFR_RNDN);
        }
    }
}

/* One n >= 2 of the 1/Zeta Dirichlet series, run as its own recurrence:
 *   u_k = t_k / n^(k+1) = (1/n) (L/n)^k / (k k!),  u_k = u_{k-1} L (k-1)/(n k^2),
 * at the precision its largest term needs (far below W for every n >= 2).
 * The multiplier is the SAME rounded L the n = 1 recurrence uses, divided by
 * n in a separate correctly rounded step: the stream starts at t_k0, so a
 * differently rounded L/n would make it the tail of a series in another
 * variable, whose truncated head no longer cancels -- an error of the size of
 * the largest term, not of R. */
typedef struct {
    long n;
    int mu;
    bool active;
    mpfr_t ur, ui;            /* the term u_k */
} RRStream;

/* The Gram series at L = Lr + i Li (Li ignored unless cplx), every term to the
 * absolute accuracy 2^sigma.  Outputs R - 1 = Sum t_k/Zeta[k+1] in (Sr, Si) and,
 * when want_d, Sum k t_k/Zeta[k+1] in (Dr, Di); all four are initialised here
 * at the work precision.  Returns false only on an absurd size request.
 *
 * Split at k0:
 *   k < k0   s_k = t_k / Zeta[k+1], mpfr_zeta_ui at the q_k bits term k needs;
 *   k >= k0  s_k = Sum_{n <= N} MoebiusMu[n] t_k / n^(k+1), where k0 is the
 *            first k with 2^(q_k/(k+1)) <= N, so n > N is below the budget
 *            for every k >= k0.  The n = 1 term is t_k itself; each n >= 2
 *            is an RRStream at its own (much lower) precision.
 * No per-k zeta (20 ms each at 10^4 bits) and no per-(k, n) power. */
static bool rr_gram_mpfr(const mpfr_t Lr, const mpfr_t Li, bool cplx, long p_bits,
                         long guard, bool want_d,
                         mpfr_t Sr, mpfr_t Si, mpfr_t Dr, mpfr_t Di) {
    double aL = cplx ? hypot(mpfr_get_d(Lr, MPFR_RNDN), mpfr_get_d(Li, MPFR_RNDN))
                     : fabs(mpfr_get_d(Lr, MPFR_RNDN));
    bool positive = !cplx && mpfr_sgn(Lr) > 0;
    if (!isfinite(aL) || aL > 1e9) return false;
    if (aL < 1e-300) aL = 1e-300;

    /* log2 of the largest |t_k|. */
    double lmax = -1e300;
    long kpeak = (long)aL;
    for (long k = (kpeak > 3 ? kpeak - 2 : 1); k <= kpeak + 3; k++) {
        double v = rr_log2_term(aL, k);
        if (v > lmax) lmax = v;
    }
    if (rr_log2_term(aL, 1) > lmax) lmax = rr_log2_term(aL, 1);
    /* Absolute target: relative to the largest term when no cancellation is
     * possible (L > 0), relative to 1 otherwise (the adaptive caller widens
     * the guard when |R| comes out small). */
    double sigma = positive ? lmax - (double)(p_bits + guard) : -(double)(p_bits + guard);
    double kmax_est = 2.0 * M_E * aL + 2.0 * (double)(p_bits + guard) + 64.0;
    double lK = log2(kmax_est + 2.0);
    double tau = sigma - lK - 8.0;                 /* per-term absolute accuracy */
    double wd = (lmax > sigma ? lmax : sigma) - tau + 16.0;
    if (wd < (double)p_bits + 32.0) wd = (double)p_bits + 32.0;
    if (wd > 1e8) return false;
    mpfr_prec_t W = (mpfr_prec_t)wd;

    /* k0: the first k with q_k / (k+1) <= log2 N  (q_k = log2|t_k| - tau). */
    const double log2N = 10.0;                     /* N = 1024 = RR_MU_MAX */
    long k0 = 1;
    while (k0 < 100000000L) {
        double q = rr_log2_term(aL, k0) + (want_d ? log2((double)k0) : 0.0) - tau + 4.0;
        if (q / (double)(k0 + 1) <= log2N) break;
        k0++;
    }

    mpfr_init2(Sr, W); mpfr_init2(Si, W); mpfr_init2(Dr, W); mpfr_init2(Di, W);
    mpfr_set_ui(Sr, 0, MPFR_RNDN); mpfr_set_ui(Si, 0, MPFR_RNDN);
    mpfr_set_ui(Dr, 0, MPFR_RNDN); mpfr_set_ui(Di, 0, MPFR_RNDN);

    mpfr_t tr, ti, t1, t2, tmp, zq;
    mpfr_inits2(W, tr, ti, t1, t2, tmp, (mpfr_ptr)0);
    mpfr_init2(zq, 64);
    mpfr_set(tr, Lr, MPFR_RNDN);
    if (cplx) mpfr_set(ti, Li, MPFR_RNDN); else mpfr_set_ui(ti, 0, MPFR_RNDN);

    RRStream* st = NULL;
    size_t nst = 0;

    for (long k = 1;; k++) {
        if (k > 1) {                                 /* t <- t L (k-1)/k^2 */
            rr_cmul(tr, ti, Lr, Li, cplx, t1, t2);
            rr_step_scalar(tr, ti, cplx, k, 1);
        }
        double lt = rr_log2abs(tr, ti, cplx);
        if ((double)k > aL + 2.0 && (k & 63) == 0) rr_trim(tr, ti, cplx, lt - tau + lK + 32.0);
        double ltd = lt + (want_d ? log2((double)k) : 0.0);
        bool streams_live = false;
        for (size_t j = 0; j < nst; j++) if (st[j].active) { streams_live = true; break; }
        if ((double)k > 2.0 * aL + 2.0 && ltd < tau && !streams_live) break;
        if (k > 100000000L) break;

        if (k < k0) {
            double q = ltd - tau + 4.0;              /* bits s_k needs */
            if (q <= 0) continue;
            mpfr_prec_t qb = (mpfr_prec_t)q + 16;
            mpfr_set_prec(zq, qb);
            mpfr_zeta_ui(zq, (unsigned long)(k + 1), MPFR_RNDN);
            mpfr_set_prec(t1, qb); mpfr_set_prec(t2, qb);
            mpfr_div(t1, tr, zq, MPFR_RNDN);
            if (cplx) mpfr_div(t2, ti, zq, MPFR_RNDN); else mpfr_set_ui(t2, 0, MPFR_RNDN);
            rr_accum(Sr, Si, Dr, Di, t1, t2, cplx, 1, k, want_d, tmp);
            mpfr_set_prec(t1, W); mpfr_set_prec(t2, W);
            continue;
        }

        if (k == k0) {
            /* Open a stream for each squarefree n <= N whose terms reach the
             * budget at some k >= k0.  Its precision is set by its largest
             * term, at k ~ |L|/n or at k0, whichever is later. */
            st = malloc(sizeof(RRStream) * (RR_MU_MAX + 1));
            if (!st) break;
            for (long n = 2; n <= RR_MU_MAX; n++) {
                if (rr_mu[n] == 0) continue;
                double ln = log2((double)n);
                long kp = (long)(aL / (double)n);
                if (kp < k0) kp = k0;
                double best = -1e300;
                for (long kk = kp; kk <= kp + 2; kk++) {
                    double v = rr_log2_term(aL, kk) - (double)(kk + 1) * ln +
                               (want_d ? log2((double)kk) : 0.0);
                    if (v > best) best = v;
                }
                double Qn = best - tau + lK + 8.0;
                if (Qn <= 0) {
                    if (kp == k0) break;     /* larger n are smaller still */
                    continue;
                }
                RRStream* s = &st[nst++];
                s->n = n; s->mu = rr_mu[n]; s->active = true;
                mpfr_prec_t qp = (mpfr_prec_t)Qn + 16;
                mpfr_inits2(qp, s->ur, s->ui, (mpfr_ptr)0);
                /* u_k0 = t_k0 / n^(k0+1) */
                mpfr_set_prec(zq, qp);
                mpfr_ui_pow_ui(zq, (unsigned long)n, (unsigned long)(k0 + 1), MPFR_RNDN);
                mpfr_div(s->ur, tr, zq, MPFR_RNDN);
                if (cplx) mpfr_div(s->ui, ti, zq, MPFR_RNDN); else mpfr_set_ui(s->ui, 0, MPFR_RNDN);
            }
        }

        /* n = 1: t_k itself (skipped once below the budget). */
        if (ltd >= tau) rr_accum(Sr, Si, Dr, Di, tr, ti, cplx, 1, k, want_d, tmp);
        for (size_t j = 0; j < nst; j++) {
            RRStream* s = &st[j];
            if (!s->active) continue;
            if (k > k0) {
                mpfr_set_prec(t1, mpfr_get_prec(s->ur));
                mpfr_set_prec(t2, mpfr_get_prec(s->ur));
                rr_cmul(s->ur, s->ui, Lr, Li, cplx, t1, t2);
                rr_step_scalar(s->ur, s->ui, cplx, k, s->n);
            }
            double lu = rr_log2abs(s->ur, s->ui, cplx) + (want_d ? log2((double)k) : 0.0);
            if ((double)k > aL / (double)s->n + 2.0 && (k & 63) == 0)
                rr_trim(s->ur, s->ui, cplx, lu - tau + lK + 32.0);
            if (lu < tau) {
                if ((double)k > 2.0 * aL / (double)s->n + 2.0) {
                    s->active = false;
                    mpfr_clears(s->ur, s->ui, (mpfr_ptr)0);
                }
                continue;
            }
            rr_accum(Sr, Si, Dr, Di, s->ur, s->ui, cplx, s->mu, k, want_d, tmp);
        }
        mpfr_set_prec(t1, W); mpfr_set_prec(t2, W);
    }
    for (size_t j = 0; j < nst; j++)
        if (st[j].active) mpfr_clears(st[j].ur, st[j].ui, (mpfr_ptr)0);
    free(st);
    mpfr_clears(tr, ti, t1, t2, tmp, zq, (mpfr_ptr)0);
    return true;
}

/* One-entry memo of the last expensive evaluation.  FindRoot asks for R(z)
 * and then R'(z) at the SAME point; at z ~ e^-34142 one Gram pass costs
 * seconds, and the derivative's sum rides along in that pass for a few
 * percent more.  Keyed on the exact MPFR point (value and width) and the bit
 * target. */
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_rr_memo_lock = PTHREAD_MUTEX_INITIALIZER;
#define RR_LOCK()   pthread_mutex_lock(&g_rr_memo_lock)
#define RR_UNLOCK() pthread_mutex_unlock(&g_rr_memo_lock)
#else
#define RR_LOCK()   ((void)0)
#define RR_UNLOCK() ((void)0)
#endif
static struct {
    bool valid, cplx;
    long p_bits;
    mpfr_prec_t out_prec;
    mpfr_t zr, zi, rre, rim, dre, dim;
} g_rr_memo;

/* Work precision (bits of |L| log2 e) above which a pass is memoised and
 * computes R' alongside R. */
#define RR_MEMO_WORK_BITS 4096

static bool rr_same_point(const mpfr_t a, const mpfr_t b) {
    return mpfr_get_prec(a) == mpfr_get_prec(b) && mpfr_equal_p(a, b);
}

/* On a hit, initialise and fill whichever of (rre, rim) / (dre, dim) are
 * non-NULL. */
static bool rr_memo_get(const mpfr_t zr, const mpfr_t zi, bool cplx, long p_bits,
                        mpfr_prec_t out_prec, mpfr_ptr rre, mpfr_ptr rim,
                        mpfr_ptr dre, mpfr_ptr dim) {
    bool hit = false;
    RR_LOCK();
    if (g_rr_memo.valid && g_rr_memo.cplx == cplx && g_rr_memo.p_bits == p_bits &&
        g_rr_memo.out_prec == out_prec && rr_same_point(g_rr_memo.zr, zr) &&
        (!cplx || rr_same_point(g_rr_memo.zi, zi))) {
        if (rre) {
            mpfr_init2(rre, out_prec); mpfr_set(rre, g_rr_memo.rre, MPFR_RNDN);
            mpfr_init2(rim, out_prec); mpfr_set(rim, g_rr_memo.rim, MPFR_RNDN);
        }
        if (dre) {
            mpfr_init2(dre, out_prec); mpfr_set(dre, g_rr_memo.dre, MPFR_RNDN);
            mpfr_init2(dim, out_prec); mpfr_set(dim, g_rr_memo.dim, MPFR_RNDN);
        }
        hit = true;
    }
    RR_UNLOCK();
    return hit;
}

static void rr_memo_put(const mpfr_t zr, const mpfr_t zi, bool cplx, long p_bits,
                        mpfr_prec_t out_prec, const mpfr_t rre, const mpfr_t rim,
                        const mpfr_t dre, const mpfr_t dim) {
    RR_LOCK();
    if (g_rr_memo.valid)
        mpfr_clears(g_rr_memo.zr, g_rr_memo.zi, g_rr_memo.rre, g_rr_memo.rim,
                    g_rr_memo.dre, g_rr_memo.dim, (mpfr_ptr)0);
    mpfr_init2(g_rr_memo.zr, mpfr_get_prec(zr)); mpfr_set(g_rr_memo.zr, zr, MPFR_RNDN);
    if (cplx) { mpfr_init2(g_rr_memo.zi, mpfr_get_prec(zi)); mpfr_set(g_rr_memo.zi, zi, MPFR_RNDN); }
    else      { mpfr_init2(g_rr_memo.zi, 16); mpfr_set_ui(g_rr_memo.zi, 0, MPFR_RNDN); }
    mpfr_init2(g_rr_memo.rre, out_prec); mpfr_set(g_rr_memo.rre, rre, MPFR_RNDN);
    mpfr_init2(g_rr_memo.rim, out_prec); mpfr_set(g_rr_memo.rim, rim, MPFR_RNDN);
    mpfr_init2(g_rr_memo.dre, out_prec); mpfr_set(g_rr_memo.dre, dre, MPFR_RNDN);
    mpfr_init2(g_rr_memo.dim, out_prec); mpfr_set(g_rr_memo.dim, dim, MPFR_RNDN);
    g_rr_memo.cplx = cplx; g_rr_memo.p_bits = p_bits; g_rr_memo.out_prec = out_prec;
    g_rr_memo.valid = true;
    RR_UNLOCK();
}

/* R'(z) = D / (z L) into (dre, dim), initialised here at out_prec. */
static void rr_finish_derivative(const mpfr_t zr, const mpfr_t zi, bool cplx,
                                 const mpfr_t Lr, const mpfr_t Li, const mpfr_t Dr,
                                 const mpfr_t Di, mpfr_prec_t out_prec, mpfr_t dre, mpfr_t dim) {
    mpfr_t a, b, c, d, den;
    mpfr_inits2(mpfr_get_prec(Dr), a, b, c, d, den, (mpfr_ptr)0);
    mpfr_init2(dre, out_prec); mpfr_init2(dim, out_prec);
    if (cplx) {
        /* z L = (zr Lr - zi Li) + i (zr Li + zi Lr) */
        mpfr_mul(a, zr, Lr, MPFR_RNDN); mpfr_mul(b, zi, Li, MPFR_RNDN);
        mpfr_sub(c, a, b, MPFR_RNDN);
        mpfr_mul(a, zr, Li, MPFR_RNDN); mpfr_mul(b, zi, Lr, MPFR_RNDN);
        mpfr_add(d, a, b, MPFR_RNDN);
        /* (Dr + i Di) / (c + i d) */
        mpfr_sqr(den, c, MPFR_RNDN); mpfr_sqr(a, d, MPFR_RNDN);
        mpfr_add(den, den, a, MPFR_RNDN);
        mpfr_mul(a, Dr, c, MPFR_RNDN); mpfr_mul(b, Di, d, MPFR_RNDN);
        mpfr_add(a, a, b, MPFR_RNDN);
        mpfr_div(dre, a, den, MPFR_RNDN);
        mpfr_mul(a, Di, c, MPFR_RNDN); mpfr_mul(b, Dr, d, MPFR_RNDN);
        mpfr_sub(a, a, b, MPFR_RNDN);
        mpfr_div(dim, a, den, MPFR_RNDN);
    } else {
        mpfr_mul(den, zr, Lr, MPFR_RNDN);
        mpfr_div(dre, Dr, den, MPFR_RNDN);
        mpfr_set_ui(dim, 0, MPFR_RNDN);
    }
    mpfr_clears(a, b, c, d, den, (mpfr_ptr)0);
}

/* R(z) into (rre, rim) and/or R'(z) into (dre, dim) -- either pair may be
 * NULL -- at an MPFR point z = zr + i zi (zi NULL: real), to p_bits bits
 * relative to the requested quantity, outputs initialised at `out_prec`.
 * Adaptive in the guard: a cancellation deeper than budgeted is recomputed.
 * Returns false when the series cannot be summed (absurd size). */
static bool rr_eval_mpfr(const mpfr_t zr, const mpfr_t zi, long p_bits, mpfr_prec_t out_prec,
                         mpfr_ptr rre, mpfr_ptr rim, mpfr_ptr dre, mpfr_ptr dim) {
    bool cplx = zi != NULL && !(mpfr_zero_p(zi) && mpfr_sgn(zr) > 0);
    if (rr_memo_get(zr, zi, cplx, p_bits, out_prec, rre, rim, dre, dim)) return true;

    /* |L| at double accuracy, for planning. */
    double aL_est;
    {
        mpfr_t tmp; mpfr_init2(tmp, 64);
        if (cplx) mpfr_hypot(tmp, zr, zi, MPFR_RNDN); else mpfr_abs(tmp, zr, MPFR_RNDN);
        mpfr_log(tmp, tmp, MPFR_RNDN);
        aL_est = fabs(mpfr_get_d(tmp, MPFR_RNDN));
        if (cplx) aL_est = hypot(aL_est, M_PI);
        mpfr_clear(tmp);
    }
    bool positive_est = !cplx && mpfr_cmp_ui(zr, 1) > 0;
    /* The first guard anticipates the smooth decay R(e^-u) ~ 21.9/u^3 of the
     * cancelling case, so the common call needs one pass, not two. */
    long guard = 24;
    if (!positive_est && aL_est > 2.0) guard += (long)(3.0 * log2(aL_est)) + 8;
    bool big = !positive_est && aL_est * 1.4427 > (double)RR_MEMO_WORK_BITS;
    bool want_d = dre != NULL || big;
    bool judge_d = dre != NULL && rre == NULL;   /* accuracy asked of R' */

    for (int attempt = 0; attempt < 6; attempt++) {
        /* L = Log z needs only the bits of the answer plus guard: the computed
         * sum is the exact Gram sum at the rounded L, and R(L') - R(L) =
         * R'(L) (L' - L) is small even where the TERMS are huge. */
        mpfr_prec_t lp = (mpfr_prec_t)(p_bits + guard + 64 + (long)log2(aL_est + 2.0));
        mpfr_t Lr, Li;
        mpfr_inits2(lp, Lr, Li, (mpfr_ptr)0);
        if (cplx) {
            mpfr_t ab; mpfr_init2(ab, lp + 16);
            mpfr_hypot(ab, zr, zi, MPFR_RNDN);
            mpfr_log(Lr, ab, MPFR_RNDN);
            mpfr_atan2(Li, zi, zr, MPFR_RNDN);
            mpfr_clear(ab);
        } else {
            mpfr_log(Lr, zr, MPFR_RNDN);
            mpfr_set_ui(Li, 0, MPFR_RNDN);
        }
        mpfr_t Sr, Si, Dr, Di;
        bool ok = rr_gram_mpfr(Lr, Li, cplx, p_bits, guard, want_d, Sr, Si, Dr, Di);
        if (!ok) { mpfr_clears(Lr, Li, (mpfr_ptr)0); return false; }
        mpfr_add_ui(Sr, Sr, 1, MPFR_RNDN);            /* R = 1 + S */

        bool positive = !cplx && mpfr_sgn(Lr) > 0;
        double lr = judge_d ? rr_log2abs(Dr, Di, cplx) : rr_log2abs(Sr, Si, cplx);
        bool enough = positive || (lr > -1e299 && lr + (double)guard >= 8.0);
        if (enough || attempt == 5) {
            mpfr_t xre, xim, xdre, xdim;
            mpfr_init2(xre, out_prec); mpfr_set(xre, Sr, MPFR_RNDN);
            mpfr_init2(xim, out_prec);
            if (cplx) mpfr_set(xim, Si, MPFR_RNDN); else mpfr_set_ui(xim, 0, MPFR_RNDN);
            if (want_d) rr_finish_derivative(zr, zi, cplx, Lr, Li, Dr, Di, out_prec, xdre, xdim);
            if (big) rr_memo_put(zr, zi, cplx, p_bits, out_prec, xre, xim, xdre, xdim);
            if (rre) {
                mpfr_init2(rre, out_prec); mpfr_swap(rre, xre);
                mpfr_init2(rim, out_prec); mpfr_swap(rim, xim);
            }
            if (dre) {
                mpfr_init2(dre, out_prec); mpfr_swap(dre, xdre);
                mpfr_init2(dim, out_prec); mpfr_swap(dim, xdim);
            }
            mpfr_clears(xre, xim, (mpfr_ptr)0);
            if (want_d) mpfr_clears(xdre, xdim, (mpfr_ptr)0);
            mpfr_clears(Sr, Si, Dr, Di, Lr, Li, (mpfr_ptr)0);
            return true;
        }
        /* |R| ~ 2^lr fell below the budget: widen the guard so that p_bits
         * survive relative to |R| (an exact zero: 64 more bits). */
        guard += (lr > -1e299) ? (long)(8.0 - lr) + 16 : 64;
        mpfr_clears(Sr, Si, Dr, Di, Lr, Li, (mpfr_ptr)0);
    }
    return false;
}

/* The MPFR value of a real leaf (Integer/Rational/Real/MPFR) at `prec`. */
static bool rr_leaf_mpfr(const Expr* e, mpfr_t out, mpfr_prec_t prec) {
    mpfr_init2(out, prec);
    if (e->type == EXPR_MPFR) { mpfr_set(out, e->data.mpfr, MPFR_RNDN); return true; }
    if (e->type == EXPR_REAL) { mpfr_set_d(out, e->data.real, MPFR_RNDN); return true; }
    if (e->type == EXPR_INTEGER) { mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); return true; }
    if (e->type == EXPR_BIGINT) { mpfr_set_z(out, e->data.bigint, MPFR_RNDN); return true; }
    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = (long)prec;
    Expr* v = numericalize((Expr*)e, spec);
    bool ok = false;
    if (v && v->type == EXPR_MPFR) { mpfr_set(out, v->data.mpfr, MPFR_RNDN); ok = true; }
    else if (v && v->type == EXPR_REAL) { mpfr_set_d(out, v->data.real, MPFR_RNDN); ok = true; }
    expr_free(v);
    if (!ok) mpfr_clear(out);
    return ok;
}

/* Classify z as an inexact number: machine (53) or MPFR at the minimum width
 * of its parts.  Returns 0 when z is not an inexact number. */
static mpfr_prec_t rr_inexact_prec(Expr* z, bool* machine, Expr** re, Expr** im) {
    *machine = false; *re = z; *im = NULL;
    Expr *cr, *ci;
    if (is_complex(z, &cr, &ci)) { *re = cr; *im = ci; }
    mpfr_prec_t p = 0;
    bool any = false;
    for (int i = 0; i < 2; i++) {
        Expr* e = i ? *im : *re;
        if (!e) continue;
        if (e->type == EXPR_REAL) { *machine = true; any = true; }
        else if (e->type == EXPR_MPFR) {
            any = true;
            mpfr_prec_t q = mpfr_get_prec(e->data.mpfr);
            if (p == 0 || q < p) p = q;
        }
    }
    if (!any) return 0;
    if (*machine) return 53;
    return p;
}

/* RiemannR (or its derivative, deriv = true) at an inexact z, at z's precision.
 * Returns NULL when z is not an inexact number. */
static Expr* rr_numeric(Expr* z, bool deriv) {
    bool machine;
    Expr *re, *im;
    mpfr_prec_t prec = rr_inexact_prec(z, &machine, &re, &im);
    if (prec == 0) return NULL;
    long p_bits = machine ? 53 : (long)prec;
    mpfr_prec_t wprec = machine ? 64 : prec;
    mpfr_t zr, zi;
    if (!rr_leaf_mpfr(re, zr, wprec + 16)) return NULL;
    bool have_im = im != NULL;
    if (have_im && !rr_leaf_mpfr(im, zi, wprec + 16)) { mpfr_clear(zr); return NULL; }
    if (!have_im) mpfr_init2(zi, 16), mpfr_set_ui(zi, 0, MPFR_RNDN);

    Expr* out = NULL;
    if (mpfr_zero_p(zr) && mpfr_zero_p(zi)) {
        /* R(0+) = 0; R'(x) ~ 1/(x L) blows up, leave the derivative alone. */
        if (!deriv) out = machine ? expr_new_real(0.0) : expr_new_mpfr_bits(prec);
    } else if (mpfr_number_p(zr) && mpfr_number_p(zi)) {
        mpfr_t rre, rim, dre, dim;
        mpfr_prec_t op = machine ? 53 : prec;
        bool cplx_in = have_im || mpfr_sgn(zr) < 0;
        if (rr_eval_mpfr(zr, cplx_in ? zi : NULL, p_bits, op,
                         deriv ? NULL : rre, deriv ? NULL : rim,
                         deriv ? dre : NULL, deriv ? dim : NULL)) {
            mpfr_ptr vr = deriv ? dre : rre, vi = deriv ? dim : rim;
            if (machine) {
                double a = mpfr_get_d(vr, MPFR_RNDN), b = mpfr_get_d(vi, MPFR_RNDN);
                if (isfinite(a) && isfinite(b))
                    out = cplx_in ? sf_machine_complex_result(a, b) : expr_new_real(a);
            } else {
                Expr* er = expr_new_mpfr_copy(vr);
                out = cplx_in ? make_complex(er, expr_new_mpfr_copy(vi)) : er;
            }
            if (deriv) mpfr_clears(dre, dim, (mpfr_ptr)0);
            else       mpfr_clears(rre, rim, (mpfr_ptr)0);
        }
    }
    mpfr_clears(zr, zi, (mpfr_ptr)0);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Builtins                                                           */
/* ------------------------------------------------------------------ */

Expr* builtin_riemannr(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("RiemannR", argc, 1, 1);
    Expr* z = res->data.function.args[0];

    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma)) {
        double vr, vi;
        if (!ma.any_complex) {
            if (sf_riemannr_r(ma.re[0], &vr)) return expr_new_real(vr);
        } else if (sf_riemannr_c(ma.re[0], ma.im[0], &vr, &vi)) {
            return sf_machine_complex_result(vr, vi);
        }
        /* declined: the MPFR path at 53 bits, rounded, below */
    }
    if (z->type == EXPR_INTEGER && z->data.integer == 1) return expr_new_integer(1);
    if (z->type == EXPR_SYMBOL && z->data.symbol.name == SYM_Infinity)
        return expr_new_symbol(SYM_Infinity);
#ifdef USE_MPFR
    return rr_numeric(z, false);
#else
    return NULL;
#endif
}

/* Derivative[1][RiemannR][z] at an inexact z: R'(z) from the same Gram pass.
 * Reached through the internal head below (see riemannr_init); symbolic z
 * keeps the inert Derivative[1][RiemannR][z], as in the reference. */
Expr* builtin_riemannr_prime(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* z = res->data.function.args[0];
    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && !ma.any_complex && ma.re[0] > 0.0 &&
        isfinite(ma.re[0]) && ma.re[0] != 1.0) {
        rr_dd L = rr_log_dd(ma.re[0]);
        double r, d;
        if (rr_gram_dd(L, &r, &d)) return expr_new_real(d / (ma.re[0] * (L.hi + L.lo)));
    }
#ifdef USE_MPFR
    return rr_numeric(z, true);
#else
    (void)z;
    return NULL;
#endif
}

/* Derivative[1][RiemannR][z] for an inexact number z (else NULL: inert). */
static Expr* rr_derivative_value(const Expr* call) {
    const Expr* dop = call->data.function.head;
    const Expr* d = dop->data.function.head;          /* Derivative[n] */
    if (d->data.function.arg_count != 1 || d->data.function.args[0]->type != EXPR_INTEGER ||
        d->data.function.args[0]->data.integer != 1 || call->data.function.arg_count != 1)
        return NULL;
    Expr* z = call->data.function.args[0];
    /* an inexact number: a Real / MPFR leaf or a Complex of them */
    bool num = z->type == EXPR_REAL
#ifdef USE_MPFR
               || z->type == EXPR_MPFR
#endif
               || (z->type == EXPR_FUNCTION && z->data.function.head->type == EXPR_SYMBOL &&
                   z->data.function.head->data.symbol.name == SYM_Complex);
    if (!num || !arg_is_inexact(z)) return NULL;
    Expr* a[1] = { expr_copy(z) };
    return eval_and_free(expr_new_function(expr_new_symbol(intern_symbol("Internal`RiemannRPrime")), a, 1));
}

void riemannr_init(void) {
#ifdef USE_MPFR
    rr_mu_sieve();
#endif
    symtab_add_builtin("RiemannR", builtin_riemannr);
    symtab_get_def("RiemannR")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_add_builtin("Internal`RiemannRPrime", builtin_riemannr_prime);
    /* An implementation detail reached only through the derivative registry,
     * not a user numeric head: no NumericFunction (so no Compile surface is
     * owed; Compile of Derivative[1][RiemannR] is the inert derivative). */
    symtab_get_def("Internal`RiemannRPrime")->attributes |= ATTR_PROTECTED;
    /* Derivative[1][RiemannR][z] at an inexact z has a value (R'(z)), while
     * D[RiemannR[x], x] stays the inert Derivative[1][RiemannR][x], as in the
     * reference.  Registered with the derivative registry (deriv.h), so it is
     * a system rule and SubValues[Derivative] stays empty. */
    deriv_register_numeric(intern_symbol("RiemannR"), rr_derivative_value, 1);
    symtab_set_docstring("RiemannR",
        "RiemannR[x]\n"
        "\tRiemann's prime-counting function R(x) = Sum_k MoebiusMu[k]/k LogIntegral[x^(1/k)],\n"
        "\tevaluated by the Gram series 1 + Sum_k Log[x]^k/(k k! Zeta[k+1]).\n"
        "RiemannR[1] = 1; other exact arguments stay unevaluated. Machine and arbitrary-\n"
        "precision real and complex x evaluate numerically (complex x and x < 0 through\n"
        "the principal Log); the working precision grows with the cancellation for\n"
        "0 < x < 1, so R near its zeros (x ~ 10^-14828) is still correct to the\n"
        "precision asked for. Derivative[1][RiemannR] evaluates at inexact points.\n"
        "Machine, packed, NDArray and Compile[] share one kernel. Listable, NumericFunction.");
}
