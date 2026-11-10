/* Machine-precision (double) kernels for the Anger-Weber and Kelvin functions.
 *
 *   sf_machine_anger_weber(0, nu, z)   AngerJ  J_nu(z) = (1/pi) Int_0^pi cos(nu t - z sin t) dt
 *   sf_machine_anger_weber(1, nu, z)   WeberE  E_nu(z) = (1/pi) Int_0^pi sin(nu t - z sin t) dt
 *   sf_machine_kelvin(0..3, nu, x)     ber, bei, ker, kei
 *
 * Each entry point is shared bit for bit by the scalar builtin (angerweber.c,
 * kelvin.c), the NDArray / packed element-wise kernels (ndkernels.c) and the
 * Compile[] VM, which calls the same NDBinaryKernel. Contract as sf_machine.h:
 * a kernel returns FALSE when it cannot deliver a machine answer (the value is
 * complex, infinite, overflows, or neither algorithm meets its error budget);
 * the caller then answers through the arb bridge.
 *
 * Real arithmetic only (no <complex.h>); the Kelvin functions reuse the raw
 * (re, im) complex Bessel kernels of sf_machine_bessel.c.
 */
#include "sf_machine.h"
#include "sf_machine_bessel.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

/* sin(pi r), cos(pi r) with EXACT values at the multiples of 1/2, so an
 * integer order gets sin(n pi) = 0 rather than 1.2e-16 (AngerJ[n, z] must be
 * BesselJ[n, z] bit for bit, and the Kelvin phase e^(-i n pi/2) exactly real or
 * imaginary). The reduction r mod 2 is exact (fmod is). */
static void aw_sincospi(double r, double* s, double* c) {
    double m = fmod(r, 2.0);
    if (m < 0.0) m += 2.0;
    double h = 2.0 * m;
    if (h == floor(h)) {
        static const double S[4] = { 0.0, 1.0, 0.0, -1.0 };
        static const double C[4] = { 1.0, 0.0, -1.0, 0.0 };
        int q = (int)h & 3;
        *s = S[q]; *c = C[q];
        return;
    }
    /* reduce to |x| <= 1/2 before scaling by pi (keeps both small values accurate) */
    if (m > 1.0) { aw_sincospi(m - 1.0, s, c); *s = -*s; *c = -*c; return; }
    if (m > 0.5) { double s2, c2; aw_sincospi(1.0 - m, &s2, &c2); *s = s2; *c = -c2; return; }
    *s = sin(M_PI * m); *c = cos(M_PI * m);
}

/* ------------------------------------------------------------------ */
/*  Anger / Weber: the S1, S2 power series                             */
/* ------------------------------------------------------------------ */
/* Sum_{k>=0} (-w)^k / (Gamma(k+p) Gamma(k+q)), w = (z/2)^2 >= 0. A term whose
 * Gamma has a pole (k+p or k+q a non-positive integer) is 0, so the sum starts
 * at the first k where both are finite. Sets the sum and the largest term
 * (the cancellation scale); false on Gamma overflow or non-convergence. */
static bool aw_sum(double p, double q, double w, double hump,
                   double* sum, double* maxabs) {
    int k0 = 0;
    if (p <= 0.0 && p == floor(p)) k0 = (int)(1.0 - p);
    if (q <= 0.0 && q == floor(q) && (int)(1.0 - q) > k0) k0 = (int)(1.0 - q);
    if (k0 > 400) return false;
    double gp = tgamma(k0 + p), gq = tgamma(k0 + q);
    if (!isfinite(gp) || !isfinite(gq) || gp == 0.0 || gq == 0.0) return false;
    double t = 1.0 / gp / gq;
    for (int j = 0; j < k0; j++) t *= -w;
    if (!isfinite(t)) return false;
    double s = t, mx = fabs(t);
    int k;
    for (k = k0 + 1; k < 2000; k++) {
        t *= -w / ((k - 1 + p) * (k - 1 + q));
        s += t;
        double at = fabs(t);
        if (at > mx) mx = at;
        if (at <= fabs(s) * 1.0e-18 && (double)k > hump) break;
        if (t == 0.0) break;
    }
    if (k >= 2000 || !isfinite(s)) return false;
    *sum = s; *maxabs = mx;
    return true;
}

/* DLMF 11.10.8-9 (S1, S2 of 11.10.10-11):
 *   J_nu = cos(nu pi/2) S1 + sin(nu pi/2) S2,  E_nu = sin(nu pi/2) S1 - cos(nu pi/2) S2,
 *   S1 = Sum (-1)^k (z/2)^(2k)   / (Gamma(k+1+nu/2) Gamma(k+1-nu/2)),
 *   S2 = Sum (-1)^k (z/2)^(2k+1) / (Gamma(k+3/2+nu/2) Gamma(k+3/2-nu/2)).
 * *abserr is the predicted ABSOLUTE error, (scale of the largest term) * eps. */
static bool aw_series(int which, double nu, double z, double* val, double* abserr) {
    double w = 0.25 * z * z, hump = 0.5 * fabs(z), a = 0.5 * nu;
    double s1, m1, s2, m2;
    if (!aw_sum(1.0 + a, 1.0 - a, w, hump, &s1, &m1)) return false;
    if (!aw_sum(1.5 + a, 1.5 - a, w, hump, &s2, &m2)) return false;
    double hz = 0.5 * z;
    s2 *= hz; m2 *= fabs(hz);
    double sn, cs;
    aw_sincospi(a, &sn, &cs);
    double v = which == 0 ? cs * s1 + sn * s2 : sn * s1 - cs * s2;
    double scale = fabs(cs) * m1 + fabs(sn) * m2;
    if (!isfinite(v)) return false;
    *val = v;
    *abserr = 2.0 * DBL_EPSILON * scale;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Anger / Weber: the Bessel + Laplace-integral form (z > 0)          */
/* ------------------------------------------------------------------ */
/* For Re z > 0 (DLMF 11.10.4-11.10.6, 11.10.13 / 11.10.14 in A_nu form):
 *   J_nu(z) = J_nu(z) + (sin(nu pi)/pi) A(nu)
 *   E_nu(z) = -Y_nu(z) - (A(-nu) + cos(nu pi) A(nu)) / pi,
 *   A(mu)   = Int_0^inf exp(-mu t - z sinh t) dt.
 * The integrand is ENTIRE in t and decays double-exponentially, so Gauss-
 * Legendre on [0, T] (T where it has fallen e^-40 below its peak) converges
 * geometrically; 32 nodes reach the double floor over the whole range used. */
#define AW_GL_N 12   /* half of the 24-point rule (symmetric nodes) */
static const double AW_GL_X[AW_GL_N] = {
    0.064056892862605626085, 0.191118867473616309159, 0.315042679696163374387,
    0.433793507626045138487, 0.545421471388839535658, 0.648093651936975569252,
    0.740124191578554364244, 0.820001985973902921954, 0.886415527004401034213,
    0.938274552002732758524, 0.974728555971309498198, 0.99518721999702136018 };
static const double AW_GL_W[AW_GL_N] = {
    0.127938195346752156974, 0.125837456346828296121, 0.121670472927803391204,
    0.115505668053725601353, 0.107444270115965634783, 0.0976186521041138882699,
    0.0861901615319532759172, 0.073346481411080305734, 0.0592985849154367807464,
    0.0442774388174198061686, 0.0285313886289336631813, 0.0123412297999871995468 };

/* Upper cut-off T for exp(-mu t - z sinh t): past the peak t_p (cosh t_p =
 * -mu/z when mu < -z, else 0), Newton on the concave decreasing
 * g(T) = -mu T - z sinh T - (g_peak - 40), which converges from the right. */
static double aw_cutoff(double mu, double z) {
    double tp = (-mu > z) ? acosh(-mu / z) : 0.0;
    double target = -mu * tp - z * sinh(tp) - 40.0;
    double T = tp + 1.0;
    for (int i = 0; i < 40; i++) {
        double g = -mu * T - z * sinh(T) - target;
        double dg = -mu - z * cosh(T);
        if (!(dg < 0.0)) { T += 1.0; continue; }
        double Tn = T - g / dg;
        if (Tn <= tp) Tn = 0.5 * (T + tp);
        if (fabs(Tn - T) < 1.0e-3 * T) { T = Tn; break; }
        T = Tn;
    }
    return T > tp ? T : tp + 1.0;
}

/* A(nu) and (when want_neg) A(-nu) on [lo, lo + 2h] by the 24-point rule.
 * The symmetric node pair t = c +- h x_i shares its exponentials: e^t =
 * e^c e^(+-h x_i) and e^(-nu t) = e^(-nu c) e^(-+nu h x_i), so a pair costs
 * four exp() calls (two of them e^(-z sinh t)) instead of six. */
static void aw_laplace_panel(double nu, double z, bool want_neg, double lo, double h,
                             double* sp, double* sm) {
    double c = lo + h;
    double ec = exp(c), gc = exp(-nu * c);
    double ap = 0.0, am = 0.0;
    for (int i = 0; i < AW_GL_N; i++) {
        double f = exp(h * AW_GL_X[i]), g = exp(-nu * h * AW_GL_X[i]);
        double etp = ec * f, etm = ec / f;               /* e^t at the pair     */
        double enp = gc * g, enm = gc / g;               /* e^(-nu t)           */
        double ep = exp(-0.5 * z * (etp - 1.0 / etp));   /* e^(-z sinh t)       */
        double em = exp(-0.5 * z * (etm - 1.0 / etm));
        double w = AW_GL_W[i];
        ap += w * (enp * ep + enm * em);
        if (want_neg) am += w * (ep / enp + em / enm);
    }
    *sp = h * ap;
    *sm = h * am;
}

/* A(nu) and (when want_neg) A(-nu), with a MEASURED quadrature error: the
 * one-panel rule on [0, T] against the two-panel rule on its halves.  The
 * 24 nodes alone are not always enough -- a negative order at small z puts
 * structure near t = 0 that one panel misses (AngerJ[-5/6, 5/6] came out
 * -0.24730 for -0.24614, against a claimed error of 2e-16, and the regime
 * choice took it over the correct series).  *ep_err / *em_err bound the
 * error of the two-panel values returned. */
static bool aw_laplace(double nu, double z, bool want_neg, double* ap, double* am,
                       double* ap_err, double* am_err) {
    double T = aw_cutoff(nu, z);
    if (want_neg) { double T2 = aw_cutoff(-nu, z); if (T2 > T) T = T2; }
    if (!isfinite(T) || T > 700.0) return false;
    double p1, m1, pa, ma, pb, mb;
    aw_laplace_panel(nu, z, want_neg, 0.0, 0.5 * T, &p1, &m1);
    aw_laplace_panel(nu, z, want_neg, 0.0, 0.25 * T, &pa, &ma);
    aw_laplace_panel(nu, z, want_neg, 0.5 * T, 0.25 * T, &pb, &mb);
    *ap = pa + pb;
    *am = ma + mb;
    *ap_err = fabs(*ap - p1);
    *am_err = fabs(*am - m1);
    return isfinite(*ap) && isfinite(*am);
}

/* *abserr is the predicted absolute error and *env the size of the components
 * (|J_nu| or |Y_nu| plus the integral terms): the function's local envelope,
 * against which an answer near one of its zeros is judged. */
static bool aw_integral(int which, double nu, double z, double* val, double* abserr,
                        double* env) {
    double sn, cs;
    aw_sincospi(nu, &sn, &cs);
    double ap = 0.0, am = 0.0, pe = 0.0, me = 0.0, b, v, scale, qerr;
    if (which == 0) {
        if (!sf_bessel_j_r(nu, z, &b)) return false;
        if (sn != 0.0 && !aw_laplace(nu, z, false, &ap, &am, &pe, &me)) return false;
        v = b + sn / M_PI * ap;
        scale = fabs(b) + fabs(sn / M_PI * ap);
        qerr = fabs(sn) / M_PI * pe;
    } else {
        if (!sf_bessel_y_r(nu, z, &b)) return false;
        if (!aw_laplace(nu, z, true, &ap, &am, &pe, &me)) return false;
        double c = (am + cs * ap) / M_PI;
        v = -b - c;
        scale = fabs(b) + (fabs(am) + fabs(cs * ap)) / M_PI;
        qerr = (me + fabs(cs) * pe) / M_PI;
    }
    if (!isfinite(v)) return false;
    *val = v;
    /* ~4 eps per component (the Bessel kernel, the quadrature sum), plus the
     * measured quadrature error */
    *abserr = 4.0 * DBL_EPSILON * scale + qerr;
    *env = scale;
    return true;
}

#define AW_ACCEPT   1.0e-13   /* budget, relative to the local envelope              */
#define AW_SERIES_Z 12.0      /* beyond this the series loses > 5 digits anyway       */

/* Regime choice works on ABSOLUTE errors: both functions oscillate through
 * zeros, where no algorithm can deliver relative accuracy, so an answer is
 * judged against the local envelope max(|value|, size of the components) --
 * the standard for an oscillatory special function. (Judging it by relative
 * error declined every element near a zero, which sent a whole packed array
 * down the List path.) */
bool sf_machine_anger_weber(int which, double nu, double z, double* out) {
    if (!isfinite(nu) || !isfinite(z)) return false;
    /* An integer order: the Anger function IS the Bessel function. */
    if (which == 0 && nu == floor(nu) && fabs(nu) < 1.0e9)
        return sf_bessel_j_r(nu, z, out);
    double sv = 0.0, se = HUGE_VAL, env = 0.0;
    bool sok = false;
    if (fabs(z) <= AW_SERIES_Z) {
        sok = aw_series(which, nu, z, &sv, &se);
        if (sok) {
            env = fabs(sv);
            if (se <= 4.0 * DBL_EPSILON * env) { *out = sv; return true; }
        }
    }
    if (z != 0.0) {
        /* parity: J_nu(-z) = J_{-nu}(z), E_nu(-z) = -E_{-nu}(z) */
        double iv, ie, ienv;
        bool neg = z < 0.0;
        if (aw_integral(which, neg ? -nu : nu, fabs(z), &iv, &ie, &ienv)) {
            if (neg && which == 1) iv = -iv;
            if (!sok || ie < se) { sv = iv; se = ie; sok = true; }
            if (ienv > env) env = ienv;
        }
    }
    if (!sok || !(se <= AW_ACCEPT * env)) return false;
    *out = sv;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Kelvin functions                                                   */
/* ------------------------------------------------------------------ */
/* ber_nu(x) + i bei_nu(x) = J_nu(x e^(3 pi i/4))              (DLMF 10.61.1)
 * ker_nu(x) + i kei_nu(x) = e^(-nu pi i/2) K_nu(x e^(pi i/4)) (DLMF 10.61.2)
 * for x > 0. ber/bei of a negative x at an integer order follow from the
 * series, ber_n(-x) = (-1)^n ber_n(x) (likewise bei); every other x < 0 value
 * is complex (or on the ker/kei cut) and declines. x = 0 is answered exactly
 * where finite: ber_0(0) = 1, ber_nu(0) = bei_nu(0) = 0 (nu > 0), bei_0(0) = 0,
 * kei_0(0) = -pi/4. */
bool sf_machine_kelvin(int which, double nu, double x, double* out) {
    if (!isfinite(nu) || !isfinite(x)) return false;
    if (x == 0.0) {
        switch (which) {
        case 0: if (nu == 0.0) { *out = 1.0; return true; }
                if (nu > 0.0) { *out = 0.0; return true; } return false;
        case 1: if (nu >= 0.0) { *out = 0.0; return true; } return false;
        case 3: if (nu == 0.0) { *out = -M_PI / 4.0; return true; } return false;
        default: return false;
        }
    }
    double rr, ri;
    if (which <= 1) {
        double flip = 1.0;
        if (x < 0.0) {
            if (nu != floor(nu) || fabs(nu) > 1.0e15) return false;
            flip = fmod(fabs(nu), 2.0) == 1.0 ? -1.0 : 1.0;
            x = -x;
        }
        if (!sf_bessel_j_c(nu, -M_SQRT1_2 * x, M_SQRT1_2 * x, &rr, &ri)) return false;
        *out = flip * (which == 0 ? rr : ri);
    } else {
        if (x < 0.0) return false;
        if (!sf_bessel_k_c(nu, M_SQRT1_2 * x, M_SQRT1_2 * x, &rr, &ri)) return false;
        double sn, cs;                         /* e^(-nu pi i/2) = cs - i sn */
        aw_sincospi(0.5 * nu, &sn, &cs);
        *out = which == 2 ? cs * rr + sn * ri : cs * ri - sn * rr;
    }
    return isfinite(*out);
}
