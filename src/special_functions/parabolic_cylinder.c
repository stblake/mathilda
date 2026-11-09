/* Mathilda -- ParabolicCylinderD[nu, z], Whittaker's parabolic cylinder function
 * D_nu(z), the solution of  w'' + (nu + 1/2 - z^2/4) w = 0  that decays along
 * the positive real axis. Entire in z for every nu, entire in nu for every z.
 *
 * Two representations, each used where it does not cancel (DLMF 12.7.14 and
 * 12.4.1 / 12.7.12 with a = -nu - 1/2, x = z^2/2):
 *
 *   (U)  D_nu(z) = 2^(nu/2) e^(-z^2/4) U(-nu/2, 1/2, z^2/2)            Re z > 0
 *   (M)  D_nu(z) = 2^(nu/2) Sqrt[Pi] e^(-z^2/4) [ M(-nu/2, 1/2, x) / Gamma((1-nu)/2)
 *                       - Sqrt[2] z M((1-nu)/2, 3/2, x) / Gamma(-nu/2) ]   any z
 *
 * (U) is the principal-branch Tricomi function: x = z^2/2 stays off its cut
 * (-inf, 0] exactly when Re z > 0, and the (z^2)^(1/2) = z that its connection
 * formula produces is then the right branch. In |arg z| < Pi/4, where D decays
 * like z^nu e^(-z^2/4), the two terms of (M) are each of size e^(+z^2/4) and
 * cancel -- (U) has no such cancellation. In Re z <= 0 the two terms of (M) add
 * (both grow like e^(z^2/4) |z|^(-nu-1) along the negative axis, with the same
 * sign), so (M) is used there. 1/Gamma vanishes at its poles, which is what
 * makes (M) collapse to the Hermite form at a non-negative integer nu.
 *
 * Numerics:
 *   machine scalar                  -> sf_machine_pcfd / _c / _cc, the shared
 *                                      double kernels (see "The real kernel"
 *                                      below; complex z rides the 1F1 / U kernels
 *                                      of sf_machine_hyper.c; declines are
 *                                      answered in-kernel by arb at 53 bits);
 *                                      identical to the NDArray / packed map and
 *                                      the Compile VM.
 *   anything else inexact           -> arb (flint_num_pcfd) at the arguments'
 *                                      precision, lowered back to machine numbers
 *                                      for a machine call -- so a value outside
 *                                      the double range comes back as an
 *                                      extended-range real, never inf / 0.
 * Exact rules: D_nu(0) = 2^(nu/2) Sqrt[Pi] / Gamma[(1-nu)/2] for every nu. An
 * exact integer order stays unevaluated, as in the reference system;
 * FunctionExpand rewrites it (Hermite form for n >= 0, Erfc recurrence for n < 0).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "parabolic_cylinder.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#include "attr.h"
#include "common.h"            /* builtin_arg_error */
#include "eval.h"              /* eval_and_free */
#include "expr.h"
#include "flint_num_bridge.h"  /* flint_num_pcfd, flint_num_call_inexact, flint_num_for_call */
#include "sf_machine.h"        /* sf_machine_hyper1f1, sf_machine_hyperu */
#include "sf_machine_hyper.h"  /* complex twins */
#include "sf_scalar.h"         /* sf_machine_scalar_args */
#include "sf_dd.h"             /* double-double re-summation */
#include "symtab.h"
#include "../sym_names.h"

/* ------------------------------------------------------------------ */
/* Machine kernels                                                     */
/* ------------------------------------------------------------------ */

/* Largest |nu| the machine kernels accept. Beyond it the 1F1 / U parameters are
 * large enough that the double kernels lose digits or decline anyway, and
 * tgamma((1 - nu)/2) leaves the double range near nu = -342. */
#define PCFD_NU_MAX 150.0

/* sqrt(pi), sqrt(2) */
#define PCFD_SQRTPI 1.77245385090551602729816748334
#define PCFD_SQRT2  1.41421356237309504880168872421

/* 1/Gamma(x) in double: exactly 0 at the poles x = 0, -1, -2, ... */
static double pcfd_rgamma(double x) {
    if (x <= 0.0 && x == floor(x)) return 0.0;
    if (x < 0.0) {
        /* reflection 1/Gamma(x) = sin(pi x) Gamma(1 - x) / pi, sin(pi x) reduced
         * exactly: tgamma of a large negative argument loses ulps */
        double r = x - 2.0 * floor(0.5 * x + 0.5);        /* r in [-1, 1) */
        if (r > 0.5) r = 1.0 - r;                          /* sin(pi r) symmetric */
        else if (r < -0.5) r = -1.0 - r;
        return sin(M_PI * r) * tgamma(1.0 - x) / M_PI;
    }
    return 1.0 / tgamma(x);
}

/* A double that is finite and not subnormal (subnormals carry < 53 bits, and a
 * zero produced by underflow is not the value). Exact zero passes only when the
 * caller says it is a genuine zero. */
static bool pcfd_normal(double v) {
    return isfinite(v) && (v == 0.0 || fabs(v) >= DBL_MIN);
}

/* Complex (vr + i vi) * exp(er + i ei). */
static void pcfd_scale_c(double vr, double vi, double er, double ei,
                         double* rr, double* ri) {
    double c = cos(ei), s = sin(ei);
    double pr = vr * c - vi * s, pi = vr * s + vi * c;
    if (fabs(er) < 700.0) {
        double m = exp(er);
        *rr = pr * m; *ri = pi * m;
    } else {
        double h = exp(0.5 * er);
        *rr = (pr * h) * h; *ri = (pi * h) * h;
    }
}

/* ------------------------------------------------------------------ */
/* The real kernel                                                     */
/* ------------------------------------------------------------------ */
/* For real nu (|nu| <= 40) and real z, in roughly 0.2 - 5 us:
 *
 *   z >= 9, expansion converges
 *                      Poincare expansion (DLMF 12.9.1)
 *                        D = z^nu e^(-z^2/4) sum_s (-1)^s (-nu)_{2s} / (s! (2 z^2)^s)
 *   1.5 < z            the recessive solution integrated BACKWARD from a point
 *                      z_A where the expansion is accurate, by Taylor steps of
 *                      the ODE  w'' = (z^2/4 - nu - 1/2) w.  Going towards 0 the
 *                      recessive solution grows (or oscillates) relative to the
 *                      dominant one, so the stepping is stable; the step lands on
 *                      exactly representable points and the coefficient
 *                      z^2/4 - nu - 1/2 is carried as hi + lo, since either
 *                      rounding would shift the solution by D'/D ulp.
 *   -12 <= z <= 1.5    the Maclaurin form (M) summed directly, in double and --
 *                      when the series itself cancels (nu > 0 makes (-nu/2)_k
 *                      alternate) -- in double-double. For z <= 0 the two terms
 *                      of (M) add except near a zero of D, so the Gamma factors
 *                      (double, by reflection for negative arguments) are not
 *                      amplified; a cancelling pair declines.
 *   z < -12            the connection formula (DLMF 12.2.15 with the expansion
 *                      of D_{-nu-1}(-i t), t = -z):
 *                        D_nu(-t) = cos(pi nu) D_nu(t)
 *                          + Sqrt[2 pi]/Gamma(-nu) t^(-nu-1) e^(t^2/4)
 *                              sum_s (nu+1)_{2s} / (s! (2 t^2)^s)
 * Every exponential is applied through an exactly reduced exponent
 * (pcfd_expscale). Every branch measures its own trouble (a non-converging
 * expansion, a cancelling sum) and declines; what remains -- a point near a zero
 * on the negative axis, |nu| > 40, and the complex arm's declines -- is answered
 * by arb at 53 bits (flint_pcfd_d) INSIDE the kernel, so the kernel only
 * declines a value outside the normal double range and every surface (scalar,
 * packed, NDArray, Compile) answers identically and stays on the buffer. Measured
 * against mpmath over 2300 points: real median 0.9 ulp, p99 18 ulp. */

/* cos(pi nu), exact at integers and half-integers. */
static double pcfd_cospi(double nu) {
    double r = nu - 2.0 * floor(0.5 * nu + 0.5);      /* r in [-1, 1) */
    if (r == 0.5 || r == -0.5) return 0.0;
    return cos(M_PI * r);
}

/* The Poincare sum S(z) and its z-derivative S'(z) of
 * D_nu(z) = z^nu e^(-z^2/4) S(z), z > 0. False when the expansion does not
 * reach ~2^-56 relative accuracy before its terms start to grow, or when its
 * terms cancel by more than 2^4. */
static bool pcfd_asym_pos(double nu, double z, double* S, double* Sp) {
    double w = 1.0 / (2.0 * z * z);
    double t = 1.0, s = 1.0, sp = 0.0, peak = 1.0, prev = 1.0;
    bool dec = false;
    for (int k = 1; k < 400; k++) {
        double f = (2.0 * k - 2.0 - nu) * (2.0 * k - 1.0 - nu);
        t *= -f * w / k;
        if (t == 0.0) { *S = s; *Sp = sp; return true; }   /* terminating */
        double at = fabs(t);
        if (at > prev && dec) return false;      /* an asymptotic series: stop at its least term */
        if (at < prev) dec = true;
        s += t;
        sp += t * (-2.0 * k / z);
        if (at > peak) peak = at;
        prev = at;
        if (at < 1.4e-17 * fabs(s)) break;
        if (k == 399) return false;
    }
    if (peak > 16.0 * fabs(s)) return false;
    *S = s; *Sp = sp;
    return true;
}

/* v * 2^E * exp(eh + el), the exponent reduced exactly (Cody-Waite ln 2), so a
 * large exponent such as z^2/4 costs no digits: exp(z^2/4) computed from a
 * rounded z^2/4 would carry a relative error of z^2/4 ulp. */
static double pcfd_expscale(double v, int E, double eh, double el) {
    const double LN2HI = 6.93147180369123816490e-01, LN2LO = 1.90821492927058770002e-10;
    if (v == 0.0) return 0.0;
    double k = floor(eh / M_LN2 + 0.5);
    if (fabs(k) > 4000.0) return eh > 0 ? INFINITY : 0.0;
    double r = (eh - k * LN2HI) - k * LN2LO + el;
    return ldexp(v * exp(r), (int)k + E);
}

/* -z^2/4 as an exact double-double (hi, lo). */
static void pcfd_mquarter_sq(double z, double* hi, double* lo) {
    sdd p = sdd_two_prod(z, z);
    *hi = -0.25 * p.hi; *lo = -0.25 * p.lo;
}

/* 1/((k+2)(k+1)), k = 0..77: the Taylor recurrence divisor, tabulated */
static const double PCFD_INV[78] = {
    0.5, 0.16666666666666666, 0.083333333333333329, 0.050000000000000003,
    0.033333333333333333, 0.023809523809523808, 0.017857142857142856, 0.013888888888888888,
    0.011111111111111112, 0.0090909090909090905, 0.007575757575757576, 0.00641025641025641,
    0.0054945054945054949, 0.0047619047619047623, 0.0041666666666666666, 0.0036764705882352941,
    0.0032679738562091504, 0.0029239766081871343, 0.002631578947368421, 0.0023809523809523812,
    0.0021645021645021645, 0.001976284584980237, 0.0018115942028985507, 0.0016666666666666668,
    0.0015384615384615385, 0.0014245014245014246, 0.0013227513227513227, 0.0012315270935960591,
    0.0011494252873563218, 0.0010752688172043011, 0.0010080645161290322, 0.000946969696969697,
    0.00089126559714795004, 0.00084033613445378156, 0.00079365079365079365, 0.00075075075075075074,
    0.00071123755334281653, 0.00067476383265856947, 0.00064102564102564103, 0.00060975609756097561,
    0.00058072009291521487, 0.00055370985603543741, 0.00052854122621564484, 0.00050505050505050505,
    0.00048309178743961351, 0.00046253469010175765, 0.00044326241134751772, 0.00042517006802721087,
    0.00040816326530612246, 0.00039215686274509802, 0.00037707390648567121, 0.00036284470246734398,
    0.00034940600978336826, 0.00033670033670033672, 0.00032467532467532468, 0.00031328320802005011,
    0.00030248033877797946, 0.00029222676797194621, 0.0002824858757062147, 0.00027322404371584699,
    0.00026441036488630354, 0.00025601638504864311, 0.000248015873015873, 0.0002403846153846154,
    0.0002331002331002331, 0.00022614201718679331, 0.00021949078138718174, 0.00021312872975277067,
    0.00020703933747412008, 0.00020120724346076458, 0.00019561815336463224, 0.00019025875190258751,
    0.00018511662347278786, 0.00018018018018018018, 0.00017543859649122806, 0.00017088174982911826,
    0.0001665001665001665, 0.00016228497241155469,
};
/* Backward integration of the recessive solution from z_A down to z (0 < z <
 * z_A). (w, wp) = (D, D') * 2^-E * exp(-(eh + el)) at z_A. */
static bool pcfd_integrate(double nu, double zA, double w, double wp, int E,
                           double eh, double el, double z, double* out) {
    double z0 = zA;
    sdd En = sdd_two_sum(0.5, nu);              /* nu + 1/2, exactly */
    double c[80], lmax = -INFINITY;
    while (z0 > z) {
        /* q0 = z0^2/4 - nu - 1/2 as hi + lo: its rounding would shift the
         * solution by ~|q0| ulp per unit length */
        sdd zz = sdd_two_prod(z0, z0);
        sdd qd = sdd_sub(sdd_add_d(sdd_from(0.25 * zz.lo), 0.25 * zz.hi), En);
        double q0 = qd.hi, q0l = qd.lo;
        /* step length in units of the local scale 1/sqrt|q|: long where the
         * solution is exponential (its Taylor terms then share one sign), short
         * where it oscillates (they alternate and cancel like e^(kappa h)) */
        double kh = (q0 > 0.0 && z0 > 0.0) ? 5.0 : 2.0;
        double h = -fmin(2.0, kh / sqrt(fabs(q0) + 0.5 * fabs(z0) + 1.0));
        /* land on an exactly representable point and step by the EXACT
         * difference: a rounded z0 + h shifts the solution by D'/D ulp */
        double zn = z0 + h;
        if (z0 > 0.0 && zn > 0.0 && zn < 0.5 * z0) zn = 0.5 * z0;   /* Sterbenz */
        if (zn < z || zn - z < 1e-3) zn = z;
        h = zn - z0;
        /* Taylor coefficients about z0: w'' = (q0 + (z0/2) t + t^2/4) w */
        c[0] = w; c[1] = wp;
        double val = w + wp * h, der = wp;
        double hp = h, mag = fabs(w) + fabs(wp * h);
        int k, small = 0;
        for (k = 0; k + 2 < 80; k++) {
            double r = q0 * c[k] + q0l * c[k];
            if (k >= 1) r += 0.5 * z0 * c[k - 1];
            if (k >= 2) r += 0.25 * c[k - 2];
            c[k + 2] = r * PCFD_INV[k];
            double dterm = (k + 2.0) * c[k + 2] * hp;   /* (k+2) c_{k+2} h^{k+1} */
            double term = c[k + 2] * hp * h;            /* c_{k+2} h^{k+2} */
            der += dterm;
            val += term;
            hp *= h;
            double at = fabs(term);
            if (at > mag) mag = at;
            if (at < 1e-18 * fabs(val) && fabs(dterm) < 1e-18 * fabs(der)) {
                if (++small >= 2) break;
            } else {
                small = 0;
            }
        }
        if (k + 2 >= 80) return false;
        /* the step's own cancellation, against the local amplitude (the value
         * alone is small near a zero of an oscillating solution) */
        double amp = fmax(fabs(val), fabs(der) / sqrt(fabs(q0) + 1.0));
        if (mag > 64.0 * amp) return false;
        w = val; wp = der; z0 = zn;
        /* keep the scaled pair in range, exactly (powers of two) */
        double m = fmax(fabs(w), fabs(wp));
        if (!isfinite(m)) return false;
        if (m > 1e150 || (m < 1e-150 && m > 0.0)) {
            int e; frexp(m, &e);
            w = ldexp(w, -e); wp = ldexp(wp, -e);
            E += e;
        }
        /* below zero the solution must stay dominant: a recessive end (an
         * order near an integer) would be swamped by the stepping errors */
        double la = log2(fmax(fabs(w), fabs(wp) / sqrt(fabs(q0) + 1.0))) + E;
        if (la > lmax) lmax = la;
        if (z0 < 0.0 && la < lmax - 12.0) return false;
    }
    double r = pcfd_expscale(w, E, eh, el);
    if (!pcfd_normal(r) || (r == 0.0 && w != 0.0)) return false;
    *out = r;
    return true;
}

/* D_nu(z) for z < z_A: the expansion directly (z >= 9), or the recessive
 * solution integrated down from z_A -- through 0 when z < 0. */
static bool pcfd_pos(double nu, double z, double* out) {
    double S, Sp, eh, el;
    if (z >= 9.0 && pcfd_asym_pos(nu, z, &S, &Sp)) {
        pcfd_mquarter_sq(z, &eh, &el);
        double r = pcfd_expscale(S * pow(z, nu), 0, eh, el);
        if (!pcfd_normal(r) || (r == 0.0 && S != 0.0)) return false;
        *out = r;
        return true;
    }
    double zA = fmax(fmax(9.0, 0.6 * fabs(nu) + 3.0), 2.0 * sqrt(fabs(nu) + 1.0) + 4.0);
    if (zA <= z) zA = z * 1.25;
    for (int tries = 0; tries < 8; tries++, zA *= 1.25) {
        if (zA > 64.0) return false;
        if (!pcfd_asym_pos(nu, zA, &S, &Sp)) continue;
        double pw = pow(zA, nu);
        double wp = ((nu / zA - 0.5 * zA) * S + Sp) * pw;
        pcfd_mquarter_sq(zA, &eh, &el);
        return pcfd_integrate(nu, zA, S * pw, wp, 0, eh, el, z, out);
    }
    return false;
}

/* sum_k (a)_k / (b)_k x^k / k!, a given exactly as a double-double (a = ah + al)
 * and x = xh + xl. Plain double first; double-double when the double sum
 * cancels by more than 2^8. */
static bool pcfd_m_series(double ah, double al, double b, double xh, double xl,
                          double* out) {
    double a = ah + al, x = xh + xl;
    double t = 1.0, s = 1.0, peak = 1.0;
    int k;
    for (k = 0; k < 3000; k++) {
        t *= (a + k) * x / ((b + k) * (k + 1.0));
        s += t;
        double at = fabs(t);
        if (at > peak) peak = at;
        if (t == 0.0 || (at < 1e-17 * fabs(s) && a + k > 0.0)) break;
    }
    if (k >= 3000 || !isfinite(s)) return false;
    if (peak <= 4.0 * fabs(s)) { *out = s; return true; }
    if (peak > 1e30 * fabs(s)) return false;         /* beyond double-double */
    /* double-double re-sum, the term ratio formed in double-double too */
    sdd A; A.hi = ah; A.lo = al;
    sdd X; X.hi = xh; X.lo = xl;
    sdd T = sdd_from(1.0), Sd = sdd_from(1.0);
    for (k = 0; k < 3000; k++) {
        sdd f = sdd_mul(sdd_add_d(A, (double)k), X);
        T = sdd_div_d(sdd_mul(T, f), (b + k) * (k + 1.0));
        Sd = sdd_add(Sd, T);
        /* the peak/|sum| cancellation measured above sets the depth needed */
        if (T.hi == 0.0 || (fabs(T.hi) * peak < 1e-19 * Sd.hi * Sd.hi && a + k > 0.0)) break;
    }
    if (k >= 3000) return false;
    *out = Sd.hi + Sd.lo;
    return true;
}

/* The Maclaurin form (M) for z <= 0 or small z > 0. */
static bool pcfd_series(double nu, double z, double* out) {
    sdd x2 = sdd_two_prod(z, z);
    double xh = 0.5 * x2.hi, xl = 0.5 * x2.lo;      /* x = z^2/2 exactly */
    double g1 = pcfd_rgamma(0.5 * (1.0 - nu));
    double g2 = pcfd_rgamma(-0.5 * nu);
    double t1 = 0.0, t2 = 0.0;
    if (g1 != 0.0) {                                 /* M(-nu/2, 1/2, x) */
        double m;
        if (!pcfd_m_series(-0.5 * nu, 0.0, 0.5, xh, xl, &m)) return false;
        t1 = g1 * m;
    }
    if (g2 != 0.0 && z != 0.0) {                     /* M((1-nu)/2, 3/2, x) */
        sdd a = sdd_two_sum(0.5, -0.5 * nu);
        double m;
        if (!pcfd_m_series(a.hi, a.lo, 1.5, xh, xl, &m)) return false;
        t2 = -PCFD_SQRT2 * z * g2 * m;
    }
    double s = t1 + t2;
    if (fabs(s) * 4.0 < fabs(t1) + fabs(t2)) return false;
    double r = pcfd_expscale(PCFD_SQRTPI * s * exp2(0.5 * nu), 0, -0.5 * xh, -0.5 * xl);
    if (!pcfd_normal(r) || (r == 0.0 && s != 0.0)) return false;
    *out = r;
    return true;
}

/* z < 0 large: the connection formula, t = -z > 0. */
static bool pcfd_neg_asym(double nu, double t, double* out) {
    double w = 1.0 / (2.0 * t * t);
    double term = 1.0, s = 1.0, prev = 1.0;
    bool dec = false;
    int k;
    for (k = 1; k < 400; k++) {
        term *= (nu + 2.0 * k - 1.0) * (nu + 2.0 * k) * w / k;
        double at = fabs(term);
        if (term == 0.0) break;
        if (at > prev && dec) return false;      /* past its least term */
        if (at < prev) dec = true;
        s += term;
        prev = at;
        if (at < 1.4e-17 * fabs(s)) break;
    }
    if (k >= 400) return false;
    double qh, ql;
    pcfd_mquarter_sq(t, &qh, &ql);                   /* -t^2/4 */
    double g = pcfd_rgamma(-nu);
    /* dominant part, without its e^(t^2/4): Sqrt[2 pi] g t^(-nu-1) s */
    double second = (g != 0.0) ? 2.50662827463100050242 * g * s * pow(t, -nu - 1.0) : 0.0;
    double cp = pcfd_cospi(nu), first = 0.0;
    if (cp != 0.0) {
        /* recessive part relative to e^(t^2/4): cos(pi nu) t^nu S e^(-t^2/2);
         * skipped when it is below the dominant part's rounding */
        double lrel = (second != 0.0) ? nu * log(t) + 2.0 * qh - log(fabs(second)) : 0.0;
        if (second == 0.0 || lrel > -45.0) {
            double S, Sp;
            if (!pcfd_asym_pos(nu, t, &S, &Sp)) return false;
            first = pcfd_expscale(cp * S * pow(t, nu), 0, 2.0 * qh, 2.0 * ql);
        }
    }
    double v = second + first;
    if (v == 0.0 || fabs(v) * 4.0 < fabs(second)) return false;
    double r = pcfd_expscale(v, 0, -qh, -ql);
    if (!pcfd_normal(r) || r == 0.0) return false;
    *out = r;
    return true;
}

/* The arb last resort, accepted only as a normal double. */
static bool pcfd_arb(double nr, double ni, double zr, double zi, double* rr, double* ri) {
    double vr, vi;
    if (!flint_pcfd_d(nr, ni, zr, zi, &vr, &vi)) return false;
    if (!isfinite(vr) || !isfinite(vi)) return false;
    double big = fmax(fabs(vr), fabs(vi));
    if (big < DBL_MIN) return false;                     /* underflow: MPFR path */
    if (fabs(vr) < DBL_MIN) vr = 0.0;
    if (fabs(vi) < DBL_MIN) vi = 0.0;
    *rr = vr; *ri = vi;
    return true;
}

bool sf_machine_pcfd(double nu, double z, double* out) {
    if (!isfinite(nu) || !isfinite(z)) return false;
    if (fabs(nu) <= 40.0) {
        if (z > 1.5) {
            if (pcfd_pos(nu, z, out)) return true;
        } else if (z >= -12.0) {
            if (pcfd_series(nu, z, out)) return true;
        } else {
            if (pcfd_neg_asym(nu, -z, out)) return true;
            if (z >= -30.0 && pcfd_series(nu, z, out)) return true;
        }
    }
    double ri;
    return pcfd_arb(nu, 0.0, z, 0.0, out, &ri) && *out != 0.0;
}

/* Complex z, real nu, through the confluent kernels: (U) for Re z > 0, the
 * Kummer-transformed (M) otherwise. */
static bool pcfd_cplx_hyper(double nu, double zr, double zi, double* rr, double* ri) {
    if (fabs(nu) > PCFD_NU_MAX) return false;
    /* x = z^2 / 2 */
    double xr = 0.5 * (zr - zi) * (zr + zi), xi = zr * zi;
    double er = 0.5 * nu * M_LN2, vr, vi;
    if (zr > 0.0) {
        double ur, ui;
        if (!sf_machine_hyperu_c(-0.5 * nu, 0.5, xr, xi, &ur, &ui)) return false;
        pcfd_scale_c(ur, ui, er - 0.5 * xr, -0.5 * xi, &vr, &vi);
        if (vr == 0.0 && vi == 0.0 && (ur != 0.0 || ui != 0.0)) return false;
    } else {
        double g1 = pcfd_rgamma(0.5 * (1.0 - nu));
        double g2 = pcfd_rgamma(-0.5 * nu);
        double t1r = 0.0, t1i = 0.0, t2r = 0.0, t2i = 0.0;
        if (g1 != 0.0) {
            double mr, mi;
            if (!sf_machine_hyper1f1_c(0.5 + 0.5 * nu, 0.0, 0.5, 0.0, -xr, -xi, &mr, &mi))
                return false;
            t1r = mr * g1; t1i = mi * g1;
        }
        if (g2 != 0.0) {
            double mr, mi;
            if (!sf_machine_hyper1f1_c(1.0 + 0.5 * nu, 0.0, 1.5, 0.0, -xr, -xi, &mr, &mi))
                return false;
            /* -Sqrt[2] z M g2 */
            double cr = -PCFD_SQRT2 * g2 * zr, ci = -PCFD_SQRT2 * g2 * zi;
            t2r = cr * mr - ci * mi;
            t2i = cr * mi + ci * mr;
        }
        double sr = t1r + t2r, si = t1i + t2i;
        if (hypot(sr, si) * 16.0 < hypot(t1r, t1i) + hypot(t2r, t2i)) return false;
        pcfd_scale_c(PCFD_SQRTPI * sr, PCFD_SQRTPI * si, er + 0.5 * xr, 0.5 * xi, &vr, &vi);
        if (vr == 0.0 && vi == 0.0 && (sr != 0.0 || si != 0.0)) return false;
    }
    if (!pcfd_normal(vr) || !pcfd_normal(vi)) {
        /* A part that underflowed beside a normal partner is only a rounding
         * residue; a part that overflowed, or both underflowing, is not
         * representable. */
        if (!isfinite(vr) || !isfinite(vi)) return false;
        double big = fmax(fabs(vr), fabs(vi));
        if (big < DBL_MIN) return false;
        if (fabs(vr) < DBL_MIN) vr = 0.0 * vr;
        if (fabs(vi) < DBL_MIN) vi = 0.0 * vi;
    }
    *rr = vr; *ri = vi;
    return true;
}

bool sf_machine_pcfd_c(double nu, double zr, double zi, double* rr, double* ri) {
    if (zi == 0.0) {
        double v;
        if (!sf_machine_pcfd(nu, zr, &v)) return false;
        *rr = v; *ri = 0.0;
        return true;
    }
    if (!isfinite(nu) || !isfinite(zr) || !isfinite(zi)) return false;
    if (pcfd_cplx_hyper(nu, zr, zi, rr, ri)) return true;
    return pcfd_arb(nu, 0.0, zr, zi, rr, ri);
}

bool sf_machine_pcfd_cc(double nr, double ni, double zr, double zi, double* rr, double* ri) {
    if (ni == 0.0) return sf_machine_pcfd_c(nr, zr, zi, rr, ri);
    return pcfd_arb(nr, ni, zr, zi, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtin                                                             */
/* ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_rat(int64_t p, int64_t q) {
    return mk_fn2("Rational", expr_new_integer(p), expr_new_integer(q));
}

/* 2^(nu/2) Sqrt[Pi] / Gamma[(1 - nu)/2], the value at z = 0 (DLMF 12.2.6). */
static Expr* pcfd_at_zero(const Expr* nu) {
    Expr* p2 = mk_fn2("Power", expr_new_integer(2),
                      mk_fn2("Times", mk_rat(1, 2), expr_copy((Expr*)nu)));
    Expr* sp = mk_fn2("Power", expr_new_symbol("Pi"), mk_rat(1, 2));
    Expr* g  = mk_fn1("Gamma",
                      mk_fn2("Plus", mk_rat(1, 2),
                             mk_fn2("Times", mk_rat(-1, 2), expr_copy((Expr*)nu))));
    Expr* rg = mk_fn2("Power", g, expr_new_integer(-1));
    return eval_and_free(expr_new_function(expr_new_symbol("Times"),
                                           (Expr*[]){ p2, sp, rg }, 3));
}

Expr* builtin_parabolic_cylinder_d(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("ParabolicCylinderD", argc, 2, 2);
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    /* Machine scalar call: the shared double kernel (the NDArray / packed map
     * and the Compile VM call the same one). */
    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && ma.n == 2) {
        double rr, ri;
        if (ma.any_complex) {
            if (sf_machine_pcfd_cc(ma.re[0], ma.im[0], ma.re[1], ma.im[1], &rr, &ri))
                return sf_machine_complex_result(rr, ri);
        } else if (sf_machine_pcfd(ma.re[0], ma.re[1], &rr)) {
            return expr_new_real(rr);
        }
        /* declined (over/underflow): the extended-range arb path below */
    }

    /* Exact zero argument, any order (symbolic included). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) return pcfd_at_zero(nu);

    /* Every other inexact numeric call: arb at the arguments' precision. */
    if (flint_num_call_inexact(res))
        return flint_num_for_call(flint_num_pcfd(nu, z), res);
    return NULL;
}

void parabolic_cylinder_init(void) {
    symtab_add_builtin("ParabolicCylinderD", builtin_parabolic_cylinder_d);
    symtab_get_def("ParabolicCylinderD")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("ParabolicCylinderD",
        "ParabolicCylinderD[nu, z]\n"
        "\tgives the parabolic cylinder function D_nu(z), the solution of\n"
        "\tw'' + (nu + 1/2 - z^2/4) w = 0 that decays as z -> +Infinity.\n"
        "Entire in nu and z. D_nu(0) = 2^(nu/2) Sqrt[Pi]/Gamma[(1-nu)/2]; an exact\n"
        "integer order stays unevaluated (FunctionExpand gives the Hermite / Erfc\n"
        "form). Numeric at machine and arbitrary precision, real and complex\n"
        "arguments; threads over packed arrays and compiles.\n"
        "D[ParabolicCylinderD[nu,z],z] = z/2 ParabolicCylinderD[nu,z] - ParabolicCylinderD[nu+1,z].");
}
