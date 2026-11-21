/* Machine-precision kernels for the zeta family -- see sf_machine_zeta.h.
 *
 * Algorithms (all in raw double / double re-im; no <complex.h>):
 *
 *  - Riemann and Hurwitz zeta: Euler-Maclaurin summation (DLMF 25.11.5 /
 *    25.2.9).  The direct sum runs to w = a + N with w >= W(s), W chosen so the
 *    Bernoulli correction series converges geometrically (its term ratio is
 *    about |s + 2k|^2 / (2 pi w)^2); the correction is truncated adaptively at
 *    2^-56 of the running total.  Integer real s returns tabulated constants
 *    (zeta(2..64), zeta(0..-60)); s >= 20 sums the Dirichlet series directly.
 *    Real s < 1/2 uses the functional equation (DLMF 25.4.1)
 *        zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s),
 *    with sin(pi s/2) reduced exactly in s (so the trivial zeros keep their
 *    relative accuracy), (2 pi)^s formed from an exact two-product, and the
 *    rounding of 1-s corrected to first order in Gamma and zeta. Complex
 *    Re s < 1/2 declines (see sf_zeta_c).
 *
 *  - Dilogarithm Li_2: the Bernoulli series in u = -log(1-z)
 *        Li_2(z) = Sum_{n>=0} B_n u^(n+1)/(n+1)!        (|u| < 2 pi)
 *    ('t Hooft & Veltman, Nucl. Phys. B153 (1979) 365; Crandall, "Note on fast
 *    polylogarithm computation", 2006), after the maps
 *        z -> 1/z   Li_2(z) = -Li_2(1/z) - pi^2/6 - log^2(-z)/2      (|z| > 1)
 *        z -> 1-z   Li_2(z) = pi^2/6 - log z log(1-z) - Li_2(1-z)   (Re z > 1/2)
 *    which land every z in |z| <= 1, Re z <= 1/2, where |u| <= 1.25 and the
 *    series needs at most ~14 terms.
 *
 *  - PolyLog Li_s, real s: the defining power series for |z| <= 1/2; the
 *    Jonquiere expansion in mu = log z (DLMF 25.12.12 for integer s, 25.12.11
 *    otherwise) for |mu| < 3.6; for integer s and larger |z| the inversion
 *    formula (DLMF 25.12.13 / Crandall eq. 1.4)
 *        Li_n(z) = (-1)^(n-1) Li_n(1/z) - (2 pi i)^n/n! B_n(1/2 + log(-z)/(2 pi i)).
 *
 *  - LerchPhi with complex z: the defining series Sum z^k (k+a)^-s, |z| < 1.
 *
 * Branch cut convention (matches the MPFR path and Mathematica): PolyLog's cut
 * is [1, Infinity), and a real z on it takes the value from below the cut.
 */
#include "sf_machine_zeta.h"
#include "sf_machine.h"   /* sf_machine_lerchphi, sf_machine_polylog */

#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* pi^2/6 split into a double and its rounding residual, so pi^2/6 - (...)
 * stays accurate to ~0.5 ulp where the subtraction is the whole answer. */
#define ZM_PI2_6_HI 1.6449340668482264
#define ZM_PI2_6_LO 3.040672350398476e-17
#define ZM_LN2PI_HI 1.8378770664093456      /* log(2 pi) */
#define ZM_LN2PI_LO (-7.756588316134483e-17)
#define ZM_EPS      1.3877787807814457e-17  /* 2^-56: correction cut-off */

/* B_2k / (2k)!, k = 1..25 (Euler-Maclaurin correction coefficients). */
static const double ZM_B2K_FACT[25] = {
    8.33333333333333333333e-2, -1.38888888888888888889e-3, 3.30687830687830687831e-5,
    -8.26719576719576719577e-7, 2.08767569878680989792e-8, -5.28419013868749318485e-10,
    1.33825365306846788328e-11, -3.38968029632258286683e-13, 8.58606205627784456414e-15,
    -2.17486869855806187304e-16, 5.5090028283602295152e-18, -1.39544646858125233407e-19,
    3.53470703962946747169e-21, -8.9535174270375468504e-23, 2.26795245233768306031e-24,
    -5.74479066887220244526e-26, 1.45517247561486490187e-27, -3.68599494066531017818e-29,
    9.33673425709504467203e-31, -2.36502241570062993456e-32, 5.99067176248213430466e-34,
    -1.51745488446829026171e-35, 3.84375812545418823223e-37, -9.73635307264669103527e-39,
    2.46624704420068095711e-40
};

/* B_2k / (2k+1)!, k = 1..16: the dilogarithm's Bernoulli series. */
static const double ZM_DILOG_C[16] = {
    2.77777777777777777778e-2, -2.77777777777777777778e-4, 4.72411186696900982615e-6,
    -9.18577307466196355085e-8, 1.8978869988970999072e-9, -4.06476164514422552681e-11,
    8.92169102045645255522e-13, -1.99392958607210756872e-14, 4.51898002961991819165e-16,
    -1.03565176121812470145e-17, 2.39521862102618674574e-19, -5.58178587432500933628e-21,
    1.30915075541832128581e-22, -3.08741980242674029324e-24, 7.31597565270220342036e-26,
    -1.74084565723400074099e-27
};

/* zeta(n), n = 2..64. */
static const double ZM_ZPOS[63] = {
    1.64493406684822643647, 1.2020569031595942854, 1.08232323371113819152,
    1.03692775514336992633, 1.01734306198444913971, 1.00834927738192282684,
    1.00407735619794433938, 1.00200839282608221442, 1.00099457512781808534,
    1.00049418860411946456, 1.0002460865533080483, 1.00012271334757848915,
    1.00006124813505870483, 1.00003058823630702049, 1.00001528225940865187,
    1.00000763719763789976, 1.00000381729326499984, 1.00000190821271655394,
    1.0000009539620338728, 1.00000047693298678781, 1.00000023845050272773,
    1.00000011921992596531, 1.00000005960818905126, 1.00000002980350351465,
    1.00000001490155482837, 1.00000000745071178984, 1.00000000372533402479,
    1.00000000186265972351, 1.00000000093132743242, 1.0000000004656629065,
    1.00000000023283118337, 1.00000000011641550173, 1.00000000005820772088,
    1.00000000002910385044, 1.00000000001455192189, 1.00000000000727595984,
    1.00000000000363797955, 1.00000000000181898965, 1.00000000000090949478,
    1.00000000000045474738, 1.00000000000022737368, 1.00000000000011368684,
    1.00000000000005684342, 1.00000000000002842171, 1.00000000000001421085,
    1.00000000000000710543, 1.00000000000000355271, 1.00000000000000177636,
    1.00000000000000088818, 1.00000000000000044409, 1.00000000000000022204,
    1.00000000000000011102, 1.00000000000000005551, 1.00000000000000002776,
    1.00000000000000001388, 1.00000000000000000694, 1.00000000000000000347,
    1.00000000000000000173, 1.00000000000000000087, 1.00000000000000000043,
    1.00000000000000000022, 1.00000000000000000011, 1.00000000000000000005
};

/* zeta(-n), n = 0..60 (= -B_{n+1}/(n+1); zero at the negative even integers). */
static const double ZM_ZNEG[61] = {
    -5.0e-1, -8.33333333333333333333e-2, 0.0, 8.33333333333333333333e-3, 0.0,
    -3.96825396825396825397e-3, 0.0, 4.16666666666666666667e-3, 0.0,
    -7.57575757575757575758e-3, 0.0, 2.10927960927960927961e-2, 0.0,
    -8.33333333333333333333e-2, 0.0, 4.43259803921568627451e-1, 0.0,
    -3.0539543302701197438, 0.0, 2.64562121212121212121e+1, 0.0,
    -2.81460144927536231884e+2, 0.0, 3.60751054639804639805e+3, 0.0,
    -5.48275833333333333333e+4, 0.0, 9.74936823850574712644e+5, 0.0,
    -2.00526957966880789461e+7, 0.0, 4.72384867721629901961e+8, 0.0,
    -1.26357247959166666667e+10, 0.0, 3.80879311252453688116e+11, 0.0,
    -1.28508504993050833333e+13, 0.0, 4.82414483548501703716e+14, 0.0,
    -2.00403106565162527381e+16, 0.0, 9.1677436031953307757e+17, 0.0,
    -4.59798883436565034904e+19, 0.0, 2.51804719214510956971e+21, 0.0,
    -1.50017334921539287337e+23, 0.0, 9.68995788746359406565e+24, 0.0,
    -6.76458823792928209909e+26, 0.0, 5.08906594686622896898e+28, 0.0,
    -4.11472887925579786977e+30, 0.0, 3.56665820953755561097e+32, 0.0
};

/* Bernoulli numbers B_0..B_40 (Bernoulli polynomials for PolyLog inversion). */
static const double ZM_BN[41] = {
    1.0, -5.0e-1, 1.66666666666666666667e-1, 0.0, -3.33333333333333333333e-2, 0.0,
    2.38095238095238095238e-2, 0.0, -3.33333333333333333333e-2, 0.0,
    7.57575757575757575758e-2, 0.0, -2.53113553113553113553e-1, 0.0,
    1.16666666666666666667, 0.0, -7.09215686274509803922, 0.0,
    5.49711779448621553885e+1, 0.0, -5.29124242424242424242e+2, 0.0,
    6.19212318840579710145e+3, 0.0, -8.65802531135531135531e+4, 0.0,
    1.42551716666666666667e+6, 0.0, -2.7298231067816091954e+7, 0.0,
    6.01580873900642368384e+8, 0.0, -1.51163157670921568627e+10, 0.0,
    4.29614643061166666667e+11, 0.0, -1.37116552050883327722e+13, 0.0,
    4.88332318973593166667e+14, 0.0, -1.92965793419400681486e+16
};


/* ------------------------------------------------------------------------ */
/* Small complex helpers (re, im pairs)                                      */
/* ------------------------------------------------------------------------ */

static void zm_mul(double ar, double ai, double br, double bi, double* r, double* i) {
    double t = ar * br - ai * bi;
    *i = ar * bi + ai * br;
    *r = t;
}

/* a / b by Smith's algorithm (no spurious overflow, no __divdc3). */
static void zm_div(double ar, double ai, double br, double bi, double* r, double* i) {
    if (fabs(br) >= fabs(bi)) {
        double q = bi / br, d = br + bi * q;
        double rr = (ar + ai * q) / d;
        *i = (ai - ar * q) / d;
        *r = rr;
    } else {
        double q = br / bi, d = bi + br * q;
        double rr = (ar * q + ai) / d;
        *i = (ai * q - ar) / d;
        *r = rr;
    }
}

static void zm_log(double xr, double xi, double* r, double* i) {
    *r = log(hypot(xr, xi));
    *i = atan2(xi, xr);
}

static void zm_exp(double xr, double xi, double* r, double* i) {
    double m = exp(xr);
    *r = m * cos(xi);
    *i = m * sin(xi);
}

/* w^-s for complex w, complex s (principal branch). A positive real w takes
 * the cheaper real-log route. */
static void zm_pow_neg(double wr, double wi, double sr, double si, double* r, double* i) {
    if (wi == 0.0 && wr > 0.0) {
        double L = log(wr);
        double m = exp(-sr * L), ph = -si * L;
        *r = m * cos(ph);
        *i = m * sin(ph);
        return;
    }
    double Lr, Li;
    zm_log(wr, wi, &Lr, &Li);
    zm_exp(-(sr * Lr - si * Li), -(sr * Li + si * Lr), r, i);
}

/* sin(pi r) and cos(pi r) with exact argument reduction in r (so that the
 * zeros at the integers are relatively accurate). */
static void zm_sincospi(double x, double* s, double* c) {
    double r = x - 2.0 * floor(0.5 * x + 0.5);     /* exact, r in [-1, 1) */
    double sg = 1.0;
    if (r < 0.0) { r = -r; sg = -1.0; }              /* sin odd, cos even */
    double cs = 1.0;
    if (r > 0.5) { r = 1.0 - r; cs = -1.0; }         /* sin(pi-x)=sin x, cos(pi-x)=-cos x */
    double t = M_PI * r;
    if (r == 0.5) { *s = sg; *c = 0.0; return; }
    *s = sg * sin(t);
    *c = cs * cos(t);
}

/* k + a rounded to w, with the exact residual: w + e == k + a (Knuth TwoSum).
 * A term w^-s then picks up the first-order correction (1 - s e / w), which
 * removes the s * (rounding of k + a) error a non-integer a otherwise puts
 * into every Hurwitz / Lerch term. */
static double zm_two_sum(double k, double a, double* e) {
    double w = k + a;
    double bb = w - k;
    *e = (k - (w - bb)) + (a - bb);
    return w;
}

/* (h + l) *= z for a double-double running power h + l (complex, parts
 * (hr, hi) + (lr, li)) and an exact double z: every product's rounding error
 * is recovered with an FMA, so z^k stays accurate to ~2^-100 relative instead
 * of drifting by k/2 ulp. A real z (zi == 0) takes the two-FMA path. */
static void zm_dd_mul(double* hr, double* hi, double* lr, double* li,
                      double zr, double zi) {
    if (zi == 0.0 && *hi == 0.0 && *li == 0.0) {
        double p = *hr * zr;
        double l = fma(*hr, zr, -p) + *lr * zr;
        double q = p + l;
        *lr = l - (q - p);
        *hr = q;
        return;
    }
    double a = *hr, b = *hi;
    double p1 = a * zr, p2 = b * zi, p3 = a * zi, p4 = b * zr;
    double e1 = fma(a, zr, -p1), e2 = fma(b, zi, -p2);
    double e3 = fma(a, zi, -p3), e4 = fma(b, zr, -p4);
    double nr = p1 - p2, ni = p3 + p4;
    double b1 = nr - p1, f1 = (p1 - (nr - b1)) + (-p2 - b1);
    double b2 = ni - p3, f2 = (p3 - (ni - b2)) + (p4 - b2);
    double nlr = (e1 - e2) + f1 + (*lr * zr - *li * zi);
    double nli = (e3 + e4) + f2 + (*lr * zi + *li * zr);
    double qr = nr + nlr, qi = ni + nli;
    *lr = nlr - (qr - nr); *li = nli - (qi - ni);
    *hr = qr; *hi = qi;
}

/* Cancellation guard: a sum whose terms are much larger than its value has
 * lost log2(ratio) bits to rounding before anything else goes wrong; past
 * this ratio the kernel declines and the arbitrary-precision path answers. */
#define ZM_CANCEL_MAX 16.0
/* PolyLog's Jonquiere expansion near an integer order, and its power series
 * for s < 0, cancel with no structure that bounds the error by the ratio
 * alone (the pieces carry the error of zeta and Gamma themselves), so they
 * decline earlier. */
#define ZM_CANCEL_PL 4.0

/* exp(s * log(2 pi)) for real s, with the product formed exactly (two-product)
 * so that a large |s| does not amplify the rounding of s * log(2 pi). */
static double zm_two_pi_pow(double s) {
    double p = s * ZM_LN2PI_HI;
    double e = fma(s, ZM_LN2PI_HI, -p) + s * ZM_LN2PI_LO;
    return exp(p) * (1.0 + e);
}

/* ------------------------------------------------------------------------ */
/* Euler-Maclaurin (Hurwitz) cores                                           */
/* ------------------------------------------------------------------------ */

/* Smallest w = a + N the EM correction series is run at. */
static double zm_em_w(double abs_s) {
    double W = 0.5 * abs_s + 4.0;
    return W < 10.0 ? 10.0 : W;
}

/* Real Hurwitz zeta, real s != 1, real a > 0. */
static bool zm_hurwitz_em_r(double s, double a, double* out) {
    double W = zm_em_w(fabs(s));
    int N = 0;
    if (a < W) {
        double n = ceil(W - a);
        if (n > 1.0e5) return false;
        N = (int)n;
    }
    double sum = 0.0;
    for (int k = N - 1; k >= 0; k--) {                  /* small terms first */
        double e, wk = zm_two_sum((double)k, a, &e);
        double f = pow(wk, -s);
        sum += e != 0.0 ? f * (1.0 - s * e / wk) : f;
    }
    double w = a + (double)N;
    double wps = pow(w, -s);
    double tail = 0.5 * wps;
    double head = w * wps / (s - 1.0);
    double tot = sum + head + tail;
    double t = s * wps / w;                         /* (s)_1 w^(-s-1) */
    double w2 = 1.0 / (w * w);
    double corr = 0.0;
    bool conv = false;
    for (int k = 1; k <= 25; k++) {
        double term = ZM_B2K_FACT[k - 1] * t;
        corr += term;
        if (fabs(term) <= ZM_EPS * fabs(tot)) { conv = true; break; }
        t *= (s + 2.0 * k - 1.0) * (s + 2.0 * k) * w2;
    }
    if (!conv) return false;
    double v = sum + (head + (tail + corr));
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* Complex Hurwitz zeta, complex s != 1, complex a with Re a > 0. */
static bool zm_hurwitz_em_c(double sr, double si, double ar, double ai,
                            double* vr, double* vi) {
    double as = hypot(sr, si);
    if (!(as < 1.0e4)) return false;
    double W = zm_em_w(as);
    int N = 0;
    if (ar < W) {
        double n = ceil(W - ar);
        if (n > 1.0e5) return false;
        N = (int)n;
    }
    double Sr = 0.0, Si = 0.0, A = 0.0;              /* A = Sum |term| */
    for (int k = N - 1; k >= 0; k--) {
        double tr, ti, e;
        double wk = zm_two_sum((double)k, ar, &e);
        zm_pow_neg(wk, ai, sr, si, &tr, &ti);
        if (e != 0.0) {
            /* term *= 1 - s e / w, w = wk + i ai */
            double qr, qi;
            zm_div(e * sr, e * si, wk, ai, &qr, &qi);
            zm_mul(tr, ti, 1.0 - qr, -qi, &tr, &ti);
        }
        Sr += tr; Si += ti;
        A += hypot(tr, ti);
    }
    double wr = ar + (double)N, wi = ai;
    double pr, pi;                                   /* w^-s */
    zm_pow_neg(wr, wi, sr, si, &pr, &pi);
    /* head = w * w^-s / (s - 1) */
    double hr, hi;
    zm_mul(wr, wi, pr, pi, &hr, &hi);
    zm_div(hr, hi, sr - 1.0, si, &hr, &hi);
    double tr = 0.5 * pr, ti = 0.5 * pi;
    double totm = hypot(Sr + hr + tr, Si + hi + ti);
    /* t = s * w^-s / w */
    double ur, ui;
    zm_mul(sr, si, pr, pi, &ur, &ui);
    zm_div(ur, ui, wr, wi, &ur, &ui);
    double w2r, w2i;                                 /* 1 / w^2 */
    zm_mul(wr, wi, wr, wi, &w2r, &w2i);
    zm_div(1.0, 0.0, w2r, w2i, &w2r, &w2i);
    double cr = 0.0, ci = 0.0;
    bool conv = false;
    for (int k = 1; k <= 25; k++) {
        double b = ZM_B2K_FACT[k - 1];
        double er = b * ur, ei = b * ui;
        cr += er; ci += ei;
        if (hypot(er, ei) <= ZM_EPS * totm) { conv = true; break; }
        /* t *= (s + 2k - 1)(s + 2k) / w^2 */
        double fr, fi;
        zm_mul(sr + 2.0 * k - 1.0, si, sr + 2.0 * k, si, &fr, &fi);
        zm_mul(fr, fi, w2r, w2i, &fr, &fi);
        zm_mul(ur, ui, fr, fi, &ur, &ui);
    }
    if (!conv) return false;
    double rr = Sr + (hr + (tr + cr));
    double ri = Si + (hi + (ti + ci));
    if (!isfinite(rr) || !isfinite(ri)) return false;
    /* A complex s / a can make the sum exponentially small against its terms
     * (large Im s with Im a != 0): decline rather than return the noise. */
    A += hypot(hr, hi) + hypot(tr, ti);
    if (A > ZM_CANCEL_MAX * hypot(rr, ri)) return false;
    *vr = rr; *vi = ri;
    return true;
}

/* ------------------------------------------------------------------------ */
/* a <= 0: the symmetric convention, shifted into (0, 1]                     */
/* ------------------------------------------------------------------------ */

/* Largest shift m the head sum below is run for; past it the
 * arbitrary-precision path answers. */
#define ZM_SHIFT_MAX 1.0e7

/* For a real a <= 0 the sums are taken in the symmetric form ((k+a)^2)^(-s/2),
 * i.e. |k+a|^-s, so with m = floor(-a) + 1 (a + m in (0, 1])
 *   Sum_k z^k |k+a|^-s = H + z^m T,   H = Sum_{k<m} z^k |k+a|^-s,
 * and T the principal-branch sum at a + m. This forms H and z^m. The k + a = 0
 * term of an integer a is 0^-s: dropped for Re s > 0, z^k at s = 0, and 0 for
 * Re s < 0 (declined at Re s = 0 != s, which has no value). |k+a| is formed by
 * TwoSum and its power corrected to first order; z^k is a double-double power;
 * the sum is Neumaier-compensated. */
static bool zm_shift_head(double zr, double zi, double sr, double si, double a,
                          double m, double* Hr, double* Hi, double* Pr, double* Pi) {
    bool real = (zi == 0.0 && si == 0.0);
    double hr = 1.0, hi = 0.0, lr = 0.0, li = 0.0;     /* z^k = h + l */
    double Sr = 0.0, Si = 0.0, cr = 0.0, ci = 0.0;
    for (double k = 0.0; k < m; k += 1.0) {
        double e, wk = zm_two_sum(k, a, &e);
        double tr = 0.0, ti = 0.0;
        if (wk == 0.0 && e == 0.0) {                   /* the 0^-s term */
            if (sr == 0.0 && si != 0.0) return false;
            if (sr == 0.0) { tr = hr + lr; ti = hi + li; }   /* s = 0: 0^0 = 1 */
        } else {
            double w = -wk, ew = -e;                   /* |k+a| = w + ew, w > 0 */
            double fr, fi;
            if (si == 0.0) {
                fr = pow(w, -sr); fi = 0.0;
                if (ew != 0.0) fr *= 1.0 - sr * ew / w;
            } else {
                zm_pow_neg(w, 0.0, sr, si, &fr, &fi);
                if (ew != 0.0) zm_mul(fr, fi, 1.0 - sr * ew / w, -si * ew / w, &fr, &fi);
            }
            if (real) { tr = (hr + lr) * fr; ti = 0.0; }
            else zm_mul(hr + lr, hi + li, fr, fi, &tr, &ti);
        }
        double nr = Sr + tr;
        cr += fabs(Sr) >= fabs(tr) ? (Sr - nr) + tr : (tr - nr) + Sr;
        Sr = nr;
        double ni = Si + ti;
        ci += fabs(Si) >= fabs(ti) ? (Si - ni) + ti : (ti - ni) + Si;
        Si = ni;
        if (zr != 1.0 || zi != 0.0) zm_dd_mul(&hr, &hi, &lr, &li, zr, zi);
    }
    *Hr = Sr + cr; *Hi = Si + ci;
    *Pr = hr + lr; *Pi = hi + li;
    return isfinite(*Hr) && isfinite(*Hi);
}

/* m = floor(-a) + 1 for a <= 0, or false past ZM_SHIFT_MAX. */
static bool zm_shift_count(double a, double* m) {
    double mf = floor(-a) + 1.0;
    if (!(mf <= ZM_SHIFT_MAX)) return false;
    *m = mf;
    return true;
}

/* ------------------------------------------------------------------------ */
/* Riemann zeta                                                              */
/* ------------------------------------------------------------------------ */

bool sf_zeta_r(double s, double* out) {
    if (!isfinite(s)) return false;
    if (s == 1.0) return false;                         /* pole */
    if (s == floor(s)) {                                /* tabulated integers */
        if (s >= 2.0 && s <= 64.0) { *out = ZM_ZPOS[(int)s - 2]; return true; }
        if (s > 64.0) { *out = 1.0; return true; }      /* 1 + 2^-65 rounds to 1 */
        if (s <= 0.0 && s >= -60.0) { *out = ZM_ZNEG[(int)(-s)]; return true; }
        if (s < 0.0 && fmod(s, 2.0) == 0.0) { *out = 0.0; return true; }
    }
    if (s >= 20.0) {
        /* The Dirichlet series converges at once: n^-s < 2^-60 for n >= 9. */
        double sum = 0.0;
        for (int n = 8; n >= 2; n--) sum += pow((double)n, -s);
        *out = 1.0 + sum;
        return true;
    }
    if (s >= 0.5) return zm_hurwitz_em_r(s, 1.0, out);
    /* Functional equation: zeta(s) = (2 pi)^s / pi sin(pi s/2) Gamma(1-s) zeta(1-s).
     * 1 - s is generally not a double: it rounds to t with residual d (an exact
     * TwoSum), and Gamma (condition ~ t psi(t), ~110 at t = 32) and zeta (pole
     * at t = 1) would turn that rounding into tens of ulps. Both factors are
     * therefore corrected to first order in d, with log-derivatives that need
     * only a few correct digits:
     *   Gamma(t + d) = Gamma(t) (1 + d psi(t)),
     *   zeta(t + d)  = zeta(t)  (1 + d zeta'/zeta(t)). */
    double t = 1.0 - s;
    double bb = t - 1.0;
    double d = (1.0 - (t - bb)) + (-s - bb);
    double z1;
    if (!sf_zeta_r(t, &z1)) return false;
    double sn, cs;
    zm_sincospi(0.5 * s, &sn, &cs);
    double g = tgamma(t);
    if (d != 0.0) {
        double tt = t < 1.0 ? t + 1.0 : t;                /* psi(t), asymptotic */
        double psi = log(tt) - 0.5 / tt - 1.0 / (12.0 * tt * tt);
        if (t < 1.0) psi -= 1.0 / t;
        double dlz;                                        /* zeta'/zeta(t) */
        if (t < 2.0) dlz = -1.0 / ((t - 1.0) * (1.0 + 0.5772156649015329 * (t - 1.0)));
        else dlz = -(0.6931471805599453 * (pow(2.0, -t) + 2.0 * pow(4.0, -t))
                     + 1.0986122886681098 * pow(3.0, -t) + 1.6094379124341003 * pow(5.0, -t));
        g *= 1.0 + d * psi;
        z1 *= 1.0 + d * dlz;
    }
    double v = zm_two_pi_pow(s) / M_PI * sn * g * z1;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool sf_zeta_c(double sr, double si, double* vr, double* vi) {
    if (!isfinite(sr) || !isfinite(si)) return false;
    if (si == 0.0) {
        double v;
        if (!sf_zeta_r(sr, &v)) return false;
        *vr = v; *vi = 0.0;
        return true;
    }
    /* Re s < 1/2 would need the functional equation, whose Gamma(1-s) factor
     * cannot be formed in double to better than ~|1-s| log|1-s| ulps (the
     * exponent of the Stirling form is large, and every ulp of it is a
     * relative ulp of the result): measured 18 ulp median, 70 max on
     * Re s in [-10, 1/2]. The arbitrary-precision path answers instead. */
    if (!(sr >= 0.5)) return false;
    if (sr >= 20.0) {
        /* The Dirichlet series, as sf_zeta_r does: n^-s has modulus n^-Re s,
         * below 2^-60 of the leading 1 for n >= 9.  Without this every
         * |s| >= 10^4 declined to the arbitrary-precision path -- 1.7 s for
         * Zeta[1.*^6 + 5. I], and a Newton iterate at Re s ~ 10^189 did not
         * return.  The phase Im s Log[n] is a rounded double, so a term is
         * used only while its phase error (|Im s Log n| ulps) times its modulus
         * stays far below an ulp of the result (~1); a term whose modulus has
         * underflowed is exactly 0 whatever its phase.  Otherwise decline. */
        double Sr = 1.0, Si = 0.0;
        for (int n = 8; n >= 2; n--) {
            double mag = pow((double)n, -sr);
            if (mag == 0.0) continue;
            double ph = si * log((double)n);
            if (!(mag * fabs(ph) * DBL_EPSILON <= 0x1p-58)) return false;
            Sr += mag * cos(ph);
            Si -= mag * sin(ph);
        }
        *vr = Sr; *vi = Si;
        return true;
    }
    return zm_hurwitz_em_c(sr, si, 1.0, 0.0, vr, vi);
}

/* ------------------------------------------------------------------------ */
/* Hurwitz zeta                                                              */
/* ------------------------------------------------------------------------ */

bool sf_hurwitz_r(double s, double a, double* out) {
    if (!isfinite(s) || !isfinite(a)) return false;
    if (s == 1.0 || !(a > 0.0)) return false;
    if (a == 1.0) return sf_zeta_r(s, out);            /* one value, one kernel */
    if (s < 0.5) return false;                         /* continuation: MPFR path */
    return zm_hurwitz_em_r(s, a, out);
}

bool sf_hurwitz_c(double sr, double si, double ar, double ai,
                  double* vr, double* vi) {
    if (!isfinite(sr) || !isfinite(si) || !isfinite(ar) || !isfinite(ai)) return false;
    if (!(ar > 0.0)) return false;
    if (si == 0.0 && ai == 0.0) {
        double v;
        if (!sf_hurwitz_r(sr, ar, &v)) return false;
        *vr = v; *vi = 0.0;
        return true;
    }
    if (ai == 0.0 && ar == 1.0) return sf_zeta_c(sr, si, vr, vi);
    if (sr < 0.5) return false;                        /* continuation: MPFR path */
    if (sr == 1.0 && si == 0.0) return false;
    return zm_hurwitz_em_c(sr, si, ar, ai, vr, vi);
}

/* Zeta[s, a]: the Hurwitz sum in the symmetric convention, which is the
 * principal one for a > 0 and, for a real a <= 0, the shifted form of
 * zm_shift_head. Shared by the Zeta builtin, its NDArray / packed / Compile[]
 * kernel and LerchPhi at z = 1. HurwitzZeta keeps the principal powers and so
 * calls sf_hurwitz_* directly. */
bool sf_zeta2_r(double s, double a, double* out) {
    if (!isfinite(s) || !isfinite(a) || s == 1.0) return false;
    if (a > 0.0) return sf_hurwitz_r(s, a, out);
    double m, H, Hi, P, Pi, T;
    if (!zm_shift_count(a, &m)) return false;
    double r = a + m;                                  /* (0, 1], exact for m >= 2 */
    if (!(r == 1.0 ? sf_zeta_r(s, &T) : sf_hurwitz_r(s, r, &T))) return false;
    if (!zm_shift_head(1.0, 0.0, s, 0.0, a, m, &H, &Hi, &P, &Pi)) return false;
    double v = H + T;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool sf_zeta2_c(double sr, double si, double ar, double ai, double* vr, double* vi) {
    if (!isfinite(sr) || !isfinite(si) || !isfinite(ar) || !isfinite(ai)) return false;
    if (ai != 0.0 || ar > 0.0) return sf_hurwitz_c(sr, si, ar, ai, vr, vi);
    if (si == 0.0) {                                   /* real s, real a <= 0 */
        double v;
        if (!sf_zeta2_r(sr, ar, &v)) return false;
        *vr = v; *vi = 0.0;
        return true;
    }
    double m, Hr, Hi, Pr, Pi, Tr, Ti;                  /* complex s, real a <= 0 */
    if (!zm_shift_count(ar, &m)) return false;
    if (!sf_hurwitz_c(sr, si, ar + m, 0.0, &Tr, &Ti)) return false;
    if (!zm_shift_head(1.0, 0.0, sr, si, ar, m, &Hr, &Hi, &Pr, &Pi)) return false;
    *vr = Hr + Tr; *vi = Hi + Ti;
    return isfinite(*vr) && isfinite(*vi);
}

/* ------------------------------------------------------------------------ */
/* LerchPhi, complex arguments                                               */
/* ------------------------------------------------------------------------ */

/* The defining series Sum_k z^k (k+a)^-s, |z| < 1, for real or complex z and
 * (s, a) real (a > 0) or complex (Re a > 0). Accuracy measures:
 *   - z^k is carried in double-double (an exact FMA residual per product), so
 *     the k-fold rounding of a running power does not grow into the terms;
 *   - k + a is formed by TwoSum and each (k+a)^-s corrected to first order;
 *   - the sum is compensated (Neumaier) and guarded against cancellation.
 * Terminates on two consecutive terms below 2^-58 of the sum once the terms
 * are guaranteed to decrease: |z| ((k+1+a)/(k+a))^(-Re s) < 1. */
static bool zm_lerch_series(double zr, double zi, double sr, double si,
                            double ar, double ai, double* vr, double* vi) {
    double zm2 = zr * zr + zi * zi;
    double logz = 0.5 * log(zm2);                      /* log |z| < 0 */
    bool real_z = (zi == 0.0);
    bool real_sa = (si == 0.0 && ai == 0.0);
    double hr = 1.0, hi = 0.0, lr = 0.0, li = 0.0;     /* z^k = h + l */
    double Sr = 0.0, Si = 0.0, cr = 0.0, ci = 0.0;     /* Neumaier sums */
    double A = 0.0;
    int small = 0;
    for (int k = 0; k < 200000; k++) {
        double e, wk = zm_two_sum((double)k, ar, &e);
        double tr, ti;
        if (real_sa) {
            double f = pow(wk, -sr);
            if (e != 0.0) f *= 1.0 - sr * e / wk;
            tr = (hr + lr) * f;
            ti = real_z ? 0.0 : (hi + li) * f;
        } else {
            double fr, fi;
            zm_pow_neg(wk, ai, sr, si, &fr, &fi);
            if (e != 0.0) {
                double qr, qi;
                zm_div(e * sr, e * si, wk, ai, &qr, &qi);
                zm_mul(fr, fi, 1.0 - qr, -qi, &fr, &fi);
            }
            zm_mul(hr + lr, hi + li, fr, fi, &tr, &ti);
        }
        /* Neumaier: S += t with the lost low part kept in c */
        double nr = Sr + tr;
        cr += fabs(Sr) >= fabs(tr) ? (Sr - nr) + tr : (tr - nr) + Sr;
        Sr = nr;
        if (!real_z || !real_sa) {
            double ni = Si + ti;
            ci += fabs(Si) >= fabs(ti) ? (Si - ni) + ti : (ti - ni) + Si;
            Si = ni;
        }
        double at = fabs(tr) + fabs(ti);
        A += at;
        if (at <= 0.25 * ZM_EPS * (fabs(Sr) + fabs(Si)) &&
            (sr >= 0.0 || logz - sr * log1p(1.0 / wk) < 0.0)) {
            if (++small >= 2) break;
        } else {
            small = 0;
        }
        zm_dd_mul(&hr, &hi, &lr, &li, zr, zi);         /* z^(k+1) */
        if (hr == 0.0 && hi == 0.0) { small = 2; break; }
    }
    if (small < 2) return false;
    Sr += cr; Si += ci;
    if (!isfinite(Sr) || !isfinite(Si)) return false;
    if (A > ZM_CANCEL_MAX * (fabs(Sr) + fabs(Si))) return false;
    *vr = Sr; *vi = Si;
    return true;
}

bool sf_lerch_r(double z, double s, double a, double* out) {
    if (!isfinite(z) || !isfinite(s) || !isfinite(a)) return false;
    if (z == 1.0) return sf_zeta2_r(s, a, out);
    if (!(fabs(z) < 1.0)) return false;
    if (!(a > 0.0)) {                 /* a <= 0: the ((k+a)^2)^(s/2) convention */
        double m, H, Hi, P, Pi, T;
        if (!zm_shift_count(a, &m) || m > 1.0e5) return false;
        if (!zm_shift_head(z, 0.0, s, 0.0, a, m, &H, &Hi, &P, &Pi)) return false;
        if (!sf_lerch_r(z, s, a + m, &T)) return false;
        double v = H + P * T;
        if (!isfinite(v)) return false;
        *out = v;
        return true;
    }
    if (z == 0.0) { *out = pow(a, -s); return isfinite(*out); }
    double vi;
    return zm_lerch_series(z, 0.0, s, 0.0, a, 0.0, out, &vi);
}

bool sf_lerch_c(double zr, double zi, double sr, double si, double ar, double ai,
                double* vr, double* vi) {
    if (!isfinite(zr) || !isfinite(zi) || !isfinite(sr) || !isfinite(si) ||
        !isfinite(ar) || !isfinite(ai)) return false;
    if (zi == 0.0 && si == 0.0 && ai == 0.0) {
        double v;
        if (!sf_machine_lerchphi(zr, sr, ar, &v)) return false;
        *vr = v; *vi = 0.0;
        return true;
    }
    /* z = 1 is Zeta[s, a], symmetric convention included. Otherwise Re a > 0:
     * then (k+a)^s is the principal power for every k, which is what both
     * conventions (the symmetric ((k+a)^2)^(s/2) and the principal one) agree
     * on. Anything else is left to the MPFR path. */
    if (zi == 0.0 && zr == 1.0) return sf_zeta2_c(sr, si, ar, ai, vr, vi);
    if (!(ar > 0.0)) return false;
    double zm2 = zr * zr + zi * zi;
    if (!(zm2 < 1.0)) return false;
    if (zm2 == 0.0) {                                  /* a^-s */
        zm_pow_neg(ar, ai, sr, si, vr, vi);
        return isfinite(*vr) && isfinite(*vi);
    }
    return zm_lerch_series(zr, zi, sr, si, ar, ai, vr, vi);
}

/* ------------------------------------------------------------------------ */
/* Dilogarithm                                                               */
/* ------------------------------------------------------------------------ */

/* Li_2 = u - u^2/4 + Sum_{k>=1} B_2k u^(2k+1)/(2k+1)!, real |u| <= log 2. */
static double zm_dilog_series_r(double u) {
    double v = u * u;
    double p = ZM_DILOG_C[8];
    for (int k = 7; k >= 0; k--) p = p * v + ZM_DILOG_C[k];
    return u - 0.25 * v + u * v * p;
}

bool sf_dilog_r(double x, double* out) {
    if (!isfinite(x) || x > 1.0) return false;
    if (x == 1.0) { *out = ZM_PI2_6_HI; return true; }
    if (x < -1.0) {
        /* Li_2(x) = -Li_2(1/x) - pi^2/6 - log^2(-x)/2 */
        double l = log(-x);
        double y = 1.0 / x;
        double v = -zm_dilog_series_r(-log1p(-y)) - 0.5 * l * l - ZM_PI2_6_HI
                   - ZM_PI2_6_LO;
        *out = v;
        return true;
    }
    if (x <= 0.5) { *out = zm_dilog_series_r(-log1p(-x)); return true; }
    /* Li_2(x) = pi^2/6 - log x log(1-x) - Li_2(1-x); for Li_2(1-x) the series
     * variable is u = -log(1-(1-x)) = -log x. */
    double lx = log(x);
    double v = ZM_PI2_6_HI - lx * log1p(-x) - zm_dilog_series_r(-lx) + ZM_PI2_6_LO;
    *out = v;
    return true;
}

/* Complex Bernoulli series in u (|u| <= ~1.3). */
static void zm_dilog_series_c(double ur, double ui, double* r, double* i) {
    double vr, vi;
    zm_mul(ur, ui, ur, ui, &vr, &vi);
    double pr = ZM_DILOG_C[15], pi = 0.0;
    for (int k = 14; k >= 0; k--) {
        zm_mul(pr, pi, vr, vi, &pr, &pi);
        pr += ZM_DILOG_C[k];
    }
    double qr, qi;
    zm_mul(vr, vi, pr, pi, &qr, &qi);
    zm_mul(ur, ui, qr, qi, &qr, &qi);
    *r = ur - 0.25 * vr + qr;
    *i = ui - 0.25 * vi + qi;
}

/* Li_2(z) for |z| <= 1 (z != 1). */
static void zm_dilog_unit(double zr, double zi, double* r, double* i) {
    if (zr <= 0.5) {
        /* u = -log(1 - z); log|1-z| = log1p(|1-z|^2 - 1)/2 */
        double t = zr * (zr - 2.0) + zi * zi;
        double ur = -0.5 * log1p(t);
        double ui = -atan2(-zi, 1.0 - zr);
        zm_dilog_series_c(ur, ui, r, i);
        return;
    }
    /* Li_2(z) = pi^2/6 - log z log(1-z) - Li_2(1-z); Li_2(1-z) uses u = -log z. */
    double lzr, lzi, l1r, l1i;
    zm_log(zr, zi, &lzr, &lzi);
    double wr = 1.0 - zr, wi = -zi;
    l1r = 0.5 * log(wr * wr + wi * wi);
    l1i = atan2(wi, wr);
    double pr, pi;
    zm_mul(lzr, lzi, l1r, l1i, &pr, &pi);
    double sr, si;
    zm_dilog_series_c(-lzr, -lzi, &sr, &si);
    *r = ZM_PI2_6_HI - pr - sr + ZM_PI2_6_LO;
    *i = -pi - si;
}

bool sf_dilog_c(double zr, double zi, double* vr, double* vi) {
    if (!isfinite(zr) || !isfinite(zi)) return false;
    if (zi == 0.0) {
        if (zr <= 1.0) {
            double v;
            if (!sf_dilog_r(zr, &v)) return false;
            *vr = v; *vi = 0.0;
            return true;
        }
        /* On the cut, from below: Re = pi^2/3 - log^2 x/2 - Li_2(1/x),
         * Im = -pi log x. */
        double l = log(zr), w;
        if (!sf_dilog_r(1.0 / zr, &w)) return false;
        *vr = 2.0 * ZM_PI2_6_HI - 0.5 * l * l - w + 2.0 * ZM_PI2_6_LO;
        *vi = -M_PI * l;
        return true;
    }
    double m2 = zr * zr + zi * zi;
    if (m2 <= 1.0) {
        zm_dilog_unit(zr, zi, vr, vi);
        return isfinite(*vr) && isfinite(*vi);
    }
    /* |z| > 1: Li_2(z) = -Li_2(1/z) - pi^2/6 - log^2(-z)/2 */
    double ir, ii;
    zm_div(1.0, 0.0, zr, zi, &ir, &ii);
    double wr, wi;
    zm_dilog_unit(ir, ii, &wr, &wi);
    double Lr, Li;
    zm_log(-zr, -zi, &Lr, &Li);
    double l2r = 0.5 * (Lr - Li) * (Lr + Li), l2i = Lr * Li;
    *vr = -wr - l2r - ZM_PI2_6_HI - ZM_PI2_6_LO;
    *vi = -wi - l2i;
    return isfinite(*vr) && isfinite(*vi);
}

/* ------------------------------------------------------------------------ */
/* PolyLog, real order, complex argument                                     */
/* ------------------------------------------------------------------------ */

/* Defining series, |z| <= 1/2. */
static bool zm_polylog_series_c(double s, double zr, double zi, double* vr, double* vi) {
    int si_ = (int)s;
    bool s_int = (s == floor(s) && s >= 0.0 && s <= 32.0);
    double Sr = 0.0, Si = 0.0, cr = 0.0, ci = 0.0, A = 0.0;
    double pr = zr, pi = zi, lr = 0.0, li = 0.0;     /* z^k, double-double */
    for (int k = 1; k < 2000; k++) {
        double ks;
        if (s_int) { ks = 1.0; for (int j = 0; j < si_; j++) ks *= (double)k; ks = 1.0 / ks; }
        else       ks = pow((double)k, -s);
        double tr = (pr + lr) * ks, ti = (pi + li) * ks;
        double nr = Sr + tr, ni = Si + ti;             /* Neumaier summation */
        cr += fabs(Sr) >= fabs(tr) ? (Sr - nr) + tr : (tr - nr) + Sr;
        ci += fabs(Si) >= fabs(ti) ? (Si - ni) + ti : (ti - ni) + Si;
        Sr = nr; Si = ni;
        A += fabs(tr) + fabs(ti);
        /* (k+1)^-s / k^-s <= 2^-s bounds the growth of the terms for s < 0, so
         * only stop once they are past their peak (|z| (1+1/k)^-s < 1). */
        if (fabs(tr) + fabs(ti) <= 0.25 * ZM_EPS * (fabs(Sr) + fabs(Si)) &&
            (s >= 0.0 || pow(1.0 + 1.0 / k, -s) * sqrt(zr * zr + zi * zi) < 1.0)) {
            Sr += cr; Si += ci;
            *vr = Sr; *vi = Si;
            if (!isfinite(Sr) || !isfinite(Si)) return false;
            return A <= (s_int ? ZM_CANCEL_MAX : ZM_CANCEL_PL) * (fabs(Sr) + fabs(Si));
        }
        zm_dd_mul(&pr, &pi, &lr, &li, zr, zi);
    }
    return false;
}

/* zeta(n) at an integer n != 1, from the tables (fallback: the kernel). */
static bool zm_zeta_int(int n, double* out) {
    if (n >= 2 && n <= 64) { *out = ZM_ZPOS[n - 2]; return true; }
    if (n <= 0 && n >= -60) { *out = ZM_ZNEG[-n]; return true; }
    return sf_zeta_r((double)n, out);
}

/* Jonquiere expansion in mu = log z (|mu| < 2 pi). */
static bool zm_polylog_jonq_c(double s, double mr, double mi, double* vr, double* vi) {
    bool s_int = (s == floor(s) && s >= 1.0 && s <= 60.0);
    int n = s_int ? (int)s : 0;
    double Hr = 0.0, Hi = 0.0;                         /* the head term */
    if (s_int) {
        /* mu^(n-1)/(n-1)! (H_{n-1} - log(-mu)) */
        double H = 0.0, f = 1.0;
        for (int j = 1; j <= n - 1; j++) { H += 1.0 / j; f *= (double)j; }
        double lr, li;
        zm_log(-mr, -mi, &lr, &li);
        double pr = 1.0, pi = 0.0;
        for (int j = 0; j < n - 1; j++) zm_mul(pr, pi, mr, mi, &pr, &pi);
        zm_mul(pr / f, pi / f, H - lr, -li, &Hr, &Hi);
    } else {
        /* Gamma(1-s) (-mu)^(s-1) */
        double g = tgamma(1.0 - s);
        double lr, li;
        zm_log(-mr, -mi, &lr, &li);
        zm_exp((s - 1.0) * lr, (s - 1.0) * li, &Hr, &Hi);
        Hr *= g; Hi *= g;
    }
    double Sr = 0.0, Si = 0.0, pr = 1.0, pi = 0.0;     /* mu^k / k! */
    double A = hypot(Hr, Hi);                          /* Sum |term| */
    int small = 0;
    for (int k = 0; k < 300; k++) {
        if (!(s_int && k == n - 1)) {
            double z;
            bool ok = s_int ? zm_zeta_int(n - k, &z) : sf_zeta_r(s - (double)k, &z);
            if (!ok) return false;
            double tr = z * pr, ti = z * pi;
            Sr += tr; Si += ti;
            A += fabs(tr) + fabs(ti);
            /* zeta vanishes at the negative even integers: a zero term is not
             * convergence, so demand two consecutive negligible terms. */
            if (k > 2 && fabs(tr) + fabs(ti) <= 0.25 * ZM_EPS * (fabs(Sr + Hr) + fabs(Si + Hi))) {
                if (++small >= 2) {
                    *vr = Hr + Sr; *vi = Hi + Si;
                    if (!isfinite(*vr) || !isfinite(*vi)) return false;
                    /* Near an integer order Gamma(1-s) and zeta(s-k) have
                     * cancelling poles: decline what the doubles cannot hold. */
                    return A <= (s_int ? ZM_CANCEL_MAX : ZM_CANCEL_PL) * hypot(*vr, *vi);
                }
            } else {
                small = 0;
            }
        }
        zm_mul(pr, pi, mr, mi, &pr, &pi);
        pr /= (double)(k + 1); pi /= (double)(k + 1);
    }
    return false;
}

/* Inversion for integer n >= 2, |z| > 1. */
static bool zm_polylog_invert_c(int n, double zr, double zi, double* vr, double* vi) {
    if (n > 40) return false;
    double ir, ii;
    zm_div(1.0, 0.0, zr, zi, &ir, &ii);
    double wr, wi;
    if (!sf_polylog_c((double)n, 0.0, ir, ii, &wr, &wi)) return false;
    /* x = 1/2 + log(-z)/(2 pi i) = 1/2 + Li/(2 pi) - i Lr/(2 pi) */
    double Lr, Li;
    zm_log(-zr, -zi, &Lr, &Li);
    double xr = 0.5 + Li / (2.0 * M_PI), xi = -Lr / (2.0 * M_PI);
    /* B_n(x) = Sum_k C(n,k) B_k x^(n-k), Horner in x */
    double br = 0.0, bi = 0.0, c = 1.0;                /* C(n, k) built downward */
    for (int k = 0; k <= n; k++) {
        /* coefficient of x^(n-k) is C(n,k) B_k; Horner from k = 0 */
        zm_mul(br, bi, xr, xi, &br, &bi);
        br += c * ZM_BN[k];
        c = c * (double)(n - k) / (double)(k + 1);
    }
    /* (2 pi i)^n / n! = ((2 pi)^n / n!) i^n */
    double m = 1.0;
    for (int j = 1; j <= n; j++) m *= 2.0 * M_PI / (double)j;
    double cr, ci;
    switch (n & 3) {
        case 0:  cr = m;    ci = 0.0;  break;
        case 1:  cr = 0.0;  ci = m;    break;
        case 2:  cr = -m;   ci = 0.0;  break;
        default: cr = 0.0;  ci = -m;   break;
    }
    double tr, ti;
    zm_mul(cr, ci, br, bi, &tr, &ti);
    double sg = (n & 1) ? 1.0 : -1.0;                  /* (-1)^(n-1) */
    *vr = sg * wr - tr;
    *vi = sg * wi - ti;
    return isfinite(*vr) && isfinite(*vi);
}

bool sf_polylog_r(double s, double x, double* out) {
    if (!isfinite(s) || !isfinite(x)) return false;
    if (s == 2.0) return sf_dilog_r(x, out);
    if (fabs(x) > 1.0) return false;                   /* complex, or the inversion */
    if (s == 1.0) {                                    /* -log(1 - x) */
        if (x == 1.0) return false;
        *out = -log1p(-x);
        return true;
    }
    if (s == 0.0) {                                    /* x / (1 - x) */
        if (x == 1.0) return false;
        *out = x / (1.0 - x);
        return true;
    }
    if (x == 0.0) { *out = 0.0; return true; }
    if (x == 1.0) return s > 1.0 ? sf_zeta_r(s, out) : false;
    if (x == -1.0) {                                   /* -(1 - 2^(1-s)) zeta(s) */
        double z;
        if (!sf_zeta_r(s, &z)) return false;
        *out = expm1((1.0 - s) * 0.6931471805599453) * z;
        return isfinite(*out);
    }
    double vr, vi;
    if (fabs(x) <= 0.5) {
        if (!zm_polylog_series_c(s, x, 0.0, &vr, &vi)) return false;
    } else {
        if (x < 0.0 && s >= 1.0) {
            /* Duplication: Li_s(-y) = 2^(1-s) Li_s(y^2) - Li_s(y). Benign for
             * s >= 1 (both pieces bounded by zeta(s)); for s < 1 Li_s(y)
             * diverges as y -> 1 and the difference cancels catastrophically. */
            double y = -x, a, b;
            if (!sf_polylog_r(s, y * y, &a) || !sf_polylog_r(s, y, &b)) return false;
            *out = pow(2.0, 1.0 - s) * a - b;
            return isfinite(*out);
        }
        /* Jonquiere in mu = log x, complex (log|x| + i pi) for negative x. */
        if (!zm_polylog_jonq_c(s, log(fabs(x)), x < 0.0 ? M_PI : 0.0, &vr, &vi)) return false;
    }
    *out = vr;
    return isfinite(vr);
}

bool sf_polylog_c(double sr, double si, double zr, double zi,
                  double* vr, double* vi) {
    if (!isfinite(sr) || !isfinite(si) || !isfinite(zr) || !isfinite(zi)) return false;
    if (si != 0.0) return false;                       /* complex order: MPFR path */
    double s = sr;
    if (s == 2.0) return sf_dilog_c(zr, zi, vr, vi);
    if (zi == 0.0) {
        if (zr >= -1.0 && zr <= 1.0) {
            double v;
            if (!sf_polylog_r(s, zr, &v)) return false;
            *vr = v; *vi = 0.0;
            return true;
        }
        if (zr > 1.0) zi = -0.0;                       /* below the cut */
    }
    if (s == 1.0) {                                    /* -log(1 - z) */
        double wr = 1.0 - zr, wi = -zi;
        *vr = -log(hypot(wr, wi));
        *vi = -atan2(wi, wr);
        return isfinite(*vr);
    }
    if (s == 0.0) {                                    /* z / (1 - z) */
        zm_div(zr, zi, 1.0 - zr, -zi, vr, vi);
        return isfinite(*vr) && isfinite(*vi);
    }
    double m2 = zr * zr + zi * zi;
    if (m2 <= 0.25) return zm_polylog_series_c(s, zr, zi, vr, vi);
    double mr, mi;
    zm_log(zr, zi, &mr, &mi);
    if (hypot(mr, mi) < 3.6) return zm_polylog_jonq_c(s, mr, mi, vr, vi);
    if (m2 > 1.0 && s == floor(s) && s >= 2.0)
        return zm_polylog_invert_c((int)s, zr, zi, vr, vi);
    return false;
}
