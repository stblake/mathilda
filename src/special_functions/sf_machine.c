/* Machine-precision (double) kernels for the exponential-integral family.
 * See sf_machine.h for why these exist separately from the MPFR modules.
 *
 * Each function uses the standard split: a power series where it converges
 * quickly and is well conditioned, a continued fraction or asymptotic form
 * beyond.  Accuracy is verified against the modules' own MPFR implementations
 * by tests/test_sf_machine.c rather than asserted here.
 */
#include "sf_machine.h"
#include "sf_machine_gamma.h"
#include "sf_machine_zeta.h"  /* zeta-family kernels */
#include "sf_machine_hyper.h"   /* sf_machine_hyperu_quad */
#include "sf_machine_bessel.h"  /* sf_bessel_i_r (StruveL asymptotic) */

#include <math.h>
#include <complex.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#define EI_EULER   0.57721566490153286061
#define EI_MAXIT   200
#define EI_EPS     1.0e-16
#define EI_FPMIN   1.0e-300

/* E_1, Ei, Si, Ci, Shi, Chi and li moved to sf_machine_expint.c, which holds
 * the whole exponential-integral family (real and complex kernels). */

bool sf_machine_sinc(double x, double* out) {
    *out = (x == 0.0) ? 1.0 : sin(x) / x;
    return isfinite(*out);
}

/* ------------------------------------------------------------------ *
 *  Further special functions                                          *
 * ------------------------------------------------------------------ */

/* erfi(x) = 2/sqrt(pi) * integral_0^x e^{t^2} dt.
 *
 * The Maclaurin series has ALL-POSITIVE terms, so unlike most series of this
 * shape it suffers no cancellation and stays accurate as far as it is
 * practical to sum it; beyond that the asymptotic form takes over.  (The
 * textbook route through Dawson's function is quicker but only good to about
 * single precision in the tabulated form.) */
bool sf_machine_erfi(double x, double* out) {
    double t = fabs(x);
    if (t == 0.0) { *out = 0.0; return true; }
    double v;
    if (t <= 6.0) {
        double term = t, sum = t;                     /* n = 0 */
        for (int n = 1; n <= 400; n++) {
            term *= t * t / (double)n;                /* x^(2n+1)/n! */
            double add = term / (2.0 * n + 1.0);
            sum += add;
            if (add < sum * EI_EPS) break;
        }
        v = 2.0 / sqrt(M_PI) * sum;
    } else {
        /* erfi(x) ~ e^{x^2}/(x sqrt(pi)) * sum (2n-1)!!/(2x^2)^n */
        double sum = 1.0, term = 1.0, y = 2.0 * t * t;
        for (int n = 1; n <= 200; n++) {
            double prev = term;
            term *= (2.0 * n - 1.0) / y;
            if (term > prev) break;                   /* asymptotic series turned */
            sum += term;
            if (term < sum * EI_EPS) break;
        }
        v = exp(t * t) / (t * sqrt(M_PI)) * sum;
    }
    *out = (x < 0.0) ? -v : v;
    return isfinite(*out);
}

/* ProductLog (Lambert W) lives in sf_machine_plog.c with its complex twin. */

/* Fresnel C and S, in the pi/2 convention Mathematica uses:
 *   C(x) = int_0^x cos(pi t^2 / 2) dt,  S(x) = int_0^x sin(pi t^2 / 2) dt.
 * Power series below the crossover, and above it the continued fraction for the
 * complex error function, whose real and imaginary parts carry both at once. */
static void frenel_double(double x, double* c, double* s) {
    const double XCROSS = 1.5;
    double ax = fabs(x);
    if (ax < 1e-150) { *c = ax; *s = 0.0; }
    else if (ax <= XCROSS) {
        double sum = 0.0, sums = 0.0, sumc = ax, sign = 1.0;
        double fact = M_PI / 2.0 * ax * ax, term = ax;
        int odd = 1, n = 3;
        for (int k = 1; k <= 400; k++) {
            term *= fact / k;
            sum += sign * term / n;
            double test = fabs(sum) * EI_EPS;
            if (odd) { sign = -sign; sums = sum; sum = sumc; }
            else     { sumc = sum; sum = sums; }
            if (term < test) break;
            odd = !odd;
            n += 2;
        }
        *s = sums; *c = sumc;
    } else {
        /* Lentz's continued fraction for the complex error function, carried in
         * real/imag pairs rather than `double _Complex`: each `1/z` here would
         * otherwise compile to libgcc's __divdc3 (Annex-G infinity-safe complex
         * divide), which dominated this branch (~13x the power-series branch).
         * The CF intermediates stay O(1)..O(n^2) for the x>1.5 range, so the plain
         * reciprocal z* / |z|^2 is safe, and the values are unchanged to ~1 ulp. */
        double pix2 = M_PI * ax * ax;
        double br = 1.0, bi = -pix2;            /* b  = 1 - i pix2        */
        double cr = 1.0 / EI_FPMIN, ci = 0.0;   /* cc = 1/FPMIN           */
        double den = br * br + bi * bi;
        double dr = br / den, di = -bi / den;   /* d  = 1/b               */
        double hr = dr, hi = di;                /* h  = d                 */
        int n = -1;
        for (int k = 2; k <= 400; k++) {
            n += 2;
            double a = -(double)n * (double)(n + 1);
            br += 4.0;                           /* b += 4 (real part)     */
            double tr = a * dr + br, ti = a * di + bi;   /* a*d + b        */
            den = tr * tr + ti * ti;
            dr = tr / den; di = -ti / den;       /* d = 1/(a*d + b)        */
            double cden = cr * cr + ci * ci;     /* cc = b + a/cc          */
            cr = br + a * cr / cden;
            ci = bi - a * ci / cden;
            double delr = cr * dr - ci * di;     /* del = cc * d           */
            double deli = cr * di + ci * dr;
            double nhr = hr * delr - hi * deli;  /* h *= del               */
            hi = hr * deli + hi * delr;
            hr = nhr;
            if (fabs(delr - 1.0) + fabs(deli) < EI_EPS) break;
        }
        double thr = ax * hr + ax * hi;          /* h = (ax - i ax) * h    */
        double thi = ax * hi - ax * hr;
        double cg = cos(0.5 * pix2), sg = sin(0.5 * pix2);
        double pr = cg * thr - sg * thi;         /* (cos + i sin) * h      */
        double pimag = cg * thi + sg * thr;
        double ur = 1.0 - pr, ui = -pimag;       /* 1 - (...)              */
        *c = 0.5 * (ur - ui);                    /* (1/2 + i/2)(ur + i ui) */
        *s = 0.5 * (ur + ui);
    }
    if (x < 0.0) { *c = -(*c); *s = -(*s); }
}

bool sf_machine_fresnel_c(double x, double* out) {
    double c, s; frenel_double(x, &c, &s); *out = c; return isfinite(c);
}
bool sf_machine_fresnel_s(double x, double* out) {
    double c, s; frenel_double(x, &c, &s); *out = s; return isfinite(s);
}

/* Digamma.  Recurrence up to where the asymptotic series is good, reflection
 * below zero; the poles at the non-positive integers decline. */
/* Digamma: the Gamma family's kernel (sf_machine_gamma.c), so PolyGamma,
 * HarmonicNumber and the hypergeometric digamma sums share one implementation. */
bool sf_machine_digamma(double x, double* out) {
    return sfg_digamma(x, out);
}

bool sf_machine_harmonic(double x, double* out) {
    double d;
    if (!sf_machine_digamma(x + 1.0, &d)) return false;
    *out = d + EI_EULER;
    return isfinite(*out);
}

/* Riemann zeta on the reals: the zeta-family kernel (sf_machine_zeta.c) --
 * tabulated integers, Euler-Maclaurin, the functional equation below 1/2. */
bool sf_machine_zeta(double s, double* out) {
    return sf_zeta_r(s, out);
}

/* Fibonacci and Lucas continued to real argument, which is what the
 * interpreter does too (Fibonacci[2.5] is a Real, not an error). */
bool sf_machine_fibonacci(double x, double* out) {
    const double PHI = 1.61803398874989484820;
    double v = (pow(PHI, x) - cos(M_PI * x) * pow(PHI, -x)) / sqrt(5.0);
    *out = v;
    return isfinite(v);
}
bool sf_machine_lucasl(double x, double* out) {
    const double PHI = 1.61803398874989484820;
    double v = pow(PHI, x) + cos(M_PI * x) * pow(PHI, -x);
    *out = v;
    return isfinite(v);
}

/* Hurwitz zeta on the reals, s >= 1/2 (s != 1), a > 0: the zeta-family kernel
 * (sf_machine_zeta.c), Euler-Maclaurin with an s-dependent shift. */
bool sf_machine_hurwitz_zeta(double s, double a, double* out) {
    return sf_hurwitz_r(s, a, out);
}

/* PolyGamma[n, x] for an integer order: the Gamma family's kernel
 * (sf_machine_gamma.c). */
bool sf_machine_polygamma(double n, double x, double* out) {
    if (n < 0.0 || n != floor(n) || n > (double)SFG_POLYGAMMA_MAX_N) return false;
    return sfg_polygamma((int)n, x, out);
}

/* Sign of Gamma(x), computed rather than taken from lgamma_r or the global
 * signgam: lgamma_r is POSIX rather than C99 and signgam is not thread-safe.
 * Gamma alternates sign on the negative axis, negative on (-1,0), positive on
 * (-2,-1), and so on — so the parity of floor(x) decides it. */
static double gamma_sign(double x) {
    if (x > 0.0) return 1.0;
    long k = (long)floor(x);
    return (k & 1L) ? -1.0 : 1.0;
}

/* Pochhammer[a, n]: the Gamma family's kernel (sf_machine_gamma.c). */
bool sf_machine_pochhammer(double a, double n, double* out) {
    return sfg_pochhammer(a, n, out);
}

/* Binomial[n, k] = Gamma(n+1) / (Gamma(k+1) Gamma(n-k+1)).  The interpreter
 * evaluates this for real arguments (Binomial[5.5, 2.] is 12.375), so a machine
 * kernel is answering the same question, not a different one. */
bool sf_machine_binomial(double n, double k, double* out) {
    double m = n - k;
    if ((k < 0.0 && k == floor(k)) || (m < 0.0 && m == floor(m))) {
        /* Gamma pole in a denominator: the coefficient is zero, unless the
         * numerator is singular too, which the MPFR path should settle. */
        if (n < 0.0 && n == floor(n)) return false;
        *out = 0.0;
        return true;
    }
    if (n < 0.0 && n == floor(n)) return false;
    double l0 = lgamma(n + 1.0), l1 = lgamma(k + 1.0), l2 = lgamma(m + 1.0);
    double v = gamma_sign(n + 1.0) * gamma_sign(k + 1.0) * gamma_sign(m + 1.0)
             * exp(l0 - l1 - l2);
    *out = v;
    return isfinite(v);
}

/* LegendreP[n, x] for a non-negative integer degree, by the three-term
 * recurrence — stable upward for this polynomial family.  A non-integer degree
 * needs the hypergeometric form and is left to the MPFR path, matching how
 * BesselJ/BesselY already decline a non-integer order. */
bool sf_machine_legendre_p(double n, double x, double* out) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    double p0 = 1.0, p1 = x;
    if (N == 0) { *out = 1.0; return true; }
    if (N == 1) { *out = x; return isfinite(x); }
    for (int k = 2; k <= N; k++) {
        double pk = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / (double)k;
        p0 = p1; p1 = pk;
    }
    *out = p1;
    return isfinite(p1);
}

/* ChebyshevT[n, x] for a non-negative integer degree, by the three-term
 * recurrence T_0 = 1, T_1 = x, T_k = 2x T_{k-1} - T_{k-2} -- stable upward for
 * this polynomial family. A non-integer, negative, or very large degree is left
 * to the MPFR / symbolic path, matching the ChebyshevU helper (the seed T_1 = x
 * is the only difference from U_1 = 2x). */
bool sf_machine_chebyshev_t(double n, double x, double* out) {
    if (n != floor(n) || n < 0.0 || !isfinite(x)) return false;
    if (n > 1000.0) {
        /* High degree: the trigonometric / hyperbolic closed form, O(1) instead
         * of O(n).  With x = cos(theta), T_n(x) = cos(n theta); a rounding of x
         * already moves T_n by |T_n'(x)| eps = n |sin(n theta)| eps / sin(theta),
         * and the formula's own error, |sin(n theta)| 2 n theta eps, is within
         * 2 theta sin(theta) <= 3.7 of that -- so it is as accurate as the
         * argument allows, on the whole interval (acos and cos are accurate
         * for every argument).  The recurrence, below, stays for n <= 1000,
         * where it is exact on small integers and dyadic arguments.  Off the
         * interval, T_n(x) = sgn(x)^n cosh(n acosh|x|); an overflow declines
         * (the MPFR path keeps the value). */
        double v;
        if (fabs(x) <= 1.0) {
            v = cos(n * acos(x));
        } else {
            v = cosh(n * acosh(fabs(x)));
            if (x < 0.0 && fmod(n, 2.0) != 0.0) v = -v;
        }
        *out = v;
        return isfinite(v);
    }
    int N = (int)n;
    if (N == 0) { *out = 1.0; return true; }
    double t0 = 1.0, t1 = x;
    if (N == 1) { *out = t1; return isfinite(t1); }
    for (int k = 2; k <= N; k++) {
        double tk = 2.0 * x * t1 - t0;
        t0 = t1; t1 = tk;
    }
    *out = t1;
    return isfinite(t1);
}

/* ChebyshevT[n, x] for a non-negative integer degree with a machine COMPLEX
 * argument x = xre + i*xim, by the same three-term recurrence carried in two
 * doubles (manual re/im, no <complex.h> -- C99-portable).  Declines for a
 * non-integer/negative/large degree or a non-finite result, mirroring the real
 * helper so the ND kernel and the scalar dispatch fall back to the MPFR path
 * exactly where the real path would. */
bool sf_machine_chebyshev_t_c(double n, double xre, double xim,
                              double* rr, double* ri) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    if (N == 0) { *rr = 1.0; *ri = 0.0; return true; }
    double t1r = xre, t1i = xim;                     /* T_1 = x */
    if (N == 1) { *rr = t1r; *ri = t1i; return isfinite(t1r) && isfinite(t1i); }
    double t0r = 1.0, t0i = 0.0;                     /* T_0 = 1 */
    const double tr = 2.0 * xre, ti = 2.0 * xim;     /* 2x */
    for (int k = 2; k <= N; k++) {
        /* tk = 2x*t1 - t0  (complex multiply, then subtract) */
        double pr = tr * t1r - ti * t1i;
        double pi = tr * t1i + ti * t1r;
        double tkr = pr - t0r, tki = pi - t0i;
        t0r = t1r; t0i = t1i;
        t1r = tkr; t1i = tki;
    }
    *rr = t1r; *ri = t1i;
    return isfinite(t1r) && isfinite(t1i);
}

/* ChebyshevU[n, x] for a non-negative integer degree, by the three-term
 * recurrence U_0 = 1, U_1 = 2x, U_k = 2x U_{k-1} - U_{k-2} -- stable upward for
 * this polynomial family. A non-integer, negative, or very large degree is left
 * to the MPFR / symbolic path, matching how sf_machine_legendre_p declines. */
bool sf_machine_chebyshev_u(double n, double x, double* out) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    if (N == 0) { *out = 1.0; return true; }
    double u0 = 1.0, u1 = 2.0 * x;
    if (N == 1) { *out = u1; return isfinite(u1); }
    for (int k = 2; k <= N; k++) {
        double uk = 2.0 * x * u1 - u0;
        u0 = u1; u1 = uk;
    }
    *out = u1;
    return isfinite(u1);
}

/* ChebyshevU[n, x] for a non-negative integer degree with a machine COMPLEX
 * argument x = xre + i*xim, by the same three-term recurrence carried in two
 * doubles (manual re/im, no <complex.h> -- C99-portable).  Declines for a
 * non-integer/negative/large degree or a non-finite result, mirroring the real
 * helper so the ND kernel and the scalar dispatch fall back to the MPFR path
 * exactly where the real path would. */
bool sf_machine_chebyshev_u_c(double n, double xre, double xim,
                              double* rr, double* ri) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    if (N == 0) { *rr = 1.0; *ri = 0.0; return true; }
    double u1r = 2.0 * xre, u1i = 2.0 * xim;        /* U_1 = 2x */
    if (N == 1) { *rr = u1r; *ri = u1i; return isfinite(u1r) && isfinite(u1i); }
    double u0r = 1.0, u0i = 0.0;                     /* U_0 = 1  */
    const double tr = 2.0 * xre, ti = 2.0 * xim;     /* 2x */
    for (int k = 2; k <= N; k++) {
        /* uk = 2x*u1 - u0  (complex multiply, then subtract) */
        double pr = tr * u1r - ti * u1i;
        double pi = tr * u1i + ti * u1r;
        double ukr = pr - u0r, uki = pi - u0i;
        u0r = u1r; u0i = u1i;
        u1r = ukr; u1i = uki;
    }
    *rr = u1r; *ri = u1i;
    return isfinite(u1r) && isfinite(u1i);
}

/* HermiteH[n, x] (physicists') for a non-negative integer degree, by the
 * three-term recurrence H_0 = 1, H_1 = 2x, H_k = 2x H_{k-1} - 2(k-1) H_{k-2} --
 * stable upward for this polynomial family. A non-integer, negative, or very
 * large degree is left to the MPFR / symbolic path, matching how
 * sf_machine_chebyshev_u declines. */
bool sf_machine_hermite_h(double n, double x, double* out) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    if (N == 0) { *out = 1.0; return true; }
    double h0 = 1.0, h1 = 2.0 * x;
    if (N == 1) { *out = h1; return isfinite(h1); }
    for (int k = 2; k <= N; k++) {
        double hk = 2.0 * x * h1 - 2.0 * (double)(k - 1) * h0;
        h0 = h1; h1 = hk;
    }
    *out = h1;
    return isfinite(h1);
}

/* HermiteH[n, x] for a non-negative integer degree with a machine COMPLEX
 * argument x = xre + i*xim, by the same three-term recurrence carried in two
 * doubles (manual re/im, no <complex.h> -- C99-portable). Declines for a
 * non-integer/negative/large degree or a non-finite result, mirroring the real
 * helper so the ND kernel and the scalar dispatch fall back to the MPFR path
 * exactly where the real path would. */
bool sf_machine_hermite_h_c(double n, double xre, double xim,
                            double* rr, double* ri) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    if (N == 0) { *rr = 1.0; *ri = 0.0; return true; }
    double h1r = 2.0 * xre, h1i = 2.0 * xim;        /* H_1 = 2x */
    if (N == 1) { *rr = h1r; *ri = h1i; return isfinite(h1r) && isfinite(h1i); }
    double h0r = 1.0, h0i = 0.0;                     /* H_0 = 1  */
    const double tr = 2.0 * xre, ti = 2.0 * xim;     /* 2x */
    for (int k = 2; k <= N; k++) {
        /* hk = 2x*h1 - 2(k-1)*h0  (complex multiply, then subtract) */
        double pr = tr * h1r - ti * h1i;
        double pi = tr * h1i + ti * h1r;
        double f = 2.0 * (double)(k - 1);
        double hkr = pr - f * h0r, hki = pi - f * h0i;
        h0r = h1r; h0i = h1i;
        h1r = hkr; h1i = hki;
    }
    *rr = h1r; *ri = h1i;
    return isfinite(h1r) && isfinite(h1i);
}

/* LaguerreL[n, a, x] for a non-negative integer degree n, by the stable forward
 * recurrence  L_0 = 1,  L_1 = 1 + a - x,
 *   (k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}
 * (a = 0 recovers the ordinary Laguerre polynomial). Evaluating the recurrence
 * rather than substituting into the explicit polynomial avoids the catastrophic
 * cancellation that ruins the polynomial form for large n. Declines for a
 * non-integer/negative/large degree or a non-finite result, matching how
 * sf_machine_chebyshev_u declines. */
bool sf_machine_laguerre_l_assoc(double n, double a, double x, double* out) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    double l0 = 1.0;                       /* L_0 = 1        */
    if (N == 0) { *out = l0; return true; }
    double l1 = 1.0 + a - x;               /* L_1 = 1 + a - x */
    if (N == 1) { *out = l1; return isfinite(l1); }
    for (int k = 1; k < N; k++) {
        double lk = ((2.0 * k + 1.0 + a - x) * l1 - ((double)k + a) * l0)
                    / (double)(k + 1);
        l0 = l1; l1 = lk;
    }
    *out = l1;
    return isfinite(l1);
}

bool sf_machine_laguerre_l(double n, double x, double* out) {
    return sf_machine_laguerre_l_assoc(n, 0.0, x, out);
}

/* LaguerreL[n, a, x] with a machine COMPLEX argument x = xre + i*xim (order n
 * and parameter a stay real), by the same recurrence carried in two doubles
 * (manual re/im, no <complex.h> -- C99-portable). */
bool sf_machine_laguerre_l_assoc_c(double n, double a, double xre, double xim,
                                   double* rr, double* ri) {
    if (n != floor(n) || n < 0.0 || n > 1000.0) return false;
    int N = (int)n;
    double l0r = 1.0, l0i = 0.0;                   /* L_0 = 1         */
    if (N == 0) { *rr = l0r; *ri = l0i; return true; }
    double l1r = 1.0 + a - xre, l1i = -xim;        /* L_1 = 1 + a - x */
    if (N == 1) { *rr = l1r; *ri = l1i; return isfinite(l1r) && isfinite(l1i); }
    for (int k = 1; k < N; k++) {
        /* c = (2k+1+a) - x  (complex);  num = c*L_k - (k+a)*L_{k-1} */
        double cr = 2.0 * k + 1.0 + a - xre, ci = -xim;
        double pr = cr * l1r - ci * l1i;
        double pi = cr * l1i + ci * l1r;
        double ka = (double)k + a;
        double inv = 1.0 / (double)(k + 1);
        double lkr = (pr - ka * l0r) * inv;
        double lki = (pi - ka * l0i) * inv;
        l0r = l1r; l0i = l1i;
        l1r = lkr; l1i = lki;
    }
    *rr = l1r; *ri = l1i;
    return isfinite(l1r) && isfinite(l1i);
}

bool sf_machine_laguerre_l_c(double n, double xre, double xim,
                             double* rr, double* ri) {
    return sf_machine_laguerre_l_assoc_c(n, 0.0, xre, xim, rr, ri);
}

/* LegendreQ (sf_machine_legendre_q) lives in sf_machine_plog.c with its complex twin. */

/* ---- Airy ---------------------------------------------------------------
 * Ai, Bi and their derivatives from the ascending series
 *     Ai = c1 f - c2 g,   Bi = sqrt(3) (c1 f + c2 g)
 * where f and g are the two solutions of y'' = x y normalised at the origin,
 * generated by the recurrence c_{n+3} = c_n / ((n+3)(n+2)) the ODE gives
 * directly.  The series converges everywhere, but for x > 0 it computes a
 * decaying Ai as a difference of growing terms and loses roughly 2*zeta/ln(10)
 * digits, so past the crossover the asymptotic expansions take over.
 *
 * Where neither form can be trusted the kernel DECLINES and the MPFR path
 * answers — an accuracy gate is the one thing a fast path must never fake. */
#define AIRY_C1 0.355028053887817239260063186004
#define AIRY_C2 0.258819403792806798405183560472
/* The two expansions do not meet in double precision, and that is a property of
 * the functions rather than of this code.  The ascending series computes a
 * decaying Ai as a difference of growing terms and loses about 2*zeta/ln(10)
 * digits; the asymptotic series has a smallest attainable error of order
 * e^(-2 zeta).  Measured against the MPFR implementation: the series holds ~1e-14
 * to |x| = 2.5 (3e-13 by |x| = 3), the asymptotic reaches ~1e-15 from |x| = 8,
 * and in between the best either can do is ~1e-5 at |x| = 3.5.
 *
 * The band between them is covered by a THIRD method: Taylor marching of the
 * defining ODE itself (see airy_march).  Nothing is approximated — the same
 * y'' = x y that generates the ascending series is integrated in short steps
 * from a point where one of the other two expansions is exact. */
#define AIRY_SERIES_MAX 2.5
#define AIRY_ASYMP_MIN  8.0

static bool airy_series(double x, double* ai, double* bi, double* aip, double* bip) {
    double x3 = x * x * x;
    double A = 1.0, B = x;                 /* f and g term values */
    double f = 1.0, g = x;
    double fp = 0.0, gp = 1.0;             /* f' and g' */
    for (int k = 0; k < 400; k++) {
        double A1 = A * x3 / ((3.0 * k + 3.0) * (3.0 * k + 2.0));
        double B1 = B * x3 / ((3.0 * k + 4.0) * (3.0 * k + 3.0));
        f += A1; g += B1;
        if (x != 0.0) {
            fp += (3.0 * k + 3.0) * A1 / x;
            gp += (3.0 * k + 4.0) * B1 / x;
        }
        A = A1; B = B1;
        if (fabs(A1) < fabs(f) * EI_EPS && fabs(B1) < (fabs(g) + EI_EPS) * EI_EPS) break;
    }
    const double S3 = 1.732050807568877293527446341506;
    *ai  = AIRY_C1 * f  - AIRY_C2 * g;
    *bi  = S3 * (AIRY_C1 * f  + AIRY_C2 * g);
    *aip = AIRY_C1 * fp - AIRY_C2 * gp;
    *bip = S3 * (AIRY_C1 * fp + AIRY_C2 * gp);
    return isfinite(*ai) && isfinite(*bi) && isfinite(*aip) && isfinite(*bip);
}

/* u_k and v_k of the asymptotic expansions, by their standard ratios. */
static void airy_uv(double* u, double* v, int n) {
    u[0] = 1.0; v[0] = 1.0;
    for (int k = 1; k <= n; k++) {
        u[k] = u[k - 1] * (6.0 * k - 5.0) * (6.0 * k - 3.0) * (6.0 * k - 1.0)
             / (216.0 * k * (2.0 * k - 1.0));
        v[k] = u[k] * (6.0 * k + 1.0) / (1.0 - 6.0 * k);
    }
}

/* Sum an asymptotic series in 1/zeta, stopping at its smallest term. */
static double airy_asum(const double* c, double zeta, int alt) {
    double sum = 0.0, p = 1.0, prev = HUGE_VAL;
    for (int k = 0; k <= 30; k++) {
        double t = c[k] * p;
        if (fabs(t) > prev) break;                 /* the series has turned */
        sum += alt && (k & 1) ? -t : t;
        prev = fabs(t);
        if (fabs(t) < fabs(sum) * EI_EPS) break;
        p /= zeta;
    }
    return sum;
}

static bool airy_asymptotic(double x, double* ai, double* bi, double* aip, double* bip) {
    double u[31], v[31];
    airy_uv(u, v, 30);
    const double SQPI = 1.772453850905516027298167483341;   /* sqrt(pi) */
    double ax = fabs(x), q = sqrt(sqrt(ax)), zeta = 2.0 / 3.0 * ax * sqrt(ax);

    if (x > 0.0) {
        double e = exp(-zeta);
        *ai  =  e / (2.0 * SQPI * q) * airy_asum(u, zeta, 1);
        *aip = -q * e / (2.0 * SQPI) * airy_asum(v, zeta, 1);
        /* Bi grows like e^{+zeta}; overflows past x ~ 104. */
        double ep = exp(zeta);
        *bi  = ep / (SQPI * q) * airy_asum(u, zeta, 0);
        *bip = q * ep / SQPI  * airy_asum(v, zeta, 0);
    } else {
        /* Oscillatory side: the even and odd parts of the same coefficients
         * ride on sin and cos of zeta + pi/4. */
        double th = zeta + M_PI / 4.0, s = sin(th), c = cos(th);
        double ue = 0.0, uo = 0.0, ve = 0.0, vo = 0.0;
        double p = 1.0, prev = HUGE_VAL;
        for (int k = 0; k <= 30; k++) {
            double tu = u[k] * p, tv = v[k] * p;
            if (fabs(tu) > prev) break;
            double sg = (k / 2) & 1 ? -1.0 : 1.0;          /* (-1)^floor(k/2) */
            if (k & 1) { uo += sg * tu; vo += sg * tv; }
            else       { ue += sg * tu; ve += sg * tv; }
            prev = fabs(tu);
            if (fabs(tu) < EI_EPS) break;
            p /= zeta;
        }
        *ai  =  (s * ue - c * uo) / (SQPI * q);
        *bi  =  (c * ue + s * uo) / (SQPI * q);
        *aip = -q * (c * ve + s * vo) / SQPI;
        *bip =  q * (s * ve - c * vo) / SQPI;
    }
    return isfinite(*ai) && isfinite(*bi) && isfinite(*aip) && isfinite(*bip);
}

/* ---- the band 2.5 < |x| < 8: Taylor marching of y'' = x y ----------------
 *
 * Neither expansion reaches double precision here, so the ODE is integrated
 * instead.  About the point x0, writing y(x0 + h) = sum c_k h^k and matching
 * y'' = (x0 + h) y term by term gives
 *
 *     (k+2)(k+1) c_{k+2} = x0 c_k + c_{k-1},        c_{-1} = 0,
 *
 * with c_0 = y(x0) and c_1 = y'(x0).  This is exact, not an approximation: Ai
 * and Bi are entire, so the series converges for every h, and a short step makes
 * it converge in a handful of terms.  No fitted coefficients, no new special
 * function, and the same recurrence shape the ascending series already uses.
 *
 * DIRECTION IS THE WHOLE PROBLEM, and getting it wrong is the classic way this
 * method quietly loses digits.  y'' = x y has one recessive and one dominant
 * solution, and rounding error in the recessive one gets amplified by the
 * dominant one's growth.  For x > 0, Ai decays and Bi grows by a factor of
 * e^(zeta(8) - zeta(2.5)) ~ 2.4e5 across the band — so marching Ai FORWARD from
 * 2.5 would let a 1e-16 seed error grow to ~1e-11.  Each solution is therefore
 * marched in the direction in which it DOMINATES:
 *
 *   Bi  forward  from  |x| = 2.5 (ascending series exact there), growing;
 *   Ai  backward from  |x| = 8   (asymptotic exact there), growing as x falls.
 *
 * For x < 0 both solutions oscillate with comparable amplitude, neither
 * dominates, and one forward march from -2.5 carries both. */
#define AIRY_MARCH_H  0.25          /* step; h^2 = 1/16 makes the series short */
#define AIRY_MARCH_N  24            /* terms per step, ample at that h */

/* Integrate y'' = x y from `xa` to `xb` (either direction), advancing the
 * solution given by (*y, *yp) in place. */
static void airy_march(double xa, double xb, double* y, double* yp) {
    double span = xb - xa;
    int nstep = (int)(fabs(span) / AIRY_MARCH_H) + 1;
    double h = span / (double)nstep;
    double x0 = xa;
    for (int s = 0; s < nstep; s++) {
        double c[AIRY_MARCH_N + 3];
        c[0] = *y; c[1] = *yp;
        for (int k = 0; k <= AIRY_MARCH_N; k++)
            c[k + 2] = (x0 * c[k] + (k >= 1 ? c[k - 1] : 0.0))
                     / (((double)k + 2.0) * ((double)k + 1.0));
        /* Horner, from the smallest term down, for both the value and the
         * derivative — summing the tail first keeps the rounding at the level of
         * the terms that contribute least. */
        double sv = 0.0, sd = 0.0;
        for (int k = AIRY_MARCH_N + 2; k >= 1; k--) {
            sv = sv * h + c[k];
            sd = sd * h + (double)k * c[k];
        }
        *y  = sv * h + c[0];
        *yp = sd;
        x0 += h;
    }
}

static bool airy_band(double x, double* ai, double* bi, double* aip, double* bip) {
    double a, b, ap, bp;
    if (x > 0.0) {
        /* Bi is dominant going up: seed from the series and march forward. */
        if (!airy_series(AIRY_SERIES_MAX, &a, &b, &ap, &bp)) return false;
        *bi = b; *bip = bp;
        airy_march(AIRY_SERIES_MAX, x, bi, bip);
        /* Ai is dominant going down: seed from the asymptotic and march back. */
        if (!airy_asymptotic(AIRY_ASYMP_MIN, &a, &b, &ap, &bp)) return false;
        *ai = a; *aip = ap;
        airy_march(AIRY_ASYMP_MIN, x, ai, aip);
    } else {
        if (!airy_series(-AIRY_SERIES_MAX, &a, &b, &ap, &bp)) return false;
        *ai = a; *aip = ap; *bi = b; *bip = bp;
        airy_march(-AIRY_SERIES_MAX, x, ai, aip);
        airy_march(-AIRY_SERIES_MAX, x, bi, bip);
    }
    return isfinite(*ai) && isfinite(*bi) && isfinite(*aip) && isfinite(*bip);
}

static bool airy_all(double x, double* ai, double* bi, double* aip, double* bip) {
    double ax = fabs(x);
    if (ax <= AIRY_SERIES_MAX) return airy_series(x, ai, bi, aip, bip);
    if (ax >= AIRY_ASYMP_MIN)  return airy_asymptotic(x, ai, bi, aip, bip);
    return airy_band(x, ai, bi, aip, bip);
}

bool sf_machine_airy_ai(double x, double* out) {
    double a, b, ap, bp; if (!airy_all(x, &a, &b, &ap, &bp)) return false;
    *out = a; return isfinite(a);
}
bool sf_machine_airy_bi(double x, double* out) {
    double a, b, ap, bp; if (!airy_all(x, &a, &b, &ap, &bp)) return false;
    *out = b; return isfinite(b);
}
bool sf_machine_airy_ai_prime(double x, double* out) {
    double a, b, ap, bp; if (!airy_all(x, &a, &b, &ap, &bp)) return false;
    *out = ap; return isfinite(ap);
}
bool sf_machine_airy_bi_prime(double x, double* out) {
    double a, b, ap, bp; if (!airy_all(x, &a, &b, &ap, &bp)) return false;
    *out = bp; return isfinite(bp);
}

/* ---- Modified Bessel I, and friends ------------------------------------- */

/* I_nu(x) from the ascending series sum (x/2)^(2k+nu) / (k! Gamma(k+nu+1)).
 * Every term is positive for x > 0, so unlike most series of this shape there is
 * no cancellation at all and the accuracy holds right up to where I_nu itself
 * overflows.  The first term is formed through lgamma so a large nu cannot
 * overflow Gamma before the ratio is taken. */
bool sf_machine_bessel_i(double nu, double x, double* out) {
    double sgn = 1.0;
    if (x < 0.0) {
        /* I_nu(-x) = (-1)^nu I_nu(x) only for integer nu; otherwise complex. */
        if (nu != floor(nu)) return false;
        if (fmod(fabs(nu), 2.0) == 1.0) sgn = -1.0;
        x = -x;
    }
    if (nu < 0.0 && nu == floor(nu)) nu = -nu;        /* I_-n = I_n */
    if (x == 0.0) { *out = (nu == 0.0) ? 1.0 : 0.0; return true; }

    double h = x / 2.0;
    double lt = nu * log(h) - lgamma(nu + 1.0);
    if (lt < -745.0) { *out = 0.0; return true; }     /* underflows to zero */
    double term = gamma_sign(nu + 1.0) * exp(lt);
    double sum = term, h2 = h * h;
    for (int k = 0; k < 4000; k++) {
        term *= h2 / ((k + 1.0) * (k + 1.0 + nu));
        sum += term;
        if (fabs(term) < fabs(sum) * EI_EPS) break;
    }
    *out = sgn * sum;
    return isfinite(*out);
}

/* 0F1, 1F1, 2F1 and pFq live in sf_machine_hyper.c. */

/* QPochhammer[a, q] = prod_{k>=0} (1 - a q^k), convergent for |q| < 1. */
bool sf_machine_qpochhammer(double a, double q, double* out) {
    if (!(fabs(q) < 1.0)) return false;
    double p = 1.0, qk = 1.0;
    for (int k = 0; k < 100000; k++) {
        double f = 1.0 - a * qk;
        p *= f;
        if (p == 0.0) break;
        qk *= q;
        if (fabs(a * qk) < EI_EPS * 0.5) break;       /* remaining factors are 1 */
    }
    *out = (p == 0.0) ? 0.0 : p;    /* a vanishing factor is an exact zero, not -0. */
    return isfinite(p);
}

/* K_nu(x), x > 0.
 *
 * Three regimes, because no single form covers the range in double:
 *   - large x: the asymptotic sqrt(pi/2x) e^-x sum a_k/x^k, which is where the
 *     I-difference below would lose everything to cancellation (I_nu and I_-nu
 *     both grow like e^x while K decays like e^-x);
 *   - nu far from an integer: K = pi/2 (I_-nu - I_nu)/sin(nu pi) directly;
 *   - nu at (or near) an integer, where that formula is 0/0: the classical
 *     log-series, whose psi terms reuse the digamma above.
 */
/* Steed-Barnett continued fraction for K_mu and K_{mu+1}, |mu| <= 1/2, x > 2.
 *
 * This is the piece that makes K usable in the middle of its range.  Both the
 * ascending log-series and the I_-nu - I_nu difference compute a decaying K from
 * quantities that grow like e^x, so they lose ~2x/ln(10) digits — measured
 * against MPFR, 3e-10 by x = 8 and no correct digits at all by x = 20 — while
 * the large-x asymptotic only reaches machine precision beyond x ~ 20.  The
 * continued fraction has neither problem. */
static bool bessk_cf2(double xmu, double x, double* kmu, double* kmu1) {
    double xmu2 = xmu * xmu;
    double a1 = 0.25 - xmu2;
    double b = 2.0 * (1.0 + x);
    double d = 1.0 / b, delh = d, h = d;
    double q1 = 0.0, q2 = 1.0;
    double c = a1, q = c, a = -a1, s = 1.0 + q * delh;
    for (int i = 2; i <= 20000; i++) {
        a -= 2.0 * (i - 1);
        c = -a * c / i;
        double qnew = (q1 - b * q2) / a;
        q1 = q2; q2 = qnew;
        q += c * qnew;
        b += 2.0;
        d = 1.0 / (b + a * d);
        delh = (b * d - 1.0) * delh;
        h += delh;
        double dels = q * delh;
        s += dels;
        if (fabs(dels / s) < EI_EPS) break;
    }
    h = a1 * h;
    *kmu  = sqrt(M_PI / (2.0 * x)) * exp(-x) / s;
    *kmu1 = *kmu * (xmu + x + 0.5 - h) / x;
    return isfinite(*kmu) && isfinite(*kmu1);
}

bool sf_machine_bessel_k(double nu, double x, double* out) {
    if (!(x > 0.0)) return false;                     /* complex for x < 0 */
    nu = fabs(nu);                                    /* K_-nu = K_nu */

    if (x > 2.0) {
        /* Continued fraction at the fractional order, then UPWARD recurrence in
         * nu, which is the stable direction for K. */
        int nl = (int)(nu + 0.5);
        double xmu = nu - (double)nl;
        double kmu, kmu1;
        if (!bessk_cf2(xmu, x, &kmu, &kmu1)) return false;
        for (int i = 1; i <= nl; i++) {
            double knew = kmu + 2.0 * (xmu + (double)i) / x * kmu1;
            kmu = kmu1; kmu1 = knew;
            if (!isfinite(kmu1)) return false;
        }
        *out = kmu;
        return isfinite(kmu) && kmu > 0.0;
    }

    double n_round = floor(nu + 0.5);
    if (fabs(nu - n_round) > 1e-6) {                  /* generic order, small x */
        double ip, im;
        if (!sf_machine_bessel_i(nu, x, &ip)) return false;
        if (!sf_machine_bessel_i(-nu, x, &im)) return false;
        double v = M_PI / 2.0 * (im - ip) / sin(nu * M_PI);
        *out = v;
        return isfinite(v);
    }

    /* Integer order n, small x: K_n(x) = 1/2 sum_{k<n} (-1)^k (n-k-1)!/k! (x/2)^(2k-n)
     *                                 + (-1)^(n+1) sum_k (x/2)^(2k+n)/(k!(n+k)!)
     *                                   * [log(x/2) - (psi(k+1) + psi(n+k+1))/2] */
    int n = (int)n_round;
    if (n > 200) return false;
    double h = x / 2.0, lh = log(h), v = 0.0;
    for (int k = 0; k < n; k++) {                     /* the finite part */
        double c = exp(lgamma((double)(n - k)) - lgamma((double)(k + 1)));
        v += 0.5 * ((k & 1) ? -1.0 : 1.0) * c * pow(h, 2.0 * k - (double)n);
    }
    double sgn = ((n + 1) & 1) ? -1.0 : 1.0;
    double term = pow(h, (double)n) / tgamma((double)n + 1.0);
    double psi1, psin;
    if (!sf_machine_digamma(1.0, &psi1)) return false;
    if (!sf_machine_digamma((double)n + 1.0, &psin)) return false;
    double sm = 0.0, pk = psi1, pnk = psin;
    for (int k = 0; k < 4000; k++) {
        double add = term * (lh - 0.5 * (pk + pnk));
        sm += add;
        if (fabs(add) < fabs(sm) * EI_EPS && k > n) break;
        term *= h * h / ((k + 1.0) * (k + 1.0 + n));
        pk  += 1.0 / (k + 1.0);                       /* psi(k+2) = psi(k+1) + 1/(k+1) */
        pnk += 1.0 / (k + 1.0 + n);
    }
    v += sgn * sm;
    *out = v;
    return isfinite(v);
}

/* Li_s(x) on the reals: the zeta-family kernel (sf_machine_zeta.c) -- the
 * dilogarithm's Bernoulli series for s = 2 (every real x <= 1), closed forms
 * for s = 0, 1, the power series for |x| <= 1/2 and the Jonquiere expansion in
 * log x (complex for negative x) up to |x| = 1. Beyond |x| = 1 the value is
 * complex (or, for x < -1, reached through the complex kernel) and this real
 * entry declines. */
bool sf_machine_polylog(double s, double x, double* out) {
    return sf_polylog_r(s, x, out);
}

/* ---- n-ary ------------------------------------------------------------- */

/* LerchPhi[z, s, a] = sum_{k>=0} z^k / (k+a)^s on the reals: the zeta-family
 * kernel (sf_machine_zeta.c) -- the defining series for |z| < 1 with a
 * double-double running power and a cancellation guard, the Hurwitz zeta at
 * z = 1. */
bool sf_machine_lerchphi(double z, double s, double a, double* out) {
    return sf_lerch_r(z, s, a, out);
}

/* WhittakerM M_{k,m}(z) = e^(-z/2) z^(m+1/2) 1F1(m-k+1/2; 2m+1; z).
 *
 * Real machine path. Declines (returns false) for z <= 0 (the z^(m+1/2) prefactor
 * has a branch cut along the negative real axis and the value is generally
 * complex -- handled by the scalar/MPFR 1F1 template instead), and relays the
 * 1F1 kernel's own declines (the 2m+1 = non-positive-integer pole, and bad
 * ascending-series cancellation at large z). A single growing PRODUCT, so no
 * cancellation beyond 1F1's own: pow overflow / exp underflow are caught by the
 * final isfinite check. Shared bit-for-bit with the NDArray kernel and Compile VM
 * (ndkernels.c ndk_WhittakerM_n) and the scalar machine fast path. */
bool sf_machine_whittaker_m(double k, double m, double z, double* out) {
    if (!(z > 0.0)) return false;              /* z <= 0 (and NaN) -> complex / template */
    double a = m - k + 0.5;
    double b = 2.0 * m + 1.0;
    double h;
    if (!sf_machine_hyper1f1(a, b, z, &h)) return false;   /* pole & cancellation declined inside */
    double e = m + 0.5;
    double v = exp(-0.5 * z) * (e == 1.0 ? z : pow(z, e)) * h;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* ---- HypergeometricU U(a,b,z) (Tricomi), real machine path for z > 0 -------- */

/* Regime 1 -- the asymptotic series U(a,b,z) ~ z^(-a) 2F0(a, a-b+1; ; -1/z),
 * summed to its least term. The series is divergent (asymptotic), so acceptance
 * is MEASURED, not assumed: we stop when a term exceeds its predecessor and keep
 * the result only when that least term is negligible against the partial sum.
 * Terminating cases (a or a-b+1 a non-positive integer) are a finite sum, not an
 * asymptotic one, and are handed to hyperu_terminating before any truncation. */
/* z^-a; a = 1 and 2 by division (<= 1 ulp, like pow, and far cheaper -- U(1,2,z)
 * = 1/z is the common terminating case). */
static double hyperu_zpow(double z, double a) {
    if (a == 1.0) return 1.0 / z;
    if (a == 2.0) return 1.0 / (z * z);
    return pow(z, -a);
}

/* Is x a non-positive integer (a pole of Gamma / a terminating Pochhammer)? */
static bool hyperu_nonpos_int(double x) { return x <= 0.0 && x == floor(x); }

/* The TERMINATING case, a or a-b+1 = -m (m a non-negative integer): U is then
 * exactly z^(-a) times the finite sum sum_{k=0}^m (a)_k (a-b+1)_k / k! (-1/z)^k
 * -- for a = -m the polynomial (-1)^m (b)_m 1F1(-m; b; z) (DLMF 13.2.7), for
 * a-b+1 = -m its Kummer image z^(1-b) U(-m, 2-b, z). It is NOT an asymptotic
 * series here, so it must not be cut at its "least term": the terms of
 * U(-1, 2, z) = z - 2 grow (1, -2/z) for z < 2, and treating that growth as
 * divergence used to discard the series and fall into the confluent log formula,
 * whose two 1/Gamma prefactors both vanish at a = -m and returned 0.0. Every
 * term is summed (compensated); the only error is cancellation, measured as
 * peak term over result and declined (-> FLINT) beyond sf_series_usable. */
static bool hyperu_terminating(double a, double a2, double z, double* out) {
    double m = hyperu_nonpos_int(a) ? -a : -a2;  /* the shorter terminating side */
    if (hyperu_nonpos_int(a) && hyperu_nonpos_int(a2) && -a2 < m) m = -a2;
    if (m > 100000.0) return false;
    int mi = (int)m;
    double term = 1.0, sum = 1.0, comp = 0.0, peak = 1.0;
    for (int k = 0; k < mi; k++) {
        term = term * ((a + k) * (a2 + k)) / ((k + 1.0) * -z);
        if (!isfinite(term)) return false;
        double at = fabs(term);
        if (at > peak) peak = at;
        double s2 = sum + term, bb = s2 - sum;
        comp += (sum - (s2 - bb)) + (term - bb);
        sum = s2;
    }
    double s = sum + comp;
    if (!sf_series_usable(peak, fabs(s))) return false;
    double v = hyperu_zpow(z, a) * s;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

static bool hyperu_asymptotic(double a, double b, double z, double* out) {
    if (!(z > 0.0)) return false;
    double a2 = a - b + 1.0;
    if (hyperu_nonpos_int(a) || hyperu_nonpos_int(a2))
        return hyperu_terminating(a, a2, z, out);
    /* Compensated (TwoSum) summation: the series alternates and runs to ~z
     * terms before its least one. Each ratio divides by -z itself rather than
     * multiplying by a once-rounded -1/z, whose error would repeat in every
     * term (a systematic k*u/2 drift). */
    double term = 1.0, sum = 1.0, comp = 0.0, last = 1.0;
    for (int k = 0; k < 400; k++) {
        double next = term * ((a + k) * (a2 + k)) / ((k + 1.0) * -z);
        double an = fabs(next);
        if (!isfinite(next)) return false;
        if (next == 0.0) {                      /* series terminated: exact */
            double v = hyperu_zpow(z, a) * (sum + comp);
            if (!isfinite(v)) return false;
            *out = v; return true;
        }
        if (an > last) {                        /* least term reached -> truncate */
            /* The least term IS the truncation error: accept it only below
             * half an ulp of the sum (it allowed 1e-13, ~450 ulp, which put
             * U(7/10, 13/10, 28.5) 334 ulp off). */
            if (!(last <= fabs(sum) * 5.5511151231257827e-17)) return false;
            double v = hyperu_zpow(z, a) * (sum + comp);
            if (!isfinite(v)) return false;
            *out = v; return true;
        }
        double s2 = sum + next, bb = s2 - sum;
        comp += (sum - (s2 - bb)) + (next - bb);
        sum = s2;
        term = next;
        last = an;
        if (an <= fabs(sum) * 1.0e-16) {        /* numerically converged */
            double v = hyperu_zpow(z, a) * (sum + comp);
            if (!isfinite(v)) return false;
            *out = v; return true;
        }
    }
    return false;                               /* never reached the least term */
}

/* 1/Gamma(x), finite everywhere: 0 at the non-positive-integer poles of Gamma,
 * where libm tgamma returns a NaN/Inf that would otherwise poison the sum. */
static double hyperu_rgamma(double x) {
    if (x <= 0.0 && x == floor(x)) return 0.0;
    return 1.0 / tgamma(x);
}

/* Regime 3 -- the Kummer connection formula for NON-integer b:
 *   U = Gamma(1-b)/Gamma(a-b+1) 1F1(a,b,z)
 *     + Gamma(b-1)/Gamma(a)     z^(1-b) 1F1(a-b+1, 2-b, z).
 * A difference of two growing terms, so a cancellation guard declines (-> FLINT)
 * when the terms dwarf their sum. */
static bool hyperu_connection(double a, double b, double z, double* out) {
    if (!(z > 0.0)) return false;
    double m1, m2;
    if (!sf_machine_hyper1f1(a, b, z, &m1)) return false;
    if (!sf_machine_hyper1f1(a - b + 1.0, 2.0 - b, z, &m2)) return false;
    double g1 = tgamma(1.0 - b) * hyperu_rgamma(a - b + 1.0);   /* = G(1-b)/G(a-b+1) */
    double g2 = tgamma(b - 1.0) * hyperu_rgamma(a);             /* = G(b-1)/G(a)     */
    double t1 = g1 * m1;
    double t2 = g2 * pow(z, 1.0 - b) * m2;
    double v = t1 + t2;
    if (!isfinite(v)) return false;
    /* Cancellation guard: the two terms' rounding is amplified by
     * (|t1|+|t2|)/|v|. This used to allow 1e13 -- U(5.5, 0.3, 3) came back
     * 1.5e-10 off, U(1.3, 5.5, 10) 60 ulp off -- so decline beyond a factor of 4
     * and let the rigorous FLINT path answer. */
    if (fabs(t1) + fabs(t2) > fabs(v) * 4.0) return false;
    *out = v;
    return true;
}

/* Regime 2 -- INTEGER b = n+1 (n >= 0), the confluent limit DLMF 13.2.9:
 *   U(a,n+1,z) = (-1)^(n+1)/(n! Gamma(a-n)) *
 *                  sum_{k>=0} (a)_k/((n+1)_k k!) z^k
 *                             [ln z + psi(a+k) - psi(k+1) - psi(n+k+1)]
 *              + (n-1)!/Gamma(a) z^(-n)
 *                  sum_{k=0}^{n-1} (a-n)_k/((1-n)_k k!) z^k         (empty for n=0).
 * The 1/Gamma factors vanish on the terminating sub-cases (a or a-n a non-positive
 * integer), so those collapse to the finite polynomial correctly -- EXCEPT a
 * itself a non-positive integer: then a-n is one too, BOTH prefactors vanish and
 * the formula degenerates to 0 (the limit needs l'Hopital; DLMF 13.2.9 excludes
 * a = 0, -1, ...). That case is the polynomial hyperu_terminating sums; reaching
 * here means it declined on cancellation, so decline too (-> FLINT). */
static bool hyperu_integer_b(double a, int n, double z, double* out) {
    if (hyperu_nonpos_int(a)) return false;
    double finite = 0.0;
    if (n >= 1) {
        double rga = hyperu_rgamma(a);                   /* 1/Gamma(a) (0 at a pole) */
        if (rga != 0.0) {
            double sum = 1.0, term = 1.0;                /* k = 0 term */
            for (int k = 1; k <= n - 1; k++) {
                /* term_k/term_{k-1} = (a-n+k-1)/(k-n) * z/k */
                term *= (a - n + (k - 1)) * z / ((double)(k - n) * (double)k);
                sum += term;
            }
            finite = tgamma((double)n) * rga * pow(z, -(double)n) * sum;  /* (n-1)! */
        }
    }
    double logpart = 0.0, peak = fabs(finite);
    double rgan = hyperu_rgamma(a - (double)n);          /* 1/Gamma(a-n) */
    if (rgan != 0.0) {
        double lnz = log(z);
        double pref = ((n & 1) ? 1.0 : -1.0) / tgamma((double)n + 1.0) * rgan;
        double term = 1.0, sum = 0.0;                    /* (a)_0/((n+1)_0) z^0/0! = 1 */
        /* The three digammas advance by the recurrence psi(x+1) = psi(x) + 1/x,
         * so seed each ONCE and step it by one division per term instead of a
         * full sf_machine_digamma call (a shift loop + log + 7-term asymptotic)
         * three times per term -- the dominant cost of this log series. The
         * asymptotic branch has already peeled off every terminating case, so a
         * is never a non-positive integer here and a+k never crosses a pole. */
        double pa, pk, pnk;
        if (!sf_machine_digamma(a, &pa)) return false;           /* psi(a)     */
        if (!sf_machine_digamma(1.0, &pk)) return false;         /* psi(1)     */
        if (!sf_machine_digamma((double)(n + 1), &pnk)) return false; /* psi(n+1) */
        /* The bracket ln z + psi(a+k) - psi(k+1) - psi(n+k+1) is itself a
         * difference: its rounding is ~u (|ln z| + |psi|s), not u |bracket|, so
         * the condition measure accumulates |term| times that scale. */
        double scale = 0.0;
        for (int k = 0; k < 200000; k++) {
            double add = term * (lnz + pa - pk - pnk);
            sum += add;
            scale += fabs(term) * (fabs(lnz) + fabs(pa) + fabs(pk) + fabs(pnk));
            if (fabs(add) < fabs(sum) * EI_EPS && k > 2) break;
            /* advance to k+1: term and the three digammas by their recurrences */
            term *= (a + k) / (double)(n + 1 + k) * z / (double)(k + 1);
            if (!isfinite(term)) return false;
            if (term == 0.0) break;                      /* (a)_k terminated it */
            pa  += 1.0 / (a + (double)k);                /* psi(a+k+1)         */
            pk  += 1.0 / (double)(k + 1);                /* psi(k+2)           */
            pnk += 1.0 / (double)(n + k + 1);            /* psi(n+k+2)         */
        }
        logpart = pref * sum;
        if (fabs(pref) * scale > peak) peak = fabs(pref) * scale;
    }
    double v = finite + logpart;
    if (!isfinite(v)) return false;
    /* Cancellation guard: this is a growing-term (1F1-like) series at large z, and
     * the finite and log parts can nearly cancel. When the largest contribution
     * dwarfs the result the double has thrown most of its digits away -- DECLINE
     * (as sf_machine_pfq does) so the caller takes the rigorous FLINT path rather
     * than return a confidently wrong value. (It allowed 1e6 -- ~10 good
     * digits -- which is not a machine-precision answer.) */
    if (peak > fabs(v) * 4.0) return false;
    *out = v;
    return true;
}

bool sf_machine_hyperu(double a, double b, double z, double* out) {
    if (!(z > 0.0)) return false;                        /* branch cut -> complex */
    /* 1. large-z asymptotic (and every exactly-terminating case). */
    if (hyperu_asymptotic(a, b, z, out)) return true;
    /* 2/3. integer vs non-integer b. */
    double rb = floor(b + 0.5);
    if (fabs(b - rb) < 1.0e-9) {                         /* b at an integer */
        if (rb >= 1.0) {
            if (hyperu_integer_b(a, (int)rb - 1, z, out)) return true;
        } else {
            /* b <= 0 integer: fold up by U(a,b,z) = z^(1-b) U(a-b+1, 2-b, z). */
            double v;
            if (!sf_machine_hyperu(a - b + 1.0, 2.0 - b, z, &v)) return false;
            double r = pow(z, 1.0 - b) * v;
            if (!isfinite(r)) return false;
            *out = r; return true;
        }
    } else if (hyperu_connection(a, b, z, out)) {        /* 3. non-integer b */
        return true;
    }
    /* 4. moderate z, where the series above cancel: the integral (a > 0), or
     *    for a <= 0 its Kummer image U = z^(1-b) U(a-b+1, 2-b, z). */
    if (a > 0.0) return sf_machine_hyperu_quad(a, b, z, out);
    if (a - b + 1.0 > 0.0) {
        double v;
        if (!sf_machine_hyperu_quad(a - b + 1.0, 2.0 - b, z, &v)) return false;
        double r = pow(z, 1.0 - b) * v;
        if (!isfinite(r)) return false;
        *out = r; return true;
    }
    return false;
}

/* WhittakerW W_{k,m}(z) = e^(-z/2) z^(m+1/2) U(m-k+1/2, 2m+1, z).
 * W_{k,m} = W_{k,-m}, so fold m -> |m| (b = 2|m|+1 >= 1). Mirrors
 * sf_machine_whittaker_m with U in place of 1F1; shared with the NDArray kernel
 * and the Compile VM (ndk_WhittakerW_n). Declines (z <= 0, cancellation, pole,
 * overflow) are relayed from sf_machine_hyperu / the final isfinite check. */
bool sf_machine_whittaker_w(double k, double m, double z, double* out) {
    if (!(z > 0.0)) return false;              /* z <= 0 (and NaN) -> complex / template */
    double am = fabs(m);                       /* even in m: fold so b = 2|m|+1 >= 1 */
    double u;
    if (!sf_machine_hyperu(am - k + 0.5, 2.0 * am + 1.0, z, &u)) return false;
    double e = am + 0.5;
    double v = exp(-0.5 * z) * (e == 1.0 ? z : pow(z, e)) * u;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* An ascending (power) series converges everywhere, but it is not USABLE
 * everywhere: wherever the value is small compared with the terms that build it
 * (which is most of the plane once |z| is large), the sum is a difference of
 * quantities far bigger than itself and the double has already thrown the answer
 * away.  Peak term over result is exactly how many bits went, so it is measured
 * rather than guessed from |z| — the loss depends on the direction of z, not
 * only its modulus.
 *
 * Above the budget the kernel DECLINES and the MPFR path answers.  A fast path
 * that is quietly wrong in half the plane would be far worse than a slow one. */
#define SF_SERIES_MAX_LOSS 1.0e3      /* ~10 bits lost; leaves ~43 for the answer */

bool sf_series_usable(double peak, double result) {
    if (!(peak > 0.0)) return true;                 /* nothing accumulated */
    if (!(result > 0.0)) return false;              /* total cancellation */
    return peak <= result * SF_SERIES_MAX_LOSS;
}

/* ------------------------------------------------------------------ */
/*  StruveH H_nu(z) -- real machine path                               */
/* ------------------------------------------------------------------ */
/* H_nu(z) = (z/2)^(nu+1) sum_{k>=0} (-1)^k (z/2)^(2k)/(Gamma(k+3/2)Gamma(k+nu+3/2))
 *
 * Two elementary regimes, each returning both a value and a PREDICTED relative
 * error, so the dispatcher can pick the more accurate at every z rather than
 * declining in a "hard band":
 *   - the convergent 1F2-form series (DLMF 11.2.1). Accurate for small/moderate
 *     |z|; its error is the measured cancellation maxterm/|sum| * eps.
 *   - the large-|z| asymptotic H_nu = Y_nu + K_nu (DLMF 11.6.1 + 10.17.4), each
 *     an optimally-truncated divergent series whose error is the size of its
 *     smallest (first omitted) term.
 * The two error curves cross near z ~ 12 (nu small), where double precision
 * bottoms out at ~2e-11 -- the intrinsic floor there (scipy's struve loses
 * accuracy in the same band), so the kernel returns the best double answer and
 * declines ONLY when even the better regime is worse than STRUVE_ACCEPT: the
 * large-nu uniform regime, a complex (z<0, non-integer nu+1) prefactor, or a
 * non-finite result. (nu = -3/2, -5/2, ..., where Gamma(nu+3/2) has a pole, is
 * NOT a decline: the series drops its vanishing 1/Gamma terms and K_nu == 0.) The caller then falls to the rigorous
 * arb bridge. Shared bit-for-bit with the NDArray kernel and the Compile VM. */

#define STRUVE_ACCEPT 1.0e-7   /* decline if even the better regime exceeds this */

/* The convergent 1F2-form series. Sets *val and the predicted relative error
 * *relerr (cancellation). Returns false only for a Gamma overflow, a complex
 * prefactor, non-convergence, or a non-finite result. */
static bool struve_series(double nu, double z, int sgn, double* val, double* relerr) {
    double w = 0.25 * z * z;                            /* (z/2)^2 >= 0           */
    /* sgn = -1: H_nu (alternating); sgn = +1: L_nu, the modified function, whose
     * series is the same with every term of one sign (DLMF 11.2.2). */
    double ws = sgn < 0 ? -w : w;
    double hump = 0.5 * fabs(z);                         /* term peak ~ |z|/2      */
    double t;
    int k0 = 0;
    double c = nu + 1.5;
    if (c <= 0.0 && c == floor(c)) {
        /* nu = -(n+1/2), n >= 1: Gamma(k+nu+3/2) = Gamma(k+1-n) has poles at
         * k = 0..n-1, so those terms are 1/Gamma = 0 and the series starts at
         * k0 = n with Gamma(1) = 1:  t_n = (-w)^n / Gamma(n+3/2), built as a
         * product (no overflowing w^n or Gamma). The resulting sum is
         * H_{-(n+1/2)} = (-1)^n J_{n+1/2} (DLMF 11.4.4). */
        k0 = (int)(-(nu + 0.5) + 0.5);
        if (k0 > 4000) return false;
        t = 2.0 / sqrt(M_PI);                           /* 1/Gamma(3/2)           */
        for (int j = 1; j <= k0; j++) t *= ws / (j + 0.5);
    } else {
        double g = tgamma(c);
        if (!isfinite(g) || g == 0.0) return false;     /* overflow / underflow   */
        t = (2.0 / sqrt(M_PI)) / g;                     /* t_0 = 1/(G(3/2)G(nu+3/2)) */
    }
    double sum = t, maxabs = fabs(t);
    int k;
    for (k = k0 + 1; k < 8000; k++) {
        double denom = (k + 0.5) * (k + nu + 0.5);
        if (denom == 0.0) return false;                 /* nu a negative half-int */
        t *= ws / denom;
        sum += t;
        double at = fabs(t);
        if (at > maxabs) maxabs = at;
        if (at <= fabs(sum) * 1.0e-18 && (double)k > hump) break;   /* converged */
    }
    if (k >= 8000) return false;                        /* did not converge       */
    double pref = pow(0.5 * z, nu + 1.0);               /* NaN if z<0 & nu+1 non-int */
    double v = pref * sum;
    if (!isfinite(v)) return false;
    *val = v;
    *relerr = (fabs(sum) > 0.0) ? (maxabs / fabs(sum)) * DBL_EPSILON : 1.0;
    return true;
}

/* The large-|z| asymptotic H_nu(z) = Y_nu(z) + K_nu(z), z > 0.
 *   Y_nu : Hankel asymptotic (DLMF 10.17.4), omega = z - (nu/2 + 1/4) pi.
 *   K_nu : H_nu - Y_nu asymptotic (DLMF 11.6.1); for a half-integer nu >= 1/2 the
 *          (nu+1/2-k) factor hits zero and the series terminates exactly.
 * Both are summed to optimal truncation (smallest term). Sets *val and the
 * predicted relative error *relerr (combined truncation floor / |H|). */
static bool struve_asymp(double nu, double z, double* val, double* relerr) {
    if (!(z > 0.0)) return false;

    double mu = 4.0 * nu * nu;
    double omega = z - (0.5 * nu + 0.25) * M_PI;
    double a = 1.0, invzk = 1.0, P = 1.0, Q = 0.0;      /* k=0 adds a_0 to P      */
    double yfloor = 1.0, prevy = 1.0;
    /* Half-integer nu: the factor mu - (2k-1)^2 vanishes at k = |nu| + 1/2, so the
     * Hankel series is a FINITE, exact sum (Y_{n+1/2} is elementary). Its terms
     * may grow before they stop (nu = -13/2 at z = 20), which is not divergence:
     * sum them all and charge only the rounding, eps * sum |a_k / z^k|. Cutting
     * at the growth instead left an 8e-10 floor where the elementary form is
     * exact to the last bit. */
    double twonu = 2.0 * nu;
    bool yterm = fabs(nu) < 390.0 && twonu == floor(twonu) && fmod(fabs(twonu), 2.0) == 1.0;
    double yabs = 1.0;
    for (int k = 1; k < 400; k++) {
        a *= (mu - (double)(2 * k - 1) * (2 * k - 1)) / ((double)k * 8.0);
        invzk /= z;
        double tk = a * invzk, atk = fabs(tk);          /* a_k / z^k              */
        if (yterm) {
            if (tk == 0.0) break;                       /* terminated: exact      */
            if (!isfinite(tk)) return false;
            yabs += atk;
        } else {
            if (atk > prevy) break;                     /* divergent: smallest term */
            prevy = atk; yfloor = atk;
        }
        int m = k / 2;
        double sgn = (m % 2 == 0) ? 1.0 : -1.0;         /* (-1)^floor(k/2)        */
        if (k % 2 == 0) P += sgn * tk;                  /* even k -> P            */
        else            Q += sgn * tk;                  /* odd k  -> Q            */
    }
    if (yterm) yfloor = 2.0 * DBL_EPSILON * yabs;      /* rounding of a finite sum */
    double scale = sqrt(2.0 / (M_PI * z));
    double Y = scale * (sin(omega) * P + cos(omega) * Q);

    double w = 0.25 * z * z;
    double c = nu + 0.5;
    double Ksum = 0.0, kfloor = 0.0;
    /* K_nu = (1/pi) sum_k Gamma(k+1/2) (z/2)^(nu-2k-1) / Gamma(nu+1/2-k): for
     * nu+1/2 = 0, -1, ... EVERY 1/Gamma(nu+1/2-k) vanishes, so K_nu == 0 exactly
     * and H_nu = Y_nu (= (-1)^n J_{n+1/2} for nu = -(n+1/2)). */
    if (!(c <= 0.0 && c == floor(c))) {
        double gk = tgamma(c);
        double term = sqrt(M_PI) * pow(0.5 * z, nu - 1.0) / gk;
        if (!isfinite(term)) return false;                  /* overflow               */
        double prevk = fabs(term), kabs = fabs(term);
        Ksum = term; kfloor = fabs(term);
        for (int k = 1; k < 400; k++) {
            term *= (k - 0.5) * (nu + 0.5 - k) / w;
            if (term == 0.0) {                              /* terminated (half-int nu) */
                kfloor = 2.0 * DBL_EPSILON * kabs;          /* exact up to rounding     */
                break;
            }
            double at = fabs(term);
            if (!yterm && at > prevk) break;                /* divergent: smallest term */
            prevk = at; kfloor = at; kabs += at;
            Ksum += term;
        }
    }
    double K = Ksum / M_PI;

    double H = Y + K;
    if (!isfinite(H)) return false;
    double abserr = scale * yfloor + kfloor / M_PI;     /* each floor in its scale */
    *val = H;
    *relerr = (fabs(H) > 0.0) ? abserr / fabs(H) : abserr;
    return true;
}

bool sf_machine_struve_h(double nu, double z, double* out) {
    if (z == 0.0) { *out = 0.0; return nu > -1.0; }

    double cv, ce; bool cok = struve_series(nu, z, -1, &cv, &ce);
    double av, ae; bool aok = struve_asymp(nu, z, &av, &ae);

    /* Prefer the regime with the smaller predicted error; accept only if that
     * error clears the budget (otherwise the arb bridge answers). */
    if (cok && (!aok || ce <= ae)) { *out = cv; return isfinite(cv) && ce < STRUVE_ACCEPT; }
    if (aok)                       { *out = av; return isfinite(av) && ae < STRUVE_ACCEPT; }
    return false;
}

/* ------------------------------------------------------------------ */
/*  StruveL L_nu(z) -- real machine path                               */
/* ------------------------------------------------------------------ */
/* L_nu(z) = (z/2)^(nu+1) sum_{k>=0} (z/2)^(2k)/(Gamma(k+3/2)Gamma(k+nu+3/2))
 * (DLMF 11.2.2): the StruveH series with every term of one sign once
 * Gamma(k+nu+3/2) > 0, so for nu > -3/2 it is cancellation-free at every real
 * z and the convergent series alone is accurate to a few ulp (the summation's
 * own rounding). For large z it switches to L_nu = I_{-nu} + the DLMF 11.6.2
 * asymptotic (struvel_asymp), which needs a handful of terms instead of ~z/2.
 * Declines where the value is complex (z < 0 with non-integer nu), the series
 * lost more than STRUVEL_ACCEPT to cancellation (nu < -3/2), Gamma overflow
 * (large nu), or the result overflows a double; the caller then answers through
 * the arb bridge. Shared bit for bit with the NDArray kernel and Compile VM. */
#define STRUVEL_ACCEPT 1.0e-14
#define STRUVEL_ASYMP_Z 30.0   /* large-z switch, see struvel_asymp */

/* L_nu(z) = I_{-nu}(z) + (1/pi) sum_{k>=0} (-1)^(k+1) Gamma(k+1/2) (z/2)^(nu-2k-1)
 *                                        / Gamma(nu+1/2-k)          (DLMF 11.6.2)
 * for z > STRUVEL_ASYMP_Z. I_{-nu} grows like e^z while the correction is
 * algebraic, so the (divergent) sum only has to reach eps relative to the
 * total; it terminates exactly at a half-integer nu. Declines (false) when the
 * smallest term is still above that floor (|nu| comparable to z), so the
 * convergent series answers instead. */
static bool struvel_asymp(double nu, double z, double* out) {
    double bi;
    if (!sf_bessel_i_r(-nu, z, &bi) || !isfinite(bi) || !(bi > 0.0)) return false;
    double c = nu + 0.5, sum = 0.0;
    if (!(c <= 0.0 && c == floor(c))) {           /* else every 1/Gamma is 0 */
        double g = tgamma(c);
        if (!isfinite(g) || g == 0.0) return false;
        double term = -sqrt(M_PI) * pow(0.5 * z, nu - 1.0) / g;
        if (!isfinite(term)) return false;
        double w = 0.25 * z * z, prev = fabs(term);
        sum = term;
        for (int k = 1; ; k++) {
            if (k > 200) return false;
            term *= -(k - 0.5) * (nu + 0.5 - k) / w;
            if (term == 0.0) break;                  /* terminated (half-int nu) */
            double at = fabs(term);
            if (at > prev) return false;             /* divergence before eps    */
            sum += term;
            if (at < 1.0e-17 * M_PI * bi) break;      /* below the rounding floor */
            prev = at;
        }
    }
    double v = bi + sum / M_PI;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool sf_machine_struve_l(double nu, double z, double* out) {
    if (z == 0.0) { *out = 0.0; return nu > -1.0; }
    if (!isfinite(nu) || !isfinite(z)) return false;
    double v, e;
    if (z > STRUVEL_ASYMP_Z && struvel_asymp(nu, z, &v)) { *out = v; return true; }
    if (!struve_series(nu, z, +1, &v, &e)) return false;
    *out = v;
    return isfinite(v) && e < STRUVEL_ACCEPT;
}
