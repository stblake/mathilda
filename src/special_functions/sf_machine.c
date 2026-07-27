/* Machine-precision (double) kernels for the exponential-integral family.
 * See sf_machine.h for why these exist separately from the MPFR modules.
 *
 * Each function uses the standard split: a power series where it converges
 * quickly and is well conditioned, a continued fraction or asymptotic form
 * beyond.  Accuracy is verified against the modules' own MPFR implementations
 * by tests/test_sf_machine.c rather than asserted here.
 */
#include "sf_machine.h"

#include <math.h>
#include <complex.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EI_EULER   0.57721566490153286061
#define EI_MAXIT   200
#define EI_EPS     1.0e-16
#define EI_FPMIN   1.0e-300

/* E_1(x) for x > 0.  Series below 1, modified-Lentz continued fraction above:
 *   E_1(x) = e^-x / (x + 1 - 1^2/(x + 3 - 2^2/(x + 5 - ...)))
 * which is the form that stays stable as x grows. */
bool sf_machine_e1(double x, double* out) {
    if (!(x > 0.0)) return false;                     /* pole at 0, complex below */
    if (x <= 1.0) {
        double sum = 0.0, fact = 1.0;
        for (int k = 1; k <= EI_MAXIT; k++) {
            fact *= -x / k;
            double term = -fact / k;
            sum += term;
            if (fabs(term) < fabs(sum) * EI_EPS) break;
        }
        double v = sum - log(x) - EI_EULER;
        *out = v;
        return isfinite(v);
    }
    double b = x + 1.0, c = 1.0 / EI_FPMIN, d = 1.0 / b, h = d;
    for (int i = 1; i <= EI_MAXIT; i++) {
        double a = -(double)i * (double)i;
        b += 2.0;
        d = 1.0 / (a * d + b);
        c = b + a / c;
        double del = c * d;
        h *= del;
        if (fabs(del - 1.0) < EI_EPS) break;
    }
    double v = h * exp(-x);
    *out = v;
    return isfinite(v);
}

/* Ei(x).  Below zero it is the principal value, which is exactly -E_1(-x); above
 * zero the power series up to where its terms stop being well scaled, then the
 * asymptotic series. */
bool sf_machine_ei(double x, double* out) {
    if (x == 0.0) return false;                        /* -Infinity */
    if (x < 0.0) {
        double e1;
        if (!sf_machine_e1(-x, &e1)) return false;
        *out = -e1;
        return isfinite(*out);
    }
    if (x < 40.0) {
        double sum = 0.0, fact = 1.0;
        for (int k = 1; k <= EI_MAXIT; k++) {
            fact *= x / k;
            double term = fact / k;
            sum += term;
            if (term < sum * EI_EPS) break;
        }
        double v = sum + log(x) + EI_EULER;
        *out = v;
        return isfinite(v);
    }
    /* Asymptotic: Ei(x) ~ e^x/x * sum k!/x^k, truncated at its smallest term. */
    double sum = 0.0, term = 1.0;
    for (int k = 1; k <= EI_MAXIT; k++) {
        double prev = term;
        term *= (double)k / x;
        if (term < EI_EPS) { sum += term; break; }
        if (term < prev) sum += term;
        else { sum -= prev; break; }                   /* series has turned */
    }
    double v = exp(x) * (1.0 + sum) / x;
    *out = v;
    return isfinite(v);
}

/* Si and Ci together: the power series for small |x|, and for larger |x| the
 * continued fraction for the complex exponential integral, whose real and
 * imaginary parts carry Ci and Si at once. */
static bool cisi_double(double x, double* si, double* ci) {
    const double TMIN = 2.0;
    double t = fabs(x);
    if (t == 0.0) { *si = 0.0; *ci = -HUGE_VAL; return false; }   /* Ci(0) = -Inf */

    if (t > TMIN) {
        double _Complex b = 1.0 + t * I;
        double _Complex c = 1.0 / EI_FPMIN;
        double _Complex d = 1.0 / b;
        double _Complex h = d;
        for (int i = 2; i <= EI_MAXIT; i++) {
            double a = -(double)(i - 1) * (double)(i - 1);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            double _Complex del = c * d;
            h *= del;
            if (fabs(creal(del) - 1.0) + fabs(cimag(del)) < EI_EPS) break;
        }
        h = (cos(t) - sin(t) * I) * h;
        *ci = -creal(h);
        *si = M_PI / 2.0 + cimag(h);
    } else {
        double sum = 0.0, sums = 0.0, sumc = 0.0, sign = 1.0, fact = 1.0;
        int odd = 1;
        for (int k = 1; k <= EI_MAXIT; k++) {
            fact *= t / k;
            double term = fact / k;
            sum += sign * term;
            double err = term / fabs(sum);
            if (odd) { sign = -sign; sums = sum; sum = sumc; }
            else      { sumc = sum; sum = sums; }
            if (err < EI_EPS) break;
            odd = !odd;
        }
        *si = sums;
        *ci = sumc + log(t) + EI_EULER;
    }
    if (x < 0.0) *si = -(*si);
    return true;
}

bool sf_machine_si(double x, double* out) {
    double si, ci;
    if (!cisi_double(x, &si, &ci)) { if (x == 0.0) { *out = 0.0; return true; } return false; }
    *out = si;
    return isfinite(si);
}

bool sf_machine_ci(double x, double* out) {
    if (x <= 0.0) return false;               /* complex for x < 0, -Inf at 0 */
    double si, ci;
    if (!cisi_double(x, &si, &ci)) return false;
    *out = ci;
    return isfinite(ci);
}

/* Shi and Chi from Ei: Chi + Shi = Ei(x) and Chi - Shi = Ei(-x), so each is a
 * half-sum.  Near zero that difference cancels the shared `EulerGamma + log x`,
 * so small |x| takes the series directly instead. */
bool sf_machine_shi(double x, double* out) {
    double t = fabs(x);
    if (t < 0.5) {
        double sum = 0.0, term = t;
        for (int k = 0; k <= EI_MAXIT; k++) {
            double add = term / (2.0 * k + 1.0);
            sum += add;
            if (fabs(add) < fabs(sum) * EI_EPS) break;
            term *= t * t / ((2.0 * k + 2.0) * (2.0 * k + 3.0));
        }
        *out = (x < 0.0) ? -sum : sum;
        return isfinite(*out);
    }
    double ep, em;
    if (!sf_machine_ei(t, &ep) || !sf_machine_ei(-t, &em)) return false;
    double v = 0.5 * (ep - em);
    *out = (x < 0.0) ? -v : v;
    return isfinite(*out);
}

bool sf_machine_chi(double x, double* out) {
    if (x <= 0.0) return false;               /* complex for x < 0, -Inf at 0 */
    double ep, em;
    if (!sf_machine_ei(x, &ep) || !sf_machine_ei(-x, &em)) return false;
    double v = 0.5 * (ep + em);
    *out = v;
    return isfinite(v);
}

/* li(x) = Ei(log x): real on [0, 1) and (1, inf), with a pole at 1. */
bool sf_machine_li(double x, double* out) {
    /* li(0) is 0 by the limit, but the interpreter answers Indeterminate there,
     * so the machine path must decline rather than be the only one to give a
     * number. */
    if (!(x > 0.0)) return false;
    if (x == 1.0) return false;               /* -Infinity */
    return sf_machine_ei(log(x), out);
}

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

/* Lambert W, principal branch.  Halley's iteration converges cubically, so a
 * decent starting guess reaches machine precision in a handful of steps. */
bool sf_machine_productlog(double x, double* out) {
    const double INV_E = 0.36787944117144232160;
    if (x < -INV_E) return false;                     /* complex below the branch point */
    if (x == 0.0) { *out = 0.0; return true; }
    if (x == -INV_E) { *out = -1.0; return true; }

    double w;
    if (x < -0.3) {                                   /* series about the branch point */
        double p = sqrt(2.0 * (M_E * x + 1.0));
        w = -1.0 + p - p * p / 3.0 + 11.0 * p * p * p / 72.0;
    } else if (x < 1.0) {
        w = x * (1.0 - x + 1.5 * x * x);
    } else {
        double l1 = log(x), l2 = log(l1);
        w = l1 - l2 + l2 / l1;
    }
    for (int i = 0; i < 60; i++) {
        double e = exp(w), we = w * e, f = we - x;
        double denom = e * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0);
        if (denom == 0.0) break;
        double dw = f / denom;
        w -= dw;
        if (fabs(dw) <= EI_EPS * (fabs(w) + 1.0)) break;
    }
    *out = w;
    return isfinite(w);
}

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
        double pix2 = M_PI * ax * ax;
        double _Complex b = 1.0 - pix2 * I;
        double _Complex cc = 1.0 / EI_FPMIN;
        double _Complex d = 1.0 / b;
        double _Complex h = d;
        int n = -1;
        for (int k = 2; k <= 400; k++) {
            n += 2;
            double a = -(double)n * (double)(n + 1);
            b += 4.0;
            d = 1.0 / (a * d + b);
            cc = b + a / cc;
            double _Complex del = cc * d;
            h *= del;
            if (fabs(creal(del) - 1.0) + fabs(cimag(del)) < EI_EPS) break;
        }
        h = (ax - ax * I) * h;
        double _Complex cs = (0.5 + 0.5 * I) * (1.0 - (cos(0.5 * pix2) + sin(0.5 * pix2) * I) * h);
        *c = creal(cs); *s = cimag(cs);
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
bool sf_machine_digamma(double x, double* out) {
    if (x <= 0.0 && x == floor(x)) return false;      /* pole */
    double v = 0.0, z = x;
    if (z < 0.0) {
        /* psi(z) = psi(1-z) - pi/tan(pi z) */
        double refl = M_PI / tan(M_PI * z);
        double up;
        if (!sf_machine_digamma(1.0 - z, &up)) return false;
        *out = up - refl;
        return isfinite(*out);
    }
    while (z < 12.0) { v -= 1.0 / z; z += 1.0; }
    double inv = 1.0 / z, inv2 = inv * inv;
    /* psi(z) ~ log z - 1/(2z) - sum B_2k / (2k z^2k) */
    v += log(z) - 0.5 * inv;
    double p = inv2;
    v -= p / 12.0;              p *= inv2;
    v += p / 120.0;             p *= inv2;
    v -= p / 252.0;             p *= inv2;
    v += p / 240.0;             p *= inv2;
    v -= p / 132.0;             p *= inv2;
    v += p * 691.0 / 32760.0;   p *= inv2;
    v -= p / 12.0;
    *out = v;
    return isfinite(v);
}

bool sf_machine_harmonic(double x, double* out) {
    double d;
    if (!sf_machine_digamma(x + 1.0, &d)) return false;
    *out = d + EI_EULER;
    return isfinite(*out);
}

/* Riemann zeta on the reals.  Euler-Maclaurin above the critical strip, the
 * functional equation below it. */
bool sf_machine_zeta(double s, double* out) {
    if (s == 1.0) return false;                       /* pole */
    if (s < 0.5) {
        /* zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s) */
        double sn = sin(M_PI * s / 2.0);
        if (sn == 0.0) { *out = 0.0; return true; }   /* trivial zeros */
        double zr;
        if (!sf_machine_zeta(1.0 - s, &zr)) return false;
        double v = pow(2.0, s) * pow(M_PI, s - 1.0) * sn * tgamma(1.0 - s) * zr;
        *out = v;
        return isfinite(v);
    }
    const int N = 16;
    double sum = 0.0;
    for (int n = 1; n < N; n++) sum += pow((double)n, -s);
    double Ns = pow((double)N, -s);
    sum += Ns * (double)N / (s - 1.0) + 0.5 * Ns;
    /* + sum_k B_2k/(2k)! * (s)_{2k-1} * N^(-s-2k+1) */
    static const double B2K_OVER_FACT[] = {          /* B_2k / (2k)! */
        1.0 / 12.0, -1.0 / 720.0, 1.0 / 30240.0, -1.0 / 1209600.0,
        1.0 / 47900160.0, -691.0 / 1307674368000.0, 1.0 / 74724249600.0,
        -3617.0 / 10670622842880000.0
    };
    double rising = s, Npow = Ns / (double)N;         /* N^(-s-1) */
    double N2 = 1.0 / ((double)N * (double)N);
    for (int k = 1; k <= 8; k++) {
        sum += B2K_OVER_FACT[k - 1] * rising * Npow;
        /* advance the rising factorial by two and N^(-s-2k+1) by N^-2 */
        rising *= (s + 2.0 * k - 1.0) * (s + 2.0 * k);
        Npow *= N2;
    }
    *out = sum;
    return isfinite(sum);
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

/* Hurwitz zeta on the reals, by the same Euler-Maclaurin expansion as the
 * Riemann zeta: sum the first N terms directly, then the tail analytically.
 * Restricted to s > 1 and a > 0, where the series is what it is defined by;
 * the analytic continuation below s = 1 is left to the MPFR path. */
bool sf_machine_hurwitz_zeta(double s, double a, double* out) {
    if (!(s > 1.0) || !(a > 0.0)) return false;
    const int N = 16;
    double sum = 0.0;
    for (int k = 0; k < N; k++) sum += pow((double)k + a, -s);
    double Na = (double)N + a;
    double Ns = pow(Na, -s);
    sum += Ns * Na / (s - 1.0) + 0.5 * Ns;
    static const double B2K_OVER_FACT[] = {
        1.0 / 12.0, -1.0 / 720.0, 1.0 / 30240.0, -1.0 / 1209600.0,
        1.0 / 47900160.0, -691.0 / 1307674368000.0, 1.0 / 74724249600.0,
        -3617.0 / 10670622842880000.0
    };
    double rising = s, Npow = Ns / Na, Na2 = 1.0 / (Na * Na);
    for (int k = 1; k <= 8; k++) {
        sum += B2K_OVER_FACT[k - 1] * rising * Npow;
        rising *= (s + 2.0 * k - 1.0) * (s + 2.0 * k);
        Npow *= Na2;
    }
    *out = sum;
    return isfinite(sum);
}

/* PolyGamma[n, x].  n = 0 is the digamma above; for n >= 1 the standard
 * identity psi^(n)(x) = (-1)^(n+1) n! zeta(n+1, x) reuses the Hurwitz zeta. */
bool sf_machine_polygamma(double n, double x, double* out) {
    if (n < 0.0 || n != floor(n) || n > 60.0) return false;
    if (n == 0.0) return sf_machine_digamma(x, out);
    if (!(x > 0.0)) return false;                 /* reflection for n >= 1 not covered */
    double hz;
    if (!sf_machine_hurwitz_zeta(n + 1.0, x, &hz)) return false;
    double v = tgamma(n + 1.0) * hz;
    if (fmod(n, 2.0) == 0.0) v = -v;              /* (-1)^(n+1) */
    *out = v;
    return isfinite(v);
}
