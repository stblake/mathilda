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
#ifndef M_E
#define M_E 2.71828182845904523536
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
    /* zeta(0) = -1/2.  The functional equation cannot reach it: sin(pi s/2) is
     * zero there while zeta(1-s) has its pole, so the product is 0 * infinity.
     * Treating that zero sine as a trivial zero returned 0 instead of -1/2 —
     * which then propagated into PolyLog, whose expansion sums zeta(n-k) and so
     * hits s = 0 for every integer order. */
    if (s == 0.0) { *out = -0.5; return true; }
    /* The trivial zeros are the NEGATIVE EVEN integers, and only those. */
    if (s < 0.0 && s == floor(s) && fmod(s, 2.0) == 0.0) { *out = 0.0; return true; }
    if (s < 0.5) {
        /* zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s) */
        double sn = sin(M_PI * s / 2.0);
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

/* Sign of Gamma(x), computed rather than taken from lgamma_r or the global
 * signgam: lgamma_r is POSIX rather than C99 and signgam is not thread-safe.
 * Gamma alternates sign on the negative axis, negative on (-1,0), positive on
 * (-2,-1), and so on — so the parity of floor(x) decides it. */
static double gamma_sign(double x) {
    if (x > 0.0) return 1.0;
    long k = (long)floor(x);
    return (k & 1L) ? -1.0 : 1.0;
}

/* Pochhammer[a, n] = Gamma(a+n)/Gamma(a), evaluated through lgamma with the
 * sign tracked separately so the poles of Gamma cancel where the ratio is
 * finite (a negative integer a with a+n still negative, say). */
bool sf_machine_pochhammer(double a, double n, double* out) {
    if (n == 0.0) { *out = 1.0; return true; }
    /* Small non-negative integer n: the product is exact and pole-free. */
    if (n == floor(n) && n > 0.0 && n <= 200.0) {
        double v = 1.0;
        for (int k = 0; k < (int)n; k++) v *= a + (double)k;
        *out = v;
        return isfinite(v);
    }
    if (a <= 0.0 && a == floor(a)) return false;          /* pole in Gamma(a) */
    if (a + n <= 0.0 && a + n == floor(a + n)) return false;
    double l1 = lgamma(a + n), l2 = lgamma(a);
    double v = gamma_sign(a + n) * gamma_sign(a) * exp(l1 - l2);
    *out = v;
    return isfinite(v);
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

/* 0F1(;a;z) = sum z^k / ((a)_k k!) — entire, so the series is the whole story. */
bool sf_machine_pfq(const double* a, size_t p, const double* b, size_t q,
                    double z, double* out);   /* defined below */

bool sf_machine_hyper0f1(double a, double z, double* out) {
    return sf_machine_pfq(NULL, 0, &a, 1, z, out);
}

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
    *out = p;
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

/* Li_s(x) on the reals, |x| <= 1.
 *
 * The defining series only converges usefully for |x| <~ 1/2, so closer to 1 the
 * Jonquiere expansion in mu = log x takes over; for negative x the duplication
 * identity folds the problem back onto positive arguments.  Beyond |x| = 1 the
 * function is complex and the kernel declines. */
static bool polylog_mu(double s, double x, double* out) {
    /* Li_s(x) = Gamma(1-s) (-mu)^(s-1) + sum_k zeta(s-k) mu^k / k!,  mu = log x */
    double mu = log(x);
    if (s == floor(s) && s >= 1.0) {
        /* Integer s: the Gamma pole and the zeta pole cancel; the standard
         * replacement for the k = s-1 term is mu^(n-1)/(n-1)! (H_{n-1} - log(-mu)). */
        int n = (int)s;
        double sum = 0.0, mk = 1.0;
        int small = 0;
        for (int k = 0; k < 200; k++) {
            if (k != n - 1) {
                double z;
                if (!sf_machine_zeta((double)(n - k), &z)) return false;
                double add = z * mk / tgamma((double)k + 1.0);
                sum += add;
                /* TWO consecutive negligible terms, not one: zeta vanishes at
                 * every negative even integer, so this series has a zero term
                 * every other step and a single-term test stops it early. */
                if (k > n && fabs(add) < fabs(sum) * EI_EPS) { if (++small >= 2) break; }
                else small = 0;
            }
            mk *= mu;
        }
        double H = 0.0;
        for (int j = 1; j <= n - 1; j++) H += 1.0 / j;
        double lm = (mu < 0.0) ? log(-mu) : log(mu);
        sum += pow(mu, (double)(n - 1)) / tgamma((double)n) * (H - lm);
        *out = sum;
        return isfinite(sum);
    }
    double sum = tgamma(1.0 - s) * pow(-mu, s - 1.0), mk = 1.0;
    int small = 0;
    for (int k = 0; k < 200; k++) {
        double z;
        if (!sf_machine_zeta(s - k, &z)) return false;
        double add = z * mk / tgamma((double)k + 1.0);
        sum += add;
        if (k > 2 && fabs(add) < fabs(sum) * EI_EPS) { if (++small >= 2) break; }
        else small = 0;
        mk *= mu;
    }
    *out = sum;
    return isfinite(sum);
}

bool sf_machine_polylog(double s, double x, double* out) {
    if (fabs(x) > 1.0) return false;                  /* complex past the cut */
    if (x == 0.0) { *out = 0.0; return true; }
    if (x == 1.0) { if (!(s > 1.0)) return false; return sf_machine_zeta(s, out); }
    if (x == -1.0) {
        double z;
        if (!sf_machine_zeta(s, &z)) return false;
        *out = -(1.0 - pow(2.0, 1.0 - s)) * z;
        return isfinite(*out);
    }
    if (fabs(x) <= 0.5) {                             /* the defining series */
        double sum = 0.0, xk = x;
        for (int k = 1; k < 2000; k++) {
            double add = xk / pow((double)k, s);
            sum += add;
            if (fabs(add) < fabs(sum) * EI_EPS) break;
            xk *= x;
        }
        *out = sum;
        return isfinite(sum);
    }
    if (x > 0.0) return polylog_mu(s, x, out);
    /* Duplication: Li_s(-y) = 2^(1-s) Li_s(y^2) - Li_s(y). */
    double y = -x, a, b;
    if (!sf_machine_polylog(s, y * y, &a)) return false;
    if (!sf_machine_polylog(s, y, &b)) return false;
    double v = pow(2.0, 1.0 - s) * a - b;
    *out = v;
    return isfinite(v);
}

/* ---- n-ary ------------------------------------------------------------- */

/* LerchPhi[z, s, a] = sum_{k>=0} z^k / (k+a)^s.
 * The defining series converges for |z| < 1 and is used directly; at z = 1 the
 * function IS the Hurwitz zeta, which has a proper expansion above. */
bool sf_machine_lerchphi(double z, double s, double a, double* out) {
    if (a <= 0.0 && a == floor(a)) return false;      /* term k = -a is a pole */
    if (z == 1.0) return sf_machine_hurwitz_zeta(s, a, out);
    if (fabs(z) >= 1.0) return false;
    double sum = 0.0, zk = 1.0;
    for (int k = 0; k < 200000; k++) {
        double d = a + k;
        if (d == 0.0) return false;
        double add = zk / pow(d, s);
        sum += add;
        if (fabs(add) < fabs(sum) * EI_EPS && k > 2) break;
        zk *= z;
        if (zk == 0.0) break;
    }
    *out = sum;
    return isfinite(sum);
}

/* 1F1(a; b; z).  For z very negative the direct series alternates and loses
 * everything to cancellation, so Kummer's transformation moves it to +z. */
bool sf_machine_hyper1f1(double a, double b, double z, double* out) {
    /* One implementation, not two: pFq is what the evaluator actually reaches
     * (Hypergeometric1F1 canonicalises to HypergeometricPFQ), and keeping a
     * second copy here meant a fix to one silently missed the other. */
    return sf_machine_pfq(&a, 1, &b, 1, z, out);
}

/* 2F1(a, b; c; z), |z| < 1 by the defining series.  Outside the disc the
 * function needs a connection formula and may be complex, so the kernel
 * declines and the MPFR path answers. */
bool sf_machine_hyper2f1(double a, double b, double cc, double z, double* out) {
    double ab[2] = { a, b };
    return sf_machine_pfq(ab, 2, &cc, 1, z, out);
}

/* pFq({a1..ap}, {b1..bq}, z) = sum_k prod (a_i)_k / prod (b_j)_k * z^k / k!.
 *
 * This is the head that actually matters: the evaluator canonicalises
 * Hypergeometric0F1, 1F1 and 2F1 all into HypergeometricPFQ before anything
 * downstream sees them.
 *
 * Convergence is decided by p against q: p <= q is entire, p = q+1 needs
 * |z| < 1, and p > q+1 diverges — the last is declined rather than truncated. */
bool sf_machine_pfq(const double* a, size_t p, const double* b, size_t q,
                    double z, double* out) {
    for (size_t j = 0; j < q; j++)
        if (b[j] <= 0.0 && b[j] == floor(b[j])) return false;   /* (b)_k hits zero */
    if (p > q + 1) return false;
    if (p == q + 1 && !(fabs(z) < 1.0)) return false;

    /* 1F1 at ANY negative z: the raw series alternates with terms larger than
     * the sum, so Kummer's transformation moves it to +z where every term is
     * positive and nothing cancels.  The condition b - a > 0 is what guarantees
     * the transformed series is positive-term; without the transform the answer
     * was wrong in the second decimal place at z = -40, and still 5e-9 at -20. */
    if (p == 1 && q == 1 && z < 0.0 && (b[0] - a[0]) > 0.0) {
        double ba = b[0] - a[0], v;
        if (!sf_machine_pfq(&ba, 1, b, 1, -z, &v)) return false;
        *out = exp(z) * v;
        return isfinite(*out);
    }

    double maxterm = 1.0;
    double term = 1.0, sum = 1.0;
    for (int k = 0; k < 200000; k++) {
        double num = 1.0, den = 1.0;
        for (size_t i = 0; i < p; i++) num *= (a[i] + k);
        for (size_t j = 0; j < q; j++) den *= (b[j] + k);
        if (den == 0.0) return false;
        term *= num * z / (den * (k + 1.0));
        if (term == 0.0) break;                 /* an upper parameter terminated it */
        sum += term;
        if (fabs(term) > maxterm) maxterm = fabs(term);
        if (!isfinite(term)) return false;
        if (fabs(term) < fabs(sum) * EI_EPS) break;
    }
    /* Cancellation guard: if the largest term dwarfed the result, the sum has
     * lost most of its digits and the kernel must DECLINE rather than hand back
     * a confidently wrong number — the MPFR path carries enough precision. */
    if (maxterm > fabs(sum) * 1.0e13) return false;
    *out = sum;
    return isfinite(sum);
}

/* The ascending series above converges everywhere, but it is not USABLE
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
