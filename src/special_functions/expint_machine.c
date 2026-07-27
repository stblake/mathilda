/* Machine-precision (double) kernels for the exponential-integral family.
 * See expint_machine.h for why these exist separately from the MPFR modules.
 *
 * Each function uses the standard split: a power series where it converges
 * quickly and is well conditioned, a continued fraction or asymptotic form
 * beyond.  Accuracy is verified against the modules' own MPFR implementations
 * by tests/test_expint_machine.c rather than asserted here.
 */
#include "expint_machine.h"

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
