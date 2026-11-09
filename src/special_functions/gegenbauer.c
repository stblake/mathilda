/* Mathilda -- GegenbauerC, the Gegenbauer (ultraspherical) polynomials and their
 * analytic continuation in the degree.
 *
 *   GegenbauerC[n, lambda, z] = C_n^lambda(z)
 *       = Gamma(n + 2 lambda) / (Gamma(2 lambda) Gamma(n + 1))
 *           * 2F1(-n, n + 2 lambda; lambda + 1/2; (1 - z)/2)          (DLMF 18.5.9)
 *   GegenbauerC[n, z]         = 2 ChebyshevT[n, z] / n   (the lambda -> 0 limit
 *                               of C_n^lambda / lambda, the reference convention)
 *
 * For a non-negative integer n the polynomial is the explicit sum
 *   C_n^lambda(z) = sum_k (-1)^k (lambda)_{n-k} (2z)^{n-2k} / (k! (n-2k)!)
 * for EVERY lambda (finite where the 2F1 form is 0 * infinity: lambda = -1/2, ...),
 * and C_n^0 = 0 -- the reference system's conventions (GegenbauerC[2, -1/2, x] =
 * 1/2 - x^2/2, GegenbauerC[n, 0, x] = 0).
 *
 * Exact rules (three arguments):
 *   lambda = 0 -> 0;  n = 0 -> 1;  lambda = 1 -> ChebyshevU[n, z];
 *   lambda = 1/2 -> LegendreP[n, z];  integer 0 <= n <= GEG_EXPAND_MAX -> the
 *   explicit polynomial (any lambda, z);  anything else exact stays unevaluated.
 * The two-argument form always rewrites to 2 ChebyshevT[n, z]/n (n = 0 gives
 * ComplexInfinity, with Power::infy, as in the reference system).
 *
 * Numerics (any inexact argument):
 *   machine -> sf_machine_gegenbauer[_c] / sf_machine_gegenbauer2_c, shared with
 *              the NDArray / packed map and Compile[]: the three-term recurrence
 *              for an integer degree (real or complex z), the machine 2F1 kernel
 *              with the Gamma prefactor otherwise, and arb at 53 bits for what
 *              the double algorithms decline (a cancelling last step, complex
 *              parameters, a pole of a Gamma factor, a huge degree).
 *   MPFR    -> flint_num_gegenbauer (arb), at the arguments' precision.
 * D[GegenbauerC[n, l, z], z] = 2 l GegenbauerC[n - 1, l + 1, z] (DLMF 18.9.19);
 * the n- and l-derivatives stay inert.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "gegenbauer.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "attr.h"
#include "common.h"            /* builtin_arg_error */
#include "eval.h"              /* eval_and_free */
#include "expr.h"
#include "flint_num_bridge.h"
#include "sf_machine.h"        /* sf_machine_hyper2f1, sf_machine_chebyshev_t[_c] */
#include "sf_machine_hyper.h"  /* sf_machine_hyper2f1_c */
#include "sf_scalar.h"
#include "sf_dd.h"
#include "symtab.h"
#include "../sym_names.h"

/* Largest integer degree run by the double recurrence (O(n)). */
#define GEG_REC_MAX 100000.0
/* Largest integer degree expanded to an explicit polynomial symbolically. */
#define GEG_EXPAND_MAX 400

/* ------------------------------------------------------------------ */
/* Machine kernels                                                     */
/* ------------------------------------------------------------------ */

/* 1/Gamma(x), exactly 0 at the poles, by reflection below 0 (tgamma of a large
 * negative argument loses ulps). 0 is also returned past the double range. */
static double geg_rgamma(double x) {
    if (x <= 0.0 && x == floor(x)) return 0.0;
    if (x < 0.0) {
        double r = x - 2.0 * floor(0.5 * x + 0.5);        /* r in [-1, 1) */
        if (r > 0.5) r = 1.0 - r;
        else if (r < -0.5) r = -1.0 - r;
        return sin(M_PI * r) * tgamma(1.0 - x) / M_PI;
    }
    return 1.0 / tgamma(x);
}

static bool geg_is_nonneg_int(double n) {
    return n >= 0.0 && n == floor(n) && n <= GEG_REC_MAX;
}

/* The three-term recurrence, real. False when the final step cancels by more
 * than 2^6 (near a zero of C_n the recurrence is only absolutely accurate). */
static bool geg_rec_r(double n, double lam, double x, double* out) {
    long N = (long)n;
    if (N == 0) { *out = 1.0; return true; }
    if (N <= 32) {
        /* plain double, kept only when no step cancels by more than 2: each
         * step then adds ~1 ulp, so a short recurrence is accurate */
        double c0 = 1.0, c1 = 2.0 * lam * x;
        bool clean = true;
        for (long k = 1; k < N && clean; k++) {
            double a = 2.0 * (k + lam) * x * c1;
            double b = (k + 2.0 * lam - 1.0) * c0;
            double c2 = (a - b) / (k + 1.0);
            if (fabs(a - b) * 2.0 < fabs(a) + fabs(b)) clean = false;
            c0 = c1; c1 = c2;
        }
        if (clean && isfinite(c1)) { *out = c1; return true; }
    }
    /* double-double: the double recurrence's rounding grows like ~4 n ulp and is
     * amplified by every cancelling step (5600 ulp at n = 1468, 668 ulp at
     * n = 8 near a zero, measured against mpmath); in double-double, with
     * k + lambda and k + 2 lambda - 1 formed exactly, it is ~1 ulp */
    sdd c0 = sdd_from(1.0), c1 = sdd_mul_d(sdd_from(2.0 * lam), x);
    sdd a = sdd_from(0.0), b = sdd_from(0.0);
    for (long k = 1; k < N; k++) {
        a = sdd_mul(sdd_mul_d(sdd_two_sum((double)k, lam), 2.0 * x), c1);
        b = sdd_mul(sdd_two_sum((double)(k - 1), 2.0 * lam), c0);
        sdd c2 = sdd_div_d(sdd_sub(a, b), k + 1.0);
        c0 = c1; c1 = c2;
        if (!isfinite(c1.hi)) return false;
    }
    double v = c1.hi + c1.lo;
    if (!isfinite(v)) return false;
    /* a cancellation the double-double carries cannot hide: 2^40 */
    if (N >= 2 && fabs(v) * 1e12 < fabs(a.hi) + fabs(b.hi)) return false;
    *out = v;
    return true;
}

/* The same recurrence for complex x (real lambda). */
static bool geg_rec_c(double n, double lam, double xr, double xi, double* rr, double* ri) {
    long N = (long)n;
    if (N == 0) { *rr = 1.0; *ri = 0.0; return true; }
    double c0r = 1.0, c0i = 0.0, c1r = 2.0 * lam * xr, c1i = 2.0 * lam * xi;
    double am = 0.0, bm = 0.0;
    for (long k = 1; k < N; k++) {
        double f = 2.0 * (k + lam);
        double ar = f * (xr * c1r - xi * c1i), ai = f * (xr * c1i + xi * c1r);
        double g = k + 2.0 * lam - 1.0;
        double br = g * c0r, bi = g * c0i;
        double c2r = (ar - br) / (k + 1.0), c2i = (ai - bi) / (k + 1.0);
        c0r = c1r; c0i = c1i; c1r = c2r; c1i = c2i;
        am = hypot(ar, ai); bm = hypot(br, bi);
    }
    if (!isfinite(c1r) || !isfinite(c1i)) return false;
    if (N >= 2 && hypot(c1r, c1i) * (double)N * 64.0 < am + bm) return false;
    *rr = c1r; *ri = c1i;
    return true;
}

/* Gamma(n + 2l) / (Gamma(2l) Gamma(n + 1)), real; false at a pole of the
 * numerator or past the double range. */
static bool geg_prefactor(double n, double lam, double* out) {
    double a = n + 2.0 * lam;
    if (a <= 0.0 && a == floor(a)) return false;
    if (fabs(a) > 170.0 || fabs(2.0 * lam) > 170.0 || fabs(n + 1.0) > 170.0) return false;
    double g = tgamma(a) * geg_rgamma(2.0 * lam) * geg_rgamma(n + 1.0);
    if (!isfinite(g)) return false;
    *out = g;
    return true;
}

/* Normal double range, with a residue imaginary part flushed. */
static bool geg_finish(double vr, double vi, double* rr, double* ri) {
    if (!isfinite(vr) || !isfinite(vi)) return false;
    double big = fmax(fabs(vr), fabs(vi));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(vr) < DBL_MIN) vr = 0.0;
    if (fabs(vi) < DBL_MIN) vi = 0.0;
    *rr = vr; *ri = vi;
    return true;
}

/* 2F1(-n, n + 2l; l + 1/2; t) for 0 <= t <= 0.6 by its power series, the
 * parameters carried exactly as double-doubles (b = n + 2l and c = l + 1/2 round
 * in double, and a rounded parameter is amplified by any cancellation among the
 * terms). Plain double first; double-double when the terms cancel by more than
 * 2^3; declines past 2^50 (the general kernel then answers). ~0.3 us where the
 * general 2F1 kernel, which must also handle every degenerate transformation,
 * costs 5-10 us here. */
static bool geg_2f1_series(double n, double lam, double x, double* out) {
    sdd Td = sdd_two_sum(1.0, -x);                        /* t = (1 - x)/2 exactly */
    Td.hi *= 0.5; Td.lo *= 0.5;
    double t = Td.hi + Td.lo;
    sdd B = sdd_add_d(sdd_two_sum(n, lam), lam);       /* n + 2 lam */
    sdd C = sdd_two_sum(lam, 0.5);                       /* lam + 1/2 */
    double a = -n, b = B.hi + B.lo, c = C.hi + C.lo;
    double term = 1.0, s = 1.0, peak = 1.0;
    int k;
    for (k = 0; k < 600; k++) {
        term *= (a + k) * (b + k) / ((c + k) * (k + 1.0)) * t;
        s += term;
        double at = fabs(term);
        if (at > peak) peak = at;
        if (at < 1e-17 * fabs(s) && a + k > 0.0) break;
    }
    if (k >= 600 || !isfinite(s)) return false;
    if (peak <= 8.0 * fabs(s)) { *out = s; return true; }
    if (peak > 1e15 * fabs(s)) return false;
    sdd T = sdd_from(1.0), S = sdd_from(1.0);
    for (k = 0; k < 600; k++) {
        sdd f = sdd_mul(sdd_two_sum(a, (double)k), sdd_add_d(B, (double)k));
        sdd g = sdd_mul_d(sdd_add_d(C, (double)k), k + 1.0);
        T = sdd_mul(sdd_div(sdd_mul(T, f), g), Td);
        S = sdd_add(S, T);
        if (fabs(T.hi) * peak < 1e-19 * S.hi * S.hi && a + k > 0.0) break;
    }
    if (k >= 600) return false;
    *out = S.hi + S.lo;
    return true;
}

/* -1 < x < -0.2, non-integer degree: the 2F1 at t = (1-x)/2 near 1 converges
 * slowly (the general kernel costs 25-250 us there), so the Gegenbauer ODE
 *     (1 - x^2) y'' - (2l + 1) x y' + n (n + 2l) y = 0
 * is stepped by Taylor series from x = 0, where y = P 2F1(...; 1/2) and
 * y' = 2l C_{n-1}^{l+1}(0) both come from the fast series. Steps stay within
 * half the distance to the singular point -1 and land on exactly representable
 * points. A step whose terms cancel (> 2^6 against the local amplitude) declines
 * to the general kernel. */
static bool geg_ode_neg(double n, double lam, double P, double x, double* out) {
    double F0, F1, P1;
    if (!geg_2f1_series(n, lam, 0.0, &F0)) return false;
    if (!geg_2f1_series(n - 1.0, lam + 1.0, 0.0, &F1)) return false;
    if (!geg_prefactor(n - 1.0, lam + 1.0, &P1)) return false;
    double y = P * F0, yp = 2.0 * lam * P1 * F1;
    double z0 = 0.0;
    const double mu = n * (n + 2.0 * lam), g = 2.0 * lam + 1.0;
    double c[64];
    for (int steps = 0; z0 > x; steps++) {
        if (steps > 400) return false;
        double h = -fmin(0.5 * (1.0 + z0), 1.5 / (fabs(n) + fabs(lam) + 1.0));
        double zn = z0 + h;
        if (zn < x || zn - x < 1e-4 * (1.0 + x)) zn = x;
        h = zn - z0;
        double p0 = (1.0 - z0) * (1.0 + z0);
        c[0] = y; c[1] = yp;
        double val = y + yp * h, der = yp, hp = h;
        double mag = fabs(y) + fabs(yp * h);
        int k, small = 0;
        for (k = 0; k + 2 < 64; k++) {
            double r = 2.0 * z0 * k * (k + 1.0) * c[k + 1] + k * (k - 1.0) * c[k]
                     + g * (z0 * (k + 1.0) * c[k + 1] + k * c[k]) - mu * c[k];
            c[k + 2] = r / (p0 * (k + 2.0) * (k + 1.0));
            double dterm = (k + 2.0) * c[k + 2] * hp;
            double term = c[k + 2] * hp * h;
            der += dterm; val += term; hp *= h;
            double at = fabs(term);
            if (at > mag) mag = at;
            if (at < 1e-18 * fabs(val) && fabs(dterm) < 1e-18 * fabs(der)) {
                if (++small >= 2) break;
            } else small = 0;
        }
        if (k + 2 >= 64) return false;
        double amp = fmax(fabs(val), fabs(der * h));
        if (mag > 64.0 * amp || !isfinite(val) || !isfinite(der)) return false;
        y = val; yp = der; z0 = zn;
    }
    *out = y;
    return isfinite(y);
}

/* The double algorithms, real or complex x, real n and lambda. */
static bool geg_double(double n, double lam, double xr, double xi, double* rr, double* ri) {
    if (lam == 0.0) { *rr = 0.0; *ri = 0.0; return true; }
    if (geg_is_nonneg_int(n)) {
        if (xi == 0.0) { *ri = 0.0; return geg_rec_r(n, lam, xr, rr); }
        return geg_rec_c(n, lam, xr, xi, rr, ri);
    }
    double c = lam + 0.5;
    if (c <= 0.0 && c == floor(c)) return false;          /* 2F1 c-pole: arb */
    double P;
    if (!geg_prefactor(n, lam, &P)) return false;
    if (P == 0.0) { *rr = 0.0; *ri = 0.0; return true; }  /* 1/Gamma(n+1) or 1/Gamma(2l) */
    double tr = 0.5 * (1.0 - xr), ti = -0.5 * xi;
    if (xi == 0.0 && tr <= 1.0) {
        double F;
        if (tr > 0.6 && xr > -1.0 && geg_ode_neg(n, lam, P, xr, rr)) { *ri = 0.0; return true; }
        if (!(tr <= 0.6 && geg_2f1_series(n, lam, xr, &F)) &&
            !sf_machine_hyper2f1(-n, n + 2.0 * lam, c, tr, &F)) return false;
        *rr = P * F; *ri = 0.0;
        return true;
    }
    /* complex value: z complex, or z < -1 (t > 1, on the 2F1 cut: the limit
     * from below, ti = -0, the reference convention) */
    if (xi == 0.0) ti = -0.0;
    double Fr, Fi;
    if (!sf_machine_hyper2f1_c(-n, 0.0, n + 2.0 * lam, 0.0, c, 0.0, tr, ti, &Fr, &Fi))
        return false;
    *rr = P * Fr; *ri = P * Fi;
    return true;
}

bool sf_machine_gegenbauer_c(double nr, double ni, double lr, double li,
                             double xr, double xi, double* rr, double* ri) {
    if (!isfinite(nr) || !isfinite(ni) || !isfinite(lr) || !isfinite(li) ||
        !isfinite(xr) || !isfinite(xi)) return false;
    double vr, vi;
    if (ni == 0.0 && li == 0.0 && geg_double(nr, lr, xr, xi, &vr, &vi) &&
        geg_finish(vr, vi, rr, ri))
        return true;
    if (!flint_gegenbauer_d(nr, ni, lr, li, xr, xi, &vr, &vi)) return false;
    return geg_finish(vr, vi, rr, ri);
}

bool sf_machine_gegenbauer(double n, double lam, double x, double* out) {
    double ri;
    if (!sf_machine_gegenbauer_c(n, 0.0, lam, 0.0, x, 0.0, out, &ri)) return false;
    return ri == 0.0;            /* a complex value (z < -1, non-integer n) declines */
}

/* T_n(x) for real n, x by the trigonometric / hyperbolic form (cos(n acos x)
 * on [-1, 1], cosh(n acosh x) above, cos(n pi) cosh(n t) + i sin(n pi) sinh(n t)
 * below, t = acosh(-x)), the recurrence for an integer n <= 1000. */
static bool geg_cheb_t(double n, double x, double* rr, double* ri) {
    if (n == floor(n) && fabs(n) <= 1000.0) {
        /* T_{k+1} = 2x T_k - T_{k-1} in double-double: the double recurrence
         * loses ~n ulp near a zero (421 ulp at n = 12 measured) */
        long N = (long)fabs(n);
        *ri = 0.0;
        if (N == 0) { *rr = 1.0; return true; }
        sdd t0 = sdd_from(1.0), t1 = sdd_from(x);
        for (long k = 1; k < N; k++) {
            sdd t2 = sdd_sub(sdd_mul_d(t1, 2.0 * x), t0);
            t0 = t1; t1 = t2;
        }
        *rr = t1.hi + t1.lo;
        return isfinite(*rr);
    }
    if (fabs(x) <= 1.0) { *rr = cos(n * acos(x)); *ri = 0.0; return true; }
    double t = acosh(fabs(x));
    if (x > 0.0) { *rr = cosh(n * t); *ri = 0.0; return isfinite(*rr); }
    /* cos(n pi) and sin(n pi) with the argument reduced exactly */
    double r = n - 2.0 * floor(0.5 * n + 0.5);
    double cp = cos(M_PI * r), sp = sin(M_PI * r);
    *rr = cp * cosh(n * t);
    *ri = sp * sinh(n * t);
    return isfinite(*rr) && isfinite(*ri);
}

bool sf_machine_gegenbauer2_c(double nr, double ni, double xr, double xi,
                              double* rr, double* ri) {
    if (!isfinite(nr) || !isfinite(ni) || !isfinite(xr) || !isfinite(xi)) return false;
    if (nr == 0.0 && ni == 0.0) return false;          /* 2/0: ComplexInfinity */
    double tr, ti;
    if (ni == 0.0 && xi == 0.0) {
        if (!geg_cheb_t(nr, xr, &tr, &ti)) return false;
    } else if (ni == 0.0 && nr == floor(nr) && nr >= 0.0 && nr <= 1000.0) {
        if (!sf_machine_chebyshev_t_c(nr, xr, xi, &tr, &ti)) return false;
    } else {
        return false;                                  /* complex T: the interpreter */
    }
    /* 2 T / n, n possibly complex */
    double s = 2.0;
    if (ni == 0.0) return geg_finish(s * tr / nr, s * ti / nr, rr, ri);
    double d = nr * nr + ni * ni;
    return geg_finish(s * (tr * nr + ti * ni) / d, s * (ti * nr - tr * ni) / d, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtin                                                             */
/* ------------------------------------------------------------------ */

static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}

static bool geg_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}
static bool geg_is_half(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Rational &&
           geg_is_int(e->data.function.args[0], 1) && geg_is_int(e->data.function.args[1], 2);
}

/* The explicit polynomial, sum_k (-1)^k (l)_{n-k} (2z)^{n-2k} / (k! (n-2k)!),
 * with the rising factorials as products (l)(l+1)... so a symbolic lambda gives
 * the reference's factored coefficients. */
static Expr* geg_expand(int64_t n, const Expr* lam, const Expr* z) {
    size_t nt = (size_t)(n / 2 + 1);
    Expr** terms = malloc(sizeof(Expr*) * nt);
    if (!terms) return NULL;
    for (int64_t k = 0; k <= n / 2; k++) {
        int64_t m = n - k;                                   /* (l)_m */
        Expr** f = malloc(sizeof(Expr*) * (size_t)(m + 2));
        if (!f) { for (int64_t j = 0; j < k; j++) expr_free(terms[j]); free(terms); return NULL; }
        size_t nf = 0;
        for (int64_t j = 0; j < m; j++)
            f[nf++] = j == 0 ? expr_copy((Expr*)lam)
                             : mk_fn2("Plus", expr_copy((Expr*)lam), expr_new_integer(j));
        /* (-1)^k 2^(n-2k) / (k! (n-2k)!) as an exact rational */
        Expr* coef = mk_fn2("Times",
            mk_fn2("Power", expr_new_integer(-1), expr_new_integer(k)),
            mk_fn2("Times",
                mk_fn2("Power", expr_new_integer(2), expr_new_integer(n - 2 * k)),
                mk_fn2("Power",
                    mk_fn2("Times",
                        expr_new_function(expr_new_symbol("Factorial"),
                                          (Expr*[]){ expr_new_integer(k) }, 1),
                        expr_new_function(expr_new_symbol("Factorial"),
                                          (Expr*[]){ expr_new_integer(n - 2 * k) }, 1)),
                    expr_new_integer(-1))));
        f[nf++] = coef;
        f[nf++] = mk_fn2("Power", expr_copy((Expr*)z), expr_new_integer(n - 2 * k));
        terms[k] = expr_new_function(expr_new_symbol("Times"), f, nf);
        free(f);
    }
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), terms, nt);
    free(terms);
    return eval_and_free(sum);
}

static Expr* geg_three(Expr* res) {
    Expr* n = res->data.function.args[0];
    Expr* lam = res->data.function.args[1];
    Expr* z = res->data.function.args[2];

    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && ma.n == 3) {
        double rr, ri;
        if (sf_machine_gegenbauer_c(ma.re[0], ma.im[0], ma.re[1], ma.im[1],
                                    ma.re[2], ma.im[2], &rr, &ri))
            return ma.any_complex || ri != 0.0 ? sf_machine_complex_result(rr, ri)
                                               : expr_new_real(rr);
        /* declined: over/underflow -> the extended-range arb path */
    }

    if (geg_is_int(lam, 0)) return expr_new_integer(0);
    if (geg_is_int(n, 0)) return expr_new_integer(1);

    if (flint_num_call_inexact(res))
        return flint_num_for_call(flint_num_gegenbauer(n, lam, z), res);

    if (geg_is_int(lam, 1))
        return eval_and_free(expr_new_function(expr_new_symbol("ChebyshevU"),
                             (Expr*[]){ expr_copy(n), expr_copy(z) }, 2));
    if (geg_is_half(lam))
        return eval_and_free(expr_new_function(expr_new_symbol("LegendreP"),
                             (Expr*[]){ expr_copy(n), expr_copy(z) }, 2));
    if (n->type == EXPR_INTEGER && n->data.integer >= 1 && n->data.integer <= GEG_EXPAND_MAX) {
        /* a numeric inexact lambda / z was handled above; here everything is
         * exact or symbolic, so the polynomial is exact */
        return geg_expand(n->data.integer, lam, z);
    }
    return NULL;
}

static Expr* geg_two(Expr* res) {
    Expr* n = res->data.function.args[0];
    Expr* z = res->data.function.args[1];
    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && ma.n == 2) {
        double rr, ri;
        if (sf_machine_gegenbauer2_c(ma.re[0], ma.im[0], ma.re[1], ma.im[1], &rr, &ri))
            return ma.any_complex || ri != 0.0 ? sf_machine_complex_result(rr, ri)
                                               : expr_new_real(rr);
    }
    /* 2 ChebyshevT[n, z] / n */
    Expr* t = expr_new_function(expr_new_symbol("ChebyshevT"),
                                (Expr*[]){ expr_copy(n), expr_copy(z) }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(2), mk_fn2("Power", expr_copy(n), expr_new_integer(-1)), t }, 3));
}

Expr* builtin_gegenbauer_c(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 3) return geg_three(res);
    if (argc == 2) return geg_two(res);
    return builtin_arg_error("GegenbauerC", argc, 2, 3);
}

void gegenbauer_init(void) {
    symtab_add_builtin("GegenbauerC", builtin_gegenbauer_c);
    symtab_get_def("GegenbauerC")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("GegenbauerC",
        "GegenbauerC[n, lambda, z]\n"
        "\tgives the Gegenbauer (ultraspherical) polynomial C_n^lambda(z), and its\n"
        "\tcontinuation Gamma[n+2lambda]/(Gamma[2lambda] Gamma[n+1]) *\n"
        "\tHypergeometric2F1[-n, n+2lambda, lambda+1/2, (1-z)/2] for any n.\n"
        "GegenbauerC[n, z]\n"
        "\tgives 2 ChebyshevT[n, z]/n, the lambda -> 0 limit of C_n^lambda/lambda.\n"
        "An integer n >= 0 expands to the explicit polynomial (for every lambda);\n"
        "GegenbauerC[n, 0, z] = 0, GegenbauerC[n, 1, z] = ChebyshevU[n, z],\n"
        "GegenbauerC[n, 1/2, z] = LegendreP[n, z]. Numeric at machine and arbitrary\n"
        "precision, real and complex; threads over packed arrays and compiles.\n"
        "D[GegenbauerC[n, l, z], z] = 2 l GegenbauerC[n-1, l+1, z].");
}
