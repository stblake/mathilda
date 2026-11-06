/*
 * besseljzero.c -- BesselJZero[n, k] and BesselYZero[n, k], the k-th positive
 * zeros of BesselJ[n, x] and BesselY[n, x].
 *
 * An exact call stays unevaluated -- which is also what lets the Hadamard-
 * product recogniser (Product`BesselZero) match the canonical infinite product
 *
 *   Product[1 - x^2/BesselJZero[n,k]^2, {k,1,Inf}] = Gamma[n+1] (2/x)^n BesselJ[n,x].
 *
 * With an inexact argument (and so under N, which hands it the integer-valued
 * real k = 2.) it evaluates, for a real order n > -1 and a positive integer k:
 *
 *   machine      sf_bessel_j_zero_r / sf_bessel_y_zero_r, the kernels the
 *                NDArray / packed map and Compile[] share: McMahon's expansion
 *                (DLMF 10.21.19, beta = (k + n/2 - 1/4) pi for J and
 *                (k + n/2 - 3/4) pi for Y) when k is large against n,
 *                otherwise a sign-change scan from the origin (consecutive
 *                zeros are more than 2 apart, so a 1/4 step never skips one),
 *                then Newton on f / f', f' = (f_{n-1} - f_{n+1})/2
 *   arbitrary    the machine zero, refined by Newton in arb at the call's
 *                precision (flint_num_bessel_j_zero / _y_zero)
 *
 * BesselYZero is defined here for a real order n >= 0 (every zero of Y_n then
 * lies beyond n, DLMF 10.21(i)); the scan starts at max(n, 0.05), where Y_n is
 * large and negative.
 *
 * Memory contract: takes ownership of res but must not free it; returns NULL
 * (leave unevaluated) or an owned result.
 */

#include "besseljzero.h"
#include "symtab.h"
#include "expr.h"
#include "attr.h"
#include "flint_num_bridge.h"
#include "numeric.h"
#include "sf_machine_bessel.h"
#include "sf_scalar.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* J_n (y == 0) or Y_n (y == 1). */
static bool bz_f(int y, double nu, double x, double* out) {
    return y ? sf_bessel_y_r(nu, x, out) : sf_bessel_j_r(nu, x, out);
}

/* f_n'(x) = (f_{n-1}(x) - f_{n+1}(x)) / 2, f = J or Y. */
static bool bjz_dj(int y, double nu, double x, double* out) {
    double a, b;
    if (!bz_f(y, nu - 1.0, x, &a) || !bz_f(y, nu + 1.0, x, &b)) return false;
    *out = 0.5 * (a - b);
    return true;
}

static bool bjz_newton(int y, double nu, double* x) {
    for (int it = 0; it < 30; it++) {
        double f, d;
        if (!bz_f(y, nu, *x, &f) || !bjz_dj(y, nu, *x, &d) || d == 0.0) return false;
        double dx = f / d;
        *x -= dx;
        if (fabs(dx) <= 4e-16 * fabs(*x)) return true;
    }
    return true;
}

static bool bz_zero_r(int y, double nu, double kd, double* out) {
    if (!isfinite(nu) || !isfinite(kd) || kd != floor(kd) || kd < 1.0 || kd > 1e6) return false;
    if (!(y ? nu >= 0.0 : nu > -1.0) || nu > 1e4) return false;
    double x;
    if (kd > nu + 20.0) {
        /* McMahon: beta = (k + n/2 - 1/4) Pi (J), (k + n/2 - 3/4) Pi (Y);
         * mu = 4 n^2. */
        double b = (kd + 0.5 * nu - (y ? 0.75 : 0.25)) * M_PI, mu = 4.0 * nu * nu, e = 8.0 * b;
        x = b - (mu - 1.0) / e
              - 4.0 * (mu - 1.0) * (7.0 * mu - 31.0) / (3.0 * e * e * e)
              - 32.0 * (mu - 1.0) * (83.0 * mu * mu - 982.0 * mu + 3779.0) /
                    (15.0 * e * e * e * e * e);
        double g = x;
        if (!bjz_newton(y, nu, &x) || fabs(x - g) > 1.0) return false;
    } else {
        /* No zero lies below n for n >= 0; J_n > 0 near the origin for
         * -1 < n < 0, and Y_n < 0 there for n >= 0. */
        int64_t k = (int64_t)kd, seen = 0;
        double h = 0.25, a = nu > 0.5 ? nu : (y ? 0.05 : 1e-3), fa;
        if (!bz_f(y, nu, a, &fa)) return false;
        double lim = (kd + 0.5 * nu + 2.0) * M_PI + 2.0 * nu + 10.0;
        for (;;) {
            double b = a + h, fb;
            if (b > lim || !bz_f(y, nu, b, &fb)) return false;
            if (fa == 0.0 || (fa < 0.0) != (fb < 0.0)) {
                if (++seen == k) {
                    /* Bisect the bracket to a Newton-safe width. */
                    double lo = a, hi = b, flo = fa;
                    for (int it = 0; it < 40; it++) {
                        double mid = 0.5 * (lo + hi), fm;
                        if (!bz_f(y, nu, mid, &fm)) return false;
                        if ((fm < 0.0) == (flo < 0.0) && fm != 0.0) { lo = mid; flo = fm; }
                        else hi = mid;
                    }
                    x = 0.5 * (lo + hi);
                    double g = x;
                    if (!bjz_newton(y, nu, &x) || fabs(x - g) > h) x = g;
                    break;
                }
            }
            a = b; fa = fb;
        }
    }
    if (!isfinite(x) || x <= 0.0) return false;
    *out = x;
    return true;
}

bool sf_bessel_j_zero_r(double nu, double kd, double* out) { return bz_zero_r(0, nu, kd, out); }
bool sf_bessel_y_zero_r(double nu, double kd, double* out) { return bz_zero_r(1, nu, kd, out); }

/* A real machine double from an inexact or exact real leaf. */
static bool bjz_real(const Expr* e, double* out, bool* inexact) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL) { *out = e->data.real; *inexact = true; return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); *inexact = true; return true; }
#endif
    return false;
}

static Expr* bz_builtin(Expr* res, int y) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n = res->data.function.args[0];
    Expr* k = res->data.function.args[1];
    double nu, kd;
    bool inexact = false;
    if (!bjz_real(n, &nu, &inexact) || !bjz_real(k, &kd, &inexact) || !inexact)
        return NULL;   /* exact (or symbolic, or complex): stays unevaluated */
    double x;
    if (!bz_zero_r(y, nu, kd, &x)) return NULL;
#ifdef USE_MPFR
    long bits = 0;
    for (int i = 0; i < 2; i++) {
        long b = numeric_min_inexact_bits(res->data.function.args[i]);
        if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    if (bits > 53) {
        Expr* r = y ? flint_num_bessel_y_zero(n, x, bits) : flint_num_bessel_j_zero(n, x, bits);
        if (r) return r;
    }
#endif
    return expr_new_real(x);
}

Expr* builtin_besseljzero(Expr* res) { return bz_builtin(res, 0); }
Expr* builtin_besselyzero(Expr* res) { return bz_builtin(res, 1); }

void besseljzero_init(void) {
    symtab_add_builtin("BesselJZero", builtin_besseljzero);
    symtab_get_def("BesselJZero")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("BesselJZero",
        "BesselJZero[n, k]\n"
        "\tThe k-th positive zero of BesselJ[n, x], for a real order n > -1 and a\n"
        "\tpositive integer k.\n"
        "\tEvaluates when an argument is inexact (so N works); an exact call\n"
        "\tstays unevaluated. Machine values via McMahon's expansion or a\n"
        "\tsign-change scan plus Newton; arbitrary precision by Newton in arb.");
    symtab_add_builtin("BesselYZero", builtin_besselyzero);
    symtab_get_def("BesselYZero")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("BesselYZero",
        "BesselYZero[n, k]\n"
        "\tThe k-th positive zero of BesselY[n, x], for a real order n >= 0 and a\n"
        "\tpositive integer k.\n"
        "\tEvaluates when an argument is inexact (so N works); an exact call\n"
        "\tstays unevaluated. Machine values via McMahon's expansion or a\n"
        "\tsign-change scan plus Newton; arbitrary precision by Newton in arb.");
}
