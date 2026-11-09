/* Mathilda -- the Hankel functions and the spherical Bessel functions.
 *
 *   HankelH1[nu, z] = BesselJ[nu, z] + I BesselY[nu, z]
 *   HankelH2[nu, z] = BesselJ[nu, z] - I BesselY[nu, z]
 *   SphericalBesselJ[n, z]   = Sqrt[Pi/(2 z)] BesselJ[n + 1/2, z]
 *   SphericalBesselY[n, z]   = Sqrt[Pi/(2 z)] BesselY[n + 1/2, z]
 *   SphericalHankelH1[n, z]  = SphericalBesselJ[n, z] + I SphericalBesselY[n, z]
 *   SphericalHankelH2[n, z]  = SphericalBesselJ[n, z] - I SphericalBesselY[n, z]
 * principal branches (cut along the negative real axis), as in the reference
 * system, which also leaves every exact call unevaluated except at z = 0.
 *
 * The sum J + I Y cancels exactly where a Hankel function is small: H1 decays
 * like e^(i z) in the upper half plane, while J and Y each grow like e^|Im z|.
 * There the Hankel functions come from the modified Bessel K, which has no such
 * cancellation (DLMF 10.27.8):
 *   H1_nu(z) = (2/Pi) e^(-i Pi (nu+1)/2) K_nu(-i z),   -Pi/2 < arg z <= Pi
 *   H2_nu(z) = (2/Pi) e^(i Pi (nu+1)/2) K_nu(i z),     -Pi < arg z <= Pi/2
 * used for Im z > 0 (H1) / Im z < 0 (H2); elsewhere J +- I Y has no cancellation
 * (on the real axis the two are simply the real and imaginary parts).
 *
 * Numerics: machine kernel sf_machine_hankel_sph (the scalar call, the NDArray /
 * packed map and Compile[]) over the shared machine Bessel kernels
 * (sf_bessel_nd), then arb at 53 bits in-kernel for whatever they decline;
 * arbitrary precision through arb (flint_num_hankel_sph).
 * D: H'_nu = (H_{nu-1} - H_{nu+1})/2; f'_n = -f_n/(2z) + (f_{n-1} - f_{n+1})/2.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "hankel_spherical.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "flint_num_bridge.h"
#include "sf_machine_bessel.h"
#include "sf_scalar.h"
#include "symtab.h"
#include "../sym_names.h"

static const char* const HS_NAME[HS_COUNT] = {
    "HankelH1", "HankelH2", "SphericalBesselJ", "SphericalBesselY",
    "SphericalHankelH1", "SphericalHankelH2" };

/* ------------------------------------------------------------------ */
/* Machine kernel                                                      */
/* ------------------------------------------------------------------ */

static bool hs_ok(double r, double i, double* rr, double* ri) {
    if (!isfinite(r) || !isfinite(i)) return false;
    double big = fmax(fabs(r), fabs(i));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(r) < DBL_MIN) r = 0.0;
    if (fabs(i) < DBL_MIN) i = 0.0;
    *rr = r; *ri = i;
    return true;
}

/* cos(pi x), sin(pi x) with the argument reduced exactly. */
static void hs_cs_pi(double x, double* c, double* s) {
    double r = x - 2.0 * floor(0.5 * x + 0.5);        /* [-1, 1) */
    *c = cos(M_PI * r);
    *s = sin(M_PI * r);
    if (r == 0.5 || r == -0.5) *c = 0.0;
    if (r == 0.0 || r == -1.0) *s = 0.0;
}

/* H^(1) (kind 1) or H^(2) (kind 2) at real order nu, complex z. */
static bool hs_hankel(int kind, double nu, double zr, double zi, double* rr, double* ri) {
    if (zr == 0.0 && zi == 0.0) return false;
    bool use_k = (kind == 1) ? (zi > 0.0) : (zi < 0.0);
    if (use_k) {
        /* H1 = (2/Pi) e^(-i Pi (nu+1)/2) K(-i z);  H2 = (2/Pi) e^(i Pi (nu+1)/2) K(i z) */
        double wr = (kind == 1) ? zi : -zi, wi = (kind == 1) ? -zr : zr;
        double kr, ki;
        if (!sf_bessel_nd(SF_BESSEL_K, nu, 0.0, wr, wi, &kr, &ki)) return false;
        double c, s;
        hs_cs_pi(0.5 * (nu + 1.0), &c, &s);
        double pr, pi;
        if (kind == 1) { pr = c * kr + s * ki; pi = c * ki - s * kr; }        /* e^{-i t} K */
        else           { pr = c * kr - s * ki; pi = c * ki + s * kr; }        /* e^{i t} K */
        return hs_ok((2.0 / M_PI) * pr, (2.0 / M_PI) * pi, rr, ri);
    }
    double jr, ji, yr, yi;
    if (!sf_bessel_nd(SF_BESSEL_J, nu, 0.0, zr, zi, &jr, &ji)) return false;
    if (!sf_bessel_nd(SF_BESSEL_Y, nu, 0.0, zr, zi, &yr, &yi)) return false;
    double sg = (kind == 1) ? 1.0 : -1.0;
    return hs_ok(jr - sg * yi, ji + sg * yr, rr, ri);
}

/* sqrt(pi/(2 z)), principal, as (re, im). */
static void hs_pref(double zr, double zi, double* pr, double* pi) {
    /* Sqrt[Pi]/Sqrt[2 z] with the PRINCIPAL Sqrt[2 z] -- the reference's form:
     * on the cut (z < 0) it differs in sign from Sqrt[Pi/(2 z)]
     * (SphericalBesselJ[2.5, -2.] = 0.114273 I there, not -0.114273 I) */
    double ar = 2.0 * zr, ai = 2.0 * zi;
    double m = hypot(ar, ai);
    double sr = sqrt(0.5 * (m + ar));
    double si = sqrt(0.5 * (m - ar));
    if (ai < 0.0) si = -si;
    double d = sr * sr + si * si;              /* Sqrt[Pi] / (sr + i si) */
    *pr = 1.77245385090551602729816748334 * sr / d;
    *pi = -1.77245385090551602729816748334 * si / d;
}

/* j_n (code 2), y_n (code 3) for an integer 0 <= n <= 200 and real x > 0: y by
 * its upward recurrence (always stable: y is dominant), j by the same
 * recurrence where it is stable (x > n), from the elementary j_0, j_1, y_0, y_1:
 *   f_{k+1} = (2k + 1)/x f_k - f_{k-1}.
 * ~20-50 ns where the half-integer Bessel kernel costs ~0.5 us. j for x <= n is
 * left to that kernel (the upward j recurrence loses digits there). */
static bool hs_sph_int(int code, int n, double x, double* out) {
    if (code == 2 && x <= (double)n) return false;
    double s = sin(x), c = cos(x), ix = 1.0 / x;
    double f0, f1;
    if (code == 2) { f0 = s * ix; f1 = (s * ix - c) * ix; }
    else           { f0 = -c * ix; f1 = (-c * ix - s) * ix; }
    if (n == 0) { *out = f0; return isfinite(f0); }
    for (int k = 1; k < n; k++) {
        double f2 = (2.0 * k + 1.0) * ix * f1 - f0;
        f0 = f1; f1 = f2;
    }
    if (!isfinite(f1)) return false;
    *out = f1;
    return true;
}

/* the spherical functions; code 2 j, 3 y, 4 h1, 5 h2 */
static bool hs_spherical(int code, double n, double zr, double zi, double* rr, double* ri) {
    double nu = n + 0.5;
    if (zr == 0.0 && zi == 0.0) {
        if (code == 2 && n == 0.0) { *rr = 1.0; *ri = 0.0; return true; }
        if (code == 2 && n > 0.0 && n == floor(n)) { *rr = 0.0; *ri = 0.0; return true; }
        return false;
    }
    /* integer order on the negative axis: the parity relations keep it real */
    if (zi == 0.0 && zr < 0.0 && n == floor(n)) {
        double sg = fmod(fabs(n), 2.0) == 0.0 ? 1.0 : -1.0;     /* (-1)^n */
        double vr, vi;
        if (code == 2 || code == 3) {
            if (!hs_spherical(code, n, -zr, 0.0, &vr, &vi)) return false;
            return hs_ok(code == 2 ? sg * vr : -sg * vr, 0.0, rr, ri);
        }
        /* h1(-x) = (-1)^n h2(x), h2(-x) = (-1)^n h1(x) */
        if (!hs_spherical(code == 4 ? 5 : 4, n, -zr, 0.0, &vr, &vi)) return false;
        return hs_ok(sg * vr, sg * vi, rr, ri);
    }
    if (zi == 0.0 && zr > 0.0 && n == floor(n) && n >= 0.0 && n <= 200.0) {
        /* the fast integer-order recurrences; h = j +- i y */
        double v, w;
        if (code == 2 || code == 3) {
            if (hs_sph_int(code, (int)n, zr, &v)) return hs_ok(v, 0.0, rr, ri);
        } else if (hs_sph_int(3, (int)n, zr, &w)) {
            double vi = 0.0;
            if (!hs_sph_int(2, (int)n, zr, &v)) {
                double jr;
                if (!sf_bessel_nd(SF_BESSEL_J, nu, 0.0, zr, 0.0, &jr, &vi)) return false;
                v = jr * 1.77245385090551602729816748334 / sqrt(2.0 * zr);
            }
            return hs_ok(v, code == 4 ? w : -w, rr, ri);
        }
    }
    double cr, ci;
    if (code == 2 || code == 3) {
        if (!sf_bessel_nd(code == 2 ? SF_BESSEL_J : SF_BESSEL_Y, nu, 0.0, zr, zi, &cr, &ci))
            return false;
    } else {
        if (!hs_hankel(code == 4 ? 1 : 2, nu, zr, zi, &cr, &ci)) return false;
    }
    double pr, pi;
    hs_pref(zr, zi, &pr, &pi);
    return hs_ok(pr * cr - pi * ci, pr * ci + pi * cr, rr, ri);
}

bool sf_machine_hankel_sph(int code, double nr, double ni, double zr, double zi,
                           double* rr, double* ri) {
    if (code < 0 || code >= HS_COUNT) return false;
    if (!isfinite(nr) || !isfinite(ni) || !isfinite(zr) || !isfinite(zi)) return false;
    double vr, vi;
    if (ni == 0.0) {
        bool ok = code <= 1 ? hs_hankel(code + 1, nr, zr, zi, &vr, &vi)
                            : hs_spherical(code, nr, zr, zi, &vr, &vi);
        if (ok) {
            /* a real input on the positive axis: j, y real */
            if (zi == 0.0 && zr > 0.0 && (code == 2 || code == 3)) vi = 0.0;
            *rr = vr; *ri = vi;
            return true;
        }
    }
    if (zr == 0.0 && zi == 0.0) return false;
    if (!flint_hankel_sph_d(code, nr, ni, zr, zi, &vr, &vi)) return false;
    return hs_ok(vr, vi, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtins                                                            */
/* ------------------------------------------------------------------ */

static bool hs_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

static Expr* hs_builtin(Expr* res, int code) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error(HS_NAME[code], argc, 2, 2);
    Expr* n = res->data.function.args[0];
    Expr* z = res->data.function.args[1];

    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && ma.n == 2) {
        double rr, ri;
        if (sf_machine_hankel_sph(code, ma.re[0], ma.im[0], ma.re[1], ma.im[1], &rr, &ri)) {
            bool real_out = (code == 2 || code == 3) && !ma.any_complex && ri == 0.0;
            return real_out ? expr_new_real(rr) : sf_machine_complex_result(rr, ri);
        }
    }
    /* z = 0 (exact) */
    if (hs_is_int(z, 0)) {
        if (code <= 1) {
            /* H(nu, 0): Y is infinite for every real order */
            if (n->type == EXPR_INTEGER || n->type == EXPR_REAL) return expr_new_symbol("ComplexInfinity");
            return NULL;
        }
        if (code == 2 && hs_is_int(n, 0)) return expr_new_integer(1);
        if (code == 2 && n->type == EXPR_INTEGER && n->data.integer > 0) return expr_new_integer(0);
        if (code != 2 && n->type == EXPR_INTEGER && n->data.integer >= 0)
            return expr_new_symbol("ComplexInfinity");
        return NULL;
    }
    if (flint_num_call_inexact(res))
        return flint_num_for_call(flint_num_hankel_sph(code, n, z), res);
    return NULL;
}

#define HS_WRAP(NAME, CODE) static Expr* builtin_##NAME(Expr* res) { return hs_builtin(res, CODE); }
HS_WRAP(hankel_h1, HS_H1) HS_WRAP(hankel_h2, HS_H2)
HS_WRAP(sph_j, HS_SJ) HS_WRAP(sph_y, HS_SY) HS_WRAP(sph_h1, HS_SH1) HS_WRAP(sph_h2, HS_SH2)

static const char* const HS_DOC[HS_COUNT] = {
    "HankelH1[nu, z]\n\tgives the Hankel function of the first kind, BesselJ[nu, z] + I BesselY[nu, z].",
    "HankelH2[nu, z]\n\tgives the Hankel function of the second kind, BesselJ[nu, z] - I BesselY[nu, z].",
    "SphericalBesselJ[n, z]\n\tgives the spherical Bessel function j_n(z) = Sqrt[Pi/(2z)] BesselJ[n+1/2, z].",
    "SphericalBesselY[n, z]\n\tgives the spherical Bessel function y_n(z) = Sqrt[Pi/(2z)] BesselY[n+1/2, z].",
    "SphericalHankelH1[n, z]\n\tgives SphericalBesselJ[n, z] + I SphericalBesselY[n, z].",
    "SphericalHankelH2[n, z]\n\tgives SphericalBesselJ[n, z] - I SphericalBesselY[n, z]."
};

void hankel_spherical_init(void) {
    Expr* (*fns[HS_COUNT])(Expr*) = { builtin_hankel_h1, builtin_hankel_h2, builtin_sph_j,
                                      builtin_sph_y, builtin_sph_h1, builtin_sph_h2 };
    for (int i = 0; i < HS_COUNT; i++) {
        symtab_add_builtin(HS_NAME[i], fns[i]);
        symtab_get_def(HS_NAME[i])->attributes |=
            (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
        char buf[600];
        snprintf(buf, sizeof buf, "%s\n%s", HS_DOC[i],
                 "Principal branches (cut along the negative real axis); exact calls stay\n"
                 "unevaluated except at z = 0. Numeric at machine and arbitrary precision,\n"
                 "real and complex (the Hankel forms through BesselK where J + I Y would\n"
                 "cancel); threads over packed arrays and compiles.");
        symtab_set_docstring(HS_NAME[i], buf);
    }
}
