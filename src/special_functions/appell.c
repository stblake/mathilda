/* Mathilda -- the Appell hypergeometric functions AppellF1, AppellF2, AppellF3,
 * AppellF4 (DLMF 16.13; definitions in appell.h).
 *
 * METHOD. Every evaluation is one outer series whose terms carry a Gauss
 * function of the other variable (DLMF 16.15.1-4 summed over one index):
 *
 *   F1 = Sum_m (a)_m (b1)_m / ((c)_m m!)  x^m 2F1(a + m, b2; c + m; y)
 *   F2 = Sum_m (a)_m (b1)_m / ((c1)_m m!) x^m 2F1(a + m, b2; c2; y)
 *   F3 = Sum_m (a1)_m (b1)_m / ((c)_m m!) x^m 2F1(a2, b2; c + m; y)
 *   F4 = Sum_m (a)_m (b)_m / ((c1)_m m!)  x^m 2F1(a + m, b + m; c2; y)
 *
 * with the inner function continued analytically (principal branch, the cut
 * y > 1 taken from below). The outer series converges for |x| < 1 (F1, F3),
 * |x| < min(1, |1 - y|) (F2) and |x| < min |1 -+ sqrt y|^2 (F4) -- for ANY y
 * off its cut, which is the analytic continuation of the double series in the
 * inner variable. Either variable may be the outer one. F1 additionally has the
 * five substitutions of its Euler integral that map [0, 1] onto itself
 * (DLMF 16.16.1-5 and their compositions), and F2 the three Pfaff-type ones
 * (DLMF 16.16.6-8); the planner (appell_plan) takes the transformation and
 * orientation with the smallest outer ratio. The F1 substitutions t -> t/(1 -
 * x + x t) move a segment of the Euler integral onto an arc, so two of them are
 * valid only on one side of an extra ray (the sector tests below; checked
 * against the integral on both sides).
 *
 * The inner functions G_m are NOT computed one by one: they satisfy a
 * three-term recurrence in m (F1: Gauss's (a, c) -> (a + 1, c + 1) relation,
 * derived from DLMF 15.5.16 / 15.5.18; F2: 15.5.11; F3: 15.5.16; F4: the
 * (a, b) -> (a + 1, b + 1) relation, derived from the hypergeometric equation
 * in theta form), run in the stable direction: forward when G_m is the dominant
 * solution, backward from two directly summed values when it is the minimal
 * one (F1 with |y| < 1, F3 with |1 - 1/y| > 1). Two Gauss evaluations and a
 * linear recurrence per call, so the machine kernel costs about as much as
 * three 2F1 evaluations.
 *
 * When the recurrence's two roots are within a factor 1.5 in modulus the G_m
 * are evaluated one by one instead (the stable direction is then only
 * polynomially stable), as is a degenerate step (F4 at a + b - c2 + 1 + 2m = 0).
 *
 * Machine precision (sf_machine_appell) runs this in double complex; a sum
 * whose terms cancel (sum of |terms| > 4 |sum|), a ratio above 0.97, a
 * denominator parameter left of 1/2, F3 with an inner variable beyond 0.6,
 * complex parameters beyond the machine 2F1, or no convergent expansion at all
 * go to arb at 53 bits INSIDE the kernel (each threshold measured against
 * mpmath), so the scalar, packed / NDArray and Compile[] answers stay
 * bit-identical. Arbitrary precision (appell_arb.c) runs the same plan in ball
 * arithmetic with an asymptotic tail bound, and for F1 the Euler integral
 * (series at the singular endpoints, acb_calc_integrate between) where no
 * expansion converges.
 *
 * Exact rules (the reference system's): x = 0 or y = 0 and the vanishing
 * parameters reduce to 2F1; F1 also x = y, a = c, c = b1 + b2, x = 1, y = 1,
 * a = 0 and a a negative integer (a polynomial). D in x and y by the
 * contiguous-parameter formulas, e.g. D[F1, x] = a b1/c F1(a+1, b1+1, b2, c+1).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "appell.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "numeric.h"
#include "flint_num_bridge.h"
#include "sf_scalar.h"
#include "sf_machine.h"
#include "sf_machine_hyper.h"
#include "symtab.h"
#include "../sym_names.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Double complex arithmetic (no <complex.h>)                          */
/* ------------------------------------------------------------------ */

typedef struct { double r, i; } cd;

static inline cd cd_(double r, double i) { cd z; z.r = r; z.i = i; return z; }
static inline cd cadd(cd a, cd b) { return cd_(a.r + b.r, a.i + b.i); }
static inline cd csub(cd a, cd b) { return cd_(a.r - b.r, a.i - b.i); }
static inline cd cmul(cd a, cd b) { return cd_(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
static inline cd caddr(cd a, double s) { return cd_(a.r + s, a.i); }
static inline cd cscal(cd a, double s) { return cd_(a.r * s, a.i * s); }
static inline double cabs_(cd a) { return hypot(a.r, a.i); }
static inline bool cis0(cd a) { return a.r == 0.0 && a.i == 0.0; }
static inline cd cconj(cd a) { return cd_(a.r, -a.i); }

static cd cdiv(cd a, cd b) {
    /* Smith's algorithm */
    if (fabs(b.r) >= fabs(b.i)) {
        if (b.r == 0.0) return cd_(a.r / 0.0, a.i / 0.0);
        double t = b.i / b.r, d = b.r + b.i * t;
        return cd_((a.r + a.i * t) / d, (a.i - a.r * t) / d);
    }
    double t = b.r / b.i, d = b.r * t + b.i;
    return cd_((a.r * t + a.i) / d, (a.i * t - a.r) / d);
}

static cd cexp_(cd a) { double e = exp(a.r); return cd_(e * cos(a.i), e * sin(a.i)); }

/* w^p on the principal branch (arg w in (-pi, pi]; a negative real w has
 * arg +pi, the value from the upper side, which is what 1 - x is for x on its
 * cut taken from below) */
static cd cpow_(cd w, cd p) {
    if (cis0(p)) return cd_(1.0, 0.0);
    if (cis0(w)) return p.r > 0.0 ? cd_(0.0, 0.0) : cd_(HUGE_VAL, 0.0);
    if (w.i == 0.0 && w.r > 0.0 && p.i == 0.0) return cd_(pow(w.r, p.r), 0.0);
    double ar = atan2(w.i, w.r);
    if (w.i == 0.0 && w.r < 0.0) ar = M_PI;
    cd l = cd_(log(cabs_(w)), ar);
    return cexp_(cmul(p, l));
}

/* sqrt on the principal branch */
static cd csqrt_(cd a) {
    if (cis0(a)) return a;
    double m = cabs_(a);
    double r = sqrt(0.5 * (m + fabs(a.r)));
    if (a.r >= 0.0) return cd_(r, a.i / (2.0 * r));
    return cd_(fabs(a.i) / (2.0 * r), a.i >= 0.0 ? r : -r);
}

static bool c_nonpos_int(cd a) { return a.i == 0.0 && a.r <= 0.0 && a.r == floor(a.r); }

/* ------------------------------------------------------------------ */
/* The transformation table                                            */
/* ------------------------------------------------------------------ */

/* F1 parameters (a, b1, b2, c) = indices 0..3; F2 (a, b1, b2, c1, c2) 0..4. */
#define ID5 { {1,0,0,0,0}, {0,1,0,0,0}, {0,0,1,0,0}, {0,0,0,1,0}, {0,0,0,0,1} }
static const AppellTr TR_F1[] = {
    /* T0: the identity */
    { 0, {0, 0}, {{0}}, ID5, AP_V_X, AP_V_Y, 0 },
    /* T3 (t -> 1 - t): (1-x)^-b1 (1-y)^-b2 F1(c-a; b1, b2; c; x/(x-1), y/(y-1)) */
    { 2, {AP_B_1X, AP_B_1Y}, {{0,-1,0,0,0}, {0,0,-1,0,0}},
      { {-1,0,0,1,0}, {0,1,0,0,0}, {0,0,1,0,0}, {0,0,0,1,0}, {0} }, AP_V_XM, AP_V_YM, 0 },
    /* T1: (1-x)^-a F1(a; c-b1-b2, b2; c; x/(x-1), (y-x)/(1-x)) */
    { 1, {AP_B_1X, 0}, {{-1,0,0,0,0}, {0}},
      { {1,0,0,0,0}, {0,-1,-1,1,0}, {0,0,1,0,0}, {0,0,0,1,0}, {0} }, AP_V_XM, AP_V_YX1, 1 },
    /* T2: (1-y)^-a F1(a; b1, c-b1-b2; c; (x-y)/(1-y), y/(y-1)) */
    { 1, {AP_B_1Y, 0}, {{-1,0,0,0,0}, {0}},
      { {1,0,0,0,0}, {0,1,0,0,0}, {0,-1,-1,1,0}, {0,0,0,1,0}, {0} }, AP_V_XY1, AP_V_YM, 2 },
    /* T5 = T3 o T1: (1-x)^(c-a-b1) (1-y)^-b2 F1(c-a; c-b1-b2, b2; c; x, (y-x)/(y-1)) */
    { 2, {AP_B_1X, AP_B_1Y}, {{-1,-1,0,1,0}, {0,0,-1,0,0}},
      { {-1,0,0,1,0}, {0,-1,-1,1,0}, {0,0,1,0,0}, {0,0,0,1,0}, {0} }, AP_V_X, AP_V_YXY, 1 },
    /* T6 = T3 o T2: (1-x)^-b1 (1-y)^(c-a-b2) F1(c-a; b1, c-b1-b2; c; (x-y)/(x-1), y) */
    { 2, {AP_B_1X, AP_B_1Y}, {{0,-1,0,0,0}, {-1,0,-1,1,0}},
      { {-1,0,0,1,0}, {0,1,0,0,0}, {0,-1,-1,1,0}, {0,0,0,1,0}, {0} }, AP_V_XYX, AP_V_Y, 2 },
};
static const AppellTr TR_F2[] = {
    { 0, {0, 0}, {{0}}, ID5, AP_V_X, AP_V_Y, 0 },
    /* (1-x)^-a F2(a; c1-b1, b2; c1, c2; x/(x-1), y/(1-x)) */
    { 1, {AP_B_1X, 0}, {{-1,0,0,0,0}, {0}},
      { {1,0,0,0,0}, {0,-1,0,1,0}, {0,0,1,0,0}, {0,0,0,1,0}, {0,0,0,0,1} }, AP_V_XM, AP_V_Y1X, 0 },
    /* (1-y)^-a F2(a; b1, c2-b2; c1, c2; x/(1-y), y/(y-1)) */
    { 1, {AP_B_1Y, 0}, {{-1,0,0,0,0}, {0}},
      { {1,0,0,0,0}, {0,1,0,0,0}, {0,0,-1,0,1}, {0,0,0,1,0}, {0,0,0,0,1} }, AP_V_X1Y, AP_V_YM, 0 },
    /* (1-x-y)^-a F2(a; c1-b1, c2-b2; c1, c2; x/(x+y-1), y/(x+y-1)) */
    { 1, {AP_B_1XY, 0}, {{-1,0,0,0,0}, {0}},
      { {1,0,0,0,0}, {0,-1,0,1,0}, {0,0,-1,0,1}, {0,0,0,1,0}, {0,0,0,0,1} }, AP_V_XS, AP_V_YS, 0 },
};
static const AppellTr TR_ID[] = { { 0, {0, 0}, {{0}}, ID5, AP_V_X, AP_V_Y, 0 } };

const AppellTr* appell_tr(int kind, int tr) {
    if (tr < 0) return NULL;
    if (kind == 1) return tr < (int)(sizeof TR_F1 / sizeof TR_F1[0]) ? &TR_F1[tr] : NULL;
    if (kind == 2) return tr < (int)(sizeof TR_F2 / sizeof TR_F2[0]) ? &TR_F2[tr] : NULL;
    return tr == 0 ? &TR_ID[0] : NULL;
}

int appell_np(int kind) { return (kind == 2 || kind == 3) ? 5 : 4; }

static cd ap_var(int v, cd x, cd y) {
    cd one = cd_(1.0, 0.0);
    switch (v) {
        case AP_V_X:   return x;
        case AP_V_Y:   return y;
        case AP_V_XM:  return cdiv(x, csub(x, one));
        case AP_V_YM:  return cdiv(y, csub(y, one));
        case AP_V_YX1: return cdiv(csub(y, x), csub(one, x));
        case AP_V_XY1: return cdiv(csub(x, y), csub(one, y));
        case AP_V_YXY: return cdiv(csub(y, x), csub(y, one));
        case AP_V_XYX: return cdiv(csub(x, y), csub(x, one));
        case AP_V_Y1X: return cdiv(y, csub(one, x));
        case AP_V_X1Y: return cdiv(x, csub(one, y));
        case AP_V_XS:  return cdiv(x, caddr(cadd(x, y), -1.0));
        case AP_V_YS:  return cdiv(y, caddr(cadd(x, y), -1.0));
    }
    return x;
}

/* `t` strictly inside the angle between 0 and `lim` (same side). The
 * arguments come from the perturbed variables, so an input on its cut lands on
 * the side its limit is taken from. */
static bool ap_between(double t, double lim) {
    if (lim > 0.0) return t > 0.0 && t < lim;
    if (lim < 0.0) return t < 0.0 && t > lim;
    return false;
}

/* The sector tests of the F1 substitutions (x, y perturbed off their cuts):
 * rule 1 (T1, T5) fails when arg(y - 1) lies between 0 and arg(1 - x); rule 2
 * (T2, T6) the same with x and y exchanged. */
static bool ap_rule_ok(int rule, cd x, cd y) {
    if (rule == 0) return true;
    cd one = cd_(1.0, 0.0);
    if (rule == 1) {
        cd u = csub(y, one), w = csub(one, x);
        return !ap_between(atan2(u.i, u.r), atan2(w.i, w.r));
    }
    cd u = csub(x, one), w = csub(one, y);
    return !ap_between(atan2(u.i, u.r), atan2(w.i, w.r));
}

/* The new parameter j of transformation `t`. */
static cd ap_newp(const AppellTr* t, int j, const cd* p, int np) {
    cd s = cd_(0.0, 0.0);
    for (int k = 0; k < np; k++)
        if (t->newp[j][k]) s = cadd(s, cscal(p[k], (double)t->newp[j][k]));
    return s;
}

void appell_engine_map(int kind, int swap, int idx[6]) {
    static const int M[5][2][6] = {
        { {0}, {0} },
        /* F1 (a, b1, b2, c) */
        { {0, 1, 3, 0, 2, 3}, {0, 2, 3, 0, 1, 3} },
        /* F2 (a, b1, b2, c1, c2) */
        { {0, 1, 3, 0, 2, 4}, {0, 2, 4, 0, 1, 3} },
        /* F3 (a1, a2, b1, b2, c) */
        { {0, 2, 4, 1, 3, 4}, {1, 3, 4, 0, 2, 4} },
        /* F4 (a, b, c1, c2) */
        { {0, 1, 2, 0, 1, 3}, {0, 1, 3, 0, 1, 2} },
    };
    for (int k = 0; k < 6; k++) idx[k] = M[kind][swap ? 1 : 0][k];
}

/* Map the (transformed) parameters q[] of `kind`, with the second variable
 * outer when `swap`, onto the engine's (n1, n2, d1; A0, B0, C0). */
static void ap_engine_params(int kind, const cd* q, int swap, cd* e) {
    int idx[6];
    appell_engine_map(kind, swap, idx);
    for (int k = 0; k < 6; k++) e[k] = q[idx[k]];
}

/* Growth rate of |G_m| (the dominant root of the inner recurrence when G is
 * the dominant solution; 1 when G_m tends to a constant). */
static double ap_growth(int kind, cd Y) {
    if (kind == 2) {
        double d = cabs_(csub(cd_(1.0, 0.0), Y));
        return d < 1.0 ? 1.0 / d : 1.0;
    }
    if (kind == 4) {
        cd s = csqrt_(Y);
        double d1 = cabs_(csub(cd_(1.0, 0.0), s)), d2 = cabs_(cadd(cd_(1.0, 0.0), s));
        double d = d1 < d2 ? d1 : d2;
        return 1.0 / (d * d);
    }
    return 1.0;
}

/* The perturbed variables: a real input beyond 1 sits on its cut and is taken
 * from below, x - i0; the planner reads sides off these. */
static cd ap_perturb(double r, double i) {
    if (i == 0.0 && r > 1.0) return cd_(r, -1e-150 * (r > 1e100 ? r : 1.0));
    return cd_(r, i);
}

bool appell_plan(int kind, const double* re, const double* im, double rmax, AppellPlan* P) {
    int np = appell_np(kind);
    cd p[AP_MAXP];
    for (int k = 0; k < np; k++) p[k] = cd_(re[k], im[k]);
    cd x = ap_perturb(re[np], im[np]), y = ap_perturb(re[np + 1], im[np + 1]);
    double best = HUGE_VAL;
    bool found = false;
    for (int tr = 0; ; tr++) {
        const AppellTr* t = appell_tr(kind, tr);
        if (!t) break;
        if (!ap_rule_ok(t->rule, x, y)) continue;
        cd q[AP_MAXP];
        for (int j = 0; j < np; j++) q[j] = ap_newp(t, j, p, np);
        cd X = ap_var(t->vx, x, y), Y = ap_var(t->vy, x, y);
        for (int swap = 0; swap < 2; swap++) {
            cd e[6];
            ap_engine_params(kind, q, swap, e);
            cd Xo = swap ? Y : X, Yi = swap ? X : Y;
            if (!isfinite(Xo.r) || !isfinite(Xo.i) || !isfinite(Yi.r) || !isfinite(Yi.i)) continue;
            double r;
            if (c_nonpos_int(e[0]) || c_nonpos_int(e[1])) r = 0.0;     /* terminates */
            else r = cabs_(Xo) * ap_growth(kind, Yi);
            /* prefer the identity (and fewer factors) on a tie */
            if (r < best - 1e-12) {
                best = r; found = true;
                P->kind = kind; P->tr = tr; P->swap = swap; P->r = r;
                P->sigY = 0;
                if (fabs(Yi.i) < 1e-100 && Yi.r > 1.0)
                    P->sigY = Yi.i > 0.0 ? 1 : -1;
            }
        }
    }
    return found && best <= rmax;
}

/* ------------------------------------------------------------------ */
/* Machine evaluation                                                  */
/* ------------------------------------------------------------------ */

/* 2F1(a, b; c; z) through the machine Gauss kernel, z on its cut from below
 * (sig <= 0) or above (sig > 0: the conjugate of the value at the conjugate
 * parameters). */
static bool ap_2f1(cd a, cd b, cd c, cd z, int sig, cd* out) {
    bool realp = a.i == 0.0 && b.i == 0.0 && c.i == 0.0;
    double rr, ri;
    if (realp && z.i == 0.0 && z.r < 1.0) {
        double v;
        if (!sf_machine_hyper2f1(a.r, b.r, c.r, z.r, &v) || !isfinite(v)) return false;
        *out = cd_(v, 0.0);
        return true;
    }
    if (z.i == 0.0 && z.r > 1.0 && sig > 0) {
        if (!sf_machine_hyper2f1_c(a.r, -a.i, b.r, -b.i, c.r, -c.i, z.r, -0.0, &rr, &ri)) return false;
        ri = -ri;
    } else if (z.i == 0.0 && z.r >= 1.0) {
        if (!sf_machine_hyper2f1_c(a.r, a.i, b.r, b.i, c.r, c.i, z.r, -0.0, &rr, &ri)) return false;
    } else {
        if (!sf_machine_hyper2f1_c(a.r, a.i, b.r, b.i, c.r, c.i, z.r, z.i, &rr, &ri)) return false;
    }
    if (!isfinite(rr) || !isfinite(ri)) return false;
    *out = cd_(rr, ri);
    return true;
}

/* The direct series of 2F1, |z| <= 0.995 (the backward recurrence's two start
 * values at large m, where the parameters are big; up to 20000 terms near
 * |z| = 1). */
static bool ap_ser2f1(cd a, cd b, cd c, cd z, cd* out) {
    if (!(cabs_(z) <= 0.995)) return false;
    cd s = cd_(1.0, 0.0), t = s;
    double sabs = 1.0;
    for (int k = 0; k < 20000; k++) {
        cd num = cmul(caddr(a, k), caddr(b, k));
        cd den = cscal(caddr(c, k), (double)(k + 1));
        t = cmul(cdiv(num, den), cmul(t, z));
        s = cadd(s, t);
        double at = cabs_(t);
        sabs += at;
        double q = cabs_(cdiv(cmul(caddr(a, k + 1), caddr(b, k + 1)),
                                cscal(caddr(c, k + 1), (double)(k + 2)))) * cabs_(z);
        if (q < 1.0 && at * q / (1.0 - q) <= 1e-18 * cabs_(s)) {
            if (!(sabs <= 1e3 * cabs_(s)) || !isfinite(s.r) || !isfinite(s.i)) return false;
            *out = s;
            return true;
        }
        if (c_nonpos_int(caddr(a, k + 1)) || c_nonpos_int(caddr(b, k + 1))) { *out = s; return true; }
    }
    return false;
}

/* The inner recurrence alpha P_m + beta P_(m+1) + gamma P_(m+2) = 0 at index m. */
static void ap_rec(int kind, cd A0, cd B0, cd C0, cd z, int m, cd* al, cd* be, cd* ga) {
    cd one = cd_(1.0, 0.0);
    if (kind == 1) {          /* G_m = 2F1(A0 + m, b; C0 + m; z) */
        cd A = caddr(A0, m), C = caddr(C0, m), b = B0;
        cd A1 = caddr(A, 1.0), C1 = caddr(C, 1.0);
        *ga = cmul(cmul(A1, csub(C1, b)), z);
        *be = cscal(cmul(C1, cadd(C, cmul(csub(A1, b), z))), -1.0);
        *al = cmul(C, C1);
    } else if (kind == 2) {   /* G_m = 2F1(A0 + m, b; c; z) */
        cd A1 = caddr(A0, m + 1.0), b = B0, c = C0;
        *al = csub(c, A1);
        *be = cadd(csub(cscal(A1, 2.0), c), cmul(csub(b, A1), z));
        *ga = cmul(A1, csub(z, one));
    } else if (kind == 3) {   /* G_m = 2F1(a, b; C0 + m; z) */
        cd a = A0, b = B0, C = caddr(C0, m), C1 = caddr(C, 1.0);
        *al = cmul(cmul(C1, C), csub(z, one));
        cd t = cmul(csub(caddr(cscal(C, 2.0), 1.0), cadd(a, b)), z);
        *be = cmul(C1, csub(C, t));
        *ga = cmul(cmul(csub(C1, a), csub(C1, b)), z);
    } else {                  /* G_m = 2F1(A0 + m, B0 + m; c; z) */
        cd A = caddr(A0, m), B = caddr(B0, m), c = C0;
        cd zm1 = csub(z, one);
        cd ABc = csub(cadd(A, B), c);
        *ga = cmul(cmul(caddr(A, 1.0), caddr(B, 1.0)), cmul(cmul(zm1, zm1), caddr(ABc, 1.0)));
        *al = cmul(cmul(caddr(csub(A, c), 1.0), caddr(csub(B, c), 1.0)), caddr(ABc, 3.0));
        /* A^2 z + 2AB - Acz - Ac + 2Az + 2A + B^2 z - Bcz - Bc + 2Bz + 2B + c^2 - cz - 3c + z + 2 */
        cd Q = cmul(cmul(A, A), z);
        Q = cadd(Q, cscal(cmul(A, B), 2.0));
        Q = csub(Q, cmul(cmul(A, c), z));
        Q = csub(Q, cmul(A, c));
        Q = cadd(Q, cscal(cmul(A, z), 2.0));
        Q = cadd(Q, cscal(A, 2.0));
        Q = cadd(Q, cmul(cmul(B, B), z));
        Q = csub(Q, cmul(cmul(B, c), z));
        Q = csub(Q, cmul(B, c));
        Q = cadd(Q, cscal(cmul(B, z), 2.0));
        Q = cadd(Q, cscal(B, 2.0));
        Q = cadd(Q, cmul(c, c));
        Q = csub(Q, cmul(c, z));
        Q = csub(Q, cscal(c, 3.0));
        Q = cadd(Q, z);
        Q = caddr(Q, 2.0);
        *be = cscal(cmul(caddr(ABc, 2.0), Q), -1.0);
    }
}

/* The inner recurrence runs backward (G minimal) for F1 with |Y| < 1 and F3
 * with |1 - 1/Y| > 1; forward otherwise. */
static bool ap_backward(int kind, cd Y) {
    if (kind == 1) return cabs_(Y) < 1.0;
    if (kind == 3) return !cis0(Y) && cabs_(csub(cd_(1.0, 0.0), cdiv(cd_(1.0, 0.0), Y))) > 1.0;
    return false;
}

/* The ratio of the moduli of the inner recurrence's two characteristic roots
 * (F1: 1, 1/Y; F2: 1, 1/(1-Y); F3: 1, 1 - 1/Y; F4: (1 -+ sqrt Y)^-2). Close to
 * 1 the stable direction is only polynomially stable and the polynomial
 * factors can win over a few hundred steps, so the G_m are then evaluated one
 * by one instead. */
static double ap_separation(int kind, cd Y) {
    cd one = cd_(1.0, 0.0);
    double r;
    if (kind == 1) r = cis0(Y) ? HUGE_VAL : 1.0 / cabs_(Y);
    else if (kind == 2) r = 1.0 / cabs_(csub(one, Y));
    else if (kind == 3) r = cis0(Y) ? HUGE_VAL : cabs_(csub(one, cdiv(one, Y)));
    else {
        cd s = csqrt_(Y);
        double d1 = cabs_(csub(one, s)), d2 = cabs_(cadd(one, s));
        r = (d2 * d2) / (d1 * d1);
    }
    if (!(r == r)) return 1.0;
    return r >= 1.0 ? r : 1.0 / r;
}
#define AP_SEP_MIN 1.5

void appell_engine_info(int kind, double yr, double yi, double* growth, int* backward) {
    cd Y = cd_(yr, yi);
    *growth = ap_growth(kind, Y);
    *backward = ap_backward(kind, Y) ? 1 : 0;
}

static void ap_shift(int kind, int* sA, int* sB, int* sC) {
    *sA = (kind == 1 || kind == 2 || kind == 4);
    *sB = (kind == 4);
    *sC = (kind == 1 || kind == 3);
}

#define AP_MAXTERMS 6000

/* S = Sum_m T_m G_m in double complex (see the header comment). False when
 * the sum cannot be certified to a few ulp. */
static bool ap_eval_d(int kind, const cd* e, cd X, cd Y, int sig, cd* out) {
    cd n1 = e[0], n2 = e[1], d1 = e[2], A0 = e[3], B0 = e[4], C0 = e[5];
    /* a denominator parameter left of 1/2: the first coefficients (and the
     * Gauss functions) pass near a pole and cancel -- measured up to 3600 ulp
     * in double -- so arb answers */
    if (d1.r < 0.5 || C0.r < 0.5) return false;
    /* F3's inner functions 2F1(a2, b2; c + m; Y) near |Y| = 1 lose up to ~140
     * ulp through the recurrence in c (measured against mpmath); arb answers */
    if (kind == 3 && (cabs_(Y) > 0.6 || C0.r < 1.0)) return false;
    int sA, sB, sC;
    ap_shift(kind, &sA, &sB, &sC);
    /* the number of outer terms */
    long M;
    double g = ap_growth(kind, Y);
    bool term = false;
    long Mt = AP_MAXTERMS + 1;
    if (c_nonpos_int(n1)) { term = true; Mt = (long)(-n1.r); }
    if (c_nonpos_int(n2) && (long)(-n2.r) < Mt) { term = true; Mt = (long)(-n2.r); }
    if (cis0(X)) { M = 0; }
    else if (term) { M = Mt; }
    else {
        double rho = cabs_(X) * g;
        if (!(rho < 1.0)) return false;
        double t = 1.0, mx = 1.0;
        M = -1;
        for (long m = 0; m < AP_MAXTERMS; m++) {
            double q = cabs_(cmul(caddr(n1, m), caddr(n2, m))) /
                       (cabs_(caddr(d1, m)) * (m + 1.0)) * rho;
            t *= q;
            if (!(t < 1e300)) return false;
            if (t > mx) mx = t;
            double qn = cabs_(cmul(caddr(n1, m + 1), caddr(n2, m + 1))) /
                        (cabs_(caddr(d1, m + 1)) * (m + 2.0)) * rho;
            if (qn < 1.0 && t * qn / (1.0 - qn) < 1e-18 * mx && m >= 1) { M = m + 1; break; }
        }
        if (M < 0) return false;
    }
    /* outer coefficients T_0..T_M */
    cd Tbuf[256];
    cd* T = M + 1 <= 256 ? Tbuf : (cd*)malloc(sizeof(cd) * (size_t)(M + 1));
    if (!T) return false;
    T[0] = cd_(1.0, 0.0);
    for (long m = 0; m < M; m++)
        T[m + 1] = cmul(T[m], cmul(cdiv(cmul(caddr(n1, m), caddr(n2, m)),
                                        cscal(caddr(d1, m), (double)(m + 1))), X));
    cd S = cd_(0.0, 0.0);
    double Sabs = 0.0;
    bool ok = true;
    if (cis0(Y)) {
        for (long m = M; m >= 0; m--) { S = cadd(S, T[m]); Sabs += cabs_(T[m]); }
    } else if (M == 0) {
        cd G;
        ok = ap_2f1(A0, B0, C0, Y, sig, &G);
        if (ok) { S = G; Sabs = cabs_(G); }
    } else if (ap_separation(kind, Y) < AP_SEP_MIN) {
        for (long m = M; ok && m >= 0; m--) {
            cd G;
            ok = ap_2f1(caddr(A0, sA * (double)m), caddr(B0, sB * (double)m),
                        caddr(C0, sC * (double)m), Y, sig, &G);
            if (!ok) break;
            cd tg = cmul(T[m], G);
            S = cadd(S, tg); Sabs += cabs_(tg);
        }
    } else if (!ap_backward(kind, Y)) {
        cd G0, G1;
        ok = ap_2f1(A0, B0, C0, Y, sig, &G0) &&
             ap_2f1(caddr(A0, sA), caddr(B0, sB), caddr(C0, sC), Y, sig, &G1);
        for (long m = 0; ok && m <= M; m++) {
            cd tg = cmul(T[m], G0);
            S = cadd(S, tg); Sabs += cabs_(tg);
            if (m == M) break;
            if (m + 1 == M) { G0 = G1; continue; }    /* G_(M+1) is not needed */
            cd al, be, ga;
            ap_rec(kind, A0, B0, C0, Y, (int)m, &al, &be, &ga);
            cd G2;
            if (cabs_(ga) <= 1e-8 * (cabs_(al) + cabs_(be))) {
                /* a degenerate step (F4: a + b - c2 + 1 + 2m = 0): directly */
                long m2 = m + 2;
                if (!ap_2f1(caddr(A0, sA * (double)m2), caddr(B0, sB * (double)m2),
                            caddr(C0, sC * (double)m2), Y, sig, &G2)) { ok = false; break; }
            } else {
                G2 = cscal(cdiv(cadd(cmul(al, G0), cmul(be, G1)), ga), -1.0);
            }
            G0 = G1; G1 = G2;
            if (!isfinite(G1.r) || !isfinite(G1.i)) ok = false;
        }
    } else {
        cd GM = cd_(0.0, 0.0), GM1 = GM;   /* G_M, G_(M+1) */
        cd aM = caddr(A0, sA * (double)M), bM = caddr(B0, sB * (double)M), cM = caddr(C0, sC * (double)M);
        cd aM1 = caddr(aM, sA), bM1 = caddr(bM, sB), cM1 = caddr(cM, sC);
        ok = (ap_ser2f1(aM, bM, cM, Y, &GM) || ap_2f1(aM, bM, cM, Y, sig, &GM)) &&
             (ap_ser2f1(aM1, bM1, cM1, Y, &GM1) || ap_2f1(aM1, bM1, cM1, Y, sig, &GM1));
        cd Gm = GM, Gm1 = GM1;
        for (long m = M; ok && m >= 0; m--) {
            cd tg = cmul(T[m], Gm);
            S = cadd(S, tg); Sabs += cabs_(tg);
            if (m == 0) break;
            cd al, be, ga;   /* at index m - 1: al G_(m-1) + be G_m + ga G_(m+1) = 0 */
            ap_rec(kind, A0, B0, C0, Y, (int)(m - 1), &al, &be, &ga);
            cd Gp;
            if (cabs_(al) <= 1e-8 * (cabs_(be) + cabs_(ga))) {
                long m1 = m - 1;
                if (!ap_2f1(caddr(A0, sA * (double)m1), caddr(B0, sB * (double)m1),
                            caddr(C0, sC * (double)m1), Y, sig, &Gp)) { ok = false; break; }
            } else {
                Gp = cscal(cdiv(cadd(cmul(be, Gm), cmul(ga, Gm1)), al), -1.0);
            }
            Gm1 = Gm; Gm = Gp;
            if (!isfinite(Gm.r) || !isfinite(Gm.i)) ok = false;
        }
    }
    if (T != Tbuf) free(T);
    if (!ok || !isfinite(S.r) || !isfinite(S.i)) return false;
    double aS = cabs_(S);
    if (!(aS > 0.0) || Sabs > 4.0 * aS) return false;
    *out = S;
    return true;
}

static bool ap_out(cd v, double* rr, double* ri) {
    if (!isfinite(v.r) || !isfinite(v.i)) return false;
    double big = fmax(fabs(v.r), fabs(v.i));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(v.r) < DBL_MIN) v.r = 0.0;
    if (fabs(v.i) < DBL_MIN) v.i = 0.0;
    *rr = v.r; *ri = v.i;
    return true;
}

/* The reductions to one Gauss function (and the F1 product at x = 1 or y = 1);
 * also the kernel's value for these arguments, so every surface agrees with the
 * builtin. 1: answered (*out set), 0: none applies, -1: no finite value. */
static int ap_special_d(int kind, const cd* p, cd x, cd y, cd* out) {
    cd one = cd_(1.0, 0.0);
    int sx = (x.i == 0.0 && x.r > 1.0) ? -1 : 0, sy = (y.i == 0.0 && y.r > 1.0) ? -1 : 0;
    switch (kind) {
        case 1: {   /* (a, b1, b2, c) */
            if (cis0(p[0])) { *out = one; return 1; }
            if (cis0(x) || cis0(p[1])) return ap_2f1(p[0], p[2], p[3], y, sy, out) ? 1 : -1;
            if (cis0(y) || cis0(p[2])) return ap_2f1(p[0], p[1], p[3], x, sx, out) ? 1 : -1;
            if (x.r == y.r && x.i == y.i) return ap_2f1(p[0], cadd(p[1], p[2]), p[3], x, sx, out) ? 1 : -1;
            if (x.r == 1.0 && x.i == 0.0) {
                cd u, v;
                if (!ap_2f1(p[0], p[1], p[3], one, -1, &u) ||
                    !ap_2f1(p[0], p[2], csub(p[3], p[1]), y, sy, &v)) return -1;
                *out = cmul(u, v); return 1;
            }
            if (y.r == 1.0 && y.i == 0.0) {
                cd u, v;
                if (!ap_2f1(p[0], p[1], csub(p[3], p[2]), x, sx, &u) ||
                    !ap_2f1(p[0], p[2], p[3], one, -1, &v)) return -1;
                *out = cmul(u, v); return 1;
            }
            return 0;
        }
        case 2:     /* (a, b1, b2, c1, c2) */
            if (cis0(p[0]) || (cis0(x) && cis0(y))) { *out = one; return 1; }
            if (cis0(x) || cis0(p[1])) return ap_2f1(p[0], p[2], p[4], y, sy, out) ? 1 : -1;
            if (cis0(y) || cis0(p[2])) return ap_2f1(p[0], p[1], p[3], x, sx, out) ? 1 : -1;
            return 0;
        case 3:     /* (a1, a2, b1, b2, c) */
            if (cis0(x) || cis0(p[0]) || cis0(p[2])) return ap_2f1(p[1], p[3], p[4], y, sy, out) ? 1 : -1;
            if (cis0(y) || cis0(p[1]) || cis0(p[3])) return ap_2f1(p[0], p[2], p[4], x, sx, out) ? 1 : -1;
            return 0;
        default:    /* (a, b, c1, c2) */
            if (cis0(p[0]) || cis0(p[1])) { *out = one; return 1; }
            if (cis0(x)) return ap_2f1(p[0], p[1], p[3], y, sy, out) ? 1 : -1;
            if (cis0(y)) return ap_2f1(p[0], p[1], p[2], x, sx, out) ? 1 : -1;
            return 0;
    }
}

/* True when a denominator parameter is a pole (a non-positive integer). */
static bool ap_pole(int kind, const cd* p) {
    if (kind == 1) return c_nonpos_int(p[3]);
    if (kind == 2) return c_nonpos_int(p[3]) || c_nonpos_int(p[4]);
    if (kind == 3) return c_nonpos_int(p[4]);
    return c_nonpos_int(p[2]) || c_nonpos_int(p[3]);
}

/* Is the value real: real arguments inside the region where the series and
 * its continuation never meet a cut. */
static bool ap_real_region(int kind, const double* re, const double* im) {
    int n = APPELL_NARGS(kind), np = appell_np(kind);
    for (int k = 0; k < n; k++) if (im[k] != 0.0) return false;
    double x = re[np], y = re[np + 1];
    switch (kind) {
        case 1: case 3: return x < 1.0 && y < 1.0;
        case 2: return x < 1.0 && y < 1.0 && x + y < 1.0;
        default:
            /* F4: inside sqrt|x| + sqrt|y| < 1 for x, y >= 0; real for
             * non-positive arguments (the region of the inner 2F1 off its cut
             * through the whole continuation) */
            if (x <= 0.0 && y <= 0.0) return true;
            return sqrt(fabs(x)) + sqrt(fabs(y)) < 1.0;
    }
}

/* Build and evaluate the plan in double. */
static bool ap_machine(int kind, const double* re, const double* im, const AppellPlan* P, cd* out) {
    int np = appell_np(kind);
    const AppellTr* t = appell_tr(kind, P->tr);
    cd p[AP_MAXP], q[AP_MAXP];
    for (int k = 0; k < np; k++) p[k] = cd_(re[k], im[k]);
    cd x = cd_(re[np], im[np]), y = cd_(re[np + 1], im[np + 1]);
    for (int j = 0; j < np; j++) q[j] = ap_newp(t, j, p, np);
    cd X = ap_var(t->vx, x, y), Y = ap_var(t->vy, x, y);
    cd e[6];
    ap_engine_params(kind, q, P->swap, e);
    cd Xo = P->swap ? Y : X, Yi = P->swap ? X : Y;
    if (P->sigY != 0) Yi.i = 0.0;       /* on its cut: the side is in sigY */
    cd S;
    if (!ap_eval_d(kind, e, Xo, Yi, P->sigY, &S)) return false;
    cd pre = cd_(1.0, 0.0);
    cd one = cd_(1.0, 0.0);
    for (int b = 0; b < t->nbase; b++) {
        cd base = t->base[b] == AP_B_1X ? csub(one, x)
                : t->base[b] == AP_B_1Y ? csub(one, y) : csub(csub(one, x), y);
        cd ex = cd_(0.0, 0.0);
        for (int k = 0; k < np; k++)
            if (t->expo[b][k]) ex = cadd(ex, cscal(p[k], (double)t->expo[b][k]));
        pre = cmul(pre, cpow_(base, ex));
    }
    *out = cmul(pre, S);
    return isfinite(out->r) && isfinite(out->i);
}

bool sf_machine_appell(int kind, const double* re, const double* im, double* rr, double* ri) {
    if (kind < 1 || kind > 4) return false;
    int n = APPELL_NARGS(kind), np = appell_np(kind);
    for (int k = 0; k < n; k++) if (!isfinite(re[k]) || !isfinite(im[k])) return false;
    cd p[AP_MAXP];
    for (int k = 0; k < np; k++) p[k] = cd_(re[k], im[k]);
    if (ap_pole(kind, p)) return false;
    cd x = cd_(re[np], im[np]), y = cd_(re[np + 1], im[np + 1]);
    bool real = ap_real_region(kind, re, im);
    cd v;
    int sp = ap_special_d(kind, p, x, y, &v);
    if (sp > 0) {
        if (real) v.i = 0.0;
        if (ap_out(v, rr, ri)) return true;
    }
    if (sp == 0) {
        AppellPlan P;
        if (appell_plan(kind, re, im, 0.97, &P) && ap_machine(kind, re, im, &P, &v)) {
            if (real) v.i = 0.0;
            if (ap_out(v, rr, ri)) return true;
        }
    }
    /* everything the double route cannot certify: arb at 53 bits */
    double vr, vi;
    if (!flint_appell_d(kind, re, im, &vr, &vi)) return false;
    if (real) vi = 0.0;
    return ap_out(cd_(vr, vi), rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtins                                                            */
/* ------------------------------------------------------------------ */

static const char* const AP_NAME[5] = { "", "AppellF1", "AppellF2", "AppellF3", "AppellF4" };

static bool ap_leaf(const Expr* e, double* re, double* im, bool* inexact, bool* cplx) {
    *cplx = false;
    if (sf_machine_leaf(e, re, inexact)) { *im = 0.0; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex &&
        e->data.function.arg_count == 2) {
        bool a, b;
        if (!sf_machine_leaf(e->data.function.args[0], re, &a) ||
            !sf_machine_leaf(e->data.function.args[1], im, &b)) return false;
        *inexact = a || b;
        *cplx = true;
        return true;
    }
    return false;
}

static bool ap_is_int(const Expr* e, int64_t v) {
    return e->type == EXPR_INTEGER && e->data.integer == v;
}

/* Apply the pure-function template `tmpl` to args and evaluate. */
static Expr* ap_apply(const char* tmpl, Expr* const* args, size_t n) {
    Expr* f = parse_expression(tmpl);
    if (!f) return NULL;
    Expr* a[8];
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
    return eval_and_free(expr_new_function(f, a, n));
}

/* e1 - e2 evaluates to exact 0 */
static bool ap_diff_zero(Expr* e1, Expr* e2) {
    Expr* d = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(e1), expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), expr_copy(e2) }, 2) }, 2));
    bool z = d && ap_is_int(d, 0);
    if (d) expr_free(d);
    return z;
}

/* (s)_k as an explicit product of k factors */
static Expr* ap_poch(Expr* s, int64_t k) {
    if (k == 0) return expr_new_integer(1);
    Expr** f = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
    for (int64_t j = 0; j < k; j++)
        f[j] = j == 0 ? expr_copy(s) : expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy(s), expr_new_integer(j) }, 2);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), f, (size_t)k);
    free(f);
    return r;
}

/* F1 with a = -K, K a positive integer: the terminating double sum */
static Expr* ap_f1_poly(Expr* const* A, int64_t K) {
    size_t nt = (size_t)((K + 1) * (K + 2) / 2), it = 0;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * nt);
    for (int64_t k = 0; k <= K; k++) {
        for (int64_t m = 0; m <= k; m++) {
            int64_t nn = k - m;
            /* (a)_k (b1)_m (b2)_n / ((c)_k m! n!) x^m y^n */
            Expr* f[7];
            f[0] = ap_poch(A[0], k);
            f[1] = ap_poch(A[1], m);
            f[2] = ap_poch(A[2], nn);
            f[3] = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ ap_poch(A[3], k), expr_new_integer(-1) }, 2);
            /* Multinomial[m, n] = k!/(m! n!): divided by k! below */
            f[4] = expr_new_function(expr_new_symbol("Multinomial"),
                       (Expr*[]){ expr_new_integer(m), expr_new_integer(nn) }, 2);
            f[5] = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(A[4]), expr_new_integer(m) }, 2);
            f[6] = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(A[5]), expr_new_integer(nn) }, 2);
            Expr* t = expr_new_function(expr_new_symbol(SYM_Times), f, 7);
            Expr* kf = expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ expr_new_function(expr_new_symbol("Factorial"),
                           (Expr*[]){ expr_new_integer(k) }, 1), expr_new_integer(-1) }, 2);
            terms[it++] = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ t, kf }, 2);
        }
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_Plus), terms, it);
    free(terms);
    return eval_and_free(r);
}

/* The exact rules (NULL: none applies). A[] in the head's order. */
static Expr* ap_exact(int kind, Expr* const* A) {
    switch (kind) {
        case 1: {
            Expr *a = A[0], *b1 = A[1], *b2 = A[2], *c = A[3], *x = A[4], *y = A[5];
            if (ap_is_int(a, 0)) return expr_new_integer(1);
            if (ap_is_int(x, 0) || ap_is_int(b1, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a, b2, c, y }, 4);
            if (ap_is_int(y, 0) || ap_is_int(b2, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a, b1, c, x }, 4);
            if (expr_eq(a, c))
                return ap_apply("((1 - #3)^(-#1) (1 - #4)^(-#2))&", (Expr*[]){ b1, b2, x, y }, 4);
            if (expr_eq(x, y))
                return ap_apply("Hypergeometric2F1[#1, #2 + #3, #4, #5]&", (Expr*[]){ a, b1, b2, c, x }, 5);
            if (a->type == EXPR_INTEGER && a->data.integer < 0 && a->data.integer >= -60)
                return ap_f1_poly(A, -a->data.integer);
            {
                Expr* s = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                        (Expr*[]){ expr_copy(b1), expr_copy(b2) }, 2));
                bool eq = s && ap_diff_zero(c, s);
                if (s) expr_free(s);
                if (eq)
                    return ap_apply("((1 - #4)^(-#1) Hypergeometric2F1[#1, #2, #2 + #3, (#5 - #4)/(1 - #4)])&",
                                    (Expr*[]){ a, b1, b2, y, x }, 5);
            }
            if (ap_is_int(x, 1))
                return ap_apply("(Hypergeometric2F1[#1, #2, #4, 1] Hypergeometric2F1[#1, #3, #4 - #2, #5])&",
                                (Expr*[]){ a, b1, b2, c, y }, 5);
            if (ap_is_int(y, 1))
                return ap_apply("(Hypergeometric2F1[#1, #2, #4 - #3, #5] Hypergeometric2F1[#1, #3, #4, 1])&",
                                (Expr*[]){ a, b1, b2, c, x }, 5);
            return NULL;
        }
        case 2: {
            Expr *a = A[0], *b1 = A[1], *b2 = A[2], *c1 = A[3], *c2 = A[4], *x = A[5], *y = A[6];
            if (ap_is_int(a, 0)) return expr_new_integer(1);
            if (ap_is_int(x, 0) || ap_is_int(b1, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a, b2, c2, y }, 4);
            if (ap_is_int(y, 0) || ap_is_int(b2, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a, b1, c1, x }, 4);
            return NULL;
        }
        case 3: {
            Expr *a1 = A[0], *a2 = A[1], *b1 = A[2], *b2 = A[3], *c = A[4], *x = A[5], *y = A[6];
            if (ap_is_int(x, 0) || ap_is_int(a1, 0) || ap_is_int(b1, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a2, b2, c, y }, 4);
            if (ap_is_int(y, 0) || ap_is_int(a2, 0) || ap_is_int(b2, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a1, b1, c, x }, 4);
            return NULL;
        }
        default: {
            Expr *a = A[0], *b = A[1], *c1 = A[2], *c2 = A[3], *x = A[4], *y = A[5];
            if (ap_is_int(a, 0) || ap_is_int(b, 0)) return expr_new_integer(1);
            if (ap_is_int(x, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a, b, c2, y }, 4);
            if (ap_is_int(y, 0))
                return ap_apply("Hypergeometric2F1[#1, #2, #3, #4]&", (Expr*[]){ a, b, c1, x }, 4);
            return NULL;
        }
    }
}

static long ap_bits(Expr* const* a, int n) {
    long bits = 0;
    for (int i = 0; i < n; i++) {
        long b = numeric_min_inexact_bits(a[i]);
        if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    return bits;
}

static Expr* ap_builtin(int kind, Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = (size_t)APPELL_NARGS(kind);
    size_t argc = res->data.function.arg_count;
    if (argc != n) return builtin_arg_error(AP_NAME[kind], argc, n, n);
    Expr** A = res->data.function.args;
    double re[7], im[7];
    bool ix, cx, any_ix = false, any_cx = false, mach = true;
    for (size_t i = 0; i < n; i++) {
        if (!ap_leaf(A[i], &re[i], &im[i], &ix, &cx)) { mach = false; break; }
        any_ix |= ix; any_cx |= cx;
    }
    if (mach && any_ix) {
        double rr, ri;
        if (sf_machine_appell(kind, re, im, &rr, &ri))
            return (!any_cx && ri == 0.0) ? expr_new_real(rr) : sf_machine_complex_result(rr, ri);
    }
    bool numeric = true;
    for (size_t i = 0; i < n; i++) numeric &= expr_is_numeric_quantity(A[i]);
    long bits = ap_bits(A, (int)n);
    if (!numeric || bits == 0) return ap_exact(kind, A);
    Expr* r = flint_num_appell(kind, (const Expr* const*)A, (int)n);
    if (r && bits <= 53) r = numeric_lower_mpfr_to_machine(r);
    return r;
}

static Expr* builtin_appellf1(Expr* res) { return ap_builtin(1, res); }
static Expr* builtin_appellf2(Expr* res) { return ap_builtin(2, res); }
static Expr* builtin_appellf3(Expr* res) { return ap_builtin(3, res); }
static Expr* builtin_appellf4(Expr* res) { return ap_builtin(4, res); }

void appell_init(void) {
    Expr* (*fn[5])(Expr*) = { NULL, builtin_appellf1, builtin_appellf2, builtin_appellf3,
                              builtin_appellf4 };
    for (int k = 1; k <= 4; k++) {
        symtab_add_builtin(AP_NAME[k], fn[k]);
        symtab_get_def(AP_NAME[k])->attributes |=
            (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    }
    symtab_set_docstring("AppellF1",
        "AppellF1[a, b1, b2, c, x, y]\n"
        "\tgives the Appell hypergeometric function of two variables\n"
        "\tSum[(a)_(m+n) (b1)_m (b2)_n/((c)_(m+n) m! n!) x^m y^n], continued\n"
        "\tanalytically from |x|, |y| < 1 (branch cuts x, y in [1, Infinity), from below).\n"
        "Numeric at machine and arbitrary precision, real and complex; threads over\n"
        "packed arrays and compiles. Reduces to Hypergeometric2F1 for x = 0, y = 0,\n"
        "x = y, a = c, c = b1 + b2; a polynomial for a negative integer a.\n"
        "D[AppellF1[a, b1, b2, c, x, y], x] = a b1/c AppellF1[a+1, b1+1, b2, c+1, x, y].");
    symtab_set_docstring("AppellF2",
        "AppellF2[a, b1, b2, c1, c2, x, y]\n"
        "\tgives the Appell hypergeometric function\n"
        "\tSum[(a)_(m+n) (b1)_m (b2)_n/((c1)_m (c2)_n m! n!) x^m y^n] (|x| + |y| < 1),\n"
        "\tcontinued analytically.\n"
        "Numeric at machine and arbitrary precision, real and complex; threads over\n"
        "packed arrays and compiles. D[AppellF2[...], x] = a b1/c1 AppellF2[a+1, b1+1, b2, c1+1, c2, x, y].");
    symtab_set_docstring("AppellF3",
        "AppellF3[a1, a2, b1, b2, c, x, y]\n"
        "\tgives the Appell hypergeometric function\n"
        "\tSum[(a1)_m (a2)_n (b1)_m (b2)_n/((c)_(m+n) m! n!) x^m y^n] (|x|, |y| < 1),\n"
        "\tcontinued analytically.\n"
        "Numeric at machine and arbitrary precision, real and complex; threads over\n"
        "packed arrays and compiles. D[AppellF3[...], x] = a1 b1/c AppellF3[a1+1, a2, b1+1, b2, c+1, x, y].");
    symtab_set_docstring("AppellF4",
        "AppellF4[a, b, c1, c2, x, y]\n"
        "\tgives the Appell hypergeometric function\n"
        "\tSum[(a)_(m+n) (b)_(m+n)/((c1)_m (c2)_n m! n!) x^m y^n] (Sqrt[|x|] + Sqrt[|y|] < 1),\n"
        "\tcontinued analytically.\n"
        "Numeric at machine and arbitrary precision, real and complex; threads over\n"
        "packed arrays and compiles. D[AppellF4[...], x] = a b/c1 AppellF4[a+1, b+1, c1+1, c2, x, y].");
}
