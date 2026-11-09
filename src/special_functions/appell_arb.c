/* Mathilda -- the Appell functions F1..F4 in ball arithmetic (arb / acb).
 *
 * The algorithm of appell.c carried out with acb: the same plan (transformation,
 * outer variable, side of the inner cut, from appell_plan on the midpoints),
 * the outer coefficients T_m exactly, the inner Gauss functions G_m from two
 * acb_hypgeom_2f1 values and the three-term recurrence in its stable direction
 * (at extra working precision, since a ball recurrence overestimates the
 * propagated error), and the truncation error of the outer series bounded by
 * |T_(M+1) G_M| q/(1 - q) times a safety factor, q a bound on the ratio
 * |T_(m+1) G_(m+1)/(T_m G_m)| for m > M from the parameters and the inner growth
 * rate. That bound is asymptotic in G (the growth of G_m is that of its
 * recurrence's dominant root, not a proven majorant), which is the one
 * non-rigorous ingredient; the precision ladder around this function then
 * certifies the digits it emits.
 *
 * Where no expansion has ratio <= 0.995 (F1 with both variables large and
 * |1 - y/x|, |1 - x/y| >= 1, a region no Euler substitution reaches), F1 is the
 * Euler integral
 *   Gamma(c)/(Gamma(a) Gamma(c-a)) Int_0^1 t^(a-1) (1-t)^(c-a-1) (1-x t)^-b1 (1-y t)^-b2 dt
 * by acb_calc_integrate (Re c > Re a > 0, x and y off their cuts). F2, F3, F4
 * there are left unevaluated.
 */
#include "appell.h"

#ifdef USE_FLINT
#include <math.h>
#include <stdlib.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_calc.h>
#include <flint/acb_poly.h>

/* rel = the largest relative radius among four balls (an exact ball adds 0) */
static void ap_rel_rad(mag_t rel, const acb_t a, const acb_t b, const acb_t c, const acb_t d) {
    const acb_struct* g[4] = { a, b, c, d };
    mag_t r, l;
    mag_init(r); mag_init(l);
    mag_zero(rel);
    for (int k = 0; k < 4; k++) {
        mag_max(r, arb_radref(acb_realref(g[k])), arb_radref(acb_imagref(g[k])));
        if (mag_is_zero(r)) continue;
        mag_mul_2exp_si(r, r, 1);
        acb_get_mag_lower(l, g[k]);
        if (mag_is_zero(l)) { mag_inf(rel); break; }
        mag_div(r, r, l);
        if (mag_cmp(r, rel) > 0) mag_set(rel, r);
    }
    mag_clear(r); mag_clear(l);
}

/* The Gauss series summed directly, |z| < 0.99, with a rigorous tail bound:
 * past term k every ratio |t_(j+1)/t_j| is at most
 *   q = |z| (1 + |a - 1|/(k + 1)) (1 + |b - c|/(k + Re c)) (or a, b exchanged),
 * so the tail is below |t_k| q/(1 - q). 0 when it cannot (no such k soon). */
static int ap_acb_ser(acb_t r, acb_srcptr a, acb_srcptr b, acb_srcptr c, acb_srcptr z,
                      double za, slong prec) {
    double am1 = hypot(arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR) - 1.0,
                       arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR)) * 1.001 + 1e-12;
    double bmc = hypot(arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR) -
                       arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR),
                       arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR) -
                       arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR)) * 1.001 + 1e-12;
    double rc = arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR) - 1e-9;
    /* the other pairing, (a + j)/(c + j) and (b + j)/(j + 1): tight when a and c
     * are both large (the shifted parameters of F1's start values) */
    double amc = hypot(arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR) -
                       arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR),
                       arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR) -
                       arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR)) * 1.001 + 1e-12;
    double bm1 = hypot(arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR) - 1.0,
                       arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR)) * 1.001 + 1e-12;
    acb_t t, u, v, s;
    acb_init(t); acb_init(u); acb_init(v); acb_init(s);
    acb_one(t); acb_one(s);
    int ok = 0;
    mag_t mt, ms, mb, err, rel;
    mag_init(mt); mag_init(ms); mag_init(mb); mag_init(err); mag_init(rel);
    ap_rel_rad(rel, a, b, c, z);
    mag_set_ui_2exp_si(mb, 8, -prec);
    mag_add(rel, rel, mb);
    for (slong k = 0; k < 100000; k++) {
        acb_add_si(u, a, k, prec); acb_add_si(v, b, k, prec); acb_mul(u, u, v, prec);
        acb_add_si(v, c, k, prec); acb_mul_si(v, v, k + 1, prec); acb_div(u, u, v, prec);
        acb_mul(u, u, z, prec);
        acb_mul(t, t, u, prec);
        /* a ball product of complex numbers overestimates the radius by up to
         * sqrt 2 a step (rectangular wrapping), so the term runs on midpoints
         * and its error -- at most (k + 2)(8 ulp + the inputs' relative radius)
         * relative -- is accumulated separately */
        acb_get_mid(t, t);
        acb_add(s, s, t, prec);
        acb_get_mag(mt, t);
        mag_mul_ui(mt, mt, (ulong)(k + 2));
        mag_mul(mt, mt, rel);
        mag_add(err, err, mt);
        if (acb_is_zero(t)) { ok = 1; break; }          /* terminated */
        slong j = k + 1;                                /* the next ratio index */
        if (j + rc <= 0.0) continue;
        double q1 = za * (1.0 + am1 / (j + 1.0)) * (1.0 + bmc / (j + rc));
        double q2 = za * (1.0 + bm1 / (j + 1.0)) * (1.0 + amc / (j + rc));
        double q = q1 < q2 ? q1 : q2;
        if (!(q < 0.999)) continue;
        acb_get_mag(mt, t);
        mag_set_d(mb, q / (1.0 - q));
        mag_mul(mt, mt, mb);
        acb_get_mag_lower(ms, s);
        mag_mul_2exp_si(ms, ms, -prec - 4);
        if (mag_cmp(mt, ms) <= 0) { acb_add_error_mag(s, mt); ok = 1; break; }
    }
    acb_add_error_mag(s, err);
    if (ok) acb_swap(r, s);
    mag_clear(mt); mag_clear(ms); mag_clear(mb); mag_clear(err); mag_clear(rel);
    acb_clear(t); acb_clear(u); acb_clear(v); acb_clear(s);
    return ok && acb_is_finite(r);
}

/* 2F1 with the inner variable on its cut from below (sig <= 0, acb's
 * convention) or above (sig > 0: conjugate of the value at conjugate
 * parameters). Inside |z| < 0.99 the direct series (ap_acb_ser): there
 * acb_hypgeom_2f1 can return an indeterminate ball for a large c (e.g.
 * 2F1(-0.0721, 1.15; 102.33; 0.514 - 0.614 i) in FLINT 3.6), and the large
 * parameters of the backward start values are exactly where the series is
 * short. */
static void ap_acb_2f1(acb_t r, acb_srcptr a, acb_srcptr b, acb_srcptr c, acb_srcptr z,
                       int sig, slong prec) {
    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 0) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 0) < 0) {
        double zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
        double zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
        double za = hypot(zr, zi);
        if (za < 0.99 && ap_acb_ser(r, a, b, c, z, za * 1.0001 + 1e-15, prec)) return;
    }
    {
        /* Pfaff, 2F1(a, b; c; z) = (1 - z)^-a 2F1(a, c - b; c; z/(z - 1)), valid
         * with principal branches off the cut, when z/(z - 1) is inside 0.99 */
        double zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
        double zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
        double w = hypot(zr, zi) / hypot(zr - 1.0, zi);
        if (w < 0.99 && zr < 0.5) {
            acb_t W, B, P;
            acb_init(W); acb_init(B); acb_init(P);
            acb_sub_si(P, z, 1, prec);
            acb_div(W, z, P, prec);
            acb_sub(B, c, b, prec);
            int ok = ap_acb_ser(r, a, B, c, W, w * 1.0001 + 1e-15, prec);
            if (ok) {
                acb_neg(P, P);                        /* 1 - z */
                acb_neg(B, a);
                acb_pow(P, P, B, prec);
                acb_mul(r, r, P, prec);
            }
            acb_clear(W); acb_clear(B); acb_clear(P);
            if (ok && acb_is_finite(r)) return;
        }
    }
    if (sig > 0) {
        acb_t A, B, C;
        acb_init(A); acb_init(B); acb_init(C);
        acb_conj(A, a); acb_conj(B, b); acb_conj(C, c);
        acb_hypgeom_2f1(r, A, B, C, z, 0, prec);
        acb_conj(r, r);
        acb_clear(A); acb_clear(B); acb_clear(C);
        return;
    }
    acb_hypgeom_2f1(r, a, b, c, z, 0, prec);
}

/* alpha P_m + beta P_(m+1) + gamma P_(m+2) = 0, the recurrences of appell.c */
static void ap_rec_acb(int kind, const acb_t A0, const acb_t B0, const acb_t C0, const acb_t z,
                       slong m, acb_t al, acb_t be, acb_t ga, slong prec) {
    acb_t A, B, C, t, u, v;
    acb_init(A); acb_init(B); acb_init(C); acb_init(t); acb_init(u); acb_init(v);
    if (kind == 1) {
        acb_add_si(A, A0, m + 1, prec);                    /* A + 1 */
        acb_add_si(C, C0, m, prec);                        /* C */
        acb_add_si(t, C, 1, prec);                         /* C + 1 */
        acb_sub(u, t, B0, prec);                           /* C + 1 - b */
        acb_mul(ga, A, u, prec); acb_mul(ga, ga, z, prec);
        acb_sub(u, A, B0, prec); acb_mul(u, u, z, prec); acb_add(u, u, C, prec);
        acb_mul(be, t, u, prec); acb_neg(be, be);
        acb_mul(al, C, t, prec);
    } else if (kind == 2) {
        acb_add_si(A, A0, m + 1, prec);                    /* A + 1 */
        acb_sub(al, C0, A, prec);
        acb_mul_2exp_si(t, A, 1); acb_sub(t, t, C0, prec);
        acb_sub(u, B0, A, prec); acb_mul(u, u, z, prec);
        acb_add(be, t, u, prec);
        acb_sub_si(u, z, 1, prec); acb_mul(ga, A, u, prec);
    } else if (kind == 3) {
        acb_add_si(C, C0, m, prec);
        acb_add_si(t, C, 1, prec);                         /* C + 1 */
        acb_sub_si(u, z, 1, prec);
        acb_mul(al, t, C, prec); acb_mul(al, al, u, prec);
        acb_mul_2exp_si(u, C, 1); acb_add_si(u, u, 1, prec);
        acb_sub(u, u, A0, prec); acb_sub(u, u, B0, prec); acb_mul(u, u, z, prec);
        acb_sub(u, C, u, prec);
        acb_mul(be, t, u, prec);
        acb_sub(u, t, A0, prec); acb_sub(v, t, B0, prec);
        acb_mul(ga, u, v, prec); acb_mul(ga, ga, z, prec);
    } else {
        acb_t ABc, zm1, Q, w;
        acb_init(ABc); acb_init(zm1); acb_init(Q); acb_init(w);
        acb_add_si(A, A0, m, prec); acb_add_si(B, B0, m, prec);
        const acb_struct* c = C0;
        acb_sub_si(zm1, z, 1, prec);
        acb_add(ABc, A, B, prec); acb_sub(ABc, ABc, c, prec);
        acb_add_si(t, A, 1, prec); acb_add_si(u, B, 1, prec);
        acb_mul(ga, t, u, prec);
        acb_mul(w, zm1, zm1, prec); acb_mul(ga, ga, w, prec);
        acb_add_si(w, ABc, 1, prec); acb_mul(ga, ga, w, prec);
        acb_sub(t, A, c, prec); acb_add_si(t, t, 1, prec);
        acb_sub(u, B, c, prec); acb_add_si(u, u, 1, prec);
        acb_mul(al, t, u, prec);
        acb_add_si(w, ABc, 3, prec); acb_mul(al, al, w, prec);
        /* Q = A^2 z + 2AB - Acz - Ac + 2Az + 2A + B^2 z - Bcz - Bc + 2Bz + 2B
         *     + c^2 - cz - 3c + z + 2
         *   = z (A^2 + B^2 - (A + B) c + 2 (A + B) - c + 1)
         *     + 2AB - (A + B) c + 2 (A + B) + c^2 - 3c + 2 */
        acb_t S;
        acb_init(S);
        acb_add(S, A, B, prec);
        acb_mul(t, A, A, prec); acb_mul(u, B, B, prec); acb_add(t, t, u, prec);
        acb_mul(u, S, c, prec); acb_sub(t, t, u, prec);
        acb_mul_2exp_si(u, S, 1); acb_add(t, t, u, prec);
        acb_sub(t, t, c, prec); acb_add_si(t, t, 1, prec);
        acb_mul(Q, t, z, prec);
        acb_mul(t, A, B, prec); acb_mul_2exp_si(t, t, 1);
        acb_mul(u, S, c, prec); acb_sub(t, t, u, prec);
        acb_mul_2exp_si(u, S, 1); acb_add(t, t, u, prec);
        acb_mul(u, c, c, prec); acb_add(t, t, u, prec);
        acb_mul_si(u, c, 3, prec); acb_sub(t, t, u, prec);
        acb_add_si(t, t, 2, prec);
        acb_add(Q, Q, t, prec);
        acb_add_si(w, ABc, 2, prec);
        acb_mul(be, w, Q, prec); acb_neg(be, be);
        acb_clear(S);
        acb_clear(ABc); acb_clear(zm1); acb_clear(Q); acb_clear(w);
    }
    acb_clear(A); acb_clear(B); acb_clear(C); acb_clear(t); acb_clear(u); acb_clear(v);
}

static double ap_absd(double r, double i) { return hypot(r, i); }

static int ap_nonpos_int_acb(const acb_t a) {
    if (!acb_is_exact(a) || !acb_is_real(a)) return 0;
    if (!arb_is_int(acb_realref(a))) return 0;
    return arf_sgn(arb_midref(acb_realref(a))) <= 0;
}

/* the variable expressions of the transformation table */
static void ap_var_acb(acb_t r, int v, const acb_t x, const acb_t y, slong prec) {
    acb_t n, d;
    acb_init(n); acb_init(d);
    switch (v) {
        case AP_V_X: acb_set(r, x); break;
        case AP_V_Y: acb_set(r, y); break;
        case AP_V_XM:  acb_sub_si(d, x, 1, prec); acb_div(r, x, d, prec); break;
        case AP_V_YM:  acb_sub_si(d, y, 1, prec); acb_div(r, y, d, prec); break;
        case AP_V_YX1: acb_sub(n, y, x, prec); acb_neg(d, x); acb_add_si(d, d, 1, prec);
                       acb_div(r, n, d, prec); break;
        case AP_V_XY1: acb_sub(n, x, y, prec); acb_neg(d, y); acb_add_si(d, d, 1, prec);
                       acb_div(r, n, d, prec); break;
        case AP_V_YXY: acb_sub(n, y, x, prec); acb_sub_si(d, y, 1, prec); acb_div(r, n, d, prec); break;
        case AP_V_XYX: acb_sub(n, x, y, prec); acb_sub_si(d, x, 1, prec); acb_div(r, n, d, prec); break;
        case AP_V_Y1X: acb_neg(d, x); acb_add_si(d, d, 1, prec); acb_div(r, y, d, prec); break;
        case AP_V_X1Y: acb_neg(d, y); acb_add_si(d, d, 1, prec); acb_div(r, x, d, prec); break;
        case AP_V_XS:  acb_add(d, x, y, prec); acb_sub_si(d, d, 1, prec); acb_div(r, x, d, prec); break;
        case AP_V_YS:  acb_add(d, x, y, prec); acb_sub_si(d, d, 1, prec); acb_div(r, y, d, prec); break;
        default: acb_indeterminate(r);
    }
    acb_clear(n); acb_clear(d);
}

/* |p| <= 1e-8 (|q| + |r|): the recurrence step dividing by p degenerates
 * (F4 at a + b - c2 + 1 + 2m = 0), and the value is computed directly */
static int ap_degenerate(acb_srcptr p, acb_srcptr q, acb_srcptr r) {
    double a = hypot(arf_get_d(arb_midref(acb_realref(p)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(p)), ARF_RND_NEAR));
    double b = hypot(arf_get_d(arb_midref(acb_realref(q)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(q)), ARF_RND_NEAR));
    double c = hypot(arf_get_d(arb_midref(acb_realref(r)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(r)), ARF_RND_NEAR));
    return a <= 1e-8 * (b + c) || acb_contains_zero(p);
}

/* r0 = the larger relative radius of two start values (infinite when a ball
 * contains zero) */
static void ap_start_rad(mag_t r0, const acb_t g0, const acb_t g1) {
    mag_t a, b;
    mag_init(a); mag_init(b);
    mag_zero(r0);
    const acb_struct* g[2] = { g0, g1 };
    for (int k = 0; k < 2; k++) {
        mag_max(a, arb_radref(acb_realref(g[k])), arb_radref(acb_imagref(g[k])));
        mag_mul_2exp_si(a, a, 1);
        if (mag_is_zero(a)) continue;
        acb_get_mag_lower(b, g[k]);
        if (mag_is_zero(b)) { mag_inf(r0); break; }
        mag_div(a, a, b);
        if (mag_cmp(a, r0) > 0) mag_set(r0, a);
    }
    mag_clear(a); mag_clear(b);
}

/* The outer series S = Sum_m T_m G_m in ball arithmetic. e[] = (n1, n2, d1,
 * A0, B0, C0); returns 0 when it cannot (too many terms). */
static int ap_eval_acb(acb_t S, int kind, acb_srcptr e, const acb_t X, const acb_t Y,
                       int sig, slong prec) {
    int sA = (kind == 1 || kind == 2 || kind == 4), sB = (kind == 4), sC = (kind == 1 || kind == 3);
    double er[6], ei[6];
    for (int k = 0; k < 6; k++) {
        er[k] = arf_get_d(arb_midref(acb_realref(e + k)), ARF_RND_NEAR);
        ei[k] = arf_get_d(arb_midref(acb_imagref(e + k)), ARF_RND_NEAR);
    }
    double Xa = ap_absd(arf_get_d(arb_midref(acb_realref(X)), ARF_RND_NEAR),
                        arf_get_d(arb_midref(acb_imagref(X)), ARF_RND_NEAR));
    double yr = arf_get_d(arb_midref(acb_realref(Y)), ARF_RND_NEAR);
    double yi = arf_get_d(arb_midref(acb_imagref(Y)), ARF_RND_NEAR);
    double g; int back;
    appell_engine_info(kind, yr, yi, &g, &back);
    /* terminating outer series */
    slong M = -1;
    int term = 0;
    for (int k = 0; k < 2; k++)
        if (ap_nonpos_int_acb(e + k)) {
            slong t = -(slong)arf_get_si(arb_midref(acb_realref(e + k)), ARF_RND_NEAR);
            if (!term || t < M) M = t;
            term = 1;
        }
    double rho = Xa * g, qtail = 0.0;
    if (acb_is_zero(X)) { M = 0; term = 1; }
    if (!term) {
        if (!(rho < 1.0)) return 0;
        double tol = ldexp(1.0, -(int)prec - 12), t = 1.0, mx = 1.0, ls = 0.0;
        for (slong m = 0; m < 200000; m++) {
            double q = ap_absd(er[0] + m, ei[0]) * ap_absd(er[1] + m, ei[1]) /
                       (ap_absd(er[2] + m, ei[2]) * (m + 1.0)) * rho;
            t *= q;
            if (t > 1e250) { t *= 1e-250; mx *= 1e-250; ls += 250.0; }
            if (t > mx) mx = t;
            double qn = ap_absd(er[0] + m + 1, ei[0]) * ap_absd(er[1] + m + 1, ei[1]) /
                        (ap_absd(er[2] + m + 1, ei[2]) * (m + 2.0)) * rho;
            if (m >= 1 && qn < 1.0 && t * qn / (1.0 - qn) < tol * mx) {
                M = m + 1; qtail = qn; break;
            }
        }
        (void)ls;
        if (M < 0) return 0;
    }
    /* The recurrence runs on MIDPOINTS (each value's radius dropped after its
     * step): it is run in its stable direction, where rounding errors grow at
     * most linearly in m, while a ball recurrence would compound the
     * overestimated radius at every step (~1 bit a step near |Y| = 1). The
     * dropped rounding error is restored below as a bound of (M + 2) ulps of
     * the working precision on each summand. Near |Y| = 1 (two roots of nearly
     * equal modulus) the stable direction is only polynomially stable, which
     * the 8 log2 M guard bits absorb. */
    slong lgM = 1;
    while ((1L << lgM) < (long)(M + 2) && lgM < 60) lgM++;
    slong wp = prec + 30 + 8 * lgM;
    acb_ptr T = _acb_vec_init(M + 2);
    acb_one(T + 0);
    acb_t t, u, al, be, ga, G0, G1, G2;
    acb_init(t); acb_init(u); acb_init(al); acb_init(be); acb_init(ga);
    acb_init(G0); acb_init(G1); acb_init(G2);
    for (slong m = 0; m <= M; m++) {
        acb_add_si(t, e + 0, m, wp); acb_add_si(u, e + 1, m, wp); acb_mul(t, t, u, wp);
        acb_add_si(u, e + 2, m, wp); acb_mul_si(u, u, m + 1, wp); acb_div(t, t, u, wp);
        acb_mul(t, t, X, wp);
        acb_mul(T + m + 1, T + m, t, wp);
        acb_get_mid(T + m + 1, T + m + 1);    /* see ap_acb_ser: no wrapping */
    }
    acb_zero(S);
    int ok = 1;
    acb_t a_, b_, c_;
    acb_init(a_); acb_init(b_); acb_init(c_);
    mag_t Sabs, mt, r0, rX;
    mag_init(Sabs); mag_init(mt); mag_init(r0); mag_init(rX);
    ap_rel_rad(rX, X, e + 0, e + 1, e + 2);
    if (acb_is_zero(Y)) {
        for (slong m = M; m >= 0; m--) acb_add(S, S, T + m, wp);
        acb_one(G0);
    } else if (!back || M == 0) {
        ap_acb_2f1(G0, e + 3, e + 4, e + 5, Y, sig, wp);
        acb_add_si(a_, e + 3, sA, wp); acb_add_si(b_, e + 4, sB, wp); acb_add_si(c_, e + 5, sC, wp);
        if (M > 0) ap_acb_2f1(G1, a_, b_, c_, Y, sig, wp); else acb_set(G1, G0);
        ap_start_rad(r0, G0, G1);
        acb_get_mid(G0, G0); acb_get_mid(G1, G1);
        for (slong m = 0; m <= M; m++) {
            acb_mul(t, T + m, G0, wp); acb_add(S, S, t, wp);
            acb_get_mag(mt, t); mag_add(Sabs, Sabs, mt);
            if (m == M) break;
            if (m + 1 == M) { acb_swap(G0, G1); continue; }   /* G_(M+1) is not needed */
            ap_rec_acb(kind, e + 3, e + 4, e + 5, Y, m, al, be, ga, wp);
            if (ap_degenerate(ga, al, be)) {
                acb_add_si(a_, e + 3, sA * (m + 2), wp); acb_add_si(b_, e + 4, sB * (m + 2), wp);
                acb_add_si(c_, e + 5, sC * (m + 2), wp);
                ap_acb_2f1(G2, a_, b_, c_, Y, sig, wp);
            } else {
                acb_mul(G2, al, G0, wp); acb_addmul(G2, be, G1, wp);
                acb_div(G2, G2, ga, wp); acb_neg(G2, G2);
            }
            acb_get_mid(G2, G2);
            acb_swap(G0, G1); acb_swap(G1, G2);
            if (!acb_is_finite(G1)) { ok = 0; break; }
        }
        /* G0 now holds G_M (the last summed) */
    } else {
        acb_add_si(a_, e + 3, sA * M, wp); acb_add_si(b_, e + 4, sB * M, wp);
        acb_add_si(c_, e + 5, sC * M, wp);
        ap_acb_2f1(G0, a_, b_, c_, Y, sig, wp);                    /* G_M */
        acb_add_si(a_, a_, sA, wp); acb_add_si(b_, b_, sB, wp); acb_add_si(c_, c_, sC, wp);
        ap_acb_2f1(G1, a_, b_, c_, Y, sig, wp);                    /* G_(M+1) */
        ap_start_rad(r0, G0, G1);
        acb_get_mid(G0, G0); acb_get_mid(G1, G1);
        acb_set(G2, G0);                                         /* keep G_M for the tail */
        for (slong m = M; m >= 0; m--) {
            acb_mul(u, T + m, G0, wp); acb_add(S, S, u, wp);
            acb_get_mag(mt, u); mag_add(Sabs, Sabs, mt);
            if (m == 0) break;
            ap_rec_acb(kind, e + 3, e + 4, e + 5, Y, m - 1, al, be, ga, wp);
            if (ap_degenerate(al, be, ga)) {
                acb_add_si(a_, e + 3, sA * (m - 1), wp); acb_add_si(b_, e + 4, sB * (m - 1), wp);
                acb_add_si(c_, e + 5, sC * (m - 1), wp);
                ap_acb_2f1(t, a_, b_, c_, Y, sig, wp);
            } else {
                acb_mul(t, be, G0, wp); acb_addmul(t, ga, G1, wp);
                acb_div(t, t, al, wp); acb_neg(t, t);
            }
            acb_get_mid(t, t);
            acb_swap(G1, G0); acb_swap(G0, t);
            if (!acb_is_finite(G0)) { ok = 0; break; }
        }
        acb_swap(G0, G2);                                          /* G_M */
    }
    if (ok && acb_is_finite(S)) {
        /* the dropped rounding error of the recurrence: (M + 2) 2^-(wp - 4)
         * relative to each summand */
        mag_mul_ui(mt, Sabs, (ulong)(M + 2));
        mag_mul_2exp_si(mt, mt, -(wp - 4));
        acb_add_error_mag(S, mt);
        /* the outer coefficients' propagated input radius, (m + 1) rX */
        mag_mul_ui(mt, Sabs, (ulong)(M + 2));
        mag_mul(mt, mt, rX);
        mag_mul_ui(mt, mt, 2);
        acb_add_error_mag(S, mt);
        /* the start values' own error, carried by the stable recurrence to a
         * bounded multiple of the same relative error in every G_m */
        mag_mul(mt, Sabs, r0);
        mag_mul_ui(mt, mt, 4);
        acb_add_error_mag(S, mt);
    }
    if (ok && !term && acb_is_finite(S)) {
        /* tail: |T_(M+1) G_M| q/(1 - q), q the ratio bound with margins */
        double qb = qtail > rho ? qtail : rho;
        double q = qb + 0.25 * (1.0 - qb);
        if (q < 1.0) {
            mag_t m1, m2;
            mag_init(m1); mag_init(m2);
            acb_get_mag(m1, T + M + 1);
            acb_get_mag(m2, G0);
            mag_mul(m1, m1, m2);
            mag_mul_ui(m1, m1, 4);
            mag_set_d(m2, q / (1.0 - q));
            mag_mul(m1, m1, m2);
            acb_add_error_mag(S, m1);
            mag_clear(m1); mag_clear(m2);
        } else {
            ok = 0;
        }
    }
    acb_clear(a_); acb_clear(b_); acb_clear(c_);
    mag_clear(Sabs); mag_clear(mt); mag_clear(r0); mag_clear(rX);
    acb_clear(t); acb_clear(u); acb_clear(al); acb_clear(be); acb_clear(ga);
    acb_clear(G0); acb_clear(G1); acb_clear(G2);
    _acb_vec_clear(T, M + 2);
    if (!ok || !acb_is_finite(S)) return 0;
    return 1;
}

/* The F1 Euler integral
 *   Gamma(c)/(Gamma(a) Gamma(c-a)) Int_0^1 t^(a-1) (1-t)^(c-a-1) (1-x t)^-b1 (1-y t)^-b2 dt,
 * Re c > Re a > 0. The integrand is a product of four factors w_j(t)^p_j,
 * w_j = u_j + v_j t. Its algebraic singular points on the path -- t = 0, t = 1
 * and t = 1/x for an x on its cut (real > 1, from below: there 1 - x t turns
 * negative with arg +pi, so beyond 1/x the factor is e^(i pi p) (x t - 1)^p)
 * -- split [0, 1] into segments. Near each singular point e the integral over
 * [e, e + delta] (or [e - delta, e]) is done exactly from the Taylor series of
 * the regular factors: Int_0^delta s^lambda Sum h_k s^k ds = Sum h_k
 * delta^(lambda+k+1)/(lambda+k+1), delta a quarter of the distance to the next
 * singularity (so the series gains two bits a term); the interior, where the
 * integrand is analytic, goes to acb_calc_integrate. */
#define AP_NFAC 4
typedef struct {
    acb_struct u[AP_NFAC], v[AP_NFAC], p[AP_NFAC];
    int flip[AP_NFAC];      /* the factor is e^(i pi p) (x t - 1)^p here */
} ap_euler_seg;

static void ap_seg_init(ap_euler_seg* S) {
    for (int j = 0; j < AP_NFAC; j++) { acb_init(S->u + j); acb_init(S->v + j); acb_init(S->p + j); }
}
static void ap_seg_clear(ap_euler_seg* S) {
    for (int j = 0; j < AP_NFAC; j++) { acb_clear(S->u + j); acb_clear(S->v + j); acb_clear(S->p + j); }
}

/* the factors on a segment: flips[k] for the x (k = 0) / y (k = 1) factors */
static void ap_seg_set(ap_euler_seg* S, acb_srcptr p, const int* flips, slong prec) {
    acb_zero(S->u + 0); acb_one(S->v + 0); acb_sub_si(S->p + 0, p + 0, 1, prec);    /* t^(a-1) */
    acb_one(S->u + 1); acb_set_si(S->v + 1, -1);                                    /* (1-t)^(c-a-1) */
    acb_sub(S->p + 1, p + 3, p + 0, prec); acb_sub_si(S->p + 1, S->p + 1, 1, prec);
    S->flip[0] = S->flip[1] = 0;
    for (int k = 0; k < 2; k++) {
        int j = 2 + k;
        S->flip[j] = flips[k];
        if (flips[k]) { acb_set_si(S->u + j, -1); acb_set(S->v + j, p + 4 + k); }
        else { acb_one(S->u + j); acb_neg(S->v + j, p + 4 + k); }
        acb_neg(S->p + j, p + 1 + k);
    }
}

static int ap_euler_f(acb_ptr out, const acb_t t, void* param, slong order, slong prec) {
    const ap_euler_seg* S = (const ap_euler_seg*)param;
    if (order > 1) flint_abort();
    acb_t w, r;
    acb_init(w); acb_init(r);
    acb_one(r);
    for (int j = 0; j < AP_NFAC; j++) {
        acb_mul(w, S->v + j, t, prec); acb_add(w, w, S->u + j, prec);
        acb_pow_analytic(w, w, S->p + j, order == 1, prec);
        acb_mul(r, r, w, prec);
    }
    acb_swap(out, r);
    acb_clear(w); acb_clear(r);
    return 0;
}

/* Int over [e, e + dir delta] of the integrand, by the series at its singular
 * point e (factor `js` vanishes there). */
static void ap_euler_end(acb_t res, const ap_euler_seg* S, int js, const acb_t e, int dir,
                         const arb_t delta, slong N, slong prec) {
    acb_poly_t h, f, g;
    acb_poly_init(h); acb_poly_init(f); acb_poly_init(g);
    acb_t w0, w1, lam, d, t, sum, k1;
    acb_init(w0); acb_init(w1); acb_init(lam); acb_init(d); acb_init(t); acb_init(sum); acb_init(k1);
    acb_poly_one(h);
    acb_t coef;
    acb_init(coef);
    acb_one(coef);
    for (int j = 0; j < AP_NFAC; j++) {
        /* w_j(e + dir s) = w_j(e) + dir v_j s */
        acb_mul(w0, S->v + j, e, prec); acb_add(w0, w0, S->u + j, prec);
        acb_set(w1, S->v + j); if (dir < 0) acb_neg(w1, w1);
        if (j == js) {
            /* (dir v s)^p, dir v > 0 */
            acb_pow(t, w1, S->p + j, prec);
            acb_mul(coef, coef, t, prec);
            acb_set(lam, S->p + j);
            continue;
        }
        acb_poly_set_coeff_acb(f, 0, w0);
        acb_poly_set_coeff_acb(f, 1, w1);
        acb_poly_pow_acb_series(g, f, S->p + j, N, prec);
        acb_poly_mullow(h, h, g, N, prec);
    }
    /* Sum_k h_k delta^(lam + k + 1)/(lam + k + 1) */
    acb_set_arb(d, delta);
    acb_add_si(t, lam, 1, prec);
    acb_pow(w0, d, t, prec);                 /* delta^(lam + 1) */
    acb_zero(sum);
    acb_one(w1);                             /* delta^k */
    mag_t last;
    mag_init(last);
    for (slong k = 0; k < N; k++) {
        acb_poly_get_coeff_acb(t, h, k);
        acb_mul(t, t, w1, prec);
        acb_add_si(k1, lam, k + 1, prec);
        acb_div(t, t, k1, prec);
        acb_add(sum, sum, t, prec);
        if (k == N - 1) acb_get_mag(last, t);
        acb_mul(w1, w1, d, prec);
    }
    /* tail: the terms fall by 4x a term (delta is a quarter of the radius) */
    mag_mul_ui(last, last, 2);
    acb_add_error_mag(sum, last);
    acb_mul(sum, sum, w0, prec);
    acb_mul(res, sum, coef, prec);
    for (int j = 0; j < AP_NFAC; j++)
        if (S->flip[j]) {
            acb_exp_pi_i(t, S->p + j, prec);       /* e^(i pi p) */
            acb_mul(res, res, t, prec);
        }
    mag_clear(last);
    acb_clear(coef);
    acb_clear(w0); acb_clear(w1); acb_clear(lam); acb_clear(d); acb_clear(t); acb_clear(sum); acb_clear(k1);
    acb_poly_clear(h); acb_poly_clear(f); acb_poly_clear(g);
}

static int ap_on_cut(const acb_t z) {
    if (!arb_is_zero(acb_imagref(z))) return 0;
    arb_t one; arb_init(one); arb_one(one);
    int r = arb_gt(acb_realref(z), one);
    arb_clear(one);
    return r;
}

static void ap_euler_acb(acb_t res, acb_srcptr p, slong prec) {
    double re[6], im[6];
    for (int k = 0; k < 6; k++) {
        re[k] = arf_get_d(arb_midref(acb_realref(p + k)), ARF_RND_NEAR);
        im[k] = arf_get_d(arb_midref(acb_imagref(p + k)), ARF_RND_NEAR);
    }
    acb_t ca;
    acb_init(ca);
    acb_sub(ca, p + 3, p + 0, prec);
    int ok = arb_is_positive(acb_realref(p + 0)) && arb_is_positive(acb_realref(ca));
    int cut[2] = { 0, 0 };
    for (int k = 0; k < 2 && ok; k++) {
        cut[k] = ap_on_cut(p + 4 + k);
        arb_t one; arb_init(one); arb_one(one);
        if (!cut[k] && arb_contains_zero(acb_imagref(p + 4 + k)) &&
            !arb_lt(acb_realref(p + 4 + k), one)) ok = 0;
        arb_clear(one);
    }
    if (cut[0] && cut[1] && re[4] == re[5]) ok = 0;
    if (!ok) { acb_indeterminate(res); acb_clear(ca); return; }
    /* singular points on the path (doubles for the geometry), with the factor
     * vanishing at each: 0 -> 0, 1 -> 1, 1/x -> 2, 1/y -> 3 */
    double pt[4]; int fac[4], np = 0;
    pt[np] = 0.0; fac[np++] = 0;
    for (int k = 0; k < 2; k++) if (cut[k]) { pt[np] = 1.0 / re[4 + k]; fac[np++] = 2 + k; }
    pt[np] = 1.0; fac[np++] = 1;
    if (np == 4 && pt[1] > pt[2]) { double t = pt[1]; pt[1] = pt[2]; pt[2] = t; int f = fac[1]; fac[1] = fac[2]; fac[2] = f; }
    /* zeros of all factors in the t-plane, for the radii */
    double zr[4], zi[4]; int nz = 0;
    zr[nz] = 0.0; zi[nz++] = 0.0;
    zr[nz] = 1.0; zi[nz++] = 0.0;
    for (int k = 0; k < 2; k++) {
        double xr = re[4 + k], xi = im[4 + k], d = xr * xr + xi * xi;
        if (d > 0.0) { zr[nz] = xr / d; zi[nz++] = -xi / d; }
    }
    slong N = prec / 2 + 20;
    ap_euler_seg S;
    ap_seg_init(&S);
    acb_t part, e, a1, b1;
    acb_init(part); acb_init(e); acb_init(a1); acb_init(b1);
    arb_t del;
    arb_init(del);
    mag_t tol;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);
    acb_zero(res);
    int flips[2] = { 0, 0 };
    for (int s = 0; s + 1 < np && ok; s++) {
        ap_seg_set(&S, p, flips, prec);
        /* exact endpoints */
        acb_t E[2];
        acb_init(E[0]); acb_init(E[1]);
        for (int side = 0; side < 2; side++) {
            int f = fac[s + side];
            if (f == 0) acb_zero(E[side]);
            else if (f == 1) acb_one(E[side]);
            else acb_inv(E[side], p + 4 + (f - 2), prec);
        }
        double len = pt[s + 1] - pt[s];
        double dl[2];
        for (int side = 0; side < 2; side++) {
            double R = HUGE_VAL, x0 = pt[s + side];
            for (int z = 0; z < nz; z++) {
                double dd = hypot(zr[z] - x0, zi[z]);
                if (dd > 1e-300 && dd < R) R = dd;
            }
            dl[side] = fmin(R / 4.0, len / 3.0);
        }
        /* [e0, e0 + d0] and [e1 - d1, e1] by series; the middle by quadrature */
        arb_set_d(del, dl[0]);
        ap_euler_end(part, &S, fac[s], E[0], +1, del, N, prec);
        acb_add(res, res, part, prec);
        arb_set_d(del, dl[1]);
        ap_euler_end(part, &S, fac[s + 1], E[1], -1, del, N, prec);
        acb_add(res, res, part, prec);
        acb_set_d(a1, dl[0]); acb_add(a1, a1, E[0], prec);
        acb_set_d(b1, -dl[1]); acb_add(b1, b1, E[1], prec);
        acb_calc_integrate(part, ap_euler_f, &S, a1, b1, prec, tol, NULL, prec);
        for (int j = 0; j < AP_NFAC; j++)
            if (S.flip[j]) { acb_exp_pi_i(e, S.p + j, prec); acb_mul(part, part, e, prec); }
        acb_add(res, res, part, prec);
        acb_clear(E[0]); acb_clear(E[1]);
        /* crossing a cut breakpoint flips its factor */
        int f = fac[s + 1];
        if (f >= 2) flips[f - 2] = 1;
    }
    acb_t g;
    acb_init(g);
    acb_gamma(g, p + 3, prec); acb_mul(res, res, g, prec);
    acb_rgamma(g, p + 0, prec); acb_mul(res, res, g, prec);
    acb_rgamma(g, ca, prec); acb_mul(res, res, g, prec);
    acb_clear(g);
    mag_clear(tol);
    arb_clear(del);
    acb_clear(part); acb_clear(e); acb_clear(a1); acb_clear(b1);
    ap_seg_clear(&S);
    acb_clear(ca);
}

/* the reductions to one Gauss function: 1 answered, 0 none */
static int ap_special_acb(acb_t res, int kind, acb_srcptr p, slong prec) {
    int np = appell_np(kind);
    acb_srcptr x = p + np, y = p + np + 1;
    int sx = 0, sy = 0;   /* on a cut: acb's 2F1 is the value from below */
    (void)sx; (void)sy;
    acb_t s, one;
    acb_init(s); acb_init(one);
    acb_one(one);
    int done = 1;
    switch (kind) {
        case 1:
            if (acb_is_zero(p + 0)) acb_one(res);
            else if (acb_is_zero(x) || acb_is_zero(p + 1)) acb_hypgeom_2f1(res, p + 0, p + 2, p + 3, y, 0, prec);
            else if (acb_is_zero(y) || acb_is_zero(p + 2)) acb_hypgeom_2f1(res, p + 0, p + 1, p + 3, x, 0, prec);
            else if (acb_is_exact(x) && acb_equal(x, y)) {
                acb_add(s, p + 1, p + 2, prec);
                acb_hypgeom_2f1(res, p + 0, s, p + 3, x, 0, prec);
            } else if (acb_is_exact(x) && acb_is_one(x)) {
                acb_hypgeom_2f1(res, p + 0, p + 1, p + 3, one, 0, prec);
                acb_sub(s, p + 3, p + 1, prec);
                acb_t v; acb_init(v);
                acb_hypgeom_2f1(v, p + 0, p + 2, s, y, 0, prec);
                acb_mul(res, res, v, prec);
                acb_clear(v);
            } else if (acb_is_exact(y) && acb_is_one(y)) {
                acb_hypgeom_2f1(res, p + 0, p + 2, p + 3, one, 0, prec);
                acb_sub(s, p + 3, p + 2, prec);
                acb_t v; acb_init(v);
                acb_hypgeom_2f1(v, p + 0, p + 1, s, x, 0, prec);
                acb_mul(res, res, v, prec);
                acb_clear(v);
            } else done = 0;
            break;
        case 2:
            if (acb_is_zero(p + 0) || (acb_is_zero(x) && acb_is_zero(y))) acb_one(res);
            else if (acb_is_zero(x) || acb_is_zero(p + 1)) acb_hypgeom_2f1(res, p + 0, p + 2, p + 4, y, 0, prec);
            else if (acb_is_zero(y) || acb_is_zero(p + 2)) acb_hypgeom_2f1(res, p + 0, p + 1, p + 3, x, 0, prec);
            else done = 0;
            break;
        case 3:
            if (acb_is_zero(x) || acb_is_zero(p + 0) || acb_is_zero(p + 2))
                acb_hypgeom_2f1(res, p + 1, p + 3, p + 4, y, 0, prec);
            else if (acb_is_zero(y) || acb_is_zero(p + 1) || acb_is_zero(p + 3))
                acb_hypgeom_2f1(res, p + 0, p + 2, p + 4, x, 0, prec);
            else done = 0;
            break;
        default:
            if (acb_is_zero(p + 0) || acb_is_zero(p + 1)) acb_one(res);
            else if (acb_is_zero(x)) acb_hypgeom_2f1(res, p + 0, p + 1, p + 3, y, 0, prec);
            else if (acb_is_zero(y)) acb_hypgeom_2f1(res, p + 0, p + 1, p + 2, x, 0, prec);
            else done = 0;
            break;
    }
    acb_clear(s); acb_clear(one);
    return done;
}

void appell_acb(acb_t res, int kind, acb_srcptr p, slong prec) {
    int n = APPELL_NARGS(kind), np = appell_np(kind);
    double re[7], im[7];
    for (int k = 0; k < n; k++) {
        if (!acb_is_finite(p + k)) { acb_indeterminate(res); return; }
        re[k] = arf_get_d(arb_midref(acb_realref(p + k)), ARF_RND_NEAR);
        im[k] = arf_get_d(arb_midref(acb_imagref(p + k)), ARF_RND_NEAR);
    }
    /* a pole: a denominator parameter a non-positive integer */
    if ((kind == 1 && ap_nonpos_int_acb(p + 3)) || (kind == 3 && ap_nonpos_int_acb(p + 4)) ||
        (kind == 2 && (ap_nonpos_int_acb(p + 3) || ap_nonpos_int_acb(p + 4))) ||
        (kind == 4 && (ap_nonpos_int_acb(p + 2) || ap_nonpos_int_acb(p + 3)))) {
        acb_indeterminate(res); return;
    }
    if (ap_special_acb(res, kind, p, prec)) return;
    AppellPlan P;
    if (!appell_plan(kind, re, im, 0.995, &P)) {
        if (kind == 1) ap_euler_acb(res, p, prec);
        else acb_indeterminate(res);
        return;
    }
    const AppellTr* t = appell_tr(kind, P.tr);
    acb_srcptr x = p + np, y = p + np + 1;
    acb_ptr q = _acb_vec_init(AP_MAXP);
    acb_ptr e = _acb_vec_init(6);
    acb_t X, Y, pre, base, ex, w, S;
    acb_init(X); acb_init(Y); acb_init(pre); acb_init(base); acb_init(ex); acb_init(w); acb_init(S);
    slong wp = prec + 10;
    for (int j = 0; j < np; j++) {
        acb_zero(q + j);
        for (int k = 0; k < np; k++)
            if (t->newp[j][k]) { acb_mul_si(w, p + k, t->newp[j][k], wp); acb_add(q + j, q + j, w, wp); }
    }
    ap_var_acb(X, t->vx, x, y, wp);
    ap_var_acb(Y, t->vy, x, y, wp);
    if (P.swap) acb_swap(X, Y);
    int idx[6];
    appell_engine_map(kind, P.swap, idx);
    for (int k = 0; k < 6; k++) acb_set(e + k, q + idx[k]);
    int ok = ap_eval_acb(S, kind, e, X, Y, P.sigY, wp);
    if (ok) {
        acb_one(pre);
        for (int b = 0; b < t->nbase; b++) {
            acb_one(base);
            if (t->base[b] == AP_B_1X) acb_sub(base, base, x, wp);
            else if (t->base[b] == AP_B_1Y) acb_sub(base, base, y, wp);
            else { acb_sub(base, base, x, wp); acb_sub(base, base, y, wp); }
            acb_zero(ex);
            for (int k = 0; k < np; k++)
                if (t->expo[b][k]) { acb_mul_si(w, p + k, t->expo[b][k], wp); acb_add(ex, ex, w, wp); }
            acb_pow(w, base, ex, wp);
            acb_mul(pre, pre, w, wp);
        }
        acb_mul(res, pre, S, prec);
    } else if (kind == 1) {
        ap_euler_acb(res, p, prec);
    } else {
        acb_indeterminate(res);
    }
    acb_clear(X); acb_clear(Y); acb_clear(pre); acb_clear(base); acb_clear(ex); acb_clear(w); acb_clear(S);
    _acb_vec_clear(q, AP_MAXP);
    _acb_vec_clear(e, 6);
}

#else
typedef int appell_arb_unused;
#endif
