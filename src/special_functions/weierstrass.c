/* Mathilda -- the Weierstrass elliptic functions (DLMF 23):
 *
 *   WeierstrassP[u, {g2, g3}]        wp(u),  wp'^2 = 4 wp^3 - g2 wp - g3
 *   WeierstrassPPrime[u, {g2, g3}]   wp'(u)
 *   WeierstrassZeta[u, {g2, g3}]     zeta(u), zeta' = -wp
 *   WeierstrassSigma[u, {g2, g3}]    sigma(u), sigma'/sigma = zeta
 *   WeierstrassHalfPeriods[{g2, g3}] {w1, w3}
 *   WeierstrassInvariants[{w1, w3}]  {g2, g3}
 *   InverseWeierstrassP[p, {g2, g3}] a u with wp(u) = p
 *
 * The half-periods follow the reference system's normalisation, recovered from
 * its outputs on 380 sampled invariants (all reproduced): the lattice basis with
 * tau = w3/w1 in the closed fundamental domain |Re tau| <= 1/2, |tau| >= 1, with
 * Re tau = +1/2 (not -1/2) and Re tau >= 0 on |tau| = 1 at the boundary, and the
 * sign Re w1 > 0 (Im w1 < 0 when w1 is imaginary). The square lattice (g3 = 0)
 * takes arg w1 in [-pi/4, pi/4), w3 = I w1; the hexagonal one (g2 = 0) arg w1 in
 * [-pi/6, pi/6) and arg w3 in [pi/3, 2 pi/3). Function values do not depend on
 * the choice.
 *
 * Machine kernel (scalar, NDArray / packed buffers, Compile[]): the lattice is
 * built once per (g2, g3) -- roots of 4t^3 - g2 t - g3 by Cardano and Newton, a
 * basis from w1 = K(m)/sqrt(e1 - e3), w3 = I K(1 - m)/sqrt(e1 - e3) (complex
 * AGM), the normalisation above, and a CERTIFICATE: the invariants recomputed
 * from the theta constants of the reduced tau must reproduce (g2, g3) to 1e-11,
 * or the next root ordering is tried (a sublattice fails it). After the
 * reduction |q| <= exp(-pi sqrt(3)/2) = 0.066, so every theta series has <= 6
 * terms. u is reduced modulo the lattice (z = u/(2 w1) = z0 + m + n tau) and
 *   wp    = (pi^2 th2^2 th3^2 th4(v)^2/th1(v)^2 - pi^2 (th2^4 + th3^4)/3) / (2w1)^2
 *   wp'   = -2 pi^3 th2(v) th3(v) th4(v) th1'(0)^2 / th1(v)^3 / (2w1)^3
 *   zeta  = (2 eta z + pi (th1'(v)/th1(v) - 2 i n)) / (2w1)
 *   sigma = 2w1 (-1)^(m+n) th1(v) / (pi th1'(0)) exp(eta z^2 - i pi tau n^2 - 2 i n v)
 * with v = pi z0, eta = pi^2 E2(tau)/6. A value the double route cannot certify
 * (cancellation beyond 2^4 near a zero, an argument whose reduction amplifies
 * rounding by more than 2^6, a far-out sigma) is answered by arb inside the
 * kernel (flint_weierstrass_d), so the buffer and Compile[] answers bit-match
 * the scalar one. A memo of eight plans serves compiled loops over u.
 *
 * Arbitrary precision: arb's acb_elliptic_p / p_jet / zeta / sigma on the
 * lattice (1, tau) built the same way in ball arithmetic (flint_num_weierstrass).
 *
 * Exact rules (the reference system's): {0, 0} gives 1/u^2, -2/u^3, 1/u, u; a
 * double root (g2^3 = 27 g3^2) the hyperbolic / trigonometric closed forms;
 * u = 0 a pole (sigma 0); parity (wp even, the others odd).
 *
 * Attributes: NumericFunction, Protected (the half-period and invariant heads
 * Protected only). The functions thread over a List first argument, as the
 * reference does.
 */
#include "weierstrass.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "numeric.h"
#include "message.h"
#include "flint_num_bridge.h"
#include "sf_scalar.h"
#include "sf_dd.h"
#include "ndarray.h"
#include "ndarray_internal.h"
#include "symtab.h"
#include "../sym_names.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Complex double helpers (no <complex.h>)                             */
/* ------------------------------------------------------------------ */

typedef struct { double r, i; } wcx;
static wcx W(double r, double i) { wcx z; z.r = r; z.i = i; return z; }
static wcx wadd(wcx a, wcx b) { return W(a.r + b.r, a.i + b.i); }
static wcx wsub(wcx a, wcx b) { return W(a.r - b.r, a.i - b.i); }
static wcx wmul(wcx a, wcx b) { return W(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
static wcx wsc(wcx a, double s) { return W(a.r * s, a.i * s); }
static wcx wneg(wcx a) { return W(-a.r, -a.i); }
static double wabs(wcx a) { return hypot(a.r, a.i); }
static wcx wdiv(wcx a, wcx b) {
    /* Smith's algorithm */
    if (fabs(b.r) >= fabs(b.i)) {
        double t = b.i / b.r, d = b.r + b.i * t;
        return W((a.r + a.i * t) / d, (a.i - a.r * t) / d);
    }
    double t = b.r / b.i, d = b.r * t + b.i;
    return W((a.r * t + a.i) / d, (a.i * t - a.r) / d);
}
static wcx wsqrt(wcx z) {
    if (z.r == 0.0 && z.i == 0.0) return W(0.0, 0.0);
    double m = wabs(z);
    double t = sqrt(0.5 * (m + fabs(z.r)));
    if (z.r >= 0.0) return W(t, z.i / (2.0 * t));
    return W(fabs(z.i) / (2.0 * t), z.i >= 0.0 ? t : -t);
}
static wcx wexp(wcx z) { double e = exp(z.r); return W(e * cos(z.i), e * sin(z.i)); }
static wcx wcbrt(wcx z) {
    double m = wabs(z);
    if (m == 0.0) return W(0.0, 0.0);
    double a = atan2(z.i, z.r) / 3.0, r = cbrt(m);
    return W(r * cos(a), r * sin(a));
}
static wcx wsin(wcx z) { return W(sin(z.r) * cosh(z.i), cos(z.r) * sinh(z.i)); }
static wcx wcos(wcx z) { return W(cos(z.r) * cosh(z.i), -sin(z.r) * sinh(z.i)); }
static wcx wsinh(wcx z) { return W(sinh(z.r) * cos(z.i), cosh(z.r) * sin(z.i)); }
static wcx wcosh(wcx z) { return W(cosh(z.r) * cos(z.i), sinh(z.r) * sin(z.i)); }
static bool wfin(wcx z) { return isfinite(z.r) && isfinite(z.i); }

/* ------------------------------------------------------------------ */
/* The lattice                                                         */
/* ------------------------------------------------------------------ */

/* Roots of 4 t^3 - g2 t - g3: Cardano on t^3 + p t + q (p = -g2/4,
 * q = -g3/4), each polished by two Newton steps. */
static void wp_roots(wcx g2, wcx g3, wcx r[3]) {
    wcx p = wsc(g2, -0.25), q = wsc(g3, -0.25);
    wcx hq = wsc(q, -0.5);                                  /* -q/2 */
    wcx p3 = wsc(p, 1.0 / 3.0);
    wcx D = wadd(wmul(hq, hq), wmul(p3, wmul(p3, p3)));
    wcx s = wsqrt(D);
    wcx A = wadd(hq, s), B = wsub(hq, s);
    if (wabs(B) > wabs(A)) A = B;                           /* no cancellation */
    wcx C = wcbrt(A);
    const wcx om[3] = { { 1.0, 0.0 }, { -0.5, 0.86602540378443864676 },
                        { -0.5, -0.86602540378443864676 } };
    for (int k = 0; k < 3; k++) {
        wcx Ck = wmul(C, om[k]);
        wcx t = (wabs(Ck) == 0.0) ? W(0.0, 0.0) : wsub(Ck, wdiv(p3, Ck));
        for (int it = 0; it < 2; it++) {
            wcx t2 = wmul(t, t);
            wcx f = wsub(wsub(wsc(wmul(t2, t), 4.0), wmul(g2, t)), g3);
            wcx fp = wsub(wsc(t2, 12.0), g2);
            if (wabs(fp) == 0.0) break;
            t = wsub(t, wdiv(f, fp));
        }
        r[k] = t;
    }
}

/* The arithmetic-geometric mean with the right choice of square root at each
 * step (|a - b| <= |a + b|): K(m) = pi / (2 AGM(1, sqrt(1 - m))), principal. */
static wcx wp_agm(wcx a, wcx b) {
    for (int it = 0; it < 64; it++) {
        if (wabs(wsub(a, b)) <= 1e-16 * wabs(a)) break;
        wcx an = wsc(wadd(a, b), 0.5);
        wcx bn = wsqrt(wmul(a, b));
        if (wabs(wsub(an, bn)) > wabs(wadd(an, bn))) bn = wneg(bn);
        a = an; b = bn;
    }
    return a;
}
static wcx wp_K(wcx m) {
    wcx M = wp_agm(W(1.0, 0.0), wsqrt(wsub(W(1.0, 0.0), m)));
    return wdiv(W(0.5 * M_PI, 0.0), M);
}

/* the basis operation (w1, w3) <- (a w1 + b w3, c w1 + d w3), accumulated */
static void wp_op(wcx* w1, wcx* w3, double M[4], double a, double b, double c, double d) {
    wcx n1 = wadd(wsc(*w1, a), wsc(*w3, b));
    wcx n3 = wadd(wsc(*w1, c), wsc(*w3, d));
    *w1 = n1; *w3 = n3;
    double m0 = a * M[0] + b * M[2], m1 = a * M[1] + b * M[3];
    double m2 = c * M[0] + d * M[2], m3 = c * M[1] + d * M[3];
    M[0] = m0; M[1] = m1; M[2] = m2; M[3] = m3;
}

#define WP_TOL 1e-12

bool sf_wp_normalize(double w1r, double w1i, double w3r, double w3i, bool sq, bool hex,
                     double M[4]) {
    wcx w1 = W(w1r, w1i), w3 = W(w3r, w3i);
    M[0] = 1.0; M[1] = 0.0; M[2] = 0.0; M[3] = 1.0;
    if (!wfin(w1) || !wfin(w3) || wabs(w1) == 0.0) return false;
    wcx t = wdiv(w3, w1);
    if (!(fabs(t.i) > 1e-300) || !wfin(t)) return false;
    if (t.i < 0.0) wp_op(&w1, &w3, M, 1, 0, 0, -1);
    int it;
    for (it = 0; it < 200; it++) {
        t = wdiv(w3, w1);
        double n = nearbyint(t.r);
        if (n != 0.0) { wp_op(&w1, &w3, M, 1, 0, -n, 1); t = wdiv(w3, w1); }
        if (t.r * t.r + t.i * t.i < 1.0 - 1e-13) wp_op(&w1, &w3, M, 0, 1, -1, 0);
        else break;
    }
    if (it == 200) return false;
    t = wdiv(w3, w1);
    if (sq) {
        /* tau = i: w1 -> i^k w1 (i w1 = w3, i w3 = -w1) into arg [-pi/4, pi/4) */
        for (int k = 0; k < 4; k++) {
            double a = atan2(w1.i, w1.r);
            if (a >= 0.25 * M_PI - WP_TOL)       wp_op(&w1, &w3, M, 0, -1, 1, 0);  /* x (-i) */
            else if (a < -0.25 * M_PI - WP_TOL)  wp_op(&w1, &w3, M, 0, 1, -1, 0);  /* x i */
            else break;
        }
        return true;
    }
    if (hex) {
        /* tau = rho = e^(i pi/3): rho w1 = w3, rho w3 = w3 - w1 */
        if (t.r < 0.0) wp_op(&w1, &w3, M, 1, 0, 1, 1);
        for (int k = 0; k < 6; k++) {
            double a = atan2(w1.i, w1.r);
            if (a >= M_PI / 6.0 - WP_TOL)        wp_op(&w1, &w3, M, 1, -1, 1, 0);  /* x rho^-1 */
            else if (a < -M_PI / 6.0 - WP_TOL)   wp_op(&w1, &w3, M, 0, 1, -1, 1);  /* x rho */
            else break;
        }
        if (atan2(w1.i, w1.r) < -WP_TOL) wp_op(&w1, &w3, M, 1, 0, -1, 1);      /* w3 = rho^2 w1 */
        return true;
    }
    if (t.r < -0.5 + WP_TOL) { wp_op(&w1, &w3, M, 1, 0, 1, 1); t = wdiv(w3, w1); }
    if (t.r * t.r + t.i * t.i < 1.0 + WP_TOL && t.r < -WP_TOL) wp_op(&w1, &w3, M, 0, 1, -1, 0);
    double aw = wabs(w1);
    if (w1.r < -WP_TOL * aw || (fabs(w1.r) <= WP_TOL * aw && w1.i > 0.0))
        wp_op(&w1, &w3, M, -1, 0, 0, -1);
    return true;
}

/* theta2(0), theta3(0), theta4(0) at the nome q (|q| small), q4 = q^(1/4) */
static void wp_theta0(wcx q, wcx q4, wcx* t2, wcx* t3, wcx* t4) {
    wcx s2 = W(1.0, 0.0), s3 = W(0.0, 0.0), s4 = W(0.0, 0.0);
    wcx q2 = wmul(q, q);
    wcx qa = W(1.0, 0.0), qstep = q2;         /* q^(n(n+1)): ratio q^(2n) */
    wcx qn = W(1.0, 0.0), qodd = q;           /* q^(n^2): ratio q^(2n-1) */
    for (int n = 1; n < 40; n++) {
        qa = wmul(qa, qstep); qstep = wmul(qstep, q2);
        qn = wmul(qn, qodd);  qodd = wmul(qodd, q2);
        s2 = wadd(s2, qa);
        s3 = wadd(s3, qn);
        s4 = (n & 1) ? wsub(s4, qn) : wadd(s4, qn);
        if (wabs(qn) < 1e-18 && wabs(qa) < 1e-18) break;
    }
    *t2 = wsc(wmul(q4, s2), 2.0);
    *t3 = wadd(W(1.0, 0.0), wsc(s3, 2.0));
    *t4 = wadd(W(1.0, 0.0), wsc(s4, 2.0));
}

/* E2(tau) = 1 - 24 Sum_n n q^(2n)/(1 - q^(2n)) */
static wcx wp_e2(wcx q) {
    wcx q2 = wmul(q, q), qk = W(1.0, 0.0), s = W(0.0, 0.0);
    for (int n = 1; n < 60; n++) {
        qk = wmul(qk, q2);
        wcx t = wsc(wdiv(qk, wsub(W(1.0, 0.0), qk)), (double)n);
        s = wadd(s, t);
        if (wabs(t) < 1e-18 * (1.0 + wabs(s))) break;
    }
    return wsub(W(1.0, 0.0), wsc(s, 24.0));
}

/* ---- double-double helpers ---- */

static scdd cd_sub(scdd a, scdd b) { return scdd_add(a, scdd_mul_d(b, -1.0)); }
static scdd cd_div(scdd a, scdd b) { return scdd_mul(a, scdd_inv(b)); }
static double cd_abs(scdd a) { return hypot(a.re.hi, a.im.hi); }
static scdd cd_from(wcx z) { return scdd_from(z.r, z.i); }
static wcx cd_d(scdd a) { return W(a.re.hi + a.re.lo, a.im.hi + a.im.lo); }
static sdd dd_pi(void) { return sdd_quick(SDD_PI_HI, SDD_PI_LO); }

/* the function constants of the plan from theta2, 3, 4 at 0 and 1/(2 w1),
 * formed in double-double and rounded once */
static void wp_consts(WpPlan* P, scdd t2, scdd t3, scdd t4, scdd i2) {
    sdd pi = dd_pi();
    scdd pi2; pi2.re = sdd_mul(pi, pi); pi2.im = sdd_from(0.0);
    scdd i22 = scdd_mul(i2, i2);
    scdd t22 = scdd_mul(t2, t2), t32 = scdd_mul(t3, t3);
    scdd A = scdd_mul(scdd_mul(pi2, scdd_mul(t22, t32)), i22);
    scdd B = scdd_add(scdd_mul(t22, t22), scdd_mul(t32, t32));
    B = scdd_mul(scdd_mul(pi2, B), i22);
    B.re = sdd_div_d(B.re, 3.0); B.im = sdd_div_d(B.im, 3.0);
    scdd t1p = scdd_mul(scdd_mul(t2, t3), t4);
    scdd pi3; pi3.re = sdd_mul(pi2.re, pi); pi3.im = sdd_from(0.0);
    scdd K = scdd_mul_d(scdd_mul(scdd_mul(pi3, scdd_mul(t1p, t1p)), scdd_mul(i22, i2)), -2.0);
    scdd pt; pt.re = sdd_mul(t1p.re, pi); pt.im = sdd_mul(t1p.im, pi);
    scdd S = scdd_inv(scdd_mul(pt, i2));
    wcx a = cd_d(A), b = cd_d(B), k = cd_d(K), sg = cd_d(S);
    P->pAr = a.r; P->pAi = a.i; P->pBr = b.r; P->pBi = b.i;
    P->kpr = k.r; P->kpi = k.i; P->scr = sg.r; P->sci = sg.i;
}

static bool wp_fill(WpPlan* P, wcx w1, wcx w3) {
    wcx tau = wdiv(w3, w1);
    if (!(tau.i > 0.0) || tau.i > 100.0) return false;
    wcx ipt = W(-M_PI * tau.i, M_PI * tau.r);                /* i pi tau */
    wcx q = wexp(ipt), q4 = wexp(wsc(ipt, 0.25));
    wcx t2, t3, t4;
    wp_theta0(q, q4, &t2, &t3, &t4);
    wcx a2 = wmul(wmul(t2, t2), wmul(t2, t2));
    wcx a3 = wmul(wmul(t3, t3), wmul(t3, t3));
    wcx a4 = wmul(wmul(t4, t4), wmul(t4, t4));
    /* the certificate: invariants of the reduced lattice */
    wcx i2 = wdiv(W(1.0, 0.0), wsc(w1, 2.0));
    wcx c = wmul(wsc(i2, M_PI), wsc(i2, M_PI));              /* (pi/(2 w1))^2 */
    wcx s8 = wadd(wadd(wmul(a2, a2), wmul(a3, a3)), wmul(a4, a4));
    wcx G2 = wsc(wmul(wmul(c, c), s8), 2.0 / 3.0);
    wcx G3 = wsc(wmul(wmul(wmul(c, c), c),
                      wmul(wmul(wadd(a3, a4), wsub(a2, a4)), wadd(a2, a3))), -4.0 / 27.0);
    wcx g2 = W(P->g2r, P->g2i), g3 = W(P->g3r, P->g3i);
    double sc = fmax(sqrt(wabs(g2)), cbrt(wabs(g3)));
    if (!(sc > 0.0) || !wfin(G2) || !wfin(G3)) return false;
    if (wabs(wsub(G2, g2)) > 1e-11 * sc * sc || wabs(wsub(G3, g3)) > 1e-11 * sc * sc * sc)
        return false;
    P->w1r = w1.r; P->w1i = w1.i; P->w3r = w3.r; P->w3i = w3.i;
    P->tr = tau.r; P->ti = tau.i;
    P->qr = q.r; P->qi = q.i; P->q4r = q4.r; P->q4i = q4.i;
    P->i2r = i2.r; P->i2i = i2.i;
    wcx et = wsc(wp_e2(q), M_PI * M_PI / 6.0);
    wp_consts(P, scdd_from(t2.r, t2.i), scdd_from(t3.r, t3.i), scdd_from(t4.r, t4.i),
              scdd_from(i2.r, i2.i));
    P->etr = et.r; P->eti = et.i;
    P->trl = P->til = P->i2rl = P->i2il = P->etrl = P->etil = 0.0;
    return true;
}


static sdd dd_sqrt(sdd a) {
    if (!(a.hi > 0.0)) return sdd_from(0.0);
    double s = sqrt(a.hi);
    sdd e = sdd_sub(a, sdd_two_prod(s, s));
    return sdd_add(sdd_from(s), sdd_div_d(e, 2.0 * s));
}
static scdd cd_sqrt(scdd z) {
    if (z.re.hi == 0.0 && z.im.hi == 0.0) return scdd_from(0.0, 0.0);
    sdd m = dd_sqrt(sdd_add(sdd_mul(z.re, z.re), sdd_mul(z.im, z.im)));
    sdd ax = z.re.hi >= 0.0 ? z.re : sdd_neg(z.re);
    sdd t = dd_sqrt(sdd_mul_d(sdd_add(m, ax), 0.5));
    scdd r;
    if (z.re.hi >= 0.0) { r.re = t; r.im = sdd_div(z.im, sdd_mul_d(t, 2.0)); }
    else {
        sdd ay = z.im.hi >= 0.0 ? z.im : sdd_neg(z.im);
        r.re = sdd_div(ay, sdd_mul_d(t, 2.0));
        r.im = z.im.hi >= 0.0 ? t : sdd_neg(t);
    }
    return r;
}
static scdd cd_agm(scdd a, scdd b) {
    for (int it = 0; it < 80; it++) {
        if (cd_abs(cd_sub(a, b)) <= 1e-31 * cd_abs(a)) break;
        scdd an = scdd_mul_d(scdd_add(a, b), 0.5);
        scdd bn = cd_sqrt(scdd_mul(a, b));
        if (cd_abs(cd_sub(an, bn)) > cd_abs(scdd_add(an, bn))) bn = scdd_mul_d(bn, -1.0);
        a = an; b = bn;
    }
    return a;
}
static scdd cd_K(scdd m) {
    scdd one = scdd_from(1.0, 0.0);
    scdd M = cd_agm(one, cd_sqrt(cd_sub(one, m)));
    scdd hp; hp.re = sdd_quick(SDD_PIO2_HI, SDD_PIO2_LO); hp.im = sdd_from(0.0);
    return cd_div(hp, M);
}
/* theta constants and the zeta constant from the (accurate) nome, in
 * double-double: theta2, 3, 4 at 0 by their q-series */
static void wp_constants(WpPlan* P, wcx q, wcx q4, scdd qd, scdd q4d, scdd i2) {
    scdd one = scdd_from(1.0, 0.0);
    scdd s2 = one, s3 = scdd_from(0.0, 0.0), s4 = scdd_from(0.0, 0.0);
    scdd q2 = scdd_mul(qd, qd), qa = one, qstep = q2, qn = one, qodd = qd;
    for (int n = 1; n < 40; n++) {
        qa = scdd_mul(qa, qstep); qstep = scdd_mul(qstep, q2);
        qn = scdd_mul(qn, qodd);  qodd = scdd_mul(qodd, q2);
        s2 = scdd_add(s2, qa);
        s3 = scdd_add(s3, qn);
        s4 = (n & 1) ? cd_sub(s4, qn) : scdd_add(s4, qn);
        if (cd_abs(qn) < 1e-34 && cd_abs(qa) < 1e-34) break;
    }
    scdd t2 = scdd_mul_d(scdd_mul(q4d, s2), 2.0);
    scdd t3 = scdd_add(one, scdd_mul_d(s3, 2.0));
    scdd t4 = scdd_add(one, scdd_mul_d(s4, 2.0));
    P->qr = q.r; P->qi = q.i; P->q4r = q4.r; P->q4i = q4.i;
    wp_consts(P, t2, t3, t4, i2);
    /* eta = pi^2 E2 / 6, E2 = 1 - 24 Sum sigma1(k) q^(2k) (the Lambert series
     * Sum n x^n/(1 - x^n) regrouped: no divisions), double-double */
    static const double SIG1[49] = { 0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24,
        31, 18, 39, 20, 42, 32, 36, 24, 60, 31, 42, 40, 56, 30, 72, 32, 63, 48, 54, 48, 91,
        38, 60, 56, 90, 42, 96, 44, 84, 78, 72, 48, 124 };
    scdd qk = scdd_from(1.0, 0.0), s = scdd_from(0.0, 0.0);
    for (int k = 1; k < 49; k++) {
        qk = scdd_mul(qk, q2);
        scdd t = scdd_mul_d(qk, SIG1[k]);
        s = scdd_add(s, t);
        if (cd_abs(t) < 1e-33 * (1.0 + cd_abs(s))) break;
    }
    scdd e2 = cd_sub(scdd_from(1.0, 0.0), scdd_mul_d(s, 24.0));
    sdd pi = dd_pi();
    sdd c = sdd_div_d(sdd_mul(pi, pi), 6.0);
    scdd et; et.re = sdd_mul(e2.re, c); et.im = sdd_mul(e2.im, c);
    P->etr = et.re.hi; P->etrl = et.re.lo; P->eti = et.im.hi; P->etil = et.im.lo;
}

static bool wp_refine(WpPlan* P, const wcx r[3], const int perm[3], const double M[4]) {
    scdd g2 = scdd_from(P->g2r, P->g2i), g3 = scdd_from(P->g3r, P->g3i);
    scdd e[3];
    for (int k = 0; k < 3; k++) {
        scdd t = cd_from(r[perm[k]]);
        for (int it = 0; it < 2; it++) {
            scdd t2 = scdd_mul(t, t);
            scdd f = cd_sub(cd_sub(scdd_mul_d(scdd_mul(t2, t), 4.0), scdd_mul(g2, t)), g3);
            wcx fp = wsub(wsc(wmul(cd_d(t), cd_d(t)), 12.0), W(P->g2r, P->g2i));
            if (wabs(fp) == 0.0) return false;
            t = cd_sub(t, cd_div(f, cd_from(fp)));
        }
        e[k] = t;
    }
    scdd d13 = cd_sub(e[0], e[2]);
    if (cd_abs(d13) == 0.0) return false;
    scdd s = cd_sqrt(d13), m = cd_div(cd_sub(e[1], e[2]), d13);
    scdd w1 = cd_div(cd_K(m), s);
    scdd k2 = cd_div(cd_K(cd_sub(scdd_from(1.0, 0.0), m)), s);
    scdd w3; w3.re = sdd_neg(k2.im); w3.im = k2.re;                /* i K(1-m)/s */
    scdd n1 = scdd_add(scdd_mul_d(w1, M[0]), scdd_mul_d(w3, M[1]));
    scdd n3 = scdd_add(scdd_mul_d(w1, M[2]), scdd_mul_d(w3, M[3]));
    if (P->real) {
        double a1 = cd_abs(n1), a3 = cd_abs(n3);
        if (fabs(n1.re.hi) < 1e-14 * a1) n1.re = sdd_from(0.0);
        if (fabs(n1.im.hi) < 1e-14 * a1) n1.im = sdd_from(0.0);
        if (fabs(n3.re.hi) < 1e-14 * a3) n3.re = sdd_from(0.0);
        if (fabs(n3.im.hi) < 1e-14 * a3) n3.im = sdd_from(0.0);
    }
    /* must agree with the certified double basis */
    wcx d1 = cd_d(n1), d3 = cd_d(n3);
    if (wabs(wsub(d1, W(P->w1r, P->w1i))) > 1e-12 * wabs(d1) ||
        wabs(wsub(d3, W(P->w3r, P->w3i))) > 1e-12 * wabs(d3))
        return false;
    scdd tau = cd_div(n3, n1);
    scdd i2 = scdd_inv(scdd_mul_d(n1, 2.0));
    sdd pi = dd_pi();
    scdd ipt; ipt.re = sdd_neg(sdd_mul(pi, tau.im)); ipt.im = sdd_mul(pi, tau.re);
    scdd q4d = scdd_exp(scdd_mul_d(ipt, 0.25));
    scdd qd = scdd_mul(q4d, q4d); qd = scdd_mul(qd, qd);              /* q = (q^(1/4))^4 */
    P->w1r = d1.r; P->w1i = d1.i; P->w3r = d3.r; P->w3i = d3.i;
    P->tr = tau.re.hi; P->trl = tau.re.lo; P->ti = tau.im.hi; P->til = tau.im.lo;
    P->i2r = i2.re.hi; P->i2rl = i2.re.lo; P->i2i = i2.im.hi; P->i2il = i2.im.lo;
    wp_constants(P, cd_d(qd), cd_d(q4d), qd, q4d, i2);
    return true;
}

bool sf_wp_plan(double g2r, double g2i, double g3r, double g3i, WpPlan* P) {
    memset(P, 0, sizeof *P);
    if (!isfinite(g2r) || !isfinite(g2i) || !isfinite(g3r) || !isfinite(g3i)) return false;
    P->g2r = g2r; P->g2i = g2i; P->g3r = g3r; P->g3i = g3i;
    P->real = (g2i == 0.0 && g3i == 0.0);
    wcx g2 = W(g2r, g2i), g3 = W(g3r, g3i);
    if (wabs(g2) == 0.0 && wabs(g3) == 0.0) { P->kind = 1; return true; }
    double s2 = wabs(g2), s3 = wabs(g3);
    if (s2 > 1e100 || s3 > 1e150 || (s2 < 1e-100 && s3 < 1e-150)) return false;
    wcx D = wsub(wmul(g2, wmul(g2, g2)), wsc(wmul(g3, g3), 27.0));
    if (D.r == 0.0 && D.i == 0.0) {
        /* a double root e = -3 g3 / (2 g2), c = sqrt(3 e) */
        wcx e = wdiv(wsc(g3, -3.0), wsc(g2, 2.0));
        wcx c = wsqrt(wsc(e, 3.0));
        P->kind = 2; P->er = e.r; P->ei = e.i; P->cr = c.r; P->ci = c.i;
        return true;
    }
    wcx r[3];
    wp_roots(g2, g3, r);
    static const int PERM[6][3] = { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} };
    for (int p = 0; p < 6; p++) {
        wcx e1 = r[PERM[p][0]], e2 = r[PERM[p][1]], e3 = r[PERM[p][2]];
        wcx d13 = wsub(e1, e3);
        if (wabs(d13) == 0.0) continue;
        wcx s = wsqrt(d13), m = wdiv(wsub(e2, e3), d13);
        wcx w1 = wdiv(wp_K(m), s);
        wcx w3 = wmul(W(0.0, 1.0), wdiv(wp_K(wsub(W(1.0, 0.0), m)), s));
        double M[4];
        if (!sf_wp_normalize(w1.r, w1.i, w3.r, w3.i, g3r == 0.0 && g3i == 0.0,
                             g2r == 0.0 && g2i == 0.0, M))
            continue;
        wcx n1 = wadd(wsc(w1, M[0]), wsc(w3, M[1]));
        wcx n3 = wadd(wsc(w1, M[2]), wsc(w3, M[3]));
        if (P->real) {
            /* a real lattice: clear the rounding residue of an exactly zero part */
            double a1 = wabs(n1), a3 = wabs(n3);
            if (fabs(n1.r) < 1e-14 * a1) n1.r = 0.0;
            if (fabs(n1.i) < 1e-14 * a1) n1.i = 0.0;
            if (fabs(n3.r) < 1e-14 * a3) n3.r = 0.0;
            if (fabs(n3.i) < 1e-14 * a3) n3.i = 0.0;
        }
        if (wp_fill(P, n1, n3)) {
            WpPlan keep = *P;
            if (!wp_refine(P, r, PERM[p], M)) *P = keep;
            return true;
        }
    }
    return false;
}

/* A small append-only-then-ring memo of plans, copied out under the lock, so a
 * compiled loop over u builds the lattice once. */
#define WP_MEMO 8
static WpPlan g_wp_memo[WP_MEMO];
static int g_wp_memo_n = 0, g_wp_memo_next = 0;
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_wp_lock = PTHREAD_MUTEX_INITIALIZER;
#define WP_LOCK()   pthread_mutex_lock(&g_wp_lock)
#define WP_UNLOCK() pthread_mutex_unlock(&g_wp_lock)
#else
#define WP_LOCK()   ((void)0)
#define WP_UNLOCK() ((void)0)
#endif

static bool wp_plan_get(double g2r, double g2i, double g3r, double g3i, WpPlan* out) {
    bool hit = false;
    WP_LOCK();
    for (int s = 0; s < g_wp_memo_n; s++) {
        const WpPlan* m = &g_wp_memo[s];
        if (m->g2r == g2r && m->g2i == g2i && m->g3r == g3r && m->g3i == g3i) {
            *out = *m; hit = true; break;
        }
    }
    WP_UNLOCK();
    if (hit) return true;
    if (!sf_wp_plan(g2r, g2i, g3r, g3i, out)) return false;
    WP_LOCK();
    g_wp_memo[g_wp_memo_next] = *out;
    g_wp_memo_next = (g_wp_memo_next + 1) % WP_MEMO;
    if (g_wp_memo_n < WP_MEMO) g_wp_memo_n++;
    WP_UNLOCK();
    return true;
}

/* ------------------------------------------------------------------ */
/* Function values                                                     */
/* ------------------------------------------------------------------ */

typedef struct { wcx t1, t1p, t2, t3, t4; double l1, l1p, l2, l3, l4; } WpTh;

/* the largest cancellation (peak term / |sum|) a machine value may carry */
#define WP_LOSS 4.0

/* theta_1..4(v), theta_1'(v) at nome q, v = vh + vl with |vl| ~ 1e-16 |vh|:
 * the series run at vh and the first-order term vl theta'(vh) restores the
 * rounding of v to a double (so a value near a zero, where the function is
 * ill-conditioned in v, is still the value at the exact argument). l* = the
 * peak term / |sum| of each series. */
static bool wp_thetas(wcx v, wcx vl, wcx q, wcx q4, WpTh* T) {
    wcx s1 = wsin(v), c1 = wcos(v);
    if (!wfin(s1) || !wfin(c1)) return false;
    wcx sk = s1, ck = c1;                      /* sin(k v), cos(k v), k = 1 */
    wcx q2 = wmul(q, q);
    wcx qa = W(1.0, 0.0), qstep = q2;          /* q^(n(n+1)) */
    wcx qn = W(1.0, 0.0), qodd = q;            /* q^(n^2) */
    wcx A1 = W(0, 0), A1p = W(0, 0), A2 = W(0, 0), A3 = W(0.5, 0), A4 = W(0.5, 0);
    wcx B1pp = W(0, 0), B2p = W(0, 0), B3p = W(0, 0), B4p = W(0, 0);   /* derivatives */
    double p1 = 0, p1p = 0, p2 = 0, p3 = 0.5, p4 = 0.5;
    for (int n = 0; n < 12; n++) {
        /* odd multiple k = 2n + 1: sk, ck hold sin/cos(k v) */
        double k = 2.0 * n + 1.0, sg = (n & 1) ? -1.0 : 1.0;
        wcx a = wmul(qa, sk), b = wmul(qa, ck);
        wcx t1 = wsc(a, sg), t1p = wsc(b, sg * k);
        A1 = wadd(A1, t1); A1p = wadd(A1p, t1p); A2 = wadd(A2, b);
        B1pp = wadd(B1pp, wsc(t1, -k * k)); B2p = wadd(B2p, wsc(a, -k));
        double m1 = wabs(t1), m1p = wabs(t1p), m2 = wabs(b);
        if (m1 > p1) p1 = m1;
        if (m1p > p1p) p1p = m1p;
        if (m2 > p2) p2 = m2;
        /* even multiple 2(n + 1) */
        wcx s2 = wadd(wmul(sk, c1), wmul(ck, s1)), c2 = wsub(wmul(ck, c1), wmul(sk, s1));
        qn = wmul(qn, qodd); qodd = wmul(qodd, q2);
        wcx e = wmul(qn, c2);
        A3 = wadd(A3, e);
        A4 = (n & 1) ? wadd(A4, e) : wsub(A4, e);
        wcx ep = wsc(wmul(qn, s2), -2.0 * (n + 1));
        B3p = wadd(B3p, ep);
        B4p = (n & 1) ? wadd(B4p, ep) : wsub(B4p, ep);
        double me = wabs(e);
        if (me > p3) p3 = me;
        if (me > p4) p4 = me;
        /* advance to the next odd multiple 2n + 3 */
        sk = wadd(wmul(s2, c1), wmul(c2, s1));
        ck = wsub(wmul(c2, c1), wmul(s2, s1));
        qa = wmul(qa, qstep); qstep = wmul(qstep, q2);
        if (!wfin(sk) || !wfin(ck) || !wfin(qa)) return false;
        /* the next terms against the largest so far (a sum that cancels below
         * its peak is caught by the loss ratio, not by more terms) */
        double tail = wabs(qa) * fmax(wabs(sk), wabs(ck)) * (k + 2.0);
        double tail2 = wabs(wmul(qn, qodd)) * 2.0 * fmax(wabs(s2), wabs(c2)) * (k + 3.0);
        if (n >= 1 && tail < 1e-18 * fmin(p1, fmin(p1p, p2)) && tail2 < 1e-18 * p3)
            break;
        if (n == 11) return false;
    }
    wcx f = wsc(q4, 2.0);
    T->t1 = wmul(f, A1); T->t1p = wmul(f, A1p); T->t2 = wmul(f, A2);
    T->t3 = wsc(A3, 2.0); T->t4 = wsc(A4, 2.0);
    if (vl.r != 0.0 || vl.i != 0.0) {
        T->t1 = wadd(T->t1, wmul(vl, T->t1p));
        T->t1p = wadd(T->t1p, wmul(vl, wmul(f, B1pp)));
        T->t2 = wadd(T->t2, wmul(vl, wmul(f, B2p)));
        T->t3 = wadd(T->t3, wmul(vl, wsc(B3p, 2.0)));
        T->t4 = wadd(T->t4, wmul(vl, wsc(B4p, 2.0)));
    }
    T->l1 = p1 / wabs(A1); T->l1p = p1p / wabs(A1p); T->l2 = p2 / wabs(A2);
    T->l3 = p3 / wabs(A3); T->l4 = p4 / wabs(A4);
    return wfin(T->t1) && wfin(T->t2) && wfin(T->t3) && wfin(T->t4) && wfin(T->t1p);
}

static bool wp_out(wcx v, double* rr, double* ri) {
    if (!wfin(v)) return false;
    double big = fmax(fabs(v.r), fabs(v.i));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(v.r) < DBL_MIN) v.r = 0.0;
    if (fabs(v.i) < DBL_MIN) v.i = 0.0;
    *rr = v.r; *ri = v.i;
    return true;
}

static bool wp_eval_degenerate(const WpPlan* P, int kind, wcx u, wcx* out) {
    if (P->kind == 1) {
        wcx iu = wdiv(W(1.0, 0.0), u);
        switch (kind) {
            case WP_P:      *out = wmul(iu, iu); break;
            case WP_PPRIME: *out = wsc(wmul(iu, wmul(iu, iu)), -2.0); break;
            case WP_ZETA:   *out = iu; break;
            default:        *out = u; break;
        }
        return true;
    }
    wcx e = W(P->er, P->ei), c = W(P->cr, P->ci);
    wcx cu = wmul(c, u), sh = wsinh(cu), ch = wcosh(cu);
    if (!wfin(sh) || !wfin(ch)) return false;
    switch (kind) {
        case WP_P:      /* e + c^2 / sinh^2(c u) */
            *out = wadd(e, wdiv(wmul(c, c), wmul(sh, sh))); break;
        case WP_PPRIME: /* -2 c^3 cosh / sinh^3 */
            *out = wsc(wdiv(wmul(wmul(c, wmul(c, c)), ch), wmul(sh, wmul(sh, sh))), -2.0); break;
        case WP_ZETA:   /* -e u + c coth(c u) */
            *out = wadd(wneg(wmul(e, u)), wdiv(wmul(c, ch), sh)); break;
        default:        /* exp(-e u^2/2) sinh(c u)/c */
            *out = wmul(wexp(wsc(wmul(e, wmul(u, u)), -0.5)), wdiv(sh, c)); break;
    }
    return true;
}

bool sf_wp_eval_plan(const WpPlan* P, int kind, double ur, double ui, double* rr, double* ri) {
    wcx u = W(ur, ui), out;
    if (!wfin(u)) return false;
    if (P->kind != 0) {
        if (ur == 0.0 && ui == 0.0) {
            if (kind != WP_SIGMA) return false;
            *rr = 0.0; *ri = 0.0; return true;
        }
        if (!wp_eval_degenerate(P, kind, u, &out)) return false;
        if (P->real && ui == 0.0) out.i = 0.0;
        return wp_out(out, rr, ri);
    }
    /* z = u/(2 w1) and its reduction z0 = z - n tau - m in double-double, so
     * z0 is accurate to ~1e-32 |z| whatever the number of periods crossed */
    scdd i2d, taud;
    i2d.re = sdd_quick(P->i2r, P->i2rl); i2d.im = sdd_quick(P->i2i, P->i2il);
    taud.re = sdd_quick(P->tr, P->trl); taud.im = sdd_quick(P->ti, P->til);
    scdd zd = scdd_mul(scdd_from(ur, ui), i2d);
    double n = nearbyint(zd.im.hi / P->ti);
    if (fabs(n) > 1e6) return false;
    scdd z1 = cd_sub(zd, scdd_mul_d(taud, n));
    double m = nearbyint(z1.re.hi);
    if (fabs(m) > 1e6) return false;
    scdd z0d = z1; z0d.re = sdd_add_d(z1.re, -m);
    wcx z = cd_d(zd), z0 = cd_d(z0d);
    double az0 = wabs(z0);
    if (az0 == 0.0 || az0 < 1e-15 * (wabs(z) + 1.0)) return false;   /* at a lattice point */
    sdd pid = dd_pi();
    sdd vre = sdd_mul(z0d.re, pid), vim = sdd_mul(z0d.im, pid);
    wcx v = W(vre.hi, vim.hi), vl = W(vre.lo, vim.lo);
    WpTh T;
    if (!wp_thetas(v, vl, W(P->qr, P->qi), W(P->q4r, P->q4i), &T)) return false;
    {
        double l = fmax(T.l1, kind == WP_P ? T.l4 : kind == WP_PPRIME ? fmax(T.l2, fmax(T.l3, T.l4))
                                        : kind == WP_ZETA ? T.l1p : 0.0);
        if (!(l <= WP_LOSS)) return false;
    }
    switch (kind) {
        case WP_P: {
            wcx r = wdiv(T.t4, T.t1);
            wcx a = wmul(W(P->pAr, P->pAi), wmul(r, r));
            wcx b = W(P->pBr, P->pBi);
            wcx d = wsub(a, b);
            if (wabs(d) * WP_LOSS < fmax(wabs(a), wabs(b))) return false;
            out = d;
            break;
        }
        case WP_PPRIME: {
            wcx num = wmul(wmul(T.t2, T.t3), T.t4);
            wcx den = wmul(T.t1, wmul(T.t1, T.t1));
            out = wmul(W(P->kpr, P->kpi), wdiv(num, den));
            break;
        }
        case WP_ZETA: {
            /* 2 eta z - 2 pi i n in double-double (the quasi-period parts cancel
             * against each other far out), + pi th1'/th1(v) */
            scdd etd; etd.re = sdd_quick(P->etr, P->etrl); etd.im = sdd_quick(P->eti, P->etil);
            scdd a = scdd_mul_d(scdd_mul(etd, zd), 2.0);
            a.im = sdd_sub(a.im, sdd_mul_d(pid, 2.0 * n));
            wcx b = wsc(wdiv(T.t1p, T.t1), M_PI);
            scdd d = scdd_add(a, scdd_from(b.r, b.i));
            if (cd_abs(d) * WP_LOSS < fmax(cd_abs(a), wabs(b))) return false;
            out = cd_d(scdd_mul(d, i2d));
            break;
        }
        default: {
            /* exponent eta z^2 - i pi tau n^2 - 2 i pi n z0, double-double */
            scdd etd; etd.re = sdd_quick(P->etr, P->etrl); etd.im = sdd_quick(P->eti, P->etil);
            scdd E = scdd_mul(etd, scdd_mul(zd, zd));
            scdd t = scdd_mul_d(taud, n * n);                         /* tau n^2 */
            scdd pt; pt.re = sdd_mul(pid, t.re); pt.im = sdd_mul(pid, t.im);
            E.re = sdd_add(E.re, pt.im);                              /* - i (pi tau n^2) */
            E.im = sdd_sub(E.im, pt.re);
            scdd pz; pz.re = sdd_mul(pid, z0d.re); pz.im = sdd_mul(pid, z0d.im);
            E.re = sdd_add(E.re, sdd_mul_d(pz.im, 2.0 * n));          /* - 2 i n (pi z0) */
            E.im = sdd_sub(E.im, sdd_mul_d(pz.re, 2.0 * n));
            if (E.re.hi > 700.0 || E.re.hi < -700.0) return false;
            wcx s = wmul(W(P->scr, P->sci), T.t1);
            if (fmod(fabs(m + n), 2.0) == 1.0) s = wneg(s);
            /* exp(hi + lo) = exp(hi) (1 + lo), |lo| ~ 1e-16 |E| */
            wcx ex = wmul(wexp(W(E.re.hi, E.im.hi)), W(1.0 + E.re.lo, E.im.lo));
            out = wmul(s, ex);
            break;
        }
    }
    if (P->real && ui == 0.0) out.i = 0.0;
    return wp_out(out, rr, ri);
}

/* The second chance for a value the double route declined (a cancellation
 * near a zero): arb on the plan's own double-double z = u/(2 w1) and tau, which
 * costs one function call at ~128 bits instead of rebuilding the lattice. */
static bool wp_eval_plan_arb(const WpPlan* P, int kind, double ur, double ui,
                             double* rr, double* ri) {
    if (P->kind != 0) return false;
    scdd i2d;
    i2d.re = sdd_quick(P->i2r, P->i2rl); i2d.im = sdd_quick(P->i2i, P->i2il);
    scdd zd = scdd_mul(scdd_from(ur, ui), i2d);
    double z[4] = { zd.re.hi, zd.re.lo, zd.im.hi, zd.im.lo };
    double t[4] = { P->tr, P->trl, P->ti, P->til };
    double vr, vi;
    if (!flint_wp_unit_d(kind, z, t, &vr, &vi)) return false;
    /* scale by (2 w1)^-2, ^-3, ^-1 or (2 w1) */
    scdd s = scdd_from(1.0, 0.0);
    if (kind == WP_P) s = scdd_mul(i2d, i2d);
    else if (kind == WP_PPRIME) s = scdd_mul(i2d, scdd_mul(i2d, i2d));
    else if (kind == WP_ZETA) s = i2d;
    else s = scdd_inv(i2d);
    wcx out = wmul(W(vr, vi), cd_d(s));
    if (P->real && ui == 0.0) out.i = 0.0;
    return wp_out(out, rr, ri);
}

bool sf_machine_weierstrass(int kind, double g2r, double g2i, double g3r, double g3i,
                            double ur, double ui, double* rr, double* ri) {
    if (kind < WP_P || kind > WP_SIGMA) return false;
    if (!isfinite(ur) || !isfinite(ui)) return false;
    WpPlan P;
    if (wp_plan_get(g2r, g2i, g3r, g3i, &P)) {
        if (sf_wp_eval_plan(&P, kind, ur, ui, rr, ri)) return true;
        if (!(ur == 0.0 && ui == 0.0) && wp_eval_plan_arb(&P, kind, ur, ui, rr, ri)) return true;
    }
    if (ur == 0.0 && ui == 0.0) {
        if (kind != WP_SIGMA) return false;
        *rr = 0.0; *ri = 0.0;
        return true;
    }
    double vr, vi;
    if (!flint_weierstrass_d(kind, ur, ui, g2r, g2i, g3r, g3i, &vr, &vi)) return false;
    if (g2i == 0.0 && g3i == 0.0 && ui == 0.0) vi = 0.0;
    return wp_out(W(vr, vi), rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtins                                                            */
/* ------------------------------------------------------------------ */

static const char* const WP_NAME[4] = {
    "WeierstrassP", "WeierstrassPPrime", "WeierstrassZeta", "WeierstrassSigma" };

static bool wp_is_list2(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List &&
           e->data.function.arg_count == 2;
}
static bool wp_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* A machine scalar (Real / exact / Complex of those): its parts and whether it
 * was inexact. */
static bool wp_leaf(const Expr* e, double* re, double* im, bool* inexact, bool* cplx) {
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

static bool wp_is_zero(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    return false;
}

/* visibly negated: a negative number, or Times[c, ...] with c < 0 rational */
static bool wp_negated(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Rational && e->data.function.arg_count == 2)
        return wp_negated(e->data.function.args[0]);
    if (h == SYM_Times && e->data.function.arg_count >= 2) {
        const Expr* c = e->data.function.args[0];
        return (c->type == EXPR_INTEGER && c->data.integer < 0) ||
               (c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
                c->data.function.head->data.symbol.name == SYM_Rational && wp_negated(c));
    }
    return false;
}

static bool wp_exact_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Rational;
}

/* apply a pure-function template (a string) to the arguments, evaluated */
static Expr* wp_apply(const char* tmpl, Expr* const* args, size_t n) {
    Expr* f = parse_expression(tmpl);
    if (!f) return NULL;
    Expr* a[4];
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
    return eval_and_free(expr_new_function(f, a, n));
}

/* the exact double-root closed forms (e = -3 g3/(2 g2) real), #1 = u, #2 = e,
 * #3 = c: c = Sqrt[3 e] (e > 0, hyperbolic) or Sqrt[-3 e] (e < 0, trigonometric) */
static const char* const WP_DEG_H[4] = {
    "(#2 + #3^2 Csch[#3 #1]^2)&",
    "(-2 #3^3 Coth[#3 #1] Csch[#3 #1]^2)&",
    "(-#2 #1 + #3 Coth[#3 #1])&",
    "(E^(-#2 #1^2/2) #1 Sinc[I #3 #1])&" };
static const char* const WP_DEG_T[4] = {
    "(-2 #2 + #3^2 Cot[#3 #1]^2)&",
    "(-2 #3^3 Cot[#3 #1] Csc[#3 #1]^2)&",
    "(-#2 #1 + #3 Cot[#3 #1])&",
    "(E^(-#2 #1^2/2) #1 Sinc[#3 #1])&" };

/* The buffer path: u an NDArray / packed List, (g2, g3) machine or exact
 * numbers; the plan is built once and the elements evaluated in parallel. */
typedef struct { const WpPlan* P; int kind; const void* in; NDType dt; double* out;
                 int any_c; } WpNdCtx;

static bool wp_nd_chunk(void* vctx, size_t lo, size_t hi) {
    WpNdCtx* c = (WpNdCtx*)vctx;
    for (size_t i = lo; i < hi; i++) {
        double xr, xi, rr, ri;
        ndt_get(c->in, i, c->dt, &xr, &xi);
        if (c->P->kind < 0 || (!sf_wp_eval_plan(c->P, c->kind, xr, xi, &rr, &ri) &&
                               !wp_eval_plan_arb(c->P, c->kind, xr, xi, &rr, &ri))) {
            if (!sf_machine_weierstrass(c->kind, c->P->g2r, c->P->g2i, c->P->g3r, c->P->g3i,
                                        xr, xi, &rr, &ri))
                return false;
        }
        c->out[2 * i] = rr; c->out[2 * i + 1] = ri;
        if (ri != 0.0) c->any_c = 1;   /* benign race: only ever set to 1 */
    }
    return true;
}

static Expr* wp_ndarray(int kind, const Expr* u, const double g[4]) {
    NDType dt = u->data.ndarray.dtype;
    if (dt == NDT_BOOL) return NULL;
    WpPlan P;
    if (!sf_wp_plan(g[0], g[1], g[2], g[3], &P)) {
        /* every element through the kernel (its arb fallback) */
        memset(&P, 0, sizeof P);
        P.kind = -1;
        P.g2r = g[0]; P.g2i = g[1]; P.g3r = g[2]; P.g3i = g[3];
    }
    size_t sz = ndarray_size(u);
    double* buf = malloc(sizeof(double) * 2 * (sz ? sz : 1));
    if (!buf) return NULL;
    WpNdCtx ctx = { &P, kind, u->data.ndarray.data, dt, buf, 0 };
    bool ok = nd_parallel_for(sz, wp_nd_chunk, &ctx);
    if (!ok) { free(buf); return NULL; }
    bool cx = ctx.any_c || ndt_is_complex(dt) || g[1] != 0.0 || g[3] != 0.0;
    NDType odt = cx ? ((dt == NDT_FLOAT32 || dt == NDT_COMPLEX32) ? NDT_COMPLEX32 : NDT_COMPLEX64)
                    : ((dt == NDT_FLOAT32) ? NDT_FLOAT32 : NDT_FLOAT64);
    void* out = malloc(ndt_elem_size(odt) * (sz ? sz : 1));
    if (!out) { free(buf); return NULL; }
    for (size_t i = 0; i < sz; i++) ndt_set(out, i, odt, buf[2 * i], buf[2 * i + 1]);
    free(buf);
    return expr_new_ndarray_like(u, u->data.ndarray.rank, u->data.ndarray.dims, out, odt);
}

static long wp_bits(const Expr* const* a, int n) {
    long bits = 0;
    for (int i = 0; i < n; i++) {
        long b = numeric_min_inexact_bits(a[i]);
        if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    return bits;
}

static Expr* wp_builtin(Expr* res, int kind) {
    const char* name = WP_NAME[kind];
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error(name, argc, 2, 2);
    Expr* u = res->data.function.args[0];
    Expr* G = res->data.function.args[1];
    if (!wp_is_list2(G)) return NULL;
    Expr* g2 = G->data.function.args[0];
    Expr* g3 = G->data.function.args[1];

    /* thread over a List first argument (the reference does) */
    if (wp_is_list(u)) {
        size_t n = u->data.function.arg_count;
        Expr** el = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!el) return NULL;
        for (size_t i = 0; i < n; i++)
            el[i] = expr_new_function(expr_new_symbol(name),
                        (Expr*[]){ expr_copy(u->data.function.args[i]), expr_copy(G) }, 2);
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), el, n);
        free(el);
        return out;
    }

    double gv[4];
    bool ix2, ix3, c2, c3;
    bool gm = wp_leaf(g2, &gv[0], &gv[1], &ix2, &c2) && wp_leaf(g3, &gv[2], &gv[3], &ix3, &c3);

    if (is_ndarray(u)) {
        Expr* out = gm ? wp_ndarray(kind, u, gv) : NULL;
        return out ? out : ndarray_delist_and_reeval(res);
    }

    /* the pole at u = 0 (sigma vanishes) */
    if (wp_is_zero(u)) {
        if (kind == WP_SIGMA) return u->type == EXPR_REAL ? expr_new_real(0.0) : expr_new_integer(0);
        mth_message("Power", "infy", "Infinite expression %s encountered.",
                    u->type == EXPR_REAL ? "1/0." : "1/0");
        return expr_new_symbol(SYM_ComplexInfinity);
    }

    double ur, ui;
    bool ixu, cu;
    if (gm && wp_leaf(u, &ur, &ui, &ixu, &cu) && (ixu || ix2 || ix3)) {
        double rr, ri;
        if (sf_machine_weierstrass(kind, gv[0], gv[1], gv[2], gv[3], ur, ui, &rr, &ri)) {
            bool real_out = !cu && !c2 && !c3 && ri == 0.0;
            return real_out ? expr_new_real(rr) : sf_machine_complex_result(rr, ri);
        }
    }

    /* g2 = g3 = 0: the degenerate lattice */
    if (g2->type == EXPR_INTEGER && g2->data.integer == 0 &&
        g3->type == EXPR_INTEGER && g3->data.integer == 0) {
        static const char* const Z[4] = { "(#1^-2)&", "(-2 #1^-3)&", "(1/#1)&", "#1&" };
        return wp_apply(Z[kind], &u, 1);
    }

    const Expr* nargs[3] = { u, g2, g3 };
    bool numeric = expr_is_numeric_quantity(u) && expr_is_numeric_quantity(g2) &&
                   expr_is_numeric_quantity(g3);
    if (numeric && wp_bits(nargs, 3) > 0) {
        Expr* r = flint_num_weierstrass(kind, u, g2, g3);
        if (r && wp_bits(nargs, 3) <= 53) r = numeric_lower_mpfr_to_machine(r);
        return r;
    }

    /* exact real invariants with a double root */
    if (wp_exact_number(g2) && wp_exact_number(g3) && !wp_is_zero(g2)) {
        Expr* D = wp_apply("(#1^3 - 27 #2^2)&", (Expr*[]){ g2, g3 }, 2);
        bool deg = D && D->type == EXPR_INTEGER && D->data.integer == 0;
        expr_free(D);
        if (deg) {
            Expr* e = wp_apply("(-3 #2/(2 #1))&", (Expr*[]){ g2, g3 }, 2);
            bool pos = e && !wp_negated(e);
            Expr* c = wp_apply(pos ? "Sqrt[3 #1]&" : "Sqrt[-3 #1]&", &e, 1);
            Expr* out = wp_apply(pos ? WP_DEG_H[kind] : WP_DEG_T[kind], (Expr*[]){ u, e, c }, 3);
            expr_free(e); expr_free(c);
            return out;
        }
    }

    if (wp_negated(u)) {
        Expr* nu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                 (Expr*[]){ expr_new_integer(-1), expr_copy(u) }, 2));
        Expr* call = expr_new_function(expr_new_symbol(name), (Expr*[]){ nu, expr_copy(G) }, 2);
        if (kind == WP_P) return eval_and_free(call);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_integer(-1), call }, 2));
    }
    return NULL;
}

static Expr* builtin_weierstrass_p(Expr* res)      { return wp_builtin(res, WP_P); }
static Expr* builtin_weierstrass_pprime(Expr* res) { return wp_builtin(res, WP_PPRIME); }
static Expr* builtin_weierstrass_zeta(Expr* res)   { return wp_builtin(res, WP_ZETA); }
static Expr* builtin_weierstrass_sigma(Expr* res)  { return wp_builtin(res, WP_SIGMA); }

/* WeierstrassHalfPeriods[{g2, g3}] */
static Expr* builtin_weierstrass_halfperiods(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("WeierstrassHalfPeriods", argc, 1, 1);
    Expr* G = res->data.function.args[0];
    if (!wp_is_list2(G)) return NULL;
    Expr* g2 = G->data.function.args[0];
    Expr* g3 = G->data.function.args[1];
    if (wp_is_zero(g2) && wp_is_zero(g3))
        return parse_expression("{Infinity, DirectedInfinity[I]}");
    double gv[4];
    bool ix2, ix3, c2, c3;
    if (wp_leaf(g2, &gv[0], &gv[1], &ix2, &c2) && wp_leaf(g3, &gv[2], &gv[3], &ix3, &c3) &&
        (ix2 || ix3)) {
        WpPlan P;
        if (sf_wp_plan(gv[0], gv[1], gv[2], gv[3], &P)) {
            if (P.kind == 2) {
                /* a double root: e > 0 {-I pi/(2 c), Infinity}, e < 0 {pi/(2 c), I Infinity} */
                if (P.ei == 0.0 && P.er != 0.0) {
                    double c = sqrt(3.0 * fabs(P.er));
                    Expr* a = P.er > 0.0 ? sf_machine_complex_result(0.0, -M_PI / (2.0 * c))
                                         : expr_new_real(M_PI / (2.0 * c));
                    Expr* b = parse_expression(P.er > 0.0 ? "Infinity" : "DirectedInfinity[I]");
                    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ a, b }, 2);
                }
                return NULL;
            }
            if (P.kind == 0) {
                Expr* a = (P.w1i == 0.0) ? expr_new_real(P.w1r) : sf_machine_complex_result(P.w1r, P.w1i);
                Expr* b = (P.w3i == 0.0) ? expr_new_real(P.w3r) : sf_machine_complex_result(P.w3r, P.w3i);
                return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ a, b }, 2);
            }
        }
    }
    /* exact invariants with a double root */
    if (wp_exact_number(g2) && wp_exact_number(g3) && !wp_is_zero(g2)) {
        Expr* D = wp_apply("(#1^3 - 27 #2^2)&", (Expr*[]){ g2, g3 }, 2);
        bool deg = D && D->type == EXPR_INTEGER && D->data.integer == 0;
        expr_free(D);
        if (deg) {
            Expr* e = wp_apply("(-3 #2/(2 #1))&", (Expr*[]){ g2, g3 }, 2);
            Expr* out = wp_apply(wp_negated(e) ? "{Pi/(2 Sqrt[-3 #1]), DirectedInfinity[I]}&"
                                               : "{-I Pi/(2 Sqrt[3 #1]), Infinity}&", &e, 1);
            expr_free(e);
            return out;
        }
    }
    const Expr* a2[2] = { g2, g3 };
    if (expr_is_numeric_quantity(g2) && expr_is_numeric_quantity(g3) && wp_bits(a2, 2) > 0) {
        Expr* r = flint_num_wp_halfperiods(g2, g3);
        if (r && wp_bits(a2, 2) <= 53) r = numeric_lower_mpfr_to_machine(r);
        return r;
    }
    return NULL;
}

/* {g2, g3} of the reduced basis (w1, w3) from the theta constants; NULL when
 * tau is too close to the real axis for the double route. A part below the
 * rounding level of the value (a rectangular or rhombic lattice) is cleared. */
static Expr* wp_invariants_machine(wcx w1, wcx w3) {
    wcx tau = wdiv(w3, w1);
    if (!wfin(tau) || !(tau.i > 0.0) || tau.i > 100.0) return NULL;
    wcx ipt = W(-M_PI * tau.i, M_PI * tau.r);
    wcx q = wexp(ipt), q4 = wexp(wsc(ipt, 0.25)), t2, t3, t4;
    wp_theta0(q, q4, &t2, &t3, &t4);
    wcx A2 = wmul(wmul(t2, t2), wmul(t2, t2)), A3 = wmul(wmul(t3, t3), wmul(t3, t3));
    wcx A4 = wmul(wmul(t4, t4), wmul(t4, t4));
    wcx i2 = wdiv(W(1.0, 0.0), wsc(w1, 2.0));
    wcx c = wmul(wsc(i2, M_PI), wsc(i2, M_PI));
    wcx s8 = wadd(wadd(wmul(A2, A2), wmul(A3, A3)), wmul(A4, A4));
    wcx G2 = wsc(wmul(wmul(c, c), s8), 2.0 / 3.0);
    wcx G3 = wsc(wmul(wmul(wmul(c, c), c),
                      wmul(wmul(wadd(A3, A4), wsub(A2, A4)), wadd(A2, A3))), -4.0 / 27.0);
    if (!wfin(G2) || !wfin(G3)) return NULL;
    double s = fmax(wabs(G2), pow(wabs(G3), 2.0 / 3.0));
    if (fabs(G2.i) < 1e-14 * s) G2.i = 0.0;
    if (fabs(G3.i) < 1e-14 * pow(s, 1.5)) G3.i = 0.0;
    Expr* x = G2.i == 0.0 ? expr_new_real(G2.r) : sf_machine_complex_result(G2.r, G2.i);
    Expr* y = G3.i == 0.0 ? expr_new_real(G3.r) : sf_machine_complex_result(G3.r, G3.i);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ x, y }, 2);
}

/* WeierstrassInvariants[{w1, w3}] */
static Expr* builtin_weierstrass_invariants(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("WeierstrassInvariants", argc, 1, 1);
    Expr* Wl = res->data.function.args[0];
    if (!wp_is_list2(Wl)) return NULL;
    Expr* w1 = Wl->data.function.args[0];
    Expr* w3 = Wl->data.function.args[1];
    double a[4];
    bool i1, i3, c1, c3;
    if (wp_leaf(w1, &a[0], &a[1], &i1, &c1) && wp_leaf(w3, &a[2], &a[3], &i3, &c3) && (i1 || i3)) {
        double M[4];
        if (!sf_wp_normalize(a[0], a[1], a[2], a[3], false, false, M)) {
            char* s = expr_to_string(Wl);
            mth_message("EllipticReducedHalfPeriods", "nind", "Half-periods %s are not independent.", s);
            free(s);
            return NULL;
        }
        wcx o1 = W(a[0], a[1]), o3 = W(a[2], a[3]);
        wcx n1 = wadd(wsc(o1, M[0]), wsc(o3, M[1])), n3 = wadd(wsc(o1, M[2]), wsc(o3, M[3]));
        Expr* r = wp_invariants_machine(n1, n3);
        if (r) return r;
    }
    const Expr* a2[2] = { w1, w3 };
    if (expr_is_numeric_quantity(w1) && expr_is_numeric_quantity(w3)) {
        if (wp_bits(a2, 2) > 0) {
            Expr* r = flint_num_wp_invariants(w1, w3);
            if (r && wp_bits(a2, 2) <= 53) r = numeric_lower_mpfr_to_machine(r);
            return r;
        }
        /* exact: the square lattice w3 = I w1 and the hexagonal w3 = E^(I Pi/3) w1 */
        Expr* t = wp_apply("Together[#2/#1]&", (Expr*[]){ w1, w3 }, 2);
        Expr* out = NULL;
        if (t) {
            Expr* sq = wp_apply("(#1 === I || #1 === -I)&", &t, 1);
            Expr* hx = wp_apply("MemberQ[{(-1)^(1/3), (-1)^(2/3), -(-1)^(1/3), -(-1)^(2/3)}, #1]&", &t, 1);
            if (sq && sq->type == EXPR_SYMBOL && sq->data.symbol.name == SYM_True)
                out = wp_apply("{Gamma[1/4]^8/(256 Pi^2 #1^4), 0}&", &w1, 1);
            else if (hx && hx->type == EXPR_SYMBOL && hx->data.symbol.name == SYM_True)
                out = wp_apply("{0, Gamma[1/3]^18/(4096 Pi^6 #1^6)}&", &w1, 1);
            expr_free(sq); expr_free(hx); expr_free(t);
        }
        return out;
    }
    return NULL;
}

void weierstrass_init(void) {
    symtab_add_builtin("WeierstrassP", builtin_weierstrass_p);
    symtab_add_builtin("WeierstrassPPrime", builtin_weierstrass_pprime);
    symtab_add_builtin("WeierstrassZeta", builtin_weierstrass_zeta);
    symtab_add_builtin("WeierstrassSigma", builtin_weierstrass_sigma);
    symtab_add_builtin("WeierstrassHalfPeriods", builtin_weierstrass_halfperiods);
    symtab_add_builtin("WeierstrassInvariants", builtin_weierstrass_invariants);
    for (int i = 0; i < 4; i++)
        symtab_get_def(WP_NAME[i])->attributes |= (ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_get_def("WeierstrassHalfPeriods")->attributes |= ATTR_PROTECTED;
    symtab_get_def("WeierstrassInvariants")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("WeierstrassP",
        "WeierstrassP[u, {g2, g3}]\n"
        "\tgives the Weierstrass elliptic function wp(u; g2, g3), the solution of\n"
        "\twp'^2 = 4 wp^3 - g2 wp - g3 with a double pole at u = 0.\n"
        "Numeric at machine and arbitrary precision, real and complex; threads over a\n"
        "List or packed array u and compiles. D[wp, u] = WeierstrassPPrime.");
    symtab_set_docstring("WeierstrassPPrime",
        "WeierstrassPPrime[u, {g2, g3}]\n"
        "\tgives the derivative of WeierstrassP[u, {g2, g3}] with respect to u.\n"
        "Numeric at machine and arbitrary precision; threads over packed u and compiles.");
    symtab_set_docstring("WeierstrassZeta",
        "WeierstrassZeta[u, {g2, g3}]\n"
        "\tgives the Weierstrass zeta function, zeta' = -WeierstrassP, zeta ~ 1/u at 0.\n"
        "Numeric at machine and arbitrary precision; threads over packed u and compiles.");
    symtab_set_docstring("WeierstrassSigma",
        "WeierstrassSigma[u, {g2, g3}]\n"
        "\tgives the Weierstrass sigma function, sigma'/sigma = WeierstrassZeta, an entire\n"
        "\tfunction with sigma ~ u at 0.\n"
        "Numeric at machine and arbitrary precision; threads over packed u and compiles.");
    symtab_set_docstring("WeierstrassHalfPeriods",
        "WeierstrassHalfPeriods[{g2, g3}]\n"
        "\tgives the half-periods {w1, w3} of the lattice with invariants g2, g3:\n"
        "\ttau = w3/w1 in the fundamental domain (Re tau = 1/2 rather than -1/2,\n"
        "\tRe tau >= 0 on |tau| = 1) and Re w1 > 0.");
    symtab_set_docstring("WeierstrassInvariants",
        "WeierstrassInvariants[{w1, w3}]\n"
        "\tgives the invariants {g2, g3} of the lattice with half-periods w1, w3.");
}
