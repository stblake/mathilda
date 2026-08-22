/*
 * solveint_conic.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "checked_int.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "linalg/hnf.h"
#include "solvethue.h"
#include "solveint_internal.h"


/* --- Conic  Y^2 == A X^2 + B X + C  with A a perfect square. ---
 *
 * When one variable Y appears only as Y^2 (coefficient +/-1) and the other, X,
 * to degree <= 2 with no cross term, and the X^2 coefficient A is a POSITIVE
 * PERFECT SQUARE p^2, the curve is a "parabola-like" hyperbola with finitely
 * many integer points.  Completing the square,
 *   (2 p Y)^2 - (2 A X + B)^2 = 4 A C - B^2 =: D,
 * so every solution comes from a factorisation  d1 * d2 = D  via
 *   Y = (d1 + d2) / (4 p),  X = (d2 - d1 - 2 B) / (4 A).
 * Enumerating the divisors of |D| is exhaustive, so an empty result is a proof
 * of no solutions.  Euler's  n^2 + n + 41 == y^2  has D = 163 (prime) and the
 * unique positive point n = 40, y = 41.  (A non-square A is a genuine Pell
 * conic -- left to the Pell continued-fraction path.) */
bool si_solve_conic(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return false;

    for (int choice = 0; choice < 2; choice++) {
        int Yv = act[choice], Xv = act[1 - choice];
        mpz_t e0, e2, e1, ec; mpz_init_set_ui(e0, 0); mpz_init_set_ui(e2, 0);
        mpz_init_set_ui(e1, 0); mpz_init_set_ui(ec, 0);
        bool shape = true;
        for (size_t t = 0; t < eq->n_terms && shape; t++) {
            const int* ex = eq->exps + t * (size_t)n;
            for (int v = 0; v < n; v++) if (v != Yv && v != Xv && ex[v] != 0) shape = false;
            if (!shape) break;
            int ey = ex[Yv], exx = ex[Xv];
            if (ey == 2 && exx == 0) mpz_add(e0, e0, eq->coefs[t]);
            else if (ey == 0 && exx == 2) mpz_add(e2, e2, eq->coefs[t]);
            else if (ey == 0 && exx == 1) mpz_add(e1, e1, eq->coefs[t]);
            else if (ey == 0 && exx == 0) mpz_add(ec, ec, eq->coefs[t]);
            else shape = false;                       /* Y linear, cross, or higher */
        }
        bool ok = shape && (mpz_cmpabs_ui(e0, 1) == 0);   /* Y^2 coefficient +/-1 */
        if (ok) {
            if (mpz_sgn(e0) < 0) { mpz_neg(e0, e0); mpz_neg(e2, e2); mpz_neg(e1, e1); mpz_neg(ec, ec); }
            /* Y^2 = A X^2 + B X + C. */
            mpz_t A, B, C, D, p, tmp; mpz_init(A); mpz_init(B); mpz_init(C);
            mpz_init(D); mpz_init(p); mpz_init(tmp);
            mpz_neg(A, e2); mpz_neg(B, e1); mpz_neg(C, ec);
            bool go = (mpz_sgn(A) > 0) && mpz_perfect_square_p(A);
            if (go) {
                mpz_sqrt(p, A);
                mpz_mul(D, A, C); mpz_mul_ui(D, D, 4);
                mpz_mul(tmp, B, B); mpz_sub(D, D, tmp);   /* D = 4AC - B^2 */
            }
            if (go && mpz_sgn(D) != 0) {
                st->max_visits = SI_MAX_NODES;
                mpz_t absD; mpz_init(absD); mpz_abs(absD, D);
                Expr* dl = divisors_ordinary(absD);
                mpz_clear(absD);
                mpz_t d1, d2, Y, X, fourp, fourA, twoB, num;
                mpz_init(d1); mpz_init(d2); mpz_init(Y); mpz_init(X);
                mpz_init(fourp); mpz_init(fourA); mpz_init(twoB); mpz_init(num);
                mpz_mul_ui(fourp, p, 4); mpz_mul_ui(fourA, A, 4); mpz_mul_ui(twoB, B, 2);
                if (dl && dl->type == EXPR_FUNCTION) {
                    for (size_t i = 0; i < dl->data.function.arg_count; i++) {
                        if (!expr_is_integer_like(dl->data.function.args[i])) continue;
                        mpz_t g; mpz_init(g); expr_to_mpz(dl->data.function.args[i], g);
                        for (int sg = 1; sg >= -1; sg -= 2) {
                            mpz_set(d1, g); if (sg < 0) mpz_neg(d1, d1);   /* d1 | D */
                            if (!mpz_divisible_p(D, d1)) continue;
                            mpz_divexact(d2, D, d1);
                            mpz_add(num, d1, d2);                          /* Y = (d1+d2)/(4p) */
                            if (!mpz_divisible_p(num, fourp)) continue;
                            mpz_divexact(Y, num, fourp);
                            mpz_sub(num, d2, d1); mpz_sub(num, num, twoB);  /* X = (d2-d1-2B)/(4A) */
                            if (!mpz_divisible_p(num, fourA)) continue;
                            mpz_divexact(X, num, fourA);
                            if (!mpz_fits_slong_p(Y) || !mpz_fits_slong_p(X)) continue;
                            int64_t vals[SI_MAX_VARS];
                            for (int k = 0; k < n; k++) vals[k] = 0;
                            vals[Yv] = mpz_get_si(Y); vals[Xv] = mpz_get_si(X);
                            if (si_verify(c, vals)) emit_full(st, vals);
                        }
                        mpz_clear(g);
                    }
                }
                if (dl) expr_free(dl);
                mpz_clear(d1); mpz_clear(d2); mpz_clear(Y); mpz_clear(X);
                mpz_clear(fourp); mpz_clear(fourA); mpz_clear(twoB); mpz_clear(num);
                mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(D); mpz_clear(p); mpz_clear(tmp);
                mpz_clear(e0); mpz_clear(e2); mpz_clear(e1); mpz_clear(ec);
                return true;                                   /* exhaustive: handled */
            }
            mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(D); mpz_clear(p); mpz_clear(tmp);
        }
        mpz_clear(e0); mpz_clear(e2); mpz_clear(e1); mpz_clear(ec);
    }
    return false;
}


/* --- Factorable binary quadratic  A x^2 + B x y + C y^2 + D x + E y + F == 0
 *     with a CROSS term and perfect-square discriminant delta = B^2 - 4 A C. ---
 *
 * When delta = k^2 > 0 the quadratic part factors into two distinct rational
 * linear forms, so the conic is a hyperbola with rational asymptotes and only
 * finitely many integer points -- Runge's simplest case.  Completing the square
 * (A != 0; swap x,y when A == 0) gives, with U = 2 A x + B y + D,
 *   4 A * eq = U^2 - k^2 y^2 + P y + Q,   P = 4 A E - 2 B D,  Q = 4 A F - D^2,
 * and multiplying by 4 k^2 with V = 2 k^2 y - P collapses it to a difference of
 * squares
 *   (2 k U)^2 - V^2 = W,   W = -(P^2 + 4 k^2 Q),
 * so every solution comes from a factorisation d1 * d2 = W via
 *   U = (d1 + d2) / (4 k),  V = (d2 - d1) / 2,
 *   y = (V + P) / (2 k^2),  x = (U - B y - D) / (2 A).
 * Enumerating the divisors of |W| is exhaustive, so an empty result is a PROOF
 * of no solutions (e.g. (x - y)(x + 2 y) == 15 -> {} by a mod-3 obstruction).
 * Generalises si_solve_conic to the cross-term case and to non-unit square
 * coefficients; the no-cross-term / Pell forms are handled before this.  W == 0
 * (two parallel/coincident lines: a positive-dimensional union) is declined. */
bool si_solve_factorable_conic(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return false;

    int Xv = act[0], Yv = act[1];
    mpz_t A, B, C, D, E, F;
    mpz_inits(A, B, C, D, E, F, NULL);
    bool shape = true;
    for (size_t t = 0; t < eq->n_terms && shape; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++)
            if (v != Xv && v != Yv && ex[v] != 0) { shape = false; break; }
        if (!shape) break;
        int ix = ex[Xv], iy = ex[Yv];
        if (ix == 2 && iy == 0) mpz_add(A, A, eq->coefs[t]);
        else if (ix == 1 && iy == 1) mpz_add(B, B, eq->coefs[t]);
        else if (ix == 0 && iy == 2) mpz_add(C, C, eq->coefs[t]);
        else if (ix == 1 && iy == 0) mpz_add(D, D, eq->coefs[t]);
        else if (ix == 0 && iy == 1) mpz_add(E, E, eq->coefs[t]);
        else if (ix == 0 && iy == 0) mpz_add(F, F, eq->coefs[t]);
        else shape = false;                               /* total degree > 2 */
    }

    /* Need A != 0; swap x <-> y when only C is nonzero.  Both zero -> bilinear
     * (handled elsewhere), decline. */
    if (shape && mpz_sgn(A) == 0 && mpz_sgn(C) != 0) {
        mpz_swap(A, C); mpz_swap(D, E); int tmp = Xv; Xv = Yv; Yv = tmp;
    }
    bool go = shape && mpz_sgn(A) != 0;

    mpz_t delta, k, P, Q, W, tmp;
    mpz_inits(delta, k, P, Q, W, tmp, NULL);
    if (go) {
        mpz_mul(delta, B, B); mpz_mul(tmp, A, C); mpz_submul_ui(delta, tmp, 4); /* B^2 - 4AC */
        go = (mpz_sgn(delta) > 0) && mpz_perfect_square_p(delta);
    }
    if (go) {
        mpz_sqrt(k, delta);
        mpz_mul(P, A, E); mpz_mul_ui(P, P, 4);
        mpz_mul(tmp, B, D); mpz_submul_ui(P, tmp, 2);      /* P = 4AE - 2BD */
        mpz_mul(Q, A, F); mpz_mul_ui(Q, Q, 4);
        mpz_mul(tmp, D, D); mpz_sub(Q, Q, tmp);            /* Q = 4AF - D^2 */
        mpz_mul(W, P, P);                                  /* W = -(P^2 + 4 k^2 Q) */
        mpz_mul(tmp, k, k); mpz_mul(tmp, tmp, Q); mpz_addmul_ui(W, tmp, 4);
        mpz_neg(W, W);
        go = (mpz_sgn(W) != 0);
    }

    bool handled = false;
    if (go) {
        handled = true;
        st->max_visits = SI_MAX_NODES;
        mpz_t absW; mpz_init(absW); mpz_abs(absW, W);
        Expr* dl = divisors_ordinary(absW);
        mpz_clear(absW);
        mpz_t d1, d2, U, V, x, y, fourk, twok2, num;
        mpz_inits(d1, d2, U, V, x, y, fourk, twok2, num, NULL);
        mpz_mul_ui(fourk, k, 4);
        mpz_mul(twok2, k, k); mpz_mul_ui(twok2, twok2, 2);
        if (dl && dl->type == EXPR_FUNCTION) {
            for (size_t i = 0; i < dl->data.function.arg_count && !st->overflow; i++) {
                if (!expr_is_integer_like(dl->data.function.args[i])) continue;
                mpz_t g; mpz_init(g); expr_to_mpz(dl->data.function.args[i], g);
                for (int sg = 1; sg >= -1; sg -= 2) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    mpz_set(d1, g); if (sg < 0) mpz_neg(d1, d1);   /* d1 | W */
                    if (!mpz_divisible_p(W, d1)) continue;
                    mpz_divexact(d2, W, d1);
                    mpz_add(num, d1, d2);                          /* U = (d1+d2)/(4k) */
                    if (!mpz_divisible_p(num, fourk)) continue;
                    mpz_divexact(U, num, fourk);
                    mpz_sub(V, d2, d1);                            /* V = (d2 - d1)/2 */
                    if (mpz_odd_p(V)) continue;
                    mpz_divexact_ui(V, V, 2);
                    mpz_add(num, V, P);                            /* y = (V+P)/(2k^2) */
                    if (!mpz_divisible_p(num, twok2)) continue;
                    mpz_divexact(y, num, twok2);
                    mpz_mul(num, B, y); mpz_sub(num, U, num); mpz_sub(num, num, D); /* x=(U-By-D)/(2A) */
                    mpz_mul_ui(tmp, A, 2);
                    if (!mpz_divisible_p(num, tmp)) continue;
                    mpz_divexact(x, num, tmp);
                    if (!mpz_fits_slong_p(x) || !mpz_fits_slong_p(y)) continue;
                    int64_t vals[SI_MAX_VARS];
                    for (int kk = 0; kk < n; kk++) vals[kk] = 0;
                    vals[Xv] = mpz_get_si(x); vals[Yv] = mpz_get_si(y);
                    if (si_verify(c, vals)) emit_full(st, vals);
                }
                mpz_clear(g);
            }
        }
        if (dl) expr_free(dl);
        mpz_clears(d1, d2, U, V, x, y, fourk, twok2, num, NULL);
    }

    mpz_clears(A, B, C, D, E, F, delta, k, P, Q, W, tmp, NULL);
    return handled;
}


/* --- Definite binary quadratic (ellipse)  A x^2 + B x y + C y^2 + D x + E y + F
 *     == 0  with discriminant  delta = B^2 - 4 A C < 0. ---
 *
 * A negative discriminant makes the quadratic part positive- (or negative-)
 * definite, so the conic is a compact ellipse with FINITELY many integer points
 * -- but a rotated one (B != 0) is not bounded by derive_bounds, which only
 * boxes sign-definite / even-only variables.  Treat the equation as a quadratic
 * in x for each fixed y:  A x^2 + (B y + D) x + (C y^2 + E y + F) == 0.  Its
 * x-discriminant  disc_x(y) = delta y^2 + (2 B D - 4 A E) y + (D^2 - 4 A F)  is a
 * DOWNWARD parabola in y (leading coefficient delta < 0), so real x exist only
 * for y in the finite interval between its roots.  Enumerate those y, solve the
 * integer quadratic in x exactly (perfect-square discriminant), and verify.
 * Exhaustive over the (bounded) ellipse, so an empty result is a PROOF. */
bool si_solve_elliptic_bqf(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return false;
    int Xv = act[0], Yv = act[1];

    mpz_t A, B, C, D, E, F; mpz_inits(A, B, C, D, E, F, NULL);
    bool shape = true;
    for (size_t t = 0; t < eq->n_terms && shape; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++)
            if (v != Xv && v != Yv && ex[v] != 0) { shape = false; break; }
        if (!shape) break;
        int ix = ex[Xv], iy = ex[Yv];
        if (ix == 2 && iy == 0) mpz_add(A, A, eq->coefs[t]);
        else if (ix == 1 && iy == 1) mpz_add(B, B, eq->coefs[t]);
        else if (ix == 0 && iy == 2) mpz_add(C, C, eq->coefs[t]);
        else if (ix == 1 && iy == 0) mpz_add(D, D, eq->coefs[t]);
        else if (ix == 0 && iy == 1) mpz_add(E, E, eq->coefs[t]);
        else if (ix == 0 && iy == 0) mpz_add(F, F, eq->coefs[t]);
        else shape = false;                                /* total degree > 2 */
    }

    mpz_t delta, tmp; mpz_inits(delta, tmp, NULL);
    if (shape) {
        mpz_mul(delta, B, B); mpz_mul(tmp, A, C); mpz_submul_ui(delta, tmp, 4);  /* B^2 - 4AC */
    }
    if (!shape || mpz_sgn(delta) >= 0) {                   /* not a definite ellipse */
        mpz_clears(A, B, C, D, E, F, delta, tmp, NULL); return false;
    }
    /* Normalise the leading coefficient positive (delta unchanged by negation). */
    if (mpz_sgn(A) < 0) {
        mpz_neg(A, A); mpz_neg(B, B); mpz_neg(C, C);
        mpz_neg(D, D); mpz_neg(E, E); mpz_neg(F, F);
    }

    /* y-range where disc_x(y) = delta y^2 + beta y + gamma >= 0. */
    mpz_t beta, gamma, inner; mpz_inits(beta, gamma, inner, NULL);
    mpz_mul(beta, B, D); mpz_mul_ui(beta, beta, 2);
    mpz_mul(tmp, A, E); mpz_submul_ui(beta, tmp, 4);       /* beta = 2BD - 4AE */
    mpz_mul(gamma, D, D); mpz_mul(tmp, A, F); mpz_submul_ui(gamma, tmp, 4); /* D^2 - 4AF */
    mpz_mul(inner, beta, beta); mpz_mul(tmp, delta, gamma);
    mpz_submul_ui(inner, tmp, 4);                          /* beta^2 - 4 delta gamma */

    bool handled = true;
    if (mpz_sgn(inner) < 0) {                              /* no real y -> {} proof */
        st->max_visits = SI_MAX_NODES;
    } else {
        double dDelta = mpz_get_d(delta), dBeta = mpz_get_d(beta), dInner = mpz_get_d(inner);
        double sq = sqrt(dInner < 0 ? 0 : dInner);
        double rA = (-dBeta + sq) / (2.0 * dDelta), rB = (-dBeta - sq) / (2.0 * dDelta);
        double dlo = (rA < rB ? rA : rB) - 2.0, dhi = (rA < rB ? rB : rA) + 2.0;
        if (!(dlo > -9e17 && dhi < 9e17) || (dhi - dlo) > (double)SI_MAX_NODES) {
            handled = false;                               /* too wide / overflow -> decline */
        } else {
            st->max_visits = SI_MAX_NODES;
            int64_t ylo = (int64_t)floor(dlo), yhi = (int64_t)ceil(dhi);
            mpz_t a1, a0, disc, root, num, twoA, yz, s2;
            mpz_inits(a1, a0, disc, root, num, twoA, yz, s2, NULL);
            mpz_mul_ui(twoA, A, 2);
            for (int64_t yy = ylo; yy <= yhi && !st->overflow; yy++) {
                mpz_set_si(yz, yy);
                /* a1 = B y + D,  a0 = C y^2 + E y + F. */
                mpz_mul(a1, B, yz); mpz_add(a1, a1, D);
                mpz_mul(a0, C, yz); mpz_add(a0, a0, E); mpz_mul(a0, a0, yz); mpz_add(a0, a0, F);
                mpz_mul(disc, a1, a1); mpz_mul(s2, A, a0); mpz_submul_ui(disc, s2, 4); /* a1^2 - 4A a0 */
                if (mpz_sgn(disc) < 0 || !mpz_perfect_square_p(disc)) continue;
                mpz_sqrt(root, disc);
                for (int sgn = 1; sgn >= -1; sgn -= 2) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    mpz_set(num, root); if (sgn < 0) mpz_neg(num, num);
                    mpz_sub(num, num, a1);                  /* -a1 +/- sqrt(disc) */
                    if (!mpz_divisible_p(num, twoA)) continue;
                    mpz_divexact(num, num, twoA);           /* x */
                    if (!mpz_fits_slong_p(num)) continue;
                    int64_t vals[SI_MAX_VARS];
                    for (int i = 0; i < n; i++) vals[i] = 0;
                    vals[Xv] = mpz_get_si(num); vals[Yv] = yy;
                    if (si_verify(c, vals)) emit_full(st, vals);
                    if (mpz_sgn(root) == 0) break;          /* double root: one x */
                }
            }
            mpz_clears(a1, a0, disc, root, num, twoA, yz, s2, NULL);
        }
    }

    mpz_clears(A, B, C, D, E, F, delta, tmp, beta, gamma, inner, NULL);
    return handled;
}
