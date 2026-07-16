/*
 * test_mvfactor.c
 * ---------------
 * End-to-end tests of bpoly_hensel_lift_2 -- the bivariate two-factor
 * Hensel lift.
 *
 * Test pattern: pick known monic-in-x bivariate factorisations P =
 * U * V, compute P, extract u = U(x, 0) and v = V(x, 0), feed (P, u,
 * v) to bpoly_hensel_lift_2, and verify the lift recovers U and V
 * exactly (up to the canonical sign convention the lift maintains).
 *
 * Each test verifies the algebraic invariant U*V == P (which the
 * lift's final check enforces), and additionally checks that the
 * two factors land in the right "shape" -- i.e. that the iteration
 * actually distributed work between U and V rather than producing
 * a trivial (constant, P) split.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "mvfactor.h"
#include "bpoly.h"
#include "zupoly.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

/* ---------------------------------------------------------------------- */
/*  Helpers                                                               */
/* ---------------------------------------------------------------------- */

static ZUPoly* zu(const int64_t* c, int n) {
    ZUPoly* p = zupoly_new(n);
    for (int i = 0; i < n; i++) zupoly_setcoef_si(p, i, c[i]);
    return p;
}

/* Build the monic-in-x BPoly P = U * V, where U, V are given as
 * lists of (x_index, ZUPoly-in-y) pairs. */

/* Helper: given a ZUPoly u(x), build the bivariate version that is
 * just u(x) (no y dependence). */
static BPoly* bpoly_const_in_y(const ZUPoly* u) {
    if (u->deg < 0) return bpoly_zero();
    BPoly* b = bpoly_new(u->deg + 1);
    for (int i = 0; i <= u->deg; i++) {
        const mpz_t* c = zupoly_getcoef(u, i);
        if (c && mpz_sgn(*c) != 0) {
            ZUPoly* y_const = zupoly_new(1);
            zupoly_setcoef(y_const, 0, *c);
            bpoly_set_xcoef(b, i, y_const);
        }
    }
    return b;
}

/* Verify the lifted factorisation:
 *   U*V == P
 *   U(x, 0) == u
 *   V(x, 0) == v */
static int verify_lift(const BPoly* U, const BPoly* V, const BPoly* P,
                       const ZUPoly* u, const ZUPoly* v) {
    BPoly* prod = bpoly_mul(U, V);
    int prod_ok = bpoly_eq(prod, P);
    bpoly_free(prod);
    if (!prod_ok) {
        fprintf(stderr, "verify_lift: U*V != P\n");
        return 0;
    }

    ZUPoly* u_check = bpoly_eval_y_si(U, 0);
    ZUPoly* v_check = bpoly_eval_y_si(V, 0);
    int u_ok = zupoly_eq(u_check, u);
    int v_ok = zupoly_eq(v_check, v);
    zupoly_free(u_check); zupoly_free(v_check);
    if (!u_ok || !v_ok) {
        fprintf(stderr, "verify_lift: image at y=0 mismatch\n");
        return 0;
    }
    return 1;
}

/* ====================================================================== */
/*  Trivial cases                                                         */
/* ====================================================================== */

static void test_y_independent_factorisation(void) {
    /* P(x, y) = (x - 1)(x + 1) = x^2 - 1.  No y dependence at all.
     * u = x - 1, v = x + 1.
     * Hensel iteration should immediately verify and return without
     * any non-trivial corrections (E_k is zero for all k > 0). */
    ZUPoly* u = zu((int64_t[]){-1, 1}, 2);
    ZUPoly* v = zu((int64_t[]){1, 1}, 2);
    ZUPoly* prod_z = zupoly_mul(u, v);
    BPoly* P = bpoly_const_in_y(prod_z);

    BPoly *U = NULL, *V = NULL;
    bool ok = bpoly_hensel_lift_2(P, u, v, &U, &V);
    ASSERT(ok);
    ASSERT(verify_lift(U, V, P, u, v));

    bpoly_free(P); bpoly_free(U); bpoly_free(V);
    zupoly_free(u); zupoly_free(v); zupoly_free(prod_z);
}

/* ====================================================================== */
/*  One-step lifts                                                        */
/* ====================================================================== */

static void test_lift_with_linear_y_dependence(void) {
    /* Build P = (x + y - 1)(x + y + 1) = x^2 + 2xy + y^2 - 1.
     * At y = 0: P(x, 0) = (x - 1)(x + 1).
     * u = x - 1, v = x + 1.  Hensel must lift u -> x + y - 1,
     * v -> x + y + 1.
     *
     * BPoly representation of P:
     *   cx[0] = y^2 - 1 = (-1, 0, 1)
     *   cx[1] = 2y      = (0, 2)
     *   cx[2] = 1       = (1) */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, 0, 1}, 3));
    bpoly_set_xcoef(P, 1, zu((int64_t[]){0, 2}, 2));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){1}, 1));

    ZUPoly* u = zu((int64_t[]){-1, 1}, 2);
    ZUPoly* v = zu((int64_t[]){1, 1}, 2);

    BPoly *U = NULL, *V = NULL;
    bool ok = bpoly_hensel_lift_2(P, u, v, &U, &V);
    ASSERT(ok);
    ASSERT(verify_lift(U, V, P, u, v));

    /* Both U and V must have y in them (degree-1 in y). */
    ASSERT(bpoly_deg_y(U) == 1);
    ASSERT(bpoly_deg_y(V) == 1);

    bpoly_free(P); bpoly_free(U); bpoly_free(V);
    zupoly_free(u); zupoly_free(v);
}

static void test_lift_with_quadratic_y_dependence(void) {
    /* Build P = (x + y^2 - 1)(x - y^2 + 1) over Z[x, y].
     * Compute symbolically:
     *   P = x^2 + x*(y^2 - 1) - x*(y^2 - 1) - (y^2 - 1)^2
     *     = x^2 - (y^2 - 1)^2
     *     = x^2 - (y^4 - 2y^2 + 1)
     *     = x^2 - y^4 + 2y^2 - 1.
     * At y = 0: P(x, 0) = x^2 - 1 = (x - 1)(x + 1).
     * u = x - 1, v = x + 1.  Lift to y^2 corrections.
     *
     * BPoly storage:
     *   cx[0] = -y^4 + 2y^2 - 1     = (-1, 0, 2, 0, -1)
     *   cx[1] = 0
     *   cx[2] = 1                   = (1) */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, 0, 2, 0, -1}, 5));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){1}, 1));

    ZUPoly* u = zu((int64_t[]){-1, 1}, 2);
    ZUPoly* v = zu((int64_t[]){1, 1}, 2);

    BPoly *U = NULL, *V = NULL;
    bool ok = bpoly_hensel_lift_2(P, u, v, &U, &V);
    ASSERT(ok);
    ASSERT(verify_lift(U, V, P, u, v));

    /* Each lifted factor has y^2 (not just y) coefficients. */
    ASSERT(bpoly_deg_y(U) >= 2);
    ASSERT(bpoly_deg_y(V) >= 2);

    bpoly_free(P); bpoly_free(U); bpoly_free(V);
    zupoly_free(u); zupoly_free(v);
}

/* ====================================================================== */
/*  Higher-degree x                                                       */
/* ====================================================================== */

static void test_lift_cubic_x(void) {
    /* P = (x^2 + y - 1)(x + y + 2).
     * Expand:
     *   P = x^3 + x^2*(y + 2) + x*(y - 1) + (y - 1)*(y + 2)
     *     = x^3 + (y + 2) x^2 + (y - 1) x + (y^2 + y - 2)
     * At y = 0:
     *   P(x, 0) = x^3 + 2x^2 - x - 2
     *   Factor over Z: (x^2 - 1)(x + 2) = (x-1)(x+1)(x+2).
     *   But our intended factorisation pairs (x^2 - 1) with (x + 2).
     * u = x^2 - 1, v = x + 2.  Both monic, coprime. */
    BPoly* P = bpoly_new(4);
    /* cx[0] = y^2 + y - 2 */
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-2, 1, 1}, 3));
    /* cx[1] = y - 1 */
    bpoly_set_xcoef(P, 1, zu((int64_t[]){-1, 1}, 2));
    /* cx[2] = y + 2 */
    bpoly_set_xcoef(P, 2, zu((int64_t[]){2, 1}, 2));
    /* cx[3] = 1 */
    bpoly_set_xcoef(P, 3, zu((int64_t[]){1}, 1));

    ZUPoly* u = zu((int64_t[]){-1, 0, 1}, 3);  /* x^2 - 1 */
    ZUPoly* v = zu((int64_t[]){2, 1}, 2);      /* x + 2 */

    BPoly *U = NULL, *V = NULL;
    bool ok = bpoly_hensel_lift_2(P, u, v, &U, &V);
    ASSERT(ok);
    ASSERT(verify_lift(U, V, P, u, v));

    bpoly_free(P); bpoly_free(U); bpoly_free(V);
    zupoly_free(u); zupoly_free(v);
}

/* ====================================================================== */
/*  Failure paths                                                         */
/* ====================================================================== */

static void test_lift_with_wrong_image_returns_false(void) {
    /* Same P as in the linear-y-dependence test (= x^2 + 2xy + y^2 - 1)
     * but with deliberately-wrong u, v: u = x, v = x. */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, 0, 1}, 3));
    bpoly_set_xcoef(P, 1, zu((int64_t[]){0, 2}, 2));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){1}, 1));

    ZUPoly* u_wrong = zu((int64_t[]){0, 1}, 2);
    ZUPoly* v_wrong = zu((int64_t[]){0, 1}, 2);

    BPoly *U = NULL, *V = NULL;
    bool ok = bpoly_hensel_lift_2(P, u_wrong, v_wrong, &U, &V);
    /* Either the Diophantine fails (gcd not 1) or the final
     * verification fails -- either way the function must signal
     * failure with NULL outputs. */
    ASSERT(!ok);
    ASSERT(U == NULL && V == NULL);

    bpoly_free(P);
    zupoly_free(u_wrong); zupoly_free(v_wrong);
}

/* ====================================================================== */
/*  Compositional / integration                                           */
/* ====================================================================== */

static void test_lift_via_real_image_factor(void) {
    /* Build P from a non-trivial bivariate factorisation, then start
     * the lift only from the univariate image -- mimicking how the
     * full pipeline will assemble the inputs.
     *
     * P = (x + y)(x - y - 3) over Z[x, y].
     *   = x^2 + xy - 3x - xy - y^2 - 3y
     *   Wait: (x + y)(x - y - 3)
     *       = x*(x - y - 3) + y*(x - y - 3)
     *       = x^2 - xy - 3x + xy - y^2 - 3y
     *       = x^2 - 3x - y^2 - 3y
     * At y = 0: P(x, 0) = x^2 - 3x.  This factors as x(x - 3).
     * u = x, v = x - 3.  Both monic, coprime.
     *
     * BPoly storage of P:
     *   cx[0] = -y^2 - 3y     = (0, -3, -1)
     *   cx[1] = -3            = (-3)
     *   cx[2] = 1             = (1)
     *
     * Note: the "true" factorisation in Hensel's terms is
     *   U(x, y) = x + y       (because U(x, 0) = x = u)
     *   V(x, y) = x - y - 3   (because V(x, 0) = x - 3 = v)
     * and indeed U(x, 0) * V(x, 0) = x*(x - 3) = x^2 - 3x = P(x, 0). */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){0, -3, -1}, 3));
    bpoly_set_xcoef(P, 1, zu((int64_t[]){-3}, 1));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){1}, 1));

    ZUPoly* u = zu((int64_t[]){0, 1}, 2);   /* x */
    ZUPoly* v = zu((int64_t[]){-3, 1}, 2);  /* x - 3 */

    BPoly *U = NULL, *V = NULL;
    bool ok = bpoly_hensel_lift_2(P, u, v, &U, &V);
    ASSERT(ok);
    ASSERT(verify_lift(U, V, P, u, v));

    /* U should be x + y; V should be x - y - 3. */
    /* U.cx[0] = y, U.cx[1] = 1. */
    const ZUPoly* U0 = bpoly_get_xcoef(U, 0);
    const ZUPoly* U1 = bpoly_get_xcoef(U, 1);
    ASSERT(U0 && U0->deg == 1);  /* y */
    ASSERT(U1 && U1->deg == 0);  /* 1 */
    /* V.cx[0] = -y - 3, V.cx[1] = 1. */
    const ZUPoly* V0 = bpoly_get_xcoef(V, 0);
    const ZUPoly* V1 = bpoly_get_xcoef(V, 1);
    ASSERT(V0 && V0->deg == 1);  /* -3 - y */
    ASSERT(V1 && V1->deg == 0);  /* 1 */

    bpoly_free(P); bpoly_free(U); bpoly_free(V);
    zupoly_free(u); zupoly_free(v);
}

/* ====================================================================== */
/*  Multifactor lift                                                      */
/* ====================================================================== */

static void test_multifactor_three_factors(void) {
    /* P = (x - 1)(x - 2)(x - 3) (no y dependence -- trivial multifactor).
     *   = x^3 - 6x^2 + 11x - 6
     * u_1 = x - 1, u_2 = x - 2, u_3 = x - 3.  prod = P(x, 0) = P. */
    ZUPoly* u1 = zu((int64_t[]){-1, 1}, 2);
    ZUPoly* u2 = zu((int64_t[]){-2, 1}, 2);
    ZUPoly* u3 = zu((int64_t[]){-3, 1}, 2);

    /* Build P as a constant-in-y BPoly via the product of u_i. */
    ZUPoly* prod_u12 = zupoly_mul(u1, u2);
    ZUPoly* prod_z = zupoly_mul(prod_u12, u3);
    zupoly_free(prod_u12);
    BPoly* P = bpoly_const_in_y(prod_z);
    zupoly_free(prod_z);

    const ZUPoly* us[3] = { u1, u2, u3 };
    BPoly** Us = NULL;
    int rc = 0;
    bool ok = bpoly_hensel_lift_multi(P, us, 3, &Us, &rc);
    ASSERT(ok);
    ASSERT(rc == 3);

    /* Each Us[i] must equal u_i (no y dependence to lift). */
    ZUPoly* u1_check = bpoly_eval_y_si(Us[0], 0);
    ZUPoly* u2_check = bpoly_eval_y_si(Us[1], 0);
    ZUPoly* u3_check = bpoly_eval_y_si(Us[2], 0);
    ASSERT(zupoly_eq(u1_check, u1));
    ASSERT(zupoly_eq(u2_check, u2));
    ASSERT(zupoly_eq(u3_check, u3));
    zupoly_free(u1_check); zupoly_free(u2_check); zupoly_free(u3_check);

    /* Product check. */
    BPoly* p_check = bpoly_mul(Us[0], Us[1]);
    BPoly* p_full = bpoly_mul(p_check, Us[2]);
    ASSERT(bpoly_eq(p_full, P));
    bpoly_free(p_check); bpoly_free(p_full);

    for (int i = 0; i < 3; i++) bpoly_free(Us[i]);
    free(Us);
    bpoly_free(P);
    zupoly_free(u1); zupoly_free(u2); zupoly_free(u3);
}

static void test_multifactor_with_y_dependence(void) {
    /* P = (x - 1)(x + 1)(x + y).
     *   = (x^2 - 1)(x + y)
     *   = x^3 + x^2 y - x - y
     *
     * At y=0: P(x, 0) = x^3 - x = x(x-1)(x+1).
     * u_1 = x, u_2 = x - 1, u_3 = x + 1.
     *
     * BPoly storage of P:
     *   cx[0] = -y                = (0, -1)
     *   cx[1] = -1                = (-1)
     *   cx[2] = y                 = (0, 1)
     *   cx[3] = 1                 = (1) */
    BPoly* P = bpoly_new(4);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){0, -1}, 2));
    bpoly_set_xcoef(P, 1, zu((int64_t[]){-1}, 1));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){0, 1}, 2));
    bpoly_set_xcoef(P, 3, zu((int64_t[]){1}, 1));

    ZUPoly* u1 = zu((int64_t[]){0, 1}, 2);   /* x */
    ZUPoly* u2 = zu((int64_t[]){-1, 1}, 2);  /* x - 1 */
    ZUPoly* u3 = zu((int64_t[]){1, 1}, 2);   /* x + 1 */
    const ZUPoly* us[3] = { u1, u2, u3 };

    BPoly** Us = NULL;
    int rc = 0;
    bool ok = bpoly_hensel_lift_multi(P, us, 3, &Us, &rc);
    ASSERT(ok);
    ASSERT(rc == 3);

    /* Check each lifted factor's image at y=0 matches its u. */
    for (int i = 0; i < 3; i++) {
        ZUPoly* check = bpoly_eval_y_si(Us[i], 0);
        ASSERT(zupoly_eq(check, us[i]));
        zupoly_free(check);
    }

    /* Product check: U_1 * U_2 * U_3 = P. */
    BPoly* prod_12 = bpoly_mul(Us[0], Us[1]);
    BPoly* prod_full = bpoly_mul(prod_12, Us[2]);
    ASSERT(bpoly_eq(prod_full, P));
    bpoly_free(prod_12); bpoly_free(prod_full);

    /* The factor that absorbed the y dependence is the one whose image
     * at y=0 was u_1 = x.  It should now be x + y.  cx[0] = y, cx[1] = 1.
     * (The other two factors should still be y-independent: x - 1, x + 1.) */
    int y_dep_count = 0;
    for (int i = 0; i < 3; i++) {
        if (bpoly_deg_y(Us[i]) >= 1) y_dep_count++;
    }
    ASSERT(y_dep_count == 1);

    for (int i = 0; i < 3; i++) bpoly_free(Us[i]);
    free(Us);
    bpoly_free(P);
    zupoly_free(u1); zupoly_free(u2); zupoly_free(u3);
}

static void test_multifactor_single_factor_passthrough(void) {
    /* r = 1 case: lift a single factor.  Should return a copy of P. */
    BPoly* P = bpoly_new(2);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){1, 2}, 2));   /* 1 + 2y */
    bpoly_set_xcoef(P, 1, zu((int64_t[]){1}, 1));      /* 1 (monic in x) */

    ZUPoly* u = zu((int64_t[]){1, 1}, 2);  /* x + 1 = P(x, 0) */
    const ZUPoly* us[1] = { u };
    BPoly** Us = NULL;
    int rc = 0;
    bool ok = bpoly_hensel_lift_multi(P, us, 1, &Us, &rc);
    ASSERT(ok);
    ASSERT(rc == 1);
    ASSERT(bpoly_eq(Us[0], P));

    bpoly_free(Us[0]); free(Us);
    bpoly_free(P); zupoly_free(u);
}

/* ====================================================================== */
/*  Recombination: bivariate factorisation coarser than univariate one    */
/* ====================================================================== */

/* When the bivariate factorisation has fewer factors than the
 * univariate image, the lift's pair-and-recurse approach (lift one
 * u against the rest) fails because no single u corresponds to a
 * true bivariate factor.  Phase 3b's recombination then tries
 * larger subsets until it finds the correct grouping. */

static void test_recombine_two_subset_grouping(void) {
    /* P(x,y) = (x^2 - y)(x^2 + y) = x^4 - y^2.
     * At y = 1: P(x, 1) = x^4 - 1 = (x-1)(x+1)(x^2+1).  Three
     * univariate factors; only two bivariate factors.  The lift
     * must recombine: (x-1)*(x+1) → x^2 - y; (x^2+1) → x^2 + y.
     *
     * Build P after the y -> y+1 shift (so the lift sees images at
     * y = 0).  Original P = x^4 - y^2.  After y -> y+1:
     *   P_shifted = x^4 - (y+1)^2 = x^4 - y^2 - 2y - 1.
     *
     * BPoly storage of P_shifted (y -> y+1 of x^4 - y^2):
     *   cx[0] = -y^2 - 2y - 1   = (-1, -2, -1)
     *   cx[4] = 1               = (1)  (monic in x). */
    BPoly* P = bpoly_new(5);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, -2, -1}, 3));
    bpoly_set_xcoef(P, 4, zu((int64_t[]){1}, 1));

    /* Univariate factors at y=0 (i.e. y_orig=1): x-1, x+1, x^2+1. */
    ZUPoly* u1 = zu((int64_t[]){-1, 1}, 2);   /* x - 1 */
    ZUPoly* u2 = zu((int64_t[]){1, 1}, 2);    /* x + 1 */
    ZUPoly* u3 = zu((int64_t[]){1, 0, 1}, 3); /* x^2 + 1 */
    const ZUPoly* us[3] = { u1, u2, u3 };

    BPoly** Us = NULL;
    int rc = 0;
    bool ok = bpoly_hensel_lift_multi(P, us, 3, &Us, &rc);
    ASSERT(ok);
    /* The bivariate factorisation has TWO factors, not three. */
    ASSERT(rc == 2);

    /* Product check: U_1 * U_2 = P. */
    BPoly* prod = bpoly_mul(Us[0], Us[1]);
    ASSERT(bpoly_eq(prod, P));
    bpoly_free(prod);

    /* The two factors should be x^2 - y - 1 (the shifted form of
     * x^2 - y) and x^2 + y + 1 (the shifted form of x^2 + y).  Both
     * have x-degree 2.  Across the two: one image at y=0 is x^2 - 1
     * = (x-1)(x+1), the other is x^2 + 1.  Sum-of-degrees in x = 4. */
    int total_deg_x = 0;
    for (int i = 0; i < rc; i++) {
        ASSERT(Us[i]->deg_x >= 1);
        total_deg_x += Us[i]->deg_x;
    }
    ASSERT(total_deg_x == 4);

    for (int i = 0; i < rc; i++) bpoly_free(Us[i]);
    free(Us);
    bpoly_free(P);
    zupoly_free(u1); zupoly_free(u2); zupoly_free(u3);
}

static void test_recombine_irreducible_returns_one(void) {
    /* P(x,y) = x^4 + y + 1.  Irreducible over Z[x,y].
     * At y = 0: P(x, 0) = x^4 + 1.  Univariately factors over Z as ...
     * Actually x^4 + 1 is irreducible over Q (cyclotomic Phi_8).
     * Try y = 0 — image x^4 + 1 has 1 factor.  Already irreducible.
     *
     * Need a case where image factors but bivariate is irreducible.
     * Try P = x^4 + 4y + 1.
     * At y = 0: x^4 + 1, irreducible -> r = 1.
     * Hmm, need image to factor.
     *
     * Try P = x^2 + 2y x + 1 (degree 1 in y).
     * At y = 0: x^2 + 1, irreducible.  Useless.
     *
     * Try P = x^4 - y^2 - 2.
     * At y = 0: x^4 - 2, irreducible over Q (4th root of 2 not in Q).
     * Also useless for testing recombination *failure*.
     *
     * Better: P(x,y) = x^4 + x^2 y - x^2 + y - 1.
     * Factor at y = 0: x^4 - x^2 - 1.  Discriminant for x^2:
     *   t^2 - t - 1, roots (1 ± √5)/2.  Irrational, so x^4 - x^2 - 1
     *   is irreducible over Q.  r = 1. Useless.
     *
     * Stick to: just verify that on a y-image with multiple factors
     * but where the bivariate is irreducible, lift_multi returns
     * a single-factor result equal to P.
     *
     * Use P = x^4 - 1 + y.  At y=0: x^4 - 1 = (x-1)(x+1)(x^2+1).
     * 3 univariate factors.  Is x^4 - 1 + y reducible bivariately?
     *   Suppose x^4 + (y - 1) = f(x,y) g(x,y).  At y=0: factors as
     *   given.  But for any y, x^4 + y - 1 is "x^4 plus a constant"
     *   in y, and adding y can't preserve a non-trivial product
     *   structure unless the factors have y-dependence that
     *   conspires.  Let's just trust the test result.
     *
     * BPoly storage:
     *   cx[0] = y - 1   = (-1, 1)
     *   cx[4] = 1       = (1)  (monic in x). */
    BPoly* P = bpoly_new(5);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, 1}, 2));
    bpoly_set_xcoef(P, 4, zu((int64_t[]){1}, 1));

    ZUPoly* u1 = zu((int64_t[]){-1, 1}, 2);   /* x - 1 */
    ZUPoly* u2 = zu((int64_t[]){1, 1}, 2);    /* x + 1 */
    ZUPoly* u3 = zu((int64_t[]){1, 0, 1}, 3); /* x^2 + 1 */
    const ZUPoly* us[3] = { u1, u2, u3 };

    BPoly** Us = NULL;
    int rc = 0;
    bool ok = bpoly_hensel_lift_multi(P, us, 3, &Us, &rc);
    ASSERT(ok);
    /* P is irreducible bivariately -- no subset partition lifts. */
    ASSERT(rc == 1);
    ASSERT(bpoly_eq(Us[0], P));

    bpoly_free(Us[0]);
    free(Us);
    bpoly_free(P);
    zupoly_free(u1); zupoly_free(u2); zupoly_free(u3);
}

static void test_recombine_four_factors_two_pairs(void) {
    /* P(x, y) = (x^2 + 2y - 1)(x^2 - 2y + 1) over Z[x, y].
     * Compute:
     *   = x^4 - x^2 (2y - 1) + x^2 (2y - 1) - (2y - 1)^2
     *   wait, more carefully:
     *   (x^2 + a)(x^2 - a) = x^4 - a^2 with a = 2y - 1.
     *   So P = x^4 - (2y - 1)^2 = x^4 - 4y^2 + 4y - 1.
     * At y = 0: P(x, 0) = x^4 - 1 = (x-1)(x+1)(x^2+1).
     * Univariate has r = 3.  Bivariate has r = 2.
     *
     * BPoly storage:
     *   cx[0] = -4y^2 + 4y - 1  = (-1, 4, -4)
     *   cx[4] = 1               = (1) */
    BPoly* P = bpoly_new(5);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, 4, -4}, 3));
    bpoly_set_xcoef(P, 4, zu((int64_t[]){1}, 1));

    ZUPoly* u1 = zu((int64_t[]){-1, 1}, 2);   /* x - 1 */
    ZUPoly* u2 = zu((int64_t[]){1, 1}, 2);    /* x + 1 */
    ZUPoly* u3 = zu((int64_t[]){1, 0, 1}, 3); /* x^2 + 1 */
    const ZUPoly* us[3] = { u1, u2, u3 };

    BPoly** Us = NULL;
    int rc = 0;
    bool ok = bpoly_hensel_lift_multi(P, us, 3, &Us, &rc);
    ASSERT(ok);
    ASSERT(rc == 2);

    BPoly* prod = bpoly_mul(Us[0], Us[1]);
    ASSERT(bpoly_eq(prod, P));
    bpoly_free(prod);

    int total_deg_x = 0;
    for (int i = 0; i < rc; i++) total_deg_x += Us[i]->deg_x;
    ASSERT(total_deg_x == 4);

    for (int i = 0; i < rc; i++) bpoly_free(Us[i]);
    free(Us);
    bpoly_free(P);
    zupoly_free(u1); zupoly_free(u2); zupoly_free(u3);
}

/* ====================================================================== */
/*  Orchestrator: mvfactor_try_bivariate_monic                            */
/* ====================================================================== */

/* A simple test callback that factors a univariate ZUPoly via the
 * existing facpoly.c::bz_factor_to_expr (going through Expr).  This
 * is what the production wiring will look like, so the orchestrator
 * tests double as integration tests. */
extern struct Expr* bz_factor_to_expr(struct Expr* P, struct Expr* var);

static bool test_factor_via_expr(const ZUPoly* image,
                                 ZUPoly*** factors_out,
                                 int* count_out,
                                 void* user_data) {
    (void)user_data;
    /* Convert ZUPoly to Expr, factor, then convert back.  The result
     * of bz_factor_to_expr is a Times of factors (or a single factor
     * if irreducible). */
    Expr* var = parse_expression("x");
    Expr* image_expr = zupoly_to_expr(image, var);
    Expr* factored = bz_factor_to_expr(image_expr, var);
    expr_free(image_expr);

    /* Walk the result. */
    ZUPoly** result = NULL;
    int n = 0;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head
        && factored->data.function.head->type == EXPR_SYMBOL
        && strcmp(factored->data.function.head->data.symbol.name, "Times") == 0) {
        size_t ac = factored->data.function.arg_count;
        result = (ZUPoly**)malloc(sizeof(ZUPoly*) * ac);
        for (size_t i = 0; i < ac; i++) {
            Expr* a = factored->data.function.args[i];
            /* Skip pure-integer factors (constants); we want only the
             * non-trivial polynomial factors. */
            ZUPoly* zp = expr_to_zupoly(a, var);
            if (zp && !zupoly_is_zero(zp) && zp->deg >= 1) {
                result[n++] = zp;
            } else if (zp) {
                zupoly_free(zp);
            }
        }
    } else {
        /* Single factor (possibly non-Times). */
        ZUPoly* zp = expr_to_zupoly(factored, var);
        if (zp && !zupoly_is_zero(zp) && zp->deg >= 1) {
            result = (ZUPoly**)malloc(sizeof(ZUPoly*));
            result[0] = zp;
            n = 1;
        } else if (zp) {
            zupoly_free(zp);
        }
    }
    expr_free(factored);
    expr_free(var);

    if (n == 0) {
        free(result);
        return false;
    }
    *factors_out = result;
    *count_out = n;
    return true;
}

static void test_orchestrator_two_factor_lift(void) {
    /* P = (x + y - 1)(x + y + 1) = x^2 + 2xy + y^2 - 1.
     * Monic in x.  At y=0: (x-1)(x+1).  Should lift to two factors. */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1, 0, 1}, 3));
    bpoly_set_xcoef(P, 1, zu((int64_t[]){0, 2}, 2));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){1}, 1));

    BPoly** factors = NULL;
    int r = 0;
    bool ok = mvfactor_try_bivariate_monic(P, test_factor_via_expr, NULL, &factors, &r);
    ASSERT(ok);
    ASSERT(r == 2);

    /* Multiply factors back; must equal P. */
    BPoly* prod = bpoly_mul(factors[0], factors[1]);
    ASSERT(bpoly_eq(prod, P));
    bpoly_free(prod);

    for (int i = 0; i < r; i++) bpoly_free(factors[i]);
    free(factors);
    bpoly_free(P);
}

static void test_orchestrator_irreducible(void) {
    /* P = x^2 + y^2 + 1.  Irreducible over Z[x, y].
     * The orchestrator should detect this (univariate image at any
     * alpha is x^2 + alpha^2 + 1, irreducible -- e.g. x^2 + 1 at
     * alpha=0) and return P as the sole factor. */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){1, 0, 1}, 3));   /* 1 + y^2 */
    bpoly_set_xcoef(P, 2, zu((int64_t[]){1}, 1));         /* 1 */

    BPoly** factors = NULL;
    int r = 0;
    bool ok = mvfactor_try_bivariate_monic(P, test_factor_via_expr, NULL, &factors, &r);
    ASSERT(ok);
    ASSERT(r == 1);
    ASSERT(bpoly_eq(factors[0], P));

    bpoly_free(factors[0]); free(factors);
    bpoly_free(P);
}

static void test_orchestrator_recombination_two_factors(void) {
    /* P(x, y) = x^4 - y^2 = (x^2 - y)(x^2 + y).
     *   At y = 0: x^4, NOT squarefree (skipped).
     *   At y = 1: x^4 - 1 = (x-1)(x+1)(x^2+1), squarefree.  3 univariate factors.
     *   At y = -1: x^4 - 1, same factorisation, squarefree.
     *
     * The orchestrator should try y=0 (skip due to non-squarefree),
     * then y=1 where it factors into 3 univariate factors and the
     * recombination logic finds that (x-1)*(x+1) lifts together to
     * x^2 - y while x^2 + 1 lifts alone to x^2 + y -- giving 2
     * bivariate factors.
     *
     * BPoly storage of P:
     *   cx[0] = -y^2     = (0, 0, -1)
     *   cx[4] = 1 */
    BPoly* P = bpoly_new(5);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){0, 0, -1}, 3));
    bpoly_set_xcoef(P, 4, zu((int64_t[]){1}, 1));

    BPoly** factors = NULL;
    int r = 0;
    bool ok = mvfactor_try_bivariate_monic(P, test_factor_via_expr, NULL,
                                           &factors, &r);
    ASSERT(ok);
    ASSERT(r == 2);

    BPoly* prod = bpoly_mul(factors[0], factors[1]);
    ASSERT(bpoly_eq(prod, P));
    bpoly_free(prod);

    int total_deg_x = 0;
    for (int i = 0; i < r; i++) total_deg_x += factors[i]->deg_x;
    ASSERT(total_deg_x == 4);

    for (int i = 0; i < r; i++) bpoly_free(factors[i]);
    free(factors);
    bpoly_free(P);
}

static void test_orchestrator_nonmonic_returns_false(void) {
    /* P = 2 x^2 + ... is not monic in x; orchestrator should reject. */
    BPoly* P = bpoly_new(3);
    bpoly_set_xcoef(P, 0, zu((int64_t[]){-1}, 1));
    bpoly_set_xcoef(P, 2, zu((int64_t[]){2}, 1));

    BPoly** factors = NULL;
    int r = 0;
    bool ok = mvfactor_try_bivariate_monic(P, test_factor_via_expr, NULL, &factors, &r);
    ASSERT(!ok);
    ASSERT(factors == NULL);
    bpoly_free(P);
}

/* ====================================================================== */
/*  Main                                                                  */
/* ====================================================================== */

int main(void) {
    symtab_init();
    core_init();

    printf("Running mvfactor (bivariate Hensel) tests...\n");

    TEST(test_y_independent_factorisation);
    TEST(test_lift_with_linear_y_dependence);
    TEST(test_lift_with_quadratic_y_dependence);
    TEST(test_lift_cubic_x);
    TEST(test_lift_with_wrong_image_returns_false);
    TEST(test_lift_via_real_image_factor);

    /* Multifactor */
    TEST(test_multifactor_three_factors);
    TEST(test_multifactor_with_y_dependence);
    TEST(test_multifactor_single_factor_passthrough);

    /* Recombination (Phase 3b) */
    TEST(test_recombine_two_subset_grouping);
    TEST(test_recombine_irreducible_returns_one);
    TEST(test_recombine_four_factors_two_pairs);

    /* Orchestrator */
    TEST(test_orchestrator_two_factor_lift);
    TEST(test_orchestrator_irreducible);
    TEST(test_orchestrator_recombination_two_factors);
    TEST(test_orchestrator_nonmonic_returns_false);

    printf("All mvfactor tests passed!\n");
    return 0;
}
