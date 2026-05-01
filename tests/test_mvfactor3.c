/*
 * test_mvfactor3.c
 * ----------------
 * Unit tests for the trivariate Hensel pipeline (Phase F2 MVP).
 *
 * Coverage:
 *   - MPoly <-> ZUPoly conversion round-trip.
 *   - Bivariate Diophantine solver (mpoly_diophantine_2).
 *   - Trivariate two-factor Hensel lift (mpoly_hensel_lift_3_2).
 *
 * The tests exercise increasingly hard cases:
 *   1. (z + xy)(z - xy)  -- z is the main variable; bivariate factors
 *      are linear in z with coefficients in (x, y).
 *   2. (z + x)(z + y)    -- linear-in-z, but each factor depends on
 *      a different second-tier variable.
 *   3. (z^2 + x)(z + y)  -- mixed degrees in z.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "mpoly.h"
#include "mvfactor3.h"
#include "zupoly.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

#define TEST(name) do { printf("Running test: %s\n", #name); name(); } while (0)

/* ---------------------------------------------------------------- */
/*  Helpers                                                         */
/* ---------------------------------------------------------------- */

/* Parse, evaluate, and convert to MPoly with the named variables.
 * Caller frees with mpoly_free. */
static MPoly* parse_to_mpoly(const char* input, Expr** vars, int n_vars) {
    Expr* e = parse_expression(input);
    Expr* ev = evaluate(e);
    MPoly* p = expr_to_mpoly(ev, vars, n_vars);
    expr_free(e); expr_free(ev);
    return p;
}

/* ---------------------------------------------------------------- */
/*  Conversion tests                                                */
/* ---------------------------------------------------------------- */

static void test_mpoly_to_zupoly_roundtrip(void) {
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { x, y };

    /* Pure-x polynomial in 2-var storage. */
    MPoly* p = parse_to_mpoly("3 x^2 - 5 x + 7", vars, 2);
    ASSERT(p != NULL);

    ZUPoly* z = mpoly_to_zupoly_in(p, 0);
    ASSERT(z != NULL);
    ASSERT(z->deg == 2);
    ASSERT(mpz_cmp_si(z->c[0],  7) == 0);
    ASSERT(mpz_cmp_si(z->c[1], -5) == 0);
    ASSERT(mpz_cmp_si(z->c[2],  3) == 0);

    MPoly* p_back = zupoly_to_mpoly_in(z, 2, 0);
    ASSERT(mpoly_eq(p, p_back));

    mpoly_free(p); mpoly_free(p_back); zupoly_free(z);
    expr_free(x); expr_free(y);
}

static void test_mpoly_to_zupoly_rejects_multivariate(void) {
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { x, y };

    /* Genuinely bivariate -- has a y term. */
    MPoly* p = parse_to_mpoly("3 x^2 + y", vars, 2);
    ASSERT(p != NULL);

    /* Asking to extract as univariate-in-x must fail. */
    ZUPoly* z = mpoly_to_zupoly_in(p, 0);
    ASSERT(z == NULL);

    mpoly_free(p);
    expr_free(x); expr_free(y);
}

/* ---------------------------------------------------------------- */
/*  Bivariate Diophantine                                           */
/* ---------------------------------------------------------------- */

static void test_diophantine_2_simple(void) {
    /* U = x - y, V = x + y, both monic in x, coprime over Q(y).
     * Solve  Δu * V + Δv * U = E  with E = 1.
     * Standard Bezout: 1/(2y) * (x+y) - 1/(2y) * (x-y) = 1.
     * Doesn't have integer solution; expect failure. */
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { x, y };

    MPoly* U = parse_to_mpoly("x - y", vars, 2);
    MPoly* V = parse_to_mpoly("x + y", vars, 2);
    MPoly* E = mpoly_from_int(2, 1);

    MPoly* du = NULL; MPoly* dv = NULL;
    bool ok = mpoly_diophantine_2(U, V, E, 0, 1, &du, &dv);
    /* Z[x,y] solution doesn't exist (Bezout coefficients are rational
     * in y).  Expect false. */
    ASSERT(!ok);

    mpoly_free(U); mpoly_free(V); mpoly_free(E);
    expr_free(x); expr_free(y);
}

static void test_diophantine_2_quadratic_e(void) {
    /* U = x - 1, V = x + 1, E = -2 (constant).
     * Δu * (x+1) + Δv * (x-1) = -2.
     * deg(Δu) < deg(U) = 1 -> Δu = a constant.
     * deg(Δv) < deg(V) = 1 -> Δv = b constant.
     * (a + b) x + (a - b) = -2.
     * a + b = 0 -> b = -a.  a - (-a) = 2a = -2 -> a = -1.
     * So Δu = -1, Δv = 1. */
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { x, y };

    MPoly* U = parse_to_mpoly("x - 1", vars, 2);
    MPoly* V = parse_to_mpoly("x + 1", vars, 2);
    MPoly* E = parse_to_mpoly("-2", vars, 2);

    MPoly* du = NULL; MPoly* dv = NULL;
    bool ok = mpoly_diophantine_2(U, V, E, 0, 1, &du, &dv);
    ASSERT(ok);

    /* Verify Δu * V + Δv * U == E. */
    MPoly* duV = mpoly_mul(du, V);
    MPoly* dvU = mpoly_mul(dv, U);
    MPoly* sum = mpoly_add(duV, dvU);
    ASSERT(mpoly_eq(sum, E));

    mpoly_free(U); mpoly_free(V); mpoly_free(E);
    mpoly_free(du); mpoly_free(dv);
    mpoly_free(duV); mpoly_free(dvU); mpoly_free(sum);
    expr_free(x); expr_free(y);
}

static void test_diophantine_2_with_y_dependence(void) {
    /* U = z - x  (var_main=z=index 0, var_y=x=index 1)
     * V = z + x
     * E = -2 x^2  (this is the actual residual from the (z-xy)(z+xy)
     *              y-lift after one iteration).
     *
     * Solve (deg_z(Δu) < 1, deg_z(Δv) < 1):
     *   Δu (z + x) + Δv (z - x) = -2 x^2
     *   z coeff:  Δu + Δv = 0          -> Δv = -Δu
     *   const:    Δu * x  - Δv * x  = -2 x^2
     *             Δu * x + Δu * x = 2 x Δu = -2 x^2
     *             Δu = -x
     * So Δu = -x, Δv = x. */
    Expr* z = parse_expression("z");
    Expr* x = parse_expression("x");
    Expr* vars[] = { z, x };  /* index 0 = z, index 1 = x */

    MPoly* U = parse_to_mpoly("z - x", vars, 2);
    MPoly* V = parse_to_mpoly("z + x", vars, 2);
    MPoly* E = parse_to_mpoly("-2 x^2", vars, 2);

    MPoly* du = NULL; MPoly* dv = NULL;
    bool ok = mpoly_diophantine_2(U, V, E, 0, 1, &du, &dv);
    ASSERT(ok);

    /* Verify identity. */
    MPoly* duV = mpoly_mul(du, V);
    MPoly* dvU = mpoly_mul(dv, U);
    MPoly* sum = mpoly_add(duV, dvU);
    ASSERT(mpoly_eq(sum, E));

    /* Expected: du = -x, dv = x. */
    MPoly* expected_du = parse_to_mpoly("-x", vars, 2);
    MPoly* expected_dv = parse_to_mpoly("x",  vars, 2);
    ASSERT(mpoly_eq(du, expected_du));
    ASSERT(mpoly_eq(dv, expected_dv));

    mpoly_free(U); mpoly_free(V); mpoly_free(E);
    mpoly_free(du); mpoly_free(dv);
    mpoly_free(duV); mpoly_free(dvU); mpoly_free(sum);
    mpoly_free(expected_du); mpoly_free(expected_dv);
    expr_free(z); expr_free(x);
}

/* ---------------------------------------------------------------- */
/*  Trivariate two-factor Hensel lift                               */
/* ---------------------------------------------------------------- */

static void test_lift_3_z_plus_xy_times_z_minus_xy(void) {
    /* P(x, y, z) = (z + xy)(z - xy) = z^2 - x^2 y^2.
     *
     * Variables: index 0 = z, index 1 = x, index 2 = y.
     * Main = z; secondaries = x, y.
     *
     * Specialise y = 1 to get a bivariate (z, x) image:
     *   P|y=1 = z^2 - x^2 = (z - x)(z + x).
     *
     * Lift y from 1 back to general y via mpoly_hensel_lift_3_2 with
     * var_y = 1 (x), var_z = 2 (y), alpha_z = 1.
     * Hmm naming: the API's "var_z" is the lift variable.  Here we are
     * lifting the y dimension, so "var_z" = y_idx = 2, alpha_z = 1. */
    Expr* z = parse_expression("z");
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { z, x, y };

    MPoly* P = parse_to_mpoly("Expand[(z + x y)(z - x y)]", vars, 3);
    ASSERT(P != NULL);

    /* Bivariate seeds: (z - x) and (z + x).  These are P|y=1, factored. */
    MPoly* U_zx = parse_to_mpoly("z - x", vars, 3);
    MPoly* V_zx = parse_to_mpoly("z + x", vars, 3);

    MPoly* U = NULL; MPoly* V = NULL;
    bool ok = mpoly_hensel_lift_3_2(P, U_zx, V_zx,
                                     /*var_main*/ 0,
                                     /*var_y   */ 1,  /* x */
                                     /*var_z   */ 2,  /* y -- the lift var */
                                     /*alpha_z */ 1,
                                     &U, &V);
    ASSERT(ok);

    /* Verify U * V == P. */
    MPoly* prod = mpoly_mul(U, V);
    ASSERT(mpoly_eq(prod, P));

    /* Verify U|y=1 == U_zx and V|y=1 == V_zx. */
    MPoly* U_at_1 = mpoly_subst_var_int(U, 2, 1);
    MPoly* V_at_1 = mpoly_subst_var_int(V, 2, 1);
    ASSERT(mpoly_eq(U_at_1, U_zx));
    ASSERT(mpoly_eq(V_at_1, V_zx));

    /* Expected exact form: U = z - xy, V = z + xy. */
    MPoly* expected_U = parse_to_mpoly("z - x y", vars, 3);
    MPoly* expected_V = parse_to_mpoly("z + x y", vars, 3);
    ASSERT(mpoly_eq(U, expected_U));
    ASSERT(mpoly_eq(V, expected_V));

    mpoly_free(P); mpoly_free(U_zx); mpoly_free(V_zx);
    mpoly_free(U); mpoly_free(V); mpoly_free(prod);
    mpoly_free(U_at_1); mpoly_free(V_at_1);
    mpoly_free(expected_U); mpoly_free(expected_V);
    expr_free(z); expr_free(x); expr_free(y);
}

static void test_lift_3_z_plus_x_times_z_plus_y(void) {
    /* P = (z + x)(z + y) = z^2 + (x+y) z + xy.
     * At y = 0: P|y=0 = z^2 + x z = z(z + x).  But the seeds are
     *   U(z, x) = z + x and V(z, x) = z (from y = 0).
     * Hmm V is just z, monic in z, deg 0 in x.
     *
     * Try y = 1: P|y=1 = z^2 + (x+1) z + x = (z + x)(z + 1).
     * Seeds: U = z + x, V = z + 1.  Both monic in z.
     *
     * Lift y from 1 back. */
    Expr* z = parse_expression("z");
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { z, x, y };

    MPoly* P = parse_to_mpoly("Expand[(z + x)(z + y)]", vars, 3);
    ASSERT(P != NULL);

    MPoly* U_seed = parse_to_mpoly("z + x", vars, 3);
    MPoly* V_seed = parse_to_mpoly("z + 1", vars, 3);

    MPoly* U = NULL; MPoly* V = NULL;
    bool ok = mpoly_hensel_lift_3_2(P, U_seed, V_seed,
                                     /*var_main*/ 0, /*var_y*/ 1, /*var_z*/ 2,
                                     /*alpha_z*/ 1, &U, &V);
    ASSERT(ok);

    /* Verify U*V = P. */
    MPoly* prod = mpoly_mul(U, V);
    ASSERT(mpoly_eq(prod, P));

    /* Expected: U = z + x, V = z + y. */
    MPoly* expected_U = parse_to_mpoly("z + x", vars, 3);
    MPoly* expected_V = parse_to_mpoly("z + y", vars, 3);
    /* The lift should preserve U_seed at y=1 and adjust the y-dependence
     * appropriately.  Our seeds ordered (U_seed = z+x, V_seed = z+1).
     * Expected lifted: U stays z+x (no y-dependence), V grows from z+1
     * to z+y. */
    ASSERT(mpoly_eq(U, expected_U));
    ASSERT(mpoly_eq(V, expected_V));

    mpoly_free(P); mpoly_free(U_seed); mpoly_free(V_seed);
    mpoly_free(U); mpoly_free(V); mpoly_free(prod);
    mpoly_free(expected_U); mpoly_free(expected_V);
    expr_free(z); expr_free(x); expr_free(y);
}

static void test_lift_3_quadratic_in_z(void) {
    /* P = (z^2 + x)(z + y).  At y = 1: (z^2 + x)(z + 1).
     *  - U = z^2 + x  (monic in z, deg 2)
     *  - V = z + 1    (monic in z, deg 1)
     * Lift y back. */
    Expr* z = parse_expression("z");
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { z, x, y };

    MPoly* P = parse_to_mpoly("Expand[(z^2 + x)(z + y)]", vars, 3);
    ASSERT(P != NULL);

    MPoly* U_seed = parse_to_mpoly("z^2 + x", vars, 3);
    MPoly* V_seed = parse_to_mpoly("z + 1", vars, 3);

    MPoly* U = NULL; MPoly* V = NULL;
    bool ok = mpoly_hensel_lift_3_2(P, U_seed, V_seed,
                                     0, 1, 2, 1, &U, &V);
    ASSERT(ok);

    MPoly* prod = mpoly_mul(U, V);
    ASSERT(mpoly_eq(prod, P));

    MPoly* expected_U = parse_to_mpoly("z^2 + x", vars, 3);
    MPoly* expected_V = parse_to_mpoly("z + y", vars, 3);
    ASSERT(mpoly_eq(U, expected_U));
    ASSERT(mpoly_eq(V, expected_V));

    mpoly_free(P); mpoly_free(U_seed); mpoly_free(V_seed);
    mpoly_free(U); mpoly_free(V); mpoly_free(prod);
    mpoly_free(expected_U); mpoly_free(expected_V);
    expr_free(z); expr_free(x); expr_free(y);
}

static void test_lift_3_irreducible_returns_false(void) {
    /* P = z^2 + x + y -- bivariate-irreducible image at any alpha,
     * cannot be split as (linear)(linear) in z. */
    Expr* z = parse_expression("z");
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { z, x, y };

    MPoly* P = parse_to_mpoly("Expand[z^2 + x + y]", vars, 3);
    ASSERT(P != NULL);

    /* Try wrong seeds: U_seed = z, V_seed = z (just to confirm the
     * lift returns false on bad input). */
    MPoly* U_seed = parse_to_mpoly("z", vars, 3);
    MPoly* V_seed = parse_to_mpoly("z", vars, 3);

    MPoly* U = NULL; MPoly* V = NULL;
    bool ok = mpoly_hensel_lift_3_2(P, U_seed, V_seed, 0, 1, 2, 0, &U, &V);
    /* The seed product is z*z = z^2, but P|y=0 = z^2 + x.  Mismatch
     * at the seed; lift_3_2 should reject. */
    ASSERT(!ok);
    ASSERT(U == NULL && V == NULL);

    mpoly_free(P); mpoly_free(U_seed); mpoly_free(V_seed);
    expr_free(z); expr_free(x); expr_free(y);
}

/* ---------------------------------------------------------------- */
/*  Main                                                            */
/* ---------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    printf("Running mvfactor3 (trivariate Hensel) tests...\n");

    /* MPoly <-> ZUPoly */
    TEST(test_mpoly_to_zupoly_roundtrip);
    TEST(test_mpoly_to_zupoly_rejects_multivariate);

    /* Bivariate Diophantine */
    TEST(test_diophantine_2_simple);
    TEST(test_diophantine_2_quadratic_e);
    TEST(test_diophantine_2_with_y_dependence);

    /* Trivariate Hensel lift */
    TEST(test_lift_3_z_plus_xy_times_z_minus_xy);
    TEST(test_lift_3_z_plus_x_times_z_plus_y);
    TEST(test_lift_3_quadratic_in_z);
    TEST(test_lift_3_irreducible_returns_false);

    printf("All mvfactor3 tests passed!\n");
    return 0;
}
