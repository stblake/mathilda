/*
 * test_dsolve_m17_stress.c — anti-overfit stress families for M17
 * (DSolve`SpecialFunctionForm: rational-coefficient hypergeometric-class 2nd-order
 * linear ODEs with two finite regular singular points, reduced to Hypergeometric2F1
 * by an affine map onto x(1-x) plus a local-exponent shift Y = s^r0 (1-s)^r1 F; and
 * the Liouville normal-form pre-pass that lets the Airy/Bessel recognisers fire on
 * an equation carrying a y' term).
 *
 * Each family is a FORWARD GENERATOR over a genuine parameter grid (Gegenbauer's
 * lambda, Jacobi's (alpha,beta), associated Legendre's order m, the shifted-interval
 * endpoints, Bessel's power), so a fix to one example is proven not to be an overfit
 * by requiring the whole grid to solve.  Verification is a numeric back-substitution
 * on the ORIGINAL equation with every symbolic parameter instantiated at a generic
 * rational (the 2F1 residual is a contiguous-relation identity zero_test cannot
 * discharge, and lie/M14-style methods reach the same 2F1); Head === List is checked
 * first so a decline cannot pass vacuously.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}
static bool lang_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  expected True: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}
#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

/* Solve `eqn`, require Head === List, then require the ORIGINAL-equation residual
 * `resid` to vanish numerically under the substitution list `subs` (parameters +
 * C[1],C[2] + x, chosen inside the mapped interval so any (x-x1)^(1-c) factor is
 * real).  A decline makes Head =!= List and the vacuous pass impossible. */
static void res_ok(const char* eqn, const char* resid, const char* subs) {
    char buf[2200];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Abs[N[(%s) /. sol[[1]] /. {%s}, 20]] < 10^-6]",
        eqn, resid, subs);
    ASSERT_TRUE(buf);
}

/* Family A — Gegenbauer / ultraspherical at symbolic degree n over a lambda grid:
 * (1-x^2) y'' - (2 lam + 1) x y' + n(n + 2 lam) y == 0.  RSPs {-1,1,oo}; the local
 * exponents at +-1 are symmetric.  Kovacic owns integer n; symbolic n is 2F1.  The
 * lambda grid avoids the half-integer values (lambda + 1/2 an integer -> the 2F1
 * lower parameter is an integer -> dependent solutions -> a correct decline to the
 * series fallback, the same limitation the base canonical-Gauss row carries). */
static void t_m17_gegenbauer(void) {
    const char* lam[] = { "1", "2", "3", "3/4", "4/3" };
    for (int i = 0; i < 5; i++) {
        char eqn[320], res[320];
        snprintf(eqn, sizeof(eqn),
            "(1 - x^2) y''[x] - (2 (%s) + 1) x y'[x] + n (n + 2 (%s)) y[x] == 0", lam[i], lam[i]);
        snprintf(res, sizeof(res),
            "(1 - x^2) y''[x] - (2 (%s) + 1) x y'[x] + n (n + 2 (%s)) y[x]", lam[i], lam[i]);
        res_ok(eqn, res, "n -> 23/10, C[1] -> 6/5, C[2] -> 4/5, x -> 3/10");
    }
}

/* Family B — Jacobi at symbolic degree n over an (alpha,beta) grid:
 * (1-x^2) y'' + (beta - alpha - (alpha + beta + 2) x) y' + n(n + alpha + beta + 1) y == 0.
 * RSPs {-1,1,oo} with ASYMMETRIC local exponents (alpha/2 at 1, beta/2 at -1) — a
 * different exponent pair than Gegenbauer, exercising the general F-homotopy.  alpha
 * is kept non-integer (the 2F1 lower parameter is alpha + 1). */
static void t_m17_jacobi(void) {
    const char* ab[][2] = { {"3/2", "1"}, {"1/2", "2"}, {"5/2", "3/2"}, {"3/4", "1/3"} };
    for (int i = 0; i < 4; i++) {
        const char* a = ab[i][0]; const char* b = ab[i][1];
        char eqn[420], res[420];
        snprintf(eqn, sizeof(eqn),
            "(1 - x^2) y''[x] + ((%s) - (%s) - ((%s) + (%s) + 2) x) y'[x] "
            "+ n (n + (%s) + (%s) + 1) y[x] == 0", b, a, a, b, a, b);
        snprintf(res, sizeof(res),
            "(1 - x^2) y''[x] + ((%s) - (%s) - ((%s) + (%s) + 2) x) y'[x] "
            "+ n (n + (%s) + (%s) + 1) y[x]", b, a, a, b, a, b);
        res_ok(eqn, res, "n -> 19/10, C[1] -> 6/5, C[2] -> 4/5, x -> 1/4");
    }
}

/* Family C — associated Legendre at symbolic degree n over a NON-INTEGER order grid
 * mu (matching the corpus, whose order is a symbolic parameter): (1-x^2) y'' - 2 x y'
 * + (n(n+1) - mu^2/(1-x^2)) y == 0.  Declined by the ordinary-Legendre row (mu != 0);
 * the affine/F-homotopy row pulls the +-mu/2 exponents at +-1 to emit verifiable 2F1.
 * Integer mu makes the 2F1 lower parameter 1-mu an integer -> a correct decline. */
static void t_m17_associated_legendre(void) {
    const char* mu[] = { "1/2", "3/2", "5/4", "7/3" };
    for (int i = 0; i < 4; i++) {
        char eqn[320], res[320];
        snprintf(eqn, sizeof(eqn),
            "(1 - x^2) y''[x] - 2 x y'[x] + (n (n + 1) - (%s)^2/(1 - x^2)) y[x] == 0", mu[i]);
        snprintf(res, sizeof(res),
            "(1 - x^2) y''[x] - 2 x y'[x] + (n (n + 1) - (%s)^2/(1 - x^2)) y[x]", mu[i]);
        res_ok(eqn, res, "n -> 23/10, C[1] -> 6/5, C[2] -> 4/5, x -> 3/10");
    }
}

/* Family D — Gauss-class equation on a SHIFTED interval {A, A+L} (RSPs not at {0,1}):
 * (x-A)(A+L-x) y'' + (L c - (a+b+1)(x-A)) y' - a b y == 0, the canonical
 * s(1-s)Y'' + (c-(a+b+1)s)Y' - ab Y == 0 pulled back through s = (x-A)/L.  Fixed
 * non-integer c = 5/4 (a=1/2, b=3/2), verified at the interval midpoint. */
static void t_m17_shifted_gauss(void) {
    /* {A, L, x-midpoint = A + L/2} */
    const char* AA[] = { "2", "1", "0", "-1" };
    const char* LL[] = { "3", "5", "5", "3" };
    const char* XM[] = { "7/2", "7/2", "5/2", "1/2" };
    for (int i = 0; i < 4; i++) {
        char eqn[420], res[420], subs[160];
        snprintf(eqn, sizeof(eqn),
            "(x - (%s))((%s) + (%s) - x) y''[x] + ((%s)(5/4) - 3 (x - (%s))) y'[x] "
            "- (3/4) y[x] == 0", AA[i], AA[i], LL[i], LL[i], AA[i]);
        snprintf(res, sizeof(res),
            "(x - (%s))((%s) + (%s) - x) y''[x] + ((%s)(5/4) - 3 (x - (%s))) y'[x] "
            "- (3/4) y[x]", AA[i], AA[i], LL[i], LL[i], AA[i]);
        snprintf(subs, sizeof(subs), "C[1] -> 11/10, C[2] -> 7/10, x -> %s", XM[i]);
        res_ok(eqn, res, subs);
    }
}

/* Family E — Liouville normal-form pre-pass: y'' + (2k/x) y' + y == 0 has no y'-free
 * form directly, but its normal form is a Bessel potential (spherical Bessel of
 * half-integer order), recovered as y = x^(...) Z(x).  Grid over k. */
static void t_m17_normalform_bessel(void) {
    for (int k = 1; k <= 4; k++) {
        char eqn[200], res[200];
        snprintf(eqn, sizeof(eqn), "y''[x] + (2 %d/x) y'[x] + y[x] == 0", k);
        snprintf(res, sizeof(res), "y''[x] + (2 %d/x) y'[x] + y[x]", k);
        res_ok(eqn, res, "C[1] -> 13/10, C[2] -> 7/10, x -> 6/5");
    }
}

/* Pinned method solves a Gegenbauer; declines a constant-coefficient ODE (not a
 * two-finite-RSP hypergeometric equation). */
static void t_m17_pinned(void) {
    /* Gegenbauer with lambda = 1 (2 lambda + 1 = 3), c = 3/2 non-integer -> 2F1 */
    ASSERT_TRUE("Head[DSolve`SpecialFunctionForm[(1 - x^2) y''[x] - 3 x y'[x] "
                "+ n (n + 2) y[x] == 0, y, x]] === List");
    /* not a 2nd-order linear homogeneous equation -> declines outright */
    ASSERT_TRUE("Head[DSolve`SpecialFunctionForm[y'[x] + y[x] == 0, y, x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m17_gegenbauer);
    TEST(t_m17_jacobi);
    TEST(t_m17_associated_legendre);
    TEST(t_m17_shifted_gauss);
    TEST(t_m17_normalform_bessel);
    TEST(t_m17_pinned);

    printf("All DSolve M17 stress tests passed.\n");
    return 0;
}
