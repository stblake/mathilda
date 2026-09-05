/*
 * test_dsolve_m12_stress.c — anti-overfit stress families for M12
 * (DSolve`SecondOrderSymmetry, the nonlinear-second-order Lie point-symmetry
 * method).
 *
 * Each family is a FORWARD GENERATOR over a coefficient grid: the ODE class is
 * one lie2 solves via a genuine symmetry (projective / scaling), so a fix to one
 * seed example is proven not to be an overfit by requiring the WHOLE grid to solve.
 * Verification is NUMERIC, not PossibleZeroQ: these solutions carry logs / radicals
 * whose residual PossibleZeroQ cannot decide (and lie2 itself relies on a numeric
 * back-substitution guard, so the symbolic residual is deliberately left
 * undecidable).  Each case asserts Head[DSolve[...]]===List first (a decline would
 * make a vacuous pass) and then |residual| < 1e-6 at two independent (C[1],C[2],x)
 * points.  The whole run stays well under the 60 s harness alarm.
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

/* Solve `eqn` automatically and require: a List result whose first branch makes
 * `resid` (= lhs - rhs) numerically vanish at two independent sample points. */
static void nl2_ok(const char* eqn, const char* resid) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, "
        "Head[sol] === List && Length[sol] >= 1 && "
        "Abs[N[(%s) /. sol[[1]] /. {C[1]->13/11, C[2]->5/7} /. x->17/10, 14]] < 10^-6 && "
        "Abs[N[(%s) /. sol[[1]] /. {C[1]->9/4, C[2]->3/8} /. x->23/10, 14]] < 10^-6]",
        eqn, resid, resid);
    ASSERT_TRUE(buf);
}

/* Family B — projective symmetry [x,y],[x^2,x y]:  x^3 y'' == a (y - x y')^2.
 * (Cheb-Terrab et al. Eq.3 / Kamke-181 class; the arbitrary-a generalization.) */
static void t_m12_projective(void) {
    const char* as[] = { "1", "2", "3", "4", "5", "6", "-1", "-2", "-3", "-4" };
    for (size_t i = 0; i < sizeof(as)/sizeof(as[0]); i++) {
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn), "x^3 y''[x] == (%s)(y[x] - x y'[x])^2", as[i]);
        snprintf(res, sizeof(res), "x^3 y''[x] - (%s)(y[x] - x y'[x])^2", as[i]);
        nl2_ok(eqn, res);
    }
}

/* Family C — projective with a linear term (same symmetry algebra):
 *   x^3 y'' == (y - x y')^2 + d x (y - x y'). */
static void t_m12_projective_linear(void) {
    const char* ds[] = { "0", "1", "2", "3", "-1", "-2" };
    for (size_t i = 0; i < sizeof(ds)/sizeof(ds[0]); i++) {
        char eqn[320], res[320];
        snprintf(eqn, sizeof(eqn),
            "x^3 y''[x] == (y[x] - x y'[x])^2 + (%s) x (y[x] - x y'[x])", ds[i]);
        snprintf(res, sizeof(res),
            "x^3 y''[x] - ((y[x] - x y'[x])^2 + (%s) x (y[x] - x y'[x]))", ds[i]);
        nl2_ok(eqn, res);
    }
}

/* Family D — separable-scaling symmetry [x,0],[0,y]:  k x^2 y y'' + y^2 == x^2 y'^2
 * (the u = y^2 -> Euler class; N5 is k=2). */
static void t_m12_scaling(void) {
    const char* ks[] = { "2", "3", "4", "5", "-1", "-2" };
    for (size_t i = 0; i < sizeof(ks)/sizeof(ks[0]); i++) {
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn),
            "(%s) x^2 y''[x] y[x] + y[x]^2 == x^2 y'[x]^2", ks[i]);
        snprintf(res, sizeof(res),
            "(%s) x^2 y''[x] y[x] + y[x]^2 - x^2 y'[x]^2", ks[i]);
        nl2_ok(eqn, res);
    }
}

/* The named seed deficiencies from the 12000.org corpus (Kamke / Murphy). */
static void t_m12_seeds(void) {
    nl2_ok("x^3 y''[x] == (-x y'[x] + y[x])^2",
           "x^3 y''[x] - (-x y'[x] + y[x])^2");                          /* N3 */
    nl2_ok("x^2 y''[x] y[x] + x^2 y'[x]^2 - 5 x y[x] y'[x] == 4 y[x]^2",
           "x^2 y''[x] y[x] + x^2 y'[x]^2 - 5 x y[x] y'[x] - 4 y[x]^2"); /* N4 */
    nl2_ok("2 x^2 y''[x] y[x] + y[x]^2 == x^2 y'[x]^2",
           "2 x^2 y''[x] y[x] + y[x]^2 - x^2 y'[x]^2");                  /* N5 */
    nl2_ok("x^2 (x + y[x]) y''[x] == (x y'[x] - y[x])^2",
           "x^2 (x + y[x]) y''[x] - (x y'[x] - y[x])^2");                /* Proj */
    nl2_ok("y''[x] == (y'[x] x - y[x])^2 / x^3",
           "y''[x] - (y'[x] x - y[x])^2 / x^3");                        /* Eq3 */
}

/* The pinned method really is the one solving these (not a neighbour). */
static void t_m12_pinned(void) {
    ASSERT_TRUE("Head[DSolve`SecondOrderSymmetry["
                "x^3 y''[x] == (y[x] - x y'[x])^2, y, x]] === List");
    ASSERT_TRUE("Abs[N[(x^3 y''[x] - (y[x]-x y'[x])^2) /. "
                "DSolve`SecondOrderSymmetry[x^3 y''[x] == (y[x]-x y'[x])^2, y, x][[1]] "
                "/. {C[1]->7/5,C[2]->3/4} /. x->13/10, 14]] < 10^-6");
    /* a linear ODE is declined by lie2 (its domain is nonlinear equations) */
    ASSERT_TRUE("Head[DSolve`SecondOrderSymmetry[y''[x] + y[x] == 0, y, x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m12_seeds);
    TEST(t_m12_projective);
    TEST(t_m12_projective_linear);
    TEST(t_m12_scaling);
    TEST(t_m12_pinned);

    printf("All DSolve M12 stress tests passed.\n");
    return 0;
}
