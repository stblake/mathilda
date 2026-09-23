/*
 * test_dsolve_m58_stress.c — anti-overfit stress families for M58
 * (DSolve`AutonomousReduction lifted to any order n >= 2: an autonomous ODE
 * y^(n) == f(y, ..., y^(n-1)) reduces by p = y'(y) through the derivative chain
 * D_{k+1} = p d/dy(D_k) to an order-(n-1) ODE in p(y), then y' == p(y) separable).
 *
 * The main regression risk is the order-2 path (now n=2 of the general chain), so the
 * primary family is a FORWARD GENERATOR over the order-2 power law y y'' == k (y')^2
 * (whose closed form y = C x-power / exponential is guaranteed for every k): the
 * general chain must reproduce it for the whole grid.  The order-3 gain
 * y y''' == y' y'' (reduced p = Sqrt[C1 y^2 + C2], elementary stage-2 quadrature) is
 * the anti-overfit witness that the lift genuinely reaches order 3.  Verification is
 * back-substitution into the ORIGINAL equation (Function form, so derivatives
 * evaluate); a decline makes Head =!= List and a vacuous pass impossible.
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

/* Solve `eqn` (2nd arg y, Function form) and require the residual `resid` back-
 * substitute to a numeric ~0 at several sample points with the generated constants
 * instantiated. */
static void auto_ok(const char* eqn, const char* resid) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Module[{r = (%s) /. sol[[1]] /. {C[1] -> 6/10, C[2] -> 13/10, C[3] -> 7/10}}, "
        "  Max[Table[Abs[N[r /. x -> xv, 20]], {xv, {7/10, 6/5, 17/10}}]] < 10^-5]]",
        eqn, resid);
    ASSERT_TRUE(buf);
}

/* Order-2 power law y y'' == k (y')^2 (n = 2 of the general chain), k = 2..5. */
static void t_m58_power_family(void) {
    for (int k = 2; k <= 5; k++) {
        char eqn[128], res[128];
        snprintf(eqn, sizeof(eqn), "y[x] y''[x] == %d y'[x]^2", k);
        snprintf(res, sizeof(res), "y[x] y''[x] - %d y'[x]^2", k);
        auto_ok(eqn, res);
    }
}

/* Order-2 Tan/Tanh law y'' == a + b (y')^2 through the general chain (a few a,b). */
static void t_m58_tantanh_family(void) {
    int ab[][2] = {{1, 1}, {2, 1}, {1, 3}, {3, 2}};
    for (size_t i = 0; i < 4; i++) {
        char eqn[128], res[128];
        snprintf(eqn, sizeof(eqn), "y''[x] == %d + %d y'[x]^2", ab[i][0], ab[i][1]);
        snprintf(res, sizeof(res), "y''[x] - %d - %d y'[x]^2", ab[i][0], ab[i][1]);
        auto_ok(eqn, res);
    }
}

/* Order-3 gain: y y''' == y' y'' (p = Sqrt[C1 y^2 + C2], elementary quadrature). */
static void t_m58_order3(void) {
    auto_ok("y[x] y'''[x] == y'[x] y''[x]", "y[x] y'''[x] - y'[x] y''[x]");
    /* pinned method reaches it and keeps all three constants */
    ASSERT_TRUE("With[{s = DSolve`AutonomousReduction[y[x] y'''[x] == y'[x] y''[x], y, x]}, "
                "Head[s] === List && Not[FreeQ[s[[1]], C[3]]]]");
}

/* The guard: a 3rd-order autonomous case with a non-elementary stage-2 quadrature
 * (2 y y''' == y', reduced p carries a Log) declines FAST rather than spinning. */
static void t_m58_guard_declines(void) {
    ASSERT_TRUE("Head[DSolve`AutonomousReduction[2 y[x] y'''[x] == y'[x], y[x], x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m58_power_family);
    TEST(t_m58_tantanh_family);
    TEST(t_m58_order3);
    TEST(t_m58_guard_declines);

    printf("All DSolve M58 stress tests passed.\n");
    return 0;
}
