/*
 * test_dsolve_m57_stress.c — anti-overfit stress families for M57
 * (DSolve`SolvableForY / DSolve`SolvableForX: the "dp" differentiation method —
 * solve F(x,y,y')=0 for y (resp. x), differentiate, recurse on the induced ODE,
 * and return the PARAMETRIC general solution {x=X(p,C), y=Y(p)}).
 *
 * Each family is a FORWARD GENERATOR over an integer exponent n: the class is one
 * whose induced first-order ODE closes for every n, so a fix to one example is
 * proven not to be an overfit by requiring the whole grid to solve.  Verification
 * is on the ORIGINAL equation: a parametric branch {x->Function[{t},X],
 * y->Function[{t},Y]} back-substitutes y'[x] = D[Y,t]/D[X,t], and the residual must
 * vanish numerically at several parameter samples (a decline makes Head=!=List, so
 * a vacuous pass is impossible).
 *
 * The flagship SolvableForY family is the 12000.org 2.1.2-347 class
 *   x^(n-1) (y')^n - n x y' + y == 0   (y = n x y' - x^(n-1) (y')^n, NOT affine in
 * x, so Lagrange declines): the induced ODE collapses to the separable
 * dx/dp = -n/(n-1) x/p.  The sign-variant 347+ is an independent family.
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

/* Solve `eqn` and require the ORIGINAL residual `resid` (in x, y[x], y'[x]) vanish
 * numerically on the parametric branch at several parameter samples t0.  `head`
 * is DSolve or a pinned DSolve`Method call spelled with y[x] as the second arg. */
static void sfy_ok(const char* head, const char* eqn, const char* resid) {
    char buf[1800];
    snprintf(buf, sizeof(buf),
        "With[{sol = %s[%s, y[x], x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Module[{xf = x /. sol[[1]], yf = y /. sol[[1]]}, "
        "  Max[Table[Abs[N[(%s) /. {"
        "     y'[x] -> (D[yf[s], s]/D[xf[s], s] /. s -> t0), "
        "     y[x] -> yf[t0], x -> xf[t0]} /. C[1] -> 3/5, 20]], "
        "  {t0, {7/10, 6/5, 19/10, 12/5}}]] < 10^-5]]",
        head, eqn, resid);
    ASSERT_TRUE(buf);
}

/* 347 family:  x^(n-1) y'^n - n x y' + y == 0,  n = 3..6  (auto DSolve). */
static void t_m57_family347(void) {
    for (int n = 3; n <= 6; n++) {
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn),
            "x^%d (y'[x])^%d - %d x y'[x] + y[x] == 0", n - 1, n, n);
        snprintf(res, sizeof(res),
            "x^%d (y'[x])^%d - %d x y'[x] + y[x]", n - 1, n, n);
        sfy_ok("DSolve", eqn, res);
    }
}

/* 347+ sign-variant:  x^(n-1) y'^n + n x y' - y == 0,  n = 3..5  (auto DSolve). */
static void t_m57_family347p(void) {
    for (int n = 3; n <= 5; n++) {
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn),
            "x^%d (y'[x])^%d + %d x y'[x] - y[x] == 0", n - 1, n, n);
        snprintf(res, sizeof(res),
            "x^%d (y'[x])^%d + %d x y'[x] - y[x]", n - 1, n, n);
        sfy_ok("DSolve", eqn, res);
    }
}

/* Pinned DSolve`SolvableForY solves 347/n=3; declines a 2nd-order ODE (max_order
 * gate) and a plainly separable ODE the specialists own. */
static void t_m57_pinned_y(void) {
    sfy_ok("DSolve`SolvableForY",
           "x^2 (y'[x])^3 - 3 x y'[x] + y[x] == 0",
           "x^2 (y'[x])^3 - 3 x y'[x] + y[x]");
    ASSERT_TRUE("Head[DSolve`SolvableForY[y''[x] + y[x] == 0, y[x], x]] =!= List");
}

/* Pinned DSolve`SolvableForX (opt-in mirror) solves the degree-1-in-x form
 * x = y y' + (y')^2  (x - y y' - y'^2 == 0), which is not caught by an earlier
 * specialist; declines a 2nd-order ODE. */
static void t_m57_pinned_x(void) {
    sfy_ok("DSolve`SolvableForX",
           "x - y[x] y'[x] - (y'[x])^2 == 0",
           "x - y[x] y'[x] - (y'[x])^2");
    ASSERT_TRUE("Head[DSolve`SolvableForX[y''[x] + y[x] == 0, y[x], x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m57_family347);
    TEST(t_m57_family347p);
    TEST(t_m57_pinned_y);
    TEST(t_m57_pinned_x);

    printf("All DSolve M57 stress tests passed.\n");
    return 0;
}
