/*
 * test_dsolve_m14_stress.c — anti-overfit stress families for M14
 * (DSolve`ChangeOfVariable: 2nd-order linear ODEs with transcendental coefficients
 * solved by a change of the independent variable t = phi(x) that rationalizes them).
 *
 * Each family is a FORWARD GENERATOR over an integer degree: the ODE class is one
 * that becomes rational (here the Legendre equation) under t = Cos[x], so a fix to
 * one example is proven not to be an overfit by requiring the whole grid to solve.
 * Verification is on the ORIGINAL transcendental-coefficient equation: the general
 * solution is a 2-parameter family, so we require the residual's C[1] and C[2]
 * coefficients to both vanish numerically (a decline would make Head=!=List and a
 * vacuous pass impossible).  The run stays well under the harness alarm.
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

/* Solve `eqn` and require both C[1],C[2] residual coefficients vanish numerically. */
static void cv_ok(const char* eqn, const char* resid) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Module[{r = Simplify[(%s) /. sol[[1]]]}, "
        "Abs[N[Coefficient[r, C[1]] /. x->13/10, 12]] < 10^-5 && "
        "Abs[N[Coefficient[r, C[2]] /. x->13/10, 12]] < 10^-5]]",
        eqn, resid);
    ASSERT_TRUE(buf);
}

/* Legendre via t = Cos[x]:  y'' + Cot[x] y' + k(k+1) y == 0. */
static void t_m14_legendre(void) {
    for (int k = 1; k <= 5; k++) {
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn), "y''[x] + Cot[x] y'[x] + %d*%d y[x] == 0", k, k+1);
        snprintf(res, sizeof(res), "y''[x] + Cot[x] y'[x] + %d*%d y[x]", k, k+1);
        cv_ok(eqn, res);
    }
}

/* The Sturm-Liouville spelling y'' Sin[x] + y' Cos[x] + m(m+1) y Sin[x] == 0. */
static void t_m14_sinform(void) {
    for (int m = 1; m <= 4; m++) {
        char eqn[300], res[300];
        snprintf(eqn, sizeof(eqn),
            "y''[x] Sin[x] + y'[x] Cos[x] + %d*%d y[x] Sin[x] == 0", m, m+1);
        snprintf(res, sizeof(res),
            "y''[x] Sin[x] + y'[x] Cos[x] + %d*%d y[x] Sin[x]", m, m+1);
        cv_ok(eqn, res);
    }
}

/* Pinned method solves; declines an already-rational ODE (not its domain). */
static void t_m14_pinned(void) {
    ASSERT_TRUE("Head[DSolve`ChangeOfVariable[y''[x] + Cot[x] y'[x] + 6 y[x] == 0, y, x]] === List");
    ASSERT_TRUE("Head[DSolve`ChangeOfVariable[y''[x] + y[x] == 0, y, x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m14_legendre);
    TEST(t_m14_sinform);
    TEST(t_m14_pinned);

    printf("All DSolve M14 stress tests passed.\n");
    return 0;
}
