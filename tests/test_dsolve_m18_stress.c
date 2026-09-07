/*
 * test_dsolve_m18_stress.c — anti-overfit stress families for M18
 * (DSolve`ReducibleIntegratingFactor: nonlinear 2nd-order ODEs solved by an
 * integrating factor mu of a restricted form; Cheb-Terrab & Roche 1999).
 *
 * Each family is a FORWARD GENERATOR over a genuine parameter grid, built so the
 * integrating factor and first integral are known by construction, and solved
 * through the PINNED method so ifactor specifically is exercised (not whatever the
 * automatic cascade would pick).  Verification is a numeric back-substitution on
 * the ORIGINAL 2nd-order residual (the solutions carry Airy / Exp-of-polynomial /
 * ArcTanh terms that zero_test cannot discharge); Head === List is checked first so
 * a decline cannot pass vacuously.
 *
 *   Family A (mu(x,y) Case A, mu == 1):  y'' + k y y' == c
 *       -> first integral y' + k y^2/2 - c x == C[1] (Riccati) -> Airy.
 *   Family B (mu(x,y) Case B, mu == 1/y):  y y'' - y'^2 + h(x) y^2 == 0
 *       -> u = Log y satisfies u'' + h == 0 -> y == Exp[C[1] x + C[2] - Int Int h].
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

/* Solve `eqn` via the PINNED reducible-integrating-factor method, require
 * Head === List, then require the ORIGINAL residual `resid` to vanish numerically
 * under the substitution list `subs`.  A decline makes Head =!= List, so the
 * vacuous pass is impossible. */
static void res_ok_if(const char* eqn, const char* resid, const char* subs) {
    char buf[2400];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve`ReducibleIntegratingFactor[%s, y, x]}, "
        "Head[sol] === List && Length[sol] >= 1 && "
        "Abs[N[(%s) /. sol[[1]] /. {%s}, 20]] < 10^-6]",
        eqn, resid, subs);
    ASSERT_TRUE(buf);
}

/* Family A: y'' + k y y' == c  (mu == 1, Case A) -> Airy. */
static void t_m18_stress_caseA_ky(void) {
    const char* ks[] = { "1", "2", "3", "5" };
    const char* cs[] = { "0", "6", "1", "4" };
    for (int i = 0; i < 4; i++) {
        char eqn[256], resid[256];
        snprintf(eqn,  sizeof eqn,  "y''[x] + %s y[x] y'[x] == %s", ks[i], cs[i]);
        snprintf(resid, sizeof resid, "y''[x] + %s y[x] y'[x] - %s", ks[i], cs[i]);
        res_ok_if(eqn, resid, "C[1] -> 6/5, C[2] -> 7/10, x -> 3/5");
    }
}

/* Family B: y y'' - y'^2 + h(x) y^2 == 0  (mu == 1/y, Case B) -> Exp. */
static void t_m18_stress_caseB_logexp(void) {
    /* h(x) with elementary double integral so the closed form exists */
    const char* hs[] = { "2 x", "6 x", "12 x^2", "Cos[x]", "Exp[x]" };
    for (int i = 0; i < 5; i++) {
        char eqn[256], resid[256];
        snprintf(eqn,  sizeof eqn,  "y[x] y''[x] - y'[x]^2 + (%s) y[x]^2 == 0", hs[i]);
        snprintf(resid, sizeof resid, "y[x] y''[x] - y'[x]^2 + (%s) y[x]^2", hs[i]);
        res_ok_if(eqn, resid, "C[1] -> 11/10, C[2] -> 3/5, x -> 7/10");
    }
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m18_stress_caseA_ky);
    TEST(t_m18_stress_caseB_logexp);

    printf("All DSolve M18 stress tests passed.\n");
    return 0;
}
