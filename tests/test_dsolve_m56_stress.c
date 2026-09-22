/*
 * test_dsolve_m56_stress.c — anti-overfit stress families for M56
 * (DSolve`ReducibleFirstIntegral: mu(x,y') integrating factors, Cheb-Terrab &
 * Roche 1999 Section 2.2, Lemma 3 Cases A/C/D; reduction-of-order first integral
 * R(x, y[x], y'[x]) == C[1]).
 *
 * Each family is a FORWARD GENERATOR over a genuine parameter grid built to stay
 * inside one Lemma-3 case, solved through the AUTOMATIC cascade (so the whole
 * DSolve dispatch is exercised — the reduction-of-order emit sits after the
 * full-solution methods).  Verification is INTRINSIC to the first integral: the
 * returned branch must be a single {R == C[1]} that genuinely contains y'[x] (a
 * real reduction of order, not a vacuous decline), and the total x-derivative
 * D[R, x], with y''[x] replaced by the ODE's Phi, must vanish numerically — i.e.
 * R is constant along every solution.  Head === List is checked first so a decline
 * cannot pass vacuously.
 *
 *   Family A (Case A, 2.36-2.40):  y'' == (a x^2 y y' + a x y^2)/y'   (Kamke 226)
 *       -> mu == y',           first integral y'^2/2 - a x^2 y^2/2 == C[1].
 *   Family C (Case C, 2.41-2.50): y'' == (1 + y'^2)/(x - y + k)      (Kamke 136)
 *       -> mu == (y'-1)/(1+y'^2), first integral Log[x-y+k] + Log[1+y'^2]/2 - ArcTan[y'].
 *   Family D (Case D, 2.57-2.64): y'' == a (c + b x + y)(1 + y'^2)^(3/2)  (Kamke 66)
 *       -> mu == (b+y')/(a(1+y'^2)^(3/2)).
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

/* Solve `eqn` (an ODE y''[x] == Phi with RHS `phi`) through the automatic cascade,
 * require the result to be a first integral {{ R == C[1] }} that contains y'[x],
 * and require D[R,x] (with y''[x] -> Phi) to vanish numerically under `subs`. */
static void fi_ok(const char* eqn, const char* phi, const char* subs) {
    char buf[3200];
    snprintf(buf, sizeof buf,
        "With[{sol = DSolve[%s, y, x]}, "
        "Head[sol] === List && Length[sol] >= 1 && "
        "MatchQ[sol[[1, 1]], _Equal] && "
        "Not[FreeQ[sol[[1, 1, 1]], Derivative[1][y][x]]] && "
        "Abs[N[((D[sol[[1, 1, 1]], x] /. Derivative[2][y][x] -> (%s)) /. {%s}), 20]] < 10^-6]",
        eqn, phi, subs);
    ASSERT_TRUE(buf);
}

/* Family A (Case A): y'' == (a x^2 y y' + a x y^2)/y' -> mu == y'. */
static void t_m56_stress_caseA(void) {
    const char* as[] = { "1", "2", "3", "5" };
    for (int i = 0; i < 4; i++) {
        char eqn[256], phi[256];
        snprintf(phi, sizeof phi, "(%s x^2 y[x] y'[x] + %s x y[x]^2)/y'[x]", as[i], as[i]);
        snprintf(eqn, sizeof eqn, "y''[x] == %s", phi);
        fi_ok(eqn, phi, "y[x] -> 9/10, Derivative[1][y][x] -> 7/5, x -> 6/5");
    }
}

/* Family C (Case C): y'' == (1 + y'^2)/(x - y + k) -> mu == (y'-1)/(1+y'^2). */
static void t_m56_stress_caseC(void) {
    const char* ks[] = { "0", "1", "2", "-1" };
    for (int i = 0; i < 4; i++) {
        char eqn[256], phi[256];
        snprintf(phi, sizeof phi, "(1 + y'[x]^2)/(x - y[x] + (%s))", ks[i]);
        snprintf(eqn, sizeof eqn, "y''[x] == %s", phi);
        /* choose a sample with x - y + k != 0 */
        fi_ok(eqn, phi, "y[x] -> 3/10, Derivative[1][y][x] -> 6/5, x -> 8/5");
    }
}

/* Family D (Case D): y'' == a (c + b x + y)(1 + y'^2)^(3/2) -> mu == (b+y')/(a(1+y'^2)^(3/2)). */
static void t_m56_stress_caseD(void) {
    const char* as[] = { "1", "2", "1", "3" };
    const char* bs[] = { "1", "1", "2", "1" };
    const char* cs[] = { "0", "1", "1", "2" };
    for (int i = 0; i < 4; i++) {
        char eqn[320], phi[256], subs[256];
        snprintf(phi, sizeof phi, "%s (%s + %s x + y[x]) (1 + y'[x]^2)^(3/2)", as[i], cs[i], bs[i]);
        snprintf(eqn, sizeof eqn, "y''[x] == %s", phi);
        snprintf(subs, sizeof subs, "y[x] -> 9/10, Derivative[1][y][x] -> 7/5, x -> 6/5");
        fi_ok(eqn, phi, subs);
    }
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m56_stress_caseA);
    TEST(t_m56_stress_caseC);
    TEST(t_m56_stress_caseD);

    printf("All DSolve M56 stress tests passed.\n");
    return 0;
}
