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

/* Solve `eqn` and require the residual vanish numerically after instantiating a
 * SYMBOLIC free parameter (regression guard for the cv_num_ok param-instantiation
 * fix: a solution carrying a symbolic degree must still verify). */
static void cv_ok_param(const char* eqn, const char* resid,
                        const char* param, const char* pval) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Module[{r = (%s) /. sol[[1]] /. {%s -> %s, C[1] -> 13/10, C[2] -> 7/10, x -> 11/10}}, "
        "Abs[N[r, 20]] < 10^-6]]",
        eqn, resid, param, pval);
    ASSERT_TRUE(buf);
}

/* Symbolic-degree Legendre via t = Cos[x]: y'' + Cot[x] y' + nu(nu+1) y == 0 with
 * nu SYMBOLIC (solution carries LegendreP[nu, Cos[x]] / LegendreQ[nu, Cos[x]]).
 * Checked at several non-integer nu — before the cv_num_ok fix these all declined
 * because the symbolic nu made the numeric verifier reject the correct transform. */
static void t_m16_legendre_symbolic(void) {
    const char* eqn = "y''[x] + Cot[x] y'[x] + nu (nu + 1) y[x] == 0";
    const char* res = "y''[x] + Cot[x] y'[x] + nu (nu + 1) y[x]";
    const char* nuv[] = { "7/3", "3/2", "12/5", "5/4" };
    for (int i = 0; i < 4; i++) cv_ok_param(eqn, res, "nu", nuv[i]);
}

/* Solve `eqn` (symbolic-parameter trig-potential family) and require the residual
 * vanish at two small generic parameter instantiations where the emitted 2F1
 * numericizes cleanly (a,b,c,p,q < 1 keep the hypergeometric argument in range). */
static void pt_verify(const char* eqn, const char* res) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Module[{r1 = (%s) /. sol[[1]] /. {a->31/100,b->22/100,c->17/100,p->42/100,q->53/100,C[1]->13/10,C[2]->7/10,x->3/5}, "
        "        r2 = (%s) /. sol[[1]] /. {a->-27/100,b->35/100,c->12/100,p->48/100,q->39/100,C[1]->9/10,C[2]->6/5,x->4/5}}, "
        "  Abs[N[r1,18]] < 10^-6 && Abs[N[r2,18]] < 10^-6]]",
        eqn, res, res);
    ASSERT_TRUE(buf);
}

/* M16 Pöschl-Teller / trigonometric-potential recognizer (SpecialFunctionForm):
 * y'' == (a + p(p-1)Csc^2 x + q(q-1)Sec^2 x) y and its Csc-only and
 * (a Cos^2+b Sin^2+c)/Sin^2 spellings -> Hypergeometric2F1 (numerically verified). */
static void t_m16_poschl_teller(void) {
    pt_verify("y''[x] == (a + p (p-1) Csc[x]^2 + q (q-1) Sec[x]^2) y[x]",
              "y''[x] - (a + p (p-1) Csc[x]^2 + q (q-1) Sec[x]^2) y[x]");
    pt_verify("y''[x] == (a + p (p-1) Csc[x]^2) y[x]",
              "y''[x] - (a + p (p-1) Csc[x]^2) y[x]");
    pt_verify("y''[x] == ((a Cos[x]^2 + b Sin[x]^2 + c) y[x])/Sin[x]^2",
              "y''[x] - ((a Cos[x]^2 + b Sin[x]^2 + c) y[x])/Sin[x]^2");
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
    TEST(t_m16_legendre_symbolic);
    TEST(t_m16_poschl_teller);
    TEST(t_m14_pinned);

    printf("All DSolve M14 stress tests passed.\n");
    return 0;
}
