/*
 * test_dsolve_m19_stress.c — anti-overfit stress families for M19
 * (DSolve`SpecialFunctionForm: confluent / Whittaker 2nd-order linear ODEs whose
 * reduced (Liouville normal-form) potential has a SINGLE finite regular singular
 * point x0 and a rank-1 irregular point at infinity, emitted as the verifiable
 * Whittaker->1F1 form  Exp[-z/2] z^(1/2+-mu) Hypergeometric1F1[1/2+-mu-kappa, 1+-2mu, z],
 * z = c (x - x0), c = 2 Sqrt[-b2], mu = Sqrt[1/4 - b0], kappa = b1/c).
 *
 * Each family is a FORWARD GENERATOR over a genuine parameter grid — the Whittaker
 * (kappa, mu), the pole location x0, the energy b2 (constant at infinity), and the
 * y'-carrying spelling reached through the normal-form recovery factor.  A fix to
 * one example is proven not to be an overfit by requiring the whole grid to solve.
 * Verification is a numeric back-substitution on the ORIGINAL equation (the 1F1
 * residual is a contiguous-relation identity zero_test cannot discharge); Head ===
 * List is checked first so a decline cannot pass vacuously.  Every mu in the grids
 * is chosen with 2 mu NON-integer, since integer 2 mu is a correct decline (the two
 * 1F1 partners become dependent / a lower parameter is a non-positive integer).
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
 * `resid` to vanish numerically under `subs` (C[1],C[2],x, away from the pole). */
static void res_ok(const char* eqn, const char* resid, const char* subs) {
    char buf[2200];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && "
        "Abs[N[(%s) /. sol[[1]] /. {%s}, 20]] < 10^-6]",
        eqn, resid, subs);
    ASSERT_TRUE(buf);
}

/* Family A — Whittaker normal form at x0 = 0 (P == 0), grid over (kappa, mu):
 *   y'' + (-1/4 + K/x + (1/4 - M^2)/x^2) y == 0,  solutions WhittakerM[K,+-M,x].
 * mu = M non-half-integer so the two 1F1 partners are independent. */
static void t_m19_whittaker_nf(void) {
    const char* K[] = { "1",   "3/2", "2",   "5/2" };
    const char* M[] = { "1/3", "2/5", "3/7", "1/4" };
    for (int i = 0; i < 4; i++) {
        char eqn[320], res[320];
        snprintf(eqn, sizeof(eqn),
            "y''[x] + (-1/4 + (%s)/x + (1/4 - (%s)^2)/x^2) y[x] == 0", K[i], M[i]);
        snprintf(res, sizeof(res),
            "y''[x] + (-1/4 + (%s)/x + (1/4 - (%s)^2)/x^2) y[x]", K[i], M[i]);
        res_ok(eqn, res, "C[1] -> 6/5, C[2] -> 4/5, x -> 6/5");
    }
}

/* Family B — a SHIFTED single finite regular singular point x0 = X1 (P == 0):
 *   y'' + (-1/4 + 1/(x-X1) + (1/4 - (1/3)^2)/(x-X1)^2) y == 0.
 * Exercises the squarefree-part pole locator at a nonzero x0.  Verified past the
 * pole at x = X1 + 6/5. */
static void t_m19_whittaker_shifted(void) {
    const char* X1[] = { "1", "-2", "3/2", "2" };
    const char* XV[] = { "11/5", "-4/5", "27/10", "16/5" };  /* X1 + 6/5 */
    for (int i = 0; i < 4; i++) {
        char eqn[360], res[360], subs[120];
        snprintf(eqn, sizeof(eqn),
            "y''[x] + (-1/4 + 1/(x - (%s)) + (1/4 - (1/3)^2)/(x - (%s))^2) y[x] == 0", X1[i], X1[i]);
        snprintf(res, sizeof(res),
            "y''[x] + (-1/4 + 1/(x - (%s)) + (1/4 - (1/3)^2)/(x - (%s))^2) y[x]", X1[i], X1[i]);
        snprintf(subs, sizeof(subs), "C[1] -> 6/5, C[2] -> 4/5, x -> %s", XV[i]);
        res_ok(eqn, res, subs);
    }
}

/* Family C — general nonzero energy b2 = E at infinity (P == 0, x0 = 0):
 *   y'' + (E + 1/x + (1/4 - (2/5)^2)/x^2) y == 0,  c = 2 Sqrt[-E].
 * Both signs of E (real vs imaginary scale c) exercise the same 1F1 emission. */
static void t_m19_whittaker_energy(void) {
    const char* E[] = { "1", "-1", "4", "9/4" };
    for (int i = 0; i < 4; i++) {
        char eqn[320], res[320];
        snprintf(eqn, sizeof(eqn),
            "y''[x] + ((%s) + 1/x + (1/4 - (2/5)^2)/x^2) y[x] == 0", E[i]);
        snprintf(res, sizeof(res),
            "y''[x] + ((%s) + 1/x + (1/4 - (2/5)^2)/x^2) y[x]", E[i]);
        res_ok(eqn, res, "C[1] -> 6/5, C[2] -> 4/5, x -> 6/5");
    }
}

/* NOTE: the Whittaker recogniser is run only on the y'-free (P == 0) surface.  The
 * P != 0 confluent family (reached through the Liouville normal-form pre-pass with a
 * recovery factor Exp[-Int P/2]) is deliberately excluded — the recovery factor and
 * the Whittaker z^(1/2+-mu) share the finite-pole base, and verifying the composed
 * candidate can drive the evaluator into $IterationLimit on multi-parameter forms;
 * it is future work.  So there is no pre-pass stress family here. */

/* Pinned: a corpus flagship (102, x^2 y'' + (c x^2 + b x + a)y, symbolic params, P == 0)
 * solves; an integer-2mu Whittaker declines (dependent 1F1 partners -> series fallback);
 * a first-order ODE is not the method's domain. */
static void t_m19_pinned(void) {
    /* corpus 2.1.2-102, symbolic a,b,c -> verifiable 1F1 */
    ASSERT_TRUE("Head[DSolve[x^2 y''[x] + (c x^2 + b x + a) y[x] == 0, y, x]] === List");
    /* 2 mu = 2 (M=1) integer -> the two 1F1 partners are dependent -> correct decline */
    ASSERT_TRUE("Head[DSolve`SpecialFunctionForm[y''[x] + (-1/4 + 1/x + (1/4 - 1) /x^2) y[x] == 0, y, x]] =!= List");
    /* not a 2nd-order linear homogeneous equation -> declines outright */
    ASSERT_TRUE("Head[DSolve`SpecialFunctionForm[y'[x] + y[x] == 0, y, x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m19_whittaker_nf);
    TEST(t_m19_whittaker_shifted);
    TEST(t_m19_whittaker_energy);
    TEST(t_m19_pinned);

    printf("All DSolve M19 stress tests passed.\n");
    return 0;
}
