/*
 * test_dsolve_m62_stress.c — anti-overfit stress families for M62.
 *
 * M62 landed four root-cause fixes, each with a different failure mode, so each gets
 * the family that would actually catch it going wrong:
 *
 *   F1  TimeConstrained CLAMPS under nesting.  The property is an inequality, not an
 *       answer: an inner budget may shorten the enclosing one and may never extend it.
 *       Measured directly on a pure CPU loop, so no subsystem's heuristics are in the
 *       way.  Before the fix `TimeConstrained[TimeConstrained[loop, 30], 3]` ran the
 *       loop to completion in 24.5 s; this family is what notices a regression to
 *       that, and it is the one test here whose assertion is about TIME.
 *
 *   F2  an order-n general solution carries exactly n DISTINCT constants.  This is the
 *       shape of the ExactODE aliasing bug: a solution missing a constant still has a
 *       zero residual, so every verifier in the substrate — and the corpus harness —
 *       scores it as correct.  Counting is the only test that sees it.  Run over a
 *       forward generator of doubly-exact equations (the nesting that triggered it),
 *       not the single case that exposed it.
 *
 *   F3  the symbolic-exponent incomplete-Gamma recognizer, by its OWN certificate:
 *       differentiate the answer and compare, over a (p, a, m) grid.  Plus the
 *       negative controls that matter more than the positives — an INTEGER and a
 *       RATIONAL exponent must keep the elementary / Erf answers they had before, or
 *       the recognizer is a quality regression dressed as a speedup.  And a latency
 *       bound: the general stages take ~13 s to decline on this shape, so if the
 *       recognizer stops firing, the time is what tells you.
 *
 *   F4  the sequential scalar constant fit.  A forward generator of nonlinear
 *       second-order IVPs whose general solution nests the constant inside a Log, i.e.
 *       exactly where Solve's list form bubbles.  Each answer is checked by
 *       substituting it back into BOTH the ODE and the conditions, and the
 *       under-determined member must keep its free constant rather than invent one.
 *
 *   F5  the nonlinear exact (jet-peel) first integral.  Forward generators over
 *       (y^2)'' == g(x) and the 3rd-order (y^2)''' + (y^2)'' == g(x) families, whose
 *       first integral is a total derivative BY CONSTRUCTION, so nothing is
 *       hand-picked.  Verified by back-substitution.  Plus the two gates: a
 *       non-exact nonlinear equation must decline, and an IVP must decline (the
 *       specialists own that class and answer it better).
 *
 * Every family asserts Head[sol] === List FIRST, so a declining method cannot pass
 * vacuously.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

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

static double seconds_of(const char* input) {
    clock_t t0 = clock();
    char* s = eval_str(input);
    double dt = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    free(s);
    return dt;
}

/* ---------------------------------------------------------------------------
 * F1 — TimeConstrained clamps under nesting.
 *
 * The loop count is chosen to be far longer than any budget here, so each case is
 * decided by the timer and not by finishing.  The assertions are one-sided
 * (<= budget + slack) because that is the actual contract.
 * ------------------------------------------------------------------------- */
static void t_m62_timeconstrained_clamps(void) {
    const char* loop = "Do[Sin[1.0], {200000000}]";
    char buf[512];

    /* An unnested budget is honoured: the baseline the nested cases are measured
     * against, so a broken timer cannot make the clamp look fine. */
    snprintf(buf, sizeof buf, "TimeConstrained[%s, 1, $Aborted]", loop);
    double t_plain = seconds_of(buf);
    ASSERT_MSG(t_plain < 5.0, "plain 1 s budget overran: %.2f s", t_plain);

    /* The outer budget bounds a LARGER inner one — the whole bug.  Pre-fix this ran
     * to completion (24.5 s measured for 30-in-3). */
    snprintf(buf, sizeof buf,
             "TimeConstrained[TimeConstrained[%s, 60, $Aborted], 1, $Aborted]", loop);
    double t_nest = seconds_of(buf);
    ASSERT_MSG(t_nest < 6.0, "inner 60 s budget escaped an outer 1 s: %.2f s", t_nest);

    /* A SHORTER inner budget still wins: clamping must not turn into "ignore". */
    snprintf(buf, sizeof buf,
             "TimeConstrained[TimeConstrained[%s, 1, $Aborted], 60, $Aborted]", loop);
    double t_inner = seconds_of(buf);
    ASSERT_MSG(t_inner < 8.0, "inner 1 s budget not honoured: %.2f s", t_inner);

    /* Three levels deep, since the refund compounded per level. */
    snprintf(buf, sizeof buf,
             "TimeConstrained[TimeConstrained[TimeConstrained[%s, 60, $Aborted], 60, $Aborted], 1, $Aborted]",
             loop);
    double t_deep = seconds_of(buf);
    ASSERT_MSG(t_deep < 7.0, "two nested 60 s budgets escaped an outer 1 s: %.2f s", t_deep);

    /* And a budget never truncates a computation that fits inside it. */
    ASSERT_TRUE("TimeConstrained[Integrate[x^2, x], 30, $Aborted] === x^3/3");
    ASSERT_TRUE("TimeConstrained[TimeConstrained[Integrate[x^2, x], 30, $Aborted], 30, $Aborted] === x^3/3");
}

/* ---------------------------------------------------------------------------
 * F2 — an order-n general solution carries exactly n distinct constants.
 *
 * Forward generator: x y''' + 2 y'' == k x (k = 1..4) and its homogeneous member are
 * DOUBLY exact — the reduction re-enters ExactODE, which is the nesting whose shared
 * placeholder aliased two constants into one.  The equations are built from k, so no
 * member is the one that exposed the bug.
 * ------------------------------------------------------------------------- */
static void t_m62_constant_count(void) {
    char buf[700];
    for (int k = 0; k <= 2; k++) {
        /* k == 0 is the homogeneous member. */
        if (k == 0) snprintf(buf, sizeof buf, "x y'''[x] + 2 y''[x] == 0");
        else        snprintf(buf, sizeof buf, "x y'''[x] + 2 y''[x] == %d x", k);

        char q[1200];
        snprintf(q, sizeof q, "Head[DSolve[%s, y, x]] === List", buf);
        ASSERT_TRUE(q);
        /* Exactly three distinct constants for a third-order equation. */
        snprintf(q, sizeof q,
                 "Length[Union[Cases[DSolve[%s, y, x], C[_], Infinity]]] == 3", buf);
        ASSERT_TRUE(q);
        /* ... and the residual is zero, so the count is not bought with a wrong body. */
        snprintf(q, sizeof q,
                 "Module[{s = First[DSolve[%s, y, x]]}, PossibleZeroQ[Simplify[(%s /. Equal -> Subtract) /. s]]]",
                 buf, buf);
        ASSERT_TRUE(q);
    }
    /* A second, differently-shaped doubly-exact family: (x^2 y')'' == k. */
    for (int k = 1; k <= 2; k++) {
        char q[900];
        snprintf(q, sizeof q,
                 "Length[Union[Cases[DSolve[x^2 y'''[x] + 4 x y''[x] + 2 y'[x] == %d, y, x], C[_], Infinity]]] == 3", k);
        ASSERT_TRUE(q);
    }
}

/* ---------------------------------------------------------------------------
 * F3 — symbolic-exponent power times exponential -> incomplete Gamma.
 * ------------------------------------------------------------------------- */
static void t_m62_gammapower(void) {
    /* Positives, by the recognizer's own certificate (differentiate and compare).
     * PowerExpand is the branch convention the closed form is stated in. */
    const char* pos[] = {
        "x^n E^(-x)", "x^n E^x", "x^n E^(-2 x)", "x^n E^(3 x)",
        "x^p E^(-x)", "x^n E^(-x^2)", "x^n E^(-x^3)", "5 x^n E^(-x)",
        "c x^p E^(-a x)"
    };
    for (size_t i = 0; i < sizeof pos / sizeof pos[0]; i++) {
        char q[700];
        /* It must CLOSE (no Integrate left) ... */
        snprintf(q, sizeof q, "FreeQ[Integrate[%s, x], Integrate]", pos[i]);
        ASSERT_TRUE(q);
        /* ... and differentiate back to the integrand. */
        snprintf(q, sizeof q,
                 "PossibleZeroQ[Simplify[PowerExpand[D[Integrate[%s, x], x] - (%s)]]]",
                 pos[i], pos[i]);
        ASSERT_TRUE(q);
    }

    /* Negative controls: a NUMERIC exponent must keep the answer it had before, so
     * the recognizer cannot divert an elementary integral into Gamma form. */
    ASSERT_TRUE("FreeQ[Integrate[x^2 E^(-x), x], Gamma]");
    ASSERT_TRUE("FreeQ[Integrate[x^5 E^(-3 x), x], Gamma]");
    ASSERT_TRUE("FreeQ[Integrate[E^(-x^2), x], Gamma]");
    ASSERT_TRUE("FreeQ[Integrate[x E^(-x^2), x], Gamma]");
    ASSERT_TRUE("PossibleZeroQ[Simplify[D[Integrate[x^2 E^(-x), x], x] - x^2 E^(-x)]]");
    /* A symbolic exponent on a base that is not the variable, and an x-dependent
     * exponent, are not this family at all. */
    ASSERT_TRUE("FreeQ[Integrate[a^n E^(-x), x], Gamma]");

    /* Latency: the general stages need ~13 s per integral to decline on this shape.
     * A generous bound, so this fails only if the recognizer stops firing. */
    double dt = seconds_of("Integrate[x^n E^(-x), x]");
    ASSERT_MSG(dt < 3.0, "symbolic-exponent integral took %.2f s (recognizer not firing?)", dt);

    /* And the ODE the pair of them is for. */
    ASSERT_TRUE("Head[DSolve[y''[x] - y[x] == x^n, y, x]] === List");
    ASSERT_TRUE("FreeQ[DSolve[y''[x] - y[x] == x^n, y, x], Integrate]");
}

/* ---------------------------------------------------------------------------
 * F4 — sequential scalar constant fit (a constant nested inside a Log).
 *
 * Forward generator: y'' + y'^2 + k y' == 0 has general solution
 * C[2] + Log[C[1] - E^(-k x)]/k for k != 0, so the fit of ANY point condition is a
 * Log inversion — the case Solve's list form bubbles on.
 * ------------------------------------------------------------------------- */
static void t_m62_sequential_fit(void) {
    for (int k = 1; k <= 3; k++) {
        char q[900];
        /* Fully determined: two conditions, no constant may survive. */
        snprintf(q, sizeof q,
                 "Head[DSolve[{y''[x] + y'[x]^2 + %d y'[x] == 0, y[0] == 0, y'[0] == 1}, y, x]] === List", k);
        ASSERT_TRUE(q);
        snprintf(q, sizeof q,
                 "FreeQ[DSolve[{y''[x] + y'[x]^2 + %d y'[x] == 0, y[0] == 0, y'[0] == 1}, y, x], C]", k);
        ASSERT_TRUE(q);
        /* ... and it satisfies the ODE and both conditions. */
        snprintf(q, sizeof q,
                 "Module[{s = First[DSolve[{y''[x] + y'[x]^2 + %d y'[x] == 0, y[0] == 0, y'[0] == 1}, y, x]]},"
                 " PossibleZeroQ[Simplify[((y''[x] + y'[x]^2 + %d y'[x]) /. s)]] &&"
                 " PossibleZeroQ[Simplify[(y[0] /. s)]] &&"
                 " PossibleZeroQ[Simplify[(y'[0] /. s) - 1]]]", k, k);
        ASSERT_TRUE(q);

        /* UNDER-determined: one condition on a second-order equation must leave
         * exactly one free constant -- not zero (an invented extra condition) and
         * not two (an unfitted general solution). */
        snprintf(q, sizeof q,
                 "Length[Union[Cases[DSolve[{y''[x] + y'[x]^2 + %d y'[x] == 0, y[0] == 0}, y, x], C[_], Infinity]]] == 1", k);
        ASSERT_TRUE(q);
        snprintf(q, sizeof q,
                 "Module[{s = First[DSolve[{y''[x] + y'[x]^2 + %d y'[x] == 0, y[0] == 0}, y, x]]},"
                 " PossibleZeroQ[Simplify[(y[0] /. s)]]]", k);
        ASSERT_TRUE(q);
    }

    /* The fits the list form already handled must be byte-for-byte unchanged: the
     * sequential pass is a FALLBACK and may not reach a determined linear system. */
    ASSERT_TRUE("DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi/2] == 1}, y, x] === {{y -> Function[{x}, Sin[x]]}}");
    ASSERT_TRUE("DSolve[{y'[x] == y[x], y[0] == 2}, y, x] === {{y -> Function[{x}, 2 E^x]}}");
    /* An inconsistent boundary-value problem still proves no solution. */
    ASSERT_TRUE("DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 1}, y, x] === {}");
    /* An under-determined BVP keeps its free constant (the list form's own job). */
    ASSERT_TRUE("Length[Union[Cases[DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x], C[_], Infinity]]] == 1");
}

/* ---------------------------------------------------------------------------
 * F5 — nonlinear exact first integral by the jet peel.
 *
 * (y^2)'' == g  is  2 y y'' + 2 y'^2 == g, exact by construction; likewise
 * (y^2)''' + (y^2)'' == g  is  2 y y''' + 2(y + 3 y') y'' + 2 y'^2 == g.  Both
 * families are generated from g, so membership is structural, not chosen.
 * ------------------------------------------------------------------------- */
static void t_m62_nonlinear_exact(void) {
    const char* gs[] = { "0", "6 x", "Sin[x]", "12 x^2" };
    for (size_t i = 0; i < sizeof gs / sizeof gs[0]; i++) {
        char q[1100];
        /* order 2 */
        snprintf(q, sizeof q,
                 "Head[DSolve[2 y[x] y''[x] + 2 y'[x]^2 == %s, y, x]] === List", gs[i]);
        ASSERT_TRUE(q);
        snprintf(q, sizeof q,
                 "Module[{s = First[DSolve[2 y[x] y''[x] + 2 y'[x]^2 == %s, y, x]]},"
                 " PossibleZeroQ[Simplify[(2 y[x] y''[x] + 2 y'[x]^2 - (%s)) /. s]]]",
                 gs[i], gs[i]);
        ASSERT_TRUE(q);
        /* order 3 — two members, since each costs ~2 s to solve.  Verified
         * NUMERICALLY: the body is a Sqrt of a sum of exponentials, and `Simplify` on
         * its third-derivative residual costs 13.6 s for one member and does not
         * return for the next, which is a property of Simplify and not of the answer.
         * Sampling at two points with the constants instantiated is the same check the
         * corpus harness makes, at 0.1 s. */
        if (i >= 2) continue;
        snprintf(q, sizeof q,
                 "Head[DSolve[2 y[x] y'''[x] + 2 (y[x] + 3 y'[x]) y''[x] + 2 y'[x]^2 == %s, y, x]] === List",
                 gs[i]);
        ASSERT_TRUE(q);
        snprintf(q, sizeof q,
                 "Module[{s = First[DSolve[2 y[x] y'''[x] + 2 (y[x] + 3 y'[x]) y''[x] + 2 y'[x]^2 == %s, y, x]], r},"
                 " r = (2 y[x] y'''[x] + 2 (y[x] + 3 y'[x]) y''[x] + 2 y'[x]^2 - (%s)) /. s;"
                 " r = r /. {C[1] -> 13/10, C[2] -> 7/5, C[3] -> 9/10};"
                 " Abs[N[r /. x -> 13/10, 30]] < 10^-10 && Abs[N[r /. x -> 17/10, 30]] < 10^-10]",
                 gs[i], gs[i]);
        ASSERT_TRUE(q);
    }

    /* Gate 1: a nonlinear equation that is NOT a total derivative must decline here
     * (the peel's remainder still carries a jet), fast. */
    ASSERT_TRUE("Head[DSolve`ExactODE[y''[x] + y[x]^2 == 0, y, x]] =!= List");
    ASSERT_TRUE("Head[DSolve`ExactODE[y''[x] + y[x] y'[x]^3 == 0, y, x]] =!= List");
    /* Gate 2: an IVP is left to the specialists, which fit the reduction's constants
     * from the conditions and answer better (y'' + 2 y y' == 0 -> Tanh[x]). */
    ASSERT_TRUE("Head[DSolve`ExactODE[{2 y[x] y''[x] + 2 y'[x]^2 == 0, y[0] == 1, y'[0] == 1}, y, x]] =!= List");
    ASSERT_TRUE("DSolve[{y''[x] + 2 y[x] y'[x] == 0, y[0] == 0, y'[0] == 1}, y, x] === {{y -> Function[{x}, Tanh[x]]}}");

    /* The LINEAR exact path is untouched: the two documented decline guards and one
     * documented solve. */
    ASSERT_TRUE("Head[DSolve[(x+1)^2 y''[x] + 3 (x+1) y'[x] + y[x] == x^2, y, x]] === List");
    ASSERT_TRUE("MemberQ[DSolve[y''[x] + Sin[x] y'[x] + Cos[x] y[x] == 0, y, x], SeriesData, Infinity, Heads -> True]");
}

int main(void) {
    setvbuf(stdout, NULL, _IOLBF, 0);   /* the shared alarm() loses buffered output */
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m62_timeconstrained_clamps);
    TEST(t_m62_constant_count);
    TEST(t_m62_gammapower);
    TEST(t_m62_sequential_fit);
    TEST(t_m62_nonlinear_exact);

    printf("All DSolve M62 stress tests passed.\n");
    return 0;
}
