/* test_integrate_line.c
 *
 * Tests for complex line / contour integration (src/calculus/integrate_line.c):
 *   Integrate[f, {x, z0, z1, ..., zn}]      (auto-dispatch, contour)
 *   Integrate`LineIntegral[f, {x, z0, ...}] (explicit entry point)
 *   Integrate`PathSingularPoints[f, {x, ...}]
 *
 * Coverage: on-path singularity detection, single-segment values, complex-plane
 * endpoint limits, divergence on a contour through a pole, branch-correct
 * closed-contour residues, polyline path independence, and confirmation that the
 * real-axis definite-integral path (Newton-Leibniz) is unaffected.
 *
 * Complex-log results that Mathilda leaves in an un-combined form (it does not
 * merge Log[I] - Log[-I] into 2 I Pi) are pinned numerically via
 * Chop[N[value - reference]] == 0, which verifies the *value* independent of the
 * surface form.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------------------------------------------------------------------------
 * On-path singularity detection.
 * ---------------------------------------------------------------------- */
static void test_detection(void) {
    /* Pole of 1/z is at 0, off the segment from 1-I to 2+3I (Re >= 1). */
    check_eq("Integrate`PathSingularPoints[1/z, {z, 1 - I, 2 + 3 I}]", "{}");
    /* The segment -1-I -> 1+I passes through the origin at its midpoint. */
    check_eq("Integrate`PathSingularPoints[1/z, {z, -1 - I, 1 + I}]", "{0}");
    /* Vertical segment 0 -> 2I hits the pole of 1/(z-I) at I. */
    check_eq("Integrate`PathSingularPoints[1/(z - I), {z, 0, 2 I}]", "{I}");
    /* Both poles +/-I lie on the imaginary-axis segment -2I -> 2I. */
    check_eq("Integrate`PathSingularPoints[1/((z - I)(z + I)), {z, -2 I, 2 I}]",
             "{-I, I}");
    /* Poles +/-I are off the real-axis segment. */
    check_eq("Integrate`PathSingularPoints[1/(z^2 + 1), {z, -2, 2}]", "{}");
    /* Regression: the horizontal segment I-1 -> I+1 passes through the pole at
     * I.  The parametrised denominator is quadratic in t with COMPLEX
     * coefficients; Solve[..., Reals] wrongly drops its real root t=1/2, so the
     * detector must solve over C and keep the real root -- otherwise the
     * divergence goes unreported. */
    check_eq("Integrate`PathSingularPoints[1/(z^2 + 1), {z, I - 1, I + 1}]",
             "{I}");
}

/* -------------------------------------------------------------------------
 * Single-segment values (segments that do not cross a branch cut).
 * ---------------------------------------------------------------------- */
static void test_segment_values(void) {
    /* Marquee example: 1/x from 1-I to 2+3I. */
    check_eq("Integrate[1/x, {x, 1 - I, 2 + 3 I}]", "Log[2 + 3*I] - Log[1 - I]");
    /* ... and its value, pinned numerically. */
    check_eq("Chop[N[Integrate[1/x, {x, 1 - I, 2 + 3 I}] "
             "- (Log[2 + 3 I] - Log[1 - I])]]", "0");

    /* Entire integrands: value depends only on the endpoints. */
    check_eq("Integrate[z, {z, 0, 1 + I}]",   "I");            /* (1+I)^2/2 */
    check_eq("Integrate[z^2, {z, 0, 1 + I}]", "-2/3 + 2/3*I"); /* (1+I)^3/3 */
    check_eq("Integrate[Cos[z], {z, 0, I}]",  "I Sinh[1]");    /* Sin[I] */
    check_eq("Simplify[Integrate[Exp[z], {z, 0, I Pi}]]", "-2"); /* e^{iPi}-1 */

    /* Rational integrand whose antiderivative is an inverse tangent. */
    check_eq("Integrate[1/(z^2 + 1), {z, 0, 1 + I}]", "ArcTan[1 + I]");
}

/* -------------------------------------------------------------------------
 * Complex-plane limits: endpoint singularities.
 * ---------------------------------------------------------------------- */
static void test_complex_limits(void) {
    /* 1/Sqrt[z] has an integrable branch-point at the endpoint 0: the one-sided
     * limit F(gamma(t)) -> 0 as t -> 0 is finite. */
    check_eq("Integrate[1/Sqrt[z], {z, 0, 1 + I}]", "2 Sqrt[1 + I]");

    /* Segment I -> 1 avoids the pole of 1/z at 0; value = Log[1] - Log[I]. */
    check_eq("Chop[N[Integrate[1/z, {z, I, 1}] + I Pi/2]]", "0");
}

/* -------------------------------------------------------------------------
 * Divergence: a genuine singularity on the contour.
 * ---------------------------------------------------------------------- */
static void test_divergence(void) {
    /* Segment through the pole at 0: contour integral is divergent -> the input
     * is returned unevaluated. */
    check_eq("Integrate[1/z, {z, -1 - I, 1 + I}]",
             "Integrate[1/z, {z, -1 - I, 1 + I}]");
    /* Pole at the endpoint 0: also divergent. */
    check_eq("Integrate[1/z, {z, I, 0}]", "Integrate[1/z, {z, I, 0}]");
    /* Segment through the pole at I of a QUADRATIC denominator (complex-
     * coefficient parametrisation): divergent -> unevaluated, not a spurious
     * finite value. */
    check_eq("Integrate[1/(z^2 + 1), {z, I - 1, I + 1}]",
             "Integrate[1/(1 + z^2), {z, -1 + I, 1 + I}]");
}

/* -------------------------------------------------------------------------
 * Polyline contours + branch-correct closed-contour residues.
 * ---------------------------------------------------------------------- */
static void test_contours(void) {
    /* Polyline through 0 -> 1 -> 1+I equals the direct segment 0 -> 1+I for an
     * entire integrand (path independence). */
    check_eq("Integrate[z^2, {z, 0, 1, 1 + I}]", "-2/3 + 2/3*I");

    /* Closed contour of an entire function integrates to 0. */
    check_eq("Integrate[z, {z, 0, 1, 1 + I, I, 0}]", "0");

    /* Winding integral of 1/z once counter-clockwise around the origin (square
     * 1 -> I -> -1 -> -I -> 1) is 2 Pi I -- recovered branch-correctly across
     * the negative-real-axis branch cut of Log. */
    check_eq("Chop[N[Integrate[1/z, {z, 1, I, -1, -I, 1}] - 2 Pi I]]", "0");
}

/* -------------------------------------------------------------------------
 * Dispatch: explicit builtin, and the real-axis path is untouched.
 * ---------------------------------------------------------------------- */
static void test_dispatch(void) {
    /* Explicit entry point matches the auto-dispatched result. */
    check_eq("Integrate`LineIntegral[1/x, {x, 1 - I, 2 + 3 I}]",
             "Log[2 + 3*I] - Log[1 - I]");
    /* Explicit builtin with real endpoints degenerates to the ordinary integral. */
    check_eq("Integrate`LineIntegral[z^2, {z, 0, 2}]", "8/3");

    /* Real definite integrals still go through Newton-Leibniz, unchanged. */
    check_eq("Integrate[1/x, {x, 1, 2}]", "Log[2]");
    check_eq("Integrate[Sin[x], {x, 0, Pi}]", "2");
    check_eq("Integrate[1/(1 + x^2), {x, 0, 1}]", "1/4 Pi");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_detection);
    TEST(test_segment_values);
    TEST(test_complex_limits);
    TEST(test_divergence);
    TEST(test_contours);
    TEST(test_dispatch);

    printf("All Integrate LineIntegral tests passed!\n");
    return 0;
}
