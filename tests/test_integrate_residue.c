/* test_integrate_residue.c
 *
 * Tests for definite integration by the residue theorem
 * (src/calculus/integrate_residue.c):
 *   Integrate[f, {x, a, b}]                       (auto-dispatch, before Newton-Leibniz)
 *   Integrate[f, {x, a, b}, Method -> "Residue"]  (strict: no NL fallback)
 *   Integrate`ContourResidue[f, {x, a, b}]        (explicit entry point)
 *
 * Coverage: rational integrands on (-Inf,Inf) including higher-order poles
 * (Family A), Fourier/Jordan integrands with Cos/Sin/Exp kernels (Family B),
 * rational-in-{Sin,Cos} integrands over a full period (Family C), the
 * principal-value half-residue (Sin[x]/x = Pi), the even half-line, dispatch
 * selection, and the negative controls that MUST fall through / stay unevaluated
 * (real-axis pole, branch point, partial period).
 *
 * Closed forms that Mathilda leaves in an equivalent-but-unsimplified surface
 * form are pinned numerically via Chop[N[value - reference]] == 0.
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
 * Family A -- rational integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */
static void test_family_rational(void) {
    check_eq("Integrate[1/(1+x^2), {x, -Infinity, Infinity}]", "Pi");
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}]", "Pi/Sqrt[2]");
    check_eq("Integrate[x^2/(1+x^4), {x, -Infinity, Infinity}]", "Pi/Sqrt[2]");
    check_eq("Integrate[1/((x^2+1)(x^2+4)), {x, -Infinity, Infinity}]", "1/6 Pi");
    /* Order-2 pole at I. */
    check_eq("Integrate[1/(1+x^2)^2, {x, -Infinity, Infinity}]", "1/2 Pi");
    /* Value-form independent numeric confirmation. */
    check_eq("Chop[N[Integrate[1/(1+x^6), {x, -Infinity, Infinity}] - 2 Pi/3]]", "0");
}

/* -------------------------------------------------------------------------
 * Family B -- Fourier / Jordan integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */
static void test_family_fourier(void) {
    check_eq("Integrate[Cos[x]/(1+x^2), {x, -Infinity, Infinity}]", "Pi/E");
    check_eq("Integrate[x Sin[x]/(1+x^2), {x, -Infinity, Infinity}]", "Pi/E");
    /* Cos[a x]/(x^2+b^2) = (Pi/b) e^{-a b}. */
    check_eq("Integrate[Cos[2 x]/(x^2+9), {x, -Infinity, Infinity}]", "(1/3 Pi)/E^6");
    /* x Sin[x]/(x^2+4) = Pi/E^2. */
    check_eq("Chop[N[Integrate[x Sin[x]/(x^2+4), {x, -Infinity, Infinity}] - Pi/E^2]]", "0");
}

/* -------------------------------------------------------------------------
 * Family C -- rational-in-{Sin,Cos} over a full period.
 * ---------------------------------------------------------------------- */
static void test_family_trig(void) {
    check_eq("Integrate[1/(2+Cos[x]), {x, 0, 2 Pi}]", "(2 Pi)/Sqrt[3]");
    check_eq("Integrate[1/(5-4 Cos[x]), {x, 0, 2 Pi}]", "2/3 Pi");
    /* (-Pi, Pi) is also a full period. */
    check_eq("Chop[N[Integrate[1/(2+Cos[x]), {x, -Pi, Pi}] - 2 Pi/Sqrt[3]]]", "0");
}

/* -------------------------------------------------------------------------
 * Principal value (simple real-axis pole, half residue) + even half-line.
 * ---------------------------------------------------------------------- */
static void test_principal_value(void) {
    /* Sin[x]/x is analytic at 0 (the kernel supplies a matching zero): the
     * ordinary integral converges to Pi via the half residue. */
    check_eq("Integrate[Sin[x]/x, {x, -Infinity, Infinity}]", "Pi");
    /* Cos[x]/x has a GENUINE pole at 0 (Cos[0] = 1 != 0): the ordinary integral
     * diverges, so the residue method must NOT return a value -- the method
     * leaves it for Newton-Leibniz (strict "Residue" stays unevaluated). */
    check_eq("Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> \"Residue\"]");
}

static void test_half_line(void) {
    check_eq("Integrate[1/(1+x^4), {x, 0, Infinity}]", "(1/2 Pi)/Sqrt[2]");
    check_eq("Chop[N[Integrate[1/(1+x^2), {x, 0, Infinity}] - Pi/2]]", "0");
}

/* -------------------------------------------------------------------------
 * Dispatch: explicit method / builtin; strict Residue has no NL fallback.
 * ---------------------------------------------------------------------- */
static void test_dispatch(void) {
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Pi/Sqrt[2]");
    check_eq("Integrate`ContourResidue[1/(2+Cos[x]), {x, 0, 2 Pi}]", "(2 Pi)/Sqrt[3]");
    /* Strict "Residue" on an integrand no family recognises: unevaluated, NO
     * Newton-Leibniz fallback. */
    check_eq("Integrate[Sin[x], {x, 0, Pi}, Method -> \"Residue\"]",
             "Integrate[Sin[x], {x, 0, Pi}, Method -> \"Residue\"]");
}

/* -------------------------------------------------------------------------
 * Negative controls + regression: ordinary definite integrals unaffected.
 * ---------------------------------------------------------------------- */
static void test_negative_controls(void) {
    /* Real-axis pole of a non-PV integrand: not a clean residue answer. */
    check_eq("Integrate[1/(1+x^3), {x, -Infinity, Infinity}]",
             "Integrate[1/(1 + x^3), {x, -Infinity, Infinity}]");
    /* Branch point (not rational): the residue method must not fire. */
    check_eq("Integrate[1/Sqrt[1+x^4], {x, -Infinity, Infinity}]",
             "Integrate[1/Sqrt[1 + x^4], {x, -Infinity, Infinity}]");
    /* Rational integrand with a real-axis pole: only a principal value exists,
     * which plain Integrate does not compute -- strict "Residue" is unevaluated. */
    check_eq("Integrate[1/(x^2-1), {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[1/(-1 + x^2), {x, -Infinity, Infinity}, Method -> \"Residue\"]");
}

static void test_regression_finite(void) {
    /* Finite definite integrals still go through Newton-Leibniz, unchanged. */
    check_eq("Integrate[1/x, {x, 1, 2}]", "Log[2]");
    check_eq("Integrate[Sin[x], {x, 0, Pi}]", "2");
    check_eq("Integrate[1/(1 + x^2), {x, 0, 1}]", "1/4 Pi");
    check_eq("Integrate[x^2, {x, 0, 1}]", "1/3");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_family_rational);
    TEST(test_family_fourier);
    TEST(test_family_trig);
    TEST(test_principal_value);
    TEST(test_half_line);
    TEST(test_dispatch);
    TEST(test_negative_controls);
    TEST(test_regression_finite);

    printf("All Integrate ContourResidue tests passed!\n");
    return 0;
}
