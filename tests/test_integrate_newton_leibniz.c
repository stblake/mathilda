/* test_integrate_newton_leibniz.c
 *
 * Unit + stress tests for definite integration by the Newton-Leibniz rule
 * (the fundamental theorem of calculus):
 *
 *   Integrate[f, {x, a, b}]                       (automatic dispatch)
 *   Integrate[f, {x, a, b}, Method -> "NewtonLeibniz"]
 *   Integrate`NewtonLeibniz[f, {x, a, b}]         (explicit entry point)
 *   Integrate`SingularPoints[expr, {x, a, b}]     (the pole detector)
 *
 * Correctness is asserted three ways:
 *   1. exact closed-form results for textbook integrals,
 *   2. the detector's interior-pole list on hand-checked rational inputs,
 *   3. a numeric differential cross-check against NIntegrate for a stress
 *      batch (the symbolic definite value must agree with the numerical one).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Evaluate `input`, assert its printed form equals `expected`. */
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

/* Assert Integrate[f, {x, a, b}] agrees numerically with NIntegrate. */
static void num_match(const char* f, const char* a, const char* b) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "N[Abs[Integrate[%s, {x, %s, %s}] - "
             "NIntegrate[%s, {x, %s, %s}]]] < 1/100000000",
             f, a, b, f, a, b);
    check_eq(buf, "True");
}

/* ------------------------------------------------------------------ */
/* 1. Pole detector: interior real roots of Denominator[Together[.]].   */
/* ------------------------------------------------------------------ */
static void test_detector(void) {
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 0, 3}]",  "{1}");
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 2, 3}]",  "{}");
    check_eq("Integrate`SingularPoints[1/((x - 1)(x - 2)), {x, 0, 3}]", "{1, 2}");
    check_eq("Integrate`SingularPoints[1/(x^2 + 1), {x, -5, 5}]", "{}");
    /* Pole exactly on an endpoint is NOT interior. */
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 1, 3}]",  "{}");
    /* No x in the denominator -> no poles. */
    check_eq("Integrate`SingularPoints[x^2 + 3 x, {x, -5, 5}]", "{}");
    /* Sorted ascending, several interior. */
    check_eq("Integrate`SingularPoints[1/((x + 1) x (x - 1)), {x, -2, 2}]",
             "{-1, 0, 1}");
}

/* ------------------------------------------------------------------ */
/* 2. Core FTC on proper integrals (continuous antiderivative).         */
/* ------------------------------------------------------------------ */
static void test_proper(void) {
    check_eq("Integrate[x^2, {x, 0, 1}]",        "1/3");
    check_eq("Integrate[Sin[x], {x, 0, Pi}]",    "2");
    check_eq("Integrate[1/(1 + x^2), {x, 0, 1}]", "1/4 Pi");
    check_eq("Integrate[E^x, {x, 0, 1}]",        "-1 + E");
    check_eq("Integrate[Cos[x], {x, 0, Pi/2}]",  "1");
    /* Pole of the integrand strictly OUTSIDE the interval: proper. */
    check_eq("Integrate[1/x^3, {x, 1, 2}]",      "3/8");
    check_eq("Integrate[1/x^2, {x, -2, -1}]",    "1/2");
    check_eq("Integrate[1/x, {x, 1, 2}]",        "Log[2]");
    check_eq("Integrate[1/x, {x, 1, E}]",        "1");
}

/* ------------------------------------------------------------------ */
/* 3. Improper integrals: one-sided limits, infinite bounds.            */
/* ------------------------------------------------------------------ */
static void test_improper(void) {
    /* Integrable endpoint singularity of the integrand. */
    check_eq("Integrate[1/Sqrt[x], {x, 0, 1}]",         "2");
    /* Infinite upper bound. */
    check_eq("Integrate[1/(1 + x^2), {x, 0, Infinity}]", "1/2 Pi");
    check_eq("Integrate[E^(-x), {x, 0, Infinity}]",      "1");
    check_eq("Integrate[1/x^2, {x, 1, Infinity}]",       "1");
}

/* ------------------------------------------------------------------ */
/* 4. Divergent integrals: interior pole -> a divergence sentinel.      */
/* ------------------------------------------------------------------ */
static void test_divergent(void) {
    /* Even-order interior pole -> +Infinity. */
    check_eq("Integrate[1/x^2, {x, -1, 1}]",     "Infinity");
    /* Odd-order interior pole -> Indeterminate (NOT a spurious finite). */
    check_eq("Integrate[1/x, {x, -1, 1}]",       "Indeterminate");
}

/* ------------------------------------------------------------------ */
/* 5. Iterated / multiple integrals (inner = last spec, done first).    */
/* ------------------------------------------------------------------ */
static void test_iterated(void) {
    check_eq("Integrate[x y, {x, 0, 1}, {y, 0, 1}]", "1/4");
    check_eq("Integrate[x, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]", "1/2");
    /* Inner bound depends on the outer variable. */
    check_eq("Integrate[1, {x, 0, 1}, {y, 0, x}]",   "1/2");
}

/* ------------------------------------------------------------------ */
/* 6. Explicit entry point, Method selection, Method pass-through.      */
/* ------------------------------------------------------------------ */
static void test_entry_points(void) {
    check_eq("Integrate`NewtonLeibniz[x^2, {x, 0, 1}]",              "1/3");
    check_eq("Integrate[x^2, {x, 0, 1}, Method -> \"NewtonLeibniz\"]", "1/3");
    /* Method pass-through to the inner indefinite integration. */
    check_eq("Integrate[x^3, {x, 0, 2}, Method -> \"BronsteinRational\"]", "4");
    /* Symbolic bounds with a pole-free antiderivative. */
    check_eq("Integrate[2 x, {x, a, b}]", "-a^2 + b^2");
}

/* ------------------------------------------------------------------ */
/* 7. Unevaluated when the antiderivative is unknown (no wrong value).  */
/* ------------------------------------------------------------------ */
static void test_unevaluated(void) {
    /* No closed-form antiderivative -> stays a definite Integrate. */
    check_eq("Head[Integrate[Sin[x^2 + Log[x]] Cos[x], {x, 0, 1}]]", "Integrate");
    /* Gaussian: neither the Erf antiderivative nor Limit[Erf, Inf] is
     * available, so this remains unevaluated rather than guessing. */
    check_eq("Head[Integrate[E^(-x^2), {x, -Infinity, Infinity}]]", "Integrate");
}

/* ------------------------------------------------------------------ */
/* 8. Stress: symbolic definite value vs NIntegrate.                    */
/* ------------------------------------------------------------------ */
static void test_stress_numeric(void) {
    num_match("x^3 - 2 x + 1", "-1", "2");
    num_match("x Exp[x]",      "0",  "1");
    num_match("Tan[x]",        "0",  "Pi/4");
    num_match("Log[x]",        "1",  "2");
    num_match("x/(1 + x^2)",   "0",  "1");
    num_match("1/Sqrt[1 - x^2]", "0", "1/2");
    num_match("Cosh[x]",       "0",  "1");
    num_match("Sinh[x]",       "0",  "2");
    num_match("x^2 Exp[-x]",   "0",  "1");
    num_match("1/(x^2 + 4)",   "0",  "2");
    num_match("1/(x^2 + 1)^2", "0",  "1");
    num_match("x Log[x]",      "1",  "2");
    /* Pole of the integrand outside the interval: still proper. */
    num_match("1/(1 - x)",     "2",  "5");
    num_match("1/x^3",         "1",  "3");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_detector);
    TEST(test_proper);
    TEST(test_improper);
    TEST(test_divergent);
    TEST(test_iterated);
    TEST(test_entry_points);
    TEST(test_unevaluated);
    TEST(test_stress_numeric);

    printf("All Integrate NewtonLeibniz tests passed!\n");
    return 0;
}
