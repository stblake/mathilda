/* test_integrate_principalvalue.c
 *
 * Tests for the Cauchy principal value of a definite integral with an interior
 * pole (src/calculus/integrate_newton_leibniz.c, integrate_newton_leibniz_try_pv):
 *   Integrate[f, {x, a, b}, PrincipalValue -> True]
 *
 * PV is valid only for odd-order interior poles (the integrand changes sign
 * across them, so the one-sided divergences cancel); an even-order pole has no
 * principal value and stays unevaluated (Integrate::idiv).  With no interior
 * pole the option is a no-op.  Values verified against the closed forms
 * (PV of 1/x over [-1,1] = 0, of 1/(x-2) over [0,3] = -Log 2, etc.).
 *
 * The value carries the antiderivative's own surface form: a Log-based
 * antiderivative reduces cleanly, an ArcTanh-based one is left as an
 * (equivalent) Re[...] expression, pinned numerically here.
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
 * Odd-order interior poles: the principal value exists and closes.
 * ---------------------------------------------------------------------- */
static void test_pv_simple(void) {
    /* Symmetric simple pole: the divergences cancel, PV = 0. */
    check_eq("Integrate[1/x, {x, -1, 1}, PrincipalValue -> True]", "0");
    /* Off-centre simple pole. */
    check_eq("Integrate[1/(x-2), {x, 0, 3}, PrincipalValue -> True]", "-Log[2]");
    check_eq("Integrate[1/(x-1), {x, 0, 3}, PrincipalValue -> True]", "Log[2]");
    /* Rational with a Log-form antiderivative. */
    check_eq("Integrate[x/(x^2-1), {x, 0, 2}, PrincipalValue -> True]", "1/2 Log[3]");
    /* Order-3 (odd) pole: PV exists, and here it is 0 by symmetry. */
    check_eq("Integrate[1/x^3, {x, -1, 1}, PrincipalValue -> True]", "0");
}

/* -------------------------------------------------------------------------
 * Even-order interior poles: no principal value -> unevaluated (idiv).
 * ---------------------------------------------------------------------- */
static void test_pv_even_pole(void) {
    check_eq("Integrate[1/x^2, {x, -1, 1}, PrincipalValue -> True]",
             "Integrate[1/x^2, {x, -1, 1}, PrincipalValue -> True]");
    check_eq("Integrate[1/(x-2)^2, {x, 0, 3}, PrincipalValue -> True]",
             "Integrate[1/(-2 + x)^2, {x, 0, 3}, PrincipalValue -> True]");
}

/* -------------------------------------------------------------------------
 * The option is off by default: a genuinely divergent integral stays
 * unevaluated without PrincipalValue -> True.
 * ---------------------------------------------------------------------- */
static void test_pv_default_off(void) {
    check_eq("Integrate[1/x, {x, -1, 1}]",
             "Integrate[1/x, {x, -1, 1}]");
    check_eq("Integrate[1/x, {x, -1, 1}, PrincipalValue -> False]",
             "Integrate[1/x, {x, -1, 1}, PrincipalValue -> False]");
}

/* -------------------------------------------------------------------------
 * No interior pole: PrincipalValue -> True is a no-op (ordinary FTC value).
 * ---------------------------------------------------------------------- */
static void test_pv_no_pole(void) {
    check_eq("Integrate[1/x, {x, 1, 2}, PrincipalValue -> True]", "Log[2]");
    check_eq("Integrate[x^2, {x, 0, 2}, PrincipalValue -> True]", "8/3");
}

/* -------------------------------------------------------------------------
 * Numeric cross-checks, incl. an ArcTanh-form antiderivative (correct value in
 * an equivalent Re[...] surface form).
 * ---------------------------------------------------------------------- */
static void test_pv_numeric(void) {
    /* Two simple interior poles: PV = -2 Log 2 (surface form Re[-4 ArcTanh[3]]). */
    check_eq("Chop[N[Integrate[1/((x-1)(x-2)), {x, 0, 3}, PrincipalValue -> True] "
             "+ 2 Log[2]]]", "0");
    /* 1/(x^2-1) over [0,2]: PV = -(1/2) Log 3. */
    check_eq("Chop[N[Integrate[1/(x^2-1), {x, 0, 2}, PrincipalValue -> True] "
             "+ 1/2 Log[3]]]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_pv_simple);
    TEST(test_pv_even_pole);
    TEST(test_pv_default_off);
    TEST(test_pv_no_pole);
    TEST(test_pv_numeric);

    printf("All integrate_principalvalue tests passed.\n");
    return 0;
}
