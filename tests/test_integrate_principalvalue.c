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
#include "message.h"

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


/* Evaluate `input`, assert its printed form is `expected` AND that no message
 * fired while it ran -- the direct guard against an internal probe's
 * diagnostic reaching an enclosing Check[] (which would then return its
 * failure branch for a successful computation). */
static void check_quiet(const char* input, const char* expected) {
    unsigned long fired0 = mth_msg_fired_count();
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
    ASSERT_STR_EQ(s, expected);
    if (mth_msg_fired_count() != fired0)
        fprintf(stderr, "FAIL: %s fired a message\n", input);
    ASSERT(mth_msg_fired_count() == fired0);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------------------------------------------------------------------------
 * Symbolic interior pole (x = a, proved interior under the assumptions): a
 * sole one needs no ordering, and a provably SIMPLE one is odd-order, so the
 * PV exists.  1/x over (-a, b): the realness probe now samples a point that
 * satisfies the assumptions, and Refine reads Re[Log[b]] as Log[b] for b > 0.
 * References: the reference system (0, Log[b/a], a (3 + Log[2]), idiv).
 * ---------------------------------------------------------------------- */
static void test_pv_symbolic_pole(void) {
    check_quiet("Integrate[1/(x-a), {x, 0, 2a}, PrincipalValue -> True, Assumptions -> a > 0]", "0");
    check_quiet("Integrate[1/x, {x, -a, b}, PrincipalValue -> True, Assumptions -> a > 0 && b > 0]",
                "-Log[a] + Log[b]");
    check_quiet("Integrate[x/(x-a), {x, 0, 3a}, PrincipalValue -> True, Assumptions -> a > 0]",
                "a (3 + Log[2])");
    /* 0 < a < b: the all-equal sample a = b = 2 used to reject this real PV.
     * The difference is an exact zero (-Log[2] - Log[5/2] + Log[25]/2), which
     * N[.., 30] cannot certify (N::meprec at the precision cap), so it is
     * tested exactly. */
    check_quiet("PossibleZeroQ[(Integrate[1/(x-a), {x, 0, b}, PrincipalValue -> True, "
                "Assumptions -> 0 < a < b] - Log[b/a - 1]) /. {a -> 2, b -> 7}]", "True");
    /* Even-order symbolic pole: no PV. */
    check_eq("Quiet[Integrate[1/(x-a)^2, {x, 0, 2a}, PrincipalValue -> True, Assumptions -> a > 0]]",
             "Integrate[1/(-a + x)^2, {x, 0, 2 a}, PrincipalValue -> True, Assumptions -> a > 0]");
    /* Ordinary integral across the symbolic pole diverges. */
    check_eq("Quiet[Integrate[1/(x-a), {x, 0, 2a}, Assumptions -> a > 0]]",
             "Integrate[1/(-a + x), {x, 0, 2 a}, Assumptions -> a > 0]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_pv_simple);
    TEST(test_pv_even_pole);
    TEST(test_pv_default_off);
    TEST(test_pv_no_pole);
    TEST(test_pv_numeric);

    TEST(test_pv_symbolic_pole);
    printf("All integrate_principalvalue tests passed.\n");
    return 0;
}
