/* test_integrate_beta.c
 *
 * Tests for definite integration by Euler-Beta reduction
 * (src/calculus/integrate_beta.c):
 *   Integrate`Beta[f, {x, 0, 1}]           x^(k-1)(1-x)^(l-1) -> Beta[k,l]
 *   Integrate`TrigPower[f, {x, 0, c}]      Sin^m Cos^n -> Beta reduction
 *   Method -> "Beta" / "TrigPower"         pinned mechanism (strict)
 *
 * Coverage: the non-elementary cases the FTC cannot reach (non-integer /
 * symbolic exponents, even powers over [0,Pi]), the Log[x]/Log[1-x] parameter
 * derivatives, the parity vanishing over [0,Pi]/[0,2Pi], the convergence-strip
 * gate (ConditionalExpression when undecided, NULL when refuted), the
 * auto-dispatch regression (integer cases still owned by Newton-Leibniz), and
 * the negative controls that MUST stay unevaluated.
 *
 * Values Mathilda leaves in an equivalent-but-unsimplified form are pinned via
 * Simplify[value - reference] == 0.
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
 * Beta on [0,1] -- explicit entry point.
 * ---------------------------------------------------------------------- */
static void test_beta_unit(void) {
    /* Integer exponents: the Beta reduction gives the same rational as FTC. */
    check_eq("Integrate`Beta[x^2 (1-x)^3, {x, 0, 1}]", "1/60");
    /* Half-integer: FTC cannot (incomplete Beta), Beta gives Pi/8. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2), {x, 0, 1}]", "1/8 Pi");
    /* Negative (but > -1) exponent of x: still convergent -> Pi/2. */
    check_eq("Integrate`Beta[x^(-1/2)(1-x)^(1/2), {x, 0, 1}]", "1/2 Pi");
    /* Symbolic exponents -> ConditionalExpression on the strip. */
    check_eq("Integrate`Beta[x^(k-1)(1-x)^(l-1), {x, 0, 1}, "
             "Assumptions -> k>0 && l>0]", "Beta[k, l]");
    check_eq("Integrate`Beta[x^(k-1)(1-x)^(l-1), {x, 0, 1}]",
             "ConditionalExpression[Beta[k, l], Re[k] > 0 && Re[l] > 0]");
}

static void test_beta_logweight(void) {
    /* Log[x] weight = d/dk Beta[k,l]. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2) Log[x], {x, 0, 1}]",
             "1/8 Pi (1/2 - 2 Log[2])");
    /* Log[1-x] weight = d/dl Beta[k,l]. */
    check_eq("Integrate`Beta[x^2 (1-x)^3 Log[1-x], {x, 0, 1}]", "-37/3600");
    /* Mixed Log[x] Log[1-x] weight = d^2/dk dl Beta. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2) Log[x] Log[1-x], {x, 0, 1}]",
             "1/48 Pi (9 - 12 Log[2] + 24 Log[2]^2 - Pi^2)");
}

static void test_beta_negative(void) {
    /* Divergent (k = -1 <= 0): the strip is refuted -> unevaluated. */
    check_eq("Integrate`Beta[x^(-2)(1-x)^(1/2), {x, 0, 1}]",
             "Integrate`Beta[Sqrt[1 - x]/x^2, {x, 0, 1}]");
    /* Not of the x^p (1-x)^q form. */
    check_eq("Integrate`Beta[1/(1+x), {x, 0, 1}]",
             "Integrate`Beta[1/(1 + x), {x, 0, 1}]");
    /* Wrong interval. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2), {x, 0, 2}]",
             "Integrate`Beta[Sqrt[x] Sqrt[1 - x], {x, 0, 2}]");
}

/* -------------------------------------------------------------------------
 * Sin^m Cos^n over the canonical trig intervals -- explicit entry point.
 * ---------------------------------------------------------------------- */
static void test_trigpower_unit(void) {
    /* [0, Pi/2] half-integer -> Gamma[3/4]^2/Sqrt[Pi] = (1/2) Beta[3/4,3/4]. */
    check_eq("Integrate`TrigPower[Sin[x]^(1/2) Cos[x]^(1/2), {x, 0, Pi/2}]",
             "Gamma[3/4]^2/Sqrt[Pi]");
    /* [0, Pi/2] with a constant factor. */
    check_eq("Integrate`TrigPower[3 Sin[x]^4 Cos[x]^2, {x, 0, Pi/2}]", "3/32 Pi");
    /* [0, Pi]: n even -> 2 x quarter. */
    check_eq("Integrate`TrigPower[Sin[x]^4 Cos[x]^2, {x, 0, Pi}]", "1/16 Pi");
    /* [0, Pi]: n even, m odd (m parity irrelevant) -> 4/15. */
    check_eq("Integrate`TrigPower[Sin[x]^3 Cos[x]^2, {x, 0, Pi}]", "4/15");
    /* [0, 2Pi]: both even -> 4 x quarter. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^2, {x, 0, 2 Pi}]", "1/4 Pi");
    /* Symbolic exponents over [0, Pi/2].  The affine strip (1+m)/2 > 0 is not
     * discharged by Simplify from m > 0 (no linear-implication reasoning), so
     * the honest result is a ConditionalExpression -- as for the residue
     * 1/(1+x^n) family. */
    check_eq("Integrate`TrigPower[Sin[x]^m Cos[x]^n, {x, 0, Pi/2}, "
             "Assumptions -> m>0 && n>0]",
             "ConditionalExpression[1/2 Beta[1/2 (1 + m), 1/2 (1 + n)], "
             "Re[1/2 (1 + m)] > 0 && Re[1/2 (1 + n)] > 0]");
    /* Over [0, Pi] only n parity matters, so m may stay symbolic (n = 2 even). */
    check_eq("Integrate`TrigPower[Sin[x]^m Cos[x]^2, {x, 0, Pi}, "
             "Assumptions -> m>0]",
             "ConditionalExpression[Beta[1/2 (1 + m), 3/2], Re[1/2 (1 + m)] > 0]");
}

static void test_trigpower_vanishing(void) {
    /* [0, Pi]: n odd -> 0. */
    check_eq("Integrate`TrigPower[Sin[x]^4 Cos[x]^3, {x, 0, Pi}]", "0");
    /* [0, 2Pi]: m odd -> 0. */
    check_eq("Integrate`TrigPower[Sin[x]^5 Cos[x]^3, {x, 0, 2 Pi}]", "0");
    /* [0, 2Pi]: n odd -> 0. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^5, {x, 0, 2 Pi}]", "0");
}

static void test_trigpower_negative(void) {
    /* Symbolic n over [0, Pi]: parity undecidable -> unevaluated. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^n, {x, 0, Pi}]",
             "Integrate`TrigPower[Sin[x]^2 Cos[x]^n, {x, 0, Pi}]");
    /* Non-canonical interval. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^2, {x, 0, Pi/3}]",
             "Integrate`TrigPower[Cos[x]^2 Sin[x]^2, {x, 0, 1/3 Pi}]");
}

/* -------------------------------------------------------------------------
 * Auto-dispatch: new cases close; the integer cases stay owned by FTC.
 * ---------------------------------------------------------------------- */
static void test_dispatch_regression(void) {
    /* NEW: previously unevaluated, now closed by the Beta mechanism. */
    check_eq("Integrate[x^(1/2)(1-x)^(1/2), {x, 0, 1}]", "1/8 Pi");
    check_eq("Integrate[Sin[x]^(1/2) Cos[x]^(1/2), {x, 0, Pi/2}]",
             "Gamma[3/4]^2/Sqrt[Pi]");
    check_eq("Integrate[Sin[x]^4 Cos[x]^2, {x, 0, Pi}]", "1/16 Pi");
    /* Regression: integer-power cases still evaluate (via FTC), same value. */
    check_eq("Integrate[x^2 (1-x)^3, {x, 0, 1}]", "1/60");
    check_eq("Integrate[Sin[x]^4 Cos[x]^2, {x, 0, Pi/2}]", "1/32 Pi");
}

static void test_method_option(void) {
    check_eq("Integrate[x^(1/2)(1-x)^(1/2), {x, 0, 1}, Method -> \"Beta\"]",
             "1/8 Pi");
    check_eq("Integrate[Sin[x]^(1/2) Cos[x]^(1/2), {x, 0, Pi/2}, "
             "Method -> \"TrigPower\"]", "Gamma[3/4]^2/Sqrt[Pi]");
    /* Strict: a non-Beta integrand under Method -> "Beta" stays unevaluated. */
    check_eq("Integrate[1/(1+x), {x, 0, 1}, Method -> \"Beta\"]",
             "Integrate[1/(1 + x), {x, 0, 1}, Method -> \"Beta\"]");
}

/* -------------------------------------------------------------------------
 * Stress: high powers, deep half-integer, double log weight.
 * ---------------------------------------------------------------------- */
static void test_stress(void) {
    check_eq("Integrate`Beta[x^10 (1-x)^7, {x, 0, 1}]", "1/350064");
    check_eq("Integrate`Beta[x^(3/2)(1-x)^(5/2), {x, 0, 1}]", "3/256 Pi");
    check_eq("Integrate`TrigPower[Sin[x]^10 Cos[x]^8, {x, 0, Pi/2}]",
             "35/131072 Pi");
    /* Deep half-integer trig closes to a clean multiple of Pi. */
    check_eq("Simplify[Integrate`TrigPower[Sin[x]^10 Cos[x]^8, {x, 0, Pi/2}] "
             "- 1/2 Beta[11/2, 9/2]]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_beta_unit);
    TEST(test_beta_logweight);
    TEST(test_beta_negative);
    TEST(test_trigpower_unit);
    TEST(test_trigpower_vanishing);
    TEST(test_trigpower_negative);
    TEST(test_dispatch_regression);
    TEST(test_method_option);
    TEST(test_stress);

    printf("All integrate_beta tests passed.\n");
    return 0;
}
