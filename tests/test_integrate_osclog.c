/* test_integrate_osclog.c
 *
 * Tests for the half-line families added to src/calculus/integrate_ramanujan.c:
 *   Integrate`SinPowerMonomial[f, {x,0,Infinity}]   Sin[r x]^k / x^m  (ssp)
 *   Integrate`OscillatoryPower[f, {x,0,Infinity}]   Cos/Sin[b x^n]    (Fresnel)
 *   Integrate`RationalLog[f, {x,0,Infinity}]        R(x) Log[x]^n     (log*rat)
 * and the matching Method -> "SinPowerMonomial" / "OscillatoryPower" /
 * "RationalLog" selectors plus the auto-dispatch (each is a pre-pass / path
 * inside the Ramanujan/Mellin engine).
 *
 * Values verified against the classical closed forms (Int sin^2/x^2 = Pi/2,
 * sin^3/x^3 = 3Pi/8, sin^4/x^4 = Pi/3; the Fresnel Gamma forms; and
 * Int Log x/((x+1)(x+2)) = (Log^2 2)/2).  Negative controls cover the divergent
 * and out-of-scope (positive-real-pole) cases that MUST stay unevaluated.
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
 * Sin[r x]^k / x^m (ssp).
 * ---------------------------------------------------------------------- */
static void test_sinpow_unit(void) {
    check_eq("Integrate`SinPowerMonomial[Sin[x]^2/x^2, {x, 0, Infinity}]", "1/2 Pi");
    check_eq("Integrate`SinPowerMonomial[Sin[x]^3/x^3, {x, 0, Infinity}]", "3/8 Pi");
    check_eq("Integrate`SinPowerMonomial[Sin[x]^4/x^2, {x, 0, Infinity}]", "1/4 Pi");
    check_eq("Integrate`SinPowerMonomial[Sin[x]^4/x^4, {x, 0, Infinity}]", "1/3 Pi");
    check_eq("Integrate`SinPowerMonomial[Sin[x]^6/x^2, {x, 0, Infinity}]", "3/16 Pi");
    /* Frequency scaling: Sin[2x]^2/x^2 = 2 * Sin[y]^2/y^2 = Pi. */
    check_eq("Integrate`SinPowerMonomial[Sin[2 x]^2/x^2, {x, 0, Infinity}]", "Pi");
}

static void test_sinpow_negative(void) {
    /* Sin^2/x diverges at infinity (~ 1/(2x)) -> unevaluated. */
    check_eq("Integrate`SinPowerMonomial[Sin[x]^2/x, {x, 0, Infinity}]",
             "Integrate`SinPowerMonomial[Sin[x]^2/x, {x, 0, Infinity}]");
    /* Not a Sin-power / monomial integrand. */
    check_eq("Integrate`SinPowerMonomial[Cos[x]/x^2, {x, 0, Infinity}]",
             "Integrate`SinPowerMonomial[Cos[x]/x^2, {x, 0, Infinity}]");
}

/* -------------------------------------------------------------------------
 * Cos/Sin[b x^n] (Fresnel-type oscillatory).
 * ---------------------------------------------------------------------- */
static void test_oscpower(void) {
    check_eq("Integrate`OscillatoryPower[Cos[x^2], {x, 0, Infinity}]",
             "(1/2 Sqrt[Pi])/Sqrt[2]");
    check_eq("Integrate`OscillatoryPower[Sin[x^2], {x, 0, Infinity}]",
             "(1/2 Sqrt[Pi])/Sqrt[2]");
    check_eq("Integrate`OscillatoryPower[Cos[2 x^2], {x, 0, Infinity}]", "1/4 Sqrt[Pi]");
    /* n = 3. */
    check_eq("Integrate`OscillatoryPower[Cos[x^3], {x, 0, Infinity}]",
             "(1/3 Pi)/Gamma[2/3]");
}

/* -------------------------------------------------------------------------
 * R(x) Log[x]^n (log*rat) with negative-real-axis poles.
 * ---------------------------------------------------------------------- */
static void test_ratlog_unit(void) {
    check_eq("Integrate`RationalLog[Log[x]/((1+x)(2+x)), {x, 0, Infinity}]",
             "1/2 Log[2]^2");
    check_eq("Integrate`RationalLog[Log[x]/((1+x)(3+x)), {x, 0, Infinity}]",
             "1/4 Log[3]^2");
    /* Second log power (n = 2). */
    check_eq("Integrate`RationalLog[Log[x]^2/((1+x)(2+x)), {x, 0, Infinity}]",
             "1/3 Log[2] (Log[2]^2 + Pi^2)");
    /* Repeated pole. */
    check_eq("Integrate`RationalLog[Log[x]/(1+x)^2, {x, 0, Infinity}]", "0");
}

static void test_ratlog_negative(void) {
    /* Pole on the positive real axis (x = 1): PV only -> unevaluated. */
    check_eq("Integrate`RationalLog[Log[x]/((1-x)(2+x)), {x, 0, Infinity}]",
             "Integrate`RationalLog[Log[x]/((2 + x) (1 - x)), {x, 0, Infinity}]");
    /* No Log weight -> not this mechanism's job. */
    check_eq("Integrate`RationalLog[1/((1+x)(2+x)), {x, 0, Infinity}]",
             "Integrate`RationalLog[1/((1 + x) (2 + x)), {x, 0, Infinity}]");
}

/* -------------------------------------------------------------------------
 * Auto-dispatch + Method selectors.
 * ---------------------------------------------------------------------- */
static void test_dispatch(void) {
    /* Auto: each family is reached through the Ramanujan pass. */
    check_eq("Integrate[Sin[x]^2/x^2, {x, 0, Infinity}]", "1/2 Pi");
    check_eq("Integrate[Log[x]/((1+x)(2+x)), {x, 0, Infinity}]", "1/2 Log[2]^2");
    /* Method selectors. */
    check_eq("Integrate[Sin[x]^3/x^3, {x, 0, Infinity}, Method -> \"SinPowerMonomial\"]",
             "3/8 Pi");
    check_eq("Integrate[Cos[x^2], {x, 0, Infinity}, Method -> \"OscillatoryPower\"]",
             "(1/2 Sqrt[Pi])/Sqrt[2]");
    check_eq("Integrate[Log[x]/((1+x)(2+x)), {x, 0, Infinity}, Method -> \"RationalLog\"]",
             "1/2 Log[2]^2");
    /* Strict: wrong family under a pinned method stays unevaluated. */
    check_eq("Integrate[Cos[x]/x^2, {x, 0, Infinity}, Method -> \"SinPowerMonomial\"]",
             "Integrate[Cos[x]/x^2, {x, 0, Infinity}, Method -> \"SinPowerMonomial\"]");
}

/* -------------------------------------------------------------------------
 * Stress: deeper powers, numeric cross-checks.
 * ---------------------------------------------------------------------- */
static void test_stress(void) {
    check_eq("Integrate`SinPowerMonomial[Sin[x]^5/x^5, {x, 0, Infinity}]", "115/384 Pi");
    check_eq("Integrate`SinPowerMonomial[Sin[x]^8/x^2, {x, 0, Infinity}]", "5/32 Pi");
    /* Distinct triple simple poles close to a clean combination of Log^2. */
    check_eq("Integrate`RationalLog[Log[x]/((1+x)(2+x)(3+x)), {x, 0, Infinity}]",
             "1/4 (-Log[3]^2 + 2 Log[2]^2)");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_sinpow_unit);
    TEST(test_sinpow_negative);
    TEST(test_oscpower);
    TEST(test_ratlog_unit);
    TEST(test_ratlog_negative);
    TEST(test_dispatch);
    TEST(test_stress);

    printf("All integrate_osclog tests passed.\n");
    return 0;
}
