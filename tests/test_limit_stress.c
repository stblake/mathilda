#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Limit assumptions STRESS corpus.
 *
 * Each case is a parametric limit whose value is fixed ONLY once the sign,
 * magnitude, or domain of a symbolic parameter is pinned by an assumption.
 * The assumption reaches Limit three ways -- the `Assumptions -> ...` option,
 * an ambient `Assuming[...]` scope, or `$Assumptions` -- exactly as it does for
 * PossibleZeroQ (see test_possiblezeroq_stress.c). Assumptions that name a
 * domain use the Element[x,dom] function form (the parser has no \[Element]).
 *
 * The corpus doubles as the development fitness function: before the
 * assumption-aware machinery lands, the parametric cases print a divergent
 * DirectedInfinity[...] husk or stay unevaluated; each phase drives the FAIL
 * count down. Cases marked "regression" already hold and must stay put -- in
 * particular a limit with NO informative assumption must print exactly what it
 * does today (a NULL AssumeCtx reproduces the legacy path byte-for-byte).
 *
 * assert_eval_eq uses libc assert() (a no-op under -DNDEBUG); run the binary
 * and grep for FAIL:. Deferred cases are listed at the bottom, NOT asserted.
 *
 * Printing notes (verified against the binary): Pi/2 prints as "1/2 Pi";
 * -Pi/2 as "-1/2 Pi"; a divergent linear limit with unknown coefficient sign
 * prints "DirectedInfinity[c]".
 */

/* Cat 1 -- monomial power x^n, exponent sign decides growth/decay. */
static void test_limit_stress_00(void) {
    assert_eval_eq("Limit[x^n, x -> Infinity, Assumptions -> n > 0]", "Infinity", 0);
    assert_eval_eq("Limit[x^n, x -> Infinity, Assumptions -> n < 0]", "0", 0);
    assert_eval_eq("Limit[x^n, x -> Infinity, Assumptions -> n > 2]", "Infinity", 0);
    assert_eval_eq("Limit[x^n, x -> 0, Direction -> \"FromAbove\", Assumptions -> n > 0]", "0", 0);
    assert_eval_eq("Limit[x^n, x -> 0, Direction -> \"FromAbove\", Assumptions -> n < 0]", "Infinity", 0);
    /* regression: concrete exponent needs no assumption. */
    assert_eval_eq("Limit[x^3, x -> Infinity]", "Infinity", 0);
}

/* Cat 2 -- exponential a^x, base magnitude/sign. Real base a>1 -> Infinity is a
 * DIFFERENT theorem from Abs[a]>1 -> ComplexInfinity (unknown phase). */
static void test_limit_stress_01(void) {
    assert_eval_eq("Limit[a^x, x -> Infinity, Assumptions -> a > 1]", "Infinity", 0);
    assert_eval_eq("Limit[a^x, x -> Infinity, Assumptions -> a > 0 && a < 1]", "0", 0);
    assert_eval_eq("Limit[a^x, x -> -Infinity, Assumptions -> a > 1]", "0", 0);
    assert_eval_eq("Limit[a^x, x -> -Infinity, Assumptions -> a > 0 && a < 1]", "Infinity", 0);
    /* regression / contrast: Abs magnitude gives ComplexInfinity, and concrete
     * bases already resolve with no assumption. */
    assert_eval_eq("Limit[a^x, x -> Infinity, Assumptions -> Abs[a] > 1]", "ComplexInfinity", 0);
    assert_eval_eq("Limit[a^x, x -> Infinity, Assumptions -> Abs[a] < 1]", "0", 0);
    assert_eval_eq("Limit[2^x, x -> Infinity]", "Infinity", 0);
    assert_eval_eq("Limit[(1/2)^x, x -> Infinity]", "0", 0);
}

/* Cat 3 -- log / exponential term sign via Sign[Log[b]] and Sign[c]. */
static void test_limit_stress_02(void) {
    assert_eval_eq("Limit[x Log[a], x -> Infinity, Assumptions -> a > 1]", "Infinity", 0);
    assert_eval_eq("Limit[x Log[a], x -> Infinity, Assumptions -> a > 0 && a < 1]", "-Infinity", 0);
    assert_eval_eq("Limit[E^(x Log[a]), x -> Infinity, Assumptions -> a > 1]", "Infinity", 0);
    assert_eval_eq("Limit[E^(c x), x -> Infinity, Assumptions -> c > 0]", "Infinity", 0);
    assert_eval_eq("Limit[E^(c x), x -> Infinity, Assumptions -> c < 0]", "0", 0);
}

/* Cat 4 -- growth-rate exponent comparison. */
static void test_limit_stress_03(void) {
    assert_eval_eq("Limit[x^n/x^m, x -> Infinity, Assumptions -> n > m]", "Infinity", 0);
    assert_eval_eq("Limit[x^n/x^m, x -> Infinity, Assumptions -> n < m]", "0", 0);
    assert_eval_eq("Limit[(x^n + x^m)/x^m, x -> Infinity, Assumptions -> n > m]", "Infinity", 0);
    assert_eval_eq("Limit[x^a/(x^a + 1), x -> Infinity, Assumptions -> a > 0]", "1", 0);
    /* regression: exp dominates any real power, so these already resolve. */
    assert_eval_eq("Limit[x^n/E^x, x -> Infinity, Assumptions -> n > 0]", "0", 0);
    assert_eval_eq("Limit[E^x/x^n, x -> Infinity, Assumptions -> n > 0]", "Infinity", 0);
}

/* Cat 5 -- Plus dominant balance with a parameter-signed leading coefficient. */
static void test_limit_stress_04(void) {
    assert_eval_eq("Limit[c x, x -> Infinity, Assumptions -> c > 0]", "Infinity", 0);
    assert_eval_eq("Limit[c x, x -> Infinity, Assumptions -> c < 0]", "-Infinity", 0);
    assert_eval_eq("Limit[c x, x -> -Infinity, Assumptions -> c > 0]", "-Infinity", 0);
    assert_eval_eq("Limit[c x^2 + x, x -> Infinity, Assumptions -> c < 0]", "-Infinity", 0);
    assert_eval_eq("Limit[a x + b, x -> Infinity, Assumptions -> a > 0]", "Infinity", 0);
    /* regression: leading term already wins / equal-degree ratio needs no sign. */
    assert_eval_eq("Limit[x^2 + c x, x -> Infinity]", "Infinity", 0);
    assert_eval_eq("Limit[(c x + 1)/(x + 1), x -> Infinity]", "c", 0);
}

/* Cat 6 -- bounded / compose-at-infinity with a parameter. */
static void test_limit_stress_05(void) {
    assert_eval_eq("Limit[Sin[x]/x^p, x -> Infinity, Assumptions -> p > 0]", "0", 0);
    assert_eval_eq("Limit[Cos[x]/x^p, x -> Infinity, Assumptions -> p > 0]", "0", 0);
    assert_eval_eq("Limit[ArcTan[c x], x -> Infinity, Assumptions -> c < 0]", "-1/2 Pi", 0);
    assert_eval_eq("Limit[Tanh[c x], x -> Infinity, Assumptions -> c > 0]", "1", 0);
    assert_eval_eq("Limit[Tanh[c x], x -> Infinity, Assumptions -> c < 0]", "-1", 0);
    /* regression: bounded/x already 0; ArcTan[c x] with c>0 already +Pi/2. */
    assert_eval_eq("Limit[Sin[a x]/x, x -> Infinity]", "0", 0);
    assert_eval_eq("Limit[Sin[x], x -> Infinity]", "Indeterminate", 0);
    assert_eval_eq("Limit[ArcTan[c x], x -> Infinity, Assumptions -> c > 0]", "1/2 Pi", 0);
}

/* Cat 7 -- integer / domain assumptions. */
static void test_limit_stress_06(void) {
    assert_eval_eq("Limit[x^n, x -> Infinity, Assumptions -> Element[n, Integers] && n > 0]", "Infinity", 0);
    assert_eval_eq("Limit[x^n, x -> Infinity, Assumptions -> Element[n, PositiveIntegers]]", "Infinity", 0);
    assert_eval_eq("Limit[x^n, x -> 0, Direction -> \"FromAbove\", Assumptions -> Element[n, PositiveIntegers]]", "0", 0);
    assert_eval_eq("Limit[x^n, x -> -Infinity, Assumptions -> Element[n, Evens] && n > 0]", "Infinity", 0);
}

/* Cat 8 -- ambient Assuming[] scope (incl. nested). */
static void test_limit_stress_07(void) {
    assert_eval_eq("Assuming[c > 0, Limit[c x, x -> Infinity]]", "Infinity", 0);
    assert_eval_eq("Assuming[a > 1, Limit[a^x, x -> Infinity]]", "Infinity", 0);
    assert_eval_eq("Assuming[Element[n, Integers] && n > 0, Limit[x^n, x -> Infinity]]", "Infinity", 0);
    assert_eval_eq("Assuming[a > 0 && a < 1, Limit[a^x, x -> Infinity]]", "0", 0);
    assert_eval_eq("Assuming[c > 0, Assuming[c < 10, Limit[c x, x -> Infinity]]]", "Infinity", 0);
}

/* Cat 9/10 -- inconsistent/trivial assumptions and directional + parameter sign. */
static void test_limit_stress_08(void) {
    /* Cat 9: an uninformative or contradictory assumption falls through. */
    assert_eval_eq("Limit[x, x -> Infinity, Assumptions -> False]", "Infinity", 0);
    assert_eval_eq("Limit[c x, x -> Infinity, Assumptions -> True]", "DirectedInfinity[c]", 0);
    /* Cat 10: one-sided pole with a signed residue. */
    assert_eval_eq("Limit[c/x, x -> 0, Direction -> \"FromAbove\", Assumptions -> c > 0]", "Infinity", 0);
    assert_eval_eq("Limit[c/x, x -> 0, Direction -> \"FromAbove\", Assumptions -> c < 0]", "-Infinity", 0);
    assert_eval_eq("Limit[c/x, x -> 0, Direction -> \"FromBelow\", Assumptions -> c > 0]", "-Infinity", 0);
    assert_eval_eq("Limit[c/x, x -> 0, Direction -> \"FromBelow\", Assumptions -> c < 0]", "Infinity", 0);
}

/* Cat 11 -- passthrough / regression. NO informative assumption => today's
 * behaviour, exactly (NULL AssumeCtx reproduces the legacy path). */
static void test_limit_stress_09(void) {
    assert_eval_eq("Limit[x^n, x -> Infinity]", "Limit[x^n, x -> Infinity]", 0);
    assert_eval_eq("Limit[a^x, x -> Infinity]", "Limit[a^x, x -> Infinity]", 0);
    assert_eval_eq("Limit[c x, x -> Infinity]", "DirectedInfinity[c]", 0);
    assert_eval_eq("Limit[1/x, x -> Infinity]", "0", 0);
    assert_eval_eq("Limit[x^2, x -> Infinity]", "Infinity", 0);
    assert_eval_eq("Limit[1/(x - 1), x -> 1, Direction -> \"FromAbove\"]", "Infinity", 0);
}

/* Deferred (documented; not asserted):
 *   Limit[x^n Log[x], x -> Infinity, Assumptions -> n > 0]        (needs product growth+log)
 *   Limit[a^x + b^x, x -> Infinity, Assumptions -> a > b > 1]     (two-term exp dominance)
 *   Limit[(1 + a/x)^x, x -> Infinity, Assumptions -> a > 0]       (parametric 1^inf)
 *   Limit[x^a - x^b, x -> Infinity, Assumptions -> a > b > 0]     (signed Plus of unbounded)
 */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_limit_stress_00);
    TEST(test_limit_stress_01);
    TEST(test_limit_stress_02);
    TEST(test_limit_stress_03);
    TEST(test_limit_stress_04);
    TEST(test_limit_stress_05);
    TEST(test_limit_stress_06);
    TEST(test_limit_stress_07);
    TEST(test_limit_stress_08);
    TEST(test_limit_stress_09);

    printf("All limit stress tests passed! (62 cases)\n");
    return 0;
}
