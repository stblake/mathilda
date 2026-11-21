/*
 * test_series_assumptions.c -- targeted, per-feature tests for the Series
 * Assumptions / Assuming[] wiring (companion to the broad test_series_stress.c,
 * mirroring test_limit_assumptions.c beside test_limit_stress.c).
 *
 * Series honours an assumption from three channels -- the `Assumptions -> ...`
 * option, an ambient `Assuming[...]` scope, and a direct `$Assumptions`
 * assignment via Block[] -- with the option overriding the ambient. A NULL /
 * uninformative / inconsistent assumption reproduces the legacy path
 * byte-for-byte. These tests pin the wiring invariants that are most likely to
 * regress silently, independent of the mathematical coverage in the stress file.
 *
 * assert_fullform uses libc assert() (a no-op under -DNDEBUG); run the binary and
 * grep for FAIL:.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static bool g_setup_done = false;
static void setup_full(void) {
    symtab_init();
    core_init();
    if (!g_setup_done) {
        const char* ups[] = { ".", "..", "../..", "../../..", "../../../..", NULL };
        for (int i = 0; ups[i]; i++) {
            char path[256];
            snprintf(path, sizeof(path), "%s/src/internal/init.m", ups[i]);
            if (access(path, F_OK) == 0) { (void)!chdir(ups[i]); break; }
        }
        g_setup_done = true;
    }
    Expr* c = parse_expression("Get[\"src/internal/init.m\"]");
    if (c) { Expr* r = evaluate(c); expr_free(c); if (r) expr_free(r); }
}
static void assert_fullform(const char* input, const char* expected) {
    assert_eval_eq(input, expected, 1);
}

/* Constants used repeatedly: the a > 0 result and the legacy (no-assumption)
 * husk for Series[Sqrt[a^2 + x], {x, 0, 1}]. */
#define A_POS "SeriesData[x, 0, List[a, Times[Rational[1, 2], Power[a, -1]]], 0, 2, 1]"
#define A_HUSK "SeriesData[x, 0, List[Power[Power[a, 2], Rational[1, 2]], " \
               "Times[Rational[1, 2], Power[a, -2], Power[Power[a, 2], Rational[1, 2]]]], 0, 2, 1]"

/* Channel 1: the Assumptions option. */
static void test_option_channel(void) {
    setup_full();
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> a > 0]", A_POS);
}

/* Channel 2: an ambient Assuming[...] scope. */
static void test_assuming_channel(void) {
    setup_full();
    assert_fullform("Assuming[a > 0, Series[Sqrt[a^2 + x], {x, 0, 1}]]", A_POS);
}

/* Channel 3: a direct $Assumptions assignment via Block[]. */
static void test_dollar_assumptions_channel(void) {
    setup_full();
    assert_fullform("Block[{$Assumptions = a > 0}, Series[Sqrt[a^2 + x], {x, 0, 1}]]", A_POS);
}

/* Precedence: the option overrides an ambient scope of the opposite sign. */
static void test_option_overrides_ambient(void) {
    setup_full();
    assert_fullform("Assuming[a < 0, Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> a > 0]]", A_POS);
}

/* NULL / uninformative / inconsistent assumptions reproduce the legacy husk
 * exactly (byte-for-byte). */
static void test_null_passthrough(void) {
    setup_full();
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}]", A_HUSK);
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> True]", A_HUSK);
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> Automatic]", A_HUSK);
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> False]", A_HUSK);
}

/* Negative controls: an assumption that leaves the parameter's sign/reality
 * unknown must not collapse Sqrt[a^2]. */
static void test_negative_controls(void) {
    setup_full();
    /* unrelated symbol */
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> b > 0]", A_HUSK);
    /* complex domain: a could be complex, so no branch decision is sound */
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> Element[a, Complexes]]", A_HUSK);
    /* real-but-unknown-sign yields Abs[a], not a */
    assert_fullform("Series[Sqrt[a^2 + x], {x, 0, 1}, Assumptions -> Element[a, Reals]]",
        "SeriesData[x, 0, List[Abs[a], Times[Rational[1, 2], Power[a, -2], Abs[a]]], 0, 2, 1]");
}

/* Expansion-variable sign: the ambient channels reach the Log-branch decision
 * (previously option-only). */
static void test_expansion_variable_sign(void) {
    setup_full();
    assert_fullform("Assuming[x < 0, Series[ExpIntegralEi[x], {x, 0, 2}]]",
        "SeriesData[x, 0, List[Plus[EulerGamma, Log[Times[-1, x]]], 1, Rational[1, 4]], 0, 3, 1]");
    assert_fullform("Block[{$Assumptions = x < 0}, Series[Abs[x], {x, 0, 2}]]",
        "SeriesData[x, 0, List[0, -1, 0], 0, 3, 1]");
}

/* Multivariate forwarding: the option reaches the inner-variable recursion. */
static void test_multivariate_forwarding(void) {
    setup_full();
    assert_fullform("Series[Sqrt[a^2] + x + y, {x, 0, 1}, {y, 0, 1}, Assumptions -> a > 0]",
        "SeriesData[x, 0, List[SeriesData[y, 0, List[a, 1], 0, 2, 1], 1], 0, 2, 1]");
}

/* Abs[g] with g nonzero at the expansion point: the branch is fixed by the
 * sign of g(x0) under the assumptions (the generic sign of x - a is undecidable
 * in x). Needs the deep Reduce bridge, since the cheap prov_* layer does not
 * sign a compound value like -a. The no-assumption path keeps the Derivative
 * husk byte-for-byte. */
static void test_abs_sign_at_point(void) {
    setup_full();
    assert_fullform("Series[Abs[x - a], {x, 0, 2}, Assumptions -> a > 0]",
        "SeriesData[x, 0, List[a, -1, 0], 0, 3, 1]");                      /* a - x */
    assert_fullform("Series[Abs[x - a], {x, 0, 2}, Assumptions -> a < 0]",
        "SeriesData[x, 0, List[Times[-1, a], 1, 0], 0, 3, 1]");           /* x - a */
    /* negative control: no assumption => legacy Derivative-of-Abs husk (its
     * value term is Abs[a]: Abs[-a] evaluates, a product's numeric factor
     * coming out of Abs) */
    assert_fullform("Series[Abs[x - a], {x, 0, 2}]",
        "SeriesData[x, 0, List[Abs[a], Derivative[1][Abs][Times[-1, a]], "
        "Times[Rational[1, 2], Derivative[2][Abs][Times[-1, a]]]], 0, 3, 1]");
}

/* Square root of an even-order zero: the leading monomial reduces as |x|^e, so
 * an odd e flips the whole series sign for x < 0. The leading exponent sits
 * above the stored nmin (1-Cos[x] has nmin=0, coefs {0,0,1/2,...}). x > 0 is
 * the principal + root (unchanged). */
static void test_even_order_root_sign(void) {
    setup_full();
    assert_fullform("Series[Sqrt[1 - Cos[x]], {x, 0, 2}, Assumptions -> x > 0]",
        "SeriesData[x, 0, List[Power[2, Rational[-1, 2]], 0], 1, 3, 1]");       /* x/Sqrt2 */
    assert_fullform("Series[Sqrt[1 - Cos[x]], {x, 0, 2}, Assumptions -> x < 0]",
        "SeriesData[x, 0, List[Times[-1, Power[2, Rational[-1, 2]]], 0], 1, 3, 1]"); /* -x/Sqrt2 */
    assert_fullform("Series[Sqrt[x^4 + a x^2], {x, 0, 3}, Assumptions -> a > 0 && x > 0]",
        "SeriesData[x, 0, List[Power[a, Rational[1, 2]], 0, "
        "Times[Rational[1, 2], Power[a, Rational[-1, 2]]]], 1, 4, 1]");
    assert_fullform("Series[Sqrt[x^4 + a x^2], {x, 0, 3}, Assumptions -> a > 0 && x < 0]",
        "SeriesData[x, 0, List[Times[-1, Power[a, Rational[1, 2]]], 0, "
        "Times[Rational[-1, 2], Power[a, Rational[-1, 2]]]], 1, 4, 1]");
}

/* Square root of a provably-negative quantity takes the imaginary principal
 * branch Sqrt[a] -> I Sqrt[-a] (via the shared apply_assumption_rules). a > 0
 * keeps the real branch. */
static void test_sqrt_negative_branch(void) {
    setup_full();
    assert_fullform("Series[Sqrt[a - x], {x, 0, 1}, Assumptions -> a > 0]",
        "SeriesData[x, 0, List[Power[a, Rational[1, 2]], "
        "Times[Rational[-1, 2], Power[a, Rational[-1, 2]]]], 0, 2, 1]");
    assert_fullform("Series[Sqrt[a - x], {x, 0, 1}, Assumptions -> a < 0]",
        "SeriesData[x, 0, List[Times[Complex[0, 1], Power[Times[-1, a], Rational[1, 2]]], "
        "Times[Complex[0, Rational[1, 2]], Power[Times[-1, a], Rational[-1, 2]]]], 0, 2, 1]");
}

/* Expansion at -Infinity by reflection: the result is the +Infinity twin with
 * the odd-exponent terms flipped, in the canonical 1/x variable. */
static void test_expansion_at_neg_infinity(void) {
    setup_full();
    /* +Infinity twin, for reference */
    assert_fullform("Series[Sqrt[x^2 + a x], {x, Infinity, 1}, Assumptions -> Element[a, Reals]]",
        "SeriesData[Power[x, -1], 0, List[1, Times[Rational[1, 2], a], "
        "Times[Rational[-1, 8], Power[a, 2]]], -1, 2, 1]");
    assert_fullform("Series[Sqrt[x^2 + a x], {x, -Infinity, 1}, Assumptions -> Element[a, Reals]]",
        "SeriesData[Power[x, -1], 0, List[-1, Times[Rational[-1, 2], a], "
        "Times[Rational[1, 8], Power[a, 2]]], -1, 2, 1]");
    assert_fullform("Series[Sqrt[a^2 x^2 + b x], {x, -Infinity, 1}, "
                    "Assumptions -> a > 0 && Element[b, Reals]]",
        "SeriesData[Power[x, -1], 0, List[Times[-1, a], "
        "Times[Rational[-1, 2], Power[a, -1], b], "
        "Times[Rational[1, 8], Power[a, -3], Power[b, 2]]], -1, 2, 1]");
}

/* Parametric (scaled) arguments at infinity: the sign of the leading coefficient
 * selects the branch / direction. ArcTan coefficient sign; saturating heads
 * (Tanh/Erf/Erfc/Fresnel) limit to a sign-dependent constant; Ei/E_1 carry the
 * scaled asymptotic. */
static void test_parametric_at_infinity(void) {
    setup_full();
    /* ArcTan[a x] honours the sign of a (and the direction) */
    assert_fullform("Series[ArcTan[a x], {x, Infinity, 1}, Assumptions -> a < 0]",
        "SeriesData[Power[x, -1], 0, List[Times[Rational[-1, 2], Pi], "
        "Times[-1, Power[a, -1]]], 0, 2, 1]");                         /* -Pi/2 - 1/(a x) */
    assert_fullform("Series[ArcTan[a x], {x, -Infinity, 1}, Assumptions -> a < 0]",
        "SeriesData[Power[x, -1], 0, List[Times[Rational[1, 2], Pi], "
        "Times[-1, Power[a, -1]]], 0, 2, 1]");                         /*  Pi/2 - 1/(a x) */
    /* saturating heads -> sign-dependent constant limit */
    assert_fullform("Series[Tanh[a x], {x, Infinity, 1}, Assumptions -> a > 0]", "1");
    assert_fullform("Series[Tanh[a x], {x, Infinity, 1}, Assumptions -> a < 0]", "-1");
    assert_fullform("Series[Erf[a x], {x, Infinity, 0}, Assumptions -> a > 0]", "1");
    assert_fullform("Series[Erf[a x], {x, Infinity, 0}, Assumptions -> a < 0]", "-1");
    assert_fullform("Series[Erfc[a x], {x, Infinity, 0}, Assumptions -> a > 0]", "0");
    assert_fullform("Series[Erfc[a x], {x, Infinity, 0}, Assumptions -> a < 0]", "2");
    assert_fullform("Series[FresnelS[a x], {x, Infinity, 0}, Assumptions -> a > 0]",
        "Rational[1, 2]");
    /* growth/decay scaled asymptotics */
    assert_fullform("Series[ExpIntegralEi[a x], {x, Infinity, 1}, Assumptions -> a > 0]",
        "SeriesData[Times[Power[a, -1], Power[x, -1]], 0, "
        "List[Power[E, Times[a, x]]], 1, 2, 1]");                      /* E^(a x)/(a x) */
    assert_fullform("Series[ExpIntegralE[1, a x], {x, Infinity, 1}, Assumptions -> a > 0]",
        "SeriesData[Times[Power[a, -1], Power[x, -1]], 0, "
        "List[Power[E, Times[-1, a, x]]], 1, 2, 1]");                  /* E^(-a x)/(a x) */
    /* bare E_1 asymptotic (new emitter) */
    assert_fullform("Series[ExpIntegralE[1, x], {x, Infinity, 1}]",
        "Times[Power[E, Times[-1, x]], SeriesData[Power[x, -1], 0, List[1], 1, 2, 1]]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_option_channel);
    TEST(test_assuming_channel);
    TEST(test_dollar_assumptions_channel);
    TEST(test_option_overrides_ambient);
    TEST(test_null_passthrough);
    TEST(test_negative_controls);
    TEST(test_expansion_variable_sign);
    TEST(test_multivariate_forwarding);
    TEST(test_abs_sign_at_point);
    TEST(test_even_order_root_sign);
    TEST(test_sqrt_negative_branch);
    TEST(test_expansion_at_neg_infinity);
    TEST(test_parametric_at_infinity);

    printf("All series assumptions tests passed!\n");
    return 0;
}
