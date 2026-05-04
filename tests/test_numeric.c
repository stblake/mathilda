/* Tests for N[] and related numeric-evaluation builtins.
 *
 * Phase 1 here covers the machine-precision path. Phase 2 tests
 * (MPFR, precision literals, Precision/Accuracy/SetPrecision/SetAccuracy)
 * are added in a follow-up section once MPFR support lands.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  Phase 1: machine-precision N[]
 * ---------------------------------------------------------------------- */

static void test_n_integer_leaf(void) {
    assert_eval_eq("N[0]",  "0.0", 0);
    assert_eval_eq("N[5]",  "5.0", 0);
    assert_eval_eq("N[-7]", "-7.0", 0);
}

static void test_n_rational_leaf(void) {
    assert_eval_eq("N[1/2]", "0.5", 0);
    assert_eval_eq("N[1/3]", "0.333333", 0);
    assert_eval_eq("N[1/7]", "0.142857", 0);
}

static void test_n_real_leaf(void) {
    /* Already inexact — N is a no-op on value. */
    assert_eval_eq("N[3.14]", "3.14", 0);
}

static void test_n_constants(void) {
    assert_eval_eq("N[Pi]",          "3.14159", 0);
    assert_eval_eq("N[E]",           "2.71828", 0);
    assert_eval_eq("N[EulerGamma]",  "0.577216", 0);
    assert_eval_eq("N[Catalan]",     "0.915966", 0);
    assert_eval_eq("N[GoldenRatio]", "1.61803", 0);
    /* Degree = Pi / 180 ≈ 0.0174533 */
    assert_eval_eq("N[Degree]",      "0.0174533", 0);
}

static void test_n_unknown_symbols_pass_through(void) {
    assert_eval_eq("N[x]",         "x", 0);
    assert_eval_eq("N[myFancy]",   "myFancy", 0);
    assert_eval_eq("N[Infinity]",  "Infinity", 0);
}

static void test_n_on_sin_cos(void) {
    assert_eval_eq("N[Sin[1]]", "0.841471", 0);
    assert_eval_eq("N[Cos[1]]", "0.540302", 0);
    /* Sin[Pi] evaluates to 0 (exact) *before* N sees it, so N[0] = 0.0. */
    assert_eval_eq("N[Sin[Pi]]", "0.0", 0);
}

static void test_n_descend_then_reevaluate(void) {
    /* N[Pi + 1] should numericalize Pi, rebuild Plus[3.14159, 1],
     * re-evaluate, and end up at 4.14159. */
    assert_eval_eq("N[Plus[Pi, 1]]", "4.14159", 0);

    /* For an unknown head, N descends into args and re-evaluates: the
     * result keeps the head. */
    assert_eval_eq("N[foo[Pi, x]]", "foo[3.14159, x]", 0);
}

static void test_n_listable(void) {
    assert_eval_eq("N[{1/3, Pi, E}]", "{0.333333, 3.14159, 2.71828}", 0);
}

static void test_n_complex(void) {
    /* Use the is_fullform=1 path to compare against FullForm output; do
     * NOT also wrap the input in FullForm[...] (that'd double-wrap). */
    assert_eval_eq("N[Complex[1, 2]]", "Complex[1.0, 2.0]", 1);
    assert_eval_eq("N[Complex[1/3, 1/7]]",
                   "Complex[0.333333, 0.142857]", 1);
}

static void test_n_preserves_hold(void) {
    /* HoldForm[Pi] is *kept* by N (not reduced to Pi numerically).
     * FullForm output reveals the structure; the plain printer elides HoldForm. */
    assert_eval_eq("N[HoldForm[Pi]]", "HoldForm[Pi]", 1);
    assert_eval_eq("N[Hold[Pi + 1]]", "Hold[Plus[Pi, 1]]", 1);
}

static void test_n_nested(void) {
    /* Nested structure with several constants and a user head. */
    assert_eval_eq("N[foo[Pi, {1/3, E}]]",
                   "foo[3.14159, List[0.333333, 2.71828]]", 1);
}

static void test_n_preserves_integer_exponents(void) {
    /* Mathematica-compatible: N[] leaves integer exponents alone, since
     * polynomial structure (and the integer power's meaning as repeated
     * multiplication) is more useful than a real-valued exponent. */
    assert_eval_eq("N[x^2]", "x^2", 0);
    assert_eval_eq("N[x^3]", "x^3", 0);
    assert_eval_eq("N[3 x^2 + 3 x + 1]", "1.0 + 3.0 x + 3.0 x^2", 0);
    /* Base is still numericalized. */
    assert_eval_eq("N[(Pi)^2]", "9.8696", 0);
    /* Non-integer exponents are also preserved when the base remains
     * symbolic after numericalization: keeping Sqrt[x] in Sqrt[x] form
     * (rather than x^0.5) preserves polynomial structure for downstream
     * Together / Cancel and matches the contagion fix in numericalize_function. */
    assert_eval_eq("N[x^(1/2)]", "Sqrt[x]", 0);
}

static void test_n_two_arg_fallback(void) {
    /* In a USE_MPFR=0 build, N[expr, n] emits a one-shot stderr warning
     * and returns a machine-precision double. With USE_MPFR enabled we
     * get arbitrary precision — that path is exercised in the Phase-2
     * tests below. */
#ifndef USE_MPFR
    assert_eval_eq("N[Pi, 20]", "3.14159", 0);
#endif
    assert_eval_eq("N[Pi, MachinePrecision]", "3.14159", 0);
}

static void test_n_bad_precision_arg(void) {
    /* A symbolic (non-numeric, non-MachinePrecision) precision argument
     * keeps the expression unevaluated. */
    assert_eval_eq("N[Pi, x]", "N[Pi, x]", 0);
}

static void test_machineprecision_is_protected(void) {
    /* Attempting to overwrite should leave the symbol untouched. */
    struct Expr* p = parse_expression("MachinePrecision = 99");
    struct Expr* r = evaluate(p);
    expr_free(p); expr_free(r);

    assert_eval_eq("MachinePrecision", "MachinePrecision", 0);
}

/* ------------------------------------------------------------------------
 *  Phase 2: arbitrary-precision N[], precision literals, Precision/Accuracy
 * ---------------------------------------------------------------------- */
static void test_inexact_contagion(void) {
    /* Inexact contagion in Plus and Times: when any summand/factor is a
     * Real, exact numeric constants (Pi, E, Sqrt[2], …) are numericalized
     * in-place instead of staying as frozen Times/Plus expressions. */
    assert_eval_eq("1. Pi",       "3.14159", 0);
    assert_eval_eq("1. + Pi",     "4.14159", 0);
    assert_eval_eq("1. E",        "2.71828", 0);
    assert_eval_eq("1. Sqrt[2]",  "1.41421", 0);
    assert_eval_eq("1. Sin[1]",   "0.841471", 0);
    assert_eval_eq("1. x Pi",     "3.14159 x", 0);
    /* No contagion when every arg is exact. */
    assert_eval_eq("2 Pi",   "2 Pi",   0);
    assert_eval_eq("1 + Pi", "1 + Pi", 0);
    /* Symbolic factors stay symbolic — Pi was the only exact numeric. */
    assert_eval_eq("1. x", "1.0 x", 0);

#ifdef USE_MPFR
    /* Precision contagion — the *lowest* precision among inexact args
     * wins. MachinePrecision is the floor; any Real collapses the
     * combination to machine even alongside higher-precision MPFR. */
    assert_eval_eq("1.0`50 Pi + 1.", "4.14159", 0);
    assert_eval_eq("Precision[1.0`50 Pi + 1.]", "MachinePrecision", 0);
    assert_eval_eq("1.0`50 + 1.", "2.0", 0);
    assert_eval_eq("Precision[1.0`50 + 1.]", "MachinePrecision", 0);

    /* Two MPFR operands — min precision wins, not max. 1.0`50 + 1.0`20
     * used to preserve 50 digits; now lands at 20. */
    assert_eval_startswith("1.0`50 Pi + 1.0`20", "4.14159265358979323846");
    assert_eval_startswith("Precision[1.0`50 Pi + 1.0`20]", "20.");
    assert_eval_startswith("Precision[1.0`30 1.0`50]",     "30.");

    /* Single MPFR (or multiple at the same precision) — no lower peer,
     * precision is preserved. */
    assert_eval_startswith("1.0`50 Pi",
                           "3.14159265358979323846264338327950288419");
    assert_eval_startswith("Precision[1.0`50 Pi]", "50.");
#endif
}

#ifdef USE_MPFR

static void test_n_prec_constants(void) {
    /* Mathematica-pinned values for Pi, E, Sqrt[2] at 40 digits. We use
     * prefix-match to tolerate last-bit rounding. */
    assert_eval_startswith("N[Pi, 40]", "3.141592653589793238462643383279502884");
    assert_eval_startswith("N[E, 40]",  "2.718281828459045235360287471352662497");
    assert_eval_startswith("N[Sqrt[2], 40]", "1.41421356237309504880168872420969807");
}

static void test_n_prec_transcendentals(void) {
    assert_eval_startswith("N[Sin[1], 40]",
                           "0.8414709848078965066525023216302989");
    assert_eval_startswith("N[Cos[1], 40]",
                           "0.5403023058681397174009366074429766");
    assert_eval_startswith("N[Log[2], 40]",
                           "0.6931471805599453094172321214581765");
    assert_eval_startswith("N[Exp[1], 40]",
                           "2.718281828459045235360287471352662");
    assert_eval_startswith("N[Sinh[1], 40]",
                           "1.175201193643801456882381850595600");
    assert_eval_startswith("N[Tanh[1], 40]",
                           "0.761594155955764888119458282604793");
    assert_eval_startswith("N[ArcTan[1], 40]",
                           "0.7853981633974483096156608458198757");
}

static void test_n_prec_rational_and_mixed(void) {
    assert_eval_startswith("N[1/7, 40]",
                           "0.1428571428571428571428571428571428571428");
    assert_eval_startswith("N[Pi + E, 40]",
                           "5.85987448204883");
    assert_eval_startswith("N[Pi^2, 30]",
                           "9.86960440108935861883449099987");
    assert_eval_startswith("N[Sqrt[Pi] + Log[2], 30]",
                           "2.46560103146546");
}

static void test_precision_literal_basic(void) {
    /* Precision literals parse and propagate through N-style walks.
     * We avoid asserting the exact printed form because the closest
     * representable MPFR value for a decimal like 3.98 shows last-digit
     * noise when rendered at its full precision — a known Phase-2 v1
     * limitation (Mathematica pads "declared" trailing zeros instead).
     * Instead we verify the Precision is what we asked for. */
    assert_eval_startswith("Precision[3.98`50]", "50.");
    assert_eval_startswith("Precision[1.5`30]",  "30.");
    assert_eval_startswith("Precision[3`30]",    "30.");
}

static void test_precision_builtin(void) {
    assert_eval_eq("Precision[3]",         "Infinity", 0);
    assert_eval_eq("Precision[1/3]",       "Infinity", 0);
    assert_eval_eq("Precision[Pi]",        "Infinity", 0);
    assert_eval_eq("Precision[3.14]",      "MachinePrecision", 0);
    /* Precision of a 50-digit literal: close to 50, not exactly (bit
     * precision rounds up). Check with a prefix that starts with "50.". */
    assert_eval_startswith("Precision[3.14`50]", "50.");
}

static void test_accuracy_builtin(void) {
    assert_eval_eq("Accuracy[0]",     "Infinity", 0);
    assert_eval_eq("Accuracy[0.0]",   "Infinity", 0);
    assert_eval_eq("Accuracy[1/2]",   "Infinity", 0);
    /* Accuracy of 0.001 ≈ 15.95 + 3 = 18.95. */
    assert_eval_startswith("Accuracy[0.001]", "18.");
    /* Accuracy of 3.14`50 ≈ 50 − log10(3.14) ≈ 49.5. */
    assert_eval_startswith("Accuracy[3.14`50]", "49.");
}

static void test_set_precision(void) {
    assert_eval_startswith("SetPrecision[Pi, 40]",
                           "3.141592653589793238462643383279502884");
    assert_eval_startswith("SetPrecision[1/3, 20]",
                           "0.33333333333333333333");
    /* Down-cast MPFR to machine. */
    assert_eval_eq("SetPrecision[Pi, MachinePrecision]", "3.14159", 0);
    /* Precision inspection of the result. */
    assert_eval_startswith("Precision[SetPrecision[Pi, 40]]", "40.");
}

static void test_set_accuracy(void) {
    /* SetAccuracy[0.1, 30] produces a value with ~30 digits of accuracy
     * (about 29 digits of precision since log10(0.1) = -1). The printed
     * result starts with 0.1… (further digits are the accumulated
     * rounding of the double 0.1 seed, which is fine — Mathematica does
     * the same). */
    assert_eval_startswith("SetAccuracy[0.1, 30]", "0.1");
}

static void test_listable_prec(void) {
    assert_eval_startswith("N[{Pi, E}, 20]",
                           "{3.141592653589793238");
}

#endif /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_n_integer_leaf);
    TEST(test_n_rational_leaf);
    TEST(test_n_real_leaf);
    TEST(test_n_constants);
    TEST(test_n_unknown_symbols_pass_through);
    TEST(test_n_on_sin_cos);
    TEST(test_n_descend_then_reevaluate);
    TEST(test_n_listable);
    TEST(test_n_complex);
    TEST(test_n_preserves_hold);
    TEST(test_n_nested);
    TEST(test_n_preserves_integer_exponents);
    TEST(test_n_two_arg_fallback);
    TEST(test_n_bad_precision_arg);
    TEST(test_machineprecision_is_protected);
    TEST(test_inexact_contagion);

#ifdef USE_MPFR
    TEST(test_n_prec_constants);
    TEST(test_n_prec_transcendentals);
    TEST(test_n_prec_rational_and_mixed);
    TEST(test_precision_literal_basic);
    TEST(test_precision_builtin);
    TEST(test_accuracy_builtin);
    TEST(test_set_precision);
    TEST(test_set_accuracy);
    TEST(test_listable_prec);
#endif

    printf("All numeric_tests passed.\n");
    return 0;
}
