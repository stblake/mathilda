/* Tests for N[] and related numeric-evaluation builtins.
 *
 * Phase 1 here covers the machine-precision path. Phase 2 tests
 * (MPFR, precision literals, Precision/Accuracy/SetPrecision/SetAccuracy)
 * are added in a follow-up section once MPFR support lands.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"
#include "message.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

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

/* Regression: a rational whose numerator or denominator overflows int64 into
 * a BigInt must still numericalize to a single real leaf, not a frozen
 * Rational[Real, Real] (e.g. the old bug N[1/10^30] -> 1.0/1.0e+30). */
static void test_n_big_rational_leaf(void) {
    assert_eval_eq("N[1/10^30]", "1e-30", 0);
    assert_eval_eq("N[(-3)/10^30]", "-3e-30", 0);
    assert_eval_eq("Head[N[1/10^30]]", "Real", 0);
    /* Out-of-range machine quotient: carried as a machine-precision MPFR
     * real (arbitrary exponent, ~16 sig figs) so it stays finite and
     * nonzero, not Rational[Real, Real]. */
    assert_eval_startswith("N[10^400/3]", "3.33333333333333");
    assert_eval_startswith("N[1/10^400]", "1.00000000000000e-400");
    /* In-range BigInt-component quotient: ordinary IEEE double (~6 sig figs). */
    assert_eval_startswith("N[(7*10^350)/(3*10^300)]", "2.33333e+50");
    /* Two-argument N computes the quotient at the requested precision. */
    assert_eval_startswith("N[(7*10^350)/(3*10^300), 30]",
                           "2.33333333333333333333333333333");
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
    /* GoldenAngle = (3 - Sqrt[5]) Pi; Glaisher A; Khinchin K. */
    assert_eval_eq("N[GoldenAngle]", "2.39996",   0);
    assert_eval_eq("N[Glaisher]",    "1.28243",   0);
    assert_eval_eq("N[Khinchin]",    "2.68545",   0);
}

/* GoldenAngle / Glaisher / Khinchin are first-class constants: Constant +
 * Protected, numeric to NumericQ, and constant to D. */
static void test_new_constants_identity(void) {
    assert_eval_eq("Attributes[GoldenAngle]", "{Constant, Protected}", 0);
    assert_eval_eq("Attributes[Glaisher]",    "{Constant, Protected}", 0);
    assert_eval_eq("Attributes[Khinchin]",    "{Constant, Protected}", 0);
    assert_eval_eq("NumericQ[GoldenAngle]", "True", 0);
    assert_eval_eq("NumericQ[Glaisher]",    "True", 0);
    assert_eval_eq("NumericQ[Khinchin]",    "True", 0);
    assert_eval_eq("D[GoldenAngle x, x]", "GoldenAngle", 0);
    assert_eval_eq("D[Glaisher x, x]",    "Glaisher",    0);
    assert_eval_eq("D[Khinchin x, x]",    "Khinchin",    0);
}

static void test_n_unknown_symbols_pass_through(void) {
    assert_eval_eq("N[x]",         "x", 0);
    assert_eval_eq("N[myFancy]",   "myFancy", 0);
    assert_eval_eq("N[Infinity]",  "Infinity", 0);
}

/* The working-precision planner's machine-integer Rational scan (no GMP):
 * a non-dyadic rational's magnitude still raises the evaluation (Sin of a
 * 15-digit argument must not inherit the argument's rounding error), a
 * dyadic one does not need to, and an inexact integer leaf still counts.
 * Expected values: 60-digit evaluations by an independent library. */
static void test_n_plan_machine_rationals(void) {
    assert_eval_eq("Abs[N[Sin[(10^15 + 1)/3]] - 0.01063131189606997605411337464807559990955] < 10^-17",
                   "True", 0);
    assert_eval_eq("Abs[N[Sin[2^60 + 1]] + 0.91732943534747918303971392200955657824] < 10^-15",
                   "True", 0);
    assert_eval_eq("N[Sin[3/4]] == Sin[0.75]", "True", 0);
    assert_eval_eq("N[{1/2, 1/3, -5/8, 7}]", "{0.5, 0.333333, -0.625, 7.0}", 0);
    assert_eval_eq("N[{1/3, 1/7}, 30] == {N[1/3, 30], N[1/7, 30]}", "True", 0);
}

/* N[E^y, p] goes straight to exp(y): values against 60-digit references,
 * precisions as for Power[e, y]. */
static void test_n_exp_power(void) {
    assert_eval_eq("Abs[N[E^(1/3), 50] - "
                   "1.39561242508608952862812531960258683759790651519940698261752`60] < 10^-48",
                   "True", 0);
    assert_eval_eq("Abs[N[E^(-81/25), 50] - "
                   "0.0391638950989870737397710903657551404259278563463088657555005`60] < 10^-50",
                   "True", 0);
    assert_eval_eq("Precision[N[E^(1/3), 30]]", "30.0", 0);
    assert_eval_eq("Abs[N[E^(-81/25)] - 0.039163895098987073739771] < 10^-17", "True", 0);
    assert_eval_eq("N[E^(1/2), 40] == N[Sqrt[E], 40]", "True", 0);
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

    /* Contagion must NOT thread N into the arguments of a non-numeric head:
     * the `1` in x[1] is an argument to x, not an arithmetic operand, so
     * `0.333 x[1]` stays `0.333 x[1]`, never `0.333 x[1.]`. Left unfixed this
     * leaves an unbound x[1.] and silently breaks NMinimize's indexed-variable
     * substitution (objective scored non-numeric → infeasible 1e+300). */
    assert_eval_eq("0.333 x[1]",   "Times[0.333, x[1]]",   1);
    assert_eval_eq("x[1] + 0.5",   "Plus[0.5, x[1]]",      1);
    assert_eval_eq("1. Log[x[2]]", "Times[1.0, Log[x[2]]]", 1);
    assert_eval_eq("1. Abs[x[3]]", "Times[1.0, Abs[x[3]]]", 1);
    assert_eval_eq("0.5 (x[1] + x[2]^2)",
                   "Times[0.5, Plus[x[1], Power[x[2], 2]]]", 1);
    /* Only the NUMERIC summands and factors take the inexact precision, as in
     * the reference (checked): a non-numeric term keeps its exact
     * coefficient, so 1. - x is not 1. - 1. x, 1. + x - x is 1., and a
     * non-numeric factor such as (2 + a) is not touched. */
    assert_eval_eq("1. + 2 x",     "Plus[1.0, Times[2, x]]", 1);
    assert_eval_eq("{1. - x, 1. + x - x, 1.5 (2 + a), 1. + Sqrt[2] x, 2.5 + Pi + x/3}",
                   "List[Plus[1.0, Times[-1, x]], 1.0, Times[1.5, Plus[2, a]], "
                   "Plus[1.0, Times[Power[2, Rational[1, 2]], x]], "
                   "Plus[5.64159, Times[Rational[1, 3], x]]]", 1);
    /* Explicit N[] is different from contagion — it threads everywhere, so
     * N[x[1]] == x[1.] exactly as in Mathematica. */
    assert_eval_eq("N[x[1]]",      "x[1.0]", 1);

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
    assert_eval_startswith("1.0`50 Pi + 1.0`20", "4.1415926535897932385");
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

static void test_n_bigint_machine_overflow(void) {
    /* Machine-mode N of an exact integer whose magnitude exceeds DBL_MAX
     * (~1.8e308) used to overflow the double conversion and print "inf.0".
     * It must instead fall back to a machine-precision MPFR real with an
     * arbitrary exponent. (Pinned against Mathematica's N[1001!].) */
    assert_eval_startswith("N[1001!]", "4.027896473371");
    assert_eval_startswith("N[1000!]", "4.023872600770");
    assert_eval_startswith("N[2^2000]", "1.148130695274");
    /* In-range integers must still produce ordinary machine doubles. */
    assert_eval_eq("N[12345]", "12345.0", 0);
    /* The arg here is already an out-of-range MPFR (Plus promotes 1.5 + 1001!
     * past DBL_MAX); N's MPFR->machine down-conversion must not re-overflow. */
    assert_eval_startswith("N[1.5 + 1001!]", "4.027896473371");
    assert_eval_startswith("N[1.5 * 1001!]", "6.041844710057");
}

static void test_n_bare_targets_machine_precision(void) {
    /* Bare N[expr] targets machine precision even when expr is already an
     * approximate arbitrary-precision number (matching Mathematica).
     * Regression: N[Pi, 100] // N (== N[N[Pi, 100]]) used to keep 100 digits
     * instead of collapsing to machine precision. */
    assert_eval_eq("N[N[Pi, 100]]", "3.14159", 0);
    assert_eval_eq("Precision[N[N[Pi, 100]]]", "MachinePrecision", 0);
    /* Other constants and a bare high-precision literal collapse too. */
    assert_eval_eq("Precision[N[N[E, 50]]]", "MachinePrecision", 0);
    assert_eval_eq("Precision[N[2.5`100]]", "MachinePrecision", 0);
    /* Exact symbols numericalize to machine precision. */
    assert_eval_eq("N[Pi]", "3.14159", 0);
    assert_eval_eq("Precision[N[Pi]]", "MachinePrecision", 0);
    /* The explicit two-argument form is a precision request: it lowers an
     * already-approximate value (N[.., 20]) but never manufactures digits it
     * does not have (N[.., 200] on a 100-digit value stays 100 digits). */
    assert_eval_startswith("Precision[N[N[Pi, 100], 20]]", "20.");
    assert_eval_startswith("Precision[N[N[Pi, 100], 200]]", "100.");
    /* Inexact contagion is unchanged: a machine real wins. */
    assert_eval_eq("Precision[1. + N[Pi, 100]]", "MachinePrecision", 0);
}

static void test_n_prec_constants(void) {
    /* Mathematica-pinned values for Pi, E, Sqrt[2] at 40 digits. We use
     * prefix-match to tolerate last-bit rounding. */
    assert_eval_startswith("N[Pi, 40]", "3.141592653589793238462643383279502884");
    assert_eval_startswith("N[E, 40]",  "2.718281828459045235360287471352662497");
    assert_eval_startswith("N[Sqrt[2], 40]", "1.41421356237309504880168872420969807");
}

static void test_n_prec_new_constants(void) {
    /* GoldenAngle = (3 - Sqrt[5]) Pi; Glaisher (Euler-Maclaurin on zeta'(2))
     * and Khinchin (BBC zeta series) — see src/numeric.c. We request more
     * digits than we assert so the checked prefix sits inside the
     * guaranteed-correct region (the system carries no display guard over
     * ceil(digits*log2 10) bits, so the very last printed digit of N[x, k]
     * is never pinned — true for every constant). References are
     * Mathematica-pinned. */
    assert_eval_startswith("N[GoldenAngle, 60]",
        "2.39996322972865332223155550663361385312499901105811");
    assert_eval_startswith("N[Glaisher, 60]",
        "1.28242712910062263687534256886979172776768892732500");
    assert_eval_startswith("N[Khinchin, 60]",
        "2.68545200106530644530971483548179569382038229399446");
    /* Deep precision exercises series convergence / Euler-Maclaurin truncation:
     * request 130 digits, assert a 110-digit prefix. */
    assert_eval_startswith("N[Glaisher, 130]",
        "1.28242712910062263687534256886979172776768892732500"
        "1192063740021740406308858826461129736491958202374394206461");
    assert_eval_startswith("N[Khinchin, 130]",
        "2.68545200106530644530971483548179569382038229399446"
        "2953051152345557218859537152002801141174931847697995153465");
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
    /* 1/7 to 40 digits rounds up its 40th: ...571428|571... -> ...571429 */
    assert_eval_startswith("N[1/7, 40]",
                           "0.1428571428571428571428571428571428571429");
    assert_eval_startswith("N[Pi + E, 40]",
                           "5.85987448204883");
    /* Pi^2 = 9.86960440108935861883449099987|615...: the 30th digit rounds up */
    assert_eval_startswith("N[Pi^2, 30]",
                           "9.86960440108935861883449099988");
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
    assert_eval_eq("Accuracy[1/2]",   "Infinity", 0);
    /* Inexact zero is finite (matches Mathematica):
     *   Accuracy[0.0] = MachinePrecision - Log10[$MinMachineNumber]
     *                 ~ 15.9546 + 307.6526 ~ 323.607. */
    assert_eval_startswith("Accuracy[0.0]", "323.");
    /* Accuracy of 0.001 ≈ 15.95 + 3 = 18.95. */
    assert_eval_startswith("Accuracy[0.001]", "18.");
    /* Accuracy of 3.14`50 ≈ 50 − log10(3.14) ≈ 49.5. */
    assert_eval_startswith("Accuracy[3.14`50]", "49.");
    /* MPFR zero @ p digits: Accuracy[0``p] = p (N[0, 50] is the exact 0) */
    assert_eval_startswith("Accuracy[SetPrecision[0, 50]]", "50.");
    assert_eval_eq("Accuracy[N[0, 50]]", "Infinity", 0);
    /* An exact cancellation x - x keeps the accuracy the operands had,
     * p - log2|x| bits rounded down: N[x, 30] is 100 bits, so 1/34 leaves
     * floor(100 + log2 34) = 105 bits (31.6 digits), 10^10 leaves
     * floor(100 - log2 10^10) = 66 bits (19.9 digits), and 1 its 100 bits.
     * It used to be the bare 100 bits whatever the operands' size. */
    assert_eval_startswith("Accuracy[N[1/34, 30] - N[1/34, 30]]", "31.6");
    assert_eval_startswith("Accuracy[N[10^10, 30] - N[10^10, 30]]", "19.8");
    assert_eval_startswith("Accuracy[N[1, 30] - N[1, 30]]", "30.0");
    assert_eval_startswith("Accuracy[N[1/34, 30] - N[1/34, 40]]", "31.6");
    /* N[x, p] does not change an approximate zero: it has no precision to
     * lower, only its accuracy (was cut to p digits, so a certified tiny
     * difference became a zero known only to 10^-p) */
    assert_eval_startswith("Accuracy[N[SetPrecision[0, 80], 20]]", "80.");
    assert_eval_startswith("Accuracy[N[SetPrecision[0, 10], 20]]", "10.");
    assert_eval_eq("N[N[1/3, 50] - N[1/3, 50], 20] < 10^-45", "True", 0);
    /* a zero known to within 2^-a, times y, is known to within |y| 2^-a:
     * accuracy a - log10|y| (rounded down to bits); 0/y gains log10|y| */
    assert_eval_startswith("Accuracy[SetPrecision[0, 80] N[10^50, 100]]", "29.");
    assert_eval_startswith("Accuracy[SetPrecision[0, 80] 10^50]", "29.");
    assert_eval_startswith("Accuracy[SetPrecision[0, 80]/10^35]", "114.");
    assert_eval_eq("SetPrecision[0, 80] 10^50 < 10^-40", "False", 0);
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

static void test_n_prec_complex_pow(void) {
    /* Regression: at MPFR precision the cpow fallback in power.c used to
     * coerce MPFR operands to zero, producing nan.0 + nan.0*I for any
     * complex-valued power. The polar-form helper in numeric.c now handles
     * these cases. The two forms below evaluate to the same primitive 8th
     * root of unity, so their difference must round to zero. */
    /* Sqrt[2]/2 = 0.70710678118654752440084436210484|90...: 32 digits round up */
    assert_eval_startswith("N[(-1)^(1/4), 32]",
                           "0.70710678118654752440084436210485");
    assert_eval_startswith("N[E^(I Pi/4), 32]",
                           "0.70710678118654752440084436210485");
    assert_eval_eq("N[E^(I Pi/4) - (-1)^(1/4), 32]",
                   "0.e-32 + 0.e-32*I", 0);
    /* Negative real base with complex exponent: i^i = exp(-pi/2). */
    assert_eval_startswith("N[I^I, 32]",
                           "0.20787957635076190854695561983");
}

/* BarnesG at arbitrary precision goes through arb: the 12-term asymptotic
 * needed ~10^(D/24) Gamma factors and overflowed at 80 digits. */
static void test_barnesg_arbitrary_precision(void) {
    assert_eval_startswith("N[BarnesG[37/10], 80]",
        "1.47004087379182024528585019115859644390187499791247432441168975852307655158787");
    assert_eval_eq("N[BarnesG[37/10], 75] - N[BarnesG[37/10], 200] == 0", "True", 0);
    assert_eval_eq("N[BarnesG[x + 2] - x Gamma[x]^2 BarnesG[x] /. x -> 17/10, 30] == 0", "True", 0);
}

/* N[expr, p] promises p good digits, so it raises the working precision until
 * two evaluations agree (up to $MaxExtraPrecision). One fixed-precision run
 * reported cancellation noise as digits: N[Cos[10^-20] - 1, 20] was an exact
 * 0.0 (it is -5*10^-41) and the exact zero below was 3*10^-52 at p = 32. */
static void test_n_adaptive_precision(void) {
    /* a genuine cancellation keeps its digits */
    /* -5*10^-41 + 4.2*10^-82: 20 digits round to -5.0000000000000000000e-41 */
    assert_eval_startswith("N[Cos[10^-20] - 1, 20]", "-5.0000000000000000000e-41");
    assert_eval_eq("Abs[N[Cos[10^-20] - 1, 32]/(-5*10^-41) - 1] < 10^-30", "True", 0);
    /* an exact zero comes out as zero at every precision */
    assert_eval_eq("N[E^(I Pi/4) - (-1)^(1/4), 50]", "0.e-50 + 0.e-50*I", 0);
    assert_eval_eq("N[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1, 40] == 0", "True", 0);
    /* genuine small values and plain cases are untouched */
    assert_eval_startswith("N[Exp[-200], 32]", "1.3838965267367375306486814569791e-87");
    assert_eval_startswith("N[Pi - 314159265358979323846/10^20, 30]", "2.64338327950288419716939937511e-21");
    /* machine inputs gain no digits under N[_, p], so their noise stays */
    assert_eval_startswith("N[0.1 + 0.2 - 0.3, 32]", "5.55112e-17");
}

/* The adaptive climb certifies a result only when it agrees with one at
 * least 32 bits lower AND the pass's own worst cancellation left p bits, so an
 * exact 0.0 produced by a cancellation never certifies itself, and a short
 * last step to the cap cannot agree on rounding noise. Expected values are
 * closed forms, not this code's output. */
/* The MPFR elementary functions run through arb's ball arithmetic and keep a
 * value only when the whole ball rounds to one number, so the answer must be
 * bit for bit the correctly rounded value MPFR itself computes. */
static void test_mpfr_elementary_bit_identical(void) {
    typedef int (*Op)(mpfr_t, const mpfr_t, mpfr_rnd_t);
    const Op ops[] = { mpfr_sin, mpfr_cos, mpfr_tan, mpfr_sinh, mpfr_cosh, mpfr_tanh,
                       mpfr_atan, mpfr_asin, mpfr_acos, mpfr_atanh, mpfr_asinh, mpfr_acosh,
                       mpfr_cot, mpfr_sec, mpfr_csc, mpfr_coth, mpfr_sech, mpfr_csch, mpfr_exp };
    const char* names[] = { "sin", "cos", "tan", "sinh", "cosh", "tanh", "atan", "asin", "acos",
                            "atanh", "asinh", "acosh", "cot", "sec", "csc", "coth", "sech", "csch",
                            "exp" };
    const long precs[] = { 65, 100, 197, 1000, 3001 };
    /* x = num/den: small, moderate, near 1 (domain edges of asin/atanh), large */
    const long nums[] = { 1, 3, 7, 999999, 1000001, 31, -5, 12345 };
    const long dens[] = { 1000000007, 10, 9, 1000000, 1000000, 3, 4, 100 };
    int checked = 0;
    for (size_t o = 0; o < sizeof(ops) / sizeof(ops[0]); o++)
        for (size_t pi = 0; pi < sizeof(precs) / sizeof(precs[0]); pi++)
            for (size_t xi = 0; xi < sizeof(nums) / sizeof(nums[0]); xi++) {
                mpfr_t x, ref;
                mpfr_init2(x, precs[pi]);
                mpfr_init2(ref, precs[pi]);
                mpfr_set_si(x, nums[xi], MPFR_RNDN);
                mpfr_div_si(x, x, dens[xi], MPFR_RNDN);
                ops[o](ref, x, MPFR_RNDN);
                Expr* in = expr_new_mpfr_copy(x);
                Expr* got = numeric_mpfr_apply_unary(in, 0, ops[o]);
                if (mpfr_nan_p(ref)) {
                    ASSERT(got == NULL);
                } else {
                    ASSERT(got != NULL && got->type == EXPR_MPFR);
                    if (!mpfr_equal_p(got->data.mpfr, ref)) {
                        printf("FAIL: %s at %ld bits of %ld/%ld differs from MPFR\n",
                               names[o], precs[pi], nums[xi], dens[xi]);
                        ASSERT(0);
                    }
                    checked++;
                }
                if (got) expr_free(got);
                expr_free(in);
                mpfr_clear(x);
                mpfr_clear(ref);
            }
    ASSERT(checked > 500);
}

/* An approximate argument cannot equal an exact table value: the inverse
 * functions answer numerically (and the exact inputs still give exact angles). */
static void test_inverse_trig_numeric_arguments(void) {
    assert_eval_eq("ArcSin[1/2] === Pi/6", "True", 0);
    assert_eval_eq("ArcTan[1] === Pi/4", "True", 0);
    assert_eval_eq("ArcCos[-1/2] === 2 Pi/3", "True", 0);
    /* the bare rationals among the table values, and radical ones */
    assert_eval_eq("{ArcSec[2], ArcCsc[-2], ArcSin[-1/2], ArcCos[Sqrt[3]/2], ArcTan[2 - Sqrt[3]]}"
                   " === {Pi/3, -Pi/6, -Pi/6, Pi/6, Pi/12}", "True", 0);
    /* a non-table argument stays as it is */
    assert_eval_eq("ArcSin[1/3] === ArcSin[1/3] && Head[ArcSin[ArcTan[1/10^8]]] === ArcSin", "True", 0);
    assert_eval_eq("Abs[ArcSin[0.5] - Pi/6] < 10^-15", "True", 0);
    assert_eval_eq("Abs[ArcCos[0.5] - Pi/3] < 10^-15", "True", 0);
    assert_eval_eq("Abs[ArcTan[1.] - Pi/4] < 10^-15", "True", 0);
    assert_eval_eq("Abs[ArcSin[N[1/2, 40]] - N[Pi/6, 45]] < 10^-39", "True", 0);
    assert_eval_eq("Abs[ArcCot[N[1, 40]] - N[Pi/4, 45]] < 10^-39", "True", 0);
    assert_eval_eq("Precision[ArcSin[N[1/2, 40]]] == 40", "True", 0);
}

/* N of a composition f[g[x]] repeats g's whole tree in every position of x in
 * f; each distinct subexpression is evaluated once per pass.  Values deep in
 * the cancellation (100 to 450 digits cancel), checked against independent
 * 3000- and 6000-digit evaluations that agree to all digits shown. */
/* Ordinary evaluation of a composition at an arbitrary-precision point: the
 * repeated numeric calls inside one sum share their value (eval.c,
 * NfMemoEntry).  The value must be the one each call gives on its own:
 * checked against an evaluation at 2000 digits (3.0603524556729171455e9,
 * which a second system also gives to 11 digits) and, term by term, against
 * the same terms evaluated as a List, whose arguments share nothing. */
static void test_eval_sum_common_numeric_calls(void) {
    assert_eval_eq(
        "fp7[x_] = Tan[Tanh[Sin[x]]] + Tanh[Sin[Tan[x]]] + Sin[Tan[Tanh[x]]] - Tan[Sin[Tanh[x]]] "
        "- Sin[Tanh[Tan[x]]] - Tanh[Tan[Sin[x]]] - Tan[Sinh[Tanh[x]]] - Sinh[Tanh[Tan[x]]] "
        "- Tanh[Tan[Sinh[x]]] + Tan[Tanh[Sinh[x]]] + Tanh[Sinh[Tan[x]]] + Sinh[Tan[Tanh[x]]]; "
        "gp7[x_] = Sinh[Tanh[Sin[x]]] + Tanh[Sin[Sinh[x]]] + Sin[Sinh[Tanh[x]]] - Sinh[Sin[Tanh[x]]] "
        "- Sin[Tanh[Sinh[x]]] - Tanh[Sinh[Sin[x]]] - Tan[Sinh[Sin[x]]] - Sinh[Sin[Tan[x]]] "
        "- Sin[Tan[Sinh[x]]] + Tan[Sin[Sinh[x]]] + Sin[Sinh[Tan[x]]] + Sinh[Tan[Sin[x]]]; "
        "Abs[N[fp7[gp7[N[10^-5, 1000]]]/gp7[fp7[N[10^-5, 1000]]], 12] / "
        "3.0603524556729171455`20*^9 - 1] < 10^-11", "True", 0);
    /* the same terms summed by Total over a List (whose arguments share
     * nothing) agree with the sum to its precision */
    assert_eval_eq("Module[{y = N[10^-5, 1000], t}, "
                   "t = {Tan[Tanh[Sin[y]]], Tanh[Sin[Tan[y]]], Sin[Tan[Tanh[y]]], -Tan[Sin[Tanh[y]]], "
                   "-Sin[Tanh[Tan[y]]], -Tanh[Tan[Sin[y]]], -Tan[Sinh[Tanh[y]]], -Sinh[Tanh[Tan[y]]], "
                   "-Tanh[Tan[Sinh[y]]], Tan[Tanh[Sinh[y]]], Tanh[Sinh[Tan[y]]], Sinh[Tan[Tanh[y]]]}; "
                   "Abs[fp7[y] - Total[t]] <= 10^-990 Max[Abs[t]]]", "True", 0);
    /* term by term: each shared call equals the call evaluated alone */
    assert_eval_eq("Module[{y = N[2/7, 300]}, "
                   "Sin[y] + Sinh[Sin[y]] + 2 Sin[y] + Cos[Sinh[Sin[y]]] - "
                   "(3 (Sin[y] + 0) + (Sinh[Sin[y]] + 0) + (Cos[Sinh[Sin[y]]] + 0)) == 0 && "
                   "Precision[Sin[y] + Sin[y]] == 300]", "True", 0);
    /* equal values at different precisions are different calls */
    assert_eval_eq("Module[{s = Sin[N[1/3, 50]] + Sin[N[1/3, 60]]}, "
                   "Precision[s] == 50 && Abs[s - 2 Sin[N[1/3, 70]]] < 10^-48]", "True", 0);
    /* a call that reports a message is not kept: both copies report it */
    unsigned long before = mth_msg_fired_count();
    assert_eval_eq("Module[{z = N[0, 30]}, Power[z, -1] + Power[z, -1]]", "Indeterminate", 0);
    ASSERT(mth_msg_fired_count() - before >= 2);
}

static void test_n_composition_common_subexpressions(void) {
    assert_eval_eq(
        "f2[x_] = Tan[Sin[x]] - Sin[Tan[x]]; g2[x_] = Tan[Tanh[x]] - Tanh[Tan[x]]; "
        "Block[{$MaxExtraPrecision = 1000}, Abs[N[f2[g2[x]]/g2[f2[x]] /. x -> 10^-6, 30] "
        "- 5.618655692684505410761`25] < 10^-20]", "True", 0);
    assert_eval_eq(
        "f3[x_] = ArcTan[ArcTanh[ArcSin[x]]] + ArcTanh[ArcSin[ArcTan[x]]] + ArcSin[ArcTan[ArcTanh[x]]] "
        "- ArcTan[ArcSin[ArcTanh[x]]] - ArcSin[ArcTanh[ArcTan[x]]] - ArcTanh[ArcTan[ArcSin[x]]] "
        "- ArcTan[ArcSinh[ArcTanh[x]]] - ArcSinh[ArcTanh[ArcTan[x]]] - ArcTanh[ArcTan[ArcSinh[x]]] "
        "+ ArcTan[ArcTanh[ArcSinh[x]]] + ArcTanh[ArcSinh[ArcTan[x]]] + ArcSinh[ArcTan[ArcTanh[x]]]; "
        "g3[x_] = ArcSinh[ArcTanh[ArcSin[x]]] + ArcTanh[ArcSin[ArcSinh[x]]] + ArcSin[ArcSinh[ArcTanh[x]]] "
        "- ArcSinh[ArcSin[ArcTanh[x]]] - ArcSin[ArcTanh[ArcSinh[x]]] - ArcTanh[ArcSinh[ArcSin[x]]] "
        "- ArcTan[ArcSinh[ArcSin[x]]] - ArcSinh[ArcSin[ArcTan[x]]] - ArcSin[ArcTan[ArcSinh[x]]] "
        "+ ArcTan[ArcSin[ArcSinh[x]]] + ArcSin[ArcSinh[ArcTan[x]]] + ArcSinh[ArcTan[ArcSin[x]]]; "
        "Block[{$MaxExtraPrecision = 1000}, Abs[N[f3[g3[x]]/g3[f3[x]] /. x -> 10^-4, 30] "
        "- 9.313225746154790314691`25] < 10^-20]", "True", 0);
    assert_eval_eq(
        "f10[x_] = Tan[Sin[Tan[Sin[x]]]] - Sin[Tan[Sin[Tan[x]]]]; g10[x_] = Sinh[Tanh[x]] - Tanh[Sinh[x]]; "
        "Block[{$MaxExtraPrecision = 1000}, Abs[N[f10[g10[x]]/g10[f10[x]] /. x -> 10^-5, 30] "
        "- 0.01562499996571180559310`25] < 10^-22]", "True", 0);
}

static void test_n_adaptive_certification(void) {
    /* E^1000 (E^(a - E^-1000) - E^a) = -E^a (1 - E^-1000/2 + ...): ~1450 bits
     * cancel. Was an exact 0.0 at every $MaxExtraPrecision. */
    assert_eval_eq("Block[{$MaxExtraPrecision = 1000}, Abs[N[E^1000 (E^(1/1000 - E^-1000) - E^(1/1000)), 30]"
                   "/N[-E^(1/1000), 40] - 1] < 10^-29]", "True", 0);
    assert_eval_eq("Block[{$MaxExtraPrecision = 1000}, Abs[N[E^1000 (E^(3/7 - E^-1000) - E^(3/7)), 25]"
                   "/N[-E^(3/7), 40] - 1] < 10^-24]", "True", 0);
    /* beyond the default 50 extra digits: N::meprec, not a silent 0.0 */
    assert_eval_eq("Check[N[E^1000 (E^(1/1000 - E^-1000) - E^(1/1000)), 30], meprecFired]",
                   "meprecFired", 0);
    /* a division by such a difference is not an answer either */
    assert_eval_eq("Check[N[1/(E^(1/1000 - E^-1000) - E^(1/1000)), 30], meprecFired]",
                   "meprecFired", 0);
    /* within the limit: E^(a - E^-100) - E^a = -E^a E^-100 (1 - E^-100/2) */
    assert_eval_eq("Abs[N[E^(1/1000 - E^-100) - E^(1/1000), 30]/N[-E^(1/1000 - 100), 40] - 1] < 10^-29",
                   "True", 0);
    /* nested differences: (E^x - 1)/x - 1 = x/2 + x^2/6 + ..., x = 10^-30,
     * 100 bits cancel twice; the climb's last step is a short one */
    assert_eval_eq("Abs[N[(Exp[10^-30] - 1) 10^30 - 1, 20]/(5*10^-31) - 1] < 10^-19", "True", 0);
    assert_eval_eq("Check[N[(Exp[10^-30] - 1) 10^30 - 1, 20]; ok, meprecFired]", "ok", 0);
    /* (Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8 = -17/2520 - 31/14175 x^2
     * + O(x^4): ~205 bits cancel after Log loses ~67 near 1. Its last step to
     * the cap (330 -> 331 bits) certified 18 wrong digits; now N::meprec at the
     * default limit, and correct with more room. 30 digits of a value of
     * size 6.7*10^-3 know it to within ~6.7*10^-33 (Accuracy 32.17), so the
     * difference from the exact value is a zero of that accuracy and cannot
     * be shown below 10^-33; 10^-31 is the bound 30 digits can certify. */
    assert_eval_eq("Check[N[(Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8 /. x -> 10^-10, 30], meprecFired]",
                   "meprecFired", 0);
    assert_eval_eq("Block[{$MaxExtraPrecision = 100}, Abs[N[(Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8 "
                   "/. x -> 10^-10, 30] - (-17/2520 - 31/14175 10^-20)] < 10^-31]", "True", 0);
    assert_eval_eq("Block[{$MaxExtraPrecision = 100}, Abs[N[(Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8 "
                   "/. x -> 10^-10, 30] - (-17/2520 - 31/14175 10^-20)] < 10^-33]", "False", 0);
    /* b^n with a huge exact n multiplies the base's rounding error by n:
     * (1 + 10^-40)^(10^40) - E = -E/2 10^-40 (1 - 11/12 10^-40 + ...); was
     * 3.5575743793 with no message */
    assert_eval_eq("Block[{$MaxExtraPrecision = 100}, Abs[N[(1 + 10^-40)^(10^40) - E, 10]"
                   "/(-E/2 10^-40) - 1] < 10^-9]", "True", 0);
    assert_eval_eq("Check[N[(1 + 10^-40)^(10^40) - E, 10], meprecFired]", "meprecFired", 0);
    assert_eval_eq("Abs[N[(1 + 10^-40)^(10^40), 30]/N[E, 40] - 1] < 10^-29", "True", 0);
    /* a large argument is not a cancellation: Sin[10^30] to 20 digits */
    assert_eval_startswith("N[Sin[10^30], 20]", "-0.0901169019121380580");
    /* an exact zero is 0 at the precision reached; it cannot be told from
     * 0 times a huge factor, so N says it is not certified */
    assert_eval_eq("Quiet[N[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1, 40]] == 0", "True", 0);
    assert_eval_eq("Check[N[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1, 40], meprecFired]", "meprecFired", 0);
    /* a part that is not cancelling is unaffected by one that is */
    assert_eval_eq("Abs[Quiet[N[Pi + (Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]), 30]] - N[Pi, 40]] < 10^-29",
                   "True", 0);
}

/* An exact zero has no digits, so N[_, p] leaves it exact, and an exact-zero
 * part stays exact through arbitrary-precision arithmetic and the kernels
 * that map an axis to an axis. Machine N keeps 0. Values are closed forms. */
static void test_n_exact_zero_parts(void) {
    assert_eval_eq("N[0, 30] === 0", "True", 0);
    assert_eval_eq("N[{0, 1/3, 0}, 20][[{1, 3}]] === {0, 0}", "True", 0);
    assert_eval_eq("N[f[0, 1/2], 20][[1]] === 0", "True", 0);
    assert_eval_eq("0 N[Pi, 20] === 0", "True", 0);
    assert_eval_eq("Times[0, N[Pi, 20], x] === 0", "True", 0);
    /* pure imaginary results: the real part is exactly 0 */
    assert_eval_eq("Re[N[2 I, 20]] === 0 && Re[N[I Pi, 20]] === 0", "True", 0);
    assert_eval_eq("Re[N[Sqrt[-2], 20]] === 0 && Abs[Im[N[Sqrt[-2], 20]] - N[Sqrt[2], 40]] < 10^-19",
                   "True", 0);
    assert_eval_eq("With[{z = Sqrt[N[-3, 25]]^3}, Re[z] === 0 && Abs[Im[z] + N[3 Sqrt[3], 40]] < 10^-23]",
                   "True", 0);
    assert_eval_eq("Complex[0, N[2, 20]]^2 === N[-4, 20]", "True", 0);
    assert_eval_eq("Complex[0, N[2, 20]] Complex[0, N[3, 20]] === N[-6, 20]", "True", 0);
    /* Erf[I y] = I Erfi[y], Erfi[I y] = I Erf[y] */
    assert_eval_eq("With[{z = N[Erf[I], 30]}, Re[z] === 0 && Abs[Im[z] - N[Erfi[1], 40]] < 10^-29]",
                   "True", 0);
    assert_eval_eq("With[{z = N[Erfi[I], 25]}, Re[z] === 0 && Abs[Im[z] - N[Erf[1], 40]] < 10^-24]",
                   "True", 0);
    /* ArcCosh[x] = I ArcCos[x] on [-1, 1); ArcCos[x] = I ArcCosh[x] for x > 1;
     * ArcCosh[1/2] had a stray real part -3.9*10^-90 at 70 digits */
    assert_eval_eq("With[{z = ArcCosh[N[1/2, 50]]}, Re[z] === 0 && Abs[Im[z] - N[Pi/3, 60]] < 10^-49]",
                   "True", 0);
    assert_eval_eq("Re[N[ArcCosh[1/2], 70]] === 0 && Re[N[ArcCosh[-1/3], 30]] === 0", "True", 0);
    assert_eval_eq("With[{z = N[ArcCos[2], 20]}, Re[z] === 0 && Abs[Im[z] - N[ArcCosh[2], 40]] < 10^-19]",
                   "True", 0);
    /* Sin[I y] = I Sinh[y], Tan[I y] = I Tanh[y] and back */
    assert_eval_eq("With[{z = Sin[Complex[0, N[2, 20]]]}, Re[z] === 0 && Abs[Im[z] - N[Sinh[2], 40]] < 10^-19]",
                   "True", 0);
    assert_eval_eq("Re[Tanh[Complex[0, N[1/3, 20]]]] === 0 && Re[Sinh[Complex[0, N[1, 20]]]] === 0",
                   "True", 0);
    /* a zero matrix stays exact, so its decomposition is the exact one */
    assert_eval_eq("QRDecomposition[N[{{0, 0}, {0, 0}}, 60]] === {{{1, 0}}, {{0, 0}}}", "True", 0);
    /* an exact 0 index or branch is now honoured: ProductLog[0, z] */
    assert_eval_eq("Abs[N[ProductLog[0, 1/3], 20] - 0.25762765304973670428`22] < 10^-19", "True", 0);
    /* machine precision is unchanged */
    assert_eval_eq("N[0]", "0.0", 0);
    assert_eval_eq("0 2.5", "0.0", 0);
    assert_eval_eq("Head[Re[N[ArcCosh[1/2]]]] === Real && Head[Re[N[Erf[I]]]] === Real", "True", 0);
}

/* An N[expr, p] that cancels nothing is evaluated exactly once. */
static void test_n_adaptive_single_pass(void) {
    const char* cases[] = {
        "N[Pi + E, 30]", "N[Sqrt[2] - 1, 50]", "N[Sin[1] + Cos[1] - 1, 30]",
        "N[Exp[-200] + Log[3], 40]", "N[{1/3 + Pi, 2 - E}, 25]", "N[Sin[10^30], 20]",
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        long before = numeric_adaptive_reevaluations();
        Expr* in = parse_expression(cases[i]);
        Expr* r = evaluate(in);
        expr_free(in);
        expr_free(r);
        long after = numeric_adaptive_reevaluations();
        if (after != before) {
            fprintf(stderr, "  %s re-evaluated %ld times\n", cases[i], after - before);
        }
        ASSERT(after == before);
    }
    /* and a cancelling one does re-evaluate */
    long before = numeric_adaptive_reevaluations();
    Expr* in = parse_expression("N[Cos[10^-20] - 1, 20]");
    Expr* r = evaluate(in);
    expr_free(in);
    expr_free(r);
    ASSERT(numeric_adaptive_reevaluations() > before);
}

/* Phase 8: end-to-end MPFR-complex corner cases. The infrastructure
 * landed in earlier phases (per-builtin wirings + numeric_mpfr_complex_pow);
 * this test fixes the behaviour at the REPL level so future regressions
 * surface here rather than in scattered builtin tests. */
static void test_mpfr_complex_corner_cases(void) {
    /* Sqrt of a negative MPFR real → Complex[0, MPFR Sqrt[|x|]].
     * The Power-with-half-exponent path routes through
     * numeric_mpfr_complex_pow's polar form. */
    assert_eval_startswith("Im[Sqrt[N[-2, 50]]]",
                           "1.41421356237309504880168872420969807856967187537");
    assert_eval_startswith("Precision[Sqrt[N[-2, 50]]]", "50.");

    /* Complex base raised to an MPFR exponent: (1+i)^3 = -2 + 2i at
     * MPFR precision. */
    assert_eval_startswith("Re[(1 + I)^N[3, 50]]",
                           "-2.0000000000000000000000000000000000000000000000");
    assert_eval_startswith("Im[(1 + I)^N[3, 50]]",
                           "2.0000000000000000000000000000000000000000000000");
    assert_eval_startswith("Precision[(1 + I)^N[3, 50]]", "50.");

    /* Sin at very high precision (100 digits) stays at full precision
     * end-to-end via Phase 4 + the helper module. */
    assert_eval_startswith("Re[Sin[Complex[N[1, 100], N[1, 100]]]]",
                           "1.298457581415977294826042365807815620313436561635");
    assert_eval_startswith("Precision[Sin[Complex[N[1, 100], N[1, 100]]]]",
                           "100.");

    /* Round-trip: Log[Exp[z]] = z at MPFR precision for principal-branch z.
     * Use a non-dyadic real part (1/3) so the round-trip exercises real
     * 50-digit output: a dyadic value like 0.5 round-trips to an exact
     * binary fraction whose trailing zeros the printer trims (%.*Rg), so
     * the printed form would be "0.5" and the prefix check vacuous. The
     * low-order digits carry the expected Exp/Log bit-level rounding, hence
     * a startswith check over the stable leading digits. */
    assert_eval_startswith(
        "Re[Log[Exp[Complex[N[1/3, 50], N[0.25, 50]]]]]",
        "0.33333333333333333333333333333333333333333");
    assert_eval_startswith(
        "Im[Log[Exp[Complex[N[1/3, 50], N[0.25, 50]]]]]",
        "0.25000000000000000000000000000000000000000");
    assert_eval_startswith(
        "Precision[Log[Exp[Complex[N[1/3, 50], N[0.25, 50]]]]]", "50.");
}

/* Guard bits: N[x, d] stores its value with NUMERIC_PREC_GUARD_BITS more
 * mantissa bits than its precision, so the d digits it shows are the true
 * value's correct rounding. Without them the stored ceil(d log2 10)-bit value
 * carried an error of up to 1.5 units of the last digit: N[Sqrt[2], 50] ended
 * ...753770 (Sqrt[2] = ...75376948...), N[FresnelC[1], 40] ...366307
 * (...3663062570...). The expected digits are computed here, independently,
 * by MPFR at 400 bits and rounded once to d significant digits. */
static void guard_expect_sqrt(int k, int d, char* out, size_t cap) {
    mpfr_t r;
    mpfr_init2(r, 400);
    mpfr_set_ui(r, (unsigned long)k, MPFR_RNDN);
    mpfr_sqrt(r, r, MPFR_RNDN);
    /* d significant digits, positional (1 <= Sqrt[k] < 10 here) */
    mpfr_snprintf(out, cap, "%.*Rf", d - 1, r);
    mpfr_clear(r);
}

static void test_n_guard_bits_display(void) {
    assert_eval_eq("N[Sqrt[2], 50]",
                   "1.4142135623730950488016887242096980785696718753769", 0);
    assert_eval_eq("N[FresnelC[1], 40]",
                   "0.7798934003768228294742064136526901366306", 0);
    /* every Sqrt[k] < 100 at 20, 30, 50 and 100 digits (the old storage
     * misrounded 68 of these 360) */
    static const int ds[] = { 20, 30, 50, 100 };
    char want[256], in[64];
    for (int k = 2; k < 100; k++) {
        int r = 1; while (r * r < k) r++;
        if (r * r == k) continue;          /* exact squares print as integers */
        for (size_t j = 0; j < sizeof ds / sizeof ds[0]; j++) {
            guard_expect_sqrt(k, ds[j], want, sizeof want);
            snprintf(in, sizeof in, "N[Sqrt[%d], %d]", k, ds[j]);
            assert_eval_eq(in, want, 0);
        }
    }
    /* the precision is the requested one, not the mantissa width: exactly d
     * for N[x, d], and the min of the operands' through Times */
    assert_eval_eq("Precision[N[Pi, 30]] === 30.", "True", 0);
    assert_eval_eq("Precision[N[Pi, 30] N[E, 40]] === 30.", "True", 0);
    assert_eval_eq("Precision[SetPrecision[2, 25]] === 25.", "True", 0);
    assert_eval_eq("Precision[1.`40] === 40.", "True", 0);
    assert_eval_eq("Precision[N[Exp[1000]]] == $MachinePrecision", "True", 0);
    assert_eval_eq("Precision[N[Sqrt[2], 11]] === 11.", "True", 0);
    /* InputForm re-reads to the identical value and precision */
    assert_eval_eq("With[{v = N[Sqrt[3], 45]}, ToExpression[ToString[v, InputForm]] === v]",
                   "True", 0);
    /* the guard is invisible to Equal: the tolerance is still that of 30
     * digits less the 7 ignored bits (a difference of 10^-31 is equal, one of
     * 10^-26 is not) */
    assert_eval_eq("N[Pi, 30] == N[Pi, 30] + 10^-31", "True", 0);
    assert_eval_eq("N[Pi, 30] == N[Pi, 30] + 10^-26", "False", 0);
}

#endif /* USE_MPFR */

/* N inside a position the head HOLDS numericalizes the structure but never
 * evaluates it: If's branches, a Function body.  N used to rebuild and
 * evaluate every argument, so N[If[c, Print[a], b]] printed and lost the
 * branch -- and Equal, which numericalizes symbolic sides, ran every held
 * side effect of `If[...] == 2`. */
static void test_n_does_not_evaluate_held_arguments(void) {
    assert_eval_eq("ClearAll[nhc, y]; nhc = 0; N[If[y == 2, nhc++, y + 1]]",
                   "If[y == 2.0, Increment[nhc], y + 1.0]", 0);
    assert_eval_eq("nhc", "0", 0);
    assert_eval_eq("nhc = 0; If[y == 2, nhc++, y + 1] == 2; nhc", "0", 0);
    /* 0 -> 1 -> 2 -> 0 (nhc 1) -> 1 -> 2 -> 1 (nhc 2): only real visits count. */
    assert_eval_eq("nhc = 0; FixedPoint[If[# == 2, nhc++, # + 1] &, 0, 6]; nhc", "2", 0);
    assert_eval_eq("N[If[y > 1, Sin[1], Cos[1]]]", "If[y > 1.0, Sin[1.0], Cos[1.0]]", 0);
    assert_eval_eq("N[Function[x, x^2 + 1/3]][2]", "4.33333", 0);
    /* Non-held positions still evaluate: the If condition is numericalized
     * and decided, and an evaluated head computes as before. */
    assert_eval_eq("N[If[True, 1/3, Pi]]", "0.333333", 0);
    assert_eval_eq("N[f[1/3, Sqrt[2]]]", "f[0.333333, 1.41421]", 0);
    assert_eval_eq("ClearAll[nhc]", "Null", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_n_integer_leaf);
    TEST(test_n_rational_leaf);
    TEST(test_n_big_rational_leaf);
    TEST(test_n_real_leaf);
    TEST(test_n_constants);
    TEST(test_new_constants_identity);
    TEST(test_n_unknown_symbols_pass_through);
    TEST(test_n_plan_machine_rationals);
    TEST(test_n_exp_power);
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
    TEST(test_n_bigint_machine_overflow);
    TEST(test_n_bare_targets_machine_precision);
    TEST(test_n_prec_constants);
    TEST(test_n_prec_new_constants);
    TEST(test_n_prec_transcendentals);
    TEST(test_n_prec_rational_and_mixed);
    TEST(test_precision_literal_basic);
    TEST(test_precision_builtin);
    TEST(test_accuracy_builtin);
    TEST(test_set_precision);
    TEST(test_set_accuracy);
    TEST(test_listable_prec);
    TEST(test_n_prec_complex_pow);
    TEST(test_n_adaptive_precision);
    TEST(test_n_adaptive_certification);
    TEST(test_n_adaptive_single_pass);
    TEST(test_n_exact_zero_parts);
    TEST(test_barnesg_arbitrary_precision);
    TEST(test_mpfr_complex_corner_cases);
    TEST(test_n_guard_bits_display);
    TEST(test_n_does_not_evaluate_held_arguments);
    TEST(test_mpfr_elementary_bit_identical);
    TEST(test_inverse_trig_numeric_arguments);
    TEST(test_n_composition_common_subexpressions);
    TEST(test_eval_sum_common_numeric_calls);
#endif

    printf("All numeric_tests passed.\n");
    return 0;
}
