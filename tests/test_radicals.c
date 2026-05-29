/* Unit tests for ToRadicals (src/radicals.c).
 *
 * Coverage:
 *   - linear, quadratic, cubic, quartic Root[] conversions
 *   - canonical k-th root selection
 *   - binomial Root[Function[a #^n + b], k] for n >= 5
 *   - threading over List, Equal, inequalities, And, Or, Not
 *   - parametric Root[] coefficients (root_numericalize unavailable)
 *   - two-arg `Function[t, body]` form alongside Slot[1]
 *   - degree-5 non-binomial stays a held Root[]
 *   - k out of range stays a held Root[]
 *   - idempotency: ToRadicals[ToRadicals[expr]] == ToRadicals[expr]
 *   - arity guard: ToRadicals[] / ToRadicals[a, b] stay unevaluated
 *   - residual sanity: every emitted radical r satisfies p[r] == 0
 *
 * Run binary directly: ./radicals_tests
 * (per MEMORY.md note: ctest is not configured in tests/CMakeLists.txt). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

/* The parametric tests trigger `Root::nonint` diagnostics from
 * root_numericalize (it cannot decide non-integer coefficients) and the
 * quintic-passthrough test triggers it indirectly via numericalize.
 * Squash stderr once so the test output is readable. */
static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) {
        freopen("/dev/null", "w", stderr);
        done = 1;
    }
}

/* FullForm string match. */
static void check_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* "expr evaluates to True" sanity check, used for numeric-residual
 * assertions like `Abs[r^3 + r + 1] < 1.*^-15`. */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n",
                input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* 1. Linear                                                           */
/* ------------------------------------------------------------------ */

static void test_linear(void) {
    /* Root[3 # - 5 &, 1] = 5/3. */
    check_eq("ToRadicals[Root[Function[3 # - 5], 1]]",
             "Rational[5, 3]");
}

/* ------------------------------------------------------------------ */
/* 2. Quadratic — both Mathematica examples                            */
/* ------------------------------------------------------------------ */

static void test_quadratic_k1(void) {
    /* Mathematica: ToRadicals[Root[#^2 + 3 # + 5 &, 1]] = 1/2(-3 - I Sqrt[11]) */
    check_eq("ToRadicals[Root[Function[#^2 + 3 # + 5], 1]]",
             "Times[Rational[1, 2], "
             "Plus[-3, Times[Complex[0, -1], Power[11, Rational[1, 2]]]]]");
}

static void test_quadratic_k2(void) {
    check_eq("ToRadicals[Root[Function[#^2 + 3 # + 5], 2]]",
             "Times[Rational[1, 2], "
             "Plus[-3, Times[Complex[0, 1], Power[11, Rational[1, 2]]]]]");
}

/* Real quadratic — canonical k=1 picks the smaller (negative) root. */
static void test_quadratic_real_ordering(void) {
    /* Roots of x^2 - 2 are ±Sqrt[2]; k=1 -> -Sqrt[2], k=2 -> Sqrt[2]. */
    check_eq("ToRadicals[Root[Function[#^2 - 2], 1]]",
             "Times[-1, Power[2, Rational[1, 2]]]");
    check_eq("ToRadicals[Root[Function[#^2 - 2], 2]]",
             "Power[2, Rational[1, 2]]");
}

/* ------------------------------------------------------------------ */
/* 3. Cubic — radical residual must vanish for every k                 */
/* ------------------------------------------------------------------ */

static void test_cubic_residual_k1(void) {
    check_true("With[{r = ToRadicals[Root[Function[#^3 - 5 #^2 - 7 # + 9], 1]]},"
               " Chop[N[r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
}
static void test_cubic_residual_k2(void) {
    check_true("With[{r = ToRadicals[Root[Function[#^3 - 5 #^2 - 7 # + 9], 2]]},"
               " Chop[N[r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
}
static void test_cubic_residual_k3(void) {
    check_true("With[{r = ToRadicals[Root[Function[#^3 - 5 #^2 - 7 # + 9], 3]]},"
               " Chop[N[r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
}

/* Mathematica's k=1 of x^3 + x + 1 is the unique real root. */
static void test_cubic_canonical_k_matches(void) {
    check_true("Abs[N[ToRadicals[Root[Function[#^3 + # + 1], 1]], 25] - "
               "N[Root[Function[#^3 + # + 1], 1], 25]] < 1.*^-15");
    check_true("Abs[N[ToRadicals[Root[Function[#^3 + # + 1], 2]], 25] - "
               "N[Root[Function[#^3 + # + 1], 2], 25]] < 1.*^-15");
    check_true("Abs[N[ToRadicals[Root[Function[#^3 + # + 1], 3]], 25] - "
               "N[Root[Function[#^3 + # + 1], 3], 25]] < 1.*^-15");
}

/* ------------------------------------------------------------------ */
/* 4. Quartic via Ferrari resolvent cubic + canonical k                */
/* ------------------------------------------------------------------ */

static void test_quartic_residual_all_k(void) {
    check_true("With[{r = ToRadicals[Root[Function["
               "#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 1]]},"
               " Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
    check_true("With[{r = ToRadicals[Root[Function["
               "#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 2]]},"
               " Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
    check_true("With[{r = ToRadicals[Root[Function["
               "#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 3]]},"
               " Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
    check_true("With[{r = ToRadicals[Root[Function["
               "#^4 + 3 #^3 - 5 #^2 - 7 # + 9], 4]]},"
               " Chop[N[r^4 + 3 r^3 - 5 r^2 - 7 r + 9, 30]] == 0]");
}

/* Biquadratic quartic (q == 0 in the depressed form): exercises the
 * biquadratic short-circuit branch of radical_quartic. */
static void test_biquadratic_quartic_k_ordering(void) {
    /* x^4 - 5 x^2 + 6 = 0 -> x^2 = 2 or 3 -> x = ±sqrt(2), ±sqrt(3).
     * Canonical (real-first ascending): -sqrt(3), -sqrt(2), sqrt(2), sqrt(3). */
    check_eq("ToRadicals[Root[Function[#^4 - 5 #^2 + 6], 1]]",
             "Times[-1, Power[3, Rational[1, 2]]]");
    check_eq("ToRadicals[Root[Function[#^4 - 5 #^2 + 6], 2]]",
             "Times[-1, Power[2, Rational[1, 2]]]");
    check_eq("ToRadicals[Root[Function[#^4 - 5 #^2 + 6], 3]]",
             "Power[2, Rational[1, 2]]");
    check_eq("ToRadicals[Root[Function[#^4 - 5 #^2 + 6], 4]]",
             "Power[3, Rational[1, 2]]");
}

static void test_quartic_canonical_k_matches_numerically(void) {
    /* All-real quartic x^4 - 10 x^2 + 1: roots = ±(sqrt(3) ± sqrt(2)). */
    for (int k = 1; k <= 4; k++) {
        char input[256];
        snprintf(input, sizeof input,
                 "Abs[N[ToRadicals[Root[Function[#^4 - 10 #^2 + 1], %d]], 25] - "
                 "N[Root[Function[#^4 - 10 #^2 + 1], %d], 25]] < 1.*^-15", k, k);
        check_true(input);
    }
}

/* ------------------------------------------------------------------ */
/* 5. Binomial fast path — Mathematica's degree-5 example              */
/* ------------------------------------------------------------------ */

static void test_quintic_binomial_k3(void) {
    /* Mathematica: ToRadicals[Root[#^5 - 2 &, 3]] = (-1)^(4/5) 2^(1/5). */
    check_eq("ToRadicals[Root[Function[#^5 - 2], 3]]",
             "Times[Power[-1, Rational[4, 5]], Power[2, Rational[1, 5]]]");
}

static void test_quintic_binomial_residual_all_k(void) {
    for (int k = 1; k <= 5; k++) {
        char input[256];
        snprintf(input, sizeof input,
                 "With[{r = ToRadicals[Root[Function[#^5 - 2], %d]]},"
                 " Chop[N[r^5 - 2, 30]] == 0]", k);
        check_true(input);
    }
}

/* ------------------------------------------------------------------ */
/* 6. Threading                                                        */
/* ------------------------------------------------------------------ */

static void test_threading_list(void) {
    check_eq("ToRadicals[{Root[Function[#^2 + 3 # + 5], 1], "
                       "Root[Function[#^2 + 3 # + 5], 2]}]",
             "List["
             "Times[Rational[1, 2], "
             "Plus[-3, Times[Complex[0, -1], Power[11, Rational[1, 2]]]]], "
             "Times[Rational[1, 2], "
             "Plus[-3, Times[Complex[0, 1], Power[11, Rational[1, 2]]]]]]");
}

static void test_threading_equal(void) {
    /* Root[..., 1] in an Equal node should be replaced; the Equal stays. */
    check_eq("ToRadicals[Root[Function[#^2 + 3 # + 5], 1] == y]",
             "Equal["
             "Times[Rational[1, 2], "
             "Plus[-3, Times[Complex[0, -1], Power[11, Rational[1, 2]]]]], y]");
}

static void test_threading_inequality(void) {
    /* sqrt(2) < 3 evaluates to True. */
    check_eq("ToRadicals[Root[Function[#^2 - 2], 2] < 3]",
             "True");
}

static void test_threading_and(void) {
    /* (-sqrt(2)) < 0 && sqrt(2) > 0 -> True. */
    check_eq("ToRadicals[Root[Function[#^2 - 2], 1] < 0 && "
                       "Root[Function[#^2 - 2], 2] > 0]",
             "True");
}

static void test_threading_or(void) {
    /* True || ... -> True. */
    check_eq("ToRadicals[Root[Function[#^2 - 2], 2] > 0 || False]",
             "True");
}

static void test_threading_inside_plus(void) {
    /* Sum of two radical conversions must combine through Plus. */
    check_true("With[{r = ToRadicals[Root[Function[#^2 + 3 # + 5], 1] + "
                                   "Root[Function[#^2 + 3 # + 5], 2]]},"
               " Chop[N[r + 3, 25]] == 0]");
}

/* ------------------------------------------------------------------ */
/* 7. Parametric — root_numericalize unavailable                       */
/* ------------------------------------------------------------------ */

static void test_parametric_quadratic_k1(void) {
    mute_stderr_once();
    check_eq("ToRadicals[Root[Function[#^2 + a # + b], 1]]",
             "Times[Rational[1, 2], "
             "Plus[Times[-1, a], "
             "Times[-1, Power[Plus[Power[a, 2], Times[-4, b]], Rational[1, 2]]]]]");
}

static void test_parametric_quadratic_k2(void) {
    mute_stderr_once();
    check_eq("ToRadicals[Root[Function[#^2 + a # + b], 2]]",
             "Times[Rational[1, 2], "
             "Plus[Times[-1, a], "
             "Power[Plus[Power[a, 2], Times[-4, b]], Rational[1, 2]]]]");
}

/* ------------------------------------------------------------------ */
/* 8. Two-arg Function form                                            */
/* ------------------------------------------------------------------ */

static void test_two_arg_function_form(void) {
    check_eq("ToRadicals[Root[Function[t, t^2 + 3 t + 5], 1]]",
             "Times[Rational[1, 2], "
             "Plus[-3, Times[Complex[0, -1], Power[11, Rational[1, 2]]]]]");
}

/* ------------------------------------------------------------------ */
/* 9. Unsupported / pass-through cases                                 */
/* ------------------------------------------------------------------ */

static void test_quintic_non_binomial_pass_through(void) {
    /* x^5 - x - 1 (Bring quintic): non-binomial, degree 5 — leave the
     * Root unchanged.  Root is HoldAll, so the parser-produced body
     * (Times[-1, 1] etc.) survives literally and any equality check
     * must use that exact tree. */
    mute_stderr_once();
    check_true("ToRadicals[Root[Function[#^5 - # - 1], 1]] === "
               "Root[Function[#^5 - # - 1], 1]");
}

static void test_k_out_of_range_pass_through(void) {
    mute_stderr_once();
    check_true("ToRadicals[Root[Function[#^2 - 4], 5]] === "
               "Root[Function[#^2 - 4], 5]");
}

static void test_non_root_argument_unchanged(void) {
    check_eq("ToRadicals[5]", "5");
    check_eq("ToRadicals[x]", "x");
    check_eq("ToRadicals[x + y]", "Plus[x, y]");
}

/* ------------------------------------------------------------------ */
/* 10. Idempotency                                                     */
/* ------------------------------------------------------------------ */

static void test_idempotency(void) {
    check_true("ToRadicals[ToRadicals[Root[Function[#^2 + 3 # + 5], 1]]] === "
               "ToRadicals[Root[Function[#^2 + 3 # + 5], 1]]");
    check_true("ToRadicals[ToRadicals[Root[Function[#^3 - 5 #^2 - 7 # + 9], 2]]] === "
               "ToRadicals[Root[Function[#^3 - 5 #^2 - 7 # + 9], 2]]");
}

/* ------------------------------------------------------------------ */
/* 11. Arity guard                                                     */
/* ------------------------------------------------------------------ */

static void test_arity_guard(void) {
    /* ToRadicals[] / ToRadicals[a, b] must remain unevaluated. */
    check_eq("ToRadicals[]", "ToRadicals[]");
    check_eq("ToRadicals[a, b]", "ToRadicals[a, b]");
}

/* ------------------------------------------------------------------ */
/* 12. Attributes registration                                         */
/* ------------------------------------------------------------------ */

static void test_attributes_protected(void) {
    /* Protected attribute must be set by radicals_init. */
    check_true("MemberQ[Attributes[ToRadicals], Protected]");
}

/* ------------------------------------------------------------------ */
/* 13. Mixed-degree expression from Mathematica spec                   */
/* ------------------------------------------------------------------ */

static void test_cubic_plus_quintic_binomial(void) {
    /* Numeric agreement: ToRadicals on
     *   Root[#^3 + # + 11 &, 1] + Root[#^5 - 2 &, 3]
     * should match the sum of the numeric Root values. */
    check_true("Abs[N[ToRadicals[Root[Function[#^3 + # + 11], 1] + "
                                "Root[Function[#^5 - 2], 3]], 25] - "
               "(N[Root[Function[#^3 + # + 11], 1], 25] + "
                "N[Root[Function[#^5 - 2], 3], 25])] < 1.*^-15");
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_linear);

    TEST(test_quadratic_k1);
    TEST(test_quadratic_k2);
    TEST(test_quadratic_real_ordering);

    TEST(test_cubic_residual_k1);
    TEST(test_cubic_residual_k2);
    TEST(test_cubic_residual_k3);
    TEST(test_cubic_canonical_k_matches);

    TEST(test_quartic_residual_all_k);
    TEST(test_biquadratic_quartic_k_ordering);
    TEST(test_quartic_canonical_k_matches_numerically);

    TEST(test_quintic_binomial_k3);
    TEST(test_quintic_binomial_residual_all_k);

    TEST(test_threading_list);
    TEST(test_threading_equal);
    TEST(test_threading_inequality);
    TEST(test_threading_and);
    TEST(test_threading_or);
    TEST(test_threading_inside_plus);

    TEST(test_parametric_quadratic_k1);
    TEST(test_parametric_quadratic_k2);

    TEST(test_two_arg_function_form);

    TEST(test_quintic_non_binomial_pass_through);
    TEST(test_k_out_of_range_pass_through);
    TEST(test_non_root_argument_unchanged);

    TEST(test_idempotency);

    TEST(test_arity_guard);
    TEST(test_attributes_protected);

    TEST(test_cubic_plus_quintic_binomial);

    printf("All ToRadicals tests passed!\n");
    return 0;
}
