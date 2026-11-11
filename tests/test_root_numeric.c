/* Unit tests for N[Root[..]] / N[Root[..], prec] in src/root_numeric.c.
 *
 * Coverage:
 *   - Wallis cubic (all-real spectrum), k=1, machine + 30 digits.
 *   - Cubic with one real + complex pair, canonical ordering.
 *   - Conjugate-pair tiebreak (smaller Im first).
 *   - All-real quartic across k = 1..4 ascending.
 *   - High-precision Bring quintic at 200 digits.
 *   - Squarefree pre-pass on (x^2 - 2)^3.
 *   - Solve round-trip: N[Solve[..]] keeps canonical k consistency.
 *   - Out-of-range k yields the original Root expression unchanged.
 *
 * Run binary directly: ./root_numeric_tests
 * (per MEMORY.md note: ctest is not configured in tests/CMakeLists.txt). */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) {
        test_mute_messages();
        done = 1;
    }
}

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

/* ------------------------------------------------------------------ */
/* 1. Wallis cubic, all-real, k=1                                      */
/* ------------------------------------------------------------------ */

static void test_wallis_machine(void) {
    /* x^3 - 2 x - 5 has one real root at ~2.0945514815...; complex pair
     * follows. Canonical k=1 is the real root. */
    check_true("Abs[N[Root[Function[#1^3 - 2 #1 - 5], 1]] - 2.0945514815423265] < 1.*^-12");
}

static void test_wallis_30_digits(void) {
    check_true("Abs[N[Root[Function[#1^3 - 2 #1 - 5], 1], 30] - "
               "2.094551481542326591482386540580] < 1.*^-25");
}

/* ------------------------------------------------------------------ */
/* 2. Cubic with one real + complex conjugate pair                     */
/* ------------------------------------------------------------------ */

static void test_one_real_first(void) {
    /* x^3 + x + 1 has one real root (~-0.6823) and a conjugate pair.
     * Canonical k=1 selects the real root. */
    check_true("Abs[N[Root[Function[#1^3 + #1 + 1], 1], 30] - "
               "(-0.6823278038280193273694837397107)] < 1.*^-25");
}

static void test_conjugate_pair_negative_im_first(void) {
    /* k=2 must be the negative-Im member of the conjugate pair.
     * Verify by polynomial residual instead of comparing against a
     * hand-typed reference (which is fragile beyond a dozen digits). */
    check_true("With[{r = N[Root[Function[#1^3 + #1 + 1], 2], 25]}, "
               "Im[r] < 0 && Abs[r^3 + r + 1] < 1.*^-20]");
}

static void test_conjugate_pair_positive_im_second(void) {
    check_true("With[{r = N[Root[Function[#1^3 + #1 + 1], 3], 25]}, "
               "Im[r] > 0 && Abs[r^3 + r + 1] < 1.*^-20]");
}

/* Roots whose real parts are EQUAL (all on the imaginary axis, or a
 * conjugate pair): the QR estimates' real parts are roundoff of either sign,
 * and the sign differs between LAPACK builds, so the order must not depend
 * on it.  49 + 18 x^2 + x^4 (all roots on the imaginary axis) orders by |Im|,
 * then the negative member first: -1.83 I, 1.83 I, -3.83 I, 3.83 I. */
static void test_equal_real_parts_order(void) {
    const char* want[] = { "-1.828", "1.828", "-3.828", "3.828" };
    for (int k = 1; k <= 4; k++) {
        char buf[256];
        snprintf(buf, sizeof buf,
            "With[{r = N[Root[Function[49 + 18 #1^2 + #1^4], %d]]},"
            " Abs[Re[r]] < 1.*^-10 && Abs[Im[r] - (%s)] < 1.*^-3]", k, want[k - 1]);
        check_true(buf);
    }
    /* x^4 + 1: Re -1/Sqrt[2] pair first, more-negative Im first in each pair */
    check_true("Chop[N[Root[Function[#1^4 + 1], #] & /@ Range[4]] - "
               "{-1 - I, -1 + I, 1 - I, 1 + I}/Sqrt[2]] == {0, 0, 0, 0}");
}

/* ------------------------------------------------------------------ */
/* 3. All-real quartic: x^4 - 10 x^2 + 1, roots = +/- (sqrt(3) +/- sqrt(2)) */
/* ------------------------------------------------------------------ */

static void test_quartic_ascending(void) {
    /* Roots = +/- (sqrt(3) +/- sqrt(2)) in ascending order:
     *   k=1: -(sqrt(3)+sqrt(2)) ~ -3.146264369941972
     *   k=2: -(sqrt(3)-sqrt(2)) ~ -0.317837245195782
     *   k=3:  (sqrt(3)-sqrt(2)) ~  0.317837245195782
     *   k=4:  (sqrt(3)+sqrt(2)) ~  3.146264369941972 */
    check_true("Abs[N[Root[Function[#1^4 - 10 #1^2 + 1], 1], 25] - "
               "(-3.146264369941972342329135)] < 1.*^-20");
    check_true("Abs[N[Root[Function[#1^4 - 10 #1^2 + 1], 2], 25] - "
               "(-0.317837245195782244725758)] < 1.*^-20");
    check_true("Abs[N[Root[Function[#1^4 - 10 #1^2 + 1], 3], 25] - "
               "0.317837245195782244725758] < 1.*^-20");
    check_true("Abs[N[Root[Function[#1^4 - 10 #1^2 + 1], 4], 25] - "
               "3.146264369941972342329135] < 1.*^-20");
}

/* ------------------------------------------------------------------ */
/* 4. High-precision stress: Bring quintic at 200 digits               */
/* ------------------------------------------------------------------ */

static void test_bring_quintic_200_digits(void) {
    /* x^5 - x - 1 (Bring quintic), real root.  Instead of comparing
     * against a hand-typed reference value, we verify by polynomial
     * residual at 200-digit precision.  Tolerance: 10^-180 (well
     * within the 200-digit working precision, leaving room for the
     * post-Newton round-down step). */
    check_true("With[{r = N[Root[Function[#1^5 - #1 - 1], 1], 200]}, "
               "Abs[r^5 - r - 1] < 1.*^-180]");
}

/* ------------------------------------------------------------------ */
/* 5. Squarefree pre-pass                                              */
/* ------------------------------------------------------------------ */

static void test_squarefree_multiplicity(void) {
    /* (x^2 - 2)^3 squarefree radical is x^2 - 2; k=1 -> -sqrt(2). */
    check_true("Abs[N[Root[Function[(#1^2 - 2)^3], 1], 20] - "
               "(-1.4142135623730950488)] < 1.*^-15");
    check_true("Abs[N[Root[Function[(#1^2 - 2)^3], 2], 20] - "
               "1.4142135623730950488] < 1.*^-15");
}

/* ------------------------------------------------------------------ */
/* 6. Solve round-trip                                                 */
/* ------------------------------------------------------------------ */

static void test_solve_roundtrip(void) {
    /* Solve emits Root[g, 1..3] for x^3 + x + 1.  N[Root[g, 1]] is the
     * canonical-1 (real) root, agreeing with the standalone Root form. */
    check_true(
        "With[{r = N[x /. First[Solve[x^3 + x + 1 == 0, x]], 25]}, "
        "Abs[Im[r]] < 1.*^-20 && Abs[r^3 + r + 1] < 1.*^-20]");
}

/* ------------------------------------------------------------------ */
/* 7. Out-of-range k                                                   */
/* ------------------------------------------------------------------ */

static void test_index_out_of_range(void) {
    /* k=7 for a degree-3 polynomial should pass through unchanged. */
    mute_stderr_once();
    check_eq("N[Root[Function[#1^3 + #1 + 1], 7], 20]",
             "Root[Function[Plus[Power[Slot[1], 3], Slot[1], 1]], 7]");
}

/* ------------------------------------------------------------------ */
/* Non-real coefficients (v0.6512).  The reference numbers the roots of  */
/* such a polynomial by real part; among equal real parts the non-real   */
/* roots come first (by imaginary part), then the real one.  Values and  */
/* order checked against the reference.                                  */

#define CLOSE(a, b, tol) \
    "Abs[(" a ") - (" b ")] < " tol

static void test_complex_coefficients_values(void) {
    check_true(CLOSE("N[Root[#^2 + I &, 1]]", "-0.7071067811865475 + 0.7071067811865476 I", "10^-14"));
    check_true(CLOSE("N[Root[#^2 + I &, 2]]", "0.7071067811865475 - 0.7071067811865476 I", "10^-14"));
    check_true(CLOSE("N[Root[#^2 - Sqrt[2] I &, 1]]", "-0.8408964152537146 - 0.8408964152537144 I", "10^-14"));
    /* 30 digits */
    check_true(CLOSE("N[Root[#^5 + I #^2 + 1 &, 3], 30]",
        "-0.196885333304452602162313488934 + 1.117463894146812028378483226312 I", "10^-28"));
    check_true("Precision[N[Root[#^5 + I #^2 + 1 &, 3], 30]] >= 29");
    /* an irreducible Gaussian cubic, 20 digits */
    check_true(CLOSE("N[Root[3 - 2 I + (2 + 2 I) #1 + (-2 - 2 I) #1^2 + #1^3 &, 3], 20]",
        "1.43526071038211799785 + 2.45426466946671164605 I", "10^-18"));
}

static void test_complex_coefficients_order(void) {
    check_true(CLOSE("N[Root[#^3 - I # + 1 &, 1]]", "-1.0046988220954716 - 0.3453793009246787 I", "10^-13"));
    check_true(CLOSE("N[Root[#^3 - I # + 1 &, 2]]", "0.22454268478871658 - 0.6903145345590079 I", "10^-13"));
    check_true(CLOSE("N[Root[#^3 - I # + 1 &, 3]]", "0.780156137306755 + 1.0356938354836867 I", "10^-13"));
    check_true(CLOSE("N[Root[#^4 + (1 + I) #^2 + 2 &, 2]]", "-0.572613580854668 - 0.8077464061219409 I", "10^-13"));
    /* real part first; a real root after the non-real ones of equal real part */
    check_true("Chop[N[Table[Root[(# - I) (# + 2) (# - 3 I) &, k], {k, 3}]] - {-2, I, 3 I}] == {0, 0, 0}");
    check_true("Chop[N[Table[Root[(# - 5) (# - I) (# + 1 - I) &, k], {k, 3}]] - {-1 + I, I, 5}] == {0, 0, 0}");
    check_true("Chop[N[Table[Root[(# - 1 - I) (# - 1 + 2 I) (# - 1) &, k], {k, 3}]] - {1 - 2 I, 1 + I, 1}] == {0, 0, 0}");
    check_true("Chop[N[Table[Root[(# - I) (# + I) (# - 2 I) (# - 1/2) &, k], {k, 4}]] - {-I, I, 2 I, 1/2}] == {0, 0, 0, 0}");
    /* multiplicity counts: a double root occupies two indices */
    check_true("Chop[N[Table[Root[(# - I)^2 (# + 3) &, k], {k, 3}]] - {-3, I, I}] == {0, 0, 0}");
}

static void test_complex_coefficients_noise_and_edges(void) {
    /* solver noise is chopped; a genuine small part is kept */
    check_eq("N[Root[(# - I) (# + 2) (# - 3 I) &, 2]]", "Complex[0.0, 1.0]");
    check_true(CLOSE("Re[N[Root[(# - (10^-10 + I)) (# - 2) &, 1]]]", "10^-10", "10^-24"));
    check_true(CLOSE("Im[N[Root[(# - (1/3 + 10^-12 I)) (# + 1) &, 2]]]", "10^-12", "10^-26"));
    /* symbolic coefficient or out-of-range index: stays a Root */
    check_true("Head[N[Root[#^3 + a I &, 1]]] === Root");
    check_true("Head[Quiet[N[Root[#^2 + I &, 3]]]] === Root");
    /* Root evaluates its index (NHoldAll, as in the reference) */
    check_true("Block[{k = 2}, Root[#^3 - 2 &, k]] === Root[#^3 - 2 &, 2]");
    check_true("Length[Union[N[Table[Root[#^3 - I &, k], {k, 3}]]]] == 3");
    /* a RootSum antiderivative over a Gaussian cubic can now be checked */
    check_true("With[{F = Integrate[(-3 + I - x)/(3 - 2 I + (2 + 2 I) x + (-2 - 2 I) x^2 + x^3), x]}, "
               "Abs[N[(D[F, x] - (-3 + I - x)/(3 - 2 I + (2 + 2 I) x + (-2 - 2 I) x^2 + x^3)) /. x -> 1/3, 20]] < 10^-15]");
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_wallis_machine);
    TEST(test_wallis_30_digits);

    TEST(test_one_real_first);
    TEST(test_conjugate_pair_negative_im_first);
    TEST(test_conjugate_pair_positive_im_second);
    TEST(test_equal_real_parts_order);

    TEST(test_quartic_ascending);

    TEST(test_bring_quintic_200_digits);

    TEST(test_squarefree_multiplicity);

    TEST(test_solve_roundtrip);

    TEST(test_index_out_of_range);

    TEST(test_complex_coefficients_values);
    TEST(test_complex_coefficients_order);
    TEST(test_complex_coefficients_noise_and_edges);

    printf("All N[Root[..]] tests passed!\n");
    return 0;
}
