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
        freopen("/dev/null", "w", stderr);
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

    TEST(test_quartic_ascending);

    TEST(test_bring_quintic_200_digits);

    TEST(test_squarefree_multiplicity);

    TEST(test_solve_roundtrip);

    TEST(test_index_out_of_range);

    printf("All N[Root[..]] tests passed!\n");
    return 0;
}
