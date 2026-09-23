/* test_field_together_cancel.c — native field-coefficient Together / Cancel.
 *
 * Together/Cancel on AlgebraicNumber[theta,{..}]-coefficient rational functions
 * (the ParallelMixedTower assembly's representation) reduce natively over K(x)
 * via gr_poly over the antic number-field ring (flint_field_together/_cancel).
 *
 * The checks are VALUE-based, not printed-form: Together's canonical form need
 * not match the generic path byte-for-byte, so each test asserts that the reduced
 * fraction equals its input at a generic rational point (to 30 digits) and that a
 * genuine common factor is actually cancelled.  Fields exercised: Q(sqrt2)
 * (degree 2), Q(i), and a degree-4 Root field (the split-specials case that
 * dominates the Charlwood suite's field Together calls).  Uses ASSERT_MSG (never
 * libc assert: Release defines NDEBUG).
 */
#include "poly.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Evaluate `input`; assert its FullForm equals `expected`. */
static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Evaluate `input`; assert it evaluates to True (a numeric/symbolic predicate). */
static void run_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, "True") != 0)
        printf("FAIL: %s\n  expected True, got %s\n", input, s);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s: expected True, got %s", input, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Together preserves value over a degree-4 Root field, Q(i), and Q(sqrt2). */
static void test_field_together_value(void) {
    /* degree-4 Root field (the split-specials case): multi-term sum of fractions. */
    run_true("With[{a = AlgebraicNumber[Root[#^4 - #^2 - 1 &, 1], {0, 1, 0, 0}], "
             "b = AlgebraicNumber[Root[#^4 - #^2 - 1 &, 1], {0, 0, 1, 0}]}, "
             "Abs[N[(Together[a/(x - a) + b/(x^2 - a)] - (a/(x - a) + b/(x^2 - a))) "
             "/. x -> 7/13, 30]] < 10^-25]");
    run_true("With[{a = AlgebraicNumber[Root[#^4 - #^2 - 1 &, 1], {0, 1, 0, 0}], "
             "b = AlgebraicNumber[Root[#^4 - #^2 - 1 &, 1], {0, 0, 1, 0}]}, "
             "Abs[N[(Together[(a x^2 + b)/(x - a) + 1/(x + a) - a/x] "
             "- ((a x^2 + b)/(x - a) + 1/(x + a) - a/x)) /. x -> 5/17, 30]] < 10^-25]");

    /* Gaussian field Q(i). */
    run_true("With[{i = AlgebraicNumber[I, {0, 1}]}, "
             "Abs[N[(Together[i/(x - i) + 1/(x + i) + 3/x] - (i/(x - i) + 1/(x + i) + 3/x)) "
             "/. x -> 9/11, 30]] < 10^-25]");

    /* Q(sqrt2). */
    run_true("With[{s = AlgebraicNumber[Sqrt[2], {0, 1}]}, "
             "Abs[N[(Together[s/(x - s) + 1/(x + s)] - (s/(x - s) + 1/(x + s))) "
             "/. x -> 4/9, 30]] < 10^-25]");
}

/* Cancel removes a genuine common factor over the field. */
static void test_field_cancel(void) {
    /* (x - a)(x - b) / (x - a) == x - b over the degree-4 Root field. */
    run_true("With[{a = AlgebraicNumber[Root[#^4 - #^2 - 1 &, 1], {0, 1, 0, 0}], "
             "b = AlgebraicNumber[Root[#^4 - #^2 - 1 &, 1], {0, 0, 1, 0}]}, "
             "Simplify[Cancel[Expand[(x - a)(x - b)]/(x - a)] - (x - b)] === 0]");
    /* Q(i): (x^2 + 1)/(x - i) == x + i. */
    run_true("With[{i = AlgebraicNumber[I, {0, 1}]}, "
             "Simplify[Cancel[(x^2 + 1)/(x - i)] - (x + i)] === 0]");
}

/* A single-fraction reduction is exact; the result is fully reduced. */
static void test_field_reduced(void) {
    /* Together of a proper single fraction whose num/den share a factor reduces it. */
    run_true("With[{s = AlgebraicNumber[Sqrt[2], {0, 1}]}, "
             "PolynomialGCD[Numerator[#], Denominator[#]] === 1 &[ "
             "Together[s (x - s)/((x - s)(x + 1))] ]] ");
}

int main(void) {
    setbuf(stdout, NULL);
    printf("Starting field_together_cancel_tests\n");
    symtab_init();
    core_init();

    TEST(test_field_together_value);
    TEST(test_field_cancel);
    TEST(test_field_reduced);

    printf("All field Together / Cancel tests passed!\n");
    return 0;
}
