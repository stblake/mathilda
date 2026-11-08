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

/* An indexed variable x[i] (or an inert user head such as Pair[p, q]) is an
 * ordinary polynomial variable on the FLINT path, so the answer is the one the
 * same fraction gets over symbols, and a 7-variable case is fast (it fell to
 * the classical multivariate Euclid: ~29 s). */
static void test_indexed_variables(void) {
    run_true("(Together[1/(1 - x[1]) - 1/(1 + x[1])] /. x[1] -> x1) === "
             "Together[1/(1 - x1) - 1/(1 + x1)]");
    run_true("Together[a/x[1] + b/x[2]] === (b x[1] + a x[2])/(x[1] x[2])");
    run_true("Cancel[(x[1]^2 - x[2]^2)/(x[1] - x[2])] === x[1] + x[2]");
    run_true("Together[1/(p[a, 1] - q) + 1/(q - p[a, 1])] === 0");
    run_true("With[{s = Sum[x[i], {i, 7}]}, "
             "First[AbsoluteTiming[Together[s^3/(s^2 - x[1]^2) - x[2]/(s + x[1])]]] < 5]");
}

/* Together combines only through the field operations: a function argument
 * and a symbolic exponent are opaque, a power with a rational exponent has its
 * base combined, and lists, equations, inequalities, logical connectives and
 * rules are threaded over (values from the reference). */
static void test_kernel_arguments_opaque(void) {
    run_true("Together[Gamma[2 - D/2]] === Gamma[2 - D/2]");
    run_true("Together[x2 FVD[p, mu] Gamma[2 - D/2]] === x2 FVD[p, mu] Gamma[2 - D/2]");
    run_true("Together[f[1/a + 1/b]] === f[1/a + 1/b]");
    run_true("Together[Abs[1/a + 1/b]] === Abs[1/a + 1/b]");
    run_true("Together[x^(1/a + 1/b)] === x^(1/a + 1/b)");
    run_true("Together[(1/a + 1/b)^n] === (1/a + 1/b)^n");
    run_true("Together[Sqrt[1/a + 1/b]] === Sqrt[(a + b)/(a b)]");
    run_true("Together[1/a + 1/b + Sin[1/a + 1/b]] === (a + b + a b Sin[1/a + 1/b])/(a b)");
    run_true("Together[{1/a + 1/b, g[1/a + 1/b]}] === {(a + b)/(a b), g[1/a + 1/b]}");
    run_true("Together[1/a + 1/b != 1/c] === ((a + b)/(a b) != 1/c)");
    run_true("Together[1/a + 1/b > 0 && x] === ((a + b)/(a b) > 0 && x)");
    run_true("Together[1/a + 1/b -> 1/c + 1/d] === ((a + b)/(a b) -> (c + d)/(c d))");
}

/* A rational exponent with a large denominator q: the power's generator
 * relation g^q - x is raised by repeated squaring, not q multiplications
 * (Together[x^(1 + 2^-40)] and an inexact exponent at WorkingPrecision,
 * p/2^52, hung -- NIntegrate of x^a e^-x at WorkingPrecision with it). */
static void test_large_denominator_exponent(void) {
    run_true("TimeConstrained[Together[x^(1 + 2^-40)], 10] === x^(1 + 2^-40)");
    run_true("TimeConstrained[Head[Together[x^SetPrecision[1.515, 60] y]], 10] === Times");
}

int main(void) {
    setbuf(stdout, NULL);
    printf("Starting field_together_cancel_tests\n");
    symtab_init();
    core_init();

    TEST(test_field_together_value);
    TEST(test_field_cancel);
    TEST(test_field_reduced);
    TEST(test_indexed_variables);
    TEST(test_kernel_arguments_opaque);
    TEST(test_large_denominator_exponent);

    printf("All field Together / Cancel tests passed!\n");
    return 0;
}
