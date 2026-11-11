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

/* ---- nested constant radicals among the coefficients (src/poly/algfield.c) ----
 *
 * D2 of tasks/polygcd_together_deficiencies: a 4639-leaf rational function of x
 * with coefficients in Q(I, Sqrt[2], q), q = Sqrt[3/4 - 1/Sqrt[2]] =
 * 1/Sqrt[2] - 1/2, so the field is Q(I, Sqrt[2]) and the value is
 * 1 + Sqrt[2] x.  Together, Cancel and Simplify ran past 20 s.  Every call is
 * under TimeConstrained so a regression fails instead of hanging. */
static const char* d2_path = "../together_nested_radical_d2.m";

#define D2_ANSWER "Plus[1, Times[Power[2, Rational[1, 2]], x]]"

static void test_nested_radical_d2(void) {
    char buf[1200];
    snprintf(buf, sizeof buf, "d2e = Get[\"%s\"]; LeafCount[d2e] == 4639", d2_path);
    run_true(buf);
    /* the input really is 1 + Sqrt[2] x */
    run_true("TimeConstrained[PossibleZeroQ[d2e - (1 + Sqrt[2] x)], 20, False]");

    /* Extension -> Automatic: over the number field, q denested */
    run_test("TimeConstrained[Together[d2e, Extension -> Automatic], 10, $Failed]", D2_ANSWER);
    run_test("TimeConstrained[Cancel[d2e, Extension -> Automatic], 10, $Failed]", D2_ANSWER);
    /* Simplify takes the field normal form */
    run_test("TimeConstrained[Simplify[d2e], 20, $Failed]", D2_ANSWER);

    /* default: the nested radical is an opaque generator (no denesting, as in
     * the reference); finishes, compact, and value-equal to the input */
    run_true("d2t = TimeConstrained[Together[d2e], 10, $Failed]; "
             "d2t =!= $Failed && LeafCount[d2t] <= 200");
    run_true("And @@ (Abs[N[(d2t - (1 + Sqrt[2] x)) /. x -> #, 30]] < 10^-20 & /@ "
             "{1/3, 7/5 + I, -2, 11/13})");
    run_true("TimeConstrained[PossibleZeroQ[d2t - d2e], 20, False]");
    run_true("d2c = TimeConstrained[Cancel[d2e], 10, $Failed]; "
             "d2c =!= $Failed && LeafCount[d2c] <= 200 && "
             "TimeConstrained[PossibleZeroQ[d2c - (1 + Sqrt[2] x)], 20, False]");
}

/* Smaller variants: a denestable nested radical, several powers of it, a
 * non-denestable one (Sqrt[1 + Sqrt[2]]), mixed with I, two variables. */
static void test_nested_radical_variants(void) {
    run_true("q = Sqrt[3/4 - 1/Sqrt[2]]; r1 = Sqrt[1 + Sqrt[2]]; True");
    /* Extension -> Automatic: the field normal form */
    run_test("Together[1/(x - q) + 1/(x + q), Extension -> Automatic]",
             "Times[8, x, Power[Plus[-3, Times[2, Power[2, Rational[1, 2]]], "
             "Times[4, Power[x, 2]]], -1]]");
    run_test("Together[(x^2 - q^2)/(x - q), Extension -> Automatic]",
             "Times[Rational[1, 2], Plus[-1, Power[2, Rational[1, 2]], Times[2, x]]]");
    run_test("Cancel[(x^2 - q^2)/(x - q), Extension -> Automatic]",
             "Times[Rational[1, 2], Plus[-1, Power[2, Rational[1, 2]], Times[2, x]]]");
    /* q - (1/Sqrt[2] - 1/2) = 0 is visible only over the field */
    run_test("Together[1/(x - q) - 1/(x - 1/Sqrt[2] + 1/2), Extension -> Automatic]", "0");
    /* q, q^3 = (3/4 - 1/Sqrt[2])^(3/2) and 1/(1 + q) together */
    run_test("Together[q^3 x + q x^2 + 1/(1 + q), Extension -> Automatic]",
             "Times[Rational[1, 8], Plus[-16, Times[16, Power[2, Rational[1, 2]]], "
             "Times[-7, x], Times[5, Power[2, Rational[1, 2]], x], Times[-4, Power[x, 2]], "
             "Times[4, Power[2, Rational[1, 2]], Power[x, 2]]]]");
    /* non-denestable: Sqrt[1 + Sqrt[2]] stays, the factor cancels */
    run_test("Together[(x^2 - 1 - Sqrt[2])/(x - r1), Extension -> Automatic]",
             "Plus[Power[Plus[1, Power[2, Rational[1, 2]]], Rational[1, 2]], x]");
    run_test("Together[(x^2 + 1 + Sqrt[2])/(x - I r1), Extension -> Automatic]",
             "Plus[Times[Complex[0, 1], Power[Plus[1, Power[2, Rational[1, 2]]], "
             "Rational[1, 2]]], x]");
    run_test("Simplify[(x^2 - 1 - Sqrt[2])/(x - r1)]",
             "Plus[Power[Plus[1, Power[2, Rational[1, 2]]], Rational[1, 2]], x]");
    /* two variables */
    run_true("With[{u = (x - y q^2)/(x - q) + y q}, "
             "With[{t = Together[u, Extension -> Automatic]}, FreeQ[t, q] && "
             "Abs[N[(t - u) /. {x -> 3/7, y -> -5/11}, 30]] < 10^-20]]");

    /* default: reference-consistent, the nested radical is kept */
    run_true("!FreeQ[Together[(x^2 - q^2)/(x - q)], Sqrt[3/4 - 1/Sqrt[2]]]");
    run_test("Together[1/(x - q) + 1/(x + q)]",
             "Times[8, Power[2, Rational[1, 2]], x, Power[Plus[4, Times[-3, "
             "Power[2, Rational[1, 2]]], Times[4, Power[2, Rational[1, 2]], Power[x, 2]]], -1]]");
    run_true("With[{u = q^3 x + q x^2 + 1/(1 + q)}, Abs[N[(Together[u] - u) /. x -> 3/7, 30]] < 10^-20]");
}

/* Inputs without a nested constant radical take the paths they always took. */
static void test_nested_radical_guards(void) {
    run_test("Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2])]",
             "Times[2, x, Power[Plus[-2, Power[x, 2]], -1]]");
    run_test("Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2]), Extension -> Automatic]",
             "Times[2, x, Power[Plus[-2, Power[x, 2]], -1]]");
    run_test("Together[(x^2 - Sqrt[2])/(x - 2^(1/4))]", "Plus[Power[2, Rational[1, 4]], x]");
    run_test("Together[1/(x - k^(1/3)) + 1/(x + k^(1/3))]",
             "Times[-2, x, Power[Plus[Power[k, Rational[2, 3]], Times[-1, Power[x, 2]]], -1]]");
    run_test("Together[a/b + c/d]",
             "Times[Power[b, -1], Power[d, -1], Plus[Times[b, c], Times[a, d]]]");
    run_test("Simplify[Sqrt[3/4 - 1/Sqrt[2]]]", "Plus[Rational[-1, 2], Power[2, Rational[-1, 2]]]");
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    printf("Starting field_together_cancel_tests\n");
    if (argc > 1) d2_path = argv[1];
    symtab_init();
    core_init();

    TEST(test_field_together_value);
    TEST(test_field_cancel);
    TEST(test_field_reduced);
    TEST(test_indexed_variables);
    TEST(test_kernel_arguments_opaque);
    TEST(test_large_denominator_exponent);
    TEST(test_nested_radical_d2);
    TEST(test_nested_radical_variants);
    TEST(test_nested_radical_guards);

    printf("All field Together / Cancel tests passed!\n");
    return 0;
}
