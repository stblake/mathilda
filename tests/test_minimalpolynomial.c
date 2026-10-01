/* test_minimalpolynomial.c -- unit tests for MinimalPolynomial.
 *
 * Coverage:
 *   - Quadratic / biquadratic radicals (Sqrt, sums of radicals).
 *   - Nested radicals.
 *   - Complex algebraic numbers (Gaussian / roots of unity).
 *   - Root[] objects shifted by a constant.
 *   - Rationals and integers (degree-1 minimal polynomials).
 *   - Cube roots / Surd.
 *   - Automatic threading over lists (Listable).
 *   - Pure-function form MinimalPolynomial[s].
 *   - The Extension -> a characteristic-polynomial form, including the
 *     non-membership case (left unevaluated).
 *   - Non-algebraic input left unevaluated (no crash).
 *
 * Output is compared in FullForm so the structure is unambiguous.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, res_str);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* Radicals of algebraic numbers. */
static void test_radicals(void) {
    run_test("MinimalPolynomial[Sqrt[3], x]",
             "Plus[-3, Power[x, 2]]");
    run_test("MinimalPolynomial[Sqrt[2] + Sqrt[3], x]",
             "Plus[1, Times[-10, Power[x, 2]], Power[x, 4]]");
    /* 2^(1/3): x^3 - 2. */
    run_test("MinimalPolynomial[2^(1/3), x]",
             "Plus[-2, Power[x, 3]]");
    /* 5^(2/3): x^3 - 25. */
    run_test("MinimalPolynomial[5^(2/3), x]",
             "Plus[-25, Power[x, 3]]");
}

/* Nested radicals. */
static void test_nested(void) {
    run_test("MinimalPolynomial[Sqrt[1 + Sqrt[3]], x]",
             "Plus[-2, Times[-2, Power[x, 2]], Power[x, 4]]");
    run_test("MinimalPolynomial[Sqrt[1 + Sqrt[2]], x]",
             "Plus[-1, Times[-2, Power[x, 2]], Power[x, 4]]");
}

/* Complex algebraic numbers. */
static void test_complex(void) {
    run_test("MinimalPolynomial[(1 + I)/Sqrt[2], x]",
             "Plus[1, Power[x, 4]]");
    /* I itself: x^2 + 1. */
    run_test("MinimalPolynomial[I, x]",
             "Plus[1, Power[x, 2]]");
}

/* Root[] objects (defining polynomial recovered exactly). */
static void test_root_objects(void) {
    run_test("MinimalPolynomial[Root[2 #1^3 - 2 #1 + 7 &, 1] + 17, x]",
             "Plus[-9785, Times[1732, x], Times[-102, Power[x, 2]], "
             "Times[2, Power[x, 3]]]");
}

/* Rationals and integers: degree-1 minimal polynomials. */
static void test_rational(void) {
    run_test("MinimalPolynomial[1/2, x]",  "Plus[-1, Times[2, x]]");
    run_test("MinimalPolynomial[5, x]",    "Plus[-5, x]");
    run_test("MinimalPolynomial[-3/4, x]", "Plus[3, Times[4, x]]");
}

/* Automatic threading over lists. */
static void test_listable(void) {
    run_test("MinimalPolynomial[{Sqrt[3], Root[-2 + #1^3 &, 2] + 1}, x]",
             "List[Plus[-3, Power[x, 2]], "
             "Plus[-3, Times[3, x], Times[-3, Power[x, 2]], Power[x, 3]]]");
}

/* Pure-function form. */
static void test_pure_function(void) {
    run_test("MinimalPolynomial[Sqrt[2] + Sqrt[3]]",
             "Function[Plus[1, Times[-10, Power[Slot[1], 2]], Power[Slot[1], 4]]]");
}

/* Extension -> a : characteristic polynomial over Q(a). */
static void test_extension(void) {
    run_test("MinimalPolynomial[Sqrt[2], x, Extension -> E^(I Pi/4)]",
             "Plus[4, Times[-4, Power[x, 2]], Power[x, 4]]");
    /* The characteristic polynomial is a power of the minimal polynomial. */
    run_test("Factor[MinimalPolynomial[Sqrt[2], x, Extension -> E^(I Pi/4)]]",
             "Power[Plus[-2, Power[x, 2]], 2]");
    /* k = 1 reproduces the ordinary minimal polynomial. */
    run_test("MinimalPolynomial[Sqrt[2], x, Extension -> Sqrt[2]]",
             "Plus[-2, Power[x, 2]]");
    /* Non-membership: Sqrt[3] is not in Q(E^(I Pi/4)) -> unevaluated. */
    run_test("MinimalPolynomial[Sqrt[3], x, Extension -> E^(I Pi/4)]",
             "MinimalPolynomial[Power[3, Rational[1, 2]], x, "
             "Rule[Extension, Power[E, Times[Complex[0, Rational[1, 4]], Pi]]]]");
}

/* Non-algebraic / malformed input is left unevaluated (no crash). */
static void test_unevaluated(void) {
    run_test("MinimalPolynomial[Pi, x]",
             "MinimalPolynomial[Pi, x]");
    run_test("MinimalPolynomial[x + Log[2], x]",
             "MinimalPolynomial[Plus[Log[2], x], x]");
}

/* A Root[] object whose NUMERIC refinement does not converge must still get its
 * minimal polynomial.  The factor among the candidates is chosen by an 80-digit
 * numeric test with a 1e-3 ceiling -- fine when there is a choice to make, but it
 * was run even with a SINGLE candidate, where G(s) = 0 holds by construction.
 * `N[this, 80]` emits three Root::conv warnings and returns the Root unevaluated,
 * so the test could not score it, and MinimalPolynomial refused outright on a
 * polynomial that IrreduciblePolynomialQ confirms irreducible.
 *
 * It refused silently, and the caller did not notice: ParallelMixed.m's
 * NontorsionDivisor fed the unevaluated result to CoefficientList, took the
 * one-element list of garbage as a number field's defining polynomial, and the
 * mod-p non-torsion certificate then ran past the whole 120 s budget on it
 * (integration corpus #109). */
static void test_root_nonconvergent_numerics(void) {
    run_test("MinimalPolynomial[Root[2298903058544965460896492004394531250000 "
             "- 250417320057749748229980468750000000 #1 "
             "+ 11606206744909286499023437500000 #1^2 "
             "- 294947862625122070312500000 #1^3 "
             "+ 4494413375854492187500 #1^4 - 42187500000000000 #1^5 "
             "+ 239062500000 #1^6 - 750000 #1^7 + #1^8 &, 3], x]",
             "Plus[2298903058544965460896492004394531250000, "
             "Times[-250417320057749748229980468750000000, x], "
             "Times[11606206744909286499023437500000, Power[x, 2]], "
             "Times[-294947862625122070312500000, Power[x, 3]], "
             "Times[4494413375854492187500, Power[x, 4]], "
             "Times[-42187500000000000, Power[x, 5]], "
             "Times[239062500000, Power[x, 6]], "
             "Times[-750000, Power[x, 7]], Power[x, 8]]");
}

int main(void) {
    symtab_init();
    core_init();

    test_radicals();
    test_nested();
    test_complex();
    test_root_objects();
    test_root_nonconvergent_numerics();
    test_rational();
    test_listable();
    test_pure_function();
    test_extension();
    test_unevaluated();

    printf("\nAll MinimalPolynomial tests passed.\n");
    return 0;
}
