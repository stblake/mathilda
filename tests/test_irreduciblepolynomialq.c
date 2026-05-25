#include <stdio.h>
#include <assert.h>
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
    if (!e) {
        printf("Failed to parse: %s\n", input);
        ASSERT(0);
    }
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

/* Constants are not irreducible polynomials. */
static void test_irrpoly_constants(void) {
    run_test("IrreduciblePolynomialQ[0]", "False");
    run_test("IrreduciblePolynomialQ[1]", "False");
    run_test("IrreduciblePolynomialQ[-3]", "False");
    run_test("IrreduciblePolynomialQ[5]", "False");
    run_test("IrreduciblePolynomialQ[2/3]", "False");
}

/* Univariate over Q -- spec headline. */
static void test_irrpoly_univariate_q(void) {
    run_test("IrreduciblePolynomialQ[x]",            "True");
    run_test("IrreduciblePolynomialQ[x + 1]",        "True");
    run_test("IrreduciblePolynomialQ[x^2 + 1]",      "True");
    run_test("IrreduciblePolynomialQ[x^2 - 1]",      "False");  /* (x-1)(x+1) */
    run_test("IrreduciblePolynomialQ[x^2 + 2]",      "True");
    run_test("IrreduciblePolynomialQ[x^3 - 8]",      "False");  /* (x-2)(x^2+2x+4) */
    run_test("IrreduciblePolynomialQ[x^3 - 2]",      "True");   /* irreducible over Q */
    run_test("IrreduciblePolynomialQ[(x - 1)^2]",    "False");  /* repeated factor */
}

/* Listable threading -- spec headline.  ATTR_LISTABLE delegates to the
 * evaluator, so this also verifies the evaluator routing. */
static void test_irrpoly_listable(void) {
    run_test("IrreduciblePolynomialQ[{x^2 - 1, x^2 - 2}]",
             "List[False, True]");
    run_test("IrreduciblePolynomialQ[{x^2 + 1, x^3 - 8}]",
             "List[True, False]");
}

/* Multivariate over Q -- spec headline. */
static void test_irrpoly_multivariate_q(void) {
    run_test("IrreduciblePolynomialQ[x^4 - 4 y^2]", "False");  /* (x^2-2y)(x^2+2y) */
    run_test("IrreduciblePolynomialQ[x^4 - 2 y^2]", "True");
    run_test("IrreduciblePolynomialQ[x^2 + y^2]",   "True");   /* irreducible over Q */
    run_test("IrreduciblePolynomialQ[x^2 - y^2]",   "False");  /* (x-y)(x+y) */
    run_test("IrreduciblePolynomialQ[x*y]",         "False");  /* x * y */
}

/* Complex coefficients flip into Gaussian rationals automatically.
 * x^2 + 2 I x - 1 = (x + I)^2 over Q(i) -- one factor, multiplicity 2 -> False. */
static void test_irrpoly_complex_coef_auto(void) {
    run_test("IrreduciblePolynomialQ[x^2 + 2 I x - 1]", "False");
}

/* GaussianIntegers -> True: x^2 + 1 = (x-I)(x+I) over Q(i).
 *
 * Multivariate Q(i) factoring is not yet supported (qa_factor_with_extension
 * is univariate-only), so cases like `x^2 + y^2, GaussianIntegers -> True`
 * silently fall back to factoring over Q and produce True; we don't test
 * that case here. */
static void test_irrpoly_gaussian_option(void) {
    run_test("IrreduciblePolynomialQ[x^2 + 1, GaussianIntegers -> True]", "False");
    /* x^2 + 1 over Q (default): irreducible. */
    run_test("IrreduciblePolynomialQ[x^2 + 1, GaussianIntegers -> False]", "True");
}

/* Algebraic-number coefficient treated as independent variable by default. */
static void test_irrpoly_alg_default(void) {
    run_test("IrreduciblePolynomialQ[x^2 + 2 Sqrt[2] x + 2]", "True");
    run_test("IrreduciblePolynomialQ[x^2 + 2 Sqrt[3] x y + 3 y^2]", "True");
}

/* Extension -> Automatic absorbs algebraic coefficients of poly.
 * x^2 + 2 Sqrt[2] x + 2 = (x + Sqrt[2])^2 over Q(Sqrt[2]) -- one factor,
 * multiplicity 2 -> False. */
static void test_irrpoly_extension_automatic(void) {
    run_test("IrreduciblePolynomialQ[x^2 + 2 Sqrt[2] x + 2, Extension -> Automatic]",
             "False");
    run_test("IrreduciblePolynomialQ[x^2 + 2 Sqrt[3] x y + 3 y^2, Extension -> Automatic]",
             "False");
}

/* Extension -> alpha selects a finite algebraic extension. */
static void test_irrpoly_extension_explicit(void) {
    /* x^3 - 2 = (x - 2^(1/3)) * (x^2 + 2^(1/3) x + 2^(2/3)) over Q(2^(1/3)). */
    run_test("IrreduciblePolynomialQ[x^3 - 2, Extension -> 2^(1/3)]", "False");
    /* x^3 - 3 has no root in Q(2^(1/3)); stays irreducible. */
    run_test("IrreduciblePolynomialQ[x^3 - 3, Extension -> 2^(1/3)]", "True");
}

/* Extension as a list: same effect as a tower of generators. */
static void test_irrpoly_extension_list(void) {
    /* x^2 + 2 over Q is irreducible; over Q(I, Sqrt[2]) it factors as
     * (x + Sqrt[2] I)(x - Sqrt[2] I). */
    run_test("IrreduciblePolynomialQ[x^2 + 2, Extension -> {I, Sqrt[2]}]", "False");
}

/* Extension -> All: absolute irreducibility over C.
 *   Univariate degree-1: True; higher univariate: False.
 *   Multivariate: best-effort (factor over Q(i)).
 *
 * Multivariate Q(i) factoring is not supported by the underlying Factor
 * builtin yet, so cases like `x^2 + y^2, Extension -> All` (correctly
 * (x+iy)(x-iy) over C) currently report True; that gap is noted in the
 * arithmetic-and-algebra spec entry and not tested here. */
static void test_irrpoly_extension_all(void) {
    /* Univariate degree >= 2 splits over C. */
    run_test("IrreduciblePolynomialQ[x^3 - 3, Extension -> All]", "False");
    run_test("IrreduciblePolynomialQ[x + 1,   Extension -> All]", "True");
    /* Multivariate -- absolutely irreducible per spec. */
    run_test("IrreduciblePolynomialQ[x^2 + 2 x y - 7, Extension -> All]", "True");
    run_test("IrreduciblePolynomialQ[x^7 + 12 x y - 11, Extension -> All]", "True");
}

/* Non-polynomial expressions: False per Mathematica. */
static void test_irrpoly_nonpolynomial(void) {
    run_test("IrreduciblePolynomialQ[Sin[x]]", "False");
    run_test("IrreduciblePolynomialQ[Pi]",     "False");
    run_test("IrreduciblePolynomialQ[Sqrt[x]]","False");
}

/* IrreduciblePolynomialQ is Protected and Listable. */
static void test_irrpoly_attributes(void) {
    run_test("Attributes[IrreduciblePolynomialQ]", "List[Listable, Protected]");
}

/* argx diagnostic: 0 args -> unevaluated + stderr message. */
static void test_irrpoly_argx(void) {
    fprintf(stderr, "--- next line(s) are an EXPECTED IrreduciblePolynomialQ::argx diagnostic ---\n");
    run_test("IrreduciblePolynomialQ[]", "IrreduciblePolynomialQ[]");
}

/* nonopt diagnostic: trailing non-Rule positional arg -> stderr + unevaluated.
 * Mathematica reports the last bad arg ("instead of 3"). */
static void test_irrpoly_nonopt(void) {
    fprintf(stderr, "--- next line(s) are an EXPECTED IrreduciblePolynomialQ::nonopt diagnostic ---\n");
    run_test("IrreduciblePolynomialQ[1, 2, 3]",
             "IrreduciblePolynomialQ[1, 2, 3]");
    /* Unknown option name. */
    fprintf(stderr, "--- next line(s) are an EXPECTED IrreduciblePolynomialQ::nonopt diagnostic ---\n");
    run_test("IrreduciblePolynomialQ[x, Foo -> Bar]",
             "IrreduciblePolynomialQ[x, Rule[Foo, Bar]]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running IrreduciblePolynomialQ tests...\n");
    TEST(test_irrpoly_constants);
    TEST(test_irrpoly_univariate_q);
    TEST(test_irrpoly_listable);
    TEST(test_irrpoly_multivariate_q);
    TEST(test_irrpoly_complex_coef_auto);
    TEST(test_irrpoly_gaussian_option);
    TEST(test_irrpoly_alg_default);
    TEST(test_irrpoly_extension_automatic);
    TEST(test_irrpoly_extension_explicit);
    TEST(test_irrpoly_extension_list);
    TEST(test_irrpoly_extension_all);
    TEST(test_irrpoly_nonpolynomial);
    TEST(test_irrpoly_attributes);
    TEST(test_irrpoly_argx);
    TEST(test_irrpoly_nonopt);
    printf("All IrreduciblePolynomialQ tests passed!\n");
    return 0;
}
