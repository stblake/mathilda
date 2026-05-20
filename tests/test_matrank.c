/*
 * tests/test_matrank.c
 *
 * Unit tests for MatrixRank[m] (and the Method / Tolerance options).
 *
 * Verification strategies:
 *
 *   1. Exact comparison via run_test for cases where the rank is
 *      unambiguous (integer / rational / symbolic matrices).
 *
 *   2. Tolerance-aware numerical tests where the input matrix has
 *      machine-precision floats; we assert the rank value directly
 *      (it is an Integer in every case).
 *
 *   3. Methods (DivisionFreeRowReduction / OneStepRowReduction /
 *      Automatic) must all produce identical rank on the same input.
 *
 * Memory: every parse_expression / evaluate result is freed.  The
 * binary is intended to run under valgrind --leak-check=full with no
 * "definitely lost" bytes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Compare FullForm representation. */
static void run_test(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, res_str);
        free(res_str);
        expr_free(res);
        expr_free(parsed);
        ASSERT(0);
    }
    printf("  PASS: %s -> %s\n", input, res_str);
    free(res_str);
    expr_free(res);
    expr_free(parsed);
}

/* Parse, evaluate, assert the result is the integer `expected`. */
static void assert_rank_eq(const char* src, int64_t expected) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* res = evaluate(parsed);
    if (res->type != EXPR_INTEGER || res->data.integer != expected) {
        char* s = expr_to_string(res);
        fprintf(stderr, "FAIL: %s\n  expected: Integer %lld\n  got:      %s\n",
                src, (long long)expected, s);
        free(s);
        expr_free(res);
        expr_free(parsed);
        ASSERT(0);
    }
    printf("  PASS: %s -> %lld\n", src, (long long)expected);
    expr_free(res);
    expr_free(parsed);
}

/* ===================== Exact-form tests ===================== */

static void test_matrank_3x3_rank2(void) {
    /* Classic rank-2 matrix. */
    run_test("MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]", "2");
}

static void test_matrank_3x3_full(void) {
    run_test("MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]", "3");
}

static void test_matrank_symbolic_full(void) {
    /* Generic symbolic matrix: assumed full rank. */
    run_test("MatrixRank[{{a, b, c}, {d, e, f}, {g, h, i}}]", "3");
}

static void test_matrank_symbolic_2x2_dep(void) {
    /* {{a, b}, {2a, 2b}} -- row 2 = 2 * row 1. */
    run_test("MatrixRank[{{a, b}, {2 a, 2 b}}]", "1");
}

static void test_matrank_symbolic_2x2_indep(void) {
    /* {{a, b}, {c, d}} -- generic, no relation. */
    run_test("MatrixRank[{{a, b}, {c, d}}]", "2");
}

static void test_matrank_3x5_rect(void) {
    /* From the spec: {{0,5,2,4,4},{2,5,0,4,0},{5,1,5,4,5}} -> 3. */
    run_test("MatrixRank[{{0, 5, 2, 4, 4}, {2, 5, 0, 4, 0}, {5, 1, 5, 4, 5}}]",
             "3");
}

static void test_matrank_2x4_rect(void) {
    /* From the spec: {{2,2,3,4},{3,2,1,3}} -> 2. */
    run_test("MatrixRank[{{2, 2, 3, 4}, {3, 2, 1, 3}}]", "2");
}

static void test_matrank_zero_matrix(void) {
    run_test("MatrixRank[{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}]", "0");
    run_test("MatrixRank[{{0, 0}, {0, 0}, {0, 0}}]", "0");
}

static void test_matrank_identity(void) {
    run_test("MatrixRank[IdentityMatrix[3]]", "3");
    run_test("MatrixRank[IdentityMatrix[5]]", "5");
    run_test("MatrixRank[IdentityMatrix[1]]", "1");
}

static void test_matrank_diag(void) {
    /* DiagonalMatrix{1,2,0,4} -> 3 (one zero on diagonal). */
    run_test("MatrixRank[DiagonalMatrix[{1, 2, 0, 4}]]", "3");
}

static void test_matrank_rational_full(void) {
    /* {{1, 1/2}, {1/3, 1/4}} -- det = 1/4 - 1/6 = 1/12 != 0 -> rank 2. */
    run_test("MatrixRank[{{1, 1/2}, {1/3, 1/4}}]", "2");
}

static void test_matrank_rational_dep(void) {
    /* {{1, 1/2, 1/3}, {2, 1, 2/3}} -- row 2 = 2 * row 1 -> rank 1. */
    run_test("MatrixRank[{{1, 1/2, 1/3}, {2, 1, 2/3}}]", "1");
}

static void test_matrank_4x4_rank3(void) {
    /* Standard {1..10..16} block with row sum dependency. */
    run_test("MatrixRank[{{1, 2, 3, 10}, {4, 5, 6, 11}, "
                       "{7, 8, 9, 12}, {13, 14, 15, 16}}]",
             "3");
}

static void test_matrank_bigint(void) {
    /* Two rows of large integers; second is 4 * first -> rank 1. */
    run_test("MatrixRank[{{100000000000000000001, 200000000000000000002}, "
             "{400000000000000000004, 800000000000000000008}}]",
             "1");
}

static void test_matrank_bigint_indep(void) {
    /* Two linearly independent rows of large integers -> rank 2. */
    run_test("MatrixRank[{{100000000000000000001, 200000000000000000002}, "
             "{300000000000000000005, 700000000000000000011}}]",
             "2");
}

static void test_matrank_1x3(void) {
    /* 1x3 non-zero row -> rank 1. */
    run_test("MatrixRank[{{1, 2, 3}}]", "1");
    /* 1x3 zero row -> rank 0. */
    run_test("MatrixRank[{{0, 0, 0}}]", "0");
}

static void test_matrank_3x1(void) {
    /* 3x1 non-zero -> rank 1. */
    run_test("MatrixRank[{{1}, {2}, {3}}]", "1");
    run_test("MatrixRank[{{0}, {0}, {0}}]", "0");
}

/* ===================== Method-option tests ===================== */

static void test_matrank_method_divfree(void) {
    run_test("MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
             "Method -> \"DivisionFreeRowReduction\"]",
             "2");
}

static void test_matrank_method_onestep(void) {
    run_test("MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
             "Method -> \"OneStepRowReduction\"]",
             "2");
}

static void test_matrank_method_cofactor(void) {
    /* Cofactor on a square non-singular matrix. */
    run_test("MatrixRank[{{1, 2}, {3, 4}}, Method -> \"CofactorExpansion\"]",
             "2");
    /* Cofactor on a singular matrix -- falls back to DivFree. */
    run_test("MatrixRank[{{1, 2}, {2, 4}}, Method -> \"CofactorExpansion\"]",
             "1");
}

static void test_matrank_method_automatic_symbol(void) {
    /* Bare Automatic symbol. */
    run_test("MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
             "Method -> Automatic]",
             "2");
}

static void test_matrank_method_automatic_string(void) {
    /* "Automatic" string. */
    run_test("MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
             "Method -> \"Automatic\"]",
             "2");
}

static void test_matrank_method_unknown(void) {
    /* Unknown method name: call left unevaluated. */
    run_test("MatrixRank[{{1, 2}}, Method -> \"NoSuchMethod\"]",
             "MatrixRank[List[List[1, 2]], Rule[Method, \"NoSuchMethod\"]]");
}

/* ===================== Tolerance-option tests ===================== */

static void test_matrank_real_full(void) {
    /* Generic real 3x3 with no near-collinearity -> rank 3. */
    assert_rank_eq("MatrixRank[{{1.25, 3.2, 3.2}, "
                              "{7.9, -1.4, 5.1}, "
                              "{1.1, 2.5, -1.5}}]", 3);
}

static void test_matrank_real_singular(void) {
    /* Row 3 = 2 * Row 1 + Row 2 -- rank 2 in exact arithmetic.  With
     * default machine tolerance on N[m] the third row vanishes. */
    assert_rank_eq("MatrixRank[{{1.5, 4.75, -3.2}, "
                              "{-1.7, 6.7, -9.3}, "
                              "{1.3, 16.2, -15.7}}]", 2);
}

static void test_matrank_complex_dep(void) {
    /* {{1+I, 2, 3-2I}, {0, 4, 5I}, {1+I, 6, 3+3I}} - row 3 = row 1 + row 2,
     * so rank is 2. */
    assert_rank_eq("MatrixRank[{{1. + I, 2, 3 - 2 I}, "
                              "{0, 4, 5 I}, "
                              "{1. + I, 6, 3 + 3 I}}]", 2);
}

/* The spec's tolerance ladder: */
static void test_matrank_tolerance_spec(void) {
    /* Exact: full rank 3. */
    assert_rank_eq("MatrixRank[{{1, 1, 1}, {0, 10^-10, 0}, {0, 0, 10^-20}}]", 3);

    /* Machine-precision default: 10^-20 falls below default tolerance,
     * so the last row is treated as zero. */
    assert_rank_eq("MatrixRank[N[{{1, 1, 1}, "
                                "{0, 10^-10, 0}, "
                                "{0, 0, 10^-20}}]]", 2);

    /* Tolerance -> 0 with machine precision: all entries count. */
    assert_rank_eq("MatrixRank[N[{{1, 1, 1}, "
                                "{0, 10^-10, 0}, "
                                "{0, 0, 10^-20}}], Tolerance -> 0]", 3);

    /* Tolerance > 10^-10: only the first row survives. */
    assert_rank_eq("MatrixRank[N[{{1, 1, 1}, "
                                "{0, 10^-10, 0}, "
                                "{0, 0, 10^-20}}], Tolerance -> 10^-8]", 1);
}

/* Exact matrix with explicit Tolerance -> still works via numerical
 * path (since the integer-valued matrix is trivially convertible). */
static void test_matrank_exact_with_tolerance(void) {
    /* {{1, 1}, {1, 1}} - rank 1 regardless of tolerance up to 1. */
    assert_rank_eq("MatrixRank[{{1, 1}, {1, 1}}, Tolerance -> 0]", 1);
    /* With Tolerance > 1: everything counts as zero -> rank 0. */
    assert_rank_eq("MatrixRank[{{1, 1}, {1, 1}}, Tolerance -> 2]", 0);
}

static void test_matrank_tolerance_rational(void) {
    /* Tolerance accepts rationals. */
    assert_rank_eq("MatrixRank[{{1, 1}, {1, 1}}, Tolerance -> 1/2]", 1);
}

static void test_matrank_tolerance_invalid(void) {
    /* Negative tolerance: call left unevaluated.  Symbolic tolerance
     * (e.g. Tolerance -> x) is also not accepted. */
    run_test("MatrixRank[{{1, 2}}, Tolerance -> -1]",
             "MatrixRank[List[List[1, 2]], Rule[Tolerance, -1]]");
    run_test("MatrixRank[{{1, 2}}, Tolerance -> y]",
             "MatrixRank[List[List[1, 2]], Rule[Tolerance, y]]");
}

static void test_matrank_method_and_tolerance(void) {
    /* Both options combined. */
    assert_rank_eq("MatrixRank[N[{{1, 1, 1}, {0, 10^-10, 0}, {0, 0, 10^-20}}], "
                  "Method -> Automatic, Tolerance -> 10^-8]", 1);
    /* Order independent. */
    assert_rank_eq("MatrixRank[N[{{1, 1, 1}, {0, 10^-10, 0}, {0, 0, 10^-20}}], "
                  "Tolerance -> 10^-8, Method -> Automatic]", 1);
}

/* ===================== Shape / arity errors ===================== */

static void test_matrank_arity(void) {
    run_test("MatrixRank[]", "MatrixRank[]");
}

static void test_matrank_non_matrix(void) {
    run_test("MatrixRank[x]", "MatrixRank[x]");
    run_test("MatrixRank[{}]", "MatrixRank[List[]]");
    run_test("MatrixRank[{1, 2, 3}]", "MatrixRank[List[1, 2, 3]]");
    /* 3-tensor (rank-3) input -- not a matrix. */
    run_test("MatrixRank[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]",
             "MatrixRank[List[List[List[1, 2], List[3, 4]], "
             "List[List[5, 6], List[7, 8]]]]");
}

/* ===================== Stress / docstring ===================== */

/* The same matrix returned by every method should give the same rank. */
static void test_matrank_method_agreement(void) {
    int64_t r_div, r_one, r_cof, r_auto;
    {
        Expr* p = parse_expression("MatrixRank[{{1, 2, 3}, {4, 5, 6}, "
                                  "{7, 8, 9}}, Method -> \"DivisionFreeRowReduction\"]");
        Expr* r = evaluate(p);
        ASSERT(r->type == EXPR_INTEGER);
        r_div = r->data.integer;
        expr_free(r);
        expr_free(p);
    }
    {
        Expr* p = parse_expression("MatrixRank[{{1, 2, 3}, {4, 5, 6}, "
                                  "{7, 8, 9}}, Method -> \"OneStepRowReduction\"]");
        Expr* r = evaluate(p);
        ASSERT(r->type == EXPR_INTEGER);
        r_one = r->data.integer;
        expr_free(r);
        expr_free(p);
    }
    {
        Expr* p = parse_expression("MatrixRank[{{1, 2, 3}, {4, 5, 6}, "
                                  "{7, 8, 9}}, Method -> \"CofactorExpansion\"]");
        Expr* r = evaluate(p);
        ASSERT(r->type == EXPR_INTEGER);
        r_cof = r->data.integer;
        expr_free(r);
        expr_free(p);
    }
    {
        Expr* p = parse_expression("MatrixRank[{{1, 2, 3}, {4, 5, 6}, "
                                  "{7, 8, 9}}]");
        Expr* r = evaluate(p);
        ASSERT(r->type == EXPR_INTEGER);
        r_auto = r->data.integer;
        expr_free(r);
        expr_free(p);
    }
    ASSERT(r_div == 2);
    ASSERT(r_one == 2);
    ASSERT(r_cof == 2);
    ASSERT(r_auto == 2);
    printf("  PASS: all four methods agree on rank=2\n");
}

/* Repeated calls -- exercise that no state accumulates and no leaks
 * occur over many invocations. */
static void test_matrank_repeated(void) {
    for (int i = 0; i < 25; i++) {
        Expr* p = parse_expression(
            "MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]");
        Expr* r = evaluate(p);
        ASSERT(r->type == EXPR_INTEGER);
        ASSERT(r->data.integer == 2);
        expr_free(r);
        expr_free(p);
    }
    printf("  PASS: 25 repeated MatrixRank calls\n");
}

/* Docstring registered via Information[MatrixRank]. */
static void test_matrank_docstring(void) {
    Expr* p = parse_expression("Information[MatrixRank]");
    Expr* r = evaluate(p);
    expr_free(r);
    expr_free(p);
    printf("  PASS: Information[MatrixRank] evaluated\n");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running MatrixRank tests...\n");

    /* Exact-form. */
    TEST(test_matrank_3x3_rank2);
    TEST(test_matrank_3x3_full);
    TEST(test_matrank_symbolic_full);
    TEST(test_matrank_symbolic_2x2_dep);
    TEST(test_matrank_symbolic_2x2_indep);
    TEST(test_matrank_3x5_rect);
    TEST(test_matrank_2x4_rect);
    TEST(test_matrank_zero_matrix);
    TEST(test_matrank_identity);
    TEST(test_matrank_diag);
    TEST(test_matrank_rational_full);
    TEST(test_matrank_rational_dep);
    TEST(test_matrank_4x4_rank3);
    TEST(test_matrank_bigint);
    TEST(test_matrank_bigint_indep);
    TEST(test_matrank_1x3);
    TEST(test_matrank_3x1);

    /* Method. */
    TEST(test_matrank_method_divfree);
    TEST(test_matrank_method_onestep);
    TEST(test_matrank_method_cofactor);
    TEST(test_matrank_method_automatic_symbol);
    TEST(test_matrank_method_automatic_string);
    TEST(test_matrank_method_unknown);

    /* Tolerance / numerical. */
    TEST(test_matrank_real_full);
    TEST(test_matrank_real_singular);
    TEST(test_matrank_complex_dep);
    TEST(test_matrank_tolerance_spec);
    TEST(test_matrank_exact_with_tolerance);
    TEST(test_matrank_tolerance_rational);
    TEST(test_matrank_tolerance_invalid);
    TEST(test_matrank_method_and_tolerance);

    /* Errors. */
    TEST(test_matrank_arity);
    TEST(test_matrank_non_matrix);

    /* Stress / docs. */
    TEST(test_matrank_method_agreement);
    TEST(test_matrank_repeated);
    TEST(test_matrank_docstring);

    printf("All MatrixRank tests passed!\n");
    symtab_clear();
    return 0;
}
