#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "Differences %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "Differences %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Spec examples ---------- */

static void test_first_differences_symbolic(void) {
    /* Differences[{a, b, c, d, e}] -> {-a + b, -b + c, -c + d, -d + e} */
    run_infix("Differences[{a, b, c, d, e}]",
              "{-a + b, -b + c, -c + d, -d + e}");
}

static void test_second_differences_symbolic(void) {
    /* Differences[{a, b, c, d, e}, 2] -> {a - 2 b + c, b - 2 c + d, c - 2 d + e} */
    run_infix("Differences[{a, b, c, d, e}, 2]",
              "{a - 2 b + c, b - 2 c + d, c - 2 d + e}");
}

static void test_first_differences_numeric(void) {
    /* Differences[{1, 4, 9, 16, 25}] -> {3, 5, 7, 9} */
    run_full("Differences[{1, 4, 9, 16, 25}]", "List[3, 5, 7, 9]");
}

static void test_linear_constant_first_difference(void) {
    /* First differences constant for a linear function (3 i + 6). */
    run_full("Differences[Table[3 i + 6, {i, 10}]]",
             "List[3, 3, 3, 3, 3, 3, 3, 3, 3]");
}

static void test_quadratic_constant_second_difference(void) {
    /* Second differences constant for a quadratic function (3 i^2 + 6). */
    run_full("Differences[Table[3 i^2 + 6, {i, 10}], 2]",
             "List[6, 6, 6, 6, 6, 6, 6, 6]");
}

static void test_degree_sequence(void) {
    /* n^5 + 2 n - 1 has constant 5th differences (= 120) and zero 6th. */
    run_full("Differences[Table[n^5 + 2 n - 1, {n, 8}], 5]",
             "List[120, 120, 120]");
    run_full("Differences[Table[n^5 + 2 n - 1, {n, 8}], 6]",
             "List[0, 0]");
}

static void test_progressively_shorter(void) {
    /* Differences[Range[10], n] is all zeros for n >= 2, length 10 - n. */
    run_full("Differences[Range[10], 1]",
             "List[1, 1, 1, 1, 1, 1, 1, 1, 1]");
    run_full("Differences[Range[10], 9]", "List[0]");
    run_full("Differences[Range[10], 10]", "List[]");
}

/* ---------- Step argument (third argument s) ---------- */

static void test_step_two(void) {
    /* Differences[list, 1, 2] subtracts elements two apart. */
    run_full("Differences[{1, 2, 4, 8, 16, 32}, 1, 2]",
             "List[3, 6, 12, 24]");
}

static void test_step_length_property(void) {
    /* Differences[list, n, s] has length l - n |s|. */
    run_full("Length[Differences[Range[20], 2, 3]]", "14");
    run_full("Length[Differences[Range[20], 1, 5]]", "15");
}

static void test_negative_step(void) {
    /* For s < 0, the i-th element is elem[i] - elem[i + |s|]; on Range this
       is the negation of the positive-step result. */
    run_full("Differences[{1, 2, 4, 8, 16}, 1, -2]",
             "List[-3, -6, -12]");
}

/* ---------- Nested lists / level specifications ---------- */

static void test_matrix_row_differences_default(void) {
    /* Differences[m] on a matrix differences successive rows (within columns). */
    run_full("Differences[{{1, 2, 3}, {4, 6, 8}, {9, 12, 15}}]",
             "List[List[3, 4, 5], List[5, 6, 7]]");
}

static void test_matrix_row_differences_levelspec(void) {
    /* Differences[m, {1, 0}] is equivalent to Differences[m] / Differences[m, 1]. */
    run_full("Differences[{{1, 2, 3}, {4, 6, 8}, {9, 12, 15}}, {1, 0}]",
             "List[List[3, 4, 5], List[5, 6, 7]]");
}

static void test_matrix_column_differences(void) {
    /* Differences[m, {0, 1}] differences columns within each row. */
    run_full("Differences[{{1, 2, 3}, {4, 6, 8}, {9, 12, 15}}, {0, 1}]",
             "List[List[1, 1], List[2, 2], List[3, 3]]");
}

static void test_matrix_mixed_levels(void) {
    /* Differences[m, {1, 1}] == Differences[Differences[m, 1], {0, 1}]. */
    run_full("Differences[{{1, 2, 3}, {4, 6, 8}, {9, 12, 15}}, {1, 1}]",
             "List[List[1, 1], List[1, 1]]");
}

static void test_matrix_symbolic_columns(void) {
    /* Column differences of a symbolic matrix (spec example, {0,1}). */
    run_infix("Differences[{{a11, a12, a13}, {a21, a22, a23}}, {0, 1}]",
              "{{-a11 + a12, -a12 + a13}, {-a21 + a22, -a22 + a23}}");
}

/* ---------- Inverse relationship with FoldList ---------- */

static void test_foldlist_inverse(void) {
    /* FoldList[Plus, x, Differences[list]] reconstructs the original list. */
    run_full("FoldList[Plus, 1, Differences[{1, 4, 9, 16, 25}]]",
             "List[1, 4, 9, 16, 25]");
}

static void test_foldlist_inverse_symbolic(void) {
    run_infix("FoldList[Plus, a, Differences[{a, b, c, d, e}]]",
              "{a, b, c, d, e}");
}

/* ---------- Head preservation & length ---------- */

static void test_length_property(void) {
    /* Differences[list, n] has length l - n. */
    run_full("Length[Differences[Range[12], 4]]", "8");
    run_full("Length[Differences[{a, b, c, d, e}]]", "4");
}

static void test_n_zero_identity(void) {
    /* Differences[list, 0] returns the list unchanged. */
    run_full("Differences[{a, b, c}, 0]", "List[a, b, c]");
}

static void test_head_preserved(void) {
    /* The result keeps the head of the input. */
    run_full("Differences[f[1, 4, 9, 16]]", "f[3, 5, 7]");
}

/* ---------- Edge cases ---------- */

static void test_singleton(void) {
    run_full("Differences[{5}]", "List[]");
}

static void test_empty(void) {
    run_full("Differences[{}]", "List[]");
}

static void test_n_exceeds_length(void) {
    run_full("Differences[{1, 2, 3}, 5]", "List[]");
}

static void test_bignum(void) {
    /* Exact arbitrary-precision arithmetic flows through. */
    run_full("Differences[{2^100, 2^101, 2^102}]",
             "List[1267650600228229401496703205376, "
             "2535301200456458802993406410752]");
}

static void test_real(void) {
    run_full("Differences[{1.5, 2.5, 4.0}]", "List[1.0, 1.5]");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_nonlist_unevaluated(void) {
    run_full("Differences[x]", "Differences[x]");
}

static void test_wrong_arg_count(void) {
    run_full("Differences[]", "Differences[]");
    run_full("Differences[{1, 2, 3}, 1, 2, 4]",
             "Differences[List[1, 2, 3], 1, 2, 4]");
}

static void test_bad_step_zero(void) {
    /* Step 0 is invalid -> unevaluated. */
    run_full("Differences[{1, 2, 3}, 1, 0]",
             "Differences[List[1, 2, 3], 1, 0]");
}

static void test_negative_order_unevaluated(void) {
    run_full("Differences[{1, 2, 3}, -1]",
             "Differences[List[1, 2, 3], -1]");
}

static void test_symbolic_order_unevaluated(void) {
    run_full("Differences[{1, 2, 3}, n]",
             "Differences[List[1, 2, 3], n]");
}

static void test_empty_levelspec_identity(void) {
    run_full("Differences[{1, 2, 3}, {}]", "List[1, 2, 3]");
}

static void test_attributes_protected(void) {
    run_full("MemberQ[Attributes[Differences], Protected]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_first_differences_symbolic);
    TEST(test_second_differences_symbolic);
    TEST(test_first_differences_numeric);
    TEST(test_linear_constant_first_difference);
    TEST(test_quadratic_constant_second_difference);
    TEST(test_degree_sequence);
    TEST(test_progressively_shorter);

    /* Step argument */
    TEST(test_step_two);
    TEST(test_step_length_property);
    TEST(test_negative_step);

    /* Nested lists / level specs */
    TEST(test_matrix_row_differences_default);
    TEST(test_matrix_row_differences_levelspec);
    TEST(test_matrix_column_differences);
    TEST(test_matrix_mixed_levels);
    TEST(test_matrix_symbolic_columns);

    /* FoldList inverse */
    TEST(test_foldlist_inverse);
    TEST(test_foldlist_inverse_symbolic);

    /* Head & length */
    TEST(test_length_property);
    TEST(test_n_zero_identity);
    TEST(test_head_preserved);

    /* Edge cases */
    TEST(test_singleton);
    TEST(test_empty);
    TEST(test_n_exceeds_length);
    TEST(test_bignum);
    TEST(test_real);

    /* Unevaluated / error cases */
    TEST(test_nonlist_unevaluated);
    TEST(test_wrong_arg_count);
    TEST(test_bad_step_zero);
    TEST(test_negative_order_unevaluated);
    TEST(test_symbolic_order_unevaluated);
    TEST(test_empty_levelspec_identity);
    TEST(test_attributes_protected);

    printf("All Differences tests passed!\n");
    return 0;
}
