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
               "Ratios %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "Ratios %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Spec examples ---------- */

static void test_first_ratios_symbolic(void) {
    /* Ratios[{a, b, c, d, e}] -> {b/a, c/b, d/c, e/d} */
    run_infix("Ratios[{a, b, c, d, e}]",
              "{b/a, c/b, d/c, e/d}");
}

static void test_second_ratios_symbolic(void) {
    /* Ratios[{a, b, c, d, e}, 2] -> {(a c)/b^2, (b d)/c^2, (c e)/d^2} */
    run_infix("Ratios[{a, b, c, d, e}, 2]",
              "{(a c)/b^2, (b d)/c^2, (c e)/d^2}");
}

static void test_third_ratios_constant(void) {
    /* Third ratios of a geometric-ish run are 1 once exhausted. */
    run_full("Ratios[{1, 2, 4, 8, 16, 32}, 3]", "List[1, 1, 1]");
}

static void test_geometric_sequence_constant(void) {
    /* First ratios are constant (= 2) for a geometric sequence 2^i. */
    run_full("Ratios[Table[2^i, {i, 0, 10}]]",
             "List[2, 2, 2, 2, 2, 2, 2, 2, 2, 2]");
}

static void test_first_ratios_rational(void) {
    /* Exact rationals flow through. */
    run_full("Ratios[{2, 3, 4}]", "List[Rational[3, 2], Rational[4, 3]]");
    run_infix("Ratios[{2, 3, 4}]", "{3/2, 4/3}");
}

/* ---------- Iterated order argument ---------- */

static void test_n_zero_identity(void) {
    /* Ratios[list, 0] returns the list unchanged. */
    run_full("Ratios[{a, b, c}, 0]", "List[a, b, c]");
}

static void test_length_property(void) {
    /* Ratios[list, n] has length l - n. */
    run_full("Length[Ratios[Range[12], 4]]", "8");
    run_full("Length[Ratios[{a, b, c, d, e}]]", "4");
    run_full("Length[Ratios[Range[10], 1]]", "9");
}

static void test_progressively_shorter(void) {
    /* On a geometric run, every iterated ratio level is all 1s. */
    run_full("Ratios[{1, 2, 4, 8, 16}, 1]", "List[2, 2, 2, 2]");
    run_full("Ratios[{1, 2, 4, 8, 16}, 2]", "List[1, 1, 1]");
    run_full("Ratios[{1, 2, 4, 8, 16}, 4]", "List[1]");
}

/* ---------- Nested lists / level specifications ---------- */

static void test_matrix_row_ratios_default(void) {
    /* Ratios[m] on a matrix takes ratios of successive rows (within columns). */
    run_infix("Ratios[{{a11, a12, a13}, {a21, a22, a23}, {a31, a32, a33}}]",
              "{{a21/a11, a22/a12, a23/a13}, {a31/a21, a32/a22, a33/a23}}");
}

static void test_matrix_row_ratios_levelspec(void) {
    /* Ratios[m, {1, 0}] is equivalent to Ratios[m] / Ratios[m, 1]. */
    run_infix("Ratios[{{a11, a12, a13}, {a21, a22, a23}}, {1, 0}]",
              "{{a21/a11, a22/a12, a23/a13}}");
}

static void test_matrix_column_ratios(void) {
    /* Ratios[m, {0, 1}] takes ratios of columns within each row. */
    run_infix("Ratios[{{a11, a12, a13}, {a21, a22, a23}, {a31, a32, a33}}, {0, 1}]",
              "{{a12/a11, a13/a12}, {a22/a21, a23/a22}, {a32/a31, a33/a32}}");
}

static void test_matrix_mixed_levels_equivalence(void) {
    /* Ratios[m, {n1, n2}] == Ratios[Ratios[m, n1], {0, n2}]. */
    run_infix("Ratios[{{a, b}, {c, d}, {e, f}}, {1, 1}]",
              "{{(a d)/(b c)}, {(c f)/(d e)}}");
    run_infix("Ratios[Ratios[{{a, b}, {c, d}, {e, f}}, 1], {0, 1}]",
              "{{(a d)/(b c)}, {(c f)/(d e)}}");
}

static void test_matrix_numeric_columns(void) {
    run_full("Ratios[{{2, 4, 8}, {3, 9, 27}}, {0, 1}]",
             "List[List[2, 2], List[3, 3]]");
}

/* ---------- Inverse relationship with FoldList ---------- */

static void test_foldlist_inverse(void) {
    /* FoldList[Times, x, Ratios[list]] reconstructs the original list. */
    run_full("FoldList[Times, 1, Ratios[{1, 2, 4, 8, 16}]]",
             "List[1, 2, 4, 8, 16]");
}

static void test_foldlist_inverse_symbolic(void) {
    run_infix("FoldList[Times, a, Ratios[{a, b, c, d, e}]]",
              "{a, b, c, d, e}");
}

/* ---------- Head preservation ---------- */

static void test_head_preserved(void) {
    /* The result keeps the head of the input. */
    run_full("Ratios[f[1, 2, 4, 8]]", "f[2, 2, 2]");
}

/* ---------- Edge cases ---------- */

static void test_singleton(void) {
    run_full("Ratios[{5}]", "List[]");
}

static void test_empty(void) {
    run_full("Ratios[{}]", "List[]");
}

static void test_n_exceeds_length(void) {
    run_full("Ratios[{1, 2, 3}, 5]", "List[]");
}

static void test_bignum(void) {
    /* Exact arbitrary-precision arithmetic flows through (powers of two -> 2). */
    run_full("Ratios[{2^100, 2^101, 2^102}]", "List[2, 2]");
}

static void test_real(void) {
    run_full("Ratios[{1.0, 2.0, 6.0}]", "List[2.0, 3.0]");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_nonlist_unevaluated(void) {
    run_full("Ratios[x]", "Ratios[x]");
}

static void test_wrong_arg_count(void) {
    run_full("Ratios[]", "Ratios[]");
    run_full("Ratios[{1, 2, 3}, 1, 2]",
             "Ratios[List[1, 2, 3], 1, 2]");
}

static void test_negative_order_unevaluated(void) {
    run_full("Ratios[{1, 2, 3}, -1]",
             "Ratios[List[1, 2, 3], -1]");
}

static void test_symbolic_order_unevaluated(void) {
    run_full("Ratios[{1, 2, 3}, n]",
             "Ratios[List[1, 2, 3], n]");
}

static void test_empty_levelspec_identity(void) {
    run_full("Ratios[{1, 2, 3}, {}]", "List[1, 2, 3]");
}

static void test_negative_level_unevaluated(void) {
    run_full("Ratios[{{1, 2}, {3, 4}}, {1, -1}]",
             "Ratios[List[List[1, 2], List[3, 4]], List[1, -1]]");
}

static void test_attributes_protected(void) {
    run_full("MemberQ[Attributes[Ratios], Protected]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_first_ratios_symbolic);
    TEST(test_second_ratios_symbolic);
    TEST(test_third_ratios_constant);
    TEST(test_geometric_sequence_constant);
    TEST(test_first_ratios_rational);

    /* Iterated order */
    TEST(test_n_zero_identity);
    TEST(test_length_property);
    TEST(test_progressively_shorter);

    /* Nested lists / level specs */
    TEST(test_matrix_row_ratios_default);
    TEST(test_matrix_row_ratios_levelspec);
    TEST(test_matrix_column_ratios);
    TEST(test_matrix_mixed_levels_equivalence);
    TEST(test_matrix_numeric_columns);

    /* FoldList inverse */
    TEST(test_foldlist_inverse);
    TEST(test_foldlist_inverse_symbolic);

    /* Head preservation */
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
    TEST(test_negative_order_unevaluated);
    TEST(test_symbolic_order_unevaluated);
    TEST(test_empty_levelspec_identity);
    TEST(test_negative_level_unevaluated);
    TEST(test_attributes_protected);

    printf("All Ratios tests passed!\n");
    return 0;
}
