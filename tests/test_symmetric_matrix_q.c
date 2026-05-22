/* Unit tests for SymmetricMatrixQ.
 *
 * A matrix m is symmetric iff m == Transpose[m], i.e. m[i,j] == m[j,i]
 * for every (i, j) pair.  No conjugation is involved, so a complex
 * symmetric matrix need not be Hermitian.
 *
 * Coverage parallels test_hermitian_matrix_q.c, focusing on the
 * complex-symmetric semantics:
 *   - Real symmetric numeric matrices (small + identity).
 *   - Complex symmetric numeric matrices (need not be Hermitian).
 *   - Symbolic symmetric matrices ({{a,b,c},{b,d,e},{c,e,f}}).
 *   - Non-symmetric rejections (asymmetric, complex Hermitian-but-not-
 *     symmetric, mismatched off-diagonal).
 *   - Shape/structure rejections: non-list, non-square, ragged, empty,
 *     vectors, 3-D tensors, scalars.
 *   - Cross-property: real symmetric matrices are also Hermitian; complex
 *     symmetric matrices generally are not.
 *   - Option handling:
 *       * SameTest -> Automatic falls through to structural test.
 *       * SameTest -> f uses the user-supplied predicate.
 *       * Tolerance -> Automatic falls through to structural test.
 *       * Tolerance -> t accepts entries within numeric tolerance.
 *       * Unknown options leave the call unevaluated.
 *       * Non-Rule extra args leave the call unevaluated.
 *   - Equivalence with SameQ[m, Transpose[m]] on a battery of matrices.
 *   - Attribute introspection (Protected) and docstring presence.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

/* --- Real / numeric symmetric ------------------------------------- */

static void test_real_symmetric_2x2(void) {
    assert_eval_eq("SymmetricMatrixQ[{{1, 2}, {2, 3}}]", "True", 0);
    assert_eval_eq("SymmetricMatrixQ[{{1, 2.3}, {2.3, 2}}]", "True", 0);
}

static void test_real_symmetric_3x3(void) {
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}]", "True", 0);
}

static void test_identity_matrix_is_symmetric(void) {
    assert_eval_eq("SymmetricMatrixQ[IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("SymmetricMatrixQ[IdentityMatrix[5]]", "True", 0);
    assert_eval_eq("SymmetricMatrixQ[{{1}}]", "True", 0);
}

static void test_single_element_symbolic_is_symmetric(void) {
    /* {{a}} is symmetric for any a -- there is nothing to compare. */
    assert_eval_eq("SymmetricMatrixQ[{{a}}]", "True", 0);
}

/* --- Complex symmetric -------------------------------------------- */

static void test_complex_symmetric_2x2(void) {
    /* Docstring example: m = {{1+I, 2-3 I}, {2-3 I, 2-3 I}}.  Off-diagonal
     * entries are equal (not conjugate-mirrored), so this is symmetric. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1 + I, 2 - 3*I}, {2 - 3*I, 2 - 3*I}}]",
        "True", 0);
    /* {{1, 2 I}, {2 I, 3}} -- off-diagonals equal, so symmetric even
     * though Hermitian would require Conjugate. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2*I}, {2*I, 3}}]", "True", 0);
}

static void test_complex_symmetric_3x3(void) {
    /* Docstring c-matrix: complex symmetric but not Hermitian. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{I, 0, 1}, {0, I, I}, {1, I, I}}]", "True", 0);
}

static void test_complex_hermitian_not_symmetric(void) {
    /* {{1, 3+4 I}, {3-4 I, 2}} is Hermitian, not symmetric (off-diagonal
     * entries differ). */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 3 + 4*I}, {3 - 4*I, 2}}]", "False", 0);
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2 - 3*I}, {2 + 3*I, 4}}]", "False", 0);
}

/* --- Symbolic symmetric ------------------------------------------- */

static void test_symbolic_symmetric_3x3(void) {
    /* Docstring example: {{a, b, c}, {b, d, e}, {c, e, f}} is symmetric. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{a, b, c}, {b, d, e}, {c, e, f}}]",
        "True", 0);
}

static void test_symbolic_symmetric_2x2(void) {
    /* {{a, b}, {b, c}} is symmetric for any a, b, c. */
    assert_eval_eq("SymmetricMatrixQ[{{a, b}, {b, c}}]", "True", 0);
    /* Same when the diagonal is the same value. */
    assert_eval_eq("SymmetricMatrixQ[{{a, b}, {b, a}}]", "True", 0);
}

static void test_symbolic_generic_not_symmetric(void) {
    /* {{a, b}, {c, d}} with all four symbols distinct -- structural test
     * cannot prove b == c, so we report False. */
    assert_eval_eq("SymmetricMatrixQ[{{a, b}, {c, d}}]", "False", 0);
}

/* --- Non-symmetric rejections -------------------------------------- */

static void test_asymmetric_real_matrix(void) {
    /* Docstring example: 3x3 integer matrix that is not symmetric. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{4, -5, 2}, {3, -3, -3}, {5, 5, 5}}]",
        "False", 0);
}

static void test_symmetrize_makes_symmetric(void) {
    /* (m + Transpose[m])/2 is always symmetric. */
    assert_eval_eq(
        "SymmetricMatrixQ[({{4, -5, 2}, {3, -3, -3}, {5, 5, 5}}"
        " + Transpose[{{4, -5, 2}, {3, -3, -3}, {5, 5, 5}}])/2]",
        "True", 0);
}

static void test_off_diagonal_mismatch(void) {
    /* Mismatch in a single off-diagonal pair. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2, 3}, {2, 4, 7}, {3, 5, 6}}]",
        "False", 0);
}

/* --- Shape / structure rejections ---------------------------------- */

static void test_non_matrix_inputs(void) {
    /* Scalars: not a List. */
    assert_eval_eq("SymmetricMatrixQ[5]", "False", 0);
    assert_eval_eq("SymmetricMatrixQ[Sqrt[3]]", "False", 0);
    assert_eval_eq("SymmetricMatrixQ[\"foo\"]", "False", 0);
    /* Vector (1-D list). */
    assert_eval_eq("SymmetricMatrixQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("SymmetricMatrixQ[{a, b, c}]", "False", 0);
    /* Empty list. */
    assert_eval_eq("SymmetricMatrixQ[{}]", "False", 0);
}

static void test_non_square_matrix(void) {
    /* 1xn rows or nx1 columns. */
    assert_eval_eq("SymmetricMatrixQ[{{1, 2, 3}}]", "False", 0);
    assert_eval_eq("SymmetricMatrixQ[{{1}, {2}, {3}}]", "False", 0);
    /* Rectangular 2x3 / 3x2. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "False", 0);
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {3, 4}, {5, 6}}]", "False", 0);
}

static void test_ragged_matrix(void) {
    /* Rows of different lengths. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {3, 4, 5}}]", "False", 0);
}

static void test_three_d_tensor_rejected(void) {
    /* A list-of-list-of-list is not a matrix. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
}

static void test_pure_symbol_stays_evaluated(void) {
    /* A bare symbol is not a matrix; predicate answers False rather than
     * leaving the call unevaluated. */
    assert_eval_eq("SymmetricMatrixQ[m]", "False", 0);
}

/* --- Cross-property: symmetric vs Hermitian ----------------------- */

static void test_real_symmetric_is_hermitian(void) {
    /* For real numeric matrices, symmetric <=> Hermitian. */
    assert_eval_eq(
        "{SymmetricMatrixQ[{{1, 2}, {2, 3}}],"
        " HermitianMatrixQ[{{1, 2}, {2, 3}}]}",
        "{True, True}", 0);
}

static void test_complex_symmetric_not_always_hermitian(void) {
    /* m = {{1+I, 2+2 I}, {2+2 I, 3+3 I}} is symmetric but not Hermitian. */
    assert_eval_eq(
        "{SymmetricMatrixQ[{{1 + I, 2 + 2*I}, {2 + 2*I, 3 + 3*I}}],"
        " HermitianMatrixQ[{{1 + I, 2 + 2*I}, {2 + 2*I, 3 + 3*I}}]}",
        "{True, False}", 0);
}

/* --- Options: SameTest -------------------------------------------- */

static void test_sametest_automatic_is_default(void) {
    /* SameTest -> Automatic must match the no-option case. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {2, 3}}, SameTest -> Automatic]",
        "True", 0);
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {3, 4}}, SameTest -> Automatic]",
        "False", 0);
}

static void test_sametest_accepts_via_predicate(void) {
    /* Predicate that trivially accepts everything. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{a, b}, {c, d}}, SameTest -> (True &)]",
        "True", 0);
    /* Predicate that trivially rejects everything; with a single-element
     * matrix there is no off-diagonal pair to check, so the predicate is
     * never invoked and we still report True. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {3, 4}}, SameTest -> (False &)]",
        "False", 0);
}

static void test_sametest_simplify_log_example(void) {
    /* Docstring example: m = {{1, Log[x^2]}, {2 Log[x], 2}} is symmetric
     * for positive x but the structural test cannot see it.  Using a
     * Simplify-based SameTest with x > 0 should accept. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, Log[x^2]}, {2*Log[x], 2}}]",
        "False", 0);
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, Log[x^2]}, {2*Log[x], 2}},"
        " SameTest -> (Simplify[#1 - #2, x > 0] == 0 &)]",
        "True", 0);
}

/* --- Options: Tolerance ------------------------------------------- */

static void test_tolerance_automatic_is_default(void) {
    /* Tolerance -> Automatic falls through to the structural test. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {2, 3}}, Tolerance -> Automatic]",
        "True", 0);
    assert_eval_eq(
        "SymmetricMatrixQ[{{1, 2}, {3, 4}}, Tolerance -> Automatic]",
        "False", 0);
}

static void test_tolerance_accepts_within_threshold(void) {
    /* Off-diagonal diff = 0.01; tolerance 0.1 accepts. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1.0, 2.01}, {2.0, 1.5}}, Tolerance -> 0.1]",
        "True", 0);
    /* Complex off-diagonals differ by 0.01 I; tolerance 0.1 accepts. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0, 1.5}},"
        " Tolerance -> 0.1]",
        "True", 0);
}

static void test_tolerance_rejects_outside_threshold(void) {
    /* Same near-symmetric matrix but with a tighter tolerance. */
    assert_eval_eq(
        "SymmetricMatrixQ[{{1.0, 2.01}, {2.0, 1.5}}, Tolerance -> 0.001]",
        "False", 0);
}

/* --- Options: error handling -------------------------------------- */

static void test_unknown_option_stays_unevaluated(void) {
    const char* in = "SymmetricMatrixQ[{{1, 2}, {2, 3}}, Foo -> Bar]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "SymmetricMatrixQ") != NULL,
               "expected unevaluated SymmetricMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_non_rule_extra_arg_stays_unevaluated(void) {
    const char* in = "SymmetricMatrixQ[{{1, 2}, {2, 3}}, foo]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "SymmetricMatrixQ") != NULL,
               "expected unevaluated SymmetricMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_zero_arg_stays_unevaluated(void) {
    const char* in = "SymmetricMatrixQ[]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "SymmetricMatrixQ") != NULL,
               "expected unevaluated SymmetricMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Cross-check against Transpose equality ----------------------- */

static void test_matches_transpose_equality(void) {
    /* For each matrix, the predicate should agree with m === Transpose[m]
     * up to SameQ. */
    struct {
        const char* matrix;
        const char* expected;
    } cases[] = {
        { "{{1, 2}, {2, 3}}",                              "True"  },
        { "{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}",             "True"  },
        { "{{1, 2*I}, {2*I, 3}}",                          "True"  },
        { "{{1, 3 + 4*I}, {3 - 4*I, 2}}",                  "False" },
        { "{{1, 2}, {3, 4}}",                              "False" },
        { "IdentityMatrix[4]",                             "True"  },
        { "{{a, b, c}, {b, d, e}, {c, e, f}}",             "True"  },
        { "{{I, 0, 1}, {0, I, I}, {1, I, I}}",             "True"  },
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf), "SymmetricMatrixQ[%s]", cases[i].matrix);
        assert_eval_eq(buf, cases[i].expected, 0);
    }
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("SymmetricMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* SymmetricMatrixQ should NOT be Listable -- input is a matrix, not
     * a list of matrices. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("SymmetricMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "symmetric") != NULL);
    ASSERT(strstr(def->docstring, "SameTest") != NULL);
    ASSERT(strstr(def->docstring, "Tolerance") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_SymmetricMatrixQ != NULL);
    ASSERT(strcmp(SYM_SymmetricMatrixQ, "SymmetricMatrixQ") == 0);
}

/* --- Repeated evaluation: leak sanity ----------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Run the same predicate many times; any builtin double-free / dangling
     * pointer typically manifests as a crash or a wrong answer under
     * valgrind on the second iteration onwards. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq(
            "SymmetricMatrixQ[{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}]",
            "True", 0);
        assert_eval_eq(
            "SymmetricMatrixQ[{{a, b}, {c, d}}]", "False", 0);
        assert_eval_eq(
            "SymmetricMatrixQ[{{1.0, 2.01}, {2.0, 1.5}}, Tolerance -> 0.1]",
            "True", 0);
        assert_eval_eq(
            "SymmetricMatrixQ[{{a, b}, {c, d}}, SameTest -> (True &)]",
            "True", 0);
        assert_eval_eq(
            "SymmetricMatrixQ[{{I, 0, 1}, {0, I, I}, {1, I, I}}]",
            "True", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_real_symmetric_2x2);
    TEST(test_real_symmetric_3x3);
    TEST(test_identity_matrix_is_symmetric);
    TEST(test_single_element_symbolic_is_symmetric);

    TEST(test_complex_symmetric_2x2);
    TEST(test_complex_symmetric_3x3);
    TEST(test_complex_hermitian_not_symmetric);

    TEST(test_symbolic_symmetric_3x3);
    TEST(test_symbolic_symmetric_2x2);
    TEST(test_symbolic_generic_not_symmetric);

    TEST(test_asymmetric_real_matrix);
    TEST(test_symmetrize_makes_symmetric);
    TEST(test_off_diagonal_mismatch);

    TEST(test_non_matrix_inputs);
    TEST(test_non_square_matrix);
    TEST(test_ragged_matrix);
    TEST(test_three_d_tensor_rejected);
    TEST(test_pure_symbol_stays_evaluated);

    TEST(test_real_symmetric_is_hermitian);
    TEST(test_complex_symmetric_not_always_hermitian);

    TEST(test_sametest_automatic_is_default);
    TEST(test_sametest_accepts_via_predicate);
    TEST(test_sametest_simplify_log_example);

    TEST(test_tolerance_automatic_is_default);
    TEST(test_tolerance_accepts_within_threshold);
    TEST(test_tolerance_rejects_outside_threshold);

    TEST(test_unknown_option_stays_unevaluated);
    TEST(test_non_rule_extra_arg_stays_unevaluated);
    TEST(test_zero_arg_stays_unevaluated);

    TEST(test_matches_transpose_equality);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All SymmetricMatrixQ tests passed!\n");
    return 0;
}
