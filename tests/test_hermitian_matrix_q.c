/* Unit tests for HermitianMatrixQ.
 *
 * A matrix m is Hermitian (self-adjoint) iff m == ConjugateTranspose[m].
 * These tests exercise the *evaluator* path (parse -> evaluate -> print)
 * so attribute application and option dispatch are all in scope.
 *
 * Coverage:
 *   - Real symmetric numeric matrices (small + identity).
 *   - Complex Hermitian numeric matrices (Cartesian and structural).
 *   - Symbolic matrices containing Conjugate[...] patterns.
 *   - Non-Hermitian rejections (asymmetric, complex symmetric, mismatched
 *     diagonal, off-diagonal sign error).
 *   - Diagonal-must-be-real rejection (e.g. {{1+I}}).
 *   - Shape/structure rejections: non-list, non-square, ragged, empty,
 *     vectors, 3-D tensors, scalars.
 *   - Option handling:
 *       * SameTest -> Automatic falls through to structural test.
 *       * SameTest -> f uses the user-supplied predicate.
 *       * Tolerance -> Automatic falls through to structural test.
 *       * Tolerance -> t accepts entries within numeric tolerance.
 *       * Unknown options leave the call unevaluated.
 *       * Non-Rule extra args leave the call unevaluated.
 *   - Equivalence with SameQ[m, ConjugateTranspose[m]] on a battery of
 *     numeric matrices.
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

/* --- Real / numeric Hermitian -------------------------------------- */

static void test_real_symmetric_2x2(void) {
    assert_eval_eq("HermitianMatrixQ[{{1, 2.3}, {2.3, 4}}]", "True", 0);
    assert_eval_eq("HermitianMatrixQ[{{1, 2}, {2, 3}}]", "True", 0);
}

static void test_real_symmetric_3x3(void) {
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}]", "True", 0);
}

static void test_identity_matrix_is_hermitian(void) {
    assert_eval_eq("HermitianMatrixQ[IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("HermitianMatrixQ[IdentityMatrix[5]]", "True", 0);
    assert_eval_eq("HermitianMatrixQ[{{1}}]", "True", 0);
}

static void test_complex_hermitian_2x2(void) {
    /* m = {{1, 3+4 I}, {3-4 I, 2}}.  Conjugate[3-4 I] = 3+4 I. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 3 + 4*I}, {3 - 4*I, 2}}]", "True", 0);
    /* Same with simpler 2 - 3 I / 2 + 3 I pair. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2 - 3*I}, {2 + 3*I, 4}}]", "True", 0);
}

static void test_complex_hermitian_3x3(void) {
    /* Diagonal real (1, 0, -1), off-diagonals are mirror-conjugates. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2 + I, 3}, {2 - I, 0, 4 - 2*I},"
        " {3, 4 + 2*I, -1}}]",
        "True", 0);
}

/* --- Symbolic Hermitian via Conjugate ------------------------------ */

static void test_symbolic_hermitian_3x3(void) {
    /* Docstring example: {{0,a,b},{Conjugate[a],1,c},
     *                     {Conjugate[b],Conjugate[c],-1}}. */
    assert_eval_eq(
        "HermitianMatrixQ[{{0, a, b}, {Conjugate[a], 1, c},"
        " {Conjugate[b], Conjugate[c], -1}}]",
        "True", 0);
}

static void test_symbolic_hermitian_2x2(void) {
    /* {{0, a}, {Conjugate[a], 0}} is Hermitian for any a. */
    assert_eval_eq(
        "HermitianMatrixQ[{{0, a}, {Conjugate[a], 0}}]", "True", 0);
    /* Swap argument order; the test must be symmetric. */
    assert_eval_eq(
        "HermitianMatrixQ[{{0, Conjugate[a]}, {a, 0}}]", "True", 0);
}

/* --- Non-Hermitian rejections -------------------------------------- */

static void test_complex_symmetric_not_hermitian(void) {
    /* A complex *symmetric* matrix is not Hermitian (off-diagonals must
     * be conjugate-mirrored, not equal). */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 3 + 4*I}, {3 + 4*I, 5}}]", "False", 0);
    /* {{1, 2 I}, {2 I, 3}}: Conjugate[2 I] = -2 I != 2 I. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2*I}, {2*I, 3}}]", "False", 0);
}

static void test_asymmetric_real_matrix(void) {
    assert_eval_eq(
        "HermitianMatrixQ[{{4, -5, 2}, {3, -3, -3}, {5, 5, 5}}]",
        "False", 0);
}

static void test_symmetrize_makes_hermitian(void) {
    /* For a real matrix, (m + Transpose[m])/2 is symmetric => Hermitian. */
    assert_eval_eq(
        "HermitianMatrixQ[({{4, -5, 2}, {3, -3, -3}, {5, 5, 5}}"
        " + Transpose[{{4, -5, 2}, {3, -3, -3}, {5, 5, 5}}])/2]",
        "True", 0);
}

static void test_symbolic_generic_not_hermitian(void) {
    /* {{a, b}, {c, d}} -- no info about a, b, c, d being conjugates. */
    assert_eval_eq("HermitianMatrixQ[{{a, b}, {c, d}}]", "False", 0);
}

static void test_complex_diagonal_not_hermitian(void) {
    /* The diagonal must be self-conjugate (purely real for numeric). */
    assert_eval_eq("HermitianMatrixQ[{{1 + I}}]", "False", 0);
    assert_eval_eq(
        "HermitianMatrixQ[{{1 + I, 2}, {2, 3}}]", "False", 0);
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2}, {2, 3 - 4*I}}]", "False", 0);
}

static void test_symbolic_diagonal_not_real(void) {
    /* For a single unknown symbol a, Conjugate[a] != a, so {{a}} is
     * not "explicitly" Hermitian. */
    assert_eval_eq("HermitianMatrixQ[{{a}}]", "False", 0);
    assert_eval_eq(
        "HermitianMatrixQ[{{a, Conjugate[b]}, {b, a}}]", "False", 0);
}

/* --- Shape / structure rejections ---------------------------------- */

static void test_non_matrix_inputs(void) {
    /* Scalars: not a List. */
    assert_eval_eq("HermitianMatrixQ[5]", "False", 0);
    assert_eval_eq("HermitianMatrixQ[\"foo\"]", "False", 0);
    /* Vector (1-D list): not a list of lists. */
    assert_eval_eq("HermitianMatrixQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("HermitianMatrixQ[{a, b, c}]", "False", 0);
    /* Empty list. */
    assert_eval_eq("HermitianMatrixQ[{}]", "False", 0);
}

static void test_non_square_matrix(void) {
    /* 1xn rows or nx1 columns -- not square. */
    assert_eval_eq("HermitianMatrixQ[{{1, 2, 3}}]", "False", 0);
    assert_eval_eq("HermitianMatrixQ[{{1}, {2}, {3}}]", "False", 0);
    /* Rectangular 2x3 / 3x2. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2, 3}, {4, 5, 6}}]", "False", 0);
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2}, {3, 4}, {5, 6}}]", "False", 0);
}

static void test_ragged_matrix(void) {
    /* Rows of different lengths. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 2}, {3, 4, 5}}]", "False", 0);
}

static void test_three_d_tensor_rejected(void) {
    /* A list-of-list-of-list is not a matrix even if outer dims square. */
    assert_eval_eq(
        "HermitianMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
}

/* --- Symbolic input that stays a function -------------------------- */

static void test_pure_symbol_stays_evaluated(void) {
    /* Pure symbol or non-list -- the predicate still answers False
     * rather than leaving the call unevaluated.  This matches MatrixQ. */
    assert_eval_eq("HermitianMatrixQ[m]", "False", 0);
}

/* --- Options: SameTest -------------------------------------------- */

static void test_sametest_automatic_is_default(void) {
    /* SameTest -> Automatic must match the no-option case. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 3 + 4*I}, {3 - 4*I, 2}},"
        " SameTest -> Automatic]",
        "True", 0);
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 3 + 4*I}, {3 + 4*I, 2}},"
        " SameTest -> Automatic]",
        "False", 0);
}

static void test_sametest_accepts_via_predicate(void) {
    /* Pair-wise predicate: trivially accept everything. */
    assert_eval_eq(
        "HermitianMatrixQ[{{a, b}, {c, d}},"
        " SameTest -> (True &)]",
        "True", 0);
    /* Pair-wise predicate: trivially reject everything. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1}}, SameTest -> (False &)]",
        "False", 0);
}

static void test_sametest_tolerance_via_predicate(void) {
    /* Custom tolerance predicate accepts a near-Hermitian matrix. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0 - 0.02*I, 1.5}},"
        " SameTest -> (Abs[#1 - #2] < 0.1 &)]",
        "True", 0);
    /* And rejects with a tighter tolerance. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0 - 0.02*I, 1.5}},"
        " SameTest -> (Abs[#1 - #2] < 0.001 &)]",
        "False", 0);
}

/* --- Options: Tolerance ------------------------------------------- */

static void test_tolerance_automatic_is_default(void) {
    /* Tolerance -> Automatic falls through to the structural test. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1, 3 + 4*I}, {3 - 4*I, 2}},"
        " Tolerance -> Automatic]",
        "True", 0);
}

static void test_tolerance_accepts_within_threshold(void) {
    /* Off-diagonal Abs[diff] = 0.01; tolerance 0.1 accepts. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0 - 0.02*I, 1.5}},"
        " Tolerance -> 0.1]",
        "True", 0);
    /* Diagonal imag part 0.005; tolerance 0.01 accepts. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1.0 + 0.005*I, 2.0}, {2.0, 1.5}},"
        " Tolerance -> 0.01]",
        "True", 0);
}

static void test_tolerance_rejects_outside_threshold(void) {
    assert_eval_eq(
        "HermitianMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0 - 0.02*I, 1.5}},"
        " Tolerance -> 0.001]",
        "False", 0);
    /* {{1, 2+0.01I},{2+0.02I, 1.5}}: diff = 0.03 I, abs = 0.03 > 0.001. */
    assert_eval_eq(
        "HermitianMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0 + 0.02*I, 1.5}},"
        " Tolerance -> 0.001]",
        "False", 0);
}

/* --- Options: error handling -------------------------------------- */

static void test_unknown_option_stays_unevaluated(void) {
    /* Unknown option name -> leave unevaluated; printer keeps the head. */
    const char* in = "HermitianMatrixQ[{{1, 2}, {2, 3}}, Foo -> Bar]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "HermitianMatrixQ") != NULL,
               "expected unevaluated HermitianMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_non_rule_extra_arg_stays_unevaluated(void) {
    /* Extra positional (non-Rule) arg -> leave unevaluated. */
    const char* in = "HermitianMatrixQ[{{1, 2}, {2, 3}}, foo]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "HermitianMatrixQ") != NULL,
               "expected unevaluated HermitianMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_zero_arg_stays_unevaluated(void) {
    /* No args -> evaluator should not crash; result keeps the head. */
    const char* in = "HermitianMatrixQ[]";
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "HermitianMatrixQ") != NULL,
               "expected unevaluated HermitianMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Cross-check against ConjugateTranspose ----------------------- */

static void test_matches_conjugate_transpose_equality(void) {
    /* For each matrix, the predicate should agree with the structural
     * test m === ConjugateTranspose[m] up to SameQ. */
    struct {
        const char* matrix;
        const char* expected;
    } cases[] = {
        { "{{1, 3 + 4*I}, {3 - 4*I, 2}}",                "True"  },
        { "{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}",           "True"  },
        { "{{1, 2*I}, {2*I, 3}}",                         "False" },
        { "{{1, 3 + 4*I}, {3 + 4*I, 5}}",                "False" },
        { "{{1, 2}, {3, 4}}",                             "False" },
        { "IdentityMatrix[4]",                            "True"  },
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        char buf[512];
        snprintf(buf, sizeof(buf), "HermitianMatrixQ[%s]", cases[i].matrix);
        assert_eval_eq(buf, cases[i].expected, 0);
    }
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("HermitianMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* HermitianMatrixQ should NOT be Listable -- the input is a matrix,
     * not a list of matrices. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("HermitianMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "Hermitian") != NULL);
    ASSERT(strstr(def->docstring, "SameTest") != NULL);
    ASSERT(strstr(def->docstring, "Tolerance") != NULL);
}

static void test_sym_pointer_interned(void) {
    /* SYM_HermitianMatrixQ must be initialised by sym_names_init() and
     * stably point to the interned name. */
    ASSERT(SYM_HermitianMatrixQ != NULL);
    ASSERT(strcmp(SYM_HermitianMatrixQ, "HermitianMatrixQ") == 0);
}

/* --- Repeated evaluation: leak sanity ----------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Run the same predicate many times; any builtin double-free / dangling
     * pointer typically manifests as a crash or a wrong answer under
     * AddressSanitizer / valgrind on the second iteration onwards. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq(
            "HermitianMatrixQ[{{1, 3 + 4*I, 2 + I},"
            " {3 - 4*I, 0, 4 - 2*I},"
            " {2 - I, 4 + 2*I, -1}}]",
            "True", 0);
        assert_eval_eq(
            "HermitianMatrixQ[{{a, b}, {c, d}}]", "False", 0);
        assert_eval_eq(
            "HermitianMatrixQ[{{1.0, 2.0 + 0.01*I}, {2.0 - 0.02*I, 1.5}},"
            " Tolerance -> 0.1]",
            "True", 0);
        assert_eval_eq(
            "HermitianMatrixQ[{{a, b}, {c, d}}, SameTest -> (True &)]",
            "True", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_real_symmetric_2x2);
    TEST(test_real_symmetric_3x3);
    TEST(test_identity_matrix_is_hermitian);
    TEST(test_complex_hermitian_2x2);
    TEST(test_complex_hermitian_3x3);

    TEST(test_symbolic_hermitian_3x3);
    TEST(test_symbolic_hermitian_2x2);

    TEST(test_complex_symmetric_not_hermitian);
    TEST(test_asymmetric_real_matrix);
    TEST(test_symmetrize_makes_hermitian);
    TEST(test_symbolic_generic_not_hermitian);
    TEST(test_complex_diagonal_not_hermitian);
    TEST(test_symbolic_diagonal_not_real);

    TEST(test_non_matrix_inputs);
    TEST(test_non_square_matrix);
    TEST(test_ragged_matrix);
    TEST(test_three_d_tensor_rejected);
    TEST(test_pure_symbol_stays_evaluated);

    TEST(test_sametest_automatic_is_default);
    TEST(test_sametest_accepts_via_predicate);
    TEST(test_sametest_tolerance_via_predicate);

    TEST(test_tolerance_automatic_is_default);
    TEST(test_tolerance_accepts_within_threshold);
    TEST(test_tolerance_rejects_outside_threshold);

    TEST(test_unknown_option_stays_unevaluated);
    TEST(test_non_rule_extra_arg_stays_unevaluated);
    TEST(test_zero_arg_stays_unevaluated);

    TEST(test_matches_conjugate_transpose_equality);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All HermitianMatrixQ tests passed!\n");
    return 0;
}
