/* Unit tests for DiagonalMatrixQ.
 *
 * `DiagonalMatrixQ[m]` returns True iff `m` is diagonal -- every entry
 * off the main diagonal is zero.  `DiagonalMatrixQ[m, k]` shifts the
 * test to the k-th diagonal (positive k = superdiagonal, negative k =
 * subdiagonal).  Works on rectangular matrices.  The option
 * `Tolerance -> t` relaxes the zero test to `Abs[e] <= t`.
 *
 * The coverage parallels test_symmetric_matrix_q.c at the shape layer
 * and exercises the additional axes specific to DiagonalMatrixQ:
 *
 *   - 1- and 2-arg surface (default k = 0; explicit k > 0 and k < 0).
 *   - Numeric / real / symbolic / complex matrices.
 *   - Identity / DiagonalMatrix-constructed / Inverse-of-diagonal.
 *   - Rectangular matrices (n x m, n != m).
 *   - Off-diagonal nonzero -> False.
 *   - Symbolic non-numeric-zero off-diagonal -> False (structural).
 *   - Empty-row matrices `{{}, {}}` are vacuously True.
 *   - Tolerance option (Automatic + explicit numeric).
 *   - Non-matrix / scalar / vector / ragged / 3-D / non-square diagonal
 *     test inputs.
 *   - argt diagnostic for zero args.
 *   - nonopt diagnostic for >= 3 positional args / non-Rule junk in the
 *     option region.
 *   - Attribute introspection (Protected; NOT Listable) and docstring
 *     presence.
 *   - sym_names.c interning of `SYM_DiagonalMatrixQ`.
 *   - Repeated evaluation under a tight loop to catch double-frees /
 *     dangling-pointer regressions under valgrind.
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

/* --- Default (main-diagonal) test --------------------------------- */

static void test_diag_3x3_symbolic_diagonal(void) {
    /* Docstring example: {{a,0,0},{0,b,0},{0,0,c}}. */
    assert_eval_eq("DiagonalMatrixQ[{{a, 0, 0}, {0, b, 0}, {0, 0, c}}]",
                   "True", 0);
}

static void test_diag_3x3_nondiagonal(void) {
    /* Docstring example: {{1,0,0},{0,0,2},{3,0,0}} -- 2 and 3 are off. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, 0, 0}, {0, 0, 2}, {3, 0, 0}}]", "False", 0);
}

static void test_diag_zero_matrix_is_diagonal(void) {
    /* The all-zeros matrix is trivially diagonal (every entry off the
     * main diagonal is zero). */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}]", "True", 0);
}

static void test_diag_1x1(void) {
    /* 1x1 matrices have only main-diagonal entries -- always diagonal. */
    assert_eval_eq("DiagonalMatrixQ[{{0}}]", "True", 0);
    assert_eval_eq("DiagonalMatrixQ[{{a}}]", "True", 0);
    assert_eval_eq("DiagonalMatrixQ[{{42}}]", "True", 0);
}

static void test_diag_2x2_pure_symbolic(void) {
    /* {{a,b},{c,d}} is not diagonal for symbolic b, c (they can't be
     * proved zero structurally). */
    assert_eval_eq("DiagonalMatrixQ[{{a, b}, {c, d}}]", "False", 0);
}

static void test_diag_2x2_diag_only(void) {
    /* Mixed symbolic-on-diagonal + zeros off-diagonal. */
    assert_eval_eq("DiagonalMatrixQ[{{a, 0}, {0, b}}]", "True", 0);
}

static void test_diag_machine_real_diagonal(void) {
    /* Real entries on the diagonal, exact zero elsewhere. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1.2, 0, 0}, {0, 2.4, 0}, {0, 0, 3.5}}]",
        "True", 0);
}

static void test_diag_complex_offdiag(void) {
    /* {{1, -3 I}, {0, 4}} -- the -3I off-diagonal entry is not zero. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, -3 I}, {0, 4}}]", "False", 0);
}

static void test_diag_identity_matrix(void) {
    /* IdentityMatrix is diagonal of every size we care about. */
    assert_eval_eq("DiagonalMatrixQ[IdentityMatrix[1]]", "True", 0);
    assert_eval_eq("DiagonalMatrixQ[IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("DiagonalMatrixQ[IdentityMatrix[5]]", "True", 0);
}

static void test_diag_diagonalmatrix_constructed(void) {
    /* DiagonalMatrix constructs a diagonal matrix; DiagonalMatrixQ
     * round-trips it.  Docstring example. */
    assert_eval_eq(
        "DiagonalMatrixQ[DiagonalMatrix[{a, b, c, d}]]", "True", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[DiagonalMatrix[{1, 2, 3, 4, 5}]]", "True", 0);
}

static void test_diag_inverse_of_diagonal(void) {
    /* Inverse[DiagonalMatrix[Range[5]]] should still be diagonal. */
    assert_eval_eq(
        "DiagonalMatrixQ[Inverse[DiagonalMatrix[Range[5]]]]", "True", 0);
}

/* --- Superdiagonals: k > 0 ---------------------------------------- */

static void test_diag_superdiagonal_k1(void) {
    /* Docstring example: {{0,a,0},{0,0,b},{0,0,0}} on k = 1. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, a, 0}, {0, 0, b}, {0, 0, 0}}, 1]",
        "True", 0);
}

static void test_diag_superdiagonal_k1_with_numbers(void) {
    /* {{0,1,0},{0,0,2},{0,0,0}} has nonzero entries only at j-i = 1. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 1, 0}, {0, 0, 2}, {0, 0, 0}}, 1]",
        "True", 0);
    /* The same matrix is NOT k = 2 diagonal. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 1, 0}, {0, 0, 2}, {0, 0, 0}}, 2]",
        "False", 0);
    /* And it is NOT a main-diagonal matrix either. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 1, 0}, {0, 0, 2}, {0, 0, 0}}]",
        "False", 0);
}

static void test_diag_superdiagonal_k2(void) {
    /* k = 2 diagonal of a 3x3: only (0, 2) is on it.  Matrix with
     * only m[0][2] nonzero is k = 2 diagonal. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 5}, {0, 0, 0}, {0, 0, 0}}, 2]",
        "True", 0);
    /* But with an extra nonzero at (1, 2) it fails on k = 2. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 5}, {0, 0, 7}, {0, 0, 0}}, 2]",
        "False", 0);
}

/* --- Subdiagonals: k < 0 ------------------------------------------ */

static void test_diag_subdiagonal_neg1(void) {
    /* Docstring example: {{0,0,0},{a,0,0},{0,b,0}} on k = -1. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {a, 0, 0}, {0, b, 0}}, -1]",
        "True", 0);
}

static void test_diag_subdiagonal_neg2_only(void) {
    /* Docstring example: {{0,0,0},{0,0,0},{3,0,0}} -- 3 is on k = -2,
     * not on k = -1. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {0, 0, 0}, {3, 0, 0}}, -2]",
        "True", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {0, 0, 0}, {3, 0, 0}}, -1]",
        "False", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {0, 0, 0}, {3, 0, 0}}]", "False", 0);
}

/* --- Rectangular matrices ----------------------------------------- */

static void test_diag_rect_2x3(void) {
    /* Docstring example: {{1,0,0},{0,2,0}} -- diagonal 2x3. */
    assert_eval_eq("DiagonalMatrixQ[{{1, 0, 0}, {0, 2, 0}}]", "True", 0);
}

static void test_diag_rect_3x2(void) {
    /* Docstring example: {{1,0},{0,2},{0,0}} -- diagonal 3x2. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, 0}, {0, 2}, {0, 0}}]", "True", 0);
}

static void test_diag_rect_with_offdiag(void) {
    /* 2x3 with off-diagonal nonzero -> False. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, 5, 0}, {0, 2, 0}}]", "False", 0);
    /* 3x2 with off-diagonal nonzero -> False. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, 0}, {3, 2}, {0, 0}}]", "False", 0);
}

static void test_diag_rect_super_k1(void) {
    /* 2x3 superdiagonal at k = 1. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 1, 0}, {0, 0, 2}}, 1]", "True", 0);
}

static void test_diag_rect_sub_neg1(void) {
    /* 3x2 subdiagonal at k = -1. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0}, {1, 0}, {0, 2}}, -1]", "True", 0);
}

static void test_diag_rect_k_out_of_range(void) {
    /* k = 5 on a 2x3 has an empty diagonal -- the whole matrix must
     * be zero.  A nonzero entry => False. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, 0, 0}, {0, 2, 0}}, 5]", "False", 0);
    /* All-zero matrix with absurd k -- True (vacuously diagonal). */
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {0, 0, 0}}, 5]", "True", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[{{0, 0, 0}, {0, 0, 0}}, -7]", "True", 0);
}

/* --- Tolerance option --------------------------------------------- */

static void test_diag_tolerance_relaxes(void) {
    /* Docstring example: a 3x3 with tiny off-diagonal entries.  Without
     * tolerance, those are non-zero -> False.  With Tolerance -> 10^-12,
     * they are absorbed and the predicate accepts. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1., 10^-12, 0}, {0, 2., 10^-13}, {0, 0, 3.}}]",
        "False", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[{{1., 10^-12, 0}, {0, 2., 10^-13}, {0, 0, 3.}},"
        " Tolerance -> 10^-12]",
        "True", 0);
}

static void test_diag_tolerance_automatic_falls_through(void) {
    /* Tolerance -> Automatic must fall through to the structural test,
     * just like SymmetricMatrixQ / HermitianMatrixQ.  A matrix that
     * fails the structural test must still fail with Automatic. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1., 0.5, 0}, {0, 2., 0}, {0, 0, 3.}},"
        " Tolerance -> Automatic]",
        "False", 0);
}

static void test_diag_tolerance_too_tight(void) {
    /* A tolerance smaller than the largest off-diagonal still rejects. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1., 10^-3, 0}, {0, 2., 0}, {0, 0, 3.}},"
        " Tolerance -> 10^-6]",
        "False", 0);
}

static void test_diag_tolerance_with_k(void) {
    /* k = 1 superdiagonal where entries on the +1 diagonal are arbitrary
     * but the rest are within tolerance. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{10^-13, 1.0, 0}, {0, 10^-14, 2.0},"
        " {0, 0, 10^-15}}, 1, Tolerance -> 10^-12]",
        "True", 0);
}

/* --- Non-matrix / non-list / shape edge cases --------------------- */

static void test_diag_scalar_inputs(void) {
    /* Docstring example: DiagonalMatrixQ[1] -> False. */
    assert_eval_eq("DiagonalMatrixQ[1]", "False", 0);
    assert_eval_eq("DiagonalMatrixQ[3.14]", "False", 0);
    assert_eval_eq("DiagonalMatrixQ[\"hello\"]", "False", 0);
    /* Pure symbol resolves to False rather than leaving unevaluated. */
    assert_eval_eq("DiagonalMatrixQ[m]", "False", 0);
}

static void test_diag_vector_inputs(void) {
    /* 1-D Lists are vectors, not matrices. */
    assert_eval_eq("DiagonalMatrixQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("DiagonalMatrixQ[{a, b}]", "False", 0);
}

static void test_diag_empty_list(void) {
    /* {} has dimension {0} -- not a matrix. */
    assert_eval_eq("DiagonalMatrixQ[{}]", "False", 0);
}

static void test_diag_n_by_0_matrix(void) {
    /* Docstring example: {{}, {}} (2x0) is vacuously diagonal. */
    assert_eval_eq("DiagonalMatrixQ[{{}, {}}]", "True", 0);
    assert_eval_eq("DiagonalMatrixQ[{{}}]", "True", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[{{}, {}, {}}]", "True", 0);
}

static void test_diag_ragged_rows(void) {
    /* Ragged row lengths -> not a matrix -> False. */
    assert_eval_eq("DiagonalMatrixQ[{{1, 0}, {0}}]", "False", 0);
    assert_eval_eq("DiagonalMatrixQ[{{1}, {0, 1}}]", "False", 0);
}

static void test_diag_three_d_tensor_rejected(void) {
    /* Higher-rank: entries are themselves Lists -> rejected. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
    assert_eval_eq(
        "DiagonalMatrixQ[{{{1}, {2}}, {{3}, {4}}}]", "False", 0);
}

/* --- Diagonality / structural-zero subtleties --------------------- */

static void test_diag_nonzero_symbolic_offdiag(void) {
    /* Off-diagonal symbol -- can't be proven zero, so reject. */
    assert_eval_eq("DiagonalMatrixQ[{{1, a}, {0, 1}}]", "False", 0);
}

static void test_diag_offdiag_explicit_real_zero(void) {
    /* Real 0.0 should also count as zero. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1.5, 0.0, 0.0}, {0.0, 2.5, 0.0},"
        " {0.0, 0.0, 3.5}}]",
        "True", 0);
}

static void test_diag_offdiag_negative_integer(void) {
    /* A negative integer off-diagonal is non-zero -> False. */
    assert_eval_eq(
        "DiagonalMatrixQ[{{1, -2, 0}, {0, 3, 0}, {0, 0, 4}}]",
        "False", 0);
}

/* --- Argument-count / option diagnostics -------------------------- */

/* Evaluate `in`, assert the printed result still contains the head
 * "DiagonalMatrixQ" (i.e. the call was left unevaluated), and free
 * everything.  Used for argt / nonopt cases. */
static void assert_diagmatrixq_unevaluated(const char* in) {
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "DiagonalMatrixQ") != NULL,
               "expected unevaluated DiagonalMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_zero_arg_stays_unevaluated(void) {
    /* DiagonalMatrixQ[] -- 0 args -> argt diagnostic + unevaluated. */
    assert_diagmatrixq_unevaluated("DiagonalMatrixQ[]");
}

static void test_four_positional_args_nonopt(void) {
    /* Docstring example: DiagonalMatrixQ[1,2,3,4] -> nonopt diagnostic
     * and the call stays unevaluated. */
    assert_diagmatrixq_unevaluated("DiagonalMatrixQ[1, 2, 3, 4]");
}

static void test_three_positional_args_nonopt(void) {
    /* m, k, then junk (not a Rule). */
    assert_diagmatrixq_unevaluated(
        "DiagonalMatrixQ[{{1, 0}, {0, 2}}, 0, foo]");
}

static void test_unknown_option_nonopt(void) {
    /* Unknown option (not Tolerance) -> nonopt. */
    assert_diagmatrixq_unevaluated(
        "DiagonalMatrixQ[{{1, 0}, {0, 2}}, Frobnicate -> 7]");
}

static void test_bad_k_nonopt(void) {
    /* A non-integer, non-Rule second arg is treated as a bad option. */
    assert_diagmatrixq_unevaluated(
        "DiagonalMatrixQ[{{1, 0}, {0, 2}}, 1.5]");
    assert_diagmatrixq_unevaluated(
        "DiagonalMatrixQ[{{1, 0}, {0, 2}}, \"foo\"]");
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("DiagonalMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* DiagonalMatrixQ is NOT Listable -- input is a single matrix. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("DiagonalMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    /* The docstring should mention "diagonal" and the two arities. */
    ASSERT(strstr(def->docstring, "diagonal") != NULL);
    ASSERT(strstr(def->docstring, "DiagonalMatrixQ[m]") != NULL);
    ASSERT(strstr(def->docstring, "DiagonalMatrixQ[m, k]") != NULL);
    ASSERT(strstr(def->docstring, "Tolerance") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_DiagonalMatrixQ != NULL);
    ASSERT(strcmp(SYM_DiagonalMatrixQ, "DiagonalMatrixQ") == 0);
}

/* --- Cross-property: DiagonalMatrix implies DiagonalMatrixQ ------- */

static void test_diag_diagonalmatrix_roundtrip(void) {
    /* Constructing with DiagonalMatrix and querying with DiagonalMatrixQ
     * must always agree -- True for any element list. */
    const char* lists[] = {
        "{1, 2, 3}",
        "{a, b, c, d}",
        "{1.0, 2.0, 3.0, 4.0, 5.0}",
        "{Sin[x], Cos[x], Tan[x]}",
    };
    for (size_t i = 0; i < sizeof(lists)/sizeof(lists[0]); i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "DiagonalMatrixQ[DiagonalMatrix[%s]]", lists[i]);
        assert_eval_eq(buf, "True", 0);
    }
}

/* --- Repeated evaluation: leak sanity ----------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Tight loop over a mix of accepted and rejected shapes / k values
     * to surface any double-free / dangling-pointer regression under
     * valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq(
            "DiagonalMatrixQ[{{a, 0, 0}, {0, b, 0}, {0, 0, c}}]",
            "True", 0);
        assert_eval_eq(
            "DiagonalMatrixQ[{{0, a, 0}, {0, 0, b}, {0, 0, 0}}, 1]",
            "True", 0);
        assert_eval_eq(
            "DiagonalMatrixQ[{{1, 0, 0}, {0, 2, 0}}]", "True", 0);
        assert_eval_eq(
            "DiagonalMatrixQ[{{1., 10^-13, 0}, {0, 2., 0}, {0, 0, 3.}},"
            " Tolerance -> 10^-12]",
            "True", 0);
        assert_eval_eq(
            "DiagonalMatrixQ[{{1, 2}, {3, 4}}]", "False", 0);
        assert_eval_eq("DiagonalMatrixQ[{1, 2, 3}]", "False", 0);
        assert_eval_eq("DiagonalMatrixQ[{}]", "False", 0);
        assert_eval_eq("DiagonalMatrixQ[{{}, {}}]", "True", 0);
        assert_eval_eq("DiagonalMatrixQ[IdentityMatrix[4]]", "True", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_diag_3x3_symbolic_diagonal);
    TEST(test_diag_3x3_nondiagonal);
    TEST(test_diag_zero_matrix_is_diagonal);
    TEST(test_diag_1x1);
    TEST(test_diag_2x2_pure_symbolic);
    TEST(test_diag_2x2_diag_only);
    TEST(test_diag_machine_real_diagonal);
    TEST(test_diag_complex_offdiag);
    TEST(test_diag_identity_matrix);
    TEST(test_diag_diagonalmatrix_constructed);
    TEST(test_diag_inverse_of_diagonal);

    TEST(test_diag_superdiagonal_k1);
    TEST(test_diag_superdiagonal_k1_with_numbers);
    TEST(test_diag_superdiagonal_k2);

    TEST(test_diag_subdiagonal_neg1);
    TEST(test_diag_subdiagonal_neg2_only);

    TEST(test_diag_rect_2x3);
    TEST(test_diag_rect_3x2);
    TEST(test_diag_rect_with_offdiag);
    TEST(test_diag_rect_super_k1);
    TEST(test_diag_rect_sub_neg1);
    TEST(test_diag_rect_k_out_of_range);

    TEST(test_diag_tolerance_relaxes);
    TEST(test_diag_tolerance_automatic_falls_through);
    TEST(test_diag_tolerance_too_tight);
    TEST(test_diag_tolerance_with_k);

    TEST(test_diag_scalar_inputs);
    TEST(test_diag_vector_inputs);
    TEST(test_diag_empty_list);
    TEST(test_diag_n_by_0_matrix);
    TEST(test_diag_ragged_rows);
    TEST(test_diag_three_d_tensor_rejected);

    TEST(test_diag_nonzero_symbolic_offdiag);
    TEST(test_diag_offdiag_explicit_real_zero);
    TEST(test_diag_offdiag_negative_integer);

    TEST(test_zero_arg_stays_unevaluated);
    TEST(test_four_positional_args_nonopt);
    TEST(test_three_positional_args_nonopt);
    TEST(test_unknown_option_nonopt);
    TEST(test_bad_k_nonopt);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_diag_diagonalmatrix_roundtrip);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All DiagonalMatrixQ tests passed!\n");
    return 0;
}
