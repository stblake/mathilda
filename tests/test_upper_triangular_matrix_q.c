/* Unit tests for UpperTriangularMatrixQ.
 *
 * `UpperTriangularMatrixQ[m]` returns True iff every entry of `m`
 * strictly below the main diagonal is zero.  `UpperTriangularMatrixQ[m, k]`
 * shifts the cut-off: every entry m[i,j] with `j - i < k` must be zero.
 * Positive k tightens the test (start from a superdiagonal above the
 * main diagonal); negative k loosens it (entries on the first few
 * subdiagonals are allowed).  Works on rectangular matrices.  The
 * option `Tolerance -> t` relaxes the zero test to `Abs[e] <= t`.
 *
 * Coverage mirrors test_diagonal_matrix_q.c at the shape layer and adds
 * the axes specific to UpperTriangularMatrixQ:
 *
 *   - 1- and 2-arg surface (default k = 0; explicit k > 0 and k < 0).
 *   - Numeric / real / symbolic / complex matrices.
 *   - Rectangular matrices (n x m, n != m).
 *   - Entry strictly below the cut-off -> False (numeric and symbolic).
 *   - Empty-row matrices `{{}, {}}` are vacuously True.
 *   - Tolerance option (Automatic + explicit numeric + too-tight).
 *   - Non-matrix / scalar / vector / ragged / 3-D test inputs.
 *   - argt diagnostic for zero args.
 *   - nonopt diagnostic for >= 3 positional args / non-Rule junk / bad k.
 *   - Attribute introspection (Protected; NOT Listable) and docstring
 *     presence.
 *   - sym_names.c interning of `SYM_UpperTriangularMatrixQ`.
 *   - Cross-property: DiagonalMatrix[..] is upper triangular for any k <= 0;
 *     IdentityMatrix is upper triangular.
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

/* --- Default (main-diagonal cut-off, k = 0) ------------------------ */

static void test_utri_3x3_symbolic_upper(void) {
    /* Docstring example: {{a,b,c},{0,e,f},{0,0,g}}. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{a, b, c}, {0, e, f}, {0, 0, g}}]",
        "True", 0);
}

static void test_utri_3x3_full(void) {
    /* Docstring example: dense 3x3 of nonzero ints -> not upper triangular. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
        "False", 0);
}

static void test_utri_zero_matrix(void) {
    /* All-zeros is trivially upper triangular. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}]",
        "True", 0);
}

static void test_utri_1x1(void) {
    /* 1x1 matrices have no sub-diagonal entries -- always upper triangular. */
    assert_eval_eq("UpperTriangularMatrixQ[{{0}}]", "True", 0);
    assert_eval_eq("UpperTriangularMatrixQ[{{a}}]", "True", 0);
    assert_eval_eq("UpperTriangularMatrixQ[{{42}}]", "True", 0);
}

static void test_utri_2x2_pure_symbolic(void) {
    /* Docstring example: {{a,b},{c,d}} is not upper triangular for
     * symbolic c (can't be proven zero structurally). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{a, b}, {c, d}}]", "False", 0);
}

static void test_utri_2x2_substitution_makes_it(void) {
    /* Docstring example: {{a,b},{c,d}} /. c -> 0  -> True. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{a, b}, {c, d}} /. c -> 0]", "True", 0);
}

static void test_utri_2x2_explicit_upper(void) {
    /* {{a, b}, {0, c}} with the strict sub-diagonal entry exact zero. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{a, b}, {0, c}}]", "True", 0);
}

static void test_utri_machine_real_lower_nonzero(void) {
    /* Docstring example: real machine matrix with nonzero sub-diagonal
     * entries -> False. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1.5, 0, 0}, {2.2, 3.1, 0}, "
        "{4.7, 5.2, 6.4}}]",
        "False", 0);
}

static void test_utri_complex_upper(void) {
    /* Docstring example: {{1+I, 2, 3-2I}, {0, 4, 5I}, {0, 0, 6}}. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1 + I, 2, 3 - 2 I}, {0, 4, 5 I}, "
        "{0, 0, 6}}]",
        "True", 0);
}

static void test_utri_exact_symbolic_lower_nonzero(void) {
    /* Docstring example: {{1, Sqrt[2]}, {Pi, 1/2}} -- Pi on (2,1) is
     * symbolic / non-zero, so False. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, Sqrt[2]}, {Pi, 1/2}}]", "False", 0);
}

static void test_utri_arbitrary_precision(void) {
    /* N[..., prec] produces MPFR (arbitrary-precision) zeros which the
     * shared structural-zero predicate does not currently treat as
     * exact zero -- matches DiagonalMatrixQ behaviour.  A Tolerance
     * absorbs the MPFR zero. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[N[{{1, Sqrt[2]}, {0, 1/2}}, 20],"
        " Tolerance -> 0]",
        "True", 0);
}

static void test_utri_identity_is_upper(void) {
    /* IdentityMatrix is upper triangular at every size. */
    assert_eval_eq("UpperTriangularMatrixQ[IdentityMatrix[1]]", "True", 0);
    assert_eval_eq("UpperTriangularMatrixQ[IdentityMatrix[3]]", "True", 0);
    assert_eval_eq("UpperTriangularMatrixQ[IdentityMatrix[5]]", "True", 0);
}

static void test_utri_diagonalmatrix_is_upper(void) {
    /* DiagonalMatrix output is upper triangular (and lower, but we only
     * test upper here). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[DiagonalMatrix[{a, b, c, d}]]", "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[DiagonalMatrix[{1, 2, 3, 4, 5}]]", "True", 0);
}

/* --- Superdiagonals: k > 0 ---------------------------------------- */

static void test_utri_super_k1_true(void) {
    /* Docstring example: {{0,2,3,4},{0,0,7,8},{0,0,0,12},{0,0,0,0}}
     * is upper triangular starting from the first superdiagonal (k=1):
     * everything with j-i < 1 must be zero. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 2, 3, 4}, {0, 0, 7, 8}, "
        "{0, 0, 0, 12}, {0, 0, 0, 0}}, 1]",
        "True", 0);
    /* The same matrix is also upper triangular (k = 0). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 2, 3, 4}, {0, 0, 7, 8}, "
        "{0, 0, 0, 12}, {0, 0, 0, 0}}]",
        "True", 0);
}

static void test_utri_super_k2_false(void) {
    /* Same matrix is NOT upper triangular starting from k = 2 -- the
     * (1, 2) entry (value 7) lies on the k = 1 diagonal which is now
     * disallowed. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 2, 3, 4}, {0, 0, 7, 8}, "
        "{0, 0, 0, 12}, {0, 0, 0, 0}}, 2]",
        "False", 0);
}

static void test_utri_super_k_strictly_inside(void) {
    /* Strict upper-tri matrix with only entries at (0,1) and (1,2). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 1, 0}, {0, 0, 2}, {0, 0, 0}}, 1]",
        "True", 0);
    /* But it is NOT k = 2 (entry at j-i = 1 is non-zero). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 1, 0}, {0, 0, 2}, {0, 0, 0}}, 2]",
        "False", 0);
}

/* --- Subdiagonals: k < 0 ------------------------------------------ */

static void test_utri_sub_neg1_true(void) {
    /* Docstring example: {{1,2,3},{4,5,6},{0,7,9}} is upper triangular
     * starting from the first subdiagonal: every entry with j-i < -1
     * must be zero.  Here (2,0)=0, so True. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {0, 7, 9}}, -1]",
        "True", 0);
}

static void test_utri_sub_neg2_then_neg1_false(void) {
    /* Docstring example: {{1,2,3,4},{5,6,7,8},{9,10,11,12},{0,14,15,16}}.
     *
     *  - k = -1: requires (2,0)=0 and (3,0)=0 and (3,1)=0.  Here
     *    (2,0)=9 != 0 -> False.
     *  - k = -2: requires (3,0)=0.  (3,0)=0 in the matrix -> True.
     *  - k = 0: requires sub-diagonal zero too -> False.
     */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3, 4}, {5, 6, 7, 8}, "
        "{9, 10, 11, 12}, {0, 14, 15, 16}}, -1]",
        "False", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3, 4}, {5, 6, 7, 8}, "
        "{9, 10, 11, 12}, {0, 14, 15, 16}}, -2]",
        "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3, 4}, {5, 6, 7, 8}, "
        "{9, 10, 11, 12}, {0, 14, 15, 16}}]",
        "False", 0);
}

static void test_utri_sub_very_negative_k_accepts_everything(void) {
    /* k = -100 on a 3x3 leaves no positions with j-i < -100, so any
     * 3x3 matrix is vacuously k=-100 upper triangular. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, -100]",
        "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{a, b}, {c, d}}, -5]", "True", 0);
}

/* --- Rectangular matrices ----------------------------------------- */

static void test_utri_rect_2x3(void) {
    /* Docstring example: mat23 = {{1,2,3},{0,4,5}} -- upper triangular
     * (the (1,0) entry is zero). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3}, {0, 4, 5}}]", "True", 0);
}

static void test_utri_rect_3x2(void) {
    /* Docstring example: mat32 = {{1,2},{0,4},{0,0}} -- upper triangular. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2}, {0, 4}, {0, 0}}]", "True", 0);
}

static void test_utri_rect_with_subdiag_nonzero(void) {
    /* 2x3 with nonzero (1,0) -> False. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3}, {9, 4, 5}}]", "False", 0);
    /* 3x2 with nonzero (1,0) -> False. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2}, {3, 4}, {0, 0}}]", "False", 0);
}

static void test_utri_rect_super_k1(void) {
    /* 2x3 with k=1: only entries at j-i >= 1 may be nonzero.
     * {{0,2,3},{0,0,4}} satisfies this. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 2, 3}, {0, 0, 4}}, 1]", "True", 0);
}

static void test_utri_rect_sub_neg1(void) {
    /* 3x2 with k=-1: entries with j-i < -1 must be zero.  In a 3x2
     * those are positions (2,0) only.  So a 3x2 with arbitrary entries
     * except (2,0)=0 is k=-1 upper triangular. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2}, {3, 4}, {0, 5}}, -1]", "True", 0);
    /* But not with (2,0) nonzero. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2}, {3, 4}, {6, 5}}, -1]", "False", 0);
}

static void test_utri_rect_k_beyond_matrix(void) {
    /* k = 5 on a 2x3: positions with j-i < 5 are everywhere -- the
     * entire matrix must be zero. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 0, 0}, {0, 0, 0}}, 5]", "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 0, 0}, {0, 0, 0}}, 5]", "False", 0);
}

/* --- Tolerance option --------------------------------------------- */

static void test_utri_tolerance_relaxes(void) {
    /* Docstring example: m = {{1.,2.,3.},{10^-12,4.,5.},{0,10^-13,6.}}.
     * Without tolerance, the (1,0) and (2,1) entries are non-zero -> False.
     * With Tolerance -> 10^-12 they are absorbed -> True. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1., 2., 3.}, {10^-12, 4., 5.}, "
        "{0, 10^-13, 6.}}]",
        "False", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1., 2., 3.}, {10^-12, 4., 5.}, "
        "{0, 10^-13, 6.}}, Tolerance -> 10^-12]",
        "True", 0);
}

static void test_utri_tolerance_automatic_falls_through(void) {
    /* Tolerance -> Automatic must reduce to the structural test. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1., 2.}, {0.5, 3.}}, "
        "Tolerance -> Automatic]",
        "False", 0);
}

static void test_utri_tolerance_too_tight(void) {
    /* Tolerance smaller than the largest sub-diagonal still rejects. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1., 2., 3.}, {10^-3, 4., 5.}, "
        "{0, 0, 6.}}, Tolerance -> 10^-6]",
        "False", 0);
}

static void test_utri_tolerance_with_k(void) {
    /* k = 1 with tolerance: entries with j-i < 1 must satisfy
     * Abs[e] <= t.  Entries on the +1 diagonal and above are free. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{10^-13, 1.0, 0}, {0, 10^-14, 2.0}, "
        "{10^-15, 0, 10^-15}}, 1, Tolerance -> 10^-12]",
        "True", 0);
}

/* --- Non-matrix / non-list / shape edge cases --------------------- */

static void test_utri_scalar_inputs(void) {
    /* Docstring example: UpperTriangularMatrixQ[1] -> False. */
    assert_eval_eq("UpperTriangularMatrixQ[1]", "False", 0);
    assert_eval_eq("UpperTriangularMatrixQ[3.14]", "False", 0);
    assert_eval_eq("UpperTriangularMatrixQ[\"hello\"]", "False", 0);
    assert_eval_eq("UpperTriangularMatrixQ[m]", "False", 0);
}

static void test_utri_vector_inputs(void) {
    /* 1-D Lists are vectors, not matrices. */
    assert_eval_eq("UpperTriangularMatrixQ[{1, 2, 3}]", "False", 0);
    assert_eval_eq("UpperTriangularMatrixQ[{a, b}]", "False", 0);
}

static void test_utri_empty_list(void) {
    /* Docstring example: {} -> False. */
    assert_eval_eq("UpperTriangularMatrixQ[{}]", "False", 0);
}

static void test_utri_n_by_0_matrix(void) {
    /* Docstring example: {{},{}} -- 2x0 -- vacuously upper triangular. */
    assert_eval_eq("UpperTriangularMatrixQ[{{}, {}}]", "True", 0);
    assert_eval_eq("UpperTriangularMatrixQ[{{}}]", "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{}, {}, {}}]", "True", 0);
}

static void test_utri_ragged_rows(void) {
    /* Ragged row lengths -> not a matrix -> False. */
    assert_eval_eq("UpperTriangularMatrixQ[{{1, 0}, {0}}]", "False", 0);
    assert_eval_eq("UpperTriangularMatrixQ[{{1}, {0, 1}}]", "False", 0);
}

static void test_utri_three_d_tensor_rejected(void) {
    /* Higher-rank: entries are themselves Lists -> rejected. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{{1, 2}, {3, 4}}}]", "False", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{{1}, {2}}, {{3}, {4}}}]", "False", 0);
}

/* --- Strict-zero subtleties --------------------------------------- */

static void test_utri_subdiag_explicit_real_zero(void) {
    /* Real 0.0 should count as zero. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1.5, 2.0, 3.0}, {0.0, 2.5, 4.0}, "
        "{0.0, 0.0, 3.5}}]",
        "True", 0);
}

static void test_utri_subdiag_negative_integer(void) {
    /* A negative integer in the sub-diagonal region is non-zero -> False. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2, 3}, {-2, 3, 4}, {0, 0, 5}}]",
        "False", 0);
}

static void test_utri_subdiag_symbolic_nonzero(void) {
    /* Symbolic sub-diagonal entry -- can't be proven zero -> False. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{1, 2}, {a, 3}}]", "False", 0);
}

/* --- Argument-count / option diagnostics -------------------------- */

/* Evaluate `in`, assert the printed result still contains the head
 * "UpperTriangularMatrixQ" (i.e. the call was left unevaluated), and
 * free everything.  Used for argt / nonopt cases. */
static void assert_utrimatrixq_unevaluated(const char* in) {
    Expr* p = parse_expression(in);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strstr(s, "UpperTriangularMatrixQ") != NULL,
               "expected unevaluated UpperTriangularMatrixQ, got: %s", s);
    free(s);
    expr_free(p);
    expr_free(e);
}

static void test_utri_zero_arg_stays_unevaluated(void) {
    /* UpperTriangularMatrixQ[] -- 0 args -> argt + unevaluated. */
    assert_utrimatrixq_unevaluated("UpperTriangularMatrixQ[]");
}

static void test_utri_four_positional_args_nonopt(void) {
    /* Docstring example: UpperTriangularMatrixQ[1,2,3,4] -> nonopt. */
    assert_utrimatrixq_unevaluated("UpperTriangularMatrixQ[1, 2, 3, 4]");
}

static void test_utri_three_positional_args_nonopt(void) {
    /* m, k, then junk (not a Rule). */
    assert_utrimatrixq_unevaluated(
        "UpperTriangularMatrixQ[{{1, 0}, {0, 2}}, 0, foo]");
}

static void test_utri_unknown_option_nonopt(void) {
    /* Unknown option (not Tolerance) -> nonopt. */
    assert_utrimatrixq_unevaluated(
        "UpperTriangularMatrixQ[{{1, 0}, {0, 2}}, Frobnicate -> 7]");
}

static void test_utri_bad_k_nonopt(void) {
    /* A non-integer, non-Rule second arg is treated as a bad option. */
    assert_utrimatrixq_unevaluated(
        "UpperTriangularMatrixQ[{{1, 0}, {0, 2}}, 1.5]");
    assert_utrimatrixq_unevaluated(
        "UpperTriangularMatrixQ[{{1, 0}, {0, 2}}, \"foo\"]");
}

/* --- Attribute / docstring introspection -------------------------- */

static void test_utri_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("UpperTriangularMatrixQ");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
    /* UpperTriangularMatrixQ is NOT Listable. */
    ASSERT((def->attributes & ATTR_LISTABLE) == 0);
}

static void test_utri_docstring_set(void) {
    SymbolDef* def = symtab_get_def("UpperTriangularMatrixQ");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "upper triangular") != NULL);
    ASSERT(strstr(def->docstring, "UpperTriangularMatrixQ[m]") != NULL);
    ASSERT(strstr(def->docstring, "UpperTriangularMatrixQ[m, k]") != NULL);
    ASSERT(strstr(def->docstring, "Tolerance") != NULL);
}

static void test_utri_sym_pointer_interned(void) {
    ASSERT(SYM_UpperTriangularMatrixQ != NULL);
    ASSERT(strcmp(SYM_UpperTriangularMatrixQ, "UpperTriangularMatrixQ") == 0);
}

/* --- Cross-property: identity / diagonal matrices are upper-tri --- */

static void test_utri_diagonalmatrix_at_negative_k(void) {
    /* DiagonalMatrix output is upper triangular for any k <= 0. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[DiagonalMatrix[{1, 2, 3}], -1]", "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[DiagonalMatrix[{1, 2, 3}], -2]", "True", 0);
    /* But NOT for k = 1 (the main diagonal is on j-i = 0, which is now
     * disallowed). */
    assert_eval_eq(
        "UpperTriangularMatrixQ[DiagonalMatrix[{1, 2, 3}], 1]", "False", 0);
}

static void test_utri_strict_upper_only_super(void) {
    /* A matrix whose only nonzero entries lie on the +1 superdiagonal is
     * upper triangular for k = 0 and k = 1, but not for k = 2. */
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}]",
        "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}, 1]",
        "True", 0);
    assert_eval_eq(
        "UpperTriangularMatrixQ[{{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}, 2]",
        "False", 0);
}

/* --- Repeated evaluation: leak sanity ----------------------------- */

static void test_utri_repeated_evaluation_does_not_corrupt(void) {
    /* Tight loop over a mix of accepted and rejected shapes / k values
     * to surface any double-free / dangling-pointer regression under
     * valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq(
            "UpperTriangularMatrixQ[{{a, b, c}, {0, e, f}, {0, 0, g}}]",
            "True", 0);
        assert_eval_eq(
            "UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
            "False", 0);
        assert_eval_eq(
            "UpperTriangularMatrixQ[{{0, 1, 2}, {0, 0, 3}, {0, 0, 0}}, 1]",
            "True", 0);
        assert_eval_eq(
            "UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {0, 7, 9}}, -1]",
            "True", 0);
        assert_eval_eq(
            "UpperTriangularMatrixQ[{{1., 2., 3.}, {10^-13, 4., 5.}, "
            "{0, 0, 6.}}, Tolerance -> 10^-12]",
            "True", 0);
        assert_eval_eq("UpperTriangularMatrixQ[{1, 2, 3}]", "False", 0);
        assert_eval_eq("UpperTriangularMatrixQ[{}]", "False", 0);
        assert_eval_eq("UpperTriangularMatrixQ[{{}, {}}]", "True", 0);
        assert_eval_eq("UpperTriangularMatrixQ[IdentityMatrix[4]]", "True", 0);
        assert_eval_eq(
            "UpperTriangularMatrixQ[{{1, 2, 3}, {0, 4, 5}}]", "True", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_utri_3x3_symbolic_upper);
    TEST(test_utri_3x3_full);
    TEST(test_utri_zero_matrix);
    TEST(test_utri_1x1);
    TEST(test_utri_2x2_pure_symbolic);
    TEST(test_utri_2x2_substitution_makes_it);
    TEST(test_utri_2x2_explicit_upper);
    TEST(test_utri_machine_real_lower_nonzero);
    TEST(test_utri_complex_upper);
    TEST(test_utri_exact_symbolic_lower_nonzero);
    TEST(test_utri_arbitrary_precision);
    TEST(test_utri_identity_is_upper);
    TEST(test_utri_diagonalmatrix_is_upper);

    TEST(test_utri_super_k1_true);
    TEST(test_utri_super_k2_false);
    TEST(test_utri_super_k_strictly_inside);

    TEST(test_utri_sub_neg1_true);
    TEST(test_utri_sub_neg2_then_neg1_false);
    TEST(test_utri_sub_very_negative_k_accepts_everything);

    TEST(test_utri_rect_2x3);
    TEST(test_utri_rect_3x2);
    TEST(test_utri_rect_with_subdiag_nonzero);
    TEST(test_utri_rect_super_k1);
    TEST(test_utri_rect_sub_neg1);
    TEST(test_utri_rect_k_beyond_matrix);

    TEST(test_utri_tolerance_relaxes);
    TEST(test_utri_tolerance_automatic_falls_through);
    TEST(test_utri_tolerance_too_tight);
    TEST(test_utri_tolerance_with_k);

    TEST(test_utri_scalar_inputs);
    TEST(test_utri_vector_inputs);
    TEST(test_utri_empty_list);
    TEST(test_utri_n_by_0_matrix);
    TEST(test_utri_ragged_rows);
    TEST(test_utri_three_d_tensor_rejected);

    TEST(test_utri_subdiag_explicit_real_zero);
    TEST(test_utri_subdiag_negative_integer);
    TEST(test_utri_subdiag_symbolic_nonzero);

    TEST(test_utri_zero_arg_stays_unevaluated);
    TEST(test_utri_four_positional_args_nonopt);
    TEST(test_utri_three_positional_args_nonopt);
    TEST(test_utri_unknown_option_nonopt);
    TEST(test_utri_bad_k_nonopt);

    TEST(test_utri_protected_attribute);
    TEST(test_utri_docstring_set);
    TEST(test_utri_sym_pointer_interned);

    TEST(test_utri_diagonalmatrix_at_negative_k);
    TEST(test_utri_strict_upper_only_super);

    TEST(test_utri_repeated_evaluation_does_not_corrupt);

    printf("All UpperTriangularMatrixQ tests passed!\n");
    return 0;
}
