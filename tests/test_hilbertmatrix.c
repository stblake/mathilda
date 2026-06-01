/* Unit tests for HilbertMatrix.
 *
 *   HilbertMatrix[n]       -- n x n Hilbert matrix, entries 1/(i + j - 1).
 *   HilbertMatrix[{m, n}]  -- m x n Hilbert matrix.
 *
 * Entries are exact Rationals by default (WorkingPrecision -> Infinity);
 * WorkingPrecision -> MachinePrecision gives machine Reals and a digit
 * count above machine precision gives MPFR Reals.
 *
 * Coverage:
 *   - Exact square / rectangular matrices, including the 1x1 and very wide
 *     / very tall edge cases.
 *   - Structural identities: symmetry of the square matrix, and that
 *     Transpose[H[{m,n}]] == H[{n,m}].
 *   - Mathematical correctness against known results: Det[H[2]] = 1/12,
 *     Det[H[3]] = 1/2160, and the exact integer Inverse[H[3]].
 *   - WorkingPrecision -> MachinePrecision and the sub-machine digit-count
 *     degrade-to-machine path.
 *   - WorkingPrecision -> Infinity round-trips to the exact default.
 *   - MPFR path: entries carry the requested Precision and approximate the
 *     true rational value.
 *   - Diagnostics that leave the call unevaluated: ::argx (0 args),
 *     ::dims (non-integer, non-positive, wrong-length / non-integer pair,
 *     symbolic), and ::nonopt (trailing non-option arguments).
 *   - Attributes (Protected), docstring, and the interned SYM_HilbertMatrix.
 *   - Repeated-evaluation loop to surface double-free / leak regressions
 *     under valgrind.
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

/* Helper: evaluate `input` and assert it is returned unevaluated (its
 * printed form is unchanged from a canonical reference string). */
static void assert_unevaluated(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "expected unevaluated %s, got: %s", expected, s);
    free(s);
    expr_free(p);
    expr_free(e);
}

/* --- Exact square matrices ------------------------------------------- */

static void test_square_1x1(void) {
    assert_eval_eq("HilbertMatrix[1]", "{{1}}", 0);
}

static void test_square_2x2(void) {
    assert_eval_eq("HilbertMatrix[2]", "{{1, 1/2}, {1/2, 1/3}}", 0);
}

static void test_square_3x3(void) {
    assert_eval_eq("HilbertMatrix[3]",
                   "{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}", 0);
}

static void test_square_4x4(void) {
    assert_eval_eq(
        "HilbertMatrix[4]",
        "{{1, 1/2, 1/3, 1/4}, {1/2, 1/3, 1/4, 1/5}, "
        "{1/3, 1/4, 1/5, 1/6}, {1/4, 1/5, 1/6, 1/7}}",
        0);
}

/* --- Exact rectangular matrices -------------------------------------- */

static void test_rect_3x5(void) {
    assert_eval_eq(
        "HilbertMatrix[{3, 5}]",
        "{{1, 1/2, 1/3, 1/4, 1/5}, {1/2, 1/3, 1/4, 1/5, 1/6}, "
        "{1/3, 1/4, 1/5, 1/6, 1/7}}",
        0);
}

static void test_rect_2x3(void) {
    assert_eval_eq("HilbertMatrix[{2, 3}]",
                   "{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}}", 0);
}

static void test_rect_tall_3x1(void) {
    /* A single column: more rows than columns. */
    assert_eval_eq("HilbertMatrix[{3, 1}]", "{{1}, {1/2}, {1/3}}", 0);
}

static void test_rect_wide_1x4(void) {
    /* A single row. */
    assert_eval_eq("HilbertMatrix[{1, 4}]", "{{1, 1/2, 1/3, 1/4}}", 0);
}

/* --- Structural identities ------------------------------------------- */

static void test_square_is_symmetric(void) {
    assert_eval_eq("HilbertMatrix[5] === Transpose[HilbertMatrix[5]]",
                   "True", 0);
    assert_eval_eq("SymmetricMatrixQ[HilbertMatrix[4]]", "True", 0);
}

static void test_transpose_swaps_dims(void) {
    /* Transpose of the m x n Hilbert matrix is the n x m Hilbert matrix,
     * since the entry 1/(i+j-1) is symmetric in i and j. */
    assert_eval_eq(
        "Transpose[HilbertMatrix[{3, 5}]] === HilbertMatrix[{5, 3}]",
        "True", 0);
}

static void test_dimensions(void) {
    assert_eval_eq("Dimensions[HilbertMatrix[{3, 5}]]", "{3, 5}", 0);
    assert_eval_eq("Dimensions[HilbertMatrix[4]]", "{4, 4}", 0);
}

/* --- Known mathematical results -------------------------------------- */

static void test_determinants(void) {
    /* Classic small Hilbert determinants. */
    assert_eval_eq("Det[HilbertMatrix[2]]", "1/12", 0);
    assert_eval_eq("Det[HilbertMatrix[3]]", "1/2160", 0);
}

static void test_inverse_exact_integers(void) {
    /* The inverse of a Hilbert matrix has exact integer entries. */
    assert_eval_eq(
        "Inverse[HilbertMatrix[3]]",
        "{{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}",
        0);
}

/* --- WorkingPrecision -> MachinePrecision ---------------------------- */

static void test_machine_precision(void) {
    assert_eval_eq(
        "HilbertMatrix[3, WorkingPrecision -> MachinePrecision]",
        "{{1.0, 0.5, 0.333333}, {0.5, 0.333333, 0.25}, "
        "{0.333333, 0.25, 0.2}}",
        0);
}

static void test_machine_entries_are_machine_numbers(void) {
    assert_eval_eq(
        "MachineNumberQ["
        "  HilbertMatrix[3, WorkingPrecision -> MachinePrecision][[1, 2]]]",
        "True", 0);
}

static void test_submachine_digits_degrade_to_machine(void) {
    /* A digit count at or below machine precision uses the doubles path. */
    assert_eval_eq("HilbertMatrix[2, WorkingPrecision -> 5]",
                   "{{1.0, 0.5}, {0.5, 0.333333}}", 0);
}

/* --- WorkingPrecision -> Infinity (explicit default) ----------------- */

static void test_infinity_is_exact_default(void) {
    assert_eval_eq("HilbertMatrix[3, WorkingPrecision -> Infinity]",
                   "{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}", 0);
}

static void test_exact_entries_have_infinite_precision(void) {
    assert_eval_eq("Precision[HilbertMatrix[3][[1, 2]]]", "Infinity", 0);
}

/* --- WorkingPrecision -> digit count (MPFR) -------------------------- */

static void test_mpfr_precision_carried(void) {
    /* 30-digit entries report Precision ~ 30.103 (= 30 / log10(2) bits
     * back to digits). */
    assert_eval_startswith(
        "Precision[HilbertMatrix[3, WorkingPrecision -> 30][[1, 1]]]",
        "30.1");
    assert_eval_startswith(
        "Precision[HilbertMatrix[3, WorkingPrecision -> 30][[2, 3]]]",
        "30.1");
}

static void test_mpfr_terminating_entries_exact(void) {
    /* Entries with terminating binary expansions are represented exactly. */
    assert_eval_eq("HilbertMatrix[2, WorkingPrecision -> 30][[1, 1]]",
                   "1.0", 0);
    assert_eval_eq("HilbertMatrix[2, WorkingPrecision -> 30][[1, 2]]",
                   "0.5", 0);
}

static void test_mpfr_recurring_entry_value(void) {
    /* 1/3 to 30 digits begins 0.3333333333333333333... */
    assert_eval_startswith(
        "HilbertMatrix[2, WorkingPrecision -> 30][[2, 2]]",
        "0.33333333333333333333");
}

/* --- Diagnostics: ::argx --------------------------------------------- */

static void test_zero_args(void) {
    assert_unevaluated("HilbertMatrix[]", "HilbertMatrix[]");
}

/* --- Diagnostics: ::dims --------------------------------------------- */

static void test_dims_real(void) {
    assert_unevaluated("HilbertMatrix[2.3]", "HilbertMatrix[2.3]");
}

static void test_dims_zero(void) {
    assert_unevaluated("HilbertMatrix[0]", "HilbertMatrix[0]");
}

static void test_dims_negative(void) {
    assert_unevaluated("HilbertMatrix[-2]", "HilbertMatrix[-2]");
}

static void test_dims_triple_list(void) {
    assert_unevaluated("HilbertMatrix[{2, 3, 4}]", "HilbertMatrix[{2, 3, 4}]");
}

static void test_dims_non_positive_pair(void) {
    assert_unevaluated("HilbertMatrix[{0, 3}]", "HilbertMatrix[{0, 3}]");
    assert_unevaluated("HilbertMatrix[{2, -1}]", "HilbertMatrix[{2, -1}]");
}

static void test_dims_symbolic(void) {
    /* A symbolic dimension cannot be resolved -- left unevaluated. */
    assert_unevaluated("HilbertMatrix[x]", "HilbertMatrix[x]");
}

/* --- Diagnostics: ::nonopt ------------------------------------------- */

static void test_nonopt_extra_args(void) {
    assert_unevaluated("HilbertMatrix[2, 34, 4, 5]",
                       "HilbertMatrix[2, 34, 4, 5]");
}

static void test_nonopt_single_bad_trailing(void) {
    assert_unevaluated("HilbertMatrix[3, foo]", "HilbertMatrix[3, foo]");
}

/* --- Attributes / docstring / symbol --------------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("HilbertMatrix");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("HilbertMatrix");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "Hilbert matrix") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_HilbertMatrix != NULL);
    ASSERT(strcmp(SYM_HilbertMatrix, "HilbertMatrix") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix of exact / machine / MPFR / error paths; allocator misuse
     * typically surfaces as a crash, wrong answer, or valgrind diagnostic
     * within a few iterations. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("HilbertMatrix[3]",
                       "{{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}", 0);
        assert_eval_eq("HilbertMatrix[{2, 4}]",
                       "{{1, 1/2, 1/3, 1/4}, {1/2, 1/3, 1/4, 1/5}}", 0);
        assert_eval_eq(
            "HilbertMatrix[2, WorkingPrecision -> MachinePrecision]",
            "{{1.0, 0.5}, {0.5, 0.333333}}", 0);
        assert_eval_startswith(
            "HilbertMatrix[2, WorkingPrecision -> 40][[2, 2]]",
            "0.3333333333333333333");
        assert_unevaluated("HilbertMatrix[-2]", "HilbertMatrix[-2]");
        assert_unevaluated("HilbertMatrix[2, 3, 4]", "HilbertMatrix[2, 3, 4]");
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_square_1x1);
    TEST(test_square_2x2);
    TEST(test_square_3x3);
    TEST(test_square_4x4);

    TEST(test_rect_3x5);
    TEST(test_rect_2x3);
    TEST(test_rect_tall_3x1);
    TEST(test_rect_wide_1x4);

    TEST(test_square_is_symmetric);
    TEST(test_transpose_swaps_dims);
    TEST(test_dimensions);

    TEST(test_determinants);
    TEST(test_inverse_exact_integers);

    TEST(test_machine_precision);
    TEST(test_machine_entries_are_machine_numbers);
    TEST(test_submachine_digits_degrade_to_machine);

    TEST(test_infinity_is_exact_default);
    TEST(test_exact_entries_have_infinite_precision);

    TEST(test_mpfr_precision_carried);
    TEST(test_mpfr_terminating_entries_exact);
    TEST(test_mpfr_recurring_entry_value);

    TEST(test_zero_args);

    TEST(test_dims_real);
    TEST(test_dims_zero);
    TEST(test_dims_negative);
    TEST(test_dims_triple_list);
    TEST(test_dims_non_positive_pair);
    TEST(test_dims_symbolic);

    TEST(test_nonopt_extra_args);
    TEST(test_nonopt_single_bad_trailing);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All HilbertMatrix tests passed!\n");
    return 0;
}
