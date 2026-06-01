/* Unit tests for HankelMatrix.
 *
 *   HankelMatrix[n]            -- n x n integer Hankel matrix (first row and
 *                                 column the integers 1..n, zeros below the
 *                                 antidiagonal).
 *   HankelMatrix[{c...}]       -- square Hankel matrix, first column given.
 *   HankelMatrix[{c...},{r...}] -- m x n Hankel matrix, first column and last
 *                                 row given.
 *
 * A Hankel matrix is constant along its antidiagonals; entry (i, j) is
 * c_{i+j-1} when i+j-1 <= m and r_{i+j-m} otherwise.  Entries are copied
 * verbatim, so symbolic, complex, exact and inexact entries flow through.
 *
 * Coverage:
 *   - Integer square form (1x1 .. 4x4) and its structural shape.
 *   - Symbolic / complex single-list (square) form.
 *   - Two-list rectangular form (wide and tall), plus the corner-sharing
 *     contract c_m == r_1.
 *   - Structural identities: square HankelMatrix[n] is symmetric; entries
 *     are constant along antidiagonals.
 *   - Precision flows from the entries: machine via N, MPFR via `20 marks,
 *     and N[..., d]; exact integers stay exact (Precision -> Infinity), and
 *     the zero fill stays an exact 0.
 *   - The ::crs warning still produces the column-element matrix.
 *   - Diagnostics / unevaluated returns: ::argb (0 args), non-list/integer
 *     spec, empty lists, and over-arity calls.
 *   - Attributes (Protected), docstring, and the interned SYM_HankelMatrix.
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
 * printed form equals a canonical reference string). */
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

/* --- Integer square form --------------------------------------------- */

static void test_int_1x1(void) {
    assert_eval_eq("HankelMatrix[1]", "{{1}}", 0);
}

static void test_int_2x2(void) {
    assert_eval_eq("HankelMatrix[2]", "{{1, 2}, {2, 0}}", 0);
}

static void test_int_3x3(void) {
    assert_eval_eq("HankelMatrix[3]",
                   "{{1, 2, 3}, {2, 3, 0}, {3, 0, 0}}", 0);
}

static void test_int_4x4(void) {
    assert_eval_eq(
        "HankelMatrix[4]",
        "{{1, 2, 3, 4}, {2, 3, 4, 0}, {3, 4, 0, 0}, {4, 0, 0, 0}}", 0);
}

/* --- Single-list (square) form --------------------------------------- */

static void test_symbolic_square(void) {
    assert_eval_eq(
        "HankelMatrix[{a, b, c, d}]",
        "{{a, b, c, d}, {b, c, d, 0}, {c, d, 0, 0}, {d, 0, 0, 0}}", 0);
}

static void test_symbolic_square_2(void) {
    assert_eval_eq("HankelMatrix[{a, b}]", "{{a, b}, {b, 0}}", 0);
}

static void test_complex_square(void) {
    assert_eval_eq(
        "HankelMatrix[{1, 1 + 2 I, 3 + 4 I}]",
        "{{1, 1 + 2*I, 3 + 4*I}, {1 + 2*I, 3 + 4*I, 0}, {3 + 4*I, 0, 0}}", 0);
}

/* --- Two-list (rectangular) form ------------------------------------- */

static void test_rect_wide(void) {
    assert_eval_eq(
        "HankelMatrix[{x, y, z}, {z, a, b, c, d}]",
        "{{x, y, z, a, b}, {y, z, a, b, c}, {z, a, b, c, d}}", 0);
}

static void test_rect_wide_integers(void) {
    assert_eval_eq(
        "HankelMatrix[{1, 2, 3, 4}, {4, 5, 6, 0, 0, 0}]",
        "{{1, 2, 3, 4, 5, 6}, {2, 3, 4, 5, 6, 0}, "
        "{3, 4, 5, 6, 0, 0}, {4, 5, 6, 0, 0, 0}}", 0);
}

static void test_rect_tall_integers(void) {
    assert_eval_eq(
        "HankelMatrix[{1, 2, 3, 4, 5, 6}, {6, 0, 0, 0}]",
        "{{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}, "
        "{4, 5, 6, 0}, {5, 6, 0, 0}, {6, 0, 0, 0}}", 0);
}

static void test_complex_rect(void) {
    assert_eval_eq(
        "HankelMatrix[{1, 1 + 2 I, 3 + 4 I}, {3 + 4 I, 5 + 6 I, 7 + 8 I}]",
        "{{1, 1 + 2*I, 3 + 4*I}, {1 + 2*I, 3 + 4*I, 5 + 6*I}, "
        "{3 + 4*I, 5 + 6*I, 7 + 8*I}}", 0);
}

/* The corner element c_m must match r_1; when it does not, the column
 * element is used (r_1 is never read by the formula) and a ::crs warning is
 * emitted -- but the matrix is still built. */
static void test_corner_mismatch_uses_column(void) {
    assert_eval_eq(
        "HankelMatrix[{1, 2, 3, 4}, {a, b, c, d}]",
        "{{1, 2, 3, 4}, {2, 3, 4, b}, {3, 4, b, c}, {4, b, c, d}}", 0);
}

/* --- Structural identities ------------------------------------------- */

static void test_square_is_symmetric(void) {
    assert_eval_eq("SymmetricMatrixQ[HankelMatrix[5]]", "True", 0);
    assert_eval_eq("HankelMatrix[5] === Transpose[HankelMatrix[5]]", "True", 0);
}

static void test_single_list_square_is_symmetric(void) {
    assert_eval_eq("SymmetricMatrixQ[HankelMatrix[{a, b, c, d}]]", "True", 0);
}

static void test_constant_along_antidiagonals(void) {
    /* H[[1,3]] == H[[2,2]] == H[[3,1]] for any Hankel matrix. */
    assert_eval_eq(
        "With[{h = HankelMatrix[{p, q, r, s}, {s, t, u, v}]}, "
        "h[[1, 3]] === h[[2, 2]] && h[[2, 2]] === h[[3, 1]]]",
        "True", 0);
}

static void test_dimensions(void) {
    assert_eval_eq("Dimensions[HankelMatrix[{1, 2, 3}, {3, 4, 5, 6, 7}]]",
                   "{3, 5}", 0);
    assert_eval_eq("Dimensions[HankelMatrix[6]]", "{6, 6}", 0);
}

/* --- Precision flows from the entries -------------------------------- */

static void test_exact_entries_are_exact(void) {
    assert_eval_eq("Precision[HankelMatrix[3][[1, 1]]]", "Infinity", 0);
}

static void test_zero_fill_is_exact_integer(void) {
    /* The antidiagonal fill is an exact integer 0, not a machine 0.0. */
    assert_eval_eq("Precision[HankelMatrix[3][[3, 3]]]", "Infinity", 0);
    assert_eval_eq("HankelMatrix[3][[3, 3]]", "0", 0);
}

static void test_machine_via_n(void) {
    assert_eval_eq(
        "N[HankelMatrix[3]]",
        "{{1.0, 2.0, 3.0}, {2.0, 3.0, 0.0}, {3.0, 0.0, 0.0}}", 0);
    assert_eval_eq("MachineNumberQ[N[HankelMatrix[3]][[1, 1]]]", "True", 0);
}

static void test_machine_entries_pass_through(void) {
    assert_eval_eq(
        "HankelMatrix[{1., 2., 3.}, {3., 2., 1.}]",
        "{{1.0, 2.0, 3.0}, {2.0, 3.0, 2.0}, {3.0, 2.0, 1.0}}", 0);
}

static void test_mpfr_precision_via_n(void) {
    /* N[..., 20] gives 20-digit entries. */
    assert_eval_startswith(
        "Precision[N[HankelMatrix[3], 20][[1, 1]]]", "20.");
}

static void test_mpfr_precision_via_backtick(void) {
    assert_eval_startswith(
        "Precision[HankelMatrix[{1`20, 2`20, 3`20}, "
        "{3`20, 2`20, 1`20}][[2, 2]]]", "20.");
}

/* --- Diagnostics / unevaluated returns ------------------------------- */

static void test_zero_args(void) {
    assert_unevaluated("HankelMatrix[]", "HankelMatrix[]");
}

static void test_non_integer_spec(void) {
    /* A bare real or symbol is neither a positive integer nor a list. */
    assert_unevaluated("HankelMatrix[2.3]", "HankelMatrix[2.3]");
    assert_unevaluated("HankelMatrix[x]", "HankelMatrix[x]");
}

static void test_non_positive_integer(void) {
    assert_unevaluated("HankelMatrix[0]", "HankelMatrix[0]");
    assert_unevaluated("HankelMatrix[-3]", "HankelMatrix[-3]");
}

static void test_empty_list(void) {
    assert_unevaluated("HankelMatrix[{}]", "HankelMatrix[{}]");
    assert_unevaluated("HankelMatrix[{1, 2}, {}]", "HankelMatrix[{1, 2}, {}]");
}

static void test_over_arity(void) {
    assert_unevaluated("HankelMatrix[{1, 2}, {2, 3}, {3, 4}]",
                       "HankelMatrix[{1, 2}, {2, 3}, {3, 4}]");
}

/* --- Attributes / docstring / symbol --------------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("HankelMatrix");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("HankelMatrix");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "Hankel matrix") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_HankelMatrix != NULL);
    ASSERT(strcmp(SYM_HankelMatrix, "HankelMatrix") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix of integer / symbolic / complex / rectangular / error paths;
     * allocator misuse typically surfaces as a crash, wrong answer, or
     * valgrind diagnostic within a few iterations. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("HankelMatrix[3]",
                       "{{1, 2, 3}, {2, 3, 0}, {3, 0, 0}}", 0);
        assert_eval_eq(
            "HankelMatrix[{a, b, c, d}]",
            "{{a, b, c, d}, {b, c, d, 0}, {c, d, 0, 0}, {d, 0, 0, 0}}", 0);
        assert_eval_eq(
            "HankelMatrix[{x, y, z}, {z, a, b, c, d}]",
            "{{x, y, z, a, b}, {y, z, a, b, c}, {z, a, b, c, d}}", 0);
        assert_eval_eq(
            "HankelMatrix[{1, 2, 3, 4}, {a, b, c, d}]",
            "{{1, 2, 3, 4}, {2, 3, 4, b}, {3, 4, b, c}, {4, b, c, d}}", 0);
        assert_unevaluated("HankelMatrix[]", "HankelMatrix[]");
        assert_unevaluated("HankelMatrix[-3]", "HankelMatrix[-3]");
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_int_1x1);
    TEST(test_int_2x2);
    TEST(test_int_3x3);
    TEST(test_int_4x4);

    TEST(test_symbolic_square);
    TEST(test_symbolic_square_2);
    TEST(test_complex_square);

    TEST(test_rect_wide);
    TEST(test_rect_wide_integers);
    TEST(test_rect_tall_integers);
    TEST(test_complex_rect);
    TEST(test_corner_mismatch_uses_column);

    TEST(test_square_is_symmetric);
    TEST(test_single_list_square_is_symmetric);
    TEST(test_constant_along_antidiagonals);
    TEST(test_dimensions);

    TEST(test_exact_entries_are_exact);
    TEST(test_zero_fill_is_exact_integer);
    TEST(test_machine_via_n);
    TEST(test_machine_entries_pass_through);
    TEST(test_mpfr_precision_via_n);
    TEST(test_mpfr_precision_via_backtick);

    TEST(test_zero_args);
    TEST(test_non_integer_spec);
    TEST(test_non_positive_integer);
    TEST(test_empty_list);
    TEST(test_over_arity);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All HankelMatrix tests passed!\n");
    return 0;
}
