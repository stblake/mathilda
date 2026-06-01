/* Unit tests for ToeplitzMatrix.
 *
 *   ToeplitzMatrix[n]             -- n x n integer Toeplitz matrix (first row
 *                                    and column the integers 1..n; entry
 *                                    (i, j) is |i - j| + 1, hence symmetric).
 *   ToeplitzMatrix[{c...}]        -- square symmetric Toeplitz matrix, first
 *                                    column (and first row) given.
 *   ToeplitzMatrix[{c...},{r...}] -- m x n Toeplitz matrix, first column and
 *                                    first row given.
 *
 * A Toeplitz matrix is constant along its diagonals; entry (i, j) is
 * c_{i-j+1} when i >= j and r_{j-i+1} otherwise.  Entries are copied verbatim
 * (plain symmetric for the single-list form -- no conjugation), so symbolic,
 * complex, exact and inexact entries flow through unchanged.
 *
 * Coverage:
 *   - Integer square form (1x1 .. 4x4) and its structural shape.
 *   - Symbolic / complex single-list (symmetric) form.
 *   - Two-list rectangular form (wide and tall), plus the corner-sharing
 *     contract c_1 == r_1.
 *   - Structural identities: square forms are symmetric; entries are constant
 *     along diagonals.
 *   - Precision flows from the entries: machine via N, MPFR via `20 marks and
 *     N[..., d]; exact integers stay exact (Precision -> Infinity).
 *   - The ::crs warning still produces the column-element matrix.
 *   - Diagnostics / unevaluated returns: ::argb (0 args), non-list/integer
 *     spec, non-positive integer, empty lists, and over-arity calls.
 *   - Attributes (Protected), docstring, and the interned SYM_ToeplitzMatrix.
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
    assert_eval_eq("ToeplitzMatrix[1]", "{{1}}", 0);
}

static void test_int_2x2(void) {
    assert_eval_eq("ToeplitzMatrix[2]", "{{1, 2}, {2, 1}}", 0);
}

static void test_int_3x3(void) {
    assert_eval_eq("ToeplitzMatrix[3]",
                   "{{1, 2, 3}, {2, 1, 2}, {3, 2, 1}}", 0);
}

static void test_int_4x4(void) {
    assert_eval_eq(
        "ToeplitzMatrix[4]",
        "{{1, 2, 3, 4}, {2, 1, 2, 3}, {3, 2, 1, 2}, {4, 3, 2, 1}}", 0);
}

/* --- Single-list (symmetric) form ------------------------------------ */

static void test_symbolic_square(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{a, b, c, d}]",
        "{{a, b, c, d}, {b, a, b, c}, {c, b, a, b}, {d, c, b, a}}", 0);
}

static void test_symbolic_square_2(void) {
    assert_eval_eq("ToeplitzMatrix[{a, b}]", "{{a, b}, {b, a}}", 0);
}

/* Single-list complex form is plain symmetric (verbatim copy, no
 * conjugation): the upper triangle equals the lower triangle. */
static void test_complex_square(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{1 + 2 I, 3 + 4 I, 5 + 6 I}]",
        "{{1 + 2*I, 3 + 4*I, 5 + 6*I}, {3 + 4*I, 1 + 2*I, 3 + 4*I}, "
        "{5 + 6*I, 3 + 4*I, 1 + 2*I}}", 0);
}

/* --- Two-list (rectangular) form ------------------------------------- */

static void test_rect_tall(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{1, 2, 3, 4, 5}, {1, 6, 7}]",
        "{{1, 6, 7}, {2, 1, 6}, {3, 2, 1}, {4, 3, 2}, {5, 4, 3}}", 0);
}

static void test_rect_wide(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{1, 2, 3}, {1, 4, 5, 6, 7}]",
        "{{1, 4, 5, 6, 7}, {2, 1, 4, 5, 6}, {3, 2, 1, 4, 5}}", 0);
}

static void test_two_list_square(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{c1, c2, c3, c4}, {c1, r2, r3, r4}]",
        "{{c1, r2, r3, r4}, {c2, c1, r2, r3}, {c3, c2, c1, r2}, "
        "{c4, c3, c2, c1}}", 0);
}

static void test_complex_rect(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{1 + 2 I, 3 + 4 I, 5 + 6 I}, {1 + 2 I, 3 + 4 I, 5 + 6 I}]",
        "{{1 + 2*I, 3 + 4*I, 5 + 6*I}, {3 + 4*I, 1 + 2*I, 3 + 4*I}, "
        "{5 + 6*I, 3 + 4*I, 1 + 2*I}}", 0);
}

/* The corner element c_1 must match r_1; when it does not, the column
 * element is used (r_1 is never read by the formula, which reads c_1 on the
 * diagonal) and a ::crs warning is emitted -- but the matrix is still built. */
static void test_corner_mismatch_uses_column(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{c1, c2, c3, c4}, {r1, r2, r3, r4}]",
        "{{c1, r2, r3, r4}, {c2, c1, r2, r3}, {c3, c2, c1, r2}, "
        "{c4, c3, c2, c1}}", 0);
}

/* --- Structural identities ------------------------------------------- */

static void test_square_is_symmetric(void) {
    assert_eval_eq("SymmetricMatrixQ[ToeplitzMatrix[5]]", "True", 0);
    assert_eval_eq("ToeplitzMatrix[5] === Transpose[ToeplitzMatrix[5]]",
                   "True", 0);
}

static void test_single_list_square_is_symmetric(void) {
    assert_eval_eq("SymmetricMatrixQ[ToeplitzMatrix[{a, b, c, d}]]", "True", 0);
}

static void test_constant_along_diagonals(void) {
    /* T[[1,1]] == T[[2,2]] == T[[3,3]] (main diagonal) and
     * T[[1,2]] == T[[2,3]] (first superdiagonal) for any Toeplitz matrix. */
    assert_eval_eq(
        "With[{t = ToeplitzMatrix[{p, q, r, s}, {p, u, v, w}]}, "
        "t[[1, 1]] === t[[2, 2]] && t[[2, 2]] === t[[3, 3]] && "
        "t[[1, 2]] === t[[2, 3]]]",
        "True", 0);
}

static void test_dimensions(void) {
    assert_eval_eq("Dimensions[ToeplitzMatrix[{1, 2, 3}, {1, 4, 5, 6, 7}]]",
                   "{3, 5}", 0);
    assert_eval_eq("Dimensions[ToeplitzMatrix[6]]", "{6, 6}", 0);
}

/* --- Precision flows from the entries -------------------------------- */

static void test_exact_entries_are_exact(void) {
    assert_eval_eq("Precision[ToeplitzMatrix[3][[1, 1]]]", "Infinity", 0);
    assert_eval_eq("Precision[ToeplitzMatrix[3][[3, 1]]]", "Infinity", 0);
}

static void test_machine_via_n(void) {
    assert_eval_eq(
        "N[ToeplitzMatrix[3]]",
        "{{1.0, 2.0, 3.0}, {2.0, 1.0, 2.0}, {3.0, 2.0, 1.0}}", 0);
    assert_eval_eq("MachineNumberQ[N[ToeplitzMatrix[3]][[1, 1]]]", "True", 0);
}

static void test_machine_entries_pass_through(void) {
    assert_eval_eq(
        "ToeplitzMatrix[{1., 2., 4., 8.}, {1., .5, .25, .125}]",
        "{{1.0, 0.5, 0.25, 0.125}, {2.0, 1.0, 0.5, 0.25}, "
        "{4.0, 2.0, 1.0, 0.5}, {8.0, 4.0, 2.0, 1.0}}", 0);
}

static void test_mpfr_precision_via_n(void) {
    /* N[..., 20] gives 20-digit entries. */
    assert_eval_startswith(
        "Precision[N[ToeplitzMatrix[3], 20][[1, 1]]]", "20.");
}

static void test_mpfr_precision_via_backtick(void) {
    assert_eval_startswith(
        "Precision[ToeplitzMatrix[{1`20, 2`20, 4`20}, "
        "{1`20, .5`20, .25`20}][[2, 2]]]", "20.");
}

/* --- Diagnostics / unevaluated returns ------------------------------- */

static void test_zero_args(void) {
    assert_unevaluated("ToeplitzMatrix[]", "ToeplitzMatrix[]");
}

static void test_non_integer_spec(void) {
    /* A bare real or symbol is neither a positive integer nor a list. */
    assert_unevaluated("ToeplitzMatrix[2.3]", "ToeplitzMatrix[2.3]");
    assert_unevaluated("ToeplitzMatrix[x]", "ToeplitzMatrix[x]");
}

static void test_non_positive_integer(void) {
    assert_unevaluated("ToeplitzMatrix[0]", "ToeplitzMatrix[0]");
    assert_unevaluated("ToeplitzMatrix[-3]", "ToeplitzMatrix[-3]");
}

static void test_empty_list(void) {
    assert_unevaluated("ToeplitzMatrix[{}]", "ToeplitzMatrix[{}]");
    assert_unevaluated("ToeplitzMatrix[{1, 2}, {}]",
                       "ToeplitzMatrix[{1, 2}, {}]");
}

static void test_over_arity(void) {
    assert_unevaluated("ToeplitzMatrix[{1, 2}, {1, 3}, {3, 4}]",
                       "ToeplitzMatrix[{1, 2}, {1, 3}, {3, 4}]");
}

/* --- Attributes / docstring / symbol --------------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("ToeplitzMatrix");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("ToeplitzMatrix");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "Toeplitz matrix") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_ToeplitzMatrix != NULL);
    ASSERT(strcmp(SYM_ToeplitzMatrix, "ToeplitzMatrix") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix of integer / symbolic / complex / rectangular / error paths;
     * allocator misuse typically surfaces as a crash, wrong answer, or
     * valgrind diagnostic within a few iterations. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("ToeplitzMatrix[3]",
                       "{{1, 2, 3}, {2, 1, 2}, {3, 2, 1}}", 0);
        assert_eval_eq(
            "ToeplitzMatrix[{a, b, c, d}]",
            "{{a, b, c, d}, {b, a, b, c}, {c, b, a, b}, {d, c, b, a}}", 0);
        assert_eval_eq(
            "ToeplitzMatrix[{1, 2, 3, 4, 5}, {1, 6, 7}]",
            "{{1, 6, 7}, {2, 1, 6}, {3, 2, 1}, {4, 3, 2}, {5, 4, 3}}", 0);
        assert_eval_eq(
            "ToeplitzMatrix[{c1, c2, c3, c4}, {r1, r2, r3, r4}]",
            "{{c1, r2, r3, r4}, {c2, c1, r2, r3}, {c3, c2, c1, r2}, "
            "{c4, c3, c2, c1}}", 0);
        assert_unevaluated("ToeplitzMatrix[]", "ToeplitzMatrix[]");
        assert_unevaluated("ToeplitzMatrix[-3]", "ToeplitzMatrix[-3]");
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

    TEST(test_rect_tall);
    TEST(test_rect_wide);
    TEST(test_two_list_square);
    TEST(test_complex_rect);
    TEST(test_corner_mismatch_uses_column);

    TEST(test_square_is_symmetric);
    TEST(test_single_list_square_is_symmetric);
    TEST(test_constant_along_diagonals);
    TEST(test_dimensions);

    TEST(test_exact_entries_are_exact);
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

    printf("All ToeplitzMatrix tests passed!\n");
    return 0;
}
