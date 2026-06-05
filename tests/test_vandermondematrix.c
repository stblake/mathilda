/* Unit tests for VandermondeMatrix.
 *
 *   VandermondeMatrix[{x1, ..., xn}]      -- n x n Vandermonde matrix; entry
 *                                            (i, j) is xi^(j-1).
 *   VandermondeMatrix[{x1, ..., xn}, k]   -- n x k Vandermonde matrix.
 *
 * The first column is all ones (xi^0, emitted as the literal 1 so 0^0 reads as
 * 1), the second column is the nodes, the third their squares, and so on.  The
 * nodes need not be numerical or distinct; symbolic nodes stay as Power
 * expressions while numeric powers fold to their value.
 *
 * Coverage:
 *   - Symbolic and numeric square forms and their structural shape.
 *   - Single node; rectangular n x k form (k < n, k = 1, k > n).
 *   - The 0^0 = 1 convention and the all-ones first column.
 *   - Structural / mathematical identities: Dimensions, the determinant equals
 *     the product of node differences, and LinearSolve recovers an
 *     interpolating polynomial's coefficients.
 *   - Precision flows from the nodes: machine via N, MPFR via `20 marks and
 *     N[..., d]; exact integer nodes stay exact (Precision -> Infinity).
 *   - Diagnostics / unevaluated returns: ::argt (0 args), empty list, the
 *     unsupported matrix-conversion form, non-list spec, non-positive k, and
 *     over-arity calls.
 *   - Attributes (Protected), docstring, and the interned SYM_VandermondeMatrix.
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

/* --- Symbolic square form -------------------------------------------- */

static void test_symbolic_square(void) {
    assert_eval_eq(
        "VandermondeMatrix[{x1, x2, x3, x4}]",
        "{{1, x1, x1^2, x1^3}, {1, x2, x2^2, x2^3}, "
        "{1, x3, x3^2, x3^3}, {1, x4, x4^2, x4^3}}", 0);
}

static void test_symbolic_square_2(void) {
    assert_eval_eq("VandermondeMatrix[{a, b}]", "{{1, a}, {1, b}}", 0);
}

/* --- Numeric square form --------------------------------------------- */

static void test_numeric_square(void) {
    assert_eval_eq("VandermondeMatrix[{2, 3, 5}]",
                   "{{1, 2, 4}, {1, 3, 9}, {1, 5, 25}}", 0);
}

static void test_single_node(void) {
    assert_eval_eq("VandermondeMatrix[{7}]", "{{1}}", 0);
}

/* --- Rectangular n x k form ------------------------------------------ */

static void test_rect_narrow(void) {
    assert_eval_eq("VandermondeMatrix[{a, b, c}, 2]",
                   "{{1, a}, {1, b}, {1, c}}", 0);
}

static void test_rect_one_column(void) {
    /* k = 1: just the all-ones column. */
    assert_eval_eq("VandermondeMatrix[{a, b, c}, 1]",
                   "{{1}, {1}, {1}}", 0);
}

static void test_rect_wide(void) {
    /* k > n: extra higher-power columns. */
    assert_eval_eq("VandermondeMatrix[{a}, 5]",
                   "{{1, a, a^2, a^3, a^4}}", 0);
}

static void test_rect_numeric(void) {
    assert_eval_eq("VandermondeMatrix[{2, 3}, 4]",
                   "{{1, 2, 4, 8}, {1, 3, 9, 27}}", 0);
}

/* --- The 0^0 = 1 convention / all-ones first column ------------------ */

static void test_zero_node_first_column_is_one(void) {
    /* A zero node must still give 1 in the first column (0^0 -> 1), not
     * Indeterminate; the rest of its row is the ordinary powers of 0. */
    assert_eval_eq("VandermondeMatrix[{0, 1, 2}]",
                   "{{1, 0, 0}, {1, 1, 1}, {1, 2, 4}}", 0);
}

static void test_first_column_all_ones(void) {
    assert_eval_eq(
        "With[{v = VandermondeMatrix[{p, q, r}]}, "
        "v[[1, 1]] === 1 && v[[2, 1]] === 1 && v[[3, 1]] === 1]",
        "True", 0);
}

static void test_second_column_is_nodes(void) {
    assert_eval_eq(
        "With[{v = VandermondeMatrix[{p, q, r}]}, "
        "{v[[1, 2]], v[[2, 2]], v[[3, 2]]}]",
        "{p, q, r}", 0);
}

/* --- Structural / mathematical identities ---------------------------- */

static void test_dimensions(void) {
    assert_eval_eq("Dimensions[VandermondeMatrix[{a, b, c, d}]]", "{4, 4}", 0);
    assert_eval_eq("Dimensions[VandermondeMatrix[{a, b, c}, 5]]", "{3, 5}", 0);
    assert_eval_eq("Dimensions[VandermondeMatrix[{a, b, c}, 1]]", "{3, 1}", 0);
}

/* The Vandermonde determinant equals the product of node differences:
 * Det = prod_{i<j} (x_j - x_i).  For {1, 2, 3}: (2-1)(3-1)(3-2) = 2. */
static void test_determinant_numeric(void) {
    assert_eval_eq("Det[VandermondeMatrix[{1, 2, 3}]]", "2", 0);
    /* {2, 4, 8}: (4-2)(8-2)(8-4) = 2*6*4 = 48. */
    assert_eval_eq("Det[VandermondeMatrix[{2, 4, 8}]]", "48", 0);
}

static void test_determinant_symbolic(void) {
    /* Det[V] - prod_{i<j}(xj - xi) expands to 0. */
    assert_eval_eq(
        "Expand[Det[VandermondeMatrix[{a, b, c}]] - "
        "(b - a) (c - a) (c - b)]",
        "0", 0);
    /* Factoring the determinant recovers the product of node differences. */
    assert_eval_eq("Factor[Det[VandermondeMatrix[{a, b, c}]]]",
                   "(-a + b) (-a + c) (-b + c)", 0);
}

/* LinearSolve[V, b] recovers the coefficients {a0, a1, ...} of the polynomial
 * p(x) = a0 + a1 x + ... interpolating the points {xi, bi}.  Nodes {1, 2, 3}
 * with values {6, 11, 18} interpolate p(x) = 3 + 2 x + x^2. */
static void test_interpolation_via_linearsolve(void) {
    assert_eval_eq(
        "LinearSolve[VandermondeMatrix[{1, 2, 3}], {6, 11, 18}]",
        "{3, 2, 1}", 0);
}

/* --- Precision flows from the nodes ---------------------------------- */

static void test_exact_nodes_are_exact(void) {
    assert_eval_eq("Precision[VandermondeMatrix[{2, 3, 5}][[3, 3]]]",
                   "Infinity", 0);
}

static void test_machine_via_n(void) {
    assert_eval_eq("N[VandermondeMatrix[{2, 3}]]",
                   "{{1.0, 2.0}, {1.0, 3.0}}", 0);
    assert_eval_eq("MachineNumberQ[N[VandermondeMatrix[{2, 3}]][[2, 2]]]",
                   "True", 0);
}

static void test_machine_nodes_pass_through(void) {
    assert_eval_eq("VandermondeMatrix[{2., 3.}]",
                   "{{1, 2.0}, {1, 3.0}}", 0);
}

static void test_mpfr_precision_via_n(void) {
    assert_eval_startswith(
        "Precision[N[VandermondeMatrix[{2, 3}], 20][[2, 2]]]", "20.");
}

static void test_mpfr_precision_via_backtick(void) {
    assert_eval_startswith(
        "Precision[VandermondeMatrix[{2`20, 3`20}][[2, 2]]]", "20.");
}

/* --- Diagnostics / unevaluated returns ------------------------------- */

static void test_zero_args(void) {
    assert_unevaluated("VandermondeMatrix[]", "VandermondeMatrix[]");
}

static void test_empty_list(void) {
    assert_unevaluated("VandermondeMatrix[{}]", "VandermondeMatrix[{}]");
}

/* A list of lists is the unsupported structured-array conversion form; it is
 * left unevaluated rather than treated as a node list. */
static void test_matrix_conversion_unevaluated(void) {
    assert_unevaluated("VandermondeMatrix[{{1, 2}, {3, 4}}]",
                       "VandermondeMatrix[{{1, 2}, {3, 4}}]");
}

static void test_non_list_spec(void) {
    assert_unevaluated("VandermondeMatrix[x]", "VandermondeMatrix[x]");
    assert_unevaluated("VandermondeMatrix[5]", "VandermondeMatrix[5]");
}

static void test_non_positive_k(void) {
    assert_unevaluated("VandermondeMatrix[{a, b}, 0]",
                       "VandermondeMatrix[{a, b}, 0]");
    assert_unevaluated("VandermondeMatrix[{a, b}, -2]",
                       "VandermondeMatrix[{a, b}, -2]");
}

static void test_non_integer_k(void) {
    assert_unevaluated("VandermondeMatrix[{a, b}, x]",
                       "VandermondeMatrix[{a, b}, x]");
}

static void test_over_arity(void) {
    assert_unevaluated("VandermondeMatrix[{a, b}, 2, 3]",
                       "VandermondeMatrix[{a, b}, 2, 3]");
}

/* --- Attributes / docstring / symbol --------------------------------- */

static void test_protected_attribute(void) {
    SymbolDef* def = symtab_get_def("VandermondeMatrix");
    ASSERT(def != NULL);
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("VandermondeMatrix");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "Vandermonde matrix") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_VandermondeMatrix != NULL);
    ASSERT(strcmp(SYM_VandermondeMatrix, "VandermondeMatrix") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix of symbolic / numeric / rectangular / error paths; allocator
     * misuse typically surfaces as a crash, wrong answer, or valgrind
     * diagnostic within a few iterations. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq(
            "VandermondeMatrix[{x1, x2, x3}]",
            "{{1, x1, x1^2}, {1, x2, x2^2}, {1, x3, x3^2}}", 0);
        assert_eval_eq("VandermondeMatrix[{2, 3, 5}]",
                       "{{1, 2, 4}, {1, 3, 9}, {1, 5, 25}}", 0);
        assert_eval_eq("VandermondeMatrix[{a, b, c}, 2]",
                       "{{1, a}, {1, b}, {1, c}}", 0);
        assert_eval_eq("Det[VandermondeMatrix[{1, 2, 3}]]", "2", 0);
        assert_unevaluated("VandermondeMatrix[]", "VandermondeMatrix[]");
        assert_unevaluated("VandermondeMatrix[{}]", "VandermondeMatrix[{}]");
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_symbolic_square);
    TEST(test_symbolic_square_2);

    TEST(test_numeric_square);
    TEST(test_single_node);

    TEST(test_rect_narrow);
    TEST(test_rect_one_column);
    TEST(test_rect_wide);
    TEST(test_rect_numeric);

    TEST(test_zero_node_first_column_is_one);
    TEST(test_first_column_all_ones);
    TEST(test_second_column_is_nodes);

    TEST(test_dimensions);
    TEST(test_determinant_numeric);
    TEST(test_determinant_symbolic);
    TEST(test_interpolation_via_linearsolve);

    TEST(test_exact_nodes_are_exact);
    TEST(test_machine_via_n);
    TEST(test_machine_nodes_pass_through);
    TEST(test_mpfr_precision_via_n);
    TEST(test_mpfr_precision_via_backtick);

    TEST(test_zero_args);
    TEST(test_empty_list);
    TEST(test_matrix_conversion_unevaluated);
    TEST(test_non_list_spec);
    TEST(test_non_positive_k);
    TEST(test_non_integer_k);
    TEST(test_over_arity);

    TEST(test_protected_attribute);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    printf("All VandermondeMatrix tests passed!\n");
    return 0;
}
