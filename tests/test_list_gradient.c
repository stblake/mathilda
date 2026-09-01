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

/* ListGradient — numerical gradient (numpy.gradient port). The expected strings
 * below are cross-checked against numpy.gradient in the plan and the smoke run. */

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "ListGradient %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Default (numpy edge_order = 1) ---------- */

static void test_default_numpy_parity(void) {
    /* numpy.gradient([1,4,9,16,25]) = [3,4,6,8,9]; interior central, edges 1st-order. */
    run_infix("ListGradient[{1, 4, 9, 16, 25}]", "{3, 4, 6, 8, 9}");
}

static void test_default_reals(void) {
    run_infix("ListGradient[{1., 4., 9., 16., 25.}]", "{3.0, 4.0, 6.0, 8.0, 9.0}");
}

static void test_linear_constant_gradient(void) {
    /* Derivative of an arithmetic ramp is the constant step, everywhere. */
    run_infix("ListGradient[{5, 7, 9, 11}]", "{2, 2, 2, 2}");
}

static void test_symbolic_default(void) {
    run_infix("ListGradient[{a, b, c, d, e}]",
              "{-a + b, -1/2 a + 1/2 c, -1/2 b + 1/2 d, -1/2 c + 1/2 e, -d + e}");
}

/* ---------- Spacing ---------- */

static void test_uniform_scalar_spacing(void) {
    run_infix("ListGradient[{1, 4, 9, 16, 25}, 2]", "{3/2, 2, 3, 4, 9/2}");
}

static void test_symbolic_spacing(void) {
    run_infix("ListGradient[{a, b, c}, h]",
              "{-a/h + b/h, -1/2 a/h + 1/2 c/h, -b/h + c/h}");
}

static void test_nonuniform_coordinates(void) {
    /* numpy.gradient([1,4,9,16,25], [0,1,3,6,10]) — non-uniform Fornberg weights. */
    run_infix("ListGradient[{1, 4, 9, 16, 25}, {0, 1, 3, 6, 10}]",
              "{3, 17/6, 73/30, 193/84, 9/4}");
}

/* ---------- Method ---------- */

static void test_forward_exact_on_quadratic(void) {
    /* 2nd-order forward is exact for a quadratic (values k^2 at unit spacing). */
    run_infix("ListGradient[{1., 4., 9., 16., 25.}, Method -> \"Forward\"]",
              "{2.0, 4.0, 6.0, 8.0, 10.0}");
}

static void test_backward_exact_on_quadratic(void) {
    run_infix("ListGradient[{1., 4., 9., 16., 25.}, Method -> \"Backward\"]",
              "{2.0, 4.0, 6.0, 8.0, 10.0}");
}

static void test_bad_method_unevaluated(void) {
    run_infix("Head[ListGradient[{1., 2., 3.}, Method -> \"Bogus\"]]", "ListGradient");
}

/* ---------- Order / WindowLength ---------- */

static void test_order4_exact_on_cubic(void) {
    /* 4th-order scheme is exact on a cubic (values k^3), edges included. */
    run_infix("ListGradient[{0, 1, 8, 27, 64, 125, 216}, DifferenceOrder -> 4]",
              "{0, 3, 12, 27, 48, 75, 108}");
}

static void test_windowlength_equals_order(void) {
    /* WindowLength -> 5 is DifferenceOrder -> 4 restated. */
    run_infix("ListGradient[{0, 1, 8, 27, 64, 125, 216}, WindowLength -> 5]",
              "{0, 3, 12, 27, 48, 75, 108}");
}

/* ---------- Multi-dimensional & Axis ---------- */

static void test_matrix_all_axes(void) {
    run_infix("ListGradient[{{1, 2, 6}, {3, 4, 5}}]",
              "{{{2, 2, -1}, {2, 2, -1}}, {{1, 5/2, 4}, {1, 1, 1}}}");
}

static void test_matrix_axis1_columns(void) {
    run_infix("ListGradient[{{1, 2, 6}, {3, 4, 5}}, Axis -> 1]",
              "{{2, 2, -1}, {2, 2, -1}}");
}

static void test_matrix_axis2_rows(void) {
    run_infix("ListGradient[{{1, 2, 6}, {3, 4, 5}}, Axis -> 2]",
              "{{1, 5/2, 4}, {1, 1, 1}}");
}

static void test_matrix_axis_negative(void) {
    /* Axis -> -1 selects the last axis (== Axis -> 2 here). */
    run_infix("ListGradient[{{1, 2, 6}, {3, 4, 5}}, Axis -> -1]",
              "{{1, 5/2, 4}, {1, 1, 1}}");
}

/* ---------- Packed / NDArray surfaces ---------- */

static void test_visible_ndarray(void) {
    run_infix("Normal[ListGradient[NDArray[{1., 4., 9., 16., 25.}]]]",
              "{3.0, 4.0, 6.0, 8.0, 9.0}");
}

static void test_packed_equals_visible(void) {
    run_infix("ListGradient[{1., 4., 9., 16., 25.}] === "
              "Normal[ListGradient[NDArray[{1., 4., 9., 16., 25.}]]]", "True");
}

static void test_ndarrayq_result(void) {
    run_infix("NDArrayQ[ListGradient[NDArray[{1., 2., 3., 4.}]]]", "True");
}

/* ---------- Number towers: complex, bignum, rational ---------- */

static void test_complex_exact(void) {
    run_infix("ListGradient[{1 + I, 2 - I, 3 + 2 I}]",
              "{1 - 2*I, 1 + 1/2*I, 1 + 3*I}");
}

static void test_complex_packed_buffer(void) {
    /* Complex machine array uses the buffer and agrees with the List path. */
    run_infix("NDArrayQ[ListGradient[NDArray[{1. + I, 2. - I, 3. + 2. I}]]]", "True");
    run_infix("Normal[ListGradient[NDArray[{1. + I, 2. - I, 3. + 2. I}]]] === "
              "ListGradient[{1. + I, 2. - I, 3. + 2. I}]", "True");
}

static void test_bignum_exact(void) {
    run_infix("ListGradient[{10^30, 3 10^30, 6 10^30}]",
              "{2000000000000000000000000000000, "
              "2500000000000000000000000000000, 3000000000000000000000000000000}");
}

static void test_rational_exact(void) {
    run_infix("ListGradient[{1/2, 1/3, 1/5}]", "{-1/6, -3/20, -2/15}");
}

/* ---------- Compile ---------- */

static void test_compiled_call(void) {
    run_infix("Compile[{{v, _Real, 1}}, ListGradient[v]][{1., 4., 9., 16., 25.}]",
              "{3.0, 4.0, 6.0, 8.0, 9.0}");
}

/* ---------- Options / attributes / edge cases ---------- */

static void test_options_defaults(void) {
    run_infix("Options[ListGradient]",
              "{Method -> \"Centered\", DifferenceOrder -> 2, WindowLength -> Automatic, "
              "Axis -> All}");
}

static void test_attributes_protected(void) {
    run_infix("MemberQ[Attributes[ListGradient], Protected]", "True");
}

static void test_nonlist_unevaluated(void) {
    run_infix("ListGradient[x]", "ListGradient[x]");
}

static void test_short_axis_unevaluated(void) {
    /* A length-1 axis has no gradient; stay unevaluated (message on stderr). */
    run_infix("ListGradient[{5}]", "ListGradient[{5}]");
}

int main(void) {
    symtab_init();
    core_init();

    /* Default numpy-parity */
    TEST(test_default_numpy_parity);
    TEST(test_default_reals);
    TEST(test_linear_constant_gradient);
    TEST(test_symbolic_default);

    /* Spacing */
    TEST(test_uniform_scalar_spacing);
    TEST(test_symbolic_spacing);
    TEST(test_nonuniform_coordinates);

    /* Method */
    TEST(test_forward_exact_on_quadratic);
    TEST(test_backward_exact_on_quadratic);
    TEST(test_bad_method_unevaluated);

    /* Order / WindowLength */
    TEST(test_order4_exact_on_cubic);
    TEST(test_windowlength_equals_order);

    /* Multi-dimensional & Axis */
    TEST(test_matrix_all_axes);
    TEST(test_matrix_axis1_columns);
    TEST(test_matrix_axis2_rows);
    TEST(test_matrix_axis_negative);

    /* Packed / NDArray surfaces */
    TEST(test_visible_ndarray);
    TEST(test_packed_equals_visible);
    TEST(test_ndarrayq_result);

    /* Number towers */
    TEST(test_complex_exact);
    TEST(test_complex_packed_buffer);
    TEST(test_bignum_exact);
    TEST(test_rational_exact);

    /* Compile */
    TEST(test_compiled_call);

    /* Options / attributes / edge cases */
    TEST(test_options_defaults);
    TEST(test_attributes_protected);
    TEST(test_nonlist_unevaluated);
    TEST(test_short_axis_unevaluated);

    printf("All ListGradient tests passed!\n");
    return 0;
}
