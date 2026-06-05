/* test_fit.c — unit tests for Fit and DesignMatrix.
 *
 * Expected output strings track the printer's canonical form; values match
 * the Wolfram Language reference results from the Fit documentation.
 */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* ---- Fit[data, funs, vars]: the L2 core ---- */

void test_fit_linear(void) {
    /* Line of best fit. */
    assert_eval_eq("Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]",
                   "0.186441 + 0.694915 x", 0);
}

void test_fit_quadratic(void) {
    assert_eval_eq("Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x, x^2}, x]",
                   "0.678392 - 0.266332 x + 0.190955 x^2", 0);
}

void test_fit_vars_as_list(void) {
    /* A bare symbol and a one-element list are equivalent. */
    assert_eval_eq("Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, {x}]",
                   "0.186441 + 0.694915 x", 0);
}

void test_fit_list_of_values(void) {
    /* {v1,...,vn} is paired with coordinates 1..n. */
    assert_eval_eq("Fit[{5,4,3,2,1,0,1,2,3,4,5}, {1, x, x^2}, x]",
                   "7.27273 - 2.0979 x + 0.174825 x^2", 0);
}

void test_fit_quartic_list(void) {
    assert_eval_eq("Fit[{5,4,3,2,1,0,1,2,3,4,5}, {1, x, x^2, x^3, x^4}, x]",
                   "4.54545 + 1.18881 x - 0.938228 x^2 + 0.13986 x^3 "
                   "- 0.00582751 x^4", 0);
}

void test_fit_plane(void) {
    /* Bivariate (plane) fit. */
    assert_eval_eq("Fit[{{0,0,0},{1,0,1},{0,1,2},{1,1,0},{1/2,1/2,1}}, "
                   "{1, x, y}, {x, y}]",
                   "0.8 - 0.5 x + 0.5 y", 0);
}

void test_fit_multivariate_quadratic(void) {
    assert_eval_eq("Fit[{{0,0,0},{1,0,1},{0,1,2},{1,1,0},{1/2,1/2,1}}, "
                   "{1, x, y, x^2, x y, y^2}, {x, y}]",
                   "0.0 + 1.75 x - 0.75 x^2 + 2.25 y - 3.0 x y - 0.25 y^2", 0);
}

void test_fit_transcendental_basis(void) {
    /* Sine basis (machine precision). */
    assert_eval_eq("Fit[N[{{-Pi,4},{-Pi/2,0},{0,1},{Pi/2,-1},{Pi,-4}}], "
                   "{Sin[x/2], Sin[x], Sin[2 x]}, x]",
                   "2.32843 Sin[x] - 4.0 Sin[0.5 x] + 0.0 Sin[2.0 x]", 0);
}

/* ---- WorkingPrecision ---- */

void test_fit_exact_infinity(void) {
    /* Exact rational arithmetic. */
    assert_eval_eq("Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x, x^2}, x, "
                   "WorkingPrecision -> Infinity]",
                   "135/199 - 53/199 x + 38/199 x^2", 0);
}

void test_fit_24_digits(void) {
    /* 24-digit MPFR: the leading digits are stable. */
    assert_eval_startswith(
        "Fit[N[{{-Pi,4},{-Pi/2,0},{0,1},{Pi/2,-1},{Pi,-4}},24], "
        "{Sin[x/2], Sin[x], Sin[2 x]}, x]",
        "2.3284271247461");
}

void test_fit_auto_keeps_input_precision(void) {
    /* WorkingPrecision -> Automatic preserves high-precision input. */
    assert_eval_startswith(
        "Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x, x^2}, x, "
        "WorkingPrecision -> 24]",
        "0.678391959798994974874");
}

/* ---- Fit[{m, v}] design-matrix form ---- */

void test_fit_design_form_exact(void) {
    assert_eval_eq("Fit[{{{1,1},{1,2},{1,3}}, {1,2,3}}, "
                   "WorkingPrecision -> Infinity]",
                   "{0, 1}", 0);
}

void test_fit_design_form_hilbert(void) {
    assert_eval_eq("Fit[{N[HilbertMatrix[4]], Range[4]}]",
                   "{-64.0, 900.0, -2520.0, 1820.0}", 0);
}

/* ---- Regularization ---- */

void test_fit_tikhonov(void) {
    assert_eval_eq("Fit[{{0.,0.},{0.001,1},{0.01,1}}, {1, x, x^2}, x, "
                   "FitRegularization -> {\"Tikhonov\", 1}]",
                   "0.499985 + 0.00549961 x + 4.9996e-05 x^2", 0);
}

void test_fit_ridge_aliases(void) {
    /* "L2" and "RidgeRegression" are aliases of "Tikhonov". */
    assert_eval_eq("Fit[{{0.,0.},{0.001,1},{0.01,1}}, {1, x, x^2}, x, "
                   "FitRegularization -> {\"L2\", 1}]",
                   "0.499985 + 0.00549961 x + 4.9996e-05 x^2", 0);
    assert_eval_eq("Fit[{{0.,0.},{0.001,1},{0.01,1}}, {1, x, x^2}, x, "
                   "FitRegularization -> {\"RidgeRegression\", 1}]",
                   "0.499985 + 0.00549961 x + 4.9996e-05 x^2", 0);
}

void test_fit_lasso_sparsity(void) {
    /* Large penalty drives all non-constant coefficients to exactly zero. */
    assert_eval_eq("Fit[{{0.,0.},{0.001,1},{0.01,1},{.1,0},{1,1}}, "
                   "{1, x, x^2, x^3, x^4}, x, FitRegularization -> {\"LASSO\", 1}]",
                   "0.5 + 0.0 x + 0.0 x^2 + 0.0 x^3 + 0.0 x^4", 0);
    /* "L1" is an alias of "LASSO". */
    assert_eval_eq("Fit[{{0.,0.},{0.001,1},{0.01,1},{.1,0},{1,1}}, "
                   "{1, x, x^2, x^3, x^4}, x, FitRegularization -> {\"L1\", 1}]",
                   "0.5 + 0.0 x + 0.0 x^2 + 0.0 x^3 + 0.0 x^4", 0);
}

void test_fit_lasso_monotone(void) {
    /* A smaller penalty restores the most significant basis function. */
    assert_eval_startswith(
        "Fit[{{0.,0.},{0.001,1},{0.01,1},{.1,0},{1,1}}, "
        "{1, x, x^2, x^3, x^4}, x, FitRegularization -> {\"LASSO\", 0.5}]",
        "0.500006");
}

/* ---- NormFunction ---- */

void test_fit_l1_norm(void) {
    /* Least absolute deviations (IRLS) reproduces the exact LAD line. */
    assert_eval_eq("Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x, "
                   "NormFunction -> Function[Norm[#, 1]]]",
                   "-1.0 + 1.0 x", 0);
}

/* ---- DesignMatrix ---- */

void test_designmatrix_univariate(void) {
    assert_eval_eq("DesignMatrix[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]",
                   "{{1, 0}, {1, 1}, {1, 3}, {1, 5}}", 0);
}

void test_designmatrix_bivariate(void) {
    assert_eval_eq("DesignMatrix[{{0,0,0},{1,0,1},{0,1,2}}, {1, x, y}, {x, y}]",
                   "{{1, 0, 0}, {1, 1, 0}, {1, 0, 1}}", 0);
}

/* ---- Attributes ---- */

void test_fit_protected(void) {
    assert_eval_eq("MemberQ[Attributes[Fit], Protected]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[DesignMatrix], Protected]", "True", 0);
}

/* ---- Error handling: malformed input stays unevaluated ---- */

void test_fit_errors(void) {
    /* Ragged data rows. */
    assert_eval_eq("Fit[{{1,2},{3}}, {1, x}, x]",
                   "Fit[{{1, 2}, {3}}, {1, x}, x]", 0);
    /* Variable count does not match the data coordinates. */
    assert_eval_eq("Fit[{{1,2,3}}, {1, x}, x]",
                   "Fit[{{1, 2, 3}}, {1, x}, x]", 0);
    /* Empty data. */
    assert_eval_eq("Fit[{}, {1, x}, x]", "Fit[{}, {1, x}, x]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_fit_linear);
    TEST(test_fit_quadratic);
    TEST(test_fit_vars_as_list);
    TEST(test_fit_list_of_values);
    TEST(test_fit_quartic_list);
    TEST(test_fit_plane);
    TEST(test_fit_multivariate_quadratic);
    TEST(test_fit_transcendental_basis);
    TEST(test_fit_exact_infinity);
    TEST(test_fit_24_digits);
    TEST(test_fit_auto_keeps_input_precision);
    TEST(test_fit_design_form_exact);
    TEST(test_fit_design_form_hilbert);
    TEST(test_fit_tikhonov);
    TEST(test_fit_ridge_aliases);
    TEST(test_fit_lasso_sparsity);
    TEST(test_fit_lasso_monotone);
    TEST(test_fit_l1_norm);
    TEST(test_designmatrix_univariate);
    TEST(test_designmatrix_bivariate);
    TEST(test_fit_protected);
    TEST(test_fit_errors);

    printf("All fit tests passed!\n");
    return 0;
}
