/*
 * test_vectoranal.c -- Grad / Div / Curl / Laplacian (vector analysis).
 *
 * Cartesian (2-arg) forms, fully general in rank, are built on D[f, {{vars}}]:
 *   Grad       scalar -> vector, vector -> Jacobian, tensor -> +1 rank
 *   Div        vector -> scalar, tensor -> contract innermost slot
 *   Curl       2-D vector -> scalar, 3-D vector -> vector, rank-2 -> scalar
 *   Laplacian  scalar, and element-wise over arrays
 * Curvilinear (3-arg) forms use orthogonal scale factors for the charts
 * "Cartesian" / "Polar" / "Cylindrical" / "Spherical"; unsupported field ranks
 * and unknown charts are left unevaluated.
 *
 * Reference outputs are Wolfram Language results; symbolic equality is checked
 * via check_zero (reduce lhs - rhs to 0), the same technique test_deriv_array.c
 * uses, since the canonical printer spelling differs harmlessly from the
 * hand-written form. Exact integer/atomic results use check (FullForm).
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Strict FullForm comparison. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, expected) == 0,
               "Mismatch for %s:\n    expected: %s\n    got:      %s",
               input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* Symbolic-equality check: reduce (expr) via Total[Flatten[Expand[Together]]]
 * to 0. Handles scalar and (nested) list valued differences alike. */
static void check_zero(const char* expr) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "Total[Flatten[Expand[Together[%s]]]]", expr);
    Expr* e = parse_expression(buf);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", buf);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, "0") == 0,
               "Expected zero for %s, got %s", expr, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* ---------------------------------------------------------------------- */
/* Grad -- Cartesian                                                       */
/* ---------------------------------------------------------------------- */
static void test_grad_cartesian(void) {
    check_zero("Grad[Sin[x^2 + y^2], {x, y}] - {2 x Cos[x^2 + y^2], 2 y Cos[x^2 + y^2]}");
    check_zero("Grad[x^2 + y^2 + z^2, {x, y, z}] - {2 x, 2 y, 2 z}");
    check_zero("Grad[x y z, {x, y, z}] - {y z, x z, x y}");
    /* A coordinate not appearing in f contributes a zero entry. */
    check("Grad[x^2, {x, y}]", "List[Times[2, x], 0]");
    /* Gradient of a constant is all zeros. */
    check("Grad[5, {x, y}]", "List[0, 0]");
    /* Grad is exactly the double-list array derivative. */
    check_zero("Grad[x y Sin[z], {x, y, z}] - D[x y Sin[z], {{x, y, z}}]");
}

static void test_grad_jacobian(void) {
    /* Vector field -> Jacobian (new slot appended at the end). */
    check_zero("Grad[{x y, y z, z x}, {x, y, z}]"
               " - {{y, x, 0}, {0, z, y}, {z, 0, x}}");
    /* Grad of a gradient = Hessian (symmetric here). */
    check_zero("Grad[Grad[x y z, {x, y, z}], {x, y, z}]"
               " - {{0, z, y}, {z, 0, x}, {y, x, 0}}");
    check_zero("Grad[x^4 + y^4 - 20 x^2 - 10 x y - 25, {x, y}]"
               " - {-40 x + 4 x^3 - 10 y, -10 x + 4 y^3}");
}

/* ---------------------------------------------------------------------- */
/* Div -- Cartesian                                                        */
/* ---------------------------------------------------------------------- */
static void test_div_cartesian(void) {
    check_zero("Div[{x^2, y^2, z^2}, {x, y, z}] - (2 x + 2 y + 2 z)");
    /* Rank-2 tensor: contract the innermost slot -> a vector. */
    check_zero("Div[{{x y, x y^2, x y^3}, {x^2 y, x^2 y^2, x^2 y^3},"
               " {x^3 y, x^3 y^2, x^3 y^3}}, {x, y, z}]"
               " - {y + 2 x y, 2 x y + 2 x^2 y, 3 x^2 y + 2 x^3 y}");
    /* Div of a gradient is the Laplacian: here Laplacian[x^2+y^2+z^2] = 6. */
    check_zero("Div[Grad[x^2 + y^2 + z^2, {x, y, z}], {x, y, z}] - 6");
}

/* ---------------------------------------------------------------------- */
/* Curl -- Cartesian                                                       */
/* ---------------------------------------------------------------------- */
static void test_curl_cartesian(void) {
    /* 2-D rotational is a scalar. */
    check("Curl[{y, -x}, {x, y}]", "-2");
    /* 3-D curl of a vector field. */
    check_zero("Curl[{x y, y z, z x}, {x, y, z}] - {-y, -z, -x}");
    /* Rank-2 tensor curl (generalized Levi-Civita) collapses to a scalar. */
    check_zero("Curl[{{x y, x y^2, x y^3}, {x^2 y, x^2 y^2, x^2 y^3},"
               " {x^3 y, x^3 y^2, x^3 y^3}}, {x, y, z}]"
               " - 1/2 (x^3 - 3 x y^2 - 3 x^2 y^2 + 2 x y^3)");
    /* Curl of a gradient vanishes identically. */
    check_zero("Curl[Grad[x^2 y z^3, {x, y, z}], {x, y, z}] - {0, 0, 0}");
}

/* ---------------------------------------------------------------------- */
/* Laplacian -- Cartesian                                                  */
/* ---------------------------------------------------------------------- */
static void test_laplacian_cartesian(void) {
    check("Laplacian[x^2 + y^2 + z^2, {x, y, z}]", "6");
    check_zero("Laplacian[Sin[x y], {x, y}] - (-(x^2 + y^2) Sin[x y])");
    /* Applies element-wise over an array (result has f's dimensions). */
    check_zero("Laplacian[{x^2, y^2}, {x, y}] - {2, 2}");
    check_zero("Laplacian[{{x y, x y^2}, {x^2 y, x^2 y^2}}, {x, y}]"
               " - {{0, 2 x}, {2 y, 2 x^2 + 2 y^2}}");
    check_zero("Laplacian[{Sin[x/y], Cos[y/x]}, {x, y}]"
               " - {(2 x y Cos[x/y] - (x^2 + y^2) Sin[x/y])/y^4,"
               " -(((x^2 + y^2) Cos[y/x] + 2 x y Sin[y/x])/x^4)}");
}

/* ---------------------------------------------------------------------- */
/* Charts: Grad / Div / Curl / Laplacian in orthogonal coordinates         */
/* ---------------------------------------------------------------------- */
static void test_grad_charts(void) {
    /* Polar: {(1/1) d/dr, (1/r) d/dtheta}. */
    check_zero("Grad[r^2 Sin[t], {r, t}, \"Polar\"] - {2 r Sin[t], r Cos[t]}");
    /* Spherical force from a potential: -Grad[k q / r] = {k q/r^2, 0, 0}. */
    check_zero("(-Grad[k q / r, {r, t, p}, \"Spherical\"]) - {k q / r^2, 0, 0}");
    /* Cartesian chart agrees with the plain 2-arg form. */
    check_zero("Grad[x^2 + y^2, {x, y}, \"Cartesian\"] - Grad[x^2 + y^2, {x, y}]");
}

static void test_div_charts(void) {
    check_zero("Div[{r Sin[t], -r Cos[t]}, {r, t}, \"Polar\"] - 3 Sin[t]");
    /* Constant field has nonzero divergence in a curvilinear chart. */
    check_zero("Div[{1, 1, 1}, {r, t, p}, \"Spherical\"] - (2/r + Cot[t]/r)");
    check_zero("Div[{r, 0, 0}, {r, t, z}, \"Cylindrical\"] - 2");
}

static void test_curl_charts(void) {
    check_zero("Curl[{r Sin[t], -r Cos[t]}, {r, t}, \"Polar\"] - (-3 Cos[t])");
    /* 3-D cylindrical curl of {0, r, 0} is {0, 0, 2}. */
    check_zero("Curl[{0, r, 0}, {r, t, z}, \"Cylindrical\"] - {0, 0, 2}");
}

static void test_laplacian_charts(void) {
    check_zero("Laplacian[Sin[r^2], {r, t}, \"Polar\"]"
               " - (4 Cos[r^2] - 4 r^2 Sin[r^2])");
    /* Cylindrical Laplacian of r^2 is the constant 4. */
    check_zero("Laplacian[r^2, {r, t, z}, \"Cylindrical\"] - 4");
}

/* ---------------------------------------------------------------------- */
/* Unsupported forms and edge cases: must be left unevaluated (NULL).       */
/* ---------------------------------------------------------------------- */
static void test_unevaluated(void) {
    /* Vector Grad in a chart needs Christoffel symbols -> unevaluated. */
    check("Head[Grad[{1, 1, 1}, {r, t, p}, \"Spherical\"]]", "Grad");
    /* Divergence of a scalar is undefined. */
    check("Head[Div[x^2, {x, y}]]", "Div");
    /* Tensor divergence in a chart is unsupported. */
    check("Head[Div[{{1, 0}, {0, 1}}, {r, t}, \"Polar\"]]", "Div");
    /* Vector Laplacian in a chart is unsupported. */
    check("Head[Laplacian[{1, 1}, {x, y}, \"Polar\"]]", "Laplacian");
    /* Unknown chart -> unevaluated (also emits Grad::chart on stderr). */
    check("Head[Grad[x, {x, y}, \"Elliptic\"]]", "Grad");
    /* Non-list variable specification -> unevaluated. */
    check("Head[Grad[x^2, x]]", "Grad");
    /* 2-D curl needs a length-2 vector; a scalar stays put. */
    check("Head[Curl[x, {x, y}]]", "Curl");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_grad_cartesian);
    TEST(test_grad_jacobian);
    TEST(test_div_cartesian);
    TEST(test_curl_cartesian);
    TEST(test_laplacian_cartesian);
    TEST(test_grad_charts);
    TEST(test_div_charts);
    TEST(test_curl_charts);
    TEST(test_laplacian_charts);
    TEST(test_unevaluated);

    printf("All vector-analysis tests passed.\n");
    symtab_clear();
    return 0;
}
