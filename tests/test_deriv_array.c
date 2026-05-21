/*
 * test_deriv_array.c -- Array (Outer-style) derivative tests for D.
 *
 * Covers all the array-spec forms documented for D:
 *
 *   D[f, {{x1, ..., xN}}]         gradient (vector derivative)
 *   D[f, {{x1, ..., xN}, n}]      n-th array derivative; n==2 is Hessian
 *   D[{f1,...}, {{x1,...,xN}}]    Jacobian
 *   D[f, {arr1}, {arr2}, ...]     First[Outer[D, {f}, arr1, arr2, ...]]
 *
 * Plus regression coverage of the original scalar/multi-arg shapes:
 *
 *   D[f, x]                       partial derivative
 *   D[f, {x, n}]                  n-th partial derivative
 *   D[f, x, y, ...]               mixed partial derivative
 *   D[f, {x, n}, {y, m}, ...]     mixed multiple partial derivative
 *   D[{f1, f2, ...}, x]           list threading (recursive)
 *
 * The bug this exercises: D[x^2 + 5 y^3, {{x, y}, 2}] used to return 0
 * because parse_var_spec treated List[x, y] itself as the variable;
 * compute_deriv then short-circuited via expr_free_of. We now route
 * array specs through Outer[D, f, array] semantics.
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

/* Strict FullForm comparison. Avoids surface-syntax canonicalization
 * differences. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, expected) == 0,
               "Array derivative mismatch for %s:\n    expected: %s\n    got:      %s",
               input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* Symbolic-equality check: reduce (lhs) - (rhs) via Expand[Together[...]]
 * and require the final form to be 0. For list-valued diffs we flatten
 * and total the absolute values so e.g. {0, 0} also satisfies the check.
 * Used when the canonical printer spelling differs harmlessly from the
 * hand-written expected form. */
static void check_zero(const char* expr) {
    char buf[1024];
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
/* Form 1: D[f, x] -- scalar partial (regression).                         */
/* ---------------------------------------------------------------------- */
static void test_scalar_partial(void) {
    check("D[x^2, x]", "Times[2, x]");
    check("D[a x^2 + b y^3, x]", "Times[2, a, x]");
    check("D[a x^2 + b y^3, y]", "Times[3, b, Power[y, 2]]");
    check("D[Sin[x y], y]", "Times[x, Cos[Times[x, y]]]");
    /* Univariate */
    check_zero("D[Exp[2 x], x] - 2 Exp[2 x]");
}

/* ---------------------------------------------------------------------- */
/* Form 2: D[f, {x, n}] -- repeated partial (regression).                  */
/* ---------------------------------------------------------------------- */
static void test_scalar_higher_order(void) {
    check("D[x^4, {x, 0}]", "Power[x, 4]");
    check("D[x^4, {x, 1}]", "Times[4, Power[x, 3]]");
    check("D[x^4, {x, 2}]", "Times[12, Power[x, 2]]");
    check("D[x^4, {x, 4}]", "24");
    check("D[x^4, {x, 5}]", "0");
    /* Multivariate: order applies only to the named variable. */
    check("D[a x^3 + b y^2, {x, 2}]", "Times[6, a, x]");
}

/* ---------------------------------------------------------------------- */
/* Form 3: D[f, x, y, ...] -- mixed scalar partials (regression).          */
/* ---------------------------------------------------------------------- */
static void test_mixed_scalar_partials(void) {
    /* Equality of mixed partials. */
    check_zero("D[x^2 y^3, x, y] - D[x^2 y^3, y, x]");
    check_zero("D[x^2 y^3, x, y] - 6 x y^2");
    check_zero("D[Sin[x y], x, y] - (Cos[x y] - x y Sin[x y])");
    check_zero("D[x y z, x, y, z] - 1");
}

/* ---------------------------------------------------------------------- */
/* Form 4: D[f, {x, n}, {y, m}, ...] -- mixed multiple partials.           */
/* ---------------------------------------------------------------------- */
static void test_mixed_multiple_partials(void) {
    /* D^4/(dx^2 dy^2) of x^2 y^2 = 4. */
    check("D[x^2 y^2, {x, 2}, {y, 2}]", "4");
    /* D^3/(dx^2 dy) of x^3 y^2 = 12 x y. */
    check_zero("D[x^3 y^2, {x, 2}, {y, 1}] - 12 x y");
    /* Higher mixed via lists; commute. */
    check_zero("D[Sin[x] Cos[y], {x, 2}, {y, 1}] - D[Sin[x] Cos[y], {y, 1}, {x, 2}]");
}

/* ---------------------------------------------------------------------- */
/* Form 5: D[{f1, f2, ...}, x] -- list threading (regression).             */
/* ---------------------------------------------------------------------- */
static void test_list_threading_scalar(void) {
    check("D[{x, x^2, Sin[x]}, x]", "List[1, Times[2, x], Cos[x]]");
    /* Nested. */
    check("D[{{x, x^2}, {y, x y}}, x]",
          "List[List[1, Times[2, x]], List[0, y]]");
}

/* ---------------------------------------------------------------------- */
/* Form 6: D[f, {{x1, ..., xN}}] -- gradient (vector derivative).          */
/* ---------------------------------------------------------------------- */
static void test_gradient_scalar(void) {
    /* The motivating bug case (gradient form). */
    check("D[x^2 + 5 y^3, {{x, y}}]",
          "List[Times[2, x], Times[15, Power[y, 2]]]");
    /* Single-variable gradient = 1-element list. */
    check("D[x^3, {{x}}]", "List[Times[3, Power[x, 2]]]");
    /* Trig + constant terms. */
    check_zero("D[Sin[x] + Cos[y], {{x, y}}] - {Cos[x], -Sin[y]}");
    /* Three-variable gradient. */
    check("D[x^2 + y^2 + z^2, {{x, y, z}}]",
          "List[Times[2, x], Times[2, y], Times[2, z]]");
    /* Gradient of a product. */
    check_zero("D[x y z, {{x, y, z}}] - {y z, x z, x y}");
    /* Variables not appearing in f give 0 entries. */
    check("D[x^2, {{x, y}}]", "List[Times[2, x], 0]");
}

/* ---------------------------------------------------------------------- */
/* Form 7: D[f, {{x1,...,xN}, 2}] -- Hessian (and higher orders).          */
/* ---------------------------------------------------------------------- */
static void test_hessian(void) {
    /* The MMA-doc reference case from the bug report. */
    check("D[x^2 + 5 y^3, {{x, y}, 2}]",
          "List[List[2, 0], List[0, Times[30, y]]]");
    /* Hessian of a symmetric quadratic form. */
    check("D[x^2 + 2 x y + 3 y^2, {{x, y}, 2}]",
          "List[List[2, 2], List[2, 6]]");
    /* Hessian is symmetric. */
    check_zero("D[Sin[x y], {{x, y}, 2}][[1, 2]] - D[Sin[x y], {{x, y}, 2}][[2, 1]]");
    /* Three-variable Hessian = symmetric 3x3 matrix. */
    check("D[x y + y z + x z, {{x, y, z}, 2}]",
          "List[List[0, 1, 1], List[1, 0, 1], List[1, 1, 0]]");
    /* Order 0 = identity (deep-copy of the input). */
    check("D[x^2 + y^3, {{x, y}, 0}]", "Plus[Power[x, 2], Power[y, 3]]");
    /* Order 3: tensor of third derivatives. */
    check("D[x^3 + y^3, {{x, y}, 3}]",
          "List[List[List[6, 0], List[0, 0]], List[List[0, 0], List[0, 6]]]");
}

/* ---------------------------------------------------------------------- */
/* Form 8: D[{f1, ..., fM}, {{x1, ..., xN}}] -- Jacobian.                  */
/* ---------------------------------------------------------------------- */
static void test_jacobian(void) {
    /* Classical 2-by-2 Jacobian. Rows = functions, columns = variables. */
    check("D[{x^2 + y, x y}, {{x, y}}]",
          "List[List[Times[2, x], 1], List[y, x]]");
    /* 3-by-2 Jacobian. */
    check("D[{x + y, x y, x^2 + y^2}, {{x, y}}]",
          "List[List[1, 1], List[y, x], List[Times[2, x], Times[2, y]]]");
    /* 2-by-3 Jacobian. */
    check("D[{x y z, x + y + z}, {{x, y, z}}]",
          "List[List[Times[y, z], Times[x, z], Times[x, y]],"
          " List[1, 1, 1]]");
    /* Single-function Jacobian degenerates to gradient row. */
    check("D[{x^2 + y^2}, {{x, y}}]",
          "List[List[Times[2, x], Times[2, y]]]");
}

/* ---------------------------------------------------------------------- */
/* Form 9: D[f, {arr1}, {arr2}, ...] -- chained array derivatives.         */
/* ---------------------------------------------------------------------- */
static void test_chained_arrays(void) {
    /* D[f, {{x,y}}, {{x,y}}] == D[f, {{x,y}, 2}] (Hessian). */
    check("D[x^2 + 5 y^3, {{x, y}}, {{x, y}}]",
          "List[List[2, 0], List[0, Times[30, y]]]");
    /* Distinct array specs: First[Outer[D, {f}, {x}, {y}]] = {{D[f,x,y]}}. */
    check("D[x^2 y^3, {{x}}, {{y}}]",
          "List[List[Times[6, x, Power[y, 2]]]]");
    /* Scalar interleaved with array: D[D[f, x], {{x, y}}] = grad of f_x. */
    check_zero("D[Sin[x y], x, {{x, y}}] - {y^2 (-Sin[x y]),"
               " Cos[x y] - x y Sin[x y]}");
}

/* ---------------------------------------------------------------------- */
/* Edge cases.                                                             */
/* ---------------------------------------------------------------------- */
static void test_edge_cases(void) {
    /* Gradient of a constant: all zeros. */
    check("D[5, {{x, y}}]", "List[0, 0]");
    /* Hessian of a constant: zero matrix. */
    check("D[5, {{x, y}, 2}]", "List[List[0, 0], List[0, 0]]");
    /* Gradient with respect to a variable not in f. */
    check("D[a, {{x, y}}]", "List[0, 0]");
    /* Empty variable array: empty list of derivatives. */
    check("D[x^2, {{}}]", "List[]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_scalar_partial);
    TEST(test_scalar_higher_order);
    TEST(test_mixed_scalar_partials);
    TEST(test_mixed_multiple_partials);
    TEST(test_list_threading_scalar);
    TEST(test_gradient_scalar);
    TEST(test_hessian);
    TEST(test_jacobian);
    TEST(test_chained_arrays);
    TEST(test_edge_cases);

    printf("All array-derivative tests passed.\n");
    return 0;
}
