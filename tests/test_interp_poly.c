/* test_interp_poly.c -- unit tests for InterpolatingPolynomial.
 *
 * Coverage:
 *   - 1-D value form {f1,f2,...} at x = 1,2,... (exact nested Newton-Horner).
 *   - 1-D pairs form {{x1,f1},...} incl. symbolic abscissae and values.
 *   - Hermite / confluent data ({...,{f,df},...}) reproducing derivatives.
 *   - Multidimensional interpolants (value, gradient, Hessian tensors) of
 *     lowest total degree, verified by Expand and the interpolation property.
 *   - Automatic-filled values/derivatives via the minimal-degree solve.
 *   - Failure diagnostics (::poised, ::noipf) and arity (::argrx): the call is
 *     left unevaluated.
 *   - Modulus->n solved exactly over Z/nZ.
 *   - Packed / NDArray data tables (int64 exact, float64 numeric).
 *   - Inexact (Real) data giving Real coefficients.
 *   - No memory leaks (every Expr freed; run under valgrind/leaks).
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Evaluate `input`; assert FullForm equals `expected` (exact / integer /
 * symbolic results). */
static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    } else printf("PASS: %s -> %s\n", input, s);
    free(s); expr_free(res); expr_free(e);
}

/* Evaluate `input`; assert its standard-form string equals `expected`. */
static void run_print(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* s = expr_to_string(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL print: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    } else printf("PASS print: %s -> %s\n", input, s);
    free(s); expr_free(res); expr_free(e);
}

/* Evaluate; assert the result is a machine Real within tol of want. */
static void run_close(const char* input, double want, double tol) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    if (res->type != EXPR_REAL) {
        char* s = expr_to_string_fullform(res);
        printf("FAIL (not a Real): %s -> %s\n", input, s); free(s); ASSERT(0);
    } else if (fabs(res->data.real - want) > tol) {
        printf("FAIL: %s -> %.10g (expected %.10g)\n", input, res->data.real, want);
        ASSERT(0);
    } else printf("PASS: %s -> %.10g\n", input, res->data.real);
    expr_free(res); expr_free(e);
}

/* Evaluate `input`; assert it is left unevaluated with head
 * InterpolatingPolynomial (a failure / arity path). */
static void run_uneval(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    int ok = res->type == EXPR_FUNCTION
        && res->data.function.head->type == EXPR_SYMBOL
        && strcmp(res->data.function.head->data.symbol.name, "InterpolatingPolynomial") == 0;
    if (!ok) {
        char* s = expr_to_string_fullform(res);
        printf("FAIL (expected unevaluated): %s -> %s\n", input, s); free(s); ASSERT(0);
    } else printf("PASS (unevaluated): %s\n", input);
    expr_free(res); expr_free(e);
}

/* --- 1-D value form: exact nested Newton-Horner, exact for exact input. --- */
static void test_value_form(void) {
    run_print("InterpolatingPolynomial[{1,4,9,16},x]", "1 + (-1 + x) (1 + x)");
    run_print("Expand[InterpolatingPolynomial[{1,4,9,16},x]]", "x^2");
    /* Mathematica's exact nested form with rational divided differences. */
    run_print("InterpolatingPolynomial[{4,7,2,8,9},x]",
              "4 + (-1 + x) (3 + (-2 + x) (-4 + (-3 + x) (19/6 - 35/24 (-4 + x))))");
    run_print("InterpolatingPolynomial[{3,7},x]", "3 + 4 (-1 + x)");
    run_full("InterpolatingPolynomial[{5},x]", "5");
    /* Interpolation property at every node. */
    run_full("InterpolatingPolynomial[{4,7,2,8,9},x] /. x->3", "2");
    run_full("InterpolatingPolynomial[{4,7,2,8,9},x] /. x->5", "9");
}

/* --- 1-D pairs form, incl. symbolic abscissae/values. --- */
static void test_pairs_form(void) {
    run_print("InterpolatingPolynomial[{{-1,4},{0,2},{1,6}},x]", "4 + (1 + x) (-2 + 3 x)");
    run_full("InterpolatingPolynomial[{{-1,4},{0,2},{1,6}},x] /. x->0", "2");
    /* symbolic value a,b,c reproduced at the nodes. */
    run_full("InterpolatingPolynomial[{{1,a},{2,b},{3,c}},x] /. x->1", "a");
    run_full("InterpolatingPolynomial[{{1,a},{2,b},{3,c}},x] /. x->2", "b");
    /* symbolic abscissae. */
    run_print("Expand[InterpolatingPolynomial[{{p,1},{q,4}},x]]",
              "1 - 3 p/(-p + q) + 3 x/(-p + q)");
}

/* --- Hermite / confluent divided differences reproduce derivatives. --- */
static void test_hermite(void) {
    run_print("InterpolatingPolynomial[{4,7,2,{8,0},9},x]",
              "4 + (-1 + x) (3 + (-2 + x) (-4 + (-3 + x) "
              "(19/6 + (-4 + x) (-107/36 + 109/72 (-4 + x)))))");
    run_full("InterpolatingPolynomial[{4,7,2,{8,0},9},x] /. x->4", "8");       /* value */
    run_full("D[InterpolatingPolynomial[{4,7,2,{8,0},9},x],x] /. x->4", "0");  /* slope */
    run_full("InterpolatingPolynomial[{4,7,2,{8,0},9},x] /. x->2", "7");
}

/* --- Multidimensional: lowest total degree, verified by Expand + property. --- */
static void test_multivariate(void) {
    const char* d = "{{{0,0},1},{{1,0},7},{{0,1},10},{{2,1},40},{{3,3},151},{{1,2},47}}";
    char buf[512];
    snprintf(buf, sizeof(buf), "Expand[InterpolatingPolynomial[%s,{x,y}]]", d);
    run_print(buf, "1 + 2 x + 4 x^2 + 3 y + 5 x y + 6 y^2");
    snprintf(buf, sizeof(buf), "InterpolatingPolynomial[%s,{x,y}] /. {x->3,y->3}", d);
    run_full(buf, "151");
    snprintf(buf, sizeof(buf), "InterpolatingPolynomial[%s,{x,y}] /. {x->1,y->2}", d);
    run_full(buf, "47");
    /* 3-D value + gradient consistent with a LINEAR polynomial (degree 1). */
    run_print("InterpolatingPolynomial[{{{1,2,3},4,{5,6,7}},{{3,2,1},0}},{x,y,z}]",
              "-34 + 5 x + 6 y + 7 z");
}

/* --- Derivative tensors: partial gradient + full Hessian. --- */
static void test_multivariate_derivs(void) {
    /* Gradient at (0,0), one gradient component Automatic -> xy solves to 0. */
    run_print("Expand[InterpolatingPolynomial["
              "{{{0,0},1,{2,3}},{{1,0},7},{{0,1},10,{Automatic,15}}},{x,y}]]",
              "1 + 2 x + 4 x^2 + 3 y + 6 y^2");
    /* Value + gradient + Hessian at one point recovers x^2 + y^2. */
    run_print("Expand[InterpolatingPolynomial[{{{0,0},0,{0,0},{{2,0},{0,2}}}},{x,y}]]",
              "x^2 + y^2");
}

/* --- Automatic-filled values (minimal-degree solve). --- */
static void test_automatic(void) {
    run_print("Expand[InterpolatingPolynomial[{{-1,Automatic,0},{0,1,1},{1,Automatic,0}},x]]",
              "1 + x - 1/3 x^3");
    /* Symbolic value where Automatic gave noipf. */
    run_print("Expand[InterpolatingPolynomial[{{-1,1},{0,a,0},{1,-1}},x]]",
              "a - a x^2 - x^3");
}

/* --- Failure diagnostics + arity: left unevaluated. --- */
static void test_failures(void) {
    /* Collinear abscissae in 2-D: not poised at total degree 1. */
    run_uneval("InterpolatingPolynomial[{{{0,0},1},{{1,1},2},{{2,2},4}},{x,y}]");
    /* Automatic drops a value: no quadratic through the remaining conditions. */
    run_uneval("InterpolatingPolynomial[{{-1,1},{0,Automatic,0},{1,-1}},x]");
    /* Arity. */
    run_uneval("InterpolatingPolynomial[]");
    run_uneval("InterpolatingPolynomial[{1,2,3}]");
}

/* --- Modulus->n: exact solve over Z/nZ. --- */
static void test_modulus(void) {
    run_print("InterpolatingPolynomial[{1,4,9,16},x,Modulus->7]", "x^2");
    /* Coefficients reduced mod 7; reproduces the data mod 7. */
    run_full("Mod[InterpolatingPolynomial[{2,3,5,1},x,Modulus->7] /. x->1, 7]", "2");
    run_full("Mod[InterpolatingPolynomial[{2,3,5,1},x,Modulus->7] /. x->2, 7]", "3");
    run_full("Mod[InterpolatingPolynomial[{2,3,5,1},x,Modulus->7] /. x->3, 7]", "5");
    run_full("Mod[InterpolatingPolynomial[{2,3,5,1},x,Modulus->7] /. x->4, 7]", "1");
    /* Modulus->0 means no modulus. */
    run_print("InterpolatingPolynomial[{1,4,9,16},x,Modulus->0]", "1 + (-1 + x) (1 + x)");
}

/* --- Packed / NDArray intake. --- */
static void test_packed_ndarray(void) {
    /* int64 data -> exact polynomial, identical to the plain-List result. */
    run_print("InterpolatingPolynomial[NDArray[{1,4,9,16},DataType->\"int64\"],x]",
              "1 + (-1 + x) (1 + x)");
    run_print("Expand[InterpolatingPolynomial[NDArray[{{-1,4},{0,2},{1,6}},"
              "DataType->\"int64\"],x]]", "2 + x + 3 x^2");
    /* float64 data -> numeric, reproduces values. */
    run_close("InterpolatingPolynomial[NDArray[{1.,4.,9.,16.}],x] /. x->3", 9.0, 1e-9);
}

/* --- Inexact (Real) data -> Real coefficients, values reproduced. --- */
static void test_inexact(void) {
    run_close("InterpolatingPolynomial[{1.,4.,9.,16.},x] /. x->2", 4.0, 1e-9);
    run_close("InterpolatingPolynomial[{1.,4.,9.,16.},x] /. x->4", 16.0, 1e-9);
    run_close("InterpolatingPolynomial[{{0.,1.},{1.,3.},{2.,7.}},x] /. x->1", 3.0, 1e-9);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_value_form);
    TEST(test_pairs_form);
    TEST(test_hermite);
    TEST(test_multivariate);
    TEST(test_multivariate_derivs);
    TEST(test_automatic);
    TEST(test_failures);
    TEST(test_modulus);
    TEST(test_packed_ndarray);
    TEST(test_inexact);

    symtab_clear();
    printf("\nAll InterpolatingPolynomial tests passed.\n");
    return 0;
}
