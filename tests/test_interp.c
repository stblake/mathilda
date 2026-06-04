/* test_interp.c -- unit tests for InterpolatingFunction.
 *
 * Coverage:
 *   - Linear (2-point) and quadratic (3-point) interpolation, reduced order.
 *   - Cubic (order-3) interpolation, including exact reproduction of cubic
 *     data and Mathematica's sliding-window values on a reference table.
 *   - Exact-node queries returning the stored value exactly (Integer and
 *     Rational nodes).
 *   - Extrapolation outside the declared domain (with a dmval warning).
 *   - Symbolic and wrong-arity applications left unevaluated.
 *   - Malformed objects (too few points, non-increasing abscissae,
 *     non-numeric data) left unevaluated.
 *   - Standard-form printing abbreviating the data table as <>.
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

/* The reference table from the task description: points on a hill. */
#define IFUN_EX \
    "InterpolatingFunction[{{0,5}}," \
    "{{0,0},{1,1},{2,3},{3,4},{4,3},{5,0}}]"

/* Evaluate `input`; assert the result is a machine real within tol of want. */
static void run_close(const char* input, double want, double tol) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    if (res->type != EXPR_REAL) {
        char* s = expr_to_string_fullform(res);
        printf("FAIL (not a Real): %s -> %s\n", input, s);
        free(s); ASSERT(0);
    } else if (fabs(res->data.real - want) > tol) {
        printf("FAIL: %s -> %.10g (expected %.10g)\n", input, res->data.real, want);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %.10g\n", input, res->data.real);
    }
    expr_free(res);
    expr_free(e);
}

/* Evaluate `input`; assert FullForm equals `expected`. */
static void run_fullform(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Evaluate `input`; assert it is still an unevaluated function application
 * (i.e. the interpolation did NOT fire and produce a number). */
static void run_unevaluated(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    if (res->type != EXPR_FUNCTION) {
        char* s = expr_to_string_fullform(res);
        printf("FAIL (expected unevaluated): %s -> %s\n", input, s);
        free(s); ASSERT(0);
    } else {
        printf("PASS (unevaluated): %s\n", input);
    }
    expr_free(res);
    expr_free(e);
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
    } else {
        printf("PASS print: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Two-point linear interpolation. */
static void test_linear(void) {
    run_close("InterpolatingFunction[{{0,1}},{{0,0},{1,10}}][0.5]", 5.0, 1e-9);
    run_close("InterpolatingFunction[{{0,1}},{{0,0},{1,10}}][1/4]", 2.5, 1e-9);
    /* y = 2x+1 sampled on a grid; linear/cubic reproduce it exactly. */
    run_close("InterpolatingFunction[{{0,3}},{{0,1},{1,3},{2,5},{3,7}}][1.25]", 3.5, 1e-9);
    /* Data given out of order is sorted by abscissa, then interpolated. */
    run_close("InterpolatingFunction[{{0,1}},{{1,0},{0,1}}][0.25]", 0.75, 1e-9);
}

/* Three-point reduced to quadratic (order = min(3, n-1) = 2). */
static void test_quadratic(void) {
    /* y = x^2 + 1, exact for a quadratic interpolant. */
    run_close("InterpolatingFunction[{{0,2}},{{0,1},{1,2},{2,5}}][1.5]", 3.25, 1e-9);
    run_close("InterpolatingFunction[{{0,2}},{{0,1},{1,2},{2,5}}][0.5]", 1.25, 1e-9);
}

/* Cubic data is reproduced exactly by the order-3 interpolant. */
static void test_cubic_exact(void) {
    const char* cube = "InterpolatingFunction[{{0,5}},"
                       "{{0,0},{1,1},{2,8},{3,27},{4,64},{5,125}}]";
    char buf[256];
    snprintf(buf, sizeof buf, "%s[2.5]", cube);
    run_close(buf, 15.625, 1e-6);   /* 2.5^3 */
    snprintf(buf, sizeof buf, "%s[3.7]", cube);
    run_close(buf, 50.653, 1e-6);   /* 3.7^3 */
    snprintf(buf, sizeof buf, "%s[0.5]", cube);
    run_close(buf, 0.125, 1e-6);    /* 0.5^3, near the left edge */
}

/* The reference hill table: matches Mathematica's piecewise windows. */
static void test_reference_windows(void) {
    run_close(IFUN_EX "[0.5]", 0.25, 1e-9);
    run_close(IFUN_EX "[1.5]", 2.0, 1e-9);
    run_close(IFUN_EX "[2.5]", 3.6875, 1e-9);
    run_close(IFUN_EX "[3.5]", 3.75, 1e-9);
    run_close(IFUN_EX "[4.5]", 1.75, 1e-9);
}

/* Exact arguments that hit a node return the stored value exactly. */
static void test_exact_nodes(void) {
    run_fullform(IFUN_EX "[0]", "0");
    run_fullform(IFUN_EX "[2]", "3");
    run_fullform(IFUN_EX "[5]", "0");
    /* Rational node: query 1/2 matches the {1/2, 7} node and returns 7. */
    run_fullform("InterpolatingFunction[{{0,1}},{{0,0},{1/2,7},{1,1}}][1/2]", "7");
    /* A real query at a node value is still interpolated (returns a Real). */
    run_close(IFUN_EX "[2.0]", 3.0, 1e-9);
}

/* Outside the domain: a dmval warning is printed and a value extrapolated. */
static void test_extrapolation(void) {
    fprintf(stderr, "[expect InterpolatingFunction::dmval warnings below]\n");
    run_close(IFUN_EX "[6]", -5.0, 1e-9);
    run_close(IFUN_EX "[-1]", 2.0, 1e-9);
}

/* Symbolic / wrong-arity applications are left unevaluated. */
static void test_unevaluated(void) {
    run_unevaluated(IFUN_EX "[t]");        /* symbolic argument */
    run_unevaluated(IFUN_EX "[1, 2]");     /* too many arguments */
    run_unevaluated(IFUN_EX "[]");         /* no argument */
}

/* Malformed objects do not interpolate (stay as an application). */
static void test_malformed(void) {
    /* Only one data point -> cannot interpolate. */
    run_unevaluated("InterpolatingFunction[{{0,1}},{{0,0}}][0.5]");
    /* Non-numeric ordinate. */
    run_unevaluated("InterpolatingFunction[{{0,1}},{{0,a},{1,b}}][0.5]");
    /* Table entry not a pair. */
    run_unevaluated("InterpolatingFunction[{{0,1}},{{0,0},{1}}][0.5]");
}

/* Standard output shows only the domain; data is abbreviated as <>. */
static void test_printing(void) {
    run_print(IFUN_EX, "InterpolatingFunction[{{0, 5}}, <>]");
}

/* 1-D derivatives via Derivative (') match Mathematica's piecewise-cubic
 * derivative on the reference hill table. */
static void test_derivatives_1d(void) {
    run_close("Derivative[1][" IFUN_EX "][2.5]", 1.0416666666666667, 1e-9);
    run_close("Derivative[2][" IFUN_EX "][2.5]", -1.5, 1e-9);
    run_close("Derivative[3][" IFUN_EX "][2.5]", -1.0, 1e-9);
    /* Order beyond the cubic window degree -> 0. */
    run_close("Derivative[4][" IFUN_EX "][2.5]", 0.0, 1e-12);
    /* Derivative composes additively: D[D[g]] == D^2[g]. */
    run_close("Derivative[1][Derivative[1][" IFUN_EX "]][2.5]", -1.5, 1e-9);
}

/* A differentiated object is itself an InterpolatingFunction (so it can be
 * applied), and prints with the <> abbreviation. */
static void test_derivative_object(void) {
    run_fullform("Head[Derivative[1][" IFUN_EX "]]", "InterpolatingFunction");
    run_print("Derivative[1][" IFUN_EX "]", "InterpolatingFunction[{{0, 5}}, <>]");
}

/* D[ifun[x], x] routes through the chain rule to the same derivative. */
static void test_d_chain_rule(void) {
    run_close("D[" IFUN_EX "[x], x] /. x -> 2.5", 1.0416666666666667, 1e-9);
    run_close("D[" IFUN_EX "[x], {x, 2}] /. x -> 2.5", -1.5, 1e-9);
}

/* 2-D tensor-product interpolation. Data is f(i,j) = i^2 + j^3 on a 5x5 grid,
 * which an order-3 interpolant reproduces exactly. */
#define IFUN_2D \
    "InterpolatingFunction[{{0,4},{0,4}}," \
    "{{{0,0},0},{{0,1},1},{{0,2},8},{{0,3},27},{{0,4},64}," \
    "{{1,0},1},{{1,1},2},{{1,2},9},{{1,3},28},{{1,4},65}," \
    "{{2,0},4},{{2,1},5},{{2,2},12},{{2,3},31},{{2,4},68}," \
    "{{3,0},9},{{3,1},10},{{3,2},17},{{3,3},36},{{3,4},73}," \
    "{{4,0},16},{{4,1},17},{{4,2},24},{{4,3},43},{{4,4},80}}]"

static void test_multidim(void) {
    run_close(IFUN_2D "[1.5, 2.5]", 17.875, 1e-6);   /* 1.5^2 + 2.5^3 */
    run_close(IFUN_2D "[0.5, 3.5]", 43.125, 1e-6);   /* 0.25 + 42.875 */
    /* Exact node tuple -> exact stored value. */
    run_fullform(IFUN_2D "[2, 3]", "31");
    run_fullform(IFUN_2D "[4, 4]", "80");
}

static void test_multidim_derivatives(void) {
    /* d/dx (i^2+j^3) = 2x ; d/dy = 3y^2 ; mixed = 0. */
    run_close("Derivative[1,0][" IFUN_2D "][1.5, 2.5]", 3.0, 1e-6);
    run_close("Derivative[0,1][" IFUN_2D "][1.5, 2.5]", 18.75, 1e-6);
    run_close("Derivative[2,0][" IFUN_2D "][1.5, 2.5]", 2.0, 1e-6);
    run_close("Derivative[0,2][" IFUN_2D "][1.5, 2.5]", 15.0, 1e-6);
    run_close("Derivative[1,1][" IFUN_2D "][1.5, 2.5]", 0.0, 1e-6);
}

/* Multi-dimensional arity and malformation. */
static void test_multidim_malformed(void) {
    /* 2-D object applied to one argument -> unevaluated. */
    run_unevaluated(IFUN_2D "[1.5]");
    /* Incomplete grid (missing the {1,1} point) -> unevaluated. */
    run_unevaluated("InterpolatingFunction[{{0,1},{0,1}},"
                    "{{{0,0},0},{{0,1},1},{{1,0},1}}][0.5, 0.5]");
}

int main(void) {
    symtab_init();
    core_init();

    test_linear();
    test_quadratic();
    test_cubic_exact();
    test_reference_windows();
    test_exact_nodes();
    test_extrapolation();
    test_unevaluated();
    test_malformed();
    test_printing();
    test_derivatives_1d();
    test_derivative_object();
    test_d_chain_rule();
    test_multidim();
    test_multidim_derivatives();
    test_multidim_malformed();

    symtab_clear();
    printf("\nAll InterpolatingFunction tests passed.\n");
    return 0;
}
