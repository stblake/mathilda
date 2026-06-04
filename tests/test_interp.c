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
#ifdef USE_MPFR
#include <mpfr.h>
#endif
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

/* Read a numeric result (EXPR_REAL or EXPR_MPFR) as a double; NAN otherwise. */
static double result_as_double(Expr* res, int* is_mpfr) {
    if (is_mpfr) *is_mpfr = 0;
    if (res->type == EXPR_REAL) return res->data.real;
#ifdef USE_MPFR
    if (res->type == EXPR_MPFR) { if (is_mpfr) *is_mpfr = 1; return mpfr_get_d(res->data.mpfr, MPFR_RNDN); }
#endif
    return NAN;
}

/* Evaluate; assert the result is numeric (Real or MPFR) within tol of want. */
static void run_val(const char* input, double want, double tol) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    double got = result_as_double(res, NULL);
    if (isnan(got)) {
        char* s = expr_to_string_fullform(res);
        printf("FAIL (not numeric): %s -> %s\n", input, s); free(s); ASSERT(0);
    } else if (fabs(got - want) > tol) {
        printf("FAIL: %s -> %.12g (expected %.12g)\n", input, got, want); ASSERT(0);
    } else printf("PASS: %s -> %.12g\n", input, got);
    expr_free(res); expr_free(e);
}

/* Evaluate; assert the result is a genuine arbitrary-precision MPFR (prec > 53)
 * within tol of want.  Falls back to a Real check when built without MPFR. */
static void run_mpfr(const char* input, double want, double tol) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
#ifdef USE_MPFR
    if (res->type != EXPR_MPFR) {
        char* s = expr_to_string_fullform(res);
        printf("FAIL (not MPFR): %s -> %s\n", input, s); free(s); ASSERT(0);
        expr_free(res); expr_free(e); return;
    }
    if (mpfr_get_prec(res->data.mpfr) <= 53) { printf("FAIL (prec<=53): %s\n", input); ASSERT(0); }
    double got = mpfr_get_d(res->data.mpfr, MPFR_RNDN);
#else
    double got = (res->type == EXPR_REAL) ? res->data.real : NAN;
#endif
    if (isnan(got) || fabs(got - want) > tol) {
        printf("FAIL mpfr: %s -> %.12g (expected %.12g)\n", input, got, want); ASSERT(0);
    } else printf("PASS mpfr: %s -> %.12g\n", input, got);
    expr_free(res); expr_free(e);
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

/* ---- Interpolation builtin ------------------------------------------- */

/* Form 1: bare function values f_i at x = 1, 2, ... */
static void test_interpolation_values(void) {
    /* Domain {{1, 6}}, abbreviated table. */
    run_print("Interpolation[{1,2,3,5,8,5}]",
              "InterpolatingFunction[{{1, 6}}, <>]");
    /* Sliding-window cubic value (matches Mathematica's 2.4375). */
    run_close("Interpolation[{1,2,3,5,8,5}][2.5]", 2.4375, 1e-9);
    /* Always passes through the data points exactly (Integer in -> Integer out). */
    run_fullform("Interpolation[{1,2,3,5,8,5}][3]", "3");
    run_fullform("Interpolation[{1,2,3,5,8,5}][1]", "1");
}

/* Form 2: {x, f} pairs with explicit abscissae. */
static void test_interpolation_xy(void) {
    /* Real abscissae give a Real domain. */
    run_print("Interpolation[{{0.,0.},{1.,1.},{2.,3.},{3.,4.},{4.,3.},{5.,0.}}]",
              "InterpolatingFunction[{{0.0, 5.0}}, <>]");
    /* Same hill data as IFUN_EX -> identical windowed value at 2.5. */
    run_close("Interpolation[{{0,0},{1,1},{2,3},{3,4},{4,3},{5,0}}][2.5]",
              3.6875, 1e-9);
    /* Non-uniform abscissae. */
    run_close("Interpolation[{{0,1},{1,3},{2,5},{3,7}}][1.25]", 3.5, 1e-9);
}

/* Form 3: multidimensional tensor-grid data. */
static void test_interpolation_multidim(void) {
    /* f(x,y) = x + y on a 2x2 grid; bilinear reproduces it exactly. */
    run_print("Interpolation[{{{0,0},0},{{0,1},1},{{1,0},1},{{1,1},2}}]",
              "InterpolatingFunction[{{0, 1}, {0, 1}}, <>]");
    run_close("Interpolation[{{{0,0},0},{{0,1},1},{{1,0},1},{{1,1},2}}][0.5,0.5]",
              1.0, 1e-9);
    /* 3x3 grid of x*y, evaluated off-node. */
    run_close("Interpolation[{{{1,1},1},{{1,2},2},{{1,3},3},"
              "{{2,1},2},{{2,2},4},{{2,3},6},"
              "{{3,1},3},{{3,2},6},{{3,3},9}}][1.5,2.5]", 3.75, 1e-9);
}

/* Interpolation[data, x] builds and immediately evaluates. */
static void test_interpolation_immediate(void) {
    run_close("Interpolation[{1,2,3,5,8,5}, 2.5]", 2.4375, 1e-9);
    run_close("Interpolation[{{0,0},{1,1},{2,3},{3,4},{4,3},{5,0}}, 2.5]",
              3.6875, 1e-9);
    /* Multidimensional immediate evaluation takes a coordinate list. */
    run_close("Interpolation[{{{0,0},0},{{0,1},1},{{1,0},1},{{1,1},2}}, {0.5,0.5}]",
              1.0, 1e-9);
}

/* InterpolationOrder option (0 = piecewise constant, 1 = linear, 3 = default). */
static void test_interpolation_order(void) {
    /* Order 0: value of the bracketing left node (x in [2,3) -> f(2)=2). */
    run_close("Interpolation[{1,2,3,5,8,5}, InterpolationOrder->0][2.5]", 2.0, 1e-9);
    /* Order 1: linear between f(2)=2 and f(3)=3. */
    run_close("Interpolation[{1,2,3,5,8,5}, InterpolationOrder->1][2.5]", 2.5, 1e-9);
    /* Order 3 (explicit) matches the default cubic value. */
    run_close("Interpolation[{1,2,3,5,8,5}, InterpolationOrder->3][2.5]", 2.4375, 1e-9);
    /* Exact node returns the stored value regardless of order. */
    run_close("Interpolation[{1,2,3,5,8,5}, InterpolationOrder->1][2.0]", 2.0, 1e-9);
}

/* Differentiation composes with an Interpolation-produced object. */
static void test_interpolation_derivative(void) {
    run_close("Interpolation[{1,2,3,5,8,5}]'[2.5]", 0.9583333333333333, 1e-9);
    run_close("D[Interpolation[{1,2,3,5,8,5}][x], x] /. x -> 2.5",
              0.9583333333333333, 1e-9);
    /* The order annotation survives differentiation (linear slope = 1). */
    run_close("Interpolation[{1,2,3,5,8,5}, InterpolationOrder->1]'[2.5]", 1.0, 1e-9);
}

/* Unsupported and malformed inputs are left unevaluated. */
static void test_interpolation_unsupported(void) {
    run_unevaluated("Interpolation[{1}]");                 /* too few points */
    run_unevaluated("Interpolation[{1,2,a}]");             /* non-numeric value */
    run_unevaluated("Interpolation[5]");                   /* not a list */
    run_unevaluated("Interpolation[{{0,0,1},{1,1,1}}]");   /* scalar-first triple (malformed) */
    run_unevaluated("Interpolation[{1,2,3}, PeriodicInterpolation->True]");
    /* Vector-valued samples are not supported yet. */
    run_unevaluated("Interpolation[{{0.,{1.,2.}},{1.,{3.,4.}}}]");
    /* An unknown Method is left unevaluated. */
    run_unevaluated("Interpolation[{1,2,3}, Method->\"Nonesuch\"]");
}

/* --- Method -> "Spline" ------------------------------------------------ */
static void test_method_spline(void) {
    /* Passes through every data node exactly. */
    run_fullform("Interpolation[{1,2,3,5,8,5}, Method->\"Spline\"][3]", "3");
    run_fullform("Interpolation[{1,2,3,5,8,5}, Method->\"Spline\"][1]", "1");
    /* n = 2 reduces to the straight line. */
    run_val("Interpolation[{{0,0},{1,10}}, Method->\"Spline\"][0.25]", 2.5, 1e-9);
    /* A spline of linear data is linear (y = 2(x-1) on x=1..4). */
    run_val("Interpolation[{0.,2.,4.,6.}, Method->\"Spline\"][1.5]", 1.0, 1e-9);
    /* Standard-form printing and method survival under differentiation. */
    run_print("Interpolation[{1,2,3,5,8,5}, Method->\"Spline\"]",
              "InterpolatingFunction[{{1, 6}}, <>]");
    run_print("Interpolation[{1,2,3,5,8,5}, Method->\"Spline\"]'",
              "InterpolatingFunction[{{1, 6}}, <>]");
    run_val("Interpolation[{0.,2.,4.,6.}, Method->\"Spline\"]'[1.5]", 2.0, 1e-9);
}

/* --- Method -> "Hermite" ---------------------------------------------- */
static void test_method_hermite(void) {
    /* Passes through nodes. */
    run_fullform("Interpolation[{1,4,9,16,25}, Method->\"Hermite\"][3]", "9");
    /* Reproduces quadratics exactly (x^2 sampled at x=1..5). */
    run_val("Interpolation[{1,4,9,16,25}, Method->\"Hermite\"][2.5]", 6.25, 1e-9);
    run_val("Interpolation[{1,4,9,16,25}, Method->\"Hermite\"][3.5]", 12.25, 1e-9);
    /* Reproduces linear data exactly. */
    run_val("Interpolation[{0.,2.,4.,6.,8.}, Method->\"Hermite\"][2.5]", 3.0, 1e-9);
    /* Derivative of x^2 interpolant: 2x at 2.5 = 5. */
    run_val("Interpolation[{1,4,9,16,25}, Method->\"Hermite\"]'[2.5]", 5.0, 1e-9);
}

/* --- 1-D derivative-supplied data (Hermite) --------------------------- */
/* f = x^3 with value + first derivative (3x^2) at x = 0..3. */
#define SUP_CUBIC "Interpolation[{{{0},0,0},{{1},1,3},{{2},8,12},{{3},27,27}}]"
/* f = x^5 with value + 1st (5x^4) + 2nd (20x^3) derivative at x = 0..3. */
#define SUP_QUINTIC "Interpolation[{{{0},0,0,0},{{1},1,5,20},{{2},32,80,160},{{3},243,405,540}}]"
static void test_supplied_1d(void) {
    /* Cubic Hermite reproduces the cubic exactly. */
    run_val(SUP_CUBIC "[1.5]", 3.375, 1e-9);
    run_val(SUP_CUBIC "[2.5]", 15.625, 1e-9);
    /* Derivative of the interpolant: 3x^2. */
    run_val(SUP_CUBIC "'[1.5]", 6.75, 1e-9);
    run_val("D[" SUP_CUBIC "[x], x] /. x -> 2.5", 18.75, 1e-9);
    /* Quintic Hermite (value + 2 derivatives) reproduces the quintic exactly. */
    run_val(SUP_QUINTIC "[1.5]", 7.59375, 1e-7);
    run_val(SUP_QUINTIC "'[1.5]", 25.3125, 1e-6);
    /* Real abscissae are accepted, with a {x} coordinate even in 1-D. */
    run_val("Interpolation[{{{0.},0.,0.},{{1.},1.,3.},{{2.},8.,12.}}][1.5]", 3.375, 1e-9);
}

/* --- n-D derivative-supplied data ------------------------------------- */
/* f = x^3 + y^3 with gradient {3x^2, 3y^2} on the 0..2 grid (cross term 0). */
#define SUP_2D_GRAD \
    "Interpolation[{{{0,0},0,{0,0}},{{0,1},1,{0,3}},{{0,2},8,{0,12}}," \
    "{{1,0},1,{3,0}},{{1,1},2,{3,3}},{{1,2},9,{3,12}}," \
    "{{2,0},8,{12,0}},{{2,1},9,{12,3}},{{2,2},16,{12,12}}}]"
/* f = x^2 + y^2 with gradient {2x,2y} and constant Hessian {{2,0},{0,2}}. */
#define SUP_2D_HESS \
    "Interpolation[{{{0,0},0,{0,0},{{2,0},{0,2}}},{{0,1},1,{0,2},{{2,0},{0,2}}}," \
    "{{0,2},4,{0,4},{{2,0},{0,2}}},{{1,0},1,{2,0},{{2,0},{0,2}}}," \
    "{{1,1},2,{2,2},{{2,0},{0,2}}},{{1,2},5,{2,4},{{2,0},{0,2}}}," \
    "{{2,0},4,{4,0},{{2,0},{0,2}}},{{2,1},5,{4,2},{{2,0},{0,2}}}," \
    "{{2,2},8,{4,4},{{2,0},{0,2}}}}]"
static void test_supplied_nd(void) {
    /* Separable cubic from gradient data: reproduced exactly. */
    run_val(SUP_2D_GRAD "[0.5, 0.5]", 0.25, 1e-9);   /* 0.125 + 0.125 */
    run_val(SUP_2D_GRAD "[1.5, 0.5]", 3.5, 1e-9);    /* 3.375 + 0.125 */
    run_val(SUP_2D_GRAD "[0.5, 1.5]", 3.5, 1e-9);
    /* Partial derivatives of the interpolant. */
    run_val("Derivative[1,0][" SUP_2D_GRAD "][1.5, 0.5]", 6.75, 1e-9);   /* 3x^2 */
    run_val("Derivative[0,1][" SUP_2D_GRAD "][0.5, 1.5]", 6.75, 1e-9);   /* 3y^2 */
    /* Exact node returns the stored value. */
    run_fullform(SUP_2D_GRAD "[1, 2]", "9");
    /* Hessian (biquintic) data reproduces x^2 + y^2. */
    run_val(SUP_2D_HESS "[0.5, 0.5]", 0.5, 1e-9);
    run_val(SUP_2D_HESS "[1.5, 0.5]", 2.5, 1e-9);    /* 2.25 + 0.25 */
    run_val("Derivative[1,0][" SUP_2D_HESS "][0.5, 0.5]", 1.0, 1e-9);   /* 2x */
    run_val("Derivative[0,1][" SUP_2D_HESS "][0.5, 0.5]", 1.0, 1e-9);   /* 2y */
}

/* --- arbitrary precision (MPFR) --------------------------------------- */
static void test_precision_mpfr(void) {
    /* High-precision argument -> MPFR result. */
    run_mpfr("Interpolation[{1,2,3,5,8,5}][N[5/2, 30]]", 2.4375, 1e-12);
    /* High-precision data -> MPFR result, agreeing with the machine value. */
    run_mpfr("Interpolation[N[{1,2,3,5,8,5}, 30]][N[5/2, 30]]", 2.4375, 1e-12);
    /* Spline / Hermite / default in MPFR all agree with the double path. */
    run_mpfr("Interpolation[N[{1,2,3,5,8,5}, 30], Method->\"Spline\"][N[5/2, 30]]",
             2.473086124401914, 1e-12);
    run_mpfr("Interpolation[N[{1,4,9,16,25}, 30], Method->\"Hermite\"][N[5/2, 30]]",
             6.25, 1e-12);
    /* Supplied-derivative cubic in MPFR. */
    run_mpfr("Interpolation[N[{{{0},0,0},{{1},1,3},{{2},8,12},{{3},27,27}}, 30]][N[3/2, 30]]",
             3.375, 1e-12);
    /* 2-D gradient data in MPFR. */
    run_mpfr("Interpolation[N[{{{0,0},0,{0,0}},{{0,1},1,{0,3}},{{0,2},8,{0,12}},"
             "{{1,0},1,{3,0}},{{1,1},2,{3,3}},{{1,2},9,{3,12}},"
             "{{2,0},8,{12,0}},{{2,1},9,{12,3}},{{2,2},16,{12,12}}}, 30]][N[1/2,30], N[1/2,30]]",
             0.25, 1e-12);
    /* Machine-precision input keeps returning a machine Real (not MPFR). */
    run_val("Interpolation[{1,2,3,5,8,5}][2.5]", 2.4375, 1e-12);
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

    test_interpolation_values();
    test_interpolation_xy();
    test_interpolation_multidim();
    test_interpolation_immediate();
    test_interpolation_order();
    test_interpolation_derivative();
    test_interpolation_unsupported();

    test_method_spline();
    test_method_hermite();
    test_supplied_1d();
    test_supplied_nd();
    test_precision_mpfr();

    symtab_clear();
    printf("\nAll InterpolatingFunction and Interpolation tests passed.\n");
    return 0;
}
