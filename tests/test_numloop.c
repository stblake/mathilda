/*
 * test_numloop.c -- the automatic numeric loop fast-path (src/numloop.c).
 *
 * Two things are checked:
 *   1. DIFFERENTIAL: for every fast-pathed loop shape, the compiled double
 *      result AGREES with the interpreter's (fast path forced off via
 *      numloop_set_enabled) to floating-point rounding. It is not bit-identical:
 *      the interpreter's Orderless Plus/Times sort operands by their runtime
 *      values before folding, so the operation order (and thus the last ULP) is
 *      data-dependent and cannot be reproduced by a static compile. Both are
 *      valid IEEE evaluations; we require a tight relative agreement so any real
 *      structural bug (wrong op/constant -> macroscopic error) is caught.
 *   2. FALLBACK: exact / symbolic loops that must NOT be fast-pathed still
 *      produce their exact/symbolic result unchanged.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "numloop.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Evaluate `input` to a machine Real and return the raw double bits. Fails the
 * test if the result is not an EXPR_REAL. */
static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    ASSERT_MSG(r && r->type == EXPR_REAL,
               "expected Real result for %s (got head %d)", input,
               r ? (int)r->type : -1);
    double v = r->data.real;
    expr_free(e);
    expr_free(r);
    return v;
}

/* Assert the fast path and the interpreter agree to floating-point rounding.
 * A genuine compile bug (wrong op, wrong constant, missing term) produces a
 * macroscopic error; ULP-level operand-order differences do not. */
static void diff(const char* input) {
    numloop_set_enabled(true);
    double fast = eval_real(input);
    numloop_set_enabled(false);
    double interp = eval_real(input);
    numloop_set_enabled(true);
    double scale = fabs(interp) > 1.0 ? fabs(interp) : 1.0;
    double rel = fabs(fast - interp) / scale;
    ASSERT_MSG(rel < 1e-9,
               "fast/interp mismatch for %s: fast=%.17g interp=%.17g rel=%.3g",
               input, fast, interp, rel);
}

/* Assert the fast path and interpreter agree element-wise for a float64 NDArray
 * result (small-array fused fast-path vs the interpreter's vectorized kernels). */
static void diff_array(const char* input) {
    numloop_set_enabled(true);
    Expr* e1 = parse_expression(input); Expr* r1 = evaluate(e1);
    numloop_set_enabled(false);
    Expr* e2 = parse_expression(input); Expr* r2 = evaluate(e2);
    numloop_set_enabled(true);
    ASSERT_MSG(r1 && r1->type == EXPR_NDARRAY && r2 && r2->type == EXPR_NDARRAY &&
               r1->data.ndarray.dtype == NDT_FLOAT64 &&
               r2->data.ndarray.dtype == NDT_FLOAT64,
               "%s: expected float64 NDArray results", input);
    size_t n1 = 1, n2 = 1;
    for (int i = 0; i < r1->data.ndarray.rank; i++) n1 *= (size_t)r1->data.ndarray.dims[i];
    for (int i = 0; i < r2->data.ndarray.rank; i++) n2 *= (size_t)r2->data.ndarray.dims[i];
    ASSERT_MSG(n1 == n2, "%s: shape mismatch %zu vs %zu", input, n1, n2);
    const double* d1 = (const double*)r1->data.ndarray.data;
    const double* d2 = (const double*)r2->data.ndarray.data;
    for (size_t k = 0; k < n1; k++) {
        double scale = fabs(d2[k]) > 1.0 ? fabs(d2[k]) : 1.0;
        ASSERT_MSG(fabs(d1[k] - d2[k]) / scale < 1e-9,
                   "%s: element %zu fast=%.17g interp=%.17g", input, k, d1[k], d2[k]);
    }
    expr_free(e1); expr_free(r1); expr_free(e2); expr_free(r2);
}

/* Assert an expression evaluates to `expected` (FullForm) with the fast path
 * enabled -- used for fallback (exact/symbolic) cases. */
static void expect_full(const char* input, const char* expected) {
    numloop_set_enabled(true);
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------------- Differential: Nest ---------------- */
static void test_diff_nest_logistic(void) {
    diff("Nest[3.5 # (1 - #)&, 1./Pi, 1000]");
}
static void test_diff_nest_cos(void) {
    diff("Nest[Cos, 1.0, 200]");
}
static void test_diff_nest_real_seed_pow(void) {
    /* integer-exponent Power via Power[#,2] */
    diff("Nest[0.5 #^2 + 0.5&, 0.3, 500]");
}
static void test_diff_nest_div(void) {
    /* #/2 parses to Times[Slot[1], Power[2,-1]] */
    diff("Nest[#/2 + 0.1&, 1.0, 300]");
}
static void test_diff_nest_transcendental(void) {
    diff("Nest[Exp[-#] + 0.1&, 0.5, 100]");
}
static void test_diff_nest_body_real_forces(void) {
    /* exact integer seed, but Real literal in body -> inexact result */
    diff("Nest[3.5 # (1 - #)&, 1/3, 400]");
}
static void test_diff_nest_sqrt(void) {
    diff("Nest[Sqrt[# + 1.0]&, 2.0, 60]");
}

/* ---------------- Differential: Do count form ---------------- */
static void test_diff_do_logistic(void) {
    diff("Module[{x = 1/Pi}, Do[x = 3.5 x (1 - x), {5000}]; x]");
}
static void test_diff_do_real_seed(void) {
    diff("Module[{x = 0.5}, Do[x = x^2 + 0.1, {100}]; x]");
}
static void test_diff_do_trig(void) {
    diff("Module[{x = 0.2}, Do[x = Sin[x] + 0.3, {500}]; x]");
}

/* ---------------- Differential: compound (multi-statement) bodies ---------------- */
static void test_diff_do_compound(void) {
    /* three statements, a temporary global y, sequential within an iteration */
    diff("Module[{x = 1/Pi}, Do[x = 3.5 x (1 - x); y = 4 x; x = y/4., {5000}]; x]");
}
static void test_diff_do_compound_multivar(void) {
    /* two persistent state variables updated with a temp */
    diff("Module[{a = 0., b = 1.}, Do[t = a + b; a = b; b = t, {40}]; a + b]");
}
static void test_diff_do_range_sumsq(void) {
    /* integer range form with a real accumulator */
    diff("Module[{s = 0.}, Do[s = s + i^2, {i, 1, 500}]; s]");
}
static void test_diff_for_compound(void) {
    diff("Module[{x = 1/Pi}, For[i = 0, i < 5000, i++, x = 3.5 x (1 - x); y = 4 x; x = y/4.]; x]");
}
static void test_diff_while_compound(void) {
    diff("Module[{s = 0.}, While[s < 100000., s = s + 1.; s = s*1.0000001]; s]");
}

/* ---------------- Differential: For ---------------- */
static void test_diff_for_logistic(void) {
    diff("Module[{x = 1/Pi}, For[i = 0, i < 5000, i++, x = 3.5 x (1 - x)]; x]");
}
static void test_diff_for_uses_counter(void) {
    /* rhs references the counter i as well as accumulator x */
    diff("Module[{x = 0.0}, For[i = 1, i <= 100, i++, x = x + 1.0/i]; x]");
}

/* ---------------- Differential: While ---------------- */
static void test_diff_while(void) {
    diff("Module[{s = 1.0}, While[s < 1000000, s = s 1.5 + 1.0]; s]");
}

/* ---------------- Differential: Fold ---------------- */
static void test_diff_fold_sumsq(void) {
    diff("Fold[#1 + #2^2 &, 0., Range[400]]");
}
static void test_diff_fold_logistic(void) {
    diff("Fold[3.5 #1 (1 - #1) + 0 #2 &, 0.2, Range[600]]");
}
static void test_diff_fold_real_elems(void) {
    /* exact seed, but Real list elements force inexactness */
    diff("Fold[#1 + #2 &, 0, Range[100.]]");
}

/* ---------------- Differential: FixedPoint ---------------- */
static void test_diff_fixedpoint_cos(void) {
    diff("FixedPoint[Cos, 1.0]");
}
static void test_diff_fixedpoint_sqrt(void) {
    /* Newton iteration for sqrt(2) */
    diff("FixedPoint[(# + 2/#)/2 &, 1.0]");
}

/* ---------------- Differential: NestWhile ---------------- */
static void test_diff_nestwhile_halve(void) {
    diff("NestWhile[#/2 &, 1000000., # > 1 &]");
}
static void test_diff_nestwhile_double(void) {
    diff("NestWhile[# 2. &, 1., # < 1000000 &]");
}
static void test_diff_nest_cos_head(void) {
    /* bare function head (not a Function[...]) is accelerated too */
    diff("Nest[Cos, 1.0, 500]");
}

/* ---------------- Differential: NDArray (fused element-wise) ---------------- */
static void test_diff_nest_array_logistic(void) {
    diff_array("Nest[3.5 # (1 - #)&, NDArray[{0.1, 0.3, 0.5, 0.7, 0.9}], 300]");
}
static void test_diff_nest_array_cos(void) {
    diff_array("Nest[Cos, NDArray[{0., 1., 2., 3.}], 50]");
}
static void test_diff_nest_array_2d(void) {
    diff_array("Nest[# + 0.5&, NDArray[{{0.1, 0.2}, {0.3, 0.4}}], 10]");
}
static void test_diff_do_array_logistic(void) {
    diff_array("Module[{a = NDArray[{0.1, 0.3, 0.5, 0.7}]}, Do[a = 3.5 a (1 - a), {300}]; a]");
}
static void test_diff_do_array_compound(void) {
    diff_array("Module[{a = NDArray[{0.1, 0.3, 0.5, 0.7}]}, "
               "Do[a = 3.5 a (1 - a); b = a + 1.; a = b - 1., {200}]; a]");
}

/* ---------------- Differential: in-place Part-assignment loops ---------------- */
static void test_diff_do_part_fill(void) {
    diff_array("Module[{v = NDArray[Table[0., {40}]]}, Do[v[[i]] = N[i]^2, {i, 1, 40}]; v]");
}
static void test_diff_do_part_reads(void) {
    /* rhs reads other elements (Fibonacci) */
    diff_array("Module[{v = NDArray[Table[0., {20}]]}, v[[1]] = 1.; v[[2]] = 1.; "
               "Do[v[[i]] = v[[i-1]] + v[[i-2]], {i, 3, 20}]; v]");
}
static void test_diff_for_part(void) {
    diff_array("Module[{v = NDArray[Table[0., {30}]]}, "
               "For[i = 1, i <= 30, i++, v[[i]] = Sin[N[i]]]; v]");
}
static void test_diff_do_part_2d(void) {
    diff_array("Module[{m = NDArray[Table[0., {5}, {5}]]}, "
               "Do[m[[i, j]] = N[i]*10 + N[j], {i, 1, 5}, {j, 1, 5}]; m]");
}
static void test_diff_do_part_2d_reads(void) {
    /* multi-index reads of other elements, mutated in place */
    diff_array("Module[{m = NDArray[Table[1., {4}, {4}]]}, "
               "Do[m[[i, j]] = m[[i, j]] + m[[j, i]], {i, 1, 4}, {j, 1, 4}]; m]");
}

/* ---------------- Map[f, list] (fast-path + fallback) ---------------- */
static void test_map_fast(void) {
    expect_full("Map[2. #&, {1., 2., 3.}]", "List[2.0, 4.0, 6.0]");
}
static void test_map_body_real(void) {
    /* exact integer elements, but a Real literal forces an inexact result */
    expect_full("Map[# + 0.5&, {1, 2, 3}]", "List[1.5, 2.5, 3.5]");
}
static void test_map_bare_head(void) {
    expect_full("Map[Sqrt, {4., 9., 16.}]", "List[2.0, 3.0, 4.0]");
}
static void test_map_fallback_exact(void) {
    expect_full("Map[#^2&, {1, 2, 3}]", "List[1, 4, 9]");
}
static void test_map_fallback_symbolic(void) {
    expect_full("Map[f, {1, 2, 3}]", "List[f[1], f[2], f[3]]");
}
static void test_map_fallback_free_symbol(void) {
    expect_full("Map[# + a&, {1., 2.}]", "List[Plus[1.0, a], Plus[2.0, a]]");
}

/* ---------------- Fallback: must NOT fast-path ---------------- */
static void test_fallback_nest_exact_int(void) {
    /* exact integer arithmetic stays exact */
    expect_full("Nest[#^2&, 2, 4]", "65536");
}
static void test_fallback_nest_symbolic(void) {
    /* Sin of an exact integer stays symbolic */
    expect_full("Nest[Sin, 2, 3]", "Sin[Sin[Sin[2]]]");
}
static void test_fallback_do_exact_int(void) {
    expect_full("Module[{x = 2}, Do[x = x + 1, {5}]; x]", "7");
}
static void test_fallback_nest_symbolic_body(void) {
    /* free symbol in body -> not numeric-closed */
    expect_full("Nest[# + a&, 0, 3]", "Times[3, a]");
}
static void test_fallback_nest_zero_times(void) {
    /* n = 0 returns the seed unchanged (still inexact here) */
    expect_full("Nest[3.5 # (1 - #)&, 0.25, 0]", "0.25");
}
static void test_fallback_fold_exact(void) {
    /* exact integer Fold stays exact */
    expect_full("Fold[#1 + #2 &, 0, {1, 2, 3, 4}]", "10");
}
static void test_fallback_nestwhile_predicate(void) {
    /* IntegerQ is not a numeric comparison -> falls back, stays exact */
    expect_full("NestWhile[#/2 &, 1024, IntegerQ]", "Rational[1, 2]");
}
static void test_fallback_nestwhile_exact_int(void) {
    /* exact integer seed + integer arithmetic -> falls back, stays exact */
    expect_full("NestWhile[# - 1 &, 10, # > 0 &]", "0");
}
static void test_fallback_fixedpoint_exact(void) {
    /* exact integer FixedPoint that terminates via the interpreter */
    expect_full("FixedPoint[Floor[#/2] &, 100]", "0");
}
static void test_fallback_do_compound_symbolic(void) {
    /* a free symbol in a compound body -> not numeric-closed, stays symbolic */
    expect_full("Module[{x = 1.0}, Do[x = x + a; x = x*2, {2}]; x]",
                "Times[2, Plus[a, Times[2, Plus[1.0, a]]]]");
}
static void test_fallback_do_compound_nonset(void) {
    /* a non-assignment statement (Print) in the body -> interpreter */
    expect_full("Module[{x = 0}, Do[x = x + 1; x, {3}]; x]", "3");
}

int main(void) {
    symtab_init();
    core_init();

    /* Differential (bit-identical fast vs interpreted) */
    TEST(test_diff_nest_logistic);
    TEST(test_diff_nest_cos);
    TEST(test_diff_nest_real_seed_pow);
    TEST(test_diff_nest_div);
    TEST(test_diff_nest_transcendental);
    TEST(test_diff_nest_body_real_forces);
    TEST(test_diff_nest_sqrt);
    TEST(test_diff_do_logistic);
    TEST(test_diff_do_real_seed);
    TEST(test_diff_do_trig);
    TEST(test_diff_for_logistic);
    TEST(test_diff_for_uses_counter);
    TEST(test_diff_while);
    TEST(test_diff_do_compound);
    TEST(test_diff_do_compound_multivar);
    TEST(test_diff_do_range_sumsq);
    TEST(test_diff_for_compound);
    TEST(test_diff_while_compound);
    TEST(test_diff_fold_sumsq);
    TEST(test_diff_fold_logistic);
    TEST(test_diff_fold_real_elems);
    TEST(test_diff_fixedpoint_cos);
    TEST(test_diff_fixedpoint_sqrt);
    TEST(test_diff_nestwhile_halve);
    TEST(test_diff_nestwhile_double);
    TEST(test_diff_nest_cos_head);
    TEST(test_diff_nest_array_logistic);
    TEST(test_diff_nest_array_cos);
    TEST(test_diff_nest_array_2d);
    TEST(test_diff_do_array_logistic);
    TEST(test_diff_do_array_compound);
    TEST(test_diff_do_part_fill);
    TEST(test_diff_do_part_reads);
    TEST(test_diff_for_part);
    TEST(test_diff_do_part_2d);
    TEST(test_diff_do_part_2d_reads);
    TEST(test_map_fast);
    TEST(test_map_body_real);
    TEST(test_map_bare_head);
    TEST(test_map_fallback_exact);
    TEST(test_map_fallback_symbolic);
    TEST(test_map_fallback_free_symbol);

    /* Fallback (exact / symbolic must be untouched) */
    TEST(test_fallback_nest_exact_int);
    TEST(test_fallback_nest_symbolic);
    TEST(test_fallback_do_exact_int);
    TEST(test_fallback_nest_symbolic_body);
    TEST(test_fallback_nest_zero_times);
    TEST(test_fallback_fold_exact);
    TEST(test_fallback_nestwhile_predicate);
    TEST(test_fallback_nestwhile_exact_int);
    TEST(test_fallback_fixedpoint_exact);
    TEST(test_fallback_do_compound_symbolic);
    TEST(test_fallback_do_compound_nonset);

    printf("All numloop tests passed!\n");
    return 0;
}
