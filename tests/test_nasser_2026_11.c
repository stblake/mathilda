/* test_nasser_2026_11.c
 *
 * Reports from Nasser's 2026-11-02 session (v0.6142), each a wrong or missing
 * answer; expected values are the reference kernel's.
 *
 *   1. CoefficientArrays of a system of EQUATIONS read each lhs == rhs as a
 *      variable-free constant (a rank-1 SparseArray of the equations), so
 *      {b, m} = CoefficientArrays[eqns, vars]; LinearSolve[m, -b] failed.
 *      lhs == rhs is the polynomial lhs - rhs.
 *   1b. The failed destructuring {b, m} = {one element} stayed unevaluated in
 *      silence; it is Set::shape and evaluates to the right side.
 *   2. Integrate[DiracDelta[x] Cos[x], {x, -Infinity, Infinity}] stayed
 *      unevaluated because Positive[Infinity] was False, so the sifting rule
 *      could not place the impulse inside an infinite interval.
 *   3. SeriesCoefficient[series, n] -- the coefficient read out of an existing
 *      SeriesData -- was not implemented.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_coefficientarrays_equations(void) {
    assert_eval_eq("Normal /@ CoefficientArrays[{2 x + 3 y + 4 z == 5, 9 x + 8 y + 7 z == 6, "
                   "x + 9 y + 2 z == 8}, {x, y, z}]",
                   "{{-5, -6, -8}, {{2, 3, 4}, {9, 8, 7}, {1, 9, 2}}}", 0);
    assert_eval_eq("Module[{b, m}, {b, m} = CoefficientArrays[{2 x + 3 y + 4 z == 5, "
                   "9 x + 8 y + 7 z == 6, x + 9 y + 2 z == 8}, {x, y, z}]; LinearSolve[m, -b]]",
                   "{-13/15, 11/15, 17/15}", 0);
    assert_eval_eq("Normal /@ CoefficientArrays[x + 2 y == 3, {x, y}]", "{-3, {1, 2}}", 0);
    assert_eval_eq("Normal /@ CoefficientArrays[{x^2 == 1, x y == 3}, {x, y}]",
                   "{{-1, -3}, {{0, 0}, {0, 0}}, {{{1, 0}, {0, 0}}, {{0, 1}, {0, 0}}}}", 0);
    /* the stored value is the arrays, not the rewritten call */
    assert_eval_eq("Module[{c}, c = CoefficientArrays[{x == 1, y == 3}, {x, y}]; Normal[c[[2]]]]",
                   "{{1, 0}, {0, 1}}", 0);
    /* an inequality is not a polynomial: ::poly, unevaluated */
    assert_eval_eq("Quiet[Head[CoefficientArrays[{x + 2 y == 3, x < 2}, {x, y}]]]",
                   "CoefficientArrays", 0);
    assert_eval_eq("Check[CoefficientArrays[{x < 2}, {x}], \"msg\"]", "\"msg\"", 0);
}

static void test_set_shape(void) {
    assert_eval_eq("Module[{u, v, r}, r = Quiet[{u, v} = {1, 2, 3}]; {r, u, v} === {{1, 2, 3}, u, v}]",
                   "True", 0);
    assert_eval_eq("Check[{s1, s2} = {1}, \"caught\"]", "\"caught\"", 0);
    assert_eval_eq("Quiet[{p1, p2} := {1}]", "$Failed", 0);
    assert_eval_eq("Check[{q1, {q2, q3}} = {1, {2}}, \"caught\"]", "\"caught\"", 0);
    /* matching shapes and threading are unchanged */
    assert_eval_eq("Module[{g1, g2}, {g1, g2} = {5, 6}; {g1, g2}]", "{5, 6}", 0);
    assert_eval_eq("Module[{h1, h2}, {h1, h2} = 7; {h1, h2}]", "{7, 7}", 0);
}

static void test_sign_predicates_infinity(void) {
    assert_eval_eq("{Positive[Infinity], Negative[-Infinity], Positive[-Infinity], "
                   "NonNegative[Infinity], NonPositive[-Infinity], Positive[I Infinity], "
                   "Negative[Infinity], Positive[2 Infinity], NonNegative[-Infinity]}",
                   "{True, True, False, True, True, False, False, True, False}", 0);
    assert_eval_eq("Positive[ComplexInfinity]", "Positive[ComplexInfinity]", 0);
    assert_eval_eq("Positive[DirectedInfinity[1]]", "True", 0);
}

static void test_dirac_infinite_interval(void) {
    assert_eval_eq("Integrate[DiracDelta[x] Cos[x], {x, -Infinity, Infinity}]", "1", 0);
    assert_eval_eq("Integrate[DiracDelta[x - 1] x^2, {x, 0, Infinity}]", "1", 0);
    assert_eval_eq("Integrate[DiracDelta[x + 1] x^2, {x, 0, Infinity}]", "0", 0);
    assert_eval_eq("Integrate[DiracDelta[2 x - 1] Exp[x], {x, -Infinity, 0}]", "0", 0);
}

static void test_seriescoefficient_seriesdata(void) {
    assert_eval_eq("SeriesCoefficient[Series[Exp[Sin[x]], {x, 0, 10}], 8]", "31/5760", 0);
    assert_eval_eq("SeriesCoefficient[Series[Exp[Sin[x]], {x, 0, 10}], 11]", "Indeterminate", 0);
    assert_eval_eq("SeriesCoefficient[Series[Exp[Sin[x]], {x, 0, 10}], -1]", "0", 0);
    assert_eval_eq("SeriesCoefficient[Series[Exp[x]/x^2, {x, 0, 4}], -2]", "1", 0);
    assert_eval_eq("SeriesCoefficient[Series[Sqrt[x] Exp[x], {x, 0, 3}], 3/2]", "1", 0);
    assert_eval_eq("SeriesCoefficient[Series[Sqrt[x] Exp[x], {x, 0, 3}], 1]", "0", 0);
    assert_eval_eq("SeriesCoefficient[Series[Exp[x y], {x, 0, 3}, {y, 0, 3}], {2, 2}]", "1/2", 0);
    /* a plain polynomial is not a series: stays unevaluated, as the reference */
    assert_eval_eq("SeriesCoefficient[x + x^2, 2]", "SeriesCoefficient[x + x^2, 2]", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    TEST(test_coefficientarrays_equations);
    TEST(test_set_shape);
    TEST(test_sign_predicates_infinity);
    TEST(test_dirac_infinite_interval);
    TEST(test_seriescoefficient_seriesdata);
    printf("All Nasser 2026-11 tests passed.\n");
    return 0;
}
