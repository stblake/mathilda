/* Tests for FunctionInterpolation.
 *
 * Every expected value is the interpolated function itself (Sin, Exp,
 * (E^x - 1)/x = Integrate[E^(x y), {y, 0, 1}], Sin[x] Cos[y]) evaluated
 * independently at points that are not sample nodes; the assertions bound
 * the interpolation error by the requested goals.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

static void assert_true(const char* input) {
    char* s = eval_str(input);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s  ->  %s", input, s);
    free(s);
}

/* Off-node probe points on [0, 1]: 0.0137 + 0.0491 k. */
#define PTS01 "Table[0.0137 + 0.0491 k, {k, 0, 20}]"

static void test_default_goals_sin(void) {
    /* default PrecisionGoal is 8 at machine precision */
    assert_true("f = FunctionInterpolation[Sin[x], {x, 0, 2 Pi}]; "
                "Head[f] === InterpolatingFunction && "
                "Max[Abs[f[#] - Sin[#]] & /@ (2 Pi " PTS01 ")] < 10^-7");
}

static void test_goal_controls_error(void) {
    /* Error tracks PrecisionGoal (with AccuracyGoal Infinity), and a tighter
     * goal needs more samples. */
    assert_true("err[pg_] := Module[{f = FunctionInterpolation[Exp[x], {x, 0, 1}, "
                "    PrecisionGoal -> pg]}, "
                "  {Max[Abs[(f[#] - Exp[#])/Exp[#]] & /@ " PTS01 "], Length[f[\"Grid\"]]}]; "
                "{e4, n4} = err[4]; {e7, n7} = err[7]; {e11, n11} = err[11]; "
                "And[e4 < 10^-4, e7 < 10^-7, e11 < 10^-11, n4 < n7 < n11]");
    /* AccuracyGoal alone (PrecisionGoal Infinity): absolute error bound */
    assert_true("f = FunctionInterpolation[100 Exp[x], {x, 0, 1}, AccuracyGoal -> 6, "
                "  PrecisionGoal -> Infinity]; "
                "Max[Abs[f[#] - 100 Exp[#]] & /@ " PTS01 "] < 10^-6");
}

static void test_interpolation_order(void) {
    assert_true("f = FunctionInterpolation[Cos[3 x], {x, 0, 1}, InterpolationOrder -> 7, "
                "  PrecisionGoal -> 12]; "
                "f[\"InterpolationOrder\"] == {7} && "
                "Max[Abs[f[#] - Cos[3 #]] & /@ " PTS01 "] < 10^-11");
}

static void test_inner_nintegrate(void) {
    /* HoldAll + Block-localised x: the inner NIntegrate sees numeric x. */
    assert_true("f = FunctionInterpolation[NIntegrate[Exp[x y], {y, 0, 1}], {x, 0, 1}]; "
                "Max[Abs[f[#] - (Exp[#] - 1)/#] & /@ " PTS01 "] < 10^-7");
}

static void test_high_interpolation_precision(void) {
    /* 30-digit values, goals beyond machine precision */
    assert_true("f = FunctionInterpolation[Exp[x], {x, 0, 1}, InterpolationPrecision -> 30, "
                "  PrecisionGoal -> 25, InterpolationOrder -> 8]; "
                "v = f[1/3]; "
                "Precision[v] >= 29 && Abs[v - Exp[1/3]] < 10^-24");
    assert_true("Max[Abs[f[#] - Exp[#]] & /@ Table[SetPrecision[1/7 + k/9, 30], {k, 0, 7}]] < 10^-24");
}

static void test_two_dimensional(void) {
    assert_true("f = FunctionInterpolation[Sin[x] Cos[y], {x, 0, 2}, {y, 0, 2}]; "
                "Max[Abs[f[#[[1]], #[[2]]] - Sin[#[[1]]] Cos[#[[2]]]] & /@ "
                "  Table[{0.0373 + 0.191 k, 1.9 - 0.177 k}, {k, 0, 9}]] < 10^-6");
}

static void test_block_semantics(void) {
    /* the variable's global value is shadowed during sampling and restored */
    assert_true("x = 7; g = FunctionInterpolation[x^2, {x, 0, 1}]; "
                "r = {g[0.5], x}; Clear[x]; Abs[r[[1]] - 0.25] < 10^-12 && r[[2]] == 7");
}

static void test_setoptions(void) {
    /* With no refinement (goals Infinity), the grid is exactly the initial
     * InterpolationPoints, whether given in the call or via SetOptions. */
    assert_true("Length[FunctionInterpolation[x, {x, 0, 1}, InterpolationPoints -> 17, "
                "  PrecisionGoal -> Infinity][\"Grid\"]] == 17");
    assert_true("old = Options[FunctionInterpolation]; "
                "SetOptions[FunctionInterpolation, InterpolationPoints -> 23, "
                "  PrecisionGoal -> Infinity]; "
                "n = Length[FunctionInterpolation[x, {x, 0, 1}][\"Grid\"]]; "
                "SetOptions[FunctionInterpolation, Sequence @@ old]; n == 23");
    assert_true("(InterpolationPoints /. Options[FunctionInterpolation]) == 11");
}

static void test_interpolant_under_precision_changes(void) {
    /* N / SetPrecision leave the object's data alone (NHoldAll), so an
     * integrand SetPrecision[k[x, y] f[y], 25] still evaluates. */
    assert_true("f = FunctionInterpolation[Exp[x], {x, 0, 1}, InterpolationPrecision -> 16]; "
                "s = SetPrecision[2 f[y], 25]; "
                "Abs[(s /. y -> 1/2) - 2 Exp[1/2]] < 10^-7 && "
                "Abs[NIntegrate[s, {y, 0, 1}, WorkingPrecision -> 25, PrecisionGoal -> 8] "
                "    - 2 (E - 1)] < 10^-7");
    /* a machine-precision table queried at a 25-digit point gives the
     * machine value (it used to stay unevaluated) */
    assert_true("g = FunctionInterpolation[Exp[x], {x, 0, 1}]; v = g[N[1/3, 25]]; "
                "MachineNumberQ[v] && Abs[v - Exp[1/3]] < 10^-7 && "
                "Abs[NIntegrate[g[y], {y, 0, 1}, WorkingPrecision -> 25, PrecisionGoal -> 8] "
                "    - (E - 1)] < 10^-7");
}

static void test_unevaluated_inputs(void) {
    test_mute_messages();
    /* non-real values, bad ranges and symbolic bounds stay unevaluated */
    assert_true("Head[FunctionInterpolation[Sqrt[x], {x, -1, 1}]] === FunctionInterpolation");
    assert_true("Head[FunctionInterpolation[x, {x, 1, 0}]] === FunctionInterpolation");
    assert_true("Head[FunctionInterpolation[x, {x, 0, a}]] === FunctionInterpolation");
    assert_true("MemberQ[Attributes[FunctionInterpolation], HoldAll]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_default_goals_sin);
    TEST(test_goal_controls_error);
    TEST(test_interpolation_order);
    TEST(test_inner_nintegrate);
    TEST(test_high_interpolation_precision);
    TEST(test_two_dimensional);
    TEST(test_block_semantics);
    TEST(test_setoptions);
    TEST(test_interpolant_under_precision_changes);
    TEST(test_unevaluated_inputs);
    printf("All funcinterp tests passed.\n");
    return 0;
}
