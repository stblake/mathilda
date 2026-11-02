/* test_rational_function_q.c
 *
 * Internal`RationalFunctionQ[expr, vars] did not exist, so FeynCalc's
 * PaVeLimitTo4 (which checks its D-dependent prefactors with it) aborted with
 * "The D-dependent prefactor must be a rational function of D".  Expected
 * values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_rational(void) {
    assert_eval_eq("Internal`RationalFunctionQ[(D-1)/(D-2), D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[(D-1)(D-2) x, D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[x Pair[a,b], D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[Sin[x], D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[(D-2)^(-2) s^(1/2), D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[1/(1+1/D), D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[(a+D)^-3, D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[1.5 D + I, D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[(D+x^(1/2))^2, D]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[x, {D, x}]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[x, {}]", "True", 0);
    assert_eval_eq("Internal`RationalFunctionQ[Sin[x]^2 + y, Sin[x]]", "True", 0);
}

static void test_not_rational(void) {
    assert_eval_eq("Internal`RationalFunctionQ[Sqrt[D], D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[2^D, D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[Sin[D], D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[f[D], D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[D^a, D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[D^(3/2), D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[D^2.5, D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[x Sin[y], {x, y}]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[{D, 1/D}, D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[D == 1, D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[f[D][x], D]", "False", 0);
    assert_eval_eq("Internal`RationalFunctionQ[x^D, D]", "False", 0);
    /* one argument: ::argr, unevaluated */
    assert_eval_eq("Quiet[Internal`RationalFunctionQ[x]]", "Internal`RationalFunctionQ[x]", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_rational);
    TEST(test_not_rational);
    printf("All RationalFunctionQ tests passed.\n");
    return 0;
}
