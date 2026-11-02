/* test_compound_variables.c
 *
 * Generalised variables in Series and GroebnerBasis.  Expected values
 * computed by the reference kernel.
 *
 *  1. Series / SeriesCoefficient accept a compound expansion variable --
 *     anything not a numeric quantity and not built by Plus, Times or Power --
 *     treated, like D does, as an independent symbol (FeynCalc's
 *     FCLoopGLIExpand expands in a GLI-valued variable).  x^2 and x + y stay
 *     Series::ivar.
 *  2. GroebnerBasis takes indexed variables x[i] and compound coefficients
 *     such as SPD[p, p] (FCLoopSingularityStructure): every maximal
 *     non-polynomial subexpression is an indeterminate.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_series_compound_variable(void) {
    assert_eval_eq("Series[f[g[q]], {g[q], 0, 2}] // InputForm",
                   "SeriesData[g[q], 0, {f[0], Derivative[1][f][0], Derivative[2][f][0]/2}, 0, 3, 1]", 0);
    assert_eval_eq("Series[Exp[a[1]], {a[1], 0, 3}] // InputForm",
                   "SeriesData[a[1], 0, {1, 1, 1/2, 1/6}, 0, 4, 1]", 0);
    assert_eval_eq("Series[h[x] + x, {h[x], 0, 1}] // InputForm",
                   "SeriesData[h[x], 0, {x, 1}, 0, 2, 1]", 0);
    assert_eval_eq("Normal[Series[Sin[Log[y]], {Log[y], 0, 3}]]", "Log[y] - 1/6 Log[y]^3", 0);
    assert_eval_eq("SeriesCoefficient[Exp[a[1]], {a[1], 0, 3}]", "1/6", 0);
    /* Plus / Times / Power and numeric quantities are not variables */
    assert_eval_eq("Quiet[Series[1/(1 - x^2), {x^2, 0, 2}]]", "Series[1/(1 - x^2), {x^2, 0, 2}]", 0);
    assert_eval_eq("Head[Quiet[Series[(x + y)^2, {x + y, 0, 2}]]]", "Series", 0);
    assert_eval_eq("Head[Quiet[Series[Exp[Sin[2]], {Sin[2], 0, 1}]]]", "Series", 0);
}

static void test_groebner_generalised_atoms(void) {
    assert_eval_eq("GroebnerBasis[{m^2 x[1]^2, 2 m^2 x[1]}, {x[1]}]", "{m^2 x[1]}", 0);
    assert_eval_eq("GroebnerBasis[{-SPD[p, p] x[1] x[2], -SPD[p, p] x[2], -SPD[p, p] x[1]}, {x[1], x[2]}]",
                   "{SPD[p, p] x[2], SPD[p, p] x[1]}", 0);
    assert_eval_eq("GroebnerBasis[{x[2]^2 x[3] + x[1]^2 x[4] + 2 x[1] x[2] x[4] + x[2]^2 x[4],"
                   " 2 x[1] x[4] + 2 x[2] x[4], 2 x[2] x[3] + 2 x[1] x[4] + 2 x[2] x[4], x[2]^2,"
                   " x[1]^2 + 2 x[1] x[2] + x[2]^2}, {x[1], x[2], x[3], x[4]}]",
                   "{x[2] x[3], x[2]^2, x[1] x[4] + x[2] x[4], x[1]^2 + 2 x[1] x[2]}", 0);
    assert_eval_eq("GroebnerBasis[{y^2 - x, x y - 1}, {x[1]}]", "{-1 + y^3, x - y^2}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_series_compound_variable);
    TEST(test_groebner_generalised_atoms);
    printf("All compound-variable tests passed.\n");
    return 0;
}
