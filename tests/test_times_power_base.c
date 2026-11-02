/* test_times_power_base.c
 *
 * Times merges a factor with a power of that same factor when the factor is
 * itself a power, as the reference does: m^2 (m^2)^(-1 + D/2) is
 * (m^2)^(D/2) (FeynCalc's FCFeynmanParametrize-ID13).  The base grouping
 * keyed m^2 as (m, 2) and (m^2)^a as (m^2, a), so the two never met.  Only
 * exact equality merges.  An inexact zero absorbs symbolic factors
 * (0. d is 0.).  Expected values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_power_base_merge(void) {
    assert_eval_eq("m^2 (m^2)^(-1 + D/2)", "(m^2)^(1/2 D)", 0);
    assert_eval_eq("(m^2)^a m^2", "(m^2)^(1 + a)", 0);
    assert_eval_eq("Sqrt[m^2] m^2", "(m^2)^(3/2)", 0);
    assert_eval_eq("x^2 (x^2)^a (x^2)^b", "(x^2)^(1 + a + b)", 0);
    /* not the same base: no merge */
    assert_eval_eq("m^4 Sqrt[m^2]", "m^4 Sqrt[m^2]", 0);
    assert_eval_eq("x^2 (x^2)^a x", "x^3 (x^2)^a", 0);
    assert_eval_eq("m^2 m^3", "m^5", 0);
}

static void test_inexact_zero_absorbs(void) {
    assert_eval_eq("0. d", "0.0", 0);
    assert_eval_eq("0. (d y)", "0.0", 0);
    assert_eval_eq("0. Sin[x]", "0.0", 0);
    assert_eval_eq("Quiet[0. Infinity]", "Indeterminate", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_power_base_merge);
    TEST(test_inexact_zero_absorbs);
    printf("All Times power-base tests passed.\n");
    return 0;
}
