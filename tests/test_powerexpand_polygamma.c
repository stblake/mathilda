/* test_powerexpand_polygamma.c
 *
 * PowerExpand of the log of a rational, and PolyGamma at integer and
 * half-integer arguments (FeynCalc's SimplifyPolyLog, 9 tests).  Expected
 * values computed by the reference kernel.
 *
 *  1. PowerExpand[Log[q]] splits a rational q into prime logs (I Pi for q < 0).
 *  2. PolyGamma auto-evaluates exactly where the reference does: integer
 *     arguments and odd orders at half-integers (closed zeta), with n! kept
 *     outside the sum.  Digamma at a half-integer stays symbolic (it was
 *     evaluated, so SimplifyPolyLog's PolyGamma rule never fired).
 *  3. FunctionExpand gives the zeta(odd) forms at integers and half-integers
 *     and the Gauss digamma closed form at every positive rational.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_powerexpand_log_rational(void) {
    assert_eval_eq("PowerExpand[Log[1048576]]", "20 Log[2]", 0);
    assert_eval_eq("PowerExpand[Log[12]]", "2 Log[2] + Log[3]", 0);
    assert_eval_eq("PowerExpand[Log[3/4]]", "-2 Log[2] + Log[3]", 0);
    assert_eval_eq("PowerExpand[Log[1/8]]", "-3 Log[2]", 0);
    assert_eval_eq("PowerExpand[Log[-8]]", "I Pi + 3 Log[2]", 0);
    assert_eval_eq("PowerExpand[Log[7]]", "Log[7]", 0);
    assert_eval_eq("PowerExpand[Log[12 x^2]]", "2 Log[2] + Log[3] + 2 Log[x]", 0);
}

static void test_polygamma_auto(void) {
    assert_eval_eq("PolyGamma[3, 2]", "6 (-1 + 1/90 Pi^4)", 0);
    assert_eval_eq("PolyGamma[1, 4]", "-49/36 + 1/6 Pi^2", 0);
    assert_eval_eq("PolyGamma[3, 1]", "1/15 Pi^4", 0);
    assert_eval_eq("PolyGamma[1, 1/2]", "1/2 Pi^2", 0);
    assert_eval_eq("PolyGamma[1, 3/2]", "-4 + 1/2 Pi^2", 0);
    assert_eval_eq("PolyGamma[3, 3/2]", "6 (-16 + 1/6 Pi^4)", 0);
    assert_eval_eq("PolyGamma[0, 3/2]", "PolyGamma[0, 3/2]", 0);
    assert_eval_eq("PolyGamma[2, 3]", "PolyGamma[2, 3]", 0);
    assert_eval_eq("PolyGamma[0, 3]", "3/2 - EulerGamma", 0);
}

static void test_polygamma_functionexpand(void) {
    assert_eval_eq("FunctionExpand[PolyGamma[2, 3]]", "2 (9/8 - Zeta[3])", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[4, 2]]", "24 (1 - Zeta[5])", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[2, 1]]", "-2 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[2, 5/2]]", "448/27 - 14 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[4, 3/2]]", "768 - 744 Zeta[5]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[0, 1/2]]", "-EulerGamma - Log[4]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[0, 7/2]]", "46/15 - EulerGamma - Log[4]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[0, 3/4]]", "-EulerGamma + 1/2 Pi - Log[8]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[0, 1/3]]",
                   "-EulerGamma - Pi/(2 Sqrt[3]) - 3/2 Log[3]", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_powerexpand_log_rational);
    TEST(test_polygamma_auto);
    TEST(test_polygamma_functionexpand);
    printf("All PowerExpand / PolyGamma tests passed.\n");
    return 0;
}
