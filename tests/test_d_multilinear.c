/* test_d_multilinear.c
 *
 * D differentiates Dot, NonCommutativeMultiply and Cross by the product rule,
 * keeping the factor order, as the reference does: D[a[x] . b[x], x] is
 * a[x] . b'[x] + a'[x] . b[x].  The generic chain rule gave
 * Derivative[1, 0][Dot][a[x], b[x]] a'[x], so FeynCalc's FunctionalD and
 * FeynRule (6 tests) came back with Derivative[..][Dot] or 0.  Cross with a
 * symbolic argument now stays unevaluated without Cross::nonn1.  Expected
 * values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_product_rule_multilinear(void) {
    assert_eval_eq("D[a[x] . b[x], x]", "a[x] . Derivative[1][b][x] + Derivative[1][a][x] . b[x]", 0);
    assert_eval_eq("D[a[x] . b[x] . c[x], x]",
                   "a[x] . b[x] . Derivative[1][c][x] + a[x] . Derivative[1][b][x] . c[x]"
                   " + Derivative[1][a][x] . b[x] . c[x]", 0);
    assert_eval_eq("D[a . b[x], x]", "a . Derivative[1][b][x]", 0);
    assert_eval_eq("D[a[x] . b[x], {x, 2}]",
                   "a[x] . Derivative[2][b][x] + 2 Derivative[1][a][x] . Derivative[1][b][x]"
                   " + Derivative[2][a][x] . b[x]", 0);
    assert_eval_eq("D[a[x] ** b[x], x]", "a[x] ** Derivative[1][b][x] + Derivative[1][a][x] ** b[x]", 0);
    assert_eval_eq("D[Cross[a[x], b[x]], x]",
                   "Cross[a[x], Derivative[1][b][x]] + Cross[Derivative[1][a][x], b[x]]", 0);
    /* a generic head keeps the chain rule */
    assert_eval_eq("D[h[a[x], b[x]], x]",
                   "Derivative[1][b][x] Derivative[0, 1][h][a[x], b[x]]"
                   " + Derivative[1][a][x] Derivative[1, 0][h][a[x], b[x]]", 0);
}

static void test_cross_symbolic_silent(void) {
    assert_eval_eq("Check[Cross[a, b], err]", "Cross[a, b]", 0);
    assert_eval_eq("Cross[{1, 2, 3}, {4, 5, 6}]", "{-3, 6, -3}", 0);
    assert_eval_eq("Quiet[Check[Cross[{1, 2}, {3, 4}], err]]", "err", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_product_rule_multilinear);
    TEST(test_cross_symbolic_silent);
    printf("All multilinear-D tests passed.\n");
    return 0;
}
