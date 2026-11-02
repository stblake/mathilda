/* test_det_expand_simplify_eqs.c
 *
 * Det of a symbolic matrix is the expanded polynomial (FeynCalc's
 * EpsContract squares an Eps through Det with no Expand after it, 3 tests),
 * and Simplify uses the equations among its assumptions.  Expected values
 * computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_det_expanded(void) {
    assert_eval_eq("Det[{{a, b, c}, {d, e, f}, {g, h, i}}]",
                   "-c e g + b f g + c d h - a f h - b d i + a e i", 0);
    assert_eval_eq("Det[{{a + b, c}, {d, e}}]", "-c d + a e + b e", 0);
    assert_eval_eq("Det[{{D, 1, 1}, {1, D, 1}, {1, 1, D}}]", "2 - 3 D + D^3", 0);
    assert_eval_eq("Det[{{x + 1, x}, {x, x - 1}}]", "-1", 0);
}

static void test_simplify_with_equations(void) {
    assert_eval_eq("Simplify[x[1] + x[2], x[1] + x[2] == 1]", "1", 0);
    assert_eval_eq("Simplify[a + b, a + b == c]", "c", 0);
    assert_eval_eq("Simplify[x^2 + y, x == y]", "y (1 + y)", 0);
    assert_eval_eq("Simplify[a x + a y, x + y == 1]", "a", 0);
    assert_eval_eq("Simplify[x y + y z, x + z == 2]", "2 y", 0);
    assert_eval_eq("Simplify[x[1] (a x[1] + b x[2]), x[1] + x[2] == 1]",
                   "x[1] (a x[1] + b x[2])", 0);
}

/* A number factor is simplified with the factors it can distribute over
 * (v0.6140): the multiplicative splitter used to lift 1/3 out of
 * (a + 3 b)/(3 b) and answer (3 + a/b)/3.  FeynCalc's FCFactorOut. */
static void test_simplify_number_factor(void) {
    assert_eval_eq("{Simplify[(a + 3 b)/(3 b)], Simplify[(a + 2 b)/(2 b)], Simplify[(a + b)/(2 b)], "
                   "Simplify[(a + 3 b)/b], Simplify[(4 (-1 + z^2))/(z (1 + z^2))]}",
                   "{1 + a/(3 b), 1 + a/(2 b), (a + b)/(2 b), 3 + a/b, (4 (-1 + z^2))/(z + z^3)}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_det_expanded);
    TEST(test_simplify_with_equations);
    TEST(test_simplify_number_factor);
    printf("All Det / equational Simplify tests passed.\n");
    return 0;
}
