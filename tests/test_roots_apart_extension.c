/* test_roots_apart_extension.c
 *
 * Solve writes roots of unity and symbolic linear roots as the reference
 * does, and Apart decomposes over radicals, in a compound variable, and
 * leaves an expression free of its variable alone (FeynCalc's Factor3 and
 * ToFCPartialFractionForm, 7 tests).  Expected values computed by the
 * reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_roots_of_unity(void) {
    assert_eval_eq("Solve[z^2 + z + 1 == 0, z]", "{{z -> -(-1)^(1/3)}, {z -> (-1)^(2/3)}}", 0);
    assert_eval_eq("Solve[z^4 + z^2 + 1 == 0, z, Complexes]",
                   "{{z -> -(-1)^(1/3)}, {z -> (-1)^(1/3)}, {z -> -(-1)^(2/3)}, {z -> (-1)^(2/3)}}", 0);
    assert_eval_eq("Solve[z^4 + 1 == 0, z]",
                   "{{z -> -(-1)^(1/4)}, {z -> (-1)^(1/4)}, {z -> -(-1)^(3/4)}, {z -> (-1)^(3/4)}}", 0);
    assert_eval_eq("Solve[z^5 - 1 == 0, z]",
                   "{{z -> 1}, {z -> -(-1)^(1/5)}, {z -> (-1)^(2/5)}, {z -> -(-1)^(3/5)}, {z -> (-1)^(4/5)}}", 0);
    assert_eval_eq("Solve[z^2 + z + 1 == 0, z, Reals]", "{}", 0);
}

static void test_symbolic_linear_roots(void) {
    assert_eval_eq("Solve[(1 - 4 x) (1 + 3 y) == 0, x, Complexes]", "{{x -> 1/4}}", 0);
    assert_eval_eq("Solve[2 (1 + y) x - 1 == 0, x]", "{{x -> 1/(2 (1 + y))}}", 0);
    assert_eval_eq("Solve[(-12 + 16 eps) x + 2 - 2 eps == 0, x]",
                   "{{x -> (-1 + eps)/(2 (-3 + 4 eps))}}", 0);
}

static void test_apart_extension(void) {
    assert_eval_eq("Together[Apart[1/((z - Sqrt[2]) (z + Sqrt[2])), z] - 1/(z^2 - 2)]", "0", 0);
    assert_eval_eq("Length[Apart[1/((z - Sqrt[2]) (z + Sqrt[2])), z]]", "2", 0);
    assert_eval_eq("Apart[g[x[3]]/(x[3]^2 - 1), x[3]]",
                   "g[x[3]]/(2 (-1 + x[3])) - g[x[3]]/(2 (1 + x[3]))", 0);
    assert_eval_eq("Apart[(a + b)/(c d), y]", "(a + b)/(c d)", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_roots_of_unity);
    TEST(test_symbolic_linear_roots);
    TEST(test_apart_extension);
    printf("All roots / Apart extension tests passed.\n");
    return 0;
}
