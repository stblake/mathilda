/* test_equal_symbolic_sides.c
 *
 * Equal/Unequal/Less numericalised a side that was not a numeric quantity:
 * N maps into the arguments of an inert call, so ff[5] == x evaluated
 * ff[5.] and fired a definition the reference never reaches.  FeynCalc's
 * GA[5.] guard printed GA::noint inside every `... == x` test (8 tests).
 * The numeric fallback now runs only on a NumericQ operand, and NumericQ
 * recognises an algebraic Root[f, k] (it was False).
 * Expected values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_no_numericalisation_of_symbolic_sides(void) {
    /* a definition on Real arguments records any leak */
    assert_eval_eq("leaks = 0; ff[r_Real] := (leaks++; 0); "
                   "{ff[5] == x, 2 ff[5] == x, ff[5] . gg[1] == x, Unequal[ff[5], x], ff[5] < 3, leaks}",
                   "{ff[5] == x, 2 ff[5] == x, ff[5] . gg[1] == x, ff[5] != x, ff[5] < 3, 0}", 0);
    assert_eval_eq("f[1/2] == 0.5", "f[1/2] == 0.5", 0);
}

static void test_numeric_sides_still_decided(void) {
    assert_eval_eq("Sqrt[2] == 1.4142135623730951", "True", 0);
    assert_eval_eq("Pi == 3.14159", "False", 0);
    assert_eval_eq("1/3 == 0.3333333333333333", "True", 0);
    assert_eval_eq("Sin[1] == Sin[1.]", "True", 0);
    assert_eval_eq("Sqrt[2] < 3/2", "True", 0);
}

static void test_algebraic_root_is_numeric(void) {
    /* NumericQ of an algebraic Root was False, so the gate above stopped
     * Equal from comparing one against a number (Minimize lost its Lagrange
     * candidates). */
    assert_eval_eq("NumericQ[Root[-1 + 6 #1 + 22 #1^2 + 16 #1^3 &, 2]]", "True", 0);
    assert_eval_eq("NumericQ[1 + Root[-1 + 6 #1 + 22 #1^2 + 16 #1^3 &, 2]^2]", "True", 0);
    assert_eval_eq("Root[-1 + 6 #1 + 22 #1^2 + 16 #1^3 &, 2] == 0.1", "False", 0);
    assert_eval_eq("NumericQ[Root[#1^2 - a &, 1]]", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_no_numericalisation_of_symbolic_sides);
    TEST(test_numeric_sides_still_decided);
    TEST(test_algebraic_root_is_numeric);
    printf("All Equal symbolic-side tests passed.\n");
    return 0;
}
