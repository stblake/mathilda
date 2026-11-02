/* test_issues_120_121.c
 *
 * GitHub #121: an inexact number inside an opaque kernel (f[1.], x[0.5]) is
 * part of an indeterminate, not a coefficient, so Factor (and the other
 * polynomial-algebra heads that round-trip inexact input through exact
 * arithmetic) must leave the polynomial's exact coefficients exact.
 * GitHub #120: CoefficientArrays with no variables answers the degree-0
 * element alone.
 * The === comparisons are the reporter's, with the reference's True.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_factor_kernel_inexact(void) {
    /* #121, verbatim */
    assert_eval_eq("Factor[f[1.] x[1] + f[1.0]] === f[1.0] (1 + x[1])", "True", 0);
    assert_eval_eq("Factor[f[2.] (-x[2] + x[3]) x[4] + f[2.]] === "
                   "-(f[2.]*(-1 + x[2]*x[4] - x[3]*x[4]))", "True", 0);
    /* the same for the siblings sharing the round trip: no machine number
     * appears outside the kernel */
    assert_eval_eq("FreeQ[FactorSquareFree[f[1.] x^2 + 2 f[1.] x + f[1.]] /. f[_] -> f, _Real]",
                   "True", 0);
    assert_eval_eq("FreeQ[Apart[f[1.]/(x (x + 1)), x] /. f[_] -> f, _Real]", "True", 0);
    assert_eval_eq("FreeQ[Cancel[(f[1.] x + f[1.])/(x + 1)] /. f[_] -> f, _Real]", "True", 0);
    /* a real coefficient still round-trips inexact */
    assert_eval_eq("FreeQ[Factor[2. x + 2.], _Real]", "False", 0);
    assert_eval_eq("FreeQ[Factor[Sin[2.] x + Sin[2.]], _Real]", "False", 0);
}

static void test_coefficientarrays_no_vars(void) {
    /* #120, verbatim */
    assert_eval_eq("CoefficientArrays[0, {}] === {0}", "True", 0);
    assert_eval_eq("CoefficientArrays[a + b, {}] === {a + b}", "True", 0);
    assert_eval_eq("CoefficientArrays[5]", "{5}", 0);
    assert_eval_eq("Normal[CoefficientArrays[{1, a}, {}]]", "{{1, a}}", 0);
    assert_eval_eq("Head[First[CoefficientArrays[{1, a}, {}]]]", "SparseArray", 0);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    TEST(test_factor_kernel_inexact);
    TEST(test_coefficientarrays_no_vars);
    printf("All issue #120/#121 tests passed.\n");
    return 0;
}
