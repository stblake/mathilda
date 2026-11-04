/* test_feyncalc_round14.c
 *
 * Kernel fixes found by FeynCalc's last failing unit tests (round 14, v0.6168+):
 *   - Factor of a product with a negative real coefficient and a sum with
 *     Gaussian coefficients decides the sign on the expansion, as the
 *     reference does: Factor[-2 (I x - y)] = 2 (-I x + y), the same as
 *     Factor[-2 I x + 2 y] (DiracTrace in the BMHV scheme left two equal
 *     traces that did not cancel).
 * Expected values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* lhs and rhs evaluate to the same expression */
static void same(const char* lhs, const char* rhs) {
    char buf[4096];
    snprintf(buf, sizeof buf, "(%s) === (%s)", lhs, rhs);
    assert_eval_eq(buf, "True", 0);
}

static void test_factor_gaussian_sign_on_expansion(void) {
    same("Factor[-2 (I x - y)]", "2*(-I*x + y)");
    same("Factor[-2 I x + 2 y]", "2*(-I*x + y)");
    same("Factor[-4 (I x - y) z]", "4*(-I*x + y)*z");
    same("Factor[-2 (I x - y) (a + b)]", "2*(a + b)*(-I*x + y)");
    same("Factor[-2 (I a b - c d) e]", "2*(-I*a*b + c*d)*e");
    /* the two BMHV traces of DiracTrace-ID46, factored, now cancel */
    same("Factor[-2 (I x - y)] - Factor[2 (-I x + y)]", "0");
    /* unchanged: a real or positive coefficient, a power of a sum */
    same("Factor[-2 (x - I y)]", "-2*(x - I*y)");
    same("Factor[(2 + 2 I) (x - y)]", "(2 + 2*I)*(x - y)");
    same("Factor[-I (x - y)]", "-I*(x - y)");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_factor_gaussian_sign_on_expansion);
    printf("All FeynCalc round-14 tests passed.\n");
    return 0;
}
