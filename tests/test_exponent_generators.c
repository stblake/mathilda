/* test_exponent_generators.c
 *
 * Two fixes found on FeynCalc's FactorList2 self-check (FCFeynmanRegularize-
 * Divergence took 3-6 s; the reference 0.06 s).  Expected values computed by
 * the reference kernel.
 *
 *  1. Together and Cancel read b^(n + k t) as b^n (b^t)^k (expsplit.c, as
 *     Factor already did), so x^(-2 e) and x^(-1 - 2 e) share the generator
 *     x^(-e): Together[2 x^(-2 e) - 2 x x^(-1 - 2 e)] is 0, not a quotient
 *     over the "LCM" x^(1 + 4 e).
 *
 *  2. Simplify's prime-power rebase (64^x -> 2^(6 x)) tested an integer base
 *     by trial division up to its square root: 33M divisions per pass for a
 *     16-digit base, which the check substitutes as random primes.  It now
 *     takes exact integer roots.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <time.h>

static void test_together_exponent_generators(void) {
    assert_eval_eq("Together[2 x^(-2 e) - 2 x x^(-1 - 2 e)]", "0", 0);
    assert_eval_eq("Together[2^n + 2^(n + 1)]", "3 2^n", 0);
    assert_eval_eq("Together[x^e/(x^(e + 1) + x^e)]", "1/(1 + x)", 0);
    assert_eval_eq("Together[2 X^(-2 e) - (2 X + 2 Y) X^(-1 - 2 e)]", "-2 X^(-1 - 2 e) Y", 0);
    assert_eval_eq("Simplify[2 (5 + 7 s)^(-2 e) - (10 + 14 s) (5 + 7 s)^(-1 - 2 e)]", "0", 0);
}

static void test_prime_rebase(void) {
    assert_eval_eq("Simplify[64^x 2^y]", "2^(6 x + y)", 0);
    assert_eval_eq("Simplify[49^x 7^y]", "7^(2 x + y)", 0);
    assert_eval_eq("Simplify[12^x]", "12^x", 0);
    clock_t t0 = clock();
    assert_eval_eq("Simplify[(4445811199353361^(3 e) - 4445811199353361^(3 e) (1 + q))/q]",
                   "-4445811199353361^(3 e)", 0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < 0.5, "Simplify with a 16-digit power base took %.2f s", secs);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_together_exponent_generators);
    TEST(test_prime_rebase);
    printf("All exponent-generator / prime-rebase tests passed.\n");
    return 0;
}
