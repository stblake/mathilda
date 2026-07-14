/* test_rt_resultant.c — Integrate`RothsteinTragerResultant[num, den, z, x].
 *
 * Cherry-substrate refactor R3: the named parametric Rothstein-Trager resultant
 *   Resultant[num - z D[den, x], den, x]
 * whose roots in z are the residues of num/den.  It is the argument generator for
 * the log / exponential-integral (Ei) parts of Cherry's special-function
 * integration (CHERRY_DESIGN.md §3.3), built on the existing `Resultant` builtin.
 *
 * We check the residue-root property against known residues (up to a nonzero
 * constant factor, so we compare Factor / roots rather than the raw scaling).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

static void test_residue_roots(void) {
    /* 1/x: single residue 1 -> root z = 1. */
    run_test("Integrate`RothsteinTragerResultant[1, x, z, x]", "Plus[1, Times[-1, z]]");
    /* 1/(x^2+1): residues ±I/2 -> 1 + 4 z^2 (roots ±I/2). */
    run_test("Factor[Integrate`RothsteinTragerResultant[1, x^2 + 1, z, x]]",
             "Plus[1, Times[4, Power[z, 2]]]");
    /* 2x/(x^2-2) = D[x^2-2]/(x^2-2): a logarithmic derivative, residue 1 at both
     * roots -> resultant has a double root at z = 1. */
    run_test("Factor[Integrate`RothsteinTragerResultant[2 x, x^2 - 2, z, x]]",
             "Times[-8, Power[Plus[-1, z], 2]]");
    /* 1/(x^2-1) = 1/2 (1/(x-1) - 1/(x+1)): residues ±1/2 -> roots ±1/2.  4z^2-1
     * factors over Q (unlike 4z^2+1), so Factor fully splits it (resultant is
     * defined up to a nonzero constant factor; the roots ±1/2 are the residues). */
    run_test("Factor[Integrate`RothsteinTragerResultant[1, x^2 - 1, z, x]]",
             "Times[-1, Plus[-1, Times[2, z]], Plus[1, Times[2, z]]]");
}

static void test_arity(void) {
    run_test("Head[Integrate`RothsteinTragerResultant[1, x, z]] "
             "=== Integrate`RothsteinTragerResultant", "True");
}

int main(void) {
    core_init();
    TEST(test_residue_roots);
    TEST(test_arity);
    printf("All RothsteinTragerResultant tests passed.\n");
    return 0;
}
