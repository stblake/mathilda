/* test_diracdelta_normal.c
 *
 * DiracDelta's sign and scale normal form (src/internal/distributions.m), as
 * the reference writes it: an argument whose terms are mostly negated has the
 * sign taken out (FeynCalc's FCFeynmanParameterJoin writes
 * DiracDelta[-1 + x1 + x2]), and a numeric factor comes out as 1/|c|.
 * Expected values computed by the reference kernel (InputForm strings).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void eq_if(const char* e, const char* expected) {
    char buf[2048], quoted[2048];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);   /* a String prints quoted */
    assert_eval_eq(buf, quoted, 0);
}

static void test_diracdelta_normal_form(void) {
    /* the reference's sign and scale normalisation (FeynCalc's
     * FCFeynmanParameterJoin writes DiracDelta[-1 + x1 + x2]) */
    eq_if("{DiracDelta[1 - x1 - x2], DiracDelta[-x], DiracDelta[1 - x], DiracDelta[1 - x + y], "
          "DiracDelta[-x + 3 y], DiracDelta[-a + b], DiracDelta[a - b], DiracDelta[1 - 2 x], "
          "DiracDelta[-Sin[x]], DiracDelta[1.5 - x], DiracDelta[1 - x, -y], HeavisideTheta[-x]}",
          "{DiracDelta[-1 + x1 + x2], DiracDelta[x], DiracDelta[-1 + x], DiracDelta[1 - x + y], "
          "DiracDelta[x - 3*y], DiracDelta[-a + b], DiracDelta[a - b], DiracDelta[1 - 2*x], "
          "DiracDelta[Sin[x]], DiracDelta[1.5 - x], DiracDelta[-1 + x, y], HeavisideTheta[-x]}");
    eq_if("{DiracDelta[-2 x], DiracDelta[x/2], DiracDelta[-x/2 + 1], DiracDelta[4 x + 6 y], "
          "DiracDelta[1 + 2 x], DiracDelta[2.5 x], DiracDelta[2 x, y], DiracDelta[Sqrt[2] x], "
          "DiracDelta[2.5 x + 2.5 y]}",
          "{DiracDelta[x]/2, 2*DiracDelta[x], 2*DiracDelta[-2 + x], DiracDelta[2*x + 3*y]/2, "
          "DiracDelta[1 + 2*x], 0.4*DiracDelta[x], DiracDelta[x, y]/2, DiracDelta[x]/Sqrt[2], "
          "DiracDelta[2.5*x + 2.5*y]}");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();          /* the rules live in src/internal/distributions.m */
    TEST(test_diracdelta_normal_form);
    printf("All DiracDelta normal-form tests passed.\n");
    return 0;
}
