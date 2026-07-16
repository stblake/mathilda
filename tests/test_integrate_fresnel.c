/* test_integrate_fresnel.c — Fresnel integrals of a Gaussian-phase trig integrand.
 *
 * K Sin[a x^2 + b x + c] / K Cos[...] -> FresnelS/FresnelC by completing the
 * square (integrate_fresnel.c).  Each answer is diff-back verified.
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
#include <stdbool.h>

static bool eval_is(const char* s, const char* expected) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    char* got = expr_to_string_fullform(r);
    bool ok = strcmp(got, expected) == 0;
    if (!ok) printf("  [%s] -> %s (expected %s)\n", s, got, expected);
    free(got); expr_free(r); expr_free(e);
    return ok;
}

/* Assert Integrate[f] closes to a Fresnel form and diff-backs exactly to f. */
static void assert_fresnel(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x]},"
        " Head[r] =!= Integrate && (!FreeQ[r, FresnelS] || !FreeQ[r, FresnelC])]", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a Fresnel form", f);
    snprintf(buf, sizeof(buf), "Simplify[D[Integrate[%s, x], x] - (%s)]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: diff-back nonzero", f);
}

static void test_fresnel(void) {
    /* pure quadratic phase */
    assert_fresnel("Sin[x^2]");
    assert_fresnel("Cos[x^2]");
    assert_fresnel("Sin[2 x^2]");
    assert_fresnel("Cos[3 x^2]");
    assert_fresnel("Sin[x^2/2]");
    /* constant prefactor */
    assert_fresnel("5 Sin[x^2]");
    assert_fresnel("Cos[3 x^2 + 1]");
    /* linear term -> complete the square */
    assert_fresnel("Sin[x^2 + x + 1]");
    assert_fresnel("Cos[x^2 - x]");
    assert_fresnel("Sin[2 x^2 + 4 x + 3]");
}

/* Non-Fresnel integrands must NOT be mis-handled: linear/cubic phase and a
 * genuinely elementary sibling. */
static void test_fresnel_declines(void) {
    /* Sin[x] is elementary (-Cos[x]); Sin[x^3] has no elementary/Fresnel form. */
    ASSERT_MSG(eval_is("FreeQ[Integrate[Sin[x], x], FresnelS]", "True"),
        "Sin[x] must not produce Fresnel");
    ASSERT_MSG(eval_is("Head[Integrate[Sin[x^3], x]] === Integrate", "True"),
        "Sin[x^3] must stay unevaluated");
}

int main(void) {
    core_init();
    TEST(test_fresnel);
    TEST(test_fresnel_declines);
    printf("All Fresnel integration tests passed.\n");
    return 0;
}
