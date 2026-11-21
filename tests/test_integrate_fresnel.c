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

/* As assert_fresnel, in the variable `v`, verified by PossibleZeroQ (the
 * verifier; an exact Simplify need not reduce the Sqrt[v] chain-rule forms),
 * and within a time bound that the ParallelMixedSpecial fall-through (~11 s
 * per call, v0.6661) cannot meet. */
static void assert_fresnel_in(const char* f, const char* v) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{tr = AbsoluteTiming[Integrate[%s, %s]]}, "
        " Head[tr[[2]]] =!= Integrate && (!FreeQ[tr[[2]], FresnelS] || !FreeQ[tr[[2]], FresnelC]) "
        " && PossibleZeroQ[D[tr[[2]], %s] - (%s)] && tr[[1]] < 1]", f, v, v, f);
    ASSERT_MSG(eval_is(buf, "True"), "Integrate[%s, %s]: expected a fast, verified Fresnel form", f, v);
}

/* v0.6662: the kernel's answer template spelled the variable as the literal
 * symbol x, so in any other variable the candidate failed its own diff-back
 * and the stage declined (Integrate[Cos[t^2], t] was ::nonelem).  And the
 * power family -- x^(half-integer) Trig[a x + c], x^(2m) Trig[a x^2 + c] --
 * reduces to the kernel by parts. */
static void test_fresnel_other_variables_and_powers(void) {
    assert_fresnel_in("Cos[t^2]", "t");
    assert_fresnel_in("Sin[3 u^2 + u - 2]", "u");
    assert_fresnel_in("Sqrt[t] Sin[t]", "t");               /* DSolve 2.2.23-2299's quadrature */
    assert_fresnel_in("t^(3/2) Cos[3 t]", "t");
    assert_fresnel_in("Sin[t]/Sqrt[t]", "t");
    assert_fresnel_in("Cos[2 t + 1]/t^(3/2)", "t");
    assert_fresnel_in("5 t^(5/2) Sin[t/2 - 3]", "t");
    assert_fresnel_in("2 t^2 Sin[t^2]", "t");
    assert_fresnel_in("t^4 Cos[3 t^2 + 1]", "t");
    assert_fresnel_in("Sin[t^2]/t^2", "t");
    assert_fresnel_in("x^(5/2) Cos[x]", "x");
}

/* Non-Fresnel integrands must NOT be mis-handled: linear/cubic phase and a
 * genuinely elementary sibling. */
static void test_fresnel_declines(void) {
    /* Sin[x] is elementary (-Cos[x]); Sin[x^3] has no elementary/Fresnel form. */
    ASSERT_MSG(eval_is("FreeQ[Integrate[Sin[x], x], FresnelS]", "True"),
        "Sin[x] must not produce Fresnel");
    /* Sin[x^3] is the GammaPower stage's (v0.6588: incomplete Gamma), never Fresnel */
    ASSERT_MSG(eval_is("With[{r = Integrate[Sin[x^3], x]}, FreeQ[r, FresnelS | FresnelC] && !FreeQ[r, Gamma]]", "True"),
        "Sin[x^3] must not produce Fresnel");
}

int main(void) {
    core_init();
    TEST(test_fresnel);
    TEST(test_fresnel_declines);
    TEST(test_fresnel_other_variables_and_powers);
    printf("All Fresnel integration tests passed.\n");
    return 0;
}
