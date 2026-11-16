/* test_integrate_gammapower.c — the incomplete-Gamma stage (integrate_gammapower.c).
 *
 * Two recognizers:
 *  - K x^p E^(a x^m), p SYMBOLIC, m a positive integer (the original stage);
 *  - K x^p T(c x^q), T in {Exp, Sin, Cos, Sinh, Cosh}, q an exact nonzero
 *    rational, s = (p+1)/q with 2s not an integer (v0.6588), through the
 *    branch-exact G(a) = -(1/q) x^(p+1) w^(-s) Gamma[s, w], w = -a x^q.
 *    Integrate[Sin[x^(3/2)], x] reached only the Meijer G tail before (3.6 s,
 *    a HypergeometricPFQ, close to that stage's 5 s budget); FriCAS answers it
 *    in Gamma functions.
 * Every answer is checked by differentiation at real points on both sides of 0
 * and at a complex point: the form needs no branch convention, so the residual
 * is zero everywhere, not only on x > 0.
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

/* Integrate[f, x] closes in under `secs` seconds to a Gamma form whose
 * derivative is f at x = 3/7, 2, -5/3 and 1 + I (to 10^-20 at 30 digits),
 * with any parameters at the given rules. */
static void assert_gamma(const char* f, const char* params, int secs) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Quiet[With[{t = AbsoluteTiming[Integrate[%s, x]]}, With[{F = t[[2]]},"
        " t[[1]] < %d && FreeQ[F, Integrate | HypergeometricPFQ] && !FreeQ[F, Gamma] &&"
        " And @@ (Abs[N[(D[F, x] - (%s)) /. %s /. x -> #, 30]] < 10^-20 & /@ {3/7, 2, -5/3, 1 + I})]]]",
        f, secs, f, params);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a verified Gamma form", f);
}

static void test_trig_power(void) {
    /* the reported case, and its family */
    assert_gamma("Sin[x^(3/2)]", "{}", 2);
    assert_gamma("Cos[x^(3/2)]", "{}", 2);
    assert_gamma("x Sin[x^(3/2)]", "{}", 2);
    assert_gamma("5 Sqrt[x] Cos[x^(5/2)]", "{}", 2);
    assert_gamma("Sin[x^3]", "{}", 2);           /* was a 3.6 s HypergeometricPFQ */
    assert_gamma("Sin[x]/x^(1/3)", "{}", 2);
    /* hyperbolic and exponential kernels */
    assert_gamma("Sinh[x^(3/2)]", "{}", 2);
    assert_gamma("Cosh[3 x^(4/3)]", "{}", 2);
    assert_gamma("E^(x^(3/2))", "{}", 2);
    assert_gamma("E^(-2 x^(5/3))", "{}", 2);
    /* symbolic coefficient and symbolic power */
    assert_gamma("Sin[a x^(3/2)]", "{a -> 7/5}", 3);
    assert_gamma("x^n Sin[x^3]", "{n -> 2/7}", 3);
}

/* 2s an integer: the elementary / Fresnel / Erf / exponential-integral forms of
 * the later stages are kept (the stage declines). */
static void test_declines(void) {
    ASSERT_MSG(eval_is("Integrate[x^2 Sin[x^3], x]", "Times[Rational[-1, 3], Cos[Power[x, 3]]]"),
               "x^2 Sin[x^3] must stay elementary");
    ASSERT_MSG(eval_is("FreeQ[Integrate[Sin[x^2], x], Gamma] && !FreeQ[Integrate[Sin[x^2], x], FresnelS]", "True"),
               "Sin[x^2] must stay a Fresnel form");
    ASSERT_MSG(eval_is("FreeQ[Integrate[Sin[Sqrt[x]], x], Gamma | Integrate]", "True"),
               "Sin[Sqrt[x]] must stay elementary");
    ASSERT_MSG(eval_is("Integrate[E^(-x^2), x]", "Times[Rational[1, 2], Power[Pi, Rational[1, 2]], Erf[x]]"),
               "E^(-x^2) must stay an Erf form");
    /* s = 0: an exponential-integral case (SinIntegral in the full kernel) */
    ASSERT_MSG(eval_is("FreeQ[Quiet[Integrate[Sin[x]/x, x]], Gamma]", "True"), "Sin[x]/x is not this stage's");
    ASSERT_MSG(eval_is("FreeQ[Integrate[E^x/Sqrt[x], x], Gamma]", "True"), "E^x/Sqrt[x] must stay Erfi");
    /* two kernels, or a non-monomial argument: not this stage's shape */
    ASSERT_MSG(eval_is("FreeQ[Integrate[Sin[x^(3/2) + x], x], Gamma[_, _]]", "True"),
               "Sin[x^(3/2) + x] is not a monomial argument");
    /* the symbolic-exponent recognizer is unchanged */
    ASSERT_MSG(eval_is("Integrate[x^n E^(-x), x]", "Times[-1, Gamma[Plus[1, n], x]]"),
               "x^n E^(-x) keeps its form");
}

int main(void) {
    core_init();
    TEST(test_trig_power);
    TEST(test_declines);
    printf("All GammaPower integration tests passed.\n");
    return 0;
}
