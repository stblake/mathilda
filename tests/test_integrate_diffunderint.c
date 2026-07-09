/* test_integrate_diffunderint.c
 *
 * Tests for definite integration by differentiation under the integral sign
 * (Integrate`DiffUnderInt, Method -> "DiffUnderInt" / "DifferentiationUnder-
 * Integral").  Correctness is asserted by Simplify[<integral> - <expected>,
 * <assumptions>] === 0 (numeric comparison would violate the project rule that
 * Integrate performs no NIntegrate crosscheck; here the check is symbolic and
 * lives in the test, not the method).
 *
 * The method computes the standard parameter-dependent families -- Laplace /
 * Fourier half-line, sinc / Frullani, and even-rational half-line -- via its own
 * closed-form evaluators, so these all return fast, clean closed forms.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Assert Simplify[<expr> - <expected>, <assum>] === 0. */
static void assert_closes(const char* integral, const char* expected,
                          const char* assum) {
    char buf[1024];
    if (assum && assum[0])
        snprintf(buf, sizeof(buf), "Simplify[(%s) - (%s), %s]", integral, expected, assum);
    else
        snprintf(buf, sizeof(buf), "Simplify[(%s) - (%s)]", integral, expected);
    assert_eval_eq(buf, "0", 0);
}

/* Assert |N[(<integral> - <expected>) /. subst]| < 1e-8.  Used where the result
 * is correct but Simplify cannot PROVE the two forms equal (e.g. it does not
 * know ArcTan[b/a] + ArcTan[a/b] == Pi/2, or combine a split Log).  This is a
 * test-side numeric check; the method itself performs no numeric integration. */
static void assert_closes_num(const char* integral, const char* expected,
                             const char* subst) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
             "Abs[N[((%s) - (%s)) /. {%s}]] < 1/100000000", integral, expected, subst);
    assert_eval_eq(buf, "True", 0);
}

/* Assert that the evaluated `input` is an unevaluated call with head `head`. */
static void assert_head_unevaluated(const char* input, const char* head) {
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* Log-in-the-denominator / power families over [0,1]. */
static void test_log_denominator(void) {
    assert_closes("Integrate[(x^a - 1)/Log[x], {x, 0, 1}, Method -> \"DiffUnderInt\"]",
                  "Log[1 + a]", "");
    assert_closes("Integrate[(x^a - x^b)/Log[x], {x, 0, 1}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> {a > -1, b > -1}]",
                  "Log[(1 + a)/(1 + b)]", "a > -1 && b > -1");
}

/* Laplace / Fourier half-line (exponential x trig). */
static void test_laplace_fourier(void) {
    assert_closes("Integrate[(Exp[-a x] - Exp[-b x])/x, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> {a > 0, b > 0}]",
                  "Log[b/a]", "a > 0 && b > 0");
    /* Correct value; Simplify cannot prove ArcTan[b/a]+ArcTan[a/b]==Pi/2. */
    assert_closes_num("Integrate[Exp[-a x] Sin[b x]/x, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> a > 0]",
                  "ArcTan[b/a]", "a -> 13/10, b -> 7/10");
    /* Correct value; Simplify cannot combine the split Log. */
    assert_closes_num("Integrate[(Exp[-a x] - Exp[-b x]) Cos[c x]/x, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", "
                  "Assumptions -> {a > 0, b > 0, Element[c, Reals]}]",
                  "1/2 Log[(b^2 + c^2)/(a^2 + c^2)]", "a -> 13/10, b -> 7/10, c -> 1/2");
    assert_closes("Integrate[Exp[-x] Sin[a x]^2/x, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\"]",
                  "1/4 Log[1 + 4 a^2]", "a > 0");
}

/* Sinc / Frullani (oscillatory /x, /x^2). */
static void test_sinc(void) {
    assert_closes("Integrate[Sin[a x]^2/x^2, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> a > 0]",
                  "Pi a/2", "a > 0");
    assert_closes("Integrate[(1 - Cos[a x])/x^2, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> a > 0]",
                  "Pi a/2", "a > 0");
    assert_closes("Integrate[(ArcTan[a x] - ArcTan[b x])/x, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> {a > 0, b > 0}]",
                  "Pi/2 Log[a/b]", "a > 0 && b > 0");
}

/* Even-rational half-line (with a Log/ArcTan outer parameter integral). */
static void test_rational_halfline(void) {
    assert_closes("Integrate[Log[1 + a^2 x^2]/(1 + x^2), {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> a > 0]",
                  "Pi Log[1 + a]", "a > 0");
    assert_closes("Integrate[ArcTan[a x]/(x (1 + x^2)), {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> a > 0]",
                  "Pi/2 Log[1 + a]", "a > 0");
    assert_closes("Integrate[Log[1 + a^2 x^2]/(x^2 (1 + x^2)), {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> a > 0]",
                  "Pi (a - Log[1 + a])", "a > 0");
}

/* Non-even (decaying) rational half-line: the differentiated inner integral is a
 * DECAYING sinc ∫₀^∞ e^{-c x} Sin[a x]/x dx = ArcTan[a/c], whose Laplace image
 * M(s) = a/((s+c)^2+a^2) is non-even -> rational_halfline_general, producing a
 * real ArcTan directly (no complex-Log reduction). */
static void test_rational_halfline_general(void) {
    assert_closes("Integrate[Exp[-c x] (1 - Cos[a x])/x^2, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\", Assumptions -> {a > 0, c > 0}]",
                  "a ArcTan[a/c] - c/2 Log[1 + a^2/c^2]", "a > 0 && c > 0");
}

/* Gaussian moment family + Erf back-integration: differentiating removes the 1/x
 * to leave the cosine moment ∫₀^∞ e^{-x^2} Cos[a x] dx = (Sqrt[Pi]/2) e^{-a^2/4},
 * whose parameter integral is an Erf the engine cannot produce. */
static void test_gaussian(void) {
    assert_closes("Integrate[Exp[-x^2] Sin[a x]/x, {x, 0, Infinity}, "
                  "Method -> \"DiffUnderInt\"]",
                  "Pi/2 Erf[a/2]", "");
}

/* Reachable three ways: Method string, the alias, and the direct builtin. */
static void test_routing(void) {
    assert_closes("Integrate[(x^a - 1)/Log[x], {x, 0, 1}, "
                  "Method -> \"DifferentiationUnderIntegral\"]",
                  "Log[1 + a]", "");
    assert_closes("Integrate`DiffUnderInt[(x^a - 1)/Log[x], {x, 0, 1}]",
                  "Log[1 + a]", "");
    /* Also picked up by the Automatic cascade (an integral only this method
     * closes): no Method option given. */
    assert_closes("Integrate[Sin[a x]^2/x^2, {x, 0, Infinity}, "
                  "Assumptions -> a > 0]",
                  "Pi a/2", "a > 0");
}

/* Safety: integrands the method cannot yet close return the input unevaluated
 * *fast* (never hang) rather than a wrong value. */
static void test_declines_cleanly(void) {
    /* Gaussian cosine moment: DiffUnderInt differentiates the parameter into a
     * Sin-Gaussian (a Dawson/Erfi form the moment family declines), so no
     * parameter closes and it returns unevaluated fast -- this integrand is
     * directly integrable by a Gaussian family, not a Feynman target. */
    assert_head_unevaluated(
        "Integrate[Exp[-x^2] Cos[2 a x], {x, 0, Infinity}, "
        "Method -> \"DiffUnderInt\"]", "Integrate");
    /* No free parameter -> not applicable. */
    assert_head_unevaluated(
        "Integrate`DiffUnderInt[Sin[x]/x, {x, 0, Infinity}]", "Integrate`DiffUnderInt");
    /* Whole-line integrand with a real axis pole (Exp[x]-1 = 0 at x = 0): the
     * integral diverges (principal value only), so DiffUnderInt declines up front
     * rather than pursue the non-terminating Integrate[x^k f] escalation that
     * used to run for minutes.  Must return fast + unevaluated. */
    assert_head_unevaluated(
        "Integrate`DiffUnderInt[Exp[a x]/(Exp[x]-1), {x, -Infinity, Infinity}, "
        "Assumptions -> 0 < a < 1]", "Integrate`DiffUnderInt");
}

void test_integrate_diffunderint(void) {
    symtab_init();
    core_init();

    TEST(test_log_denominator);
    TEST(test_laplace_fourier);
    TEST(test_sinc);
    TEST(test_rational_halfline);
    TEST(test_rational_halfline_general);
    TEST(test_gaussian);
    TEST(test_routing);
    TEST(test_declines_cleanly);

    printf("All Integrate DiffUnderInt tests passed!\n");
}

int main(void) {
    test_integrate_diffunderint();
    return 0;
}
