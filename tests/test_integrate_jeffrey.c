/* test_integrate_jeffrey.c
 *
 * Tests for the Jeffrey-Rich continuous Weierstrass-substitution integrator
 * (Integrate`Weierstrass, Method -> "Weierstrass"; added 2026-06-09).
 *
 * Reference: D. J. Jeffrey and A. D. Rich, "The Evaluation of Trigonometric
 * Integrals Avoiding Spurious Discontinuities", ACM TOMS 20(1), 1994.
 *
 * Correctness is asserted by the universal predicate
 *   PossibleZeroQ[D[Integrate[f, x] /. Floor[_] -> 0, x] - f]
 * The "/. Floor[_] -> 0" strips the secular continuity term K Floor[(x-b)/p]
 * (whose derivative is 0 almost everywhere) before differentiating the smooth
 * antiderivative, so the check survives surface-form changes and the Floor term.
 * PossibleZeroQ (numeric two-phase sampler) is used rather than Simplify[.]===0
 * because the back-substituted Tan[x/2]-rational residues defeat Simplify on
 * some integrands (e.g. 1/(1 + Sin[x]^2), whose antiderivative carries
 * Tan[x/2]^3) and would otherwise time out.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Assert that the evaluated `input` is an unevaluated call with head `head`. */
static void assert_head_unevaluated(const char* input, const char* head) {
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol.name, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* Assert Simplify[D[Integrate[f,x] /. Floor[_]->0, x] - f] === 0 for the given
 * integrand string `f` (var x), routed through the Automatic cascade. */
static void assert_integral_ok(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[D[Integrate[%s, x] /. Floor[_] -> 0, x] - (%s)]", f, f);
    assert_eval_eq(buf, "True", 0);
}

/* Same, but forcing Method -> "Weierstrass". */
static void assert_integral_ok_forced(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[D[Integrate[%s, x, Method -> \"Weierstrass\"]"
        " /. Floor[_] -> 0, x] - (%s)]", f, f);
    assert_eval_eq(buf, "True", 0);
}

/* Assert the integral carries a (continuity) Floor term, i.e. its value is not
 * free of Floor.  Confirms the secular correction was actually emitted. */
static void assert_has_floor(const char* input, bool expected) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "FreeQ[%s, Floor]", input);
    assert_eval_eq(buf, expected ? "False" : "True", 0);
}

/* The Automatic cascade intercepts rational trig integrands with a kernel in a
 * denominator and returns a continuous antiderivative. */
static void test_automatic_trig(void) {
    /* Paper eq. (1)/(10): 3/(5 - 4 Cos[x]). */
    assert_integral_ok("3/(5 - 4 Cos[x])");
    /* Classic spurious-discontinuity example 1/(2 + Cos[x]). */
    assert_integral_ok("1/(2 + Cos[x])");
    /* Sin in the denominator. */
    assert_integral_ok("1/(5 + 4 Sin[x])");
    /* Both Sin and Cos. */
    assert_integral_ok("1/(2 + Sin[x] + Cos[x])");
    /* Even power of a kernel in the denominator. */
    assert_integral_ok("1/(1 + Sin[x]^2)");
    /* Numerator that is a derivative of the denominator (clean log). */
    assert_integral_ok("Sin[x]/(2 + Cos[x])");
}

/* Continuity correction: the headline examples must carry the Floor secular
 * term (the whole point of the Jeffrey-Rich algorithm). */
static void test_continuity_correction(void) {
    assert_has_floor("Integrate[3/(5 - 4 Cos[x]), x]", true);
    assert_has_floor("Integrate[1/(2 + Cos[x]), x]", true);
}

/* Hyperbolic generalisation via u = Tanh[x/2]: tanh(x/2) is a monotone
 * bijection, so the result is continuous and carries NO Floor term. */
static void test_hyperbolic(void) {
    /* Previously unevaluated by the cascade. */
    assert_integral_ok("1/(2 + Cosh[x])");
    assert_integral_ok("1/(3 + Sinh[x])");
    assert_integral_ok("1/(5 + 3 Cosh[x])");
    /* No spurious discontinuity -> no Floor term. */
    assert_has_floor("Integrate[1/(2 + Cosh[x]), x]", false);
}

/* Method plumbing: the option string and the explicit package head both route
 * to the routine and close the integrand. */
static void test_method_plumbing(void) {
    assert_integral_ok_forced("1/(5 + 4 Sin[x])");
    assert_eval_eq(
        "Simplify[D[Integrate`Weierstrass[1/(2 + Cos[x]), x] /. Floor[_] -> 0, x]"
        " - 1/(2 + Cos[x])]", "0", 0);
    /* Explicit method handles even polynomial trig (no denominator gate). */
    assert_eval_eq(
        "Simplify[D[Integrate`Weierstrass[Sin[x], x] /. Floor[_] -> 0, x]"
        " - Sin[x]]", "0", 0);
    /* Multiple-angle arguments are reduced by the TrigExpand pre-pass:
     * Sin[2 x], Cosh[x] Cosh[2 x] (= Cosh[x]^3 + Cosh[x] Sinh[x]^2). */
    assert_integral_ok_forced("Sin[2 x]/(2 + Cos[x])");
    assert_integral_ok_forced("Cosh[x] Cosh[2 x]");
}

/* Strict: the forced method/package head returns unevaluated when the integrand
 * is not a rational function of the trig/hyperbolic kernels of x. */
static void test_strict_no_match(void) {
    /* Polynomial in x outside any kernel. */
    assert_head_unevaluated(
        "Integrate[x + Sin[x], x, Method -> \"Weierstrass\"]", "Integrate");
    /* Kernel argument is a nonlinear function of x (TrigExpand cannot reduce it
     * to a kernel of the bare variable). */
    assert_head_unevaluated(
        "Integrate`Weierstrass[Sin[x^2], x]", "Integrate`Weierstrass");
    /* A non-trig transcendental factor of x. */
    assert_head_unevaluated(
        "Integrate`Weierstrass[Exp[x]/(2 + Cos[x]), x]", "Integrate`Weierstrass");
    /* Mixed trig + hyperbolic: not handled. */
    assert_head_unevaluated(
        "Integrate`Weierstrass[Sin[x] + Cosh[x], x]", "Integrate`Weierstrass");
    /* No trig/hyperbolic kernel at all. */
    assert_head_unevaluated(
        "Integrate`Weierstrass[x^2, x]", "Integrate`Weierstrass");
    /* A radical of a kernel is not a rational function of it. */
    assert_head_unevaluated(
        "Integrate`Weierstrass[Sqrt[1 + Cos[x]], x]", "Integrate`Weierstrass");
}

/* The Automatic cascade must NOT intercept polynomial trig: those keep their
 * clean table/Risch forms.  Integrate[Sin[x], x] -> -Cos[x], free of Floor. */
static void test_polynomial_trig_left_clean(void) {
    assert_eval_eq("Integrate[Sin[x], x]", "-Cos[x]", 0);
    assert_has_floor("Integrate[Sin[x], x]", false);
    assert_has_floor("Integrate[Sin[x] Cos[x], x]", false);
}

void test_integrate_jeffrey(void) {
    symtab_init();
    core_init();

    TEST(test_automatic_trig);
    TEST(test_continuity_correction);
    TEST(test_hyperbolic);
    TEST(test_method_plumbing);
    TEST(test_strict_no_match);
    TEST(test_polynomial_trig_left_clean);

    printf("All Integrate Weierstrass tests passed!\n");
}

int main(void) {
    test_integrate_jeffrey();
    return 0;
}
