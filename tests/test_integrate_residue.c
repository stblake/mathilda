/* test_integrate_residue.c
 *
 * Tests for definite integration by the residue theorem
 * (src/calculus/integrate_residue.c):
 *   Integrate[f, {x, a, b}]                       (auto-dispatch, before Newton-Leibniz)
 *   Integrate[f, {x, a, b}, Method -> "Residue"]  (strict: no NL fallback)
 *   Integrate`ContourResidue[f, {x, a, b}]        (explicit entry point)
 *
 * Coverage: rational integrands on (-Inf,Inf) including higher-order poles
 * (Family A), Fourier/Jordan integrands with Cos/Sin/Exp kernels (Family B),
 * rational-in-{Sin,Cos} integrands over a full period (Family C), the
 * principal-value half-residue (Sin[x]/x = Pi), the even half-line, dispatch
 * selection, and the negative controls that MUST fall through / stay unevaluated
 * (real-axis pole, branch point, partial period).
 *
 * Closed forms that Mathilda leaves in an equivalent-but-unsimplified surface
 * form are pinned numerically via Chop[N[value - reference]] == 0.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------------------------------------------------------------------------
 * Family A -- rational integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */
static void test_family_rational(void) {
    check_eq("Integrate[1/(1+x^2), {x, -Infinity, Infinity}]", "Pi");
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}]", "Pi/Sqrt[2]");
    check_eq("Integrate[x^2/(1+x^4), {x, -Infinity, Infinity}]", "Pi/Sqrt[2]");
    check_eq("Integrate[1/((x^2+1)(x^2+4)), {x, -Infinity, Infinity}]", "1/6 Pi");
    /* Order-2 pole at I. */
    check_eq("Integrate[1/(1+x^2)^2, {x, -Infinity, Infinity}]", "1/2 Pi");
    /* Value-form independent numeric confirmation. */
    check_eq("Chop[N[Integrate[1/(1+x^6), {x, -Infinity, Infinity}] - 2 Pi/3]]", "0");
}

/* -------------------------------------------------------------------------
 * Family B -- Fourier / Jordan integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */
static void test_family_fourier(void) {
    check_eq("Integrate[Cos[x]/(1+x^2), {x, -Infinity, Infinity}]", "Pi/E");
    check_eq("Integrate[x Sin[x]/(1+x^2), {x, -Infinity, Infinity}]", "Pi/E");
    /* Cos[a x]/(x^2+b^2) = (Pi/b) e^{-a b}. */
    check_eq("Integrate[Cos[2 x]/(x^2+9), {x, -Infinity, Infinity}]", "(1/3 Pi)/E^6");
    /* x Sin[x]/(x^2+4) = Pi/E^2. */
    check_eq("Chop[N[Integrate[x Sin[x]/(x^2+4), {x, -Infinity, Infinity}] - Pi/E^2]]", "0");
    /* Higher-order (double) complex poles: Cos[a x]/(x^2+1)^2 = (Pi/2)(1+a) e^{-a}.
     * Regression for the Laurent-pad bug that under-expanded the Taylor factor at
     * a complex pole (I = Complex[0,1]) and dropped the residue's product-rule
     * cross term -- gave (Pi/2) a e^{-a} instead of (Pi/2)(1+a) e^{-a}. */
    check_eq("Integrate[Cos[x]/(x^2+1)^2, {x, -Infinity, Infinity}]", "Pi/E");
    check_eq("Integrate[Cos[3 x]/(x^2+1)^2, {x, -Infinity, Infinity}]", "(2 Pi)/E^3");
    check_eq("Chop[N[Integrate[Cos[2 x]/(x^2+1)^2, {x, -Infinity, Infinity}] "
             "- (Pi/2)(1+2)/E^2]]", "0");
    /* Triple pole (order-3): Cos[x]/(x^2+1)^3 = 7 Pi/(8 E) (matches NIntegrate). */
    check_eq("Integrate[Cos[x]/(x^2+1)^3, {x, -Infinity, Infinity}]", "(7/8 Pi)/E");

    /* Bare complex-exponential kernel Exp[I k x], which the evaluator normalises
     * to Power[E, .] rather than an Exp[.] head. Both spellings must match. */
    check_eq("Integrate[Exp[I x]/(x^2+1), {x, -Infinity, Infinity}]", "Pi/E");
    check_eq("Integrate[E^(I x)/(x^2+1), {x, -Infinity, Infinity}]", "Pi/E");
    check_eq("Integrate[Exp[2 I x]/(x^2+1), {x, -Infinity, Infinity}]", "Pi/E^2");
    /* Lower-half-plane closure (negative frequency). */
    check_eq("Integrate[Exp[-I x]/(x^2+1), {x, -Infinity, Infinity}]", "Pi/E");
    /* Double pole through the Exp kernel path. */
    check_eq("Integrate[Exp[I x]/(x^2+1)^2, {x, -Infinity, Infinity}]", "Pi/E");
    /* Complex-valued result: Re part is odd (0), Im part = Pi/E. */
    check_eq("Integrate[x Exp[I x]/(x^2+1), {x, -Infinity, Infinity}]", "(I Pi)/E");
}

/* -------------------------------------------------------------------------
 * Family C -- rational-in-{Sin,Cos} over a full period.
 * ---------------------------------------------------------------------- */
static void test_family_trig(void) {
    check_eq("Integrate[1/(2+Cos[x]), {x, 0, 2 Pi}]", "(2 Pi)/Sqrt[3]");
    check_eq("Integrate[1/(5-4 Cos[x]), {x, 0, 2 Pi}]", "2/3 Pi");
    /* (-Pi, Pi) is also a full period. */
    check_eq("Chop[N[Integrate[1/(2+Cos[x]), {x, -Pi, Pi}] - 2 Pi/Sqrt[3]]]", "0");
}

/* -------------------------------------------------------------------------
 * Principal value (simple real-axis pole, half residue) + even half-line.
 * ---------------------------------------------------------------------- */
static void test_principal_value(void) {
    /* Sin[x]/x is analytic at 0 (the kernel supplies a matching zero): the
     * ordinary integral converges to Pi via the half residue. */
    check_eq("Integrate[Sin[x]/x, {x, -Infinity, Infinity}]", "Pi");
    /* Cos[x]/x has a GENUINE pole at 0 (Cos[0] = 1 != 0): the ordinary integral
     * diverges, so the residue method must NOT return a value -- the method
     * leaves it for Newton-Leibniz (strict "Residue" stays unevaluated). */
    check_eq("Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> \"Residue\"]");
}

static void test_half_line(void) {
    check_eq("Integrate[1/(1+x^4), {x, 0, Infinity}]", "(1/2 Pi)/Sqrt[2]");
    check_eq("Chop[N[Integrate[1/(1+x^2), {x, 0, Infinity}] - Pi/2]]", "0");
}

/* -------------------------------------------------------------------------
 * Dispatch: explicit method / builtin; strict Residue has no NL fallback.
 * ---------------------------------------------------------------------- */
static void test_dispatch(void) {
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Pi/Sqrt[2]");
    check_eq("Integrate`ContourResidue[1/(2+Cos[x]), {x, 0, 2 Pi}]", "(2 Pi)/Sqrt[3]");
    /* Strict "Residue" on an integrand no family recognises: unevaluated, NO
     * Newton-Leibniz fallback. */
    check_eq("Integrate[Sin[x], {x, 0, Pi}, Method -> \"Residue\"]",
             "Integrate[Sin[x], {x, 0, Pi}, Method -> \"Residue\"]");
}

/* -------------------------------------------------------------------------
 * Negative controls + regression: ordinary definite integrals unaffected.
 * ---------------------------------------------------------------------- */
static void test_negative_controls(void) {
    /* Real-axis pole of a non-PV integrand: not a clean residue answer. */
    check_eq("Integrate[1/(1+x^3), {x, -Infinity, Infinity}]",
             "Integrate[1/(1 + x^3), {x, -Infinity, Infinity}]");
    /* Branch point (not rational): the residue method must not fire.  (Under
     * strict Method -> "Residue" it stays unevaluated; the auto-dispatch now
     * closes the even integrand via the symmetry reduction -> half-line, giving
     * Gamma[1/4]^2/(2 Sqrt[Pi]).) */
    check_eq("Integrate[1/Sqrt[1+x^4], {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[1/Sqrt[1 + x^4], {x, -Infinity, Infinity}, Method -> \"Residue\"]");
    /* Rational integrand with a real-axis pole (Family A): the ordinary integral
     * genuinely diverges (only a principal value would exist).  The residue
     * method now flags this conclusively, so the dispatcher emits Integrate::idiv
     * (to stderr) and leaves the integral unevaluated for both the explicit
     * "Residue" method and the Automatic path. */
    check_eq("Integrate[1/(x^2-1), {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[1/(-1 + x^2), {x, -Infinity, Infinity}, Method -> \"Residue\"]");
    check_eq("Integrate[1/(x^2-1), {x, -Infinity, Infinity}]",
             "Integrate[1/(-1 + x^2), {x, -Infinity, Infinity}]");
    /* Trig integrand whose denominator vanishes on the real axis (Family C):
     * the periodic integral diverges -> unevaluated. */
    check_eq("Integrate[1/(1-Cos[x]), {x, 0, 2 Pi}, Method -> \"Residue\"]",
             "Integrate[1/(1 - Cos[x]), {x, 0, 2 Pi}, Method -> \"Residue\"]");
}

static void test_regression_finite(void) {
    /* Finite definite integrals still go through Newton-Leibniz, unchanged. */
    check_eq("Integrate[1/x, {x, 1, 2}]", "Log[2]");
    check_eq("Integrate[Sin[x], {x, 0, Pi}]", "2");
    check_eq("Integrate[1/(1 + x^2), {x, 0, 1}]", "1/4 Pi");
    check_eq("Integrate[x^2, {x, 0, 1}]", "1/3");
}

/* =========================================================================
 * Assumptions-driven / new contour families.  A parametric closed form is
 * pinned numerically via a rational substitution + Chop[N[value - reference]]
 * (exact string matching a symbolic surface form is brittle); purely numeric
 * closed forms are matched exactly.
 * ====================================================================== */

/* -------------------------------------------------------------------------
 * Option plumbing: Integrate accepts `Assumptions -> ...` (no Integrate::method
 * error), in any order, alongside Method, on definite and indefinite forms.
 * ---------------------------------------------------------------------- */
static void test_assumptions_option(void) {
    /* Assumptions accepted on a numeric definite integral (option ignored, value
     * unchanged) -- and combined with Method, in either order. */
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}, Assumptions -> a > 0]",
             "Pi/Sqrt[2]");
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}, "
             "Method -> \"Residue\", Assumptions -> a > 0]", "Pi/Sqrt[2]");
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0, Method -> \"Residue\"]", "Pi/Sqrt[2]");
    /* Indefinite form accepts (and ignores) Assumptions rather than mis-reading
     * it as a Method value. */
    check_eq("Integrate[x^2, x, Assumptions -> a > 0]", "1/3 x^3");
    /* A genuinely unrecognised trailing option is still rejected (unevaluated). */
    check_eq("Integrate[1/(1+x^4), {x, -Infinity, Infinity}, Bogus -> 3]",
             "Integrate[1/(1 + x^4), {x, -Infinity, Infinity}, Bogus -> 3]");
}

/* -------------------------------------------------------------------------
 * Family B with symbolic parameters (Fourier/Jordan under Assumptions).
 * ---------------------------------------------------------------------- */
static void test_fourier_symbolic(void) {
    /* Integrate[Cos[k x]/(x^2+a^2)] = (Pi/a) E^{-a k},  a>0, k>0. */
    check_eq("Integrate[Cos[k x]/(x^2+a^2), {x, -Infinity, Infinity}, "
             "Assumptions -> {a > 0, k > 0}]", "(Pi E^(-a k))/a");
    /* Sin is odd -> 0. */
    check_eq("Integrate[Sin[k x]/(x^2+a^2), {x, -Infinity, Infinity}, "
             "Assumptions -> {a > 0, k > 0}]", "0");
    /* Numeric confirmation at a generic rational point. */
    check_eq("Chop[N[(Integrate[Cos[k x]/(x^2+a^2), {x, -Infinity, Infinity}, "
             "Assumptions -> {a > 0, k > 0}] - Pi E^(-a k)/a) /. {a -> 13/10, k -> 7/10}]]",
             "0");
    /* A different denominator scale. */
    check_eq("Chop[N[(Integrate[Cos[k x]/(x^2+b^2), {x, -Infinity, Infinity}, "
             "Assumptions -> {b > 0, k > 0}] - Pi E^(-b k)/b) /. {b -> 9/5, k -> 2}]]",
             "0");
    /* Bare complex-exponential kernel Exp[I k x] with symbolic parameters:
     * Integrate[Exp[I k x]/(x^2+a^2)] = (Pi/a) E^{-a k},  a>0, k>0. */
    check_eq("Integrate[Exp[I k x]/(x^2+a^2), {x, -Infinity, Infinity}, "
             "Assumptions -> {a > 0, k > 0}]", "(Pi E^(-a k))/a");
    check_eq("Chop[N[(Integrate[Exp[I k x]/(x^2+a^2), {x, -Infinity, Infinity}, "
             "Assumptions -> {a > 0, k > 0}] - Pi E^(-a k)/a) /. {a -> 12/10, k -> 9/10}]]",
             "0");
}

/* -------------------------------------------------------------------------
 * Rectangular contour: quasi-periodic Exp[c x] R(Exp[x]) on (-Inf, Inf).
 * ---------------------------------------------------------------------- */
static void test_rectangular(void) {
    /* Integrate[Exp[a x]/(Exp[x]+1)] = Pi/Sin[Pi a],  0 < a < 1. */
    check_eq("Integrate[Exp[a x]/(Exp[x]+1), {x, -Infinity, Infinity}, "
             "Assumptions -> 0 < a < 1]", "Pi Csc[Pi a]");
    check_eq("Chop[N[(Integrate[Exp[a x]/(Exp[x]+1), {x, -Infinity, Infinity}, "
             "Assumptions -> 0 < a < 1] - Pi/Sin[Pi a]) /. a -> 37/100]]", "0");
    /* Exp[a x]/(Exp[x]-1) style denominator is a genuine axis pole (Exp[x]=1 at
     * x=0): no clean value -> stays unevaluated. */
    check_eq("Integrate[Exp[a x]/(Exp[x]-1), {x, -Infinity, Infinity}, "
             "Assumptions -> 0 < a < 1]",
             "Integrate[E^(a x)/(-1 + E^x), {x, -Infinity, Infinity}, Assumptions -> 0 < a < 1]");
}

/* -------------------------------------------------------------------------
 * Keyhole / Mellin: branch power x^p R(x) on (0, Inf).
 * ---------------------------------------------------------------------- */
static void test_mellin(void) {
    check_eq("Integrate[x^(1/3)/(x^2+1), {x, 0, Infinity}]", "Pi/Sqrt[3]");
    check_eq("Integrate[Sqrt[x]/(x^2+1), {x, 0, Infinity}]", "Pi/Sqrt[2]");
    /* x^{-1/2}/(1+x) = Pi. */
    check_eq("Chop[N[Integrate[1/(Sqrt[x] (1+x)), {x, 0, Infinity}] - Pi]]", "0");
    /* Divergent branch power (s = 7/2 exceeds the decay order 2): unevaluated. */
    check_eq("Integrate[x^(5/2)/(x^2+1), {x, 0, Infinity}]",
             "Integrate[x^(5/2)/(1 + x^2), {x, 0, Infinity}]");
    /* Higher-order poles: the keyhole sum is over residues of the FULL integrand
     * x^(s-1) R(x), not x_k^(s-1) Res(R).  A pure double pole has Res(R) = 0, so
     * the old formula silently returned 0 -- these guard that regression.
     * B(3/2,1/2) = Pi/2, B(3/2,3/2) = Pi/8, B(4/3,2/3) = 2 Pi/(3 Sqrt[3]). */
    check_eq("Chop[N[Integrate[Sqrt[x]/(1+x)^2, {x, 0, Infinity}] - Pi/2]]", "0");
    check_eq("Chop[N[Integrate[Sqrt[x]/(1+x)^3, {x, 0, Infinity}] - Pi/8]]", "0");
    check_eq("Chop[N[Integrate[x^(1/3)/(1+x)^2, {x, 0, Infinity}] - 2 Pi/(3 Sqrt[3])]]", "0");
}

/* -------------------------------------------------------------------------
 * Sector contour: x^m/(c + x^n), symbolic exponent n.
 * ---------------------------------------------------------------------- */
static void test_sector(void) {
    /* Integrate[1/(1+x^n)] = (Pi/n) Csc[Pi/n],  n > 1. */
    check_eq("Integrate[1/(1+x^n), {x, 0, Infinity}, Assumptions -> n > 1]",
             "(Pi Csc[Pi/n])/n");
    check_eq("Chop[N[(Integrate[1/(1+x^n), {x, 0, Infinity}, Assumptions -> n > 1] "
             "- Pi/(n Sin[Pi/n])) /. n -> 23/10]]", "0");
    /* Numeric n, monomial numerator. */
    check_eq("Chop[N[Integrate[x/(1+x^4), {x, 0, Infinity}] - Pi/4]]", "0");
    check_eq("Chop[N[Integrate[1/(1+x^3), {x, 0, Infinity}] - 2 Pi/(3 Sqrt[3])]]", "0");
    /* Under-constrained exponent (only n > 0, so n <= 1 possible -> divergence
     * not excluded): the sector residue method must not fire.  The Mellin /
     * Ramanujan method, which runs later in the Automatic cascade, does close it
     * -- as a ConditionalExpression that states the missing convergence bound
     * (n > 1), matching Wolfram.  (Beta[1/n, 1-1/n]/n = (Pi/n) Csc[Pi/n].) */
    check_eq("Integrate[1/(1+x^n), {x, 0, Infinity}, Assumptions -> n > 0]",
             "ConditionalExpression[Beta[1/n, 1 - 1/n]/n, 1/n > 0 && 1/n < 1]");
}

/* -------------------------------------------------------------------------
 * Negative controls specific to the symbolic-parameter path.
 * ---------------------------------------------------------------------- */
static void test_symbolic_negative_controls(void) {
    /* Strict Method -> "Residue" (no Newton-Leibniz fallback) isolates the
     * residue method's own decision.  A free parameter left two-sided unbounded
     * by the assumptions does not determine the pole sign -> refuse. */
    check_eq("Integrate[Cos[k x]/(x^2+a^2), {x, -Infinity, Infinity}, "
             "Method -> \"Residue\", Assumptions -> k > 0]",
             "Integrate[Cos[k x]/(a^2 + x^2), {x, -Infinity, Infinity}, "
             "Method -> \"Residue\", Assumptions -> k > 0]");
    /* No assumptions at all: symbolic poles are undecidable -> refuse. */
    check_eq("Integrate[Cos[k x]/(x^2+a^2), {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[Cos[k x]/(a^2 + x^2), {x, -Infinity, Infinity}, Method -> \"Residue\"]");
}

/* -------------------------------------------------------------------------
 * Family A (rational) also fires for symbolic parameters under Assumptions and
 * closes to a clean rational form (not a Sqrt[-4 a^2] surface).
 * ---------------------------------------------------------------------- */
static void test_rational_symbolic(void) {
    check_eq("Integrate[1/(x^2+a^2), {x, -Infinity, Infinity}, Assumptions -> a > 0]",
             "Pi/a");
    check_eq("Chop[N[(Integrate[1/(x^2+a^2)^2, {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0] - Pi/(2 a^3)) /. a -> 7/5]]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_family_rational);
    TEST(test_family_fourier);
    TEST(test_family_trig);
    TEST(test_principal_value);
    TEST(test_half_line);
    TEST(test_dispatch);
    TEST(test_negative_controls);
    TEST(test_regression_finite);
    TEST(test_assumptions_option);
    TEST(test_fourier_symbolic);
    TEST(test_rectangular);
    TEST(test_mellin);
    TEST(test_sector);
    TEST(test_rational_symbolic);
    TEST(test_symbolic_negative_controls);

    printf("All Integrate ContourResidue tests passed!\n");
    return 0;
}
