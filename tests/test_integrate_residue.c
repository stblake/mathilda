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
    /* Higher-order NUMERIC pole: 1/(2+Cos[x])^2 = 4 Pi/(3 Sqrt[3]). */
    check_eq("Chop[N[Integrate[1/(2+Cos[x])^2, {x, 0, 2 Pi}] - 4 Pi/(3 Sqrt[3])]]", "0");
}

/* -------------------------------------------------------------------------
 * Unit-circle trig with SYMBOLIC parameters (Case 16): the in-disk pole is
 * classified at a FindInstance representative point of the assumption region
 * (here a coupled a > b > 0), and its residue is taken by the analytic-part
 * derivative so an order-n parametric radical pole stays fast.
 * ---------------------------------------------------------------------- */
static void test_trig_symbolic(void) {
    /* Order 1: Integrate[1/(a+b Cos[x]), a>b>0] = 2 Pi/Sqrt[a^2-b^2]. */
    check_eq("With[{r = Integrate[1/(a + b Cos[x]), {x, 0, 2 Pi}, Assumptions -> a > b > 0, "
             "Method -> \"Residue\"]}, {FreeQ[r, Integrate], "
             "Simplify[r - 2 Pi/Sqrt[a^2 - b^2]] == 0}]", "{True, True}");
    /* Order 3 (Case 16): = Pi (2a^2+b^2)/(a^2-b^2)^(5/2). */
    check_eq("With[{r = Integrate[1/(a + b Cos[x])^3, {x, 0, 2 Pi}, Assumptions -> a > b > 0, "
             "Method -> \"Residue\"]}, {FreeQ[r, Integrate], "
             "Simplify[r - Pi (2 a^2 + b^2)/(a^2 - b^2)^(5/2)] == 0}]", "{True, True}");

    /* Case In[11]: both Cos and Sin present, 1/(a + b Cos + c Sin)^2 =
     * 2 Pi a/(a^2-b^2-c^2)^(3/2), a > Sqrt[b^2+c^2].  Formerly leaked Power::infy
     * / Infinity::indet and declined: FindInstance picked the degenerate b=c=0,
     * collapsing the pole quadratic's leading coefficient.  Now a generic
     * (all-nonzero) representative is chosen and the internal N-probe runs under
     * Quiet.  (Value correct; the 4's are not pulled out of the ^(3/2), so pin
     * numerically.) */
    check_eq("With[{r = Integrate[1/(a + b Cos[theta] + c Sin[theta])^2, {theta, 0, 2 Pi}, "
             "Assumptions -> a > Sqrt[b^2 + c^2], Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[(r - 2 Pi a/(a^2 - b^2 - c^2)^(3/2)) "
             "/. {a -> 3, b -> 1, c -> 1}]]}]", "{True, 0}");
}

/* -------------------------------------------------------------------------
 * Finite interval (-1,1) with the Chebyshev weight 1/Sqrt[1-x^2] (Case 8):
 * x = Cos[t] -> (1/2) of the full-period trig integral.  Running before
 * Newton-Leibniz, it also pre-empts the FTC branch's WRONG-SIGN antiderivative.
 * ---------------------------------------------------------------------- */
static void test_chebyshev_weight(void) {
    /* Case 8: Integrate[1/((1+x^2) Sqrt[1-x^2]), {x,-1,1}] = Pi/Sqrt[2]
     * (Newton-Leibniz returns the wrong sign -Pi/Sqrt[2] here). */
    check_eq("Integrate[1/((1 + x^2) Sqrt[1 - x^2]), {x, -1, 1}]", "Pi/Sqrt[2]");
    check_eq("Integrate[1/((1 + x^2) Sqrt[1 - x^2]), {x, -1, 1}, Method -> \"Residue\"]",
             "Pi/Sqrt[2]");
    /* Weight only: Integrate[1/Sqrt[1-x^2], {x,-1,1}] = Pi. */
    check_eq("Integrate[1/Sqrt[1 - x^2], {x, -1, 1}, Method -> \"Residue\"]", "Pi");
    /* Another rational numerator: 1/((2+x^2) Sqrt[1-x^2]) = Pi/Sqrt[6]. */
    check_eq("Chop[N[Integrate[1/((2 + x^2) Sqrt[1 - x^2]), {x, -1, 1}, Method -> \"Residue\"] "
             "- Pi/Sqrt[6]]]", "0");
    /* No Chebyshev weight: the family declines, Newton-Leibniz owns it. */
    check_eq("Integrate[1/(1 + x^2), {x, -1, 1}, Method -> \"Residue\"]",
             "Integrate[1/(1 + x^2), {x, -1, 1}, Method -> \"Residue\"]");
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
    /* Higher-order pole with a SYMBOLIC parameter (regression). The order-2 pole
     * at z = I a formerly under-padded its Laurent series and dropped the
     * product-rule cross term, so Integrate[Cos[x]/(x^2+a^2)^2] returned the
     * wrong Pi E^-a/(2 a^2) instead of Pi (1+a) E^-a/(2 a^3). */
    check_eq("Integrate[Cos[x]/(x^2+a^2)^2, {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0, Method -> \"Residue\"]",
             "(1/2 Pi (1 + a) E^(-a))/a^3");
    check_eq("Chop[N[(Integrate[Cos[x]/(x^2+a^2)^2, {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0, Method -> \"Residue\"] "
             "- Pi (1 + a) E^(-a)/(2 a^3)) /. a -> 7/5]]", "0");

    /* Case In[14] -- a removable axis pole at x=0 coexisting with an enclosed
     * pole, symbolic parameters.  This FORMERLY RETURNED A WRONG 0: the Sin
     * extraction conjugated the contour value via ReplaceAll[I -> -I] (which
     * leaves a stored Complex[0,-1] atom untouched), so the enclosed-pole term
     * cancelled against itself; and FactorTerms fabricated a 0 numerical
     * content from a Together that left a b^-1 in the numerator.  The closing
     * now uses ComplexExpand[Im[...]] and FactorTerms never returns a 0 content
     * for a nonzero input.  Correct value (Pi/b^2)(1 - E^(-a b)). */
    check_eq("Integrate[Sin[a x]/(x (x^2+b^2)), {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]",
             "(Pi (1 - E^(-a b)))/b^2");
    check_eq("With[{r = Integrate[Sin[a x]/(x (x^2+b^2)), {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[(r - Pi (1 - E^(-a b))/b^2) /. {a -> 1, b -> 1}]]}]",
             "{True, 0}");

    /* Case In[15] -- complex-exponential kernel over a quadratic with two
     * upper-half-plane poles at I b +- Sqrt[c-b^2].  Value
     * -2 Pi E^(-a b) Sin[a Sqrt[c-b^2]]/Sqrt[c-b^2] (non-vacuous numeric pin). */
    check_eq("With[{r = Integrate[Exp[I a x]/(x^2 - 2 I b x - c), {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0 && b > 0 && c > b^2, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], "
             "Chop[N[(r + 2 Pi Exp[-a b] Sin[a Sqrt[c - b^2]]/Sqrt[c - b^2]) "
             "/. {a -> 1, b -> 1, c -> 4}]]}]", "{True, 0}");

    /* Case In[13] -- multi-frequency Fourier: a DIFFERENCE of two cosines whose
     * individual integrals diverge, so the problem is not separable.  The sum
     * Cos[a x] - Cos[b x] lifts to Re[(E^(I a x) - E^(I b x))], and the x=0
     * double pole of 1/x^2 is reduced to a simple one by the cancellation
     * (numerator -> 0 at x=0), admitted only because Limit[f, x->0] is finite.
     * Half-line Pi (b-a)/2, whole-line Pi (b-a). */
    check_eq("Integrate[(Cos[a x] - Cos[b x])/x^2, {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]", "1/2 Pi (-a + b)");
    check_eq("Integrate[(Cos[a x] - Cos[b x])/x^2, {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]", "Pi (-a + b)");
    check_eq("With[{r = Integrate[(Cos[a x] - Cos[b x])/x^2, {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[(r - Pi (b - a)/2) /. {a -> 2, b -> 5}]]}]",
             "{True, 0}");
    /* The single Cos[a x]/x^2 (no cancellation -> genuine x=0 double pole, the
     * integral diverges) must DECLINE, not fabricate a finite value. */
    check_eq("Integrate[Cos[a x]/x^2, {x, 0, Infinity}, Assumptions -> a > 0, "
             "Method -> \"Residue\"]",
             "Integrate[Cos[a x]/x^2, {x, 0, Infinity}, Assumptions -> a > 0, "
             "Method -> \"Residue\"]");
    /* Sibling Sin difference: both pieces integrate to Pi/2, difference 0. */
    check_eq("Integrate[(Sin[a x] - Sin[b x])/x, {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]", "0");
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
 * Symbolic branch-power exponent: x^a R(x) with a a free parameter.  The
 * convergence interval is taken from the exponent (s = a+1), so a symbolic a is
 * no longer wrongly refused.
 * ---------------------------------------------------------------------- */
static void test_mellin_symbolic_exponent(void) {
    /* Integrate[x^a/(x+1)^3, -1<a<2] = Pi a (1-a)/(2 Sin[Pi a]). */
    check_eq("Integrate[x^a/(x+1)^3, {x, 0, Infinity}, Assumptions -> -1 < a < 2, "
             "Method -> \"Residue\"]", "1/2 Pi a (-1 + a) Csc[Pi (1 + a)]");
    check_eq("Chop[N[(Integrate[x^a/(x+1)^3, {x, 0, Infinity}, Assumptions -> -1 < a < 2, "
             "Method -> \"Residue\"] - Pi a (1 - a)/(2 Sin[Pi a])) /. a -> 1/2]]", "0");
}

/* -------------------------------------------------------------------------
 * Keyhole with a logarithm: Integrate[x^p (Log x)^m R(x), {x,0,Inf}], p a
 * non-negative integer, m >= 1.  The (Log z)^(k+1) contour gives a triangular
 * system solved for I_0..I_m from residues alone.
 * ---------------------------------------------------------------------- */
static void test_keyhole_log(void) {
    /* Each pin asserts BOTH that the residue method actually evaluated
     * (FreeQ[r, Integrate] -- so a silent decline + N[] NIntegrate fallback
     * cannot pass vacuously) AND that the closed form is numerically correct. */

    /* Integrate[Log[x]/(1+x^6)] = -Sqrt[3] Pi^2/18 (m=1, six simple poles on the
     * sixth roots of -1 -- needs exact Arg on those roots to collapse). */
    check_eq("With[{r = Integrate[Log[x]/(1+x^6), {x, 0, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r + Sqrt[3] Pi^2/18]]}]", "{True, 0}");
    /* Integrate[(Log x)^2/(x^2+x+1)] = 16 Pi^3/(81 Sqrt[3]) (m=2, primitive cube
     * roots). */
    check_eq("With[{r = Integrate[Log[x]^2/(x^2+x+1), {x, 0, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - 16 Pi^3/(81 Sqrt[3])]]}]", "{True, 0}");
    /* Higher-order pole with a log: Integrate[Log[x]/(1+x^2)^2] = -Pi/4. */
    check_eq("Integrate[Log[x]/(1+x^2)^2, {x, 0, Infinity}, Method -> \"Residue\"]", "-1/4 Pi");
    /* Integer power folded in (p=2): Integrate[x^2 Log[x]/(1+x^6)] = 0 by the
     * x -> 1/x symmetry of x^2/(1+x^6). */
    check_eq("With[{r = Integrate[x^2 Log[x]/(1+x^6), {x, 0, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r]]}]", "{True, 0}");
    /* A pure log integral: Integrate[Log[x]/(1+x^4)] = -Sqrt[2] Pi^2/16 (fourth
     * roots of -1, carrying an I factor). */
    check_eq("With[{r = Integrate[Log[x]/(1+x^4), {x, 0, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r + Sqrt[2] Pi^2/16]]}]", "{True, 0}");

    /* Case 14 -- NON-integer branch power with a log: the (1 - e^(2 Pi i a))
     * keyhole.  Integrate[Sqrt[x] Log[x]/(x^2+1)^2] = Pi(Pi-4)/(8 Sqrt[2])
     * (a = 1/2, m = 1, double poles at +/- i). */
    check_eq("With[{r = Integrate[Sqrt[x] Log[x]/(x^2+1)^2, {x, 0, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - Pi (Pi - 4)/(8 Sqrt[2])]]}]", "{True, 0}");
}

/* -------------------------------------------------------------------------
 * Keyhole-log with a pole on the branch cut (0, Inf): the Cauchy principal
 * value, each axis pole contributing the average of its two keyhole-branch
 * residues.  A simple pole at z = 1 under a Log factor is REMOVABLE (Log 1 = 0),
 * so it is admitted without the option; a genuine axis pole needs
 * PrincipalValue -> True.
 * ---------------------------------------------------------------------- */
static void test_keyhole_pv(void) {
    /* Case 18: Integrate[Log[x]/(x^3-1)] = 4 Pi^2/27.  Removable pole at x = 1
     * (the integrand is in fact continuous there), admitted with no option. */
    check_eq("With[{r = Integrate[Log[x]/(x^3-1), {x, 0, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - 4 Pi^2/27]]}]", "{True, 0}");
    check_eq("Integrate[Log[x]/(x^3-1), {x, 0, Infinity}, Method -> \"Residue\"]", "4/27 Pi^2");
    /* A genuine (non-removable) axis pole at x = 2: declines without the option,
     * returns the principal value Pi^2/8 with it. */
    check_eq("Integrate[Log[x]/(x^2-4), {x, 0, Infinity}, Method -> \"Residue\"]",
             "Integrate[Log[x]/(-4 + x^2), {x, 0, Infinity}, Method -> \"Residue\"]");
    check_eq("With[{r = Integrate[Log[x]/(x^2-4), {x, 0, Infinity}, Method -> \"Residue\", "
             "PrincipalValue -> True]}, {FreeQ[r, Integrate], Chop[N[r - Pi^2/8]]}]", "{True, 0}");

    /* Case 7: Integrate[x/Sinh[x], {x, -Inf, Inf}] = Pi^2/2.  The rectangular
     * substitution w = Exp[x] turns x/Sinh x into the keyhole-log integrand
     * 2 Log[w]/(w^2-1) with a removable pole at w = 1 (x = 0). */
    check_eq("With[{r = Integrate[x/Sinh[x], {x, -Infinity, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - Pi^2/2]]}]", "{True, 0}");
    check_eq("Integrate[x/Sinh[x], {x, -Infinity, Infinity}, Method -> \"Residue\"]", "1/2 Pi^2");
}

/* -------------------------------------------------------------------------
 * Gaussian / shifted-rectangle on (-Inf, Inf): an entire E^(quadratic) kernel,
 * TrigToExp'd to a sum of pure Gaussians and closed by completing the square.
 * ---------------------------------------------------------------------- */
static void test_gaussian(void) {
    /* Case 15: Integrate[Exp[-x^2] Cos[2 a x]] = Sqrt[Pi] E^(-a^2).  Symbolic a
     * (formerly HUNG); concrete a = 1 declined.  Both close now. */
    check_eq("Integrate[Exp[-x^2] Cos[2 x], {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Sqrt[Pi]/E");
    check_eq("Integrate[Exp[-x^2] Cos[2 a x], {x, -Infinity, Infinity}, Method -> \"Residue\", "
             "Assumptions -> a > 0]", "Sqrt[Pi] E^(-a^2)");
    /* Shifted exponent, Exp kernel, and a scaled leading coefficient. */
    check_eq("Integrate[Exp[-x^2 + 3 x], {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "E^(9/4) Sqrt[Pi]");
    check_eq("With[{r = Integrate[Exp[-x^2 + I x], {x, -Infinity, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - Sqrt[Pi] Exp[-1/4]]]}]", "{True, 0}");
    check_eq("With[{r = Integrate[Exp[-2 x^2] Cos[x], {x, -Infinity, Infinity}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - Sqrt[Pi/2] Exp[-1/8]]]}]", "{True, 0}");
    /* Odd kernel integrates to 0. */
    check_eq("Integrate[Exp[-x^2] Sin[2 x], {x, -Infinity, Infinity}, Method -> \"Residue\"]", "0");
    /* Positive leading coefficient diverges: the family must decline. */
    check_eq("Integrate[Exp[x^2], {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "Integrate[E^x^2, {x, -Infinity, Infinity}, Method -> \"Residue\"]");

    /* Case In[5]: SYMBOLIC leading coefficient Exp[-a x^2] Cos[b x].  The
     * convergence gate now proves Re A = -a < 0 via Refine (the interval test
     * saw the compound -a as unbounded), and Exponent handles the non-flat
     * Times from the complex Fourier frequency.  Whole line = Sqrt[Pi/a]
     * E^(-b^2/(4a)); even half-line = half of it. */
    check_eq("Integrate[Exp[-a x^2] Cos[b x], {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]",
             "Sqrt[Pi/a] E^(-(1/4 b^2)/a)");
    check_eq("Integrate[Exp[-a x^2] Cos[b x], {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]",
             "1/2 Sqrt[Pi/a] E^(-(1/4 b^2)/a)");
    check_eq("With[{r = Integrate[Exp[-a x^2] Cos[b x], {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[(r - Sqrt[Pi/a] Exp[-b^2/(4 a)]/2) "
             "/. {a -> 1, b -> 3}]]}]", "{True, 0}");
    /* The dropped Element[b, Reals] assumption must not block it: the Refine gate
     * works on the FindInstance representative path too. */
    check_eq("Integrate[Exp[-a x^2] Cos[b x], {x, 0, Infinity}, "
             "Assumptions -> a > 0 && Element[b, Reals], Method -> \"Residue\"]",
             "1/2 Sqrt[Pi/a] E^(-(1/4 b^2)/a)");
}

/* -------------------------------------------------------------------------
 * Honest declines (deliberate, documented): integrands outside the residue
 * repertoire stay UNEVALUATED under strict Method -> "Residue" -- never a wrong
 * or forced value.  Pinned so a future accidental "answer" is caught.
 * ---------------------------------------------------------------------- */
static void test_honest_declines(void) {
    /* Case 17: Integrate[Log[x]/Cosh[x]] -- the Mellin transform of Sech
     * differentiated at s = 1 (a Dirichlet-beta derivative / Gamma[1/4] constant),
     * not a residue sum; the w = Exp[x] reduction fails (Log x -> Log[Log w]). */
    check_eq("Integrate[Log[x]/Cosh[x], {x, 0, Infinity}, Method -> \"Residue\"]",
             "Integrate[Log[x] Sech[x], {x, 0, Infinity}, Method -> \"Residue\"]");
    /* Case 21: a Hankel-type contour fragment, ambiguous / divergent in the
     * ordinary sense. */
    check_eq("Integrate[Exp[x] x^(-s), {x, 1, -Infinity}, Method -> \"Residue\"]",
             "Integrate[E^x x^(-s), {x, 1, -Infinity}, Method -> \"Residue\"]");
    /* Case In[20]: ArcTan[a x]/(x (1+b^2 x^2)) on {0, Infinity} = (Pi/2) Log[1+a/b].
     * The integrand has a BRANCH CUT (ArcTan), not isolated poles; the value is
     * reached by parametric differentiation (Feynman: d/da of the integral is a
     * rational contour integral Pi/(2(a+b)), then integrate in a), which is not
     * the residue theorem.  A deliberate, honest decline -- no Feynman engine is
     * built on the residue path.  (Mathematica returns the closed form via a
     * different route.) */
    check_eq("Integrate[ArcTan[a x]/(x (1 + b^2 x^2)), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]",
             "Integrate[ArcTan[a x]/(x (1 + b^2 x^2)), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]");
}

/* -------------------------------------------------------------------------
 * Parametrized contour on (0, 2Pi): Integrate[g(c Exp[I t]) (I c Exp[I t])] =
 * Contour[g, |z|=c] = 2 Pi i Sum Res, including essential singularities.  Each
 * pin asserts FreeQ[r, Integrate] so a decline cannot pass via a numeric
 * fallback; the value is genuinely complex (2 Pi i), not real.
 * ---------------------------------------------------------------------- */
static void test_contour_param(void) {
    /* Case 9/19: Contour[Exp[1/z] Sin[1/z], |z|=1] -- essential singularity at 0,
     * Res = 1 (w^1 coeff of Exp[w] Sin[w]).  = 2 Pi i. */
    check_eq("With[{r = Integrate[Exp[Exp[-I t]] Sin[Exp[-I t]] I Exp[I t], {t, 0, 2 Pi}, "
             "Method -> \"Residue\"]}, {FreeQ[r, Integrate], Chop[N[r - 2 Pi I]]}]", "{True, 0}");
    /* Case 20: Contour[Exp[2z]/(z^4 (z-1/2)), |z|=1] -- order-4 pole at 0 (Res
     * -128/3) + simple pole at 1/2 (Res 16 E).  = 2 Pi i (16 E - 128/3). */
    check_eq("With[{r = Integrate[(Exp[2 Exp[I t]]/(Exp[4 I t] (Exp[I t] - 1/2))) I Exp[I t], "
             "{t, 0, 2 Pi}, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[r - 2 Pi I (16 E - 128/3)]]}]", "{True, 0}");
    /* Case 22: Contour[z^5 Cos[1/z^2], |z|=2] -- essential singularity, Res = 0
     * (no z^-1 term in z^5 Sum (-1)^k z^(-4k)/(2k)!).  = 0. */
    check_eq("With[{r = Integrate[(2 Exp[I t])^5 Cos[1/(2 Exp[I t])^2] 2 I Exp[I t], "
             "{t, 0, 2 Pi}, Method -> \"Residue\"]}, {FreeQ[r, Integrate], Chop[N[r]]}]",
             "{True, 0}");
    /* The ordinary rational-in-{Cos,Sin} trig family is unaffected. */
    check_eq("Integrate[1/(2 + Cos[t]), {t, 0, 2 Pi}, Method -> \"Residue\"]", "(2 Pi)/Sqrt[3]");
}

/* -------------------------------------------------------------------------
 * Mellin-Barnes / Bromwich vertical line: Integrate[F(s) x^-s, {s, c-iInf,
 * c+iInf}] = 2 Pi i Sum over F's left poles.  The limit c +/- I Infinity must
 * survive as a directed infinity (plus.c classify fix) for the line to be seen.
 * ---------------------------------------------------------------------- */
static void test_mellin_barnes(void) {
    /* Case 11: inverse Mellin of Gamma[s] is e^-x; = 2 Pi i e^-x (was a silent 0). */
    check_eq("With[{r = Integrate[Gamma[s] x^(-s), {s, 1/2 - I Infinity, 1/2 + I Infinity}]}, "
             "{FreeQ[r, Integrate], Simplify[r - 2 Pi I Exp[-x]] == 0}]", "{True, True}");
    /* Gamma[2 s] x^-s = I Pi e^-Sqrt[x] (B = 2 scales the pole ladder). */
    check_eq("With[{r = Integrate[Gamma[2 s] x^(-s), {s, 1/2 - I Infinity, 1/2 + I Infinity}]}, "
             "{FreeQ[r, Integrate], Simplify[r - I Pi Exp[-Sqrt[x]]] == 0}]", "{True, True}");
    /* Orientation: the downward line is the negative. */
    check_eq("With[{r = Integrate[Gamma[s] x^(-s), {s, 1/2 + I Infinity, 1/2 - I Infinity}]}, "
             "{FreeQ[r, Integrate], Simplify[r + 2 Pi I Exp[-x]] == 0}]", "{True, True}");
    /* Negative control: no Gamma-pole ladder -> declines (stays unevaluated). */
    check_eq("Head[Integrate[1/s, {s, 1/2 - I Infinity, 1/2 + I Infinity}, Method -> \"Residue\"]]",
             "Integrate");
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

    /* Case In[19]: SYMBOLIC numerator and denominator exponents,
     * x^(2m)/(1 + x^(2n)) = (Pi/(2n)) Csc[Pi(2m+1)/(2n)].  The convergence
     * s = 2m+1 < 2n must be provable: state it as n >= m+1 (the integer form of
     * n > m; Refine lacks the integer-gap step n > m && Integers => n >= m+1). */
    check_eq("Integrate[x^(2 m)/(1 + x^(2 n)), {x, 0, Infinity}, "
             "Assumptions -> Element[m, Integers] && Element[n, Integers] && n >= m + 1 && m >= 0, "
             "Method -> \"Residue\"]", "(1/2 Pi Csc[(1/2 Pi (1 + 2 m))/n])/n");
    check_eq("Chop[N[(Integrate[x^(2 m)/(1 + x^(2 n)), {x, 0, Infinity}, "
             "Assumptions -> Element[m, Integers] && Element[n, Integers] && n >= m + 1 && m >= 0, "
             "Method -> \"Residue\"] - Pi/(2 n) Csc[Pi(2 m+1)/(2 n)]) /. {m -> 1, n -> 2}]]", "0");
}

/* Generalized Beta on (0, Inf): x^a/(x+b)^c, non-integer c (branch point at -b). */
static void test_beta(void) {
    /* Case In[16]: x^a/(x+b)^c = b^(a+1-c) Gamma[a+1] Gamma[c-a-1]/Gamma[c]. */
    check_eq("Integrate[x^a/(x+b)^c, {x, 0, Infinity}, "
             "Assumptions -> b > 0 && c > a + 1 && a > -1, Method -> \"Residue\"]",
             "(Gamma[1 + a] Gamma[-1 - a + c] b^(1 + a - c))/Gamma[c]");
    check_eq("Chop[N[(Integrate[x^a/(x+b)^c, {x, 0, Infinity}, "
             "Assumptions -> b > 0 && c > a + 1 && a > -1, Method -> \"Residue\"] "
             "- b^(a+1-c) Beta[a+1, c-a-1]) /. {a -> 1/2, b -> 1, c -> 3}]]", "0");
    check_eq("Chop[N[(Integrate[x^a/(x+b)^c, {x, 0, Infinity}, "
             "Assumptions -> b > 0 && c > a + 1 && a > -1, Method -> \"Residue\"] "
             "- b^(a+1-c) Beta[a+1, c-a-1]) /. {a -> 1, b -> 2, c -> 4}]]", "0");
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
/* -------------------------------------------------------------------------
 * Mellin after a power substitution u = x^nu: x^(mu-1) G(kappa x^nu) ->
 * (1/nu) * Mellin[G](mu/nu), G in {Exp, Sin, Cos} (the Gamma / generalized
 * Fresnel family).
 * ---------------------------------------------------------------------- */
static void test_mellin_power(void) {
    /* Case In[9]: x^(s-1) e^(-a x^2) = (1/2) a^(-s/2) Gamma[s/2]  (convergence
     * stated as s > 0, not the Refine-unfriendly complex spelling Re[s] > 0). */
    check_eq("Integrate[x^(s-1) Exp[-a x^2], {x, 0, Infinity}, "
             "Assumptions -> s > 0 && a > 0, Method -> \"Residue\"]",
             "1/2 Gamma[1/2 s] a^(-1/2 s)");
    check_eq("Chop[N[(Integrate[x^(s-1) Exp[-a x^2], {x,0,Infinity}, "
             "Assumptions -> s>0 && a>0, Method -> \"Residue\"] - a^(-s/2) Gamma[s/2]/2) "
             "/. {s -> 3/2, a -> 1}]]", "0");
    /* Case In[10]: x^p Sin[x^2] = (1/2) Gamma[(p+1)/2] Sin[Pi(p+1)/4]. */
    check_eq("Integrate[x^p Sin[x^2], {x, 0, Infinity}, "
             "Assumptions -> -1 < p < 1, Method -> \"Residue\"]",
             "1/2 Gamma[1/2 (1 + p)] Sin[1/4 Pi (1 + p)]");
    /* The Cos sibling. */
    check_eq("Integrate[x^p Cos[x^2], {x, 0, Infinity}, "
             "Assumptions -> -1 < p < 1, Method -> \"Residue\"]",
             "1/2 Gamma[1/2 (1 + p)] Cos[1/4 Pi (1 + p)]");
    /* A plain generalized-Fresnel (mu = 1): Sin[x^3] = Gamma[4/3] Sin[Pi/6] = Gamma[4/3]/2. */
    check_eq("Chop[N[Integrate[Sin[x^3], {x, 0, Infinity}, Method -> \"Residue\"] "
             "- Gamma[4/3]/2]]", "0");
}

/* -------------------------------------------------------------------------
 * Periodic-strip (rectangle) contour: N(x)/Cosh[b x] on (-Inf, Inf) via the
 * quasi-period Cosh[b(x+i Pi/b)] = -Cosh[b x], and the scale-normalisation
 * Sinh[a x] -> Sinh[x] that reduces x/Sinh[a x] to the a=1 w=Exp[x] route.
 * ---------------------------------------------------------------------- */
static void test_hyperbolic_strip(void) {
    /* Case In[4]: E^(a x)/Cosh[Pi x] = Sec[a/2]  (b = Pi concrete). */
    check_eq("Integrate[Exp[a x]/Cosh[Pi x], {x, -Infinity, Infinity}, "
             "Assumptions -> -Pi < a < Pi, Method -> \"Residue\"]", "Sec[1/2 a]");
    /* Case In[21]: Cosh[a x]/Cosh[b x] = (Pi/b) Sec[Pi a/(2 b)]. */
    check_eq("Integrate[Cosh[a x]/Cosh[b x], {x, -Infinity, Infinity}, "
             "Assumptions -> b > a > 0, Method -> \"Residue\"]",
             "(Pi Sec[(1/2 Pi a)/b])/b");
    /* Case In[18]: E^(I a x)/Cosh[b x] = (Pi/b) Sech[Pi a/(2 b)]. */
    check_eq("Integrate[Exp[I a x]/Cosh[b x], {x, -Infinity, Infinity}, "
             "Assumptions -> Element[a, Reals] && b > 0, Method -> \"Residue\"]",
             "(Pi Sech[(1/2 Pi a)/b])/b");
    /* Bare Sech integrates to Pi/b (the alpha = 0 term). */
    check_eq("Integrate[Sech[Pi x], {x, -Infinity, Infinity}, Method -> \"Residue\"]",
             "1");
    /* Odd numerator over Cosh integrates to 0. */
    check_eq("Integrate[Sinh[a x]/Cosh[b x], {x, -Infinity, Infinity}, "
             "Assumptions -> b > a > 0, Method -> \"Residue\"]", "0");
    /* Numeric confirmations. */
    check_eq("Chop[N[(Integrate[Exp[a x]/Cosh[Pi x], {x,-Infinity,Infinity}, "
             "Assumptions -> -Pi<a<Pi, Method -> \"Residue\"] - Sec[a/2]) /. a -> 1]]", "0");
    check_eq("Chop[N[(Integrate[Cosh[a x]/Cosh[b x], {x,-Infinity,Infinity}, "
             "Assumptions -> b>a>0, Method -> \"Residue\"] - Pi Sec[Pi a/(2 b)]/b) "
             "/. {a -> 1, b -> 2}]]", "0");

    /* Case In[12]: x/Sinh[a x] on [0,Inf) = Pi^2/(4 a^2).  Even reduction ->
     * whole line; the scale-normalisation u = a x reduces Sinh[a x] to Sinh[x]
     * so the w=Exp[x] keyhole route closes.  The symbolic-a closed form carries
     * an (unsimplified, but correct) keyhole Arg surface, so pin numerically. */
    check_eq("With[{r = Integrate[x/Sinh[a x], {x, 0, Infinity}, "
             "Assumptions -> a > 0, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[(r - Pi^2/(4 a^2)) /. a -> 3/2]]}]",
             "{True, 0}");
    /* Concrete scale is fully simplified. */
    check_eq("Integrate[x/Sinh[2 x], {x, 0, Infinity}, Method -> \"Residue\"]",
             "1/16 Pi^2");
}

static void test_rational_symbolic(void) {
    check_eq("Integrate[1/(x^2+a^2), {x, -Infinity, Infinity}, Assumptions -> a > 0]",
             "Pi/a");
    check_eq("Chop[N[(Integrate[1/(x^2+a^2)^2, {x, -Infinity, Infinity}, "
             "Assumptions -> a > 0] - Pi/(2 a^3)) /. a -> 7/5]]", "0");

    /* Case In[7]: two conjugate pole-pairs on [0, Inf).  solve_roots returns the
     * quartic roots unfactored, and RootReduce canonicalises the residue sum into
     * nested surds (Sqrt[(a+b)^2 (a-b)^2]) that plain Simplify leaves standing;
     * the all-positive PowerExpand/Refine close in close_algebraic collapses them
     * to the clean Pi/(2 a b (a+b)).  (Both with and without a != b.) */
    check_eq("Integrate[1/((x^2+a^2)(x^2+b^2)), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0 && a != b, Method -> \"Residue\"]",
             "(1/2 Pi)/(a^2 b + a b^2)");
    check_eq("With[{r = Integrate[1/((x^2+a^2)(x^2+b^2)), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"Residue\"]}, "
             "{FreeQ[r, Integrate], Chop[N[(r - Pi/(2 a b (a+b))) /. {a -> 1, b -> 2}]]}]",
             "{True, 0}");
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
    TEST(test_mellin_symbolic_exponent);
    TEST(test_keyhole_log);
    TEST(test_keyhole_pv);
    TEST(test_gaussian);
    TEST(test_honest_declines);
    TEST(test_contour_param);
    TEST(test_mellin_barnes);
    TEST(test_sector);
    TEST(test_beta);
    TEST(test_trig_symbolic);
    TEST(test_chebyshev_weight);
    TEST(test_hyperbolic_strip);
    TEST(test_mellin_power);
    TEST(test_rational_symbolic);
    TEST(test_symbolic_negative_controls);

    printf("All Integrate ContourResidue tests passed!\n");
    return 0;
}
