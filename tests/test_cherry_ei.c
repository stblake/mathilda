/* test_cherry_ei.c — Cherry's rational exponential integral: ExpIntegralEi engine.
 *
 * Pins for cherry_ei.c (CHERRY_PLAN.md §3): the base-field ei examples from
 * G. W. Cherry, "An Analysis of the Rational Exponential Integral" (1989), plus
 * the Macsyma d8/d11 examples from the 1986 companion, and the decline-safety
 * cases (algebraic ei-argument constants are a later phase and must NOT emit a
 * partial/wrong form).  Each positive pin is diff-back verified and confirmed to
 * contain ExpIntegralEi (so it is the ei engine that closed it, not a fluke).
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

static Expr* eval_str(const char* s) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool eval_is(const char* s, const char* expected) {
    Expr* r = eval_str(s);
    char* got = expr_to_string_fullform(r);
    bool ok = strcmp(got, expected) == 0;
    if (!ok) printf("  [%s] -> %s (expected %s)\n", s, got, expected);
    free(got);
    expr_free(r);
    return ok;
}

/* Integrate f via RischTranscendental; assert it (a) closes to an ExpIntegralEi
 * form and (b) diff-backs exactly to f. */
static void assert_ei(const char* f) {
    char buf[1400];
    /* (a) result contains ExpIntegralEi and is not left unevaluated */
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x, Method -> \"RischTranscendental\"]},"
        " Head[r] =!= Integrate && !FreeQ[r, ExpIntegralEi]]", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a closed ExpIntegralEi form", f);
    /* (b) exact diff-back */
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate[%s, x, Method -> \"RischTranscendental\"], x] - (%s)]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: diff-back nonzero", f);
}

/* Assert f closes to an Erfi (error-function) form and diff-backs exactly. */
static void assert_erf(const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate[%s, x, Method -> \"RischTranscendental\"]},"
        " Head[r] =!= Integrate && !FreeQ[r, Erfi]]", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a closed Erfi form", f);
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate[%s, x, Method -> \"RischTranscendental\"], x] - (%s)]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: diff-back nonzero", f);
}

/* Assert f closes (not unevaluated) and diff-backs exactly — WITHOUT requiring an
 * ExpIntegralEi in the result (for integrands whose ei terms cancel to an
 * elementary closed form, or pure-elementary siblings). */
static void assert_closes(const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "Head[Integrate[%s, x, Method -> \"RischTranscendental\"]] =!= Integrate", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected to close", f);
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate[%s, x, Method -> \"RischTranscendental\"], x] - (%s)]", f, f);
    ASSERT_MSG(eval_is(buf, "0"), "%s: diff-back nonzero", f);
}

/* Assert f is left unevaluated by the transcendental engine (clean decline). */
static void assert_declines(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Head[Integrate`RischTranscendental[%s, x]] === Integrate`RischTranscendental", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a clean decline", f);
}

/* Cherry 1989 base-field ei examples. */
static void test_cherry_1989_ei(void) {
    /* Ex 5.1: INT e^(1/x) dx = x E^(1/x) - ei(1/x)  (P2 q-side term, alpha = 0). */
    assert_ei("E^(1/x)");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[E^(1/x), x, Method -> \"RischTranscendental\"]"
        " - (x E^(1/x) - ExpIntegralEi[1/x])]", "0"),
        "Ex 5.1 exact form");

    /* e^x/x^2 = ei(x) - e^x/x  (P1: double resultant root alpha = 0). */
    assert_ei("E^x/x^2");
}

/* Cherry 1986 companion Macsyma examples d8, d11 (rational-constant ei). */
static void test_cherry_1986_macsyma(void) {
    /* d8: INT e^x/(x+1)^2 dx = e^(-1) ei(x+1) - e^x/(x+1). */
    assert_ei("E^x/(x + 1)^2");
    /* d11: INT (x^2+3) e^x/(x^2+3x+2) dx = -7 e^(-2) ei(x+2) + 4 e^(-1) ei(x+1) + e^x.
     * Multi-term ei via the P1 resultant roots alpha = 1, 2. */
    assert_ei("(x^2 + 3) E^x/(x^2 + 3 x + 2)");
    /* A shifted single-pole sibling: e^x/(x+2)^2. */
    assert_ei("E^x/(x + 2)^2");
}

/* Complex-quadratic ei: a single irreducible-quadratic denominator with COMPLEX
 * conjugate roots (f a polynomial, q = 1) closes with a complex-conjugate pair of
 * ExpIntegralEi terms — the diff-back over Q(i) is exact and fast. */
static void test_complex_ei(void) {
    /* Q(i) conjugate pairs (roots p + q I, q rational: NO radical) close fast —
     * the coefficient solve and diff-back over Q(i) do not trip the generic GCD. */
    assert_ei("E^x/(x^2 + 1)");        /* alpha = +-I         (Q(i))    */
    assert_ei("E^x/(x^2 + 9)");        /* alpha = +-3 I       (Q(i))    */
    assert_ei("E^(2 x)/(x^2 + 1)");    /* scaled exponent               */
    assert_ei("E^x/(x^2 + 2 x + 2)");  /* alpha = -1 +- I     (Q(i))    */
    assert_ei("E^x/(x^2 + 2 x + 5)");  /* alpha = -1 -+ 2 I   (Q(i))    */

    /* Q(i sqrt d) conjugate pairs (roots carry Sqrt[d], d not a perfect square):
     * DECLINE cleanly.  The coefficient solve / diff-back for these route through
     * Simplify/Together over Q(i sqrt d) with the E^x kernel, where the generic
     * multivariate GCD blows up in exact_poly_div (a pre-existing hang).  The
     * engine gates these complex-radical roots out (cherry_ei.c expr_has_radical)
     * so Integrate returns unevaluated rather than hanging.  Restoring them to a
     * closed form needs the FLINT number-field GCD (CHERRY_BLOCKERS A1, deferred). */
    assert_declines("E^x/(x^2 + x + 1)");        /* (1 -+ I sqrt3)/2  Q(i sqrt3) */
    assert_declines("E^x/(x^2 + 3)");            /* -+ I sqrt3        Q(i sqrt3) */
    assert_declines("(2 x + 1) E^x/(x^2 + x + 3)"); /*                Q(i sqrt11) */
    assert_declines("x E^x/(x^2 + x + 1)");      /*                   Q(i sqrt3) */
    assert_declines("(x^2 + 1) E^x/(x^2 + x + 1)"); /* d12            Q(i sqrt3) */
    assert_declines("E^x/(x^2 + 2 x + 3)");      /*                   Q(i sqrt2) */
    assert_declines("E^x/(x^2 + 2 x + 7)");      /*                   Q(i sqrt6) */
    assert_declines("E^x/(x^2 - 2 x + 3)");
    assert_declines("(3 x + 1) E^x/(x^2 + 2 x + 3)");
}

/* Constant exponent offset: E^(c + h(x)) = E^c E^(h(x)) — the constant folds into the
 * cofactor rather than inflating deg(p) and defeating P2. */
static void test_const_offset(void) {
    assert_ei("E^(1/x + 2)");            /* E^2 (x E^(1/x) - ei(1/x)) */
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[E^(1/x + 2), x, Method -> \"RischTranscendental\"]"
        " - E^2 (x E^(1/x) - ExpIntegralEi[1/x])]", "0"), "E^(1/x+2) exact form");
    assert_ei("E^((x - 1)/x)");          /* = E^(1 - 1/x) */
    assert_ei("E^(1/x - 1)");
}

/* Nonlinear exponents f = p/q close too — the engine is NOT limited to linear f. */
static void test_nonlinear_exponent(void) {
    /* v = x^2:  INT e^(x^2)/x dx = (1/2) ei(x^2)  (P1 root alpha = 0). */
    assert_ei("E^(x^2)/x");
}

/* Algebraic-constant layer (C2 §7): REAL algebraic ei-argument constants
 * (Cherry's algebraically-closed constant field).  The resultant roots alpha_i
 * are irrational reals (sqrt d, golden ratio, ...); the matching system is solved
 * over the extension they generate and the answer carries e^(-alpha_i) prefactors
 * with algebraic exponents. */
static void test_cherry_algebraic_ei(void) {
    /* Cherry p.894: INT e^x/(x^2-2) dx
     *   = (1/(2 sqrt2)) (e^sqrt2 ei(x-sqrt2) - e^(-sqrt2) ei(x+sqrt2)). */
    assert_ei("E^x/(x^2 - 2)");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[E^x/(x^2-2), x, Method -> \"RischTranscendental\"]"
        " - (E^Sqrt[2] ExpIntegralEi[x-Sqrt[2]] - E^(-Sqrt[2]) ExpIntegralEi[x+Sqrt[2]])"
        "   /(2 Sqrt[2])]", "0"),
        "p.894 exact form");
    assert_ei("E^x/(x^2 - 3)");         /* sqrt3 */
    assert_ei("E^x/(x^2 - x - 1)");     /* golden ratio (1 +- sqrt5)/2 */
    assert_ei("E^(2 x)/(x^2 - 2)");     /* scaled exponent, sqrt2 */
    assert_ei("E^x/(x^2 - 5)");         /* sqrt5 */
}

/* The ei answers are non-elementary: ElementaryIntegralQ stays False even though
 * a closed ei form is produced. */
static void test_elementaryq_false(void) {
    ASSERT_MSG(eval_is("Risch`ElementaryIntegralQ[E^(1/x), x]", "False"),
        "ElementaryIntegralQ[E^(1/x)] should be False");
    ASSERT_MSG(eval_is("Risch`ElementaryIntegralQ[E^x/(x + 1)^2, x]", "False"),
        "ElementaryIntegralQ[E^x/(x+1)^2] should be False");
}

/* ================= EXTENSIVE C1 STRESS TEST =================
 * A broad battery over the base-field rational exponential integral, all
 * empirically verified (closes + exact diff-back).  Each group varies one
 * structural axis of gamma = g E^f. */

/* f = a x + b (linear exponent): single- and multi-pole denominators, repeated
 * poles, improper cofactors, and scaled/negated exponents. */
static void test_stress_linear_exp(void) {
    /* single simple pole -> one ei */
    assert_ei("E^x/(x - 2)");
    assert_ei("E^x/(x + 5)");
    /* repeated poles (den(y) NOT squarefree: multiplicity m -> y-pole m-1) */
    assert_ei("E^x/(x + 1)^2");
    assert_ei("E^x/(x + 1)^3");
    assert_ei("E^x/(x + 3)^4");
    assert_ei("E^x/(x - 1)^3");
    /* alpha = 0 (pole at the origin) at several multiplicities */
    assert_ei("E^x/x^2");
    assert_ei("E^x/x^3");
    assert_ei("E^x/x^4");
    /* several distinct rational poles -> several ei terms */
    assert_ei("E^x/((x + 1)(x + 2))");
    assert_ei("E^x/((x - 1)(x + 2)(x + 3))");
    assert_ei("E^x/((x + 1)(x + 2)(x + 3)(x + 4))");
    assert_ei("E^x/((x + 1)^2 (x + 2))");
    assert_ei("E^x/((x + 1)^2 (x + 2)^2)");
    /* quadratic denominators that split over Q */
    assert_ei("E^x/(x^2 - 1)");
    assert_ei("E^x/(x^2 - 4)");
    assert_ei("E^x/(x^2 - 5 x + 6)");
    /* improper cofactor -> polynomial part + ei */
    assert_ei("(x^3 + 1) E^x/(x^2 + 3 x + 2)");
    assert_ei("(2 x - 3) E^x/((x - 1)(x - 4))");
    /* scaled / negated exponents */
    assert_ei("E^(-x)/(x - 1)");
    assert_ei("E^(-x)/(x + 1)^2");
    assert_ei("E^(-x)/((x + 1)(x - 1))");
    assert_ei("E^(2 x)/(x + 1)^2");
    assert_ei("E^(2 x)/(x + 1)^3");
    assert_ei("E^(3 x)/(x - 1)");
    assert_ei("E^(x/2)/(2 x + 1)");
}

/* f = 1/x, 1/x^2 (reciprocal exponents): the P2 q-side term and pole-at-origin
 * structure. */
static void test_stress_reciprocal_exp(void) {
    assert_ei("E^(1/x)");            /* Cherry Ex 5.1 */
    assert_ei("E^(1/x)/x");
    assert_ei("E^(1/x)/(x (x + 1))");
    assert_ei("E^(1/x^2)/x");
    assert_ei("x E^(1/x)");
    assert_ei("x^2 E^(1/x)");
    /* ei cancels to elementary at these higher reciprocal powers */
    assert_closes("E^(1/x)/x^3");
    assert_closes("E^(1/x^2)/x^3");
}

/* f = x^2, x^3 (polynomial exponents of degree > 1). */
static void test_stress_poly_exp(void) {
    assert_ei("E^(x^2)/x");
    assert_ei("E^(x^2)/x^3");
    assert_ei("E^(x^3)/x");
}

/* ei terms cancel to a purely ELEMENTARY closed form (still must close + diff0). */
static void test_stress_elementary_cancel(void) {
    assert_closes("x E^x/(x + 1)^2");   /* -> E^x/(x+1) */
    assert_closes("E^(1/x)/x^2");       /* -> -E^(1/x)  */
    assert_closes("(x + 1) E^(1/x)/x^4");
    /* pure-elementary siblings (polynomial cofactor / exact derivative) */
    assert_closes("x^2 E^(x^3)");       /* -> E^(x^3)/3 */
    assert_closes("(-2/x^3) E^(1/x^2)");/* -> E^(1/x^2)  */
    assert_closes("(x^3 + 2 x) E^x");
    assert_closes("(x^5 + x) E^x");
}

/* Cherry 1989 erf (completing-square, §3): q = s^2 admits error functions.
 * (Cherry's erf(u) = INT u' e^(u^2) du = (Sqrt[Pi]/2) Erfi(u), so answers carry a
 * Sqrt[Pi] factor against Mathilda's Erfi normalization.) */
static void test_cherry_erf(void) {
    /* Ex 5.2: INT (1/x + 1/x^2) e^(1/x^2) dx = -1/2 ei(1/x^2) - erf(1/x). */
    assert_erf("(1/x + 1/x^2) E^(1/x^2)");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[(1/x+1/x^2) E^(1/x^2), x, Method -> \"RischTranscendental\"]"
        " - (-1/2 ExpIntegralEi[1/x^2] - (Sqrt[Pi]/2) Erfi[1/x])]", "0"),
        "Ex 5.2 exact form");
    /* Pure error function (no ei part). */
    assert_erf("E^(1/x^2)");
    /* Two-erf, completing-square with a concrete algebraic constant (a = 4 -> +-2). */
    assert_erf("E^((x^4 + 4)/x^2)");
    assert_erf("E^((x^4 + 9)/x^2)");
    assert_erf("E^((4 x^4 + 1)/x^2)");
    assert_erf("(x^4 - 1)/x^4 E^(1/x^2)");
}

/* Cherry Ex 5.4: two error functions with a SYMBOLIC parameter a (constants
 * quadratic over Q(a), i.e. +-2 Sqrt[a]).  The coefficient solve keeps a a
 * parameter (Solve over the ansatz unknowns, not SolveAlways). */
static void test_cherry_erf_symbolic(void) {
    assert_erf("E^((x^4 + a)/x^2)");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[E^((x^4+a)/x^2), x, Method -> \"RischTranscendental\"]"
        " - (Sqrt[Pi]/4)(E^(-2 Sqrt[a]) Erfi[(x^2+Sqrt[a])/x]"
        "                + E^(2 Sqrt[a]) Erfi[(x^2-Sqrt[a])/x])]", "0"),
        "Ex 5.4 exact form (symbolic a)");
}

/* Extensive C2 stress: more real-algebraic ei constants and erf cases. */
static void test_stress_c2(void) {
    /* algebraic ei over Q(sqrt d): assorted radicands, scaled/negated exponents,
     * non-monic and improper cofactors. */
    assert_ei("E^x/(x^2 - 6)");
    assert_ei("E^x/(x^2 - 7)");
    assert_ei("E^x/(x^2 - 10)");
    assert_ei("E^x/(x^2 - 3 x + 1)");     /* roots (3 +- sqrt5)/2 */
    assert_ei("E^(-x)/(x^2 - 3)");
    assert_ei("(x + 1) E^x/(x^2 - 3)");
    assert_ei("E^x/(2 x^2 - 3)");         /* non-monic denominator */
    /* erf: multi-term and symbolic-parameter completing-square. */
    assert_erf("E^((9 x^4 + 16)/x^2)");
    assert_erf("E^((x^4 + b^2)/x^2)");    /* symbolic b -> beta = +-2 b */
    assert_erf("(x^4 - 1)/x^4 E^(1/x^2)");

    /* KNOWN GAPS (shared-function limits, NOT Cherry logic) — currently decline
     * cleanly; pinned so a future fix to the algebraic normaliser trips here:
     *  - Together/Solve over Q(sqrt d) does not reduce conjugate linear factors
     *    (x - 1 - sqrt2)(x - 1 + sqrt2) back to x^2 - 2 x - 1, over-inflating the
     *    linear system (a rational-quadratic denominator with irrational roots);
     *  - PolynomialSqrt fails on a polynomial with NUMERIC radical coefficients
     *    (x^4 + 2 sqrt2 x^2 + 2), though it succeeds for a symbolic radical. */
    assert_declines("E^x/(x^2 - 2 x - 1)");
    assert_declines("E^x/((x^2 - 2)(x + 1))");
    assert_declines("E^((x^4 + 2)/x^2)");
}

/* Decline-safety: algebraic / complex ei-argument constants (Cherry §7, a later
 * phase) MUST leave the integral unevaluated rather than emit a partial/wrong form. */
static void test_stress_declines(void) {
    /* q = x^3 is NOT a perfect square -> no erf; and no ei -> clean decline. */
    assert_declines("E^((x^6 + 8)/x^3)");
    /* e^(1/x^4): outside the ei/erf class (matching identity has no solution). */
    assert_declines("E^(1/x^4)");
    /* real-quadratic exponent with a NEGATIVE completing-square constant -> the
     * erf constant beta = +-2 I is complex; deferred, declines cleanly. */
    assert_declines("E^((x^4 - 1)/x^2)");
    /* A lone conjugate pair over Q(i)/Q(i sqrt3) now closes (see test_complex_ei);
     * what still defers is a complex pair MIXED with another pole — a P2/reciprocal
     * term.  E^(1/x)/(x^2+1) has a linear-exponent denominator q = x (degq != 0),
     * so the numeric-complex gate keeps the +-I pair out (mixing it with the P2
     * q-side term would blow up Together over Q(i)); clean decline. */
    assert_declines("E^(1/x)/(x^2 + 1)");
    /* Irreducible cubic: only the real root is admitted, so the (complex) rest
     * of the pole set is unmatched -> clean decline (not a wrong partial form). */
    assert_declines("E^x/(x^3 - 2)");
}

/* Cherry Thm 5.4 case b (flat multi-term exponential): f rational in a single
 * kernel E^w with SEVERAL commensurate Laurent terms Sum_i p_i E^(i w).  The
 * single-shape rt_cherry_ei declines (its cofactor would carry a residual
 * exponential); rt_cherry_exp_multiterm integrates term-by-term and sums. */
static void test_cherry_multiterm(void) {
    /* two essential terms, shared denominator */
    assert_ei("(E^x + E^(2 x))/(x - 1)");
    ASSERT_MSG(eval_is(
        "Simplify[Integrate[(E^x + E^(2 x))/(x - 1), x, Method -> \"RischTranscendental\"]"
        " - (E ExpIntegralEi[x - 1] + E^2 ExpIntegralEi[2 x - 2])]", "0"),
        "multi-term (E^x+E^(2x))/(x-1) exact form");
    /* product form E^x (E^x + 1) = E^(2x) + E^x over a quadratic denominator */
    assert_ei("E^x (E^x + 1)/((x - 1)(x - 2))");
    assert_ei("(E^x + E^(2 x))/((x - 1)(x - 2))");
    /* three commensurate terms */
    assert_ei("(E^x + E^(2 x) + E^(3 x))/(x - 1)");
    /* a term whose per-term integral needs a P1 resultant (quadratic pole) */
    assert_ei("(E^x + E^(2 x))/(x^2 - 2)");
    /* single-term integrands still go through rt_cherry_ei, unaffected */
    assert_ei("E^(2 x)/(x - 1)");
}

int main(void) {
    core_init();
    TEST(test_cherry_1989_ei);
    TEST(test_cherry_1986_macsyma);
    TEST(test_nonlinear_exponent);
    TEST(test_cherry_algebraic_ei);
    TEST(test_cherry_erf);
    TEST(test_cherry_erf_symbolic);
    TEST(test_complex_ei);
    TEST(test_const_offset);
    TEST(test_elementaryq_false);
    TEST(test_cherry_multiterm);
    /* extensive stress battery */
    TEST(test_stress_linear_exp);
    TEST(test_stress_reciprocal_exp);
    TEST(test_stress_poly_exp);
    TEST(test_stress_elementary_cancel);
    TEST(test_stress_c2);
    TEST(test_stress_declines);
    printf("All Cherry ExpIntegralEi tests passed.\n");
    return 0;
}
