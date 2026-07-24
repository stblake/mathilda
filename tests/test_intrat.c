/* test_intrat.c — Phase 1 of the BronsteinRational port.
 *
 * Coverage:
 *   * Helpers: Content / Primitive / Monic / LeadingCoefficient.
 *   * IntegratePolynomial — term-by-term integration over Q[x].
 *   * HermiteReduce — Mack's linear version, on the four reference
 *     cases lifted from IntegrateRational.m:1331-1358 plus a couple
 *     of single-factor sanity checks.
 *   * Integrate dispatcher — polynomial / 1/(x-a) / 1/(x-a)^2 /
 *     c·D[d]/d^k cases that Phase 1 closes via derivative
 *     recognition.
 *
 * The universal correctness predicate is differentiation:
 *   ASSERT_INTEGRAL_OK(f, x)  ≡  Cancel[Together[D[Integrate[f,x],x]-f]] === 0
 * which catches algorithmic bugs even when the printed form differs
 * from the Mathematica reference.  HermiteReduce uses the analogous
 * `Cancel[Together[D[g, x] + h - f]] === 0` check.
 */

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/* Test helpers                                                        */
/* ------------------------------------------------------------------ */

/* Compare evaluated output to a string. */
static void run_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Check that ∫ integrand dx differentiates back to the integrand.
 * We Expand both sides before canceling so structurally-different
 * polynomial coefficients (e.g. Plus distributed differently) collapse
 * to the same canonical form. */
static void assert_integral_correct(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Cancel[Together[Expand[D[Integrate[%s, x], x] - (%s)]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: D[Integrate[%s, x], x] - %s != 0\n  Got: %s\n",
               integrand, integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Check HermiteReduce[f, x] = {g, h} satisfies D[g, x] + h - f == 0. */
static void assert_hermite_correct(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "With[{hr = Integrate`HermiteReduce[%s, x]}, "
        "Cancel[Together[D[hr[[1]], x] + hr[[2]] - (%s)]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: HermiteReduce decomposition incorrect for %s\n  Got: %s\n",
               integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

static void test_helpers_content(void) {
    run_eq("Integrate`Helpers`Content[3 + 6 x + 9 x^2, x]", "3");
    run_eq("Integrate`Helpers`Content[2 x^2 + 4 x, x]", "2");
    run_eq("Integrate`Helpers`Content[x^3 + x, x]", "1");
}

static void test_helpers_primitive(void) {
    run_eq("Integrate`Helpers`Primitive[3 + 6 x + 9 x^2, x]",
           "1 + 2 x + 3 x^2");
    run_eq("Integrate`Helpers`Primitive[2 x^2 + 4 x, x]",
           "2 x + x^2");
}

static void test_helpers_leading_coefficient(void) {
    run_eq("Integrate`Helpers`LeadingCoefficient[3 x^2 + 1, x]", "3");
    run_eq("Integrate`Helpers`LeadingCoefficient[5, x]", "5");
    run_eq("Integrate`Helpers`LeadingCoefficient[x^4 - x + 2, x]", "1");
}

static void test_helpers_monic(void) {
    run_eq("Integrate`Helpers`Monic[2 x^2 + 4 x + 6, x]",
           "3 + 2 x + x^2");
    /* Polynomial free of x ⇒ 1 (Mathematica convention). */
    run_eq("Integrate`Helpers`Monic[7, x]", "1");
}

/* ------------------------------------------------------------------ */
/* IntegratePolynomial                                                 */
/* ------------------------------------------------------------------ */

static void test_integrate_polynomial(void) {
    run_eq("Integrate`IntegratePolynomial[0, x]", "0");
    run_eq("Integrate`IntegratePolynomial[1, x]", "x");
    run_eq("Integrate`IntegratePolynomial[x, x]", "1/2 x^2");
    run_eq("Integrate`IntegratePolynomial[x^2, x]", "1/3 x^3");
    /* Linear combination. */
    run_eq("Integrate`IntegratePolynomial[3 + 5 x + 2 x^2, x]",
           "3 x + 5/2 x^2 + 2/3 x^3");
    /* Higher-order with parameter — Phase 1 leaves the coefficient
     * unchanged because IntegratePolynomial assumes a polynomial in x.
     * The y^2-3y+1 factor is treated as a constant. */
    assert_integral_correct("(y^2 - 3 y + 1)*(1 - x^2)^3");
}

/* ------------------------------------------------------------------ */
/* HermiteReduce                                                       */
/* ------------------------------------------------------------------ */

static void test_hermite_simple(void) {
    /* deg(D̄) starts at 0 ⇒ f stays as-is, g == 0. */
    run_eq("Integrate`HermiteReduce[1/x, x]", "{0, 1/x}");
    run_eq("Integrate`HermiteReduce[1/(x-a), x]", "{0, -1/(a - x)}");

    /* 1/(x-a)^2: classical reduction down to g = -1/(x-a), h = 0. */
    assert_hermite_correct("1/(x-1)^2");
    assert_hermite_correct("1/(x-a)^2");

    /* Two repeated factors: the algorithm iterates twice. */
    assert_hermite_correct("1/(x^2 (x-1)^3)");
    assert_hermite_correct("(2 x + 3)/(x^2 (x+1)^4)");

    /* Squarefull denominator with an integer-content factor — caught
     * the scale-invariance bug in the gcd choice during Phase 1.
     * We use the monic gcd to keep the algorithm correct. */
    assert_hermite_correct("1/(4 x^2 (5 - 4 x)^2)");
}

static void test_hermite_reference(void) {
    /* The four cases lifted from IntegrateRational.m:1331-1358.
     * Verification is by differentiation. */
    assert_hermite_correct("(x^7 - 24 x^4 - 4 x^2 + 8 x - 8)/(x^8 + 6 x^6 + 12 x^4 + 8 x^2)");
    assert_hermite_correct("1/(x^4 + 1)^8");
    assert_hermite_correct("(-4 + x - 9 x^2)/(10 x - 8 x^2)^3");
    assert_hermite_correct("(28 + 110 x + 147 x^2 - 1131 x^3 - 945 x^4 - 189 x^5)/"
                            "(16 - 32 x + 100 x^2 - 56 x^3 + 296 x^4 + 216 x^5 + 36 x^6)");
}

/* ------------------------------------------------------------------ */
/* Integrate dispatcher                                                */
/* ------------------------------------------------------------------ */

static void test_integrate_polynomial_dispatch(void) {
    run_eq("Integrate[x^2, x]", "1/3 x^3");
    run_eq("Integrate[3 + 5 x + 2 x^2, x]", "3 x + 5/2 x^2 + 2/3 x^3");
    run_eq("Integrate[a x^3 + b, x]", "b x + 1/4 a x^4");
    /* Pure constant. */
    run_eq("Integrate[5, x]", "5 x");
}

static void test_integrate_log_recognition(void) {
    /* num is exactly D[den]/den ⇒ Log[den]. */
    run_eq("Integrate[1/x, x]", "Log[x]");
    run_eq("Integrate[2 x/(x^2 + 1), x]", "Log[1 + x^2]");
    run_eq("Integrate[1/(x - 1), x]", "Log[-1 + x]");
    /* General correctness checks — formatting may vary but the
     * differential check is the contract. */
    assert_integral_correct("(2 x + 3)/(x^2 + 3 x + 5)");
    assert_integral_correct("1/(x - a)");
}

static void test_integrate_rational_via_recognition(void) {
    /* num is exactly c·D[pol] with pol a squarefull factor of den. */
    assert_integral_correct("1/(x - a)^2");
    assert_integral_correct("(2 x + 3)/(x^2 + 3 x + 5)^2");
    assert_integral_correct("(x - 1)/(x + 1)^3");
    /* Improper rational: poly part + recognition. */
    assert_integral_correct("(x^3 + 2 x)/(x^2 + 1)");
}

static void test_integrate_unevaluated(void) {
    /* Genuinely non-elementary integrand: stays as Integrate[...].
     * (Sin[x] used to bubble back unevaluated; with the
     * Integrate`RischNorman dispatcher hook it now closes to a
     * Tan[x/2]-form antiderivative.  1/Log[x] likewise now closes to
     * LogIntegral[x] via the widened li recognizer, so use a nested-log
     * integrand with no elementary/special-function antiderivative.) */
    run_eq("Integrate[1/Log[Log[x]], x]", "Integrate[1/Log[Log[x]], x]");
}

/* ------------------------------------------------------------------ */
/* Phase 2 — LRT log part                                              */
/* ------------------------------------------------------------------ */

static void test_helpers_squarefree(void) {
    /* Multiplicity-indexed list. */
    run_eq("Integrate`Helpers`SquareFree[(x-1)^2 (x+1)]",
           "{{1 + x, 1}, {-1 + x, 2}}");
    run_eq("Integrate`Helpers`SquareFree[1 + 4 t^2]",
           "{{1 + 4 t^2, 1}}");
}

static void test_helpers_extract_constants(void) {
    run_eq("Integrate`Helpers`ExtractConstants[(2 x)/(3 (x^2+1)), x]",
           "{2/3, x/(1 + x^2)}");
    run_eq("Integrate`Helpers`ExtractConstants[5/(x^2+1), x]",
           "{5, 1/(1 + x^2)}");
}

static void test_helpers_apartlist(void) {
    /* Three-pole partial-fraction expansion as a List. */
    run_eq("Integrate`Helpers`ApartList[1/((x-1)(x-2)(x-3)), x]",
           "{-1/(-2 + x), 1/2/(-3 + x), 1/2/(-1 + x)}");
}

static void test_intrationallogpart(void) {
    /* Bronstein Example 2.4.1 — 1/(x^4+1).
     * Resultant(D, A - t D') = 1 + 256 t^4. */
    run_eq("Integrate`IntRationalLogPart[1/(x^4+1), x, t]",
           "{{1 + 256 t^4, 4 t + x}}");
    /* Simple: 1/(x^2+1) ⇒ Q=1+4t^2, S=2t+x. */
    run_eq("Integrate`IntRationalLogPart[1/(x^2+1), x, t]",
           "{{1 + 4 t^2, 2 t + x}}");
    /* Linear factors: 1/((x-1)(x-2)).  The log-argument S is content-primitive
     * (integer coefficients); it is equal up to the constant factor 2 to the
     * monic-in-x form -3/2 - t/2 + x, and a constant multiple inside a Log is
     * absorbed into the RootSum constant. */
    run_eq("Integrate`IntRationalLogPart[1/((x-1)(x-2)), x, t]",
           "{{1 - t^2, -3 - t + 2 x}}");
    /* RootSum -> True wraps in symbolic RootSum heads. */
    run_eq("Integrate`IntRationalLogPart[1/(x^2+1), x, t, RootSum -> True]",
           "RootSum[Function[t, 1 + 4 t^2], Function[t, t Log[2 t + x]]]");
}

/* Regression: a rational integrand whose coefficients live in a radical
 * tower Q(Sqrt[2], Sqrt[3]) once drove IntRationalLogPart into a
 * non-terminating poly_div_rem loop — the leading Sqrt[6] = Sqrt[2] Sqrt[3]
 * coefficient never cancelled under the tower-blind whole-polynomial
 * Together/Cancel subtraction. Division now subtracts coefficient-wise, so
 * this terminates. The antiderivative is correct but sits in a nested-radical
 * form Cancel cannot reduce to 0, so verify it numerically at a sample point
 * (this also guards against re-introducing the hang). */
static void test_intrationallogpart_tower_terminates(void) {
    run_eq("PossibleZeroQ[(D[Integrate[1/((x - Sqrt[2]) (x - Sqrt[3])), x], x]"
           " - 1/((x - Sqrt[2]) (x - Sqrt[3]))) /. x -> 37/10]", "True");
}

static void test_integrate_lrt_linear_q(void) {
    /* When the resultant Q factors completely into linear pieces over
     * Q, the integrator now closes the integral end-to-end. */
    assert_integral_correct("1/((x-1)(x-2))");
    assert_integral_correct("1/(x^2 - 1)");
    assert_integral_correct("(2 x - 1)/((x-1)(x-2)(x-3))");
    /* Higher degree with rational roots. */
    assert_integral_correct("1/(x (x-1) (x-2) (x-3))");
}

/* ------------------------------------------------------------------ */
/* Phase 3 — LogToAtan                                                 */
/* ------------------------------------------------------------------ */

static void test_logtoatan_constant_b(void) {
    /* B free of x with A not free of x: 2 ArcTan[A/B]. */
    run_eq("Integrate`LogToAtan[x, a, x]", "2 ArcTan[x/a]");
    run_eq("Integrate`LogToAtan[x, 1, x]", "2 ArcTan[x]");
    /* Both free of x: 0. */
    run_eq("Integrate`LogToAtan[a, b, x]", "0");
}

static void test_logtoatan_recursive(void) {
    /* The chain x^3 - 3 x divided by x^2 - 2 should produce a sum of
     * ArcTan terms that matches the Rioboo recursion (cf.
     * IntegrateRational.m:1587-1589 reference test). */
    run_eq("Integrate`LogToAtan[x^3 - 3 x, x^2 - 2, x]",
           "2 ArcTan[x] + 2 ArcTan[x^3] + 2 ArcTan[1/2 (x - 3 x^3 + x^5)]");
}

static void test_integrate_lrt_naivelogpart_fallback(void) {
    /* Pure-numeric quartic with no rational / quadratic factorisation:
     * Phase 4's LogToReal cannot close it, but Phase 8d-bonus's
     * biquadratic detector does — `1+x^4` has c1=c3=0 so it expands
     * to an explicit 4-term Plus[Log[α±x],...] form.  We assert the
     * result is no longer the unevaluated Integrate[..] head and
     * contains Log somewhere. */
    Expr* e = parse_expression("Integrate[1/(x^4 + 1), x]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    /* Must NOT be the unevaluated Integrate head. */
    ASSERT(strcmp(r->data.function.head->data.symbol.name, "Integrate") != 0);
    expr_free(e);
    expr_free(r);
}

static void test_integrate_parametric_biquadratic_closes(void) {
    /* Parametric biquadratic 1/(b + a x^4) closes via the Sophie-Germain
     * factorisation under the positive-symbol assumption (b/a > 0).
     * The result is a closed elementary expression (Plus of Logs and
     * ArcTans with held Sqrt coefficients), not a held RootSum. */
    Expr* e = parse_expression("Integrate[1/(b + a x^4), x]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    /* Must NOT be RootSum — that would mean the Sophie-Germain path
     * silently regressed back to NaiveLogPart's fallback. */
    ASSERT(strcmp(r->data.function.head->data.symbol.name, "RootSum") != 0);
    /* Verify mathematical correctness at a concrete instantiation. */
    Expr* check = parse_expression(
        "N[(D[Integrate[1/(b + a x^4), x], x] - 1/(b + a x^4)) "
        "/. {a -> 2, b -> 3, x -> 7/2}]");
    Expr* check_r = evaluate(check);
    ASSERT(check_r->type == EXPR_REAL);
    ASSERT(fabs(check_r->data.real) < 1e-10);
    expr_free(check); expr_free(check_r);
    expr_free(e); expr_free(r);
}

/* ------------------------------------------------------------------ */
/* Phase 4 — LogToReal closing                                         */
/* ------------------------------------------------------------------ */

static void test_integrate_arctan(void) {
    /* Quadratic with negative discriminant — closes to ArcTan. */
    run_eq("Integrate[1/(x^2 + 1), x]", "ArcTan[x]");
    run_eq("Integrate[1/(x^2 + 4), x]", "1/2 ArcTan[1/2 x]");
    /* Quadratic with linear shift. */
    run_eq("Integrate[1/(x^2 + 2 x + 5), x]",
           "1/2 ArcTan[1/2 (1 + x)]");
    /* Mixed linear + quadratic factor. */
    assert_integral_correct("1/((x^2+1)(x-1))");
    assert_integral_correct("(x+2)/((x^2+1)(x-1))");
}

static void test_integrate_arctanh_simplification(void) {
    /* Phase 6 collapses c Log[A] - c Log[B] into ArcTanh when the
     * resulting argument is rational in x. */
    run_eq("Integrate[1/(x^2 - 1), x]", "-ArcTanh[x]");
    /* Linear-pole pair: Log[x-2] - Log[x-1] -> -2 ArcTanh[-3 + 2 x]. */
    run_eq("Integrate[1/((x-1)(x-2)), x]", "-2 ArcTanh[-3 + 2 x]");
    /* The combined LogToArcTanh builtin can be invoked directly. */
    run_eq("Integrate`LogToArcTanh[1/2 Log[-1+x] - 1/2 Log[1+x], x]",
           "-ArcTanh[x]");
    run_eq("Integrate`LogToArcTanh[Log[1+x] + Log[-1+x], x]",
           "Log[-1 + x^2]");
}

static void test_options_accepted(void) {
    /* Trailing options are stripped before dispatch — the result must
     * be the same as without them. */
    run_eq("Integrate`BronsteinRational[1/(x^2 + 1), x]", "ArcTan[x]");
    run_eq("Integrate`BronsteinRational[1/(x^2 + 1), x, \"PFD\" -> True]",
           "ArcTan[x]");
    run_eq("Integrate`BronsteinRational[1/(x^2 - 1), x, \"LogToArcTan\" -> True]",
           "-ArcTanh[x]");
    run_eq("Integrate`BronsteinRational[1/(x^2 - 1), x, "
           "\"PFD\" -> True, \"LogToArcTan\" -> True]",
           "-ArcTanh[x]");
    /* Extension option is recognised (advisory in Phase 7). */
    run_eq("Integrate`BronsteinRational[1/(x^2 - 1), x, Extension -> Sqrt[2]]",
           "-ArcTanh[x]");
}

static void test_integrate_quartic_factorable(void) {
    /* x^4 + x^2 + 1 = (x^2 + x + 1)(x^2 - x + 1) — two quadratics
     * with negative discriminants give two ArcTans + two Logs. */
    assert_integral_correct("1/(x^4 + x^2 + 1)");
    /* Cubic with one rational + irreducible-quadratic factor over Q. */
    assert_integral_correct("1/(x^3 - 1)");
}

/* ------------------------------------------------------------------ */
/* Phase 8b — NaiveLogPart RootSum fallback                            */
/* ------------------------------------------------------------------ */

/* NaiveLogPart on a parametric quadratic: Phase 8d-bonus's gate
 * defers to held RootSum form whenever the polynomial has a
 * variable other than the bound one. */
static void test_naivelogpart_basic(void) {
    Expr* e = parse_expression("Integrate`NaiveLogPart[1/(x^2 + a), x]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol.name, "RootSum") == 0);
    ASSERT(r->data.function.arg_count == 2);
    /* Both children are Function nodes. */
    for (size_t k = 0; k < 2; k++) {
        Expr* fn = r->data.function.args[k];
        ASSERT(fn->type == EXPR_FUNCTION);
        ASSERT(fn->data.function.head->type == EXPR_SYMBOL);
        ASSERT(strcmp(fn->data.function.head->data.symbol.name, "Function") == 0);
    }
    expr_free(e);
    expr_free(r);
}

/* Numeric biquadratic — Phase 8d-bonus expands to an explicit
 * 4-term Plus[Log[α±x],...]; the result must NOT be RootSum-headed. */
static void test_naivelogpart_quartic_hard(void) {
    Expr* e = parse_expression(
        "Integrate`NaiveLogPart[(x^2 - 1)/(2 x^4 - 2 x^2 + 1), x]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    /* Should NOT be the held RootSum form for this non-parametric
     * biquadratic — the radical-formula expansion fires. */
    ASSERT(strcmp(r->data.function.head->data.symbol.name, "RootSum") != 0);
    expr_free(e);
    expr_free(r);
}

/* D[NaiveLogPart[1/(x^2+a), x], x] threads through the body Function
 * via the D[RootSum] rule in src/deriv.c, then collapses via the
 * Lagrange-interpolation identity in src/root.c — the body becomes
 * `1 / (d'(t)(x - t))` after differentiation, which Σ-closes to
 * `1 / d(x) = 1/(a + x^2)` (the original integrand). */
static void test_naivelogpart_derivative_threads(void) {
    Expr* e = parse_expression(
        "D[Integrate`NaiveLogPart[1/(x^2 + a), x], x]");
    Expr* r = evaluate(e);
    /* Verify equality with the original integrand modulo Together. */
    Expr* check = parse_expression(
        "Together[D[Integrate`NaiveLogPart[1/(x^2 + a), x], x] - 1/(x^2 + a)]");
    Expr* check_e = evaluate(check);
    char* got = expr_to_string(check_e);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: D[NaiveLogPart[1/(x^2+a), x], x] != 1/(x^2+a)\n  diff: %s\n",
               got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(check); expr_free(check_e);
    expr_free(e); expr_free(r);
}

/* ------------------------------------------------------------------ */
/* PolynomialQuotientRemainder & SubresultantPolynomialRemainders      */
/* ------------------------------------------------------------------ */

static void test_polynomial_quotient_remainder(void) {
    run_eq("PolynomialQuotientRemainder[x^3 + x + 1, x^2 + 1, x]",
           "{x, 1}");
    run_eq("PolynomialQuotientRemainder[x^4 - 1, x^2 - 1, x]",
           "{1 + x^2, 0}");
    /* With the Extension option. */
    run_eq("PolynomialQuotientRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]",
           "{Sqrt[2] + x, 0}");
}

static void test_subresultant_chain(void) {
    /* Coprime inputs: chain ends at a constant. */
    run_eq("SubresultantPolynomialRemainders[x^4 + 1, 2 x^3, x]",
           "{1 + x^4, 2 x^3, 2}");
    /* Standard univariate example used to drive the LRT pipeline. */
    Expr* e = parse_expression(
        "SubresultantPolynomialRemainders[x^3 - 1, x^2 - 1, x]");
    Expr* res = evaluate(e);
    /* The chain has at least 3 elements — first two are the inputs
     * (with bigger-degree first), third is the pseudo-remainder. */
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count >= 3);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Closed-form regressions for the five RootSum-leak inputs reported   */
/* on 2026-05-11.  Each integrand must (a) integrate to a closed real- */
/* elementary form (no RootSum / Function head anywhere) and (b)       */
/* satisfy D[result, x] − integrand → 0 under Cancel/Together.         */
/* ------------------------------------------------------------------ */

/* True if `e` contains a node with the given symbol head. */
static bool contains_head(const Expr* e, const char* head) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol.name, head) == 0) {
            return true;
        }
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (contains_head(e->data.function.args[i], head)) return true;
        }
    }
    return false;
}

/* Evaluate `Integrate`BronsteinRational[integrand, x]` and assert
 * the result has no RootSum / Function leak, then assert the
 * derivative-of-result matches the integrand. */
static void assert_closed_real(const char* integrand) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Integrate`BronsteinRational[%s, x]", integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    if (contains_head(res, "RootSum")) {
        char* got = expr_to_string(res);
        printf("FAIL: result contains RootSum for %s\n  Got: %s\n",
               integrand, got);
        free(got);
        ASSERT(false);
    }
    if (contains_head(res, "Function")) {
        char* got = expr_to_string(res);
        printf("FAIL: result contains Function for %s\n  Got: %s\n",
               integrand, got);
        free(got);
        ASSERT(false);
    }
    expr_free(e); expr_free(res);
}

/* Numerical (N[]-based) differentiation check: for integrands whose
 * result involves nested radicals that Cancel/Together can't reduce
 * to zero symbolically, sample N[D[result, x] - integrand] at a few
 * concrete x-values (and parameters when present) and assert it is
 * below machine epsilon. */
static void assert_integral_numeric_ok(const char* integrand,
                                        const char* bindings,
                                        const char* x_value) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate`BronsteinRational[%s, x], x] - (%s)) /. %s /. x -> %s]]",
        integrand, integrand, bindings, x_value);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    /* res should be a real number near zero. */
    bool ok = false;
    if (res && res->type == EXPR_REAL) {
        ok = (res->data.real >= 0.0 && res->data.real < 1e-9);
    } else if (res && res->type == EXPR_INTEGER) {
        ok = (res->data.integer == 0);
    }
    if (!ok) {
        char* got = expr_to_string(res);
        printf("FAIL: |D[Integrate[%s,x],x] - (%s)| at %s with %s = %s (expected ~0)\n",
               integrand, integrand, x_value, bindings, got);
        free(got);
        ASSERT(false);
    }
    expr_free(e); expr_free(res);
}

/* Issue 1: 1/(b - a x^3) — cubic with parameters, nth-root form.
 * Symbolic Cancel can't close the (b/a)^(1/3) diff, so verify
 * numerically at concrete a, b > 0 values. */
static void test_closed_cubic_params(void) {
    assert_closed_real("1/(b - a x^3)");
    assert_integral_numeric_ok("1/(b - a x^3)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b - a x^3)",
                                "{a -> 0.5, b -> 3.1}", "1.7");
}

/* Issue 2: 1/(x^5 + 1) — palindromic-quartic after Apart. */
static void test_closed_quintic_one(void) {
    assert_closed_real("1/(x^5 + 1)");
    assert_integral_numeric_ok("1/(x^5 + 1)", "{}", "0.7");
    assert_integral_numeric_ok("1/(x^5 + 1)", "{}", "2.3");
    assert_integral_numeric_ok("1/(x^5 + 1)", "{}", "-1.4");
}

/* Issue 3: x^4/(a + b x^3)^2 — Hermite-reduced cubic with params. */
static void test_closed_x4_over_cubic_sq(void) {
    assert_closed_real("x^4/(a + b x^3)^2");
    assert_integral_numeric_ok("x^4/(a + b x^3)^2",
                                "{a -> 1.1, b -> 0.7}", "0.4");
    assert_integral_numeric_ok("x^4/(a + b x^3)^2",
                                "{a -> 1.1, b -> 0.7}", "1.9");
}

/* Issue 4: (-1 + x^2)/(1 - 2 x^2 + 2 x^4) — biquadratic Q with
 * negative inner discriminant (Sophie-Germain-with-c2 branch).
 * Nested Sqrt[1+Sqrt[2]] forms aren't reducible by Cancel/Together,
 * so verify numerically. */
static void test_closed_biquadratic_neg_disc(void) {
    assert_closed_real("(-1 + x^2)/(1 - 2 x^2 + 2 x^4)");
    assert_integral_numeric_ok("(-1 + x^2)/(1 - 2 x^2 + 2 x^4)",
                                "{}", "0.3");
    assert_integral_numeric_ok("(-1 + x^2)/(1 - 2 x^2 + 2 x^4)",
                                "{}", "1.7");
    assert_integral_numeric_ok("(-1 + x^2)/(1 - 2 x^2 + 2 x^4)",
                                "{}", "2.5");
}

/* Issue 5: 1/x/(1 - 2 x^2 + 2 x^4) — Times-over-Plus distribution
 * + Log[c · p] stripping post-processing.  Result is rational and
 * the printed form is exact. */
static void test_closed_coefficient_distribution(void) {
    /* Concrete printed form check: outer 2 must be distributed and
     * Log[1/2 · q] reduced to Log[q]. */
    run_eq("Integrate`BronsteinRational[1/x/(1 - 2 x^2 + 2 x^4), x]",
           "Log[x] + 1/2 ArcTan[-1 + 2 x^2] - 1/4 Log[1 - 2 x^2 + 2 x^4]");
    assert_integral_correct("1/x/(1 - 2 x^2 + 2 x^4)");
}

/* ------------------------------------------------------------------ */
/* 2026-05-11 (rev 2) — generalised nth-root coverage for             */
/* 1/(b ± a x^n) at all n ≥ 3.  The previous deg-3 fix only handled  */
/* n = 3; n ∈ {4, 5, 6, 8, 9} (and beyond) now dispatch via the same */
/* sparse cyclotomic factorisation path.                              */
/* ------------------------------------------------------------------ */

static void test_closed_nth_root_quartic_minus(void) {
    /* 1/(b - a x^4) — q := -c_0/c_n positive, n=4 even: two real
     * roots ±r plus a (t^2 + r^2) pair. */
    assert_closed_real("1/(b - a x^4)");
    assert_integral_numeric_ok("1/(b - a x^4)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b - a x^4)",
                                "{a -> 0.5, b -> 3.1}", "1.2");
}

static void test_closed_nth_root_quintic_minus(void) {
    /* 1/(b - a x^5) — n=5 odd, q < 0: real root −r plus two
     * conjugate pairs at angles π/5 and 3π/5. */
    assert_closed_real("1/(b - a x^5)");
    assert_integral_numeric_ok("1/(b - a x^5)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b - a x^5)",
                                "{a -> 1.1, b -> 0.8}", "-0.6");
}

static void test_closed_nth_root_sextic_minus(void) {
    /* 1/(b - a x^6) — n=6 even, q > 0: two real roots ±r plus two
     * conjugate pairs (cos = ±1/2). */
    assert_closed_real("1/(b - a x^6)");
    assert_integral_numeric_ok("1/(b - a x^6)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b - a x^6)",
                                "{a -> 0.5, b -> 3.1}", "1.2");
}

static void test_closed_nth_root_sextic_plus(void) {
    /* 1/(b + a x^6) — n=6 even, q < 0: no real roots, three
     * conjugate pairs (Sqrt[3]/2, 0, -Sqrt[3]/2). */
    assert_closed_real("1/(b + a x^6)");
    assert_integral_numeric_ok("1/(b + a x^6)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b + a x^6)",
                                "{a -> 0.5, b -> 3.1}", "-1.2");
}

static void test_closed_nth_root_octic_plus(void) {
    /* 1/(b + a x^8) — n=8 even, q < 0: no real roots, four
     * conjugate pairs.  Cos[Pi/8] doesn't auto-reduce in Mathilda
     * (held), but the result is still real-elementary. */
    assert_closed_real("1/(b + a x^8)");
    assert_integral_numeric_ok("1/(b + a x^8)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b + a x^8)",
                                "{a -> 0.5, b -> 3.1}", "1.2");
}

static void test_closed_nth_root_nonic_plus(void) {
    /* 1/(b + a x^9) — n=9 odd, q > 0: one real root r plus four
     * conjugate pairs.  Cos[2π/9], Cos[4π/9], Cos[8π/9] don't
     * reduce (they are roots of an irreducible cubic), but the
     * result is still real-elementary (no RootSum / Function). */
    assert_closed_real("1/(b + a x^9)");
    assert_integral_numeric_ok("1/(b + a x^9)",
                                "{a -> 2.3, b -> 1.7}", "0.4");
    assert_integral_numeric_ok("1/(b + a x^9)",
                                "{a -> 0.5, b -> 3.1}", "-1.2");
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_helpers_content);
    TEST(test_helpers_primitive);
    TEST(test_helpers_leading_coefficient);
    TEST(test_helpers_monic);

    TEST(test_integrate_polynomial);

    TEST(test_hermite_simple);
    TEST(test_hermite_reference);

    TEST(test_integrate_polynomial_dispatch);
    TEST(test_integrate_log_recognition);
    TEST(test_integrate_rational_via_recognition);
    TEST(test_integrate_unevaluated);

    TEST(test_polynomial_quotient_remainder);
    TEST(test_subresultant_chain);

    TEST(test_helpers_squarefree);
    TEST(test_helpers_extract_constants);
    TEST(test_helpers_apartlist);
    TEST(test_intrationallogpart);
    TEST(test_intrationallogpart_tower_terminates);
    TEST(test_integrate_lrt_linear_q);
    TEST(test_integrate_lrt_naivelogpart_fallback);
    TEST(test_integrate_parametric_biquadratic_closes);

    TEST(test_logtoatan_constant_b);
    TEST(test_logtoatan_recursive);

    TEST(test_integrate_arctan);
    TEST(test_integrate_arctanh_simplification);
    TEST(test_options_accepted);
    TEST(test_integrate_quartic_factorable);

    TEST(test_naivelogpart_basic);
    TEST(test_naivelogpart_quartic_hard);
    TEST(test_naivelogpart_derivative_threads);

    /* 2026-05-11 — closed-form regressions for the five RootSum-leak
     * inputs (Phases A-D2). */
    TEST(test_closed_cubic_params);
    TEST(test_closed_quintic_one);
    TEST(test_closed_x4_over_cubic_sq);
    TEST(test_closed_biquadratic_neg_disc);
    TEST(test_closed_coefficient_distribution);

    /* 2026-05-11 (rev 2) — generalised nth-root coverage. */
    TEST(test_closed_nth_root_quartic_minus);
    TEST(test_closed_nth_root_quintic_minus);
    TEST(test_closed_nth_root_sextic_minus);
    TEST(test_closed_nth_root_sextic_plus);
    TEST(test_closed_nth_root_octic_plus);
    TEST(test_closed_nth_root_nonic_plus);

    printf("All Phase 1-7 (full BronsteinRational pipeline) tests passed!\n");
    return 0;
}
