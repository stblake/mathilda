/* test_integrate_beta.c
 *
 * Tests for definite integration by Euler-Beta reduction
 * (src/calculus/integrate_beta.c):
 *   Integrate`Beta[f, {x, 0, 1}]           x^(k-1)(1-x)^(l-1) -> Beta[k,l]
 *   Integrate`TrigPower[f, {x, 0, c}]      Sin^m Cos^n -> Beta reduction
 *   Method -> "Beta" / "TrigPower"         pinned mechanism (strict)
 *
 * Coverage: the non-elementary cases the FTC cannot reach (non-integer /
 * symbolic exponents, even powers over [0,Pi]), the Log[x]/Log[1-x] parameter
 * derivatives, the parity vanishing over [0,Pi]/[0,2Pi], the convergence-strip
 * gate (ConditionalExpression when undecided, NULL when refuted), the
 * auto-dispatch regression (integer cases still owned by Newton-Leibniz), and
 * the negative controls that MUST stay unevaluated.
 *
 * Values Mathilda leaves in an equivalent-but-unsimplified form are pinned via
 * Simplify[value - reference] == 0.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "message.h"

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
 * Beta on [0,1] -- explicit entry point.
 * ---------------------------------------------------------------------- */
static void test_beta_unit(void) {
    /* Integer exponents: the Beta reduction gives the same rational as FTC. */
    check_eq("Integrate`Beta[x^2 (1-x)^3, {x, 0, 1}]", "1/60");
    /* Half-integer: FTC cannot (incomplete Beta), Beta gives Pi/8. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2), {x, 0, 1}]", "1/8 Pi");
    /* Negative (but > -1) exponent of x: still convergent -> Pi/2. */
    check_eq("Integrate`Beta[x^(-1/2)(1-x)^(1/2), {x, 0, 1}]", "1/2 Pi");
    /* Symbolic exponents -> ConditionalExpression on the strip. */
    check_eq("Integrate`Beta[x^(k-1)(1-x)^(l-1), {x, 0, 1}, "
             "Assumptions -> k>0 && l>0]", "Beta[k, l]");
    check_eq("Integrate`Beta[x^(k-1)(1-x)^(l-1), {x, 0, 1}]",
             "ConditionalExpression[Beta[k, l], Re[k] > 0 && Re[l] > 0]");
}

static void test_beta_logweight(void) {
    /* Log[x] weight = d/dk Beta[k,l]. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2) Log[x], {x, 0, 1}]",
             "1/16 Pi (1 - 4 Log[2])");
    /* Log[1-x] weight = d/dl Beta[k,l]. */
    check_eq("Integrate`Beta[x^2 (1-x)^3 Log[1-x], {x, 0, 1}]", "-37/3600");
    /* Mixed Log[x] Log[1-x] weight = d^2/dk dl Beta. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2) Log[x] Log[1-x], {x, 0, 1}]",
             "1/48 Pi (9 - Pi^2 - 12 Log[2] + 24 Log[2]^2)");
}

static void test_beta_negative(void) {
    /* Divergent (k = -1 <= 0): the strip is refuted -> unevaluated. */
    check_eq("Integrate`Beta[x^(-2)(1-x)^(1/2), {x, 0, 1}]",
             "Integrate`Beta[Sqrt[1 - x]/x^2, {x, 0, 1}]");
    /* Not of the x^p (1-x)^q form. */
    check_eq("Integrate`Beta[1/(1+x), {x, 0, 1}]",
             "Integrate`Beta[1/(1 + x), {x, 0, 1}]");
    /* Wrong interval. */
    check_eq("Integrate`Beta[x^(1/2)(1-x)^(1/2), {x, 0, 2}]",
             "Integrate`Beta[Sqrt[1 - x] Sqrt[x], {x, 0, 2}]");
}

/* -------------------------------------------------------------------------
 * Sin^m Cos^n over the canonical trig intervals -- explicit entry point.
 * ---------------------------------------------------------------------- */
static void test_trigpower_unit(void) {
    /* [0, Pi/2] half-integer -> Gamma[3/4]^2/Sqrt[Pi] = (1/2) Beta[3/4,3/4]. */
    check_eq("Integrate`TrigPower[Sin[x]^(1/2) Cos[x]^(1/2), {x, 0, Pi/2}]",
             "Gamma[3/4]^2/Sqrt[Pi]");
    /* [0, Pi/2] with a constant factor. */
    check_eq("Integrate`TrigPower[3 Sin[x]^4 Cos[x]^2, {x, 0, Pi/2}]", "3/32 Pi");
    /* [0, Pi]: n even -> 2 x quarter. */
    check_eq("Integrate`TrigPower[Sin[x]^4 Cos[x]^2, {x, 0, Pi}]", "1/16 Pi");
    /* [0, Pi]: n even, m odd (m parity irrelevant) -> 4/15. */
    check_eq("Integrate`TrigPower[Sin[x]^3 Cos[x]^2, {x, 0, Pi}]", "4/15");
    /* [0, 2Pi]: both even -> 4 x quarter. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^2, {x, 0, 2 Pi}]", "1/4 Pi");
    /* Symbolic exponents over [0, Pi/2].  The affine strip (1+m)/2 > 0 is not
     * discharged by Simplify from m > 0 (no linear-implication reasoning); the
     * gate's Refine fallback proves it, so the value is bare. */
    check_eq("Integrate`TrigPower[Sin[x]^m Cos[x]^n, {x, 0, Pi/2}, "
             "Assumptions -> m>0 && n>0]",
             "1/2 Beta[1/2 (1 + m), 1/2 (1 + n)]");
    /* undecided without assumptions: still conditioned */
    check_eq("Integrate`TrigPower[Sin[x]^m Cos[x]^n, {x, 0, Pi/2}]",
             "ConditionalExpression[1/2 Beta[1/2 (1 + m), 1/2 (1 + n)], "
             "Re[1/2 (1 + m)] > 0 && Re[1/2 (1 + n)] > 0]");
    /* Over [0, Pi] only n parity matters, so m may stay symbolic (n = 2 even). */
    check_eq("Integrate`TrigPower[Sin[x]^m Cos[x]^2, {x, 0, Pi}, "
             "Assumptions -> m>0]",
             "Beta[1/2 (1 + m), 3/2]");
}

static void test_trigpower_vanishing(void) {
    /* [0, Pi]: n odd -> 0. */
    check_eq("Integrate`TrigPower[Sin[x]^4 Cos[x]^3, {x, 0, Pi}]", "0");
    /* [0, 2Pi]: m odd -> 0. */
    check_eq("Integrate`TrigPower[Sin[x]^5 Cos[x]^3, {x, 0, 2 Pi}]", "0");
    /* [0, 2Pi]: n odd -> 0. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^5, {x, 0, 2 Pi}]", "0");
}

static void test_trigpower_negative(void) {
    /* Symbolic n over [0, Pi]: parity undecidable -> unevaluated. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^n, {x, 0, Pi}]",
             "Integrate`TrigPower[Cos[x]^n Sin[x]^2, {x, 0, Pi}]");
    /* Non-canonical interval. */
    check_eq("Integrate`TrigPower[Sin[x]^2 Cos[x]^2, {x, 0, Pi/3}]",
             "Integrate`TrigPower[Cos[x]^2 Sin[x]^2, {x, 0, 1/3 Pi}]");
}

/* -------------------------------------------------------------------------
 * Auto-dispatch: new cases close; the integer cases stay owned by FTC.
 * ---------------------------------------------------------------------- */
static void test_dispatch_regression(void) {
    /* NEW: previously unevaluated, now closed by the Beta mechanism. */
    check_eq("Integrate[x^(1/2)(1-x)^(1/2), {x, 0, 1}]", "1/8 Pi");
    check_eq("Integrate[Sin[x]^(1/2) Cos[x]^(1/2), {x, 0, Pi/2}]",
             "Gamma[3/4]^2/Sqrt[Pi]");
    check_eq("Integrate[Sin[x]^4 Cos[x]^2, {x, 0, Pi}]", "1/16 Pi");
    /* Regression: integer-power cases still evaluate (via FTC), same value. */
    check_eq("Integrate[x^2 (1-x)^3, {x, 0, 1}]", "1/60");
    check_eq("Integrate[Sin[x]^4 Cos[x]^2, {x, 0, Pi/2}]", "1/32 Pi");
}

static void test_method_option(void) {
    check_eq("Integrate[x^(1/2)(1-x)^(1/2), {x, 0, 1}, Method -> \"Beta\"]",
             "1/8 Pi");
    check_eq("Integrate[Sin[x]^(1/2) Cos[x]^(1/2), {x, 0, Pi/2}, "
             "Method -> \"TrigPower\"]", "Gamma[3/4]^2/Sqrt[Pi]");
    /* The Euler extension (one extra linear factor) now closes 1/(1+x) on [0,1]
     * under the Beta mechanism -- (1+x)^(-1) is the a=b=1 Euler integrand, value
     * Log[2].  (Under Automatic, Newton-Leibniz computed the same value anyway.) */
    check_eq("Integrate[1/(1+x), {x, 0, 1}, Method -> \"Beta\"]", "Log[2]");
    /* Strict: a genuinely non-Beta/Euler integrand (a non-linear factor base)
     * stays unevaluated under Method -> "Beta". */
    check_eq("Integrate[1/(1+x^2), {x, 0, 1}, Method -> \"Beta\"]",
             "Integrate[1/(1 + x^2), {x, 0, 1}, Method -> \"Beta\"]");
}

/* -------------------------------------------------------------------------
 * Stress: high powers, deep half-integer, double log weight.
 * ---------------------------------------------------------------------- */
static void test_stress(void) {
    check_eq("Integrate`Beta[x^10 (1-x)^7, {x, 0, 1}]", "1/350064");
    check_eq("Integrate`Beta[x^(3/2)(1-x)^(5/2), {x, 0, 1}]", "3/256 Pi");
    check_eq("Integrate`TrigPower[Sin[x]^10 Cos[x]^8, {x, 0, Pi/2}]",
             "35/131072 Pi");
    /* Deep half-integer trig closes to a clean multiple of Pi. */
    check_eq("Simplify[Integrate`TrigPower[Sin[x]^10 Cos[x]^8, {x, 0, Pi/2}] "
             "- 1/2 Beta[11/2, 9/2]]", "0");
}


/* Evaluate `input`, assert its printed form is `expected` AND that no message
 * fired while it ran -- the direct guard against an internal probe's
 * diagnostic reaching an enclosing Check[] (which would then return its
 * failure branch for a successful computation). */
static void check_quiet(const char* input, const char* expected) {
    unsigned long fired0 = mth_msg_fired_count();
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
    ASSERT_STR_EQ(s, expected);
    if (mth_msg_fired_count() != fired0)
        fprintf(stderr, "FAIL: %s fired a message\n", input);
    ASSERT(mth_msg_fired_count() == fired0);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------------------------------------------------------------------------
 * Widened Beta strip: a Log[1-x]^j weight vanishes like x^j at 0 (Log[x]^i
 * like (1-x)^i at 1), so the integral converges for A > -j, B > -i and its
 * value is the removable limit of the Beta parameter derivative.  Before,
 * the boundary substitution Beta[0,1] PolyGamma[...] raised Infinity::indet
 * (seen by Check[]) and the stage declined.  References: the reference system
 * (Zeta[3], -Pi^2/6, -Pi^4/45, -Pi^4/180, Pi^2/3).
 * ---------------------------------------------------------------------- */
static void test_beta_widened_strip(void) {
    check_quiet("Integrate[Log[x] Log[1-x]/x, {x, 0, 1}]", "Zeta[3]");
    check_quiet("Integrate[Log[1-x]/x, {x, 0, 1}]", "-1/6 Pi^2");
    check_quiet("Integrate[Log[x]^2 Log[1-x]/x, {x, 0, 1}]", "-1/45 Pi^4");
    check_quiet("Integrate[Log[x] Log[1-x]^2/x, {x, 0, 1}]", "-1/180 Pi^4");
    check_quiet("Integrate[Log[x] Log[1-x]/(1-x), {x, 0, 1}]", "Zeta[3]");
    check_quiet("Integrate[Log[1-x]^2/x^2, {x, 0, 1}]", "1/3 Pi^2");
    check_quiet("Check[Integrate[Log[x] Log[1-x]/x, {x, 0, 1}], \"msg\"]", "Zeta[3]");
}

/* King's reflection x -> c - x exists only for a finite c.  On a half-line the
 * stage used to substitute x -> Infinity - x and raise Infinity::indet, which
 * made Check[] reject a Gaussian-Bessel integral another stage then closed. */
static void test_kingreflect_halfline_gate(void) {
    check_quiet("Integrate[x Exp[-p x^2] BesselJ[0, b x], {x, 0, Infinity}, "
                "Assumptions -> p > 0 && b > 0]", "1/(2 E^(b^2/(4 p)) p)");
    check_quiet("Check[Integrate[x Exp[-p x^2] BesselJ[0, b x], {x, 0, Infinity}, "
                "Assumptions -> p > 0 && b > 0], \"msg\"]", "1/(2 E^(b^2/(4 p)) p)");
}

/* -------------------------------------------------------------------------
 * Reductions into the Beta family (integrate_tanbeta_try,
 * integrate_invtrig_beta_try): Roach (2008) #9 and #7.  With a symbolic
 * exponent the antiderivative is non-elementary, and before these stages the
 * default cascade reached a closing stage only after Newton-Leibniz's
 * indefinite search (55-66 s, ParallelMixedTower 45 s of it): #9 then declined
 * (67 s), #7 closed after a further 125 s in DiffUnderInt.  Each value is
 * checked against NIntegrate at a generic parameter, and the call is timed.
 * ---------------------------------------------------------------------- */

/* Integrate[f, {x, lo, hi}, opts] returns `expected` within `secs` seconds,
 * and agrees with NIntegrate at the parameter point `pt`. */
static void check_timed_value(const char* f, const char* spec, const char* opts,
                              const char* expected, const char* pt, double secs) {
    char buf[2048];
    snprintf(buf, sizeof buf,
             "With[{t = AbsoluteTiming[Integrate[%s, %s%s]]}, "
             "{t[[1]] < %g, ToString[t[[2]]] === \"%s\", "
             "Abs[N[(t[[2]] /. ConditionalExpression[v_, _] :> v) /. %s, 20] - "
             "NIntegrate[(%s) /. %s, %s, WorkingPrecision -> 30]] < 10^-15}]",
             f, spec, opts, secs, expected, pt, f, pt, spec);
    check_eq(buf, "{True, True, True}");
}

static void test_tanbeta(void) {
    /* #9: t = Tan[x] gives (1 - t^2)^a on [0, 1], a Half-Beta, in Gamma form */
    check_timed_value("Cos[x]^(-2 a - 2) Cos[2 x]^a", "{x, 0, Pi/4}",
                      ", Assumptions -> a > 0",
                      "(Sqrt[Pi] Gamma[1 + a])/(2 Gamma[3/2 + a])", "a -> 2/3", 10);
    /* the strip undecided without assumptions */
    check_eq("Integrate[Cos[x]^(-2 a - 2) Cos[2 x]^a, {x, 0, Pi/4}]",
             "ConditionalExpression[(Sqrt[Pi] Gamma[1 + a])/(2 Gamma[3/2 + a]), "
             "Re[1 + a] > 0]");
    /* t^a: a plain Beta */
    check_timed_value("Sin[x]^a/Cos[x]^(a + 2)", "{x, 0, Pi/4}",
                      ", Assumptions -> a > 0", "1/(1 + a)", "a -> 5/7", 10);
    /* past Pi/4 the factor 1 - t^2 changes sign: declined (pinned, so the
     * decline is visible rather than followed by the rest of the cascade) */
    check_eq("Head[Integrate[Cos[x]^(-2 a - 2) Cos[2 x]^a, {x, 0, Pi/3}, "
             "Assumptions -> a > 0, Method -> \"Beta\"]]", "Integrate");
}

static void test_invtrig_beta(void) {
    /* #7: by parts, Pi/(4 n) minus a Half-Beta remainder */
    check_timed_value("ArcSin[x] x^(2 n - 1)", "{x, 0, 1}", ", Assumptions -> n > 0",
                      "(Pi - (Sqrt[Pi] Gamma[1/2 + n])/Gamma[1 + n])/(4 n)", "n -> 2/3", 10);
    check_timed_value("ArcCos[x] x^p", "{x, 0, 1}", ", Assumptions -> p > -1",
                      "(Sqrt[Pi] Gamma[1/2 (2 + p)])/((2 + 2 p) Gamma[1/2 (3 + p)])",
                      "p -> -1/3", 10);
    /* convergence (Re p > -1) undecided: conditioned; refuted: declined */
    check_eq("Head[Integrate[ArcSin[x] x^(2 n - 1), {x, 0, 1}]]", "ConditionalExpression");
    check_eq("Head[Integrate[ArcSin[x] x^p, {x, 0, 1}, Assumptions -> p < -2, "
             "Method -> \"Beta\"]]", "Integrate");
    /* a numeric exponent stays with Newton-Leibniz (elementary) */
    check_eq("Integrate[ArcSin[x] x^2, {x, 0, 1}]", "-2/9 + 1/6 Pi");
}

/* The strip of a Half-Beta is proved by Refine when Simplify leaves it
 * unreduced (1/2 + n > 0 under n > 0): no needless ConditionalExpression. */
static void test_beta_gate_refine(void) {
    check_eq("Integrate[x^(2 n)/Sqrt[1 - x^2], {x, 0, 1}, Assumptions -> n > 0]",
             "1/2 Beta[1/2 + n, 1/2]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_tanbeta);
    TEST(test_invtrig_beta);
    TEST(test_beta_gate_refine);

    TEST(test_beta_unit);
    TEST(test_beta_logweight);
    TEST(test_beta_negative);
    TEST(test_trigpower_unit);
    TEST(test_trigpower_vanishing);
    TEST(test_trigpower_negative);
    TEST(test_dispatch_regression);
    TEST(test_method_option);
    TEST(test_stress);

    TEST(test_beta_widened_strip);
    TEST(test_kingreflect_halfline_gate);
    printf("All integrate_beta tests passed.\n");
    return 0;
}
