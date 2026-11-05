/* test_integrate_newton_leibniz.c
 *
 * Unit + stress tests for definite integration by the Newton-Leibniz rule
 * (the fundamental theorem of calculus):
 *
 *   Integrate[f, {x, a, b}]                       (automatic dispatch)
 *   Integrate[f, {x, a, b}, Method -> "NewtonLeibniz"]
 *   Integrate`NewtonLeibniz[f, {x, a, b}]         (explicit entry point)
 *   Integrate`SingularPoints[expr, {x, a, b}]     (the pole detector)
 *
 * Correctness is asserted three ways:
 *   1. exact closed-form results for textbook integrals,
 *   2. the detector's interior-pole list on hand-checked rational inputs,
 *   3. a numeric differential cross-check against NIntegrate for a stress
 *      batch (the symbolic definite value must agree with the numerical one).
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

/* Evaluate `input`, assert its printed form equals `expected`. */
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

/* Assert Integrate[f, {x, a, b}] agrees numerically with NIntegrate. */
static void num_match(const char* f, const char* a, const char* b) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "N[Abs[Integrate[%s, {x, %s, %s}] - "
             "NIntegrate[%s, {x, %s, %s}]]] < 1/100000000",
             f, a, b, f, a, b);
    check_eq(buf, "True");
}

/* ------------------------------------------------------------------ */
/* 1. Pole detector: interior real roots of Denominator[Together[.]].   */
/* ------------------------------------------------------------------ */
static void test_detector(void) {
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 0, 3}]",  "{1}");
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 2, 3}]",  "{}");
    check_eq("Integrate`SingularPoints[1/((x - 1)(x - 2)), {x, 0, 3}]", "{1, 2}");
    check_eq("Integrate`SingularPoints[1/(x^2 + 1), {x, -5, 5}]", "{}");
    /* Pole exactly on an endpoint is NOT interior. */
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 1, 3}]",  "{}");
    /* No x in the denominator -> no poles. */
    check_eq("Integrate`SingularPoints[x^2 + 3 x, {x, -5, 5}]", "{}");
    /* Sorted ascending, several interior. */
    check_eq("Integrate`SingularPoints[1/((x + 1) x (x - 1)), {x, -2, 2}]",
             "{-1, 0, 1}");
    /* Reversed limits (a > b): interior classification is orientation-
     * independent, so the pole between the bounds is still detected. */
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 3, 0}]",  "{1}");
    check_eq("Integrate`SingularPoints[1/((x - 1)(x - 2)), {x, 3, 0}]", "{1, 2}");
    check_eq("Integrate`SingularPoints[1/(x - 1), {x, 3, 2}]",  "{}");
}

/* ------------------------------------------------------------------ */
/* 2. Core FTC on proper integrals (continuous antiderivative).         */
/* ------------------------------------------------------------------ */
static void test_proper(void) {
    check_eq("Integrate[x^2, {x, 0, 1}]",        "1/3");
    check_eq("Integrate[Sin[x], {x, 0, Pi}]",    "2");
    check_eq("Integrate[1/(1 + x^2), {x, 0, 1}]", "1/4 Pi");
    check_eq("Integrate[E^x, {x, 0, 1}]",        "-1 + E");
    check_eq("Integrate[Cos[x], {x, 0, Pi/2}]",  "1");
    /* Pole of the integrand strictly OUTSIDE the interval: proper. */
    check_eq("Integrate[1/x^3, {x, 1, 2}]",      "3/8");
    check_eq("Integrate[1/x^2, {x, -2, -1}]",    "1/2");
    check_eq("Integrate[1/x, {x, 1, 2}]",        "Log[2]");
    check_eq("Integrate[1/x, {x, 1, E}]",        "1");
}

/* ------------------------------------------------------------------ */
/* 3. Improper integrals: one-sided limits, infinite bounds.            */
/* ------------------------------------------------------------------ */
static void test_improper(void) {
    /* Integrable endpoint singularity of the integrand. */
    check_eq("Integrate[1/Sqrt[x], {x, 0, 1}]",         "2");
    /* Infinite upper bound. */
    check_eq("Integrate[1/(1 + x^2), {x, 0, Infinity}]", "1/2 Pi");
    check_eq("Integrate[E^(-x), {x, 0, Infinity}]",      "1");
    check_eq("Integrate[1/x^2, {x, 1, Infinity}]",       "1");
    /* Issue #41: endpoint singularity (x=1) AND infinite upper bound at once.
     * Antiderivative is Sqrt[2] ArcTan[Sqrt[-1+x]/Sqrt[2]]; the endpoint value
     * Limit[Sqrt[2] ArcTan[Sqrt[-1+x]/Sqrt[2]], x -> Infinity] = Pi/Sqrt[2]
     * used to hang the whole integral. */
    check_eq("Integrate[1/((x + 1) Sqrt[x - 1]), {x, 1, Infinity}]",
             "Pi/Sqrt[2]");
}

/* ------------------------------------------------------------------ */
/* 4. Divergent integrals: interior pole -> Integrate::idiv + unevaluated. */
/* ------------------------------------------------------------------ */
static void test_divergent(void) {
    /* Genuine divergence emits Integrate::idiv (to stderr) and leaves the
     * integral unevaluated, matching Mathematica.  The warning does not affect
     * the printed result compared here. */
    check_eq("Integrate[1/x^2, {x, -1, 1}]", "Integrate[1/x^2, {x, -1, 1}]");
    check_eq("Integrate[1/x, {x, -1, 1}]",   "Integrate[1/x, {x, -1, 1}]");
    /* Reversed limits (a > b) with an interior pole are equally divergent:
     * the pole at 0 lies between the bounds regardless of orientation, so this
     * must NOT collapse to a spurious finite branch value (e.g. I Pi). */
    check_eq("Integrate[1/x, {x, 1, -1}]",    "Integrate[1/x, {x, 1, -1}]");
    check_eq("Integrate[1/x, {x, 2, -2}]",    "Integrate[1/x, {x, 2, -2}]");
    check_eq("Integrate[1/(x - 3), {x, 5, 1}]",
             "Integrate[1/(-3 + x), {x, 5, 1}]");
    /* Interior poles at x = +/-1 with F = -ArcTanh[x]: the one-sided boundary
     * values are +/-ArcTanh[1] = +/-Infinity.  These must NOT cancel to a
     * spurious finite branch residue (regression: previously returned -I Pi). */
    check_eq("Integrate[1/(x^2 - 1), {x, -Infinity, Infinity}]",
             "Integrate[1/(-1 + x^2), {x, -Infinity, Infinity}]");
    check_eq("Integrate[1/(x^2 - 1), {x, -Infinity, Infinity}, "
             "Method -> \"NewtonLeibniz\"]",
             "Integrate[1/(-1 + x^2), {x, -Infinity, Infinity}, "
             "Method -> \"NewtonLeibniz\"]");
    check_eq("Integrate[1/(1 - x^2), {x, 0, 2}]",
             "Integrate[1/(1 - x^2), {x, 0, 2}]");
}

/* ------------------------------------------------------------------ */
/* 4c. Reversed limits without an interior pole: Integrate[.,{x,a,b}]    */
/*     with a > b is the negation of the forward integral.               */
/* ------------------------------------------------------------------ */
static void test_reversed_limits(void) {
    check_eq("Integrate[x^2, {x, 1, 0}]",     "-1/3");
    check_eq("Integrate[Sin[x], {x, Pi, 0}]", "-2");
    check_eq("Integrate[1/x, {x, E, 1}]",     "-1");
    /* Reversed with a pole strictly OUTSIDE [min,max]: still proper. */
    check_eq("Integrate[1/x^2, {x, 2, 1}]",   "-1/2");
}

/* ------------------------------------------------------------------ */
/* 4b. Continuous branch-form antiderivatives (Weierstrass): the Floor-  */
/*     corrected antiderivative is continuous, so F(b)-F(a) is exact once */
/*     spurious F-denominator poles are ignored and the value is          */
/*     cross-checked against NIntegrate.                                  */
/* ------------------------------------------------------------------ */
static void test_continuous_periodic(void) {
    check_eq("Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]",   "(2 Pi)/Sqrt[3]");
    check_eq("Integrate[1/(5 + 4 Cos[x]), {x, 0, 2 Pi}]", "2/3 Pi");
    check_eq("Integrate[1/(3 + 2 Sin[x]), {x, 0, 2 Pi}]", "(2 Pi)/Sqrt[5]");
    check_eq("Integrate[Cos[x]^2, {x, 0, 2 Pi}]",         "Pi");
    check_eq("Integrate[Sin[x]^2, {x, 0, 2 Pi}]",         "Pi");
}

/* ------------------------------------------------------------------ */
/* 5. Iterated / multiple integrals (inner = last spec, done first).    */
/* ------------------------------------------------------------------ */
static void test_iterated(void) {
    check_eq("Integrate[x y, {x, 0, 1}, {y, 0, 1}]", "1/4");
    check_eq("Integrate[x, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]", "1/2");
    /* Inner bound depends on the outer variable. */
    check_eq("Integrate[1, {x, 0, 1}, {y, 0, x}]",   "1/2");
}

/* ------------------------------------------------------------------ */
/* 6. Explicit entry point, Method selection, Method pass-through.      */
/* ------------------------------------------------------------------ */
static void test_entry_points(void) {
    check_eq("Integrate`NewtonLeibniz[x^2, {x, 0, 1}]",              "1/3");
    check_eq("Integrate[x^2, {x, 0, 1}, Method -> \"NewtonLeibniz\"]", "1/3");
    /* Method pass-through to the inner indefinite integration. */
    check_eq("Integrate[x^3, {x, 0, 2}, Method -> \"BronsteinRational\"]", "4");
    /* Symbolic bounds with a pole-free antiderivative. */
    check_eq("Integrate[2 x, {x, a, b}]", "-a^2 + b^2");
}

/* ------------------------------------------------------------------ */
/* 7. Unevaluated when the antiderivative is unknown (no wrong value).  */
/* ------------------------------------------------------------------ */
static void test_unevaluated(void) {
    /* No closed-form antiderivative -> stays a definite Integrate. */
    check_eq("Head[Integrate[Sin[x^2 + Log[x]] Cos[x], {x, 0, 1}]]", "Integrate");
    /* Gaussian: the strict Newton-Leibniz method now closes it directly. The
     * antiderivative is (Sqrt[Pi]/2) Erf[x], and constant-factor linearity in
     * Limit resolves Limit[(Sqrt[Pi]/2) Erf[x], x -> Infinity] = Sqrt[Pi]/2,
     * so the boundary difference evaluates to Sqrt[Pi]. (The auto-dispatch
     * also closes it via the even-symmetry -> half-line Ramanujan route.) */
    check_eq("Integrate[E^(-x^2), {x, -Infinity, Infinity}, "
             "Method -> \"NewtonLeibniz\"]", "Sqrt[Pi]");
    check_eq("Integrate[E^(-x^2), {x, -Infinity, Infinity}]", "Sqrt[Pi]");
}

/* ------------------------------------------------------------------ */
/* 8. Stress: symbolic definite value vs NIntegrate.                    */
/* ------------------------------------------------------------------ */
static void test_stress_numeric(void) {
    num_match("x^3 - 2 x + 1", "-1", "2");
    num_match("x Exp[x]",      "0",  "1");
    num_match("Tan[x]",        "0",  "Pi/4");
    num_match("Log[x]",        "1",  "2");
    num_match("x/(1 + x^2)",   "0",  "1");
    num_match("1/Sqrt[1 - x^2]", "0", "1/2");
    num_match("Cosh[x]",       "0",  "1");
    num_match("Sinh[x]",       "0",  "2");
    num_match("x^2 Exp[-x]",   "0",  "1");
    num_match("1/(x^2 + 4)",   "0",  "2");
    num_match("1/(x^2 + 1)^2", "0",  "1");
    num_match("x Log[x]",      "1",  "2");
    /* Pole of the integrand outside the interval: still proper. */
    num_match("1/(1 - x)",     "2",  "5");
    num_match("1/x^3",         "1",  "3");
    /* Continuous branch-form (Weierstrass) antiderivatives. */
    num_match("1/(2 + Cos[x])",   "0", "2 Pi");
    num_match("1/(5 + 4 Cos[x])", "0", "2 Pi");
    num_match("1/(3 + 2 Sin[x])", "0", "2 Pi");
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
 * A never-vanishing factor in the denominator.  Together[E^-t/t^2] has the
 * denominator t^2 E^t; Solve cannot treat it as a unit, so the pole set went
 * undecidable and every improper e^{-t}/t^n integral declined -- and with a
 * parameter the cascade fell through to DiffUnderInt, which ground for a
 * minute probing bases outside the assumed domain (x = -1 under x > 0).
 * References (the reference system): 1/E + ExpIntegralEi[-1];
 * E^-x - x Gamma[0, x] (= E^-x + x ExpIntegralEi[-x] for x > 0).
 * ---------------------------------------------------------------------- */
static void test_nl_exponential_denominator(void) {
    check_quiet("Integrate[Exp[-t]/t^2, {t, 1, Infinity}]", "1/E + ExpIntegralEi[-1]");
    check_quiet("Integrate[Exp[-x t]/t^2, {t, 1, Infinity}, Assumptions -> x > 0]",
                "E^(-x) + x ExpIntegralEi[-x]");
    /* Equal to the reference form numerically, to 30 digits. */
    check_quiet("Abs[N[(E^(-x) + x ExpIntegralEi[-x]) - (E^-x - x Gamma[0, x]) /. x -> 7/3, 30]] < 10^-28",
                "True");
    /* DiffUnderInt pinned: declines, silently and quickly (no idiv from its
     * probes, no base outside x > 0). */
    check_quiet("Integrate`DiffUnderInt[Exp[-x t]/t^2, {t, 1, Infinity}, Assumptions -> x > 0]",
                "Integrate`DiffUnderInt[1/(E^(t x) t^2), {t, 1, Infinity}, Assumptions -> x > 0]");
}

/* The antiderivative of Log[t]/(1 + t^2)^n carries ArcTanh[1 +/- 2 I t]; its
 * Series at the t = 0 branch point dropped the Log of the inner slope
 * (-Log[2 I]/2), so every n >= 2 was off by a multiple of Pi (n = 2 by
 * -Pi/4).  And Limit at Infinity of the rational antiderivative's
 * Log[1 + x + x^2] - Log[1 - x + x^2], scaled, split a Log across a sign
 * change and came back -15/16 I Pi, so Newton-Leibniz called the convergent
 * x^4/(1 + x^2 + x^4)^3 divergent.  Reference values. */
static void test_nl_branch_point_antiderivatives(void) {
    check_eq("Integrate[Log[t]/(1 + t^2)^2, {t, 0, 1}]", "1/8 (-4 Catalan - Pi)");
    check_eq("Integrate[Log[t]/(1 + t^2)^3, {t, 0, 1}]", "1/16 (-1 - 6 Catalan - 2 Pi)");
    /* Mathar M8659, both instantiations. */
    check_eq("Abs[N[Integrate[Log[t]/(1 + t^2)^7, {t, 0, 7/4}] + 0.6655494122101222]] < 10^-10",
             "True");
    check_eq("Abs[N[Integrate[Log[t]/(1 + t^2)^3, {t, 0, 1/5}] + 0.5070095992828864]] < 10^-10",
             "True");
    check_eq("Integrate[x^4/(1 + x^2 + x^4)^3, {x, 0, Infinity}, Method -> \"NewtonLeibniz\"]",
             "Pi/(48 Sqrt[3])");
    check_eq("Limit[3/32 (Log[1 + x + x^2] - Log[1 - x + x^2]), x -> Infinity]", "0");
}

/* Mathar M14298 (symbolic form) recursed without end in the power-generator
 * abstraction: a Together nested inside another read the outer table's
 * ExpSplit`g0 as its own.  It must return (here unevaluated), not crash. */
static void test_nl_nested_expsplit_no_crash(void) {
    check_eq("Head[Quiet[Integrate[s^(1 + mu) (1 - s^2)^nu BesselJ[mu, s t], {s, 0, 1}, "
             "Assumptions -> mu > 0 && nu > 0 && t > 0]]] =!= Null", "True");
}

/* An elliptic antiderivative is cross-checked.  For Mathar M6771 the indefinite
 * cascade returns EllipticE/EllipticF of ArcSin[Sqrt[9/(8 (x - 5/8))]], whose
 * argument exceeds 1 on [5/6, 5/4]: there it is an antiderivative of -f, and it
 * telescoped to -0.12793.  The cross-check rejects it and the elliptic stage
 * answers.  (Reached once DerivativeDivides stopped grinding on the elliptic
 * kernel; before, the integral timed out.)  Reference value. */
static void test_nl_elliptic_antiderivative_crosschecked(void) {
    check_eq("Abs[N[Integrate[Sqrt[((5/4 - x) (-5/8 + x))/(7/4 - x)], {x, 5/6, 5/4}] - "
             "0.1279340024731467]] < 10^-10", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_nl_elliptic_antiderivative_crosschecked);

    TEST(test_nl_branch_point_antiderivatives);
    TEST(test_nl_nested_expsplit_no_crash);

    TEST(test_detector);
    TEST(test_proper);
    TEST(test_improper);
    TEST(test_divergent);
    TEST(test_reversed_limits);
    TEST(test_continuous_periodic);
    TEST(test_iterated);
    TEST(test_entry_points);
    TEST(test_unevaluated);
    TEST(test_stress_numeric);

    TEST(test_nl_exponential_denominator);
    printf("All Integrate NewtonLeibniz tests passed!\n");
    return 0;
}
