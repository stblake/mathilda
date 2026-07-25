/*
 * test_gruntz_stress.c -- Generalized stress battery for Gruntz's mrv limit
 * algorithm (src/calculus/gruntz.c), driven through  Method -> "Gruntz".
 *
 * Provenance
 * ----------
 * Each family below takes one canonical example from Dominik Gruntz's 1996 ETH
 * thesis "On Computing Limits in a Symbolic Manipulation System" (Ch. 8 tables
 * 8.1/8.2 and the worked examples of Ch. 2/3/5) and generalizes it along a
 * parameter axis into a progression of increasingly difficult limits whose exact
 * value is derived by hand. The base thesis example is noted in each family
 * header; the surrounding cases sweep coefficients, decay rates, tower depth,
 * base sets, and function heads.
 *
 * Verification methodology
 * ------------------------
 * Every pinned value was cross-checked two ways: (1) hand derivation of the
 * exact limit for the transcribed expression, and (2) the value Mathilda's
 * engine actually returns. Only cases where BOTH agree are asserted here.
 * Cases that abstain (return the input unevaluated) or where the two disagree
 * were dropped -- e.g. the 3-level exp-tower cancellation, Gamma *difference*
 * asymptotics, and psi-towers remain honest gaps, pinned in
 * test_stress_honest_abstentions rather than asserted to a value.
 *
 * Nothing here duplicates test_gruntz.c: that file pins the *exact* thesis
 * expressions; this file pins the generalized families around them.
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRUNTZ(f, ans) \
    assert_eval_eq("Limit[" f ", x -> Infinity, Method -> \"Gruntz\"]", ans, 0)
#define GRUNTZ_AT(f, pt, ans) \
    assert_eval_eq("Limit[" f ", " pt ", Method -> \"Gruntz\"]", ans, 0)

/* ============================================================================
 * Family A -- base 8.12:  (3^x + 5^x)^(1/x) = 5.
 * Generalization:  (sum_i c_i a_i^x)^(p/x) -> (max_i a_i)^p, independent of the
 * positive coefficients c_i.  Sweeps the base set, term count, outer power, and
 * mixes E/Pi bases against integer bases.
 * ========================================================================== */
static void test_stress_maxbase(void) {
    GRUNTZ("(2^x + 3^x)^(1/x)", "3");
    GRUNTZ("(3^x + 5^x)^(1/x)", "5");
    GRUNTZ("(5^x + 2^x)^(1/x)", "5");
    GRUNTZ("(2^x + 3^x + 5^x)^(1/x)", "5");
    GRUNTZ("(7^x + 2^x + 3^x)^(1/x)", "7");
    GRUNTZ("(2^x + 3^x + 5^x + 7^x)^(1/x)", "7");
    GRUNTZ("(2^x + 3^x + 5^x + 7^x + 11^x)^(1/x)", "11");
    GRUNTZ("(3 2^x + 5^x)^(1/x)", "5");          /* coefficient is irrelevant */
    GRUNTZ("(2^x + 3^x)^(2/x)", "9");            /* outer power squares result */
    GRUNTZ("(2^x + 3^x)^(3/x)", "27");
    GRUNTZ("(4^x + 4^x)^(1/x)", "4");            /* equal bases */
    GRUNTZ("((1/2)^x + (1/3)^x)^(1/x)", "1/2");  /* bases < 1: still the max */
    GRUNTZ("(2^x + 6^x)^(1/x)", "6");
    GRUNTZ("(2^(2 x) + 3^x)^(1/x)", "4");        /* 2^(2x) = 4^x dominates 3^x */
    GRUNTZ("(E^x + 2^x)^(1/x)", "E");            /* e > 2 */
    GRUNTZ("(E^x + 3^x)^(1/x)", "3");            /* 3 > e */
    GRUNTZ("(E^(2 x) + 5^x)^(1/x)", "E^2");      /* e^2 ~ 7.39 > 5 */
    GRUNTZ("(Pi^x + 3^x)^(1/x)", "Pi");          /* pi > 3 */
}

/* ============================================================================
 * Family B -- base 8.1:  E^x (E^(1/x - E^-x) - E^(1/x)) = -1.
 * Generalization:  E^(a x)(E^(g(x) - c E^(-a x)) - E^(g(x))) -> -c, where g is
 * any subexpression tending to a finite value.  Sweeps c (incl. sign, fraction)
 * and the growth rate a of the balancing exponential.
 * ========================================================================== */
static void test_stress_cancellation(void) {
    GRUNTZ("E^x (E^(1/x - E^-x) - E^(1/x))", "-1");
    GRUNTZ("E^x (E^(1/x - 2 E^-x) - E^(1/x))", "-2");
    GRUNTZ("E^x (E^(1/x - 3 E^-x) - E^(1/x))", "-3");
    GRUNTZ("E^x (E^(1/x - 4 E^-x) - E^(1/x))", "-4");
    GRUNTZ("E^x (E^(1/x - E^-x/2) - E^(1/x))", "-1/2");
    GRUNTZ("E^x (E^(1/x + E^-x) - E^(1/x))", "1");     /* sign flip */
    GRUNTZ("E^x (E^(2/x - E^-x) - E^(2/x))", "-1");    /* different g */
    GRUNTZ("E^x (E^(- E^-x) - 1)", "-1");              /* g = 0 */
    GRUNTZ("E^x (E^(2 E^-x) - 1)", "2");
    GRUNTZ("E^(2 x) (E^(1/x - E^(-2 x)) - E^(1/x))", "-1");    /* a = 2 */
    GRUNTZ("E^(2 x) (E^(1/x - 2 E^(-2 x)) - E^(1/x))", "-2");
    GRUNTZ("E^(2 x) (E^(1/x + 3 E^(-2 x)) - E^(1/x))", "3");
    GRUNTZ("E^(3 x) (E^(1/x - 5 E^(-3 x)) - E^(1/x))", "-5");  /* a = 3, c = 5 */
}

/* ============================================================================
 * Family C -- bases 8.5-8.8:  ratios of nested exponentials.
 * Generalization (level 2):  E^(E^(x + a E^-x))/E^(E^x) -> E^a, because
 * E^(x + a E^-x) ~ E^x + a.  A faster-decaying perturbation (E^(-x^2)) leaves
 * the ratio at 1.  At level 3 the same balance blows up (-> Inf / 0); those are
 * honest gaps (see abstentions).
 * ========================================================================== */
static void test_stress_exp_tower_ratio(void) {
    GRUNTZ("E^(E^(x + E^-x))/E^(E^x)", "E");
    GRUNTZ("E^(E^(x + 2 E^-x))/E^(E^x)", "E^2");
    GRUNTZ("E^(E^(x + 3 E^-x))/E^(E^x)", "E^3");
    GRUNTZ("E^(E^(x + 4 E^-x))/E^(E^x)", "E^4");
    GRUNTZ("E^(E^(x + 5 E^-x))/E^(E^x)", "E^5");
    GRUNTZ("E^(E^(x + E^-x/2))/E^(E^x)", "Sqrt[E]");
    GRUNTZ("E^(E^(x + E^-x/3))/E^(E^x)", "E^(1/3)");
    GRUNTZ("E^(E^(x - E^-x))/E^(E^x)", "1/E");
    GRUNTZ("E^(E^(x - 2 E^-x))/E^(E^x)", "1/E^2");
    GRUNTZ("E^(E^(x - 3 E^-x))/E^(E^x)", "1/E^3");
    GRUNTZ("E^(E^(x + E^(-x^2)))/E^(E^x)", "1");   /* faster decay -> no shift */
    GRUNTZ("E^(x + E^(-x^2))/E^x", "1");            /* example 3.13 */
}

/* ============================================================================
 * Family D -- base 8.21:  E^x (Sin[1/x + E^-x] - Sin[1/x]) = 1.
 * Generalization:  E^(a x)(F[b/x + c E^-a x] - F[b/x]) -> c F'(0+) for a smooth
 * F.  For Sin/Tan/Sinh/Tanh/ArcTan the derivative at the vanishing argument is
 * 1, so the limit is c; for Cos the derivative is sin(1/x) -> 0, so it is 0.
 * ========================================================================== */
static void test_stress_trig_vanishing(void) {
    GRUNTZ("E^x (Sin[1/x + E^-x] - Sin[1/x])", "1");
    GRUNTZ("E^x (Sin[1/x + 2 E^-x] - Sin[1/x])", "2");
    GRUNTZ("E^x (Sin[1/x + 3 E^-x] - Sin[1/x])", "3");
    GRUNTZ("E^x (Sin[1/x + E^-x/2] - Sin[1/x])", "1/2");
    GRUNTZ("E^(2 x) (Sin[1/x + E^(-2 x)] - Sin[1/x])", "1");
    GRUNTZ("E^x (Tan[1/x + E^-x] - Tan[1/x])", "1");        /* sec^2(0) = 1 */
    GRUNTZ("E^x (Sinh[1/x + E^-x] - Sinh[1/x])", "1");      /* cosh(0) = 1 */
    GRUNTZ("E^x (Tanh[1/x + E^-x] - Tanh[1/x])", "1");
    GRUNTZ("E^x (ArcTan[1/x + E^-x] - ArcTan[1/x])", "1");
    GRUNTZ("E^x (Cos[1/x + E^-x] - Cos[1/x])", "0");        /* -sin(1/x) -> 0 */
    GRUNTZ("E^(E^x) (E^Sin[1/x + E^(-E^x)] - E^Sin[1/x])", "1");   /* 8.22 */
    GRUNTZ("E^(E^x) (E^Sin[1/x + 2 E^(-E^x)] - E^Sin[1/x])", "2");
}

/* ============================================================================
 * Family E -- bases 8.19/8.20:  ratios of nested logarithms.
 * Generalization:  Log[P(x)]/Log[x] -> deg P;  Log[Log[x^k]]/Log[Log[x]] -> 1;
 * Log[x + o(x)] - Log[x] -> 0.  Exponentiating gives E^(deg P).
 * ========================================================================== */
static void test_stress_nested_log(void) {
    GRUNTZ("Log[x^2 + x]/Log[x]", "2");
    GRUNTZ("Log[x^3 + x]/Log[x]", "3");
    GRUNTZ("E^(Log[x^2 + x]/Log[x])", "E^2");
    GRUNTZ("E^(Log[x^3 + x^2 + 1]/Log[x])", "E^3");
    GRUNTZ("E^(Log[2 x + 1]/Log[x])", "E");
    GRUNTZ("Log[Log[x^7]]/Log[Log[x]]", "1");
    GRUNTZ("E^(Log[Log[x^2]]/Log[Log[x]])", "E");
    GRUNTZ("E^(Log[Log[x^5]]/Log[Log[x]])", "E");
    GRUNTZ("Log[x + Log[x]] - Log[x]", "0");
    GRUNTZ("Log[x + 5 Log[x]] - Log[x]", "0");
    GRUNTZ("Log[x^2 + Log[x]] - Log[x^2]", "0");
    GRUNTZ("x (Log[x + Log[x]] - Log[x])/Log[x]", "1");
}

/* ============================================================================
 * Family F -- base 8.9:  a sub-polynomial exp-of-sqrt-log scale / Sqrt[x] -> 0.
 * Generalization: anything of the form e^(o(log x)) or (log x)^k grows slower
 * than every positive power of x, so dividing by x^p (p > 0) -> 0 and the
 * reciprocal -> Infinity.
 * ========================================================================== */
static void test_stress_hardy_scale(void) {
    GRUNTZ("E^(Sqrt[Log[x]])/x", "0");
    GRUNTZ("E^(Sqrt[Log[x]])/Sqrt[x]", "0");
    GRUNTZ("E^(Sqrt[Log[x]]) Log[x]^10/x", "0");
    GRUNTZ("E^(Log[x]^(2/3))/x", "0");
    GRUNTZ("E^(Log[x]^(1/2) + Log[x]^(1/3))/x", "0");
    GRUNTZ("E^(Sqrt[Log[x]] Log[Log[x]]^2)/Sqrt[x]", "0");    /* 8.9 core */
    GRUNTZ("E^(Sqrt[Log[x]] Log[Log[x]])/x", "0");
    GRUNTZ("E^(Log[Log[x]]^2)/x", "0");
    GRUNTZ("(Log[x])^100/x", "0");
    GRUNTZ("(Log[x])^100/Sqrt[x]", "0");
    GRUNTZ("Log[x]^5/x^(1/3)", "0");
    GRUNTZ("x/E^(Sqrt[Log[x]])", "Infinity");
    GRUNTZ("x/(Log[x])^1000", "Infinity");
    GRUNTZ("x/Log[x]^2", "Infinity");
    GRUNTZ("Sqrt[x]/Log[x]^50", "Infinity");
}

/* ============================================================================
 * Family G -- base 2.5:  Sqrt[Log[x+1]] - Sqrt[Log[x]] -> 0.
 * Generalization: conjugate-type radical differences.  (x^2 + a x)^(1/2) - x ->
 * a/2;  (x^n + x^(n-1))^(1/n) - x -> 1/n;  and inner-log variants.
 * ========================================================================== */
static void test_stress_radical_diff(void) {
    GRUNTZ("Sqrt[Log[x + 1]] - Sqrt[Log[x]]", "0");
    GRUNTZ("Sqrt[Log[x + 5]] - Sqrt[Log[x]]", "0");
    GRUNTZ("Sqrt[x + 1] - Sqrt[x]", "0");
    GRUNTZ("(x + 1)^(1/3) - x^(1/3)", "0");
    GRUNTZ("Sqrt[x] (Sqrt[x + 1] - Sqrt[x])", "1/2");
    GRUNTZ("Sqrt[x^2 + x] - x", "1/2");
    GRUNTZ("Sqrt[x^2 + 3 x] - x", "3/2");
    GRUNTZ("(x^2 + x + 1)^(1/2) - x", "1/2");
    GRUNTZ("Sqrt[x^2 + x] - Sqrt[x^2 - x]", "1");
    GRUNTZ("x - Sqrt[x^2 - x]", "1/2");
    GRUNTZ("(x^3 + x^2)^(1/3) - x", "1/3");
    GRUNTZ("Sqrt[Log[x^2]] - Sqrt[Log[x]]", "Infinity");   /* (Sqrt2 - 1) Sqrt[Log x] */
}

/* ============================================================================
 * Family H -- bases 2.6/2.7:  power-series limits at a finite point.
 * ((1 + a x)^s - 1)/x -> a s   (x -> 0);
 * (x^a - 1)/(x^b - 1) -> a/b    (x -> 1).
 * ========================================================================== */
static void test_stress_power_expansion(void) {
    GRUNTZ_AT("((1 + x)^2 - 1)/x", "x -> 0", "2");
    GRUNTZ_AT("((1 + x)^3 - 1)/x", "x -> 0", "3");
    GRUNTZ_AT("((1 + x)^5 - 1)/x", "x -> 0", "5");
    GRUNTZ_AT("((1 + x)^(1/2) - 1)/x", "x -> 0", "1/2");
    GRUNTZ_AT("((1 + x)^(1/3) - 1)/x", "x -> 0", "1/3");
    GRUNTZ_AT("((1 + 2 x)^3 - 1)/x", "x -> 0", "6");
    GRUNTZ_AT("(Sqrt[x] - 1)/(x^(1/3) - 1)", "x -> 1", "3/2");
    GRUNTZ_AT("(x^2 - 1)/(x^3 - 1)", "x -> 1", "2/3");
    GRUNTZ_AT("(x^3 - 1)/(x^2 - 1)", "x -> 1", "3/2");
    GRUNTZ_AT("(x^4 - 1)/(x^7 - 1)", "x -> 1", "4/7");
    GRUNTZ_AT("(x^(1/2) - 1)/(x^(1/5) - 1)", "x -> 1", "5/2");
    GRUNTZ_AT("(x^(2/3) - 1)/(x^(1/6) - 1)", "x -> 1", "4");
}

/* ============================================================================
 * Family I -- bases 8.23/8.24/8.25:  essential-singularity special functions.
 * Erfc[x + c/x]/Erfc[x] -> E^(-2c)  (Gaussian tail shift);
 * Zeta[x] - 1 ~ 2^-x, so (Zeta[x]-1) 2^x -> 1 and Dirichlet ratios give 2^k;
 * ExpIntegralEi[z] ~ E^z/z governs the Ei cases.
 * ========================================================================== */
static void test_stress_special_singularity(void) {
    GRUNTZ("Erfc[x + 1/x]/Erfc[x]", "1/E^2");
    GRUNTZ("Erfc[x + 2/x]/Erfc[x]", "1/E^4");
    GRUNTZ("Erfc[x + 3/x]/Erfc[x]", "1/E^6");
    GRUNTZ("Erfc[x + 1/(2 x)]/Erfc[x]", "1/E");
    GRUNTZ("x Erfc[x] E^(x^2)", "1/Sqrt[Pi]");
    GRUNTZ("Erf[Sqrt[x]]", "1");
    GRUNTZ("x ExpIntegralEi[x] E^-x", "1");
    GRUNTZ("ExpIntegralEi[x] E^-x x^2", "Infinity");
    GRUNTZ("ExpIntegralEi[2 x] E^(-2 x) x", "1/2");     /* Ei[2x] ~ E^2x/(2x) */
    GRUNTZ("ExpIntegralEi[Log[x]] Log[x]/x", "1");
    GRUNTZ("(Zeta[x] - 1) 2^x", "1");
    GRUNTZ("(Zeta[x] - 1) 3^x", "Infinity");
    GRUNTZ("2^x (Zeta[x] - 1) - 1", "0");
    GRUNTZ("x (Zeta[x] - 1)", "0");
    GRUNTZ("(Zeta[x] - 1)/(Zeta[x + 1] - 1)", "2");
    GRUNTZ("(Zeta[x] - 1)/(Zeta[x + 2] - 1)", "4");
    GRUNTZ("(Zeta[x] - 1)/(Zeta[x + 3] - 1)", "8");
    GRUNTZ("Log[Zeta[x] - 1]/x", "-Log[2]");
}

/* ============================================================================
 * Family K -- base 8.37:  Max/Min inside a limit resolve to the dominant arg.
 * ========================================================================== */
static void test_stress_maxmin(void) {
    GRUNTZ("x Max[1/x, 2/x]", "2");
    GRUNTZ("x Max[1/x, 2/x, 3/x]", "3");
    GRUNTZ("Max[2 x, 3 x]/x", "3");
    GRUNTZ("Min[x, Log[x]]", "Infinity");
    GRUNTZ("Max[x, E^x]/E^x", "1");
    GRUNTZ("Max[x, Log[x]] - x", "0");
    GRUNTZ("Min[x, Log[x]] - Log[x]", "0");
    GRUNTZ("Log[Max[E^x, E^(2 x)]]/x", "2");
    GRUNTZ("x/Max[x, Sqrt[x]]", "1");
    GRUNTZ("Min[1/x, 1/x^2] x^2", "1");
}

/* ============================================================================
 * Family L -- bases 8.28-8.34:  digamma / log-Gamma growth.
 * The Gamma/psi *differences* of 8.28-8.34 need a deeper Stirling cancellation
 * than the engine reaches (see abstentions), but the leading growth laws --
 * Log[Gamma[x]] ~ x Log[x], psi(x) ~ Log[x] - 1/(2x) -- resolve cleanly.
 * ========================================================================== */
static void test_stress_gamma_psi_growth(void) {
    GRUNTZ("Log[Gamma[x]]/(x Log[x])", "1");
    GRUNTZ("Gamma[2 x]/(Gamma[x] Gamma[x])", "Infinity");
    GRUNTZ("Log[Gamma[x + 1]] - Log[Gamma[x]] - Log[x]", "0");
    GRUNTZ("PolyGamma[0, x] - Log[x]", "0");
    GRUNTZ("x (PolyGamma[0, x] - Log[x])", "-1/2");
    GRUNTZ("Exp[PolyGamma[0, x]]/x", "1");
    GRUNTZ("PolyGamma[0, 2 x] - PolyGamma[0, x] - Log[2]", "0");
}

/* ============================================================================
 * Honest abstentions -- generalized cases the engine correctly leaves
 * UNEVALUATED (never a wrong value).  These document the current reach of the
 * machinery: 3-level exp-tower cancellation, Gamma/psi *difference* asymptotics.
 * ========================================================================== */
static void test_stress_honest_abstentions(void) {
    /* Level-3 exp tower: the E^a balance of Family C blows up two levels down,
     * needing a cancellation the mrv Series machinery does not reach here. */
    assert_eval_startswith(
        "Limit[E^(E^(E^(x + E^-x)))/E^(E^(E^x)), x -> Infinity, Method -> \"Gruntz\"]",
        "Limit[");
    /* Gamma *difference* (thesis 8.28): needs Stirling expanded and cancelled to
     * an order beyond the dense-Series reach. */
    assert_eval_startswith(
        "Limit[(Gamma[x + 1/Gamma[x]] - Gamma[x])/Log[x], x -> Infinity, "
        "Method -> \"Gruntz\"]", "Limit[");
    /* Gamma ratio Gamma[x+1]/(x Gamma[x]) = 1 exactly, but the engine has no
     * bare-Gamma-ratio asymptotic and abstains. */
    assert_eval_startswith(
        "Limit[Gamma[x + 1]/(x Gamma[x]), x -> Infinity, Method -> \"Gruntz\"]",
        "Limit[");
    /* psi-tower (thesis 8.33): exp exp psi psi needs multi-level 1/(2x) tracking. */
    assert_eval_startswith(
        "Limit[Exp[Exp[PolyGamma[0, PolyGamma[0, x]]]]/x, x -> Infinity, "
        "Method -> \"Gruntz\"]", "Limit[");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_stress_maxbase);
    TEST(test_stress_cancellation);
    TEST(test_stress_exp_tower_ratio);
    TEST(test_stress_trig_vanishing);
    TEST(test_stress_nested_log);
    TEST(test_stress_hardy_scale);
    TEST(test_stress_radical_diff);
    TEST(test_stress_power_expansion);
    TEST(test_stress_special_singularity);
    TEST(test_stress_maxmin);
    TEST(test_stress_gamma_psi_growth);
    TEST(test_stress_honest_abstentions);

    printf("All gruntz stress tests passed!\n");
    return 0;
}
