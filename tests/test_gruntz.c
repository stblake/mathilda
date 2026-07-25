/*
 * test_gruntz.c -- Gruntz's mrv limit algorithm (src/calculus/gruntz.c).
 *
 * Exercises the engine both via the explicit  Method -> "Gruntz"  option and
 * via the Automatic cascade fallback, against:
 *   - the benchmark problems from Dominik Gruntz's 1996 ETH PhD thesis
 *     "On Computing Limits in a Symbolic Manipulation System" (Ch. 8 tables
 *     8.1/8.2, plus the worked examples of Ch. 3 and 5);
 *   - a large battery of home-made hard limits that stress the same
 *     machinery (nested exponentials with cancellation, log-towers, radical
 *     differences, one-sided/-Infinity/finite reductions);
 *   - Phase 2 essential-singularity isolation (thesis 5.2) for the
 *     semi-tractable special functions Erf, Erfc, ExpIntegralEi.
 *
 * Notes on the pinned answers:
 *   - Every asserted value has been checked to be mathematically correct for
 *     the exact expression as transcribed (a couple differ from the thesis'
 *     printed answer because the thesis uses a slightly different nesting).
 *   - Genuinely out-of-scope cases (PolyGamma/Zeta Stirling-type expansions we
 *     have no Series-at-Infinity for, the 8.31 Stirling difference, symbolic-sign
 *     limits, bare oscillation under an exclusive Gruntz method) must be left
 *     UNEVALUATED -- an honest gap, never a wrong finite value.
 *     test_honest_abstentions pins exactly that. (Deep log-tower cancellation,
 *     thesis 8.19, now resolves: see test_log_tower.)
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRUNTZ(f, ans) \
    assert_eval_eq("Limit[" f ", x -> Infinity, Method -> \"Gruntz\"]", ans, 0)

/* Basic limits and the finite / one-sided / -Infinity reductions. */
static void test_basic(void) {
    GRUNTZ("(1 + 1/x)^x", "E");
    GRUNTZ("Log[x]", "Infinity");
    GRUNTZ("x/Log[x]", "Infinity");
    GRUNTZ("Log[x]/x", "0");
    assert_eval_eq("Limit[Sin[x]/x, x -> 0, Method -> \"Gruntz\"]", "1", 0);
    assert_eval_eq("Limit[(E^x - 1)/x, x -> 0, Method -> \"Gruntz\"]", "1", 0);
    assert_eval_eq("Limit[(1 - Cos[x])/x^2, x -> 0, Method -> \"Gruntz\"]", "1/2", 0);
    assert_eval_eq("Limit[Log[1 + x]/x, x -> 0, Method -> \"Gruntz\"]", "1", 0);
    assert_eval_eq("Limit[x E^x, x -> -Infinity, Method -> \"Gruntz\"]", "0", 0);
}

/* The cancellation-heavy exp-log benchmark (thesis Table 8.1). */
static void test_thesis_exp_log(void) {
    /* 8.1 */
    GRUNTZ("E^x (E^(1/x - E^-x) - E^(1/x))", "-1");
    /* 8.5 */
    GRUNTZ("E^(E^(E^x + E^-x))/E^(E^(E^x))", "Infinity");
    /* 8.6 (= e for this exact nesting) */
    GRUNTZ("E^(E^(E^x))/E^(E^(E^x - E^(-E^x)))", "E");
    /* 8.7 */
    GRUNTZ("E^(E^(E^x))/E^(E^(E^x - E^(-E^(E^x))))", "1");
    /* 8.8 */
    GRUNTZ("E^(E^x)/E^(E^x - E^(-E^(E^x)))", "1");
    /* 8.9 */
    GRUNTZ("(Log[x])^2 E^(Sqrt[Log[x]] (Log[Log[x]])^2) "
           "E^(Sqrt[Log[Log[x]]] (Log[Log[Log[x]]])^3)/Sqrt[x]", "0");
    /* 8.11 */
    GRUNTZ("(E^(x E^-x/(E^-x + E^(-2 x^2/(x + 1)))) - E^x)/x", "-E^2");
    /* 8.12 */
    GRUNTZ("(3^x + 5^x)^(1/x)", "5");
    /* 8.13 */
    GRUNTZ("x/Log[x^(Log[x]^(2/Log[x]))]", "Infinity");
    /* 8.17 */
    GRUNTZ("E^(x E^-x/(E^-x + E^(-2 x^2/(x + 1))))/E^x", "1");
    /* 8.20 (= e) */
    GRUNTZ("E^(Log[Log[x + E^(Log[x] Log[Log[x]])]]/"
           "Log[Log[Log[E^x + x + Log[x]]]])", "E");
}

/* Tractable trig at a vanishing argument (thesis 8.21, 8.22 / example 5.1). */
static void test_thesis_trig(void) {
    /* 5.1 / 8.21 */
    GRUNTZ("E^x (Sin[1/x + E^-x] - Sin[1/x])", "1");
    /* 8.22 */
    GRUNTZ("E^(E^x) (E^Sin[1/x + E^(-E^x)] - E^Sin[1/x])", "1");
}

/* Worked examples from the algorithm chapters. */
static void test_worked_examples(void) {
    /* Example 3.15 */
    GRUNTZ("(E^(1/x - E^-x) - E^(1/x))/E^-x", "-1");
    /* Example 3.13's mrv driver, end-to-end */
    GRUNTZ("E^(x + E^(-x^2))/E^x", "1");
    /* Example 5.4 (ExpIntegralEi, essential-singularity isolation) */
    GRUNTZ("ExpIntegralEi[x + E^-x] E^-x x", "1");
}

/* PolyGamma at Infinity (Phase-2 isolation, DLMF 5.11.2). Digamma
 * PolyGamma[0,x] ~ Log[x] - 1/(2x) - ...; PolyGamma[m,x] (m>=1) ~
 * (-1)^(m-1)(m-1)!/x^m + .... The isolation branch rewrites PolyGamma[m,g]
 * (g -> +oo) as its Stirling series so the mrv engine can cancel it. */
static void test_polygamma_at_infinity(void) {
    /* Digamma growth rate and its Log[x] head. */
    GRUNTZ("PolyGamma[x]/Log[x]", "1");
    GRUNTZ("PolyGamma[0, x]/Log[x]", "1");
    GRUNTZ("PolyGamma[x] - Log[x]", "0");
    /* Next-order term: psi(x) - Log[x] ~ -1/(2x), so x(...) -> -1/2. */
    GRUNTZ("x (PolyGamma[x] - Log[x])", "-1/2");
    /* Polygamma decay rates: psi^(m)(x) ~ (-1)^(m-1)(m-1)!/x^m. */
    GRUNTZ("x PolyGamma[1, x]", "1");
    GRUNTZ("x^2 PolyGamma[2, x]", "-1");
    GRUNTZ("x^3 PolyGamma[3, x]", "2");
    GRUNTZ("x^4 PolyGamma[4, x]", "-6");
    /* Ratios that exercise the argument-scaling of the series. */
    GRUNTZ("PolyGamma[1, x]/PolyGamma[1, 2 x]", "2");   /* (2x/x)^1 = 2   */
    GRUNTZ("PolyGamma[2, x]/PolyGamma[2, 3 x]", "9");   /* (3x/x)^2 = 9   */
    /* Nested digamma: psi(psi(x)) ~ Log[psi(x)] ~ Log[Log[x]]. */
    GRUNTZ("PolyGamma[0, PolyGamma[0, x]] - Log[Log[x]]", "0");
    /* Exp of digamma: E^psi(x) ~ x E^(-1/(2x)) ~ x. */
    GRUNTZ("Exp[PolyGamma[0, x]]/x", "1");
    /* Plain Automatic must resolve the headline cases too. */
    assert_eval_eq("Limit[PolyGamma[x] - Log[x], x -> Infinity]", "0", 0);
    assert_eval_eq("Limit[x PolyGamma[1, x], x -> Infinity]", "1", 0);
}

/* Max/Min at Infinity, resolved by eventual dominance. Max[a, b] eventually
 * equals whichever argument is larger for large x -- decided by the
 * leading-term sign of a - b (so Max[1/x, 2/x] = 2/x though both -> 0), and the
 * rewrite recurses so nested Max/Min and Max inside a larger expression work. */
static void test_maxmin_at_infinity(void) {
    GRUNTZ("Max[x, x^2]", "Infinity");
    GRUNTZ("Min[x, x^2]", "Infinity");               /* Min = x -> Infinity */
    /* Both arguments -> 0, but 2/x is eventually larger. */
    GRUNTZ("x Max[1/x, 2/x]", "2");
    GRUNTZ("x Min[1/x, 2/x]", "1");
    GRUNTZ("x Max[1/x, 2/x, 3/x]", "3");             /* n-ary */
    /* Converging-from-below: 1 - 1/x > 1 - 2/x eventually. */
    GRUNTZ("Max[1 - 1/x, 1 - 2/x]", "1");
    /* Decaying exponentials: Max keeps the slower-decaying Exp[-x]. */
    GRUNTZ("Exp[x] Max[Exp[-x], Exp[-2 x]]", "1");
    GRUNTZ("Exp[2 x] Min[Exp[-x], Exp[-2 x]]", "1"); /* Min keeps Exp[-2x] */
    GRUNTZ("Min[x, Log[x]]", "Infinity");            /* Log[x] is the smaller */
    GRUNTZ("Max[x + 1, x]", "Infinity");             /* tie-break by constant */
    /* Nested: inner Max = x^2, then Min[x^2, x^3] = x^2. */
    GRUNTZ("Min[Max[x, x^2], x^3]", "Infinity");
    /* Max and Min of the same pair combine cleanly: (2/x - 1/x) x = 1. */
    GRUNTZ("x (Max[1/x, 2/x] - Min[1/x, 2/x])", "1");
    /* Automatic cascade resolves them via the Gruntz layer too. */
    assert_eval_eq("Limit[x Max[1/x, 2/x], x -> Infinity]", "2", 0);
}

/* Modified Bessel functions at Infinity (monotonic, mrv-tractable).
 * BesselK[nu, x] ~ Sqrt[Pi/(2x)] E^-x (decay); BesselI[nu, x] ~ E^x/Sqrt[2 Pi x]
 * (growth, plus an exponentially-dominated I E^-x subdominant term the engine
 * discards). The leading envelope is nu-independent, so a symbolic order works.
 * (BesselJ/BesselY are OSCILLATORY and stay honest abstentions -- see below.) */
static void test_bessel_at_infinity(void) {
    /* BesselK: decay Exp[-x] Sqrt[Pi/(2x)]. */
    GRUNTZ("BesselK[0, x]", "0");
    GRUNTZ("Exp[x] Sqrt[x] BesselK[0, x]", "Sqrt[1/2 Pi]");
    GRUNTZ("x Exp[x] BesselK[0, x]", "Infinity");
    GRUNTZ("BesselK[0, x]/BesselK[0, 2 x]", "Infinity");   /* E^-x / E^-2x */
    GRUNTZ("BesselK[2, x]/BesselK[0, x]", "1");            /* same envelope */
    /* Symbolic order: leading envelope is order-independent. */
    GRUNTZ("Exp[x] Sqrt[x] BesselK[nu, x]", "Sqrt[1/2 Pi]");
    /* BesselI: growth Exp[x]/Sqrt[2 Pi x]. */
    GRUNTZ("BesselI[0, x] Exp[-x] Sqrt[x]", "1/Sqrt[2 Pi]");
    GRUNTZ("BesselI[0, x] Exp[-x]", "0");
    GRUNTZ("BesselI[0, x]/BesselI[0, 2 x]", "0");          /* E^x / E^2x */
    /* Automatic cascade resolves them too. */
    assert_eval_eq("Limit[Exp[x] Sqrt[x] BesselK[0, x], x -> Infinity]", "Sqrt[1/2 Pi]", 0);
    assert_eval_eq("Limit[BesselI[0, x] Exp[-x] Sqrt[x], x -> Infinity]", "1/Sqrt[2 Pi]", 0);
}

/* Zeta at Infinity. Zeta[x] = 1 + 2^-x + 3^-x + ... collapses onto its
 * Dirichlet head as x -> +oo; the isolation replaces Zeta[g] by that truncated
 * exp-log sum so the mrv engine cancels it. The dominant correction is 2^-x. */
static void test_zeta_at_infinity(void) {
    GRUNTZ("Zeta[x]", "1");
    /* Leading correction 2^-x. */
    GRUNTZ("(Zeta[x] - 1) 2^x", "1");
    GRUNTZ("x (Zeta[x] - 1)", "0");                 /* 2^-x beats any power */
    GRUNTZ("Zeta[x]^x", "1");                        /* x 2^-x -> 0 */
    GRUNTZ("Log[Zeta[x] - 1]/x", "-Log[2]");        /* log of the 2^-x scale */
    /* Second correction 3^-x, exposed after subtracting 2^-x. */
    GRUNTZ("(Zeta[x] - 1 - 2^(-x)) 3^x", "1");
    /* Shifted argument in the isolation (fresh-dummy expansion handles it). */
    GRUNTZ("(Zeta[x - 1] - 1) 2^x", "2");
    GRUNTZ("(Zeta[2 x] - 1) 4^x", "1");
    /* Different mrv classes (2^-x vs 4^-x) -> clean domination. */
    GRUNTZ("(Zeta[x] - 1)/(Zeta[2 x] - 1)", "Infinity");
    /* Automatic cascade resolves the headline case too. */
    assert_eval_eq("Limit[(Zeta[x] - 1) 2^x, x -> Infinity]", "1", 0);
}

/* ArcTan at Infinity. ArcTan[x] -> Pi/2 needs the at-infinity Series hook
 * (Pi/2 - 1/x + 1/(3x^3) - ...); the rate limit x(Pi/2 - ArcTan[x]) -> 1 needs
 * that series to compose inside a product (kernel identity
 * ArcTan[u] = Pi/2 - ArcTan[1/u] when the argument blows up). */
static void test_arctan_at_infinity(void) {
    GRUNTZ("ArcTan[x]", "1/2 Pi");
    GRUNTZ("x (Pi/2 - ArcTan[x])", "1");
    GRUNTZ("x^2 (Pi/2 - ArcTan[x]) - x", "0");
    /* Plain Automatic must also resolve the rate limit. */
    assert_eval_eq("Limit[x (Pi/2 - ArcTan[x]), x -> Infinity]", "1", 0);
}

/* Deep log-tower limits with leading-order cancellation (thesis 8.19). Two
 * nested logs of nearly-equal scale cancel, leaving a lower-order tail; the
 * engine surfaces it by factoring the w-pole out of each Log before Series
 * (avoiding the branch-contaminated Log[-a]+Log[-1/a] split) and running the
 * expansion in the positive log-scale P = -Log[w] so Log[-Log[w]] stays real. */
static void test_log_tower(void) {
    /* Building-block cancellations. */
    GRUNTZ("Log[x + Log[x]] - Log[x]", "0");
    GRUNTZ("(Log[x + Log[x]] - Log[x]) x/Log[x]", "1");
    GRUNTZ("x (Log[x + Log[x]] - Log[x])", "Infinity");
    GRUNTZ("Log[Log[x] + Log[Log[x]]] - Log[Log[x]]", "0");
    /* Thesis 8.19, end-to-end. */
    GRUNTZ("(Log[Log[x] + Log[Log[x]]] - Log[Log[x]])/"
           "Log[Log[x] + Log[Log[Log[x]]]] Log[x]", "1");
}

/* Robustness of the log-tower path: the mrv engine must NEVER hang, even on a
 * deep tower it cannot resolve. A 3+-level tower re-introduces a Log[-Log[w]]
 * branch (I Pi) at a depth the single positive-scale substitution cannot fix;
 * feeding that to Series would blow up, so a Complex[] in the frozen-scale
 * input is rejected and the engine abstains FAST instead of looping. These
 * cases previously made plain Limit (Automatic) and Method->"Series"/"Asymptotic"
 * hang because the cascade recurses through the Gruntz layer on deep sub-limits.
 * Each assertion here simply completing (no timeout) is the regression check. */
static void test_log_tower_no_hang(void) {
    /* Plain Limit (Automatic) must resolve 8.19 -- the asymptotic/series cascade
     * layers abstain, then the top-level Gruntz layer closes it. */
    assert_eval_eq("Limit[(Log[Log[x] + Log[Log[x]]] - Log[Log[x]])/"
                   "Log[Log[x] + Log[Log[Log[x]]]] Log[x], x -> Infinity]", "1", 0);
    /* Under the slower cascade methods, 8.19 must terminate (abstain is fine). */
    assert_eval_startswith("Limit[(Log[Log[x] + Log[Log[x]]] - Log[Log[x]])/"
        "Log[Log[x] + Log[Log[Log[x]]]] Log[x], x -> Infinity, Method -> \"Asymptotic\"]",
        "Limit[");
    assert_eval_startswith("Limit[(Log[Log[x] + Log[Log[x]]] - Log[Log[x]])/"
        "Log[Log[x] + Log[Log[Log[x]]]] Log[x], x -> Infinity, Method -> \"Series\"]",
        "Limit[");
    /* (A still-deeper 3-level tower the engine cannot resolve also abstains
     * without hanging -- exercised by the stress battery, kept out of the fast
     * suite because its bounded fallback search costs ~15s.) */
}

/* ---- Home-made hard limits, batch 1: elementary but non-trivial. ------- */
static void test_stress_elementary(void) {
    GRUNTZ("x^(1/x)", "1");
    GRUNTZ("Log[x + 1]/Log[x]", "1");
    GRUNTZ("(Log[x])^(1/x)", "1");
    GRUNTZ("x (Log[x + 1] - Log[x])", "1");
    GRUNTZ("x (E^(1/x) - 1)", "1");
    GRUNTZ("(E^x + x)^(1/x)", "E");
    GRUNTZ("(E^x + E^(2 x))^(1/x)", "E^2");
    GRUNTZ("x^(1/Log[x])", "E");
    GRUNTZ("(1 + 1/x)^(x^2)", "Infinity");
    GRUNTZ("Log[Log[x]]/Log[x]", "0");
    GRUNTZ("Sqrt[x^2 + x] - x", "1/2");
    GRUNTZ("Sqrt[x^2 + x] - Sqrt[x^2 - x]", "1");
    GRUNTZ("Cosh[x]^(1/x)", "E");
    GRUNTZ("Log[E^x + E^(2 x)]/x", "2");
    GRUNTZ("x Sin[1/x]", "1");
    GRUNTZ("x^2 (1 - Cos[1/x])", "1/2");
    GRUNTZ("Sqrt[x] (Sqrt[x + 1] - Sqrt[x])", "1/2");
    GRUNTZ("x - x^2 Log[1 + 1/x]", "1/2");
    GRUNTZ("(2^x + 3^x + 5^x)^(1/x)", "5");
    GRUNTZ("Log[Log[Log[x]]]/Log[Log[x]]", "0");
}

/* ---- Home-made hard limits, batch 2: nested-exp cancellation. ---------- */
static void test_stress_nested_exp(void) {
    GRUNTZ("E^(E^x + E^(-x))/E^(E^x)", "1");
    GRUNTZ("E^(E^x - x)/E^(E^x)", "0");
    GRUNTZ("x (E^(1/x) - E^(1/(x + 1)))", "0");
    GRUNTZ("(E^(E^(E^x)) - E^(E^(E^x - E^-x)))/E^(E^(E^x) - E^x)", "Infinity");
    GRUNTZ("Log[x] Log[Log[x]]/Log[x + Log[x]]", "Infinity");
    GRUNTZ("x ((1 + 1/x)^x - E)", "-1/2 E");
    GRUNTZ("E^(-x) (Cosh[x] + Sinh[x])", "1");
    GRUNTZ("(E^(1/x) + E^(-1/x))^x", "Infinity");
    GRUNTZ("Log[1 + E^x]/x", "1");
    GRUNTZ("E^(x + E^(x + E^x))/E^(E^(E^x))", "Infinity");
    GRUNTZ("(E^(x^2) - E^(x^2 - Log[x]))/E^(x^2 - Log[x])", "Infinity");
    GRUNTZ("E^(E^x)/(E^(E^x) + x)", "1");
    GRUNTZ("(7^x - 5^x)/(7^x + 3^x)", "1");
    GRUNTZ("(7^x - 5^x)/(3^x - 2^x)", "Infinity");
    GRUNTZ("x^100/E^x", "0");
    GRUNTZ("E^x/x^1000", "Infinity");
    GRUNTZ("(Log[x])^x/x^(Log[x])", "Infinity");
    GRUNTZ("E^(Sqrt[Log[x]])/x", "0");
    GRUNTZ("x/E^(Sqrt[Log[x]])", "Infinity");
    GRUNTZ("(x + 1)^(x + 1)/(x^x E^x x)", "0");
    GRUNTZ("(1 + 1/x + 1/x^2)^x", "E");
    GRUNTZ("E^(x + 2 x^2)/E^(3 x^2)", "0");
    GRUNTZ("E^(3 x^2)/E^(x + 2 x^2)", "Infinity");
    GRUNTZ("Log[Log[x^(E^x)]]/x", "1");
    GRUNTZ("(2 + Cos[1/x])^x/2^x", "Infinity");
    GRUNTZ("x^(Log[x])/E^(Log[x]^2)", "1");
}

/* ---- Home-made hard limits, batch 3: finite point / one-sided / -oo. --- */
static void test_stress_reductions(void) {
    GRUNTZ("(x^x)^(1/x)/x", "1");
    GRUNTZ("x^(1/x^2)", "1");
    GRUNTZ("Log[1 + E^x]/x", "1");     /* +Infinity */
    assert_eval_eq("Limit[Log[1 + E^x]/x, x -> -Infinity, Method -> \"Gruntz\"]", "0", 0);
    assert_eval_eq("Limit[Tan[x]/x, x -> 0, Method -> \"Gruntz\"]", "1", 0);
    assert_eval_eq("Limit[(1 - Cos[x])/(x Sin[x]), x -> 0, Method -> \"Gruntz\"]", "1/2", 0);
    assert_eval_eq("Limit[(1/Sin[x] - 1/x), x -> 0, Method -> \"Gruntz\"]", "0", 0);
    assert_eval_eq("Limit[x^x, x -> 0, Direction -> \"FromAbove\", Method -> \"Gruntz\"]", "1", 0);
    assert_eval_eq("Limit[(E^x - 1 - x)/x^2, x -> 0, Method -> \"Gruntz\"]", "1/2", 0);
    assert_eval_eq("Limit[E^(-1/x^2)/x^100, x -> 0, Direction -> \"FromAbove\", "
                   "Method -> \"Gruntz\"]", "0", 0);
    assert_eval_eq("Limit[x (1 + 1/x)^x, x -> Infinity, Method -> \"Gruntz\"]", "Infinity", 0);
}

/* ---- Phase 2: essential-singularity isolation (Erf / Erfc / Ei). ------- */
static void test_thesis_special(void) {
    /* 8.23 (Erf-difference form; the two leading 1's cancel). */
    GRUNTZ("(Erf[x - E^-x] - Erf[x]) E^x E^(x^2)", "-2/Sqrt[Pi]");
    /* Same via the numerically-stable Erfc spelling. */
    GRUNTZ("(Erfc[x] - Erfc[x - E^-x]) E^x E^(x^2)", "-2/Sqrt[Pi]");
    /* Direct singular values. */
    GRUNTZ("Erf[x]", "1");
    GRUNTZ("Erfc[x]", "0");
    GRUNTZ("Erf[Sqrt[x]]", "1");
    /* Erfc tail: x Erfc[x] E^(x^2) -> 1/Sqrt[Pi]. */
    GRUNTZ("x Erfc[x] E^(x^2)", "1/Sqrt[Pi]");
    GRUNTZ("(1 - Erf[x]) E^(x^2) x Sqrt[Pi]", "1");
    /* Erf at -Infinity via the reflection Erf[g] = -Erf[-g]. */
    assert_eval_eq("Limit[Erf[x], x -> -Infinity, Method -> \"Gruntz\"]", "-1", 0);
    /* ExpIntegralEi tail: Ei[x] E^-x x -> 1. */
    GRUNTZ("ExpIntegralEi[x] E^-x x", "1");
    GRUNTZ("ExpIntegralEi[x^2 + x] E^(-x^2 - x) x", "0");
    GRUNTZ("ExpIntegralEi[Log[x]] Log[x]/x", "1");
    /* Ratio of two Erfc's -> E^-2. */
    GRUNTZ("Erfc[x + 1/x]/Erfc[x]", "1/E^2");
}

/* Phase 2, LogGamma: the Stirling series is an additive exp-log head plus a
 * Bernoulli Laurent tail, so isolation stays cheap (unlike Gamma = Exp[LogGamma],
 * whose x^x-scale towers blow up the Series machinery -- see honest_abstentions).
 * Relies on the new Series[LogGamma, {x, Infinity, n}] hook. */
static void test_thesis_loggamma(void) {
    GRUNTZ("LogGamma[x]/(x Log[x])", "1");
    GRUNTZ("LogGamma[2 x]/(x Log[x])", "2");
    /* Stirling constant: LogGamma[x] - (x-1/2)Log[x] + x -> Log[2 Pi]/2. */
    GRUNTZ("LogGamma[x] - (x - 1/2) Log[x] + x", "1/2 Log[2 Pi]");
    /* Next-order Stirling: (LogGamma[x] - x Log[x] + x)/Log[x] -> -1/2. */
    GRUNTZ("(LogGamma[x] - x Log[x] + x)/Log[x]", "-1/2");
    /* Functional-equation residual: x (LogGamma[x+1]-LogGamma[x]-Log[x]) -> 0. */
    GRUNTZ("x (LogGamma[x + 1] - LogGamma[x] - Log[x])", "0");
}

/* Phase 2, Gamma: Gamma[g] = Exp[LogGamma[g]] — the isolation branch builds the
 * x^x-scale Stirling tower from LogGamma's cheap additive series and wraps it in
 * Exp, handing the mrv engine an explicit exp-log head. These resolved only once
 * the FLINT-backed is_zero_poly removed the exponential-cost zero test that the
 * deep tower cancellations previously triggered. (The 8.31 Stirling-difference
 * needs a still-deeper cancellation and remains an honest abstention.) */
static void test_thesis_gamma(void) {
    /* Stirling ratio: Gamma[x]/(Sqrt[2Pi] x^(x-1/2) E^-x) -> 1. */
    GRUNTZ("Gamma[x]/(Sqrt[2 Pi] x^(x - 1/2) E^-x)", "1");
    /* Thesis 5.5: Log[Gamma[Gamma[x]]]/E^x -> Infinity. */
    GRUNTZ("Log[Gamma[Gamma[x]]]/E^x", "Infinity");
}

/* The engine is wired into Limit's Automatic cascade as a fallback: these
 * hard limits must resolve WITHOUT an explicit Method. */
static void test_automatic_fallback(void) {
    assert_eval_eq("Limit[E^x (E^(1/x - E^-x) - E^(1/x)), x -> Infinity]", "-1", 0);
    assert_eval_eq("Limit[(3^x + 5^x)^(1/x), x -> Infinity]", "5", 0);
    assert_eval_eq("Limit[x/Log[x^(Log[x]^(2/Log[x]))], x -> Infinity]", "Infinity", 0);
    /* Phase-2 special functions also resolve via Automatic. */
    assert_eval_eq("Limit[(Erf[x - E^-x] - Erf[x]) E^x E^(x^2), x -> Infinity]", "-2/Sqrt[Pi]", 0);
    assert_eval_eq("Limit[ExpIntegralEi[x + E^-x] E^-x x, x -> Infinity]", "1", 0);
}

/* Documented gaps and honest abstentions: these must be left UNEVALUATED
 * (printed head "Limit["), never collapsed to a wrong finite value.
 *   - the 8.31 Gamma Stirling-difference needs a deeper x^x-tower cancellation
 *     than the Series machinery reaches.
 *   - symbolic-sign limits are genuinely undecidable without assumptions.
 *   - under an *exclusive* Gruntz method, bare oscillation (Sin[x] at oo)
 *     has no mrv expansion (the Automatic cascade handles those elsewhere). */
static void test_honest_abstentions(void) {
    /* 8.31 Stirling-difference: Gamma[x+1]/Sqrt[2Pi] minus its own first two
     * asymptotic terms is a true zero, but proving it needs the x^x-scale
     * Stirling tower expanded and cancelled to an order the dense-polynomial
     * Series machinery does not reach here, so the engine abstains (never a
     * wrong value). The Stirling *ratio* and thesis 5.5 — which do NOT require
     * that deep a cancellation — now resolve; see test_thesis_gamma. */
    assert_eval_startswith(
        "Limit[Gamma[x + 1]/Sqrt[2 Pi] - E^-x (x^(x + 1/2) + x^(x - 1/2)/12), "
        "x -> Infinity, Method -> \"Gruntz\"]",
        "Limit[");
    /* Symbolic-sign: value depends on Sign[s] / which of a,b dominates. */
    assert_eval_startswith(
        "Limit[Log[x]/x^s, x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    assert_eval_startswith(
        "Limit[(a^x + b^x)^(1/x), x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* Bare oscillation under an exclusive Gruntz method. */
    assert_eval_startswith(
        "Limit[Sin[x], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* PolyGamma with a symbolic order has no known asymptotic form. */
    assert_eval_startswith(
        "Limit[PolyGamma[n, x], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* PolyGamma at -Infinity runs into the pole lattice (reflection needed). */
    assert_eval_startswith(
        "Limit[PolyGamma[0, x], x -> -Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* Zeta at -Infinity diverges through the trivial zeros -- no asymptotic. */
    assert_eval_startswith(
        "Limit[Zeta[x], x -> -Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* Same-mrv-class ratio of multi-term exp sums with a base-shifted argument:
     * a *general* mrv-engine gap (the same limit with Zeta stripped out,
     * (2^-x + 3^-x)/(2^-(x+1) + 3^-(x+1)), also abstains), not Zeta-specific. */
    assert_eval_startswith(
        "Limit[(Zeta[x] - 1)/(Zeta[x + 1] - 1), x -> Infinity, Method -> \"Gruntz\"]",
        "Limit[");
    /* Oscillatory Bessel functions (Cos/Sin[x] envelopes) have no mrv
     * expansion: bare J/Y stay unevaluated (their decaying-to-0 case needs a
     * bounded-oscillation squeeze the mrv engine does not do), and amplitude
     * factors that expose the oscillation (Sqrt[x] BesselJ) genuinely have no
     * limit. Never a wrong value. */
    assert_eval_startswith(
        "Limit[BesselJ[0, x], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    assert_eval_startswith(
        "Limit[BesselY[0, x], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    assert_eval_startswith(
        "Limit[Sqrt[x] BesselJ[0, x], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* BesselK at -Infinity enters the complex/oscillatory negative-argument
     * regime -- left unevaluated. */
    assert_eval_startswith(
        "Limit[BesselK[0, x], x -> -Infinity, Method -> \"Gruntz\"]", "Limit[");
    /* Max/Min whose argument ordering hinges on bounded oscillation cannot be
     * decided by the leading-term sign oracle (Sin[x] has no eventual sign),
     * so they abstain -- even though Max[Sin[x], 2] is mathematically 2. */
    assert_eval_startswith(
        "Limit[Max[x, x + Sin[x]], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
    assert_eval_startswith(
        "Limit[Max[Sin[x], 2], x -> Infinity, Method -> \"Gruntz\"]", "Limit[");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_thesis_exp_log);
    TEST(test_thesis_trig);
    TEST(test_worked_examples);
    TEST(test_polygamma_at_infinity);
    TEST(test_zeta_at_infinity);
    TEST(test_bessel_at_infinity);
    TEST(test_maxmin_at_infinity);
    TEST(test_arctan_at_infinity);
    TEST(test_log_tower);
    TEST(test_log_tower_no_hang);
    TEST(test_stress_elementary);
    TEST(test_stress_nested_exp);
    TEST(test_stress_reductions);
    TEST(test_thesis_special);
    TEST(test_thesis_loggamma);
    TEST(test_thesis_gamma);
    TEST(test_automatic_fallback);
    TEST(test_honest_abstentions);

    printf("All gruntz tests passed!\n");
    return 0;
}
