/* test_integrate_cascade_trager_bronstein.c -- the Trager and Bronstein stages
 * of the Automatic Integrate cascade (BRONSTEIN_PLAN.md Phase 8).
 *
 *   Integrate`CascadeGateQ["Trager" | "Bronstein", f, x]   the structural gates
 *   $TragerCascadeBudget / $BronsteinCascadeBudget         the budgets
 *
 * The gates are tested directly on accept/reject integrands; the cascade then
 * closes Trager-only and Bronstein integrands with plain Integrate[f, x]; a
 * certificate announces Integrate::nonelem once and leaves the integral
 * unevaluated (or to the special-function stage); a budget expiry declines and
 * the cascade continues.
 */
#include "bronstein_test_utils.h"
#include "symtab.h"

static const char* const trager_accept[] = {
    "1/(x Sqrt[1 - x^3])", "(x^4 + 1)^(-1/4)", "Sqrt[x^3 + 1]/x", "1/Sqrt[x^2 + 1]",
    "x^(1/3)/(1 + x)", "(x + 1)^(2/3) x", "1/((x - 1) Sqrt[x^3 + 1])", "Sqrt[x] + x^(1/3)",
    "1/(x^2 Sqrt[x^4 - 1])", "(2 x^3 + 3)^(1/3)/(5 x^6 + 7)^(1/3)", "Sqrt[(x + 1)/(x - 1)]",
    "x/(1 + x^2)^(3/2)", "1/(1 + Sqrt[x])", "(x^2 + x)^(1/4)/x", "1/Sqrt[x^3 - x]",
    /* symbolic parameters: admitted at the Trager stage by default (Phase 4.4
       step 6, Integrate`$CascadeParameters = {"Trager"}) */
    "1/((x - a) Sqrt[x^2 + 1])", "1/Sqrt[x^3 + a]", NULL };
static const char* const trager_reject[] = {
    "x^2 + 1", "1/(x^2 + 1)", "Log[x] Sqrt[x]", "Sin[x] Sqrt[x]", "Sqrt[2] x^(1/2)",
    "E^x Sqrt[x]", "Sqrt[x + Log[x]]", "1.5 Sqrt[x]", "x^Pi", "Sqrt[x]^Sqrt[2]",
    "ArcTan[Sqrt[x]]", "Sqrt[2] + x", "Pi Sqrt[x]", NULL };
static const char* const bronstein_accept[] = {
    "Sqrt[x + Log[x]]", "(x + E^x)^(1/3)", "Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x])",
    "1/(x Sqrt[Log[x]^2 + 1])", "E^x/Sqrt[E^(2 x) + 1]", "1/Sqrt[E^x + 1]",
    "Log[x] Sqrt[Log[x] + x]", "Sqrt[Log[x]]", "E^x Sqrt[1 + E^x]", "1/Sqrt[E^(x^2) + x]",
    "Sqrt[Log[Log[x]]]", "x Sqrt[E^x - 1]", "(Log[x] + 1)^(2/3)/x", "Sqrt[E^x + Log[x]]",
    "E^(2 x)/Sqrt[1 - E^(2 x)]",
    /* Shape B: a radical of x with a Log/Exp kernel (v0.6527-28) */
    "Log[x] Sqrt[x^2 + 1]", "E^x/Sqrt[x^2 + 1]", "Log[x + Sqrt[x^3 + 1]]", "E^(x Sqrt[x^3 + 1])/x",
    "Sqrt[x] E^x", NULL };
static const char* const bronstein_reject[] = {
    "Sqrt[x]", "Log[x]", "E^x", "Sqrt[2] Log[x]", "Log[2] Sqrt[x]",
    "Sin[x] Sqrt[Log[x]]", "Sqrt[Log[x]] + Tan[x]", "Sqrt[2] Sqrt[Log[x]]",
    "Erf[x] Sqrt[Log[x]]", "1.5 Sqrt[Log[x]]", "Sqrt[Log[x]]^Sqrt[3]", NULL };

static void test_gates(void) {
    char buf[512];
    for (int i = 0; trager_accept[i]; i++) {
        snprintf(buf, sizeof buf, "Integrate`CascadeGateQ[\"Trager\", %s, x]", trager_accept[i]);
        BR_TRUE(buf);
    }
    for (int i = 0; trager_reject[i]; i++) {
        snprintf(buf, sizeof buf, "Integrate`CascadeGateQ[\"Trager\", %s, x]", trager_reject[i]);
        BR_FALSE(buf);
    }
    for (int i = 0; bronstein_accept[i]; i++) {
        snprintf(buf, sizeof buf, "Integrate`CascadeGateQ[\"Bronstein\", %s, x]", bronstein_accept[i]);
        BR_TRUE(buf);
    }
    for (int i = 0; bronstein_reject[i]; i++) {
        snprintf(buf, sizeof buf, "Integrate`CascadeGateQ[\"Bronstein\", %s, x]", bronstein_reject[i]);
        BR_FALSE(buf);
    }
    BR_EQ("Head[Integrate`CascadeGateQ[\"Other\", Sqrt[x], x]]", "Integrate`CascadeGateQ");
    BR_EQ("Head[Integrate`CascadeGateQ[\"Trager\", Sqrt[x], 2]]", "Integrate`CascadeGateQ");
    BR_TRUE("MemberQ[Attributes[Integrate`CascadeGateQ], Protected]");
    br_assertions++;
    if (!symtab_get_docstring("Integrate`CascadeGateQ") || !symtab_get_docstring("$TragerCascadeBudget") ||
        !symtab_get_docstring("$BronsteinCascadeBudget")) {
        fprintf(stderr, "FAIL: a cascade symbol has no docstring\n");
        exit(1);
    }
    BR_EQ("{$TragerCascadeBudget, $BronsteinCascadeBudget}", "{10, 10}");
    /* symbolic parameters: admitted at the stages Integrate`$CascadeParameters
       names (Phase 4.4 step 6: default {"Trager"}); Pi stays out of Trager's
       gate, a Constant is in Bronstein's */
    BR_EQ("Integrate`$CascadeParameters", "{\"Trager\"}");
    BR_TRUE("Integrate`CascadeGateQ[\"Trager\", 1/Sqrt[x^3 + a], x]");
    BR_FALSE("Integrate`CascadeGateQ[\"Bronstein\", Sqrt[a + Log[x]]/x, x]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = False}, Integrate`CascadeGateQ[\"Trager\", 1/Sqrt[x^3 + a], x]]");
    BR_TRUE("Block[{Integrate`$CascadeParameters = True}, Integrate`CascadeGateQ[\"Trager\", 1/Sqrt[x^3 + a], x]]");
    BR_TRUE("Block[{Integrate`$CascadeParameters = True}, Integrate`CascadeGateQ[\"Trager\", 1/((x - a) Sqrt[x^2 + b]), x]]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = True}, Integrate`CascadeGateQ[\"Trager\", Pi Sqrt[x], x]]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = True}, Integrate`CascadeGateQ[\"Trager\", Sin Sqrt[x], x]]");
    BR_TRUE("Block[{Integrate`$CascadeParameters = True}, Integrate`CascadeGateQ[\"Bronstein\", Sqrt[a + Log[x]]/x, x]]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = True}, Integrate`CascadeGateQ[\"Bronstein\", Sqrt[a + Log[x]] + Sin[x], x]]");
    BR_TRUE("Integrate`CascadeGateQ[\"Bronstein\", Pi Sqrt[x + Log[x]], x]");
    /* a list of stage names admits the parameters at those stages only */
    BR_TRUE("Block[{Integrate`$CascadeParameters = {\"Trager\"}}, Integrate`CascadeGateQ[\"Trager\", 1/Sqrt[x^3 + a], x]]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = {\"Trager\"}}, Integrate`CascadeGateQ[\"Bronstein\", Sqrt[a + Log[x]]/x, x]]");
    BR_TRUE("Block[{Integrate`$CascadeParameters = {\"Bronstein\"}}, Integrate`CascadeGateQ[\"Bronstein\", Sqrt[a + Log[x]]/x, x]]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = {\"Bronstein\"}}, Integrate`CascadeGateQ[\"Trager\", 1/Sqrt[x^3 + a], x]]");
    BR_FALSE("Block[{Integrate`$CascadeParameters = {}}, Integrate`CascadeGateQ[\"Trager\", 1/Sqrt[x^3 + a], x]]");
    br_assertions++;
    if (!symtab_get_docstring("Integrate`$CascadeParameters")) {
        fprintf(stderr, "FAIL: Integrate`$CascadeParameters has no docstring\n");
        exit(1);
    }
}

/* Plain Integrate now closes the Trager-only and Bronstein integrands. */
static void test_cascade_closes(void) {
    BR_ANTIDERIV("1/(x Sqrt[1 - x^3])", "Automatic");
    BR_ANTIDERIV("(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]])", "Automatic");
    BR_ANTIDERIV("1/(x Sqrt[Log[x]^2 + 1])", "Automatic");
    BR_ANTIDERIV("(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3))", "Automatic");
    BR_ANTIDERIV("(((3 x + 1) Log[x] - x^3 + x^2) Sqrt[x + Log[x]] - (2 x^2 - x - 1) Log[x] - 2 x^3 + x^2 + x)/"
                 "(x Log[x]^2 - (x^3 - 2 x^2) Log[x] - x^4 + x^3)", "Automatic");
    BR_EQ("Integrate[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x]", "-2/Sqrt[x + Log[x]]");
}

/* A certificate: Integrate::nonelem once; no elementary-only stage runs after it,
 * the integral stays unevaluated (no special-function form exists here). */
static void test_cascade_certificates(void) {
    BR_EQ("Check[Integrate[Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), x], \"msg\"]", "\"msg\"");
    BR_EQ("Quiet[Head[Integrate[Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), x]]]", "Integrate");
    /* Trager certifies it non-elementary; the special-function stage still
     * runs and answers with its (verified) elliptic form */
    /* (the cascade HOLDS Integrate::nonelem and drops it when a later stage
     * answers, so a special-function answer comes without the message) */
    BR_EQ("Check[Integrate[5/((x - 3) Sqrt[x^3 - 2]), x]; \"ok\", \"msg\"]", "\"ok\"");
    BR_TRUE("Quiet[FreeQ[Integrate[5/((x - 3) Sqrt[x^3 - 2]), x], Integrate]]");
    /* a certified non-elementary integral still reaches the special-function
     * stage: Sqrt[Log[x]] has an Erfi form */
    BR_TRUE("Quiet[FreeQ[Integrate[Sqrt[Log[x]], x], Integrate]]");
}

/* FriCAS comparison C1: the cascade loses what its methods find.
 * C1.1 -- ParallelMixedSpecial returned Inactive[Integrate][f, x] (no progress
 *         at all) as a "complete" answer and hid Trager's certificate (35
 *         algebraic integrands of benchmark 256).
 * C1.2 -- a large tower integrand reaches the Risch stage before Linearity and
 *         DerivativeDivides grind (early probe), and a speculative Linearity
 *         split runs its pieces without the special-function tail; a Risch
 *         certificate skips the elementary-only PMT/Goursat stages. */
static void test_cascade_keeps_method_results(void) {
    BR_EQ("Check[Integrate[x^3/Sqrt[2 + x^6], x], \"msg\", {Integrate::nonelem}]", "\"msg\"");
    BR_EQ("Quiet[Head[Integrate[1/Sqrt[-1 + x^6], x]]]", "Integrate");
    BR_TRUE("Quiet[FreeQ[Integrate[x^3/Sqrt[2 + x^6], x], Inactive]]");
    /* 55 s in Linearity before (elementary; Risch closes it in 10 ms) */
    BR_TRUE("TimeConstrained[With[{f = (-2 (1 + E^x) Log[E^x + x] Log[1 + Tanh[x]^2])/(E^x + x) - "
            "(2 Log[E^x + x]^2 Sech[x]^2 Tanh[x])/(1 + Tanh[x]^2)}, With[{F = Integrate[f, x]}, "
            "FreeQ[F, Integrate] && Abs[N[(D[F, x] - f) /. x -> 13/10, 20]] < 10^-10]], 20, False]");
    /* 55 s in a speculative Linearity split before (non-elementary) */
    BR_EQ("TimeConstrained[Check[Integrate[x Log[Tan[x]] Tan[x] (1 + Tan[x]^2), x]; \"none\", \"msg\", "
          "{Integrate::nonelem}], 20, \"timeout\"]", "\"msg\"");
}

/* With the parameters admitted the Automatic cascade closes a parametric
 * Bronstein integrand; by default (Trager only) the Trager one. */
static void test_cascade_parameters(void) {
    BR_TRUE("With[{F = Integrate[1/((x - a) Sqrt[x^2 + 1]), x]},"
            " FreeQ[F, Integrate] && PossibleZeroQ[D[F, x] - 1/((x - a) Sqrt[x^2 + 1])]]");
    /* not elementary over Q(a): certified by the default cascade */
    /* (1/Sqrt[x^3 + a] is certified too, but the special-function stage then
     * answers with a verified HypergeometricPFQ form and the held message
     * drops; this one has no such form) */
    BR_EQ("Check[Integrate[1/((x - 2) Sqrt[x^3 + a]), x]; \"ok\", \"msg\", {Integrate::nonelem}]", "\"msg\"");
    /* a parametric Kummer piece inside ParallelMixedTower: its own integrator
     * first (3 s), so the compact three-log answer of the parameter-free run
     * comes back, not Trager's 257-leaf RootSum rendering (benchmark 258, cM-2-197) */
    BR_TRUE("With[{F = Integrate[2/(2 + a + E^(2 x))^(1/3), x]}, FreeQ[F, Integrate] && LeafCount[F] <= 130 &&"
            " PossibleZeroQ[D[F, x] - 2/(2 + a + E^(2 x))^(1/3)]]");
    /* ParallelMixed's log combine over Q(params): gathered by argument when
     * smaller; without parameters the identity */
    BR_TRUE("ParallelMixed`Private`pmtParamLogCombine[a Log[1 + x]/(2 + a) + 2 Log[1 + x]/(2 + a) + Log[x], x] === Log[x] + Log[1 + x]");
    BR_TRUE("ParallelMixed`Private`pmtParamLogCombine[Log[1 + x] + Log[1 + x], x] === 2 Log[1 + x]");
    BR_TRUE("Block[{Integrate`$CascadeParameters = True}, With[{F = Integrate[Sqrt[a + Log[x]]/x, x]},"
            " FreeQ[F, Integrate] && PossibleZeroQ[D[F, x] - Sqrt[a + Log[x]]/x]]]");
    BR_TRUE("Block[{Integrate`$CascadeParameters = True}, With[{F = Integrate[1/((x - a) Sqrt[x^2 + 1]), x]},"
            " FreeQ[F, Integrate] && PossibleZeroQ[D[F, x] - 1/((x - a) Sqrt[x^2 + 1])]]]");
}

/* A budget expiry declines: the cascade continues to the later stages. */
static void test_budgets(void) {
    /* a timeout inside a NESTED Integrate (the DerivativeDivides candidate bound)
     * once left the cascade depth raised for the session, silencing every later
     * Integrate::nonelem: the depth is now restored across the unwind */
    br_run("TimeConstrained[Integrate[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]) + Sqrt[Sin[x]]/(Sin[x]^2 + 1), x], 1/2, 0]; True");
    BR_EQ("Check[Integrate[Sqrt[Log[x]]/(1 + Log[x]), x], \"msg\"]", "\"msg\"");
    BR_EQ("Block[{$TragerCascadeBudget = 1/1000}, Head[Integrate[1/(x Sqrt[1 - x^3]), x]]] =!= Hold", "True");
    BR_TRUE("Block[{$BronsteinCascadeBudget = 1/1000},"
            " Module[{F = Integrate[1/(x Sqrt[Log[x]^2 + 1]), x]},"
            "  FreeQ[F, Integrate] || Head[F] === Integrate]]");
    BR_TRUE("Block[{$BronsteinCascadeBudget = Infinity},"
            " PossibleZeroQ[D[Integrate[1/(x Sqrt[Log[x]]), x], x] - 1/(x Sqrt[Log[x]])]]");
}

int main(void) {
    br_init();
    test_gates();
    test_cascade_closes();
    test_cascade_certificates();
    test_cascade_keeps_method_results();
    test_cascade_parameters();
    test_budgets();
    BR_DONE("integrate_cascade_trager_bronstein_tests");
    return 0;
}
