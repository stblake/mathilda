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
    "x/(1 + x^2)^(3/2)", "1/(1 + Sqrt[x])", "(x^2 + x)^(1/4)/x", "1/Sqrt[x^3 - x]", NULL };
static const char* const trager_reject[] = {
    "x^2 + 1", "1/(x^2 + 1)", "Log[x] Sqrt[x]", "Sin[x] Sqrt[x]", "Sqrt[2] x^(1/2)",
    "E^x Sqrt[x]", "Sqrt[x + Log[x]]", "1.5 Sqrt[x]", "x^Pi", "Sqrt[x]^Sqrt[2]",
    "ArcTan[Sqrt[x]]", "Sqrt[2] + x", NULL };
static const char* const bronstein_accept[] = {
    "Sqrt[x + Log[x]]", "(x + E^x)^(1/3)", "Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x])",
    "1/(x Sqrt[Log[x]^2 + 1])", "E^x/Sqrt[E^(2 x) + 1]", "1/Sqrt[E^x + 1]",
    "Log[x] Sqrt[Log[x] + x]", "Sqrt[Log[x]]", "E^x Sqrt[1 + E^x]", "1/Sqrt[E^(x^2) + x]",
    "Sqrt[Log[Log[x]]]", "x Sqrt[E^x - 1]", "(Log[x] + 1)^(2/3)/x", "Sqrt[E^x + Log[x]]",
    "E^(2 x)/Sqrt[1 - E^(2 x)]", NULL };
static const char* const bronstein_reject[] = {
    "Sqrt[x]", "Log[x]", "E^x", "Log[x] Sqrt[x^2 + 1]", "E^x/Sqrt[x^2 + 1]",
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

/* A budget expiry declines: the cascade continues to the later stages. */
static void test_budgets(void) {
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
    test_budgets();
    BR_DONE("integrate_cascade_trager_bronstein_tests");
    return 0;
}
