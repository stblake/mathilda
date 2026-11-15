/* test_bronstein_params_bronstein.c -- Phase 4.4 step 5 of the Bronstein
 * completion campaign: Bronstein's algorithm over Q(a) (symbolic parameters,
 * transcendental constants with D = 0; theory note E4.4).
 *
 *  - the coefficient fields carry the parameters (brWithParams, brFn,
 *    brFieldOver over a residue field of Q(a)), with an explicit root in
 *    radicals for presentation (brParamRoot);
 *  - the certified Risch DE (Risch`RischDECertified's five-argument form) and
 *    the linear constant systems (brConstantEquations) are over Q(a);
 *  - a kappa in Q(a) builds its special places (brGammaForm,
 *    brMthPowerCofactor);
 *  - special-place torsion specialises the parameters
 *    (brParamSpecialSpecialization, brSpecialTorsionOrderParams);
 *  - end to end, every verdict agrees with the specialisations a = 2, -5/7.
 *
 * The parameter is spelled Global`a: with Bronstein`Private` on the harness's
 * $ContextPath, a bare `a` resolves to the package symbol of that name once
 * one exists (a kernel resolution defect, tasks/lessons.md), whose pattern
 * variables a_ would then capture it.
 */
#include "bronstein_test_utils.h"

/* Step 5: Bronstein over Q(a). */
static void test_bronstein_field(void) {
    /* the coefficient fields carry the parameters */
    BR_EQ("brWithParams[{x, t}]", "{x, t}");
    BR_EQ("Block[{tragerParams = {Global`a, b}}, brWithParams[{x, b}]]", "{x, b, a}");
    BR_TRUE("Block[{tragerParams = {Global`a}}, brFn[{x}] === fieldFn[{x, Global`a}]]");
    BR_TRUE("brFn[{x}] === fieldFn[{x}]");
    BR_TRUE("With[{F = brFn[{x}, t, t^2 - 2]}, F[\"vars\"] === {x} && F[\"mu\"] === t^2 - 2 && N[F[\"theta\"]^2] == 2.]");
    BR_TRUE("Block[{tragerParams = {Global`a}}, brFn[{x}, t, 2 t^2 - Global`a] === fieldFn[{x, Global`a}, {t, t^2 - Global`a/2, Sqrt[Global`a/2]}]]");
    BR_TRUE("Block[{tragerParams = {Global`a}}, brFn[{x}, t, t^3 + Global`a t + 1][\"theta\"] === None]");
    /* an explicit root over Q(a) */
    BR_TRUE("brParamRoot[t - Global`a/2, t] === Global`a/2");
    BR_EQ("brParamRoot[t^3 - Global`a, t]", "a^(1/3)");
    BR_TRUE("Together[(brParamRoot[t^2 + Global`a t + 1, t] /. t -> 0)^2 + Global`a brParamRoot[t^2 + Global`a t + 1, t] + 1] === 0 ||"
            " PossibleZeroQ[brParamRoot[t^2 + Global`a t + 1, t]^2 + Global`a brParamRoot[t^2 + Global`a t + 1, t] + 1]");
    BR_EQ("brParamRoot[t^3 + Global`a t + 1, t]", "None");
    /* a residue field over Q(a) extended by k: its parameters and generator travel */
    BR_TRUE("Block[{tragerParams = {Global`a}}, brFieldOver[{x}, fieldFn[{Global`a}, {t, t^2 - Global`a, None}]] ==="
            " fieldFn[{x, Global`a}, {t, t^2 - Global`a, Sqrt[Global`a]}]]");
    BR_TRUE("Block[{tragerParams = {Global`a}}, brFieldOver[{x}, fieldFn[{Global`a}]] === fieldFn[{x, Global`a}]]");
    BR_TRUE("brFieldOver[{x}, fieldQ[]] === fieldFn[{x}]");
    /* rational in the tower over Q(a) */
    BR_FALSE("brRationalInQ[Global`a + x, {x}]");
    BR_TRUE("Block[{tragerParams = {Global`a}}, brRationalInQ[(Global`a + x)/(Global`a^2 - x), {x}]]");
    BR_FALSE("Block[{tragerParams = {Global`a}}, brRationalInQ[Sqrt[Global`a] + x, {x}]]");
    /* constants of the field: Q(a) elements */
    BR_TRUE("Block[{tragerParams = {Global`a}}, brConstantFieldElementQ[1/(Global`a + 8)]]");
    BR_TRUE("brConstantFieldElementQ[3/4]");
    BR_FALSE("brConstantFieldElementQ[1/(Global`a + 8)]");
    BR_FALSE("Block[{tragerParams = {Global`a}}, brConstantFieldElementQ[Sqrt[Global`a]]]");
    /* the m-th power cofactor: a factor in the parameters alone is constant */
    BR_EQ("Block[{tragerParams = {Global`a}}, brMthPowerCofactor[Global`a x^2, 2]]", "x");
    BR_EQ("brMthPowerCofactor[Global`a x^2, 2]", "None");
    BR_EQ("Block[{tragerParams = {Global`a}}, brMthPowerCofactor[(Global`a + x) x^2, 2]]", "None");
    /* linear constant systems over Q(a): exact zero test, nonzero Q(a) constants are False */
    BR_TRUE("Block[{tragerParams = {Global`a}}, brConstantEquations[{Global`a/(Global`a + 1) - 1 + 1/(Global`a + 1) == 0, Global`a + 1 == 0,"
          " (Global`a - 1) k1 == Global`a^2 - 1}, {k1}]] === {True, False, Expand[Numerator[Together[(Global`a - 1) k1 - (Global`a^2 - 1)]]] == 0}");
    BR_EQ("brConstantEquations[{k1 == 1}, {k1}]", "{k1 == 1}");
    BR_TRUE("Block[{tragerParams = {Global`a}}, brLinearConstantSystemQ[{Global`a k1 + 1/(Global`a + 1) == 0}, {k1}]]");
    BR_FALSE("brLinearConstantSystemQ[{Global`a k1 == 0}, {k1}]");
    BR_FALSE("Block[{tragerParams = {Global`a}}, brLinearConstantSystemQ[{Global`a k1^2 == 1}, {k1}]]");
    /* the constant field Q(i)(a) is not yet supported */
    BR_THROWS("Block[{tragerParams = {Global`a}, tragerConstantI = True}, brParamFieldCheck[]]", "Unimplemented",
              "symbolic parameters together with Q(i) or radical constants");
    BR_EQ("Block[{tragerParams = {Global`a}}, brParamFieldCheck[]]", "Null");
    /* an element of Q(a)[t]/(mu) without a root in radicals is not presented by value */
    BR_THROWS("Block[{tragerParams = {Global`a}}, brTermExpression[brFn[{x}, t, t^3 + Global`a t + 1], t x]]", "Unimplemented",
              "log terms over an algebraic extension of Q(params) outside the per-orbit construction");
    BR_EQ("Block[{tragerParams = {Global`a}}, brTermExpression[brFn[{x}, t, t^2 - Global`a], t x]]", "Sqrt[a] x");
    /* the per-orbit construction first over Q(a) */
    BR_TRUE("Block[{tragerParams = {Global`a}}, brOrbitFirstQ[<|\"P\" -> z^2 - Global`a, \"z\" -> z|>]]");
    BR_FALSE("brOrbitFirstQ[<|\"P\" -> z^2 - 2, \"z\" -> z|>]");
    /* the C Risch DE over Q(a): the five-argument form declares the parameters */
    BR_EQ("Risch`RischDECertified[1/x, Global`a/x, x, {}, {Global`a}]", "a");
    BR_EQ("Risch`RischDECertified[1/x, Global`a/x, x]", "$Failed");
    BR_EQ("Risch`RischDECertified[-2 x, Global`a, x, {}, {Global`a}]", "\"NoSolution\"");
    BR_EQ("Head[Risch`RischDECertified[-2 x, Global`a, x, {}, {x}]]", "Risch`RischDECertified");
    BR_EQ("Block[{tragerParams = {Global`a}}, brRischDE[-2 x, Global`a, x, {}]]", "{\"NoSolution\"}");
    /* the cascade's Bronstein stage stays closed to parameters (step 6 A/B gate) */
    BR_FALSE("Integrate`CascadeGateQ[\"Bronstein\", Sqrt[Global`a + Log[x]]/x, x]");
    BR_TRUE("Integrate`CascadeGateQ[\"Bronstein\", Pi Sqrt[x + Log[x]], x]");
}

/* Step 5: special-place torsion over Q(a) by specialising the parameters,
 * on y^2 = theta^3 + a theta + 1, theta = E^x (1/Sqrt[E^(3x) + a E^x + 1]):
 * residues at theta = 0 only, genus 1. */
static void test_special_torsion_params(void) {
    /* the tower first: brT1 exists (in Bronstein`Private`) only once it is
       built.  The parameter is held in pa: with Bronstein`Private` on the
       harness's $ContextPath, a later bare `a` can resolve to a package
       symbol of that name once one exists (a kernel resolution defect: the
       Global` probe misses the bare-stored Global`a, tasks/lessons.md) */
    br_run("pa = Global`a; Block[{tragerParams = {pa}}, TP = brTower[1/Sqrt[E^(3 x) + pa E^x + 1], x]]; True");
    br_run("Block[{tragerParams = {pa}}, cvP = brCurve[TP, 2, brT1^3 + pa brT1 + 1];"
           " cP = 1/(brT1^3 + pa brT1 + 1); nrP = brNormalResidues[TP, cvP, cP];"
           " ssP = brSpecialSetup[TP, cvP, cP, nrP, brExpSpecialResidues[TP, cvP, cP]]]; True");
    BR_EQ("{ssP[\"genus\"], ssP[\"fr\"][\"a\"], ssP[\"F0\"]}", "{1, 1, <|\"Type\" -> \"Fn\", \"vars\" -> {a}|>}");
    /* a good point: over Q, the curve specialised */
    BR_TRUE("Block[{tragerParams = {pa}}, Module[{r = brParamSpecialSpecialization[TP, cvP, ssP[\"rsF\"], ssP[\"fr\"],"
            " ssP[\"locs\"], ssP[\"rhoLists\"], {pa -> 2}]}, r[\"rs\"][\"F0\"] === fieldQ[] && r[\"cv\"][\"G\"] === 1 + 2 brT1 + brT1^3 &&"
            " r[\"rhoLists\"] === {{1, -1}, {}}]]");
    /* a bad point: theta = 1 (the frame point) is a root of G at a = -2 */
    BR_EQ("Block[{tragerParams = {pa}}, brParamSpecialSpecialization[TP, cvP, ssP[\"rsF\"], ssP[\"fr\"],"
          " ssP[\"locs\"], ssP[\"rhoLists\"], {pa -> -2}]]", "$Failed");
    /* the order: infinite at the specialisation, a certificate over Q(a) */
    BR_THROWS("Block[{tragerParams = {pa}}, brSpecialTorsionOrderParams[TP, cvP, ssP[\"rsF\"], ssP[\"fr\"], ssP[\"locs\"],"
              " ssP[\"rhoLists\"], ssP[\"qb\"][\"coeffs\"][[All, 1]]]]", "NonElementary",
              "residue divisor of infinite order: its orders modulo 7 and 13 disagree");
    /* brSpecialTorsionOrder dispatches to it; torsionOrderWith refuses parameters */
    BR_THROWS("Block[{tragerParams = {pa}}, brSpecialTorsionOrder[TP, cvP, ssP[\"rsF\"], ssP[\"fr\"], ssP[\"locs\"],"
              " ssP[\"rhoLists\"], ssP[\"qb\"][\"coeffs\"][[All, 1]]]]", "NonElementary",
              "residue divisor of infinite order: its orders modulo 7 and 13 disagree");
    BR_THROWS("Block[{tragerParams = {pa}}, torsionOrderWith[<|\"genus\" -> 1|>, Null]]", "Internal", "parameterTorsionWith");
}

/* Step 5 end to end: every verdict agrees with the specialisations at
 * a = 2, -5/7 (tools: E3 examples with a parameter, E^x -> a E^x families). */
static void test_bronstein_end_to_end(void) {
    const char* elem[] = {
        "Sqrt[Global`a + Log[x]]/x", "1/(x Sqrt[Log[x]^2 + Global`a])", "E^x/Sqrt[E^(2 x) + Global`a]", "Sqrt[E^x + Global`a]",
        "1/(x Log[x] Sqrt[Log[x]^3 + Global`a])", "1/Sqrt[E^(3 x) + Global`a]", "1/(x (Global`a + Log[x]) Sqrt[Log[x]^2 + 1])",
        "E^x/(Global`a + E^(3 x))^(1/3)", "1/(x (Global`a + Log[x]^3)^(1/3))",
        "Together[D[x + Log[-E^x + (Global`a + E^(3 x))^(1/3)], x]]" };
    for (size_t i = 0; i < sizeof elem / sizeof elem[0]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf,
                 "With[{f = %s}, With[{F = Quiet[Integrate`Bronstein[f, x]]}, !ListQ[F] && PossibleZeroQ[D[F, x] - f]]]", elem[i]);
        BR_TRUE(buf);
    }
    BR_EQ("Integrate[Sqrt[Global`a + Log[x]]/x, x, Method -> \"Bronstein\"]", "1/3 (2 a Sqrt[a + Log[x]] + 2 Log[x] Sqrt[a + Log[x]])");
    const char* nonelem[] = {
        "Sqrt[Global`a + Log[x]]", "Sqrt[Log[x] + Global`a x]", "Log[x]/Sqrt[x^2 + Global`a]", "1/(x Sqrt[Log[x]^3 + Global`a])",
        /* elementary at a = 1 only: the exceptional set */
        "(x + Global`a)/(x Sqrt[x + Log[x]])",
        "1/(x (Log[x] - 2) Sqrt[Log[x]^3 + Global`a])", "1/Sqrt[1 + Global`a E^x + E^(3 x)]",
        "(1 + Log[x])/(x (Global`a + Log[x] + Log[x]^3)^(1/3))" };
    for (size_t i = 0; i < sizeof nonelem / sizeof nonelem[0]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "MatchQ[Quiet[Integrate`Bronstein[%s, x]], {\"not elementary\", _String}]", nonelem[i]);
        BR_TRUE(buf);
    }
    /* not yet: parameters with Q(i) */
    BR_TRUE("MatchQ[Quiet[Integrate`Bronstein[Sqrt[Global`a + Log[x]] + I Sqrt[Global`a + Log[x]]/x, x]], {\"implementation incomplete\", _String}]");
}

int main(void) {
    br_init();
    test_bronstein_field();
    test_special_torsion_params();
    test_bronstein_end_to_end();
    BR_DONE("bronstein_params_bronstein_tests");
    return 0;
}
