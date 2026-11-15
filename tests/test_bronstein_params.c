/* test_bronstein_params.c -- Phase 4.4 of the Bronstein completion campaign:
 * symbolic parameters (transcendental constants, D = 0; theory note E4.4).
 *
 * Step 1, the soundness guards that must hold before any parameter refusal is
 * lifted -- each fed the bad input it protects against:
 *  - gfMod: a parameter never reaches GF(p) (it would be read as one more
 *    variable);
 *  - residualCanonicalQ: a radical of a non-rational base (Sqrt[a^2 + 1])
 *    prevents the orbit residual's certificate;
 *  - brParamConstantQ / brBCLimited: a constant of Q(a) is a constant, not
 *    "no solution";
 *  - brGammaForm: a kappa in Q(a) does not take the radical-gamma branch;
 *  - brAlgReduceY: the variable is given, not read off G (which may hold a).
 */
#include "bronstein_test_utils.h"

static void test_guards(void) {
    /* gfMod */
    BR_EQ("gfMod[3 x + 1/2, 7]", "4 + 3 x");
    BR_THROWS("Block[{tragerParams = {a}}, gfMod[a x + 1, 7]]", "Internal", "gfModParameter");
    BR_EQ("Block[{tragerParams = {a}}, gfMod[2 x + 9, 7]]", "2 + 2 x");
    /* residualCanonicalQ */
    BR_TRUE("residualCanonicalQ[{x 2^(1/3) + 2^(2/3)/x}]");
    BR_FALSE("residualCanonicalQ[{x Sqrt[a^2 + 1]}]");
    BR_FALSE("residualCanonicalQ[{(1 + a)^(1/3)}]");
    BR_TRUE("residualCanonicalQ[{a/(x + a), 0}]");
    /* a whitelist: any atom that is not a field element blocks it (E4.4e) */
    BR_FALSE("residualCanonicalQ[{x Log[a]}]");
    BR_FALSE("residualCanonicalQ[{x E^a}]");
    BR_FALSE("residualCanonicalQ[{x^a}]");
    BR_FALSE("residualCanonicalQ[{Abs[a] x}]");
    BR_FALSE("residualCanonicalQ[{Pi x}]");
    BR_FALSE("residualCanonicalQ[{x Root[#^3 - 2 &, 1]}]");
    BR_FALSE("residualCanonicalQ[{Sqrt[2] x + Sqrt[3]}]");
    BR_TRUE("residualCanonicalQ[{(1 + 2 I) x/(x^2 + a)}]");
    /* brParamConstantQ */
    BR_TRUE("brParamConstantQ[3/2 + Sqrt[2]]");
    BR_FALSE("brParamConstantQ[a]");
    BR_TRUE("Block[{tragerParams = {a, b}}, brParamConstantQ[(a^2 + 1)/b]]");
    BR_FALSE("Block[{tragerParams = {a}}, brParamConstantQ[a x]]");
    BR_FALSE("Block[{tragerParams = {a}}, brParamConstantQ[a + Log[x]]]");
    /* brGammaForm: a kappa in Q(a) is refused (step 5 builds its places) */
    BR_EQ("brGammaForm[3 x^2, 2]", "{3, x}");
    BR_EQ("brGammaForm[a x^2, 2]", "None");
    BR_THROWS("Block[{tragerParams = {a}}, brGammaForm[a x^2, 2]]", "Unimplemented",
              "special place with a constant kappa depending on a symbolic parameter");
    /* brAlgReduceY: x explicit; with G = x^3 + a the old First[Variables[G]] took a */
    BR_TRUE("Together[brAlgReduceY[1/(1 + Y), Y, 2, x^3 + a, x] - (Y - 1)/(x^3 + a - 1)] === 0");
    BR_TRUE("Together[brAlgReduceY[1/(1 + Y), Y, 2, x^3 + 1, x] - (Y - 1)/x^3] === 0");
}

/* Step 2: the parameter-aware constant field of the Trager core. */
static void test_core_field(void) {
    /* tragerParamsOf: the symbols other than x, not System` constants/heads */
    BR_EQ("tragerParamsOf[(x + b)/Sqrt[x^2 + a] + Pi + Sin[c], x]", "{a, b, c}");
    BR_EQ("tragerParamsOf[1/Sqrt[x^3 + 1], x]", "{}");
    /* constantFieldNumberQ: polynomials in the parameters, rational coefficients */
    BR_TRUE("Block[{tragerParams = {a}}, constantFieldNumberQ[a^2/3 + 1]]");
    BR_FALSE("Block[{tragerParams = {a}}, constantFieldNumberQ[Sqrt[2] a]]");
    BR_FALSE("constantFieldNumberQ[a]");
    /* baseField / extField */
    BR_TRUE("baseField[] === fieldQ[]");
    BR_TRUE("Block[{tragerParams = {a}}, baseField[] === fieldFn[{a}]]");
    BR_TRUE("Block[{tragerParams = {a}}, extField[t, 2 t^2 - a] === fieldFn[{a}, {t, t^2 - a/2, None}]]");
    BR_TRUE("extField[t, t^2 - 2] === fieldNF[t, t^2 - 2, None]");
    /* rc over Q(a)[t]/(t^2 - a): a canonical zero test */
    BR_TRUE("Block[{tragerParams = {a}}, rc[extField[t, t^2 - a], t^3 - a t] === 0]");
    /* qBasisParams: Q-span of residues in Q(a) */
    BR_TRUE("Block[{tragerParams = {a}}, qBasis[baseField[], {a/2, a, 1/(a + 1)}]] ==="
            " <|\"basis\" -> {a/2, 1/(1 + a)}, \"coeffs\" -> {{1, 0}, {2, 0}, {0, 1}}|>");
    BR_TRUE("Block[{tragerParams = {a}}, qBasis[baseField[], {a, 2 a, 3 a}]] ==="
            " <|\"basis\" -> {a}, \"coeffs\" -> {{1}, {2}, {3}}|>");
    /* the guards of the parameter path */
    BR_EQ("Trager`TragerIntegrate[1/(x Sqrt[x^3 + a x + 1]), x]",
          "{\"not elementary\", \"residue divisor is not a torsion point of the Jacobian\"}");
    BR_TRUE("StringQ[Trager`TragerIntegrate[a/Sqrt[x^2 + Sqrt[2]], x][[1]]]");
}

/* dropConstantContent: the content free of x (parameters, their radicals) */
static void test_presentation(void) {
    BR_EQ("dropConstantContent[2 (1 + a^2) x + 4 (1 + a^2) Sqrt[x^2 + 1], x]", "x + 2 Sqrt[1 + x^2]");
    BR_EQ("dropConstantContent[a/Sqrt[1 + a^2] + (a x)/Sqrt[1 + a^2], x]", "1 + x");
    BR_EQ("dropConstantContent[x/(x + 1), x]", "x/(1 + x)");
    BR_TRUE("Block[{tragerParams = {a}}, tidyLog[(a + a x)/(3 a^2), x] === 1 + x]");
    BR_TRUE("tidyLog[(2 + 2 x)/3, x] === 1 + x");
}

/* End to end over Q(a): verified antiderivatives, an exact certificate */
static void test_end_to_end(void) {
    BR_TRUE("StringQ[First[Trager`TragerIntegrate[1/Sqrt[x^3 + a], x]]] && StringContainsQ[First[Trager`TragerIntegrate[1/Sqrt[x^3 + a], x]], \"not elementary\"]");
    const char* elem[] = {
        "1/((x - a) Sqrt[x^2 + 1])", "x/Sqrt[x^4 + a]", "Sqrt[x + a]/x",
        "1/((x + b) Sqrt[x^2 + a])",
        /* a RootSum over Q(a): PossibleZeroQ needed the RootSum precision fix */
        "(x + a)/((x^2 - a) Sqrt[x^2 + a])" };
    for (size_t i = 0; i < sizeof elem / sizeof elem[0]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf,
                 "With[{f = %s}, With[{F = Trager`TragerIntegrate[f, x]}, !ListQ[F] && PossibleZeroQ[D[F, x] - f]]]", elem[i]);
        BR_TRUE(buf);
    }
    /* presentation: no parameter content left in a log argument */
    BR_TRUE("LeafCount[Trager`TragerIntegrate[1/((x - a) Sqrt[x^2 + 1]), x]] < 150");
}

/* Step 3: torsion by specialising the parameters (E4.4d).  Hand-made
 * components on y^2 = x^3 + a:
 *   ctxM: r = Y/(x G) (f = 1/(x y)), residues +-t/a over Q(a)(t), t^2 = a,
 *         at the places (0, +-t); delta = P - (-P) = 2P with P = (0, Sqrt[a])
 *         a flex: order exactly 3 for every a != 0;
 *   ctxQ: r = Y/((x - 1) G), residues +-t/(1 + a), t^2 = 1 + a. */
static void test_torsion_specialisation(void) {
    br_run("tsG = x^3 + a; "
           "ctxM = <|\"B\" -> curveBasis[2, tsG, x], \"genus\" -> 1, \"F\" -> fieldFn[{a}, {t, t^2 - a, None}],"
           " \"residues\" -> {t/a, -t/a}, \"N\" -> 1, \"E\" -> Expand[x tsG], \"Ep\" -> D[Expand[x tsG], x], \"Bp\" -> x|>; "
           "specM = <|\"res\" -> {t/a, -t/a}, \"b\" -> {1, -1}, \"ref\" -> None|>; "
           "ctxQ = <|\"B\" -> curveBasis[2, tsG, x], \"genus\" -> 1, \"F\" -> fieldFn[{a}, {t, t^2 - 1 - a, None}],"
           " \"residues\" -> {t/(1 + a), -t/(1 + a)}, \"N\" -> 1, \"E\" -> Expand[(x - 1) tsG],"
           " \"Ep\" -> D[Expand[(x - 1) tsG], x], \"Bp\" -> x - 1|>; "
           "specQ = <|\"res\" -> {t/(1 + a), -t/(1 + a)}, \"b\" -> {1, -1}, \"ref\" -> None|>; True");
    /* candidates: rationals, deterministic, growing height */
    BR_TRUE("MatchQ[paramSpecializationPoint[{a, b}, 1], {a -> _?(Head[#] === Integer || Head[#] === Rational &), b -> _}]");
    BR_TRUE("paramSpecializationPoint[{a}, 5] === paramSpecializationPoint[{a}, 5]");
    BR_TRUE("Length[DeleteDuplicates[Table[paramSpecializationPoint[{a}, k], {k, 60}]]] > 30");
    /* values: a vanishing denominator is $Failed */
    BR_EQ("paramSpecializeValue[(x + a)/(a - 2), {a -> 3}]", "3 + x");
    BR_EQ("paramSpecializeValue[(x + a)/(a - 2), {a -> 2}]", "$Failed");
    BR_EQ("paramSpecializeValue[{a, 1/a}, {a -> 0}]", "$Failed");
    /* squarefree profile and place degree */
    BR_EQ("paramSquarefreeProfile[(x - a)^2 (x^3 + a), x]", "{{1, 2}, {3, 1}}");
    BR_EQ("paramSquarefreeProfile[x^3 + a /. a -> 0, x]", "{{1, 3}}");
    BR_EQ("paramPlaceDegree[ctxM[\"F\"], x, 1, tsG, ctxM[\"Ep\"], t/a, 2, x]", "1");
    /* the residue field: degree and irreducibility kept */
    BR_TRUE("paramSpecializedField[fieldFn[{a}], {a -> 2}] === fieldQ[]");
    BR_TRUE("paramSpecializedField[fieldFn[{a}, {t, t^2 - a, None}], {a -> 2}] === fieldNF[t, t^2 - 2, None]");
    BR_EQ("paramSpecializedField[fieldFn[{a}, {t, t^2 - a, None}], {a -> 4}]", "$Failed");
    BR_EQ("paramSpecializedField[fieldFn[{a}, {t, t^2 - 1/a, None}], {a -> 0}]", "$Failed");
    /* goodness: the curve profile (a0 = 0), a vanishing denominator (a0 = -1),
       a reducible mu (a0 = 3); a good point (a0 = 1) */
    BR_EQ("paramGoodSpecialization[ctxQ, specQ, {a -> 0}]", "$Failed");
    BR_EQ("paramGoodSpecialization[ctxQ, specQ, {a -> -1}]", "$Failed");
    BR_EQ("paramGoodSpecialization[ctxQ, specQ, {a -> 3}]", "$Failed");
    BR_TRUE("Module[{r = paramGoodSpecialization[ctxQ, specQ, {a -> 1}]}, ListQ[r] && r[[1, \"F\"]] === fieldNF[t, t^2 - 2, None] &&"
            " r[[1, \"residues\"]] === {t/2, -t/2} && r[[2, \"res\"]] === {t/2, -t/2} && r[[1, \"genus\"]] === 1]");
    /* the search returns a good point; the order at it is the generic order */
    BR_TRUE("Block[{tragerParams = {a}}, ListQ[paramSpecialization[ctxM, specM]]]");
    BR_EQ("Block[{tragerParams = {a}}, torsionOrderSpec[ctxM, specM]]", "3");
    BR_EQ("Block[{tragerParams = {a}}, paramTorsionOrderSpec[ctxM, specM]]", "3");
    /* the screen: 3 delta passes, 2 delta is refuted (a proof over Q(a)) */
    BR_TRUE("Block[{tragerParams = {a}}, torsionScreen[ctxM, specM, 3]]");
    BR_FALSE("Block[{tragerParams = {a}}, paramTorsionScreen[ctxM, specM, 2]]");
    /* infinite order at the specialisation: a certificate over Q(a) */
    BR_TRUE("Block[{tragerParams = {a}}, MatchQ[brFailTag[brCatch[torsionOrderSpec[ctxQ, specQ]]], {\"NonElementary\", _String}]]");
}

/* Step 3 end to end: genus 1 over Q(a), both ways */
static void test_torsion_end_to_end(void) {
    /* torsion generically: verified antiderivatives */
    const char* elem[] = {
        "1/(x Sqrt[x^3 + a])",                                    /* 2P, P a flex: order 3 */
        "Together[D[Log[Sqrt[x^3 + a x + 1] - 1], x]]",
        "Together[D[Log[x Sqrt[x^3 + a] - a + x^2], x]]",
        "Together[D[(x + a) Sqrt[x^3 + a x + 1]/(x - 1) + Log[Sqrt[x^3 + a x + 1] - 1], x]]" };
    for (size_t i = 0; i < sizeof elem / sizeof elem[0]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf,
                 "With[{f = %s}, With[{F = Trager`TragerIntegrate[f, x]}, !ListQ[F] && PossibleZeroQ[D[F, x] - f]]]", elem[i]);
        BR_TRUE(buf);
    }
    /* infinite order over Q(a) (both at a = 3, -5/7, 2; elementary at a = 0,
       the countable exceptional set of E4.4b) */
    BR_TRUE("MatchQ[Trager`TragerIntegrate[1/((x - 2) Sqrt[x^3 + a]), x], {\"not elementary\", _String}]");
    BR_TRUE("!ListQ[Trager`TragerIntegrate[1/(x Sqrt[x^3 + 1]), x]]");
}

int main(void) {
    br_init();
    test_guards();
    test_core_field();
    test_presentation();
    test_end_to_end();
    test_torsion_specialisation();
    test_torsion_end_to_end();
    BR_DONE("bronstein_params_tests");
    return 0;
}
