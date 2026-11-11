/* test_bronstein_residues.c -- expansions at the special places (Puiseux.m),
 * reduction at the special places of an exponential theta (Special.m) and the
 * residues (Residues.m), BRONSTEIN_PLAN.md Phase 4.
 *
 * Pinned: [B] p.139 fails the resultant criterion (a certificate); the
 * component of [B] p.147 (the paper's theta-form integrand) has residues +-1
 * (its share log((x+y)/(x-y)) of the answer) and none at infinity; [B] p.158 is
 * reduced to 0 at the special places (its algebraic part 3x); arcsinh-type
 * integrands carry their residues at the places above infinity.
 * Properties: series identities, and an exact derivative D(b Y) has no residue
 * anywhere and reduces to 0 at the special places.
 */
#include "bronstein_test_utils.h"

static void setup(void) {
    br_run("TL1 = brTower[Sqrt[Log[x]^2 + 1], x];"
           " TE1 = brTower[Sqrt[E^x + 1], x];"
           " TE2 = brTower[Sqrt[E^(2 x) + 1], x];"
           " T139 = brTower[Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), x];"
           " T158 = brTower[(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3)), x];"
           " T147 = brTower[Sqrt[x + Log[x]], x]; True");
    br_run("cvL1 = brCurve[TL1, 2, brT1^2 + 1]; cvE1 = brCurve[TE1, 2, brT1 + 1];"
           " cvE2 = brCurve[TE2, 2, brT1^2 + 1]; cv139 = brCurve[T139, 3, brT2];"
           " cv158 = brCurve[T158, 3, (brT1 + x)^2]; cv147 = brCurve[T147, 2, brT1 + x]; True");
}

static void test_series(void) {
    BR_EQ("brThetaOrder[brT1^3 + x brT1^2, brT1]", "2");
    BR_EQ("brThetaOrder[x + brT1, brT1]", "0");
    BR_EQ("brThetaOrder[5, brT1]", "0");
    BR_THROWS("brThetaOrder[0, brT1]", "Internal", "brThetaOrderZero");
    BR_EQ("brSeriesCoefficients[1/(1 - u), u, 3]", "{1, 1, 1, 1}");
    BR_EQ("brSeriesCoefficients[(1 + x u)/(1 + u), u, 2]", "{1, -1 + x, 1 - x}");
    BR_THROWS("brSeriesCoefficients[1/u, u, 2]", "Internal", "brSeriesCoefficientsPole");
    BR_EQ("brBinomialSeries[u, 1/2, u, 3]", "{1, 1/2, -1/8, 1/16}");
    BR_EQ("brBinomialSeries[{x}, -1, u, 2]", "{1, -x, x^2}");
    BR_EQ("brBinomialSeries[0, 1/3, u, 2]", "{1, 0, 0}");
    BR_THROWS("brBinomialSeries[1 + u, 1/2, u, 2]", "Internal", "brBinomialSeriesConstantTerm");
    BR_EQ("brLaurentInf[TL1, (brT1^2 + x)/(brT1 - 1), -1]", "{{1, 1}, {0, 1}, {-1, 1 + x}}");
    BR_EQ("brLaurentInf[TL1, 0, 0]", "{}");
    BR_EQ("brLaurentInf[TL1, 1/brT1, 0]", "{}");
    BR_EQ("brLaurentZero[TE1, (1 + brT1)/(brT1^2 x), 0]", "{{-2, 1/x}, {-1, 1/x}, {0, 0}}");
    BR_EQ("brLaurentZero[TE1, brT1^3, 1]", "{}");
    /* properties: (Sum rho_i u^i)^m == Ghat mod u^(n+1); Laurent truncation */
    BR_TRUE("SeedRandom[71]; And @@ Table[Module[{a = Sum[RandomInteger[{-3, 3}] u^i, {i, 1, 3}], m = RandomInteger[{2, 4}], n = 5, s},"
            "   s = brBinomialSeries[a, 1/m, u, n];"
            "   Expand[Normal[Series[Sum[s[[i + 1]] u^i, {i, 0, n}]^m - (1 + a), {u, 0, n}]]] === 0], {20}]");
    BR_TRUE("SeedRandom[72]; And @@ Table[Module[{e = (brT1^RandomInteger[{1, 3}] + x)/(brT1 + RandomInteger[{1, 3}]), l, d},"
            "   l = brLaurentInf[TL1, e, -3]; d = Together[e - Total[(#[[2]] brT1^#[[1]]) & /@ l]];"
            "   Exponent[Numerator[d], brT1] - Exponent[Denominator[d], brT1] < -3], {20}]");
}

static void test_place_data(void) {
    BR_EQ("brPlaceDataInf[TL1, cvL1]", "<|\"g\" -> 2, \"lc\" -> 1, \"count\" -> 2, \"r\" -> 1, \"unramified\" -> True|>");
    BR_EQ("brPlaceDataInf[TE1, cvE1][\"r\"]", "2");
    BR_EQ("brPlaceDataInf[T158, cv158]", "<|\"g\" -> 2, \"lc\" -> 1, \"count\" -> 1, \"r\" -> 3, \"unramified\" -> False|>");
    BR_EQ("brPlaceDataZero[TE1, cvE1]", "<|\"e\" -> 0, \"G0\" -> 1, \"count\" -> 2, \"r\" -> 1, \"unramified\" -> True|>");
    BR_EQ("brPlaceDataZero[T139, cv139]", "<|\"e\" -> 1, \"G0\" -> 1, \"count\" -> 1, \"r\" -> 3, \"unramified\" -> False|>");
    BR_EQ("brPlaceDataZero[T158, cv158][\"G0\"]", "x^2");
    BR_EQ("brRootFactorInfAt[TL1, cvL1, 3]", "{1, 0, 1/2, 0}");
    BR_EQ("brRootFactorInfAt[T147, cv147, 2]", "{1, 1/2 x, -1/8 x^2}");
    BR_EQ("brSValueZero[TE1, cvE1]", "0");
    BR_EQ("brSValueInf[TE1, cvE1]", "1/2");
    BR_EQ("brSValueInf[T158, cv158]", "2/3");
    BR_EQ("brSValueZero[T158, cv158]", "2/(3 x)");
}

static void test_special_reduce(void) {
    BR_EQ("brSpecialOrderZero[TE1, (1 + brT1)/brT1^2]", "-2");
    BR_EQ("brSpecialOrderZero[TE1, 0]", "Infinity");
    BR_EQ("brSpecialDegreeInf[TE1, (1 + brT1^3)/brT1^2]", "1");
    BR_EQ("brSpecialDegreeInf[TE1, 0]", "-Infinity");
    BR_EQ("brSpecialPoleQ[TE1, cvE1, 1/brT1]", "{True, False}");
    BR_EQ("brSpecialPoleQ[TE1, cvE1, 1/(brT1 + 1)]", "{False, False}");
    BR_EQ("brSpecialPoleQ[T158, cv158, 1]", "{False, True}");
    BR_EQ("brSpecialPoleQ[TE1, cvE1, brT1]", "{False, True}");
    /* [B] p.158: the special reduction finds the whole algebraic part 3x */
    BR_EQ("brSpecialReduce[T158, cv158, (3 brT1 + 5 x + 2 brT1 x)/(brT1 + x)]",
          "<|\"Algebraic\" -> 3 x, \"Remainder\" -> 0|>");
    BR_EQ("brSpecialReduce[T147, cv147, 1/x]", "<|\"Algebraic\" -> 0, \"Remainder\" -> 1/x|>");
    /* E^(x^2) Sqrt[E^x + 1]-type pole at infinity with no reduction: Int E^(x/2)... */
    br_run("TE3 = brTower[E^(x^2) Sqrt[E^(x^2) + 1], x]; True");
    br_run("cvE3 = brCurve[TE3, 2, brT1 + 1]; True");
    BR_THROWS("brSpecialReduce[TE3, cvE3, brT1]", "NonElementary", "no reduction at theta = infinity ([B] Prop 2.18)");
    BR_THROWS("brSpecialReduce[TE1, cvE1, 1/(x brT1)]", "NonElementary", "no reduction at theta = 0 ([B] Prop 2.18)");
    /* property: an exact derivative reduces to 0 at the special places */
    BR_TRUE("SeedRandom[73]; And @@ Table[Module[{b = (RandomInteger[{1, 3}] x + RandomInteger[{-2, 2}])"
            "   brT1^RandomInteger[{-2, 2}], c, r},"
            "   c = Together[brD[TE1, b] + cvE1[\"s\"] b]; r = brSpecialReduce[TE1, cvE1, c];"
            "   Together[r[\"Remainder\"]] === 0 && Together[r[\"Algebraic\"] - b] === 0], {20}]");
}

static void test_normal_residues(void) {
    /* [B] p.139: the resultant criterion certifies non-elementarity */
    BR_EQ("brResultantCriterion[T139, brResultantR[T139, cv139, 1/(1 + brT2), z], z][[1]]", "\"NonConstant\"");
    BR_THROWS("brNormalResidues[T139, cv139, 1/(1 + brT2)]", "NonElementary",
              "resultant criterion ([B] Prop 2.9): R(z) has non-constant coefficients");
    /* [B] p.147 (theta-form integrand), the y-component: residues +-1 */
    br_run("c147 = (brT1 + 3 brT1 x + x^2 - x^3)/(x (brT1^2 + 2 brT1 x + x^2 - brT1 x^2 - x^3)); True");
    BR_TRUE("Module[{r = brNormalResidues[T147, cv147, c147]}, r[\"P\"] === r[\"z\"]^2 - 1]");
    BR_TRUE("Module[{r = brNormalResidues[TL1, cvL1, 1/(x brT1)]}, r[\"P\"] === r[\"z\"]^2 - 1]");
    BR_EQ("brResultantR[TL1, cvL1, 1/((1 + brT1^2) x), z]", "1");
    BR_EQ("brResultantR[TL1, cvL1, brT1^2, z]", "1");
    BR_EQ("brResultantCriterion[TL1, 1, z]", "{\"Constant\", 1}");
    BR_EQ("brResultantCriterion[TL1, 3 x z^2 - 3 x, z]", "{\"Constant\", -1 + z^2}");
    BR_EQ("brResultantCriterion[TL1, z^2 - x, z][[1]]", "\"NonConstant\"");
    BR_EQ("brNonzeroResiduePolynomial[z^3 (z - 2), z]", "-2 + z");
    BR_EQ("brNonzeroResiduePolynomial[z^3, z]", "1");
    BR_EQ("brNonzeroResiduePolynomial[(z^2 - 2)^2 z, z]", "-2 + z^2");
    BR_EQ("brNonzeroResiduePolynomial[1, z]", "1");
    /* property: an exact derivative D(b Y) has no residue at a normal place */
    BR_TRUE("SeedRandom[74]; And @@ Table[Module[{b = (x + RandomInteger[{1, 3}] brT1)/(brT1 + RandomInteger[{1, 4}] x)^RandomInteger[{1, 2}],"
            "   c, hr, r},"
            "   c = Together[brD[T147, b] + cv147[\"s\"] b]; hr = brHermiteReduce[T147, cv147, c];"
            "   r = brNormalResidues[T147, cv147, hr[\"Remainder\"]]; r[\"P\"] === 1], {20}]");
}

static void test_special_residues(void) {
    BR_EQ("brGammaForm[4 x^2, 2]", "{4, x}");
    BR_EQ("brGammaForm[x, 2]", "None");
    BR_EQ("brGammaForm[3, 3]", "{3, 1}");
    BR_EQ("brCoefficientAtInf[TL1, cvL1, 1/((1 + brT1^2) x), -1]", "1/x");
    BR_EQ("brCoefficientAtInf[TL1, cvL1, brT1, 0]", "1/2");
    BR_EQ("brCoefficientAtInf[TL1, cvL1, 1, 0]", "0");
    BR_EQ("brCoefficientAtZero[TE1, cvE1, 1/(1 + brT1), 0]", "1");
    BR_EQ("brCoefficientAtZero[TE1, cvE1, brT1, 0]", "0");
    BR_EQ("brCoefficientAtZero[TE1, cvE1, (2 + brT1)/(1 + brT1), 1]", "0");
    /* log: d/dx Log[theta + Sqrt[theta^2 + 1]]: lambda = 1 at infinity; Sqrt[theta^2 + 1]/(x theta): 0 */
    BR_EQ("brLogInfinityResidues[TL1, cvL1, 1/((1 + brT1^2) x)]",
          "<|\"ramified\" -> False, \"kappa\" -> 1, \"u\" -> 1, \"lambda\" -> 1|>");
    BR_EQ("brLogInfinityResidues[TL1, cvL1, 1/(x brT1)][\"lambda\"]", "0");
    BR_EQ("brLogInfinityResidues[T147, cv147, c147]", "<|\"ramified\" -> True, \"tau\" -> 0|>");
    BR_THROWS("brLogInfinityResidues[TL1, cvL1, 1/((1 + brT1^2) (x^2 + 1))]", "Decline",
              "residue at infinity: integral in k with a non-rational logarithm");
    br_run("TL2 = brTower[Sqrt[Log[x]^2 + x], x]; True");
    br_run("cvL2 = brCurve[TL2, 2, brT1^2 + x]; True");
    BR_THROWS("brLogInfinityResidues[TL2, cvL2, 1/((brT1^2 + x) (x + 1))]", "NonElementary",
              "residue at infinity undefined ([B] Def 2.10, Prop 2.11)");
    /* exp: values and their integrals */
    BR_EQ("brExpValueIntegral[TE1, 1]", "{\"Solution\", 1}");
    BR_EQ("brExpValueIntegral[TE1, 0]", "{\"Solution\", 0}");
    BR_EQ("brExpValueIntegral[TE1, 1/x]", "{\"Solution\", 0}");
    BR_EQ("brExpValueIntegral[TE1, x]", "{\"NoSolution\"}");
    BR_EQ("brExpSpecialResidues[TE1, cvE1, 1/(1 + brT1)]",
          "<|\"zero\" -> <|\"ramified\" -> False, \"kappa\" -> 1, \"u\" -> 1, \"c\" -> 1, \"w\" -> 1|>, \"inf\" -> <|\"ramified\" -> True|>|>");
    BR_EQ("brExpSpecialResidues[TE2, cvE2, brT1/(1 + brT1^2)][\"inf\"][\"c\"]", "1/2");
    BR_THROWS("brExpSpecialResidues[TE1, cvE1, x/(1 + brT1)]", "NonElementary",
              "residue at a special place undefined ([B] Def 2.21, Prop 2.22)");
    BR_TRUE("brSpecialResiduesZeroQ[T147, <|\"ramified\" -> True, \"tau\" -> 0|>]");
    BR_FALSE("brSpecialResiduesZeroQ[TL1, <|\"ramified\" -> False, \"lambda\" -> 1|>]");
    BR_TRUE("brSpecialResiduesZeroQ[TE1, <|\"zero\" -> <|\"ramified\" -> True|>, \"inf\" -> <|\"ramified\" -> False, \"c\" -> 0|>|>]");
    BR_FALSE("brSpecialResiduesZeroQ[TE1, brExpSpecialResidues[TE1, cvE1, 1/(1 + brT1)]]");
}

int main(void) {
    br_init();
    setup();
    test_series();
    test_place_data();
    test_special_reduce();
    test_normal_residues();
    test_special_residues();
    BR_DONE("bronstein_residues_tests");
    return 0;
}
