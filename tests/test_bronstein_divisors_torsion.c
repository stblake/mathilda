/* test_bronstein_divisors_torsion.c -- residue divisors over k(t)
 * (Divisors.m), torsion orders by specialization (Specialize.m) and the log
 * part with its weight-1 projection (LogPart.m), BRONSTEIN_PLAN.md Phases 5-6.
 *
 * Pinned: [B] p.147's component (the paper's theta-form integrand) has one log
 * term whose generator is a k-multiple of x + y, remainder (1+x)/(x(theta+x));
 * divisors at the places above infinity (arcsinh family, genus 0); genus-1
 * certificates lifted from Trager's Nagell-Lutz examples (theta = Log x);
 * a genus-1 torsion divisor of order 3 over Q(x) (y^2 = theta^3 + x^2).
 * Properties: for random logands w, the component of D log w has a log part
 * whose derivative's coordinate 1 equals the component (remainder 0).
 */
#include "bronstein_test_utils.h"

static void setup(void) {
    br_run("T147 = brTower[Sqrt[x + Log[x]], x]; TL1 = brTower[Sqrt[Log[x]^2 + 1], x];"
           " TE1 = brTower[Sqrt[E^x + 1], x]; TE2 = brTower[Sqrt[E^(2 x) + 1], x];"
           " TG = brTower[Sqrt[Log[x]^3 - 2], x]; TQ = brTower[Sqrt[Log[x]^3 + x^2], x]; True");
    br_run("cv147 = brCurve[T147, 2, brT1 + x]; cvL1 = brCurve[TL1, 2, brT1^2 + 1];"
           " cvE1 = brCurve[TE1, 2, brT1 + 1]; cvE2 = brCurve[TE2, 2, brT1^2 + 1];"
           " cvG = brCurve[TG, 2, brT1^3 - 2]; cvQ = brCurve[TQ, 2, brT1^3 + x^2];"
           " c147 = (brT1 + 3 brT1 x + x^2 - x^3)/(x (brT1^2 + 2 brT1 x + x^2 - brT1 x^2 - x^3));"
           " nr147 = brNormalResidues[T147, cv147, c147]; rs147 = brResidueSetup[T147, cv147, c147, nr147]; True");
}

static void test_setup_and_ideals(void) {
    BR_EQ("rs147[\"F0\"]", "<|\"Type\" -> \"Q\"|>");
    BR_EQ("Sort[rs147[\"roots\"]]", "{-1, 1}");
    BR_EQ("{Length[rs147[\"basis\"]], Abs[rs147[\"coeffs\"]]}", "{1, {{1}, {1}}}");
    BR_EQ("rs147[\"Bp\"]", "brT1 x + x^2 - x^3");   /* x (theta - x^2 + x): the x-content stays */
    BR_EQ("rs147[\"Fk\"][\"vars\"]", "{x}");
    BR_EQ("brPlaceDegree[rs147, 1]", "1");
    BR_EQ("brPlaceDegree[rs147, -1]", "1");
    BR_EQ("brPlaceDegree[rs147, 2]", "0");
    BR_TRUE("Module[{I = brDivisorIdeal[rs147, {1, 0}]}, KeyExistsQ[I, \"M\"] && Length[I[\"M\"]] === 2]");
    BR_TRUE("Module[{I = brDivisorIdeal[rs147, {0, 0}]}, I === idealUnit[rs147[\"B\"]]]");
    /* the residue divisor of p.147 is principal: a generator, k-multiple of x + y */
    BR_TRUE("Module[{u = brPrincipalGenerator[rs147, rs147[\"coeffs\"][[All, 1]], 1]},"
            " u =!= None && Together[u[[2]]/u[[1]] - 2 x/(x^2 + brT1 + x)] === 0]");   /* (x + y)^2 up to k */
    BR_EQ("brPrincipalGenerator[rs147, {1, 1}, 1]", "None");   /* degree 2: not principal */
}

static void test_specialization(void) {
    BR_EQ("brSpecializationPoint[{x}, 1]", "{x -> 64/5}");
    BR_TRUE("Length[brSpecializationPoint[{x, y}, 3]] === 2 &&"
            " brSpecializationPoint[{x, y}, 3] =!= brSpecializationPoint[{x, y}, 4]");
    BR_EQ("brSquarefreeProfile[(brT1 + x)^2 (brT1 - 1) x, brT1]", "{{1, 1}, {1, 2}}");
    BR_EQ("brSquarefreeProfile[x, brT1]", "{}");
    BR_EQ("brSpecializeValue[(x + 1)/(x - 2), {x -> 3}]", "4");
    BR_EQ("brSpecializeValue[(x + 1)/(x - 2), {x -> 2}]", "$Failed");
    BR_TRUE("brGoodSpecializationQ[T147, cv147, rs147, {x -> 5}]");
    BR_FALSE("brGoodSpecializationQ[T147, cv147, rs147, {x -> 0}]");
    BR_FALSE("brGoodSpecializationQ[T147, cv147, rs147, {x -> 1}]");  /* Bp loses its theta-term */
    BR_TRUE("brGoodSpecializationQ[T147, cv147, rs147, brGoodSpecialization[T147, cv147, rs147]]");
    BR_TRUE("Module[{ctx = brSpecializedContext[T147, cv147, rs147, {x -> 5}]},"
            " ctx[\"genus\"] === 0 && FreeQ[{ctx[\"N\"], ctx[\"Bp\"]}, x] && ctx[\"F\"] === rs147[\"F0\"]]");
    BR_EQ("brTorsionOrder[T147, cv147, rs147, {1, -1}]", "1");
}

static void test_genus_one(void) {
    /* lifted Nagell-Lutz: (3, 5) on y^2 = theta^3 - 2 has infinite order */
    BR_THROWS("brComponentLogPart[TG, cvG, 5/(x (brT1 - 3))]", "NonElementary",
              "residue divisor of infinite order: its orders modulo 11 and 17 disagree");
    /* y^2 = theta^3 + x^2, P = (0, x), Q = (0, -x): P - Q has order 3 over Q(x) */
    /* D Log[(y - x)/(y + x)] = (3 - 2 theta) y/(theta (theta^3 + x^2)): order 3 */
    br_run("cQ = (3 - 2 brT1)/(brT1 (brT1^3 + x^2)); rsQ = brResidueSetup[TQ, cvQ, cQ, brNormalResidues[TQ, cvQ, cQ]]; True");
    BR_EQ("brTorsionOrder[TQ, cvQ, rsQ, rsQ[\"coeffs\"][[All, 1]]]", "3");
    BR_TRUE("Module[{lp = brComponentLogPart[TQ, cvQ, cQ]}, Together[lp[\"remainder\"]] === 0]");
}

static void test_special_divisors(void) {
    BR_EQ("brSpecialLocations[TL1, cvL1, <|\"ramified\" -> False, \"kappa\" -> 1, \"u\" -> 1, \"lambda\" -> 1|>]",
          "{<|\"where\" -> \"inf\", \"ramified\" -> False, \"kappa\" -> 1, \"u\" -> 1, \"coef\" -> 1|>}");
    BR_EQ("brSpecialLocations[TL1, cvL1, <|\"ramified\" -> True, \"tau\" -> 0|>]", "{}");
    BR_EQ("Length[brSpecialLocations[TE1, cvE1, brExpSpecialResidues[TE1, cvE1, 1/(1 + brT1)]]]", "2");
    br_run("frL = brShiftFrame[TL1, cvL1, <|\"N\" -> 1, \"Ep\" -> 1/x, \"Bp\" -> 1|>, 1, uu]; True");
    BR_EQ("{frL[\"k\"], frL[\"Bu\"][\"G\"]}", "{1, 1 + 2 uu + 2 uu^2}");
    BR_TRUE("Module[{I = brPointPlaceIdeal[fieldFn[{x}], frL[\"Bu\"], uu, 1]}, Length[I[\"M\"]] === 2]");
    BR_THROWS("brPointPlaceIdeal[fieldFn[{x}], frL[\"Bu\"], uu, 2]", "Internal", "brPointPlaceNotOnCurve");
    BR_EQ("Length[brSpecialPlaceIdeals[fieldFn[{x}], frL, <|\"where\" -> \"inf\", \"ramified\" -> False, \"u\" -> 1|>, {1, -1}]]", "2");
    BR_EQ("brSpecialTaus[TL1, cvL1, {<|\"where\" -> \"inf\", \"ramified\" -> False, \"coef\" -> 1|>}, {{1, -1}}]", "{-1, 1}");
    BR_EQ("brSpecialTaus[TE1, cvE1, {<|\"where\" -> \"zero\", \"ramified\" -> False, \"coef\" -> 1|>,"
          " <|\"where\" -> \"inf\", \"ramified\" -> True, \"coef\" -> 0|>}, {{1, -1}, {}}]", "{0, -2, 2}");
    /* arcsinh family: the log part closes with remainder 0 */
    BR_TRUE("Together[brComponentLogPart[TL1, cvL1, 1/((1 + brT1^2) x)][\"remainder\"]] === 0");
    BR_TRUE("Together[brComponentLogPart[TE2, cvE2, brT1/(1 + brT1^2)][\"remainder\"]] === 0");
    BR_TRUE("Together[brComponentLogPart[TE1, cvE1, 1/(1 + brT1)][\"remainder\"]] === 0");
    BR_THROWS("brComponentLogPartSpecial[TG, cvG, 1/x, <|\"P\" -> 1, \"z\" -> z|>, <|\"ramified\" -> False, \"lambda\" -> 1|>]",
              "Decline", "residues at a special place on a curve of positive genus (not yet supported)");
}

static void test_log_part(void) {
    BR_EQ("brAfD[T147, cv147, fieldFn[{x}], cv147[\"B\"], {x, 0}]", "{1, 0}");
    BR_EQ("brAfD[T147, cv147, fieldFn[{x}], cv147[\"B\"], {0, 1}]", "{0, (1 + x)/(2 x (brT1 + x))}");
    BR_EQ("brLogDerivative[T147, cv147, fieldFn[{x}], cv147[\"B\"], {x, 0}]", "{1/x, 0}");
    BR_EQ("brLogPartDerivative[T147, cv147, fieldFn[{x}], cv147[\"B\"], {{2, {x, 0}}}]", "{2/x, 0}");
    BR_EQ("brLogPartDerivative[T147, cv147, fieldFn[{x}], cv147[\"B\"], {}]", "{0, 0}");
    /* [B] p.147 component: one term, remainder (1 + x)/(x (theta + x)) */
    BR_TRUE("Module[{lp = brComponentLogPart[T147, cv147, c147]},"
            " Length[lp[\"terms\"]] === 1 && Together[lp[\"remainder\"] - (1 + x)/(x (brT1 + x))] === 0]");
    BR_EQ("brComponentLogPart[T147, cv147, 0][\"terms\"]", "{}");
    BR_TRUE("Together[brComponentLogPart[T147, cv147, 1/x][\"remainder\"] - 1/x] === 0");
    BR_THROWS("brComponentLogPart[T147, cv147, 1/(x (brT1 + x + 1))]", "NonElementary",
              "resultant criterion ([B] Prop 2.9): R(z) has non-constant coefficients");
    BR_TRUE("Module[{lp = brComponentLogPart[T147, cv147, c147], F, B},"
            " F = lp[\"F\"]; B = lp[\"B\"];"
            " Head[brLogPartFinish[T147, cv147, c147, F, B, lp[\"terms\"]]] === Association]");
    /* brNormalizeGenerator: the k(theta)-content of the coordinates is divided out */
    BR_EQ("brNormalizeGenerator[fieldFn[{x}], cv147[\"B\"], {2 x (x + brT1), 2 x}, {x, brT1}]", "{brT1 + x, 1}");
    BR_EQ("brNormalizeGenerator[fieldFn[{x}], cv147[\"B\"], {(x + brT1)/x, 1/x}, {x, brT1}]", "{brT1 + x, 1}");
    BR_EQ("brNormalizeGenerator[fieldFn[{x}], cv147[\"B\"], {0, 3/x}, {x, brT1}]", "{0, 1}");
    BR_EQ("brNormalizeGenerator[fieldFn[{x}], cv147[\"B\"], {0, 0}, {x, brT1}]", "{0, 0}");
    BR_EQ("brNormalizeGenerator[fieldFn[{x}], cv147[\"B\"], {x + brT1, 1}, {x, brT1}]", "{brT1 + x, 1}");
    /* property: the ratio to the input is free of the curve (an element of k(theta)) and
     * the weight-1 log term's derivative is unchanged */
    BR_TRUE("SeedRandom[82]; And @@ Table[Module[{u = {RandomInteger[{1, 3}] x + RandomInteger[{-2, 2}] brT1 + 1, RandomInteger[{1, 2}]},"
            "   g = (RandomInteger[{1, 3}] x + brT1)/(x + RandomInteger[{1, 4}]), v},"
            "   v = brNormalizeGenerator[fieldFn[{x}], cv147[\"B\"], g u, {x, brT1}];"
            "   Together[v[[1]] u[[2]] - v[[2]] u[[1]]] === 0],"
            " {20}]");
    /* property: the y-component of D Log[p + q y] has a log part leaving remainder 0 */
    BR_TRUE("SeedRandom[81]; And @@ Table[Module[{p = RandomInteger[{1, 3}] x + RandomInteger[{-2, 2}] brT1 + RandomInteger[{1, 3}],"
            "   q = RandomInteger[{1, 2}], f, cY},"
            "   f = Together[D[Log[(p + q Sqrt[brT1 + x]) /. brT1 -> Log[x]], x] /. Log[x] -> brT1];"
            "   cY = Together[(f - (f /. Sqrt[brT1 + x] -> -Sqrt[brT1 + x]))/(2 Sqrt[brT1 + x])];"
            "   Together[brComponentLogPart[T147, cv147, brHermiteReduce[T147, cv147, cY][\"Remainder\"]][\"remainder\"]] === 0],"
            " {12}]");
}

int main(void) {
    br_init();
    setup();
    test_setup_and_ideals();
    test_specialization();
    test_genus_one();
    test_special_divisors();
    test_log_part();
    BR_DONE("bronstein_divisors_torsion_tests");
    return 0;
}
