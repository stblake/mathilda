/* test_bronstein_hermite_special.c -- Hermite reduction at the normal places
 * (Hermite.m), integration in k (LimitedIntegrate.m) and the algebraic part
 * D b + s b = c (AlgebraicPart.m), BRONSTEIN_PLAN.md Phase 3.
 *
 * Pinned: [B] p.134 closes by Hermite reduction alone (b = -2/(theta + x),
 * remainder 0); [B] p.158's algebraic part is 3x; the special-place quantity of
 * [B] p.147 integrates to x - x^2 + 1*Log[x] (lambda = 1, tau = -2).
 * Properties over seeded families: the Hermite identity c = D b + s b + r with
 * a squarefree normal remainder; the algebraic part recovers the unique b from
 * c = D b + s b (C path and the package's own log-case solver separately);
 * limited integration recovers lambda.
 */
#include "bronstein_test_utils.h"

static void setup(void) {
    /* The towers first: they intern Bronstein`Private`brT1, ...  A statement
     * that mentions brT1 must be PARSED after that, or the parser creates
     * Global`brT1 (a different symbol, earlier on $ContextPath). */
    br_run("T134 = brTower[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x];"
           " T158 = brTower[(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3)), x];"
           " TE = brTower[E^x Sqrt[Log[x]], x]; TS = brTower[Sqrt[Log[x]], x]; TX = brTower[Sqrt[x Log[x] + 1], x];"
           " True");
    BR_EQ("Context[brT1]", "\"Bronstein`Private`\"");
    br_run("cv134 = brCurve[T134, 2, brT1 + x]; cv158 = brCurve[T158, 3, (brT1 + x)^2];"
           " TL134 = brLowerTower[T134]; TLE = brLowerTower[TE];"
           " cvS = brCurve[TS, 2, brT1]; cvX = brCurve[TX, 2, x brT1 + 1]; True");
}

static void test_tower_additions(void) {
    BR_EQ("brFromSurface[T134, Log[x]^2 + x]", "brT1^2 + x");
    BR_EQ("brFromSurface[T158, E^(2 x) + E^(x + 1)]", "brT1^2 + brT1 E");
    BR_EQ("brFromSurface[T134, Log[x + 1]]", "Log[1 + x]");
    BR_EQ("TL134[\"syms\"]", "{}");
    BR_EQ("{TLE[\"syms\"], TLE[\"kernels\"], TLE[\"shape\"]}", "{{brT1}, {E^x}, \"Lower\"}");
    BR_EQ("brToSurface[TLE, brT1^2]", "E^(2 x)");
}

static void test_hermite(void) {
    BR_EQ("brCollectTheta[T134, (x brT1^2 + brT1)/x]", "brT1^2 + brT1/x");
    BR_EQ("brCollectTheta[T134, 0]", "0");
    BR_THROWS("brCollectTheta[T134, 1/brT1]", "Internal", "brCollectThetaNotPolynomial");
    BR_EQ("brNormalMultipleFactors[T134, x (brT1 + x)^2 (brT1 - 1)]", "{{brT1 + x, 2}}");
    BR_EQ("brNormalMultipleFactors[T158, brT1^3 (brT1 + x)]", "{}");
    BR_EQ("brNormalMultipleFactors[T158, brT1^3 (brT1 + x)^3]", "{{brT1 + x, 3}}");
    /* the [B] p.134 congruence: 0 == -(x+1) B1, x + 1 == -(x+1)/2 B2 (mod theta + x) */
    BR_EQ("brHermiteSolve[T134, x + 1, Together[(brT1 + x) cv134[\"s\"] - brD[T134, brT1 + x]], brT1 + x]", "-2 x");
    /* with U = x (D0 = x (theta + x)^2), the paper's B2 = -2 */
    BR_EQ("brHermiteSolve[T134, x + 1, Together[x ((brT1 + x) cv134[\"s\"] - brD[T134, brT1 + x])], brT1 + x]", "-2");
    BR_THROWS("brHermiteSolve[T134, 1, brT1 + x, brT1 + x]", "Internal", "brHermiteNonInvertibleCongruence");
    BR_THROWS("brHermiteSolve[T134, 1, 1/(brT1 + x), brT1 + x]", "Internal", "brHermiteDenominatorNotCoprime");
    BR_EQ("brHermiteReduce[T134, cv134, (1 + x)/(x (brT1 + x)^2)]",
          "<|\"Algebraic\" -> -2/(brT1 + x), \"Remainder\" -> 0|>");
    BR_EQ("brHermiteReduce[T134, cv134, 1/(brT1 - 1)]",
          "<|\"Algebraic\" -> 0, \"Remainder\" -> 1/(-1 + brT1)|>");
    BR_TRUE("brNormalDenominatorSquarefreeQ[T134, 1/(brT1 - 1)]");
    BR_FALSE("brNormalDenominatorSquarefreeQ[T134, 1/(brT1 - 1)^2]");
    BR_TRUE("brNormalDenominatorSquarefreeQ[T158, 1/brT1^3]");
    /* exponential: the special power of theta is not reduced */
    BR_EQ("brHermiteReduce[T158, cv158, 1/brT1^2][\"Algebraic\"]", "0");
    /* property: c = D b + s b + r, r normal-squarefree (log and exp) */
    BR_TRUE("SeedRandom[61]; And @@ Table[Module[{T = RandomChoice[{{T134, cv134}, {T158, cv158}}], V, c, hr},"
            "   V = brT1 + RandomInteger[{1, 3}] x + RandomInteger[{-2, 2}];"
            "   c = (RandomInteger[{1, 3}] brT1 + x)/(V^RandomInteger[{2, 4}] (brT1 + 2)^RandomInteger[{1, 2}]);"
            "   hr = brHermiteReduce[T[[1]], T[[2]], c];"
            "   Together[brD[T[[1]], hr[\"Algebraic\"]] + T[[2]][\"s\"] hr[\"Algebraic\"] + hr[\"Remainder\"] - c] === 0 &&"
            "   brNormalDenominatorSquarefreeQ[T[[1]], hr[\"Remainder\"]]], {20}]");
}

static void test_integrate_in_k(void) {
    BR_TRUE("MatchQ[ExpIntegralEi[x], _?(!FreeQ[#, brNonElementaryHeads[]] &)]");
    BR_EQ("brIntegrateInK[TL134, -2 x + 1 + 1/x]", "{\"Elementary\", x - x^2 + Log[x]}");
    BR_EQ("brIntegrateInK[TL134, 0]", "{\"Elementary\", 0}");
    BR_EQ("brIntegrateInK[TLE, brT1/x]", "{\"NonElementary\"}");
    BR_EQ("brIntegrateInK[TLE, brT1^2]", "{\"Elementary\", 1/2 brT1^2}");
    BR_EQ("brIntegrateInK[TL134, 1/(x^2 + 1)][[1]]", "\"Unknown\"");
    BR_EQ("brSplitLogs[TL134, x^2 + 3 Log[x] - Log[x + 1]/2]", "{x^2, {{3, x}, {-1/2, 1 + x}}}");
    BR_EQ("brSplitLogs[TL134, x^2]", "{x^2, {}}");
    BR_EQ("brSplitLogs[TL134, Log[x]^2]", "None");
    BR_EQ("brSplitLogs[TL134, x Log[x]]", "None");
    /* basis reduction: Log[x^2] and Log[2 x] reduce to Log[x] at the level of derivatives */
    BR_EQ("brLogBasisReduce[TL134, {x, x^2, 2 x, x + 1}][[1]]", "{x, 1 + x}");
    BR_EQ("brLogBasisReduce[TL134, {x, x^2, 2 x, x + 1}][[2]]", "{{1, 0}, {2, 0}, {1, 0}, {0, 1}}");
    BR_EQ("brLogBasisReduce[TLE, {brT1 + 1}][[1]]", "{1 + brT1}");
    BR_EQ("brLogBasisReduce[TLE, {brT1^3}]", "{{}, {{}}, {3 x}}");
    BR_EQ("brLogBasisReduce[TL134, {}]", "{{}, {}, {}}");
}

static void test_param_limited(void) {
    BR_EQ("brParamIntegrate[TL134, 1/x, {1/(x + 1)}]", "{\"NoSolution\"}");
    /* logs Log x (1 + kappa) and Log(x + 1) (2): 2 == 0 is inconsistent */
    BR_EQ("brParamIntegrate[TL134, 1/x + 2/(x + 1), {1/x}]", "{\"NoSolution\"}");
    BR_EQ("Module[{r = brParamIntegrate[TL134, 2 x + 1/x, {1/x}]}, Solve[r[[4]], r[[5]]]]", "{{brKappa1 -> -1}}");
    BR_EQ("brParamIntegrate[TL134, 2 x, {}]", "{\"Solution\", x^2, {}, {}, {}}");
    BR_EQ("brParamIntegrate[TLE, brT1/x, {1}]", "{\"NoSolution\"}");
    BR_EQ("brParamIntegrate[TLE, 1, {brT1/x}][[1]]", "\"Unknown\"");
    BR_EQ("brLimitedIntegrate[TL134, -2 x + 1 + 1/x, 1/x]", "{\"Solution\", x - x^2, 1}");
    BR_EQ("brLimitedIntegrate[TL134, 3 x^2, 1/x]", "{\"Solution\", x^3, 0}");
    BR_EQ("brLimitedIntegrate[TL134, 1/(x + 1), 1/x]", "{\"NoSolution\"}");
    BR_EQ("brLimitedIntegrate[TLE, brT1/x, 1/x]", "{\"NoSolution\"}");
    /* property: a = v' + lambda w recovers lambda and v up to a constant */
    BR_TRUE("SeedRandom[62]; And @@ Table[Module[{v = (RandomInteger[{-3, 3}] + x^2)/(x + RandomInteger[{1, 4}]),"
            "   lam = RandomInteger[{-4, 4}], w = 1/(x + 5), r},"
            "   r = brLimitedIntegrate[TL134, Together[D[v, x] + lam w], w];"
            "   r[[1]] === \"Solution\" && r[[3]] === lam && FreeQ[Together[D[r[[2]] - v, x]], x] &&"
            "   Together[D[r[[2]] - v, x]] === 0], {20}]");
}

static void test_algebraic_part(void) {
    BR_EQ("brMthPowerCofactor[4 x^2, 2]", "x");
    BR_EQ("brMthPowerCofactor[x, 2]", "None");
    BR_EQ("brMthPowerCofactor[3, 2]", "1");
    BR_EQ("brMthPowerCofactor[x^3/(x + 1)^3, 3]", "x/(1 + x)");
    BR_EQ("brMthPowerCofactor[brT1^2 x^4, 2]", "brT1 x^2");
    BR_EQ("brIdenticallyZero[a x + b, {x}]", "{a == 0, b == 0}");
    BR_EQ("brIdenticallyZero[0, {x}]", "{}");
    BR_EQ("brIdenticallyZero[(a - 1)/x, {x}]", "{-1 + a == 0}");
    BR_TRUE("brAlgebraicResidualZeroQ[T158, cv158[\"s\"], (3 brT1 + 5 x + 2 brT1 x)/(brT1 + x), 3 x]");
    BR_FALSE("brAlgebraicResidualZeroQ[T158, cv158[\"s\"], (3 brT1 + 5 x + 2 brT1 x)/(brT1 + x), 2 x]");
    BR_TRUE("brAlgebraicResidualZeroQ[T134, cv134[\"s\"], 0, 0]");
    BR_EQ("brKSolveRDE[TLE, 1, brT1 (x + 1)]", "{\"Solution\", 1/4 (brT1 + 2 brT1 x)}");
    BR_EQ("brKSolveRDE[TL134, 0, 1/x]", "{\"NoSolution\"}");
    BR_EQ("brKSolveRDE[TL134, 1, x]", "{\"Solution\", -1 + x}");
    /* [B] p.158 algebraic part; [B] p.139 component and Sqrt[Log x]: proven none */
    BR_EQ("brAlgebraicPart[T158, cv158, (3 brT1 + 5 x + 2 brT1 x)/(brT1 + x)]", "{\"Solution\", 3 x}");
    BR_EQ("brAlgebraicPart[T134, cv134, 0]", "{\"Solution\", 0}");
    BR_EQ("brAlgebraicPart[TS, cvS, 1]", "{\"NoSolution\"}");
    BR_EQ("brPolyRDELog[TS, cvS, 1]", "{\"NoSolution\"}");
    BR_EQ("brPolyRDELog[T134, cv134, 1/x]", "{\"NoSolution\"}");
    /* the own log-case solver recovers b = x theta, 1 + theta^2 */
    BR_EQ("brPolyRDELog[T134, cv134, Together[brD[T134, x brT1] + cv134[\"s\"] x brT1]]", "{\"Solution\", brT1 x}");
    BR_EQ("brPolyRDELog[TS, cvS, Together[brD[TS, 1 + brT1^2] + cvS[\"s\"] (1 + brT1^2)]]", "{\"Solution\", 1 + brT1^2}");
    /* property: c = D b + s b for a random polynomial b; both paths recover b */
    BR_TRUE("SeedRandom[63]; And @@ Table[Module[{b = Sum[RandomInteger[{-3, 3}] x^RandomInteger[{0, 2}] brT1^i,"
            "   {i, 0, RandomInteger[{0, 2}]}], c, r1, r2},"
            "   c = Together[brD[T134, b] + cv134[\"s\"] b];"
            "   r1 = brAlgebraicPart[T134, cv134, c]; r2 = brPolyRDELog[T134, cv134, c];"
            "   r1[[1]] === \"Solution\" && Together[r1[[2]] - b] === 0 &&"
            "   r2[[1]] === \"Solution\" && Together[r2[[2]] - b] === 0], {20}]");
    /* a non-constant leading coefficient (beta != 0, no homogeneous solution in k) */
    BR_TRUE("SeedRandom[64]; And @@ Table[Module[{b = Sum[RandomInteger[{-3, 3}] (x + RandomInteger[{0, 2}]) brT1^i,"
            "   {i, 0, RandomInteger[{0, 2}]}], c, r},"
            "   c = Together[brD[TX, b] + cvX[\"s\"] b]; r = brPolyRDELog[TX, cvX, c];"
            "   r[[1]] === \"Solution\" && Together[r[[2]] - b] === 0], {20}]");
    /* exponential: the C path recovers b in k(theta) */
    BR_TRUE("SeedRandom[65]; And @@ Table[Module[{b = (RandomInteger[{1, 3}] x + RandomInteger[{-2, 2}] brT1)/"
            "   (1 + RandomInteger[{0, 1}] x), c, r},"
            "   c = Together[brD[T158, b] + cv158[\"s\"] b]; r = brAlgebraicPart[T158, cv158, c];"
            "   r[[1]] === \"Solution\" && Together[r[[2]] - b] === 0], {20}]");
}

/* Temporary end-to-end driver (BRONSTEIN_PLAN.md Phase 3): [B] p.134 and
 * p.158 assembled from the coordinate-0 integral and the components'
 * algebraic parts, verified by PossibleZeroQ against the integrand. */
static void test_end_to_end(void) {
    br_run("brTmpIntegrate[f_] := Module[{T = brTower[f, x], cf, tot = 0, cc, cv, hr, ap},"
           "  cf = brCanonicalForm[T];"
           "  tot = Integrate`RischTranscendental[brToSurface[T, cf[\"coeffs\"][[1]]], x];"
           "  Do[If[cf[\"coeffs\"][[i + 1]] =!= 0,"
           "    cc = brComponentCurve[T, cf[\"n\"], cf[\"g\"], i]; cv = brCurve[T, cc[\"m\"], cc[\"G\"]];"
           "    hr = brHermiteReduce[T, cv, Together[cf[\"coeffs\"][[i + 1]] cc[\"s\"]]];"
           "    ap = brAlgebraicPart[T, cv, hr[\"Remainder\"]];"
           "    tot += brToSurface[T, (hr[\"Algebraic\"] + ap[[2]]) (cf[\"y\"]^i/cc[\"s\"]) /. cf[\"y\"] -> cf[\"yexpr\"]]],"
           "   {i, 1, cf[\"n\"] - 1}]; tot]; True");
    BR_EQ("brTmpIntegrate[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]])]", "-2/Sqrt[x + Log[x]]");
    BR_TRUE("Module[{f = (3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3)), F},"
            " F = brTmpIntegrate[f]; PossibleZeroQ[F - (3 x (x + E^x)^(2/3) + 3 Log[x])] &&"
            " PossibleZeroQ[D[F, x] - f]]");
}

int main(void) {
    br_init();
    setup();
    test_tower_additions();
    test_hermite();
    test_integrate_in_k();
    test_param_limited();
    test_algebraic_part();
    test_end_to_end();
    BR_DONE("bronstein_hermite_special_tests");
    return 0;
}
