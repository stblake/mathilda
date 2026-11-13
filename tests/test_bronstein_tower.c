/* test_bronstein_tower.c -- the front end (src/internal/bronstein/Tower.m) and
 * the k[theta] layer (Poly.m), BRONSTEIN_PLAN.md Phase 2.
 *
 * Every routine is called directly on hand-checked inputs, including its
 * decline/throw paths, and the tower is checked by two properties over a seeded
 * random family: the back substitution reproduces the integrand, and the
 * derivation of each monomial agrees with D of its surface kernel.
 */
#include "bronstein_test_utils.h"

static void test_helpers(void) {
    BR_TRUE("brCheckHeads[Log[x] Sqrt[x], x]");
    BR_THROWS("brCheckHeads[Sin[x] Log[x], x]", "Unimplemented", "trigonometric or hyperbolic function (TrigToExp front end over Q(i))");
    BR_THROWS("brCheckHeads[1.5 Log[x], x]", "Decline", "approximate number in the integrand");
    BR_THROWS("brCheckHeads[I Log[x], x]", "Unimplemented", "complex constant in the integrand (constant field Q(i))");
    BR_THROWS("brCheckHeads[ArcTan[x], x]", "Unimplemented", "trigonometric or hyperbolic function (TrigToExp front end over Q(i))");
    BR_THROWS("brCheckHeads[Erf[x] Log[x], x]", "Decline", "non-elementary function in the integrand");
    BR_THROWS("brCheckHeads[Abs[x] Sqrt[x], x]", "Decline", "non-elementary function in the integrand");
    BR_THROWS("brCheckHeads[Root[#^5 - # - x &, 1], x]", "Incomplete", "a general algebraic function (Root/RootSum) in the integrand");
    BR_TRUE("MatchQ[Log[x], brKernelPattern[]] && MatchQ[E^x, brKernelPattern[]]");
    BR_FALSE("MatchQ[x^2, brKernelPattern[]]");
    BR_TRUE("MatchQ[Sin[x], _?(!FreeQ[#, brUnsupportedHeads[]] &)]");
    BR_TRUE("brDependsQ[x + t, {t}]");
    BR_FALSE("brDependsQ[x^2, {t}]");
    BR_TRUE("brDependsQ[Log[x], {x, t}]");
    BR_EQ("brInnermostKernels[Log[Log[x]] + E^x, {x}]", "{E^x, Log[x]}");
    BR_EQ("brInnermostKernels[Log[t1 + 1] + Log[2], {x, t1}]", "{Log[1 + t1]}");
    BR_EQ("brInnermostKernels[x^2, {x}]", "{}");
    BR_TRUE("brRadicalQ[Sqrt[x + 1], {x}]");
    BR_FALSE("brRadicalQ[Sqrt[2] x, {x}]");
    BR_TRUE("brRadicalQ[(t + x)^(1/3), {t}]");
    BR_EQ("brRules[x, {}]", "{x -> 1}");
    BR_EQ("brRules[x, {<|\"sym\" -> t1, \"D\" -> 1/x|>}]", "{x -> 1, t1 -> 1/x}");
    BR_EQ("brMonoSpec[{<|\"sym\" -> t1, \"kind\" -> \"Log\", \"D\" -> 1/x|>}]", "{{t1, \"Log\", 1/x}}");
    BR_EQ("brNewSymbol[3]", "brT3");
    BR_EQ("Context[Evaluate[brNewSymbol[1]]]", "\"Bronstein`Private`\"");
    BR_TRUE("brRationalInQ[(x + t)/(x - 1), {x, t}]");
    BR_FALSE("brRationalInQ[Sqrt[x], {x}]");
    BR_FALSE("brRationalInQ[x + y, {x}]");
    BR_FALSE("brRationalInQ[Log[x], {x}]");
    BR_TRUE("brRationalInQ[3/7, {x}]");
}

static void test_monomials(void) {
    br_run("mL = First[brAddMonomial[x, {}, \"Log\", x + 1, Log[x + 1]]];"
           " mE = First[brAddMonomial[x, {}, \"Exp\", x^2, E^(x^2)]]; True");
    BR_EQ("mL[[1]][\"D\"]", "1/(1 + x)");
    BR_EQ("mE[[1]][\"D\"]", "2 brT1 x");
    BR_EQ("Last[brAddMonomial[x, mL, \"Log\", brT1, Log[Log[x + 1]]]]", "brT2");
    BR_EQ("First[brAddMonomial[x, mL, \"Log\", brT1, Log[Log[x + 1]]]][[2]][\"D\"]", "1/(brT1 (1 + x))");
    BR_THROWS("brAddMonomial[x, {}, \"Log\", 3, Log[3]]", "Internal", "brAddMonomialConstant");
    /* exponential relations */
    br_run("mX = First[brAddMonomial[x, {}, \"Exp\", x, E^x]]; True");
    BR_EQ("brExpRewrite[x, {}, x]", "{\"new\"}");
    BR_EQ("brExpRewrite[x, mX, 2 x + 1]", "{\"rewrite\", brT1^2 E}");
    BR_EQ("brExpRewrite[x, mX, x^2]", "{\"new\"}");
    BR_EQ("brExpRewrite[x, mX, x/3]", "{\"rebase\", 1, 3}");
    BR_EQ("brExpRewrite[x, mL, 2 brT1]", "{\"rewrite\", (1 + x)^2}");
    /* E^(Log[x + 1]/2) is the radical Sqrt[x + 1] (pointwise) */
    BR_EQ("brExpRewrite[x, mL, brT1/2]", "{\"rewrite\", Sqrt[1 + x]}");
    /* re-basing E^x as (E^(x/3))^3 */
    br_run("rb = brRebaseExp[x, Join[mX, First[brAddMonomial[x, mX, \"Log\", 1 + brT1, Log[1 + E^x]]][[{2}]]], 1, 3]; True");
    BR_EQ("rb[[2]]", "brT1 -> brT1r3^3");
    BR_EQ("rb[[1]][[1]][\"D\"]", "1/3 brT1r3");
    BR_EQ("rb[[1]][[1]][\"kernel\"]", "E^(1/3 x)");
    BR_EQ("rb[[1]][[2]][\"arg\"]", "1 + brT1r3^3");
    BR_TRUE("Together[rb[[1]][[2]][\"D\"] - 3 brT1r3^2 (1/3 brT1r3)/(1 + brT1r3^3)] === 0");
    BR_EQ("brBack[mL]", "{brT1 -> Log[1 + x]}");
}

static void test_kernels_to_symbols(void) {
    BR_EQ("brKernelsToSymbols[Log[x] + E^x, x][[1]]", "brT1 + brT2");
    BR_EQ("#[\"kernel\"] & /@ brKernelsToSymbols[Log[Log[x]], x][[2]]", "{Log[x], Log[Log[x]]}");
    BR_EQ("brKernelsToSymbols[E^(2 x) + E^x, x][[1]]", "brT1 + brT1^2");
    BR_EQ("brKernelsToSymbols[E^(x/2) + E^x, x][[1]]", "brT1 + brT1^2");
    BR_EQ("brKernelsToSymbols[x^2, x]", "{x^2, {}}");
    BR_THROWS("brKernelsToSymbols[Log[x] + Log[x^2], x]", "Unimplemented", "reducible logarithm (rewrite with an explicit constant)");
    /* a radical inside a kernel: admitted as the single curve monomial (Shape B,
     * test_bronstein_shapeb.c), declined over the tower */
    BR_THROWS("brKernelsToSymbols[Log[Sqrt[Log[x]] + 1], x]", "Unimplemented", "radical inside a logarithm or exponential above other monomials (transcendental Risch over an algebraic base)");
}

static void test_detect_radical(void) {
    BR_EQ("brDetectRadical[Sqrt[t + x] + (t + x)^(1/3), {x, t}]", "<|\"b\" -> t + x, \"n\" -> 6|>");
    BR_EQ("brDetectRadical[1/Sqrt[1/(t + x)] + Sqrt[t + x], {x, t}]", "<|\"b\" -> 1/(t + x), \"n\" -> 2|>");
    BR_EQ("brDetectRadical[x + t, {x, t}]", "None");
    BR_THROWS("brDetectRadical[Sqrt[x + Sqrt[t]], {x, t}]", "Incomplete", "nested radicals");
    BR_THROWS("brDetectRadical[Sqrt[x] + Sqrt[t], {x, t}]", "Incomplete", "several radicals generating a non-cyclic Kummer extension (not one simple radical)");
    BR_THROWS("brDetectRadical[Sqrt[x^3] + Sqrt[x], {x, t}]", "Unimplemented", "several radicands generating one simple radical extension up to a locally constant unit or a constant field extension");
    BR_THROWS("brDetectRadical[Sqrt[x] + (t + 1)^(1/3), {x, t}]", "Unimplemented", "several radicands generating one simple radical extension up to a locally constant unit or a constant field extension");
    BR_THROWS("brDetectRadical[Sqrt[2] Sqrt[x], {x, t}]", "Unimplemented", "algebraic constant in the integrand (a number-field constant field)");
}

/* brTower on the paper examples and on every shape and decline. */
static void test_tower(void) {
    br_run("T134 = brTower[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x];"
           " T139 = brTower[Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), x];"
           " T158 = brTower[(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3)), x]; True");
    BR_EQ("{T134[\"shape\"], T134[\"theta\"], T134[\"kind\"], T134[\"rules\"]}", "{\"A\", brT1, \"Log\", {x -> 1, brT1 -> 1/x}}");
    BR_EQ("T134[\"radical\"]", "<|\"b\" -> brT1 + x, \"n\" -> 2|>");
    BR_EQ("{T139[\"theta\"], T139[\"lower\"], T139[\"rules\"]}", "{brT2, {brT1}, {x -> 1, brT1 -> brT1, brT2 -> brT1/(1 + brT1)}}");
    BR_EQ("{T158[\"theta\"], T158[\"kind\"], T158[\"radical\"][\"n\"]}", "{brT1, \"Exp\", 3}");
    BR_EQ("brTower[x^2 + 1/x, x][\"shape\"]", "\"Rational\"");
    BR_EQ("brTower[E^(x + 1) Log[x], x][\"shape\"]", "\"Transcendental\"");
    BR_EQ("brTower[Sqrt[x^3 + 1]/x, x][\"shape\"]", "\"Algebraic\"");
    BR_EQ("brTower[Log[x] Sqrt[x^2 + 1], x][\"shape\"]", "\"B\"");
    /* theta is the highest monomial of the radicand and moves to the top */
    BR_EQ("With[{T = brTower[E^x Sqrt[Log[x]], x]}, {T[\"theta\"], T[\"lower\"], T[\"kernels\"]}]",
          "{brT2, {brT1}, {E^x, Log[x]}}");
    BR_EQ("With[{T = brTower[Sqrt[Log[x]] + E^x, x]}, {T[\"theta\"], T[\"kernels\"]}]",
          "{brT2, {E^x, Log[x]}}");
    BR_THROWS("brTower[Sqrt[x + Log[x]] Log[1 + Log[x]], x]", "Unimplemented", "a monomial above the radicand's top monomial depends on it (radical in the middle of the tower)");
    BR_THROWS("brTower[Log[x] + Sqrt[Log[x] + Sqrt[Log[x]]], x]", "Incomplete", "nested radicals");
    BR_THROWS("brTower[Sqrt[x] Sqrt[x + 1] Log[x], x]", "Incomplete", "several radicals generating a non-cyclic Kummer extension (not one simple radical)");
    BR_THROWS("brTower[Tan[x] Sqrt[Log[x]], x]", "Unimplemented", "trigonometric or hyperbolic function (TrigToExp front end over Q(i))");
    BR_EQ("brTower[Exp[x] Sqrt[Exp[x] + 1], x][\"shape\"]", "\"A\"");
    /* brTowerD, brEta, brToSurface */
    BR_EQ("brTowerD[T139, brT2^2]", "(2 brT1 brT2)/(1 + brT1)");
    BR_EQ("brEta[T134]", "1/x");
    BR_EQ("brEta[T158]", "1");
    BR_EQ("brToSurface[T139, brT2 + brT1]", "E^x + Log[1 + E^x]");
}

/* Properties over a seeded random family of surface integrands. */
static void test_tower_properties(void) {
    br_run("SeedRandom[41]; brK = {Log[x], E^x, Log[1 + x], E^(x^2), Log[x + 2], E^(2 x)};"
           " brR[] := RandomInteger[{-3, 3}];"
           " brGen[] := Module[{k1 = RandomChoice[brK], k2 = RandomChoice[brK], p},"
           "   p = brR[] + RandomInteger[{1, 3}] k1 + brR[] x + RandomChoice[{0, 1}] k2;"
           "   (x + brR[] k1)/(1 + x^2 + RandomInteger[{0, 1}] k2^2) * RandomChoice[{1, Sqrt[p], p^(1/3), 1/Sqrt[p]}]]; True");
    /* back substitution reproduces the integrand; kernels' derivations agree */
    BR_TRUE("And @@ Table[Module[{f = brGen[], T}, T = brCatch[brTower[f, x]];"
            "   brFailQ[T] && brFailTag[T][[1]] === \"Decline\" ||"
            "   (PossibleZeroQ[brToSurface[T, T[\"f\"]] - f] &&"
            "    And @@ MapThread[PossibleZeroQ[D[#1, x] - brToSurface[T, #2]] &,"
            "      {T[\"kernels\"], Rest[T[\"rules\"]][[All, 2]]}])], {30}]");
    /* Shape A: theta is last, no monomial depends on theta, radicand involves theta */
    BR_TRUE("And @@ Table[Module[{f = brGen[], T}, T = brCatch[brTower[f, x]];"
            "   brFailQ[T] || T[\"shape\"] =!= \"A\" ||"
            "   (Last[T[\"syms\"]] === T[\"theta\"] && !FreeQ[T[\"radical\"][\"b\"], T[\"theta\"]] &&"
            "    FreeQ[Most[Rest[T[\"rules\"]]][[All, 2]], T[\"theta\"]])], {30}]");
}

static void test_flatten(void) {
    BR_TRUE("brMoebiusRadicandQ[x + 1, x]");
    BR_TRUE("brMoebiusRadicandQ[(2 x + 1)/(x - 3), x]");
    BR_FALSE("brMoebiusRadicandQ[x^2 + 1, x]");
    BR_FALSE("brMoebiusRadicandQ[(x + 1)/(2 x + 2), x]");
    BR_FALSE("brMoebiusRadicandQ[3, x]");
    BR_EQ("brFlattenLinear[Log[x]/Sqrt[x + 1], x, x + 1, 2, u]",
          "<|\"g\" -> 2 Log[-1 + u^2], \"u\" -> u, \"x\" -> -1 + u^2, \"back\" -> u -> Sqrt[1 + x]|>");
    BR_EQ("brFlattenLinear[Sqrt[x^2 + 1], x, x^2 + 1, 2, u]", "None");
    BR_TRUE("Module[{r = brFlattenLinear[E^x ((x - 1)/(x + 2))^(1/3), x, (x - 1)/(x + 2), 3, u]},"
            " r =!= None && FreeQ[r[\"g\"], Power[_, 1/3]]]");
    /* property: g(u)/x'(u) at u = b^(1/n) is the integrand */
    BR_TRUE("SeedRandom[42]; And @@ Table[Module[{a = RandomInteger[{1, 4}], c = RandomInteger[{-3, 3}],"
            "   n = RandomInteger[{2, 4}], f, r},"
            "   f = Log[x + 2] (a x + c)^(1/n) + E^x/(a x + c)^(1/n);"
            "   r = brFlattenLinear[f, x, a x + c, n, u];"
            "   PossibleZeroQ[(r[\"g\"]/D[r[\"x\"], u] /. u -> (a x + c)^(1/n)) - f]], {20}]");
}

static void test_poly(void) {
    BR_EQ("brD[T134, brT1^2 + x]", "1 + (2 brT1)/x");
    BR_EQ("brThetaContent[T134, x^2 brT1^2 + 2 x brT1]", "{x, 2 brT1 + brT1^2 x}");
    BR_EQ("brThetaContent[T134, 0]", "{0, 0}");
    BR_EQ("brThetaContent[T134, (brT1 + 1)/x]", "{1/x, 1 + brT1}");
    BR_EQ("brThetaContent[T134, -3 brT1]", "{-3, brT1}");
    BR_THROWS("brThetaContent[T134, 1/brT1]", "Internal", "brThetaContentNotPolynomial");
    BR_EQ("brNumericalLC[-x^2 + 1]", "-1");
    BR_EQ("brNumericalLC[3]", "1");
    BR_EQ("brSplit[T134, brT1 (brT1 + x)]", "{brT1 (brT1 + x), 1}");
    BR_EQ("brSplit[T158, brT1^2 (brT1 + x)]", "{brT1 + x, brT1^2}");
    BR_EQ("brSplitSquarefree[T158, brT1^3 (brT1 + x)^2 (brT1 - 1)]", "{{-1 + brT1, brT1 + x, 1}, {1, 1, brT1}}");
    BR_TRUE("brNormalQ[T134, brT1 + x]");
    BR_FALSE("brNormalQ[T134, (brT1 + x)^2]");
    BR_TRUE("brSpecialQ[T158, brT1^2]");
    BR_FALSE("brSpecialQ[T158, brT1 + 1]");
    BR_TRUE("brSpecialQ[T158, x + 1]");
    BR_EQ("brSpecialPlaces[T134]", "{\"Infinity\"}");
    BR_EQ("brSpecialPlaces[T158]", "{0, \"Infinity\"}");
    /* properties: c pp == p; pn ps == p with ps special, pn normal (squarefree part) */
    BR_TRUE("SeedRandom[43]; And @@ Table[Module[{p = Expand[(RandomInteger[{1, 3}] x + 1) (brT1^2 + RandomInteger[{-2, 2}] x brT1 + 1)], cp},"
            "   cp = brThetaContent[T134, p]; Together[cp[[1]] cp[[2]] - p] === 0], {20}]");
    BR_TRUE("SeedRandom[44]; And @@ Table[Module[{p = Expand[brT1^RandomInteger[{0, 3}] (brT1 + RandomInteger[{1, 3}] x)^RandomInteger[{1, 2}]], s},"
            "   s = brSplit[T158, p]; Together[s[[1]] s[[2]] - p] === 0 && brSpecialQ[T158, s[[2]]]], {20}]");
}

int main(void) {
    br_init();
    test_helpers();
    test_monomials();
    test_kernels_to_symbols();
    test_detect_radical();
    test_tower();
    test_tower_properties();
    test_flatten();
    test_poly();
    BR_DONE("bronstein_tower_tests");
    return 0;
}
