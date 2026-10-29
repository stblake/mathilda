/* test_trager_jacobian.c -- reduced divisor classes (src/internal/trager/
 * Jacobian.m) and the code built on them: the spec-based torsion search
 * (TorsionModP.m: specFrameIdeal, specIdealX, placeIdealAt, modPFrame,
 * orderModPSpec, torsionOrderSpec) and the generator routes of
 * TragerLogTerms.m (generatorTerms, combineGenerator, rationalElementInverse,
 * logSumResidual).
 *
 * Every expectation is independent of the code under test:
 *   - the frame transform is checked by substituting back (curve equation,
 *     differential), and deg Gt = k m - e prime to m by hand;
 *   - orders at the single place over t = Infinity by hand from
 *     ord(t) = -m, ord(Yt) = -deg Gt;
 *   - class arithmetic against the ideal-based principality test of
 *     PrincipalGen.m (an independent route), on Bronstein's Ex. 6
 *     (Int dx/(x Sqrt[1-x^3]), residue divisor of order exactly 3) and on the
 *     non-torsion point (3, 5) of y^2 = x^3 - 2;
 *   - baby-step giant-step against linear stepping, and canonical keys against
 *     two different addition orders;
 *   - generators by the logarithmic-derivative identity (logSumResidual = 0).
 */
#include "trager_test_utils.h"

/* Context of a remainder c Y dx on Y^m = G in the working frame (regular at
 * infinity), with the residue data of Residues.m. */
static void setup(void) {
    char* r = tr_eval(
        "tjCtx[c_, mm_, GG_] := Module[{B1 = curveBasis[mm, GG, x], cr = Together[c], num, den, rs},"
        "  num = Numerator[cr]; den = Denominator[cr];"
        "  rs = componentResidues[num, den, mm, GG, x];"
        "  withFrameAlpha[Join[rs, <|\"B\" -> B1, \"genus\" -> curveGenus[mm, GG, x], \"N\" -> num,"
        "    \"E\" -> den, \"Ep\" -> Expand[D[den, x]], \"Ysym\" -> Unique[\"Yj\"]|>]]];"
        /* Bronstein Ex. 6 (regular at infinity in x): residues +-1, order 3 */
        "tjEx6 = tjCtx[1/(x (1 - x^3)), 2, 1 - x^3];"
        "tjS6 = <|\"res\" -> tjEx6[\"residues\"], \"b\" -> {1, -1}, \"ref\" -> None|>;"
        /* (3, 5) on y^2 = x^3 - 2 has infinite order (Nagell-Lutz) */
        "tjP35 = tjCtx[5/((x - 3) (x^3 - 2)), 2, x^3 - 2];"
        "tjSP = <|\"res\" -> tjP35[\"residues\"], \"b\" -> {1, -1}, \"ref\" -> None|>;"
        "True");
    free(r);
}

/* ------------------------------------------------------------ frame basics */
static void test_frame_basics(void) {
    /* rootMultiplicity, over Q and GF(5) */
    TR_EQ("rootMultiplicity[fieldQ[], x^2 (x + 1)^3, 0, x]", "2");
    TR_EQ("rootMultiplicity[fieldQ[], x^2 (x + 1)^3, -1, x]", "3");
    TR_EQ("rootMultiplicity[fieldQ[], x^2 (x + 1)^3, 1, x]", "0");
    TR_EQ("rootMultiplicity[fieldGF[5], (x - 1)^2 (x + 3), 1, x]", "2");
    TR_EQ("rootMultiplicity[fieldGF[5], (x - 1)^2 (x + 3), 2, x]", "1");   /* -3 = 2 mod 5 */
    /* frameRoots: rational roots of multiplicity prime to m */
    TR_EQ("frameRoots[fieldQ[], (x - 1) (x^2 + 1), 2, x]", "{1}");
    TR_EQ("frameRoots[fieldQ[], x^2 (x^2 + 3), 2, x]", "{}");            /* e = 2 not prime to 2 */
    TR_EQ("frameRoots[fieldQ[], x^2 (x^2 + 3), 3, x]", "{0}");           /* e = 2 prime to 3 */
    TR_EQ("Sort[frameRoots[fieldGF[5], x^2 + 1, 2, x]]", "{2, 3}");
    TR_EQ("frameRoots[fieldGF[7], x^2 + 1, 2, x]", "{}");                /* -1 not a square mod 7 */
    /* curveBasisOver GF(p) agrees with the reduction of the Q basis (good
     * reduction: 1 - x^3 is squarefree mod 7) */
    TR_TRUE("Module[{b1 = curveBasisOver[fieldGF[7], 2, 1 - x^3, x], b2 = basisMod[curveBasis[2, 1 - x^3, x], 7]},"
            " b1[\"d\"] === b2[\"d\"] && b1[\"mult\"] === b2[\"mult\"] && b1[\"delta\"] === b2[\"delta\"]]");
    TR_TRUE("curveBasisOver[fieldQ[], 3, x^2 + 1, x] === curveBasis[3, x^2 + 1, x]");
}

/* ------------------------------------------------------------- jacFrame */
static void test_jac_frame(void) {
    /* Ex. 6 centred at a = 1 (simple root of 1 - x^3): M = 3, k = 2,
     * deg Gt = k m - e = 4 - 1 = 3, prime to m = 2 */
    char* r = tr_eval("tjFr = jacFrame[fieldQ[], 2, 1 - x^3, 1/(x (1 - x^3)), x, 1]; True");
    free(r);
    TR_EQ("{tjFr[\"k\"], tjFr[\"a\"], tjFr[\"B\"][\"deg\"]}", "{2, 1, 3}");
    TR_EQ("GCD[tjFr[\"B\"][\"deg\"], 2]", "1");
    /* curve: Yt^2 = Gt(t) is the image of Y^2 = 1 - x^3 under x = 1 + 1/t,
     * Yt = Y t^2: Gt = t^4 (1 - (1 + 1/t)^3) */
    TR_TRUE("PossibleZeroQ[tjFr[\"B\"][\"G\"] - Expand[t0^4 (1 - (1 + 1/t0)^3)] /. t0 -> tjFr[\"t\"]]");
    /* differential: c_t Yt dt = c Y dx, i.e. (N/E)(t) Yt = c(x) Y dx/dt with
     * dx/dt = -1/t^2 and Y = Yt/t^2 */
    TR_TRUE("Module[{t = tjFr[\"t\"]}, PossibleZeroQ[tjFr[\"N\"]/tjFr[\"E\"] -"
            " (1/(x (1 - x^3)) /. x -> 1 + 1/t) t^-2 (-1/t^2)]]");
    /* E is monic, Ep its derivative, Bp its branch-free part */
    TR_TRUE("Module[{t = tjFr[\"t\"]}, Coefficient[tjFr[\"E\"], t, Exponent[tjFr[\"E\"], t]] === 1 &&"
            " tjFr[\"Ep\"] === Expand[D[tjFr[\"E\"], t]] &&"
            " PolynomialGCD[tjFr[\"Bp\"], tjFr[\"B\"][\"G\"]] === 1]");
    /* frameToX inverts the substitution */
    TR_TRUE("Module[{t = tjFr[\"t\"]}, Together[frameToX[tjFr, t^3 + Yt0 t, Yt0, Yx] -"
            " (1/(x - 1)^3 + Yx/(x - 1)^2/(x - 1))] === 0]");
    /* over GF(5): the same frame reduced */
    TR_TRUE("Module[{f5 = jacFrame[fieldGF[5], 2, 1 - x^3, 1/(x (1 - x^3)), x, 1]},"
            " f5[\"B\"][\"G\"] === gfMod[tjFr[\"B\"][\"G\"] /. tjFr[\"t\"] -> f5[\"t\"], 5]]");
}

/* ------------------------------------------- orders at P, Popov reduction */
static void test_place_order_reduced_rows(void) {
    char* r = tr_eval("tjBt = tjFr[\"B\"]; tjT = tjFr[\"t\"]; True");
    free(r);
    /* ord(t) = -m = -2, ord(Yt) = -deg Gt = -3, ord(1) = 0, ord(t Yt) = -5;
     * the order of a sum is the least order of its terms */
    TR_EQ("placeOrder[fieldQ[], tjBt, {tjT, 0}, 1]", "-2");
    TR_EQ("placeOrder[fieldQ[], tjBt, {0, 1}, 1]", "-3");
    TR_EQ("placeOrder[fieldQ[], tjBt, {1, 0}, 1]", "0");
    TR_EQ("placeOrder[fieldQ[], tjBt, {0, tjT}, 1]", "-5");
    TR_EQ("placeOrder[fieldQ[], tjBt, {tjT^2, 1}, 1]", "-4");
    TR_EQ("placeOrder[fieldQ[], tjBt, {1, 1}, tjT^2]", "1");             /* (1 + Yt)/t^2 */
    TR_EQ("placeOrder[fieldQ[], tjBt, {0, 0}, 1]", "Infinity");
    TR_EQ("leadingColumn[fieldQ[], tjBt, {tjT^2, 1}]", "0");
    TR_EQ("leadingColumn[fieldQ[], tjBt, {tjT, 1}]", "1");
    /* reducedRows: two rows leading in the same column are combined; the
     * leading columns end distinct and the module is unchanged */
    TR_TRUE("Module[{rows = {{tjT^2, 1}, {tjT, 0}}, rr},"
            " rr = reducedRows[fieldQ[], tjBt, rows];"
            " DuplicateFreeQ[leadingColumn[fieldQ[], tjBt, #] & /@ rr] &&"
            " hnf[fieldQ[], rr, tjT] === hnf[fieldQ[], rows, tjT]]");
    TR_TRUE("Module[{rows = {{1 + tjT^3, tjT}, {tjT^2, 1 + tjT}}, rr},"
            " rr = reducedRows[fieldQ[], tjBt, rows];"
            " DuplicateFreeQ[leadingColumn[fieldQ[], tjBt, #] & /@ rr] &&"
            " hnf[fieldQ[], rr, tjT] === hnf[fieldQ[], rows, tjT]]");
    /* idealNormDegree: unit ideal 0, (t) ideal = m, fractional (1/t) = -m */
    TR_EQ("idealNormDegree[fieldQ[], tjBt, idealUnit[tjBt]]", "0");
    TR_EQ("idealNormDegree[fieldQ[], tjBt, <|\"M\" -> tjT IdentityMatrix[2], \"den\" -> 1|>]", "2");
    TR_EQ("idealNormDegree[fieldQ[], tjBt, <|\"M\" -> IdentityMatrix[2], \"den\" -> tjT|>]", "-2");
    /* elementInverse: h (num/den) = 1 */
    TR_TRUE("Module[{h = {1 + tjT, tjT^2}, inv},"
            " inv = elementInverse[fieldQ[], tjBt, h];"
            " afTimes[fieldQ[], tjBt, h, rc[fieldQ[], #/inv[[2]]] & /@ inv[[1]]] === {1, 0}]");
    TR_TRUE("Module[{h = {2, 3}, inv, F5 = fieldGF[5], B5 = curveBasisOver[fieldGF[5], 2, tjBt[\"G\"], tjT]},"
            " inv = elementInverse[F5, B5, h];"
            " afTimes[F5, B5, h, rc[F5, #/inv[[2]]] & /@ inv[[1]]] === {1, 0}]");
}

/* ---------------------------------------------------- reduction and classes */
static void test_classes_char0(void) {
    char* r = tr_eval(
        "tjJn = specFrameIdeal[fieldQ[], tjFr, tjS6, Identity];"
        "tjC1 = jacClass[fieldQ[], tjBt, tjJn[[1]], tjJn[[2]], True]; True");
    free(r);
    /* genus 1: reduced degree <= 1; delta is not principal, so d = 1 */
    TR_EQ("tjC1[\"d\"]", "1");
    TR_FALSE("jacZeroQ[tjC1]");
    /* the reduction relation: Jnew = h^(-1) J, i.e. Jnew (num) = den J for
     * h = num/den */
    TR_TRUE("Module[{rr = jacReduce[fieldQ[], tjBt, tjJn[[1]], tjJn[[2]]], h, num, den},"
            " h = rr[[3]]; den = Fold[PolynomialLCM, 1, Denominator[Together[#]] & /@ h];"
            " num = Together[# den] & /@ h;"
            " idealEqualQ[fieldQ[], tjBt, idealTimes[fieldQ[], tjBt, rr[[1]], idealFromGenerators[fieldQ[], tjBt, {num}]],"
            "  <|\"M\" -> Map[Expand[den #] &, tjJn[[1]][\"M\"], {2}], \"den\" -> tjJn[[1]][\"den\"]|>]]");
    /* 2 delta not principal, 3 delta principal (order exactly 3) */
    TR_FALSE("jacZeroQ[jacMultiple[fieldQ[], tjBt, tjC1, 2]]");
    TR_TRUE("jacZeroQ[jacMultiple[fieldQ[], tjBt, tjC1, 3]]");
    TR_FALSE("jacZeroQ[jacMultiple[fieldQ[], tjBt, tjC1, 4]]");
    TR_TRUE("jacZeroQ[jacMultiple[fieldQ[], tjBt, tjC1, 6]]");
    /* independent route: the ideal of 3 delta in the working frame is
     * principal, that of 2 delta is not (PrincipalGen.m) */
    TR_TRUE("principalQ[fieldQ[], tjEx6[\"B\"], idealPower[fieldQ[], tjEx6[\"B\"], specIdealX[tjEx6, tjS6], 3]]");
    TR_FALSE("principalQ[fieldQ[], tjEx6[\"B\"], idealPower[fieldQ[], tjEx6[\"B\"], specIdealX[tjEx6, tjS6], 2]]");
    /* addition is commutative and keys are canonical: c + 2c = 2c + c = 3c */
    TR_TRUE("Module[{c2 = jacAdd[fieldQ[], tjBt, tjC1, tjC1]},"
            " jacKey[jacAdd[fieldQ[], tjBt, tjC1, c2]] === jacKey[jacAdd[fieldQ[], tjBt, c2, tjC1]]]");
    TR_TRUE("jacKey[jacMultiple[fieldQ[], tjBt, tjC1, 4]] === jacKey[tjC1]");      /* 4c = c */
    /* jacOrder (linear then BSGS) */
    TR_EQ("jacOrder[fieldQ[], tjBt, tjC1, 50]", "3");
    TR_EQ("jacMinimalOrder[fieldQ[], tjBt, tjC1, 12]", "3");
    TR_TRUE("MatchQ[catchTrager[jacOrder[fieldQ[], tjBt, tjC1, 2]], tragerFail[\"Decline\", _, _]]");
    /* the compact generator H: 3 delta + div(Prod h^e) = 0 */
    TR_TRUE("Module[{c3 = jacMultiple[fieldQ[], tjBt, tjC1, 3, True]}, ListQ[c3[\"H\"]] && Length[c3[\"H\"]] >= 1]");
    TR_TRUE("mergeFactors[{{a1, 1}, {a2, 2}, {a1, -1}, {a3, 1}}] === {{a2, 2}, {a3, 1}}");
}

/* ------------------------------------------- GF(p): order search, two gens */
static void test_classes_modp(void) {
    char* r = tr_eval(
        "tjGp = goodPrime[tjEx6, 7]; {tjFr7, tjRed7} = modPFrame[tjEx6, tjGp];"
        "tjF7 = fieldGF[7]; tjJ7 = specFrameIdeal[tjF7, tjFr7, tjS6, tjRed7];"
        "tjC7 = jacClass[tjF7, tjFr7[\"B\"], tjJ7[[1]], tjJ7[[2]]]; True");
    free(r);
    /* the order mod 7 divides... equals the characteristic-0 order 3 here
     * (3 is prime to 7, so it survives reduction) */
    TR_EQ("jacOrder[tjF7, tjFr7[\"B\"], tjC7, 20]", "3");
    TR_EQ("orderModPSpec[tjEx6, tjS6, tjGp]", "3");
    /* linear stepping and BSGS agree on the non-torsion point mod p */
    TR_TRUE("Module[{gp = goodPrime[tjP35, 11], fr, red, F, J, c, lin, M = 1},"
            " {fr, red} = modPFrame[tjP35, gp]; F = fieldGF[11];"
            " J = specFrameIdeal[F, fr, tjSP, red]; c = jacClass[F, fr[\"B\"], J[[1]], J[[2]]];"
            " lin = c; While[!jacZeroQ[lin], lin = jacAdd[F, fr[\"B\"], lin, c]; M++];"
            " M === jacOrderGens[F, fr[\"B\"], jacWithGenerators[F, fr[\"B\"], c], 40]]");
    /* the two-generator form multiplies like the full ideal */
    TR_TRUE("Module[{c2 = jacAdd[tjF7, tjFr7[\"B\"], tjC7, tjC7], cg = jacWithGenerators[tjF7, tjFr7[\"B\"], tjC7]},"
            " KeyExistsQ[cg, \"gens\"] && jacKey[jacAdd[tjF7, tjFr7[\"B\"], c2, cg]] === jacKey[jacAdd[tjF7, tjFr7[\"B\"], c2, tjC7]]]");
    TR_TRUE("Module[{g = twoGenerators[tjF7, tjFr7[\"B\"], tjC7[\"J\"]]},"
            " g =!= None && kc[tjF7, g[[2]] - kMonic[tjF7, Times @@ Diagonal[tjC7[\"J\"][\"M\"]], tjFr7[\"t\"]]] === 0]");
    /* torsionOrderSpec: Ex. 6 has order 3; (3, 5) is a certificate */
    TR_EQ("torsionOrderSpec[tjEx6, tjS6]", "3");
    TR_TRUE("MatchQ[catchTrager[torsionOrderSpec[tjP35, tjSP]], tragerFail[\"NonElementary\", _String, _]]");
    /* the legacy entry points delegate to the spec versions */
    TR_EQ("torsionOrder[tjEx6, {1, -1}]", "3");
    TR_EQ("orderModP[tjEx6, {1, -1}, goodPrime[tjEx6, 13]]", "3");
}

/* ---------------------------------------------------- the place over alpha */
static void test_place_ideal(void) {
    /* 1 - x^3 has the simple root 1: P^2 = (x - 1) O (totally ramified, m = 2) */
    TR_TRUE("Module[{P = placeIdealAt[fieldQ[], tjEx6[\"B\"], 1]},"
            " idealNormDegree[fieldQ[], tjEx6[\"B\"], P] === 1 &&"
            " idealEqualQ[fieldQ[], tjEx6[\"B\"], idealPower[fieldQ[], tjEx6[\"B\"], P, 2],"
            "  idealFromGenerators[fieldQ[], tjEx6[\"B\"], {{x - 1, 0}}]]]");
    /* a double root of G on a cube-root curve: Y^3 = x^2 (x + 1), place over 0
     * (e = 2, s = 2, r = 1: uniformiser Y^2/x); P^3 = (x) */
    TR_TRUE("Module[{B3 = curveBasis[3, x^2 (x + 1), x], P},"
            " P = placeIdealAt[fieldQ[], B3, 0];"
            " idealNormDegree[fieldQ[], B3, P] === 1 &&"
            " idealEqualQ[fieldQ[], B3, idealPower[fieldQ[], B3, P, 3], idealFromGenerators[fieldQ[], B3, {{x, 0, 0}}]]]");
    TR_TRUE("MatchQ[catchTrager[placeIdealAt[fieldQ[], curveBasis[2, x^2 + 1, x], 0]], tragerFail[\"Internal\", _, _]]");
    /* specIdealX with a Pinf term is Delta - kappa P in the working frame:
     * genus 0 (Y^3 = (3x + 2)^2), so it is principal */
    TR_TRUE("Module[{c0, K, g, I0, k0, sp, B3 = curveBasis[3, Expand[(3 x + 2)^2], x], cr = Together[(-2)/(-x^2 + 2 x - 3)]},"
            " c0 = withFrameAlpha[Join[componentResiduePolynomial[Numerator[cr], Denominator[cr], 3, B3[\"G\"], x],"
            "  <|\"B\" -> B3, \"genus\" -> 0, \"N\" -> Numerator[cr], \"E\" -> Denominator[cr],"
            "    \"Ep\" -> Expand[D[Denominator[cr], x]], \"Ysym\" -> Yj0|>]];"
            " K = fieldNF[Zq, c0[\"R\"] /. c0[\"Z\"] -> Zq, None];"
            " c0 = Join[c0, <|\"F\" -> K|>];"
            " I0 = residueIdeal[K, c0[\"B\"], Zq, c0[\"N\"], c0[\"Ep\"], c0[\"Bp\"]];"
            " k0 = idealDegree[K, c0[\"B\"], I0];"
            " sp = <|\"res\" -> {Zq}, \"b\" -> {1}, \"ref\" -> {\"Pinf\", -k0}|>;"
            " c0[\"alpha\"] =!= None && principalQ[K, c0[\"B\"], specIdealX[c0, sp]]]");
}

/* ------------------------------------------------- generators and residuals */
static void test_generators(void) {
    /* Ex. 6 through both generator routes: direct (N <= DirectPowerMaxOrder)
     * and reduced classes (forced by lowering the threshold); each a valid
     * logarithmic part: (1/3) Sum e Log[u] has the residues of the remainder */
    char* r = tr_eval("tjRem = {0, 1/(x (1 - x^3))}; True");
    free(r);
    TR_TRUE("Module[{gt = generatorTerms[tjEx6, tjS6, 3]}, ListQ[gt] &&"
            " afZeroQ[fieldQ[], logSumResidual[fieldQ[], tjEx6[\"B\"], tjRem, {1/3 #[[1]], #[[2]]} & /@ gt]]]");
    TR_TRUE("Module[{old = tragerOption[\"DirectPowerMaxOrder\"], gt},"
            " tragerOption[\"DirectPowerMaxOrder\"] = 1;"
            " gt = generatorTerms[tjEx6, tjS6, 3];"
            " tragerOption[\"DirectPowerMaxOrder\"] = old;"
            " ListQ[gt] && afZeroQ[fieldQ[], logSumResidual[fieldQ[], tjEx6[\"B\"], tjRem, {1/3 #[[1]], #[[2]]} & /@ gt]]]");
    /* not principal at N = 2: $Failed on both routes */
    TR_EQ("generatorTerms[tjEx6, tjS6, 2]", "$Failed");
    TR_TRUE("Module[{old = tragerOption[\"DirectPowerMaxOrder\"], gt},"
            " tragerOption[\"DirectPowerMaxOrder\"] = 1; gt = generatorTerms[tjEx6, tjS6, 2];"
            " tragerOption[\"DirectPowerMaxOrder\"] = old; gt === $Failed]");
    /* logSumResidual detects a wrong coefficient */
    TR_FALSE("Module[{gt = generatorTerms[tjEx6, tjS6, 3]},"
             " afZeroQ[fieldQ[], logSumResidual[fieldQ[], tjEx6[\"B\"], tjRem, {1/2 #[[1]], #[[2]]} & /@ gt]]]");
    /* logSumResidual and logPartResidual agree on zero-ness */
    TR_TRUE("Module[{lg = tragerLogTerms[tjEx6]},"
            " afZeroQ[fieldQ[], logSumResidual[fieldQ[], tjEx6[\"B\"], tjRem, lg]] &&"
            " afZeroQ[fieldQ[], logPartResidual[fieldQ[], tjEx6[\"B\"], tjRem, lg]]]");
    /* rationalElementInverse: u (1/u) = 1 for rational coordinates, over Q and Q(i) */
    TR_TRUE("Module[{u = {(1 + x)/x, 1/(x - 1)}},"
            " afTimes[fieldQ[], tjEx6[\"B\"], u, rationalElementInverse[fieldQ[], tjEx6[\"B\"], u]] === {1, 0}]");
    TR_TRUE("Module[{K = fieldNF[ti, ti^2 + 1, I], B2 = curveBasis[2, x^2 + 1, x], u},"
            " u = {(x + ti)/(x - 2), 1};"
            " afTimes[K, B2, u, rationalElementInverse[K, B2, u]] === {1, 0}]");
    /* combineGenerator: a small product is multiplied out; a large one kept */
    TR_TRUE("Module[{u = {1 + x, 1}, v = {x, 0}, c},"
            " c = combineGenerator[fieldQ[], tjEx6[\"B\"], {{2, u}, {-1, v}}];"
            " Length[c] === 1 && c[[1, 1]] === 1 &&"
            " afTimes[fieldQ[], tjEx6[\"B\"], c[[1, 2]], v] === afTimes[fieldQ[], tjEx6[\"B\"], u, u]]");
    TR_EQ("Length[combineGenerator[fieldQ[], tjEx6[\"B\"], {{5, {1 + x, 1}}, {-4, {x, 0}}}]]", "2");
}

/* ------------------------------------------------- Hensel-lifted residue powers */
static void test_hensel_powers(void) {
    /* k Delta and -k Delta from the lifted two generators equal the k-th power
     * of the residue ideal and of its inverse (idealPower: an independent route) */
    TR_TRUE("Module[{B = tjEx6[\"B\"], N0 = tjEx6[\"N\"], Ep = tjEx6[\"Ep\"], Bp = tjEx6[\"Bp\"], ok = True},"
            " Do[ok = ok && idealEqualQ[fieldQ[], B, divisorIdeal[fieldQ[], B, {1}, {k}, N0, Ep, Bp],"
            "   idealPower[fieldQ[], B, residueIdeal[fieldQ[], B, 1, N0, Ep, Bp], k]], {k, {1, 2, 3, 5, -1, -2, -4}}]; ok]");
    /* over Q(Sqrt[6]) on y^2 = x^2 + 1 (residues +-Sqrt[6]/4) */
    TR_TRUE("Module[{c6 = tjCtx[1/(x^2 - 2), 2, x^2 + 1], ok = True},"
            " Do[ok = ok && idealEqualQ[c6[\"F\"], c6[\"B\"],"
            "   divisorIdeal[c6[\"F\"], c6[\"B\"], {c6[\"residues\"][[1]]}, {k}, c6[\"N\"], c6[\"Ep\"], c6[\"Bp\"]],"
            "   idealPower[c6[\"F\"], c6[\"B\"], residueIdeal[c6[\"F\"], c6[\"B\"], c6[\"residues\"][[1]], c6[\"N\"], c6[\"Ep\"], c6[\"Bp\"]], k]],"
            "  {k, {2, 3, -3}}]; ok]");
    /* henselRoot: phi_k^m = G modulo A^k, phi_k = phi modulo A */
    TR_TRUE("Module[{ad = residueData[fieldQ[], tjEx6[\"B\"], 1, tjEx6[\"N\"], tjEx6[\"Ep\"], tjEx6[\"Bp\"]], ph},"
            " ph = henselRoot[fieldQ[], 1 - x^3, 2, ad[[2]], ad[[1]], 5, x];"
            " PolynomialRemainder[Expand[ph^2 - (1 - x^3)], Expand[ad[[1]]^5], x] === 0 &&"
            " PolynomialRemainder[ph - ad[[2]], ad[[1]], x] === 0]");
    /* scaleSpec */
    TR_TRUE("scaleSpec[<|\"res\" -> {a}, \"b\" -> {1, -2}, \"ref\" -> {\"Pinf\", -3}|>, 4] ==="
            " <|\"res\" -> {a}, \"b\" -> {4, -8}, \"ref\" -> {\"Pinf\", -12}|>");
    TR_TRUE("scaleSpec[<|\"res\" -> {a}, \"b\" -> {1}, \"ref\" -> {\"R0\", 2, -1}|>, 3][\"ref\"] === {\"R0\", 2, -3}");
}

int main(void) {
    tr_init();
    setup();
    test_frame_basics();
    test_jac_frame();
    test_place_order_reduced_rows();
    test_classes_char0();
    test_classes_modp();
    test_place_ideal();
    test_generators();
    test_hensel_powers();
    printf("trager_jacobian: %d assertions passed\n", tr_assertions);
    return 0;
}
