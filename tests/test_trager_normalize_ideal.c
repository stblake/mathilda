/* test_trager_normalize_ideal.c
 *
 * Unit tests for src/internal/trager/Normalize.m, Ideal.m and PrincipalGen.m
 * -- Hermite normal form over K[x] ([T] Ch 2 p.25), normalisation of a module
 * basis at infinity ([T] Ch 2 Sec 3 pp.27-29), fractional ideals of the
 * integral closure as divisors ([T] Ch 5 Sec 3 pp.59-62, Ch 2 Sec 2) and the
 * principality test ([T] Ch 6 Sec 1 pp.64-65), over the three coefficient
 * fields Q, Q(Sqrt[2]) = fieldNF[t, t^2 - 2, Sqrt[2]] and GF(p).
 *
 * Coverage, function by function:
 *   hnf                  upper triangular, monic diagonal, entries above a
 *                        pivot of lower degree, SAME module: every input row
 *                        lies in the span of the HNF rows (exact
 *                        back-substitution) and det(HNF) equals the monic gcd
 *                        of the maximal minors of the input (so the modules
 *                        coincide); uniqueness: a known HNF scrambled by a
 *                        unimodular matrix comes back exactly (Q, Q(t), GF(7));
 *                        zero rows, redundant rows, m = 1..5; rank-deficient
 *                        throws
 *   rowInfData           hand-computed x-adic orders and leading vectors
 *   normalizeAtInfinity  a hand example, already-normal input unchanged; the
 *                        result spans the same module (containment + det up to
 *                        a constant), leading vectors independent (kRank = m),
 *                        k_i recomputed from rowInfData, and the defining
 *                        property: ord(Sum p_i row_i) = min(k_i - deg p_i)
 *   idealFromGenerators  principal ideals: det of the HNF = monic norm
 *                        (Resultant_Y, independent of afNorm)
 *   residueIdeal         hand HNFs for Bronstein Ex. 6 and over Q(Sqrt[2])
 *   idealTimes           P + Q = div(x) on Ex. 6; commutative, associative,
 *                        unit
 *   idealPower           k = 0, 1, 2, 3, -1, -2; 3P = (Y - 1) on Ex. 6
 *   idealUnit, idealInverse   I I^-1 = O over Q, Q(t), GF(p); inverse of (x)
 *   polyAdjugate         A adj(A) = det(A) I over Q, Q(t), GF(p); 1x1
 *   idealReduce          common factors removed, value preserved
 *   idealEqualQ          equal under different denominators; unequal
 *   principalGenerator / principalQ
 *                        genus 0: always principal, div(u) = delta checked as
 *                        (u) O = I(delta) and u regular at infinity, plus the
 *                        hand generator (1 + Y)/z of the shifted 1/Sqrt[x^2+1];
 *                        Bronstein Ex. 6 (genus 1): delta, 2 delta not
 *                        principal, 3 delta = div((Y - 1)/(Y + 1)); Ex. 4
 *                        (genus 3): order 4, generator x^4/(1 + Y); the
 *                        non-torsion divisor (3, 5) - (3, -5) on Y^2 = x^3 - 2
 *                        (rank-1 Mordell curve); over GF(p), principal at the
 *                        group order #E(F_p) counted by brute force
 */
#include "trager_test_utils.h"

/* Evaluate a statement for its side effects. */
static void tr_do(const char* in) { free(tr_eval(in)); }

static void define_helpers(void) {
    tr_do(
        "tQ = fieldQ[]; tN = fieldNF[t, t^2 - 2, Sqrt[2]]; tG7 = fieldGF[7]; tG5 = fieldGF[5];"
        /* HNF shape: upper triangular, monic diagonal, reduced above pivots. */
        "tShapeQ[tF_, tH_] := Module[{tn = Length[tH]},"
        "  AllTrue[Flatten[Table[tH[[ti, tj]], {ti, tn}, {tj, ti - 1}]], # === 0 &] &&"
        "  AllTrue[Range[tn], kLC[tF, tH[[#, #]], x] === 1 &] &&"
        "  AllTrue[Flatten[Table[{ti, tj}, {tj, tn}, {ti, tj - 1}], 1],"
        "    kDeg[tF, tH[[#[[1]], #[[2]]]], x] < kDeg[tF, tH[[#[[2]], #[[2]]]], x] &]];"
        /* row in the F[x]-span of the rows of the upper-triangular tH:
         * exact back-substitution with polynomial quotients. */
        "tInSpanQ[tF_, tH_, trow_] := Module[{tr = kc[tF, #] & /@ trow, tok = True, tqr},"
        "  Do[If[tok, tqr = kDivRem[tF, tr[[tj]], tH[[tj, tj]], x];"
        "      If[tqr[[2]] =!= 0, tok = False,"
        "        tr = MapThread[kc[tF, #1 - tqr[[1]] #2] &, {tr, tH[[tj]]}]]], {tj, Length[tH]}];"
        "  tok && AllTrue[tr, # === 0 &]];"
        /* gcd of the maximal minors, monic */
        "tMinorGcd[tF_, tA_] := Module[{tms},"
        "  tms = Select[kc[tF, Det[tA[[#]]]] & /@ Subsets[Range[Length[tA]], {Length[First[tA]]}], # =!= 0 &];"
        "  kMonic[tF, Fold[kGCD[tF, #1, #2, x] &, tms], x]];"
        "tHnfOK[tF_, tA_] := Module[{tH = hnf[tF, tA, x]},"
        "  tShapeQ[tF, tH] && AllTrue[tA, tInSpanQ[tF, tH, #] &] &&"
        "  kc[tF, Det[tH]] === tMinorGcd[tF, tA]];"
        /* normalizeAtInfinity contract */
        "tNormOK[tF_, tB_, tM_, tden_] := Module[{tres = normalizeAtInfinity[tF, tB, tM, tden, x], tR, tH, trat, tdat},"
        "  tR = tres[\"rows\"]; tH = hnf[tF, tM, x];"
        "  trat = rc[tF, Det[tR]/Det[tM]];"
        "  tdat = rowInfData[tF, tB, #, tden] & /@ tR;"
        "  AllTrue[tR, tInSpanQ[tF, tH, #] &] && FreeQ[trat, x] && trat =!= 0 &&"
        "  kRank[tF, Last /@ tdat] === Length[tM] && (First /@ tdat) === tres[\"k\"]];"
        /* the normal-at-infinity property on a polynomial combination tp */
        "tNormProp[tF_, tB_, tR_, tks_, tden_, tp_] := Module[{tsum},"
        "  tsum = kc[tF, #] & /@ Total[tp tR];"
        "  rowInfData[tF, tB, tsum, tden][[1]] === Min[tks - (kDeg[tF, #, x] & /@ tp)]];"
        /* norm of an element by a resultant (independent of afNorm) */
        "tNormY[tB_, tu_] := Module[{tY},"
        "  Expand[Together[Resultant[tY^tB[\"m\"] - tB[\"G\"], Sum[tu[[ti + 1]] tY^ti/tB[\"d\"][[ti + 1]], {ti, 0, tB[\"m\"] - 1}], tY]]]];"
        /* (u) O as a fractional ideal */
        "tPrincipalIdeal[tF_, tB_, tu_] := Module[{tden = PolynomialLCM @@ Join[{1}, Denominator[Together[#]] & /@ tu], trow},"
        "  trow = kc[tF, Together[# tden]] & /@ tu;"
        "  <|\"M\" -> idealFromGenerators[tF, tB, {trow}][\"M\"], \"den\" -> tden|>];"
        /* div(u) = delta: (u) O = I(delta) and u has no pole at infinity
         * (deg delta = 0 then forces no zero there either). */
        "tGenQ[tF_, tB_, tI_, tu_] := ListQ[tu] && !afZeroQ[tF, tu] &&"
        "  idealEqualQ[tF, tB, tPrincipalIdeal[tF, tB, tu], tI] && integralAtInfinityQ[tF, tB, tu];"
        "tUnitQ[tF_, tB_, tI_] := idealEqualQ[tF, tB, tI, idealUnit[tB]];"
        "tInvOK[tF_, tB_, tI_] := tUnitQ[tF, tB, idealTimes[tF, tB, tI, idealInverse[tF, tB, tI]]];"
        "tAdjOK[tF_, tA_] := Module[{tn = Length[tA]},"
        "  Map[kc[tF, #] &, tA . polyAdjugate[tF, tA, x], {2}] === kc[tF, Det[tA]] IdentityMatrix[tn]];"
        /* proportional coordinate lists over F */
        "tPropQ[tF_, tu_, tv_] := Module[{ti = SelectFirst[Range[Length[tv]], tv[[#]] =!= 0 &], tl},"
        "  tl = rc[tF, tu[[ti]]/tv[[ti]]]; FreeQ[tl, x] && tl =!= 0 && afZeroQ[tF, tu - tl tv]];"
        /* #E(F_p) for Y^2 = g(x) by brute force */
        "tCount[tg_, tp_] := 1 + Sum[1 + JacobiSymbol[Mod[tg /. x -> ta, tp], tp], {ta, 0, tp - 1}];"
        "True");
}

/* ------------------------------------------------------------- Normalize.m */

static void test_hnf_Q(void) {
    /* A reduced HNF scrambled by a unimodular U comes back exactly. */
    tr_do("tH0 = {{x^2 + 1, 2, 5/2}, {0, x - 2, 1/2}, {0, 0, x^3}};"
          "tU = {{1, x, 0}, {0, 1, 0}, {0, 0, 1}} . {{1, 0, 0}, {x^2 - 1, 1, 0}, {0, 0, 1}} ."
          "     {{1, 0, 0}, {0, 1, 0}, {3, x + 2, 1}} . {{0, 1, 0}, {1, 0, 0}, {0, 0, 1}};"
          "tA = Expand[tU . tH0]; True");
    TR_TRUE("tShapeQ[tQ, tH0]");
    TR_TRUE("hnf[tQ, tA, x] === tH0");
    TR_TRUE("hnf[tQ, Join[tA, {Expand[x tA[[1]] + tA[[2]]]}], x] === tH0");     /* redundant row */
    TR_TRUE("hnf[tQ, Join[{{0, 0, 0}}, tA], x] === tH0");                       /* zero row */
    TR_TRUE("hnf[tQ, 3 tA, x] === tH0");                                       /* unit scaling */
    TR_TRUE("tHnfOK[tQ, tA]");
    /* An unreduced upper-triangular input is reduced above the pivots:
     * x - (x - 2) = 2, 3 - 1/2 = 5/2. */
    TR_TRUE("hnf[tQ, {{x^2 + 1, x, 3}, {0, x - 2, 1/2}, {0, 0, x^3}}, x] === tH0");
    TR_EQ("hnf[tQ, {{2 x + 4}}, x]", "{{2 + x}}");
    TR_EQ("hnf[tQ, {{2 x + 4}, {x^2 - 4}}, x]", "{{2 + x}}");
    TR_EQ("hnf[tQ, {{2, 3}, {1, 1}}, x]", "{{1, 0}, {0, 1}}");
    TR_EQ("hnf[tQ, {{0, 0}, {x, 1}, {0, x}}, x]", "{{x, 1}, {0, x}}");
    TR_EQ("hnf[tQ, {{0, x}, {x, 1}}, x]", "{{x, 1}, {0, x}}");
    /* non-square inputs: same module via the minors */
    TR_TRUE("tHnfOK[tQ, {{x^2, x}, {x, 1 + x}, {x^3, 0}, {1, x^2}}]");
    TR_TRUE("tHnfOK[tQ, {{x^2 - 1, x + 1}, {x^2 + x, x}, {x^3, x^2 - 1}}]");
    /* m = 4, 5 */
    TR_TRUE("tHnfOK[tQ, {{x, 1, 0, 2}, {x^2, x + 1, 1, 0}, {1, 0, x, x^2}, {0, x^3, 1, x - 1}}]");
    TR_TRUE("tHnfOK[tQ, {{x, 1, 0, 2, 0}, {x^2, x + 1, 1, 0, 1}, {1, 0, x, x^2, 0},"
            " {0, x^3, 1, x - 1, x}, {x^2 + 1, 0, 0, 0, x^2}, {x, x, x, x, x}}]");
    TR_TRUE("tHnfOK[tQ, {{x^2 (x - 1), 0}, {0, x^2 (x - 1)}, {x (x - 1), x (x - 1)}}]");
    /* rank deficiency */
    TR_TRUE("catchTrager[hnf[tQ, {{x, 1}, {x^2, x}}, x]][[1 ;; 2]] === tragerFail[\"Internal\", \"hnfRankDeficient\"]");
    TR_TRUE("catchTrager[hnf[tQ, {{x, 1}}, x]][[1 ;; 2]] === tragerFail[\"Internal\", \"hnfRankDeficient\"]");
    TR_TRUE("catchTrager[hnf[tQ, {{0, 0}, {0, 0}}, x]][[1 ;; 2]] === tragerFail[\"Internal\", \"hnfRankDeficient\"]");
}

static void test_hnf_NF_GF(void) {
    /* Q(Sqrt[2]): reduced HNF with t-coefficients, scrambled by U(t x). */
    tr_do("tH0n = {{x^2 + t, 4, 3}, {0, x - 2 t, 1/2}, {0, 0, x^3 + t}};"
          "tAn = Map[kc[tN, #] &, Expand[(tU /. x -> t x) . tH0n], {2}]; True");
    TR_TRUE("tShapeQ[tN, tH0n]");
    TR_TRUE("hnf[tN, tAn, x] === tH0n");
    TR_TRUE("tHnfOK[tN, tAn]");
    /* x^2 - 2 = (x - t)(x + t): the module generated by x - t and x^2 - 2 */
    TR_TRUE("hnf[tN, {{x^2 - 2}, {x + t}}, x] === {{x + t}}");
    TR_TRUE("tHnfOK[tN, {{x - t, 1}, {x^2 - 2, x + t}, {t, x}}]");
    TR_TRUE("tHnfOK[tN, {{t x, 1, 0}, {0, x - t, t}, {x^2, 0, 1 + t x}}]");
    /* a t-dependent leading coefficient is made monic: (t x + 1) -> x + t/2 */
    TR_TRUE("hnf[tN, {{t x + 1}}, x] === {{x + t/2}}");
    /* GF(7) */
    tr_do("tH0g = {{x^2 + 1, 2, 6}, {0, x + 5, 4}, {0, 0, x^3}};"
          "tAg = gfMod[#, 7] & /@ Expand[tU . tH0g]; True");
    TR_TRUE("tShapeQ[tG7, tH0g]");
    TR_TRUE("hnf[tG7, tAg, x] === tH0g");
    TR_TRUE("hnf[tG7, {{x^2 + 1, x, 3}, {0, x + 5, 4}, {0, 0, x^3}}, x] === tH0g");
    TR_TRUE("tHnfOK[tG7, tAg]");
    TR_TRUE("tHnfOK[tG5, {{x^2, x}, {x, 1 + x}, {x^3, 0}, {1, x^2}}]");
    TR_TRUE("tHnfOK[tG5, {{3 x + 1, 2, 0}, {x^2, 4 x, 1}, {0, 1, x^4 + 2}}]");
    /* over GF(2), x^2 + 1 = (x + 1)^2 */
    TR_TRUE("hnf[fieldGF[2], {{x^2 + 1}, {x^3 + x^2 + x + 1}}, x] === {{1 + x^2}}");
    TR_TRUE("hnf[fieldGF[2], {{x + 1, 1}, {x^2 + 1, 0}}, x] === {{1 + x, 1}, {0, 1 + x}}");
}

static void test_rowInfData(void) {
    /* Y^2 = x^3 + 1: delta = {0, 2}, v_1 = Y/x^2. */
    tr_do("tBe = curveBasis[2, x^3 + 1, x]; True");
    TR_EQ("tBe[\"delta\"]", "{0, 2}");
    TR_EQ("rowInfData[tQ, tBe, {x^2, x}, 1]", "{-3, {0, 1}}");
    TR_EQ("rowInfData[tQ, tBe, {x^2, x}, x]", "{-2, {0, 1}}");
    TR_EQ("rowInfData[tQ, tBe, {x^2, 1}, 1]", "{-2, {1, 1}}");
    TR_EQ("rowInfData[tQ, tBe, {x^3 + 1, 2 x}, 3 x^2]", "{-1, {1/3, 2/3}}");
    TR_EQ("rowInfData[tQ, tBe, {1, 0}, x^2]", "{2, {1, 0}}");
    TR_EQ("rowInfData[tQ, tBe, {0, 0}, 1]", "{Infinity, {0, 0}}");
    /* Y^3 = x^2 + 1: delta = {0, 1, 2}. */
    tr_do("tB32 = curveBasis[3, x^2 + 1, x]; True");
    TR_EQ("tB32[\"delta\"]", "{0, 1, 2}");
    TR_EQ("rowInfData[tQ, tB32, {x, 1, x^2}, x]", "{-3, {0, 0, 1}}");
    TR_EQ("rowInfData[tQ, tB32, {2 x, 3, 0}, x]", "{0, {2, 3, 0}}");
    /* Y^3 = x^2 (x + 1): d = {1, 1, x}, delta = {0, 1, 2 - 1} = {0, 1, 1}. */
    TR_EQ("curveBasis[3, x^2 (x + 1), x][\"delta\"]", "{0, 1, 1}");
    /* over Q(Sqrt[2]): leading vector entries divided by lc(den) = 2 */
    TR_TRUE("rowInfData[tN, tBe, {t x^2, x}, 2 x] === {-2, {0, 1/2}}");
    TR_TRUE("rowInfData[tN, tBe, {t x^2, 1}, 2 x] === {-1, {t/2, 1/2}}");
    TR_TRUE("rowInfData[tN, tBe, {t x^3, 1}, 2 x] === {-2, {t/2, 0}}");
    /* over GF(7): 1/3 = 5 */
    TR_TRUE("rowInfData[tG7, basisMod[tBe, 7], {x^3 + 1, 2 x}, 3 x^2] === {-1, {5, 3}}");
}

static void test_normalizeAtInfinity(void) {
    /* Hand example on Y^2 = x^3+1.  Rows {x^2, 1} and {x^3, x + 1}: k = -2, -3,
     * both leading vectors {1, 1}; c = {1, -1}, i0 = 2 (least k), new row
     * x row1 - row2 = {0, -1} (up to sign) with k = -2. */
    tr_do("tNr = normalizeAtInfinity[tQ, tBe, {{x^2, 1}, {x^3, x + 1}}, 1, x]; True");
    TR_TRUE("MatchQ[tNr[\"rows\"], {{x^2, 1}, {0, 1} | {0, -1}}]");   /* c is unique up to a scalar */
    TR_EQ("tNr[\"k\"]", "{-2, -2}");
    TR_TRUE("tNormOK[tQ, tBe, {{x^2, 1}, {x^3, x + 1}}, 1]");
    /* the original rows violate the normal property, the result satisfies it */
    TR_FALSE("tNormProp[tQ, tBe, {{x^2, 1}, {x^3, x + 1}}, {-2, -3}, 1, {x, -1}]");
    TR_TRUE("tNormProp[tQ, tBe, tNr[\"rows\"], tNr[\"k\"], 1, {x, 1}]");
    TR_TRUE("tNormProp[tQ, tBe, tNr[\"rows\"], tNr[\"k\"], 1, {x^2 + 1, x^3}]");
    TR_TRUE("tNormProp[tQ, tBe, tNr[\"rows\"], tNr[\"k\"], 1, {1, x^2}]");
    /* already normal: unchanged */
    TR_TRUE("normalizeAtInfinity[tQ, tBe, {{1, 0}, {x, 1}}, 1, x] === <|\"rows\" -> {{1, 0}, {x, 1}}, \"k\" -> {0, -2}|>");
    /* invariants, m = 2..4, Q / Q(t) / GF(p) */
    TR_TRUE("tNormOK[tQ, tBe, {{x, 1}, {x^2 + 1, x^3}}, x^2]");
    TR_TRUE("tNormOK[tQ, tB32, {{x, 1, 0}, {x^2, x + 1, 1}, {1, 0, x}}, x^2]");
    TR_TRUE("tNormOK[tQ, tB32, {{x^3, x^2 + 1, x}, {x^2, x, 1}, {0, 0, x^5 + 1}}, 1]");
    TR_TRUE("tNormOK[tQ, curveBasis[4, x^3 + 1, x], {{x, 1, 0, 0}, {x^2, x + 1, 1, 0}, {1, 0, x, x^2}, {0, x^3, 1, x - 1}}, x]");
    TR_TRUE("tNormOK[tQ, curveBasis[5, x^3 + 2, x], IdentityMatrix[5] + x^2 RotateLeft[IdentityMatrix[5]], 1]");
    TR_TRUE("tNormOK[tN, tBe, {{x^2, t}, {x^3 + t, x + 1}}, 1]");
    TR_TRUE("tNormOK[tN, tB32, {{t x, 1, 0}, {x^2, x + t, 1}, {1, 0, t x}}, x]");
    TR_TRUE("tNormOK[tG7, basisMod[tBe, 7], {{x^2, 1}, {x^3, x + 1}}, 1]");
    TR_TRUE("tNormOK[tG5, basisMod[tB32, 5], {{x^3, x^2 + 1, x}, {x^2, x, 1}, {0, 0, x^5 + 1}}, 1]");
    /* the normal property after a multi-step normalisation */
    TR_TRUE("Module[{tr = normalizeAtInfinity[tQ, tB32, {{x^3, x^2 + 1, x}, {x^2, x, 1}, {0, 0, x^5 + 1}}, 1, x]},"
            " AllTrue[{{1, 1, 1}, {x, 1, 0}, {x^2, x - 1, 3}, {0, x, 1}},"
            "   tNormProp[tQ, tB32, tr[\"rows\"], tr[\"k\"], 1, #] &]]");
    /* a zero row (rank-deficient module) is an internal error */
    TR_TRUE("catchTrager[normalizeAtInfinity[tQ, tBe, {{x^2, 1}, {x^3, x}}, 1, x]][[1 ;; 2]] ==="
            " tragerFail[\"Internal\", \"normalizeAtInfinityZeroRow\"]");
}

/* ----------------------------------------------------------------- Ideal.m */

static void test_ideal_basic(void) {
    tr_do("tB6 = curveBasis[2, 1 - x^3, x]; tB3 = curveBasis[3, x^2 (x + 1), x]; True");
    TR_TRUE("idealUnit[tB6] === <|\"M\" -> {{1, 0}, {0, 1}}, \"den\" -> 1|>");
    TR_TRUE("idealUnit[tB3][\"M\"] === IdentityMatrix[3]");
    /* principal ideals from integral generators */
    TR_TRUE("idealFromGenerators[tQ, tB6, {{x, 0}}][\"M\"] === {{x, 0}, {0, x}}");
    TR_TRUE("idealFromGenerators[tQ, curveBasis[2, x, x], {{0, 1}}][\"M\"] === {{x, 0}, {0, 1}}");  /* (Y), Y^2 = x */
    TR_TRUE("idealFromGenerators[tQ, tB3, {{x - 1, 0, 0}}][\"M\"] === (x - 1) IdentityMatrix[3]");
    TR_TRUE("idealFromGenerators[tQ, tB6, {{1, 0}}][\"M\"] === IdentityMatrix[2]");
    /* det(HNF of (u)) = monic norm of u, by a resultant */
    const char* gens[][2] = {
        {"tB3", "{0, 1, 0}"}, {"tB3", "{0, 0, 1}"}, {"tB3", "{x - 1, 1, 0}"}, {"tB3", "{2, 0, x + 3}"},
        {"tB6", "{0, 1}"}, {"tB6", "{-1, 1}"}, {"tB6", "{x^2, x + 1}"},
        {"curveBasis[4, x^3 + 1, x]", "{1, 0, 0, 1}"}, {"curveBasis[5, x^2 + 1, x]", "{0, 0, 1, 0, x}"},
        {"curveBasis[4, (x^2 + 1)^2 (x - 3), x]", "{0, 1, 0, 1}"},
    };
    char buf[1024];
    for (size_t i = 0; i < sizeof gens / sizeof gens[0]; i++) {
        snprintf(buf, sizeof buf,
                 "Module[{tB = %s, tu = %s, tI}, tI = idealFromGenerators[tQ, tB, {tu}];"
                 " tShapeQ[tQ, tI[\"M\"]] && Expand[Det[tI[\"M\"]] - kMonic[tQ, tNormY[tB, tu], x]] === 0]",
                 gens[i][0], gens[i][1]);
        TR_TRUE(buf);
    }
    TR_TRUE("Expand[Det[idealFromGenerators[tQ, tB3, {{0, 1, 0}}][\"M\"]] - x^2 (x + 1)] === 0");      /* N(Y) = G */
    TR_TRUE("Expand[Det[idealFromGenerators[tQ, tB3, {{0, 0, 1}}][\"M\"]] - x (x + 1)^2] === 0");      /* N(Y^2/x) */
    /* two generators: (x, Y) on Y^2 = x is (Y) */
    TR_TRUE("idealEqualQ[tQ, curveBasis[2, x, x], idealFromGenerators[tQ, curveBasis[2, x, x], {{x, 0}, {0, 1}}],"
            " idealFromGenerators[tQ, curveBasis[2, x, x], {{0, 1}}]]");
}

static void test_residueIdeal(void) {
    /* Bronstein Ex. 6: Y^2 = 1 - x^3, r = Y/(x - x^4), E' = 1 - 4x^3, Bp = x.
     * (x, Y - (1 - 4x^3)) = (x, Y - 1): rows {-1, 1}, {x, 0} -> {{1, -1}, {0, x}}. */
    tr_do("tE6 = Expand[D[x - x^4, x]];"
          "tIp = residueIdeal[tQ, tB6, 1, 1, tE6, x]; tIm = residueIdeal[tQ, tB6, -1, 1, tE6, x]; True");
    TR_TRUE("tIp === <|\"M\" -> {{1, -1}, {0, x}}, \"den\" -> 1|>");
    TR_TRUE("tIm === <|\"M\" -> {{1, 1}, {0, x}}, \"den\" -> 1|>");
    /* a residue that does not occur gives the unit ideal (no place) */
    TR_TRUE("tUnitQ[tQ, tB6, residueIdeal[tQ, tB6, 2, 1, tE6, x]]");
    /* Q(Sqrt[2]): Y^2 = x^2 + 1, N = 1, E = x - 1: residue t at (1, Sqrt[2]).
     * 1 - (t/2) Y = -(t/2)(Y - t) lies in the ideal: HNF {{1, -t/2}, {0, x - 1}}. */
    tr_do("tB2 = curveBasis[2, x^2 + 1, x];"
          "tIa = residueIdeal[tN, tB2, t, 1, 1, x - 1]; tIb = residueIdeal[tN, tB2, -t, 1, 1, x - 1]; True");
    TR_TRUE("tIa === <|\"M\" -> {{1, -t/2}, {0, x - 1}}, \"den\" -> 1|>");
    TR_TRUE("tIb === <|\"M\" -> {{1, t/2}, {0, x - 1}}, \"den\" -> 1|>");
    /* Ex. 4 on Y^2 = x^8 + 1: one place each, norm x */
    tr_do("tB4 = curveBasis[2, x^8 + 1, x]; tE4 = Expand[D[x^9 + x, x]];"
          "tJp = residueIdeal[tQ, tB4, 1, 1, tE4, x]; tJm = residueIdeal[tQ, tB4, -1, 1, tE4, x]; True");
    TR_EQ("Det[tJp[\"M\"]]", "x");
    TR_TRUE("tJp[\"M\"] === {{1, -1}, {0, x}}");
}

static void test_ideal_arith(void) {
    /* P + Q = div(x) on Ex. 6 (x = 0 unramified, two places) */
    TR_TRUE("idealTimes[tQ, tB6, tIp, tIm] === <|\"M\" -> {{x, 0}, {0, x}}, \"den\" -> 1|>");
    TR_TRUE("idealEqualQ[tQ, tB6, idealTimes[tQ, tB6, tIp, tIm], idealTimes[tQ, tB6, tIm, tIp]]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealTimes[tQ, tB6, tIp, idealUnit[tB6]], tIp]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealTimes[tQ, tB6, idealTimes[tQ, tB6, tIp, tIp], tIm],"
            " idealTimes[tQ, tB6, tIp, idealTimes[tQ, tB6, tIp, tIm]]]");
    /* powers */
    TR_TRUE("tUnitQ[tQ, tB6, idealPower[tQ, tB6, tIp, 0]]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, 1], tIp]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, 2], idealTimes[tQ, tB6, tIp, tIp]]");
    /* 3P = div(Y - 1) restricted to the finite places: Y - 1 ~ -x^3/2 at P,
     * its only finite zero; so P^3 = (Y - 1) and Q^3 = (Y + 1). */
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, 3], idealFromGenerators[tQ, tB6, {{-1, 1}}]]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIm, 3], idealFromGenerators[tQ, tB6, {{1, 1}}]]");
    TR_FALSE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, 3], idealFromGenerators[tQ, tB6, {{1, 1}}]]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, -1], idealInverse[tQ, tB6, tIp]]");
    TR_TRUE("tUnitQ[tQ, tB6, idealTimes[tQ, tB6, idealPower[tQ, tB6, tIp, -2], idealPower[tQ, tB6, tIp, 2]]]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, -3], idealInverse[tQ, tB6, idealPower[tQ, tB6, tIp, 3]]]");
    TR_TRUE("idealEqualQ[tQ, tB6, idealPower[tQ, tB6, tIp, 5], idealTimes[tQ, tB6, idealPower[tQ, tB6, tIp, 2], idealPower[tQ, tB6, tIp, 3]]]");
    /* inverses */
    TR_TRUE("tInvOK[tQ, tB6, tIp]");
    TR_TRUE("tInvOK[tQ, tB6, tIm]");
    TR_TRUE("tInvOK[tQ, tB6, idealTimes[tQ, tB6, tIp, idealInverse[tQ, tB6, tIm]]]");
    TR_TRUE("tInvOK[tQ, tB4, tJp]");
    TR_TRUE("tInvOK[tQ, tB3, idealFromGenerators[tQ, tB3, {{x, 1, 0}, {x^2, 0, 1}}]]");
    TR_TRUE("tInvOK[tQ, tB3, idealFromGenerators[tQ, tB3, {{0, 1, 0}}]]");
    TR_TRUE("tInvOK[tQ, curveBasis[4, x^3 + 1, x], idealFromGenerators[tQ, curveBasis[4, x^3 + 1, x], {{x + 1, 1, 0, 0}, {x^2 - 1, 0, 0, 0}}]]");
    TR_TRUE("tInvOK[tN, tB2, tIa]");
    TR_TRUE("tInvOK[tN, tB2, idealTimes[tN, tB2, tIa, idealInverse[tN, tB2, tIb]]]");
    TR_TRUE("tInvOK[tG7, basisMod[tB6, 7], residueIdeal[tG7, basisMod[tB6, 7], 1, 1, gfMod[tE6, 7], x]]");
    TR_TRUE("tInvOK[tG5, basisMod[curveBasis[3, x^3 + 2, x], 5], idealFromGenerators[tG5, basisMod[curveBasis[3, x^3 + 2, x], 5], {{x + 1, 2, 0}, {0, x, 3}}]]");
    /* (x)^-1 = (1/x): identity rows over den x */
    TR_TRUE("idealEqualQ[tQ, tB6, idealInverse[tQ, tB6, idealFromGenerators[tQ, tB6, {{x, 0}}]],"
            " <|\"M\" -> IdentityMatrix[2], \"den\" -> x|>]");
    TR_TRUE("tUnitQ[tQ, tB6, idealInverse[tQ, tB6, idealUnit[tB6]]]");
    /* the inverse of the residue ideal over Q(t) has a prime-field denominator */
    TR_TRUE("FreeQ[idealInverse[tN, tB2, tIa][\"den\"], t]");
}

static void test_polyAdjugate(void) {
    TR_EQ("polyAdjugate[tQ, {{x + 1}}, x]", "{{1}}");
    TR_TRUE("tAdjOK[tQ, {{x, 1}, {x^2, 3}}]");
    TR_TRUE("tAdjOK[tQ, {{x, 1, 2}, {0, x^2, 3}, {1, x, 1}}]");
    TR_TRUE("tAdjOK[tQ, {{x, 1, 2, 0}, {0, x^2, 3, 1}, {1, x, 1, x}, {x^3, 0, 0, 1}}]");
    TR_TRUE("polyAdjugate[tQ, {{a1, b1}, {c1, d1}}, x] === {{d1, -b1}, {-c1, a1}}");
    TR_TRUE("tAdjOK[tN, {{t x, 1}, {x^2, t + 1}}]");
    TR_TRUE("tAdjOK[tN, {{t x, 1, 0}, {0, x - t, t}, {x^2, 0, 1 + t x}}]");
    TR_TRUE("tAdjOK[tG7, {{3 x, 1, 5}, {0, x^2, 6}, {1, 4 x, 1}}]");
}

static void test_idealReduce_equal(void) {
    /* common factor x(x - 1) of den and all entries removed */
    TR_TRUE("idealReduce[tQ, tB2, <|\"M\" -> {{x (x - 1), x (x - 1)}, {0, x^2 (x - 1)}}, \"den\" -> x (x - 1)|>] ==="
            " <|\"M\" -> {{1, 1}, {0, x}}, \"den\" -> 1|>");
    /* partial: x divides everything once, x - 1 does not divide the entries */
    TR_TRUE("idealReduce[tQ, tB2, <|\"M\" -> {{x, 0}, {0, x^2}}, \"den\" -> x^2 (x - 1)|>] ==="
            " <|\"M\" -> {{1, 0}, {0, x}}, \"den\" -> Expand[x (x - 1)]|>");
    /* nothing to remove: den made monic, the ideal unchanged */
    tr_do("tIr = idealReduce[tQ, tB2, <|\"M\" -> {{x, 1}, {0, x^2}}, \"den\" -> 2 x|>]; True");
    TR_TRUE("tIr[\"den\"] === x");
    TR_TRUE("idealEqualQ[tQ, tB2, tIr, <|\"M\" -> {{x, 1}, {0, x^2}}, \"den\" -> 2 x|>]");
    TR_TRUE("tUnitQ[tQ, tB2, idealReduce[tQ, tB2, <|\"M\" -> {{3, 0}, {0, 3}}, \"den\" -> 7|>]]");
    TR_TRUE("idealReduce[tG7, basisMod[tB2, 7], <|\"M\" -> {{x + 1, 0}, {0, x^2 + x}}, \"den\" -> x + 1|>] ==="
            " <|\"M\" -> {{1, 0}, {0, x}}, \"den\" -> 1|>");
    /* idealEqualQ */
    TR_TRUE("idealEqualQ[tQ, tB6, <|\"M\" -> {{x, 0}, {0, x}}, \"den\" -> x|>, idealUnit[tB6]]");
    TR_TRUE("idealEqualQ[tQ, tB6, <|\"M\" -> {{x, 0}, {0, x}}, \"den\" -> x^2|>, <|\"M\" -> IdentityMatrix[2], \"den\" -> x|>]");
    TR_FALSE("idealEqualQ[tQ, tB6, tIp, tIm]");
    TR_FALSE("idealEqualQ[tQ, tB6, tIp, idealUnit[tB6]]");
    TR_FALSE("idealEqualQ[tQ, tB6, <|\"M\" -> IdentityMatrix[2], \"den\" -> x|>, idealUnit[tB6]]");
    TR_FALSE("idealEqualQ[tN, tB2, tIa, tIb]");
}

/* ---------------------------------------------------------- PrincipalGen.m */

static void test_principal_genus0(void) {
    /* 1/Sqrt[x^2+1] shifted (a = 0): r = -Y1/(z (1 + z^2)) dz on Y1^2 = 1 + z^2.
     * Residue N(0) Y(P)/E'(0) = -Y(P): +1 at Q = (0, -1), -1 at P = (0, 1).
     * (1 + Y)/z: zero of order 1 at Q (1 + Y = -z^2/2 + ...), pole at P,
     * finite and nonzero at both places over z = Infinity: div = Q - P. */
    tr_do("tBz = curveBasis[2, 1 + z^2, z];"
          "tK1 = residueIdeal[tQ, tBz, 1, -1, 1 + 3 z^2, z]; tK2 = residueIdeal[tQ, tBz, -1, -1, 1 + 3 z^2, z];"
          "tDz = idealTimes[tQ, tBz, tK1, idealInverse[tQ, tBz, tK2]];"
          "tuz = principalGenerator[tQ, tBz, tDz]; True");
    TR_TRUE("ListQ[tuz]");
    TR_TRUE("Module[{tl = Together[tuz[[1]] z]}, FreeQ[tl, z] && tl =!= 0 && Together[tuz[[2]] z] === tl]");
    TR_TRUE("principalQ[tQ, tBz, tDz]");
    TR_TRUE("tGenQ[tQ, tBz, tDz, tuz]");
    /* Q(Sqrt[2]) on Y^2 = x^2 + 1: delta = (1, Sqrt[2]) - (1, -Sqrt[2]).
     * By hand (2 Y - t (1 + x))/(x - 1) vanishes at (1, Sqrt[2]), has a
     * pole at (1, -Sqrt[2]) and tends to 2 - Sqrt[2] or -2 - Sqrt[2] at infinity. */
    tr_do("tDn = idealTimes[tN, tB2, tIa, idealInverse[tN, tB2, tIb]]; tun = principalGenerator[tN, tB2, tDn]; True");
    TR_TRUE("tGenQ[tN, tB2, tDn, tun]");
    TR_TRUE("tPropQ[tN, tun, {-t (1 + x)/(x - 1), 2/(x - 1)}]");
    TR_TRUE("principalQ[tN, tB2, tDn]");
    TR_TRUE("principalQ[tN, tB2, idealPower[tN, tB2, tDn, -3]]");
    /* genus 0, m = 3: Y^3 = x + 1, E = x - 7: residues are Y(P) = 2, 2w, 2w^2 */
    tr_do("tB31 = curveBasis[3, x + 1, x]; tcr = componentResidues[1, x - 7, 3, x + 1, x]; tF3 = tcr[\"F\"];"
          "tL = residueIdeal[tF3, tB31, #, 1, 1, x - 7] & /@ tcr[\"residues\"]; True");
    TR_EQ("curveGenus[3, x + 1, x]", "0");
    TR_EQ("Length[tL]", "3");
    TR_TRUE("AllTrue[Subsets[Range[3], {2}], Function[tp, Module[{tI = idealTimes[tF3, tB31, tL[[tp[[1]]]], idealInverse[tF3, tB31, tL[[tp[[2]]]]]]},"
            " tGenQ[tF3, tB31, tI, principalGenerator[tF3, tB31, tI]]]]]");
    /* the unit ideal (delta = 0): a nonzero constant */
    TR_TRUE("Module[{tu = principalGenerator[tQ, tB6, idealUnit[tB6]]}, tu[[2]] === 0 && NumberQ[tu[[1]]] && tu[[1]] =!= 0]");
    /* genus 0 divisors of the form 2 delta - delta' */
    TR_TRUE("Module[{tI = idealTimes[tQ, tBz, idealPower[tQ, tBz, tK1, 2], idealPower[tQ, tBz, tK2, -2]]},"
            " tGenQ[tQ, tBz, tI, principalGenerator[tQ, tBz, tI]]]");
}

static void test_principal_genus_positive(void) {
    /* Bronstein Ex. 6, genus 1: delta = P - Q; P = (0, 1) is 3-torsion on
     * Y^2 = 1 - x^3 (the tangent Y = 1 meets with multiplicity 3), so
     * delta ~ 2P has order exactly 3:  3 delta = div((Y - 1)/(Y + 1)), and
     * (Y - 1)/(Y + 1) = (Y - 1)^2/(-x^3) = ((x^3 - 2) + 2 Y)/x^3. */
    tr_do("tD6 = idealTimes[tQ, tB6, tIp, idealInverse[tQ, tB6, tIm]]; True");
    TR_EQ("curveGenus[2, 1 - x^3, x]", "1");
    TR_TRUE("principalGenerator[tQ, tB6, tD6] === None");
    TR_FALSE("principalQ[tQ, tB6, tD6]");
    TR_FALSE("principalQ[tQ, tB6, idealPower[tQ, tB6, tD6, 2]]");
    TR_FALSE("principalQ[tQ, tB6, idealPower[tQ, tB6, tD6, -1]]");
    tr_do("tu6 = principalGenerator[tQ, tB6, idealPower[tQ, tB6, tD6, 3]]; True");
    TR_TRUE("tPropQ[tQ, tu6, {(x^3 - 2)/x^3, 2/x^3}]");
    TR_TRUE("tGenQ[tQ, tB6, idealPower[tQ, tB6, tD6, 3], tu6]");
    TR_TRUE("principalQ[tQ, tB6, idealPower[tQ, tB6, tD6, 6]]");
    TR_TRUE("principalQ[tQ, tB6, idealPower[tQ, tB6, tD6, -3]]");
    TR_FALSE("principalQ[tQ, tB6, idealPower[tQ, tB6, tD6, 4]]");
    /* Bronstein Ex. 4, genus 3: delta = P - Q over x = 0 on Y^2 = x^8 + 1;
     * 4 delta = div(x^4/(1 + Y)) and x^4/(1 + Y) = (Y - 1)/x^4 up to sign. */
    tr_do("tD4 = idealTimes[tQ, tB4, tJp, idealInverse[tQ, tB4, tJm]]; True");
    TR_EQ("curveGenus[2, x^8 + 1, x]", "3");
    TR_EQ("Table[principalQ[tQ, tB4, idealPower[tQ, tB4, tD4, k]], {k, 1, 3}]", "{False, False, False}");
    tr_do("tu4 = principalGenerator[tQ, tB4, idealPower[tQ, tB4, tD4, 4]]; True");
    TR_TRUE("tPropQ[tQ, tu4, {-1/x^4, 1/x^4}]");
    TR_TRUE("tGenQ[tQ, tB4, idealPower[tQ, tB4, tD4, 4], tu4]");
    /* The log-derivative identity of Ex. 4: 1/4 u'/u = Y/(x (x^8 + 1)), i.e.
     * u' = 4 (Y/(x(x^8+1))) u in K(x)(Y). */
    TR_TRUE("afZeroQ[tQ, afMinus[tQ, afD[tQ, tB4, tu4], afTimes[tQ, tB4, {0, 4/(x (x^8 + 1))}, tu4]]]");
    /* ... and of Ex. 6: 1/3 u'/u = Y/(x - x^4). */
    TR_TRUE("afZeroQ[tQ, afMinus[tQ, afD[tQ, tB6, tu6], afTimes[tQ, tB6, {0, 3/(x - x^4)}, tu6]]]");
    /* Non-torsion: (3, 5) generates the rank-1 Mordell curve Y^2 = x^3 - 2;
     * delta = (3, 5) - (3, -5) ~ 2 (3, 5) has infinite order. */
    tr_do("tBm = curveBasis[2, x^3 - 2, x];"
          "tDm = idealTimes[tQ, tBm, residueIdeal[tQ, tBm, 5, 1, 1, x - 3], idealInverse[tQ, tBm, residueIdeal[tQ, tBm, -5, 1, 1, x - 3]]]; True");
    TR_EQ("Table[principalQ[tQ, tBm, idealPower[tQ, tBm, tDm, k]], {k, 1, 6}]",
          "{False, False, False, False, False, False}");
}

static void test_principal_GF(void) {
    /* Reduction mod p (good, p not 2, 3) is injective on prime-to-p torsion:
     * Ex. 6 keeps order 3 and Ex. 4 order 4. */
    tr_do("tOrd6[tp_] := Module[{tBp = basisMod[tB6, tp], tF = fieldGF[tp], tI1, tI2, tD},"
          "  tI1 = residueIdeal[tF, tBp, 1, 1, gfMod[tE6, tp], x]; tI2 = residueIdeal[tF, tBp, tp - 1, 1, gfMod[tE6, tp], x];"
          "  tD = idealTimes[tF, tBp, tI1, idealInverse[tF, tBp, tI2]];"
          "  Table[principalQ[tF, tBp, idealPower[tF, tBp, tD, k]], {k, 1, 6}]];"
          "tOrd4[tp_] := Module[{tBp = basisMod[tB4, tp], tF = fieldGF[tp], tI1, tI2, tD},"
          "  tI1 = residueIdeal[tF, tBp, 1, 1, gfMod[tE4, tp], x]; tI2 = residueIdeal[tF, tBp, tp - 1, 1, gfMod[tE4, tp], x];"
          "  tD = idealTimes[tF, tBp, tI1, idealInverse[tF, tBp, tI2]];"
          "  Table[principalQ[tF, tBp, idealPower[tF, tBp, tD, k]], {k, 1, 8}]];"
          /* Mordell curve mod p: delta is principal at #E(F_p) (Lagrange) */
          "tMordell[tp_] := Module[{tBp = basisMod[tBm, tp], tF = fieldGF[tp], tI1, tI2, tD, tn = tCount[x^3 - 2, tp]},"
          "  tI1 = residueIdeal[tF, tBp, 5, 1, 1, x - 3]; tI2 = residueIdeal[tF, tBp, tp - 5, 1, 1, x - 3];"
          "  tD = idealTimes[tF, tBp, tI1, idealInverse[tF, tBp, tI2]];"
          "  {principalQ[tF, tBp, idealPower[tF, tBp, tD, tn]],"
          "   tGenQ[tF, tBp, idealPower[tF, tBp, tD, tn], principalGenerator[tF, tBp, idealPower[tF, tBp, tD, tn]]],"
          "   Module[{to = SelectFirst[Range[tn], principalQ[tF, tBp, idealPower[tF, tBp, tD, #]] &]}, Mod[tn, to] === 0]}];"
          "True");
    TR_EQ("tOrd6[5]", "{False, False, True, False, False, True}");
    TR_EQ("tOrd6[7]", "{False, False, True, False, False, True}");
    TR_EQ("tOrd6[11]", "{False, False, True, False, False, True}");
    TR_EQ("Table[curveGenus[2, x^8 + 1, x, tp], {tp, {5, 7}}]", "{3, 3}");
    TR_EQ("tOrd4[5]", "{False, False, False, True, False, False, False, True}");
    TR_EQ("tOrd4[7]", "{False, False, False, True, False, False, False, True}");
    TR_EQ("tMordell[7]", "{True, True, True}");
    TR_EQ("tMordell[11]", "{True, True, True}");
    TR_EQ("tMordell[13]", "{True, True, True}");
}

/* The closed-form residue ideal, the two-generator divisor ideals and the
 * Hermite form modulo the determinant (added for speed) must agree EXACTLY
 * with the generic constructions they replace: Hermite forms are unique. */
static void test_fast_ideal_machinery(void) {
    char* r = tr_eval(
        "tB6 = curveBasis[6, 1 + x, x]; tKc = fieldNF[t, 36 + 6 t + t^2, None];"
        "tc6 = (-6 - 7 x)/(6 (-1 - x + x^6 + 2 x^7 + x^8)); tN6 = Numerator[tc6];"
        "tE6 = Expand[D[Denominator[tc6], x]]; tBp6 = -6 + 6 x^6 + 6 x^7;"
        "tB2x = curveBasis[2, 1 - x^3, x]; True");
    free(r);
    /* residueIdeal (closed form) == residueIdealGeneric (Hermite form of the 2m
     * generator products), over Q (Bronstein ex. 6), Q(Sqrt[-3]) with m = 6,
     * and a residue that occurs nowhere (unit ideal). */
    TR_TRUE("residueIdeal[tQ, tB2x, 1, 1, D[x (1 - x^3), x], -x] === residueIdealGeneric[tQ, tB2x, 1, 1, D[x (1 - x^3), x], -x]");
    TR_TRUE("residueIdeal[tQ, tB2x, -1, 1, D[x (1 - x^3), x], -x] === residueIdealGeneric[tQ, tB2x, -1, 1, D[x (1 - x^3), x], -x]");
    TR_TRUE("residueIdeal[tKc, tB6, 1/6, tN6, tE6, tBp6] === residueIdealGeneric[tKc, tB6, 1/6, tN6, tE6, tBp6]");
    TR_TRUE("residueIdeal[tKc, tB6, -1/6, tN6, tE6, tBp6] === residueIdealGeneric[tKc, tB6, -1/6, tN6, tE6, tBp6]");
    TR_TRUE("residueIdeal[tQ, tB2x, 5, 1, D[x (1 - x^3), x], -x] === idealUnit[tB2x]");
    /* residueData: A carries the places, phi^m == G mod A, N phi == alpha E' mod A */
    TR_TRUE("Module[{ad = residueData[tKc, tB6, 1/6, tN6, tE6, tBp6]},"
            " kRem[tKc, ad[[2]]^6 - (1 + x), ad[[1]], x] === 0 &&"
            " kRem[tKc, tN6 ad[[2]] - 1/6 tE6, ad[[1]], x] === 0 && kDeg[tKc, ad[[1]], x] > 0]");
    TR_EQ("residueData[tQ, tB2x, 5, 1, D[x (1 - x^3), x], -x]", "{1, 0}");
    /* residueGenerators: the +/- generator ideals are Delta and Delta^-1 */
    TR_TRUE("Module[{g = residueGenerators[tKc, tB6, 1/6, tN6, tE6, tBp6, 1]},"
            " idealEqualQ[tKc, tB6, idealTimesGenerators[tKc, tB6, idealUnit[tB6], Most[g], Last[g]],"
            " residueIdeal[tKc, tB6, 1/6, tN6, tE6, tBp6]]]");
    TR_TRUE("Module[{g = residueGenerators[tKc, tB6, 1/6, tN6, tE6, tBp6, -1]},"
            " idealEqualQ[tKc, tB6, idealTimesGenerators[tKc, tB6, idealUnit[tB6], Most[g], Last[g]],"
            " idealInverse[tKc, tB6, residueIdeal[tKc, tB6, 1/6, tN6, tE6, tBp6]]]]");
    /* the determinant reported with the generators is the module's determinant */
    TR_TRUE("Module[{g = residueGenerators[tKc, tB6, -1/6, tN6, tE6, tBp6, -1], H},"
            " H = hnf[tKc, Flatten[Table[afTimes[tKc, tB6, afUnit[tB6, k], gg], {gg, First[g]}, {k, 0, 5}], 1], x];"
            " kc[tKc, Times @@ Diagonal[H] - Last[g]] === 0]");
    /* divisorIdeal == products of residue ideals and inverses */
    TR_TRUE("idealEqualQ[tKc, tB6, divisorIdeal[tKc, tB6, {1/6, -1/6}, {1, -1}, tN6, tE6, tBp6],"
            " idealTimes[tKc, tB6, residueIdeal[tKc, tB6, 1/6, tN6, tE6, tBp6],"
            " idealInverse[tKc, tB6, residueIdeal[tKc, tB6, -1/6, tN6, tE6, tBp6]]]]");
    TR_TRUE("idealEqualQ[tKc, tB6, divisorIdeal[tKc, tB6, {1/6}, {-2}, tN6, tE6, tBp6],"
            " idealPower[tKc, tB6, residueIdeal[tKc, tB6, 1/6, tN6, tE6, tBp6], -2]]");
    TR_TRUE("idealEqualQ[tQ, tB2x, divisorIdeal[tQ, tB2x, {1, -1}, {3, 0}, 1, D[x (1 - x^3), x], -x],"
            " idealPower[tQ, tB2x, residueIdeal[tQ, tB2x, 1, 1, D[x (1 - x^3), x], -x], 3]]");
    TR_TRUE("divisorIdeal[tQ, tB2x, {1}, {0}, 1, D[x (1 - x^3), x], -x] === idealUnit[tB2x]");
    /* hnfModDet == hnf whenever D is a multiple of the determinant */
    TR_TRUE("Module[{rw = {{x^2 - 1, x + 3}, {x + 1, x^3}, {2 x, 1 + x^2}}, H},"
            " H = hnf[tQ, rw, x];"
            " hnfModDet[tQ, rw, x, Times @@ Diagonal[H]] === H &&"
            " hnfModDet[tQ, rw, x, Expand[(x^5 + 2) Times @@ Diagonal[H]]] === H]");
    TR_TRUE("Module[{rw = {{x - t, x^2}, {1 + t x, x + t}, {x^3, 2}}, H},"
            " H = hnf[tN, rw, x]; hnfModDet[tN, rw, x, kc[tN, Times @@ Diagonal[H]]] === H]");
    TR_TRUE("Module[{rw = {{x + 3, 2 x}, {x^2 + 1, 5}, {1, x^4}}, H},"
            " H = hnf[tG7, rw, x]; hnfModDet[tG7, rw, x, kc[tG7, Times @@ Diagonal[H]]] === H]");
    /* a constant modulus falls back to the plain Hermite form */
    TR_TRUE("hnfModDet[tQ, {{1, x}, {0, x + 1}}, x, 1] === hnf[tQ, {{1, x}, {0, x + 1}}, x]");
}

int main(void) {
    tr_init();
    define_helpers();
    test_hnf_Q();
    test_hnf_NF_GF();
    test_rowInfData();
    test_normalizeAtInfinity();
    test_ideal_basic();
    test_residueIdeal();
    test_ideal_arith();
    test_polyAdjugate();
    test_idealReduce_equal();
    test_principal_genus0();
    test_principal_genus_positive();
    test_principal_GF();
    test_fast_ideal_machinery();
    printf("trager_normalize_ideal: %d assertions passed\n", tr_assertions);
    return 0;
}
