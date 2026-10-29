/* test_trager_residues.c
 *
 * Unit tests for src/internal/trager/Residues.m -- the residues of a Hermite
 * remainder (N/E) Y dx on Y^m = G, the field they generate and an integer
 * Q-basis of their span ([T] Ch 5 Sec 2 pp.56-59; [B] Sec 2.5, Examples 4-6
 * pp.19-20).
 *
 * Every residue expectation is derived by hand from [T] Theorem 2 p.57: at an
 * unramified place P over a simple root x0 of E (G(x0) != 0),
 *     Res_P (N/E) Y dx = N(x0) Y(P) / E'(x0),   Y(P)^m = G(x0),
 * and the residue at a branch place is zero.  Independently of that formula,
 * tPlaceRes checks that the residue at the place with Y(P) = the principal
 * root, computed as Limit-free substitution (x - x0) N/E at x0 times G(x0)^(1/m),
 * is a root of the computed R.
 *
 * Coverage, function by function:
 *   branchFreeDenominator  E / gcd(E, G) up to a constant
 *   residuePolynomial      Bronstein Ex. 4 (Z^2 - 1), Ex. 5 (1: no residues),
 *                          Ex. 6 (Z^2 - 1); m = 3, 4, 5; repeated factor in G;
 *                          zero residues dropped; branch-only denominator;
 *                          N = 0; nontrivial N and E'; monic and squarefree
 *   residueField           rational roots (Q); quadratic Q(Sqrt[2]); mixed
 *                          linear + quadratic; Q(i); Q(Sqrt[2], Sqrt[3]) of
 *                          degree 4; the degree-6 splitting field of Z^3 - 2;
 *                          the a-priori "too large" decline; the roots satisfy
 *                          R in F exactly and numerically under theta
 *   rootObjects            count, Root heads, each a root, pairwise distinct
 *   minimalPolynomialIn    exact values
 *   algebraicToPoly        AlgebraicNumber, rational, the internal throw
 *   qBasis                 alpha_i == Sum b_ij beta_j exactly, integer b,
 *                          independent beta, rank by hand; Galois relations
 *                          r1 + r2 + r3 = 0; the empty list
 *   componentResidues      all keys, consistency of every field, no-residue case
 */
#include "trager_test_utils.h"

/* Evaluate a statement for its side effects. */
static void tr_do(const char* in) { free(tr_eval(in)); }

static void define_helpers(void) {
    tr_do(
        "tQ = fieldQ[];"
        /* Same polynomial up to a nonzero rational constant. */
        "tPropQ[tp_, tq_, tv_] := Module[{tr = Together[tp/tq]}, FreeQ[tr, tv] && tr =!= 0];"
        /* R monic and squarefree in Z. */
        "tMonicSqfQ[tR_, tZ_] := Coefficient[tR, tZ, Exponent[tR, tZ]] === 1 &&"
        "  Exponent[PolynomialGCD[tR, D[tR, tZ]], tZ] === 0;"
        /* The residue at the place over x0 with Y(P) = G(x0)^(1/m) (principal)
         * is a root of R: (x - x0) N/E -> N(x0)/E'(x0). */
        "tPlaceRes[tN_, tE_, tm_, tG_, tx0_, tR_] := Module[{tres},"
        "  tres = (Together[(x - tx0) tN/tE] /. x -> tx0) (tG /. x -> tx0)^(1/tm);"
        "  Abs[N[tR /. Z -> tres, 40]] < 10^-30];"
        /* Field consistency of a residueField result. */
        "tRootsInF[trf_, tR_] := AllTrue[trf[\"roots\"], kc[trf[\"F\"], tR /. Z -> #] === 0 &];"
        "tRootsNumeric[trf_, tR_] := Module[{tF = trf[\"F\"]},"
        "  If[tF[\"Type\"] === \"Q\", True,"
        "    Abs[N[tF[\"mu\"] /. tF[\"t\"] -> tF[\"theta\"], 40]] < 10^-30 &&"
        "    AllTrue[trf[\"roots\"], Abs[N[tR /. Z -> (# /. tF[\"t\"] -> tF[\"theta\"]), 40]] < 10^-30 &]]];"
        "tRootsDistinct[trf_] := Module[{tF = trf[\"F\"], trs = trf[\"roots\"]},"
        "  AllTrue[Subsets[trs, {2}], kc[tF, #[[1]] - #[[2]]] =!= 0 &]];"
        /* qBasis contract: alpha_i == Sum_j b_ij beta_j in F, integer b. */
        "tQBasisOK[tF_, tal_] := Module[{tqb = qBasis[tF, tal], tb, tbe},"
        "  tb = tqb[\"coeffs\"]; tbe = tqb[\"basis\"];"
        "  AllTrue[Flatten[tb], IntegerQ] &&"
        "  AllTrue[Range[Length[tal]], kc[tF, tal[[#]] - tb[[#]] . tbe] === 0 &]];"
        /* rank over Q of a list of field elements (coefficient vectors) */
        "tQRank[tF_, tal_] := If[tF[\"Type\"] === \"Q\", MatrixRank[List /@ tal],"
        "  MatrixRank[PadRight[CoefficientList[kc[tF, #], tF[\"t\"]], Exponent[tF[\"mu\"], tF[\"t\"]]] & /@ tal]];"
        "True");
}

static void test_branchFreeDenominator(void) {
    TR_TRUE("tPropQ[branchFreeDenominator[x (x^2 + 1) (x - 1), (x^2 + 1) (x + 2), x], x (x - 1), x]");
    TR_TRUE("tPropQ[branchFreeDenominator[x - x^4, 1 - x^3, x], x, x]");          /* Ex. 6 */
    TR_TRUE("tPropQ[branchFreeDenominator[x^9 + x, x^8 + 1, x], x, x]");          /* Ex. 4 */
    TR_TRUE("FreeQ[branchFreeDenominator[x^2 + 1, (x^2 + 1) (x + 2), x], x]");     /* all branch */
    TR_TRUE("FreeQ[branchFreeDenominator[1 - x^3, 1 - x^3, x], x]");              /* Ex. 5 */
    TR_TRUE("Expand[branchFreeDenominator[(x - 1) (x - 2), x + 1, x] - (x - 1) (x - 2)] === 0"); /* coprime */
    TR_TRUE("tPropQ[branchFreeDenominator[x (x - 1), x^2 (x + 1), x], x - 1, x]");  /* repeated factor of G */
}

static void test_residuePolynomial(void) {
    /* Bronstein Ex. 4: Y^2 = x^8+1, r = Y/(x (x^8+1)); residues +-1 at (0, +-1). */
    TR_EQ("residuePolynomial[1, x^9 + x, 2, x^8 + 1, x, Z]", "-1 + Z^2");
    /* Ex. 5: x Y/(1 - x^3): E divides G, every place is a branch place -> 1. */
    TR_EQ("residuePolynomial[x, 1 - x^3, 2, 1 - x^3, x, Z]", "1");
    /* Ex. 6: Y/(x - x^4) on Y^2 = 1 - x^3: E'(0) = 1, Y(0) = +-1. */
    TR_EQ("residuePolynomial[1, x - x^4, 2, 1 - x^3, x, Z]", "-1 + Z^2");
    /* m = 4: Y^4 = x^2+1, E = x: residues are the 4th roots of 1. */
    TR_EQ("residuePolynomial[1, x, 4, x^2 + 1, x, Z]", "-1 + Z^4");
    /* m = 5: Y^5 = x + 2, E = x - 30, N = 3: residues 3 * 2 zeta^k: Z^5 = 3^5 32. */
    TR_EQ("residuePolynomial[3, x - 30, 5, x + 2, x, Z]", "-7776 + Z^5");
    /* zero residue at x0 = 0 (N(0) = 0) dropped; at x0 = 1: Y = +-2. */
    TR_EQ("residuePolynomial[x, x (x - 1), 2, x + 3, x, Z]", "-4 + Z^2");
    /* repeated factor in G: Y^3 = x^2(x+1); x = 0 is a branch point, at x0 = 1 Y^3 = 2. */
    TR_EQ("residuePolynomial[1, x (x - 1), 3, x^2 (x + 1), x, Z]", "-2 + Z^3");
    TR_EQ("residuePolynomial[1, x - 1, 3, x^2 (x + 1), x, Z]", "-2 + Z^3");
    /* irreducible E: Y^3 = x, E = x^2+1: Z^3 = x0/(2 x0)^3 = 1/(8 x0^2) = -1/8 at both roots. */
    TR_EQ("residuePolynomial[1, x^2 + 1, 3, x, x, Z]", "1/8 + Z^3");
    /* Y^2 = x^2+1, E = x^3 - x: +-1 over 0, +-Sqrt[2]/2 over +-1. */
    TR_EQ("residuePolynomial[1, x^3 - x, 2, x^2 + 1, x, Z]", "1/2 - 3/2 Z^2 + Z^4");
    /* N = x + 2, E = x^2 - 1 on Y^2 = x^2+1: Z^2 = 9*2/4 over 1, 1*2/4 over -1. */
    TR_TRUE("Expand[residuePolynomial[x + 2, x^2 - 1, 2, x^2 + 1, x, Z] - (Z^2 - 9/2) (Z^2 - 1/2)] === 0");
    /* non-monic E and constant N: E = 2x - 2, Y(1) = +-2 on Y^2 = x + 3:
     * residues 1 * (+-2)/2 = +-1; N = 6, E = 3x on Y^2 = x + 4: 6 (+-2)/3 = +-4. */
    TR_EQ("residuePolynomial[1, 2 x - 2, 2, x + 3, x, Z]", "-1 + Z^2");
    TR_EQ("residuePolynomial[6, 3 x, 2, x + 4, x, Z]", "-16 + Z^2");
    /* non-monic G: Y^3 = 2 x^2 + 6, E = x - 1: Z^3 = 8 */
    TR_EQ("residuePolynomial[1, x - 1, 3, 2 x^2 + 6, x, Z]", "-8 + Z^3");
    /* branch places only / N = 0 */
    TR_EQ("residuePolynomial[1, x, 2, x (x + 1), x, Z]", "1");
    TR_EQ("residuePolynomial[0, x, 2, x + 1, x, Z]", "1");
    TR_EQ("residuePolynomial[5, 7, 2, x + 1, x, Z]", "1");
    /* Bp = x + 1 (x - 1 is a branch factor): E = x^2 - 1 on Y^2 = x - 1; Y(-1)^2 = -2, E'(-1) = -2. */
    TR_EQ("residuePolynomial[1, x^2 - 1, 2, x - 1, x, Z]", "1/2 + Z^2");
    /* Independent check at the principal place over each root of Bp. */
    TR_TRUE("tPlaceRes[x + 2, x^2 - 1, 2, x^2 + 1, 1, residuePolynomial[x + 2, x^2 - 1, 2, x^2 + 1, x, Z]]");
    TR_TRUE("tPlaceRes[x + 2, x^2 - 1, 2, x^2 + 1, -1, residuePolynomial[x + 2, x^2 - 1, 2, x^2 + 1, x, Z]]");
    TR_TRUE("tPlaceRes[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, 2, residuePolynomial[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, x, Z]]");
    TR_TRUE("tPlaceRes[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, 5, residuePolynomial[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, x, Z]]");
    TR_TRUE("tPlaceRes[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, -1, residuePolynomial[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, x, Z]]");
    TR_TRUE("tPlaceRes[1, x (x - 3), 4, x^3 + 1, 3, residuePolynomial[1, x (x - 3), 4, x^3 + 1, x, Z]]");
    TR_TRUE("tPlaceRes[2 x, (x - 1) (x - 2), 5, x^2 + x + 1, 2, residuePolynomial[2 x, (x - 1) (x - 2), 5, x^2 + x + 1, x, Z]]");
    /* (x^2+3)/((x-2)(x-5)(x+1)) Y on Y^3 = x^2+7: degree 3 * 3 distinct x0. */
    TR_EQ("Exponent[residuePolynomial[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, x, Z], Z]", "9");
    /* monic squarefree in general */
    TR_TRUE("tMonicSqfQ[residuePolynomial[x^2 + 3, (x - 2) (x - 5) (x + 1), 3, x^2 + 7, x, Z], Z]");
    TR_TRUE("tMonicSqfQ[residuePolynomial[3 x - 1, x^3 + x + 1, 2, x^5 + 2, x, Z], Z]");
    TR_TRUE("tMonicSqfQ[residuePolynomial[1, x^3 - x, 2, x^2 + 1, x, Z], Z]");
    /* a zero root never survives */
    TR_TRUE("Coefficient[residuePolynomial[x (x - 2), x (x - 1) (x - 2) (x - 3), 2, x + 7, x, Z], Z, 0] =!= 0");
}

static void test_rootObjects(void) {
    TR_EQ("Length[rootObjects[Z^3 - 2, Z]]", "3");
    TR_TRUE("AllTrue[rootObjects[Z^3 - 2, Z], Head[#] === Root &]");
    TR_TRUE("AllTrue[rootObjects[Z^3 - 2, Z], Abs[N[#^3 - 2, 30]] < 10^-25 &]");
    TR_TRUE("Length[Union[Round[10^8 N[rootObjects[Z^3 - 2, Z], 20]]]] === 3");
    TR_EQ("Length[rootObjects[Z^2 + Z + 1, Z]]", "2");
    TR_TRUE("AllTrue[rootObjects[Z^2 + Z + 1, Z], Abs[N[#^2 + # + 1, 30]] < 10^-25 &]");
    TR_TRUE("AllTrue[rootObjects[Z^5 - Z - 1, Z], Abs[N[#^5 - # - 1, 30]] < 10^-25 &]");
    TR_EQ("Length[rootObjects[Z^5 - Z - 1, Z]]", "5");
}

static void test_minimalPolynomialIn(void) {
    TR_EQ("minimalPolynomialIn[Sqrt[2] + Sqrt[3], t]", "1 - 10 t^2 + t^4");
    TR_EQ("minimalPolynomialIn[3/2, t]", "-3/2 + t");
    TR_EQ("minimalPolynomialIn[Root[#^3 - 2 &, 1]/2, t]", "-1/4 + t^3");
    TR_EQ("minimalPolynomialIn[Sqrt[2]/3, t]", "-2/9 + t^2");
    TR_EQ("minimalPolynomialIn[I, t]", "1 + t^2");
}

static void test_algebraicToPoly(void) {
    TR_EQ("algebraicToPoly[AlgebraicNumber[Sqrt[2], {1, 3}], Sqrt[2], t]", "1 + 3 t");
    TR_EQ("algebraicToPoly[AlgebraicNumber[Root[#^3 - 2 &, 1], {0, 0, -1/2}], Root[#^3 - 2 &, 1], t]", "-1/2 t^2");
    TR_EQ("algebraicToPoly[-7/3, Sqrt[2], t]", "-7/3");
    TR_EQ("algebraicToPoly[5, Sqrt[2], t]", "5");
    TR_TRUE("catchTrager[algebraicToPoly[Sqrt[2], Sqrt[2], t]][[1 ;; 2]] === tragerFail[\"Internal\", \"AlgebraicToPolyUnexpected\"]");
}

static void test_residueField(void) {
    /* rational residues */
    tr_do("tRF1 = residueField[(Z - 1) (Z + 2), Z, t]; True");
    TR_TRUE("tRF1[\"F\"] === fieldQ[]");
    TR_EQ("Sort[tRF1[\"roots\"]]", "{-2, 1}");
    TR_EQ("Sort[residueField[Z, Z, t][\"roots\"]]", "{0}");
    /* quadratic */
    tr_do("tRF2 = residueField[Z^2 - 2, Z, t]; True");
    TR_TRUE("tRF2[\"F\"][\"Type\"] === \"NF\"");
    TR_EQ("fieldDegree[tRF2[\"F\"]]", "2");
    TR_EQ("Length[tRF2[\"roots\"]]", "2");
    TR_TRUE("tRootsInF[tRF2, Z^2 - 2]");
    TR_TRUE("tRootsNumeric[tRF2, Z^2 - 2]");
    TR_TRUE("tRootsDistinct[tRF2]");
    TR_TRUE("kc[tRF2[\"F\"], Total[tRF2[\"roots\"]]] === 0");
    /* mixed linear + quadratic */
    tr_do("tRF3 = residueField[(Z - 1) (Z^2 - 2), Z, t]; True");
    TR_EQ("fieldDegree[tRF3[\"F\"]]", "2");
    TR_EQ("Length[tRF3[\"roots\"]]", "3");
    TR_TRUE("tRootsInF[tRF3, (Z - 1) (Z^2 - 2)]");
    TR_TRUE("tRootsNumeric[tRF3, (Z - 1) (Z^2 - 2)]");
    TR_TRUE("tRootsDistinct[tRF3]");
    TR_TRUE("MemberQ[tRF3[\"roots\"], 1]");
    /* Q(i) */
    tr_do("tRF4 = residueField[Z^4 - 1, Z, t]; True");
    TR_EQ("fieldDegree[tRF4[\"F\"]]", "2");
    TR_TRUE("tRootsInF[tRF4, Z^4 - 1] && tRootsNumeric[tRF4, Z^4 - 1] && tRootsDistinct[tRF4]");
    TR_EQ("Length[tRF4[\"roots\"]]", "4");
    /* Q(Sqrt[2], Sqrt[3]), degree 4 */
    tr_do("tRF5 = residueField[(Z^2 - 2) (Z^2 - 3), Z, t]; True");
    TR_EQ("fieldDegree[tRF5[\"F\"]]", "4");
    TR_TRUE("tRootsInF[tRF5, (Z^2 - 2) (Z^2 - 3)]");
    TR_TRUE("tRootsNumeric[tRF5, (Z^2 - 2) (Z^2 - 3)]");
    TR_TRUE("tRootsDistinct[tRF5]");
    TR_TRUE("Exponent[tRF5[\"F\"][\"mu\"], tRF5[\"F\"][\"t\"]] === 4 && Coefficient[tRF5[\"F\"][\"mu\"], t, 4] === 1");
    /* the splitting field of Z^3 - 2 has degree 6 */
    tr_do("tRF6 = residueField[Z^3 - 2, Z, t]; True");
    TR_EQ("fieldDegree[tRF6[\"F\"]]", "6");
    TR_TRUE("tRootsInF[tRF6, Z^3 - 2] && tRootsNumeric[tRF6, Z^3 - 2] && tRootsDistinct[tRF6]");
    TR_TRUE("kc[tRF6[\"F\"], Total[tRF6[\"roots\"]]] === 0");              /* r1 + r2 + r3 = 0 */
    TR_TRUE("kc[tRF6[\"F\"], Times @@ tRF6[\"roots\"]] === 2");            /* r1 r2 r3 = 2 */
    /* cube roots of unity */
    tr_do("tRF7 = residueField[Z^3 - 1, Z, t]; True");
    TR_EQ("fieldDegree[tRF7[\"F\"]]", "2");
    TR_TRUE("tRootsInF[tRF7, Z^3 - 1] && tRootsNumeric[tRF7, Z^3 - 1] && tRootsDistinct[tRF7]");
    /* no degree limit: the splitting field of an S4 quartic (degree 24) by
     * the algebraic tower, every root a root, all distinct */
    tr_do("tRF8 = residueFieldTower[Z^4 - Z - 1, Z, t]; True");
    TR_EQ("fieldDegree[tRF8[\"F\"]]", "24");
    TR_TRUE("tRootsInF[tRF8, Z^4 - Z - 1] && tRootsNumeric[tRF8, Z^4 - Z - 1] && tRootsDistinct[tRF8]");
    /* the octic of 1/((x^2 + 2)(x^3 + 1)^(1/4)), where ToNumberField gives up:
     * Q(gamma, i, gamma') of degree 16 */
    tr_do("tRF9 = residueField[Z^8 - Z^4/288 + 1/36864, Z, t]; True");
    TR_EQ("fieldDegree[tRF9[\"F\"]]", "16");
    TR_TRUE("tRootsInF[tRF9, Z^8 - Z^4/288 + 1/36864] && tRootsDistinct[tRF9]");
    /* tower helpers: factoring over K (Trager's norm method), extension */
    TR_TRUE("Module[{K = fieldNF[u, u^3 - 2, None], fs},"
            " fs = factorOverField[K, Z^3 - 2, Z];"
            " Sort[Exponent[#, Z] & /@ fs] === {1, 2} && kc[K, Times @@ fs - (Z^3 - 2)] === 0]");
    TR_TRUE("Module[{K = fieldNF[u, u^2 + 3, None], fs},"
            " fs = factorOverField[K, Z^2 + Z + 1, Z];"
            " Exponent[#, Z] & /@ fs === {1, 1} && kc[K, Times @@ fs - (Z^2 + Z + 1)] === 0]");
    TR_TRUE("Module[{K = fieldNF[u, u^3 - 2, None], h, K2},"
            " h = SelectFirst[factorOverField[K, Z^3 - 2, Z], Exponent[#, Z] == 2 &];"
            " K2 = extendFieldByRoot[K, h, Z];"
            " fieldDegree[K2] === 6 && kc[K2, K2[\"embed\"]^3 - 2] === 0]");
}

static void test_qBasis(void) {
    /* over Q: span of {1/2, 1/3, -5/6} is Q (1/6) with coordinates 3, 2, -5 */
    TR_TRUE("qBasis[tQ, {1/2, 1/3, -5/6}] === <|\"basis\" -> {1/6}, \"coeffs\" -> {{3}, {2}, {-5}}|>");
    TR_TRUE("qBasis[tQ, {1, -1}] === <|\"basis\" -> {1}, \"coeffs\" -> {{1}, {-1}}|>");
    TR_TRUE("qBasis[tQ, {}] === <|\"basis\" -> {}, \"coeffs\" -> {}|>");
    TR_TRUE("tQBasisOK[tQ, {2/3, -4/9, 7}]");
    TR_EQ("Length[qBasis[tQ, {2/3, -4/9, 7}][\"basis\"]]", "1");
    /* number fields: ranks by hand */
    const char* nf[][2] = {
        {"tRF2", "1"},   /* +-Sqrt[2] */
        {"tRF3", "2"},   /* 1, +-Sqrt[2] */
        {"tRF4", "2"},   /* +-1, +-i */
        {"tRF5", "2"},   /* +-Sqrt[2], +-Sqrt[3] */
        {"tRF6", "2"},   /* cube roots of 2: one relation */
        {"tRF7", "2"},   /* 1, w, w^2 with 1 + w + w^2 = 0 */
    };
    char buf[512];
    for (size_t i = 0; i < sizeof nf / sizeof nf[0]; i++) {
        snprintf(buf, sizeof buf, "tQBasisOK[%s[\"F\"], %s[\"roots\"]]", nf[i][0], nf[i][0]);
        TR_TRUE(buf);
        snprintf(buf, sizeof buf, "Length[qBasis[%s[\"F\"], %s[\"roots\"]][\"basis\"]]", nf[i][0], nf[i][0]);
        TR_EQ(buf, nf[i][1]);
        snprintf(buf, sizeof buf, "Module[{tqb = qBasis[%s[\"F\"], %s[\"roots\"]]}, tQRank[%s[\"F\"], tqb[\"basis\"]] === Length[tqb[\"basis\"]]]",
                 nf[i][0], nf[i][0], nf[i][0]);
        TR_TRUE(buf);
        snprintf(buf, sizeof buf, "tQRank[%s[\"F\"], %s[\"roots\"]]", nf[i][0], nf[i][0]);
        TR_EQ(buf, nf[i][1]);
    }
    /* +-Sqrt[2], +-2 Sqrt[2]: one Q-direction, integer coordinates +-1, +-2 */
    TR_TRUE("Module[{trf = residueField[(Z^2 - 2) (Z^2 - 8), Z, t], tqb}, tqb = qBasis[trf[\"F\"], trf[\"roots\"]];"
            " fieldDegree[trf[\"F\"]] === 2 && Length[tqb[\"basis\"]] === 1 &&"
            " Sort[Abs[Flatten[tqb[\"coeffs\"]]]] === {1, 1, 2, 2} && tQBasisOK[trf[\"F\"], trf[\"roots\"]]]");
    /* Y^2 = x, E = x^2 - 2: Z^2 = 1/(4 x0), so Z^4 = 1/32 */
    TR_EQ("residuePolynomial[1, x^2 - 2, 2, x, x, Z]", "-1/32 + Z^4");
    /* Galois relations: the residues sum to zero, so the coefficient rows do */
    TR_TRUE("Total[qBasis[tRF6[\"F\"], tRF6[\"roots\"]][\"coeffs\"]] === {0, 0}");
    TR_TRUE("Total[qBasis[tRF7[\"F\"], tRF7[\"roots\"]][\"coeffs\"]] === {0, 0}");
    TR_TRUE("Total[qBasis[tRF5[\"F\"], tRF5[\"roots\"]][\"coeffs\"]] === {0, 0}");
    /* scaled residues: alpha = Sqrt[2]/4 and 3 Sqrt[2]/2 -> one beta, integer coords */
    TR_TRUE("Module[{tF = tRF2[\"F\"], ts = tRF2[\"roots\"][[1]]}, tQBasisOK[tF, {ts/4, 3 ts/2, -ts/6}] &&"
            " Length[qBasis[tF, {ts/4, 3 ts/2, -ts/6}][\"basis\"]] === 1]");
    TR_TRUE("Module[{tF = tRF5[\"F\"], trs = tRF5[\"roots\"]}, tQBasisOK[tF, Join[trs/3, {1/2}]] &&"
            " Length[qBasis[tF, Join[trs/3, {1/2}]][\"basis\"]] === 3]");
    /* the basis elements are rational multiples of residues */
    TR_TRUE("Module[{tF = tRF6[\"F\"], tqb = qBasis[tRF6[\"F\"], tRF6[\"roots\"]]},"
            " AllTrue[tqb[\"basis\"], Function[tb, AnyTrue[tRF6[\"roots\"], tQRank[tF, {tb, #}] === 1 &]]]]");
}

static void test_componentResidues(void) {
    /* Bronstein Ex. 6 */
    tr_do("tCR6 = componentResidues[1, x - x^4, 2, 1 - x^3, x]; True");
    TR_TRUE("Sort[Keys[tCR6]] === Sort[{\"R\", \"Bp\", \"F\", \"residues\", \"basis\", \"coeffs\", \"Z\", \"t\"}]");
    TR_TRUE("Expand[tCR6[\"R\"] - (tCR6[\"Z\"]^2 - 1)] === 0");
    TR_TRUE("tPropQ[tCR6[\"Bp\"], x, x]");
    TR_TRUE("tCR6[\"F\"] === fieldQ[]");
    TR_EQ("Sort[tCR6[\"residues\"]]", "{-1, 1}");
    TR_EQ("Length[tCR6[\"basis\"]]", "1");
    TR_TRUE("AllTrue[Range[2], tCR6[\"residues\"][[#]] === tCR6[\"coeffs\"][[#]] . tCR6[\"basis\"] &]");
    /* Bronstein Ex. 4 */
    tr_do("tCR4 = componentResidues[1, x^9 + x, 2, x^8 + 1, x]; True");
    TR_EQ("Sort[tCR4[\"residues\"]]", "{-1, 1}");
    TR_TRUE("tPropQ[tCR4[\"Bp\"], x, x]");
    /* Bronstein Ex. 5: no residues */
    tr_do("tCR5 = componentResidues[x, 1 - x^3, 2, 1 - x^3, x]; True");
    TR_EQ("tCR5[\"R\"]", "1");
    TR_EQ("tCR5[\"residues\"]", "{}");
    TR_EQ("tCR5[\"basis\"]", "{}");
    TR_TRUE("tCR5[\"F\"] === fieldQ[]");
    /* residues +-Sqrt[2], +-Sqrt[3] over (x-1)(x-2) on Y^2 = x + 1 */
    tr_do("tCRq = componentResidues[1, (x - 1) (x - 2), 2, x + 1, x]; True");
    TR_TRUE("Expand[tCRq[\"R\"] - ((tCRq[\"Z\"]^2 - 2) (tCRq[\"Z\"]^2 - 3))] === 0");
    TR_EQ("fieldDegree[tCRq[\"F\"]]", "4");
    TR_EQ("Length[tCRq[\"residues\"]]", "4");
    TR_TRUE("AllTrue[tCRq[\"residues\"], kc[tCRq[\"F\"], tCRq[\"R\"] /. tCRq[\"Z\"] -> #] === 0 &]");
    TR_TRUE("AllTrue[Range[4], kc[tCRq[\"F\"], tCRq[\"residues\"][[#]] - tCRq[\"coeffs\"][[#]] . tCRq[\"basis\"]] === 0 &]");
    TR_TRUE("tCRq[\"F\"][\"t\"] === tCRq[\"t\"]");
    TR_TRUE("Expand[tCRq[\"Bp\"] - (x - 1) (x - 2)] === 0");
    /* m = 3: cube roots of 2 over x = 0 on Y^3 = x^3 + 2 */
    tr_do("tCR3 = componentResidues[1, x, 3, x^3 + 2, x]; True");
    TR_EQ("fieldDegree[tCR3[\"F\"]]", "6");
    TR_TRUE("Total[tCR3[\"coeffs\"]] === {0, 0}");
    TR_TRUE("AllTrue[Range[3], kc[tCR3[\"F\"], tCR3[\"residues\"][[#]] - tCR3[\"coeffs\"][[#]] . tCR3[\"basis\"]] === 0 &]");
    /* the Z, t symbols are fresh each call */
    TR_TRUE("componentResidues[1, x, 3, x^3 + 2, x][\"Z\"] =!= tCR3[\"Z\"]");
}

int main(void) {
    tr_init();
    define_helpers();
    test_branchFreeDenominator();
    test_residuePolynomial();
    test_rootObjects();
    test_minimalPolynomialIn();
    test_algebraicToPoly();
    test_residueField();
    test_qBasis();
    test_componentResidues();
    printf("trager_residues: %d assertions passed\n", tr_assertions);
    return 0;
}
