/* test_trager_basis_arith.c
 *
 * Unit tests for src/internal/trager/Basis.m (curveBasis: the integral basis
 * w_i = Y^i/d_i of Q(x)[Y]/(Y^m - G), its logarithmic derivatives s_i, the
 * exponents delta_i at infinity, the multiplication table) and
 * src/internal/trager/Arithmetic.m (element arithmetic on that basis: afZero,
 * afOne, afUnit, afCanon, afZeroQ, afPlus, afMinus, afScale, afTimes,
 * afPower, afD, afFromY, afToY, afMultMatrix, afNorm, afInverse,
 * afPolynomialQ, basisMod).
 *
 * Every later stage (Hermite reduction, residues, ideals, the principality
 * test) computes on this representation, so the suite pins down the basis
 * from its DEFINING properties rather than its formula:
 *   - w_i = Y^i/d_i is integral (d_i^m | G^i) and d_i is maximal (for every
 *     prime factor p of G, (p d_i)^m does not divide G^i), d_i | d_(i+1);
 *   - the multiplication factors are polynomials and w_i w_j = factor w_k
 *     holds as an identity in Q(x)[Y]/(Y^m - G);
 *   - s_i is the logarithmic derivative of w_i under Y' = G' Y/(m G);
 *   - delta_i - 1 < i deg G/m - deg d_i <= delta_i (Y^i/d_i has degree
 *     i deg G/m - deg d_i at infinity).
 * The arithmetic is checked by ring and derivation laws over seven curves
 * (m = 2..6, constant factors and repeated factors): y * y^(m-1) = G,
 * commutativity and associativity of afTimes, the Leibniz rule, afD against
 * the derivative of afToY, multiplicativity of the norm and N(Y) =
 * (-1)^(m-1) G, a * afInverse[a] = 1, afPower against repeated afTimes, and
 * basisMod reducing products and derivatives consistently with reduction of
 * the results computed over Q.  Hand-worked values for y^2 = x^3 + 1,
 * y^3 = x^2 (x + 1) and y^6 = x^3 (x + 1)^2 anchor the exact form.
 *
 * Top-level names carry a `t` prefix so they cannot collide with
 * Trager`Private` symbols, which are on $ContextPath.  The symbol tYv stands
 * for Y in afToY / afFromY.
 */
#include "trager_test_utils.h"

#define NCURVES 7

static void setup_basis(void) {
    char* s = tr_eval(
        "tQ = fieldQ[]; tN2 = fieldNF[t, t^2 - 2, Sqrt[2]]; tG5 = fieldGF[5]; tG7 = fieldGF[7];"
        "tCurves = {{2, x^3 + 1}, {2, x (x - 1) (x + 2)}, {3, x^2 (x + 1)}, {3, 2 x^2 (x - 1)},"
        "  {4, x^3 (x^2 + 1)}, {5, x^2 (x - 1)^3}, {6, x^3 (x + 1)^2}};"
        "tBs = curveBasis[#[[1]], #[[2]], x] & /@ tCurves;"
        /* reduce an expression in tYv modulo tYv^m - G */
        "tRedY[tB_, e_] := Module[{tt = Together[e]},"
        "  Together[PolynomialRemainder[Expand[Numerator[tt]], tYv^tB[\"m\"] - tB[\"G\"], tYv]/Denominator[tt]]];"
        /* the derivation of Q(x)(Y): Y' = G' Y/(m G) */
        "tDY[tB_, e_] := D[e, x] + D[e, tYv] D[tB[\"G\"], x] tYv/(tB[\"m\"] tB[\"G\"]);"
        /* deterministic test elements with rational-function coordinates */
        "tElt[tB_, k_] := Table[(i + k + 1 + x^Mod[i k + 1, 3])/(x + i + k + 3), {i, 0, tB[\"m\"] - 1}];"
        "tPolyElt[tB_, k_] := Table[x^Mod[i + k, 3] + i - k, {i, 0, tB[\"m\"] - 1}];"
        "tSame[F_, ta_, tb_] := afZeroQ[F, afMinus[F, ta, tb]];"
        "True");
    free(s);
}

/* ---- Basis.m: hand values -------------------------------------------------- */
static void test_basis_values(void) {
    /* y^2 = x^3 + 1: d = {1, 1}, delta_1 = Ceiling[3/2] = 2 */
    TR_TRUE("Module[{tb = curveBasis[2, x^3 + 1, x]}, tb[\"d\"] === {1, 1} && tb[\"delta\"] === {0, 2} &&"
            " tb[\"deg\"] === 3 && tb[\"m\"] === 2 && tb[\"x\"] === x &&"
            " Together[tb[\"s\"] - {0, 3 x^2/(2 (1 + x^3))}] === {0, 0} &&"
            " tb[\"mult\"] === {{{0, 1}, {1, 1}}, {{1, 1}, {0, 1 + x^3}}}]");
    TR_TRUE("Sort[curveBasis[2, x^3 + 1, x][\"factors\"]] === Sort[{{1 + x, 1}, {1 - x + x^2, 1}}]");
    /* y^3 = x^2 (x + 1): d = {1, 1, x}; delta = {0, 1, 1};
       s_1 = (3x + 2)/(3x(x + 1)), s_2 = (3x + 1)/(3x(x + 1));
       w1 w1 = x w2, w1 w2 = x(x + 1) w0, w2 w2 = (x + 1) w1 */
    TR_TRUE("Module[{tb = curveBasis[3, x^2 (x + 1), x]}, tb[\"d\"] === {1, 1, x} && tb[\"delta\"] === {0, 1, 1} &&"
            " Together[tb[\"s\"] - {0, (3 x + 2)/(3 x (x + 1)), (3 x + 1)/(3 x (x + 1))}] === {0, 0, 0} &&"
            " tb[\"mult\"][[2, 2]] === {2, x} && tb[\"mult\"][[2, 3, 1]] === 0 &&"
            " Expand[tb[\"mult\"][[2, 3, 2]]] === x + x^2 && tb[\"mult\"][[3, 3, 1]] === 1 &&"
            " Expand[tb[\"mult\"][[3, 3, 2]]] === 1 + x]");
    /* y^6 = x^3 (x + 1)^2: d_i = x^Floor[i/2] (x + 1)^Floor[i/3] */
    TR_TRUE("Module[{tb = curveBasis[6, x^3 (x + 1)^2, x]},"
            " Expand[tb[\"d\"] - {1, 1, x, x (x + 1), x^2 (x + 1), x^2 (x + 1)}] === {0, 0, 0, 0, 0, 0} &&"
            " tb[\"delta\"] === {0, 1, 1, 1, 1, 2} && tb[\"deg\"] === 5 &&"
            " Together[tb[\"s\"][[2]] - (5 x + 3)/(6 x (x + 1))] === 0]");
    /* constant factor: y^3 = 2x^2(x - 1) -> d = {1, 1, x} and G stored expanded */
    TR_TRUE("Module[{tb = curveBasis[3, 2 x^2 (x - 1), x]}, tb[\"d\"] === {1, 1, x} &&"
            " tb[\"G\"] === -2 x^2 + 2 x^3 && Sort[tb[\"factors\"]] === Sort[{{x, 2}, {-1 + x, 1}}]]");
    /* y^4 = x^3 (x^2 + 1): d_i = x^Floor[3i/4] -> {1, 1, x, x^2}; M = 5:
       delta = {0, 2 - 0, 3 - 1, 4 - 2} */
    TR_TRUE("Module[{tb = curveBasis[4, x^3 (x^2 + 1), x]}, tb[\"d\"] === {1, 1, x, x^2} && tb[\"delta\"] === {0, 2, 2, 2}]");
    /* y^5 = x^2 (x - 1)^3: d_i = x^Floor[2i/5] (x-1)^Floor[3i/5]:
       {1, 1, x - 1, x (x - 1), x (x - 1)^2}; M = 5: delta_i = i - deg d_i */
    TR_TRUE("Module[{tb = curveBasis[5, x^2 (x - 1)^3, x]},"
            " Expand[tb[\"d\"] - {1, 1, x - 1, x (x - 1), x (x - 1)^2}] === {0, 0, 0, 0, 0} &&"
            " tb[\"delta\"] === {0, 1, 1, 1, 1}]");
}

/* ---- Basis.m: defining properties, per curve -------------------------- */
static void test_basis_properties(void) {
    char buf[2048];
    for (int c = 1; c <= NCURVES; c++) {
#define CHK(fmt) do { snprintf(buf, sizeof buf, "With[{tB = tBs[[%d]]}, " fmt "]", c); TR_TRUE(buf); } while (0)
        /* d_0 = d_1 = 1 */
        CHK("tB[\"d\"][[1]] === 1 && tB[\"d\"][[2]] === 1");
        /* integrality: d_i^m | G^i */
        CHK("AllTrue[Range[0, tB[\"m\"] - 1], PolynomialRemainder[Expand[tB[\"G\"]^#], Expand[tB[\"d\"][[# + 1]]^tB[\"m\"]], x] === 0 &]");
        /* maximality: (p d_i)^m does not divide G^i for any prime factor p */
        CHK("AllTrue[Range[0, tB[\"m\"] - 1], Function[ti, AllTrue[tB[\"factors\"],"
            " PolynomialRemainder[Expand[tB[\"G\"]^ti], Expand[(#[[1]] tB[\"d\"][[ti + 1]])^tB[\"m\"]], x] =!= 0 &]]]");
        /* d_i | d_(i+1) */
        CHK("AllTrue[Range[1, tB[\"m\"] - 1], PolynomialRemainder[tB[\"d\"][[# + 1]], tB[\"d\"][[#]], x] === 0 &]");
        /* every multiplication factor is a polynomial, index k = (i + j) mod m */
        CHK("AllTrue[Flatten[Table[{tB[\"mult\"][[i + 1, j + 1, 1]] === Mod[i + j, tB[\"m\"]],"
            " PolynomialQ[Together[tB[\"mult\"][[i + 1, j + 1, 2]]], x]},"
            " {i, 0, tB[\"m\"] - 1}, {j, 0, tB[\"m\"] - 1}]], TrueQ]");
        /* the table is symmetric */
        CHK("tB[\"mult\"] === Transpose[tB[\"mult\"]]");
        /* w_i w_j = factor w_k as an identity modulo Y^m - G */
        CHK("AllTrue[Flatten[Table[With[{tk = tB[\"mult\"][[i + 1, j + 1]]},"
            " tRedY[tB, tYv^i/tB[\"d\"][[i + 1]] tYv^j/tB[\"d\"][[j + 1]] - tk[[2]] tYv^tk[[1]]/tB[\"d\"][[tk[[1]] + 1]]] === 0],"
            " {i, 0, tB[\"m\"] - 1}, {j, 0, tB[\"m\"] - 1}]], TrueQ]");
        /* s_i is the logarithmic derivative of w_i */
        CHK("AllTrue[Range[0, tB[\"m\"] - 1], Together[tDY[tB, tYv^#/tB[\"d\"][[# + 1]]]/(tYv^#/tB[\"d\"][[# + 1]]) - tB[\"s\"][[# + 1]]] === 0 &]");
        /* delta_i bounds the degree at infinity of w_i */
        CHK("AllTrue[Range[0, tB[\"m\"] - 1], With[{tv = # tB[\"deg\"]/tB[\"m\"] - Exponent[tB[\"d\"][[# + 1]], x]},"
            " tB[\"delta\"][[# + 1]] - 1 < tv <= tB[\"delta\"][[# + 1]]] &]");
        /* the factor list multiplies back to G up to a constant */
        CHK("FreeQ[Together[tB[\"G\"]/Times @@ (#[[1]]^#[[2]] & /@ tB[\"factors\"])], x] && AllTrue[tB[\"factors\"], 0 < #[[2]] < tB[\"m\"] &]");
#undef CHK
    }
}

/* ---- Arithmetic.m: constructors and linear operations -------------------- */
static void test_linear(void) {
    TR_TRUE("afZero[tBs[[7]]] === {0, 0, 0, 0, 0, 0}");
    TR_TRUE("afOne[tBs[[3]]] === {1, 0, 0}");
    TR_TRUE("afUnit[tBs[[3]], 2] === {0, 0, 1}");
    TR_TRUE("afUnit[tBs[[1]], 0] === afOne[tBs[[1]]]");
    TR_TRUE("afCanon[tQ, {x/x, (x^2 - 1)/(x - 1)}] === {1, 1 + x}");
    TR_TRUE("afCanon[tN2, {t^2, 1/(x - t)}] === {2, (t + x)/(-2 + x^2)}");
    TR_TRUE("afCanon[tG5, {5 x, 6 x + 1/2}] === {0, 3 + x}");

    TR_TRUE("afZeroQ[tQ, {0, 0, 0}]");
    TR_TRUE("afZeroQ[tQ, {x/x - 1, 0}]");
    TR_FALSE("afZeroQ[tQ, {0, 1/x}]");
    TR_TRUE("afZeroQ[tN2, {t^2 - 2, (t^2 - 2)/x}]");
    TR_TRUE("afZeroQ[tG5, {5, 10 x}]");
    TR_FALSE("afZeroQ[tG5, {5, 11 x}]");

    TR_TRUE("afPlus[tQ, {1, x}, {1/x, -x}] === {(1 + x)/x, 0}");
    TR_TRUE("afMinus[tQ, {1, x}, {1, x}] === {0, 0}");
    TR_TRUE("afPlus[tG5, {3, 4 x}, {2, x}] === {0, 0}");
    TR_TRUE("afPlus[tN2, {t, 0}, {t, 1}] === {2 t, 1}");
    TR_TRUE("afScale[tQ, 1/x, {x, x^2, 0}] === {1, x, 0}");
    TR_TRUE("afScale[tN2, t, {t, x}] === {2, t x}");
    TR_TRUE("afScale[tG7, 1/2, {2, 4 x}] === {1, 2 x}");
    TR_TRUE("Module[{ta = tElt[tBs[[5]], 1], tb = tElt[tBs[[5]], 2]},"
            " tSame[tQ, afMinus[tQ, afPlus[tQ, ta, tb], tb], ta] && tSame[tQ, afPlus[tQ, ta, tb], afPlus[tQ, tb, ta]]]");
    TR_TRUE("Module[{ta = tElt[tBs[[4]], 3]}, tSame[tQ, afScale[tQ, 2, ta], afPlus[tQ, ta, ta]]]");

    /* afFromY / afToY */
    TR_TRUE("afFromY[tQ, tBs[[7]], tYv^3/x + tYv, tYv] === {0, 1, 0, 1 + x, 0, 0}");
    TR_TRUE("Together[afToY[tBs[[7]], {1, 0, 0, 1, 0, 0}, tYv] - (1 + tYv^3/(x + x^2))] === 0");
    TR_TRUE("afFromY[tQ, tBs[[1]], 3, tYv] === {3, 0}");
    TR_TRUE("afFromY[tQ, tBs[[1]], (tYv + 1)/(x + 1), tYv] === {1/(1 + x), 1/(1 + x)}");
    TR_TRUE("afFromY[tN2, tBs[[1]], t^3 tYv, tYv] === {0, 2 t}");
    TR_TRUE("catchTrager[afFromY[tQ, tBs[[1]], 1/tYv, tYv]][[2]] === \"afFromYDenominatorInY\"");
    /* round trips on every curve */
    TR_TRUE("AllTrue[tBs, Function[tB, tSame[tQ, afFromY[tQ, tB, afToY[tB, tElt[tB, 2], tYv], tYv], tElt[tB, 2]]]]");
    TR_TRUE("AllTrue[tBs, Function[tB, Together[afToY[tB, afFromY[tQ, tB, Sum[(k + x) tYv^k, {k, 0, tB[\"m\"] - 1}], tYv], tYv]"
            " - Sum[(k + x) tYv^k, {k, 0, tB[\"m\"] - 1}]] === 0]]");
    /* Y = w_1 (d_1 = 1) */
    TR_TRUE("AllTrue[tBs, afFromY[tQ, #, tYv, tYv] === afUnit[#, 1] &]");
}

/* ---- Arithmetic.m: multiplication, power, norm, inverse --------------- */
static void test_ring(void) {
    char buf[2048];
    /* hand values on y^2 = x^3 + 1 */
    TR_TRUE("afTimes[tQ, tBs[[1]], {0, 1}, {0, 1}] === {1 + x^3, 0}");
    /* (1 + y)(1 - y) = 1 - (x^3 + 1) = -x^3 */
    TR_TRUE("afTimes[tQ, tBs[[1]], {1, 1}, {1, -1}] === {-x^3, 0}");
    TR_TRUE("afNorm[tQ, tBs[[1]], {1, 1}] === -x^3");
    TR_TRUE("afInverse[tBs[[1]], {1, 1}] === {-1/x^3, 1/x^3}");
    /* y^3 = x^2 (x + 1): w1 * w2 = Y^3/x = x (x + 1) */
    TR_TRUE("Expand[afTimes[tQ, tBs[[3]], {0, 1, 0}, {0, 0, 1}]] === {x + x^2, 0, 0}");
    TR_TRUE("afPower[tQ, tBs[[3]], {0, 1, 0}, 2] === {0, 0, x}");
    TR_TRUE("afPower[tQ, tBs[[3]], {0, 1, 0}, 0] === {1, 0, 0}");
    TR_TRUE("Expand[afPower[tQ, tBs[[3]], {0, 1, 0}, 3]] === {x^2 + x^3, 0, 0}");
    /* over Q(Sqrt[2]): (t + y)(t - y) = 2 - (x^3 + 1) = 1 - x^3 */
    TR_TRUE("afTimes[tN2, tBs[[1]], {t, 1}, {t, -1}] === {1 - x^3, 0}");
    TR_TRUE("afTimes[tN2, tBs[[1]], {t, 1}, {-t, 1}] === {-1 + x^3, 0}");
    TR_TRUE("afNorm[tN2, tBs[[1]], {t, 1}] === 1 - x^3");
    /* over GF(5): y * y = x^3 + 1 */
    TR_TRUE("afTimes[tG5, basisMod[tBs[[1]], 5], {0, 3}, {0, 2}] === {1 + x^3, 0}");

    for (int c = 1; c <= NCURVES; c++) {
#define CHK(fmt) do { snprintf(buf, sizeof buf, "With[{tB = tBs[[%d]], tm = tBs[[%d]][\"m\"]}, " fmt "]", c, c); TR_TRUE(buf); } while (0)
        /* y * y^(m-1) = G */
        CHK("tSame[tQ, afTimes[tQ, tB, afFromY[tQ, tB, tYv, tYv], afFromY[tQ, tB, tYv^(tm - 1), tYv]], afFromY[tQ, tB, tB[\"G\"], tYv]]");
        /* afPower[Y, m] = G, afPower[Y, m + 1] = G Y */
        CHK("tSame[tQ, afPower[tQ, tB, afUnit[tB, 1], tm], afScale[tQ, tB[\"G\"], afOne[tB]]] &&"
            " tSame[tQ, afPower[tQ, tB, afUnit[tB, 1], tm + 1], afScale[tQ, tB[\"G\"], afUnit[tB, 1]]]");
        /* afTimes agrees with multiplying the Y-expressions modulo Y^m - G */
        CHK("Module[{ta = tElt[tB, 1], tb = tElt[tB, 2]}, Together[afToY[tB, afTimes[tQ, tB, ta, tb], tYv]"
            " - tRedY[tB, afToY[tB, ta, tYv] afToY[tB, tb, tYv]]] === 0]");
        /* commutativity, associativity, distributivity, identity */
        CHK("Module[{ta = tElt[tB, 1], tb = tElt[tB, 2]}, tSame[tQ, afTimes[tQ, tB, ta, tb], afTimes[tQ, tB, tb, ta]]]");
        CHK("Module[{ta = tElt[tB, 1], tb = tElt[tB, 2], tc = tPolyElt[tB, 3]},"
            " tSame[tQ, afTimes[tQ, tB, afTimes[tQ, tB, ta, tb], tc], afTimes[tQ, tB, ta, afTimes[tQ, tB, tb, tc]]]]");
        CHK("Module[{ta = tElt[tB, 1], tb = tElt[tB, 2], tc = tPolyElt[tB, 3]},"
            " tSame[tQ, afTimes[tQ, tB, ta, afPlus[tQ, tb, tc]], afPlus[tQ, afTimes[tQ, tB, ta, tb], afTimes[tQ, tB, ta, tc]]]]");
        CHK("Module[{ta = tElt[tB, 4]}, tSame[tQ, afTimes[tQ, tB, ta, afOne[tB]], ta] && afZeroQ[tQ, afTimes[tQ, tB, ta, afZero[tB]]]]");
        /* afPower against repeated afTimes */
        CHK("Module[{ta = tPolyElt[tB, 1]}, tSame[tQ, afPower[tQ, tB, ta, 5],"
            " afTimes[tQ, tB, ta, afTimes[tQ, tB, ta, afTimes[tQ, tB, ta, afTimes[tQ, tB, ta, ta]]]]]]");
        /* the product of integral elements is integral (polynomial coordinates) */
        CHK("afPolynomialQ[tQ, afTimes[tQ, tB, tPolyElt[tB, 1], tPolyElt[tB, 2]]]");
        /* multiplication matrix: row k is w_k a; its determinant is the norm */
        CHK("Module[{ta = tElt[tB, 1], tM}, tM = afMultMatrix[tQ, tB, ta]; Length[tM] === tm &&"
            " AllTrue[Range[0, tm - 1], tSame[tQ, tM[[# + 1]], afTimes[tQ, tB, afUnit[tB, #], ta]] &] &&"
            " Together[Det[tM] - afNorm[tQ, tB, ta]] === 0]");
        /* N(Y) = (-1)^(m-1) G, N(r) = r^m for r in Q(x) */
        CHK("Together[afNorm[tQ, tB, afUnit[tB, 1]] - (-1)^(tm - 1) tB[\"G\"]] === 0");
        CHK("Together[afNorm[tQ, tB, afScale[tQ, (x + 2)/(x - 3), afOne[tB]]] - ((x + 2)/(x - 3))^tm] === 0");
        /* the norm is multiplicative */
        CHK("Module[{ta = tElt[tB, 1], tb = tPolyElt[tB, 2]}, Together[afNorm[tQ, tB, afTimes[tQ, tB, ta, tb]]"
            " - afNorm[tQ, tB, ta] afNorm[tQ, tB, tb]] === 0]");
        /* a * a^-1 = 1.  With rational-function coordinates afInverse is
           only exercised for m <= 3: at m = 4 its LinearSolve over Q(x) does
           not finish in 100 s (see the report); polynomial coordinates are
           tested on every curve. */
        if (c <= 4)
            CHK("Module[{ta = tElt[tB, 1]}, tSame[tQ, afTimes[tQ, tB, ta, afInverse[tB, ta]], afOne[tB]]]");
        CHK("Module[{ta = tPolyElt[tB, 2]}, tSame[tQ, afTimes[tQ, tB, afInverse[tB, ta], ta], afOne[tB]]]");
        /* 1/Y = Y^(m-1)/G */
        CHK("tSame[tQ, afInverse[tB, afUnit[tB, 1]], afFromY[tQ, tB, tYv^(tm - 1)/tB[\"G\"], tYv]]");
        /* the same laws over Q(Sqrt[2]) */
        CHK("Module[{ta = tElt[tB, 1] + t, tb = t tPolyElt[tB, 2]}, afTimes[tN2, tB, ta, tb] === afTimes[tN2, tB, tb, ta] &&"
            " Together[afNorm[tN2, tB, afTimes[tN2, tB, ta, tb]] - rc[tN2, afNorm[tN2, tB, ta] afNorm[tN2, tB, tb]]] === 0]");
#undef CHK
    }
    /* the inverse of zero is an internal error */
    TR_TRUE("catchTrager[afInverse[tBs[[3]], {0, 0, 0}]][[2]] === \"afInverseZeroDivisor\"");
}

/* ---- Arithmetic.m: derivation --------------------------------------------- */
static void test_derivation(void) {
    char buf[2048];
    /* y' = 3x^2/(2(x^3 + 1)) y on y^2 = x^3 + 1 */
    TR_TRUE("afD[tQ, tBs[[1]], {0, 1}] === {0, (3 x^2)/(2 (1 + x^3))}");
    TR_TRUE("afD[tQ, tBs[[1]], {x^2, 0}] === {2 x, 0}");
    TR_TRUE("afZeroQ[tQ, afD[tQ, tBs[[5]], {7, 0, 0, 0}]]");
    /* over GF(5) with the reduced basis: d/dx (x^5 + x) = 1 */
    TR_TRUE("afD[tG5, basisMod[tBs[[1]], 5], {x^5 + x, 0}] === {1, 0}");
    for (int c = 1; c <= NCURVES; c++) {
#define CHK(fmt) do { snprintf(buf, sizeof buf, "With[{tB = tBs[[%d]], tm = tBs[[%d]][\"m\"]}, " fmt "]", c, c); TR_TRUE(buf); } while (0)
        /* afD agrees with the derivative of afToY under Y' = G' Y/(m G) */
        CHK("Module[{ta = tElt[tB, 1]}, tSame[tQ, afD[tQ, tB, ta], afFromY[tQ, tB, Together[tDY[tB, afToY[tB, ta, tYv]]], tYv]]]");
        /* the Leibniz rule */
        CHK("Module[{ta = tElt[tB, 1], tb = tPolyElt[tB, 2]}, tSame[tQ, afD[tQ, tB, afTimes[tQ, tB, ta, tb]],"
            " afPlus[tQ, afTimes[tQ, tB, afD[tQ, tB, ta], tb], afTimes[tQ, tB, ta, afD[tQ, tB, tb]]]]]");
        /* d/dx Y^m = d/dx G */
        CHK("tSame[tQ, afD[tQ, tB, afPower[tQ, tB, afUnit[tB, 1], tm]], afScale[tQ, D[tB[\"G\"], x], afOne[tB]]]");
        /* linearity */
        CHK("Module[{ta = tElt[tB, 1], tb = tElt[tB, 3]}, tSame[tQ, afD[tQ, tB, afPlus[tQ, ta, afScale[tQ, 3, tb]]],"
            " afPlus[tQ, afD[tQ, tB, ta], afScale[tQ, 3, afD[tQ, tB, tb]]]]]");
        /* (1/a)' = -a'/a^2 */
        CHK("Module[{ta = tPolyElt[tB, 1], ti}, ti = afInverse[tB, ta]; tSame[tQ, afD[tQ, tB, ti],"
            " afScale[tQ, -1, afTimes[tQ, tB, afD[tQ, tB, ta], afTimes[tQ, tB, ti, ti]]]]]");
        /* Leibniz over Q(Sqrt[2]) */
        CHK("Module[{ta = t tElt[tB, 1], tb = tPolyElt[tB, 2] + t}, tSame[tN2, afD[tN2, tB, afTimes[tN2, tB, ta, tb]],"
            " afPlus[tN2, afTimes[tN2, tB, afD[tN2, tB, ta], tb], afTimes[tN2, tB, ta, afD[tN2, tB, tb]]]]]");
#undef CHK
    }
}

/* ---- Arithmetic.m: integrality, reduction mod p --------------------------- */
static void test_integrality_mod(void) {
    char buf[2048];
    TR_TRUE("afPolynomialQ[tQ, {x, 1 + x^2}]");
    TR_FALSE("afPolynomialQ[tQ, {x, 1/x}]");
    TR_TRUE("afPolynomialQ[tQ, {0, 0}]");
    TR_TRUE("afPolynomialQ[tQ, {x^2/(x), 3/2}]");
    TR_TRUE("afPolynomialQ[tN2, {t x, (t^2 - 2)/x}]");
    TR_FALSE("afPolynomialQ[tN2, {1/(x - t), 0}]");
    TR_TRUE("afPolynomialQ[tG5, {5/x + x, 1/2}]");
    TR_FALSE("afPolynomialQ[tG5, {1/(x + 5), 0}]");
    /* Y^i/d_i is integral, Y^i/(p d_i) is not (p a prime factor of G) */
    TR_TRUE("AllTrue[tBs, Function[tB, AllTrue[Range[0, tB[\"m\"] - 1],"
            " afPolynomialQ[tQ, afFromY[tQ, tB, tYv^#/tB[\"d\"][[# + 1]], tYv]] &&"
            " !afPolynomialQ[tQ, afFromY[tQ, tB, tYv^#/(tB[\"d\"][[# + 1]] tB[\"factors\"][[1, 1]]), tYv]] &]]]");

    /* basisMod: hand values for y^2 = x^3 + 1 mod 3, where x^3 + 1 = (x + 1)^3 */
    TR_TRUE("Module[{tb = basisMod[tBs[[1]], 3]}, tb[\"p\"] === 3 && tb[\"G\"] === 1 + x^3 &&"
            " tb[\"factors\"] === {{1 + x, 3}} && tb[\"s\"] === {0, 0} && tb[\"m\"] === 2]");
    /* y^3 = 2x^2(x - 1) mod 5: G = 2x^3 - 2x^2 -> 2x^3 + 3x^2 */
    TR_TRUE("basisMod[tBs[[4]], 5][\"G\"] === 3 x^2 + 2 x^3");
    for (int c = 1; c <= NCURVES; c++) {
        const char* ps[] = {"7", "11"};
        for (int k = 0; k < 2; k++) {
#define CHK(fmt) do { snprintf(buf, sizeof buf, "With[{tB = tBs[[%d]], tp = %s}, With[{tBp = basisMod[tB, tp], tF = fieldGF[tp]}, " fmt "]]", c, ps[k]); TR_TRUE(buf); } while (0)
            /* every component of the descriptor is reduced */
            CHK("tBp[\"G\"] === gfMod[tB[\"G\"], tp] && tBp[\"d\"] === (gfMod[#, tp] & /@ tB[\"d\"]) &&"
                " tBp[\"mult\"] === Map[{#[[1]], gfMod[#[[2]], tp]} &, tB[\"mult\"], {2}] && tBp[\"delta\"] === tB[\"delta\"]");
            /* s_i reduced as rational functions */
            CHK("AllTrue[Range[tB[\"m\"]], rcZeroQ[tF, tBp[\"s\"][[#]] - tB[\"s\"][[#]]] &]");
            /* the reduced factors multiply back to G mod p */
            CHK("Module[{tpr = Times @@ (#[[1]]^#[[2]] & /@ tBp[\"factors\"])}, kZeroQ[tF, tBp[\"G\"] - kLC[tF, tBp[\"G\"], x] tpr]]");
            /* products reduce consistently: (a b mod p) = (a mod p)(b mod p) */
            CHK("Module[{ta = tPolyElt[tB, 1], tb = tPolyElt[tB, 2]}, afTimes[tF, tBp, ta, tb] === afCanon[tF, afTimes[tQ, tB, ta, tb]]]");
            CHK("Module[{ta = tElt[tB, 1], tb = tPolyElt[tB, 2]}, afTimes[tF, tBp, ta, tb] === afCanon[tF, afTimes[tQ, tB, ta, tb]]]");
            /* derivatives reduce consistently (p does not divide m) */
            CHK("Module[{ta = tElt[tB, 2]}, afD[tF, tBp, ta] === afCanon[tF, afD[tQ, tB, ta]]]");
            /* the norm reduces consistently */
            CHK("Module[{ta = tPolyElt[tB, 1]}, afNorm[tF, tBp, ta] === rc[tF, afNorm[tQ, tB, ta]]]");
            /* y * y^(m-1) = G over GF(p) */
            CHK("afTimes[tF, tBp, afUnit[tBp, 1], afFromY[tF, tBp, tYv^(tBp[\"m\"] - 1), tYv]] === afFromY[tF, tBp, tBp[\"G\"], tYv]");
#undef CHK
        }
    }
}

int main(void) {
    tr_init();
    setup_basis();
    test_basis_values();
    test_basis_properties();
    test_linear();
    test_ring();
    test_derivation();
    test_integrality_mod();
    printf("trager_basis_arith: %d assertions passed\n", tr_assertions);
    return 0;
}
