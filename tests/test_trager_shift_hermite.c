/* test_trager_shift_hermite.c
 *
 * Unit tests for src/internal/trager/Shift.m and Hermite.m -- the behaviour
 * of a differential at x = Infinity, the change of variable x = a + 1/z that
 * moves infinity to a finite point ([T] Ch 4 Sec 1, Sec 3 pp.48-50; [B] Sec
 * 2.4), and algebraic Hermite reduction ([T] Ch 4 Sec 2 pp.44-48; [B] Sec
 * 2.1-2.2, Example 3 p.15).
 *
 * Coverage, function by function:
 *   infScaledOrder      hand-computed lam(c, j) = (m/r) ord_P(c w_j) on curves
 *                       with m = 2..5, unramified / ramified infinity, d_j != 1
 *   regularAtInfinityQ  c Y^j dx at the threshold m + gcd(m, M), derived from
 *                       ord_P(x) = -r, ord_P(Y) = -r M/m, ord_P(dx) = -r - 1
 *   integralAtInfinityQ the same with threshold 0
 *   shiftCurve          the chosen a avoids the roots of G and of the
 *                       denominator; k = Ceiling[M/m]; Gz and cz satisfy the
 *                       substitution identities x = a + 1/z, Y = Y1 z^-k,
 *                       dx = -dz/z^2; Gz is reduced (reducedCurveQ) and of the
 *                       SAME genus (curveGenus); a 30-digit numeric check of
 *                       c(x) G(x)^(1/m) dx/dz = cz(z) Gz(z)^(1/m) at a real
 *                       point where both principal roots are real; the
 *                       unbounded shift-point search
 *   unshift             exact round trip against the forward substitution,
 *                       and the differential identity read backwards
 *   hermiteReduce       the invariant afD[v] + r == f exactly, squarefree
 *                       remainder denominator; Bronstein Example 3; hand
 *                       antiderivatives (x^2+1)^(-3/2), x^(-3/2), x^(-4/3),
 *                       x^(-7/3)(x+1)^(1/3); branch-point V sharing a factor
 *                       with G; m = 2..5; repeated factors in G; f = v0'
 *                       (v - v0 must be integral)
 *   commonDenominator   LCM of coordinate denominators
 *   hermiteSolve        B C == A (mod V), deg B < deg V; both internal throws
 *   squarefreeDenominatorQ
 *
 * Expected values are derived by hand (see the comments) or are invariants
 * computed with Mathilda builtins independently of the code under test.
 * Assertions are hard failures (trager_test_utils.h).
 */
#include "trager_test_utils.h"

/* Evaluate a statement for its side effects. */
static void tr_do(const char* in) { free(tr_eval(in)); }

/* Helpers in Global` (t-prefixed so they never collide with package names). */
static void define_helpers(void) {
    char* s = tr_eval(
        "tQ = fieldQ[];"
        /* Hermite invariant: afD[v] + r == f exactly and r squarefree. */
        "tHermOK[tm_, tG_, tf_] := Module[{tb = curveBasis[tm, tG, x], th},"
        "  th = hermiteReduce[tb, tf];"
        "  afZeroQ[tQ, afMinus[tQ, afPlus[tQ, afD[tQ, tb, th[\"Algebraic\"]], th[\"Remainder\"]], tf]] &&"
        "  squarefreeDenominatorQ[th[\"Remainder\"], x]];"
        "tHerm[tm_, tG_, tf_] := hermiteReduce[curveBasis[tm, tG, x], tf];"
        /* f = v0': the algebraic part differs from v0 by an integral element. */
        "tHermDeriv[tm_, tG_, tv0_] := Module[{tb = curveBasis[tm, tG, x], th, tf},"
        "  tf = afD[tQ, tb, tv0]; th = hermiteReduce[tb, tf];"
        "  afPolynomialQ[tQ, afMinus[tQ, tv0, th[\"Algebraic\"]]] &&"
        "  afZeroQ[tQ, afMinus[tQ, afPlus[tQ, afD[tQ, tb, th[\"Algebraic\"]], th[\"Remainder\"]], tf]] &&"
        "  squarefreeDenominatorQ[th[\"Remainder\"], x]];"
        /* shiftCurve: every property of the contract. */
        "tShiftOK[tc_, tm_, tG_] := Module[{ts = shiftCurve[tc, tm, tG, x, z], ta, tk, tM = Exponent[tG, x], tz0 = 1/3, tx0, tnum},"
        "  ta = ts[\"a\"]; tk = ts[\"k\"]; tx0 = ta + 3;"
        "  tnum = N[(tc /. x -> tx0) (tG /. x -> tx0)^(1/tm) (-1/tz0^2) - (ts[\"c\"] /. z -> tz0) (ts[\"G\"] /. z -> tz0)^(1/tm), 30];"
        "  {(tG /. x -> ta) =!= 0, (Denominator[Together[tc]] /. x -> ta) =!= 0,"
        "   tk === Ceiling[tM/tm],"
        "   Together[ts[\"G\"] - z^(tk tm) (tG /. x -> ta + 1/z)] === 0,"
        "   PolynomialQ[ts[\"G\"], z],"
        "   Together[(tc /. x -> ta + 1/z) z^(-tk) (-1/z^2) - ts[\"c\"]] === 0,"
        "   reducedCurveQ[tm, ts[\"G\"], z],"
        "   curveGenus[tm, ts[\"G\"], z] === curveGenus[tm, tG, x],"
        "   Abs[tnum] < 10^-20}];"
        /* Earlier candidates in the search order must all be excluded. */
        "tShiftFirstQ[tc_, tm_, tG_] := Module[{ta = shiftCurve[tc, tm, tG, x, z][\"a\"], tc0, tpos},"
        "  tc0 = Join[{0}, Flatten[Table[{i, -i}, {i, 64}]]]; tpos = First[Flatten[Position[tc0, ta]]];"
        "  AllTrue[Take[tc0, tpos - 1], (tG /. x -> #) === 0 || (Denominator[Together[tc]] /. x -> #) === 0 &]];"
        /* unshift round trip */
        "tRoundTrip[ts_, te_] := Together[unshift[te /. {x -> ts[\"a\"] + 1/z, Y -> Y z^(-ts[\"k\"])}, ts, x, z, Y] - te] === 0;"
        "True");
    free(s);
}

/* ---------------------------------------------------------------- Shift.m */

/* infScaledOrder: lam(c, j) = -m deg c + m deg d_j - j M. */
static void test_infScaledOrder(void) {
    /* Y^2 = x^2+1 (M = 2, d = {1, 1}). */
    tr_do("tB2 = curveBasis[2, x^2 + 1, x]; True");
    TR_EQ("infScaledOrder[tQ, tB2, 1, 1]", "-2");             /* Y ~ x */
    TR_EQ("infScaledOrder[tQ, tB2, 1/(x^2 + 1), 1]", "2");
    TR_EQ("infScaledOrder[tQ, tB2, x, 0]", "-2");
    TR_EQ("infScaledOrder[tQ, tB2, 1/x^3, 0]", "6");
    TR_EQ("infScaledOrder[tQ, tB2, 0, 1]", "Infinity");
    TR_EQ("infScaledOrder[tQ, tB2, 0, 0]", "Infinity");
    TR_EQ("infScaledOrder[tQ, tB2, (x + 1)/(x^3 + 2), 1]", "2");
    /* Y^3 = x^2 (x+1): d = {1, 1, x}, M = 3, r = 1.  w_2 = Y^2/x ~ x. */
    tr_do("tB3 = curveBasis[3, x^2 (x + 1), x]; True");
    TR_EQ("tB3[\"d\"]", "{1, 1, x}");
    TR_EQ("infScaledOrder[tQ, tB3, 1, 2]", "-3");
    TR_EQ("infScaledOrder[tQ, tB3, 1, 1]", "-3");
    TR_EQ("infScaledOrder[tQ, tB3, 1/x^2, 2]", "3");
    TR_EQ("infScaledOrder[tQ, tB3, x/(x + 1), 1]", "-3");
    /* Y^3 = x^2+1: M = 2, r = 3; ord_P(Y) = -2. */
    tr_do("tB32 = curveBasis[3, x^2 + 1, x]; True");
    TR_EQ("infScaledOrder[tQ, tB32, 1, 1]", "-2");
    TR_EQ("infScaledOrder[tQ, tB32, 1, 2]", "-4");
    TR_EQ("infScaledOrder[tQ, tB32, 1/x, 1]", "1");
    /* Y^5 = x^3 + 2: M = 3, r = 5. */
    tr_do("tB5 = curveBasis[5, x^3 + 2, x]; True");
    TR_EQ("infScaledOrder[tQ, tB5, 1/x^3, 2]", "9");
    TR_EQ("infScaledOrder[tQ, tB5, 1/x^2, 2]", "4");
    TR_EQ("infScaledOrder[tQ, tB5, x, 4]", "-17");
}

/* regularAtInfinityQ: threshold m + gcd(m, M). */
static void test_regularAtInfinityQ(void) {
    /* 1/Sqrt[x^2+1] dx has a log pole at infinity; Y/x^3 dx ~ dx/x^2 does not. */
    TR_FALSE("regularAtInfinityQ[tQ, tB2, {0, 1/(x^2 + 1)}]");
    TR_TRUE("regularAtInfinityQ[tQ, tB2, {0, 1/(x^3 + 1)}]");
    TR_TRUE("regularAtInfinityQ[tQ, tB2, {1/x^2, 0}]");
    TR_FALSE("regularAtInfinityQ[tQ, tB2, {1/x, 0}]");                  /* dx/x */
    TR_FALSE("regularAtInfinityQ[tQ, tB2, {1/x^2, 1/x^2}]");            /* one bad coordinate */
    TR_TRUE("regularAtInfinityQ[tQ, tB2, {0, 0}]");
    /* Y^2 = 1 - x^3 (Bronstein Ex. 6): Y/(x - x^4) dx is regular (thr 3, lam 5). */
    tr_do("tB6 = curveBasis[2, 1 - x^3, x]; True");
    TR_TRUE("regularAtInfinityQ[tQ, tB6, {0, 1/(x - x^4)}]");
    /* Ex. 5: x Y/(1 - x^3) dx: lam = 4 - 3 = 1 < 3 -> pole at infinity. */
    TR_FALSE("regularAtInfinityQ[tQ, tB6, {0, x/(1 - x^3)}]");
    /* Y^3 = x^2+1: ord(Y x^-k dx) = -2 + 3k - 4 >= 0 iff k >= 2. */
    TR_FALSE("regularAtInfinityQ[tQ, tB32, {0, 1/x, 0}]");
    TR_TRUE("regularAtInfinityQ[tQ, tB32, {0, 1/x^2, 0}]");
    /* Y^4 = x^2+1: r = 2, ord(Y x^-k dx) = -1 + 2k - 3 >= 0 iff k >= 2. */
    tr_do("tB42 = curveBasis[4, x^2 + 1, x]; True");
    TR_FALSE("regularAtInfinityQ[tQ, tB42, {0, 1/x, 0, 0}]");
    TR_TRUE("regularAtInfinityQ[tQ, tB42, {0, 1/x^2, 0, 0}]");
    /* Y^5 = x^3+2: Y^2 c dx regular iff deg c <= -3. */
    TR_FALSE("regularAtInfinityQ[tQ, tB5, {0, 0, 1/x^2, 0, 0}]");
    TR_TRUE("regularAtInfinityQ[tQ, tB5, {0, 0, 1/x^3, 0, 0}]");
}

static void test_integralAtInfinityQ(void) {
    TR_TRUE("integralAtInfinityQ[tQ, tB2, {0, 1/x}]");     /* Y/x -> +-1 */
    TR_FALSE("integralAtInfinityQ[tQ, tB2, {0, 1}]");
    TR_FALSE("integralAtInfinityQ[tQ, tB2, {x, 0}]");
    TR_TRUE("integralAtInfinityQ[tQ, tB2, {1, 0}]");
    TR_TRUE("integralAtInfinityQ[tQ, tB2, {1/x, 1/x^2}]");
    TR_TRUE("integralAtInfinityQ[tQ, tB2, {0, 0}]");
    /* Y^3 = x^2+1: Y/x has order 1; Y^2/x ~ x^(1/3) has a pole; Y^2/x^2 not. */
    TR_TRUE("integralAtInfinityQ[tQ, tB32, {0, 1/x, 0}]");
    TR_FALSE("integralAtInfinityQ[tQ, tB32, {0, 0, 1/x}]");
    TR_TRUE("integralAtInfinityQ[tQ, tB32, {0, 0, 1/x^2}]");
    /* Y^3 = x^2(x+1): w_2 = Y^2/x ~ x has a pole, w_2/x does not. */
    TR_FALSE("integralAtInfinityQ[tQ, tB3, {0, 0, 1}]");
    TR_TRUE("integralAtInfinityQ[tQ, tB3, {0, 0, 1/x}]");
    /* over Q(Sqrt[2]): only the degree of each coordinate matters */
    TR_TRUE("integralAtInfinityQ[fieldNF[t, t^2 - 2, Sqrt[2]], tB2, {t, t/x}]");
    TR_FALSE("integralAtInfinityQ[fieldNF[t, t^2 - 2, Sqrt[2]], tB2, {t x + 1, 0}]");
    TR_TRUE("regularAtInfinityQ[fieldNF[t, t^2 - 2, Sqrt[2]], tB2, {0, t/(x^3 + t)}]");
}

static void test_shiftCurve(void) {
    /* 1/Sqrt[x^2+1]: a = 0, k = 1, Gz = z^2 + 1, cz = -1/(z (1+z^2)) by hand. */
    tr_do("tS1 = shiftCurve[1/(x^2 + 1), 2, x^2 + 1, x, z]; True");
    TR_EQ("tS1[\"a\"]", "0");
    TR_EQ("tS1[\"k\"]", "1");
    TR_TRUE("Expand[tS1[\"G\"] - (z^2 + 1)] === 0");
    TR_TRUE("Together[tS1[\"c\"] + 1/(z (1 + z^2))] === 0");
    TR_EQ("tShiftOK[1/(x^2 + 1), 2, x^2 + 1]", "{True, True, True, True, True, True, True, True, True}");
    /* G(0) = 0 -> a = 1; Gz = z^3 G(1 + 1/z) = 2z^3 + 4z^2 + 3z + 1 (k = 1). */
    tr_do("tS2 = shiftCurve[1/x, 3, x^3 + x, x, z]; True");
    TR_EQ("tS2[\"a\"]", "1");
    TR_TRUE("Expand[tS2[\"G\"] - (2 z^3 + 4 z^2 + 3 z + 1)] === 0");
    /* Roots of G at 0, +-1 and of the denominator at 2 force a = -2;
     * k = Ceiling[3/2] = 2. */
    tr_do("tS3 = shiftCurve[1/(x - 2), 2, x^3 - x, x, z]; True");
    TR_EQ("tS3[\"a\"]", "-2");
    TR_EQ("tS3[\"k\"]", "2");
    TR_TRUE("tShiftFirstQ[1/(x - 2), 2, x^3 - x]");
    TR_TRUE("tShiftFirstQ[1/((x - 1) (x + 1) x (x - 2)), 3, x^2 + 2]");
    TR_EQ("shiftCurve[1/((x - 1) (x + 1) x (x - 2)), 3, x^2 + 2, x, z][\"a\"]", "-2");
    TR_EQ("shiftCurve[1, 2, x (x + 1) (x - 1) (x + 2) (x - 2) + 7, x, z][\"a\"]", "0");
    /* The full property list across m = 2..5, M < m, M = m, M > m, repeated
     * factors in G (Y^3 = x^2(x+1), Y^4 = (x^2+1)^2 (x-3)). */
    const char* cases[] = {
        "tShiftOK[1/x, 3, x^3 + x]",
        "tShiftOK[1/(x - 2), 2, x^3 - x]",
        "tShiftOK[1, 3, x^2 + 1]",
        "tShiftOK[x, 3, x^4 + 2]",
        "tShiftOK[1/x, 4, x^3 + x]",
        "tShiftOK[1/(x - 1)^2, 5, x^3 + 2]",
        "tShiftOK[x^2, 2, x^5 + 3]",
        "tShiftOK[1/x, 3, x^2 (x + 1)]",
        "tShiftOK[1, 4, (x^2 + 1)^2 (x - 3)]",
        "tShiftOK[(2 x^8 + 1)/(x^17 + 2 x^9 + x), 2, x^8 + 1]",
        "tShiftOK[x/(1 - x^3), 2, 1 - x^3]",
        "tShiftOK[(x^2 + 3)/((x + 1)^2 (x - 5)), 5, x^7 + x + 1]",
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++)
        TR_EQ(cases[i], "{True, True, True, True, True, True, True, True, True}");
    /* The z-factor of Gz is z^(k m - M), multiplicity < m. */
    TR_TRUE("Exponent[shiftCurve[x, 3, x^4 + 2, x, z][\"G\"], z, Min] === 2");
    TR_TRUE("Exponent[shiftCurve[1, 2, x^2 + 1, x, z][\"G\"], z, Min] === 0");
    /* No bound on the shift point: every |a| <= 64 is a pole of c, so a = 65. */
    TR_EQ("shiftCurve[1/Product[x - i, {i, -64, 64}], 2, x^2 + 1, x, z][\"a\"]", "65");
    TR_EQ("firstOrdinaryInteger[{x (x - 1) (x + 1) (x - 2)}, x]", "-2");
}

static void test_unshift(void) {
    TR_TRUE("Together[unshift[z Y, tS1, x, z, Y] - Y/x^2] === 0");  /* a = 0, k = 1 */
    TR_TRUE("tRoundTrip[tS1, x^3 Y + 1/(x - 7)]");
    TR_TRUE("tRoundTrip[tS2, (x^2 + Y^2)/(x - 1) + Y x]");
    TR_TRUE("tRoundTrip[tS3, Y^2 x/(x + 2)^3 - 5]");
    TR_TRUE("tRoundTrip[shiftCurve[1, 4, (x^2 + 1)^2 (x - 3), x, z], Y^3/(x^2 + 1) + Y]");
    /* Gz back in x: (x - a)^(-k m) G(x). */
    TR_TRUE("Together[unshift[tS3[\"G\"], tS3, x, z, Y] - (x + 2)^(-4) (x^3 - x)] === 0");
    /* The differential read backwards: cz(z) Y1 dz = c(x) Y dx with
     * dz/dx = -1/(x - a)^2. */
    TR_TRUE("Together[unshift[tS2[\"c\"] Y, tS2, x, z, Y] (-1/(x - 1)^2) - Y/x] === 0");
    TR_TRUE("Together[unshift[tS3[\"c\"] Y, tS3, x, z, Y] (-1/(x + 2)^2) - Y/(x - 2)] === 0");
    /* Y1^m - Gz maps to (x-a)^(-km) (Y^m - G): the curve is preserved. */
    TR_TRUE("Together[unshift[Y^2 - tS3[\"G\"], tS3, x, z, Y] - (Y^2 - (x^3 - x))/(x + 2)^4] === 0");
}

/* --------------------------------------------------------------- Hermite.m */

static void test_hermiteReduce_bronstein(void) {
    /* Bronstein Example 3 p.15:
     *   Int (2x^8+1) Sqrt[x^8+1]/(x^17+2x^9+x)
     *     = -Sqrt[x^8+1]/(4(x^8+1)) + Int Sqrt[x^8+1]/(x(x^8+1)). */
    tr_do("tH3 = tHerm[2, x^8 + 1, {0, (2 x^8 + 1)/(x^17 + 2 x^9 + x)}]; True");
    TR_TRUE("Together[tH3[\"Algebraic\"] - {0, -1/(4 (x^8 + 1))}] === {0, 0}");
    TR_TRUE("Together[tH3[\"Remainder\"] - {0, 1/(x (x^8 + 1))}] === {0, 0}");
    TR_TRUE("tHermOK[2, x^8 + 1, {0, (2 x^8 + 1)/(x^17 + 2 x^9 + x)}]");
    TR_TRUE("Sort[Keys[tH3]] === {\"Algebraic\", \"Remainder\"}");
}

static void test_hermiteReduce_hand(void) {
    /* (x^2+1)^(-3/2) = d/dx [x/Sqrt[x^2+1]]: v = x Y/(x^2+1), r = 0. */
    tr_do("tHa = tHerm[2, x^2 + 1, {0, 1/(x^2 + 1)^2}]; True");
    TR_TRUE("Together[tHa[\"Algebraic\"] - {0, x/(x^2 + 1)}] === {0, 0}");
    TR_EQ("tHa[\"Remainder\"]", "{0, 0}");
    /* Branch point (V = x divides G = x): x^(-3/2) = d/dx[-2 x^(-1/2)] = d/dx[-2Y/x]. */
    tr_do("tHb = tHerm[2, x, {0, 1/x^2}]; True");
    TR_TRUE("Together[tHb[\"Algebraic\"] - {0, -2/x}] === {0, 0}");
    TR_EQ("tHb[\"Remainder\"]", "{0, 0}");
    /* m = 3 branch point: x^(-4/3) = Y^2/x^2 on Y^3 = x; integral -3 x^(-1/3) = -3 Y^2/x. */
    tr_do("tHc = tHerm[3, x, {0, 0, 1/x^2}]; True");
    TR_TRUE("Together[tHc[\"Algebraic\"] - {0, 0, -3/x}] === {0, 0, 0}");
    TR_EQ("tHc[\"Remainder\"]", "{0, 0, 0}");
    /* Y^3 = x^2(x+1), Y/x^3 = x^(-7/3)(x+1)^(1/3); integral
     * -3/4 x^(-4/3)(x+1)^(4/3) = -3(x+1)/(4x^2) Y (differentiate to check). */
    tr_do("tHd = tHerm[3, x^2 (x + 1), {0, 1/x^3, 0}]; True");
    TR_TRUE("Together[tHd[\"Algebraic\"] - {0, -3 (x + 1)/(4 x^2), 0}] === {0, 0, 0}");
    TR_EQ("tHd[\"Remainder\"]", "{0, 0, 0}");
    /* Rational coordinate: 1/(x-1)^2 integrates to -1/(x-1). */
    tr_do("tHe = tHerm[2, x^3 + 1, {1/(x - 1)^2, 0}]; True");
    TR_TRUE("Together[tHe[\"Algebraic\"] - {-1/(x - 1), 0}] === {0, 0}");
    TR_EQ("tHe[\"Remainder\"]", "{0, 0}");
    /* Already squarefree: nothing to do. */
    tr_do("tHf = tHerm[2, x^8 + 1, {0, 1/(x (x^8 + 1))}]; True");
    TR_EQ("tHf[\"Algebraic\"]", "{0, 0}");
    TR_TRUE("Together[tHf[\"Remainder\"] - {0, 1/(x (x^8 + 1))}] === {0, 0}");
    TR_EQ("tHerm[2, x^8 + 1, {0, 0}][\"Remainder\"]", "{0, 0}");
    /* Polynomial coordinates have denominator 1: untouched. */
    TR_TRUE("tHerm[3, x^2 + 1, {x^2, 3, x}][\"Algebraic\"] === {0, 0, 0}");
}

static void test_hermiteReduce_invariants(void) {
    const char* cases[] = {
        "tHermOK[2, x^2 + 1, {0, 1/(x^2 + 1)^2}]",
        "tHermOK[2, x, {0, 1/x^2}]",
        "tHermOK[3, x, {0, 0, 1/x^2}]",
        "tHermOK[3, x^2 (x + 1), {0, 1/x^3, 0}]",
        /* m = 4: V = x - 1 (ordinary) and x^2 + 1 (ordinary), multiplicities 3 and 2. */
        "tHermOK[4, x^3 + 1, {0, 1/(x - 1)^3, 0, x/(x^2 + 1)^2}]",
        /* m = 5, mixed multiplicities, w_4 polynomial. */
        "tHermOK[5, x^2 + x + 1, {0, 0, 1/(x^5 (x - 2)^2), 0, 1}]",
        /* both coordinates, V = x + 1 a branch factor of G = x^3 + 1 */
        "tHermOK[2, x^3 + 1, {1/(x - 1)^2, x/(x + 1)^3}]",
        /* branch-point V: x^2 + 1 and x - 3 divide G */
        "tHermOK[2, (x^2 + 1) (x - 3), {0, 1/((x - 3)^3 (x^2 + 1)^2)}]",
        /* repeated factor in G, m = 4: G = (x^2+1)^2 (x-3), d = {1, 1, x^2+1, x^2+1} */
        "tHermOK[4, (x^2 + 1)^2 (x - 3), {0, 1/(x^2 + 1)^2, 1/(x^2 + 1)^3, x/(x - 3)^2}]",
        /* m = 3, branch factor of multiplicity 2 in G */
        "tHermOK[3, x^2 (x + 1), {1/x^4, 1/(x^3 (x + 1)^2), (x + 5)/(x^2 (x - 1)^3)}]",
        /* high multiplicity */
        "tHermOK[2, x^2 + 1, {0, 1/(x^2 + 1)^5}]",
        "tHermOK[2, x^3 + x + 1, {0, (x^4 + 1)/((x - 1)^4 (x + 2)^2 x)}]",
        "tHermOK[3, x^3 + 2, {0, 1/(x - 1)^2, 1/(x^2 - 2)^3}]",
        "tHermOK[5, x^3 + 2, {0, 1/x^2, 0, 1/(x + 1)^3, 1/(x^3 + 2)^2}]",
        /* m = 5 with every coordinate polar at both branch factors of G */
        "tHermOK[5, x^2 (x - 1), {1/x^3, 1/x^4, 1/(x^2 (x - 1)^2), 1/(x - 1)^3, x/(x^2 (x - 1))^2}]",
        /* non-monic G and content in the denominators */
        "tHermOK[2, 2 x^3 + 3, {0, 1/(2 x^3 + 3)^2}]",
        "tHermOK[3, 5 (x^2 + 1)^2 (x + 4), {0, 1/(x^2 + 1)^3, 1/((x + 4)^2 (x^2 + 1)^4)}]",
        "tHermOK[4, x^3 (x + 1), {0, 1/x^5, 1/x^3, 1/(x^4 (x + 1)^2)}]",
        "tHermOK[2, x^2 + 1, {0, 1/(3 (x^2 + 1)^2)}]",
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) TR_TRUE(cases[i]);
    /* f = v0' exactly: the algebraic part equals v0 up to an integral element. */
    TR_TRUE("tHermDeriv[4, x^3 + 1, {0, 1/(x - 1)^2, 0, x/(x^2 + 1)}]");
    TR_TRUE("tHermDeriv[2, x^2 + 1, {1/x, x/(x^2 + 1)^2}]");
    TR_TRUE("tHermDeriv[3, x^2 (x + 1), {0, 1/x^2, (x + 2)/(x - 3)^3}]");
    TR_TRUE("tHermDeriv[2, x, {0, 1/x^3}]");
    TR_TRUE("tHermDeriv[5, x^4 + 1, {0, 0, 1/(x^4 + 1), 0, 1/(x + 3)^2}]");
    /* the m = 4 case: remainder exactly zero, v = v0 */
    TR_TRUE("Module[{tb = curveBasis[4, x^3 + 1, x], tv0 = {0, 1/(x - 1)^2, 0, x/(x^2 + 1)}, th},"
            " th = hermiteReduce[tb, afD[tQ, tb, tv0]];"
            " th[\"Remainder\"] === {0, 0, 0, 0} && Together[th[\"Algebraic\"] - tv0] === {0, 0, 0, 0}]");
}

static void test_commonDenominator(void) {
    TR_TRUE("Expand[commonDenominator[{1/x^2, 1/(x (x + 1)), 3}, x] - x^2 (x + 1)] === 0");
    TR_EQ("commonDenominator[{0, 0}, x]", "1");
    TR_EQ("commonDenominator[{x^2, 5}, x]", "1");
    TR_TRUE("Expand[commonDenominator[{1/(x^2 - 1), 1/(x + 1)^2}, x] - (x - 1) (x + 1)^2] === 0");
    /* the denominator clears every coordinate */
    TR_TRUE("Module[{tf = {x/2, 3/(4 x), 1/(x^2 + 1)}, tD}, tD = commonDenominator[tf, x];"
            " AllTrue[tf, PolynomialQ[Together[# tD], x] &] && Exponent[tD, x] === 3]");
}

static void test_hermiteSolve(void) {
    /* B x == 1 (mod x^2+1): B = -x. */
    TR_EQ("hermiteSolve[1, x, x^2 + 1, x]", "-x");
    TR_EQ("hermiteSolve[0, x, x^2 + 1, x]", "0");
    /* B 2 == 6 (mod x - 1) -> 3. */
    TR_EQ("hermiteSolve[6, 2, x - 1, x]", "3");
    /* rational C: B (x^2+3)/(x+5) == x (mod x^3 - 2) */
    TR_TRUE("Module[{tb = hermiteSolve[x, (x^2 + 3)/(x + 5), x^3 - 2, x]},"
            " Exponent[tb, x] < 3 && PolynomialRemainder[Expand[tb (x^2 + 3) - x (x + 5)], x^3 - 2, x] === 0]");
    TR_TRUE("Module[{tb = hermiteSolve[x^4 + 1, (2 x + 1)/(x^2 + 7), (x - 1) (x^2 + 1), x]},"
            " Exponent[tb, x] < 3 && PolynomialRemainder[Expand[tb (2 x + 1) - (x^4 + 1) (x^2 + 7)], (x - 1) (x^2 + 1), x] === 0]");
    /* non-invertible congruence / denominator sharing a factor with V */
    TR_TRUE("catchTrager[hermiteSolve[1, x - 1, x^2 - 1, x]][[1 ;; 2]] === tragerFail[\"Internal\", \"HermiteNonInvertibleCongruence\"]");
    TR_TRUE("catchTrager[hermiteSolve[1, 1/(x - 1), x^2 - 1, x]][[1 ;; 2]] === tragerFail[\"Internal\", \"HermiteDenominatorNotCoprime\"]");
}

static void test_squarefreeDenominatorQ(void) {
    TR_TRUE("squarefreeDenominatorQ[{1/x, 1/(x (x + 1))}, x]");
    TR_FALSE("squarefreeDenominatorQ[{1/x, 1/(x^2 (x + 1))}, x]");
    TR_TRUE("squarefreeDenominatorQ[{x, 3}, x]");
    TR_TRUE("squarefreeDenominatorQ[{0, 0}, x]");
    TR_FALSE("squarefreeDenominatorQ[{1/(x^2 + 1)^2, 0}, x]");
    /* x and x+1 separately squarefree, the LCM too */
    TR_TRUE("squarefreeDenominatorQ[{1/x, 1/(x + 1), 1/(x^2 + 1)}, x]");
    /* a repeated factor that appears only across two coordinates is not repeated */
    TR_TRUE("squarefreeDenominatorQ[{1/(x - 1), 1/(x - 1)}, x]");
    TR_FALSE("squarefreeDenominatorQ[{1/(x - 1), 1/(x^2 - 2 x + 1)}, x]");
}

int main(void) {
    tr_init();
    define_helpers();
    test_infScaledOrder();
    test_regularAtInfinityQ();
    test_integralAtInfinityQ();
    test_shiftCurve();
    test_unshift();
    test_hermiteReduce_bronstein();
    test_hermiteReduce_hand();
    test_hermiteReduce_invariants();
    test_commonDenominator();
    test_hermiteSolve();
    test_squarefreeDenominatorQ();
    printf("trager_shift_hermite: %d assertions passed\n", tr_assertions);
    return 0;
}
