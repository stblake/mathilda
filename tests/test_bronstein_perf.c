/* test_bronstein_perf.c -- batch E of the gap register (performance), and the
 * soundness bug its profiling found (tests/bronstein_paper/open_gaps.m):
 *
 *  - G-adv-t2: Resultant over Gaussian rationals returned an uncancelled
 *    fraction, and Trager's logDerivativeNF then divided it as a polynomial
 *    forever; the norm is now taken on the determinant's numerator;
 *  - G-adv-t3: reduceRadicalPowers ran Together over hundreds of copies of a
 *    degree-8 Root; the Roots are opaque symbols there;
 *  - G-perf-1: powerBasisInverse inverts a denominator free of the radical (or
 *    a monomial in it) directly, not by an n x n system (degree-30 Kummer curve);
 *  - G-perf-4 (partial): the norm in gaussianRootObjects writes the radicals of
 *    one base as powers of one generator; primitiveFieldByNorm builds the
 *    absolute residue field when ToNumberField gives up; fieldConstRep reads
 *    the constants of such a field;
 *  - G-adv-b5 and a false certificate: the C limited integration folds a new
 *    logarithm or ArcTan by its derivative into the tower, and the ArcTan
 *    obstruction certifies only over a real tower
 *    (D[Log[(u - I)/(u + I)]^2] was certified not elementary).
 */
#include "bronstein_test_utils.h"

static void test_resultant_gaussian(void) {
    /* the resultant of polynomials is a polynomial */
    BR_TRUE("PolynomialQ[Resultant[-1 + x^2 + (4 - 4 I) Z^2 + (4 - 4 I) x Z^2,"
            " -1/32 + I/32 + (-1/4 - I/4) Z^2 + Z^4, Z], x]");
    BR_TRUE("Expand[Resultant[-1 + x^2 + (4 - 4 I) Z^2 + (4 - 4 I) x Z^2, -1/32 + I/32 + (-1/4 - I/4) Z^2 + Z^4, Z]"
            " - (-1 - 4 x + (-6 + 2 I) x^2 + (-4 + 8 I) x^3 + 12 I x^4 + (4 + 8 I) x^5 + (6 + 2 I) x^6 + 4 x^7 + x^8)] === 0");
    BR_TRUE("PolynomialQ[Resultant[(1 + x) Z^2 + I, Z^2 - I Z + 1/2, Z], x]");
    BR_EQ("Resultant[x^2 - 2, x - 1, x]", "-1");
}

static void test_trager_gaussian(void) {
    /* G-adv-t2, G-adv-t3 */
    BR_ANTIDERIV("1/((x^2 + I) Sqrt[x + 1])", "Trager");
    BR_ANTIDERIV("(1 + (1 + I) x^4)^(3/4)", "Trager");
    /* reduceRadicalPowers with Root objects: an identity */
    BR_TRUE("With[{r = Root[#^3 - # - 1 &, 1]}, Together[reduceRadicalPowers[(r Y^3 + x)/(Y^2 + r), Y, 2, x + 1]"
            " - (r Y^3 + x)/(Y^2 + r) /. Y^2 -> x + 1 /. Y^3 -> Y (x + 1)] === 0]");
    BR_TRUE("FreeQ[reduceRadicalPowers[Root[#^2 - 2 &, 2] Y^3, Y, 2, x], Y^3]");
}

static void test_power_basis_inverse(void) {
    BR_EQ("powerBasisInverse[x + 1, Y, 3, x^2, x]", "{1/(1 + x), 0, 0}");
    /* 1/(2 Y) = Y^2/(2 g) on Y^3 = g */
    BR_EQ("powerBasisInverse[2 Y, Y, 3, x^2, x]", "{0, 0, 1/(2 x^2)}");
    BR_TRUE("With[{c = powerBasisInverse[1 + Y, Y, 3, x, x]}, Together[Expand[(1 + Y) c . {1, Y, Y^2}] /. {Y^3 -> x, Y^4 -> x Y}] === 1]");
    /* G-perf-1 */
    BR_EQ("BronsteinIntegrate[Sqrt[x] (x + 1)^(1/5) (x + 2)^(1/3), x][[1]]", "\"not elementary\"");
    BR_ANTIDERIV("D[x^2 Sqrt[x] (x + 1)^(1/5) (x + 2)^(1/3), x]", "Bronstein");
}

static void test_norm_fields(void) {
    /* one generator per base: the norm of a degree-6 polynomial over Q(2^(1/3)) has degree 18 */
    BR_EQ("Length[gaussianRootObjects[8/1493721 - 1/(1493721 2^(2/3)) - (64 2^(2/3))/1493721 + (2 Z^4)/6147"
          " + (256 2^(1/3) Z^4)/6147 - (16 2^(2/3) Z^4)/6147 + Z^6, Z]]", "6");
    BR_TRUE("With[{r = gaussianRootObjects[Z^2 - 2^(1/3), Z]}, Max[Abs[N[(#^2 - 2^(1/3)) & /@ r, 30]]] < 10^-20]");
    /* Trager's primitive element: the field of a root of Z^2 - 2^(1/3) Z + 1 */
    BR_TRUE("Module[{q = Z^2 - 2^(1/3) Z + 1, r1, pf, F, g, c},"
            " r1 = First[gaussianRootObjects[q, Z]];"
            " pf = primitiveFieldByNorm[q, Z, t, r1, {2^(1/3), 2^(2/3)}];"
            " F = pf[\"F\"]; g = pf[\"gamma\"]; c = 2^(1/3) /. F[\"consts\"];"
            " Exponent[F[\"mu\"], t] === 6 &&"
            " PolynomialRemainder[Expand[(q /. 2^(1/3) -> c /. 2^(2/3) -> c^2) /. Z -> g], F[\"mu\"], t] === 0 &&"
            " PolynomialRemainder[Expand[c^3 - 2], F[\"mu\"], t] === 0]");
    BR_EQ("primitiveFieldByNorm[Z^2 - Sqrt[2] Z + Sqrt[3], Z, t, 1, {Sqrt[2], Sqrt[3]}]", "None");   /* two bases */
    /* fieldConstRep on such a field: powers and a binomial Root */
    BR_TRUE("Module[{q = Z^2 - 2^(1/3) Z + 1, pf, F, c},"
            " pf = primitiveFieldByNorm[q, Z, t, First[gaussianRootObjects[q, Z]], {2^(1/3)}]; F = pf[\"F\"];"
            " c = fieldConstRep[F, 2^(1/3)];"
            " PolynomialRemainder[Expand[fieldConstRep[F, 2^(2/3)] - c^2], F[\"mu\"], t] === 0 &&"
            " PolynomialRemainder[Expand[fieldConstRep[F, 2^(-1/3)] c - 1], F[\"mu\"], t] === 0 &&"
            " PolynomialRemainder[Expand[fieldConstRep[F, Root[-16 + #^3 &, 1]] - 2 c], F[\"mu\"], t] === 0]");
    /* G-perf-4 neighbours that answer now */
    BR_ANTIDERIV("1/((x^2 - 2^(2/3)) Sqrt[x + 1])", "Trager");
    BR_ANTIDERIV("1/((x - 2^(1/3)) Sqrt[x^2 + 1])", "Trager");
}

static void test_c_log_fold(void) {
    /* the false certificate (pre-existing): D[t^2] with t = Log[(u - I)/(u + I)] */
    BR_TRUE("Risch`ElementaryIntegralQ[Together[D[Log[(u - I)/(u + I)]^2, u]], u]");
    BR_TRUE("Risch`ElementaryIntegralQ[Together[D[Log[(u - I)/(u + I)] Log[1 + u^2], u]], u]");
    BR_TRUE("Risch`ElementaryIntegralQ[Together[D[Log[(u - I)/(u + I)] Log[u^2 - 1], u]], u]");
    BR_TRUE("Integrate`RischTranscendental[Together[D[Log[(u - I)/(u + I)]^2, u]], u] === Log[(u - I)/(u + I)]^2");
    /* the real dilogarithm obstruction still certifies */
    BR_EQ("Risch`ElementaryIntegralQ[Log[1 + x^2]/(1 + x^2), x]", "False");
    BR_EQ("Risch`ElementaryIntegralQ[(3 x + 1) Log[x^2 + 1]^5, x]", "False");
    BR_EQ("Risch`ElementaryIntegralQ[Log[x]/(1 + x^2), x]", "False");
    /* products of logarithms of quotients (timed out in the flat ansatz) */
    BR_TRUE("TimeConstrained[With[{G = Log[u/(1 + u)] Log[(u - 2)/(u + 2)]},"
            " PossibleZeroQ[D[Integrate`RischTranscendental[Together[D[G, u]], u] - G, u]]], 20, False]");
    BR_TRUE("TimeConstrained[With[{G = Log[u/(1 + u^2)] Log[(u - I)/(u + I)]},"
            " PossibleZeroQ[D[Integrate`RischTranscendental[Together[D[G, u]], u] - G, u], Assumptions -> Element[u, Complexes]]], 20, False]");
    /* G-adv-b5 */
    BR_ANTIDERIV("Together[D[ArcSin[x] Log[x], x]]", "Bronstein");
}

int main(void) {
    br_init();
    test_resultant_gaussian();
    test_trager_gaussian();
    test_power_basis_inverse();
    test_norm_fields();
    test_c_log_fold();
    BR_DONE("bronstein_perf_tests");
    return 0;
}
