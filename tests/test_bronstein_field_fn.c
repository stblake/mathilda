/* test_bronstein_field_fn.c -- the rational-function coefficient field
 * fieldFn (src/internal/trager/Field.m) and the Trager ideal machinery over it
 * (BRONSTEIN_PLAN.md Phase 1a).
 *
 * The curve of a Bronstein Shape-A integrand lives over k = Q(x, t_1..) in the
 * top monomial theta, so Trager's Normalize/Ideal/PrincipalGen run with a
 * coefficient field whose scalars are rational functions.  Every routine added
 * or generalised for that is tested here directly: hand-checked values, the
 * edge and throw paths, and seeded property tests (division identity, Bezout
 * identity, inverse, norm multiplicativity, ideal products, principality).
 */
#include "bronstein_test_utils.h"

/* fieldFn records, fieldExtQ, fieldScalarSymbols, fieldScalarQ, fieldBase,
 * fieldDegree, nativeQ[F]. */
static void test_records(void) {
    br_run("F = fieldFn[{x}]; F2 = fieldFn[{x}, {t, t^2 - 2, Sqrt[2]}];"
           " F3 = fieldFn[{x, u}, {t, t^3 - t - 1, Root[#^3 - # - 1 &, 1]}]; True");
    BR_EQ("F[\"Type\"]", "\"Fn\"");
    BR_EQ("fieldFn[{x}, None] === F", "True");
    BR_FALSE("fieldExtQ[F]");
    BR_TRUE("fieldExtQ[F2]");
    BR_TRUE("fieldExtQ[fieldNF[t, t^2 - 2, Sqrt[2]]]");
    BR_FALSE("fieldExtQ[fieldQ[]]");
    BR_EQ("fieldScalarSymbols[F]", "{x}");
    BR_EQ("fieldScalarSymbols[F3]", "{x, u, t}");
    BR_EQ("fieldScalarSymbols[fieldQ[]]", "{}");
    BR_EQ("fieldScalarSymbols[fieldNF[t, t^2 - 2, Sqrt[2]]]", "{t}");
    BR_TRUE("fieldScalarQ[F, (x^2 + 1)/x]");
    BR_FALSE("fieldScalarQ[F, th x]");
    BR_TRUE("fieldScalarQ[F3, t u/(x + t)]");
    BR_TRUE("fieldScalarQ[F, 7]");
    BR_EQ("fieldBase[F2] === F", "True");
    BR_EQ("fieldBase[fieldNF[t, t^2 - 2, Sqrt[2]]] === fieldQ[]", "True");
    BR_EQ("fieldBase[fieldGF[7]] === fieldGF[7]", "True");
    BR_EQ("{fieldDegree[F], fieldDegree[F2], fieldDegree[F3]}", "{1, 2, 3}");
    /* the FLINT kernels never see a rational-function field */
    BR_FALSE("Block[{}, nativeQ[F]]");
    BR_FALSE("nativeQ[F2]");
    BR_EQ("nativeQ[fieldQ[]] === nativeQ[]", "True");
}

/* kc over Fn (kcFn, fnCollect, fnInvertT), including the throw on a
 * non-polynomial and canonical zero. */
static void test_kc(void) {
    BR_EQ("kc[F, (x th^2 + th)/x + th^2]", "2 th^2 + th/x");
    BR_EQ("kc[F, (th - x)(th + x) - th^2 + x^2]", "0");
    BR_EQ("kc[F, 3/x]", "3/x");
    BR_EQ("kc[F, 0]", "0");
    BR_TRUE("Coefficient[kc[F, (1 + x) th^2/x + th/(x - 1)], th, 2] === (1 + x)/x");
    BR_THROWS("kc[F, 1/th]", "Internal", "kcFnNotPolynomial");
    BR_THROWS("kc[F, x/(th + x)]", "Internal", "kcFnNotPolynomial");
    /* extension: t-denominators are inverted, numerators reduced mod mu */
    BR_EQ("kc[F2, t^2]", "2");
    BR_TRUE("Together[kc[F2, (t th + 1)/(t + x)] (t + x) - (t th + 1)] /. t -> Sqrt[2] // Together // PossibleZeroQ");
    BR_TRUE("FreeQ[Denominator[Together[kc[F2, 1/(t + x)]]], t]");
    BR_EQ("kc[F2, t^3 - 2 t]", "0");
    BR_TRUE("Exponent[kc[F3, t^5], t] <= 2");
    /* fnInvertT directly */
    BR_EQ("fnInvertT[F2, 1, 3]", "{1, 3}");
    BR_TRUE("Module[{r = fnInvertT[F2, 1, t + 1]}, FreeQ[r[[2]], t] &&"
            " kc[F2, r[[1]] (t + 1)/r[[2]]] === 1]");
    BR_THROWS("fnInvertT[fieldFn[{x}, {t, t^2 - 1, 1}], 1, t - 1]", "Internal", "fnInvertTNotUnit");
    /* fnCollect directly */
    BR_EQ("fnCollect[F, th^2 x + th, x]", "th^2 + th/x");
    BR_EQ("fnCollect[F, 0, x]", "0");
    BR_EQ("fnCollect[F2, t^2 th, 1]", "2 th");
}

/* kInv, kDeg, kLC, kMonic over Fn. */
static void test_scalars(void) {
    BR_EQ("kInv[F, x/(x + 1)]", "(1 + x)/x");
    BR_EQ("kInv[F, 5]", "1/5");
    BR_THROWS("kInv[F, 0]", "Internal", "kInvZero");
    BR_THROWS("kInv[F, th]", "Internal", "kInvNotScalar");
    BR_EQ("kc[F2, (t + x) kInv[F2, t + x]]", "1");
    BR_EQ("kc[F3, (t^2 + u t + x) kInv[F3, t^2 + u t + x]]", "1");
    BR_EQ("kDeg[F, th^3/x + th, th]", "3");
    BR_EQ("kDeg[F, 0, th]", "-Infinity");
    BR_EQ("kLC[F, (1 + x) th^2/x + th, th]", "(1 + x)/x");
    BR_EQ("kMonic[F, x th^2 + th, th]", "th^2 + th/x");
    BR_EQ("kMonic[F, 0, th]", "0");
}

/* kDivRem / kQuo / kRem / kExactQuo over Fn: hand values and the division
 * identity a == q b + r, deg r < deg b on a seeded family. */
static void test_division(void) {
    BR_EQ("kDivRem[F, th^3 + x, x th + 1, th]", "{1/x^3 - th/x^2 + th^2/x, (-1 + x^4)/x^3}");
    BR_EQ("kDivRem[F, th + 1, th^2, th]", "{0, 1 + th}");
    BR_EQ("kExactQuo[F, (th - x)(th + 1/x), th - x, th]", "th + 1/x");
    BR_THROWS("kExactQuo[F, th^2 + 1, th - x, th]", "Internal", "kExactQuoRemainder");
    BR_THROWS("kDivRem[F, th, 0, th]", "Internal", "kDivRemByZero");
    BR_TRUE("Module[{qr = kDivRem[F2, th^2 + t, t th + x, th]},"
            " kc[F2, qr[[1]] (t th + x) + qr[[2]] - th^2 - t] === 0 && kDeg[F2, qr[[2]], th] < 1]");
    BR_TRUE("SeedRandom[11]; rp[d_] := Sum[(RandomInteger[{-3, 3}] + RandomInteger[{-2, 2}] x)/"
            "(1 + RandomInteger[{0, 2}] x^2) th^i, {i, 0, d}];"
            " And @@ Table[Module[{a = rp[RandomInteger[{0, 5}]], b = rp[RandomInteger[{1, 3}]] + th^4, qr},"
            "   qr = kDivRem[F, a, b, th];"
            "   kc[F, qr[[1]] b + qr[[2]] - a] === 0 && kDeg[F, qr[[2]], th] < kDeg[F, b, th]], {25}]");
}

/* kExtGCD (Euclid over Fn) and kGCD (the Q(vars)[x] PolynomialGCD fast path,
 * and Euclid over an extension): Bezout identity, monic gcd, agreement. */
static void test_gcd(void) {
    BR_EQ("kExtGCD[F, (th - x)(th + 1), (th - x)(th - 2), th]", "{th - x, 1/3, -1/3}");
    BR_EQ("kExtGCD[F, 0, 0, th]", "{0, 0, 0}");
    BR_EQ("kExtGCD[F, x th, 0, th]", "{th, 1/x, 0}");
    BR_EQ("kGCD[F, (th - x)(th + 1) x, (th - x)(th - 2), th]", "th - x");
    BR_EQ("kGCD[F, 0, x th + 1, th]", "th + 1/x");
    BR_EQ("kGCD[F, th^2 + x, th + 1, th]", "1");
    BR_EQ("kGCD[F2, (th - t)(th + x), (th - t)(th + 1), th]", "-t + th");
    BR_TRUE("SeedRandom[12]; rp[d_] := Sum[(RandomInteger[{-3, 3}] + RandomInteger[{-1, 1}] x) th^i, {i, 0, d}] + th^(d + 1);"
            " And @@ Table[Module[{c = rp[1], a, b, eg},"
            "   a = kc[F, c rp[RandomInteger[{0, 2}]]]; b = kc[F, c rp[RandomInteger[{0, 2}]]];"
            "   eg = kExtGCD[F, a, b, th];"
            "   kc[F, eg[[2]] a + eg[[3]] b - eg[[1]]] === 0 && kLC[F, eg[[1]], th] === 1 &&"
            "   kRem[F, a, eg[[1]], th] === 0 && kRem[F, b, eg[[1]], th] === 0 &&"
            "   kc[F, kGCD[F, a, b, th] - eg[[1]]] === 0], {20}]");
    BR_TRUE("SeedRandom[13]; And @@ Table[Module[{c = th + RandomInteger[{-3, 3}] t + x, a, b, eg},"
            "   a = kc[F2, c (th + RandomInteger[{1, 4}])]; b = kc[F2, c (th - RandomInteger[{1, 4}] t)];"
            "   eg = kExtGCD[F2, a, b, th];"
            "   kc[F2, eg[[2]] a + eg[[3]] b - eg[[1]]] === 0 && kc[F2, eg[[1]] - c] === 0], {20}]");
}

/* kTrace, kNorm over Fn. */
static void test_trace_norm(void) {
    BR_EQ("kTrace[F, x th]", "th x");
    BR_EQ("kTrace[F2, t]", "0");
    BR_EQ("kTrace[F2, 3 + x t]", "6");
    BR_EQ("kNorm[F2, t th + x]", "-2 th^2 + x^2");
    BR_EQ("kNorm[F, x th + 1]", "1 + th x");
    BR_TRUE("Together[kNorm[F2, (t th + x)/(x + 1)] - (x^2 - 2 th^2)/(x + 1)^2] === 0");
    /* the norm is multiplicative */
    BR_TRUE("SeedRandom[14]; And @@ Table[Module[{a = th + RandomInteger[{-3, 3}] t + x,"
            "   b = RandomInteger[{1, 3}] t th + RandomInteger[{-2, 2}]},"
            "   Together[kNorm[F2, kc[F2, a b]] - kNorm[F2, a] kNorm[F2, b]] === 0], {20}]");
}

/* rc / rcFn and the rational-function helpers built on it; afPolynomialQ. */
static void test_rc(void) {
    BR_EQ("rc[F, (th^2 - x^2)/(th - x)]", "th + x");
    BR_EQ("rc[F, th/x - th/x]", "0");
    BR_TRUE("rcZeroQ[F, 1/(th + x) - x/(x th + x^2)]");
    BR_FALSE("rcZeroQ[F, 1/th]");
    BR_EQ("rc[F2, 1/(t th + 1)]", "(-1 + t th)/(-1 + 2 th^2)");
    BR_TRUE("FreeQ[Denominator[rc[F2, (th + 1)/(t x + th)]], t]");
    BR_TRUE("rcZeroQ[F2, rc[F2, 1/(t + th)] (t + th) - 1]");
    BR_EQ("rcDeg[F, (th^3 + x)/(x th + 1), th]", "2");
    BR_EQ("rcLeading[F, (x th^3 + 1)/((x + 1) th^3), th]", "x/(1 + x)");
    BR_TRUE("afPolynomialQ[F, {th/x, (1 + x)/(x - 1)}]");
    BR_FALSE("afPolynomialQ[F, {1/th, 0}]");
    BR_TRUE("afPolynomialQ[fieldQ[], {x^2, 1}]");
    BR_FALSE("afPolynomialQ[fieldQ[], {1/x, 1}]");
}

/* Ideals over Fn (Phase 1a spike): HNF, inverse, products, normalisation at
 * infinity and the principality test, on curves Y^m = G(th) over Q(x). */
static void test_ideals(void) {
    /* genus 0: y^2 = th^2 + x; u1 = y + th + 1, u2 = y + th + 2 */
    br_run("B0 = curveBasis[2, th^2 + x, th]; u1 = {th + 1, 1}; u2 = {th + 2, 1};"
           " I1 = idealFromGenerators[F, B0, {u1}]; I2 = idealFromGenerators[F, B0, {u2}]; True");
    BR_EQ("I1", "<|\"M\" -> {{1, 2/(1 + x)}, {0, th + 1/2 (1 - x)}}, \"den\" -> 1|>");
    BR_TRUE("idealEqualQ[F, B0, idealTimes[F, B0, I2, idealInverse[F, B0, I2]], idealUnit[B0]]");
    BR_TRUE("idealEqualQ[F, B0, idealTimes[F, B0, I1, I2], idealTimes[F, B0, I2, I1]]");
    br_run("J0 = idealTimes[F, B0, I1, idealInverse[F, B0, I2]]; g0 = principalGenerator[F, B0, J0]; True");
    BR_FALSE("g0 === None");
    /* the generator is u1/u2 up to a constant of k */
    BR_TRUE("Module[{r = afTimes[F, B0, afTimes[F, B0, g0, u2], afInverse[B0, u1]]},"
            " r[[2]] === 0 && fieldScalarQ[F, r[[1]]] && r[[1]] =!= 0]");
    BR_TRUE("principalQ[F, B0, idealUnit[B0]]");
    /* genus 1: y^2 = th^3 + x^2, P = (0, x), Q = (0, -x); P - Q is not principal */
    br_run("B1 = curveBasis[2, th^3 + x^2, th];"
           " IP = idealFromGenerators[F, B1, {{th, 0}, {-x, 1}}];"
           " IQ = idealFromGenerators[F, B1, {{th, 0}, {x, 1}}]; True");
    BR_EQ("IP", "<|\"M\" -> {{1, -1/x}, {0, th}}, \"den\" -> 1|>");
    BR_EQ("principalGenerator[F, B1, idealTimes[F, B1, IP, idealInverse[F, B1, IQ]]]", "None");
    /* idealReduce removes a common factor of den over Fn */
    BR_TRUE("Module[{I = idealReduce[F, B0, <|\"M\" -> (th + x) IdentityMatrix[2], \"den\" -> th + x|>]},"
            " I[\"den\"] === 1 && I[\"M\"] === IdentityMatrix[2]]");
    /* the same over Q(sqrt 2)(x): residue ideals there are over an extension */
    br_run("IA = idealFromGenerators[F2, B0, {{th + t, 1}}];"
           " IB = idealFromGenerators[F2, B0, {{th + 2 t, 1}}]; True");
    BR_TRUE("idealEqualQ[F2, B0, idealTimes[F2, B0, IB, idealInverse[F2, B0, IB]], idealUnit[B0]]");
    BR_TRUE("Module[{g = principalGenerator[F2, B0, idealTimes[F2, B0, IA, idealInverse[F2, B0, IB]]], r},"
            " g =!= None && (r = afTimes[F2, B0, afTimes[F2, B0, g, {th + 2 t, 1}],"
            "   afInverse[B0, {th + t, 1}] /. t -> Sqrt[2]] /. t -> Sqrt[2];"
            "  PossibleZeroQ[r[[2]]] && FreeQ[Together[r[[1]]], th])]");
    /* commonDenFactorFn directly */
    BR_EQ("commonDenFactorFn[F2, th + x, {{t (th + x), 0}, {0, (th + x) th}}, th]", "th + x");
    BR_EQ("commonDenFactorFn[F2, th + x, {{t, 0}, {0, th + x}}, th]", "1");
}

/* Property: on random genus-0 curves y^2 = th^2 + a(x) and random affine
 * u1 = y + th + c1, u2 = y + th + c2, the divisor (u1) - (u2) at the finite
 * places is principal with generator a k-multiple of u1/u2. */
static void test_principal_property(void) {
    BR_TRUE("SeedRandom[21]; And @@ Table[Module[{a = RandomInteger[{1, 5}] x + RandomInteger[{1, 5}],"
            "   c1 = RandomInteger[{1, 6}], c2, B, J, g, r},"
            "   c2 = c1 + RandomInteger[{1, 4}];"
            "   B = curveBasis[2, th^2 + a, th];"
            "   J = idealTimes[F, B, idealFromGenerators[F, B, {{th + c1, 1}}],"
            "     idealInverse[F, B, idealFromGenerators[F, B, {{th + c2, 1}}]]];"
            "   g = principalGenerator[F, B, J];"
            "   g =!= None && (r = afTimes[F, B, afTimes[F, B, g, {th + c2, 1}], afInverse[B, {th + c1, 1}]];"
            "     rcZeroQ[F, r[[2]]] && fieldScalarQ[F, r[[1]]])], {20}]");
}

int main(void) {
    br_init();
    test_records();
    test_kc();
    test_scalars();
    test_division();
    test_gcd();
    test_trace_norm();
    test_rc();
    test_ideals();
    test_principal_property();
    BR_DONE("bronstein_field_fn_tests");
    return 0;
}
