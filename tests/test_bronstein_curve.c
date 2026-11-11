/* test_bronstein_curve.c -- the integrand on the curve over k(theta) and the
 * per-component curves (src/internal/bronstein/Curve.m), BRONSTEIN_PLAN.md
 * Phase 2.  Hand-checked canonical forms of the paper examples, every decline
 * path, and seeded properties: the canonical form reproduces the integrand,
 * and each component generator Y = yh^i/s satisfies Y^m = G exactly.
 */
#include "bronstein_test_utils.h"

static void test_substitute(void) {
    BR_EQ("brSubstituteRadical[Sqrt[t + x] + (t + x)^(-3/2), {x, t}, t + x, 2, y]", "1/y^3 + y");
    BR_EQ("brSubstituteRadical[(1/(t + x))^(1/2), {x, t}, t + x, 2, y]", "1/y");
    BR_EQ("brSubstituteRadical[Sqrt[x] + Sqrt[t], {x, t}, t, 2, y]", "Sqrt[x] + y");
}

static void test_canonical(void) {
    br_run("T134 = brTower[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x]; c134 = brCanonicalForm[T134];"
           " T158 = brTower[(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3)), x];"
           " c158 = brCanonicalForm[T158];"
           " T139 = brTower[Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), x]; c139 = brCanonicalForm[T139]; True");
    BR_EQ("{c134[\"n\"], c134[\"g\"], c134[\"coeffs\"]}", "{2, brT1 + x, {0, (1 + x)/(x (brT1 + x)^2)}}");
    BR_EQ("c134[\"yexpr\"]", "Sqrt[brT1 + x]");
    BR_EQ("c158[\"coeffs\"]", "{3/x, 0, (3 brT1 + 5 x + 2 brT1 x)/(brT1 + x)}");
    BR_EQ("{c139[\"g\"], c139[\"coeffs\"]}", "{brT2, {0, 1/(1 + brT2), 0}}");
    /* a radicand with a denominator: yh = q b^(1/n), g = q^(n-1) p */
    BR_EQ("With[{c = brCanonicalForm[brTower[Sqrt[Log[x]/(x + 1)], x]]}, {c[\"q\"], c[\"g\"], c[\"coeffs\"]}]",
          "{1 + x, brT1 + brT1 x, {0, 1/(1 + x)}}");
    BR_THROWS("brCanonicalForm[brTower[E^x + Log[x], x]]", "Internal", "brCanonicalFormNoRadical");
    BR_THROWS("brCanonicalForm[brTower[1/(Sqrt[Log[x]^2] - Log[x]), x]]", "Decline",
              "the denominator is a zero divisor: the radical is reducible");
}

static void test_components(void) {
    BR_EQ("brComponentCurve[T134, 2, brT1 + x, 1]", "<|\"m\" -> 2, \"G\" -> brT1 + x, \"s\" -> 1, \"exact\" -> True|>");
    BR_EQ("brComponentCurve[T158, 3, brT1 + x, 2]", "<|\"m\" -> 3, \"G\" -> brT1^2 + 2 brT1 x + x^2, \"s\" -> 1, \"exact\" -> True|>");
    /* yh^3 = theta^4 (theta + x): component 1 is theta * Y, Y^3 = theta (theta + x) */
    BR_EQ("brComponentCurve[T158, 3, brT1^4 (brT1 + x), 1]", "<|\"m\" -> 3, \"G\" -> brT1^2 + brT1 x, \"s\" -> brT1, \"exact\" -> True|>");
    /* d > 1 with a rational content: the locally constant rho is allowed */
    BR_EQ("brComponentCurve[T158, 4, 4 brT1^2, 1][\"exact\"]", "False");
    BR_THROWS("brComponentCurve[T158, 4, x brT1^2, 1]", "Decline",
              "component over a radical extension of the constant field (non-constant content)");
    BR_THROWS("brComponentCurve[T158, 2, x brT1^2, 1]", "Decline",
              "component over a radical extension of the constant field (non-constant content)");
    /* brCurve: s = (1/m) DG/G, content and primitive part, Trager basis over k[theta] */
    br_run("cv = brCurve[T134, 2, x^2 (brT1 + x)]; True");
    BR_EQ("{cv[\"s\"], cv[\"content\"], cv[\"primitive\"]}", "{(1 + 2 brT1 + 3 x)/(2 x (brT1 + x)), x^2, brT1 + x}");
    BR_EQ("cv[\"B\"][\"d\"]", "{1, 1}");
    BR_EQ("brCurve[T158, 3, brT1 + x][\"s\"]", "(1 + brT1)/(3 (brT1 + x))");
    BR_TRUE("brReducedCurveQ[T134, 2, brT1 + x]");
    BR_TRUE("brReducedCurveQ[T158, 3, (brT1 + x)^2]");
    BR_FALSE("brReducedCurveQ[T158, 3, (brT1 + x)^3]");
    BR_FALSE("brReducedCurveQ[T158, 4, (brT1 + x)^2]");
    BR_FALSE("brReducedCurveQ[T158, 1, brT1 + x]");
    BR_FALSE("brReducedCurveQ[T158, 2, x + 1]");
}

/* Properties over a seeded random family of Shape A integrands. */
static void test_properties(void) {
    br_run("SeedRandom[51]; brR[] := RandomInteger[{-3, 3}];"
           " brA[] := Module[{th = RandomChoice[{Log[x], E^x, Log[x + 1]}], n = RandomInteger[{2, 4}], p},"
           "   p = th + brR[] x + RandomInteger[{1, 3}];"
           "   (1 + brR[] th)/(x + th) p^(RandomInteger[{1, n - 1}]/n) + x/p^(1/n)]; True");
    /* the canonical form reproduces the integrand */
    BR_TRUE("And @@ Table[Module[{f = brA[], T, c}, T = brTower[f, x]; c = brCanonicalForm[T];"
            "   PossibleZeroQ[brToSurface[T, (c[\"coeffs\"] . (c[\"y\"]^Range[0, c[\"n\"] - 1])) /. c[\"y\"] -> c[\"yexpr\"]] - f]], {25}]");
    /* every component generator satisfies Y^m = G (mod yh^n = g), exactly */
    BR_TRUE("And @@ Table[Module[{f = brA[], T, c, ok = True}, T = brTower[f, x]; c = brCanonicalForm[T];"
            "   Do[With[{cc = brComponentCurve[T, c[\"n\"], c[\"g\"], i]},"
            "     If[cc[\"exact\"], ok = ok && Together[reduceRadicalPowersProbe[(c[\"y\"]^i/cc[\"s\"])^cc[\"m\"], c] - cc[\"G\"]] === 0]],"
            "     {i, 1, c[\"n\"] - 1}]; ok], {20}]");
}

int main(void) {
    br_init();
    /* (yh^i/s)^m with powers of yh reduced modulo yh^n = g */
    br_run("reduceRadicalPowersProbe[e_, c_] := Module[{t = Together[e]},"
           " Together[reduceYDegree[Numerator[t], c[\"y\"], c[\"n\"], c[\"g\"]]/"
           "   reduceYDegree[Denominator[t], c[\"y\"], c[\"n\"], c[\"g\"]]]]; True");
    test_substitute();
    test_canonical();
    test_components();
    test_properties();
    BR_DONE("bronstein_curve_tests");
    return 0;
}
