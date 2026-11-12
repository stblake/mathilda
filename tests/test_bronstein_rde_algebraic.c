/* test_bronstein_rde_algebraic.c -- the Risch differential equation
 * z' + f z = g over Q(x)(Y), Y^m = G(x) (src/internal/bronstein/RDEAlgebraic.m;
 * Bronstein 1990 Part 3 with k = Q, theta = x; ISSAC'91).
 *
 * Every routine directly; the ISSAC'91 equation; solutions on square- and
 * cube-root curves (each matched against the known z); certified "no
 * solution" cases (each with a mathematical reason: an elliptic integral, a
 * non-elementary exponential integral); and the bound-validity test: widening
 * every bound (brAlgRDE[.., widen]) never changes solvability.
 */
#include "bronstein_test_utils.h"

static void test_valuations(void) {
    BR_TRUE("brAlgNuP[x^2 (x + 1)/(x - 1)^3, x - 1, x] === -3");
    BR_TRUE("brAlgNuP[x^2 (x + 1)/(x - 1)^3, x, x] === 2");
    BR_TRUE("brAlgNuP[0, x, x] === Infinity");
    BR_TRUE("brAlgDegree[(x^3 + 1)/x, x] === 2");
    BR_TRUE("brAlgDegree[0, x] === -Infinity");
    BR_TRUE("Sort[brAlgPoleFactors[{{1/x, 1/(x^2 + 1)}, {x, 1/(x (x + 1))}}, x]] === Sort[{x, 1 + x^2, 1 + x}]");
    BR_TRUE("brAlgRamification[2, x^3 + 1, x + 1, x] === 2");
    BR_TRUE("brAlgRamification[2, x^3 + 1, x, x] === 1");
    BR_TRUE("brAlgRamification[3, x^2 + 1, x^2 + 1, x] === 3");
}

static void test_bounds(void) {
    br_run("B3 = curveBasis[2, 1 + x^3, x];"
           " fe = afFromY[fieldQ[], B3, (5 x^3 + 2)/(2 (1 + x^3)) Y, Y]; True");
    /* tau: the value of r f P/P' at the places over P; only positive integers count */
    BR_TRUE("brAlgTauCandidates[B3, {-2/x, 0}, x, 1] === {}");
    BR_TRUE("brAlgTauCandidates[B3, {1/x, 0}, x, 1] === {1}");
    BR_TRUE("brAlgInfinityData[B3] === <|\"r\" -> 2, \"q\" -> 3, \"lc\" -> 1, \"d\" -> 3|>");
    BR_TRUE("brAlgInfinityCandidates[B3, {-3/x, 0}] === {-3}");
    BR_TRUE("brAlgInfinityCandidates[B3, fe] === {0}");
    BR_TRUE("brAlgGammaInfinity[B3, {x, 1}] === -3");
    BR_TRUE("brAlgGammaInfinity[B3, {0, 0}] === Infinity");
    BR_TRUE("brAlgBounds[B3, fe, {1, 1}] === <|\"D\" -> 1, \"Bdeg\" -> {2, 1}|>");
    BR_TRUE("brAlgDLinear[B3, {c1 x, c2}] === {c1, Together[3 c2 x^2/(2 (1 + x^3))]}");
    BR_TRUE("brAlgTimesLinear[B3, {0, 1}, {c1, c2}] === {c2 (1 + x^3), c1}");
    BR_TRUE("Together[brAlgReduceY[1/(1 + Y), Y, 2, x^3 + 1] - (Y - 1)/x^3] === 0");
}

static void test_solve(void) {
    /* dY[e]: the derivation on Q(x)(Y), Y' = G'/(m G) Y */
    br_run("dY[e_, m_, G_] := D[e, x] + D[e, Y] D[G, x]/(m G) Y;"
           " rdeMatch[f_, z_, m_, G_] := Module[{g = Together[dY[z, m, G] + f z], r},"
           "   r = brRDEAlgebraic[f, g, x, m, G, Y];"
           "   MatchQ[r, {\"Solution\", _}] && Together[brAlgReduceY[r[[2]] - z, Y, m, G]] === 0]; True");
    /* the ISSAC'91 equation: z' + (5x^3 + 2)/(2 Sqrt[x^3 + 1]) z = .., z = (2/x)(1 + Sqrt[x^3 + 1]) */
    BR_TRUE("rdeMatch[(5 x^3 + 2)/(2 Y), (2/x)(1 + Y), 2, x^3 + 1]");
    BR_TRUE("rdeMatch[2 x, x Y, 2, x^3 + 1]");
    BR_TRUE("rdeMatch[0, Y, 2, x^3 + 1]");
    BR_TRUE("rdeMatch[x Y, Y^2/x + x, 3, x^2 + 1]");
    BR_TRUE("rdeMatch[1/x, Y/(x - 1)^2 + 3, 2, x^3 + 1]");
    BR_TRUE("rdeMatch[-(3 x^2/(2 (1 + x^3))) + 1/(x + 2), Y^2/(x + 2), 3, x^2 + 1]");
    BR_TRUE("brRDEAlgebraic[x, 0, x, 2, x^3 + 1, Y] === {\"Solution\", 0}");
    /* certified "no solution" */
    BR_TRUE("brRDEAlgebraic[0, 1/Y, x, 2, x^3 + 1, Y] === {\"NoSolution\"}");        /* elliptic integral */
    BR_TRUE("brRDEAlgebraic[1, Y, x, 2, x^3 + 1, Y] === {\"NoSolution\"}");          /* Int E^x Y not elementary */
    BR_TRUE("brRDEAlgebraic[Y, 1, x, 2, x^3 + 1, Y] === {\"NoSolution\"}");
    BR_TRUE("brAlgRDE[B3, fe, afFromY[fieldQ[], B3, 1, Y]] === {\"NoSolution\"}");
}

static void test_bound_validity(void) {
    /* widening every bound never changes solvability (a brute-force check of the bounds) */
    br_run("widenAgree[f_, g_, m_, G_] := Module[{B = curveBasis[m, G, x], fe, ge, r0, r3},"
           "   fe = afFromY[fieldQ[], B, brAlgReduceY[f, Y, m, G], Y]; ge = afFromY[fieldQ[], B, brAlgReduceY[g, Y, m, G], Y];"
           "   r0 = brAlgRDE[B, fe, ge]; r3 = brAlgRDE[B, fe, ge, 3]; First[r0] === First[r3]]; True");
    BR_TRUE("widenAgree[(5 x^3 + 2)/(2 Y), Together[dY[(2/x)(1 + Y), 2, x^3 + 1] + (5 x^3 + 2)/(2 Y) (2/x)(1 + Y)], 2, x^3 + 1]");
    BR_TRUE("widenAgree[0, 1/Y, 2, x^3 + 1]");
    BR_TRUE("widenAgree[1, Y, 2, x^3 + 1]");
    BR_TRUE("widenAgree[1/x, Y/x, 2, x^3 + 1]");
    BR_TRUE("widenAgree[x Y, Y/(x + 1), 3, x^2 + 1]");
    /* seeded: z random on random curves, f random; the solution is always found */
    br_run("SeedRandom[1112]; rr[] := RandomInteger[{-3, 3}];"
           " curves = {{2, x^3 + 1}, {2, x^4 + 2}, {3, x^2 + 1}, {2, x^3 - x + 1}}; True");
    BR_TRUE("And @@ Table[Module[{cG = RandomChoice[curves], z, f},"
            "   z = (rr[] + x + Y)/(x - 1 - rr[]^2) + rr[] x^2; f = rr[] + Y/(x + 3);"
            "   rdeMatch[f, z, cG[[1]], cG[[2]]]], {8}]");
}

int main(void) {
    br_init();
    test_valuations();
    test_bounds();
    test_solve();
    test_bound_validity();
    BR_DONE("bronstein_rde_algebraic_tests");
    return 0;
}
