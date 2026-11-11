/* test_bronstein_rde.c -- the tri-state Risch DE (BRONSTEIN_PLAN.md Phase 1c):
 *
 *   Risch`RischDECertified[f, g, x (, kernels)]   (C, integrate_risch_transcendental.c)
 *   Risch`RischDECertifiedReason[]                (C)
 *   brRischDE / brRDEResidualZeroQ / brRDENoSolutionQ / brRDESolution
 *                                                 (src/internal/bronstein/RDE.m)
 *
 * The decisive property is SOUNDNESS of "NoSolution": for a random y in the
 * tower and g = D[y] + f y the answer must never be "NoSolution" (and a returned
 * solution must satisfy the identity).  The suite also pins hand-checked
 * certificates (Q(x) and tower), each undecided path by its reason, agreement
 * with Risch`RischDE wherever the latter solves, and the regression of the
 * non-canonical leading coefficient in the polynomial loops (a false "no
 * solution" for y = E^x + 1/(1+x), f = E^x).
 */
#include "bronstein_test_utils.h"

/* Solutions: base field and tower, with and without extra kernels. */
static void test_solutions(void) {
    BR_EQ("Risch`RischDECertified[1, x, x]", "-1 + x");
    BR_EQ("Risch`RischDECertified[0, 0, x]", "0");
    BR_EQ("Risch`RischDECertified[0, 1/x, x, {Log[x]}]", "Log[x]");
    BR_EQ("Risch`RischDECertified[0, Log[x], x]", "-x + x Log[x]");
    BR_EQ("Risch`RischDECertified[-2 x, E^(x^2), x]", "E^x^2 x");
    BR_EQ("Risch`RischDECertified[1/(2 Log[x] x), 1/x, x]", "2/3 Log[x]");
    /* [B] p.134: the algebraic part of (x+1)/((x Log x + x^2) Sqrt[x + Log x]) */
    BR_EQ("Risch`RischDECertified[(1 + 1/x)/(2 (x + Log[x])), (x + 1)/(x (x + Log[x])), x]", "2");
    /* regression: the polynomial loops canonicalize c (was a false NoSolution) */
    BR_TRUE("Module[{y = E^x + 1/(1 + x), g}, g = Together[D[y, x] + E^x y];"
            " Together[Risch`RischDECertified[E^x, g, x, {E^x}] - y] === 0]");
    BR_TRUE("Module[{y = (-2 + E^(2 x) + 3 x)/(1 + x), g}, g = Together[D[y, x] + E^x y];"
            " Together[Risch`RischDE[E^x, g, x] - y] === 0]");
}

/* Certified "NoSolution": every one is a theorem (checked by hand). */
static void test_certificates(void) {
    BR_EQ("Risch`RischDECertified[0, 1/x, x]", "\"NoSolution\"");            /* Log x not in Q(x) */
    BR_EQ("Risch`RischDECertified[-2 x, 1, x]", "\"NoSolution\"");           /* Erf */
    BR_EQ("Risch`RischDECertified[2 x, 1, x]", "\"NoSolution\"");            /* Erfi */
    BR_EQ("Risch`RischDECertified[-1/x, 1, x]", "\"NoSolution\"");           /* x Log x */
    BR_EQ("Risch`RischDECertified[0, 1/Log[x], x]", "\"NoSolution\"");       /* li */
    BR_EQ("Risch`RischDECertifiedReason[]", "None");
}

/* Undecided paths, each by its recorded reason. */
static void test_undecided(void) {
    BR_EQ("Risch`RischDECertified[1, Sqrt[x], x]", "$Failed");
    BR_EQ("Risch`RischDECertifiedReason[]", "\"base: coefficients not rational in x\"");
    BR_EQ("Risch`RischDECertified[1, E^x/x, x]", "$Failed");
    BR_EQ("Risch`RischDECertifiedReason[]", "\"tower: antidifferentiation D H = c\"");
    BR_EQ("Risch`RischDECertified[1 + 1/x, Log[x], x]", "$Failed");
    BR_EQ("Risch`RischDECertifiedReason[]", "\"RdeBoundDegreePrim: cancellation case\"");
    /* f with a normal simple pole of positive-integer residue (1 at Log x = 0):
     * RdeNormalDenominator's bound is unproven */
    BR_EQ("Risch`RischDECertified[1/(x Log[x]), x, x]", "$Failed");
    BR_EQ("Risch`RischDECertifiedReason[]", "\"tower: f not weakly normalized\"");
    BR_EQ("Risch`RischDECertified[1, Sin[x], x]", "$Failed");
    /* malformed calls stay unevaluated */
    BR_EQ("Head[Risch`RischDECertified[1, x]]", "Risch`RischDECertified");
    BR_EQ("Head[Risch`RischDECertified[1, x, 2]]", "Risch`RischDECertified");
    BR_EQ("Head[Risch`RischDECertified[1, x, x, Log[x]]]", "Risch`RischDECertified");
    BR_EQ("Head[Risch`RischDECertifiedReason[1]]", "Risch`RischDECertifiedReason");
    /* attributes and documentation */
    BR_TRUE("MemberQ[Attributes[Risch`RischDECertified], Protected]");
    br_assertions++;
    if (!symtab_get_docstring("Risch`RischDECertified") ||
        !symtab_get_docstring("Risch`RischDECertifiedReason")) {
        fprintf(stderr, "FAIL: Risch`RischDECertified(Reason) has no docstring\n");
        exit(1);
    }
    BR_TRUE("MemberQ[Attributes[Risch`RischDECertifiedReason], Protected]");
}

/* Soundness and agreement on seeded random families over several towers. */
static void test_properties(void) {
    br_run("brks = {{Log[x]}, {E^x}, {Log[x], E^x}, {Log[1 + x]}, {E^(x^2)}, {E^x, Log[E^x + 1]}};"
           " rnd[] := RandomInteger[{-3, 3}];"
           " brY[t_] := (rnd[] + rnd[] x + rnd[] t + rnd[] t^2)/"
           "   (1 + RandomInteger[{0, 2}] x + RandomChoice[{0, 1}] t^RandomInteger[{1, 2}]);"
           " brF[t_] := RandomChoice[{0, 1, x, 1/x, rnd[] D[t, x]/t, 1/(x t), (1 + x)/(t + x), rnd[] + D[t, x]}];"
           " True");
    /* never a false NoSolution; every solution satisfies the identity */
    BR_TRUE("SeedRandom[31]; And @@ Table[Module[{k = RandomChoice[brks], t, y, f, g, r},"
            "   t = Last[k]; y = brY[t]; f = brF[t]; g = Together[D[y, x] + f y];"
            "   r = Risch`RischDECertified[f, g, x, k];"
            "   r =!= \"NoSolution\" && (r === $Failed || PossibleZeroQ[D[r, x] + f r - g])], {80}]");
    /* agreement: where RischDE solves, RischDECertified solves with the same y
     * (the solution is unique unless the homogeneous equation has one) */
    BR_TRUE("SeedRandom[32]; And @@ Table[Module[{k = RandomChoice[brks], t, y, f, g, r0, r1},"
            "   t = Last[k]; y = brY[t]; f = brF[t]; g = Together[D[y, x] + f y];"
            "   r0 = Risch`RischDE[f, g, x]; r1 = Risch`RischDECertified[f, g, x, k];"
            "   Head[r0] === Risch`RischDE || r1 === $Failed ||"
            "   PossibleZeroQ[D[r1 - r0, x] + f (r1 - r0)]], {40}]");
    /* g = 1/t-type right-hand sides with no elementary antiderivative: never a
     * solution; a certificate when decided */
    BR_TRUE("SeedRandom[33]; And @@ Table[Module[{a = RandomInteger[{1, 4}], r},"
            "   r = Risch`RischDECertified[0, a/Log[x], x, {Log[x]}];"
            "   r === \"NoSolution\" || r === $Failed], {20}]");
}

/* RDE.m wrappers. */
static void test_wrappers(void) {
    BR_EQ("brRischDE[1, x, x, {}]", "{\"Solution\", -1 + x}");
    BR_EQ("brRischDE[0, 1/x, x, {}]", "{\"NoSolution\"}");
    BR_EQ("brRischDE[0, 1/x, x, {Log[x]}]", "{\"Solution\", Log[x]}");
    BR_EQ("brRischDE[1, Sqrt[x], x, {}]", "{\"Unknown\", \"base: coefficients not rational in x\"}");
    BR_EQ("brRischDE[1 + 1/x, Log[x], x, {Log[x]}]", "{\"Unknown\", \"RdeBoundDegreePrim: cancellation case\"}");
    BR_EQ("brRischDE[1/x, 1/x + Log[x], x, {Log[x]}]", "{\"Solution\", 1/4 (4 - x + 2 x Log[x])}");
    BR_TRUE("brRDEResidualZeroQ[1, x, x, x - 1]");
    BR_FALSE("brRDEResidualZeroQ[1, x, x, x]");
    BR_TRUE("brRDEResidualZeroQ[0, 1/x, x, Log[x]]");
    BR_TRUE("brRDENoSolutionQ[{\"NoSolution\"}]");
    BR_FALSE("brRDENoSolutionQ[{\"Unknown\", \"x\"}]");
    BR_FALSE("brRDENoSolutionQ[{\"Solution\", 1}]");
    BR_EQ("brRDESolution[{\"Solution\", x^2}]", "x^2");
    BR_EQ("brRDESolution[{\"NoSolution\"}]", "None");
    BR_EQ("brRDESolution[{\"Unknown\", \"r\"}]", "None");
}

int main(void) {
    br_init();
    test_solutions();
    test_certificates();
    test_undecided();
    test_properties();
    test_wrappers();
    BR_DONE("bronstein_rde_tests");
    return 0;
}
