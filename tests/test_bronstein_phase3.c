/* test_bronstein_phase3.c -- Phase 3 of the Bronstein completion campaign
 * (BRONSTEIN_PLAN.md, "2026-11-14 -- Completion campaign"): the Shape A front
 * end and log-part gaps closed one site id at a time.
 *
 *  - Tower:cyclic-compositum: several radicals that are one simple radical in
 *    disguise, rewritten pointwise (brKummerVectors, brCompositumPlan,
 *    brCompositumRewrite, brCompositumBack) and integrated;
 *  - Tower:exp-of-log-radical: E^(r Log[a]) = a^r, and the log-linear part of
 *    an exponent split off (brSplitLogExponent: E^(x + Log[x]/2) = Sqrt[x] E^x);
 *    monomials rewritten away pruned (brPruneMonomials);
 *  - Tower:reducible-log: Log[a] = L0 + kappa with kappa locally constant,
 *    f = Sum K^j g_j (brReducibleLogSplit, brIntegrateKappaSplit,
 *    brNonzeroSomewhereQ); the C dependent-logarithm retry checked against the
 *    original integrand (Int Log[x^2]/x was Log[x]^2).
 *  - LogPart:special-genus (#9): residues at special places on a curve of
 *    positive genus -- torsion by specialization with the special places
 *    included (brSpecialSetup, brSpecialDivisorIdeal, brSpecialLocationsGoodQ,
 *    brSpecialTorsionOrder, brLowOrder; bronstein_extensions.tex E3).
 */
#include "bronstein_test_utils.h"

static void test_kummer_vectors(void) {
    BR_EQ("brKummerVectors[{{x, 1/2}, {x + 1, 1/3}}, {x}]", "{{{1/2, 0}, {0, 1/3}}, 6}");
    BR_EQ("brKummerVectors[{{x^3, 1/2}, {x, 1/2}}, {x}]", "{{{1/2}, {1/2}}, 2}");
    BR_EQ("brKummerVectors[{{2 x^2, 1/2}}, {x}]", "{{{0}}, 1}");               /* constants dropped */
    BR_EQ("brKummerVectors[{{1/x, 1/2}, {-x, 1/2}}, {x}]", "{{{1/2}, {1/2}}, 2}"); /* -1/2 = 1/2 mod 1; associates */
}

static void test_compositum_plan(void) {
    /* Z/2 x Z/3: y = Sqrt[x] (x+1)^(1/3), y^6 = x^3 (x+1)^2 */
    BR_EQ("brCompositumPlan[{{x, 2}, {x + 1, 3}}, {x}][\"N\"]", "6");
    BR_EQ("Expand[brCompositumPlan[{{x, 2}, {x + 1, 3}}, {x}][\"B\"] - x^3 (x + 1)^2]", "0");
    /* every rho_j = R_j y^k_j holds pointwise: a numeric check off the real axis */
    BR_TRUE("Module[{p = brCompositumPlan[{{x, 2}, {x + 1, 3}}, {x}], y, rho, v},"
            " rho = {Sqrt[x], (x + 1)^(1/3)}; y = Sqrt[x]^p[\"c\"][[1]] ((x + 1)^(1/3))^p[\"c\"][[2]];"
            " v = Table[N[rho[[j]] - p[\"R\"][[j]] y^p[\"k\"][[j]] /. x -> z], {j, 2}, {z, {-3 + I, -1/2, 2 - 5 I, 7}}];"
            " Max[Abs[Flatten[v]]] < 10^-10]");
    /* Z/6 from one radicand at two orders */
    BR_EQ("brCompositumPlan[{{x, 2}, {x, 3}}, {x}][\"N\"]", "6");
    /* Z/2 x Z/5 x Z/3 */
    BR_EQ("brCompositumPlan[{{x, 2}, {x + 1, 5}, {x + 2, 3}}, {x}][\"N\"]", "30");
    /* one radical only up to a locally constant sign: Sqrt[x^2 (x+1)] = +-x Sqrt[x+1] */
    BR_EQ("brCompositumPlan[{{x^2 (x + 1), 2}, {x + 1, 2}}, {x}]", "None");
    /* over a larger constant field: Sqrt[2 x] = Sqrt[2] Sqrt[x] */
    BR_EQ("brCompositumPlan[{{2 x, 2}, {x, 2}}, {x}]", "None");
    /* property: for coprime a, b the plan exists, has order a b, and its R_j are rational */
    BR_TRUE("And @@ Flatten[Table[With[{p = brCompositumPlan[{{x, a}, {x + 1, b}}, {x}]},"
            " p =!= None && p[\"N\"] === a b && FreeQ[p[\"R\"], Power[_, _Rational]]],"
            " {a, 2, 5}, {b, 2, 5}, {i, If[GCD[a, b] == 1, 1, 0]}]]");
}

static void test_compositum_rewrite(void) {
    BR_EQ("brCompositumRewrite[Sqrt[x] + Sqrt[x + 1], x]", "None");             /* not cyclic */
    BR_EQ("brCompositumRewrite[Sqrt[x], x]", "None");                          /* one radical */
    BR_EQ("brCompositumRewrite[Sqrt[x^2 (x + 1)] + Sqrt[x + 1], x]", "None");  /* sign */
    BR_EQ("brCompositumRewrite[Sqrt[x] (x + 1)^(1/3), x][\"N\"]", "6");
    BR_EQ("brCompositumRewrite[Sqrt[Log[x]] (Log[x] + 1)^(1/3), x][\"N\"]", "6");
    /* the rewrite is an identity of functions */
    BR_TRUE("Module[{rw = brCompositumRewrite[Sqrt[x] + (x + 1)^(1/3) + x^(3/2)/(x + 1)^(2/3), x]},"
            " Max[Abs[Table[N[(Sqrt[x] + (x + 1)^(1/3) + x^(3/2)/(x + 1)^(2/3)) - brCompositumBack[rw, rw[\"g\"]] /. x -> z],"
            " {z, {-3 + I, -1/2, 2 - 5 I, 7}}]]] < 10^-10]");
}

static void test_compositum_integrate(void) {
    BR_ANTIDERIV("(Sqrt[x] + (x + 1)^(1/3))/(x (x + 1))", "Bronstein");
    BR_ANTIDERIV("D[x^2 Sqrt[x] (x + 1)^(1/3), x]", "Bronstein");
    BR_ANTIDERIV("D[Sqrt[Log[x]] (Log[x] + 1)^(1/3), x]", "Bronstein");
    /* Chebyshev: x^(1/2) (1 + x)^(1/3) is not elementary; certified through
       the single-radical presentation (the same differential field) */
    BR_EQ("BronsteinIntegrate[Sqrt[x] (x + 1)^(1/3), x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[Sqrt[Log[x]] (Log[x] + 1)^(1/3)/x, x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinScope[Sqrt[x] (x + 1)^(1/3), x]", "\"InScope\"");
    /* the sign case stays in scope but unimplemented, with its site */
    BR_EQ("First[BronsteinScope[Sqrt[x^2 (x + 1)] + Sqrt[x + 1], x]]", "\"implementation incomplete\"");
}

static void test_exp_of_log(void) {
    BR_EQ("brSplitLogExponent[x, {<|\"sym\" -> t1, \"kind\" -> \"Log\", \"arg\" -> x|>}, x + t1/2]", "{x, Sqrt[x]}");
    BR_EQ("brSplitLogExponent[x, {<|\"sym\" -> t1, \"kind\" -> \"Log\", \"arg\" -> x + 1|>}, x^2 - 3 t1]", "{x^2, 1/(1 + x)^3}");
    BR_EQ("brSplitLogExponent[x, {<|\"sym\" -> t1, \"kind\" -> \"Log\", \"arg\" -> x|>}, t1/2 + 3]", "None");   /* constant rest: ExpReducible */
    BR_EQ("brSplitLogExponent[x, {<|\"sym\" -> t1, \"kind\" -> \"Log\", \"arg\" -> x|>}, x t1]", "None");       /* not a constant multiple */
    BR_EQ("brSplitLogExponent[x, {<|\"sym\" -> t1, \"kind\" -> \"Exp\", \"arg\" -> x|>}, x + t1/2]", "None");   /* only logarithms */
    /* a monomial rewritten away is pruned: Log[Log[x]] after E^(x + Log[Log[x]]/2) = Sqrt[Log[x]] E^x */
    BR_EQ("Length[brPruneMonomials[t3 Sqrt[t1], {<|\"sym\" -> t1, \"arg\" -> x, \"D\" -> 1/x|>,"
          " <|\"sym\" -> t2, \"arg\" -> t1, \"D\" -> 1/(x t1)|>, <|\"sym\" -> t3, \"arg\" -> x, \"D\" -> t3|>}]]", "2");
    BR_EQ("#[\"sym\"] & /@ brPruneMonomials[t2, {<|\"sym\" -> t1, \"arg\" -> x, \"D\" -> 1/x|>,"
          " <|\"sym\" -> t2, \"arg\" -> t1, \"D\" -> 1/(x t1)|>, <|\"sym\" -> t3, \"arg\" -> x, \"D\" -> t3|>}]", "{t1, t2}");
    BR_EQ("brPruneMonomials[x, {<|\"sym\" -> t1, \"arg\" -> x, \"D\" -> 1/x|>}]", "{}");
    /* the integrand is rewritten pointwise and decided */
    BR_ANTIDERIV("E^(x + Log[x]/2) (1 + 1/(2 x))", "Bronstein");
    BR_ANTIDERIV("D[E^(x^2 + Log[x + 1]/3), x]", "Bronstein");
    BR_ANTIDERIV("D[E^(Log[Log[x]]/2 + x) , x]", "Bronstein");
    BR_EQ("BronsteinIntegrate[E^(x + Log[x]/2)/x, x][[1]]", "\"not elementary\"");            /* E^x/Sqrt[x] */
    BR_EQ("BronsteinIntegrate[E^(x + Log[x + 1]/3) Log[x], x][[1]]", "\"not elementary\"");
    /* a radical of a logarithm's argument written as an exponential (ExpReducible path) */
    BR_EQ("Head[brCatch[brTower[E^(Log[x + Log[x]]/2), x]]]", "Association");
}

static void test_reducible_log(void) {
    BR_EQ("brReducibleLogSplit[Sqrt[x] Log[x], x]", "None");                  /* nothing reducible */
    BR_EQ("Length[brReducibleLogSplit[Log[x^2] Log[x] Sqrt[Log[x]], x][\"terms\"]]", "2");
    BR_EQ("brReducibleLogSplit[Log[x^2] Log[x] Sqrt[Log[x]], x][\"kappas\"][[1, 2]]", "-2 Log[x] + Log[x^2]");
    BR_EQ("brReducibleLogSplit[Log[x]/(1 + Log[x^2]), x]", "None");          /* K in a denominator: a parameter */
    BR_EQ("brReducibleLogSplit[Log[x] Log[Log[x^2]], x]", "None");          /* a monomial above depends on it */
    BR_EQ("brReducibleLogSplit[Log[x] E^(Log[x^2]^2), x]", "None");
    /* review finding: PossibleZeroQ[Log[E^x] - x] is True, so a zero test may
       never stand for "kappa = 0"; the answer carries kappa and is valid everywhere */
    BR_EQ("BronsteinIntegrate[(Log[E^x] - x) Sqrt[Log[x]] + 1/x, x][[1]]", "\"not elementary\"");
    BR_TRUE("brNonzeroSomewhereQ[Log[E^x] - x, x]");
    BR_TRUE("With[{F = BronsteinIntegrate[(Log[E^x] - x) Sqrt[Log[x]]/x + 1/x, x]}, !ListQ[F] && !FreeQ[F, Log[E^x]]]");
    BR_EQ("BronsteinIntegrate[(Log[x^2] - 2 Log[x]) E^(x^2), x][[1]]", "\"not elementary\"");
    /* the kappa mode of the tower builder: K for kappa, recorded in brKappas */
    BR_FALSE("brKappaMode");
    BR_EQ("brKappas", "{}");
    BR_EQ("Block[{brKappaMode = True, brKappas = {}}, brKernelsToSymbols[Log[x^2] Log[x], x]; Length[brKappas]]", "1");
    BR_EQ("Head[brCatch[brKernelsToSymbols[Log[x^2] Log[x], x]]]", "tragerFail");   /* outside the mode */
    /* brIntegrateKappaSplit directly */
    BR_TRUE("With[{F = brIntegrateKappaSplit[Log[x^2] Sqrt[Log[x]]/x, x, brReducibleLogSplit[Log[x^2] Sqrt[Log[x]]/x, x]]},"
            " PossibleZeroQ[D[F, x] - Log[x^2] Sqrt[Log[x]]/x]]");
    BR_EQ("brFailTag[brCatch[brIntegrateKappaSplit[(Log[x^2] - 2 Log[x]) Sqrt[Log[x]], x,"
          " brReducibleLogSplit[(Log[x^2] - 2 Log[x]) Sqrt[Log[x]], x]]]][[1]]", "\"NonElementary\"");
    BR_TRUE("brNonzeroSomewhereQ[Log[x^2] - 2 Log[x], x]");
    BR_FALSE("brNonzeroSomewhereQ[Log[2 x] - Log[2] - Log[x], x]");
    /* answers: the branch difference is carried explicitly */
    BR_ANTIDERIV("Log[x^2] Sqrt[Log[x]]/x", "Bronstein");
    BR_ANTIDERIV("Log[x (x + 1)] Sqrt[Log[x] + Log[x + 1]] (1/x + 1/(x + 1))", "Bronstein");
    BR_ANTIDERIV("Log[2 x] Sqrt[x]", "Bronstein");
    BR_ANTIDERIV("(Log[x^2] - 2 Log[x]) Sqrt[Log[x]]/x", "Bronstein");
    /* certificates: g_0 = 2 Log[x] Sqrt[Log[x]] is not elementary */
    BR_EQ("BronsteinIntegrate[Log[x^2] Sqrt[Log[x]] Log[x], x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[Sqrt[Log[x]] Log[x^2], x][[1]]", "\"not elementary\"");
    /* only g_1 = Sqrt[Log[x]] is not elementary: kappa != 0 somewhere, transcendental there */
    BR_EQ("BronsteinIntegrate[(Log[x^2] - 2 Log[x]) Sqrt[Log[x]], x][[1]]", "\"not elementary\"");
    /* the C dependent-logarithm retry no longer answers Log[x]^2 for Int Log[x^2]/x */
    BR_TRUE("With[{r = Integrate`RischTranscendental[Log[x^2]/x, x]}, Head[r] === Integrate`RischTranscendental || PossibleZeroQ[D[r, x] - Log[x^2]/x]]");
    BR_EQ("Integrate`RischTranscendental[Log[2 x]/x, x]", "1/2 (2 Log[2] Log[x] + Log[x]^2)");
    BR_TRUE("PossibleZeroQ[D[Integrate[Log[x^2]/x, x], x] - Log[x^2]/x]");
}

static void test_special_genus(void) {
    BR_EQ("brLowOrder[t^2 + t^3, t]", "2");
    BR_EQ("brLowOrder[1 + t, t]", "0");
    BR_EQ("brLowOrder[0, t]", "Infinity");
    /* the component of D Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)] on y^3 = theta^3 + theta + 1 (genus 1) */
    br_run("T9 = brTower[Together[D[Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)], x]], x]; cf9 = brCanonicalForm[T9];"
           " cc9 = brComponentCurve[T9, cf9[\"n\"], cf9[\"g\"], 1]; cv9 = brCurve[T9, cc9[\"m\"], cc9[\"G\"]];"
           " c9 = brSpecialReduce[T9, cv9, brHermiteReduce[T9, cv9, Together[cf9[\"coeffs\"][[2]] cc9[\"s\"]]][\"Remainder\"]][\"Remainder\"];"
           " ss9 = brSpecialSetup[T9, cv9, c9, brNormalResidues[T9, cv9, c9], brLogInfinityResidues[T9, cv9, c9]]; True");
    BR_EQ("ss9[\"genus\"]", "1");
    BR_EQ("Length[ss9[\"finRoots\"]]", "3");
    BR_TRUE("brSpecialLocationsGoodQ[T9, cv9, ss9[\"rsF\"], ss9[\"locs\"], {x -> 3}]");
    BR_TRUE("KeyExistsQ[brSpecialDivisorIdeal[ss9[\"Fk\"], ss9[\"fr\"], ss9[\"finRoots\"], ss9[\"locs\"], ss9[\"rhoLists\"],"
            " ss9[\"qb\"][\"coeffs\"][[All, 1]]], \"M\"]");
    BR_EQ("brSpecialTorsionOrder[T9, cv9, ss9[\"rsF\"], ss9[\"fr\"], ss9[\"locs\"], ss9[\"rhoLists\"], ss9[\"qb\"][\"coeffs\"][[All, 1]]]", "1");
    /* a frame-degenerate special location is caught in characteristic 0 */
    BR_FALSE("brSpecialLocationsGoodQ[T9, cv9, ss9[\"rsF\"], {<|\"where\" -> \"inf\", \"u\" -> 1/x|>}, {x -> 0}]");
    /* review finding (2026-11-14): x -> 64/5 puts the residue place theta = x - 59/5
       on the frame point theta = a = 1, so the divisor lost a place and a false
       "orders disagree" certificate came out; the frame degrees now gate it */
    br_run("Rv = Together[D[Log[(Log[x]^3 + Log[x] - x + 59/5)^(1/3) - Log[x]], x]]; Tv = brTower[Rv, x]; cfv = brCanonicalForm[Tv];"
           " ccv = brComponentCurve[Tv, cfv[\"n\"], cfv[\"g\"], 1]; cvv = brCurve[Tv, ccv[\"m\"], ccv[\"G\"]];"
           " cv0 = brSpecialReduce[Tv, cvv, brHermiteReduce[Tv, cvv, Together[cfv[\"coeffs\"][[2]] ccv[\"s\"]]][\"Remainder\"]][\"Remainder\"];"
           " ssv = brSpecialSetup[Tv, cvv, cv0, brNormalResidues[Tv, cvv, cv0], brLogInfinityResidues[Tv, cvv, cv0]]; True");
    BR_FALSE("brFrameDegreesKeptQ[ssv[\"fr\"], {x -> 64/5}]");
    BR_TRUE("brFrameDegreesKeptQ[ssv[\"fr\"], {x -> 3}]");
    BR_EQ("Table[brSpecialTorsionOrder[Tv, cvv, ssv[\"rsF\"], ssv[\"fr\"], ssv[\"locs\"], ssv[\"rhoLists\"], ssv[\"qb\"][\"coeffs\"][[All, j]]],"
          " {j, Length[ssv[\"qb\"][\"basis\"]]}]", "{1, 1}");
    BR_EQ("brFrameIdealDegree[idealUnit[ss9[\"fr\"][\"Bu\"]], ss9[\"fr\"][\"u\"], 3]", "0");
    BR_EQ("brFrameIdealDegree[brSpecialDivisorIdeal[ss9[\"Fk\"], ss9[\"fr\"], ss9[\"finRoots\"], ss9[\"locs\"], ss9[\"rhoLists\"],"
          " ss9[\"qb\"][\"coeffs\"][[All, 1]]], ss9[\"fr\"][\"u\"], 3]", "0");
    /* integrals that stopped at LogPart:special-genus: elementary (Chebyshev, derivatives) */
    BR_ANTIDERIV("1/(x (1 + Log[x]^3)^(1/3))", "Bronstein");
    BR_ANTIDERIV("E^x/(1 + E^(3 x))^(1/3)", "Bronstein");
    BR_ANTIDERIV("1/(x (2 + Log[x]^3)^(1/3))", "Bronstein");
    BR_ANTIDERIV("Together[D[Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[x + Log[-E^x + (1 + E^(3 x))^(1/3)], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[(1 + Log[x] + Log[x]^3)^(1/3) + Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)], x]]", "Bronstein");
    /* and a certificate: the residue divisor has infinite order */
    BR_EQ("BronsteinIntegrate[(1 + Log[x])/(x (1 + Log[x] + Log[x]^3)^(1/3)), x][[1]]", "\"not elementary\"");
}

int main(void) {
    br_init();
    test_kummer_vectors();
    test_compositum_plan();
    test_compositum_rewrite();
    test_compositum_integrate();
    test_exp_of_log();
    test_reducible_log();
    test_special_genus();
    BR_DONE("bronstein_phase3_tests");
    return 0;
}
