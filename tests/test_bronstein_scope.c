/* test_bronstein_scope.c -- batch D of the gap register (scope gaps closed in
 * the Bronstein completion campaign, tests/bronstein_paper/open_gaps.m):
 *
 *  - G-trig-1 / G-nf-11 (E6d.1): square roots equal to one radical up to a
 *    locally constant sign, merged in the front end (brSignMerge) and
 *    integrated by the sign split, then presented over the user's radicals
 *    (brSignPresent): ArcCosh[x]/Sqrt[x^2 - 1], Sqrt[-2 x] + Sqrt[x];
 *  - E6d-5 (G-sign-1 stays open): a sign seen only on a cut never certifies
 *    (brSignPatternSeenQ);
 *  - G-red-1 (E6a.3): several reducible logarithms whose branch differences are
 *    rational multiples of 2 Pi I (brKappaUnitOrder, brKappaRegionCertificate);
 *  - G-qi-1: integrals in k with Gaussian coefficients and ArcTan over Q(i)
 *    (brIntegrateInK, brSplitLogs);
 *  - G-adv-t4: a rational Q(i) integrand reaches Trager's rational integrator
 *    (algebraicNumberQ);
 *  - G-core-4: the C limited integration folds a dependent logarithm by its
 *    derivative (Int Log[x^2]/x = Log[x^2]^2/4).
 */
#include "bronstein_test_utils.h"

static void test_sign_merge(void) {
    /* (a) a product of square roots inside one Times */
    BR_EQ("Length[brSignMerge[Log[x + Sqrt[x - 1] Sqrt[x + 1]]/Sqrt[x^2 - 1], x][\"signs\"]]", "1");
    BR_TRUE("Module[{sm = brSignMerge[Log[x + Sqrt[x - 1] Sqrt[x + 1]]/Sqrt[x^2 - 1], x]},"
            " Length[DeleteDuplicates[Cases[sm[\"g\"], Power[b_, 1/2 | -1/2] :> Together[b], {0, Infinity}]]] === 1]");
    /* (b) a negative rational multiple of another radicand: the constant field gains I */
    BR_TRUE("Module[{sm = brSignMerge[Sqrt[-x] + Sqrt[x], x]}, !FreeQ[sm[\"g\"], Complex] && Length[sm[\"signs\"]] === 1]");
    /* the rewrite is an identity once the sign symbol is its value */
    BR_TRUE("Module[{f = Sqrt[-x] Log[x] + Sqrt[x], sm},"
            " sm = brSignMerge[f, x];"
            " Max[Abs[Table[N[(f - (sm[\"g\"] /. (Rule @@@ sm[\"signs\"]))) /. x -> z, 30], {z, {-3 + I, -1/2, 2 - 5 I, 7, -4 - I/3}}]]] < 10^-20]");
    BR_TRUE("Module[{f = 1/(Sqrt[x] Sqrt[x + 1]) + Sqrt[x^2 + x], sm},"
            " sm = brSignMerge[f, x];"
            " Max[Abs[Table[N[(f - (sm[\"g\"] /. (Rule @@@ sm[\"signs\"]))) /. x -> z, 30], {z, {-3 + I, -1/2, 2 - 5 I, 7, -4 - I/3}}]]] < 10^-20]");
    /* no merge: one radical, a non-cyclic sum, positive contents (already one radical) */
    BR_EQ("brSignMerge[Sqrt[x], x]", "None");
    BR_EQ("brSignMerge[Sqrt[x] + Sqrt[x + 1], x]", "None");
    BR_EQ("brSignMerge[Sqrt[x] + Sqrt[x + 1] + Sqrt[x] Sqrt[x + 1], x]", "None");
}

static void test_sign_merge_integrate(void) {
    BR_ANTIDERIV("ArcCosh[x]/Sqrt[x^2 - 1]", "Bronstein");                          /* G-trig-1 */
    BR_ANTIDERIV("Sqrt[-2 x] + Sqrt[x]", "Bronstein");                              /* G-nf-11 */
    BR_ANTIDERIV("1/(Sqrt[x] Sqrt[x + 1])", "Bronstein");
    BR_ANTIDERIV("Sqrt[x - 1] Sqrt[x + 1]", "Bronstein");
    BR_ANTIDERIV("Sqrt[-x] Log[x]", "Bronstein");
    BR_ANTIDERIV("E^x Sqrt[-x] Sqrt[x]", "Bronstein");
    BR_EQ("BronsteinScope[ArcCosh[x]/Sqrt[x^2 - 1], x]", "\"InScope\"");
    BR_EQ("BronsteinScope[Sqrt[-2 x] + Sqrt[x], x]", "\"InScope\"");
    /* Sqrt[t^2 + t], t = Log[x], up to a sign: like Sqrt[1 + Log[x]^2], not elementary */
    BR_EQ("BronsteinIntegrate[Sqrt[Log[x]] Sqrt[Log[x] + 1], x][[1]]", "\"not elementary\"");
    /* f = eps E^(x^2): both sign patterns are non-elementary and seen on a
       region, so the sign split certifies f */
    BR_EQ("BronsteinIntegrate[E^(x^2) Sqrt[x - 1] Sqrt[x + 1]/Sqrt[x^2 - 1], x][[1]]", "\"not elementary\"");
}

static void test_sign_present(void) {
    /* presented over the user's radicals: eps ArcCosh-form, no sign indicators */
    BR_TRUE("FreeQ[BronsteinIntegrate[ArcCosh[x]/Sqrt[x^2 - 1], x], Plus[1, Times[-1, _]]]");
    BR_TRUE("BronsteinIntegrate[1/(Sqrt[x] Sqrt[x + 1]), x] === Log[1 + 2 x + 2 Sqrt[x] Sqrt[1 + x]]");
    /* brSignPresent keeps F when there is nothing to write back */
    BR_EQ("brSignPresent[1, F0, {}, <|\"signs\" -> {}, \"back\" -> {}|>, x]", "F0");
}

static void test_sign_cut_only(void) {
    /* G-sign-1 / E6d-5: eps = -1 only on the cut x < -1, never on an open set:
       no certificate (the real interval on a cut is not a region) */
    BR_EQ("brSignPatternSeenQ[{Sqrt[x + 1] Sqrt[1/(x + 1)]}, {-1}, x]", "False");
    BR_TRUE("brSignPatternSeenQ[{Sqrt[x + 1] Sqrt[1/(x + 1)]}, {1}, x]");
    BR_EQ("brSignPatternSeenQ[{Sqrt[x^2] /x}, {2}, x]", "False");
    BR_EQ("BronsteinIntegrate[(1 + E^(x^2))/2 + (1 - E^(x^2))/2 Sqrt[x + 1] Sqrt[1/(x + 1)], x][[1]]",
          "\"implementation incomplete\"");
}

static void test_kappa_unit_order(void) {
    BR_EQ("brKappaUnitOrder[Log[x^2] - 2 Log[x]]", "1");
    BR_EQ("brKappaUnitOrder[Log[x^3] - 3 Log[x]]", "1");
    BR_EQ("brKappaUnitOrder[Log[-x] - Log[x]]", "2");                 /* E^kappa = -1 */
    BR_EQ("brKappaUnitOrder[Log[x]/2 - Log[Sqrt[x]]]", "2");
    BR_EQ("brKappaUnitOrder[Log[2 x] - Log[x]]", "None");            /* Log[2] + 2 Pi I m */
    BR_EQ("brKappaUnitOrder[Log[E^x] - x]", "1");
    /* property: M kappa is in 2 Pi I Z at random points */
    BR_TRUE("And @@ Flatten[Table[With[{k = Log[x^n] - n Log[x], M = brKappaUnitOrder[Log[x^n] - n Log[x]]},"
            " Table[With[{v = N[M k /. x -> z, 30]/(2 Pi I)}, Abs[v - Round[v]] < 10^-20],"
            " {z, {-3 + I, -1/2 - 4 I, 2 - 5 I, -7 + I/9}}]], {n, 2, 6}]]");
}

static void test_kappa_region(void) {
    /* G-red-1 */
    BR_EQ("BronsteinIntegrate[(Log[x^2] - 2 Log[x]) (Log[x^3] - 3 Log[x]) Sqrt[Log[x]], x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[(Log[x^2] - 2 Log[x]) Sqrt[Log[x]] + (Log[x^3] - 3 Log[x]) Sqrt[Log[x]], x][[1]]", "\"not elementary\"");
    BR_ANTIDERIV("(Log[x^2] - 2 Log[x])/x + (Log[x^3] - 3 Log[x]) Sqrt[Log[x]]/x", "Bronstein");
    BR_TRUE("Module[{rl = brReducibleLogSplit[(Log[x^2] - 2 Log[x]) (Log[x^3] - 3 Log[x]) Sqrt[Log[x]], x]},"
            " MatchQ[brKappaRegionCertificate[rl[\"kappas\"], rl[\"terms\"], x], {{_, _}, 2}]]");
    /* every coefficient elementary: no certificate */
    BR_TRUE("Module[{rl = brReducibleLogSplit[(Log[x^2] - 2 Log[x]) (Log[x^3] - 3 Log[x])/x, x]},"
            " brKappaRegionCertificate[rl[\"kappas\"], rl[\"terms\"], x] === None]");
}

static void test_qi_limited(void) {
    /* G-qi-1 */
    BR_EQ("BronsteinIntegrate[Sqrt[I + Log[x]^2], x][[1]]", "\"not elementary\"");
    BR_ANTIDERIV("Log[x]/(x Sqrt[I + Log[x]^2])", "Bronstein");
    BR_ANTIDERIV("Sqrt[I + Log[x]]/x", "Bronstein");
    BR_TRUE("Block[{tragerConstantI = True}, brSplitLogs[<|\"syms\" -> {}, \"x\" -> x|>, I/2 x + (2 + I) Log[x + I]] =!= None]");
    BR_EQ("Block[{tragerConstantI = False}, brSplitLogs[<|\"syms\" -> {}, \"x\" -> x|>, x + I Log[x]]]", "None");
}

static void test_trager_qi_rational(void) {
    /* G-adv-t4 */
    BR_TRUE("algebraicNumberQ[1 + 2 I]");
    BR_TRUE("algebraicNumberQ[Sqrt[2] + I/3]");
    BR_EQ("algebraicNumberQ[1.5 I]", "False");
    BR_TRUE("algebraicRationalFunctionQ[1/(x^2 + I), x]");
    BR_ANTIDERIV("1/(x^2 + I)", "Trager");
    BR_ANTIDERIV("x/(x^2 + I)^2", "Trager");
}

static void test_c_dependent_log(void) {
    /* G-core-4 */
    BR_TRUE("Integrate`RischTranscendental[Log[x^2]/x, x] === Log[x^2]^2/4");
    BR_TRUE("Integrate`RischTranscendental[Log[x^2]^2/x, x] === Log[x^2]^3/6");
    BR_TRUE("Integrate`RischTranscendental[Log[x^3]/x, x] === Log[x^3]^2/6");
    BR_TRUE("With[{F = Integrate`RischTranscendental[Log[(x + 1)^2]/(x + 1), x]},"
            " PossibleZeroQ[D[F, x] - Log[(x + 1)^2]/(x + 1)]]");
}

/* Adversarial review of batch D (2026-11-15): false certificates, all on
 * integrands identically 0 or 1 off the cuts, and one wrong C answer. */
static const char* review_not_certified[] = {
    /* class A1: a merged radicand with a square factor kept a hidden sign */
    "E^(x^2) (Sqrt[x] Sqrt[(1 + x)/x]/Sqrt[1 + x] - 1)",
    "(1 + E^(x^2))/2 + (1 - E^(x^2))/2 Sqrt[x] Sqrt[(x + 1)/x]/Sqrt[x + 1]",
    "(1 + E^(x^2))/2 + (1 - E^(x^2))/2 Sqrt[(x - 1)/(x + 1)] Sqrt[x + 1]/Sqrt[x - 1]",
    /* class A2: a region-local inner certificate lifted through the sign split */
    "(Sqrt[x^2 - 1] + Sqrt[x - 1] Sqrt[x + 1]) (Log[x^2] - 2 Log[x]) E^(x^2)",
    "(1 + Sqrt[x - 1] Sqrt[x + 1]/Sqrt[x^2 - 1]) (Log[x^2] - 2 Log[x]) Sqrt[Log[x]]",
    "(1 + Sqrt[x - 1] Sqrt[x + 1]/Sqrt[x^2 - 1]) (Log[x^2] - 2 Log[x]) (Log[x^3] - 3 Log[x]) Sqrt[Log[x]]",
    "(Sqrt[x^2 + 2 x] + Sqrt[x] Sqrt[x + 2]) (Log[(x + 1)^2] - 2 Log[x + 1]) E^(x^2)",
    "(Sqrt[x] + I Sqrt[-x]) (Log[-x] - Log[x] + I Pi) E^(x^2)",
    /* class B: a transcendental constant in a kappa coefficient */
    "((Log[-x] - Log[x])^2 + Pi^2) Sqrt[Log[x]]",
    "((Log[-x] - Log[x])^2 + Pi^2) (1 + Log[x^2] - 2 Log[x]) Sqrt[Log[x]]",
    "(Log[x^2] - 2 Log[x]) (Log[x^2] - 2 Log[x] + 2 Pi I) (Log[x^2] - 2 Log[x] - 2 Pi I) Sqrt[Log[x]]",
    NULL
};

static void test_review_regressions(void) {
    char buf[1024];
    for (int i = 0; review_not_certified[i]; i++) {
        snprintf(buf, sizeof buf, "First[BronsteinIntegrate[%s, x]] =!= \"not elementary\"", review_not_certified[i]);
        BR_TRUE(buf);
    }
    BR_TRUE("brAlgebraicConstantsOnlyQ[Sqrt[2] E^x Sqrt[Log[x]], x]");
    BR_EQ("brAlgebraicConstantsOnlyQ[Pi^2 Sqrt[Log[x]], x]", "False");
    BR_EQ("brAlgebraicConstantsOnlyQ[Log[2] Sqrt[Log[x]], x]", "False");
    BR_EQ("brAlgebraicConstantsOnlyQ[E Sqrt[Log[x]], x]", "False");
    BR_EQ("brInRegion", "False");
    /* the review's second round: power regained */
    BR_TRUE("brSquareFreeSqrt[x^2 (1 + x)^2, x] === x (1 + x)");
    BR_TRUE("With[{c = brSquareFreeSqrt[(x^2 - 1)^3/x^2, x]}, Together[c^2 - (x^2 - 1)^3/x^2] === 0]");
    BR_ANTIDERIV("Sqrt[-1 - x] Sqrt[-x] Sqrt[x] Sqrt[1 + x]", "Bronstein");
    BR_ANTIDERIV("Sqrt[I x] Sqrt[x]", "Bronstein");
    /* a constant I Pi is (2 Pi I)/2 in the grading of E6a.3, not a refusal */
    BR_EQ("BronsteinIntegrate[(Log[-x] - Log[x]) (Log[x^2] - 2 Log[x]) Sqrt[Log[x]], x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[E^(x^2) (2 Pi I + Log[x^2] - 2 Log[x]) (Log[x^3] - 3 Log[x]), x][[1]]", "\"not elementary\"");
    BR_TRUE("First[BronsteinIntegrate[(Log[x^2] - 2 Log[x]) (Log[x^2] - 2 Log[x] + 2 Pi I) (Log[x^2] - 2 Log[x] - 2 Pi I) (Log[x^3] - 3 Log[x] + 1) Sqrt[Log[x]], x]] =!= \"not elementary\"");
    /* class C: the C dilogarithm path (PowerExpand check) answered wrongly for x < 0 */
    BR_TRUE("With[{r = Integrate`RischTranscendental[Log[x^2]/(1 + x), x]},"
            " Head[r] === Integrate`RischTranscendental || PossibleZeroQ[D[r, x] - Log[x^2]/(1 + x)]]");
    BR_TRUE("PossibleZeroQ[D[Integrate`RischTranscendental[Log[x]/(1 + x), x], x] - Log[x]/(1 + x)]");
}

int main(void) {
    br_init();
    test_sign_merge();
    test_sign_merge_integrate();
    test_sign_present();
    test_sign_cut_only();
    test_kappa_unit_order();
    test_kappa_region();
    test_qi_limited();
    test_trager_qi_rational();
    test_c_dependent_log();
    test_review_regressions();
    BR_DONE("bronstein_scope_tests");
    return 0;
}
