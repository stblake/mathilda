/* test_bronstein_trig.c -- Phase 4.1b and 4.2 of the Bronstein completion
 * campaign: I modulo p, Gaussian residue polynomials and kappa, and the
 * trigonometric front end (split from test_bronstein_phase4.c: the integrals
 * here are heavier than one suite's 120 s alarm allows together).
 */
#include "bronstein_test_utils.h"

static void test_qi_41b(void) {
    /* a residue polynomial with Gaussian coefficients: Root objects of its norm */
    BR_EQ("Length[gaussianRootObjects[z^2 - I, z]]", "2");
    BR_TRUE("AllTrue[gaussianRootObjects[z^2 - I, z], Abs[N[#^2 - I, 30]] < 10^-20 &]");
    BR_TRUE("Block[{tragerConstantI = True}, KeyExistsQ[residueField[z^2 - I, z, tt][\"F\"], \"I\"]]");
    /* torsion over Q(i) (genus 1): derivatives of logarithms are elementary */
    BR_ANTIDERIV("Together[D[Log[Sqrt[x^3 + 1] + I x], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[Sqrt[x^3 + I] + Log[Sqrt[x^3 + I] + x + I], x]]", "Bronstein");
    BR_EQ("BronsteinIntegrate[1/((x - I) Sqrt[x^3 + 1]), x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[Sqrt[x^3 + I], x][[1]]", "\"not elementary\"");
    /* a Gaussian kappa (special places), the review cases of v0.6551, now answered */
    BR_ANTIDERIV("Sqrt[E^x + I]", "Bronstein");
    BR_ANTIDERIV("(E^x + I)^(1/3)", "Bronstein");
    BR_ANTIDERIV("E^x/Sqrt[1 + I E^(2 x)]", "Bronstein");
    /* the identity check over a residue field containing i */
    BR_ANTIDERIV("Together[D[Log[Log[x] + (Log[x]^3 + I)^(1/3)], x]]", "Bronstein");
}

static void test_trig_front_end(void) {
    BR_TRUE("brTrigFrontEnd[Sqrt[x]] === Sqrt[x]");
    BR_TRUE("FreeQ[brTrigFrontEnd[Sqrt[Tan[x]]], Tan] && !FreeQ[brTrigFrontEnd[Sqrt[Tan[x]]], E]");
    BR_TRUE("brTrigPresentation[Sin[x], -(E^(I x) + E^(-I x))/2, x] === -Cos[x]");
    BR_EQ("BronsteinScope[Sqrt[Tan[x]], x]", "\"InScope\"");
    BR_ANTIDERIV("Sqrt[Tan[x]]", "Bronstein");
    BR_ANTIDERIV("Tan[x]^(1/3)", "Bronstein");
    BR_ANTIDERIV("Sin[x] Sqrt[Cos[x]]", "Bronstein");
    BR_ANTIDERIV("Sqrt[1 + Sin[x]]", "Bronstein");
    BR_ANTIDERIV("ArcSin[x]/Sqrt[1 - x^2]", "Bronstein");
    BR_ANTIDERIV("x ArcTan[x]/Sqrt[1 + x^2]", "Bronstein");
    BR_ANTIDERIV("Sinh[x] Sqrt[Cosh[x]]", "Bronstein");
    /* classical non-elementary: elliptic, and Sec^(1/3) */
    BR_EQ("BronsteinIntegrate[Sqrt[Sin[x]], x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[Sqrt[Cos[x]], x][[1]]", "\"not elementary\"");
    BR_EQ("BronsteinIntegrate[Sec[x]^(1/3), x][[1]]", "\"not elementary\"");
}

static void test_qi_review_fixes(void) {
    /* review of E4.1 (2026-11-14): kappa = I was read as non-constant
       (FreeQ[I, _Symbol] is False), the special residue taken as 0, and these
       elementary integrals certified NOT elementary.  Now refused. */
    BR_EQ("brGammaForm[x^2, 2]", "{1, x}");
    BR_EQ("brFailTag[brCatch[brGammaForm[I x^2, 2]]][[1]]", "\"Unimplemented\"");   /* outside Q(i) */
    BR_EQ("Block[{tragerConstantI = True}, brGammaForm[I x^2, 2]]", "{I, x}");             /* inside Q(i) (4.1b) */
    BR_EQ("brGammaForm[x, 2]", "None");
    BR_FALSE("MatchQ[BronsteinIntegrate[Sqrt[E^x + I], x], {\"not elementary\", _}]");
    BR_FALSE("MatchQ[BronsteinIntegrate[(E^x + I)^(1/3), x], {\"not elementary\", _}]");
    BR_FALSE("MatchQ[BronsteinIntegrate[Sqrt[1 + I Log[x]^2]/x, x], {\"not elementary\", _}]");
    BR_FALSE("MatchQ[BronsteinIntegrate[E^x/Sqrt[1 + I E^(2 x)], x], {\"not elementary\", _}]");
    /* a RootSum answer over Q(i) keeps no private symbol (complexSplit's ii leaked) */
    BR_TRUE("With[{F = BronsteinIntegrate[Sqrt[x]/(x^2 - I), x]}, !ListQ[F] &&"
            " FreeQ[F, s_Symbol /; StringContainsQ[Context[s], \"Private\"]] && PossibleZeroQ[D[F, x] - Sqrt[x]/(x^2 - I)]]");
    BR_ANTIDERIV("Sqrt[x]/(x^4 + 1) + I/x", "Bronstein");
    /* gfMod: I -> the least square root of -1 mod p = 1 (4); refused otherwise */
    BR_EQ("gfSqrtMinusOne[5]", "2");
    BR_EQ("gfSqrtMinusOne[13]", "5");
    BR_EQ("gfMod[x + I, 5]", "2 + x");
    BR_EQ("gfMod[(x + I) (x - I) - x^2 - 1, 13]", "0");                        /* a ring homomorphism */
    BR_THROWS("gfMod[x + I, 7]", "Unimplemented", "reduction of complex constants modulo a prime p != 1 mod 4");
    BR_EQ("gfMod[x/2 + 3, 5]", "3 + 3 x");
    BR_TRUE("denominatorsPrimeToQ[x + I/3, 5] && !denominatorsPrimeToQ[x + I/5, 5]");
}

static void test_trig_review(void) {
    /* review of E4.1b/E4.2 (2026-11-14): no false certificate; these fixes */
    /* residues all in Q(i) (fieldQ shortcut treated Q(i) as one-dimensional) */
    BR_TRUE("With[{F = BronsteinIntegrate[Sec[x] Sqrt[Sin[x]], x]}, !ListQ[F] &&"
            " PossibleZeroQ[D[F, x] - Sec[x] Sqrt[Sin[x]], Assumptions -> Element[x, Complexes]]]");
    BR_TRUE("With[{F = BronsteinIntegrate[Sqrt[Sin[x]] Tan[x], x]}, !ListQ[F]]");
    /* hyperbolic answers stay exponential: Cosh - Sinh cancels at machine precision */
    BR_TRUE("With[{F = BronsteinIntegrate[Sqrt[Tanh[x]], x]}, !ListQ[F] && FreeQ[F, Sinh | Cosh] && NumericQ[N[(F /. x -> 12) - (F /. x -> 10)]]]");
    /* a Gaussian resultant left as uncancelled fractions by Expand */
    BR_ANTIDERIV("Together[D[Log[1 + E^x Sqrt[E^(3 x) + I]], x]]", "Bronstein");
}

int main(void) {
    br_init();
    test_qi_review_fixes();
    test_qi_41b();
    test_trig_front_end();
    test_trig_review();
    BR_DONE("bronstein_trig_tests");
    return 0;
}
