/* test_bronstein_nf.c -- Phase 4.3 of the Bronstein completion campaign:
 * number-field constants (radicals of rationals such as Sqrt[2], 2^(1/3)).
 *
 *  - the constant field from the integrand (tragerConstantsOf, brAlgConstantsOf,
 *    tragerExtension, tragerFactorList over Q(constants), constantFieldNumberQ);
 *  - residue fields that contain the constants (residueField "consts",
 *    fieldSubI / fieldConstRep, norms of algebraic coefficients in
 *    gaussianRootObjects);
 *  - a radical of a constant inside a radicand is not nesting;
 *  - certificates over a number field are demoted (nfDemoteCertificate,
 *    NF:certificate) until the base-field arithmetic is field-aware (4.3b);
 *    a constant factor is pulled out by linearity, so Sqrt[2] Sqrt[Log[x]] is
 *    decided over Q;
 *  - reduction mod p refuses radical constants (TorsionModP:algebraic-constant);
 *  - regression: a nested Return in residueField dropped the "I" key and gave a
 *    false certificate on Sqrt[Tan[x]] (found while adding 4.3).
 */
#include "bronstein_test_utils.h"

static void test_nf_primitives(void) {
    BR_EQ("Sort[tragerConstantsOf[Sqrt[x + Sqrt[2]] + 2^(1/3) x]]", "{2^(1/3), Sqrt[2]}");
    BR_EQ("tragerConstantsOf[Sqrt[x] + I]", "{}");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, tragerExtension[]]", "{Sqrt[2]}");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}, tragerConstantI = True}, tragerExtension[]]", "{I, Sqrt[2]}");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, Length[tragerFactorList[z^2 - 2]]]", "3");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, constantFieldNumberQ[1 + Sqrt[2]/3]]");
    BR_FALSE("Block[{tragerAlgConstants = {Sqrt[2]}}, constantFieldNumberQ[Pi]]");
    BR_FALSE("constantFieldNumberQ[Sqrt[2]]");                                  /* outside the mode */
    br_run("rfN = Block[{tragerAlgConstants = {Sqrt[2]}}, residueField[z^2 - 3, z, tt]]; True");
    BR_TRUE("KeyExistsQ[rfN[\"F\"], \"consts\"]");
    BR_EQ("kc[rfN[\"F\"], Sqrt[2]^2 - 2]", "0");
    BR_EQ("kc[rfN[\"F\"], fieldConstRep[rfN[\"F\"], Sqrt[2]]^2 - 2]", "0");
    BR_EQ("kc[rfN[\"F\"], fieldSubI[rfN[\"F\"], Sqrt[6]]^2 - 6]", "0");          /* a product the arithmetic formed */
    BR_EQ("Length[gaussianRootObjects[z^2 - Sqrt[2], z]]", "2");
    /* outside a declared constant field a radical constant is refused */
    BR_THROWS("gfMod[x + Sqrt[2], 7]", "Unimplemented",
              "reduction mod p of a radical constant outside the declared constant field");
    /* Phase 4.3b: a degree-1 prime of Q(Sqrt[2]) above 7 (2 = 3^2 mod 7), none above 5 */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{r = gfMod[Sqrt[2], 7]}, Mod[r^2 - 2, 7] === 0]]");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, gfConstRoot[5]]", "None");
    /* one homomorphism: the images of Sqrt[2], Sqrt[3] and Sqrt[6] multiply */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2], Sqrt[3]}}, Mod[gfMod[Sqrt[2], 23] gfMod[Sqrt[3], 23] - gfMod[Sqrt[6], 23], 23] === 0]");
    BR_TRUE("Block[{tragerAlgConstants = {2^(1/3)}}, With[{r = gfMod[2^(1/3), 31]}, Mod[r^3 - 2, 31] === 0 && gfMod[2^(2/3), 31] === Mod[r^2, 31]]]");
    /* I and a radical together, under the same map */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}, tragerConstantI = True}, With[{i = gfImageI[17], s = gfMod[Sqrt[2], 17]}, Mod[i^2 + 1, 17] === 0 && Mod[s^2 - 2, 17] === 0]]");
    BR_EQ("gfModQ[x/2 + 3, 7]", "3 + 4 x");
    BR_EQ("Keys[gfConstField[{Sqrt[2]}]]", "{\"g\", \"mu\", \"theta\", \"reps\"}");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, Mod[gfConstImage[Sqrt[2], 7]^2, 7] === 2]");
    BR_EQ("brAlgConstantsOf[Sqrt[E^x + Sqrt[3]]]", "{Sqrt[3]}");
}

/* Phase 4.3b: base-field arithmetic over Q(constants) */
static void test_nf_base_field(void) {
    BR_EQ("Block[{tragerAlgConstants = {2^(1/3)}}, tragerSquareFreeList[z^2 - 2^(1/3) z^3]]",
          "{{1, 1}, {z, 2}, {1 - 2^(1/3) z, 1}}");
    BR_EQ("Block[{tragerAlgConstants = {}}, tragerSquareFreeList[z^2 - 2 z^3]]", "{{-1, 1}, {z, 2}, {-1 + 2 z, 1}}");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{r = tragerTogether[3/(2 z^3 + 2 Sqrt[2] z^6) + "
            "(-3 + 3 z + z^2 - 3 Sqrt[2] z^3 - 2 Sqrt[2] z^5)/(2 z^3 - 2 z^4 + 4 Sqrt[2] z^6 - 2 Sqrt[2] z^7 + 4 z^9)]}, "
            "Head[r] =!= Plus && Exponent[Denominator[r], z, Min] === 1]]");
    BR_EQ("Block[{tragerAlgConstants = {}}, tragerTogether[1/x + 1/(x + 1)]]", "(1 + 2 x)/(x (1 + x))");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, tragerGCD[x^2 - 2, x - Sqrt[2]]]", "-Sqrt[2] + x");
    BR_EQ("Block[{tragerAlgConstants = {}}, tragerGCD[x^2 - 1, x - 1]]", "-1 + x");
    /* former false verdicts (Hermite over Q missed a repeated factor / a cancellation) */
    BR_ANTIDERIV("(x - 2^(1/3))^(-2/3)", "Bronstein");
    BR_FALSE("MatchQ[Block[{nfDemoteCertificate = Function[r, r]}, BronsteinIntegrate[Together[D[Sqrt[x^3 + Sqrt[2]] + "
             "Log[Sqrt[x^3 + Sqrt[2]] + x], x]], x]], {\"not elementary\", _}]");
}

static void test_nf_integrals(void) {
    BR_ANTIDERIV("Sqrt[2] Sqrt[Log[x]]/x", "Bronstein");
    BR_ANTIDERIV("Sqrt[x + Sqrt[2]]/x", "Bronstein");
    BR_ANTIDERIV("Sqrt[Sqrt[2] + Log[x]]/x", "Bronstein");
    BR_ANTIDERIV("Sqrt[x + 2^(1/3)]/x", "Bronstein");
    BR_ANTIDERIV("D[Log[Sqrt[x + Sqrt[2]] + 1], x]", "Bronstein");
    BR_ANTIDERIV("x/Sqrt[x^4 + Sqrt[5]]", "Bronstein");
    /* by linearity the constant leaves: a certificate over Q */
    BR_EQ("BronsteinIntegrate[Sqrt[2] Sqrt[Log[x]], x][[1]]", "\"not elementary\"");
    /* otherwise a verdict over a number field is not yet a proof */
    BR_EQ("BronsteinIntegrate[Sqrt[Sqrt[2] + Log[x]], x]",
          "{\"implementation incomplete\", \"non-elementarity over a number-field constant field is not yet certified "
          "(no algebraic part: the Risch DE D b + s b = c has no solution ([B] Prop 2.15/2.18))\"}");
    BR_EQ("First[BronsteinIntegrate[1/((x - Sqrt[2]) Sqrt[x^3 + 1]), x]]", "\"implementation incomplete\"");
    BR_EQ("BronsteinScope[Sqrt[x + Sqrt[2]], x]", "\"InScope\"");               /* not nesting */
    /* review of 4.3: the linearity rule must not certify a zero factor */
    BR_EQ("BronsteinIntegrate[((1 + Sqrt[2])^2 - 3 - 2 Sqrt[2]) Sqrt[Log[x]], x]", "0");
    BR_EQ("BronsteinIntegrate[((Sqrt[2] - 1) (Sqrt[2] + 1) - 1) E^(x^2), x]", "0");
    BR_FALSE("MatchQ[BronsteinIntegrate[(Log[4] - 2 Log[2]) Sqrt[Log[x]], x], {\"not elementary\", _}]");
    BR_TRUE("brAlgebraicConstantQ[1 + Sqrt[3]] && !brAlgebraicConstantQ[Log[4] - 2 Log[2]] && !brAlgebraicConstantQ[Pi]");
    /* PMT: no certificate with radical constants (its arithmetic is over Q) */
    BR_FALSE("MatchQ[Quiet[Integrate`ParallelMixedTower[1/(Sqrt[2] + x^2)^(3/2), x]], {\"not elementary\", ___}]");
    /* regression: Sqrt[Tan[x]] (the nested Return in residueField) */
    BR_TRUE("!ListQ[BronsteinIntegrate[Sqrt[Tan[x]], x]]");
}

int main(void) {
    br_init();
    test_nf_primitives();
    test_nf_base_field();
    test_nf_integrals();
    BR_DONE("bronstein_nf_tests");
    return 0;
}
