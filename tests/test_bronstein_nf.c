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
#include <string.h>

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

/* Phase 4.3b C3: residue and log fields above the constant field */
static void test_nf_residue_fields(void) {
    /* without constants in use: the relative presentation, generator kept */
    BR_EQ("Block[{tragerAlgConstants = {}}, fieldOfRoot[Z^2 - 2, Z, Z][\"gamma\"]]", "Z");
    BR_EQ("Block[{tragerAlgConstants = {}}, fieldOfRoot[Z^2 - Sqrt[2], Z, Z][\"F\"][\"mu\"]]", "-Sqrt[2] + Z^2");
    /* with Sqrt[2] in use: the absolute field Q(2^(1/4)), Sqrt[2] recorded */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{fr = fieldOfRoot[Z^2 - Sqrt[2], Z, t]}, "
            "Exponent[fr[\"F\"][\"mu\"], t] === 4 && kc[fr[\"F\"], fr[\"gamma\"]^2 - Sqrt[2]] === 0]]");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, KeyExistsQ[fieldOfRoot[Z^2 - Sqrt[2], Z, t][\"F\"], \"consts\"]]");
    /* the constant field as a number field, and the canonical form over it */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{F = constNF[]}, Exponent[F[\"mu\"], F[\"t\"]] === 2]]");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{d = -1/3*1/Sqrt[2] + (2*Sqrt[2])/27}, "
          "kcConstQ[1 + 1/(3 Sqrt[2] d) - 2 Sqrt[2]/(27 d)]]]", "0");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{r = kcConstQ[x/(1 + Sqrt[2]) + 3]}, "
            "FreeQ[Denominator[Together[r]], Sqrt[2]] && PossibleZeroQ[r - x/(1 + Sqrt[2]) - 3]]]");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, kc[fieldQ[], (x + 1)^2]]", "1 + 2 x + x^2");
    /* integrands whose residue orbits carry the constant */
    BR_ANTIDERIV("Together[D[Log[Sqrt[x^2 + 1] + x + Sqrt[2]], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[Log[Sqrt[x^3 + 1] + Sqrt[2] x + 1], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[Log[Sqrt[x^4 + 1] + Sqrt[2] x], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[Log[Sqrt[E^x + 2] + Sqrt[2]], x]]", "Bronstein");
    BR_ANTIDERIV("Sqrt[E^x + Sqrt[3]]", "Bronstein");
    BR_ANTIDERIV("Together[D[2^(1/3) Log[(x^3 + 1)^(1/3) + 2^(1/3) x], x]]", "Bronstein");
    /* the generator-rename trap: an orbit over Q(i) whose factor carries Sqrt[2] */
    BR_ANTIDERIV("(1 + I) Sqrt[x]/(x^2 + I)", "Bronstein");
    /* trace shortcut: no verdict over a number field of radical constants */
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, traceCertificate[<|\"Z\" -> z, \"R\" -> z^2 - 2|>]]", "Null");
}

/* Phase 4.3b C4: monic factors in canonical form over the constant field */
static void test_nf_monic(void) {
    BR_EQ("Block[{tragerAlgConstants = {}}, tragerMonic[2 z^2 + 4, z]]", "2 + z^2");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{m = tragerMonic[(1 + Sqrt[2]) z + 1, z]}, "
            "Coefficient[m, z, 1] === 1 && FreeQ[m, Power[_Plus, -1]] && PossibleZeroQ[(m /. z -> 0) - 1/(1 + Sqrt[2])]]]");
    /* two shapes of the same factor normalise to one canonical polynomial */
    BR_TRUE("Block[{tragerAlgConstants = {2^(1/3)}}, tragerMonic[2^(1/3) z - 1, z] === tragerMonic[2 z - 2^(2/3), z]]");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, brNonzeroResiduePolynomial[z^2 (Sqrt[2] z - 1)^2, z]]",
          "-1/Sqrt[2] + z");
    BR_ANTIDERIV("Together[D[Log[Sqrt[2] x + Sqrt[x^2 + 1]], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[Log[(1 + Sqrt[2]) x^2 + Sqrt[x^4 + 1]], x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[Log[2^(1/3) x + (x^3 + 1)^(1/3)], x]]", "Bronstein");
    BR_ANTIDERIV("1/((Sqrt[2] x^2 + 1) Sqrt[x^2 + 1])", "Bronstein");
    BR_ANTIDERIV("x/((Sqrt[2] x^3 - 1) Sqrt[x^2 + 2])", "Bronstein");
}

/* Phase 4.3b C6: Bronstein's k arithmetic over Q(constants) */
static void test_nf_k_arithmetic(void) {
    BR_TRUE("T6 = brTower[Log[x] Sqrt[Log[x] + 1], x]; th6 = T6[\"theta\"]; True");
    BR_TRUE("Block[{tragerAlgConstants = {2^(1/3)}}, brConstantCoefficientsQ[th6 - 2^(2/3)]]");
    BR_FALSE("Block[{tragerAlgConstants = {}}, brConstantCoefficientsQ[th6 - 2^(2/3)]]");
    BR_FALSE("Block[{tragerAlgConstants = {Sqrt[2]}}, brConstantCoefficientsQ[th6 - 2]]");
    /* the C builtin lost the polynomial over Q(2^(1/3)) ({{}, {}}); over K0 it splits */
    BR_EQ("Block[{tragerAlgConstants = {2^(1/3)}}, brSplitSquarefreeNF[T6, Expand[(th6 - 2^(2/3))^2 (th6 + 1)]]]",
          "{{1 + brT1, 2^(2/3) - brT1}, {1, 1}}");
    BR_EQ("Block[{tragerAlgConstants = {2^(1/3)}}, brSplitSquarefree[T6, Expand[(th6 - 2^(2/3))^2 (th6 + 1)]]]",
          "{{1 + brT1, 2^(2/3) - brT1}, {1, 1}}");
    /* an exponential theta is special */
    BR_TRUE("T7 = brTower[E^x Sqrt[E^x + 1], x]; th7 = T7[\"theta\"]; "
            "Block[{tragerAlgConstants = {Sqrt[2]}}, brSpecialPartNF[T7, th7] === th7 && brSpecialPartNF[T7, th7 - Sqrt[2]] === 1]");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, brNormalQ[T7, th7 - Sqrt[2]] && !brNormalQ[T7, th7^2 - Sqrt[2] th7]]");
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, brSpecialQ[T7, Sqrt[2] th7] && !brSpecialQ[T7, th7 - Sqrt[2]]]");
    /* the coefficient field of Hermite over K0, and the map back */
    BR_TRUE("Block[{tragerAlgConstants = {2^(1/3)}}, With[{K = brHermiteK[T6]}, K[\"Type\"] === \"Fn\" && KeyExistsQ[K, \"consts\"] && K[\"theta\"] === 2^(1/3)]]");
    BR_EQ("Block[{tragerAlgConstants = {2^(1/3)}}, With[{K = brHermiteK[T6]}, "
          "brHermiteBack[K, kQuo[K, kc[K, Expand[(th6 - 2^(2/3))^2 (th6 + 1)]], kc[K, th6 - 2^(2/3)], th6]] - "
          "Expand[(th6 - 2^(2/3)) (th6 + 1)]]]", "0");
    /* former FALSE certificates (demoted) over Q(2^(1/3)): now verified answers */
    BR_ANTIDERIV("Together[D[Sqrt[2^(1/3) + Log[x]]/(-2^(2/3) + Log[x]), x]]", "Bronstein");
    BR_ANTIDERIV("Together[D[(Log[x] Sqrt[2^(1/3) + Log[x]])/(-2^(1/3) + Log[x]^2), x]]", "Bronstein");
    /* a sum Together over the extension left uncombined: read as not algebraic */
    BR_ANTIDERIV("Together[D[(Sqrt[2] + Log[x])^(1/3)/(-(Sqrt[2]/(1 + Sqrt[2])) + Log[x]), x]]", "Bronstein");
    /* exact cancellation over k(K0): the factor divides only modulo mu */
    BR_EQ("Block[{tragerAlgConstants = {2^(1/3)}}, brCancelK[T6, Expand[(th6 - 2^(2/3)) (th6 + 1)]/Expand[(th6 - 2^(2/3))^2]]]",
          "(1 + brT1)/(-2^(2/3) + brT1)");
    BR_EQ("Block[{tragerAlgConstants = {}}, brCancelK[T6, (th6^2 - 1)/(th6 - 1)]]", "1 + brT1");
    BR_ANTIDERIV("Together[D[((3 - Sqrt[2] + Log[x])^(2/3) (Sqrt[2] + Log[x]))/(1 - Sqrt[2] + Log[x]), x]]", "Bronstein");
    /* the combination helpers */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, !tragerCombinedQ[(1/(x - y) + 2)/(x + 1)] && tragerCombinedQ[(x + 1)/(x - Sqrt[2])]]");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, tragerCombineConstNF[1/(x - Sqrt[2]) - 1/(x - Sqrt[2])]]", "0");
}

/* Phase 4.3b C7: radicals equal over a larger constant field */
static void test_nf_radical_content(void) {
    BR_EQ("tragerRadicalContent[2 x]", "2");
    BR_EQ("tragerRadicalContent[-6 x/5]", "6/5");
    BR_EQ("tragerRadicalContent[x + 1]", "1");
    BR_EQ("tragerExtractRadicalContent[Sqrt[2 x] + Sqrt[x], x]", "Sqrt[x] + Sqrt[2] Sqrt[x]");
    BR_EQ("tragerExtractRadicalContent[Sqrt[-2 x] + Sqrt[-x], x]", "Sqrt[-x] + Sqrt[2] Sqrt[-x]");
    /* only when it merges radicands: a lone radical keeps its constant field Q */
    BR_EQ("tragerExtractRadicalContent[Sqrt[-2 x], x]", "Sqrt[-2 x]");
    BR_EQ("tragerExtractRadicalContent[(3 x)^(-4/5) (x + 1), x]", "(1 + x)/(3 x)^(4/5)");
    BR_ANTIDERIV("(-3 x + 2)/(2 x - 2) (3 x)^(-4/5)", "Trager");
    BR_EQ("tragerExtractRadicalContent[(8 x)^(1/3), x]", "2 x^(1/3)");
    BR_EQ("tragerExtractRadicalContent[Sqrt[x + 1], x]", "Sqrt[1 + x]");
    /* an identity on every branch: the value at a negative and a complex point */
    BR_TRUE("With[{e = Sqrt[6 x/5] + (16 x^2 - 8)^(1/3)}, With[{d = tragerExtractRadicalContent[e, x] - e},"
            " Abs[N[d /. x -> -7/3, 30]] < 10^-25 && Abs[N[d /. x -> 1/2 + 3 I, 30]] < 10^-25]]");
    BR_ANTIDERIV("Sqrt[2 x] + Sqrt[x]", "Bronstein");
    BR_ANTIDERIV("Sqrt[2 x]/(1 + Sqrt[x])", "Bronstein");
    BR_ANTIDERIV("Sqrt[3 Log[x]] Sqrt[Log[x]]/x", "Bronstein");
    BR_ANTIDERIV("Sqrt[2 x + 2] Sqrt[x + 1] Log[x]", "Bronstein");
    BR_ANTIDERIV("Sqrt[2 x] + Sqrt[x]", "Trager");
    /* genuinely two radicals: still Incomplete */
    BR_EQ("First[BronsteinIntegrate[1/(Sqrt[2 x] Sqrt[x^3 + 1]) + Sqrt[x^3 + 1]/Sqrt[x], x]]", "\"implementation incomplete\"");
}

/* Adversarial review of 4.3b (2026-11-15): 14 false certificates from the
 * per-orbit log part over Q(constants), hidden by nfDemoteCertificate.  Run
 * with the demotion OFF: none may come back "not elementary". */
static void test_nf_review_orbits(void) {
    /* the relative structure of an absolute orbit field over K0 = Q(Sqrt[2]) */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, fr8 = fieldOfRoot[Z^2 - Sqrt[2], Z, Z]; True]");
    BR_FALSE("fr8[\"F\"][\"t\"] === Z");                                 /* never the root symbol itself */
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, fr8[\"F\"][\"rel\"][\"d\"]]", "2");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, fr8[\"F\"][\"rel\"][\"m0\"]]", "2");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, relToZ[fr8[\"F\"], fr8[\"gamma\"]]]", "Z");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, tragerTogether[relToZ[fr8[\"F\"], fr8[\"gamma\"]^2] - Sqrt[2]]]", "0");
    /* the relative trace: Tr(gamma) = 0, Tr(gamma^2) = 2 Sqrt[2] (an absolute trace would give 0) */
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, relTrace[fr8[\"F\"], fr8[\"gamma\"]]]", "0");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, relTrace[fr8[\"F\"], fr8[\"gamma\"]^2]]", "2 Sqrt[2]");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, Length[First[relCoeffs[fr8[\"F\"], fr8[\"gamma\"]]]]]", "4");
    BR_EQ("Block[{tragerAlgConstants = {Sqrt[2]}}, Keys[relativeStructure[fr8[\"F\"], fr8[\"gamma\"], Z^2 - Sqrt[2], Z, {Sqrt[2]}]]]",
          "{\"Z\", \"q\", \"d\", \"theta0\", \"m0\", \"Minv\", \"s\"}");
    /* the constant field's generator is written in the constants, never a Root object */
    BR_EQ("gfConstField[{I, Sqrt[2]}][\"theta\"]", "I + Sqrt[2]");
    BR_EQ("gfConstField[{Sqrt[2], Sqrt[3]}][\"theta\"]", "Sqrt[2] + Sqrt[3]");
    /* an algebraic constant denominator whose image vanishes makes the prime bad:
     * 3 - Sqrt[2] has norm 7, and at p = 7 the image of Sqrt[2] is 3 or 4 */
    BR_TRUE("Block[{tragerAlgConstants = {Sqrt[2]}}, With[{r = gfConstImage[Sqrt[2], 7]}, "
            "denominatorsPrimeToQ[x/(3 - Sqrt[2]), 7] === (Mod[3 - r, 7] =!= 0)]]");
    const char* cases[] = {
        "1/((x^2 - Sqrt[2]) Sqrt[x + 1])", "x/Sqrt[Expand[(x^2 - Sqrt[2])^2 (x + 1)]]",
        "1/((x - I) Sqrt[x - Sqrt[2]])", "1/((x^2 + 1) Sqrt[x - Sqrt[2]])", "1/((x^2 + 1) Sqrt[x + Sqrt[2]])",
        "1/((x^2 + 1) Sqrt[x + 1 + Sqrt[2]])", "1/((x^2 - 3) Sqrt[x + Sqrt[2]])", "1/((x^2 - Sqrt[2]) Sqrt[x])",
        "1/((x^2 + Sqrt[2]) Sqrt[x + 1])", "x/((x^2 + x + 1) Sqrt[x - Sqrt[2]])",
        "1/(x (Log[x]^2 + 1) Sqrt[Log[x] + Sqrt[2]])", "1/(x (Log[x]^2 - Sqrt[2]) Sqrt[Log[x]])",
        "E^x/((E^(2 x) + 1) Sqrt[E^x + Sqrt[2]])", NULL };
    for (int i = 0; cases[i]; i++) {
        char buf[600];
        snprintf(buf, sizeof buf,
                 "!MatchQ[Block[{nfDemoteCertificate = Function[r, r]}, BronsteinIntegrate[%s, x]], {\"not elementary\", ___}]",
                 cases[i]);
        BR_TRUE(buf);
        if (!strstr(cases[i], "Log") && !strstr(cases[i], "E^")) {
            snprintf(buf, sizeof buf,
                     "!MatchQ[Block[{nfDemoteCertificate = Function[r, r]}, Trager`TragerIntegrate[%s, x]], {\"not elementary\", ___}]",
                     cases[i]);
            BR_TRUE(buf);
        }
    }
    BR_ANTIDERIV("1/((x^2 - Sqrt[2]) Sqrt[x + 1])", "Bronstein");
    BR_ANTIDERIV("1/((x - I) Sqrt[x - Sqrt[2]])", "Bronstein");
    BR_ANTIDERIV("1/(x (Log[x]^2 + 1) Sqrt[Log[x] + Sqrt[2]])", "Bronstein");
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
    test_nf_residue_fields();
    test_nf_monic();
    test_nf_k_arithmetic();
    test_nf_radical_content();
    test_nf_review_orbits();
    test_nf_integrals();
    BR_DONE("bronstein_nf_tests");
    return 0;
}
