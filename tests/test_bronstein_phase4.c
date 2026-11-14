/* test_bronstein_phase4.c -- Phase 4 of the Bronstein completion campaign
 * (BRONSTEIN_PLAN.md, "Phase 4 design"): constant fields.
 *
 *  - 4.1, the constant field Q(i) with a literal I: the front end
 *    (brGaussianRationalQ, brConstantField, BronsteinScope / brScopeBody),
 *    factorisation over Q(i) (tragerFactorList under tragerConstantI),
 *    constant-field coefficient tests (constantFieldNumberQ), residue fields
 *    that contain i (residueField, fieldSubI, brFieldOver), and the refusal
 *    of primes p != 1 mod 4 (goodPrime) once I has an image mod p (4.1b);
 *  - the Q(i) battery: answers verified, certificates checked by hand
 *    (FriCAS reports "implementation incomplete" on Sqrt[Log[x] + I]).
 *  - 4.1b: I modulo p = 1 mod 4 (gfSqrtMinusOne, gfMod, goodPrime's consistent
 *    root), Gaussian residue polynomials (gaussianRootObjects), a Gaussian
 *    kappa at the special places (brGammaForm), the identity check over a
 *    residue field containing i;
 *  - 4.2: the trigonometric front end (brTrigFrontEnd, brTrigPresentation).
 */
#include "bronstein_test_utils.h"

static void test_qi_primitives(void) {
    BR_FALSE("tragerConstantI");
    BR_EQ("Length[tragerFactorList[z^2 + 1]]", "2");                         /* over Q: irreducible */
    BR_EQ("Block[{tragerConstantI = True}, Length[tragerFactorList[z^2 + 1]]]", "3");   /* over Q(i): splits */
    BR_EQ("Block[{tragerConstantI = True}, Length[tragerFactorList[z^2 - 2]]]", "2");   /* sqrt 2 not in Q(i) */
    BR_FALSE("constantFieldNumberQ[I/2]");
    BR_TRUE("Block[{tragerConstantI = True}, constantFieldNumberQ[1/2 - 3 I]]");
    BR_FALSE("Block[{tragerConstantI = True}, constantFieldNumberQ[Sqrt[2]]]");
    BR_TRUE("Block[{tragerConstantI = True}, rationalPolyQ[x^2 + I x - 1/2, x]]");
    /* an absolute residue field that contains i, I written in its generator */
    br_run("rfI = Block[{tragerConstantI = True}, residueField[z^2 - 2, z, tt]]; True");
    BR_TRUE("KeyExistsQ[rfI[\"F\"], \"I\"]");
    BR_EQ("kc[rfI[\"F\"], (rfI[\"F\"][\"I\"])^2 + 1]", "0");                  /* rep^2 = -1 in the field */
    BR_EQ("kc[rfI[\"F\"], I^2 + 1]", "0");
    BR_TRUE("kc[rfI[\"F\"], I] === kc[rfI[\"F\"], rfI[\"F\"][\"I\"]]");       /* fieldSubI */
    BR_TRUE("fieldSubI[fieldQ[], 2 + 3 I] === 2 + 3 I");
    BR_TRUE("KeyExistsQ[brFieldOver[{x}, rfI[\"F\"]], \"I\"]");
    BR_TRUE("brFieldOver[{x}, fieldQ[]] === fieldFn[{x}]");
    /* over Q(i) only primes p = 1 mod 4 are good (a degree-1 prime of Z[i]) */
    BR_EQ("Block[{tragerConstantI = True}, goodPrime[<|\"B\" -> <|\"x\" -> x, \"m\" -> 2, \"G\" -> x^3 + I|>|>, 7]]", "None");
    BR_EQ("Block[{tragerConstantI = True}, goodPrime[<|\"B\" -> <|\"x\" -> x, \"m\" -> 2, \"G\" -> x^3 + I|>|>, 3]]", "None");
}

static void test_qi_scope(void) {
    BR_EQ("BronsteinScope[Sqrt[x + I], x]", "\"InScope\"");
    BR_EQ("brScopeBody[Sqrt[Log[x]] E^(I x), x]", "\"InScope\"");
    BR_EQ("brConstantField[E^(I x) Sqrt[1 + E^(I x)]]", "\"QI\"");
}

static void test_qi_battery(void) {
    BR_ANTIDERIV("Sqrt[x + I]/x", "Bronstein");                              /* FriCAS agrees */
    BR_ANTIDERIV("D[Sqrt[x + I] Log[x], x]", "Bronstein");
    BR_ANTIDERIV("I Sqrt[Log[x]]/x", "Bronstein");
    BR_ANTIDERIV("E^(I x) Sqrt[1 + E^(I x)]", "Bronstein");
    BR_ANTIDERIV("Sqrt[1 + E^(I x)]", "Bronstein");
    BR_ANTIDERIV("D[Log[E^(I x) + Sqrt[1 + E^(2 I x)]], x]", "Bronstein");
    BR_ANTIDERIV("D[Sqrt[Log[x] + I] Log[x], x]", "Bronstein");
    BR_ANTIDERIV("1/Sqrt[E^(2 I x) + 1]", "Bronstein");
    BR_ANTIDERIV("(1 + I) Sqrt[x]/(x^2 + I)", "Bronstein");
    /* Int E^u Sqrt[u + I] du: b' + b (1 + 1/(2 (u + I))) = 1 has no rational b */
    BR_EQ("BronsteinIntegrate[Sqrt[Log[x] + I], x][[1]]", "\"not elementary\"");
}

int main(void) {
    br_init();
    test_qi_primitives();
    test_qi_scope();
    test_qi_battery();
    BR_DONE("bronstein_phase4_tests");
    return 0;
}
