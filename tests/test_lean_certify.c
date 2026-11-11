/*
 * test_lean_certify.c -- Certify / ProofObject certificate GENERATION
 * (src/internal/lean/Lean.m).  Needs no Lean installation: every case runs
 * with Verify -> False (or the heuristic Engine -> "Mathilda") and inspects
 * the generated statement, source and properties.  The Lean-judged cases are
 * in test_lean_verify.c.
 */
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* input must evaluate to True. */
static void check_true(const char* input) {
    assert_eval_eq(input, "True", 0);
}

/* Evaluate input; it must give a String equal to expected. */
static void check_str(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL: does not parse: %s\n", input); abort(); }
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    if (!v || v->type != EXPR_STRING || strcmp(v->data.string, expected) != 0) {
        char* got = v ? expr_to_string(v) : NULL;
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected,
                got ? got : "(null)");
        free(got);
        abort();
    }
    expr_free(v);
}

static void test_identity_source(void) {
    assert_eval_eq("leanO = Certify[(x+1)^3 == x^3+3x^2+3x+1, Verify -> False]; Null", "Null", 0);
    check_true("leanO[\"Verified\"] === Missing[\"NotChecked\"]");
    check_true("StringContainsQ[leanO[\"LeanSource\"], \"theorem mathilda_cert_1\"]");
    check_true("StringContainsQ[leanO[\"LeanSource\"], \"#print axioms mathilda_cert_1\"]");
    check_true("StringContainsQ[leanO[\"LeanSource\"], \"namespace Mathilda\"]");
    check_true("StringContainsQ[leanO[\"LeanSource\"], \"set_option maxHeartbeats 2400000 in\"]");
    check_true("!StringContainsQ[leanO[\"LeanSource\"], \"sorry\"]");
    check_str("leanO[\"LeanTheorem\"]",
              "\xE2\x88\x80 x : Int, (1 + x) ^ 3 = 1 + 3 * x + 3 * x ^ 2 + x ^ 3");
    check_str("leanO[\"ProofScript\"]", "by\n  intros\n  grind");
    check_true("leanO[\"Properties\"] === {\"Axioms\", \"Cached\", \"Certificate\", \"Domain\", "
               "\"Engine\", \"Hypotheses\", \"Input\", \"Kind\", \"LeanSource\", \"LeanTheorem\", "
               "\"LeanVersion\", \"Logs\", \"ProofScript\", \"Result\", \"Theorem\", \"Tier\", "
               "\"Timing\", \"TrustedSteps\", \"Verified\"}");
    check_true("leanO[\"Kind\"] === \"Identity\" && leanO[\"Engine\"] === \"Lean\" && "
               "leanO[\"Tier\"] === 0 && leanO[\"Domain\"] === Integers");
    check_true("leanO[\"TrustedSteps\"] === {} && leanO[\"Hypotheses\"] === {}");
    check_true("leanO[\"Cached\"] === False");
    check_true("leanO[\"Nope\"] === Missing[\"UnknownProperty\", \"Nope\"]");
    check_str("ToString[leanO]",
              "ProofObject[\xC2\xABIdentity\xC2\xBB Verified: Missing[NotChecked] | Lean Tier 0 | 0.0 s]");
    /* Deterministic: the same claim gives the same source. */
    check_true("Certify[(x+1)^3 == x^3+3x^2+3x+1, Verify -> False][\"LeanSource\"] === "
               "leanO[\"LeanSource\"]");
}

static void test_computations(void) {
    check_true("ProofObject[Factor[x^2-1], Engine -> \"Lean\", Verify -> False][\"Kind\"] == "
               "\"Factorization\"");
    assert_eval_eq("leanF = Certify[Factor[x^10 - 1], Verify -> False]; Null", "Null", 0);
    check_true("leanF[\"TrustedSteps\"] === {\"irreducibility of the factors (not certified)\"}");
    check_true("StringStartsQ[leanF[\"LeanTheorem\"], \"\xE2\x88\x80 x : Int, x ^ 10 - 1 = \"]");
    check_true("leanF[\"Input\"] === HoldForm[Factor[x^10 - 1]]");
    check_true("leanF[\"Result\"] === Factor[x^10 - 1]");
    check_true("Certify[Expand[(x+y)^2], Verify -> False][\"Kind\"] === \"Identity\"");
    check_true("Certify[Collect[(x+1)^2 + x y, x], Verify -> False][\"Kind\"] === \"Identity\"");
    /* A rational coefficient selects Rationals. */
    check_true("Certify[Expand[(x/2+1)^2], Verify -> False][\"Domain\"] === Rationals");
    check_true("StringContainsQ[Certify[Expand[(x/2+1)^2], Verify -> False][\"LeanTheorem\"], "
               "\"(1 / 4 : Rat)\"]");
}

static void test_rational_identity(void) {
    assert_eval_eq("leanT = Certify[Together[1/(x-1) - 1/(x+1)], Verify -> False]; Null", "Null", 0);
    check_true("leanT[\"Certificate\"][\"Form\"] === \"ClearedDenominators\"");
    check_true("Length[leanT[\"Hypotheses\"]] === 2");
    check_true("Length[leanT[\"TrustedSteps\"]] === 1 && "
               "StringStartsQ[First[leanT[\"TrustedSteps\"]], \"field algebra: the cross-multiplied\"]");
    check_true("!StringContainsQ[leanT[\"LeanTheorem\"], \"/\"]");
    check_true("Certify[Apart[(x^2+1)/(x^2-1)], Verify -> False][\"Certificate\"][\"Form\"] === "
               "\"ClearedDenominators\"");
}

static void test_options(void) {
    /* Hypotheses: the automatic domain is every field of characteristic 0
     * (over Int or Rat an implication can be vacuous: x^2 + 4 == 0 has no
     * rational solution, but x = 2 I is one in C); an order relation makes it
     * every ordered field.  Domain -> Integers is honoured and labelled. */
    check_str("Certify[x^2 == x, Assumptions -> {x == 1}, Verify -> False][\"LeanTheorem\"]",
              "\xE2\x88\x80 {K : Type} [Lean.Grind.Field K] [Lean.Grind.IsCharP K 0], \xE2\x88\x80 x : K, x = 1 \xE2\x86\x92 x ^ 2 = x");
    check_true("Certify[x == 0, Assumptions -> 2 x == 1, Verify -> False][\"Domain\"] === "
               "Complexes");
    check_true("Certify[x == 0, Assumptions -> 2 x == 1, Domain -> Rationals, Verify -> False]"
               "[\"Domain\"] === Rationals");
    check_true("Certify[x == 0, Assumptions -> 2 x == 1, Domain -> Integers, Verify -> False]"
               "[\"Domain\"] === Integers");
    check_true("Certify[x + 1 == 1 + x, Verify -> False][\"Domain\"] === Integers");
    /* Variables that occur only in the assumptions are bound too. */
    check_str("Certify[x == x, Assumptions -> {y > 0, Element[x, Reals]}, Verify -> False]"
              "[\"LeanTheorem\"]",
              "\xE2\x88\x80 {K : Type} [Lean.Grind.Field K] [LE K] [LT K] [Std.LawfulOrderLT K] "
              "[Std.IsLinearOrder K] [Lean.Grind.OrderedRing K], \xE2\x88\x80 x y : K, y > 0 \xE2\x86\x92 x = x");
    check_true("Certify[x == x, Assumptions -> {y > 0, Element[x, Reals]}, Verify -> False]"
               "[\"Domain\"] === Reals");
    /* An order relation makes only ITS variables real (review 3, finding B). */
    assert_eval_eq("Check[Certify[x == x, Assumptions -> y > 0, Verify -> False], \"realvars\"]",
                   "\"realvars\"", 0);
    /* The theorem takes the generic binders as parameters; the footer states
     * the full statement. */
    check_true("StringContainsQ[Certify[x^2 == x, Assumptions -> {x == 1}, Verify -> False]"
               "[\"LeanSource\"], \"theorem mathilda_cert_1 {K : Type} [Lean.Grind.Field K] \"]");
    check_true("StringContainsQ[Certify[x == x, TimeConstraint -> 5, Verify -> False]"
               "[\"LeanSource\"], \"set_option maxHeartbeats 200000 in\"]");
    check_str("Certify[x == x, Method -> \"omega\", Verify -> False][\"ProofScript\"]",
              "by\n  intros\n  omega");
    check_true("Certify[x + 1 == 1 + x, Domain -> Rationals, Verify -> False][\"Domain\"] === "
               "Rationals");
    check_true("Options[Certify][[1]] === (Engine -> \"Lean\")");
    /* A symbol holding a claim is evaluated, then dispatched. */
    check_true("leanC = (x + 1)^2 == x^2 + 2 x + 1; Certify[leanC, Verify -> False][\"Kind\"] === "
               "\"Identity\"");
}

static void test_errors(void) {
    test_mute_messages();
    assert_eval_eq("Certify[Expand[(x+1)^2], Engine -> \"Foo\"]", "$Failed", 0);
    assert_eval_eq("Check[Certify[Expand[(x+1)^2], Engine -> \"Foo\"], \"engine\"]",
                   "\"engine\"", 0);
    assert_eval_eq("Check[Certify[Sin[x]^2 + Cos[x]^2 == 1, Verify -> False], \"nocert\"]",
                   "\"nocert\"", 0);
    assert_eval_eq("Check[Certify[x == x, Tier -> 1, LeanProject -> None], \"nomathlib\"]",
                   "\"nomathlib\"", 0);
    assert_eval_eq("Check[Certify[x == x, Domain -> Reals, LeanProject -> None], \"tier\"]",
                   "\"tier\"", 0);
    /* Not in Q(x): no Tier 0 antiderivative certificate. */
    assert_eval_eq("Certify[Integrate[Sin[x], x], Verify -> False]", "$Failed", 0);
    assert_eval_eq("Check[Certify[Solve[x^2 == 2, x], Verify -> False], \"nocert\"]",
                   "\"nocert\"", 0);
    assert_eval_eq("Certify[Inverse[{{1, a}, {0, 1}}], Verify -> False]", "$Failed", 0);
    assert_eval_eq("Certify[Det[IdentityMatrix[7]], Verify -> False]", "$Failed", 0);
}

/* Row 4: the cleared chain-rule identity for an antiderivative. */
static void test_antiderivative(void) {
    assert_eval_eq("leanR = Certify[Integrate[(x^4-3x^2+6)/(x^6-5x^4+5x^2+4), x], "
                   "Verify -> False]; Null", "Null", 0);
    check_true("leanR[\"Kind\"] === \"Antiderivative\" && leanR[\"Tier\"] === 0 && "
               "leanR[\"Domain\"] === Rationals");
    check_true("leanR[\"Result\"] === ArcTan[x] + ArcTan[x^3] + ArcTan[(x - 3 x^3 + x^5)/2]");
    check_true("Length[leanR[\"TrustedSteps\"]] === 4");
    check_true("leanR[\"Certificate\"][\"PreCheck\"] === True");
    check_true("leanR[\"Certificate\"][\"Form\"] === \"ClearedDerivative\"");
    check_true("Length[leanR[\"Hypotheses\"]] === 4");
    /* The structured pieces are visible: q^2 + p^2 of the third ArcTan and
     * p' q - p q' of its argument. */
    check_true("StringContainsQ[leanR[\"LeanTheorem\"], "
               "\"(2 ^ 2 + (x - 3 * x ^ 3 + x ^ 5) ^ 2)\"]");
    check_true("StringContainsQ[leanR[\"LeanTheorem\"], "
               "\"((1 - 9 * x ^ 2 + 5 * x ^ 4) * 2 - (x - 3 * x ^ 3 + x ^ 5) * 0)\"]");
    check_true("StringEndsQ[leanR[\"LeanTheorem\"], \" = 0\"]");
    /* Surds: Sqrt[2] becomes s1 and the identity holds modulo s1^2 - 2 as an
     * unconditional polynomial identity -- never a vacuous s1^2 = 2 -> ... */
    assert_eval_eq("leanS = Certify[Integrate[1/(x^2+2), x], Verify -> False]; Null", "Null", 0);
    check_true("StringStartsQ[leanS[\"LeanTheorem\"], \"\xE2\x88\x80 x s1 : Rat, \"]");
    check_true("StringContainsQ[leanS[\"LeanTheorem\"], \"(s1 ^ 2 - 2) * \"]");
    check_true("!StringContainsQ[leanS[\"LeanSource\"], \"s1 ^ 2 = 2 \xE2\x86\x92\"]");
    check_true("leanS[\"Certificate\"][\"Surds\"][[1, 2]] === Sqrt[2]");
    /* Prime surd generators: Sqrt[3/2] and Sqrt[6] share s1 = Sqrt[2], s2 = Sqrt[3]. */
    check_true("Certify[Integrate[1/(3x^2-2), x], Verify -> False][\"Certificate\"]"
               "[\"Surds\"][[All, 2]] === {Sqrt[2], Sqrt[3]}");
    check_true("Certify[Integrate[1/(x^2-1), x], Verify -> False][\"Certificate\"]"
               "[\"PreCheck\"] === True");
    /* Test hook: a wrong antiderivative fails the local pre-check but is still
     * emitted, so that Lean (not Mathilda) refutes it. */
    assert_eval_eq("Lean`Private`$antiderivOverride = ArcTan[x] + ArcTan[x^3]; "
                   "leanW = Certify[Integrate[(x^4-3x^2+6)/(x^6-5x^4+5x^2+4), x], "
                   "Verify -> False]; Lean`Private`$antiderivOverride = None; Null", "Null", 0);
    check_true("leanW[\"Certificate\"][\"PreCheck\"] === False");
    check_true("StringContainsQ[leanW[\"LeanSource\"], \"theorem mathilda_cert_1\"]");
}

/* Rows 5-7: Bezout, ideal membership, Solve. */
static void test_algebra(void) {
    assert_eval_eq("leanB = Certify[PolynomialExtendedGCD[x^3-1, x^2-1, x], Verify -> False]; "
                   "Null", "Null", 0);
    check_true("leanB[\"Kind\"] === \"Bezout\" && leanB[\"TrustedSteps\"] === {}");
    check_true("leanB[\"Certificate\"][\"Cofactors\"] === {1, -x}");
    check_true("leanB[\"Certificate\"][\"PreCheck\"] === True");
    check_true("StringCount[leanB[\"LeanTheorem\"], \"\xE2\x88\xA7\"] === 2");
    check_true("Certify[PolynomialReduce[x^2 y + x y^2 + y^2, {x y - 1, y^2 - 1}, {x, y}], "
               "Verify -> False][\"Kind\"] === \"IdealMembership\"");
    /* Ideal membership: PolynomialReduce leaves a remainder, so the cofactors
     * come from the degree-bounded ansatz. */
    assert_eval_eq("leanI = Certify[Implies[x^2+y^2 == 1 && x == y, 2 y^2 == 1], "
                   "Verify -> False]; Null", "Null", 0);
    check_true("Expand[leanI[\"Certificate\"][\"Cofactors\"] . {x^2+y^2-1, x-y} - "
               "(2 y^2 - 1)] === 0");
    check_true("StringContainsQ[leanI[\"ProofScript\"], \"rw [mth_k1, mth_k2] at mth_hc\"]");
    check_true("Lean`Private`cofactors[x^2 - y^2, {x - y}, {x, y}] === {x + y}");
    assert_eval_eq("leanQ = Certify[Solve[x^3-6x^2+11x-6 == 0, x], Verify -> False]; Null",
                   "Null", 0);
    check_true("leanQ[\"Kind\"] === \"Roots\" && leanQ[\"Certificate\"][\"Completeness\"]");
    check_true("StringContainsQ[leanQ[\"LeanTheorem\"], \"\xE2\x86\x94\"]");
    check_true("Certify[Solve[(x-1)^2 (2x+1) == 0, x], Verify -> False][\"Certificate\"]"
               "[\"Multiplicities\"] === {1, 2}");
    assert_eval_eq("leanL = Certify[Solve[{2x+y == 3, x+3y == 5}, {x, y}], Verify -> False]; "
                   "Null", "Null", 0);
    check_true("leanL[\"Certificate\"][\"Completeness\"] === False");
    check_true("leanL[\"ProofScript\"] === \"by\\n  decide +kernel\"");
    check_true("Length[leanL[\"TrustedSteps\"]] === 1");
}

/* Rows 8-10: matrices, checked by kernel evaluation of the prelude. */
static void test_matrices(void) {
    /* Every matrix statement states the shapes (matMul / matVec zipWith would
     * silently truncate ragged input). */
    check_str("Certify[Det[{{2,1},{1,1}}], Verify -> False][\"LeanTheorem\"]",
              "Mathilda.HasShape [[2, 1], [1, 1]] 2 2 \xE2\x88\xA7 "
              "Mathilda.det [[2, 1], [1, 1]] = 1");
    check_str("Certify[LinearSolve[{{2,1},{1,3}},{3,5}], Verify -> False][\"LeanTheorem\"]",
              "Mathilda.HasShape [[2, 1], [1, 3]] 2 2 \xE2\x88\xA7 "
              "Mathilda.HasLength [(4 / 5 : Rat), (7 / 5 : Rat)] 2 \xE2\x88\xA7 "
              "Mathilda.matVec [[2, 1], [1, 3]] [(4 / 5 : Rat), (7 / 5 : Rat)] = [3, 5]");
    check_true("StringStartsQ[Certify[Inverse[{{2,1},{1,1}}], Verify -> False]"
               "[\"LeanTheorem\"], \"Mathilda.HasShape [[2, 1], [1, 1]] 2 2 \xE2\x88\xA7 "
               "Mathilda.HasShape [[1, -1], [-1, 2]] 2 2 \xE2\x88\xA7 "
               "Mathilda.matMul [[2, 1], [1, 1]] [[1, -1], [-1, 2]] = "
               "Mathilda.idMat 2 \xE2\x88\xA7 \"]");
    check_true("StringCount[Certify[NullSpace[{{1,2,3},{2,4,6}}], Verify -> False]"
               "[\"LeanTheorem\"], \"Mathilda.HasLength\"] === 2");
    check_true("Certify[NullSpace[{{1,2,3},{2,4,6}}], Verify -> False][\"Kind\"] === "
               "\"NullSpace\"");
    check_true("Certify[Inverse[{{2,1},{1,1}}], Verify -> False][\"ProofScript\"] === "
               "\"by\\n  decide +kernel\"");
    /* The prelude placeholders never leak a symbol onto $ContextPath. */
    check_true("Context[matMul] === \"Global`\"");
}

/* Row 11: batches (generation only; Lean runs are in test_lean_verify.c). */
static void test_batch(void) {
    assert_eval_eq("leanBa = Certify[{Factor[x^6-1], Inverse[{{2,1},{1,1}}], Sin[x] == 1}, "
                   "Verify -> False]; Null", "Null", 0);
    check_true("Length[leanBa] === 3 && leanBa[[3]] === $Failed");
    check_true("StringContainsQ[leanBa[[2]][\"LeanSource\"], \"theorem mathilda_cert_1\"] && "
               "StringContainsQ[leanBa[[2]][\"LeanSource\"], \"theorem mathilda_cert_2\"]");
    check_true("leanBa[[1]][\"LeanSource\"] === leanBa[[2]][\"LeanSource\"]");
    check_true("Certify[{}] === {}");
    check_true("Certify[{x == x, Integrate[1/x, x]}, Engine -> \"Mathilda\"][[2]]"
               "[\"Verified\"] === Missing[\"Heuristic\", True]");
}

static void test_heuristic_engine(void) {
    /* The PossibleZeroQ false positive: the heuristic engine says True. The
     * claim is false (the difference is x/10^25); test_lean_verify.c checks
     * that Lean says False. */
    assert_eval_eq("leanH = Certify[(1+x)^20 == Expand[(1+x)^20] - 10^-25 x, "
                   "Engine -> \"Mathilda\"]; Null", "Null", 0);
    /* "Verified" is never True without Lean: the heuristic answer is carried
     * as Missing["Heuristic", True | False]. */
    check_true("leanH[\"Verified\"] === Missing[\"Heuristic\", True] && "
               "leanH[\"Tier\"] === \"Heuristic\" && leanH[\"Engine\"] === \"Mathilda\"");
    check_true("leanH[\"Domain\"] === Rationals");
    check_true("Certify[(x+1)^2 == x^2 + 1, Engine -> \"Mathilda\"][\"Verified\"] === "
               "Missing[\"Heuristic\", False]");
    check_true("StringContainsQ[ToString[leanH], \"Mathilda Heuristic\"]");
}

/* Row 12: the exact certificate behind Reduce's SOS stage. */
static void test_sos_certificate(void) {
#ifdef USE_LAPACK
    assert_eval_eq("leanSC = Reduce`SOSCertificate[x^4 + y^4 - x^2 y^2 + 1, {}, {x, y}]; "
                   "ListQ[leanSC]", "True", 0);
    /* The blocks {mult, basis, Gram} reproduce q exactly. */
    check_true("Expand[Total[(#1 (#2 . #3 . #2)) & @@@ leanSC] - "
               "(x^4 + y^4 - x^2 y^2 + 1)] === 0");
    check_true("Lean`Private`ldl[leanSC[[1, 3]]] =!= $Failed");
#else
    assert_eval_eq("Reduce`SOSCertificate[x^4 + y^4 - x^2 y^2 + 1, {}, {x, y}]", "$Failed", 0);
#endif
    assert_eval_eq("Reduce`SOSCertificate[x^2 - 1, {}, {x}]", "$Failed", 0);
    assert_eval_eq("Attributes[Reduce`SOSCertificate]", "{Protected}", 0);
    /* Exact LDL^T: PSD decided exactly, zero pivots allowed. */
    check_true("Lean`Private`ldl[{{4, 2}, {2, 2}}] === {{{1, 0}, {1/2, 1}}, {4, 1}}");
    check_true("Lean`Private`ldl[{{1, 2}, {2, 1}}] === $Failed");
    check_true("Lean`Private`ldl[{{0, 0}, {0, 1}}] === {{{1, 0}, {0, 1}}, {0, 1}}");
    check_true("Lean`Private`ldl[{{0, 1}, {1, 1}}] === $Failed");
}

/* Rows 13-14: nonnegativity and optimality certificates. */
static void test_sos_optimality(void) {
    assert_eval_eq("leanP = Certify[ForAll[{x,y}, x^4-2x^2+1+y^2 >= 0], Verify -> False]; Null",
                   "Null", 0);
    check_true("leanP[\"Kind\"] === \"Nonnegativity\" && leanP[\"Domain\"] === Reals");
    check_true("StringContainsQ[leanP[\"ProofScript\"], \":= Mathilda.sqnn\"]");
    check_true("StringContainsQ[leanP[\"ProofScript\"], \"have mth_e : \"]");
    /* Stated over every ordered field: no Rat -> Real transfer step. */
    check_true("leanP[\"TrustedSteps\"] === {}");
    /* The quadratic path never calls Reduce`SOSCertificate. */
    assert_eval_eq("Lean`Private`$sosCalls = 0; "
                   "leanQ2 = Certify[ForAll[{x,y}, x^2+y^2 >= 2 x y], Verify -> False]; "
                   "Lean`Private`$sosCalls", "0", 0);
    check_true("leanQ2[\"Certificate\"][\"Method\"] === \"Quadratic\"");
    check_true("leanQ2[\"Certificate\"][\"Squares\"] === {{1, x - y, 1}}");
    check_str("leanQ2[\"LeanTheorem\"]", "\xE2\x88\x80 {K : Type} [Lean.Grind.Field K] [LE K] [LT K] [Std.LawfulOrderLT K] "
              "[Std.IsLinearOrder K] [Lean.Grind.OrderedRing K], \xE2\x88\x80 x y : K, x ^ 2 + y ^ 2 \xE2\x89\xA5 2 * x * y");
    /* Putinar: a multiplier block uses the hypothesis 0 <= g. */
    check_true("StringContainsQ[Certify[ForAll[{x}, x >= 0, x^3 + x >= 0], Verify -> False]"
               "[\"ProofScript\"], \"Mathilda.mulnn Mathilda.sqnn mth_g1\"]");
    assert_eval_eq("leanM = Certify[Minimize[{x^2+y^2, x+y >= 1}, {x,y}], Verify -> False]; Null",
                   "Null", 0);
    check_true("leanM[\"Kind\"] === \"Optimality\" && leanM[\"Certificate\"][\"Multipliers\"] === {1}");
    check_true("leanM[\"Result\"] === {1/2, {x -> 1/2, y -> 1/2}}");
    check_true("leanM[\"Certificate\"][\"PreCheck\"] === True");
    check_true("leanM[\"Domain\"] === Reals && leanM[\"TrustedSteps\"] === {}");
    assert_eval_eq("leanLP = Certify[Maximize[{3x+2y, x+y<=4, x+3y<=6, x>=0, y>=0}, {x,y}], "
                   "Verify -> False]; Null", "Null", 0);
    check_true("leanLP[\"Certificate\"][\"Multipliers\"] === {3, 0, 0, 1} && "
               "leanLP[\"Certificate\"][\"LinearProgram\"] === True");
    check_true("StringContainsQ[leanLP[\"LeanTheorem\"], \"3 * x + 2 * y \xE2\x89\xA4 12\"]");
    check_true("StringContainsQ[Certify[Maximize[{x y, x+y == 2}, {x,y}], Verify -> False]"
               "[\"ProofScript\"], \"subst mth_v1\"]");
    /* (x^2 - 1)^2 + y^2 has two zeros: the exact Gram search, not the numeric engine. */
    check_true("Certify[Minimize[x^4-2x^2+1+y^2, {x,y}], Verify -> False][\"Certificate\"]"
               "[\"Method\"] === \"ExactGram\"");
    test_mute_messages();
    assert_eval_eq("Check[Certify[ForAll[x, x^3 >= 0], Verify -> False], \"nocert\"]",
                   "\"nocert\"", 0);
    assert_eval_eq("Certify[Minimize[{x, x > 0}, x], Verify -> False]", "$Failed", 0);
}

/* Row 15: Pratt certificates. */
static void test_pratt(void) {
    check_true("Lean`Private`pratt[274177] === {{3, 2, {{2, 1}}}, {7, 3, {{2, 1}, {3, 1}}}, "
               "{17, 3, {{2, 4}}}, {274177, 5, {{2, 8}, {3, 2}, {7, 1}, {17, 1}}}}");
    assert_eval_eq("leanFI = Certify[FactorInteger[2^64+1], Verify -> False]; Null", "Null", 0);
    check_true("StringCount[leanFI[\"LeanSource\"], \"theorem mathilda_cert_1_pratt_\"] === 2");
    check_true("StringStartsQ[leanFI[\"LeanTheorem\"], "
               "\"18446744073709551617 = 274177 * 67280421310721 \xE2\x88\xA7 \"]");
    check_true("leanFI[\"Kind\"] === \"Primality\" && Length[leanFI[\"TrustedSteps\"]] === 1");
    check_true("!StringContainsQ[leanFI[\"LeanSource\"], \"@NAME@\"]");
    check_true("Certify[PrimeQ[1000001], Verify -> False][\"Theorem\"] === "
               "HoldForm[!PrimeQ[1000001]]");
    /* Not prime: stated with the definition Mathilda.IsPrime, proved from a factor
     * (composites) or from 2 <= n (n < 2). */
    check_str("Certify[PrimeQ[1000001], Verify -> False][\"LeanTheorem\"]",
              "\xC2\xAC Mathilda.IsPrime 1000001");
    check_str("Certify[PrimeQ[1], Verify -> False][\"LeanTheorem\"]", "\xC2\xAC Mathilda.IsPrime 1");
    check_str("Certify[PrimeQ[2], Verify -> False][\"LeanTheorem\"]", "Mathilda.IsPrime 2");
    check_true("Certify[PrimeQ[0], Verify -> False][\"TrustedSteps\"] === {}");
    check_true("StringStartsQ[Certify[PrimeQ[1000003], Verify -> False][\"LeanTheorem\"], "
               "\"Mathilda.PrattCertifies 1000003 [\"]");
}

/* Row 16: closed-form sums. */
static void test_sums(void) {
    assert_eval_eq("leanSu = Certify[Sum[k 2^k, {k, 0, n}], Verify -> False]; Null", "Null", 0);
    check_true("leanSu[\"Kind\"] === \"ClosedFormSum\" && leanSu[\"Result\"] === 2 + 2^(1 + n) (-1 + n)");
    check_str("leanSu[\"LeanTheorem\"]",
              "\xE2\x88\x80 m1 : Nat, mathilda_cert_1_S m1 = 2 + 2 ^ m1 * 2 * ((m1 : Rat) - 1)");
    check_true("StringContainsQ[leanSu[\"LeanSource\"], "
               "\"def mathilda_cert_1_t (k1 : Nat) : Rat := 2 ^ k1 * (k1 : Rat)\"]");
    check_true("StringContainsQ[leanSu[\"ProofScript\"], \"induction m1 with\"]");
    assert_eval_eq("leanS3 = Certify[Sum[k^3, {k, 1, n}], Verify -> False]; Null", "Null", 0);
    check_true("leanS3[\"Input\"] === HoldForm[Sum[k^3, {k, 1, n}]]");
    check_true("StringContainsQ[leanS3[\"LeanSource\"], \"| m1 + 1 => mathilda_cert_1_S m1 + "
               "mathilda_cert_1_t (m1 + 1 + 1)\"]");
    /* A claimed closed form is what gets proved. */
    check_true("Certify[Sum[k^3,{k,1,n}] == n^2 (n+1)^2/4 + 1, Verify -> False]"
               "[\"Certificate\"][\"PreCheck\"] === False");
    test_mute_messages();
    assert_eval_eq("Certify[Sum[1/k, {k, 1, n}], Verify -> False]", "$Failed", 0);
}

/* Row 17: Tier 1 HasDerivAt generation (no Lean, no Mathlib needed). */
static void test_tier1_generation(void) {
    assert_eval_eq("leanT1 = Certify[Integrate[(1+x)/(1+x^2), x], Tier -> 1, "
                   "LeanProject -> \"/nonexistent\", Verify -> False]; Null", "Null", 0);
    check_true("leanT1[\"Tier\"] === 1 && leanT1[\"TrustedSteps\"] === {} && "
               "leanT1[\"Domain\"] === Reals");
    check_true("StringStartsQ[leanT1[\"LeanSource\"], \"import Mathlib\"]");
    check_true("StringStartsQ[leanT1[\"LeanTheorem\"], \"\xE2\x88\x80 x : Real, \"]");
    check_true("StringContainsQ[leanT1[\"LeanTheorem\"], \"HasDerivAt (fun x : Real => \"]");
    check_true("StringContainsQ[leanT1[\"ProofScript\"], \"refine mth_d.congr_deriv ?_\"]");
    check_true("leanT1[\"Hypotheses\"] === {HoldForm[1 + x^2 != 0]}");
    /* ArcTanh is written through Log (Mathlib has no Real.artanh). */
    check_true("StringContainsQ[Certify[Integrate[1/(x^2-1), x], Tier -> 1, LeanProject -> "
               "\"/nonexistent\", Verify -> False][\"LeanTheorem\"], \"Real.log (1 - x)\"]");
    test_mute_messages();
    assert_eval_eq("Check[Certify[Integrate[1/(x^2+1), x], Tier -> 1, LeanProject -> None], "
                   "\"nomathlib\"]", "\"nomathlib\"", 0);
}

/* Lean`Oracle: the server side of the lean/ tactics (plan section 6).  Answers
 * are plain data (numerals, v<i>, + - * / ^, lists, tuples), which the Lean
 * side parses with its own parser and never elaborates as Lean syntax. */
static void test_oracle(void) {
    check_str("Lean`Oracle[\"LinearCombination\", {2 v2^2, 1}, {{v1^2 + v2^2, 1}, {v1, v2}}, "
              "{v1, v2}]",
              "[1, -v1 - v2]");
    /* Every cofactor list really is a certificate: a - b == Sum c_i (l_i - r_i). */
    check_true("Expand[(2 v2^2 - 1) - {1, -v1 - v2} . {v1^2 + v2^2 - 1, v1 - v2}] === 0");
    check_str("Lean`Oracle[\"LinearCombination\", {v1, 2}, {{v1, 1}}, {v1}]",
              "MATHILDA-FAIL: the goal is not a polynomial combination of the hypotheses "
              "(no cofactors up to degree 4)");
    check_str("Lean`Oracle[\"LinearCombination\", {v1 + 1, 1 + v1}, {}, {v1}]",
              "[]");
    check_str("Lean`Oracle[\"SOS\", {2 v1 v2, v1^2 + v2^2}, {}, {v1, v2}]",
              "[(1, v1 - v2, 0)]");
    check_str("Lean`Oracle[\"SOS\", {0, v1^3}, {v1}, {v1}]",
              "[(1, v1, 1)]");
    check_str("Lean`Oracle[\"SOS\", {0, v1^4 - 2 v1^2 + 1 + v2^2}, {}, {v1, v2}]",
              "[(1, 1 - v1 ^ 2, 0), (1, v2, 0)]");
    check_str("Lean`Oracle[\"SOS\", {v1, 1}, {}, {v1}]",
              "MATHILDA-FAIL: no sum-of-squares certificate found");
    check_str("Lean`Oracle[\"Pratt\", 274177]",
              "[(3, 2, [(2, 1)]), (7, 3, [(2, 1), (3, 1)]), (17, 3, [(2, 4)]), "
              "(274177, 5, [(2, 8), (3, 2), (7, 1), (17, 1)])]");
    check_str("Lean`Oracle[\"Pratt\", 274179]", "MATHILDA-FAIL: 274179 is not prime");
    check_true("StringStartsQ[Lean`Oracle[\"Pratt\", 2], \"MATHILDA-FAIL:\"]");
    check_str("Lean`Oracle[\"Eval\", \"Factor[x^4 - 1]\"]",
              "(x - 1) * (1 + x) * (1 + x ^ 2)");
    check_true("StringStartsQ[Lean`Oracle[\"Eval\", \"Factor[\"], \"MATHILDA-FAIL:\"]");
    check_true("StringStartsQ[Lean`Oracle[\"Nope\", 1], \"MATHILDA-FAIL: unsupported query\"]");
    assert_eval_eq("Context[Oracle]", "\"Lean`\"", 0);
}

/* Soundness hardening (generation side; the Lean verdicts are in
 * test_lean_verify.c). */
static void test_hardening(void) {
    test_mute_messages();
    /* Method: the reviewer's macro injection and other command-level text are
     * refused before any file is written; a plain multi-line script is
     * indented under `by`. */
    assert_eval_eq("Certify[{x + 1 == 1 + x, x == x + 1}, Method -> \"intros\\n  grind\\n"
                   "set_option hygiene false in macro_rules | `(command| set_option "
                   "maxHeartbeats $n:num in $_:command) => `(def mathilda_cert_2 : Nat := 0)\"]",
                   "$Failed", 0);
    assert_eval_eq("Check[Certify[x == x, Method -> \"grind\\nopen Nat in\\ntheorem t : True "
                   ":= trivial\", Verify -> False], \"method\"]", "\"method\"", 0);
    assert_eval_eq("Certify[x == x, Method -> \"run_tac pure ()\", Verify -> False]",
                   "$Failed", 0);
    assert_eval_eq("Certify[x == x, Method -> \"grind -- \\u00bb\", Verify -> False]",
                   "$Failed", 0);
    assert_eval_eq("Certify[x == x, Method -> 3, Verify -> False]", "$Failed", 0);
    check_str("Certify[x == x, Method -> \"intros\\ngrind\", Verify -> False][\"ProofScript\"]",
              "by\n  intros\n  grind");
    check_true("Certify[x == x, Method -> \"intros\\ngrind\", Verify -> False]"
               "[\"LeanSource\"] // StringContainsQ[\"theorem mathilda_cert_1_footer :\\n    "
               "\\[ForAll] x : Int, x = x :=\\n  mathilda_cert_1\\n"
               "#print axioms mathilda_cert_1_footer\"]");
    /* Every generated file disables auto-bound implicits. */
    check_true("StringContainsQ[Certify[x == x, Verify -> False][\"LeanSource\"], "
               "\"set_option autoImplicit false\"]");
    /* Closed statements: a free symbol is refused (it could resolve to a Lean
     * global such as Inhabited.default), a bound one is fine. */
    assert_eval_eq("Lean`Internal`ToLean[HoldForm[ForAll[x, x == default]], \"Int\", True, {}, "
                   "True]", "$Failed", 0);
    check_str("Lean`Internal`ToLean[HoldForm[ForAll[x, x == x]], \"Int\", True, {}, True]",
              "\xE2\x88\x80 x : Int, x = x");
    check_str("Lean`Internal`ToLean[HoldForm[x == default], \"Int\", True, {}]",
              "x = default");
    assert_eval_eq("Certify[Solve[x^2 == 1 && default x == 0, x], Verify -> False]",
                   "$Failed", 0);
    /* Solve: every constraint is in the soundness statement, and every
     * denominator is stated nonzero at the solution. */
    check_true("StringContainsQ[Certify[Solve[x^2 == 1 && x > 0, x], Verify -> False]"
               "[\"LeanTheorem\"], \"(1 : Rat) > 0\"]");
    check_true("StringContainsQ[Certify[Solve[x^2 == 1 && x != 1, x], Verify -> False]"
               "[\"LeanTheorem\"], \"(-1 : Rat) \xE2\x89\xA0 1\"]");
    check_true("Quiet[Certify[Solve[(x-1)/(x^2-1) == (x-1)/(2 x - 2), x], Verify -> False]]"
               "[\"Certificate\"][\"Denominators\"] === {-1 + x^2, -2 + 2 x}");
    check_true("StringStartsQ[Quiet[Certify[Solve[(x-1)/(x^2-1) == (x-1)/(2 x - 2), x], "
               "Verify -> False]][\"LeanTheorem\"], \"(1 ^ 2 - 1 : Rat) \xE2\x89\xA0 0 "
               "\xE2\x88\xA7 (2 * 1 - 2 : Rat) \xE2\x89\xA0 0 \xE2\x88\xA7 \"]");
    assert_eval_eq("Certify[Solve[x^2 == 1 || x == 3, x], Verify -> False]", "$Failed", 0);
    check_true("Certify[Solve[x^2 == 1 && x > 0, x], Verify -> False][\"Certificate\"]"
               "[\"Completeness\"] === False");
    /* Sums whose term (or closed form) has a denominator in the variable. */
    assert_eval_eq("Certify[Sum[(k^2 - k)/(k (k - 1)), {k, 1, n}] == n - 1, Verify -> False]",
                   "$Failed", 0);
    assert_eval_eq("Certify[Sum[1/(k (k + 1)), {k, 1, n}], Verify -> False]", "$Failed", 0);
    /* The cache key is the theorem text, not a hash of it. */
    check_true("ListQ[Lean`Private`cacheKey[<|\"Statement\" -> \"True\", \"Proof\" -> "
               "\"by trivial\", \"TierUsed\" -> 0|>, 0, \"lean\", None]]");
}

static void test_attributes(void) {
    assert_eval_eq("Attributes[Certify]", "{HoldFirst, Protected}", 0);
    assert_eval_eq("Attributes[ProofObject]", "{HoldFirst, Protected}", 0);
    assert_eval_eq("Attributes[LeanVerify]", "{Protected}", 0);
    assert_eval_eq("Context[Certify]", "\"Lean`\"", 0);
}


/* Third review round (v0.6520): no Lean needed. */

/* Finding A: LeanVerify's allowlist grammar (src/lean_check.c). */
static void test_review3_checksource(void) {
    test_mute_messages();
#define CK(src) "Lean`Internal`CheckSource[\"" src "\", False, False, None]"
#define CKM(src) "Lean`Internal`CheckSource[\"" src "\", True, False, None]"
    /* Accepted: theorem / example / namespace / section / open / variable /
     * universe / numeric set_option / #print axioms, comments, strings. */
    check_true(CK("theorem t : 1 + 1 = 2 := rfl\\n") "[\"Theorems\"] === {\"t\"}");
    check_true(CK("namespace A\\nexample : True := trivial\\nsection\\nvariable (n : Nat)\\n"
                  "theorem u : n = n := rfl\\nend\\nend A\\n") "[\"Theorems\"] === "
               "{\"mathilda_example_1\", \"A.u\"}");
    check_true("StringStartsQ[" CK("example : True := trivial\\n") "[\"Source\"], "
               "\"noncomputable def _root_.mathilda_example_1 : True\"]");
    check_true(CK("universe u\\nopen Nat in\\nset_option maxHeartbeats 1000 in\\n"
                  "theorem t : True := by\\n  trivial -- sorry in a comment\\n"
                  "/- def /- nested -/ -/\\n#print axioms t\\n") "[\"OK\"]");
    check_true(CK("theorem t : \\\"#eval\\\" = \\\"#eval\\\" := rfl\\n") "[\"OK\"]");
    check_true(CKM("import Mathlib\\nimport Mathlib.Tactic\\nlemma t : True := trivial\\n") "[\"OK\"]");
    /* The trusted prelude, byte-identical, is skipped (it has defs). */
    check_true("Lean`Internal`CheckSource[Lean`Private`prelude[] <> \"\\ntheorem t : True := "
               "trivial\\n\", False, False, Lean`Private`prelude[]][\"OK\"]");
    check_true("Lean`Internal`CheckSource[StringReplace[Lean`Private`prelude[], \"powMod\" -> "
               "\"powMad\"], False, False, Lean`Private`prelude[]][\"Kind\"] === \"cmd\"");
    /* Imports: none in Tier 0; only Mathlib with a project. */
    check_true(CK("import Lean\\ntheorem t : True := trivial\\n") "[\"Kind\"] === \"import\"");
    check_true(CKM("import Lean\\ntheorem t : True := trivial\\n") "[\"Kind\"] === \"import\"");
    check_true(CK("import Mathlib\\ntheorem t : True := trivial\\n") "[\"Kind\"] === \"import\"");
    check_true(CKM("import Mathilda\\n") "[\"Kind\"] === \"import\"");
    check_true(CKM("import Batteries\\n") "[\"Kind\"] === \"import\"");
    check_true(CKM("import Mathlib Lean\\n") "[\"Kind\"] === \"import\"");
    /* Commands outside the grammar. */
    check_true(CK("def f : Nat := 1\\n") "[\"What\"] === \"def\"");
    check_true(CK("abbrev f : Nat := 1\\n") "[\"What\"] === \"abbrev\"");
    check_true(CK("instance : Inhabited Nat := \\u27e81\\u27e9\\n") "[\"What\"] === \"instance\"");
    check_true(CK("axiom bad : False\\n") "[\"What\"] === \"axiom\"");
    check_true(CK("macro \\\"x\\\" : term => pure 1\\n") "[\"What\"] === \"macro\"");
    check_true(CK("attribute [simp] Nat.add_comm\\n") "[\"What\"] === \"attribute\"");
    check_true(CK("#check 1\\n") "[\"What\"] === \"#check\"");
    check_true(CK("module\\n") "[\"What\"] === \"module\"");
    check_true(CK("lemma t : True := trivial\\n") "[\"What\"] === \"lemma\"");
    check_true(CK("set_option debug.skipKernelTC true\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("set_option maxHeartbeats true\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("open Lean.Elab\\n") "[\"What\"] === \"Elab\"");
    check_true(CK("theorem t : True := trivial\\n:= rfl\\n") "[\"What\"] === \":\"");
    check_true(CK("  theorem t : True := trivial\\n") "[\"Kind\"] === \"cmd\"");
    /* Escapes inside an accepted command, also on a continuation line. */
    check_true(CK("theorem t : 1 = 1 := by\\n  native_decide\\n") "[\"What\"] === \"native_decide\"");
    check_true(CK("theorem t : 1 = 1 := by\\n  decide +native\\n") "[\"What\"] === \"native\"");
    check_true(CK("theorem t : 1 = 1 := by\\n  sorry\\n") "[\"What\"] === \"sorry\"");
    check_true(CK("theorem t : 1 = 1 := by\\n  run_tac pure ()\\n") "[\"What\"] === \"run_tac\"");
    check_true(CK("theorem t : 1 = 1 := by_elab pure ()\\n") "[\"What\"] === \"by_elab\"");
    check_true(CK("theorem t : 1 = 1 := rfl\\n  #eval 1\\n") "[\"What\"] === \"#eval\"");
    check_true(CK("theorem t : 1 = 1 := rfl\\n  def f := 1\\n") "[\"What\"] === \"def\"");
    check_true(CK("theorem t : True := by\\n  set_option debug.x true in trivial\\n")
               "[\"What\"] === \"set_option\"");
    check_true(CK("theorem t : True := eval% trivial\\n") "[\"What\"] === \"eval%\"");
    check_true(CK("theorem t : True := (by exact (unsafeBaseIO (pure trivial)))\\n")
               "[\"What\"] === \"unsafeBaseIO\"");
    check_true(CK("@[simp] theorem t : True := trivial\\n") "[\"What\"] === \"@[ (an attribute)\"");
    check_true(CK("theorem t : `x = `x := rfl\\n") "[\"Kind\"] === \"token\"");
    check_true(CK("theorem t : s!\\\"{1}\\\" = \\\"1\\\" := rfl\\n") "[\"Kind\"] === \"token\"");
    /* Lexical ambiguity. */
    check_true(CK("theorem t : True := trivial /- open\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : \\\"x = 1 := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : True := trivial\\n\\tdef f := 1\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : True := trivial\\r\\n") "[\"Kind\"] === \"lex\"");
    check_true("Lean`Internal`CheckSource[\"theorem t : \\[LeftGuillemet]sorry\\[RightGuillemet] "
               "= 1 := rfl\\n\", False, False, None][\"Kind\"] === \"token\"");
    /* A char literal holding a quote does not open a string. */
    check_true(CK("theorem t : '\\\"' = '\\\"' := rfl\\n") "[\"OK\"]");
    /* Certify mode: the generator auxiliaries only. */
    check_true("Lean`Internal`CheckSource[\"def mathilda_cert_3_S : Nat := 1\\n\", False, True, "
               "None][\"OK\"]");
    check_true("Lean`Internal`CheckSource[\"def mathilda_cert_x : Nat := 1\\n\", False, True, "
               "None][\"OK\"] === False");
#undef CK
#undef CKM
    /* Every generated Certify file passes the check (the internal assertion). */
    check_true("AllTrue[{Certify[Factor[x^10 - 1], Verify -> False], "
               "Certify[Sum[k 2^k, {k, 0, n}], Verify -> False], "
               "Certify[FactorInteger[2^64 + 1], Verify -> False], "
               "Certify[Integrate[1/(x^2 + 2), x], Verify -> False], "
               "Certify[Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}], "
               "Verify -> False], Certify[Inverse[{{2, 1}, {1, 1}}], Verify -> False], "
               "Certify[Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1], Verify -> False], "
               "Certify[Integrate[(1 + x)/(1 + x^2), x], Tier -> 1, LeanProject -> \"/nonexistent\", "
               "Verify -> False]}, "
               "TrueQ[Lean`Private`checkSource[#[\"LeanSource\"], #[\"Tier\"] === 1, True][\"OK\"]] &]");
}

/* Findings B, C, D. */
static void test_review3_claims(void) {
    test_mute_messages();
    /* B: an order relation makes only its own variables real. */
    assert_eval_eq("Check[Certify[x == 7, Assumptions -> x^2 + 4 == 0 && y > 0, Verify -> False], "
                   "\"realvars\"]", "\"realvars\"", 0);
    assert_eval_eq("Check[Certify[x == 7, Assumptions -> x^2 + 4 == 0 && y >= y, Verify -> False], "
                   "\"realvars\"]", "\"realvars\"", 0);
    assert_eval_eq("Check[Certify[x^2 == -4, Assumptions -> y > 0, Verify -> False], "
                   "\"realvars\"]", "\"realvars\"", 0);
    check_true("Certify[x == 7, Assumptions -> x^2 + 4 == 0 && Element[x, Reals] && y > 0, "
               "Verify -> False][\"Domain\"] === Reals");
    check_true("Certify[x == y, Assumptions -> {x == y, Element[{x, y}, Reals]}, Verify -> False]"
               "[\"Domain\"] === Reals");
    check_true("Certify[x == 1, Assumptions -> 0 < x < 2, Verify -> False][\"Domain\"] === Reals");
    check_true("Certify[x == 7, Assumptions -> x^2 + 4 == 0, Verify -> False][\"Domain\"] === "
               "Complexes");
    /* C: a computation inside a claim is refused, not evaluated. */
    const char* evals[] = {
        "Simplify[(x^2 - 1)/(x - 1)] == x + 1", "Refine[Sqrt[x^2], x < 0] == -x",
        "Integrate[Sin[x], {x, 0, Pi}] == 2", "Integrate[x, {x, 0, 1}] == 1/2",
        "Exp[Log[x]] == x", "Mod[7, 3] == 1", "Limit[x, x -> 0] == 0", "Sqrt[4] == 2",
        "Sin[0] == 0", "Det[{{a, b}, {c, d}}] == a d - b c", "N[1/3] == 1/3",
        "Simplify[Sin[x]^2 + Cos[x]^2] == 1", "PowerExpand[Sqrt[x^2]] == x",
        "Simplify[x^0] == 1", "f[x] == f[x]", "D[Sin[x], x] == Cos[x]",
        "Implies[x == 1, Simplify[x^2] == 1]", "ForAll[{x}, Abs[x]^2 >= 0]",
        "Expand[Integrate[x, x]]", NULL };
    for (int i = 0; evals[i]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "Check[Certify[%s, Verify -> False], \"eval\"]", evals[i]);
        assert_eval_eq(buf, "\"eval\"", 0);
    }
    /* ... and the value-preserving rational-function heads are read through. */
    check_str("Certify[Factor[x^2 - 1] == (x - 1) (x + 1), Verify -> False][\"LeanTheorem\"]",
              "\xE2\x88\x80 x : Int, x ^ 2 - 1 = (x - 1) * (1 + x)");
    check_true("With[{c = Certify[Factor[x^2 - 1] == (x - 1) (x + 1), Verify -> False]}, "
               "c[\"Kind\"] === \"Factorization\" && "
               "MemberQ[c[\"TrustedSteps\"], \"irreducibility of the factors (not certified)\"]]");
    check_str("Certify[Expand[(x + 1)^2] == x^2 + 2 x + 1, Verify -> False][\"LeanTheorem\"]",
              "\xE2\x88\x80 x : Int, (1 + x) ^ 2 = 1 + 2 * x + x ^ 2");
    check_true("Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, Verify -> False]"
               "[\"Input\"] === HoldForm[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x]");
    /* x/x is not cancelled to 1 before Lean sees it. */
    check_true("With[{c = Certify[x/x == 1, Verify -> False]}, c[\"LeanTheorem\"] === "
               "\"\\[ForAll] x : Int, x = x\" && c[\"Hypotheses\"] === {HoldForm[x != 0]}]");
    /* D[F, x] of a polynomial: Mathilda's derivative is a trusted step. */
    check_true("AnyTrue[Certify[D[x^3, x] == 3 x^2, Verify -> False][\"TrustedSteps\"], "
               "StringStartsQ[#, \"Mathilda's differentiation\"] &]");
    /* Integrate[f, x] == F checks the CLAIMED F. */
    assert_eval_eq("leanIC = Certify[Integrate[x^2, x] == x^3/2, Verify -> False]; Null", "Null", 0);
    check_true("leanIC[\"Kind\"] === \"Antiderivative\" && leanIC[\"Result\"] === x^3/2 && "
               "leanIC[\"Certificate\"][\"PreCheck\"] === False");
    check_true("Certify[x^3/3 == Integrate[x^2, x], Verify -> False][\"Certificate\"][\"PreCheck\"]"
               " === True");
    /* D: Sum's upper bound must be a variable; ignored options are reported. */
    assert_eval_eq("Check[Certify[Sum[1, {k, 1, Infinity}] == Infinity, Verify -> False], "
                   "\"nocert\"]", "\"nocert\"", 0);
    assert_eval_eq("Check[Certify[Sum[1/2^k, {k, 1, Infinity}], Verify -> False], \"nocert\"]",
                   "\"nocert\"", 0);
    assert_eval_eq("Check[Certify[Sum[k, {k, 1, n}] == n (n + 1)/2, Assumptions -> n < 0, "
                   "Verify -> False], \"optignored\"]", "\"optignored\"", 0);
    assert_eval_eq("Check[Certify[Inverse[{{2, 1}, {1, 1}}], Domain -> Integers, Verify -> False], "
                   "\"optignored\"]", "\"optignored\"", 0);
    check_true("Certify[Sum[k, {k, 1, n}] == n (n + 1)/2, Verify -> False][\"Kind\"] === "
               "\"ClosedFormSum\"");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_identity_source);
    TEST(test_computations);
    TEST(test_rational_identity);
    TEST(test_options);
    TEST(test_errors);
    TEST(test_antiderivative);
    TEST(test_algebra);
    TEST(test_matrices);
    TEST(test_batch);
    TEST(test_heuristic_engine);
    TEST(test_sos_certificate);
    TEST(test_sos_optimality);
    TEST(test_pratt);
    TEST(test_sums);
    TEST(test_tier1_generation);
    TEST(test_oracle);
    TEST(test_hardening);
    TEST(test_attributes);
    TEST(test_review3_checksource);
    TEST(test_review3_claims);
    printf("All Lean certificate generation tests passed.\n");
    return 0;
}
