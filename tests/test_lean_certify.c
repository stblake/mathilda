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

#include "lean_check_corpus.h"

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
    check_true("Options[Certify][[1]] === (\"Engine\" -> \"Lean\")");
    /* A symbol holding a claim has a value: refused (v0.6521), unless the
     * caller asks for the value explicitly with Evaluate. */
    assert_eval_eq("leanC = (x + 1)^2 == x^2 + 2 x + 1; Check[Certify[leanC, Verify -> False], "
                   "\"value\"]", "\"value\"", 0);
    check_true("Certify[Evaluate[leanC], Verify -> False][\"Kind\"] === \"Identity\"");
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
    /* Domain -> Reals is the Tier 0 ordered-field form (v0.6522; it used to
     * need Mathlib); an unknown Domain is refused. */
    check_true("Certify[x == x, Domain -> Reals, LeanProject -> None, Verify -> False]"
               "[\"Tier\"] === 0");
    assert_eval_eq("Check[Certify[x == x, Domain -> Algebraics, Verify -> False], \"domain\"]",
                   "\"domain\"", 0);
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

/* Lean`Internal`Oracle: the server side of the lean/ tactics (plan section 6).  Answers
 * are plain data (numerals, v<i>, + - * / ^, lists, tuples), which the Lean
 * side parses with its own parser and never elaborates as Lean syntax. */
static void test_oracle(void) {
    check_str("Lean`Internal`Oracle[\"LinearCombination\", {2 v2^2, 1}, {{v1^2 + v2^2, 1}, {v1, v2}}, "
              "{v1, v2}]",
              "[1, -v1 - v2]");
    /* Every cofactor list really is a certificate: a - b == Sum c_i (l_i - r_i). */
    check_true("Expand[(2 v2^2 - 1) - {1, -v1 - v2} . {v1^2 + v2^2 - 1, v1 - v2}] === 0");
    check_str("Lean`Internal`Oracle[\"LinearCombination\", {v1, 2}, {{v1, 1}}, {v1}]",
              "MATHILDA-FAIL: the goal is not a polynomial combination of the hypotheses "
              "(no cofactors up to degree 4)");
    check_str("Lean`Internal`Oracle[\"LinearCombination\", {v1 + 1, 1 + v1}, {}, {v1}]",
              "[]");
    check_str("Lean`Internal`Oracle[\"SOS\", {2 v1 v2, v1^2 + v2^2}, {}, {v1, v2}]",
              "[(1, v1 - v2, 0)]");
    check_str("Lean`Internal`Oracle[\"SOS\", {0, v1^3}, {v1}, {v1}]",
              "[(1, v1, 1)]");
    check_str("Lean`Internal`Oracle[\"SOS\", {0, v1^4 - 2 v1^2 + 1 + v2^2}, {}, {v1, v2}]",
              "[(1, 1 - v1 ^ 2, 0), (1, v2, 0)]");
    check_str("Lean`Internal`Oracle[\"SOS\", {v1, 1}, {}, {v1}]",
              "MATHILDA-FAIL: no sum-of-squares certificate found");
    check_str("Lean`Internal`Oracle[\"Pratt\", 274177]",
              "[(3, 2, [(2, 1)]), (7, 3, [(2, 1), (3, 1)]), (17, 3, [(2, 4)]), "
              "(274177, 5, [(2, 8), (3, 2), (7, 1), (17, 1)])]");
    check_str("Lean`Internal`Oracle[\"Pratt\", 274179]", "MATHILDA-FAIL: 274179 is not prime");
    check_true("StringStartsQ[Lean`Internal`Oracle[\"Pratt\", 2], \"MATHILDA-FAIL:\"]");
    check_str("Lean`Internal`Oracle[\"Eval\", \"Factor[x^4 - 1]\"]",
              "(x - 1) * (1 + x) * (1 + x ^ 2)");
    check_true("StringStartsQ[Lean`Internal`Oracle[\"Eval\", \"Factor[\"], \"MATHILDA-FAIL:\"]");
    check_true("StringStartsQ[Lean`Internal`Oracle[\"Nope\", 1], \"MATHILDA-FAIL: unsupported query\"]");
    /* The oracle is reached by its full name only: Lean`Internal` is not on
     * $ContextPath, so the short name Oracle stays the user's. */
    assert_eval_eq("Context[Oracle]", "\"Global`\"", 0);
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
    assert_eval_eq("Attributes[Certify]", "{HoldAll, Protected}", 0);
    assert_eval_eq("Attributes[ProofObject]", "{HoldAll, Protected}", 0);
    assert_eval_eq("Attributes[LeanVerify]", "{Protected}", 0);
    assert_eval_eq("Context[Certify]", "\"Lean`\"", 0);
}


/* Third review round (v0.6520): no Lean needed. */

/* Finding A: LeanVerify's allowlist grammar (src/lean_check.c). */
static void test_review3_checksource(void) {
    test_mute_messages();
#define CK(src) "Lean`Internal`CheckSource[\"" src "\", False, False, None]"
#define CKM(src) "Lean`Internal`CheckSource[\"" src "\", True, False, None]"
    /* Accepted: theorem / example / section / variable / universe / numeric
     * set_option (namespace and open are refused since v0.6522). */
    check_true(CK("theorem t : 1 + 1 = 2 := rfl\\n") "[\"Theorems\"] === {\"t\"}");
    check_true(CK("example : True := trivial\\nsection\\nvariable (n : Nat)\\n"
                  "theorem u : n = n := rfl\\nend\\n") "[\"Theorems\"] === "
               "{\"mathilda_example_1\", \"u\"}");
    check_true("StringStartsQ[" CK("example : True := trivial\\n") "[\"Source\"], "
               "\"noncomputable def _root_.mathilda_example_1 : True\"]");
    check_true(CK("universe u\\nset_option maxHeartbeats 1000 in\\n"
                  "theorem t : True := by\\n  trivial\\n") "[\"OK\"]");
    check_true(CKM("import Mathlib\\nimport Mathlib.Tactic\\nlemma t : True := trivial\\n") "[\"OK\"]");
    /* The trusted prelude, byte-identical, is skipped (it has defs). */
    check_true("Lean`Internal`CheckSource[Lean`Private`prelude[] <> \"\\ntheorem t : True := "
               "trivial\\n\", False, False, Lean`Private`prelude[]][\"OK\"]");
    check_true("Lean`Internal`CheckSource[StringReplace[Lean`Private`prelude[], \"powMod\" -> "
               "\"powMad\"], False, False, Lean`Private`prelude[]][\"Kind\"] === \"lex\"");
    /* ... and must itself be free of strings and comments (v0.6521). */
    check_true("Lean`Internal`CheckSource[\"theorem t : True := trivial\\n\", False, False, "
               "\"/- a comment -/\\n\"][\"Kind\"] === \"lex\"");
    check_true("!StringContainsQ[Lean`Private`prelude[], \"/-\" | \"--\" | \"\\\"\"]");
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
    check_true(CK("attribute [simp] Nat.add_comm\\n") "[\"What\"] === \"attribute\"");
    check_true(CK("module\\n") "[\"What\"] === \"module\"");
    check_true(CK("lemma t : True := trivial\\n") "[\"What\"] === \"lemma\"");
    check_true(CK("set_option debug.skipKernelTC true\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("set_option maxHeartbeats true\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("open Lean.Elab\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("theorem t : True := trivial\\n:= rfl\\n") "[\"What\"] === \":\"");
    check_true(CK("  theorem t : True := trivial\\n") "[\"Kind\"] === \"cmd\"");
    /* open: refused in every form (v0.6522). */
    check_true(CK("open Nat in\\ntheorem t : True := trivial\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("open scoped Nat\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("open Nat hiding add\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("open Nat (add)\\n") "[\"Kind\"] === \"cmd\"");
    /* A starter glued to another character is another token. */
    check_true(CK("variable? (n : Nat)\\n") "[\"Kind\"] === \"cmd\"");
    check_true(CK("theorem' t : True := trivial\\n") "[\"Kind\"] === \"cmd\"");
    /* Escapes inside an accepted command, also on a continuation line. */
    check_true(CK("theorem t : 1 = 1 := by\\n  native_decide\\n") "[\"What\"] === \"native_decide\"");
    check_true(CK("theorem t : 1 = 1 := by\\n  decide +native\\n") "[\"What\"] === \"native\"");
    check_true(CK("theorem t : 1 = 1 := by\\n  sorry\\n") "[\"What\"] === \"sorry\"");
    check_true(CK("theorem t : 1 = 1 := by\\n  run_tac pure ()\\n") "[\"What\"] === \"run_tac\"");
    check_true(CK("theorem t : 1 = 1 := by_elab pure ()\\n") "[\"What\"] === \"by_elab\"");
    check_true(CK("theorem t : 1 = 1 := rfl\\n  def f := 1\\n") "[\"What\"] === \"def\"");
    check_true(CK("theorem t : True := by\\n  set_option debug.x true in trivial\\n")
               "[\"What\"] === \"set_option\"");
    check_true(CK("theorem t : True := (by exact (unsafeBaseIO (pure trivial)))\\n")
               "[\"What\"] === \"unsafeBaseIO\"");
    /* A keyword after a numeric literal (Lean splits 2theorem into 2, theorem). */
    check_true(CK("theorem t : True := trivial\\n  2theorem u : False := sorry\\n")
               "[\"What\"] === \"theorem\"");
    /* The minimal lexical surface: all of these are "lex" now. */
    check_true(CK("theorem t : 1 = 1 := rfl\\n  #eval 1\\n") "[\"What\"] === \"#eval\"");
    check_true(CK("#check 1\\n") "[\"What\"] === \"#check\"");
    check_true(CK("theorem t : True := eval% trivial\\n") "[\"What\"] === \"eval%\"");
    check_true(CK("@[simp] theorem t : True := trivial\\n") "[\"What\"] === \"@[ (an attribute)\"");
    check_true(CK("theorem t : `x = `x := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : s!\\\"{1}\\\" = \\\"1\\\" := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : \\\"#eval\\\" = \\\"#eval\\\" := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : True := trivial -- c\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : True := trivial /- open\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : True := trivial\\n\\tdef f := 1\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : True := trivial\\r\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : (fun x => x) $ True := trivial\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : 7 % 2 = 1 := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : (\\\\x => x) 1 = 1 := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true("Lean`Internal`CheckSource[\"theorem t : \\[LeftGuillemet]x\\[RightGuillemet] "
               "= 1 := rfl\\n\", False, False, None][\"Kind\"] === \"lex\"");
    check_true(CK("theorem t : 'a' = 'a' := rfl\\n") "[\"Kind\"] === \"lex\"");
    check_true(CK("theorem t (h' : True) : True := h'\\n") "[\"OK\"]");
    check_true(CK("theorem t : \xE2\x8A\x95 = 1 := rfl\\n") "[\"What\"] === \"the character U+2295\"");
    check_true("Lean`Internal`CheckSource[\"theorem t : \\[Alpha] = \\[Beta] := rfl\\n\", False, "
               "False, None][\"OK\"]");
    /* Certify mode: the generator auxiliaries and footer lines only. */
    check_true("Lean`Internal`CheckSource[\"def mathilda_cert_3_S : Nat := 1\\n\", False, True, "
               "None][\"OK\"]");
    check_true("Lean`Internal`CheckSource[\"def mathilda_cert_x : Nat := 1\\n\", False, True, "
               "None][\"OK\"] === False");
    check_true("Lean`Internal`CheckSource[\"theorem t : True := trivial\\n#print axioms "
               "mathilda_cert_2_footer\\n\", False, True, None][\"OK\"]");
    check_true("Lean`Internal`CheckSource[\"theorem t : True := trivial\\n#print axioms "
               "mathilda_cert_2_footer\\n\", False, False, None][\"Kind\"] === \"lex\"");
    check_true("Lean`Internal`CheckSource[\"theorem t : True := trivial\\n#print axioms "
               "mathilda_cert_2_footer x\\n\", False, True, None][\"Kind\"] === \"lex\"");
    check_true("Lean`Internal`CheckSource[\"theorem t : True := trivial\\n  #print axioms "
               "mathilda_cert_2_footer\\n\", False, True, None][\"Kind\"] === \"lex\"");
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

/* Fourth review round (v0.6521): no Lean needed. */

/* CheckSource[src, mathlib, False, None]["OK"] for the bytes src. */
static bool review4_ok(const char* src, bool mathlib) {
    char* lit = lean_mstring(src);
    size_t n = strlen(lit) + 96;
    char* in = malloc(n);
    if (!in) abort();
    snprintf(in, n, "TrueQ[Lean`Internal`CheckSource[%s, %s, False, None][\"OK\"]]", lit,
             mathlib ? "True" : "False");
    Expr* e = parse_expression(in);
    Expr* r = e ? evaluate(e) : NULL;
    bool ok = r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    if (r) expr_free(r);
    if (e) expr_free(e);
    free(in);
    free(lit);
    return ok;
}

/* Fifth review: `open X` and `namespace X` activate X's scoped extensions
 * (scoped command syntax included, e.g. Lean.PostprocessTraces), which the
 * parser-table-derived keyword list cannot see; both are refused outright. */
static void test_review5_scoped(void) {
    const char* srcs[] = {
        "open Lean.PostprocessTraces\ntheorem t : True := trivial\n",
        "open Nat\ntheorem c22 : succ 1 = 2 := rfl\n",
        "namespace Lean.PostprocessTraces\ntheorem t : True := trivial\nend Lean.PostprocessTraces\n",
        "namespace Foo\ntheorem bar : True := trivial\nend Foo\n",
        "theorem t : True := trivial\nopen Lean\ntheorem u : True := trivial\n",
        NULL };
    for (int i = 0; srcs[i]; i++) {
        ASSERT(!review4_ok(srcs[i], false));
        ASSERT(!review4_ok(srcs[i], true));
    }
    /* The qualified-name rewrites stay accepted. */
    ASSERT(review4_ok("theorem c22 : Nat.succ 1 = 2 := rfl\n", false));
    ASSERT(review4_ok("section Foo\ntheorem bar : True := trivial\nend Foo\n", false));
}

/* Findings 1-3: every reproducer hides a column-0 #eval inside what the
 * v0.6520 checker took for a string / character literal / comment; with the
 * minimal lexical surface each is refused before Lean runs, at the very
 * construct that did the hiding. */
static void test_review4_reproducers(void) {
    const char* xs[] = { REVIEW4_X1, REVIEW4_X2, REVIEW4_X3, REVIEW4_X4, NULL };
    for (int i = 0; xs[i]; i++) {
        ASSERT(!review4_ok(xs[i], false));
        ASSERT(!review4_ok(xs[i], true));
    }
    /* The constructs themselves, alone. */
    const char* hide[] = {
        "theorem t : True := trivial\n  1'\"'\n",              /* (a) quote after a digit */
        "theorem t : True := (fun (_ : Nat) \xE2\x86\xA6'x') 0\n", /* (a) after non-ASCII */
        "theorem t : True := trivial\n  x!'y'\n",             /* (b) identifier with ! */
        "theorem t : True := trivial\n  x?'y'\n",             /* (b) identifier with ? */
        "theorem t : True := trivial\n/--/ x\n-/\n",           /* (c) doc comment opener */
        "theorem t : True := trivial\n  h.2'x'\n",            /* field index then quote */
        NULL };
    for (int i = 0; hide[i]; i++) ASSERT(!review4_ok(hide[i], false));
    /* An identifier prime is still accepted. */
    ASSERT(review4_ok("theorem t (h' h'' : True) : True := h''\n", false));
    ASSERT(review4_ok("theorem t (x1' : Nat) : x1' = x1' := rfl\n", false));
}

/* Property: every corpus source is accepted; inserting any refused
 * construct at any of several pseudo-random character boundaries makes the
 * checker refuse it. */
static void test_review4_property(void) {
    unsigned long seed = 12345;
    int tried = 0;
    for (int i = 0; LEAN_CORPUS[i]; i++) {
        const char* src = LEAN_CORPUS[i];
        size_t n = strlen(src);
        if (!review4_ok(src, false)) {
            printf("FAIL: corpus source %d refused:\n%s", i, src);
            ASSERT(0);
        }
        ASSERT(review4_ok(src, true));
        for (int j = 0; LEAN_REFUSED[j]; j++) {
            size_t m = strlen(LEAN_REFUSED[j]);
            for (int k = 0; k < 3; k++) {
                size_t pos;
                do {
                    seed = seed * 6364136223846793005UL + 1442695040888963407UL;
                    pos = (size_t)((seed >> 33) % (n + 1));
                } while (pos < n && ((unsigned char)src[pos] & 0xC0) == 0x80);
                char* t = malloc(n + m + 1);
                if (!t) abort();
                memcpy(t, src, pos);
                memcpy(t + pos, LEAN_REFUSED[j], m);
                memcpy(t + pos + m, src + pos, n - pos + 1);
                if (review4_ok(t, false)) {
                    printf("FAIL: accepted after inserting \"%s\" at %zu:\n%s", LEAN_REFUSED[j],
                           pos, t);
                    ASSERT(0);
                }
                free(t);
                tried++;
            }
        }
    }
    ASSERT(tried > 1000);
    /* Strict UTF-8: overlong, surrogate, stray and truncated sequences, and
     * codepoints outside the allowlist (a 4-byte letter, a zero-width space). */
    const char* bad8[] = { "\xC0\x80", "\xE0\x80\xAF", "\xED\xA0\x80", "\xFF", "\x80",
                           "\xE2\x88", "\xF0\x9D\x94\xB8", "\xE2\x80\x8B", NULL };
    for (int i = 0; bad8[i]; i++) {
        char t[96];
        snprintf(t, sizeof t, "theorem t : True := trivial %s\n", bad8[i]);
        ASSERT(!review4_ok(t, false));
    }
    ASSERT(review4_ok("theorem t (\xCE\xB1 \xCE\xA9 : Nat) : \xCE\xB1 + \xCE\xA9 = \xCE\xA9 + \xCE\xB1 := "
                      "Nat.add_comm _ _\n", false));
}

/* Findings 4-5: Evaluate and symbols with values in a claim. */
static void test_review4_claims(void) {
    test_mute_messages();
    /* 4: Evaluate inside a claim is refused, not evaluated. */
    assert_eval_eq("Check[ProofObject[Evaluate[Simplify[(x^2 - 1)/(x - 1)]] == x + 1, "
                   "Verify -> False], \"eval\"]", "\"eval\"", 0);
    assert_eval_eq("Check[Certify[x + 1 == Evaluate[Simplify[(x^2 - 1)/(x - 1)]], "
                   "Verify -> False], \"eval\"]", "\"eval\"", 0);
    assert_eval_eq("Check[Certify[{Evaluate[Expand[(x + 1)^2]] == x^2 + 2 x + 1}, "
                   "Verify -> False], \"eval\"]", "\"eval\"", 0);
    assert_eval_eq("Check[Certify[Factor[Evaluate[x^2 - 1]], Verify -> False], \"eval\"]",
                   "\"eval\"", 0);
    /* ... the raw-fraction reader is HoldAllComplete as well. */
    check_true("MemberQ[Attributes[Lean`Private`rawH], HoldAllComplete] && "
               "MemberQ[Attributes[Lean`Private`cw], HoldAllComplete]");
    /* 5: a symbol with a value is refused (Certify::value), not substituted. */
    assert_eval_eq("leanX4 = 5; Check[ProofObject[leanX4^2 == 25, Verify -> False], "
                   "\"value\"]", "\"value\"", 0);
    assert_eval_eq("leanE4 = Simplify[(x^2 - 1)/(x - 1)]; Check[ProofObject[leanE4 == x + 1, "
                   "Verify -> False], \"value\"]", "\"value\"", 0);
    assert_eval_eq("leanG4 = 0; Check[Certify[(a^2 - leanG4^2)/(a - leanG4) == a + leanG4, "
                   "Verify -> False], \"value\"]", "\"value\"", 0);
    /* ... also in Assumptions, which are checked held before they evaluate. */
    assert_eval_eq("leanA4 = 0; Check[Certify[y == 7, Assumptions -> y leanA4 == 1, "
                   "Verify -> False], \"value\"]", "\"value\"", 0);
    assert_eval_eq("Check[Certify[y == 7, Assumptions -> Evaluate[y == 7], Verify -> False], "
                   "\"eval\"]", "\"eval\"", 0);
    /* An explicit top-level Evaluate certifies the value, as asked. */
    check_true("StringContainsQ[Certify[Evaluate[(leanX4 + x)^2 == x^2 + 10 x + 25], "
               "Verify -> False][\"LeanTheorem\"], \"(5 + x) ^ 2\"]");
    /* Constants are not values. */
    check_true("Head[Certify[x + 0 == x, Verify -> False]] === ProofObject");
    /* Cleared, the symbol is a variable again. */
    check_true("Clear[leanX4]; StringContainsQ[Certify[(leanX4 + 1)^2 == leanX4^2 + 2 leanX4 + 1, "
               "Verify -> False][\"LeanTheorem\"], \"(1 + leanX4) ^ 2\"]");
    /* Method scripts satisfy the same source rules. */
    assert_eval_eq("Check[Certify[x == x, Method -> \"simp -- c\", Verify -> False], \"m\"]",
                   "\"m\"", 0);
    assert_eval_eq("Check[Certify[x == x, Method -> \"exact (x 2theorem)\", Verify -> False], "
                   "\"m\"]", "\"m\"", 0);
    assert_eval_eq("Check[Certify[x == x, Method -> \"show 7 % 2 = 1\", Verify -> False], \"m\"]",
                   "\"m\"", 0);
    check_str("Certify[x == x, Method -> \"intro x\\nrfl\", Verify -> False][\"ProofScript\"]",
              "by\n  intro x\n  rfl");
    /* Every generated Certify file passes the new check (the internal
     * assertion): no comment in the header, the comment-free prelude, the
     * footer lines. */
    check_true("!StringContainsQ[Certify[x == x, Verify -> False][\"LeanSource\"], \"/-\"]");
}

/* Sixth review (v0.6522): no Lean needed. */

#define BIND_CHAR0 "\xE2\x88\x80 {K : Type} [Lean.Grind.Field K] [Lean.Grind.IsCharP K 0], "
#define BIND_ORD "\xE2\x88\x80 {K : Type} [Lean.Grind.Field K] [LE K] [LT K] " \
    "[Std.LawfulOrderLT K] [Std.IsLinearOrder K] [Lean.Grind.OrderedRing K], "

/* Finding 1: every Domain value has one meaning, any other is refused. */
static void test_review6_domain(void) {
    test_mute_messages();
    /* The reproducer: Complexes with a hypothesis was stated over Rat. */
    check_str("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Complexes, Verify -> False]"
              "[\"LeanTheorem\"]", BIND_CHAR0 "\xE2\x88\x80 x : K, x ^ 2 = -4 \xE2\x86\x92 x = 7");
    check_true("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Complexes, Verify -> False]"
               "[\"Domain\"] === Complexes");
    /* Complexes without hypotheses: the generic form too. */
    check_str("Certify[(x + 1)^2 == x^2 + 2 x + 1, Domain -> Complexes, Verify -> False]"
              "[\"LeanTheorem\"]",
              BIND_CHAR0 "\xE2\x88\x80 x : K, (1 + x) ^ 2 = 1 + 2 * x + x ^ 2");
    /* Complexes: an order relation needs every variable real. */
    assert_eval_eq("Check[Certify[x == 7, Assumptions -> x^2 + 4 == 0 && y > 0, "
                   "Domain -> Complexes, Verify -> False], \"realvars\"]", "\"realvars\"", 0);
    /* Reals: every ordered field, hypotheses or not (Tier 0, no Mathlib). */
    check_str("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Reals, Verify -> False]"
              "[\"LeanTheorem\"]", BIND_ORD "\xE2\x88\x80 x : K, x ^ 2 = -4 \xE2\x86\x92 x = 7");
    check_true("With[{c = Certify[x == x, Domain -> Reals, Verify -> False]}, "
               "c[\"Domain\"] === Reals && c[\"Tier\"] === 0]");
    /* Rationals / Integers: Rat / Int, labelled. */
    check_str("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Rationals, Verify -> False]"
              "[\"LeanTheorem\"]", "\xE2\x88\x80 x : Rat, x ^ 2 = -4 \xE2\x86\x92 x = 7");
    check_str("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Integers, Verify -> False]"
              "[\"LeanTheorem\"]", "\xE2\x88\x80 x : Int, x ^ 2 = -4 \xE2\x86\x92 x = 7");
    /* Automatic, hypothesis-free: Int / Rat (equivalent), unchanged. */
    check_true("Certify[x + 1 == 1 + x, Verify -> False][\"Domain\"] === Integers");
    /* Anything else is refused, before any generator runs. */
    const char* bad[] = { "Algebraics", "\"Complexes\"", "Booleans", "Primes", "_",
                          "{Reals}", NULL };
    for (int i = 0; bad[i]; i++) {
        char buf[256];
        snprintf(buf, sizeof buf, "Check[Certify[x == 7, Assumptions -> x^2 == -4, "
                 "Domain -> %s, Verify -> False], \"domain\"]", bad[i]);
        assert_eval_eq(buf, "\"domain\"", 0);
    }
    assert_eval_eq("Check[Certify[Inverse[{{2, 1}, {1, 1}}], Domain -> Algebraics, "
                   "Verify -> False], \"domain\"]", "\"domain\"", 0);
    /* Element declarations are read the same way on every path. */
    check_true("Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Reals], "
               "Domain -> Complexes, Verify -> False][\"Domain\"] === Reals");
    check_true("Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Reals], "
               "Verify -> False][\"Domain\"] === Reals");
    check_str("Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Reals], "
              "Domain -> Rationals, Verify -> False][\"LeanTheorem\"]",
              "\xE2\x88\x80 x : Rat, x ^ 3 = 1 \xE2\x86\x92 x = 1");
    check_str("Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Integers], "
              "Domain -> Integers, Verify -> False][\"LeanTheorem\"]",
              "\xE2\x88\x80 x : Int, x ^ 3 = 1 \xE2\x86\x92 x = 1");
    assert_eval_eq("Check[Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Integers], "
                   "Verify -> False], \"element\"]", "\"element\"", 0);
    assert_eval_eq("Check[Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Integers], "
                   "Domain -> Reals, Verify -> False], \"element\"]", "\"element\"", 0);
    /* The generic form's K cannot capture a variable named K or Rat. */
    assert_eval_eq("Check[Certify[K == Rat, Assumptions -> K^2 == Rat^2, Verify -> False], "
                   "\"nocert\"]", "\"nocert\"", 0);
}

/* Finding 2: Integrate[f, x] == F is certified as D[F, x] == f, and says so. */
static void test_review6_integrate_claim(void) {
    test_mute_messages();
    assert_eval_eq("Check[Certify[Integrate[2 x, x] == x^2 + 5, Verify -> False], "
                   "\"antideriv\"]", "\"antideriv\"", 0);
    assert_eval_eq("leanI6 = Quiet[Certify[Integrate[2 x, x] == x^2 + 5, Verify -> False]]; Null",
                   "Null", 0);
    check_true("leanI6[\"Input\"] === HoldForm[Integrate[2 x, x] == x^2 + 5]");
    check_true("leanI6[\"Theorem\"] === HoldForm[ForAll[x, D[x^2 + 5, x] == 2 x]]");
    check_true("AnyTrue[leanI6[\"TrustedSteps\"], StringStartsQ[#, \"the claim Integrate[f, x] "
               "== F is certified as D[F, x] == f\"] && StringContainsQ[#, \"locally constant\"] &]");
    assert_eval_eq("leanI7 = Quiet[Certify[Integrate[-1/(1 + x^2), x] == ArcTan[1/x], "
                   "Verify -> False]]; Null", "Null", 0);
    check_true("leanI7[\"Theorem\"] === HoldForm[ForAll[x, D[ArcTan[1/x], x] == -1/(1 + x^2)]]");
    /* Certify[Integrate[f, x]] itself (Mathilda's own F) is unchanged. */
    check_true("!AnyTrue[Certify[Integrate[2 x, x], Verify -> False][\"TrustedSteps\"], "
               "StringContainsQ[#, \"locally constant\"] &]");
}

/* Finding 3: a symbol must be a plain variable (no value of any kind). */
static void test_review6_values(void) {
    test_mute_messages();
    const char* cases[] = {
        "leanP1 = Simplify[(x^2 - 1)/(x - 1)]; Protect[leanP1]; Certify[leanP1 == x + 1, "
            "Verify -> False]",
        "leanP2 /: leanP2^2 = -1; Certify[leanP2^2 == -1, Verify -> False]",
        "leanP3 /: leanP3 + 1 = 0; Certify[leanP3 + 1 == 0, Verify -> False]",
        "leanP4 = 5; Protect[leanP4]; Certify[Integrate[2 w, w] == w^2 + leanP4, Verify -> False]",
        "leanP5 = 3; Certify[Integrate[2 leanP5, leanP5] == leanP5^2, Verify -> False]",
        "leanP6[y_] := y; Certify[leanP6 == leanP6, Verify -> False]",
        "leanP7[a_][b_] := a; Certify[leanP7 + 1 == 1 + leanP7, Verify -> False]",
        "N[leanP8] = 2; Certify[leanP8 == leanP8, Verify -> False]",
        "SetAttributes[leanP9, Constant]; Certify[leanP9 + 1 == 1 + leanP9, Verify -> False]",
        "leanP10 = 2; Certify[Sum[leanP10, {leanP10, 1, n}] == n (n + 1)/2, Verify -> False]",
        NULL };
    for (int i = 0; cases[i]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "Check[%s, \"value\"]", cases[i]);
        assert_eval_eq(buf, "\"value\"", 0);
    }
    /* The message names the kind of value. */
    check_true("Clear[leanM6]; leanM6 /: leanM6^2 = -1; "
               "StringContainsQ[ToString[Lean`Private`symbolProblem[leanM6]], \"UpValues\"]");
    /* Protected alone is harmless; System constants are not values. */
    check_true("Protect[leanQ6]; Head[Certify[leanQ6 + 1 == 1 + leanQ6, Verify -> False]] === "
               "ProofObject");
    check_str("Certify[I^2 == -1, Verify -> False][\"LeanTheorem\"]", "(-1 : Int) = -1");
}

/* Finding 4: generator-private symbols cannot be named in a claim. */
static void test_review6_private(void) {
    test_mute_messages();
    assert_eval_eq("Check[Certify[Integrate[2, x] == Sqrt[2] Lean`Private`s1 x, Verify -> False], "
                   "\"context\"]", "\"context\"", 0);
    assert_eval_eq("Check[Certify[x == Lean`Internal`foo, Verify -> False], \"context\"]",
                   "\"context\"", 0);
    assert_eval_eq("Check[Certify[x == 1, Assumptions -> Lean`Private`s1 == 1, "
                   "Verify -> False], \"context\"]", "\"context\"", 0);
    /* The surd variables are fresh per call (printed s1, s2, ... in Lean). */
    check_true("With[{a = Keys[Certify[Integrate[1/(x^2 + 2), x], Verify -> False]"
               "[\"Certificate\"][\"Surds\"]], b = Keys[Certify[Integrate[1/(x^2 + 2), x], "
               "Verify -> False][\"Certificate\"][\"Surds\"]]}, a =!= b && "
               "Context @@ a === \"Lean`Private`\"]");
    check_true("StringStartsQ[Certify[Integrate[1/(x^2 + 2), x], Verify -> False]"
               "[\"LeanTheorem\"], \"\\[ForAll] x s1 : Rat, \"]");
}

/* Finding 5: hand-made ProofObjects are not authenticated. */
static void test_review6_proofobject(void) {
    check_true("ProofObject[<|\"Verified\" -> True, \"Input\" -> HoldForm[1 == 2]|>]"
               "[\"Verified\"] === Missing[\"Unauthenticated\"]");
    check_true("With[{p = Certify[x == x, Verify -> False]}, p[\"Verified\"] === "
               "Missing[\"NotChecked\"]]");
    /* Editing a genuine object's data unauthenticates it. */
    check_true("With[{a = Append[First[Certify[x == x, Verify -> False]], \"Verified\" -> True]}, "
               "ProofObject[a]][\"Verified\"] === Missing[\"Unauthenticated\"]");
    check_true("StringContainsQ[ToString[ProofObject[<|\"Verified\" -> True, \"Kind\" -> "
               "\"Identity\", \"Engine\" -> \"Lean\", \"Tier\" -> 0, \"Timing\" -> 0.|>]], "
               "\"Unauthenticated\"]");
}

/* Finding 6: the named trusted steps. */
static void test_review6_trusted(void) {
    check_true("MemberQ[Certify[Integrate[1/x + 1/(x + 1), x], Verify -> False]"
               "[\"TrustedSteps\"], s_String /; StringContainsQ[s, \"rational part R = r/w\"]]");
    check_true("MemberQ[Certify[Sum[k, {k, 1, n}], Verify -> False][\"TrustedSteps\"], "
               "s_String /; StringStartsQ[s, \"the closed form is proved for integer n >= 1 "
               "only\"]]");
}

/* Finding 7: Assumptions are walked held; SetOptions defaults are refused. */
static void test_review6_assumptions(void) {
    test_mute_messages();
    const char* evals[] = {
        "Reduce[x^2 == 1 && x > 0, x]", "Simplify[x^2 == 1 && x > 0]", "x == leanF6[2]",
        "D[x^2, x] == 2 x", "x == y + Simplify[(1 + z)^20 - Expand[(1 + z)^20] + 10^-25 z]",
        "Evaluate[x == 1]", NULL };
    assert_eval_eq("leanF6[y_] := y^2; Null", "Null", 0);
    for (int i = 0; evals[i]; i++) {
        char buf[512];
        snprintf(buf, sizeof buf, "Check[Certify[x == 1, Assumptions -> %s, Verify -> False], "
                 "\"eval\"]", evals[i]);
        assert_eval_eq(buf, "\"eval\"", 0);
    }
    /* Relations, connectives and Element pass, as written. */
    check_str("Certify[x == 1, Assumptions -> x == 1 || x == 2, Verify -> False]"
              "[\"LeanTheorem\"]",
              BIND_CHAR0 "\xE2\x88\x80 x : K, x = 1 \xE2\x88\xA8 x = 2 \xE2\x86\x92 x = 1");
    check_true("Certify[x == 1, Assumptions -> 0 < x < 2, Verify -> False][\"Domain\"] === Reals");
    /* A default set with SetOptions was evaluated by SetOptions: refused. */
    assert_eval_eq("SetOptions[Certify, Assumptions -> x > 0]; "
                   "leanD6 = Check[Certify[x == x, Verify -> False], \"optdefault\"]; "
                   "SetOptions[Certify, Assumptions -> {}]; leanD6", "\"optdefault\"", 0);
    assert_eval_eq("SetOptions[Certify, Assumptions :> Reduce[x^2 == 1 && x > 0, x]]; "
                   "leanD7 = Check[Certify[x == 1, Verify -> False], \"optdefault\"]; "
                   "SetOptions[Certify, Assumptions -> {}]; leanD7", "\"optdefault\"", 0);
    /* ... unless the call gives its own Assumptions. */
    check_true("SetOptions[Certify, Assumptions -> x > 0]; "
               "leanD8 = Head[Certify[x == x, Assumptions -> {}, Verify -> False]]; "
               "SetOptions[Certify, Assumptions -> {}]; leanD8 === ProofObject");
}

/* Refusal messages name what was written (v0.6523).  Sum[t, {k, n}] is read
 * as Sum[t, {k, 1, n}], as Sum defines it. */
static void test_refusals(void) {
    test_mute_messages();
    check_true("leanSn = Certify[Sum[k^2, {k, n}], Verify -> False]; "
               "leanSn[\"Kind\"] === \"ClosedFormSum\" && "
               "leanSn[\"Certificate\"][\"LowerLimit\"] === 1");
    check_true("Certify[Sum[k^2, {k, n}] == n (n + 1) (2 n + 1)/6, Verify -> False]"
               "[\"Certificate\"][\"PreCheck\"] === True");
    assert_eval_eq("Check[Certify[Integrate[1/Sqrt[x], x], Verify -> False], \"fracpow\"]",
                   "\"fracpow\"", 0);
    check_true("Certify[Solve[x^3 + x == 0, x, Reals], Verify -> False] === $Failed");
    check_true("Certify[x == x, Domain -> \"Complexes\", Verify -> False] === $Failed");
}

/* Option names (v0.6523): Certify's options are matched by name, so the
 * package creates no Engine / Tier / Domain / Verify / ... symbols of its own
 * (in Lean`, on $ContextPath, they captured the user's symbols, which then
 * survived ClearAll["Global`*"]).  A symbol or a string of that name works;
 * an unknown name, or a name symbol that has a value, is refused.  Runs last:
 * it clears Global`. */
static void test_option_names(void) {
    test_mute_messages();
    check_true("Context /@ {Engine, Tier, Domain, Verify, KeepFiles, LeanExecutable, "
               "LeanProject, Oracle, Counterexample, SearchTime} === ConstantArray[\"Global`\", 10]");
    check_true("Sort[Names[\"Lean`*\"]] === Sort[{\"$LeanExecutable\", \"$LeanProject\", \"$LeanVersion\", "
               "\"Certify\", \"LeanVerify\", \"ProofObject\", \"Refute\"}]");
    check_true("Certify[x + 1 == 1 + x, Engine -> \"Mathilda\"][\"Engine\"] === \"Mathilda\"");
    check_true("Certify[x + 1 == 1 + x, \"Engine\" -> \"Mathilda\"][\"Engine\"] === \"Mathilda\"");
    check_true("Certify[x + 1 == 1 + x, Verify -> False][\"Verified\"] === Missing[\"NotChecked\"]");
    check_true("Certify[x + 1 == 1 + x, \"Verify\" -> False, \"Domain\" -> Rationals]"
               "[\"Domain\"] === Rationals");
    assert_eval_eq("Check[Certify[x == x, Teir -> 1, Verify -> False], \"optx\"]", "\"optx\"", 0);
    assert_eval_eq("Check[LeanVerify[\"theorem t : True := trivial\", Keepfiles -> True], "
                   "\"optx\"]", "\"optx\"", 0);
    /* A name symbol with a value: refused while the option is still held. */
    assert_eval_eq("Engine = \"zzz\"; leanE = Check[Certify[x == x, Engine -> \"Lean\", "
                   "Verify -> False], \"optvalue\"]; Clear[Engine]; leanE", "\"optvalue\"", 0);
    check_true("Engine = \"zzz\"; leanE = Certify[x == x, \"Engine\" -> \"Mathilda\"]"
               "[\"Engine\"]; Clear[Engine]; leanE === \"Mathilda\"");
    /* The user's definitions of these names are ordinary Global` ones. */
    assert_eval_eq("Domain = 7; Verify[x_] := x + 1; Tier = 2; ClearAll[\"Global`*\"]; "
                   "{Domain, Verify[1], Tier}", "{Domain, Verify[1], Tier}", 0);
    check_true("Certify[x + 1 == 1 + x, Engine -> \"Mathilda\"][\"Engine\"] === \"Mathilda\"");
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
    TEST(test_review4_reproducers);
    TEST(test_review4_property);
    TEST(test_review5_scoped);
    TEST(test_review4_claims);
    TEST(test_review6_domain);
    TEST(test_review6_integrate_claim);
    TEST(test_review6_values);
    TEST(test_review6_private);
    TEST(test_review6_proofobject);
    TEST(test_review6_trusted);
    TEST(test_review6_assumptions);
    TEST(test_refusals);
    TEST(test_option_names);
    printf("All Lean certificate generation tests passed.\n");
    return 0;
}
