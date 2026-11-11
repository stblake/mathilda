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
    check_true("leanT[\"TrustedSteps\"] === {\"Together normal form of each side "
               "(num/den used to clear denominators)\"}");
    check_true("!StringContainsQ[leanT[\"LeanTheorem\"], \"/\"]");
    check_true("Certify[Apart[(x^2+1)/(x^2-1)], Verify -> False][\"Certificate\"][\"Form\"] === "
               "\"ClearedDenominators\"");
}

static void test_options(void) {
    check_str("Certify[x^2 == x, Assumptions -> {x == 1}, Verify -> False][\"LeanTheorem\"]",
              "\xE2\x88\x80 x : Int, x = 1 \xE2\x86\x92 x ^ 2 = x");
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
    assert_eval_eq("Certify[Integrate[x, x], Verify -> False]", "$Failed", 0);
}

static void test_heuristic_engine(void) {
    /* The PossibleZeroQ false positive: the heuristic engine says True. The
     * claim is false (the difference is x/10^25); test_lean_verify.c checks
     * that Lean says False. */
    assert_eval_eq("leanH = Certify[(1+x)^20 == Expand[(1+x)^20] - 10^-25 x, "
                   "Engine -> \"Mathilda\"]; Null", "Null", 0);
    check_true("leanH[\"Verified\"] === True && leanH[\"Tier\"] === \"Heuristic\" && "
               "leanH[\"Engine\"] === \"Mathilda\"");
    check_true("leanH[\"Domain\"] === Rationals");
    check_true("Certify[(x+1)^2 == x^2 + 1, Engine -> \"Mathilda\"][\"Verified\"] === False");
}

static void test_attributes(void) {
    assert_eval_eq("Attributes[Certify]", "{HoldFirst, Protected}", 0);
    assert_eval_eq("Attributes[ProofObject]", "{HoldFirst, Protected}", 0);
    assert_eval_eq("Attributes[LeanVerify]", "{Protected}", 0);
    assert_eval_eq("Context[Certify]", "\"Lean`\"", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_identity_source);
    TEST(test_computations);
    TEST(test_rational_identity);
    TEST(test_options);
    TEST(test_errors);
    TEST(test_heuristic_engine);
    TEST(test_attributes);
    printf("All Lean certificate generation tests passed.\n");
    return 0;
}
