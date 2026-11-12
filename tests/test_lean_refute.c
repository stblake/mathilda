/*
 * test_lean_refute.c -- counterexample search, refutation statements and
 * Refute (src/internal/lean/Refute.m).  Needs no Lean installation: every
 * call runs with "Verify" -> False (the refutation's Lean source is built
 * and passed through the source check, not run) or with Engine -> "Exact".
 *
 * Faithfulness is the point: a witness must lie in the statement's domain
 * (FindInstance's domain argument is not trusted), satisfy its hypotheses,
 * falsify its conclusion and keep every denominator of the claim nonzero.
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

static void check_true(const char* input) {
    assert_eval_eq(input, "True", 0);
}

/* Refute[<call>, "Verify" -> False]["Counterexample"] === witness. */
static void check_witness(const char* call, const char* witness) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Refute[%s, \"Verify\" -> False][\"Counterexample\"] === %s",
             call, witness);
    check_true(buf);
}

static void test_search(void) {
    /* Int (auto): the grid's first point by height, -1 before 1. */
    check_witness("(x+1)^3 == x^3+3x^2+2x+1", "{x -> -1}");
    /* Ordered field: rationals, grid order. */
    check_witness("ForAll[{x, y}, x^2 + y^2 >= 2 x y + 1/10]", "{x -> 0, y -> 0}");
    /* Domain -> Integers: an integer witness. */
    check_witness("x^3 + y^3 == (x + y)^3, \"Domain\" -> Integers", "{x -> -1, y -> -1}");
    check_true("AllTrue[Refute[x^3 + y^3 == (x + y)^3, \"Domain\" -> Integers, "
               "\"Verify\" -> False][\"Counterexample\"][[All, 2]], IntegerQ]");
    /* Denominators: the cleared statement is false at the pole x = 0, which is
     * no counterexample to the rational claim. */
    check_witness("1/x == 1", "{x -> -1}");
    check_true("Refute[(x^2 - 1)/(x - 1) == x + 2, \"Verify\" -> False][\"Counterexample\"] "
               "=!= {x -> 1}");
    /* Sums: the first n >= a where the closed form is wrong. */
    check_witness("Sum[k^3, {k, 1, n}] == n^2 (n+1)^2/4 + 1", "{n -> 1}");
    /* Hypotheses: an antecedent that holds at the witness. */
    check_witness("x^2 == x, Assumptions -> x^2 == 2 x", "{x -> 2}");
    check_witness("Implies[x^2 + y^2 == 2 && x == y, x == 1]", "{x -> -1, y -> -1}");
    /* Determinism: the same witness on every call. */
    check_true("Table[Refute[ForAll[{x, y, z}, x^2 + y^2 + z^2 >= x y + y z + 1], "
               "\"Verify\" -> False][\"Counterexample\"], {3}] === "
               "ConstantArray[{x -> 0, y -> 0, z -> 0}, 3]");
}

static void test_exact_only(void) {
    test_mute_messages();
    /* No rational witness: FindInstance over the rationals returns -2 I, which
     * is refused for the kernel; the char-0 statement is refuted over C, by
     * exact arithmetic only. */
    assert_eval_eq("leanR1 = Refute[Implies[x^2 + 4 == 0, x == 7]]; Null", "Null", 0);
    check_true("MemberQ[{{x -> -2 I}, {x -> 2 I}}, leanR1[\"Counterexample\"]]");
    check_true("leanR1[\"Verified\"] === Missing[\"Exact\", True]");
    check_true("leanR1[\"Engine\"] === \"Exact\" && leanR1[\"Kind\"] === \"Refutation\"");
    check_true("MemberQ[leanR1[\"TrustedSteps\"], s_String /; StringContainsQ[s, \"algebraic "
               "witness\"]]");
    /* The domain check itself: -2 I is not a rational, nor a real. */
    check_true("!Lean`Private`valueOK[\"Rat\", -2 I] && !Lean`Private`valueOK[\"OrderedField\", "
               "-2 I] && Lean`Private`valueOK[\"CharZeroField\", -2 I] && "
               "Lean`Private`valueOK[\"OrderedField\", -Sqrt[2]]");
    /* A claimed antiderivative: D[F, x] - f is a nonzero rational at x = 0. */
    assert_eval_eq("leanR2 = Refute[Integrate[2 x, x] == x^2 + x]; Null", "Null", 0);
    check_true("leanR2[\"Counterexample\"] === {x -> 0} && "
               "leanR2[\"Verified\"] === Missing[\"Exact\", True]");
    /* Engine -> "Exact" on a kernel-refutable claim. */
    assert_eval_eq("leanR3 = Refute[(x+1)^3 == x^3+3x^2+2x+1, \"Engine\" -> \"Exact\"]; Null",
                   "Null", 0);
    check_true("leanR3[\"Verified\"] === Missing[\"Exact\", True]");
    check_true("leanR3[\"Properties\"] === {\"Cached\", \"Certificate\", \"Checks\", "
               "\"Counterexample\", \"Domain\", \"Engine\", \"Hypotheses\", \"Input\", \"Kind\", "
               "\"RefutedStatement\", \"Result\", \"Theorem\", \"Tier\", \"Timing\", "
               "\"TrustedSteps\", \"Verified\"}");
    check_true("MatchQ[leanR3[\"Checks\"], {_String -> True}]");
    check_true("StringContainsQ[ToString[leanR3], \"\\[LeftGuillemet]Refutation\\[RightGuillemet] "
               "Verified: Missing[Exact, True] | Exact |\"]");
    check_true("Refute[(x+1)^3 == x^3+3x^2+2x+1, \"Engine\" -> \"Exact\", \"Verify\" -> False]"
               "[\"Verified\"] === Missing[\"NotChecked\"]");
}

static void test_lean_text(void) {
    /* Int identity. */
    assert_eval_eq("leanT1 = Refute[(x+1)^3 == x^3+3x^2+2x+1, \"Verify\" -> False]; Null",
                   "Null", 0);
    check_true("leanT1[\"Properties\"] === {\"Axioms\", \"Cached\", \"Certificate\", "
               "\"Counterexample\", \"Domain\", \"Engine\", \"Hypotheses\", \"Input\", \"Kind\", "
               "\"LeanSource\", \"LeanTheorem\", \"LeanVersion\", \"Logs\", \"ProofScript\", "
               "\"RefutedStatement\", \"Result\", \"Theorem\", \"Tier\", \"Timing\", "
               "\"TrustedSteps\", \"Verified\"}");
    check_true("leanT1[\"Kind\"] === \"Refutation\" && leanT1[\"Verified\"] === "
               "Missing[\"NotChecked\"]");
    check_true("leanT1[\"LeanTheorem\"] === "
               "\"\\[Not] (\\[ForAll] x : Int, (1 + x) ^ 3 = 1 + 2 * x + 3 * x ^ 2 + x ^ 3)\"");
    check_true("leanT1[\"ProofScript\"] === \"fun h => absurd (h (-1 : Int)) (by decide +kernel)\"");
    check_true("leanT1[\"Theorem\"] === HoldForm[!ForAll[{x}, (1 + x)^3 == 1 + 2 x + 3 x^2 + x^3]]");
    check_true("leanT1[\"RefutedStatement\"] === HoldForm[ForAll[{x}, (1 + x)^3 == 1 + 2 x + "
               "3 x^2 + x^3]]");
    check_true("leanT1[\"Input\"] === HoldForm[(x+1)^3 == x^3+3x^2+2x+1]");
    /* Generic ordered field: (K := Rat). */
    assert_eval_eq("leanT2 = Refute[ForAll[{x, y}, x^2 + y^2 >= 2 x y + 1/10], "
                   "\"Verify\" -> False]; Null", "Null", 0);
    check_true("leanT2[\"ProofScript\"] === \"fun h => absurd (h (K := Rat) (0 : Rat) (0 : Rat)) "
               "(by decide +kernel)\"");
    check_true("StringStartsQ[leanT2[\"LeanTheorem\"], \"\\[Not] (\\[ForAll] {K : Type}\"]");
    /* Generic char 0 with an antecedent. */
    assert_eval_eq("leanT3 = Refute[x^2 == x, Assumptions -> x^2 == 2 x, \"Verify\" -> False]; "
                   "Null", "Null", 0);
    check_true("leanT3[\"ProofScript\"] === \"fun h => absurd (h (K := Rat) (2 : Rat) "
               "(by decide +kernel)) (by decide +kernel)\"");
    /* A sum: the preamble's definitions, the witness as m = n - a. */
    assert_eval_eq("leanT4 = Refute[Sum[k^3, {k, 1, n}] == n^2 (n+1)^2/4 + 1, "
                   "\"Verify\" -> False]; Null", "Null", 0);
    check_true("leanT4[\"ProofScript\"] === \"fun h => absurd (h 0) (by decide +kernel)\"");
    check_true("StringContainsQ[leanT4[\"LeanSource\"], \"def mathilda_cert_1_S\"]");
    check_true("leanT4[\"Theorem\"] === HoldForm[!ForAll[n, n >= 1, Sum[k^3, {k, 1, n}] == "
               "1 + 1/4 n^2 (1 + n)^2]]");
    /* Every source passes the source check. */
    check_true("AllTrue[{leanT1, leanT2, leanT3, leanT4}, "
               "TrueQ[Lean`Private`checkSource[#[\"LeanSource\"], False, True][\"OK\"]] &]");
}

static void test_refusals(void) {
    test_mute_messages();
    assert_eval_eq("Check[Refute[PrimeQ[7]], \"nostmt\", Refute::nostmt]", "\"nostmt\"", 0);
    assert_eval_eq("Check[Refute[Integrate[x, x]], \"nostmt\", Refute::nostmt]", "\"nostmt\"", 0);
    assert_eval_eq("Check[Refute[Factor[x^2 - 1]], \"nostmt\", Refute::nostmt]", "\"nostmt\"", 0);
    /* A true claim: no counterexample (Motzkin's polynomial is nonnegative). */
    assert_eval_eq("Check[Refute[ForAll[{x, y}, x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1 >= 0], "
                   "\"SearchTime\" -> 2], \"none\", Refute::none]", "\"none\"", 0);
    check_true("Quiet[Refute[ForAll[{x, y}, x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1 >= 0], "
               "\"SearchTime\" -> 2]] === $Failed");
    check_true("Quiet[Refute[(x+1)^2 == x^2 + 2 x + 1, \"Verify\" -> False]] === $Failed");
    /* Options. */
    assert_eval_eq("Check[Refute[x == x + 1, Foo -> 1], \"optx\", Refute::optx]", "\"optx\"", 0);
    assert_eval_eq("Check[Refute[x == x + 1, \"Engine\" -> \"Z3\"], \"engine\", Refute::engine]",
                   "\"engine\"", 0);
    assert_eval_eq("Check[Refute[x == x + 1, \"Domain\" -> Algebraics], \"domain\", "
                   "Certify::domain]", "\"domain\"", 0);
    /* The claim and the assumptions are walked held, as in Certify. */
    assert_eval_eq("Check[Refute[Simplify[x] == x + 1], \"eval\", Certify::eval]", "\"eval\"", 0);
    assert_eval_eq("leanRv = 3; leanRr = Check[Refute[leanRv == 2], \"value\", Certify::value]; "
                   "Clear[leanRv]; leanRr", "\"value\"", 0);
    /* Failure injection: a non-witness from the search hook is refused by the
     * exact check, so no refutation is built. */
    assert_eval_eq("Lean`Private`$witnessOverride = {x -> 0}; "
                   "leanRw = Quiet[Refute[(x+1)^3 == x^3+3x^2+2x+1, \"Verify\" -> False]]; "
                   "Lean`Private`$witnessOverride = None; leanRw", "$Failed", 0);
    assert_eval_eq("Lean`Private`$witnessOverride = {x -> 1/2}; "
                   "leanRw = Quiet[Refute[(x+1)^3 == x^3+3x^2+2x+1, \"Verify\" -> False]]; "
                   "Lean`Private`$witnessOverride = None; leanRw", "$Failed", 0);
}

/* A Domain / Assumptions that the refuted statement omits (only identity
 * claims state them) must restrict the witness: a refutation is never
 * reported for a claim other than the one stated.  x^2 >= x is true over
 * the integers, so there is no counterexample; x -> 1/2 refutes it for
 * x > 0. */
static void test_domain_assumptions(void) {
    test_mute_messages();
    assert_eval_eq("Check[Refute[x^2 >= x, Assumptions -> Element[x, Integers], "
                   "\"Verify\" -> False, \"SearchTime\" -> 2], \"none\", Refute::none]",
                   "\"none\"", 0);
    check_true("Quiet[Refute[x^2 >= x, \"Domain\" -> Integers, \"Verify\" -> False, "
               "\"SearchTime\" -> 2]] === $Failed");
    check_witness("x^2 >= x", "{x -> 1/2}");
    check_witness("x^2 >= x, Assumptions -> x > 0", "{x -> 1/2}");
    check_true("MemberQ[Refute[x^2 >= x, Assumptions -> x > 0, \"Verify\" -> False]"
               "[\"TrustedSteps\"], s_String /; StringContainsQ[s, \"the witness satisfies x > 0\"]]");
    /* No Certify::optignored any more: the options are honoured. */
    check_true("Block[{$MessageList = {}}, Quiet[Refute[x^2 >= x, Assumptions -> x > 0, "
               "\"Verify\" -> False]]; $MessageList] === {}");
    /* Implication of equations: x -> -1 violates x > 0; x == 1 holds otherwise. */
    check_true("Quiet[Refute[Implies[x^2 == 1, x == 1], Assumptions -> x > 0, "
               "\"Verify\" -> False, \"SearchTime\" -> 2]] === $Failed");
    /* Sums: the first n >= a that satisfies the assumption. */
    check_witness("Sum[k, {k, 1, n}] == n^2", "{n -> 2}");
    check_witness("Sum[k, {k, 1, n}] == n^2, Assumptions -> n > 3", "{n -> 4}");
    /* The helper: identity claims state the options themselves (no "Assume"). */
    check_true("!KeyExistsQ[Lean`Private`withAssume[<|\"Form\" -> \"Forall\", \"Vars\" -> {x}|>, "
               "<|\"Kind\" -> \"Identity\"|>, <|\"Domain\" -> Integers, \"Assumptions\" -> x > 0|>], "
               "\"Assume\"]");
    check_true("Lean`Private`withAssume[<|\"Form\" -> \"Forall\", \"Vars\" -> {x, y}|>, "
               "<|\"Kind\" -> \"Nonnegativity\"|>, <|\"Domain\" -> Integers, "
               "\"Assumptions\" -> {x > 0}|>][\"Assume\"] === "
               "Hold[And[Element[x, Integers], Element[y, Integers], x > 0]]");
}

/* Certify searches only after a kernel REJECTION: no Lean (or Verify ->
 * False) means no search and the usual 19 keys. */
static void test_certify_no_lean(void) {
    test_mute_messages();
    check_true("Length[Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1, \"LeanExecutable\" -> "
               "\"/nonexistent/lean\"]][\"Properties\"]] === 19");
    check_true("!KeyExistsQ[First[Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1, "
               "\"LeanExecutable\" -> \"/nonexistent/lean\"]]], \"Counterexample\"]");
    check_true("Length[Certify[(x+1)^3 == x^3+3x^2+2x+1, Verify -> False][\"Properties\"]] === 19");
    check_true("!KeyExistsQ[First[Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1, Engine -> \"Exact\"]]], "
               "\"Counterexample\"]");
    check_true("Options[Certify][[1]] === (\"Engine\" -> \"Lean\") && "
               "Options[Certify][[-1]] === (\"Counterexample\" -> Automatic)");
    check_true("Lean`Private`searchBudget[Automatic] === 2 && "
               "Lean`Private`searchBudget[0.5] === 0.5");
}

/* Adversarial review (v0.6528). */
static void test_review(void) {
    test_mute_messages();
    /* Finding 1/2: a Sum with an iterator the sum generator does not take
     * is not stated by evaluating it (that refuted Mathilda's evaluation of
     * the claim, not the claim: n -> -3 where the claim holds, n -> -1/2, a
     * non-integer upper limit).  No certificate, no refutation. */
    assert_eval_eq("Check[Refute[Sum[k, {k, 1, n + 1}] == (n + 1) (n + 2)/2 + "
                   "n (n + 1) (n - 1) (n + 2) (n - 2)/120, \"Verify\" -> False], \"nocert\", "
                   "Certify::nocert]", "\"nocert\"", 0);
    assert_eval_eq("Check[Refute[Sum[k, {k, 1, 2 n}] == n (2 n + 1) + n^3, "
                   "\"Verify\" -> False], \"nocert\", Certify::nocert]", "\"nocert\"", 0);
    check_true("Quiet[Refute[Sum[k, {k, n, 2 n}] == 0, \"Engine\" -> \"Exact\"]] === $Failed");
    /* Certify: the same claims have no certificate (any engine). */
    assert_eval_eq("Check[Certify[Sum[k, {k, 1, n + 1}] == (n + 1) (n + 2)/2, "
                   "\"Verify\" -> False], \"nocert\", Certify::nocert]", "\"nocert\"", 0);
    check_true("AllTrue[{Hold[Sum[k, {k, 1, n + 1}] == (n + 1) (n + 2)/2], "
               "Hold[Sum[k, {k, 1, n - 1}] == n (n - 1)/2], Hold[Sum[k, {k, 1, 2 n}] == n (2 n + 1)], "
               "Hold[Sum[k, {k, n, 2 n}] == 3 n (n + 1)/2], Hold[Sum[k, {k, 1, n - 2}] == "
               "(n - 2) (n - 1)/2]}, Function[h, Quiet[ReleaseHold[h /. Equal[a_, b_] :> "
               "Certify[Equal[a, b], \"Engine\" -> \"Exact\"]]] === $Failed]]");
    /* The supported shape is unchanged. */
    check_witness("Sum[k, {k, 1, n}] == n^2", "{n -> 2}");
    /* Finding 4: a claim variable named K is not captured by the generic
     * type variable (K1 here). */
    assert_eval_eq("leanT5 = Refute[ForAll[K, K^2 >= K + 1/10], \"Verify\" -> False]; Null",
                   "Null", 0);
    check_true("leanT5[\"ProofScript\"] === \"fun h => absurd (h (K1 := Rat) (0 : Rat)) "
               "(by decide +kernel)\"");
    check_true("StringStartsQ[leanT5[\"LeanTheorem\"], \"\\[Not] (\\[ForAll] {K1 : Type} "
               "[Lean.Grind.Field K1]\"] && StringContainsQ[leanT5[\"LeanTheorem\"], "
               "\"\\[ForAll] K : K1, K ^ 2\"]");
    check_true("TrueQ[Lean`Private`checkSource[leanT5[\"LeanSource\"], False, True][\"OK\"]]");
    check_true("Lean`Private`genericVar[\"x + K\", \"K1\"] === \"K2\" && "
               "Lean`Private`genericVar[\"Kx + xK + K.f\"] === \"K\"");
    check_true("Lean`Private`splitGeneric[\"\\[ForAll] \" <> "
               "Lean`Private`genericBinders[\"CharZero\", \"K3\"] <> \", body\"] === "
               "{\" \" <> Lean`Private`genericBinders[\"CharZero\", \"K3\"], \"body\"}");
    /* Finding 5: only a top-level And of the Assumptions is split. */
    check_witness("ForAll[{x, y}, x^2 >= y], Assumptions -> Not[y > 0 && y > 1/2]",
                  "{x -> 0, y -> 1/2}");
    check_true("Lean`Private`topConjuncts[{a && b, Not[c && d] && (e || (f && g))}] === "
               "{a, b, Not[c && d], e || (f && g)}");
}

/* The pieces of the search on their own. */
static void test_helpers(void) {
    check_true("Lean`Private`refZero[(-2 I)^2 + 4] === True");
    check_true("Lean`Private`refZero[-2 I - 7] === False");
    check_true("Lean`Private`refZero[(1 + Sqrt[2])^2 - 3 - 2 Sqrt[2]] === True");
    check_true("Lean`Private`refZero[-1 - Sqrt[2]] === False");
    check_true("Lean`Private`refSign[Sqrt[2] - 1] === Lean`Private`refUndecided");
    check_true("Lean`Private`implDepth[Hold[Implies[a, Implies[b, c]]]] === 2");
    check_true("Lean`Private`algebraicQ[1/2 + Sqrt[3] I] && !Lean`Private`algebraicQ[Pi] && "
               "!Lean`Private`algebraicQ[1.5]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_search);
    TEST(test_exact_only);
    TEST(test_lean_text);
    TEST(test_refusals);
    TEST(test_domain_assumptions);
    TEST(test_certify_no_lean);
    TEST(test_helpers);
    TEST(test_review);
    printf("All Lean refutation tests passed.\n");
    return 0;
}
