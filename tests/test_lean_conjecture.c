/*
 * test_lean_conjecture.c -- Conjecture (src/internal/lean/Conjecture.m) and
 * Fibonacci sum certificates.  Needs no Lean installation: every call runs
 * with Engine -> "Exact" (Mathilda's exact re-check, status "Checked") or
 * with "Verify" -> False (the Lean source is built and passed through the
 * source check, not run).
 *
 * The point: a guess is only a candidate until it is certified, the claim is
 * certified as written (the user's iterator name survives), and a guess
 * with no certificate kind is reported "Tested", never more.
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

static void test_polynomial(void) {
    assert_eval_eq("leanJ1 = Conjecture[Sum[k^4, {k, 1, n}], n, \"Engine\" -> \"Exact\"]; Null",
                   "Null", 0);
    check_true("leanJ1[\"Form\"] === \"Polynomial\"");
    check_true("Expand[leanJ1[\"Guess\"] - n (n + 1) (2 n + 1) (3 n^2 + 3 n - 1)/30] === 0");
    check_true("leanJ1[\"Status\"] === \"Checked\" && "
               "leanJ1[\"Verified\"] === Missing[\"Exact\", True]");
    check_true("leanJ1[\"FittedOn\"] === {1, 12} && leanJ1[\"TestedOn\"] === {1, 52}");
    check_true("leanJ1[\"Properties\"] === {\"Engine\", \"FittedOn\", \"Form\", \"Guess\", "
               "\"Input\", \"Kind\", \"Proof\", \"Status\", \"TestedOn\", \"Theorem\", \"Timing\", "
               "\"Variable\", \"Verified\"}");
    check_true("leanJ1[\"Kind\"] === \"Conjecture\" && leanJ1[\"Variable\"] === n && "
               "leanJ1[\"Engine\"] === \"Exact\"");
    /* The claim is certified as written: the iterator k survives. */
    check_true("leanJ1[\"Input\"] === HoldForm[Sum[k^4, {k, 1, n}]]");
    check_true("MatchQ[leanJ1[\"Proof\"][\"Input\"], HoldForm[Sum[k^4, {k, 1, n}] == _]]");
    check_true("MatchQ[leanJ1[\"Theorem\"], HoldForm[ForAll[n, n >= 1, Sum[k^4, {k, 1, n}] == _]]]");
    check_true("StringStartsQ[ToString[leanJ1], \"ProofObject[\\[LeftGuillemet]Conjecture"
               "\\[RightGuillemet] Checked: Sum[k^4, {k, 1, n}] == \"] && "
               "StringContainsQ[ToString[leanJ1], \"| Exact |\"]");
    /* {k, n} starts at 1. */
    check_true("Conjecture[Sum[k^2, {k, n}], n, \"Engine\" -> \"Exact\"][\"FittedOn\"] === {1, 12}");
}

static void test_forms(void) {
    assert_eval_eq("leanJ2 = Conjecture[Sum[k^2 2^k, {k, 0, n}], n, \"Engine\" -> \"Exact\"]; Null",
                   "Null", 0);
    check_true("leanJ2[\"Form\"] === \"Geometric\" && leanJ2[\"Status\"] === \"Checked\"");
    check_true("Expand[leanJ2[\"Guess\"] - (2^(n + 1) (n^2 - 2 n + 3) - 6)] === 0");
    assert_eval_eq("leanJ3 = Conjecture[Sum[Fibonacci[k]^2, {k, 1, n}], n, \"Engine\" -> \"Exact\"]; "
                   "Null", "Null", 0);
    check_true("leanJ3[\"Form\"] === \"Fibonacci\" && leanJ3[\"Status\"] === \"Checked\"");
    check_true("leanJ3[\"Guess\"] === Fibonacci[n] Fibonacci[n + 1]");
    /* No certificate kind for Binomial sums: tested, not proved; no message. */
    assert_eval_eq("leanJ4 = Block[{$MessageList = {}}, {Conjecture[Sum[Binomial[n, k]^2, "
                   "{k, 0, n}], n, \"Engine\" -> \"Exact\"], $MessageList}]; Null", "Null", 0);
    check_true("leanJ4[[2]] === {}");
    check_true("leanJ4[[1]][\"Form\"] === \"Hypergeometric\" && "
               "leanJ4[[1]][\"Guess\"] === Binomial[2 n, n]");
    check_true("leanJ4[[1]][\"Status\"] === \"Tested\" && "
               "leanJ4[[1]][\"Verified\"] === Missing[\"NoCertificate\", \"Certify::eval\"] && "
               "leanJ4[[1]][\"Proof\"] === Missing[\"NoCertificate\"]");
    check_true("leanJ4[[1]][\"Theorem\"] === HoldForm[ForAll[n, n >= 0, "
               "Sum[Binomial[n, k]^2, {k, 0, n}] == Binomial[2 n, n]]]");
    check_true("StringContainsQ[ToString[leanJ4[[1]]], \"Tested: \"] && "
               "StringContainsQ[ToString[leanJ4[[1]]], \"| no certificate |\"]");
    /* Rational and named hypergeometric forms. */
    check_true("Conjecture[Sum[1/(k (k + 1)), {k, 1, n}], n, \"Engine\" -> \"Exact\"][\"Guess\"] "
               "=== n/(n + 1)");
    check_true("Conjecture[Product[(2 k - 1)/(2 k), {k, 1, n}], n, \"Engine\" -> \"Exact\"]"
               "[\"Guess\"] === Binomial[2 n, n]/4^n");
    /* "Forms" restricts the guessers. */
    check_true("Conjecture[Sum[k^3, {k, 1, n}], n, \"Engine\" -> \"Exact\", \"Forms\" -> "
               "{\"Rational\"}][\"Form\"] === \"Rational\"");
}

static void test_verify_false(void) {
    assert_eval_eq("leanJ5 = {Conjecture[Sum[k^4, {k, 1, n}], n, \"Verify\" -> False], "
                   "Conjecture[Sum[Fibonacci[k]^2, {k, 1, n}], n, \"Verify\" -> False]}; Null",
                   "Null", 0);
    check_true("AllTrue[leanJ5, #[\"Status\"] === \"Tested\" && "
               "Head[#[\"Proof\"]] === ProofObject && "
               "#[\"Proof\"][\"Verified\"] === Missing[\"NotChecked\"] &]");
    check_true("AllTrue[leanJ5, #[\"Engine\"] === {\"Lean\", \"Exact\"} &]");
    check_true("StringContainsQ[leanJ5[[2]][\"Proof\"][\"LeanTheorem\"], \"mathilda_cert_1_F\"]");
    check_true("StringContainsQ[leanJ5[[2]][\"Proof\"][\"LeanSource\"], "
               "\"def mathilda_cert_1_F : Nat \\[RightArrow] Rat\\n  | 0 => 0\\n  | 1 => 1\\n"
               "  | n + 2 => mathilda_cert_1_F n + mathilda_cert_1_F (n + 1)\\n\"]");
    check_true("AllTrue[leanJ5, TrueQ[Lean`Private`checkSource[#[\"Proof\"][\"LeanSource\"], "
               "False, True][\"OK\"]] &]");
    check_true("MemberQ[leanJ5[[2]][\"Proof\"][\"TrustedSteps\"], s_String /; "
               "StringContainsQ[s, \"Fibonacci on the naturals\"]]");
    /* Without Fibonacci nothing of it appears. */
    check_true("!StringContainsQ[leanJ5[[1]][\"Proof\"][\"LeanSource\"], \"_F\"]");
}

static void test_refusals(void) {
    test_mute_messages();
    assert_eval_eq("Check[Conjecture[Sum[1.5^k, {k, 0, n}], n], \"nonexact\", "
                   "Conjecture::nonexact]", "\"nonexact\"", 0);
    assert_eval_eq("Check[Conjecture[Sum[Prime[k], {k, 1, n}], n], \"noguess\", "
                   "Conjecture::noguess]", "\"noguess\"", 0);
    assert_eval_eq("n = 3; leanJr = Check[Conjecture[Sum[k, {k, 1, n}], n], \"var\", "
                   "Conjecture::var]; Clear[n]; leanJr", "\"var\"", 0);
    assert_eval_eq("Check[Conjecture[Sum[k, {k, 1, n}], 2 n], \"var\", Conjecture::var]",
                   "\"var\"", 0);
    assert_eval_eq("Check[Conjecture[Sum[k, {k, 1, n}], n, Foo -> 1], \"optx\", Conjecture::optx]",
                   "\"optx\"", 0);
    assert_eval_eq("Check[Conjecture[Sum[k, {k, 1, n}], n, \"FitPoints\" -> 2], \"badopt\", "
                   "Conjecture::badopt]", "\"badopt\"", 0);
    assert_eval_eq("Check[Conjecture[Sum[k, {k, 1, n}], n, \"Forms\" -> {\"Magic\"}], \"badopt\", "
                   "Conjecture::badopt]", "\"badopt\"", 0);
    check_true("MemberQ[Attributes[Conjecture], HoldAll] && "
               "MemberQ[Attributes[Conjecture], Protected]");
    check_true("Options[Conjecture][[1]] === (\"Engine\" -> {\"Lean\", \"Exact\"})");
    /* An edited object is not authenticated. */
    check_true("With[{a = Append[First[Conjecture[Sum[k, {k, 1, n}], n, \"Engine\" -> \"Exact\"]], "
               "\"Status\" -> \"Proved\"]}, StringContainsQ[ToString[ProofObject[a]], "
               "\"Unauthenticated\"]]");
}

/* Fibonacci sum certificates on their own (v0.6527). */
static void test_fibonacci_sums(void) {
    test_mute_messages();
    check_true("Certify[Sum[Fibonacci[k], {k, 0, n}] == Fibonacci[n + 2] - 1, "
               "Engine -> \"Exact\"][\"Verified\"] === Missing[\"Exact\", True]");
    check_true("Quiet[Certify[Sum[Fibonacci[k], {k, 0, n}] == Fibonacci[n + 2], "
               "Engine -> \"Exact\"]][\"Verified\"] === Missing[\"Exact\", False]");
    /* A negative shift in the term, and Fibonacci[2 k], have no certificate. */
    assert_eval_eq("Check[Certify[Sum[Fibonacci[k - 1], {k, 1, n}] == Fibonacci[n + 1] - 1, "
                   "Verify -> False], \"nocert\", Certify::nocert]", "\"nocert\"", 0);
    assert_eval_eq("Check[Certify[Sum[Fibonacci[2 k], {k, 1, n}] == Fibonacci[2 n + 1] - 1, "
                   "Verify -> False], \"nocert\", Certify::nocert]", "\"nocert\"", 0);
    /* Fibonacci is read only inside a sum claim. */
    assert_eval_eq("Check[Certify[Fibonacci[n + 2] == Fibonacci[n] + Fibonacci[n + 1], "
                   "Verify -> False], \"eval\", Certify::eval]", "\"eval\"", 0);
    check_true("Lean`Private`exFibLin[Fibonacci[n + 3] - 2 Fibonacci[n + 1] - Fibonacci[n], n] "
               "=== 0");
}

/* Adversarial review (v0.6528): a Sum whose iterator the closed-form sum
 * generator does not take ({k, 1, n - 2}, {k, 1, n + 1}, {k, 1, 2 n},
 * {k, n, 2 n}) used to fall to the identity generator, which evaluated both
 * sides: Mathilda's own Sum output was certified against itself (Lean:
 * "Proved", the theorem a tautology, no trusted step; false at n = 0 for
 * {k, 1, n - 2}).  Such a claim has no certificate: "Tested". */
static void test_review_unsupported_sums(void) {
    test_mute_messages();
    assert_eval_eq("leanJu = {Conjecture[Sum[k, {k, 1, n - 2}], n, \"Start\" -> 2, "
                   "\"Engine\" -> \"Exact\"], Conjecture[Sum[k, {k, 1, n + 1}], n, "
                   "\"Engine\" -> \"Exact\"], Conjecture[Sum[k, {k, 1, 2 n}], n, "
                   "\"Engine\" -> \"Exact\"], Conjecture[Sum[k, {k, n, 2 n}], n, "
                   "\"Engine\" -> \"Exact\"], Conjecture[Sum[k, {k, 1, n - 1}], n, "
                   "\"Start\" -> 1, \"Engine\" -> \"Exact\"]}; Null", "Null", 0);
    check_true("AllTrue[leanJu, #[\"Status\"] === \"Tested\" && "
               "#[\"Verified\"] === Missing[\"NoCertificate\", \"Certify::nocert\"] && "
               "#[\"Proof\"] === Missing[\"NoCertificate\"] &]");
    /* The theorem shown is the claim as written, not a tautology. */
    check_true("leanJu[[1]][\"Theorem\"] === HoldForm[ForAll[n, n >= 2, "
               "Sum[k, {k, 1, n - 2}] == 1/2 (-2 + n) (-1 + n)]]");
    /* The supported shape is still certified. */
    check_true("Conjecture[Sum[k, {k, 1, n}], n, \"Engine\" -> \"Exact\"][\"Status\"] === "
               "\"Checked\"");
}

/* Adversarial review (v0.6528): every verdict-bearing key is authenticated,
 * not only "Verified". */
static void test_review_verdict_keys(void) {
    test_mute_messages();
    assert_eval_eq("leanJa = Conjecture[Sum[k, {k, 1, n}], n, \"Engine\" -> \"Exact\"]; Null",
                   "Null", 0);
    check_true("leanJa[\"Status\"] === \"Checked\"");
    check_true("With[{a = Append[First[leanJa], \"Timing\" -> 0]}, "
               "ProofObject[a][\"Status\"]] === Missing[\"Unauthenticated\"]");
    check_true("With[{a = Append[First[leanJa], \"Status\" -> \"Proved\"]}, "
               "ProofObject[a][\"Status\"]] === Missing[\"Unauthenticated\"]");
    assert_eval_eq("leanJb = Certify[(x+1)^2 == x^2+2x+1, \"Engine\" -> \"Exact\"]; Null",
                   "Null", 0);
    check_true("With[{a = Join[First[leanJb], <|\"Refuted\" -> True, "
               "\"Counterexample\" -> {x -> 1}, \"RefutationTheorem\" -> HoldForm[False]|>]}, "
               "{ProofObject[a][\"Refuted\"], ProofObject[a][\"Counterexample\"], "
               "ProofObject[a][\"RefutationTheorem\"]}] === "
               "ConstantArray[Missing[\"Unauthenticated\"], 3]");
    /* Keys that carry no verdict are still readable. */
    check_true("With[{a = Append[First[leanJa], \"Timing\" -> 0]}, "
               "ProofObject[a][\"Guess\"]] === 1/2 n (1 + n)");
    /* An authentic object reports them. */
    check_true("Refute[x^2 == x, \"Engine\" -> \"Exact\"][\"Counterexample\"] === {x -> -1}");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_polynomial);
    TEST(test_forms);
    TEST(test_verify_false);
    TEST(test_refusals);
    TEST(test_fibonacci_sums);
    TEST(test_review_unsupported_sums);
    TEST(test_review_verdict_keys);
    printf("All Lean conjecture tests passed.\n");
    return 0;
}
