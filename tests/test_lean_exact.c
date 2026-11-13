/*
 * test_lean_exact.c -- the "Exact" backend of Certify and multi-engine
 * consensus (src/internal/lean/ExactBackend.m, Engines.m).  Needs no Lean
 * installation: Exact re-checks every certificate with Mathilda's own exact
 * arithmetic, and the Lean engine is only listed with Verify -> False.
 *
 * Exact never reports "Verified" -> True: its verdicts are
 * Missing["Exact", True | False] or Missing["NotChecked"].
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

/* Certify[<call>, Engine -> "Exact"]["Verified"] must be `verdict`. */
static void check_exact(const char* call, const char* verdict) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Quiet[Certify[%s, Engine -> \"Exact\"]][\"Verified\"] === %s",
             call, verdict);
    check_true(buf);
}

static void test_exact_object(void) {
    assert_eval_eq("leanX = Certify[(x+1)^3 == x^3+3x^2+3x+1, Engine -> \"Exact\"]; Null",
                   "Null", 0);
    check_true("leanX[\"Properties\"] === {\"Cached\", \"Certificate\", \"Checks\", \"Domain\", "
               "\"Engine\", \"Hypotheses\", \"Input\", \"Kind\", \"Result\", \"Theorem\", "
               "\"Tier\", \"Timing\", \"TrustedSteps\", \"Verified\"}");
    check_true("leanX[\"Verified\"] === Missing[\"Exact\", True]");
    check_true("leanX[\"Tier\"] === \"Exact\" && leanX[\"Engine\"] === \"Exact\"");
    check_true("leanX[\"Kind\"] === \"Identity\" && leanX[\"Domain\"] === Integers");
    check_true("StringContainsQ[ToString[leanX], \"| Exact |\"]");
    check_true("MatchQ[leanX[\"Checks\"], {_String -> True}]");
    check_true("StringStartsQ[First[leanX[\"TrustedSteps\"]], \"Mathilda's exact rational "
               "arithmetic\"]");
    check_true("leanX[\"Cached\"] === False");
}

static void test_exact_true(void) {
    static const char* items[] = {
        /* identities: Int, Rat, cleared denominators, equality assumptions */
        "(x+1)^3 == x^3+3x^2+3x+1",
        "Expand[(x/2 + 1)^2]",
        "Together[1/(x-1) - 1/(x+1)]",
        "Factor[x^10 - 1]",
        "x^2 == 1, Assumptions -> x == 1",
        "x^3 == x, Assumptions -> {x^2 == 1, y == 2}",
        /* Bezout, PolynomialReduce */
        "PolynomialExtendedGCD[x^4 - 1, x^3 + 2 x^2 + 2 x + 1, x]",
        "PolynomialReduce[x^3 y - x y^2, {x^2 - y, x y - 1}, {x, y}]",
        /* Tier 0 antiderivatives (one with a surd), the claim form */
        "Integrate[1/(x^2 + 2), x]",
        "Integrate[(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4), x]",
        "Integrate[2 x, x] == x^2 + 5",
        /* membership */
        "Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1]",
        /* Solve completeness (multiple root) and soundness */
        "Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x]",
        "Solve[(x - 1)^2 (x + 2) == 0, x]",
        "Solve[{x + y == 3, x - y == 1}, {x, y}]",
        /* matrices */
        "Inverse[{{2, 1, 0}, {1, 3, 1}, {0, 1, 4}}]",
        "Det[{{2, -1, 0, 3}, {1, 4, 2, 0}, {0, 5, -2, 1}, {3, 0, 1, 1}}]",
        "LinearSolve[{{2, 1}, {1, 3}}, {1, 2}]",
        "NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
        /* SOS: Quadratic, SOSCertificate, ExactGram; a hypothesis multiplier */
        "ForAll[{x, y}, x^2 + y^2 >= 2 x y]",
        "ForAll[{x, y}, x^4 + y^4 - x^2 y^2 + 1 >= 0]",
        "ForAll[{x, y}, x^4 - 2 x^2 + 1 + y^2 >= 0]",
        "Implies[x >= 0, x^3 >= 0]",
        /* optimality: QP, LP, an eliminated equality */
        "Minimize[{x^2 + y^2, x + y >= 1}, {x, y}]",
        "Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}]",
        "Minimize[{x^2 + y^2, x + y == 1}, {x, y}]",
        /* primality */
        "FactorInteger[2^64 + 1]",
        "PrimeQ[2]", "PrimeQ[274177]", "PrimeQ[274179]", "PrimeQ[1]",
        /* sums */
        "Sum[k^3, {k, 1, n}]",
        "Sum[k 2^k, {k, 0, n}]",
        "Sum[(-1)^k k, {k, 0, n}]",
        /* the Tier 1 shape (HasDerivAt): D[F, x] == f by Together, no Mathlib */
        "Integrate[(1 + x)/(1 + x^2), x]",
        "Integrate[(1 + x)/(1 + x^2), x], \"Tier\" -> 1",
        NULL
    };
    for (int i = 0; items[i]; i++) check_exact(items[i], "Missing[\"Exact\", True]");
    check_true("Quiet[Certify[Integrate[(1 + x)/(1 + x^2), x], Engine -> \"Exact\", "
               "\"Tier\" -> 1]][\"Certificate\"][\"Form\"] === \"HasDerivAt\"");
    check_true("Certify[PrimeQ[274177], Engine -> \"Exact\"][\"TrustedSteps\"][[-1]] === "
               "\"Lucas's primality theorem (Mathilda replays the Pratt certificate with "
               "PowerMod; 2 is the base of the scheme)\"");
}

static void test_exact_false(void) {
    test_mute_messages();
    /* Showcase row 17: the PossibleZeroQ false positive. */
    assert_eval_eq("Check[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, "
                   "Engine -> \"Exact\"], \"exactfail\"]", "\"exactfail\"", 0);
    check_exact("(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x", "Missing[\"Exact\", False]");
    check_true("MemberQ[Quiet[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, "
               "Engine -> \"Exact\"]][\"Checks\"][[All, 2]], False]");
    /* A wrong antiderivative through the test hook. */
    assert_eval_eq("Lean`Private`$antiderivOverride = ArcTan[x] + ArcTan[x^3]; "
                   "leanXW = Check[Certify[Integrate[(x^4-3x^2+6)/(x^6-5x^4+5x^2+4), x], "
                   "Engine -> \"Exact\"], \"exactfail\"]; "
                   "Lean`Private`$antiderivOverride = None; leanXW", "\"exactfail\"", 0);
    assert_eval_eq("Lean`Private`$antiderivOverride = ArcTan[x] + ArcTan[x^3]; "
                   "leanXW = Quiet[Certify[Integrate[(x^4-3x^2+6)/(x^6-5x^4+5x^2+4), x], "
                   "Engine -> \"Exact\"]][\"Verified\"]; "
                   "Lean`Private`$antiderivOverride = None; leanXW",
                   "Missing[\"Exact\", False]", 0);
    /* Concrete failures of other kinds. */
    check_exact("x^2 == x", "Missing[\"Exact\", False]");
    check_exact("Sum[k, {k, 1, n}] == n (n + 1)/2 + 1", "Missing[\"Exact\", False]");
    check_exact("Solve[(x - 1)/(x^2 - 1) == (x - 1)/(2 x - 2), x]", "Missing[\"Exact\", False]");
    /* An inequality hypothesis: no decision procedure, not False. */
    check_exact("x^2 == x, Assumptions -> {x > 0}", "Missing[\"NotChecked\"]");
    check_true("StringContainsQ[First[First[Certify[x^2 == x, Assumptions -> {x > 0}, "
               "Engine -> \"Exact\"][\"Checks\"]]], \"no exact decision procedure\"]");
    /* Verify -> False runs nothing. */
    check_true("Certify[x + 1 == 1 + x, Engine -> \"Exact\", Verify -> False][\"Verified\"] === "
               "Missing[\"NotChecked\"]");
    /* The exponential-polynomial zero test. */
    check_true("Lean`Private`expPolyZeroQ[2^(n + 1) - 2 2^n, n]");
    check_true("Lean`Private`expPolyZeroQ[(-1)^(n + 2) - (-1)^n + 4^n - 2^(2 n), n]");
    check_true("!Lean`Private`expPolyZeroQ[3^n - 2^n, n]");
    check_true("Lean`Private`exactLaplaceDet[{{2, -1, 0, 3}, {1, 4, 2, 0}, {0, 5, -2, 1}, "
               "{3, 0, 1, 1}}] === Det[{{2, -1, 0, 3}, {1, 4, 2, 0}, {0, 5, -2, 1}, {3, 0, 1, 1}}]");
}

static void test_multi_engine(void) {
    assert_eval_eq("leanM = Certify[(x+1)^3 == x^3+3x^2+3x+1, Engine -> {\"Lean\", \"Exact\"}, "
                   "Verify -> False]; Null", "Null", 0);
    check_true("leanM[\"Verified\"] === Missing[\"NotChecked\"]");
    check_true("leanM[\"Engine\"] === {\"Lean\", \"Exact\"}");
    check_true("leanM[\"Engines\"] === <|\"Lean\" -> Missing[\"NotChecked\"], "
               "\"Exact\" -> Missing[\"NotChecked\"]|>");
    check_true("leanM[\"Properties\"] === {\"Axioms\", \"Cached\", \"Certificate\", \"Domain\", "
               "\"Engine\", \"EngineDetails\", \"Engines\", \"Hypotheses\", \"Input\", \"Kind\", "
               "\"LeanSource\", \"LeanTheorem\", \"LeanVersion\", \"Logs\", \"ProofScript\", "
               "\"Result\", \"Theorem\", \"Tier\", \"Timing\", \"TrustedSteps\", \"Verified\"}");
    check_true("leanM[\"Tier\"] === 0 && Sort[Keys[leanM[\"EngineDetails\"]]] === "
               "{\"Exact\", \"Lean\"}");
    check_true("KeyExistsQ[leanM[\"EngineDetails\"][\"Exact\"], \"Checks\"]");
    check_true("StringContainsQ[ToString[leanM], \"| Lean Tier 0 + Exact |\"]");
    /* No kernel: consensus of the exact check (the heuristic only adds its entry). */
    assert_eval_eq("leanM2 = Certify[(x+1)^3 == x^3+3x^2+3x+1, Engine -> {\"Exact\", "
                   "\"Mathilda\"}]; Null", "Null", 0);
    check_true("leanM2[\"Verified\"] === Missing[\"Consensus\", True]");
    check_true("leanM2[\"Engines\"] === <|\"Exact\" -> Missing[\"Exact\", True], "
               "\"Mathilda\" -> Missing[\"Heuristic\", True]|>");
    check_true("leanM2[\"Tier\"] === \"Exact\" && KeyExistsQ[First[leanM2], \"Checks\"]");
    check_true("MemberQ[leanM2[\"TrustedSteps\"], "
               "\"Mathilda's PossibleZeroQ / Simplify (a heuristic check, not a proof)\"]");
    /* A one-element list gives the multi-engine shape. */
    check_true("KeyExistsQ[First[Certify[x == x, Engine -> {\"Exact\"}]], \"Engines\"]");
    /* Row 17's claim: the exact engine fails, the multi-engine verdict is False. */
    check_true("Quiet[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, "
               "Engine -> {\"Exact\", \"Mathilda\"}]][\"Verified\"] === False");
    /* Engine values: a list with a duplicate or an unknown name is refused. */
    test_mute_messages();
    assert_eval_eq("Check[Certify[x == x, Engine -> {\"Lean\", \"Lean\"}], \"engine\"]",
                   "\"engine\"", 0);
    assert_eval_eq("Check[Certify[x == x, Engine -> {}], \"engine\"]", "\"engine\"", 0);
    assert_eval_eq("Check[Certify[x == x, Engine -> {\"Exact\", \"Z3\"}], \"engine\"]",
                   "\"engine\"", 0);
}

static void test_consensus_rules(void) {
    test_mute_messages();
    assert_eval_eq("leanG = <|\"Input\" -> HoldForm[a == b]|>; Null", "Null", 0);
#define CONS(vs, engines) "Check[Lean`Private`consensus[leanG, " vs ", " engines "], \"disagree\"]"
    /* Rule 1: a False wins; no disagreement when the kernel agrees. */
    check_true(CONS("{False, Missing[\"Exact\", False]}", "{\"Lean\", \"Exact\"}") " === False");
    check_true(CONS("{Missing[\"NotAvailable\"], Missing[\"Exact\", False]}",
                    "{\"Lean\", \"Exact\"}") " === False");
    /* ... and both disagreement triggers. */
    check_true(CONS("{True, Missing[\"Exact\", False]}", "{\"Lean\", \"Exact\"}")
               " === \"disagree\"");
    check_true(CONS("{False, Missing[\"Exact\", True]}", "{\"Lean\", \"Exact\"}")
               " === \"disagree\"");
    check_true("Quiet[Lean`Private`consensus[leanG, {True, Missing[\"Exact\", False]}, "
               "{\"Lean\", \"Exact\"}]] === False");
    /* Rule 2: every kernel True. */
    check_true(CONS("{True, Missing[\"Exact\", True]}", "{\"Lean\", \"Exact\"}") " === True");
    check_true(CONS("{True, Missing[\"NotChecked\"]}", "{\"Lean\", \"Exact\"}") " === True");
    /* Rule 3: no kernel. */
    check_true(CONS("{Missing[\"Exact\", True]}", "{\"Exact\"}")
               " === Missing[\"Consensus\", True]");
    check_true(CONS("{Missing[\"NotChecked\"]}", "{\"Exact\"}") " === Missing[\"NotChecked\"]");
    /* Rule 4: the first missing kernel verdict. */
    check_true(CONS("{Missing[\"NotAvailable\"], Missing[\"Exact\", True]}",
                    "{\"Lean\", \"Exact\"}") " === Missing[\"NotAvailable\"]");
    /* The heuristic only contributes its entry. */
    check_true(CONS("{False, Missing[\"Heuristic\", True]}", "{\"Lean\", \"Mathilda\"}")
               " === False");
    check_true(CONS("{True, Missing[\"Heuristic\", False]}", "{\"Lean\", \"Mathilda\"}")
               " === True");
#undef CONS
}

static void test_guard_and_options(void) {
    test_mute_messages();
    /* A backend that is not a kernel cannot report True. */
    assert_eval_eq("Lean`Private`registerBackend[\"ZZTest\", <|\"Name\" -> \"ZZTest\", "
                   "\"Trust\" -> \"Heuristic\", \"Encode\" -> (<||> &), "
                   "\"Run\" -> (Table[<|\"Verified\" -> True, \"Timing\" -> 0., \"Logs\" -> {}|>, "
                   "{Length[#1]}] &), \"Present\" -> Function[{item, v, o}, "
                   "With[{a = <|\"Verified\" -> v[\"Verified\"], \"Timing\" -> 0.|>}, "
                   "ProofObject[a]]], \"Label\" -> (\"ZZ\" &)|>]", "\"ZZTest\"", 0);
    check_true("First[Certify[x == x, Engine -> \"ZZTest\"]][\"Verified\"] === "
               "Missing[\"ZZTest\", True]");
    check_true("Certify[x == x, Engine -> {\"ZZTest\"}][\"Engines\"] === "
               "<|\"ZZTest\" -> Missing[\"ZZTest\", True]|>");
    check_true("Certify[x == x, Engine -> {\"ZZTest\"}][\"Verified\"] === "
               "Missing[\"ZZTest\", True]");
    assert_eval_eq("Lean`Private`$backends = KeyDrop[Lean`Private`$backends, \"ZZTest\"]; "
                   "KeyExistsQ[Lean`Private`$backends, \"ZZTest\"]", "False", 0);
    /* Engine-only options. */
    assert_eval_eq("Check[Certify[x == x, Engine -> \"Exact\", Method -> \"omega\"], "
                   "\"engineopt\"]", "\"engineopt\"", 0);
    assert_eval_eq("Check[Certify[x == x, Engine -> \"Exact\", \"Tier\" -> 1], \"engineopt\"]",
                   "\"engineopt\"", 0);
    check_true("Check[Certify[x == x, Engine -> \"Exact\", \"LeanExecutable\" -> "
               "\"/nonexistent\", TimeConstraint -> 5][\"Verified\"], \"msg\"] === "
               "Missing[\"Exact\", True]");
    check_true("Check[Certify[x == x, Engine -> {\"Lean\", \"Exact\"}, Method -> \"omega\", "
               "Verify -> False][\"Kind\"], \"msg\"] === \"Identity\"");
    /* An Exact-only call issues no Lean message (no Certify::tier for a Tier 1 shape). */
    check_true("Check[Certify[Integrate[(1 + x)/(1 + x^2), x], Engine -> \"Exact\"]"
               "[\"Verified\"], \"msg\"] === Missing[\"Exact\", True]");
}

static void test_timing(void) {
    check_true("leanXT = AbsoluteTiming[Quiet[Certify[{(x+1)^3 == x^3+3x^2+3x+1, "
               "Together[1/(x-1) - 1/(x+1)], "
               "PolynomialExtendedGCD[x^4 - 1, x^3 + 2 x^2 + 2 x + 1, x], "
               "PolynomialReduce[x^3 y - x y^2, {x^2 - y, x y - 1}, {x, y}], "
               "Integrate[1/(x^2 + 2), x], Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1], "
               "Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x], Solve[{x + y == 3, x - y == 1}, {x, y}], "
               "Inverse[{{2, 1, 0}, {1, 3, 1}, {0, 1, 4}}], "
               "Det[{{2, -1, 0, 3}, {1, 4, 2, 0}, {0, 5, -2, 1}, {3, 0, 1, 1}}], "
               "LinearSolve[{{2, 1}, {1, 3}}, {1, 2}], NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}], "
               "ForAll[{x, y}, x^2 + y^2 >= 2 x y], ForAll[{x, y}, x^4 - 2 x^2 + 1 + y^2 >= 0], "
               "Minimize[{x^2 + y^2, x + y >= 1}, {x, y}], "
               "Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}], "
               "FactorInteger[2^64 + 1], PrimeQ[274177], PrimeQ[274179], Sum[k^3, {k, 1, n}], "
               "Sum[k 2^k, {k, 0, n}], Integrate[(1 + x)/(1 + x^2), x]}, Engine -> \"Exact\"]]]; "
               "AllTrue[leanXT[[2]], #[\"Verified\"] === Missing[\"Exact\", True] &]");
    /* Wall clock, hang guard under ctest -j load (about 0.1 s idle). */
    check_true("leanXT[[1]] < 2");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_exact_object);
    TEST(test_exact_true);
    TEST(test_exact_false);
    TEST(test_multi_engine);
    TEST(test_consensus_rules);
    TEST(test_guard_and_options);
    TEST(test_timing);
    printf("All Lean Exact-engine tests passed.\n");
    return 0;
}
