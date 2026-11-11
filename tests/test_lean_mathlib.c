/*
 * test_lean_mathlib.c -- Tier 1 (Mathlib) certificates judged by Lean.
 *
 * Gated: exits 77 (ctest SKIP_RETURN_CODE) unless the environment variable
 * MATHILDA_LEAN_PROJECT names a Lake project with Mathlib built for the same
 * toolchain as the Lean executable, and a Lean 4 >= 4.35 is found.  Nothing
 * here depends on a particular project path.  `import Mathlib` costs about
 * 10-40 s per Lean run, so the cases are certified as one batch.
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

static int eval_true(const char* input) {
    Expr* p = parse_expression(input);
    Expr* v = evaluate(p);
    expr_free(p);
    int ok = v && v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "True") == 0;
    expr_free(v);
    return ok;
}

/* Tier 1 HasDerivAt antiderivatives: one batch, one `import Mathlib`. */
static void test_tier1_antiderivatives(void) {
    assert_eval_eq("leanT1 = Certify[{Integrate[(1+x)/(1+x^2), x], Integrate[1/(x^2-1), x], "
                   "Integrate[1/(x^2+2), x]}, Tier -> 1]; Null", "Null", 0);
    check_true("Length[leanT1] === 3 && AllTrue[leanT1, #[\"Verified\"] === True &]");
    check_true("AllTrue[leanT1, #[\"TrustedSteps\"] === {} && #[\"Tier\"] === 1 &]");
    check_true("SubsetQ[{\"propext\", \"Classical.choice\", \"Quot.sound\"}, "
               "leanT1[[1]][\"Axioms\"]]");
    /* Hypotheses x != -1 and x != 1 (through the Log rewrite of ArcTanh). */
    check_true("SubsetQ[leanT1[[2]][\"Hypotheses\"], {HoldForm[1 + x != 0], HoldForm[1 - x != 0]}]");
}

/* A wrong antiderivative is refuted at Tier 1 too. */
static void test_tier1_negative(void) {
    test_mute_messages();
    assert_eval_eq("Lean`Private`$antiderivOverride = ArcTan[x] + Log[1 + x^2]; "
                   "leanT1N = Certify[Integrate[(1+x)/(1+x^2), x], Tier -> 1]; "
                   "Lean`Private`$antiderivOverride = None; Null", "Null", 0);
    check_true("leanT1N[\"Verified\"] === False");
}

/* Third review (v0.6520): with a project, LeanVerify accepts `import Mathlib`
 * (and `lemma`), and still refuses every other import. */
static void test_review3_mathlib_grammar(void) {
    test_mute_messages();
    assert_eval_eq("leanM = LeanVerify[\"import Mathlib\\nlemma t : (2 : \xE2\x84\x9D) + 2 = 4 := by "
                   "norm_num\\n\"]; Null", "Null", 0);
    check_true("leanM[\"Success\"] === True && leanM[\"Replay\"] === True");
    check_true("LeanVerify[\"import Mathlib\\nimport Lean\\ntheorem t : True := trivial\\n\"]"
               "[\"Rejected\"] === {\"import Lean\"}");
    check_true("LeanVerify[\"import Mathlib\\nimport Aesop\\ntheorem t : True := trivial\\n\"]"
               "[\"Rejected\"] === {\"import Aesop\"}");
}

/* Fourth review (v0.6521): the minimal surface with Mathlib.  Ordinary
 * Mathlib proofs verify; Mathlib's own command words and a keyword glued to
 * one of its notation tokens (U+2192 `o`, `[MOD`, ...) are refused unrun. */
static void test_review4_mathlib_surface(void) {
    test_mute_messages();
    assert_eval_eq("leanM4 = LeanVerify[\"import Mathlib\\nlemma m1 (x : \xE2\x84\x9D) : 0 \xE2\x89\xA4 x ^ 2 := by "
                   "positivity\\ntheorem m2 (a b : \xE2\x84\x95) (h : a \xE2\x88\xA3 b) (hb : 0 < b) : a \xE2\x89\xA4 b := "
                   "Nat.le_of_dvd hb h\\n\"]; Null", "Null", 0);
    check_true("leanM4[\"Success\"] === True && leanM4[\"Replay\"] === True");
    check_true("LeanVerify[\"import Mathlib\\nlemma t : True := trivial\\n  x \xE2\x86\x92orun_cmd y\\n\"]"
               "[\"ExitCode\"] === Missing[\"NotRun\"]");
    check_true("LeanVerify[\"import Mathlib\\nlemma t : True := trivial\\n  x [MODtheorem y\\n\"]"
               "[\"ExitCode\"] === Missing[\"NotRun\"]");
    check_true("LeanVerify[\"import Mathlib\\nlemma t : True := trivial\\n  sudo set_option x\\n\"]"
               "[\"ExitCode\"] === Missing[\"NotRun\"]");
    check_true("LeanVerify[\"import Mathlib\\nlemma t : True := by\\n  plausible\\n\"]"
               "[\"ExitCode\"] === Missing[\"NotRun\"]");
}

int main(void) {
#ifdef __GNUC__
    /* The shared alarm(120) in test_utils.h is too short for this suite: each
     * Tier 1 run imports Mathlib (10-38 s) and is then replayed by a fresh
     * kernel under `lake env` (about 9.5 s); the suite takes about 110 s alone. */
    alarm(900);
#endif
    const char* proj = getenv("MATHILDA_LEAN_PROJECT");
    if (!proj || !*proj) {
        printf("SKIP: MATHILDA_LEAN_PROJECT is not set (Tier 1 needs a Mathlib project)\n");
        return 77;
    }
    symtab_init();
    core_init();
    if (!eval_true("With[{v = $LeanVersion}, StringQ[v] && With[{n = StringCases[v, "
                   "RegularExpression[\"version (\\\\d+)\\\\.(\\\\d+)\"] -> {\"$1\", \"$2\"}]}, "
                   "n =!= {} && OrderedQ[{{4, 35}, ToExpression[First[n]]}]]]")) {
        printf("SKIP: no Lean 4 >= 4.35 found\n");
        return 77;
    }
    TEST(test_tier1_antiderivatives);
    TEST(test_tier1_negative);
    TEST(test_review3_mathlib_grammar);
    TEST(test_review4_mathlib_surface);
    printf("All Lean Tier 1 (Mathlib) tests passed.\n");
    return 0;
}
