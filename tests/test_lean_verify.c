/*
 * test_lean_verify.c -- Certify / LeanVerify judged by a real Lean 4 kernel.
 *
 * Gated: exits 77 (ctest SKIP_RETURN_CODE) when no Lean 4 executable is
 * found ($LeanVersion is Missing) or its version is below 4.35 (core grind).
 * Discovery is Lean.m's: $LeanExecutable, MATHILDA_LEAN, the newest
 * ~/.elan toolchain, then lean on PATH.
 *
 * The negative cases are the point: a false identity, the PossibleZeroQ
 * false positive, a sorry and a native_decide proof must all be rejected.
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
#ifndef _WIN32
#include <sys/stat.h>
#include <sys/types.h>
#endif

static void check_true(const char* input) {
    assert_eval_eq(input, "True", 0);
}

/* 1 when a Lean >= 4.35 is available. */
static int lean_available(void) {
    Expr* p = parse_expression(
        "With[{v = $LeanVersion}, StringQ[v] && With[{n = StringCases[v, "
        "RegularExpression[\"version (\\\\d+)\\\\.(\\\\d+)\"] -> {\"$1\", \"$2\"}]}, "
        "n =!= {} && OrderedQ[{{4, 35}, ToExpression[First[n]]}]]]");
    Expr* v = evaluate(p);
    expr_free(p);
    int ok = v && v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "True") == 0;
    expr_free(v);
    return ok;
}

static void test_positive(void) {
    check_true("Certify[Factor[x^10 - 1]][\"Verified\"] === True");
    assert_eval_eq("leanV = Certify[(x+1)^3 == x^3+3x^2+3x+1]; Null", "Null", 0);
    check_true("leanV[\"Verified\"] === True");
    check_true("SubsetQ[{\"propext\", \"Classical.choice\", \"Quot.sound\"}, leanV[\"Axioms\"]]");
    check_true("NumericQ[leanV[\"Timing\"]] && leanV[\"Timing\"] < 10");
    check_true("StringQ[leanV[\"LeanVersion\"]]");
    check_true("Certify[Together[1/(x-1) - 1/(x+1)]][\"Verified\"] === True");
    check_true("Certify[Expand[(x/2 + y)^3]][\"Verified\"] === True");
    check_true("ProofObject[Factor[x^4 - 1], Engine -> \"Lean\"][\"Verified\"] === True");
}

static void test_negative(void) {
    test_mute_messages();
    assert_eval_eq("leanN = Certify[(x+1)^3 == x^3+3x^2+3x+2]; Null", "Null", 0);
    check_true("leanN[\"Verified\"] === False");
    check_true("Length[leanN[\"Logs\"]] > 0");
    check_true("MemberQ[leanN[\"Axioms\"], \"sorryAx\"]");
    assert_eval_eq("Check[Certify[(x+1)^2 == x^2 + 1], \"rejected\"]", "\"rejected\"", 0);
    /* Lean catches the CAS: PossibleZeroQ calls this claim true; it is false
     * (the two sides differ by x/10^25). */
    check_true("Certify[(1+x)^20 == Expand[(1+x)^20] - 10^-25 x, Engine -> \"Mathilda\"]"
               "[\"Verified\"] === Missing[\"Heuristic\", True]");
    check_true("Certify[(1+x)^20 == Expand[(1+x)^20] - 10^-25 x, Engine -> \"Lean\"]"
               "[\"Verified\"] === False");
}

static void test_leanverify(void) {
    test_mute_messages();
    check_true("LeanVerify[\"theorem ok : 2 + 2 = 4 := by decide\"][\"Success\"] === True");
    check_true("LeanVerify[\"theorem ok : 2 + 2 = 4 := by decide\"][\"Axioms\"] === "
               "<|\"ok\" -> {}|>");
    check_true("LeanVerify[\"theorem s : 1 + 1 = 3 := by sorry\"][\"Success\"] === False");
    check_true("LeanVerify[\"theorem nd : 2 + 2 = 4 := by native_decide\"][\"Success\"] === False");
    check_true("LeanVerify[\"theorem bad : 2 + 2 = 5 := by decide\"][\"Success\"] === False");
    check_true("With[{r = LeanVerify[\"theorem bad : 2 + 2 = 5 := by decide\"]}, "
               "AnyTrue[r[\"Messages\"], #[\"Severity\"] === \"error\" && #[\"Line\"] === 1 &]]");
    check_true("Lean`Internal`RunProcess[{\"/nonexistent/lean-binary\", \"--version\"}, None]"
               "[\"ExitCode\"] =!= 0");
}

/* Phase 2 kinds: antiderivatives, Bezout, ideal membership, Solve, matrices. */
static void test_phase2_positive(void) {
    assert_eval_eq("leanR = Certify[Integrate[(x^4-3x^2+6)/(x^6-5x^4+5x^2+4), x]]; Null",
                   "Null", 0);
    check_true("leanR[\"Verified\"] === True && leanR[\"Tier\"] === 0");
    check_true("leanR[\"Result\"] === ArcTan[x] + ArcTan[x^3] + ArcTan[(x - 3 x^3 + x^5)/2]");
    check_true("leanR[\"TrustedSteps\"] =!= {} && leanR[\"Timing\"] < 5");
    check_true("Certify[Integrate[1/(x^2-1), x]][\"Verified\"] === True");
    check_true("Certify[Integrate[1/(x^2+2), x]][\"Verified\"] === True");
    check_true("Certify[Integrate[1/(3x^2-2), x]][\"Verified\"] === True");
    check_true("Certify[Integrate[(2x+3)/(x^2+3x+2), x]][\"Verified\"] === True");
    assert_eval_eq("leanB = Certify[PolynomialExtendedGCD[x^3-1, x^2-1, x]]; Null", "Null", 0);
    check_true("leanB[\"Verified\"] === True && leanB[\"Kind\"] === \"Bezout\"");
    check_true("leanB[\"Certificate\"][\"Cofactors\"] == {1, -x}");
    check_true("Certify[PolynomialReduce[x^2 y + x y^2 + y^2, {x y - 1, y^2 - 1}, {x, y}]]"
               "[\"Verified\"] === True");
    check_true("Certify[Implies[x^2+y^2 == 1 && x == y, 2 y^2 == 1]][\"Verified\"] === True");
    assert_eval_eq("leanQ = Certify[Solve[x^3-6x^2+11x-6 == 0, x]]; Null", "Null", 0);
    check_true("leanQ[\"Verified\"] === True && "
               "StringContainsQ[leanQ[\"LeanTheorem\"], \"\xE2\x86\x94\"]");
    check_true("Certify[Solve[(x-1)^2 (2x+1) == 0, x]][\"Verified\"] === True");
    check_true("Certify[Solve[{2x+y == 3, x+3y == 5}, {x, y}]][\"Verified\"] === True");
    check_true("Certify[Inverse[{{1,2,0},{0,1,3},{4,0,1}}]][\"Verified\"] === True");
    check_true("Certify[Det[{{2,1},{1,1}}]][\"Verified\"] === True");
    check_true("Certify[LinearSolve[{{2,1},{1,3}},{3,5}]][\"Verified\"] === True");
    check_true("Certify[NullSpace[{{1,2,3},{2,4,6}}]][\"Verified\"] === True");
}

static void test_phase2_negative(void) {
    test_mute_messages();
    /* A wrong antiderivative (test hook): Mathilda's pre-check fails, the
     * theorem is still emitted, and Lean refutes it. */
    assert_eval_eq("Lean`Private`$antiderivOverride = ArcTan[x] + ArcTan[x^3]; "
                   "leanW = Certify[Integrate[(x^4-3x^2+6)/(x^6-5x^4+5x^2+4), x]]; "
                   "Lean`Private`$antiderivOverride = None; Null", "Null", 0);
    check_true("leanW[\"Verified\"] === False && leanW[\"Certificate\"][\"PreCheck\"] === False");
    /* A perturbed inverse entry. */
    check_true("LeanVerify[StringReplace[Certify[Inverse[{{1,2,0},{0,1,3},{4,0,1}}], "
               "Verify -> False][\"LeanSource\"], \"(12 / 25 : Rat)\" -> \"(13 / 25 : Rat)\"]]"
               "[\"Success\"] === False");
}

/* Phase 3, Tier 0: SOS, optimality, Pratt, sums. */
static void test_phase3_positive(void) {
    check_true("Certify[ForAll[{x,y}, x^4-2x^2+1+y^2 >= 0]][\"Verified\"] === True");
    check_true("Certify[ForAll[{x,y}, x^2+y^2 >= 2 x y]][\"Verified\"] === True");
    check_true("Certify[ForAll[{x}, x >= 0, x^3 + x >= 0]][\"Verified\"] === True");
    assert_eval_eq("leanM = Certify[Minimize[{x^2+y^2, x+y >= 1}, {x,y}]]; Null", "Null", 0);
    check_true("leanM[\"Verified\"] === True && leanM[\"Certificate\"][\"Multipliers\"] === {1}");
    assert_eval_eq("leanLP = Certify[Maximize[{3x+2y, x+y<=4, x+3y<=6, x>=0, y>=0}, {x,y}]]; "
                   "Null", "Null", 0);
    check_true("leanLP[\"Verified\"] === True && leanLP[\"Kind\"] === \"Optimality\" && "
               "AllTrue[leanLP[\"Certificate\"][\"Multipliers\"], # >= 0 &]");
    check_true("Certify[Maximize[{x y, x+y == 2}, {x,y}]][\"Verified\"] === True");
    check_true("Certify[Minimize[x^4-2x^2+1+y^2, {x,y}]][\"Verified\"] === True");
    assert_eval_eq("leanFI = Certify[FactorInteger[2^64+1]]; Null", "Null", 0);
    check_true("leanFI[\"Verified\"] === True && leanFI[\"Timing\"] < 5");
    check_true("StringCount[leanFI[\"LeanSource\"], \"theorem mathilda_cert_1_pratt_\"] === 2");
    check_true("Certify[PrimeQ[1000003]][\"Verified\"] === True");
    assert_eval_eq("leanPC = Certify[PrimeQ[1000001]]; Null", "Null", 0);
    check_true("leanPC[\"Verified\"] === True && leanPC[\"Kind\"] === \"Primality\" && "
               "leanPC[\"Theorem\"] === HoldForm[!PrimeQ[1000001]]");
    check_true("Certify[PrimeQ[2]][\"Verified\"] === True");
    check_true("Certify[FactorInteger[360]][\"Verified\"] === True");
    check_true("Certify[Sum[k 2^k, {k, 0, n}]][\"Verified\"] === True");
    check_true("Certify[Sum[k^3, {k, 1, n}]][\"Verified\"] === True");
    check_true("Certify[Sum[k^2 (1/3)^k, {k, 0, n}]][\"Verified\"] === True");
    check_true("Certify[Sum[2^k, {k, 2, n}]][\"Verified\"] === True");
}

static void test_phase3_negative(void) {
    test_mute_messages();
    /* A Pratt certificate whose witness is a square (never a generator). */
    check_true("With[{s = Certify[PrimeQ[1000003], Verify -> False][\"LeanSource\"], "
               "a = Lean`Private`pratt[1000003][[-1, 2]]}, "
               "LeanVerify[StringReplace[s, \"(1000003, \" <> ToString[a] <> \",\" -> "
               "\"(1000003, \" <> ToString[Mod[a^2, 1000003]] <> \",\"]][\"Success\"]] === False");
    assert_eval_eq("leanSN = Certify[Sum[k^3,{k,1,n}] == n^2 (n+1)^2/4 + 1]; Null", "Null", 0);
    check_true("leanSN[\"Verified\"] === False && Length[leanSN[\"Logs\"]] > 0");
}

/* Row 11: one Lean run for a batch, the session cache, line-span verdicts. */
static void test_batch(void) {
    test_mute_messages();
    assert_eval_eq("Lean`Private`$certCache = <||>; "
                   "leanT1 = First[AbsoluteTiming[Certify[Factor[x^8 - 1]]]]; "
                   "{leanTB, leanBa} = AbsoluteTiming[Certify[{Factor[x^6-1], "
                   "Inverse[{{2,1},{1,1}}], PolynomialExtendedGCD[x^3-1,x^2-1,x]}]]; Null",
                   "Null", 0);
    check_true("Length[leanBa] === 3 && AllTrue[leanBa, #[\"Verified\"] === True &]");
    check_true("NoneTrue[leanBa, #[\"Cached\"] &]");
    check_true("leanTB < 1.5 leanT1");
    assert_eval_eq("leanBb = Certify[{Factor[x^6-1], Inverse[{{2,1},{1,1}}], "
                   "PolynomialExtendedGCD[x^3-1,x^2-1,x]}]; Null", "Null", 0);
    check_true("AllTrue[leanBb, #[\"Cached\"] && #[\"Verified\"] === True &]");
    check_true("leanBb[[1]][\"Timing\"] === leanBa[[1]][\"Timing\"]");
    /* One bad member: only it is False. */
    assert_eval_eq("leanBc = Certify[{Factor[x^5-1], (x+1)^2 == x^2 + 1, Det[{{3,1},{1,1}}]}]; "
                   "Null", "Null", 0);
    check_true("(#[\"Verified\"] & /@ leanBc) === {True, False, True}");
    check_true("Length[leanBc[[2]][\"Logs\"]] > 0 && "
               "NoneTrue[leanBc[[1]][\"Logs\"], #[\"Severity\"] === \"error\" &]");
}

/* Regressions for the adversarial review's false "Verified -> True" results:
 * every reproducer must now be refused or come back False, while the
 * legitimate certificates stay True. */
static void test_review_method_injection(void) {
    test_mute_messages();
    /* Finding 1: a Method whose second line replaced a later theorem by a
     * def (so #print axioms reported nothing) -- now refused outright. */
    assert_eval_eq("Certify[{x + 1 == 1 + x, x == x + 1}, Method -> \"intros\\n  grind\\n"
                   "set_option hygiene false in macro_rules | `(command| set_option "
                   "maxHeartbeats $n:num in $_:command) => `(def mathilda_cert_2 : Nat := 0)\"]",
                   "$Failed", 0);
    /* A plain multi-line script still works; each user-Method item runs in its
     * own file and the false one stays False. */
    assert_eval_eq("leanM = Certify[{x + 1 == 1 + x, x == x + 1}, Method -> \"intros\\ngrind\"]; "
                   "Null", "Null", 0);
    check_true("leanM[[1]][\"Verified\"] === True && leanM[[2]][\"Verified\"] === False");
    check_true("leanM[[1]][\"LeanSource\"] =!= leanM[[2]][\"LeanSource\"]");
    /* The trusted footer: a chunk that declares something else under the
     * theorem's name (here a def of type Nat) fails the statement re-check. */
    check_true("Lean`Private`runBatch[{<|\"Name\" -> \"mathilda_cert_1\", "
               "\"Statement\" -> \"(1 : Nat) = 2\", "
               "\"Text\" -> \"def mathilda_cert_1 : Nat := 0\\n\"|>}, 0, "
               "Lean`Private`leanExe[Automatic], None][[1, \"Verified\"]] === False");
    /* ... and the same chunk with an honest theorem is accepted. */
    check_true("Lean`Private`runBatch[{<|\"Name\" -> \"mathilda_cert_1\", "
               "\"Statement\" -> \"(1 : Nat) = 1\", "
               "\"Text\" -> \"theorem mathilda_cert_1 : (1 : Nat) = 1 := rfl\\n\"|>}, 0, "
               "Lean`Private`leanExe[Automatic], None][[1, \"Verified\"]] === True");
}

static void test_review_domain(void) {
    test_mute_messages();
    /* Finding 2: Assumptions over the automatic Int domain were vacuous. */
    check_true("Certify[x == 0, Assumptions -> 2 x == 1][\"Verified\"] === False");
    check_true("Certify[x == 1, Assumptions -> 0 < x < 2][\"Verified\"] === False");
    check_true("Certify[x^2 == 1, Assumptions -> x^2 <= 1 && x != 0][\"Verified\"] === False");
    check_true("Certify[x^2 == x, Assumptions -> {x == 1}][\"Verified\"] === True");
    assert_eval_eq("leanIM = Certify[Implies[x^2+y^2 == 1 && x == y, 2 y^2 == 1]]; Null",
                   "Null", 0);
    check_true("leanIM[\"Verified\"] === True && leanIM[\"Domain\"] === Complexes");
}

static void test_review_solve(void) {
    test_mute_messages();
    /* Finding 3: Solve's spurious root x = 1 zeroes a denominator; Lean's
     * x / 0 = 0 accepted it. */
    check_true("Certify[Solve[(x-1)/(x^2-1) == (x-1)/(2 x - 2), x]][\"Verified\"] === False");
    check_true("Certify[Solve[(x-1)/(x^2-1) == 1/2, x]][\"Verified\"] === False");
    check_true("Certify[Solve[x^2 + 1/(x-1) == 1 + 1/(x-1), x]][\"Verified\"] === True");
    /* Finding 7: the inequalities / Unequal are part of the statement. */
    assert_eval_eq("leanS7 = Certify[Solve[x^2 == 1 && x > 0, x]]; Null", "Null", 0);
    check_true("leanS7[\"Verified\"] === True && "
               "StringContainsQ[leanS7[\"LeanTheorem\"], \"(1 : Rat) > 0\"]");
    check_true("Certify[Solve[x^2 == 1 && x != 1, x]][\"Verified\"] === True");
    /* A wrong solution set for a constrained system is caught by the
     * inequality conjunct (the equation alone holds at x = -1). */
    check_true("LeanVerify[StringReplace[Certify[Solve[x^2 == 1 && x > 0, x], Verify -> False]"
               "[\"LeanSource\"], \"(1 ^ 2 : Rat) = 1 \xE2\x88\xA7 (1 : Rat) > 0\" -> "
               "\"((-1) ^ 2 : Rat) = 1 \xE2\x88\xA7 (-1 : Rat) > 0\"]][\"Success\"] === False");
}

static void test_review_pratt(void) {
    test_mute_messages();
    /* Finding 4: powMod had 256 bits of fuel, so the composite Fermat number
     * F8 = 2^256 + 1 passed as prime with the certificate [(2, 256)]. */
    check_true("LeanVerify[Lean`Private`prelude[] <> \"\\ntheorem fake_prime : "
               "Mathilda.PrattCertifies (2^256 + 1) [(2^256 + 1, 3, [(2, 256)])] := "
               "by decide +kernel\\n\"][\"Success\"] === False");
    /* Genuine certificates still verify, above 2^256 too (51 * 2^257 + 1 is
     * prime with a three-step certificate; exponents of 258 bits). */
    check_true("Certify[FactorInteger[2^64 + 1]][\"Verified\"] === True");
    check_true("Certify[PrimeQ[1000003]][\"Verified\"] === True");
    check_true("Certify[PrimeQ[51 2^257 + 1]][\"Verified\"] === True");
}

static void test_review_leanverify(void) {
    test_mute_messages();
    /* Finding 5: an axiom of False used from an example / def, and
     * native_decide in an example, reported Success -> True. */
    check_true("LeanVerify[\"axiom bad : False\\nexample : 1 = 2 := bad.elim\\n"
               "def t : 2 + 2 = 5 := bad.elim\\n\"][\"Success\"] === False");
    check_true("LeanVerify[\"example : 10^9 + 7 = 1000000007 := by native_decide\\n\"]"
               "[\"Success\"] === False");
    check_true("LeanVerify[\"theorem a1 : (1:Nat) = 1 := rfl\\nset_option debug.skipKernelTC "
               "true\\n\"][\"Success\"] === False");
    check_true("LeanVerify[\"#eval IO.println \\\"'t' does not depend on any axioms\\\"\\n\"]"
               "[\"Success\"] === False");
    /* Examples are named, so their axioms are checked; a clean one passes. */
    assert_eval_eq("leanLV = LeanVerify[\"example : 1 + 1 = 2 := rfl\\n\"]; Null", "Null", 0);
    check_true("leanLV[\"Success\"] === True && "
               "leanLV[\"Axioms\"] === <|\"mathilda_example_1\" -> {}|>");
    check_true("LeanVerify[\"example : 1 + 1 = 3 := by sorry\\n\"][\"Success\"] === False");
}

static void test_review_tempfiles(void) {
    test_mute_messages();
#ifndef _WIN32
    /* Finding 6: every run uses a private mkdtemp directory (mode 0700). */
    Expr* p = parse_expression(
        "LeanVerify[\"theorem ok : 2 + 2 = 4 := by decide\", KeepFiles -> True][\"SourceFile\"]");
    Expr* v = evaluate(p);
    expr_free(p);
    ASSERT(v && v->type == EXPR_STRING);
    const char* path = v->data.string;
    ASSERT(strstr(path, "mathilda-lean-") != NULL);
    char* dir = malloc(strlen(path) + 1);
    ASSERT(dir != NULL);
    strcpy(dir, path);
    char* slash = strrchr(dir, '/');
    ASSERT(slash != NULL);
    *slash = '\0';
    struct stat st;
    ASSERT(stat(dir, &st) == 0);
    ASSERT(S_ISDIR(st.st_mode));
    ASSERT((st.st_mode & 0777) == 0700);
    ASSERT(remove(path) == 0);
    ASSERT(remove(dir) == 0);
    free(dir);
    expr_free(v);
#endif
}

static void test_review_shapes(void) {
    test_mute_messages();
    /* Finding 8: zipWith truncation -- a too-short vector "solves" a system
     * unless the statement also fixes the shapes. */
    check_true("LeanVerify[Lean`Private`prelude[] <> \"\\ntheorem t : Mathilda.matVec "
               "[[1, 1]] [1] = [1] := by decide +kernel\\n\"][\"Success\"] === True");
    check_true("LeanVerify[Lean`Private`prelude[] <> \"\\ntheorem t : Mathilda.HasShape "
               "[[1, 1]] 1 2 \xE2\x88\xA7 Mathilda.HasLength [1] 2 \xE2\x88\xA7 Mathilda.matVec "
               "[[1, 1]] [1] = [1] := by decide +kernel\\n\"][\"Success\"] === False");
    check_true("Certify[NullSpace[{{1,2,3},{2,4,6}}]][\"Verified\"] === True");
    check_true("Certify[Det[{{1,2,0},{0,1,3},{4,0,1}}]][\"Verified\"] === True");
}


/* Second review round.  Each reproducer must now fail closed. */

/* A source that defeats LeanVerify's textual deny-list: a tactic
 * elaborator, attached with `attribute [tactic ...]`, runs a command built
 * from strings that adds a theorem `1 = 0` under debug.skipKernelTC.  Lean
 * itself reports no error and "does not depend on any axioms"; only the
 * kernel replay of the compiled module catches it. */
static const char* const REPLAY_EXPLOIT =
    "import Lean\n"
    "open Lean Elab Tactic Meta\n"
    "\n"
    "def evilDecide : Tactic := fun _ => do\n"
    "  let g ← getMainGoal\n"
    "  let src := \"set_option de\" ++ \"bug.skip\" ++ \"KernelTC true in\\nrun_\" ++ \"cmd Lean.Elab.Command.liftCoreM <| Lean.add\" ++\n"
    "    \"Decl (.thmDecl { name := `one_eq_zero.evil_aux, levelParams := [], \" ++\n"
    "    \"type := Lean.mkApp3 (Lean.mkConst ``Eq [Lean.Level.one]) (Lean.mkConst ``Nat) \" ++\n"
    "    \"(Lean.mkNatLit 1) (Lean.mkNatLit 0), value := Lean.mkConst ``True.intro })\"\n"
    "  let env ← getEnv\n"
    "  match Parser.runParserCategory env `command src with\n"
    "  | .ok stx => liftCommandElabM (Command.elabCommand stx)\n"
    "  | .error e => throwError e\n"
    "  g.assign (mkConst `one_eq_zero.evil_aux)\n"
    "\n"
    "attribute [tactic Lean.Parser.Tactic.decide] evilDecide\n"
    "\n"
    "theorem one_eq_zero : (1 : Nat) = 0 := by decide\n";

/* The reviewer's lv2.lean: the same attack through @[tactic ...] and addDecl. */
static const char* const REVIEW2_LV2 =
    "import Lean\n"
    "open Lean Elab Tactic Meta\n"
    "\n"
    "@[tactic Lean.Parser.Tactic.decide] def evilDecide : Tactic := fun _ => do\n"
    "  let g ← getMainGoal\n"
    "  let ty ← instantiateMVars (← g.getType)\n"
    "  let nm := (← Term.getDeclName?).getD `anon ++ `aux\n"
    "  let o := Name.mkStr (Name.mkStr .anonymous (\"de\" ++ \"bug\")) \"skipKernelTC\"\n"
    "  withOptions (fun opts => opts.setBool o true) do\n"
    "    addDecl (.thmDecl { name := nm, levelParams := [], type := ty, value := mkConst ``True.intro })\n"
    "  g.assign (mkConst nm)\n"
    "\n"
    "theorem one_eq_zero : (1 : Nat) = 0 := by decide\n";

/* Write text to a fixed file name in the working directory; 1 on success. */
static int write_fixture(const char* path, const char* text) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    fputs(text, f);
    return fclose(f) == 0;
}

static void test_review2_replay(void) {
    test_mute_messages();
    /* Finding A: the guarantee is the kernel replay. */
    ASSERT(write_fixture("lean_review2_replay.lean", REPLAY_EXPLOIT));
    ASSERT(write_fixture("lean_review2_lv2.lean", REVIEW2_LV2));
    assert_eval_eq("leanRP = LeanVerify[ReadString[\"lean_review2_replay.lean\"]]; Null",
                   "Null", 0);
    check_true("!KeyExistsQ[leanRP, \"Rejected\"]");          /* the deny-list missed it */
    check_true("leanRP[\"ExitCode\"] === 0");                  /* Lean's own run was clean */
    check_true("leanRP[\"Success\"] === False");
    check_true("StringQ[leanRP[\"Replay\"]] && "
               "StringContainsQ[leanRP[\"Replay\"], \"declaration type mismatch\"]");
    /* lv2.lean is refused early by the deny-list ... */
    check_true("LeanVerify[ReadString[\"lean_review2_lv2.lean\"]][\"Success\"] === False");
    /* ... and, run past the deny-list, by the replay. */
    assert_eval_eq("leanRP2 = Lean`Internal`VerifySource[ReadString[\"lean_review2_lv2.lean\"], "
                   "Lean`Private`leanExe[Automatic], None, False, Lean`Private`replayScript[]]; Null",
                   "Null", 0);
    check_true("leanRP2[\"ExitCode\"] === 0 && leanRP2[\"Success\"] === False && "
               "StringQ[leanRP2[\"Replay\"]]");
    remove("lean_review2_replay.lean");
    remove("lean_review2_lv2.lean");
    /* The replay's own axiom report: an unused axiom or a sorry anywhere in
     * the module fails, an honest file passes with Replay -> True. */
    assert_eval_eq("leanRP3 = LeanVerify[\"theorem t : 2 + 2 = 4 := by decide\\n\"]; Null",
                   "Null", 0);
    check_true("leanRP3[\"Success\"] === True && leanRP3[\"Replay\"] === True && "
               "leanRP3[\"ReplayAxioms\"][\"t\"] === {}");
    /* Finding A (Mathilda route): the tactic package and its options are refused. */
    check_true("MemberQ[LeanVerify[\"import Mathilda\\ntheorem t : 1 = 1 := rfl\\n\"]"
               "[\"Rejected\"], \"import Mathilda\"]");
    check_true("LeanVerify[\"set_option mathilda.path \\\"/bin/sh\\\"\\ntheorem t : 1 = 1 := rfl\\n\"]"
               "[\"Success\"] === False");
    /* Certify: every Lean run is replayed too; a batch with a rejected chunk
     * (whose file Lean does not compile) re-runs the clean chunks. */
    assert_eval_eq("leanRB = Certify[{x + 0 == x, x == x + 1, 2 x == x + x}]; Null", "Null", 0);
    check_true("(#[\"Verified\"] & /@ leanRB) === {True, False, True}");
}

static void test_review2_domain(void) {
    test_mute_messages();
    /* Finding C: hypotheses no rational satisfies made the statement vacuous
     * over Rat.  Over every char-0 field, x = 2 I (resp. Sqrt[2]) refutes them. */
    check_true("Certify[x == 7, Assumptions -> x^2 + 4 == 0, Method -> \"intro x h\\n"
               "have h2 : 0 \\[LessEqual] x ^ 2 := Mathilda.sqnn\\ngrind\"][\"Verified\"] === False");
    check_true("Certify[x == 0, Assumptions -> x^2 == -1, Method -> \"intro x h\\n"
               "have h2 : 0 \\[LessEqual] x ^ 2 := Mathilda.sqnn\\ngrind\"][\"Verified\"] === False");
    check_true("Certify[x == 7, Assumptions -> x^2 == 2][\"Verified\"] === False");
    /* Honest statements stay True in the generic forms. */
    check_true("Certify[x^2 == 4, Assumptions -> x == 2][\"Verified\"] === True");
    assert_eval_eq("leanGD = Certify[ForAll[{x, y}, x^2 + y^2 >= 2 x y]]; Null", "Null", 0);
    check_true("leanGD[\"Verified\"] === True && leanGD[\"Domain\"] === Reals && "
               "leanGD[\"TrustedSteps\"] === {}");
    check_true("Certify[Minimize[{x^2 + y^2, x + y >= 1}, {x, y}]][\"Verified\"] === True");
}

static void test_review2_low(void) {
    test_mute_messages();
    /* Finding D(i): the unreduced denominators are cross-multiplied, so Lean
     * checks (x^2 - 1) = (x - 1)(x + 1), not 1 + x = 1 + x. */
    assert_eval_eq("leanD1 = Certify[(x^2 - 1)/(x - 1) == x + 1]; Null", "Null", 0);
    check_true("leanD1[\"Verified\"] === True && "
               "leanD1[\"LeanTheorem\"] === \"\\[ForAll] x : Int, x ^ 2 - 1 = (x - 1) * (1 + x)\"");
    check_true("leanD1[\"Hypotheses\"] === {HoldForm[-1 + x != 0]}");
    check_true("Certify[(x^3 - 1)/(x - 1) == x^2 + x + 2][\"Verified\"] === False");
    /* Finding D(ii): not prime / prime by the definition Mathilda.IsPrime. */
    assert_eval_eq("leanD2 = Certify[PrimeQ[0]]; Null", "Null", 0);
    check_true("leanD2[\"Verified\"] === True && leanD2[\"Result\"] === False && "
               "leanD2[\"LeanTheorem\"] === \"\\[Not] Mathilda.IsPrime 0\"");
    check_true("Certify[PrimeQ[1]][\"Verified\"] === True");
    check_true("Certify[PrimeQ[2]][\"Verified\"] === True && "
               "Certify[PrimeQ[2]][\"TrustedSteps\"] === {}");
    check_true("Certify[PrimeQ[1000001]][\"Verified\"] === True");
}

int main(void) {
#ifdef __GNUC__
    /* The shared alarm(120) in test_utils.h is too short for this suite: every
     * Lean run is now also replayed by a fresh kernel (about 1.1 s each, v0.6519),
     * for about 90 runs. */
    alarm(600);
#endif
    symtab_init();
    core_init();
    if (!lean_available()) {
        printf("SKIP: no Lean 4 >= 4.35 found (set $LeanExecutable or MATHILDA_LEAN)\n");
        return 77;
    }
    TEST(test_positive);
    TEST(test_negative);
    TEST(test_leanverify);
    TEST(test_phase2_positive);
    TEST(test_phase2_negative);
    TEST(test_batch);
    TEST(test_phase3_positive);
    TEST(test_phase3_negative);
    TEST(test_review_method_injection);
    TEST(test_review_domain);
    TEST(test_review_solve);
    TEST(test_review_pratt);
    TEST(test_review_leanverify);
    TEST(test_review_tempfiles);
    TEST(test_review_shapes);
    TEST(test_review2_replay);
    TEST(test_review2_domain);
    TEST(test_review2_low);
    printf("All Lean verification tests passed.\n");
    return 0;
}
