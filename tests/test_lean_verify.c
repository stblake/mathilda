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
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "lean_check_corpus.h"

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
    check_true("NumericQ[leanV[\"Timing\"]] && leanV[\"Timing\"] < 60"); /* hang guard only: wall clock under ctest -j load */
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
    check_true("leanR[\"TrustedSteps\"] =!= {} && leanR[\"Timing\"] < 60"); /* hang guard only */
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
    check_true("leanRan[LeanVerify[leanNoFooter[StringReplace[Certify[Inverse[{{1,2,0},{0,1,3},"
               "{4,0,1}}], Verify -> False][\"LeanSource\"], \"(12 / 25 : Rat)\" -> "
               "\"(13 / 25 : Rat)\"]]]][\"Success\"] === False");
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
    check_true("leanFI[\"Verified\"] === True && leanFI[\"Timing\"] < 60"); /* hang guard only */
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
               "leanRan[LeanVerify[leanNoFooter[StringReplace[s, \"(1000003, \" <> ToString[a] <> "
               "\",\" -> \"(1000003, \" <> ToString[Mod[a^2, 1000003]] <> \",\"]]]][\"Success\"]] "
               "=== False");
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
    check_true("leanRan[LeanVerify[leanNoFooter[StringReplace[Certify[Solve[x^2 == 1 && x > 0, x], "
               "Verify -> False][\"LeanSource\"], \"(1 ^ 2 : Rat) = 1 \xE2\x88\xA7 (1 : Rat) > 0\" -> "
               "\"((-1) ^ 2 : Rat) = 1 \xE2\x88\xA7 (-1 : Rat) > 0\"]]]][\"Success\"] === False");
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
    /* Since v0.6520 the source is refused before Lean runs: `import Lean`
     * is outside LeanVerify's allowlist grammar (third review). */
    assert_eval_eq("leanRP = LeanVerify[ReadString[\"lean_review2_replay.lean\"]]; Null",
                   "Null", 0);
    check_true("leanRP[\"Rejected\"] === {\"import Lean\"}");
    check_true("leanRP[\"Success\"] === False && leanRP[\"ExitCode\"] === Missing[\"NotRun\"]");
    check_true("LeanVerify[ReadString[\"lean_review2_lv2.lean\"]][\"Success\"] === False");
    /* Run past the source check (straight to the runner), the kernel replay
     * still rejects both. */
    assert_eval_eq("leanRP1 = Lean`Internal`VerifySource[ReadString[\"lean_review2_replay.lean\"], "
                   "Lean`Private`leanExe[Automatic], None, False, Lean`Private`replayScript[]]; Null",
                   "Null", 0);
    check_true("leanRP1[\"ExitCode\"] === 0 && leanRP1[\"Success\"] === False && "
               "StringQ[leanRP1[\"Replay\"]] && "
               "StringContainsQ[leanRP1[\"Replay\"], \"declaration type mismatch\"]");
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



/* Third review round.  Each reproducer must now fail closed. */

/* The reviewer's p1/pwn.lean: a user tactic elaborator adds `bogus` with
 * addDeclCore (doCheck := false) and writes forged .olean.server /
 * .olean.private files next to the module, which the v0.6519 replay read
 * instead of the real one: Success -> True for (1 : Nat) = 0. */
static const char* const REVIEW3_PWN =
    "import Lean\n"
    "open Lean Meta Elab Tactic\n"
    "set_option Elab.async false\n"
    "\n"
    "def pwn : Tactic := fun _ => do\n"
    "  let g ← getMainGoal\n"
    "  let ty ← g.getType\n"
    "  let nm := `bogus.aux\n"
    "  let decl := Declaration.thmDecl { name := nm, levelParams := [], type := ty, value := mkConst ``True.intro }\n"
    "  let env ← getEnv\n"
    "  match env.addDeclCore 0 1000 decl none (doCheck := false) with\n"
    "  | .ok env' => modifyThe Core.State fun s => { s with env := env' }\n"
    "  | .error _ => throwError \"x\"\n"
    "  -- plant the module parts the replay checker reads in place of the real one\n"
    "  let fn ← getFileName\n"
    "  let base : System.FilePath := (System.FilePath.mk fn).withExtension \"olean\"\n"
    "  let okThm : ConstantInfo := .thmInfo { name := `bogus, levelParams := [], type := mkConst ``True, value := mkConst ``True.intro }\n"
    "  let okDef : ConstantInfo := .defnInfo { name := `pwn, levelParams := [], type := mkConst ``True, value := mkConst ``True.intro, hints := .opaque, safety := .safe }\n"
    "  let md : ModuleData := { isModule := false, imports := #[], constNames := #[`pwn, `bogus], constants := #[okDef, okThm], extraConstNames := #[], entries := #[] }\n"
    "  saveModuleData (base.addExtension \"server\") `X md\n"
    "  saveModuleData (base.addExtension \"private\") `Y md\n"
    "  g.assign (mkConst nm)\n"
    "attribute [tactic Lean.Parser.Tactic.decide] pwn\n"
    "\n"
    "theorem bogus : (1 : Nat) = 0 := by decide\n";

/* The reviewer's p2/race.lean: the same tactic spawns a background shell
 * that overwrites the .olean after Lean wrote it and before the replay
 * read it.  It must be refused before Lean runs: never run it here. */
static const char* const REVIEW3_RACE =
    "import Lean\n"
    "open Lean Meta Elab Tactic\n"
    "set_option Elab.async false\n"
    "\n"
    "def pwn : Tactic := fun _ => do\n"
    "  let g ← getMainGoal\n"
    "  let ty ← g.getType\n"
    "  let nm := `bogus.aux\n"
    "  let decl := Declaration.thmDecl { name := nm, levelParams := [], type := ty, value := mkConst ``True.intro }\n"
    "  match (← getEnv).addDeclCore 0 1000 decl none (doCheck := false) with\n"
    "  | .ok env' => modifyThe Core.State fun s => { s with env := env' }\n"
    "  | .error _ => throwError \"x\"\n"
    "  let fn ← getFileName\n"
    "  let base : System.FilePath := (System.FilePath.mk fn).withExtension \"olean\"\n"
    "  let forged : System.FilePath := (System.FilePath.mk fn).withExtension \"forged\"\n"
    "  let okThm : ConstantInfo := .thmInfo { name := `bogus, levelParams := [], type := mkConst ``True, value := mkConst ``True.intro }\n"
    "  let okDef : ConstantInfo := .defnInfo { name := `pwn, levelParams := [], type := mkConst ``True, value := mkConst ``True.intro, hints := .opaque, safety := .safe }\n"
    "  let md : ModuleData := { isModule := false, imports := #[{ module := `Init }], constNames := #[`pwn, `bogus], constants := #[okDef, okThm], extraConstNames := #[], entries := #[] }\n"
    "  saveModuleData forged `Z md\n"
    "  let cmd := s!\"while [ ! -s '{base}' ]; do sleep 0.02; done; sleep 0.3; cp '{forged}' '{base}'\"\n"
    "  let _ ← IO.Process.spawn { cmd := \"/bin/sh\", args := #[\"-c\", cmd], stdin := .null, stdout := .null, stderr := .null }\n"
    "  g.assign (mkConst nm)\n"
    "attribute [tactic Lean.Parser.Tactic.decide] pwn\n"
    "\n"
    "theorem bogus : (1 : Nat) = 0 := by decide\n";

static void test_review3_allowlist(void) {
    test_mute_messages();
    /* Finding A: both attacks are refused by the allowlist grammar, unrun
     * (so the race's background process is never spawned). */
    ASSERT(write_fixture("lean_review3_pwn.lean", REVIEW3_PWN));
    ASSERT(write_fixture("lean_review3_race.lean", REVIEW3_RACE));
    assert_eval_eq("leanP1 = LeanVerify[ReadString[\"lean_review3_pwn.lean\"]]; Null", "Null", 0);
    check_true("leanP1[\"Success\"] === False && leanP1[\"ExitCode\"] === Missing[\"NotRun\"] && "
               "leanP1[\"Rejected\"] === {\"import Lean\"}");
    assert_eval_eq("Check[LeanVerify[ReadString[\"lean_review3_pwn.lean\"]], \"import\"]",
                   "\"import\"", 0);
    assert_eval_eq("leanP2 = LeanVerify[ReadString[\"lean_review3_race.lean\"]]; Null", "Null", 0);
    check_true("leanP2[\"Success\"] === False && leanP2[\"ExitCode\"] === Missing[\"NotRun\"]");
    /* Even with a project, only Mathlib may be imported. */
    check_true("LeanVerify[ReadString[\"lean_review3_pwn.lean\"], LeanProject -> "
               "\"/nonexistent\"][\"Rejected\"] === {\"import Lean\"}");
    /* Defence in depth: run pwn.lean straight through the runner (it spawns
     * nothing).  Its forged module parts appear in the private directory,
     * so the replay is refused; and the directory is removed afterwards,
     * forged files and all (finding E). */
    assert_eval_eq("leanP3 = Lean`Internal`VerifySource[ReadString[\"lean_review3_pwn.lean\"], "
                   "Lean`Private`leanExe[Automatic], None, True, Lean`Private`replayScript[]]; Null",
                   "Null", 0);
    check_true("leanP3[\"Success\"] === False && StringQ[leanP3[\"Replay\"]] && "
               "StringContainsQ[leanP3[\"Replay\"], \"unexpected file\"]");
#ifndef _WIN32
    {
        Expr* p = parse_expression("leanP3[\"SourceFile\"]");
        Expr* v = evaluate(p);
        expr_free(p);
        ASSERT(v && v->type == EXPR_STRING);
        char* dir = malloc(strlen(v->data.string) + 1);
        ASSERT(dir != NULL);
        strcpy(dir, v->data.string);
        char* slash = strrchr(dir, '/');
        ASSERT(slash != NULL);
        *slash = '\0';
        /* KeepFiles keeps exactly the source; every other file is gone. */
        DIR* d = opendir(dir);
        ASSERT(d != NULL);
        int others = 0;
        struct dirent* de;
        while ((de = readdir(d)) != NULL)
            if (strcmp(de->d_name, ".") && strcmp(de->d_name, "..") &&
                strcmp(de->d_name, "MathildaCert.lean"))
                others++;
        closedir(d);
        ASSERT(others == 0);
        ASSERT(remove(v->data.string) == 0);
        ASSERT(remove(dir) == 0);
        free(dir);
        expr_free(v);
    }
#endif
    remove("lean_review3_pwn.lean");
    remove("lean_review3_race.lean");
    /* Accepted sources still verify: examples in a section, numeric
     * set_option ... in (open ... in and #print axioms are refused since
     * v0.6521, open and namespace since v0.6522). */
    check_true("LeanVerify[\"section A\\nexample : 1 + 1 = 2 := rfl\\n"
               "set_option maxHeartbeats 400000 in\\ntheorem t : True := by\\n  trivial\\n"
               "end A\\n\"][\"Success\"] === True");
    /* Refused, unrun: user definitions, attributes, commands and escapes. */
    check_true("LeanVerify[\"def f : Nat := 1\\ntheorem t : f = 1 := rfl\\n\"][\"Rejected\"] === "
               "{\"def\"}");
    check_true("LeanVerify[\"theorem t : 1 = 1 := by\\n  decide +native\\n\"][\"Rejected\"] === "
               "{\"native\"}");
    check_true("LeanVerify[\"theorem t : 1 = 1 := rfl\\n  #eval 1\\n\"][\"Rejected\"] === "
               "{\"#eval\"}");
    check_true("LeanVerify[\"module\\ntheorem t : 1 + 1 = 2 := rfl\\n\"][\"Rejected\"] === "
               "{\"module\"}");
    check_true("LeanVerify[\"import Mathlib\\ntheorem t : 1 = 1 := rfl\\n\", LeanProject -> None]"
               "[\"Rejected\"] === {\"import Mathlib\"}");
}

static void test_review3_domains(void) {
    test_mute_messages();
    /* Finding B: an order relation makes only its own variables real. */
    check_true("Certify[x == 7, Assumptions -> x^2 + 4 == 0 && y > 0, Method -> \"intro x y h\\n"
               "have h2 : 0 \\[LessEqual] x ^ 2 := Mathilda.sqnn\\ngrind\"] === $Failed");
    /* Declared real, the statement is the faithful ordered-field one (and
     * true: no real x has x^2 + 4 = 0). */
    assert_eval_eq("leanB = Certify[x == 7, Assumptions -> x^2 + 4 == 0 && Element[x, Reals], "
                   "Method -> \"intro x h\\nhave h2 : 0 \\[LessEqual] x ^ 2 := Mathilda.sqnn\\n"
                   "grind\"]; Null", "Null", 0);
    check_true("leanB[\"Verified\"] === True && leanB[\"Domain\"] === Reals");
    /* Finding C: the claim is translated as written. */
    check_true("Certify[Simplify[(x^2 - 1)/(x - 1)] == x + 1] === $Failed");
    check_true("Certify[Factor[x^2 - 1] == (x - 1) (x + 1)][\"Verified\"] === True");
    check_true("Certify[Integrate[x^2, x] == x^3/3][\"Verified\"] === True");
    check_true("Certify[Integrate[x^2, x] == x^3/2][\"Verified\"] === False");
    check_true("Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x][\"Verified\"] === False");
}

/* Fourth review round (v0.6521). */

/* Evaluate `fmt` with %s replaced by a Mathilda string literal of src. */
static Expr* eval_with_source(const char* fmt, const char* src) {
    char* lit = lean_mstring(src);
    size_t n = strlen(lit) + strlen(fmt) + 1;
    char* in = malloc(n);
    if (!in) abort();
    snprintf(in, n, fmt, lit);
    Expr* e = parse_expression(in);
    Expr* r = e ? evaluate(e) : NULL;
    if (e) expr_free(e);
    free(in);
    free(lit);
    return r;
}

static bool eval_true_with_source(const char* fmt, const char* src) {
    Expr* r = eval_with_source(fmt, src);
    bool ok = r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    if (r) expr_free(r);
    return ok;
}

static void test_review4_verify(void) {
    test_mute_messages();
    /* Findings 1-3: the reproducers (x2 gave "Success" -> True for
     * `theorem bad : False` in v0.6520) are refused before Lean runs. */
    const char* xs[] = { REVIEW4_X1, REVIEW4_X2, REVIEW4_X3, REVIEW4_X4, NULL };
    for (int i = 0; xs[i]; i++)
        ASSERT(eval_true_with_source("With[{r = LeanVerify[%s]}, r[\"Success\"] === False && "
                                     "r[\"ExitCode\"] === Missing[\"NotRun\"] && "
                                     "Length[r[\"Rejected\"]] === 1]", xs[i]));
    /* The whole corpus, as one source, is proved by Lean (and each source on
     * its own is accepted by the checker: test_lean_certify.c). */
    size_t len = 1;
    for (int i = 0; LEAN_CORPUS[i]; i++) len += strlen(LEAN_CORPUS[i]);
    char* all = malloc(len);
    ASSERT(all != NULL);
    all[0] = '\0';
    for (int i = 0; LEAN_CORPUS[i]; i++) strcat(all, LEAN_CORPUS[i]);
    ASSERT(eval_true_with_source("With[{r = LeanVerify[%s]}, r[\"Success\"] === True && "
                                 "r[\"Replay\"] === True && Length[Keys[r[\"Axioms\"]]] >= 34]",
                                 all));
    /* Refused constructs inserted into an accepted source: not run. */
    int nref = 0;
    while (LEAN_REFUSED[nref]) nref++;
    for (int j = 0; j < nref; j += 7) {
        size_t n = strlen(LEAN_CORPUS[14]), m = strlen(LEAN_REFUSED[j]);
        char* t = malloc(n + m + 1);
        ASSERT(t != NULL);
        memcpy(t, LEAN_CORPUS[14], n - 1);
        memcpy(t + n - 1, LEAN_REFUSED[j], m);
        strcpy(t + n - 1 + m, "\n");
        if (!eval_true_with_source("LeanVerify[%s][\"ExitCode\"] === Missing[\"NotRun\"]", t)) {
            Expr* r = eval_with_source("LeanVerify[%s]", t);
            char* p = r ? expr_to_string(r) : NULL;
            printf("FAIL: ran after inserting \"%s\": %s\n", LEAN_REFUSED[j], p ? p : "(null)");
            free(p);
            if (r) expr_free(r);
            ASSERT(0);
        }
        free(t);
    }
    free(all);
    /* A false theorem in the accepted surface is refuted by Lean. */
    check_true("leanRan[LeanVerify[\"theorem f1 : 2 + 2 = 5 := by decide\\n\"]][\"Success\"] === "
               "False");
    check_true("leanRan[LeanVerify[\"theorem f2 (n : Nat) : n + 1 = n := by omega\\n\"]]"
               "[\"Success\"] === False");
    /* Findings 4-5 end to end: refused, so nothing reaches Lean. */
    check_true("ProofObject[Evaluate[Simplify[(x^2 - 1)/(x - 1)]] == x + 1] === $Failed");
    check_true("leanV4 = 5; With[{r = ProofObject[leanV4^2 == 25]}, Clear[leanV4]; r] === $Failed");
    /* A Method script inside the surface still verifies. */
    check_true("Certify[x + 0 == x, Method -> \"intro x\\ngrind\"][\"Verified\"] === True");
}

/* Sixth review round (v0.6522), judged by Lean. */
static void test_review6_verify(void) {
    test_mute_messages();
    /* Finding 1: Domain -> Complexes is stated over every char-0 field, so the
     * vacuous-over-Q claim (x = 2 I is a counterexample) is refuted. */
    check_true("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Complexes, Method -> "
               "\"intro x h\\nhave h2 : 0 \\[LessEqual] x ^ 2 := Mathilda.sqnn\\ngrind\"]"
               "[\"Verified\"] === False");
    check_true("Certify[x == 7, Assumptions -> x^2 + 4 == 0, Domain -> Complexes][\"Verified\"] "
               "=== False");
    /* ... and true generic statements still verify, for each domain. */
    check_true("Certify[(x + 1)^2 == x^2 + 2 x + 1, Domain -> Complexes][\"Verified\"] === True");
    check_true("Certify[(x + 1)^2 == x^2 + 2 x + 1, Domain -> Reals][\"Verified\"] === True");
    check_true("Certify[x^2 == 1, Assumptions -> x == 1, Domain -> Reals][\"Verified\"] === True");
    check_true("Certify[x/2 + x/2 == x, Domain -> Complexes][\"Verified\"] === True");
    check_true("Certify[x == 1, Assumptions -> x^3 == 1 && Element[x, Reals], "
               "Domain -> Complexes] =!= $Failed");
    check_true("Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Algebraics] === $Failed");
    /* Finding 2: the claim is certified as D[F, x] == f, and the theorem says so. */
    assert_eval_eq("leanI6 = Quiet[Certify[Integrate[2 x, x] == x^2 + 5]]; Null", "Null", 0);
    check_true("leanI6[\"Verified\"] === True && "
               "leanI6[\"Theorem\"] === HoldForm[ForAll[x, D[x^2 + 5, x] == 2 x]] && "
               "AnyTrue[leanI6[\"TrustedSteps\"], StringContainsQ[#, \"locally constant\"] &]");
    /* Findings 3, 4, 7: refused before Lean runs. */
    check_true("leanE6 = Simplify[(x^2 - 1)/(x - 1)]; Protect[leanE6]; "
               "With[{r = Certify[leanE6 == x + 1]}, Unprotect[leanE6]; Clear[leanE6]; r] "
               "=== $Failed");
    check_true("leanT6 /: leanT6^2 = -1; With[{r = Certify[leanT6^2 == -1]}, Clear[leanT6]; r] "
               "=== $Failed");
    check_true("Certify[Integrate[2, x] == Sqrt[2] Lean`Private`s1 x] === $Failed");
    check_true("Certify[x == 1, Assumptions -> Reduce[x^2 == 1 && x > 0, x]] === $Failed");
    /* Finding 5: a genuine verdict is authenticated, a minted one is not. */
    check_true("With[{p = Certify[x + 1 == 1 + x]}, p[\"Verified\"] === True && "
               "With[{a = First[p]}, ProofObject[a]][\"Verified\"] === True]");
    check_true("ProofObject[<|\"Verified\" -> True|>][\"Verified\"] === "
               "Missing[\"Unauthenticated\"]");
}

/* Engine -> {"Lean", "Exact"} (v0.6525): True only when the kernel accepted;
 * on a false claim both engines say False and they do not disagree. */
static void test_multi_engine_verify(void) {
    test_mute_messages();
    assert_eval_eq("leanME = Certify[(x+1)^3 == x^3+3x^2+3x+1, Engine -> {\"Lean\", \"Exact\"}]; "
                   "Null", "Null", 0);
    check_true("leanME[\"Verified\"] === True");
    check_true("leanME[\"Engines\"] === <|\"Lean\" -> True, \"Exact\" -> Missing[\"Exact\", True]|>");
    check_true("StringContainsQ[ToString[leanME], \"Verified: True | Lean Tier 0 + Exact |\"]");
    check_true("Check[leanMF = Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, "
               "Engine -> {\"Lean\", \"Exact\"}]; \"none\", \"disagree\", Certify::disagree] === "
               "\"none\"");
    check_true("leanMF[\"Verified\"] === False");
    check_true("leanMF[\"Engines\"] === <|\"Lean\" -> False, \"Exact\" -> Missing[\"Exact\", False]|>");
}

/* Counterexamples (v0.6526): a kernel rejection triggers a bounded search;
 * a rational witness gives a kernel-checked refutation not (S). */
static void test_refute_verify(void) {
    test_mute_messages();
    assert_eval_eq("leanCX = Check[Certify[(x+1)^3 == x^3+3x^2+2x+1], \"refuted\", "
                   "Certify::refuted]; Null", "Null", 0);
    check_true("leanCX === \"refuted\"");
    assert_eval_eq("leanCX = Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1]]; Null", "Null", 0);
    check_true("leanCX[\"Verified\"] === False && leanCX[\"Refuted\"] === True");
    check_true("leanCX[\"Counterexample\"] === {x -> -1}");
    check_true("leanCX[\"RefutationTheorem\"] === HoldForm[!ForAll[{x}, (1 + x)^3 == "
               "1 + 2 x + 3 x^2 + x^3]]");
    check_true("leanCX[\"Properties\"] === Sort[{\"Axioms\", \"Cached\", \"Certificate\", "
               "\"Counterexample\", \"Domain\", \"Engine\", \"Hypotheses\", \"Input\", \"Kind\", "
               "\"LeanSource\", \"LeanTheorem\", \"LeanVersion\", \"Logs\", \"ProofScript\", "
               "\"RefutationSource\", \"RefutationTheorem\", \"Refuted\", \"Result\", \"Theorem\", "
               "\"Tier\", \"Timing\", \"TrustedSteps\", \"Verified\"}]");
    check_true("StringStartsQ[ToString[leanCX], \"ProofObject[\\[LeftGuillemet]Identity"
               "\\[RightGuillemet] Refuted: {x -> -1} | Lean Tier 0 |\"]");
    /* The refutation source is accepted by LeanVerify on its own. */
    check_true("LeanVerify[leanNoFooter[leanCX[\"RefutationSource\"]]][\"Success\"] === True");
    check_true("Length[Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1, \"Counterexample\" -> False]]"
               "[\"Properties\"]] === 19");
    /* Showcase row 17's claim. */
    assert_eval_eq("leanC17 = Quiet[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x]]; Null",
                   "Null", 0);
    check_true("leanC17[\"Refuted\"] === True && leanC17[\"Counterexample\"] === {x -> -1}");
    /* Refute with Lean: each form kernel-checked, standard axioms. */
    assert_eval_eq("leanRL = Quiet[{Refute[(x+1)^3 == x^3+3x^2+2x+1], "
                   "Refute[ForAll[{x, y}, x^2 + y^2 >= 2 x y + 1/10]], "
                   "Refute[x^2 == x, Assumptions -> x^2 == 2 x], "
                   "Refute[Implies[x^2 + y^2 == 2 && x == y, x == 1]], "
                   "Refute[x^3 + y^3 == (x + y)^3, \"Domain\" -> Integers], "
                   "Refute[1/x == 1], "
                   "Refute[Sum[k^3, {k, 1, n}] == n^2 (n+1)^2/4 + 1]}]; Null", "Null", 0);
    check_true("AllTrue[leanRL, #[\"Verified\"] === True &]");
    check_true("AllTrue[leanRL, SubsetQ[{\"propext\", \"Classical.choice\", \"Quot.sound\"}, "
               "#[\"Axioms\"]] &]");
    check_true("StringContainsQ[ToString[First[leanRL]], \"Refutation\\[RightGuillemet] "
               "Verified: True | Lean Tier 0 |\"]");
    /* An algebraic witness has no kernel proof: Exact-only. */
    check_true("Quiet[Refute[Implies[x^2 + 4 == 0, x == 7]]][\"Verified\"] === "
               "Missing[\"Exact\", True]");
    /* A closed fact (Solve's spurious root zeroes a denominator). */
    assert_eval_eq("leanCS = Quiet[Certify[Solve[(x - 1)/(x^2 - 1) == (x - 1)/(2 x - 2), x]]]; "
                   "Null", "Null", 0);
    check_true("leanCS[\"Verified\"] === False && leanCS[\"Refuted\"] === True && "
               "leanCS[\"Counterexample\"] === {}");
    /* A wrong antiderivative claim: Exact-only refutation, no Lean source. */
    assert_eval_eq("leanCA = Quiet[Certify[Integrate[2 x, x] == x^2 + x]]; Null", "Null", 0);
    check_true("leanCA[\"Refuted\"] === Missing[\"Exact\", True] && "
               "!KeyExistsQ[First[leanCA], \"RefutationSource\"]");
    /* With the Exact engine listed. */
    check_true("Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1, Engine -> {\"Lean\", \"Exact\"}]]"
               "[\"Refuted\"] === True");
    /* The search honours an Assumptions the sum statement omits: n -> 2
     * would refute the statement but not the claim for n > 3. */
    assert_eval_eq("leanCN = Quiet[Certify[Sum[k, {k, 1, n}] == n^2, Assumptions -> n > 3]]; "
                   "Null", "Null", 0);
    check_true("leanCN[\"Verified\"] === False && leanCN[\"Counterexample\"] === {n -> 4}");
    /* Failure injection: a non-witness is refused, so nothing is attached. */
    assert_eval_eq("Lean`Private`$witnessOverride = {x -> 0}; "
                   "leanCW = Quiet[Certify[(x+1)^3 == x^3+3x^2+2x+1]]; "
                   "Lean`Private`$witnessOverride = None; Null", "Null", 0);
    check_true("leanCW[\"Verified\"] === False && !KeyExistsQ[First[leanCW], \"Refuted\"] && "
               "!KeyExistsQ[First[leanCW], \"Counterexample\"]");
}

/* Conjecture (v0.6527): guessed, tested and proved by induction in Lean. */
static void test_conjecture_verify(void) {
    test_mute_messages();
    assert_eval_eq("leanJ = {Conjecture[Sum[k^4, {k, 1, n}], n], "
                   "Conjecture[Sum[Fibonacci[k]^2, {k, 1, n}], n]}; Null", "Null", 0);
    check_true("AllTrue[leanJ, #[\"Status\"] === \"Proved\" && #[\"Verified\"] === True &]");
    check_true("AllTrue[leanJ, #[\"Proof\"][\"Verified\"] === True && "
               "#[\"Proof\"][\"Engines\"][\"Lean\"] === True &]");
    check_true("leanJ[[2]][\"Form\"] === \"Fibonacci\" && "
               "SubsetQ[{\"propext\", \"Classical.choice\", \"Quot.sound\"}, "
               "leanJ[[2]][\"Proof\"][\"Axioms\"]]");
    check_true("StringStartsQ[ToString[leanJ[[1]]], \"ProofObject[\\[LeftGuillemet]Conjecture"
               "\\[RightGuillemet] Proved: Sum[k^4, {k, 1, n}] == \"] && "
               "StringContainsQ[ToString[leanJ[[1]]], \"| Lean Tier 0 + Exact |\"]");
    check_true("Certify[Sum[Fibonacci[k], {k, 0, n}] == Fibonacci[n + 2] - 1][\"Verified\"] === True");
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
    /* Certify sources end with the trusted footer's #print axioms lines,
     * which only Certify's own check accepts: tamper tests give LeanVerify the
     * source without them, and require that Lean actually ran. */
    assert_eval_eq("leanNoFooter[s_String] := StringReplace[s, "
                   "RegularExpression[\"#print axioms [^\\n]*\\n\"] -> \"\"]; "
                   "leanRan[r_] := If[r[\"ExitCode\"] === Missing[\"NotRun\"], "
                   "<|\"Success\" -> \"not run\"|>, r]; Null", "Null", 0);
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
    TEST(test_review3_allowlist);
    TEST(test_review3_domains);
    TEST(test_review4_verify);
    TEST(test_review6_verify);
    TEST(test_multi_engine_verify);
    TEST(test_refute_verify);
    TEST(test_conjecture_verify);
    printf("All Lean verification tests passed.\n");
    return 0;
}
