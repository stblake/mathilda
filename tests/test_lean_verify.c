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
               "[\"Verified\"] === True");
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

int main(void) {
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
    printf("All Lean verification tests passed.\n");
    return 0;
}
