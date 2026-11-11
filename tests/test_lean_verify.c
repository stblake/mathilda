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
    printf("All Lean verification tests passed.\n");
    return 0;
}
