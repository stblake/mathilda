/*
 * test_lean_form.c -- LeanForm display text and the strict Lean`Internal`ToLean
 * translator (src/print_lean.c).  Needs no Lean installation.
 *
 * Every expected string is the Lean 4 syntax the printer is specified to
 * produce (print_lean.h); the strict-form cases pin the refusals that keep a
 * certificate's statement faithful (inexact numbers, division over Int,
 * Mathlib-only heads under core_only, unknown heads).
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

/* ToString[e, LeanForm] must be expected. */
static void lf(const char* e, const char* expected) {
    char buf[1024];
    snprintf(buf, sizeof buf, "ToString[%s, LeanForm]", e);
    check_str(buf, expected);
}

static void test_plan_table(void) {
    lf("(x+y)^2 == x^2+2x y+y^2", "(x + y) ^ 2 = x ^ 2 + 2 * x * y + y ^ 2");
    lf("-x^2", "-(x ^ 2)");
    lf("x - 3", "x - 3");
    lf("3/4 x", "(3 / 4) * x");
    lf("ForAll[{x,y}, x^2+y^2 >= 2 x y]", "\xE2\x88\x80 x y, x ^ 2 + y ^ 2 \xE2\x89\xA5 2 * x * y");
    lf("Function[{x}, x+1]", "fun x => x + 1");
    lf("{1, 2}", "[1, 2]");
    lf("Log[x]", "Real.log x");
    lf("f[x, y]", "f x y");
    lf("a`b", "\xC2\xAB" "a`b" "\xC2\xBB");
}

static void test_arithmetic(void) {
    lf("-x y", "-(x * y)");
    lf("-2 x", "-(2 * x)");
    lf("-x", "-x");
    lf("a - 2 b", "a - 2 * b");
    lf("x/(y z)", "x / (y * z)");
    lf("x^-1", "1 / x");
    lf("1/x^2 + 1", "1 + 1 / x ^ 2");
    lf("Sqrt[2] x", "Real.sqrt 2 * x");
    lf("-1/2", "-(1 / 2)");
    lf("2^70", "1180591620717411303424");
    lf("(-2)^x", "(-2) ^ x");
    lf("HoldForm[(x y)^2]", "(x * y) ^ 2");
    lf("x^y^z", "x ^ y ^ z");
    lf("f[g[x], -3, 1/2]", "f (g x) (-3) (1 / 2)");
    lf("1.5", "1.5");
    lf("\"s\\\"q\"", "\"s\\\"q\"");
}

static void test_logic_and_symbols(void) {
    lf("a < b <= c", "a < b \xE2\x88\xA7 b \xE2\x89\xA4 c");
    lf("Implies[a && b, c || !d]", "a \xE2\x88\xA7 b \xE2\x86\x92 c \xE2\x88\xA8 \xC2\xAC" "d");
    lf("Equivalent[a, b]", "a \xE2\x86\x94 b");
    lf("a != b", "a \xE2\x89\xA0 b");
    lf("(a || b) && c", "(a \xE2\x88\xA8 b) \xE2\x88\xA7 c");
    lf("Exists[x, x > 0, x^2 == 2]", "\xE2\x88\x83 x, x > 0 \xE2\x88\xA7 x ^ 2 = 2");
    lf("fun + $x", "\xC2\xAB" "fun" "\xC2\xBB + \xC2\xAB$x\xC2\xBB");
    lf("{True, False}", "[True, False]");
    lf("Pi + E", "Real.exp 1 + Real.pi");
    /* A binder as an operand is parenthesised. */
    lf("p && ForAll[x, x == x]", "p \xE2\x88\xA7 (\xE2\x88\x80 x, x = x)");
    check_str("ToString[LeanForm[x^2]]", "x ^ 2");
}

static void test_format_rules(void) {
    assert_eval_eq("Format[leanFmtF[x_], LeanForm] := leanFmtG[x, 1]; Null", "Null", 0);
    lf("leanFmtF[y] + 1", "1 + leanFmtG y 1");
    /* The strict translator ignores display rules. */
    check_str("Lean`Internal`ToLean[leanFmtF[y], None, True, {\"leanFmtF\"}]", "leanFmtF y");
}

static void test_strict(void) {
    check_str("Lean`Internal`ToLean[HoldForm[ForAll[{x}, x/2 + x/2 == x]], \"Rat\", True, {}]",
              "\xE2\x88\x80 x : Rat, x / 2 + x / 2 = x");
    check_str("Lean`Internal`ToLean[HoldForm[3/4 x == y], \"Rat\", True, {}]",
              "(3 / 4 : Rat) * x = y");
    /* A closed relation ascribes its left side so numerals are not Nat. */
    check_str("Lean`Internal`ToLean[HoldForm[2 + 2 == 4], \"Int\", True, {}]", "(2 + 2 : Int) = 4");
    check_str("Lean`Internal`ToLean[HoldForm[S[n] == n^2], None, True, {\"S\"}]", "S n = n ^ 2");
    check_str("Lean`Internal`ToLean[Log[x], \"Real\", False, {}]", "Real.log x");
    /* Refusals give $Failed and LeanForm::untrans. */
    test_mute_messages();
    assert_eval_eq("Lean`Internal`ToLean[Log[x], \"Rat\", True, {}]", "$Failed", 0);
    assert_eval_eq("Check[Lean`Internal`ToLean[Log[x], \"Rat\", True, {}], \"untrans\"]",
                   "\"untrans\"", 0);
    assert_eval_eq("Lean`Internal`ToLean[1.5, \"Rat\", True, {}]", "$Failed", 0);
    assert_eval_eq("Lean`Internal`ToLean[x/2, \"Int\", True, {}]", "$Failed", 0);
    assert_eval_eq("Lean`Internal`ToLean[x^(1/3), \"Rat\", True, {}]", "$Failed", 0);
    assert_eval_eq("Lean`Internal`ToLean[Pi, \"Real\", True, {}]", "$Failed", 0);
    assert_eval_eq("Lean`Internal`ToLean[Infinity, \"Real\", False, {}]", "$Failed", 0);
    assert_eval_eq("Lean`Internal`ToLean[f[x], \"Rat\", True, {}]", "$Failed", 0);
    assert_eval_eq("Lean`Internal`ToLean[#^2 &, \"Rat\", True, {}]", "$Failed", 0);
    /* Malformed calls stay unevaluated. */
    assert_eval_eq("Head[Lean`Internal`ToLean[x, \"Complex\", True, {}]] === Lean`Internal`ToLean", "True", 0);
}

static void test_attributes(void) {
    assert_eval_eq("Attributes[LeanForm]", "{Protected}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_plan_table);
    TEST(test_arithmetic);
    TEST(test_logic_and_symbols);
    TEST(test_format_rules);
    TEST(test_strict);
    TEST(test_attributes);
    printf("All LeanForm tests passed.\n");
    return 0;
}
