/* Unit tests for Goto / Label (Wolfram-Language imperative control flow).
 *
 * Covers forward and backward jumps within a CompoundExpression, the canonical
 * Newton-iteration loop from the spec, jumps whose Label is the last statement,
 * a Goto fired inside a nested call (If branch) reaching the enclosing
 * CompoundExpression, inner->outer propagation when the Label lives in an
 * enclosing CompoundExpression, and the inert cases (bare Label, unmatched
 * Goto). See src/eval.c (eval_is_inflight_goto + arg-loop short-circuit +
 * eval_report_uncaught_goto), src/core.c (builtin_compoundexpression
 * backward-jump loop) and src/funcprog.c (builtin_goto/builtin_label).
 *
 * NOTE: the unmatched-Goto cases deliberately emit "Goto::nolabel" on stderr;
 * that is expected output, not a failure. */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ----- Forward jump: statements between Goto and its Label are skipped ----- */
void test_forward_jump() {
    /* a=2 is skipped; result is the value set before the jump. */
    assert_eval_eq("Module[{a}, a = 1; Goto[skip]; a = 2; Label[skip]; a]", "1", 0);
    /* Explicit CompoundExpression head, integer tag. */
    assert_eval_eq("Module[{a}, a = 10; Goto[2]; a = 20; Label[2]; a]", "10", 0);
}

/* ----- Backward jump forms a loop that terminates on a condition ----- */
void test_backward_loop() {
    /* Sum 1..5 via a Label/Goto loop: 1+2+3+4+5 = 15. */
    assert_eval_eq(
        "Module[{i = 0, s = 0}, Label[top]; i = i + 1; s = s + i; "
        "If[i < 5, Goto[top]]; s]", "15", 0);
    /* Count-down loop landing exactly on the exit. */
    assert_eval_eq(
        "Module[{n = 3, c = 0}, Label[again]; c = c + 1; n = n - 1; "
        "If[n > 0, Goto[again]]; c]", "3", 0);
}

/* ----- Canonical Newton sqrt iteration from the spec ----- */
void test_newton_example() {
    assert_eval_eq(
        "f[a_] := Module[{x = 1., xp}, Label[begin]; "
        "If[Abs[xp - x] < 10^-8, Goto[end]]; xp = x; x = (x + a/x)/2; "
        "Goto[begin]; Label[end]; x]", "Null", 0);
    assert_eval_startswith("f[2]", "1.41421");
    assert_eval_eq("Clear[f]", "Null", 0);
}

/* ----- Jump target Label as the final statement yields Null ----- */
void test_label_last() {
    assert_eval_eq("Module[{a}, a = 1; Goto[done]; a = 2; Label[done]]", "Null", 0);
}

/* ----- Goto fired inside a nested call (If branch) reaches the enclosing CE ----- */
void test_nested_in_if() {
    /* If[True, Goto[skip]] surfaces the Goto as the If statement's value; the
     * CompoundExpression then jumps, so a stays 0 (a=1 is skipped). */
    assert_eval_eq(
        "Module[{a}, a = 0; If[True, Goto[skip]]; a = 1; Label[skip]; a]", "0", 0);
    /* Goto fired from a DownValue body still propagates out to the CE. */
    assert_eval_eq("jmp[] := Goto[out]", "Null", 0);
    assert_eval_eq(
        "Module[{a}, a = 5; jmp[]; a = 6; Label[out]; a]", "5", 0);
    assert_eval_eq("Clear[jmp]", "Null", 0);
}

/* ----- Inner CompoundExpression Goto whose Label is in the enclosing one ----- */
void test_inner_outer_propagation() {
    /* The inner (parenthesised) CompoundExpression has no Label[out], so the
     * sentinel propagates to the outer CE which owns the Label. a is set to 1
     * by the inner CE, a=2 is skipped. */
    assert_eval_eq(
        "Module[{a}, a = 0; (a = 1; Goto[out]); a = 2; Label[out]; a]", "1", 0);
}

/* ----- Inert cases (unmatched Goto emits Goto::nolabel on stderr) ----- */
void test_inert_cases() {
    assert_eval_eq("Label[foo]", "Null", 0);          /* bare Label is a no-op */
    assert_eval_eq("Goto[foo]", "Goto[foo]", 0);      /* no enclosing Label: inert */
    /* A CompoundExpression with no matching Label returns the Goto unevaluated;
     * reaching the top level uncaught emits Goto::nolabel (stderr). */
    assert_eval_eq("(aa; Goto[missing]; bb)", "Goto[missing]", 0);
    /* An unmatched Goto inside a Module also surfaces unchanged (+ message). */
    assert_eval_eq("Module[{x = 1}, Goto[gone]; x]", "Goto[gone]", 0);
    /* Wrong arity is not a Goto sentinel: left unevaluated, no message. */
    assert_eval_eq("Goto[a, b]", "Goto[a, b]", 0);
    assert_eval_eq("Label[a, b]", "Label[a, b]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_forward_jump);
    TEST(test_backward_loop);
    TEST(test_newton_example);
    TEST(test_label_last);
    TEST(test_nested_in_if);
    TEST(test_inner_outer_propagation);
    TEST(test_inert_cases);

    printf("All goto/label tests passed!\n");
    symtab_clear();
    return 0;
}
