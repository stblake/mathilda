/* Rational and Complex numbers are ATOMS.
 *
 * Both are stored as two-argument EXPR_FUNCTION nodes (Rational[3, 4],
 * Complex[2, 3]), so a structural builtin that only tests the node type
 * reaches inside them: Length[3/4] was 2, Apply[g, 3/4] was g[3, 4],
 * Append[3/4, 1] was Rational[3, 4, 1]. In the reference every one of these
 * treats the number as an atom: no parts, ::normal for the builtins that need
 * an expression with parts. Every expected output below was produced by the
 * reference (wolframscript).
 *
 * Pattern matching is NOT part of this: MatchQ[1/2, h_[___]] is True in the
 * reference as well (Rational[a_, b_] matches 1/2), and the matcher keeps that
 * -- pinned below so it is not "fixed" by mistake. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

static void test_atomq_and_matching_unchanged(void) {
    assert_eval_eq("AtomQ[3/4]", "True", 0);
    assert_eval_eq("AtomQ[2 + 3 I]", "True", 0);
    assert_eval_eq("MatchQ[1/2, h_[___]]", "True", 0);
    assert_eval_eq("MatchQ[1/2, Rational[a_, b_]]", "True", 0);
    assert_eval_eq("MatchQ[1 + 2 I, Complex[a_, b_]]", "True", 0);
    assert_eval_eq("Cases[{1/2, f[x]}, h_[___] :> h]", "{Rational, f}", 0);
    assert_eval_eq("1/2 /. Rational[a_, b_] :> {a, b}", "{1, 2}", 0);
}

static void test_length_dimensions(void) {
    assert_eval_eq("Length[3/4]", "0", 0);
    assert_eval_eq("Length[2 + 3 I]", "0", 0);
    assert_eval_eq("Length[{3/4, 1}]", "2", 0);
    assert_eval_eq("Dimensions[3/4]", "{}", 0);
    assert_eval_eq("Dimensions[2 + 3 I]", "{}", 0);
    assert_eval_eq("Dimensions[{3/4, 1}]", "{2}", 0);
    assert_eval_eq("Dimensions[{{3/4}, {1}}]", "{2, 1}", 0);
    assert_eval_eq("ArrayDepth[3/4]", "0", 0);
    assert_eval_eq("ArrayDepth[2 + 3 I]", "0", 0);
}

static void test_apply_replacepart(void) {
    assert_eval_eq("g @@ (3/4)", "3/4", 0);
    assert_eval_eq("g @@ (2 + 3 I)", "2 + 3*I", 0);
    assert_eval_eq("Apply[g, 3/4, {0}]", "3/4", 0);
    assert_eval_eq("g @@@ {3/4}", "{3/4}", 0);
    assert_eval_eq("Plus @@ (3/4)", "3/4", 0);
    assert_eval_eq("g @@ {3/4}", "g[3/4]", 0);
    assert_eval_eq("ReplacePart[3/4, 1 -> 5]", "3/4", 0);
    assert_eval_eq("ReplacePart[3/4, {1} -> 5]", "3/4", 0);
    assert_eval_eq("ReplacePart[3/4, {} -> 5]", "3/4", 0);
    assert_eval_eq("ReplacePart[2 + 3 I, 1 -> 5]", "2 + 3*I", 0);
    assert_eval_eq("ReplacePart[{3/4, 1}, 1 -> 5]", "{5, 1}", 0);
}

/* ::normal and the call left unevaluated */
static void test_normal_message(void) {
    assert_eval_eq("Quiet[Take[3/4, 1]]", "Take[3/4, 1]", 0);
    assert_eval_eq("Quiet[Drop[3/4, 1]]", "Drop[3/4, 1]", 0);
    assert_eval_eq("Quiet[Append[3/4, 1]]", "Append[3/4, 1]", 0);
    assert_eval_eq("Quiet[Prepend[3/4, 1]]", "Prepend[3/4, 1]", 0);
    assert_eval_eq("Quiet[Sort[3/4]]", "Sort[3/4]", 0);
    assert_eval_eq("Quiet[SortBy[3/4, Minus]]", "SortBy[3/4, Minus]", 0);
    assert_eval_eq("Quiet[Reverse[3/4]]", "Reverse[3/4]", 0);
    assert_eval_eq("Quiet[Flatten[3/4]]", "Flatten[3/4]", 0);
    assert_eval_eq("Quiet[Select[3/4, True &]]", "Select[3/4, True &]", 0);
    assert_eval_eq("Quiet[Take[2 + 3 I, 1]]", "Take[2 + 3*I, 1]", 0);
    assert_eval_eq("Quiet[Append[2 + 3 I, 1]]", "Append[2 + 3*I, 1]", 0);
    assert_eval_eq("Quiet[Sort[2 + 3 I]]", "Sort[2 + 3*I]", 0);
    /* the message fires (Check sees it) and names the builtin */
    assert_eval_eq("Check[Take[3/4, 1], msg, Take::normal]", "msg", 0);
    assert_eval_eq("Check[Append[3/4, 1], msg, Append::normal]", "msg", 0);
    assert_eval_eq("Check[Reverse[3/4], msg, Reverse::normal]", "msg", 0);
    assert_eval_eq("Check[Flatten[2 + 3 I], msg, Flatten::normal]", "msg", 0);
    assert_eval_eq("Check[Select[3/4, True &], msg, Select::normal]", "msg", 0);
    /* lists of numbers are untouched */
    assert_eval_eq("Take[{3/4, 1}, 1]", "{3/4}", 0);
    assert_eval_eq("Reverse[{3/4, 2 + 3 I}]", "{2 + 3*I, 3/4}", 0);
    assert_eval_eq("Flatten[{{3/4}, {1}}]", "{3/4, 1}", 0);
    assert_eval_eq("Sort[{3/4, 1/2}]", "{1/2, 3/4}", 0);
    assert_eval_eq("Append[{3/4}, 1]", "{3/4, 1}", 0);
}

static void test_delete_insert(void) {
    assert_eval_eq("Quiet[Delete[3/4, 1]]", "Delete[3/4, 1]", 0);
    assert_eval_eq("Quiet[Insert[3/4, 1, 1]]", "Insert[3/4, 1, 1]", 0);
    assert_eval_eq("Quiet[Delete[2 + 3 I, 1]]", "Delete[2 + 3*I, 1]", 0);
    assert_eval_eq("Check[Delete[3/4, 1], msg, Delete::partw]", "msg", 0);
    assert_eval_eq("Check[Insert[3/4, 1, 1], msg, Insert::ins]", "msg", 0);
    assert_eval_eq("Delete[3/4, {}]", "3/4", 0);
    assert_eval_eq("Insert[3/4, 1, {}]", "3/4", 0);
    assert_eval_eq("Delete[{3/4, 1}, 1]", "{1}", 0);
}

int main() {
    symtab_init();
    core_init();
    TEST(test_atomq_and_matching_unchanged);
    TEST(test_length_dimensions);
    TEST(test_apply_replacepart);
    TEST(test_normal_message);
    TEST(test_delete_insert);
    printf("All number atomicity tests passed.\n");
    return 0;
}
