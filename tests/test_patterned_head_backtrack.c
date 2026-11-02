/* test_patterned_head_backtrack.c
 *
 * A call whose head carries patterns of its own (a SubValue LHS such as
 * FeynCalc's R[nN_, k_][i___Integer, mm:(_Integer)..][a_List] := ... /; cond)
 * matched its head with no continuation, so the head committed to its first
 * sequence split: when the guard failed for that split, the match failed
 * instead of trying i = {0, 0}, mm = {1, 1}.  PaVeReduce left R[2, 1][0, 0, 1,
 * 1][{0, m1^2, m2^2}] unevaluated (20 tests).  match.c now passes the
 * argument level (and the guard) as the head match's continuation.
 * Expected values computed by the reference kernel.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_subvalue_guard_backtracks(void) {
    assert_eval_eq("g1[i___Integer, mm:(_Integer)..] := {{i}, {mm}} /; {mm}[[1]] === 1; g1[0, 0, 1, 1]",
                   "{{0, 0}, {1, 1}}", 0);
    assert_eval_eq("h1[n_][i___Integer, mm:(_Integer)..][a_List] := {{i}, {mm}} /; {mm}[[1]] === 1;"
                   " h1[2][0, 0, 1, 1][{x}]", "{{0, 0}, {1, 1}}", 0);
    assert_eval_eq("q1[a_] := 1; k1[n_][i___Integer, mm:(_Integer)..][a_List] :="
                   " Block[{q}, q = 1; {{i}, {mm}}] /; ({mm}[[1]] === q1[a]); k1[2][0, 0, 1, 1][{x}]",
                   "{{0, 0}, {1, 1}}", 0);
    /* the guard fails for every split: unevaluated */
    assert_eval_eq("h2[n_][i___Integer, mm:(_Integer)..][a_List] := {i} /; {mm}[[1]] === 7;"
                   " h2[2][0, 0, 1, 1][{x}]", "h2[2][0, 0, 1, 1][{x}]", 0);
}

static void test_patterned_head_in_matchq(void) {
    assert_eval_eq("MatchQ[f[1][a, b, c], f[n_][x__, y__] /; Length[{x}] == 2]", "True", 0);
    assert_eval_eq("Cases[{f[1][a, b, c]}, f[n_][x__, y__] /; Length[{y}] == 2 :> {x}]", "{{a}}", 0);
    assert_eval_eq("MatchQ[f[a, b][c], f[x__, y_][z_] /; Length[{x}] == 1]", "True", 0);
}

static void test_holdpattern_sequence_width(void) {
    /* HoldPattern is transparent to an argument's width (v0.6137), as in the
     * reference: FeynCalc's FCMellinJoin checks {HoldPattern[(...)..]}. */
    assert_eval_eq("MatchQ[{{1}, {1}}, {HoldPattern[({1} | {-1}) ..]}]", "True", 0);
    assert_eval_eq("MatchQ[{1, 1}, {HoldPattern[1 ..]}]", "True", 0);
    assert_eval_eq("MatchQ[{1, 2}, {HoldPattern[__]}]", "True", 0);
    assert_eval_eq("MatchQ[{1, 2}, {HoldPattern[_]}]", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_subvalue_guard_backtracks);
    TEST(test_patterned_head_in_matchq);
    TEST(test_holdpattern_sequence_width);
    printf("All patterned-head backtracking tests passed.\n");
    return 0;
}
