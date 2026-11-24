/* Builtin-audit regressions: pattern matching (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* PatternSequence#1 (high): "Condition on a PatternSequence never matches". */
static void test_condition_on_patternsequence(void) {
    audit_run("ClearAll[psh, x, y, s]");
    audit_chk("PatternSequence#1", "MatchQ[f[1, 2], f[PatternSequence[x_, y_] /; x < y]]", "True");
    audit_chk("PatternSequence#1", "MatchQ[f[2, 1], f[PatternSequence[x_, y_] /; x < y]]", "False");
    audit_run("psh[s : PatternSequence[x_, y_] /; x < y] := {s}");
    audit_chk("PatternSequence#1", "{psh[1, 2], psh[2, 1]}", "{{1, 2}, psh[2, 1]}");
    audit_chk("PatternSequence#1", "f[1, 2, 3] /. f[PatternSequence[x_, y_] /; x > 1, z_] :> {x, y, z}", "f[1, 2, 3]");
    audit_chk("PatternSequence#1", "f[3, 2, 1] /. f[PatternSequence[x_, y_] /; x > 1, z_] :> {x, y, z}", "{3, 2, 1}");
    /* sanity */
    audit_chk("PatternSequence#1", "MatchQ[f[1, 2], f[PatternSequence[x_, y_]]]", "True");
    audit_chk("PatternSequence#1", "MatchQ[f[1, 2], f[x__ /; Length[{x}] == 2]]", "True");
}

/* S14 (Wave 3): one level-spec reading for the pattern heads. */
static void test_s14_level_specs(void) {
    audit_run("ClearAll[a, b, c, x]");
    /* Count#1 / DeleteCases#1: -n is levels 1 through -n, as Level reads it */
    audit_chk("Count#1", "Count[{a, {b}}, _, -1]", "3");
    audit_true("Count#1", "Count[{a, {b}}, _, -1] == Length[Level[{a, {b}}, -1]]");
    audit_chk("Count#1", "Count[{a, {b}}, _, {-1}]", "2");
    audit_chk("Count#1", "Count[5, _?NumberQ, -1]", "0");
    audit_true("DeleteCases#1", "DeleteCases[{a, {b}}, _List, -1] === {a}");
    audit_true("DeleteCases#1", "DeleteCases[{a, {b}}, _List, {-1}] === {a, {b}}");
    /* Count#2: a BigInt level is Infinity */
    audit_chk("Count#2", "Count[{1, {2}}, _, 2^70]", "3");
    audit_chk("Count#2", "Count[{1, {2}}, _, {1, 2^70}]", "3");
    /* Cases#2: ::level for a bad spec, ::innf for a bad count */
    audit_chk("Cases#2", "Check[Cases[{1}, _, x], fired, Cases::level]", "fired");
    audit_chk("Cases#2", "Check[Cases[{1, 2}, _, {1, 2, 3}], fired, Cases::level]", "fired");
    audit_chk("Cases#2", "Check[Cases[{1, 2}, _, 1, -1], fired, Cases::innf]", "fired");
    audit_true("Cases#2", "Quiet[Head[Cases[{1, 2}, _, 1, -1]]] === Cases");
    audit_true("Cases#2", "Cases[{1, 2, 3}, _, 1, Infinity] === {1, 2, 3}");
    audit_true("Cases#2", "Cases[{1, 2, 3}, _, 1, 2] === {1, 2}");
    /* FirstPosition#2: the message is FirstPosition's, the call unevaluated */
    audit_chk("FirstPosition#2", "Check[FirstPosition[{1, 2}, 1, None, x], fired, FirstPosition::level]",
              "fired");
    audit_true("FirstPosition#2", "Quiet[Head[FirstPosition[{1, 2}, 1, None, x]]] === FirstPosition");
    audit_true("FirstPosition#2", "FirstPosition[{1, {2}}, 2, None, {2}] === {2, 1}");
    /* Optional#1: a whole-subject Optional matches what its pattern matches */
    audit_chk("Optional#1", "MatchQ[2, x_:0]", "True");
    audit_chk("Optional#1", "MatchQ[2, x_.]", "True");
    audit_true("Optional#1", "Cases[{1, a, 2}, _Integer:0] === {1, 2}");
    audit_chk("Optional#1", "Replace[2, x_:0 :> x + 1]", "3");
    audit_chk("Optional#1", "MatchQ[a, _Integer:0]", "False");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Count#3", "Quiet[Check[Count[{1, 2}, 1, \"a\"], fired, Count::level]]", "fired");
    audit_chk("Count#3", "Quiet[Check[Count[{1, 2}, 1, {1}, 5], fired, Count::nonopt]]", "fired");
    audit_chk("Count#3", "Quiet[Check[Count[], fired, Count::argb]]", "fired");
    audit_chk("Count#3", "Count[{a -> b, c}, a -> b]", "1");
    audit_chk("Count#3", "Count[a -> b][{a -> b}]", "1");
    audit_chk("MatchQ#3", "Quiet[Check[MatchQ[1, 2, 3], fired, MatchQ::argt]]", "fired");
    audit_chk("MatchQ#3", "Quiet[MatchQ[1, 2, 3]]", "MatchQ[1, 2, 3]");
    audit_chk("FirstPosition#1", "Quiet[Check[FirstPosition[{a}], fired, FirstPosition::argbu]]", "fired");
    audit_chk("FirstPosition#1", "Quiet[FirstPosition[{a}]]", "FirstPosition[{a}]");
    audit_chk("ReplaceList#3", "Quiet[Check[ReplaceList[a, b, c, d], fired, ReplaceList::argt]]", "fired");
    audit_chk("ReplaceList#3", "Quiet[ReplaceList[a, b, c, d]]", "ReplaceList[a, b, c, d]");
}

int main(void) {
    audit_init();
    TEST(test_condition_on_patternsequence);
    TEST(test_s14_level_specs);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_pattern_matching");
}
