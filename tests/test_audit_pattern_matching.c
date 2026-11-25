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

/* Wave 4 G4: Attributes[s] = list, MatchQ operator form and Unevaluated, Longest
 * around PatternSequence, Position's count, Default for Plus/Times/Power, and
 * a Condition around OrderlessPatternSequence. */
static void test_wave4_g4_patterns(void) {
    audit_run("ClearAll[pat1, f, x, y, a, b, c]");
    audit_chk("Attributes#1", "Attributes[pat1] = {Listable, Orderless}; Attributes[pat1]", "{Listable, Orderless}");
    audit_chk("Attributes#1", "pat1[{1, 2}]", "{pat1[1], pat1[2]}");
    audit_chk("Attributes#1", "pat1[b, a]", "pat1[a, b]");
    audit_chk("MatchQ#1", "MatchQ[_Integer][5]", "True");
    audit_chk("MatchQ#1", "Select[{1, a, 2}, MatchQ[_Integer]]", "{1, 2}");
    audit_chk("MatchQ#2", "MatchQ[Unevaluated[1 + 1], _Plus]", "True");
    audit_chk("MatchQ#2", "MatchQ[Unevaluated[Sin[Pi]], _Sin]", "True");
    audit_chk("MatchQ#2", "MatchQ[1 + 1, _Plus]", "False");
    audit_chk("PatternSequence#2", "MatchQ[f[a, b], f[Longest[PatternSequence[a, b]]]]", "True");
    audit_chk("PatternSequence#2", "f[a, b, c] /. f[Longest[x:PatternSequence[__]], y___] :> {{x}, {y}}", "{{a, b, c}, {}}");
    audit_chk("PatternSequence#2", "f[a, b, c] /. f[x:PatternSequence[__], y___] :> {{x}, {y}}", "{{a}, {b, c}}");
    audit_chk("Position#2", "Check[Position[{a, a}, a, {1}, -1], fired, Position::innf]", "fired");
    audit_chk("Position#2", "Quiet[Position[{a, a}, a, {1}, x]]", "Position[{a, a}, a, {1}, x]");
    audit_chk("Position#2", "Quiet[Position[{a, b}, a, {1}, 2.5]]", "Position[{a, b}, a, {1}, 2.5]");
    audit_chk("Position#2", "{Position[{a, a}, a, {1}, 1], Position[{a, a}, a, {1}, Infinity]}", "{{{1}}, {{1}, {2}}}");
    audit_chk("Default#1", "{Default[Plus], Default[Times], Default[Power, 2]}", "{0, 1, 1}");
    audit_chk("Default#1", "{Default[Plus, 1], Default[Power]}", "{0, Default[Power]}");
    audit_chk("OrderlessPatternSequence#1", "MatchQ[f[1, 2], f[OrderlessPatternSequence[x_, y_] /; x < y]]", "True");
    audit_chk("OrderlessPatternSequence#1", "MatchQ[f[1, 2], f[OrderlessPatternSequence[x_, y_]] /; x < y]", "True");
    audit_chk("OrderlessPatternSequence#1", "f[1, 2] /. f[OrderlessPatternSequence[x_, y_] /; x > y] :> {x, y}", "{2, 1}");
    audit_chk("OrderlessPatternSequence#1", "MatchQ[f[1, 2], f[OrderlessPatternSequence[x_, y_] /; x > 5]]", "False");
}

/* Wave 5 G4: Attributes argument checks, Position operator form and arity, DefaultValues
 * and a Protected Default, Heads :> v, SparseArray FirstPosition, ReplaceList checks. */
static void test_wave5_g4_patterns(void) {
    audit_chk("Attributes#2", "Check[Attributes[1], fired, Attributes::ssle]", "fired");
    audit_chk("Attributes#2", "Check[Attributes[Hold[w5x]], fired, Attributes::ssle]", "fired");
    audit_chk("Attributes#2", "Quiet[Attributes[1]]", "Attributes[1]");
    audit_chk("Attributes#2", "Attributes[HoldPattern[Plus]]", "{Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}");
    audit_chk("Attributes#3", "Check[Attributes[], fired, Attributes::argx]", "fired");
    audit_chk("Attributes#3", "Check[Attributes[w5a, w5b], fired, Attributes::argx]", "fired");
    audit_chk("Attributes#3", "Quiet[Attributes[w5a, w5b]]", "Attributes[w5a, w5b]");
    audit_chk("Position#3", "Check[Position[{a}, a, 1, 1, 1], fired, Position::nonopt]", "fired");
    audit_chk("Position#3", "Quiet[Position[{a}, a, 1, 1, 1]]", "Position[{a}, a, 1, 1, 1]");
    audit_chk("Position#3", "Check[Position[], fired, Position::argb]", "fired");
    audit_chk("Position#4", "Position[{a}]", "Position[{a}]");
    audit_chk("Position#4", "Head[Position[{a}]]", "Position");
    audit_chk("Position#4", "Position[a][{b, a, c}]", "{{2}}");
    audit_run("ClearAll[w5df]");
    audit_chk("Default#3", "Attributes[Default]", "{Protected}");
    audit_chk("Default#4", "Default[w5df] = 7; DefaultValues[w5df]", "{HoldPattern[Default[w5df]] :> 7}");
    audit_chk("Default#4", "w5dg[x_, y_.] := {x, y}; Default[w5dg] = 5; w5dg[1]", "{1, 5}");
    audit_chk("Default#4", "Default[w5df] =.; DefaultValues[w5df]", "{}");
    audit_chk("DeleteCases#3", "DeleteCases[{f[a], g[b]}, f, Infinity, Heads :> True]", "{a, g[b]}");
    audit_chk("FirstPosition#3", "FirstPosition[SparseArray[{0, 5}], 5]", "{2}");
    audit_chk("FirstPosition#3", "FirstPosition[SparseArray[{{1, 0}, {0, 7}}], 7]", "{2, 2}");
    audit_chk("ReplaceList#1", "Check[ReplaceList[{a, b, c}, {___, x_, ___} :> x, -1], fired, ReplaceList::innf]", "fired");
    audit_chk("ReplaceList#1", "Quiet[ReplaceList[{a, b, c}, {___, x_, ___} :> x, n]]", "ReplaceList[{a, b, c}, {___, x_, ___} :> x, n]");
    audit_chk("ReplaceList#1", "ReplaceList[{a, b, c}, {___, x_, ___} :> x, 2]", "{a, b}");
    audit_chk("ReplaceList#2", "Check[ReplaceList[x, 1], fired, ReplaceList::reps]", "fired");
    audit_chk("ReplaceList#2", "Quiet[ReplaceList[x, 1]]", "ReplaceList[x, 1]");
    audit_chk("Verbatim#1", "StringMatchQ[\"a*\", Verbatim[\"a*\"]]", "True");
    audit_chk("Verbatim#1", "StringMatchQ[\"abc\", Verbatim[\"a*\"]]", "False");
    audit_chk("Verbatim#1", "StringMatchQ[\"abc\", \"a*\"]", "True");
}

int main(void) {
    audit_init();
    TEST(test_condition_on_patternsequence);
    TEST(test_s14_level_specs);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_g4_patterns);
    TEST(test_wave5_g4_patterns);
    return audit_finish("audit_pattern_matching");
}
