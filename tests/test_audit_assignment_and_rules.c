/* Builtin-audit regressions: assignment and rules (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* DivideBy#1 / TimesBy#1 (high): "Part index with side effects is evaluated
 * twice: reads one element, writes another". */
static void test_op_assign_part_once(void) {
    audit_run("ClearAll[ai, al]; ai = 1; al = {8, 16, 32}");
    audit_chk("DivideBy#1", "{al[[ai++]] /= 2, ai, al}", "{4, 2, {4, 16, 32}}");
    audit_run("ai = 1; al = {1, 2, 3}");
    audit_chk("TimesBy#1", "al[[ai++]] *= 10; {al, ai}", "{{10, 2, 3}, 2}");
    audit_run("ai = 1; al = {1, 2, 3}");
    audit_chk("TimesBy#1", "al[[ai++]] += 10; {al, ai}", "{{11, 2, 3}, 2}");
    audit_run("ai = 1; al = {1, 2, 3}");
    audit_chk("TimesBy#1", "al[[ai++]] -= 10; {al, ai}", "{{-9, 2, 3}, 2}");
    /* sanity */
    audit_run("al = {1, 2, 3}");
    audit_chk("TimesBy#1", "al[[3]] *= 2; al", "{1, 2, 6}");
    audit_run("ai = 3");
    audit_chk("TimesBy#1", "ai /= 2; ai", "3/2");
}

/* ReplacePart#1 (high): "When several rules match the same position, the
 * last rule wins instead of the first". */
static void test_replacepart_first_rule(void) {
    audit_chk("ReplacePart#1", "ReplacePart[{a, b, c, d}, {1 -> x, _ -> y}]", "{x, y, y, y}");
    audit_chk("ReplacePart#1", "ReplacePart[{a, b, c}, {1 -> x, 1 -> y}]", "{x, b, c}");
    audit_chk("ReplacePart#1", "ReplacePart[{a, b, c}, {{1}, {3}} -> x]", "{x, b, x}");
    audit_chk("ReplacePart#1", "ReplacePart[{a, b, c}, {2 -> x, 3 -> y}]", "{a, x, y}");
}

/* ReplaceAt#1 (high): "RuleDelayed RHS is evaluated before insertion". */
static void test_replaceat_delayed_unevaluated(void) {
    audit_chk("ReplaceAt#1", "ReplaceAt[Hold[1 + 1], x_Plus :> Length[x], 1]", "Hold[Length[1 + 1]]");
    audit_chk("ReplaceAt#1", "ReplaceAt[Hold[a, b], a :> 1 + 1, 1]", "Hold[1 + 1, b]");
    /* sanity: outside a Hold the result evaluates as usual */
    audit_chk("ReplaceAt#1", "ReplaceAt[{a, b}, a :> 1 + 1, 1]", "{2, b}");
    audit_chk("ReplaceAt#1", "ReplaceAt[{1, 2}, x_ :> x^2, 2]", "{1, 4}");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("ReplaceAll#3", "Quiet[Check[ReplaceAll[a, b, c], fired, ReplaceAll::argt]]", "fired");
    audit_chk("ReplaceAll#3", "Quiet[ReplaceAll[a, b, c]]", "ReplaceAll[a, b, c]");
    audit_chk("SetDelayed#3", "Quiet[Check[SetDelayed[a], fired, SetDelayed::argr]]", "fired");
    audit_chk("SetDelayed#3", "Quiet[SetDelayed[a]]", "SetDelayed[a]");
    audit_chk("AddTo#3", "Quiet[Check[AddTo[z, 1, 2], fired, AddTo::argrx]]", "fired");
    audit_chk("AddTo#3", "Quiet[AddTo[z, 1, 2]]", "AddTo[z, 1, 2]");
    audit_chk("Dispatch#1", "Quiet[Check[Dispatch[{a -> 1}, {b -> 2}], fired, Dispatch::argx]]", "fired");
    audit_chk("Dispatch#1", "Quiet[Dispatch[{a -> 1}, {b -> 2}]]", "Dispatch[{a -> 1}, {b -> 2}]");
    audit_chk("OptionQ#1", "Quiet[Check[OptionQ[Sequence[a -> 1, b -> 2]], fired, OptionQ::argx]]", "fired");
    audit_chk("OptionQ#1", "Quiet[OptionQ[Sequence[a -> 1, b -> 2]]]", "OptionQ[a -> 1, b -> 2]");
    audit_chk("SubtractFrom#3", "Quiet[Check[SubtractFrom[z], fired, SubtractFrom::argr]]", "fired");
    audit_chk("SubtractFrom#3", "Quiet[SubtractFrom[z]]", "SubtractFrom[z]");
    audit_chk("TimesBy#4", "Quiet[Check[TimesBy[q], fired, TimesBy::argr]]", "fired");
    audit_chk("TimesBy#4", "Quiet[TimesBy[q]]", "TimesBy[q]");
    audit_chk("Replace#4", "Quiet[Check[Replace[{x}, x -> 1, {1}, Foo -> True], fired, Replace::optx]]", "fired");
    audit_chk("Replace#4", "Quiet[Replace[{x}, x -> 1, {1}, Foo -> True]]", "Replace[{x}, x -> 1, {1}, Foo -> True]");
    audit_chk("SetDelayed#1", "Quiet[Check[Sin[s12y_] := 7, fired, SetDelayed::write]]", "fired");
    audit_chk("Replace#4", "Replace[f[x], f -> g, {0, 1}, Heads -> True]", "g[x]");
    audit_chk("SetDelayed#1", "Quiet[Pi := 3]", "$Failed");
    audit_chk("SetDelayed#1", "Quiet[Sin[s12y_] := 7]", "$Failed");
    audit_chk("SetDelayed#1", "Quiet[Pi = 3]", "3");
}

/* Wave 4 G4: ClearAll/Remove drop messages and options, string patterns in
 * Protect/Unprotect, ::reps for non-rule replacement specs, more operator
 * forms, Attributes[s] = list, ReplacePart through a head and with Heads,
 * op-assignment to indexed l-values, Format[...] =., OptionValue lists,
 * all-or-nothing SetOptions, and Unset over a list. */
static void test_wave4_g4_assignment(void) {
    audit_run("ClearAll[fco, mco, xa1, xa2, fa1, fa2, f, x, y, a, b, c, g3]");
    audit_run("Options[fco] = {a -> 1}; ClearAll[fco]");
    audit_chk("ClearAll#1", "Options[fco]", "{}");
    audit_chk("Options#2", "Options[fco]", "{}");
    audit_run("mco::err = \"bad\"; ClearAll[mco]");
    audit_chk("ClearAll#2", "mco::err", "mco::err");
    audit_run("xa1 = 1; xa2 = 2");
    audit_chk("Protect#1", "Protect[\"xa*\"]", "{\"xa1\", \"xa2\"}");
    audit_chk("Protect#1", "{Attributes[xa1], Names[\"xa*\"]}", "{{Protected}, {\"xa1\", \"xa2\"}}");
    audit_run("Protect[fa1, fa2]");
    audit_chk("Unprotect#1", "Unprotect[\"fa*\"]", "{\"fa1\", \"fa2\"}");
    audit_chk("Unprotect#1", "{Attributes[fa1], Names[\"fa*\"]}", "{{}, {\"fa1\", \"fa2\"}}");
    audit_chk("Replace#1", "Check[Replace[x, 1], fired, Replace::reps]", "fired");
    audit_chk("Replace#1", "Quiet[{Replace[{a, b}, {a -> 1, b}], Replace[x, {x -> 1, {x -> 2}}]}]", "{Replace[{a, b}, {a -> 1, b}], Replace[x, {x -> 1, {x -> 2}}]}");
    audit_chk("Replace#1", "Replace[x, {{x -> 1}, {x -> 2}}]", "{1, 2}");
    audit_chk("ReplaceAll#1", "Check[{x, y} /. 5, fired, ReplaceAll::reps]", "fired");
    audit_chk("ReplaceAll#1", "Quiet[{{x, y} /. {x -> 1, 7}, x /. a}]", "{ReplaceAll[{x, y}, {x -> 1, 7}], ReplaceAll[x, a]}");
    audit_chk("ReplaceAll#1", "{x /. {}, x /. {{x -> 1}, {}}}", "{x, {1, x}}");
    audit_chk("ReplaceRepeated#1", "ReplaceRepeated[x -> y][x]", "y");
    audit_chk("ReplaceRepeated#1", "ReplaceRepeated[{x -> a, a -> b}][{x}]", "{b}");
    audit_chk("Set#2", "Attributes[g3] = {Listable}", "{Listable}");
    audit_chk("Set#2", "g3[{1, 2}]", "{g3[1], g3[2]}");
    audit_chk("ReplacePart#2", "ReplacePart[f[a][b], {0, 1} -> z]", "f[z][b]");
    audit_chk("ReplacePart#2", "ReplacePart[f[g][x], {0, 0} -> h]", "h[g][x]");
    audit_chk("ReplacePart#3", "ReplacePart[{a, b, c}, 2 -> x, Heads -> True]", "{a, x, c}");
    audit_run("ClearAll[ffg4, cg4, asg4]; ffg4[1] = 5");
    audit_chk("AddTo#2", "ffg4[1] += 2", "7");
    audit_chk("AddTo#2", "ffg4[1]", "7");
    audit_run("asg4 = <|\"a\" -> 1|>");
    audit_chk("AddTo#2", "{asg4[\"a\"] += 5, asg4}", "{6, <|\"a\" -> 6|>}");
    audit_chk("AddTo#2", "Check[ffg4[9] += 1, fired, AddTo::rvalue]", "fired");
    audit_chk("SubtractFrom#2", "ffg4[1] -= 2", "5");
    audit_chk("TimesBy#2", "{ffg4[1] *= 2, ffg4[1]}", "{10, 10}");
    audit_run("cg4[_] = 0; cg4[1] = 6");
    audit_chk("DivideBy#2", "{cg4[1] /= 2, cg4[1], cg4[7]}", "{3, 3, 0}");
    audit_run("ClearAll[u3g4, kg4]; Format[u3g4[x_]] := \"U\"");
    audit_chk("Format#2", "Unset[Format[u3g4[x_]]]; ToString[u3g4[1]]", "\"u3g4[1]\"");
    audit_run("Format[kg4[x_]] := x + 1");
    audit_chk("FormatValues#1", "Format[kg4[x_]] =.", "Null");
    audit_chk("FormatValues#1", "FormatValues[kg4]", "{}");
    audit_run("ClearAll[fov4]; Options[fov4] = {a -> 1, b -> 2}");
    audit_chk("OptionValue#1", "OptionValue[fov4, {a, b}]", "{1, 2}");
    audit_chk("OptionValue#1", "OptionValue[fov4, {a, nope}]", "OptionValue[fov4, {a, nope}]");
    audit_run("ClearAll[opg4]; opg4::usage = \"hi\"; Remove[opg4]");
    audit_chk("Remove#1", "opg4::usage", "opg4::usage");
    audit_chk("ReplaceAt#2", "ReplaceAt[a -> xx, 2][{a, a, a}]", "{a, xx, a}");
    audit_chk("ReplaceAt#3", "Check[ReplaceAt[{a, b, c}, foo, 1], fired, ReplaceAt::reps]", "fired");
    audit_chk("ReplaceAt#3", "Quiet[ReplaceAt[{a, b}, \"s\", 1]]", "ReplaceAt[{a, b}, \"s\", 1]");
    audit_run("ClearAll[gso4, hso4]; Options[gso4] = {a -> 1}; Options[hso4] = {b -> 1}");
    audit_chk("SetOptions#1", "Quiet[SetOptions[{gso4, hso4}, a -> 99]]", "SetOptions[{gso4, hso4}, a -> 99]");
    audit_chk("SetOptions#1", "Options[gso4]", "{a -> 1}");
    audit_run("ClearAll[pg4, qg4]; pg4 = 1; qg4 = 2");
    audit_chk("Unset#1", "{pg4, qg4} =.", "{Null, Null}");
    audit_chk("Unset#1", "{pg4, qg4}", "{pg4, qg4}");
    audit_run("ClearAll[lg4]; lg4 = {1, 2, 3}");
    audit_chk("Set#1", "Check[lg4[[5]] = 1, fired, Set::partw]", "fired");
    audit_chk("Set#1", "Check[lg4[[1, 1]] = 9, fired, Set::partd]", "fired");
    audit_chk("Set#1", "{Quiet[lg4[[5]] = 1], lg4}", "{1, {1, 2, 3}}");
    audit_run("ClearAll[e2g4]; Format[e2g4[x_]] := \"pos\" /; x > 0");
    audit_chk("Format#1", "{ToString[e2g4[1]], ToString[e2g4[-1]]}", "{\"pos\", \"e2g4[-1]\"}");
    audit_run("ClearAll[hmb4]; hmb4 /: MakeBoxes[hmb4[x_], fmt_] := \"H\"");
    audit_chk("MakeBoxes#1", "MakeBoxes[hmb4[2], StandardForm]", "\"H\"");
}

int main(void) {
    audit_init();
    TEST(test_op_assign_part_once);
    TEST(test_replacepart_first_rule);
    TEST(test_replaceat_delayed_unevaluated);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_g4_assignment);
    return audit_finish("audit_assignment_and_rules");
}
