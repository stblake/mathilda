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

int main(void) {
    audit_init();
    TEST(test_op_assign_part_once);
    TEST(test_replacepart_first_rule);
    TEST(test_replaceat_delayed_unevaluated);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_assignment_and_rules");
}
