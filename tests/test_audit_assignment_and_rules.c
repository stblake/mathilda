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

int main(void) {
    audit_init();
    TEST(test_op_assign_part_once);
    TEST(test_replacepart_first_rule);
    TEST(test_replaceat_delayed_unevaluated);
    return audit_finish("audit_assignment_and_rules");
}
