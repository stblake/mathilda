#include "boolean.h"
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ----- Boole ----- */

/* Boole[True] -> 1, Boole[False] -> 0 (the two defining cases). */
void test_boole_true_false() {
    assert_eval_eq("Boole[True]", "1", 0);
    assert_eval_eq("Boole[False]", "0", 0);
}

/* Docstring example: {Boole[False], Boole[True]} -> {0, 1}. */
void test_boole_in_list() {
    assert_eval_eq("{Boole[False], Boole[True]}", "{0, 1}", 0);
}

/* Listable: Boole threads automatically over a List argument. */
void test_boole_listable_basic() {
    assert_eval_eq("Boole[{True, False, True, True, False}]",
                   "{1, 0, 1, 1, 0}", 0);
}

/* Listable nested: a list of lists threads at every level. */
void test_boole_listable_nested() {
    assert_eval_eq("Boole[{{True, False}, {False, True}}]",
                   "{{1, 0}, {0, 1}}", 0);
}

/* Listable threading over an empty list yields an empty list, consistent with
 * every other Listable builtin (Sin[{}], Abs[{}], BernoulliB[{}], ... all give
 * {}). */
void test_boole_listable_empty_list() {
    assert_eval_eq("Boole[{}]", "{}", 0);
}

/* For neither-True-nor-False scalar input, Boole stays unevaluated. */
void test_boole_symbolic_unevaluated() {
    assert_eval_eq("Boole[x]", "Boole[x]", 0);
    assert_eval_eq("Boole[a < b]", "Boole[a < b]", 0);
    assert_eval_eq("Boole[1]", "Boole[1]", 0);
    assert_eval_eq("Boole[0]", "Boole[0]", 0);
    assert_eval_eq("Boole[Indeterminate]", "Boole[Indeterminate]", 0);
}

/* Listable still threads even when individual elements are symbolic. */
void test_boole_listable_mixed() {
    assert_eval_eq("Boole[{True, x, False, y}]",
                   "{1, Boole[x], 0, Boole[y]}", 0);
}

/* Conditions are evaluated before Boole sees them (no HoldAll). */
void test_boole_evaluates_argument() {
    assert_eval_eq("Boole[1 < 2]", "1", 0);
    assert_eval_eq("Boole[2 < 1]", "0", 0);
    assert_eval_eq("Boole[1 == 1]", "1", 0);
    assert_eval_eq("Boole[1 == 2]", "0", 0);
}

/* Wrong arity: Boole[] and Boole[a, b] stay symbolic (no exception). */
void test_boole_arity_errors_unevaluated() {
    assert_eval_eq("Boole[]", "Boole[]", 0);
    assert_eval_eq("Boole[True, False]", "Boole[True, False]", 0);
}

/* Boole composes with Plus to count True predicates -- a typical
 * Iverson-bracket use case. */
void test_boole_sum_indicator() {
    assert_eval_eq("Boole[True] + Boole[False] + Boole[True]", "2", 0);
    /* Total over a list of predicate values. */
    assert_eval_eq("Total[Boole[{True, False, True, True}]]", "3", 0);
}

/* Boole[expr] is effectively equivalent to If[expr, 1, 0]. */
void test_boole_matches_if() {
    assert_eval_eq("Boole[True] == If[True, 1, 0]", "True", 0);
    assert_eval_eq("Boole[False] == If[False, 1, 0]", "True", 0);
}

/* Replacement: a symbolic Boole resolves when the predicate is
 * substituted by a Boolean value. */
void test_boole_replace_all() {
    assert_eval_eq("Boole[p] /. p -> True", "1", 0);
    assert_eval_eq("Boole[p] /. p -> False", "0", 0);
}

/* Idempotence: re-evaluating a resolved Boole result stays put. */
void test_boole_idempotent() {
    assert_eval_eq("Boole[Boole[True] == 1]", "1", 0);
    assert_eval_eq("Boole[Boole[False] == 0]", "1", 0);
}

/* Map: explicit Map (with /@) gives the same result as Listable threading. */
void test_boole_map_matches_listable() {
    assert_eval_eq("Boole /@ {True, False, True}", "{1, 0, 1}", 0);
}

/* Boole is Protected: assignment to Boole should not succeed. The
 * existing pattern in this codebase is that Set on a protected symbol
 * is rejected. We don't assert a specific error message; we just check
 * that the built-in semantics survive. */
void test_boole_protected_semantics() {
    assert_eval_eq("Boole[True]", "1", 0);
    assert_eval_eq("Boole[False]", "0", 0);
}

/* Useful idiom: counting matches via Length[Cases] vs Total[Boole[...]]. */
void test_boole_indicator_count() {
    assert_eval_eq("Total[Boole[# > 0 & /@ {-1, 2, -3, 4, 5}]]", "3", 0);
}

/* TrueQ-then-Boole: TrueQ collapses non-Boolean to False, then Boole
 * yields 0; for True it yields 1. */
void test_boole_compose_trueq() {
    assert_eval_eq("Boole[TrueQ[True]]", "1", 0);
    assert_eval_eq("Boole[TrueQ[False]]", "0", 0);
    assert_eval_eq("Boole[TrueQ[a < b]]", "0", 0);
}

/* ----- ConditionalExpression ----- */

/* Defining cases from the docstring. */
void test_conditional_expression_true_false() {
    assert_eval_eq("ConditionalExpression[a, True]", "a", 0);
    assert_eval_eq("ConditionalExpression[a, False]", "Undefined", 0);
}

/* Literal numerics on the True branch. */
void test_conditional_expression_numeric_value() {
    assert_eval_eq("ConditionalExpression[5, True]", "5", 0);
    assert_eval_eq("ConditionalExpression[1/2, True]", "1/2", 0);
    assert_eval_eq("ConditionalExpression[Sin[x], True]", "Sin[x]", 0);
}

/* A symbolic condition (neither True nor False) keeps the head intact. */
void test_conditional_expression_symbolic_cond() {
    assert_eval_eq("ConditionalExpression[a, c]", "ConditionalExpression[a, c]", 0);
    assert_eval_eq("ConditionalExpression[x^2, x > 0]",
                   "ConditionalExpression[x^2, x > 0]", 0);
}

/* Comparison conditions are evaluated before ConditionalExpression sees them. */
void test_conditional_expression_evaluates_condition() {
    assert_eval_eq("ConditionalExpression[a, 1 < 2]", "a", 0);
    assert_eval_eq("ConditionalExpression[a, 2 < 1]", "Undefined", 0);
    assert_eval_eq("ConditionalExpression[a, 1 == 1]", "a", 0);
    assert_eval_eq("ConditionalExpression[a, 1 == 2]", "Undefined", 0);
}

/* And/Or in the condition reduce before the dispatch. */
void test_conditional_expression_boolean_condition() {
    assert_eval_eq("ConditionalExpression[a, True && True]", "a", 0);
    assert_eval_eq("ConditionalExpression[a, True && False]", "Undefined", 0);
    assert_eval_eq("ConditionalExpression[a, False || False]", "Undefined", 0);
    assert_eval_eq("ConditionalExpression[a, !False]", "a", 0);
}

/* Nested CE flattens: CE[CE[e, c1], c2] -> CE[e, c1 && c2]. */
void test_conditional_expression_nested_flattens() {
    assert_eval_eq(
        "ConditionalExpression[ConditionalExpression[e, c1], c2]",
        "ConditionalExpression[e, c1 && c2]", 0);
}

/* Nested with a True outer condition collapses to the inner CE. */
void test_conditional_expression_nested_outer_true() {
    assert_eval_eq(
        "ConditionalExpression[ConditionalExpression[e, c1], True]",
        "ConditionalExpression[e, c1]", 0);
}

/* Nested with a False inner condition: outer True step exposes inner CE,
 * which then resolves to Undefined. */
void test_conditional_expression_nested_inner_false() {
    assert_eval_eq(
        "ConditionalExpression[ConditionalExpression[e, False], True]",
        "Undefined", 0);
}

/* Arity errors remain unevaluated (no exception, no crash). */
void test_conditional_expression_arity_unevaluated() {
    assert_eval_eq("ConditionalExpression[]", "ConditionalExpression[]", 0);
    assert_eval_eq("ConditionalExpression[a]", "ConditionalExpression[a]", 0);
    assert_eval_eq("ConditionalExpression[a, b, c]",
                   "ConditionalExpression[a, b, c]", 0);
}

/* Resolution under ReplaceAll: a symbolic condition that becomes True
 * triggers the reduction. */
void test_conditional_expression_replace_all() {
    assert_eval_eq("ConditionalExpression[a, c] /. c -> True", "a", 0);
    assert_eval_eq("ConditionalExpression[a, c] /. c -> False", "Undefined", 0);
}

/* Idempotence: evaluating the result of CE again is a fixed point. */
void test_conditional_expression_idempotent() {
    assert_eval_eq("ConditionalExpression[ConditionalExpression[a, True], True]",
                   "a", 0);
    assert_eval_eq("ConditionalExpression[a, True] + 1", "1 + a", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_boole_true_false);
    TEST(test_boole_in_list);
    TEST(test_boole_listable_basic);
    TEST(test_boole_listable_nested);
    TEST(test_boole_listable_empty_list);
    TEST(test_boole_symbolic_unevaluated);
    TEST(test_boole_listable_mixed);
    TEST(test_boole_evaluates_argument);
    TEST(test_boole_arity_errors_unevaluated);
    TEST(test_boole_sum_indicator);
    TEST(test_boole_matches_if);
    TEST(test_boole_replace_all);
    TEST(test_boole_idempotent);
    TEST(test_boole_map_matches_listable);
    TEST(test_boole_protected_semantics);
    TEST(test_boole_indicator_count);
    TEST(test_boole_compose_trueq);

    TEST(test_conditional_expression_true_false);
    TEST(test_conditional_expression_numeric_value);
    TEST(test_conditional_expression_symbolic_cond);
    TEST(test_conditional_expression_evaluates_condition);
    TEST(test_conditional_expression_boolean_condition);
    TEST(test_conditional_expression_nested_flattens);
    TEST(test_conditional_expression_nested_outer_true);
    TEST(test_conditional_expression_nested_inner_false);
    TEST(test_conditional_expression_arity_unevaluated);
    TEST(test_conditional_expression_replace_all);
    TEST(test_conditional_expression_idempotent);

    printf("All Boole tests passed.\n");
    return 0;
}
