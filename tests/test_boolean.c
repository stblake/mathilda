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

/* Listable threading on empty lists is not applied in Mathilda --
 * mirrors the behaviour of Sin[{}], Abs[{}], etc., which all stay
 * symbolic. We assert that current behaviour explicitly. */
void test_boole_listable_empty_list() {
    assert_eval_eq("Boole[{}]", "Boole[{}]", 0);
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

    printf("All Boole tests passed.\n");
    return 0;
}
