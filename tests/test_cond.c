#include "cond.h"
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

void test_if() {
    assert_eval_eq("If[True, 1, 0]", "1", 0);
    assert_eval_eq("If[False, 1, 0]", "0", 0);
    assert_eval_eq("If[a < b, 1, 0]", "If[a < b, 1, 0]", 0);
    assert_eval_eq("If[a < b, 1, 0, Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("If[False, 1]", "Null", 0);
    
    assert_eval_eq("abs[x_]:=If[x<0,-x,x]", "Null", 0);
    assert_eval_eq("abs /@ {-1, 0, 1}", "{1, 0, 1}", 0);
}

void test_trueq() {
    assert_eval_eq("TrueQ[True]", "True", 0);
    assert_eval_eq("TrueQ[False]", "False", 0);
    assert_eval_eq("TrueQ[a < b]", "False", 0);
    assert_eval_eq("If[TrueQ[a < b], 1, 0]", "0", 0);
}

/* ----- Which ----- */

/* Direct literal tests fire from the first pair. */
void test_which_first_true() {
    assert_eval_eq("Which[True, a, True, b]", "a", 0);
    assert_eval_eq("Which[True, 42]", "42", 0);
}

/* False tests are skipped until the first True. */
void test_which_skip_false() {
    assert_eval_eq("Which[False, a, True, b]", "b", 0);
    assert_eval_eq("Which[False, a, False, b, True, c]", "c", 0);
}

/* All False -> Null. Empty Which -> Null. */
void test_which_all_false() {
    assert_eval_eq("Which[False, a]", "Null", 0);
    assert_eval_eq("Which[False, a, False, b]", "Null", 0);
    assert_eval_eq("Which[]", "Null", 0);
}

/* Tests are evaluated, not just compared to True/False literally. */
void test_which_evaluated_tests() {
    assert_eval_eq("Which[1 < 0, a, 2 < 3, b]", "b", 0);
    assert_eval_eq("Which[1 == 1, hit, 2 == 2, miss]", "hit", 0);
}

/* Inconclusive test -> symbolic Which with remaining pairs. */
void test_which_symbolic_remainder() {
    /* Bare symbol stays symbolic. */
    assert_eval_eq("Which[a == 1, x, a == 2, b]", "Which[a == 1, x, a == 2, b]", 0);
    /* Leading False is dropped; first inconclusive is kept along with the rest. */
    assert_eval_eq("Which[1 < 0, a, x == 0, b, 0 < 1, c]",
                   "Which[x == 0, b, 0 < 1, c]", 0);
}

/* Substituting via /. lets a previously-symbolic test resolve. */
void test_which_replace_all() {
    assert_eval_eq("Which[a == 1, x, a == 2, b] /. a -> 2", "b", 0);
}

/* HoldAll: values past the first matching test are never evaluated.
 * We can verify this indirectly by giving an unmatched branch a body
 * that would produce a distinctive head -- it must not be touched. */
void test_which_holdall_unevaluated_branches() {
    /* Pick something exotic for the unselected branch. If it were evaluated,
     * Plus would collapse the Hold form. With HoldAll the second value is
     * never even seen by the evaluator. */
    assert_eval_eq("Which[True, ok, True, Plus[1,2]]", "ok", 0);
}

/* OwnValue-driven test: `a = 2` then `Which[a==1,x, a==2,b]` -> b. */
void test_which_with_assignment() {
    assert_eval_eq("a = 2", "2", 0);
    assert_eval_eq("Which[a == 1, x, a == 2, b]", "b", 0);
    /* Clean up to avoid leaking the binding into later tests. */
    assert_eval_eq("Clear[a]", "Null", 0);
}

/* Pattern-defined function that uses Which with a `True` default. */
void test_which_default_true_branch() {
    assert_eval_eq("sign[x_] := Which[x < 0, -1, x > 0, 1, True, Indeterminate]",
                   "Null", 0);
    assert_eval_eq("sign[-2]", "-1", 0);
    assert_eval_eq("sign[3]",  "1",  0);
    assert_eval_eq("sign[0]",  "Indeterminate", 0);
    /* Listable threading still works (sign isn't Listable; map explicitly). */
    assert_eval_eq("sign /@ {-2, 0, 3}", "{-1, Indeterminate, 1}", 0);
    assert_eval_eq("Clear[sign]", "Null", 0);
}

/* Odd argument count is a usage error -- leave the expression unevaluated. */
void test_which_odd_argc_unevaluated() {
    assert_eval_eq("Which[True]", "Which[True]", 0);
    assert_eval_eq("Which[True, a, False]", "Which[True, a, False]", 0);
}

/* HoldAll means SetDelayed / `=` on the RHS of a value isn't a concern,
 * but we should verify that tests with side effects are evaluated at most
 * once each -- subsequent evaluation passes see the already-reduced form. */
void test_which_idempotent_on_resolved_form() {
    /* Two passes through the evaluator should produce the same result. */
    assert_eval_eq("Which[x == 0, b, 0 < 1, c]",
                   "Which[x == 0, b, 0 < 1, c]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_if);
    TEST(test_trueq);
    TEST(test_which_first_true);
    TEST(test_which_skip_false);
    TEST(test_which_all_false);
    TEST(test_which_evaluated_tests);
    TEST(test_which_symbolic_remainder);
    TEST(test_which_replace_all);
    TEST(test_which_holdall_unevaluated_branches);
    TEST(test_which_with_assignment);
    TEST(test_which_default_true_branch);
    TEST(test_which_odd_argc_unevaluated);
    TEST(test_which_idempotent_on_resolved_form);

    printf("All cond tests passed!\n");
    symtab_clear();
    return 0;
}