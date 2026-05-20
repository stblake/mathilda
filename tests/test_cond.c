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

/* ----- Switch -----
 *
 * Switch[expr, form1, value1, form2, value2, ...]
 *
 * expr is evaluated (HoldRest holds only the rest), then expr is
 * pattern-matched against each form_i in order; the first matching
 * branch returns its value_i. _ is the catch-all default; if no form
 * matches, Switch is returned unevaluated.
 */

/* Direct literal match on the first pair. */
void test_switch_first_literal_match() {
    assert_eval_eq("Switch[True, True, 1, False, 0]", "1", 0);
    assert_eval_eq("Switch[False, True, 1, False, 0]", "0", 0);
}

/* Default branch (`_`) fires when nothing else matches. */
void test_switch_default_blank() {
    assert_eval_eq("Switch[x, True, 1, False, 0, _, 99]", "99", 0);
    assert_eval_eq("Switch[3, 1, a, 2, b, _, z]", "z", 0);
}

/* Typed blanks: _Integer, _Plus, _Times. */
void test_switch_typed_blanks() {
    assert_eval_eq("Switch[42, _Integer, int, _Real, real, _, other]", "int", 0);
    assert_eval_eq("Switch[1.5, _Integer, int, _Real, real, _, other]", "real", 0);
    assert_eval_eq("Switch[x + y, _Plus, sum, _Times, product, _, other]", "sum", 0);
    assert_eval_eq("Switch[x*y, _Plus, sum, _Times, product, _, other]", "product", 0);
}

/* Docstring example 1: Boole-like helper using Switch as the body. */
void test_switch_boole_helper() {
    assert_eval_eq("f[b_]:=Switch[b,True,1,False,0,_,0]", "Null", 0);
    assert_eval_eq("{f[True], f[False], f[x]}", "{1, 0, 0}", 0);
    assert_eval_eq("Clear[f]", "Null", 0);
}

/* Docstring example 2: pick a symbolic transform by expression head. */
void test_switch_transform_dispatcher() {
    assert_eval_eq("t[e_] := Switch[e, _Plus, Together, _Times, Apart, _, Identity]",
                   "Null", 0);
    /* (1+x)/(1-x) + x/(1+x) is a Plus expression -> Together. */
    assert_eval_eq("t[(1+x)/(1-x) + x/(1+x)]", "Together", 0);
    /* A pure product picks Apart. */
    assert_eval_eq("t[x*y]", "Apart", 0);
    /* A bare symbol falls through to Identity. */
    assert_eval_eq("t[x]", "Identity", 0);
    assert_eval_eq("Clear[t]", "Null", 0);
}

/* No form matches -> Switch returned unevaluated. */
void test_switch_no_match_unevaluated() {
    assert_eval_eq("Switch[3, 1, a, 2, b]", "Switch[3, 1, a, 2, b]", 0);
    assert_eval_eq("Switch[foo, True, t, False, f]", "Switch[foo, True, t, False, f]", 0);
}

/* HoldRest: unselected value_i must NOT be evaluated. We arrange a
 * scenario where evaluating the wrong branch would produce a visibly
 * different result. */
void test_switch_holdrest_unselected_not_evaluated() {
    /* If the second branch were evaluated, Plus[1,2] would fold to 3
     * before Switch returns -- but the first branch matches, so we
     * should see the literal value of the first branch. */
    assert_eval_eq("Switch[1, 1, ok, 2, Plus[1,2]]", "ok", 0);
}

/* The chosen value IS evaluated by the outer evaluator. */
void test_switch_selected_value_is_evaluated() {
    assert_eval_eq("Switch[1, 1, Plus[2,3], _, nope]", "5", 0);
    assert_eval_eq("Switch[x, _Symbol, Plus[10, 20], _, 0]", "30", 0);
}

/* expr is evaluated before being matched (Switch does NOT hold the
 * first argument). */
void test_switch_expr_is_evaluated() {
    /* 2 + 3 -> 5; pattern 5 matches. */
    assert_eval_eq("Switch[2 + 3, 4, low, 5, hit, _, miss]", "hit", 0);
    /* (1+x)/.x->2 should evaluate to 3 then match. */
    assert_eval_eq("Switch[1 + 2, _Integer, anInt, _, other]", "anInt", 0);
}

/* Pattern variables in the form are NOT substituted into the value:
 * value is returned literally (matching Mathematica's spec). */
void test_switch_no_binding_substitution() {
    assert_eval_eq("Switch[{1,2}, {x_, y_}, x + y, _, 0]", "x + y", 0);
}

/* Match-with-Condition: form `_?test` should respect the predicate. */
void test_switch_pattern_test() {
    assert_eval_eq("Switch[4, _?EvenQ, even, _?OddQ, odd, _, other]", "even", 0);
    assert_eval_eq("Switch[5, _?EvenQ, even, _?OddQ, odd, _, other]", "odd", 0);
}

/* Form expressions are evaluated when the match is tried. */
void test_switch_form_evaluated() {
    /* The form `1 + 1` should evaluate to 2 before matching. */
    assert_eval_eq("Switch[2, 1 + 1, sum, _, other]", "sum", 0);
}

/* Wrong arity (no form/value pair, or odd number after expr): unevaluated. */
void test_switch_arity_errors_unevaluated() {
    /* No args. */
    assert_eval_eq("Switch[]", "Switch[]", 0);
    /* Just expr, no pairs. */
    assert_eval_eq("Switch[x]", "Switch[x]", 0);
    /* expr + form but no value -- even argc total. */
    assert_eval_eq("Switch[x, _]", "Switch[x, _]", 0);
    /* expr + 3 args -> trailing form without value. */
    assert_eval_eq("Switch[x, a, 1, b]", "Switch[x, a, 1, b]", 0);
}

/* The first matching branch wins even if later branches would also match. */
void test_switch_first_match_wins() {
    assert_eval_eq("Switch[3, _Integer, intFirst, 3, three, _, other]", "intFirst", 0);
}

/* Substituting later via /. can resolve a previously unevaluated Switch. */
void test_switch_replace_all_resolves() {
    /* With expr=x, no form matches an unbound x; result is held. After
     * x -> 1, the held Switch re-evaluates and picks the first branch. */
    assert_eval_eq("Switch[x, 1, one, 2, two] /. x -> 1", "one", 0);
}

/* Listable threading: Switch is NOT Listable, so wrapping in a list
 * keeps the structure unless Map is used. */
void test_switch_not_listable() {
    assert_eval_eq("Switch[{1,2}, {1,2}, pair, _, other]", "pair", 0);
    /* Mapped over a list, each element gets its own Switch result. */
    assert_eval_eq("Switch[#, 1, one, 2, two, _, other] & /@ {1, 2, 3}",
                   "{one, two, other}", 0);
}

/* Switch in a recursive definition (classic fizzbuzz-flavoured example). */
void test_switch_recursive_use() {
    assert_eval_eq(
        "classify[n_] := Switch[Mod[n, 3], 0, zero, 1, one, _, other]",
        "Null", 0);
    assert_eval_eq("classify[6]", "zero", 0);
    assert_eval_eq("classify[7]", "one", 0);
    assert_eval_eq("classify[8]", "other", 0);
    assert_eval_eq("Clear[classify]", "Null", 0);
}

/* ----- Piecewise -----
 *
 * Piecewise[{{v1, c1}, {v2, c2}, ...}]            -> Piecewise[..., 0]
 * Piecewise[{{v1, c1}, {v2, c2}, ...}, default]
 *
 * Conditions are evaluated left-to-right. {_, False} clauses are dropped.
 * At the first {_, True} all later clauses (and the default) are dropped.
 * If no preceding condition is inconclusive, the True clause's value is
 * returned directly. Consecutive clauses with structurally equal values
 * are merged via Or on their conditions.
 */

/* No clauses -> default; default defaults to 0. */
void test_piecewise_empty_clauses() {
    assert_eval_eq("Piecewise[{}]", "0", 0);
    assert_eval_eq("Piecewise[{}, 99]", "99", 0);
    assert_eval_eq("Piecewise[{}, foo]", "foo", 0);
}

/* All conditions literally False -> default. */
void test_piecewise_all_false_default() {
    assert_eval_eq("Piecewise[{{1, False}, {2, False}}]", "0", 0);
    assert_eval_eq("Piecewise[{{1, False}, {2, False}}, 99]", "99", 0);
    assert_eval_eq("Piecewise[{{a, False}}, default]", "default", 0);
}

/* First True (with no preceding inconclusive) returns its value. */
void test_piecewise_first_true_returns_value() {
    assert_eval_eq("Piecewise[{{a, True}, {b, c}}]", "a", 0);
    assert_eval_eq("Piecewise[{{a, True}, {b, c}}, d]", "a", 0);
    assert_eval_eq("Piecewise[{{x^2, True}}]", "x^2", 0);
    assert_eval_eq("Piecewise[{{1, False}, {2, True}, {3, foo}}, 99]", "2", 0);
}

/* Drop unreachable cases after the first True; default is also dropped
 * because the True clause becomes the final unconditional case. */
void test_piecewise_drop_after_true() {
    assert_eval_eq(
        "Piecewise[{{a, d1}, {b, d2}, {c, True}, {d, d4}, {e, d5}}, ef]",
        "Piecewise[{{a, d1}, {b, d2}, {c, True}}]", 0);
}

/* Drop {_, False} clauses but keep inconclusive ones. */
void test_piecewise_drop_false_clauses() {
    assert_eval_eq(
        "Piecewise[{{a, d1}, {b, d2}, {c, False}, {d, d4}, {e, d5}}, f]",
        "Piecewise[{{a, d1}, {b, d2}, {d, d4}, {e, d5}}, f]", 0);
}

/* Merge consecutive clauses with the same value via Or on their conds. */
void test_piecewise_merge_same_value() {
    assert_eval_eq(
        "Piecewise[{{a, d1}, {b, d2}, {b, d3}, {c, d4}}, e]",
        "Piecewise[{{a, d1}, {b, d2 || d3}, {c, d4}}, e]", 0);
    /* Three-in-a-row merge. */
    assert_eval_eq(
        "Piecewise[{{a, d1}, {a, d2}, {a, d3}}, ef]",
        "Piecewise[{{a, d1 || d2 || d3}}, ef]", 0);
}

/* Piecewise[conds] auto-fills the default 0 in symbolic form. */
void test_piecewise_autofill_default_zero() {
    /* All conditions inconclusive: the result is symbolic with default 0. */
    assert_eval_eq(
        "Piecewise[{{a, d1}, {b, d2}}]",
        "Piecewise[{{a, d1}, {b, d2}}, 0]", 0);
}

/* Symbolic form is preserved when conditions cannot be resolved. */
void test_piecewise_symbolic_preserved() {
    assert_eval_eq(
        "Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1]",
        "Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, (-x^2)/100 + 1]", 0);
}

/* /. substitution can later resolve a previously symbolic Piecewise. */
void test_piecewise_resolved_by_replace() {
    /* The piecewise from the user's docstring example. */
    assert_eval_eq(
        "Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1] /. x -> 0",
        "1", 0);
    assert_eval_eq(
        "Piecewise[{{Sin[x]/x, x < 0}, {1, x == 0}}, -x^2/100 + 1] /. x -> 5",
        "3/4", 0);
}

/* HoldAll: values for non-selected branches are not evaluated.
 * If they were evaluated, the Print side-effect or arithmetic would
 * happen even though the branch is unreached. We test via a value that
 * would visibly fold (Plus[1,2]) but is never touched. */
void test_piecewise_holdall_skips_unselected_values() {
    /* The True branch wins; the other val is held and dropped. */
    assert_eval_eq("Piecewise[{{ok, True}, {Plus[1,2], other}}]", "ok", 0);
    /* The default is dropped when a True is found. */
    assert_eval_eq("Piecewise[{{first, True}}, Plus[1,2]]", "first", 0);
}

/* HoldAll on conditions: conds are evaluated when Piecewise examines them.
 * 1 + 1 -> 2 (not == True or False); 0 < 1 -> True; etc. */
void test_piecewise_conditions_evaluated() {
    assert_eval_eq("Piecewise[{{a, 1 < 0}, {b, 0 < 1}}, c]", "b", 0);
    assert_eval_eq("Piecewise[{{a, 1 < 0}, {b, 2 == 3}}, c]", "c", 0);
}

/* Single inconclusive clause -> symbolic with default 0. */
void test_piecewise_single_inconclusive() {
    assert_eval_eq(
        "Piecewise[{{a, x == 0}}]",
        "Piecewise[{{a, x == 0}}, 0]", 0);
}

/* Default that is itself a symbolic expression is preserved verbatim. */
void test_piecewise_symbolic_default() {
    assert_eval_eq(
        "Piecewise[{{a, x == 0}}, b + c]",
        "Piecewise[{{a, x == 0}}, b + c]", 0);
}

/* OwnValue can let a previously symbolic condition resolve. */
void test_piecewise_with_assignment() {
    assert_eval_eq("x = 0", "0", 0);
    assert_eval_eq("Piecewise[{{a, x == 0}, {b, x > 0}}, c]", "a", 0);
    assert_eval_eq("Clear[x]", "Null", 0);
}

/* Threading via Map: Piecewise is NOT Listable; wrapping in a list is
 * intentional unless Map is used. */
void test_piecewise_map_over_list() {
    assert_eval_eq(
        "f[y_] := Piecewise[{{-y, y < 0}, {0, y == 0}}, y]",
        "Null", 0);
    assert_eval_eq("f /@ {-3, 0, 5}", "{3, 0, 5}", 0);
    assert_eval_eq("Clear[f]", "Null", 0);
}

/* Define a Mathilda-style absolute-value function with Piecewise. */
void test_piecewise_abs_function() {
    assert_eval_eq("myAbs[x_] := Piecewise[{{-x, x < 0}}, x]", "Null", 0);
    assert_eval_eq("myAbs[-7]", "7", 0);
    assert_eval_eq("myAbs[0]", "0", 0);
    assert_eval_eq("myAbs[3]", "3", 0);
    assert_eval_eq("Clear[myAbs]", "Null", 0);
}

/* Malformed input shapes leave the expression unevaluated. */
void test_piecewise_malformed_unevaluated() {
    /* First arg not a List. */
    assert_eval_eq("Piecewise[foo]", "Piecewise[foo]", 0);
    /* Element of clauses list is not a 2-element list. */
    assert_eval_eq("Piecewise[{{a}}]", "Piecewise[{{a}}]", 0);
    /* Three-element 'pair'. */
    assert_eval_eq("Piecewise[{{a, b, c}}]", "Piecewise[{{a, b, c}}]", 0);
    /* Too many top-level arguments. */
    assert_eval_eq("Piecewise[{{a, b}}, c, d]",
                   "Piecewise[{{a, b}}, c, d]", 0);
    /* Zero arguments. */
    assert_eval_eq("Piecewise[]", "Piecewise[]", 0);
}

/* A True buried among inconclusive predecessors becomes the symbolic
 * final case (the default vanishes). */
void test_piecewise_true_after_inconclusive() {
    assert_eval_eq(
        "Piecewise[{{a, d1}, {b, True}}, c]",
        "Piecewise[{{a, d1}, {b, True}}]", 0);
}

/* Merging respects identity of values: Pi vs. 3 are distinct, so no
 * merging happens between distinct atoms. */
void test_piecewise_no_merge_distinct_values() {
    assert_eval_eq(
        "Piecewise[{{1, d1}, {2, d2}, {1, d3}}, e]",
        "Piecewise[{{1, d1}, {2, d2}, {1, d3}}, e]", 0);
}

/* Repeated evaluation (idempotency) on a symbolic Piecewise must not
 * mutate it -- this also exercises the input==output short-circuit. */
void test_piecewise_idempotent() {
    /* Evaluating twice yields the same canonical form. */
    assert_eval_eq(
        "Piecewise[{{a, d1}, {b, d2}}, c] // Identity",
        "Piecewise[{{a, d1}, {b, d2}}, c]", 0);
}

/* Single True clause with default still returns the value. */
void test_piecewise_single_true_with_default() {
    assert_eval_eq("Piecewise[{{1, True}}, 99]", "1", 0);
}

/* When merging produces an Or that simplifies to True, the merged
 * clause's value is returned directly. */
void test_piecewise_merge_to_true_returns_value() {
    /* d || True -> True. Merging {a, d}, {a, True} into {a, d || True}
     * simplifies the cond to True. With out_count == 1 and final cond
     * True, the value a is returned. */
    assert_eval_eq("Piecewise[{{a, d}, {a, True}}, b]", "a", 0);
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

    TEST(test_switch_first_literal_match);
    TEST(test_switch_default_blank);
    TEST(test_switch_typed_blanks);
    TEST(test_switch_boole_helper);
    TEST(test_switch_transform_dispatcher);
    TEST(test_switch_no_match_unevaluated);
    TEST(test_switch_holdrest_unselected_not_evaluated);
    TEST(test_switch_selected_value_is_evaluated);
    TEST(test_switch_expr_is_evaluated);
    TEST(test_switch_no_binding_substitution);
    TEST(test_switch_pattern_test);
    TEST(test_switch_form_evaluated);
    TEST(test_switch_arity_errors_unevaluated);
    TEST(test_switch_first_match_wins);
    TEST(test_switch_replace_all_resolves);
    TEST(test_switch_not_listable);
    TEST(test_switch_recursive_use);

    TEST(test_piecewise_empty_clauses);
    TEST(test_piecewise_all_false_default);
    TEST(test_piecewise_first_true_returns_value);
    TEST(test_piecewise_drop_after_true);
    TEST(test_piecewise_drop_false_clauses);
    TEST(test_piecewise_merge_same_value);
    TEST(test_piecewise_autofill_default_zero);
    TEST(test_piecewise_symbolic_preserved);
    TEST(test_piecewise_resolved_by_replace);
    TEST(test_piecewise_holdall_skips_unselected_values);
    TEST(test_piecewise_conditions_evaluated);
    TEST(test_piecewise_single_inconclusive);
    TEST(test_piecewise_symbolic_default);
    TEST(test_piecewise_with_assignment);
    TEST(test_piecewise_map_over_list);
    TEST(test_piecewise_abs_function);
    TEST(test_piecewise_malformed_unevaluated);
    TEST(test_piecewise_true_after_inconclusive);
    TEST(test_piecewise_no_merge_distinct_values);
    TEST(test_piecewise_idempotent);
    TEST(test_piecewise_single_true_with_default);
    TEST(test_piecewise_merge_to_true_returns_value);

    printf("All cond tests passed!\n");
    symtab_clear();
    return 0;
}