/*
 * test_return.c -- coverage for Return[] / Return[expr] / Return[expr, h]
 *
 * Mathematica's Return is a marker that propagates outward through the
 * current evaluation tree until it hits a "scope boundary" -- the body
 * of a Function, a Module / Block / With, or a Do / For / While. The
 * boundary strips the marker and yields the wrapped value as its own
 * result. Two argument forms exist:
 *
 *   Return[expr]    -- yield expr from the nearest boundary.
 *   Return[]        -- yield Null from the nearest boundary.
 *   Return[expr, h] -- yield expr from the nearest boundary whose head
 *                      is the symbol h. Boundaries with a different
 *                      head propagate the marker outward unchanged.
 *
 * These tests exercise every supported boundary head (Function, Module,
 * Block, With, Do, For, While), the propagation paths through nested
 * boundaries, the propagation through CompoundExpression and
 * If / Which (which do *not* trap), the top-level no-boundary case, and
 * a handful of edge cases (Return inside argument lists, Return whose
 * value is itself an expression, Return[expr, h] for a non-existent
 * boundary).
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "parse.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Replace assert-based comparisons with hard-fail because the build
 * may set -DNDEBUG (assert becomes a no-op under Release). */
static void check_eq(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) {
        fprintf(stderr, "FAIL: parse failed for %s\n", input);
        exit(1);
    }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, str);
        free(str);
        expr_free(evaluated);
        exit(1);
    }
    free(str);
    expr_free(evaluated);
}

/* -----------------------------------------------------------------
 * Return-as-marker: at top level, with no enclosing boundary, the
 * marker survives evaluation unchanged.
 * ----------------------------------------------------------------- */
static void test_return_at_top_level_is_inert(void) {
    check_eq("Return[5]", "Return[5]");
    check_eq("Return[]", "Return[]");
    check_eq("Return[a + b]", "Return[a + b]");
    /* Args are evaluated (Return has no Hold attribute). */
    check_eq("Return[1 + 2 + 3]", "Return[6]");
}

/* Top-level Return[v, h] also stays inert because no boundary is
 * actually entered. */
static void test_return_two_arg_at_top_is_inert(void) {
    check_eq("Return[5, Module]", "Return[5, Module]");
    check_eq("Return[5, Function]", "Return[5, Function]");
}

/* -----------------------------------------------------------------
 * Function boundary (1-arg form)
 * ----------------------------------------------------------------- */
static void test_return_in_function_simple(void) {
    /* Function[x, Return[x]][5] -> 5 */
    check_eq("Function[x, Return[x]][5]", "5");
    /* Slot form */
    check_eq("(Return[#] &)[7]", "7");
    /* Multi-param form */
    check_eq("Function[{x, y}, Return[x + y]][3, 4]", "7");
    /* Return[] -> Null */
    check_eq("Function[x, Return[]][9]", "Null");
}

static void test_return_in_function_skips_following_statements(void) {
    /* The "after" branch of CompoundExpression must not be evaluated
     * once Return has fired. We use a side-effect to verify. */
    check_eq("Clear[r]; r = 0; "
                   "Function[x, Return[x]; r = 999][42]", "42");
    check_eq("r", "0");
}

static void test_return_in_function_inside_if(void) {
    /* Common pattern: Return inside an If branch. The If does not trap
     * Return; the Function body does. */
    check_eq("Function[x, If[x > 0, Return[positive], Return[negative]]][3]",
                   "positive");
    check_eq("Function[x, If[x > 0, Return[positive], Return[negative]]][-3]",
                   "negative");
}

/* -----------------------------------------------------------------
 * Module / Block / With boundaries (1-arg form)
 * ----------------------------------------------------------------- */
static void test_return_in_module(void) {
    check_eq("Module[{}, Return[42]]", "42");
    check_eq("Module[{}, Return[]]", "Null");
    check_eq("Module[{x = 10}, Return[x + 1]]", "11");
}

static void test_return_in_block(void) {
    check_eq("Block[{}, Return[42]]", "42");
    check_eq("Block[{}, Return[]]", "Null");
    check_eq("Block[{x = 10}, Return[x + 1]]", "11");
}

static void test_return_in_with(void) {
    check_eq("With[{x = 10}, Return[x + 1]]", "11");
    /* With value substitution still happens before Return fires. */
    check_eq("With[{x = 5}, Return[x*x]]", "25");
}

static void test_return_in_module_skips_subsequent_statements(void) {
    check_eq("Clear[r]; r = 0; "
                   "Module[{}, Return[7]; r = 999]", "7");
    check_eq("r", "0");
}

/* -----------------------------------------------------------------
 * Do / For / While boundaries (1-arg form, already supported but
 * re-verified after the refactor).
 * ----------------------------------------------------------------- */
static void test_return_in_do_n_times(void) {
    check_eq("Do[Return[\"hit\"], {3}]", "\"hit\"");
    check_eq("Do[Return[], {3}]", "Null");
    /* Conditional Return inside a Do. */
    check_eq("Do[If[i == 3, Return[i*10]], {i, 5}]", "30");
}

static void test_return_in_do_list_iter(void) {
    check_eq("Do[If[i == b, Return[found]], {i, {a, b, c}}]", "found");
}

static void test_return_in_do_arithmetic_progression(void) {
    /* Continue + Return in the same arithmetic-progression Do.
     * Continue must skip even values; Return triggers on i==5. */
    check_eq("Do[If[Mod[i,2] == 0, Continue[]]; "
                   "If[i == 5, Return[i]], {i, 10}]", "5");
}

static void test_return_in_for(void) {
    check_eq("For[i = 1, i < 10, i = i + 1, If[i == 4, Return[i*100]]]",
                   "400");
}

static void test_return_in_while(void) {
    check_eq("Clear[n]; n = 0; While[True, n = n + 1; If[n == 5, Return[n]]]",
                   "5");
    check_eq("n", "5");
}

/* -----------------------------------------------------------------
 * CompoundExpression does NOT trap Return: it short-circuits the
 * remaining statements but lets the marker propagate outward.
 * ----------------------------------------------------------------- */
static void test_compound_expression_propagates_return(void) {
    /* No enclosing boundary -> Return survives. */
    check_eq("(1; Return[2]; 3)", "Return[2]");
    /* The Return[2] short-circuits both the inner and outer
     * CompoundExpression -- CompoundExpression does not trap Return.
     * We verify in two parts:
     *   (a) the whole line yields Return[2];
     *   (b) by reading `r` *after* that, that `r = 999` was skipped. */
    check_eq("Clear[r]; r = 0; (Return[2]; r = 999)", "Return[2]");
    check_eq("r", "0");
}

static void test_compound_expression_inside_module_lets_return_reach(void) {
    /* The CompoundExpression bubbles Return up to the surrounding
     * Module, which strips it. */
    check_eq("Module[{}, 1; Return[2]; 3]", "2");
    check_eq("Module[{}, 1; Return[]; 3]", "Null");
}

/* -----------------------------------------------------------------
 * If / Which / Switch do not trap Return -- it propagates outward.
 * ----------------------------------------------------------------- */
static void test_if_does_not_trap_return(void) {
    /* No enclosing boundary -> Return survives the If. */
    check_eq("If[True, Return[5], Return[6]]", "Return[5]");
    check_eq("If[False, Return[5], Return[6]]", "Return[6]");
    /* Wrap with Module to confirm the Return does pierce the If. */
    check_eq("Module[{}, If[True, Return[5]]; 999]", "5");
}

/* -----------------------------------------------------------------
 * Two-arg form: Return[expr, h] targets the nearest boundary whose
 * head is h. Boundaries with a different head must propagate.
 * ----------------------------------------------------------------- */
static void test_return_targeted_at_module(void) {
    /* The inner Do does NOT consume Return[5, Module]; the outer
     * Module does. */
    check_eq("Module[{}, Do[Return[5, Module], {3}]]", "5");
}

static void test_return_targeted_at_function(void) {
    check_eq("Function[x, Do[Return[x*100, Function], {3}]][7]",
                   "700");
}

static void test_return_targeted_at_block(void) {
    check_eq("Block[{}, Do[Return[hit, Block], {3}]]", "hit");
}

static void test_return_targeted_at_for(void) {
    /* Outer For traps Return[v, For]; inner Do propagates. */
    check_eq("For[i = 1, i < 10, i = i + 1, "
                   "Do[Return[outer, For], {2}]]", "outer");
}

static void test_return_targeted_at_while(void) {
    check_eq("Clear[n]; n = 0; "
                   "While[True, n = n + 1; "
                   "Do[Return[n*1000, While], {2}]]", "1000");
}

static void test_return_targeted_at_do(void) {
    /* Inner For propagates Return[..., Do]; the outer Do consumes it. */
    check_eq("Do[For[k = 1, k < 5, k = k + 1, Return[caught, Do]], {3}]",
                   "caught");
}

/* -----------------------------------------------------------------
 * Two-arg form, mismatched target: no boundary in the call stack has
 * the requested head, so Return[expr, h] propagates all the way out
 * and is visible at top level.
 * ----------------------------------------------------------------- */
static void test_return_targeted_no_match_propagates_to_top(void) {
    /* No Block in scope -> Return[5, Block] survives the Do AND Module. */
    check_eq("Module[{}, Do[Return[5, Block], {3}]]",
                   "Return[5, Block]");
}

/* -----------------------------------------------------------------
 * Nested boundaries: 1-arg Return is consumed by the *innermost*
 * boundary (Mathematica semantics).
 * ----------------------------------------------------------------- */
static void test_return_one_arg_innermost_boundary(void) {
    /* Inner Do consumes Return[5]; the outer Module sees Null from
     * the Do (since the CompoundExpression yields the Do's value,
     * which is 5 -- not a Return). The Module then yields whatever
     * the body's last statement produced; here we add a sentinel. */
    check_eq("Module[{}, Do[Return[5], {3}]; 99]", "99");
    /* The Do's own value is recoverable. */
    check_eq("Module[{}, Do[Return[5], {3}]]", "5");
}

/* -----------------------------------------------------------------
 * Edge cases.
 * ----------------------------------------------------------------- */
static void test_return_value_is_a_list(void) {
    check_eq("Module[{}, Return[{1, 2, 3}]]", "{1, 2, 3}");
    check_eq("Function[x, Return[{x, x+1}]][7]", "{7, 8}");
}

static void test_return_value_is_unevaluated_symbolic(void) {
    /* Return[a + b] where a and b are unbound: symbolic propagates. */
    check_eq("Module[{}, Return[a + b]]", "a + b");
}

static void test_return_arg_is_evaluated_before_unwind(void) {
    /* Return[v]'s payload is fully evaluated before the marker exits.
     * Here, the inner Plus must be reduced to 6. */
    check_eq("Module[{}, Return[1 + 2 + 3]]", "6");
}

static void test_return_inside_pure_function_inner(void) {
    /* The innermost Function consumes Return[x*2]; the outer Function
     * just yields whatever the inner Function produced. */
    check_eq("Function[y, Function[x, Return[x*2]][y + 1]][4]",
                   "10");
}

static void test_return_targeting_non_existent_head_at_top(void) {
    /* No Foo boundary anywhere -> marker survives. */
    check_eq("Module[{}, Return[5, Foo]]",
                   "Return[5, Foo]");
}

static void test_return_inside_user_defined_function(void) {
    /* User-level f[x_] := body. Without a wrapping Function/Module, the
     * Return marker survives at top level (matches Mathematica). */
    check_eq("Clear[f]; f[x_] := Return[x + 1]; f[5]",
                   "Return[6]");
    /* Wrapping the body in a Module traps it. */
    check_eq("Clear[g]; g[x_] := Module[{}, Return[x + 1]]; g[5]",
                   "6");
}

static void test_return_zero_arg_in_each_boundary(void) {
    check_eq("Function[x, Return[]][1]", "Null");
    check_eq("Module[{}, Return[]]", "Null");
    check_eq("Block[{}, Return[]]", "Null");
    check_eq("With[{x = 1}, Return[]]", "Null");
    check_eq("Do[Return[], {3}]", "Null");
    check_eq("For[i = 1, i < 5, i = i + 1, Return[]]", "Null");
    check_eq("While[True, Return[]]", "Null");
}

static void test_return_two_arg_each_boundary(void) {
    check_eq("Function[x, Return[100, Function]][1]", "100");
    check_eq("Module[{}, Return[100, Module]]", "100");
    check_eq("Block[{}, Return[100, Block]]", "100");
    check_eq("With[{x = 1}, Return[100, With]]", "100");
    check_eq("Do[Return[100, Do], {3}]", "100");
    check_eq("For[i = 1, i < 5, i = i + 1, Return[100, For]]", "100");
    check_eq("While[True, Return[100, While]]", "100");
}

/* -----------------------------------------------------------------
 * Stress: deeply nested mixed boundaries.
 * ----------------------------------------------------------------- */
static void test_return_through_deep_nesting(void) {
    /* Function -> Module -> Do -> For -> While: 1-arg Return is caught
     * by the innermost (While), which yields a*2. The For body is now
     * a non-Return value (12), so For keeps looping until i hits 99 and
     * eventually returns Null. Do then sees Null, iterates 3 times,
     * yields Null. Module yields Null, Function yields Null. So the
     * whole expression yields Null -- matching Mathematica's "Return
     * exits the innermost trapping construct" rule. */
    check_eq("Function[a, "
                   "Module[{}, "
                   "Do["
                   "For[i = 1, i < 99, i = i + 1, "
                   "While[True, Return[a*2]]"
                   "], {3}]]][6]",
                   "Null");

    /* Without For/Do swallowing the value, the While yields 12 and
     * Module simply evaluates to its body's last value (12). The
     * Function yields the same. */
    check_eq("Function[a, Module[{}, While[True, Return[a*2]]]][6]",
             "12");

    /* Same nesting as the first case but Return targets Function
     * explicitly: the marker must skip While, For, Do, Module and
     * arrive at Function, yielding 6*3 = 18. */
    check_eq("Function[a, "
                   "Module[{}, "
                   "Do["
                   "For[i = 1, i < 99, i = i + 1, "
                   "While[True, Return[a*3, Function]]"
                   "], {3}]]][6]",
                   "18");

    /* Targeting Module specifically: While catches if 1-arg, but the
     * 2-arg form with target=Module skips While, For, Do until Module
     * traps it. */
    check_eq("Function[a, "
                   "Module[{}, "
                   "While[True, Return[a*5, Module]]"
                   "]][4]",
                   "20");
}

/* -----------------------------------------------------------------
 * Sanity: Return symbol exists in the symbol table with a docstring
 * (so `?Return` works) and is Protected.
 * ----------------------------------------------------------------- */
static void test_return_symbol_metadata(void) {
    check_eq("Attributes[Return]", "{Protected}");
    /* Docstring should be set so Information / ?Return is non-empty. */
    const char* doc = symtab_get_docstring("Return");
    ASSERT(doc != NULL);
    ASSERT(strstr(doc, "Return[expr]") != NULL);
    ASSERT(strstr(doc, "Return[]") != NULL);
    ASSERT(strstr(doc, "h") != NULL); /* mentions the 2-arg form */
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_return_at_top_level_is_inert);
    TEST(test_return_two_arg_at_top_is_inert);

    TEST(test_return_in_function_simple);
    TEST(test_return_in_function_skips_following_statements);
    TEST(test_return_in_function_inside_if);

    TEST(test_return_in_module);
    TEST(test_return_in_block);
    TEST(test_return_in_with);
    TEST(test_return_in_module_skips_subsequent_statements);

    TEST(test_return_in_do_n_times);
    TEST(test_return_in_do_list_iter);
    TEST(test_return_in_do_arithmetic_progression);
    TEST(test_return_in_for);
    TEST(test_return_in_while);

    TEST(test_compound_expression_propagates_return);
    TEST(test_compound_expression_inside_module_lets_return_reach);

    TEST(test_if_does_not_trap_return);

    TEST(test_return_targeted_at_module);
    TEST(test_return_targeted_at_function);
    TEST(test_return_targeted_at_block);
    TEST(test_return_targeted_at_for);
    TEST(test_return_targeted_at_while);
    TEST(test_return_targeted_at_do);

    TEST(test_return_targeted_no_match_propagates_to_top);

    TEST(test_return_one_arg_innermost_boundary);

    TEST(test_return_value_is_a_list);
    TEST(test_return_value_is_unevaluated_symbolic);
    TEST(test_return_arg_is_evaluated_before_unwind);
    TEST(test_return_inside_pure_function_inner);
    TEST(test_return_targeting_non_existent_head_at_top);
    TEST(test_return_inside_user_defined_function);

    TEST(test_return_zero_arg_in_each_boundary);
    TEST(test_return_two_arg_each_boundary);
    TEST(test_return_through_deep_nesting);

    TEST(test_return_symbol_metadata);

    printf("All Return tests passed!\n");
    return 0;
}
