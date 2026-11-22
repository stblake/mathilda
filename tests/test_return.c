/*
 * test_return.c -- coverage for Return[] / Return[expr] / Return[expr, h]
 *
 * The reference semantics, every expected value below checked against the
 * reference kernel's interactive top level (2026-11-22):
 *
 *   Return[expr] / Return[] is a VALUE that CompoundExpression, While and
 *   For stop at and pass up. Module, Block, With and pure Functions pass
 *   it through as well. Only three things consume it: a definition's rule
 *   boundary (f[x_] := ..., OwnValues, SubValues, UpValues), Do and Scan.
 *   At the top level it stays as it is: Module[{}, Return[1]; 2] at the
 *   prompt is Return[1], and Map[Return, {1}] is {Return[1]}.
 *
 *   Return[expr, h] is a genuine non-local exit to the nearest enclosing
 *   h. It leaves Table, Map, Select and function arguments on the way. The
 *   reference supports it for the heads that are frames on its stack
 *   (Module, Block, Do, While, For, Scan). Mathilda also accepts
 *   Return[v, With], Return[v, Function] and Return[v, f] for a
 *   definition's head f, where the reference says Return::nofunc.
 *
 * Before v0.6695 Mathilda let Module, Block, With, Function, While and For
 * consume a plain Return. A Return inside a With or a While therefore
 * left only that construct and fell through to the rest of the function.
 * In a decision procedure that turned "undecided" into a false answer.
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
    check_eq("Function[x, Return[x]][5]", "Return[5]");
    /* Slot form */
    check_eq("(Return[#] &)[7]", "Return[7]");
    /* Multi-param form */
    check_eq("Function[{x, y}, Return[x + y]][3, 4]", "Return[7]");
    /* Return[] -> Null */
    check_eq("Function[x, Return[]][9]", "Return[]");
}

static void test_return_in_function_skips_following_statements(void) {
    /* The "after" branch of CompoundExpression must not be evaluated
     * once Return has fired. We use a side-effect to verify. */
    check_eq("Clear[r]; r = 0; "
                   "Function[x, Return[x]; r = 999][42]", "Return[42]");
    check_eq("r", "0");
}

static void test_return_in_function_inside_if(void) {
    /* Return inside an If branch: neither the If nor the Function traps it;
     * with no definition around, the marker reaches the top level. */
    check_eq("Function[x, If[x > 0, Return[positive], Return[negative]]][3]",
                   "Return[positive]");
    check_eq("Function[x, If[x > 0, Return[positive], Return[negative]]][-3]",
                   "Return[negative]");
}

/* -----------------------------------------------------------------
 * Module / Block / With boundaries (1-arg form)
 * ----------------------------------------------------------------- */
static void test_return_in_module(void) {
    check_eq("Module[{}, Return[42]]", "Return[42]");
    check_eq("Module[{}, Return[]]", "Return[]");
    check_eq("Module[{x = 10}, Return[x + 1]]", "Return[11]");
}

static void test_return_in_block(void) {
    check_eq("Block[{}, Return[42]]", "Return[42]");
    check_eq("Block[{}, Return[]]", "Return[]");
    check_eq("Block[{x = 10}, Return[x + 1]]", "Return[11]");
}

static void test_return_in_with(void) {
    check_eq("With[{x = 10}, Return[x + 1]]", "Return[11]");
    /* With value substitution still happens before Return fires. */
    check_eq("With[{x = 5}, Return[x*x]]", "Return[25]");
}

static void test_return_in_module_skips_subsequent_statements(void) {
    check_eq("Clear[r]; r = 0; "
                   "Module[{}, Return[7]; r = 999]", "Return[7]");
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
                   "Return[400]");
}

static void test_return_in_while(void) {
    check_eq("Clear[n]; n = 0; While[True, n = n + 1; If[n == 5, Return[n]]]",
                   "Return[5]");
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
    /* The CompoundExpression stops at the Return and the Module passes it
     * on: only a definition (or Do / Scan) strips it. */
    check_eq("Module[{}, 1; Return[2]; 3]", "Return[2]");
    check_eq("Module[{}, 1; Return[]; 3]", "Return[]");
}

/* -----------------------------------------------------------------
 * If / Which / Switch do not trap Return -- it propagates outward.
 * ----------------------------------------------------------------- */
static void test_if_does_not_trap_return(void) {
    /* No enclosing boundary -> Return survives the If. */
    check_eq("If[True, Return[5], Return[6]]", "Return[5]");
    check_eq("If[False, Return[5], Return[6]]", "Return[6]");
    /* Wrap with Module to confirm the Return does pierce the If. */
    check_eq("Module[{}, If[True, Return[5]]; 999]", "Return[5]");
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
    check_eq("Module[{}, Return[{1, 2, 3}]]", "Return[{1, 2, 3}]");
    check_eq("Function[x, Return[{x, x+1}]][7]", "Return[{7, 8}]");
}

static void test_return_value_is_unevaluated_symbolic(void) {
    /* Return[a + b] where a and b are unbound: symbolic propagates. */
    check_eq("Module[{}, Return[a + b]]", "Return[a + b]");
}

static void test_return_arg_is_evaluated_before_unwind(void) {
    /* Return[v]'s payload is fully evaluated before the marker exits.
     * Here, the inner Plus must be reduced to 6. */
    check_eq("Module[{}, Return[1 + 2 + 3]]", "Return[6]");
}

static void test_return_inside_pure_function_inner(void) {
    /* Neither Function consumes Return[x*2]: it reaches the top level. */
    check_eq("Function[y, Function[x, Return[x*2]][y + 1]][4]",
                   "Return[10]");
}

static void test_return_targeting_non_existent_head_at_top(void) {
    /* No Foo boundary anywhere -> marker survives. */
    check_eq("Module[{}, Return[5, Foo]]",
                   "Return[5, Foo]");
}

static void test_return_inside_user_defined_function(void) {
    /* User-level f[x_] := body: the definition itself is a boundary, so the
     * Return yields its value (the reference gives 6). */
    check_eq("Clear[f]; f[x_] := Return[x + 1]; f[5]",
                   "6");
    /* The Module passes the Return on to the definition, which strips it. */
    check_eq("Clear[g]; g[x_] := Module[{}, Return[x + 1]]; g[5]",
                   "6");
}

static void test_return_zero_arg_in_each_boundary(void) {
    check_eq("Function[x, Return[]][1]", "Return[]");
    check_eq("Module[{}, Return[]]", "Return[]");
    check_eq("Block[{}, Return[]]", "Return[]");
    check_eq("With[{x = 1}, Return[]]", "Return[]");
    check_eq("Do[Return[], {3}]", "Null");
    check_eq("For[i = 1, i < 5, i = i + 1, Return[]]", "Return[]");
    check_eq("While[True, Return[]]", "Return[]");
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
    /* Function -> Module -> Do -> For -> While: the While and the For stop
     * at Return[12] and pass it up; the Do consumes it and yields 12, which
     * the Module and the Function return (the reference gives 12). */
    check_eq("Function[a, "
                   "Module[{}, "
                   "Do["
                   "For[i = 1, i < 99, i = i + 1, "
                   "While[True, Return[a*2]]"
                   "], {3}]]][6]",
                   "12");

    /* With no Do or definition around, nothing consumes it. */
    check_eq("Function[a, Module[{}, While[True, Return[a*2]]]][6]",
             "Return[12]");

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

    /* Return[v, Module] skips the While and stops at the Module. */
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

/* A definition is a boundary too: a rule whose evaluated right-hand side is
 * Return[v] yields v, for own, down, sub and up values; only one top-level
 * Return is stripped. Every expected value is the reference's. FeynCalc's
 * cartesianIntegralQ[ex_] := Which[..., Return[True], ...] leaked Return[True],
 * so TrueQ of it was False. */
static void test_return_from_definition(void) {
    check_eq("rd1[x_] := Return[x]; rd1[1]", "1");
    check_eq("rd2[x_] := Which[x > 0, Return[True], True, Return[False]]; {rd2[1], rd2[-1], TrueQ[rd2[1]]}",
             "{True, False, True}");
    check_eq("rd3[x_] := If[x > 0, Return[1]]; {rd3[1], rd3[-1]}", "{1, Null}");
    check_eq("rd4[x_] := (Do[Return[x], {3}]; 7); rd4[1]", "7");
    check_eq("rd5[x_] := rd6[x]; rd6[x_] := Return[x + 1]; rd5[1]", "2");
    check_eq("rd7[x_] := Return[]; rd7[1]", "Null");
    check_eq("rd8[x_][y_] := Return[x y]; rd8[2][3]", "6");
    check_eq("rh9 /: rk9[rh9[x_]] := Return[x]; rk9[rh9[5]]", "5");
    check_eq("rx12 := Return[3]; rx12", "3");
    /* only the top level, and only once */
    check_eq("rd11[x_] := Return[Return[x]]; rd11[2]", "Return[2]");
    check_eq("rd14[x_] := Return[x] + 1; rd14[2]", "1 + Return[2]");
    check_eq("rd15[x_] := {Return[x]}; rd15[2]", "{Return[2]}");
    check_eq("rd16[x_] := Hold[Return[x]]; rd16[1]", "Hold[Return[1]]");
    check_eq("rfac[n_] := If[n <= 1, Return[1], Return[n rfac[n - 1]]]; rfac[20]", "2432902008176640000");
}


/* -----------------------------------------------------------------
 * The reference table (2026-11-22). Each case is a definition, so the
 * Return has a function boundary to leave, and the value shows which
 * construct consumed it. Every expected value is the reference kernel's.
 * Before v0.6695 Mathilda gave 2 / 99 for every row marked (*): With,
 * Block, a nested Module, While, For and a pure Function consumed the
 * Return and the function carried on.
 * ----------------------------------------------------------------- */
static void test_return_reference_table(void) {
    /* CompoundExpression, If, Which, Switch: transparent */
    check_eq("rt1[] := (1; (Return[1]; 3); 2); rt1[]", "1");
    check_eq("rt2[] := Module[{}, If[True, If[True, Return[1]]]; 2]; rt2[]", "1");
    check_eq("rt3[] := Module[{}, Which[True, Return[1]]; 2]; rt3[]", "1");
    check_eq("rt4[] := Module[{}, Switch[1, 1, Return[1]]; 2]; rt4[]", "1");
    /* Module, Block, With: transparent (*) */
    check_eq("rt5[] := (Module[{}, Return[1]]; 2); rt5[]", "1");
    check_eq("rt6[] := Module[{}, Module[{}, If[True, Return[1]]; 3]; 2]; rt6[]", "1");
    check_eq("rt7[] := Module[{}, Block[{}, Return[1]]; 2]; rt7[]", "1");
    check_eq("rt8[] := Module[{}, With[{a = 1}, If[True, Return[a]]]; 2]; rt8[]", "1");
    check_eq("rt9[] := Module[{}, If[True, With[{a = 5}, Return[a]]]; 2]; rt9[]", "5");
    /* While, For: stop and pass up (*) */
    check_eq("rt10[] := Module[{i = 0}, While[i < 3, i++; If[i == 2, Return[i]]]; 99]; rt10[]", "2");
    check_eq("rt11[] := Module[{}, For[i = 0, i < 3, i++, If[i == 1, Return[i]]]; 99]; rt11[]", "1");
    /* Do, Scan: consume */
    check_eq("rt12[x_] := (Do[If[x > 5, Return[a]], {3}]; x); rt12[6]", "6");
    check_eq("rt13[] := Module[{}, Do[Do[If[j == 2, Return[j]], {j, 3}], {2}]; 99]; rt13[]", "99");
    check_eq("rt14[] := Module[{}, Scan[If[# == 2, Return[#]] &, {1, 2, 3}]; 99]; rt14[]", "99");
    check_eq("rt15[] := Module[{}, Scan[If[# == 2, Return[#]] &, {1, 2, 3}]]; rt15[]", "2");
    /* pure Function: transparent (*); Table and Map keep it as a value */
    check_eq("rt16[] := Module[{}, (Return[#]; 5) &[3]; 99]; rt16[]", "3");
    check_eq("rt17[] := Module[{}, Function[x, Return[x]; 5][3]; 99]; rt17[]", "3");
    check_eq("rt18[] := Module[{}, Table[If[i == 2, Return[i], i], {i, 3}]]; rt18[]", "{1, Return[2], 3}");
    check_eq("rt19[] := Module[{}, Map[If[# == 2, Return[#], #] &, {1, 2, 3}]]; rt19[]", "{1, Return[2], 3}");
    check_eq("rt20[] := Module[{}, Select[{1, 2, 3}, (If[# == 2, Return[True]]; False) &]]; rt20[]", "{}");
    /* Catch: transparent; Throw still goes to it */
    check_eq("rt21[] := Module[{}, Catch[If[True, Return[1]]; 3]; 2]; rt21[]", "1");
    check_eq("rt22[x_] := Catch[Do[If[x > 5, Throw[a]], {3}]; x]; rt22[6]", "a");
    /* Check / TimeConstrained / Quiet pass it to the definition */
    check_eq("rt23[] := Check[Module[{}, Return[9]; 1], fail]; rt23[]", "9");
    check_eq("rt24[] := TimeConstrained[Module[{}, Return[9]; 1], 10, fail]; rt24[]", "9");
    check_eq("rt25[] := (Quiet[Module[{}, Return[9]; 1]]; 5); rt25[]", "9");
    /* a holding wrapper's own definition is a boundary */
    check_eq("SetAttributes[rtw, HoldAll]; rtw[b_] := Block[{}, b]; "
             "rt26[] := (rtw[Module[{}, Return[6]; 7]]; 8); rt26[]", "8");
    /* memoisation stores the Return, which the stored rule strips again */
    check_eq("rtm[x_] := rtm[x] = Module[{}, If[x > 0, Return[x]]; 0]; {rtm[3], rtm[3]}", "{3, 3}");
    /* Set passes it: the CompoundExpression stops at the assignment */
    check_eq("rt27[] := Module[{y}, y = Module[{}, Return[1]]; {y}]; rt27[]", "1");
}

/* Return[v, h] leaves Table, Map, Select, Sum and function arguments on its
 * way to h (the reference unwinds through them). Before v0.6695 the marker
 * stayed in the list and the function went on with it:
 * f[] := Module[{}, Table[If[i == 2, Return[i, Module], i], {i, 3}]; 99]
 * was 99 and {..., Return[None, Module], ...} reached callers. */
static void test_return_targeted_unwinds_builtins(void) {
    check_eq("ru1[] := Module[{}, Table[If[i == 2, Return[i, Module], i], {i, 3}]; 99]; ru1[]", "2");
    check_eq("ru2[] := Module[{}, Map[If[# == 2, Return[#, Module], #] &, {1, 2, 3}]; 99]; ru2[]", "2");
    check_eq("ru3[] := Module[{}, Select[{1, 2, 3}, If[# == 2, Return[#, Module], True] &]; 99]; ru3[]", "2");
    check_eq("ru4[] := Module[{}, Sum[If[k == 2, Return[k, Module], k], {k, 3}]; 99]; ru4[]", "2");
    check_eq("ru5[] := Module[{}, h5[Return[1, Module]]; 2]; ru5[]", "1");
    check_eq("ru6[] := Module[{}, Hold[Return[1, Module]]; 2]; ru6[]", "2");
    check_eq("ru7[] := Module[{}, Scan[Return[#, Module] &, {1, 2}]; 99]; ru7[]", "1");
    check_eq("ru8[] := Module[{}, Do[Do[If[j == 2, Return[j, Do]], {j, 3}]; k = 5, {2}]; 99]; ru8[]", "99");
    check_eq("ru9[] := Module[{i = 0}, While[i < 3, i++; Do[Return[i, While], {2}]]]; ru9[]", "1");
    check_eq("ru10[] := Module[{}, Block[{}, Do[Return[1, Block], {3}]; 3]; 2]; ru10[]", "2");
    /* the shape of ParallelMixed.m's ResInfAll: an early exit from a Table */
    check_eq("ru11[s_] := Module[{}, Table[If[j == 1, Return[None, Module]]; j, {j, 0, s}]]; "
             "{ru11[0], ru11[3]}", "{{0}, None}");
    /* a Return with three arguments is not a Return: it stays as it is */
    check_eq("ru12[] := Module[{}, Return[1, Module, 3]; 2]; ru12[]", "Return[1, Module, 3]");
}

/* Do[body, s1, s2, ...] is ONE loop in the reference: Break[] and Return
 * leave all of it, Continue[] advances the innermost iterator. Mathilda
 * runs it as nested levels, and before v0.6695 the inner level consumed
 * them. Break left one level only, and a Return's value was lost (Null). */
static void test_return_multi_iterator_do(void) {
    check_eq("rm1[] := Module[{}, Do[Return[{i, j}], {i, 2}, {j, 2}]]; rm1[]", "{1, 1}");
    check_eq("rm2[] := Module[{}, Do[If[j == 2, Return[{i, j}, Do]], {i, 2}, {j, 2}]]; rm2[]", "{1, 2}");
    check_eq("rm3[] := Module[{}, Do[Return[1], {i, 2}, {j, 2}]; 99]; rm3[]", "99");
    check_eq("r = {}; Do[If[j == 2, Break[]]; AppendTo[r, {i, j}], {i, 2}, {j, 3}]; r", "{{1, 1}}");
    check_eq("r = {}; Do[If[j == 2, Continue[]]; AppendTo[r, {i, j}], {i, 2}, {j, 3}]; r",
             "{{1, 1}, {1, 3}, {2, 1}, {2, 3}}");
    check_eq("r = {}; Do[If[k == 2, Break[]]; AppendTo[r, {i, j, k}], {i, 2}, {j, 2}, {k, 2}]; r",
             "{{1, 1, 1}}");
    check_eq("s = 0; Do[s += i j, {i, 3}, {j, 3}]; s", "36");
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
    TEST(test_return_from_definition);

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
    TEST(test_return_reference_table);
    TEST(test_return_targeted_unwinds_builtins);
    TEST(test_return_multi_iterator_do);

    printf("All Return tests passed!\n");
    return 0;
}
