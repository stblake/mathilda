#include "iter.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>
#include <string.h>

void test_do() {
    assert_eval_eq("Do[x, 3]", "Null", 0);
    assert_eval_eq("Do[x, {i, 2}]", "Null", 0);

    // Testing variable substitution
    assert_eval_eq("x = 0; Do[x = x + i, {i, 5}]; x", "15", 0);

    // Testing multiple loops
    assert_eval_eq("x = 0; Do[x = x + i * j, {i, 3}, {j, 2}]; x", "18", 0);

    // Testing Control Flow (Break)
    assert_eval_eq("x = 0; Do[If[i == 3, Break[]]; x = x + i, {i, 5}]; x", "3", 0);

    // Testing Control Flow (Continue)
    assert_eval_eq("x = 0; Do[If[i == 3, Continue[]]; x = x + i, {i, 5}]; x", "12", 0);

    // Testing Control Flow (Return)
    assert_eval_eq("Do[If[i == 3, Return[i * 10]], {i, 5}]", "30", 0);

    // Testing infinite loop bounded by Break
    assert_eval_eq("x = 0; i = 1; Do[If[i > 5, Break[]]; x = x + i; i = i + 1, Infinity]; x", "15", 0);
}

void test_for() {
    assert_eval_eq("For[i=0, i<4, i=i+1, x]", "Null", 0);
    assert_eval_eq("x = 0; For[i=0, i<4, i=i+1, x = x + i]; x", "6", 0);
    assert_eval_eq("t=1; For[k=1, k<=5, k=k+1, t=t*k; If[k<2, Continue[]]; t=t+2]; t", "292", 0);
    assert_eval_eq("For[i=1, i<1000, i=i+1, If[i>10, Break[]]]; i", "11", 0);
    assert_eval_eq("For[a < b, 1, 0]", "Null", 0);
}

/* ===================================================================== */
/*  Break / Continue                                                     */
/* ===================================================================== */

/* Break[] exits the nearest enclosing Do/For/While; the loop then yields
 * Null. Continue[] skips the remainder of the body and moves to the next
 * iteration. These cases mirror the Wolfram Language reference transcripts. */
void test_break_in_all_loops() {
    // Do: sum 1+2 then Break when i>2  -> 3
    assert_eval_eq("s=0; Do[If[i>2,Break[]]; s=s+i, {i,10}]; s", "3", 0);
    // For: exit as soon as i>2  -> i left at 3
    assert_eval_eq("For[i=1,i<=10,i++,If[i>2,Break[]]]; i", "3", 0);
    // While: exit as soon as i>2  -> i left at 3
    assert_eval_eq("i=1; While[i<=10, If[i>2,Break[]]; i++]; i", "3", 0);
}

void test_continue_in_all_loops() {
    // Skip even i, sum the odds 1..10  -> 25 in each loop kind.
    assert_eval_eq("r=0; Do[If[EvenQ[i],Continue[]]; r+=i, {i,10}]; r", "25", 0);
    assert_eval_eq("r=0; For[i=1,i<=10,i++, If[EvenQ[i],Continue[]]; r+=i]; r", "25", 0);
    assert_eval_eq("r=0; i=0; While[i<10, i++; If[EvenQ[i],Continue[]]; r+=i]; r", "25", 0);
}

// Continue must still advance a Do arithmetic-progression counter (guards the
// iter.c range-advance branch) -- sum odds in {i,1,6,2} skipping 3 -> 1+5 = 6.
void test_continue_advances_do_range_counter() {
    assert_eval_eq("r=0; Do[If[i==3,Continue[]]; r+=i, {i,1,6,2}]; r", "6", 0);
}

// A loop exited via Break yields Null.
void test_break_loop_returns_null() {
    assert_eval_eq("Do[Break[], {i,5}]", "Null", 0);
    assert_eval_eq("Clear[i]; i=0; While[True, Break[]]", "Null", 0);
}

// Break escapes only the innermost loop, not the enclosing one.
void test_break_is_local_to_innermost_loop() {
    // Inner Do increments s once (j==2 breaks after one add) per outer i.
    assert_eval_eq("s=0; Do[Do[If[j==2,Break[]]; s=s+1, {j,5}], {i,3}]; s", "3", 0);
}

// Break[]/Continue[] with no enclosing loop: emit <head>::nofwd (stderr) and
// return Hold[...] so the marker is rendered inert.
void test_break_continue_out_of_loop() {
    assert_eval_eq("Break[]", "Hold[Break[]]", 0);
    assert_eval_eq("Continue[]", "Hold[Continue[]]", 0);
    // Reached at the end of a CompoundExpression, still uncaught at top level.
    assert_eval_eq("1; Break[]", "Hold[Break[]]", 0);
}

// Wrong arity: Break/Continue take no arguments; the call stays unevaluated
// (an argx message is emitted to stderr).
void test_break_continue_arity() {
    assert_eval_eq("Break[1]", "Break[1]", 0);
    assert_eval_eq("Continue[1, 2]", "Continue[1, 2]", 0);
}

// Both are Protected builtins and cannot be redefined.
void test_break_continue_protected() {
    assert_eval_eq("MemberQ[Attributes[Break], Protected]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[Continue], Protected]", "True", 0);
}

// Docstrings are attached and describe the loop behaviour.
void test_break_continue_docstrings() {
    struct Expr* eb = parse_expression("Information[\"Break\"]");
    struct Expr* rb = evaluate(eb);
    ASSERT(rb != NULL && rb->type == EXPR_STRING);
    ASSERT(strstr(rb->data.string, "Break") != NULL);
    ASSERT(strstr(rb->data.string, "loop") != NULL);
    expr_free(eb);
    expr_free(rb);

    struct Expr* ec = parse_expression("Information[\"Continue\"]");
    struct Expr* rc = evaluate(ec);
    ASSERT(rc != NULL && rc->type == EXPR_STRING);
    ASSERT(strstr(rc->data.string, "Continue") != NULL);
    ASSERT(strstr(rc->data.string, "loop") != NULL);
    expr_free(ec);
    expr_free(rc);
}

/* ===================================================================== */
/*  While                                                                */
/* ===================================================================== */

// While with False test -- body must never execute.
void test_while_false_test_never_executes_body() {
    assert_eval_eq("Clear[x]; x = 0; While[False, x = x + 1]; x", "0", 0);
}

// Basic counted loop: accumulate sum of 1..3.
void test_while_basic_accumulate() {
    assert_eval_eq("Clear[n, s]; n = 1; s = 0; While[n < 4, s = s + n; n = n + 1]; {n, s}",
                   "{4, 6}", 0);
}

// Loop that exits because test immediately fails.
void test_while_test_fails_first_iteration() {
    assert_eval_eq("Clear[n]; n = 10; While[n < 4, n = n + 1]; n", "10", 0);
}

// While[test] -- body implicitly Null; test has a side effect. Manually
// step a counter from within the test using CompoundExpression to emulate
// `While[++n < 4]`.
void test_while_null_body() {
    assert_eval_eq("Clear[n]; n = 1; While[n = n + 1; n < 4]; n", "4", 0);
}

// GCD via While -- canonical Euclid's-algorithm example from the spec.
void test_while_gcd_via_euclid() {
    assert_eval_eq("Clear[a, b, t]; a = 27; b = 6; While[b != 0, t = b; b = Mod[a, b]; a = t]; a",
                   "3", 0);
}

// Break[] must exit While and leave the enclosing counter at the right value.
void test_while_break() {
    assert_eval_eq("Clear[n]; n = 1; While[True, If[n > 10, Break[]]; n = n + 1]; n",
                   "11", 0);
}

// Continue[] must skip the rest of the body and re-enter the test loop.
// Accumulate 1..10 but skip even numbers: expect 1+3+5+7+9 = 25.
void test_while_continue() {
    assert_eval_eq(
        "Clear[n, s]; n = 0; s = 0;"
        " While[n < 10,"
        "   n = n + 1;"
        "   If[Mod[n, 2] == 0, Continue[]];"
        "   s = s + n];"
        " s", "25", 0);
}

// Return[val] inside the body makes While yield val (rather than Null).
void test_while_return_value() {
    assert_eval_eq("Clear[n]; n = 0; While[True, n = n + 1; If[n == 5, Return[n * n]]]",
                   "25", 0);
}

// While returning Null when no Return fires.
void test_while_returns_null() {
    assert_eval_eq("Clear[n]; n = 0; While[n < 3, n = n + 1]",
                   "Null", 0);
}

// Nested While loops: compute the sum of the products i*j for
// 1 <= i < 3 and 1 <= j < 3  --> 1+2+2+4 = 9.
void test_while_nested() {
    assert_eval_eq(
        "Clear[i, j, s]; i = 1; s = 0;"
        " While[i < 3,"
        "   j = 1;"
        "   While[j < 3,"
        "     s = s + i * j;"
        "     j = j + 1];"
        "   i = i + 1];"
        " s", "9", 0);
}

// Break must only escape the innermost While, not the outer loop.
void test_while_break_is_local_to_inner_loop() {
    assert_eval_eq(
        "Clear[i, j, count]; i = 0; count = 0;"
        " While[i < 3,"
        "   i = i + 1;"
        "   j = 0;"
        "   While[True, j = j + 1; If[j >= 2, Break[]]];"
        "   count = count + j];"
        " {i, count}", "{3, 6}", 0);
}

// A False body-result (not Break) does NOT exit the loop --
// the loop only exits when TEST is not True.
void test_while_body_returning_false_does_not_exit() {
    assert_eval_eq("Clear[n]; n = 0; While[n < 5, n = n + 1; False]; n", "5", 0);
}

// While should not accept 0 or more than 2 arguments -- it must remain
// unevaluated (original head preserved) in those cases.
void test_while_wrong_argc_stays_unevaluated() {
    assert_eval_eq("Head[While[]]", "While", 0);
    assert_eval_eq("Head[While[True, 1, 2]]", "While", 0);
}

// Arguments are held: a symbolic test that doesn't evaluate to True simply
// stops the loop -- it must not be evaluated eagerly before the loop starts.
void test_while_holdall_symbolic_test() {
    // `cond` is an unbound symbol so `cond` doesn't become True; loop exits.
    assert_eval_eq("Clear[cond, x]; x = 0; While[cond, x = x + 1]; x", "0", 0);
}

// Test changes state that eventually terminates the loop.
void test_while_test_becomes_false() {
    assert_eval_eq("Clear[n]; n = 0; While[n < 100, n = n + 7]; n", "105", 0);
}

// Return exits only the innermost While (matching Do/For semantics in
// Mathilda): the inner While yields the Return's argument as its value,
// and the outer While continues normally.
void test_while_return_escapes_innermost_loop() {
    assert_eval_eq(
        "Clear[i, j, outer]; i = 0; outer = 0;"
        " While[i < 3,"
        "   j = 0;"
        "   While[True, j = j + 1; If[j == 2, Return[stopped]]];"
        "   outer = outer + 1;"
        "   i = i + 1];"
        " {i, outer}", "{3, 3}", 0);
}

// Docstring is attached and mentions the key behaviour.
void test_while_has_docstring() {
    struct Expr* e = parse_expression("Information[\"While\"]");
    struct Expr* r = evaluate(e);
    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_STRING);
    ASSERT(strstr(r->data.string, "While[test, body]") != NULL);
    ASSERT(strstr(r->data.string, "Break") != NULL);
    ASSERT(strstr(r->data.string, "Continue") != NULL);
    expr_free(e);
    expr_free(r);
}

/* =====================================================================
 *  Iterator termination near the int64 boundary and past the old caps
 *  (regression suite for GitHub issue #52).
 *
 *  The loops used to terminate on a double comparison (`val <= max_val`);
 *  near 2^63 consecutive int64 values collapse to the same double, so Do
 *  never stopped, Sum/Product ran to their term cap, and Table ran to a
 *  1000000 cap and silently truncated. The fix compares the exact running
 *  value against the exact bound (iter_range_continue), raises Table's cap
 *  to match Sum/Product (100000000, decline instead of truncate), and makes
 *  the auto-compiled loops (numloop_do_range) increment overflow-safe.
 * ===================================================================== */

/* Table / Sum / Do terminate correctly when the bounds sit at the very top of
 * the int64 range — the exact case in issue #52. */
void test_iter_int64_boundary_ascending() {
    assert_eval_eq("Length[Table[i, {i, 9223372036854775805, 9223372036854775807}]]", "3", 0);
    assert_eval_eq("Table[i, {i, 9223372036854775805, 9223372036854775807}]",
                   "{9223372036854775805, 9223372036854775806, 9223372036854775807}", 0);
    /* step > 1 landing exactly on the bound */
    assert_eval_eq("Table[i, {i, 9223372036854775803, 9223372036854775807, 2}]",
                   "{9223372036854775803, 9223372036854775805, 9223372036854775807}", 0);
    /* Sum/Do over the same top-of-range span */
    assert_eval_eq("Sum[i, {i, 9223372036854775805, 9223372036854775807}]",
                   "27670116110564327418", 0);
    assert_eval_eq("Sum[1, {i, 9223372036854775805, 9223372036854775807}]", "3", 0);
    assert_eval_eq("Product[1, {i, 9223372036854775805, 9223372036854775807}]", "1", 0);
}

/* Descending near the top, and ascending off the very bottom (INT64_MIN). */
void test_iter_int64_boundary_descending_and_min() {
    assert_eval_eq("Table[i, {i, 9223372036854775807, 9223372036854775805, -1}]",
                   "{9223372036854775807, 9223372036854775806, 9223372036854775805}", 0);
    assert_eval_eq("Table[i, {i, -9223372036854775808, -9223372036854775806}]",
                   "{-9223372036854775808, -9223372036854775807, -9223372036854775806}", 0);
    assert_eval_eq("Length[Table[i, {i, -9223372036854775808, -9223372036854775806}]]", "3", 0);
}

/* The auto-compiled Do path (numloop_do_range) must also terminate at the
 * boundary rather than overflow int64 into a ~2^63-iteration wrap. Bodies that
 * assign a machine-numeric accumulator take that path. */
void test_iter_autocompiled_boundary() {
    /* integer accumulator */
    assert_eval_eq("Module[{s=0}, Do[s=s+1, {i, 9223372036854775805, 9223372036854775807}]; s]",
                   "3", 0);
    /* the summed values promote out of int64 exactly */
    assert_eval_eq("Module[{s=0}, Do[s=s+i, {i, 9223372036854775805, 9223372036854775807}]; s]",
                   "27670116110564327418", 0);
    /* real accumulator (the double-block loop) */
    assert_eval_eq("Module[{s=0.}, Do[s=s+1., {i, 9223372036854775805, 9223372036854775807}]; s]",
                   "3.", 0);
    /* in-place Part-assignment loop (partloop) */
    assert_eval_eq("a=ConstantArray[0,3]; Do[a[[i-9223372036854775804]]=i, "
                   "{i, 9223372036854775805, 9223372036854775807}]; a",
                   "{9223372036854775805, 9223372036854775806, 9223372036854775807}", 0);
}

/* Table must not truncate a legitimate range longer than the old 1000000 cap.
 * `Table[i, {i, 1, 2000000}]` used to come back with 1000001 elements. */
void test_table_no_million_truncation() {
    assert_eval_eq("Length[Table[i, {i, 1, 2000000}]]", "2000000", 0);
    assert_eval_eq("Last[Table[i, {i, 1, 1500000}]]", "1500000", 0);
    assert_eval_eq("Table[i, {i, 999999, 1000002}]", "{999999, 1000000, 1000001, 1000002}", 0);
    /* Do never truncated but must stay correct at that size */
    assert_eval_eq("Module[{s=0}, Do[s=s+1, {i, 1, 2000000}]; s]", "2000000", 0);
    assert_eval_eq("Sum[1, {i, 1, 2000000}]", "2000000", 0);
}

/* A range whose exact element count exceeds the backstop (100000000) returns
 * Table[...] unevaluated — and does so at once (O(1) length check), not after
 * allocating a hundred million elements. */
void test_table_overcap_declines() {
    assert_eval_eq("Head[Table[i, {i, 1, 200000000}]]", "Table", 0);
    /* exactly at a workable size below the cap still evaluates to a list */
    assert_eval_eq("Head[Table[i, {i, 1, 5}]]", "List", 0);
}

/* The exact/inexact semantics of the running value are unchanged by the fix:
 * mixed Integer/Rational stays exact, Real steps stay Real, descending works. */
void test_iter_exactness_preserved() {
    assert_eval_eq("Table[i, {i, 0, 1, 1/3}]", "{0, 1/3, 2/3, 1}", 0);
    assert_eval_eq("Table[i, {i, 0., 1., 0.25}]", "{0., 0.25, 0.5, 0.75, 1.}", 0);
    assert_eval_eq("Table[i, {i, 5, 1, -1}]", "{5, 4, 3, 2, 1}", 0);
    assert_eval_eq("Table[i, {i, 1, 2, 0.5}]", "{1., 1.5, 2.}", 0);
    assert_eval_eq("Sum[1/i, {i, 1, 5}]", "137/60", 0);
    assert_eval_eq("Sum[i, {i, 0, 1, 1/4}]", "5/2", 0);
    assert_eval_eq("Product[i, {i, 1, 6}]", "720", 0);
    /* The exact running value promotes to BigInt on the last advance (807 + 1)
     * while the bound is still int64, exercising the helper's GMP compare and
     * stopping at the bound rather than continuing into the BigInt. */
    assert_eval_eq("Table[i, {i, 9223372036854775806, 9223372036854775807}]",
                   "{9223372036854775806, 9223372036854775807}", 0);
}

/* With/Module/Block must substitute a bound name into a length-1 Table iteration
 * COUNT {n} — it is a count, not a binding. Only {i, ...} binds i as an iterator
 * variable. Regression: With[{n=3}, Table[x, {n}]] used to stay unevaluated. */
void test_scoping_count_iterator() {
    assert_eval_eq("With[{n = 3}, Table[7, {n}]]", "{7, 7, 7}", 0);
    assert_eval_eq("Module[{n = 4}, Table[0, {n}]]", "{0, 0, 0, 0}", 0);
    assert_eval_eq("Block[{n = 2}, Table[5, {n}]]", "{5, 5}", 0);
    /* the count expression is substituted, then evaluated */
    assert_eval_eq("With[{n = 2}, Length[Table[0, {2 n}]]]", "4", 0);
    /* an iterator VARIABLE {i, ...} is still shadowed, not substituted */
    assert_eval_eq("With[{i = 99}, Table[i, {i, 3}]]", "{1, 2, 3}", 0);
    /* the range bound {j, n} still gets the With value */
    assert_eval_eq("With[{n = 3}, Table[j, {j, n}]]", "{1, 2, 3}", 0);
    /* nested counts */
    assert_eval_eq("With[{n = 2}, Table[Table[0, {n}], {n}]]", "{{0, 0}, {0, 0}}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_do);
    TEST(test_for);

    TEST(test_break_in_all_loops);
    TEST(test_continue_in_all_loops);
    TEST(test_continue_advances_do_range_counter);
    TEST(test_break_loop_returns_null);
    TEST(test_break_is_local_to_innermost_loop);
    TEST(test_break_continue_out_of_loop);
    TEST(test_break_continue_arity);
    TEST(test_break_continue_protected);
    TEST(test_break_continue_docstrings);

    TEST(test_while_false_test_never_executes_body);
    TEST(test_while_basic_accumulate);
    TEST(test_while_test_fails_first_iteration);
    TEST(test_while_null_body);
    TEST(test_while_gcd_via_euclid);
    TEST(test_while_break);
    TEST(test_while_continue);
    TEST(test_while_return_value);
    TEST(test_while_returns_null);
    TEST(test_while_nested);
    TEST(test_while_break_is_local_to_inner_loop);
    TEST(test_while_body_returning_false_does_not_exit);
    TEST(test_while_wrong_argc_stays_unevaluated);
    TEST(test_while_holdall_symbolic_test);
    TEST(test_while_test_becomes_false);
    TEST(test_while_return_escapes_innermost_loop);
    TEST(test_while_has_docstring);

    /* Issue #52: int64-boundary termination + past-the-old-cap ranges. */
    TEST(test_iter_int64_boundary_ascending);
    TEST(test_iter_int64_boundary_descending_and_min);
    TEST(test_iter_autocompiled_boundary);
    TEST(test_table_no_million_truncation);
    TEST(test_table_overcap_declines);
    TEST(test_iter_exactness_preserved);
    TEST(test_scoping_count_iterator);

    printf("All iter tests passed!\n");
    symtab_clear();
    return 0;
}
