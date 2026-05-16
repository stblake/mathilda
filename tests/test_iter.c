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

int main() {
    symtab_init();
    core_init();

    TEST(test_do);
    TEST(test_for);

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

    printf("All iter tests passed!\n");
    symtab_clear();
    return 0;
}
