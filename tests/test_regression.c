#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "match.h"
#include "core.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_regression_flat_orderless_eval() {
    symtab_init();
    core_init();
    
    // a + b + a -> Plus[b, Times[2, a]] or Plus[Times[2, a], b]
    Expr* e1 = parse_expression("Plus[a, b, a]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_FUNCTION);
    ASSERT(res1->data.function.arg_count == 2);
    // Since Plus is Orderless, the result depends on comparison. 
    // Times[2, a] vs b. 
    // Let's just check the arguments exist.
    bool found_times = false;
    bool found_b = false;
    for (size_t i = 0; i < 2; i++) {
        Expr* arg = res1->data.function.args[i];
        if (arg->type == EXPR_SYMBOL && strcmp(arg->data.symbol.name, "b") == 0) found_b = true;
        if (arg->type == EXPR_FUNCTION && strcmp(arg->data.function.head->data.symbol.name, "Times") == 0) found_times = true;
    }
    ASSERT(found_b && found_times);
    expr_free(e1); expr_free(res1);
}

void test_regression_set_evaluation() {
    symtab_init();
    core_init();
    
    // 1. Set own value
    Expr* e1 = parse_expression("Set[xx, 100]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_INTEGER && res1->data.integer == 100);
    expr_free(e1); expr_free(res1);
    
    // 2. Evaluate own value
    Expr* e2 = parse_expression("xx");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_INTEGER && res2->data.integer == 100);
    expr_free(e2); expr_free(res2);
    
    // 3. Set down value
    Expr* e3 = parse_expression("SetDelayed[f[Pattern[y, Blank[]]], Plus[y, 5]]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol.name, "Null") == 0);
    expr_free(e3); expr_free(res3);
    
    // 4. Evaluate down value with evaluated argument
    // Now that Plus is implemented, f[xx] -> 100 + 5 -> 105
    Expr* e4 = parse_expression("f[xx]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4->type == EXPR_INTEGER && res4->data.integer == 105);
    expr_free(e4); expr_free(res4);
    }

    void test_regression_infinite_eval() {
    symtab_init();
    core_init();

    // Set is HoldFirst, so the LHS is not evaluated before the assignment
    // installs a DownValue on the head symbol. Use an unprotected symbol
    // (ufoo) so the Protected check does not reject the assignment.
    Expr* e1 = parse_expression("Set[ufoo[x], 50]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_INTEGER && res1->data.integer == 50);
    expr_free(e1); expr_free(res1);

    // ufoo[x] should now evaluate to 50 via the freshly installed DownValue.
    Expr* e2 = parse_expression("ufoo[x]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_INTEGER && res2->data.integer == 50);
    expr_free(e2); expr_free(res2);
}

void test_regression_nested_replace() {
    symtab_init();
    core_init();
    // Test that bindings deep in an expression tree are properly replaced.
    Expr* e1 = parse_expression("SetDelayed[g[Pattern[a, Blank[]], Pattern[b, Blank[]]], List[b, a, b]]");
    Expr* res1 = evaluate(e1);
    expr_free(e1); expr_free(res1);
    
    Expr* e2 = parse_expression("g[10, 20]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_FUNCTION);
    ASSERT(strcmp(res2->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(res2->data.function.arg_count == 3);
    ASSERT(res2->data.function.args[0]->data.integer == 20);
    ASSERT(res2->data.function.args[1]->data.integer == 10);
    ASSERT(res2->data.function.args[2]->data.integer == 20);
    expr_free(e2); expr_free(res2);
}

void test_regression_clear() {
    symtab_init();
    core_init();
    
    // Set x = 42
    Expr* e1 = parse_expression("Set[x, 42]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_INTEGER && res1->data.integer == 42);
    
    // Evaluate x
    Expr* e2 = parse_expression("x");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_INTEGER && res2->data.integer == 42);
    
    // Clear[x]
    Expr* e3 = parse_expression("Clear[x]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol.name, "Null") == 0);
    
    // Evaluate x again (should be unbound, so evaluates to symbol 'x')
    Expr* e4 = parse_expression("x");
    Expr* res4 = evaluate(e4);
    ASSERT(res4->type == EXPR_SYMBOL && strcmp(res4->data.symbol.name, "x") == 0);
    
    // SetDelayed f[x_] := x + 1
    Expr* e5 = parse_expression("SetDelayed[f[Pattern[y, Blank[]]], Plus[y, 1]]");
    Expr* res5 = evaluate(e5);
    
    // Evaluate f[10]
    Expr* e6 = parse_expression("f[10]");
    Expr* res6 = evaluate(e6);
    // Now that Plus is implemented, f[10] -> 10 + 1 -> 11
    ASSERT(res6->type == EXPR_INTEGER && res6->data.integer == 11);
    
    // Clear[f]
    Expr* e7 = parse_expression("Clear[f]");
    Expr* res7 = evaluate(e7);
    
    // Evaluate f[10] again (should be unbound, so evaluates to f[10])
    Expr* e8 = parse_expression("f[10]");
    Expr* res8 = evaluate(e8);
    ASSERT(res8->type == EXPR_FUNCTION && strcmp(res8->data.function.head->data.symbol.name, "f") == 0);
    
    expr_free(e1); expr_free(res1);
    expr_free(e2); expr_free(res2);
    expr_free(e3); expr_free(res3);
    expr_free(e4); expr_free(res4);
    expr_free(e5); expr_free(res5);
    expr_free(e6); expr_free(res6);
    expr_free(e7); expr_free(res7);
    expr_free(e8); expr_free(res8);
}

void test_regression_condition_rhs_setdelayed() {
    symtab_init();
    core_init();

    /* Condition on the RHS of SetDelayed should be moved to the LHS pattern.
     * f[x_] := body /; test  is equivalent to  f[x_] /; test := body */
    assert_eval_eq("f[x_] := ppp[x] /; x > 0", "Null", 0);
    assert_eval_eq("f[1]", "ppp[1]", 0);
    assert_eval_eq("f[-1]", "f[-1]", 0);

    /* Multi-argument case */
    assert_eval_eq("h[x_, y_] := x + y /; x > y", "Null", 0);
    assert_eval_eq("h[5, 3]", "8", 0);
    assert_eval_eq("h[2, 7]", "h[2, 7]", 0);

    /* Condition on the LHS should still work */
    assert_eval_eq("g[x_] /; x > 0 := x^2", "Null", 0);
    assert_eval_eq("g[3]", "9", 0);
    assert_eval_eq("g[-2]", "g[-2]", 0);
}

int main() {
    printf("Running extensive regression tests...\n");
    TEST(test_regression_flat_orderless_eval);
    TEST(test_regression_set_evaluation);
    TEST(test_regression_infinite_eval);
    TEST(test_regression_nested_replace);
    TEST(test_regression_clear);
    TEST(test_regression_condition_rhs_setdelayed);
    printf("All regression tests passed!\n");
    return 0;
}