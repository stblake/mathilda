#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_hold() {
    // Hold[1] should evaluate to Hold[1]
    Expr* e = parse_expression("Hold[1]");
    Expr* res = evaluate(e);
    
    ASSERT(res != NULL);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.head->data.symbol, "Hold") == 0);
    
    expr_free(e);
    expr_free(res);
}

void test_orderless() {
    // Plus[b, a, c] should evaluate to Plus[a, b, c]
    Expr* e = parse_expression("Plus[b, a, c]");
    Expr* res = evaluate(e);
    
    ASSERT(res != NULL);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count == 3);
    
    ASSERT(res->data.function.args[0]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[0]->data.symbol, "a") == 0);
    
    ASSERT(res->data.function.args[1]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[1]->data.symbol, "b") == 0);
    
    ASSERT(res->data.function.args[2]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[2]->data.symbol, "c") == 0);
    
    expr_free(e);
    expr_free(res);
}

void test_flat() {
    // Plus[Plus[a, b], c] should evaluate to Plus[a, b, c]
    Expr* e = parse_expression("Plus[Plus[a, b], c]");
    Expr* res = evaluate(e);
    
    ASSERT(res != NULL);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count == 3);
    
    ASSERT(res->data.function.args[0]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[0]->data.symbol, "a") == 0);
    
    ASSERT(res->data.function.args[1]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[1]->data.symbol, "b") == 0);
    
    ASSERT(res->data.function.args[2]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[2]->data.symbol, "c") == 0);
    
    expr_free(e);
    expr_free(res);

    // Test SetAttributes[f, Flat] and f[f[x]] -> f[x]
    Expr* e_attrs = parse_expression("SetAttributes[f, Flat]");
    Expr* res_attrs = evaluate(e_attrs);
    expr_free(e_attrs);
    expr_free(res_attrs);

    Expr* e1 = parse_expression("f[f[x]]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_FUNCTION && res1->data.function.arg_count == 1);
    ASSERT(strcmp(res1->data.function.head->data.symbol, "f") == 0);
    ASSERT(res1->data.function.args[0]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res1->data.function.args[0]->data.symbol, "x") == 0);
    expr_free(e1);
    expr_free(res1);

    Expr* e2 = parse_expression("f[f[x], f[y]]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 2);
    ASSERT(strcmp(res2->data.function.head->data.symbol, "f") == 0);
    expr_free(e2);
    expr_free(res2);
}

void test_flat_and_orderless() {
    // Plus[c, Plus[b, a]] should evaluate to Plus[a, b, c]
    Expr* e = parse_expression("Plus[c, Plus[b, a]]");
    Expr* res = evaluate(e);
    
    ASSERT(res != NULL);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count == 3);
    
    ASSERT(res->data.function.args[0]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[0]->data.symbol, "a") == 0);
    
    ASSERT(res->data.function.args[1]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[1]->data.symbol, "b") == 0);
    
    ASSERT(res->data.function.args[2]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[2]->data.symbol, "c") == 0);
    
    expr_free(e);
    expr_free(res);
}

void test_set_and_ownvalues() {
    // Test Set[x, 42] and evaluating x
    Expr* e1 = parse_expression("Set[x, 42]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 != NULL && res1->type == EXPR_INTEGER && res1->data.integer == 42);
    
    Expr* e2 = parse_expression("x");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 != NULL && res2->type == EXPR_INTEGER && res2->data.integer == 42);
    
    expr_free(e1); expr_free(res1);
    expr_free(e2); expr_free(res2);
}

void test_setdelayed_and_downvalues() {
    // Test SetDelayed[f[Pattern[y, Blank[]]], Plus[y, 1]]
    Expr* e1 = parse_expression("SetDelayed[f[Pattern[y, Blank[]]], Plus[y, 1]]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 != NULL && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "Null") == 0);
    
    // Evaluate f[10] -> 11
    Expr* e2 = parse_expression("f[10]");
    Expr* res2 = evaluate(e2);
    // Now that Plus is implemented, f[10] -> 10 + 1 -> 11
    ASSERT(res2 != NULL && res2->type == EXPR_INTEGER && res2->data.integer == 11);
    
    expr_free(e1); expr_free(res1);
    expr_free(e2); expr_free(res2);
}

void test_length() {
    Expr* e1 = parse_expression("Length[{1, 2, 3}]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_INTEGER && res1->data.integer == 3);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Length[5]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_INTEGER && res2->data.integer == 0);
    expr_free(e2); expr_free(res2);
}

void test_dimensions() {
    Expr* e1 = parse_expression("Dimensions[{{1, 2}, {3, 4}}]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION);
    ASSERT(strcmp(res1->data.function.head->data.symbol, "List") == 0);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT(res1->data.function.args[0]->type == EXPR_INTEGER && res1->data.function.args[0]->data.integer == 2);
    ASSERT(res1->data.function.args[1]->type == EXPR_INTEGER && res1->data.function.args[1]->data.integer == 2);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Dimensions[{{1, 2}, {3}}]"); // Ragged
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION);
    ASSERT(res2->data.function.arg_count == 1);
    ASSERT(res2->data.function.args[0]->data.integer == 2);
    expr_free(e2); expr_free(res2);

    /* Two-argument form: limit the result to n levels. */
    assert_eval_eq("Dimensions[{{{{a, b}}}}, 2]", "{1, 1}", 0);
    assert_eval_eq("Dimensions[{{{{a, b}}}}, 4]", "{1, 1, 1, 2}", 0);
    assert_eval_eq("Dimensions[{{1, 2}, {3, 4}}, 1]", "{2}", 0);
    assert_eval_eq("Dimensions[{{1, 2}, {3, 4}}, 5]", "{2, 2}", 0);
    assert_eval_eq("Dimensions[{1, 2, 3}, 0]", "{}", 0);
    assert_eval_eq("Dimensions[5, 3]", "{}", 0);

    /* Edge cases from the WL spec. */
    assert_eval_eq("Dimensions[1]", "{}", 0);
    assert_eval_eq("Dimensions[{1, 2}]", "{2}", 0);
    assert_eval_eq("Dimensions[{{1, 2}}]", "{1, 2}", 0);
    assert_eval_eq("Dimensions[{}]", "{0}", 0);
}

void test_part() {
    Expr* e1 = parse_expression("{a, b, c}[[2]]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "b") == 0);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("{a, b, c}[[-1]]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol, "c") == 0);
    expr_free(e2); expr_free(res2);

    Expr* e3 = parse_expression("{a, b, c}[[0]]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol, "List") == 0);
    expr_free(e3); expr_free(res3);
}

void test_builtin_head() {
    Expr* e1 = parse_expression("Head[g[x]]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "g") == 0);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Head[123]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol, "Integer") == 0);
    expr_free(e2); expr_free(res2);

    /* Atomic heads. Use fresh symbols since earlier tests bind x, y, etc. */
    assert_eval_eq("Head[3.14]", "Real", 0);
    assert_eval_eq("Head[headTestSym1]", "Symbol", 0);
    assert_eval_eq("Head[\"foo\"]", "String", 0);
    assert_eval_eq("Head[headTestSym1 + headTestSym2]", "Plus", 0);
    assert_eval_eq("Head[{headTestSym1, headTestSym2, headTestSym3}]", "List", 0);

    /* Two-argument form: Head[expr, h] -> h[Head[expr]]. */
    assert_eval_eq("Head[headTestSym1 + headTestSym2, headWrap]", "headWrap[Plus]", 0);
    assert_eval_eq("Head[{1, 2}, headWrap]", "headWrap[List]", 0);
    /* The wrapper is itself evaluated: Identity[x] -> x. */
    assert_eval_eq("Head[headTestSym1 + headTestSym2, Identity]", "Plus", 0);
}

void test_first_last_most_rest() {
    Expr* e1 = parse_expression("First[{a, b, c}]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "a") == 0);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Last[{a, b, c}]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol, "c") == 0);
    expr_free(e2); expr_free(res2);

    Expr* e3 = parse_expression("Most[{a, b, c}]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && res3->data.function.arg_count == 2);
    ASSERT(strcmp(res3->data.function.args[0]->data.symbol, "a") == 0);
    ASSERT(strcmp(res3->data.function.args[1]->data.symbol, "b") == 0);
    expr_free(e3); expr_free(res3);

    Expr* e4 = parse_expression("Rest[{a, b, c}]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && res4->data.function.arg_count == 2);
    ASSERT(strcmp(res4->data.function.args[0]->data.symbol, "b") == 0);
    ASSERT(strcmp(res4->data.function.args[1]->data.symbol, "c") == 0);
    expr_free(e4); expr_free(res4);
}

void test_append_prepend() {
    Expr* e1 = parse_expression("Append[{1, 2}, 3]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && res1->data.function.arg_count == 3);
    ASSERT(res1->data.function.args[2]->type == EXPR_INTEGER && res1->data.function.args[2]->data.integer == 3);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Prepend[{1, 2}, 0]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 3);
    ASSERT(res2->data.function.args[0]->type == EXPR_INTEGER && res2->data.function.args[0]->data.integer == 0);
    expr_free(e2); expr_free(res2);
}

void test_append_prepend_to() {
    // xlist = {1, 2}
    Expr* e1 = parse_expression("xlist = {1, 2}");
    Expr* res1 = evaluate(e1);
    expr_free(e1); expr_free(res1);

    // AppendTo[xlist, 3]
    Expr* e2 = parse_expression("AppendTo[xlist, 3]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 3);
    ASSERT(res2->data.function.args[2]->data.integer == 3);
    expr_free(e2); expr_free(res2);

    // Check xlist
    Expr* e3 = parse_expression("xlist");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && res3->data.function.arg_count == 3);
    expr_free(e3); expr_free(res3);

    // PrependTo[xlist, 0]
    Expr* e4 = parse_expression("PrependTo[xlist, 0]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && res4->data.function.arg_count == 4);
    ASSERT(res4->data.function.args[0]->data.integer == 0);
    expr_free(e4); expr_free(res4);
}

void test_attributes() {
    Expr* e1 = parse_expression("Attributes[Plus]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol, "List") == 0);
    // Should contain Flat, Listable, NumericFunction, OneIdentity, Orderless
    bool found_flat = false;
    for (size_t i = 0; i < res1->data.function.arg_count; i++) {
        if (strcmp(res1->data.function.args[i]->data.symbol, "Flat") == 0) found_flat = true;
    }
    ASSERT(found_flat);
    expr_free(e1); expr_free(res1);
}

void test_listable() {
    // Plus[{1, 2}, {3, 4}] -> {4, 6}
    Expr* e1 = parse_expression("Plus[{1, 2}, {3, 4}]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol, "List") == 0);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT(res1->data.function.args[0]->type == EXPR_INTEGER && res1->data.function.args[0]->data.integer == 4);
    ASSERT(res1->data.function.args[1]->type == EXPR_INTEGER && res1->data.function.args[1]->data.integer == 6);
    expr_free(e1); expr_free(res1);

    // {a,b,c} + 1 -> {Plus[1, a], Plus[1, b], Plus[1, c]}
    Expr* e2 = parse_expression("{a, b, c} + 1");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 3);
    // Elements are sorted because Plus is Orderless: Plus[1, a]
    ASSERT(res2->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res2->data.function.args[0]->data.function.head->data.symbol, "Plus") == 0);
    expr_free(e2); expr_free(res2);

    // {a,b,c} + {1,2,3} -> {Plus[1, a], Plus[2, b], Plus[3, c]}
    Expr* e3 = parse_expression("{a, b, c} + {1, 2, 3}");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && res3->data.function.arg_count == 3);
    expr_free(e3); expr_free(res3);

    // {{a,b},c} + 1 -> {{Plus[1, a], Plus[1, b]}, Plus[1, c]}
    Expr* e4 = parse_expression("{{a, b}, c} + 1");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && res4->data.function.arg_count == 2);
    ASSERT(res4->data.function.args[0]->type == EXPR_FUNCTION && res4->data.function.args[0]->data.function.arg_count == 2);
    expr_free(e4); expr_free(res4);

    // {a,b,c} * 2 -> {2a, 2b, 2c}
    Expr* e5 = parse_expression("{a, b, c} * 2");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && res5->data.function.arg_count == 3);
    ASSERT(res5->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res5->data.function.args[0]->data.function.head->data.symbol, "Times") == 0);
    expr_free(e5); expr_free(res5);

    // Test Mod
    Expr* e6 = parse_expression("Mod[{7, 8, 9}, 3]");
    Expr* res6 = evaluate(e6);
    char* s6_full = expr_to_string_fullform(res6);
    ASSERT(strcmp(s6_full, "List[1, 2, 0]") == 0);
    expr_free(e6);
    expr_free(res6);
    free(s6_full);

    // Test Quotient
    Expr* e7 = parse_expression("Quotient[{10, 11, 12}, 3]");
    Expr* res7 = evaluate(e7);
    char* s7_full = expr_to_string_fullform(res7);
    ASSERT(strcmp(s7_full, "List[3, 3, 4]") == 0);
    expr_free(e7);
    expr_free(res7);
    free(s7_full);

    // Test QuotientRemainder
    Expr* e8 = parse_expression("QuotientRemainder[{10, 11, 12}, 3]");
    Expr* res8 = evaluate(e8);
    char* s8_full = expr_to_string_fullform(res8);
    ASSERT(strcmp(s8_full, "List[List[3, 1], List[3, 2], List[4, 0]]") == 0);
    expr_free(e8);
    expr_free(res8);
    free(s8_full);
}

void test_oneidentity() {
    // Plus[unbound] -> unbound
    Expr* e1 = parse_expression("Plus[unbound]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "unbound") == 0);
    expr_free(e1); expr_free(res1);
}

void test_list_assignment() {
    // {a, b} = {1, 2}
    Expr* e1 = parse_expression("{a, b} = {1, 2}");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && res1->data.function.arg_count == 2);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("a");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_INTEGER && res2->data.integer == 1);
    expr_free(e2); expr_free(res2);

    Expr* e3 = parse_expression("b");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_INTEGER && res3->data.integer == 2);
    expr_free(e3); expr_free(res3);
}

void test_plus() {
    // Clear x to avoid interference from previous tests
    symtab_clear_symbol("x");

    // Plus[1, 2] -> 3
    Expr* e1 = parse_expression("1 + 2");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_INTEGER && res1->data.integer == 3);
    expr_free(e1); expr_free(res1);

    // Plus[1, 2.5] -> 3.5
    Expr* e2 = parse_expression("1 + 2.5");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_REAL && res2->data.real == 3.5);
    expr_free(e2); expr_free(res2);

    // Plus[x, x] -> 2*x
    Expr* e3 = parse_expression("x + x");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION);
    ASSERT(strcmp(res3->data.function.head->data.symbol, "Times") == 0);
    ASSERT(res3->data.function.args[0]->data.integer == 2);
    expr_free(e3); expr_free(res3);

    // Plus[2*x^2, 3*x^2] -> 5*x^2
    // Note: Power is symbolic for now
    Expr* e4 = parse_expression("2*x^2 + 3*x^2");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION);
    ASSERT(strcmp(res4->data.function.head->data.symbol, "Times") == 0);
    ASSERT(res4->data.function.args[0]->data.integer == 5);
    expr_free(e4); expr_free(res4);

    // Plus[1, 2*x, x^2, 4*x] -> Plus[1, 6*x, x^2]
    Expr* e6 = parse_expression("1 + 2*x + x^2 + 4*x");
    Expr* res6 = evaluate(e6);
    ASSERT(res6 && res6->type == EXPR_FUNCTION);
    ASSERT(strcmp(res6->data.function.head->data.symbol, "Plus") == 0);
    ASSERT(res6->data.function.arg_count == 3);
    // Check canonical order: 1, 6*x, x^2
    ASSERT(res6->data.function.args[0]->type == EXPR_INTEGER && res6->data.function.args[0]->data.integer == 1);
    ASSERT(res6->data.function.args[1]->type == EXPR_FUNCTION && strcmp(res6->data.function.args[1]->data.function.head->data.symbol, "Times") == 0);
    ASSERT(res6->data.function.args[2]->type == EXPR_FUNCTION && strcmp(res6->data.function.args[2]->data.function.head->data.symbol, "Power") == 0);
    expr_free(e6); expr_free(res6);

    // Plus[1, x, -1] -> x
    Expr* e5 = parse_expression("1 + x - 1");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_SYMBOL && strcmp(res5->data.symbol, "x") == 0);
    expr_free(e5); expr_free(res5);
}

void test_divide() {
    // Divide[1, 2] -> Rational[1, 2]
    Expr* e1 = parse_expression("Divide[1, 2]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol, "Rational") == 0);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT(res1->data.function.args[0]->data.integer == 1);
    ASSERT(res1->data.function.args[1]->data.integer == 2);
    expr_free(e1); expr_free(res1);

    // Divide[1., 2] -> 0.5
    Expr* e2 = parse_expression("Divide[1., 2]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_REAL && res2->data.real == 0.5);
    expr_free(e2); expr_free(res2);

    // Divide[77, 11] -> 7
    Expr* e3 = parse_expression("Divide[77, 11]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_INTEGER && res3->data.integer == 7);
    expr_free(e3); expr_free(res3);

    // Divide[105, 55] -> Rational[21, 11]
    Expr* e4 = parse_expression("Divide[105, 55]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && strcmp(res4->data.function.head->data.symbol, "Rational") == 0);
    ASSERT(res4->data.function.args[0]->data.integer == 21);
    ASSERT(res4->data.function.args[1]->data.integer == 11);
    expr_free(e4); expr_free(res4);

    // Divide[x, y] -> Times[x, Power[y, -1]]
    Expr* e5 = parse_expression("Divide[x, y]");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && strcmp(res5->data.function.head->data.symbol, "Times") == 0);
    ASSERT(res5->data.function.arg_count == 2);
    ASSERT(strcmp(res5->data.function.args[1]->data.function.head->data.symbol, "Power") == 0);
    expr_free(e5); expr_free(res5);
}

void test_rational_arithmetic() {
    // 1 + 1/2 -> 3/2
    Expr* e1 = parse_expression("1 + 1/2");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol, "Rational") == 0);
    ASSERT(res1->data.function.args[0]->data.integer == 3);
    ASSERT(res1->data.function.args[1]->data.integer == 2);
    expr_free(e1); expr_free(res1);

    // 3/2 * 1/2 -> 3/4
    Expr* e2 = parse_expression("3/2 * 1/2");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && strcmp(res2->data.function.head->data.symbol, "Rational") == 0);
    ASSERT(res2->data.function.args[0]->data.integer == 3);
    ASSERT(res2->data.function.args[1]->data.integer == 4);
    expr_free(e2); expr_free(res2);

    // 1/2 + 1/2 -> 1
    Expr* e3 = parse_expression("1/2 + 1/2");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_INTEGER && res3->data.integer == 1);
    expr_free(e3); expr_free(res3);

    // 1/2 + 0.5 -> 1.0
    Expr* e4 = parse_expression("1/2 + 0.5");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_REAL && res4->data.real == 1.0);
    expr_free(e4); expr_free(res4);
}

void test_power() {
    symtab_clear_symbol("a");
    symtab_clear_symbol("b");
    symtab_clear_symbol("x");
    symtab_clear_symbol("y");
    symtab_clear_symbol("z");

    // Power[2, 3] -> 8
    Expr* e1 = parse_expression("2^3");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_INTEGER && res1->data.integer == 8);
    expr_free(e1); expr_free(res1);

    // Power[4, 1/2] -> 2
    Expr* e2 = parse_expression("4^(1/2)");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_INTEGER && res2->data.integer == 2);
    expr_free(e2); expr_free(res2);

    // Power[8, 1/3] -> 2
    Expr* e3 = parse_expression("8^(1/3)");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_INTEGER && res3->data.integer == 2);
    expr_free(e3); expr_free(res3);

    // (a*b)^2 -> a^2 * b^2
    Expr* e4 = parse_expression("(a*b)^2");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && strcmp(res4->data.function.head->data.symbol, "Times") == 0);
    // Elements should be Power[a, 2] and Power[b, 2]
    expr_free(e4); expr_free(res4);

    // (a^2)^3 -> a^6
    Expr* e5 = parse_expression("(a^2)^3");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && strcmp(res5->data.function.head->data.symbol, "Power") == 0);
    ASSERT(res5->data.function.args[1]->type == EXPR_INTEGER && res5->data.function.args[1]->data.integer == 6);
    expr_free(e5); expr_free(res5);

    // (-1)^(1/2) -> Complex[0, 1]
    Expr* e6 = parse_expression("(-1)^(1/2)");
    Expr* res6 = evaluate(e6);
    ASSERT(res6 && res6->type == EXPR_FUNCTION && strcmp(res6->data.function.head->data.symbol, "Complex") == 0);
    ASSERT(res6->data.function.args[1]->data.integer == 1);
    expr_free(e6); expr_free(res6);
}

void test_sqrt() {
    // Sqrt[4] -> 2
    Expr* e1 = parse_expression("Sqrt[4]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_INTEGER && res1->data.integer == 2);
    expr_free(e1); expr_free(res1);

    // Sqrt[2] -> Power[2, 1/2]
    Expr* e2 = parse_expression("Sqrt[2]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && strcmp(res2->data.function.head->data.symbol, "Power") == 0);
    expr_free(e2); expr_free(res2);

    // Sqrt[{1, 4, 9}] -> {1, 2, 3}
    Expr* e3 = parse_expression("Sqrt[{1, 4, 9}]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && strcmp(res3->data.function.head->data.symbol, "List") == 0);
    ASSERT(res3->data.function.arg_count == 3);
    ASSERT(res3->data.function.args[0]->data.integer == 1);
    ASSERT(res3->data.function.args[1]->data.integer == 2);
    ASSERT(res3->data.function.args[2]->data.integer == 3);
    expr_free(e3); expr_free(res3);

    // Sqrt[x^2] -> Power[Power[x, 2], 1/2]
    Expr* e4 = parse_expression("Sqrt[x^2]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && strcmp(res4->data.function.head->data.symbol, "Power") == 0);
    // Should NOT automatically simplify to x
    expr_free(e4); expr_free(res4);

    // Sqrt[45] -> 3*Sqrt[5]
    Expr* e5 = parse_expression("Sqrt[45]");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && strcmp(res5->data.function.head->data.symbol, "Times") == 0);
    ASSERT(res5->data.function.arg_count == 2);
    ASSERT(res5->data.function.args[0]->type == EXPR_INTEGER && res5->data.function.args[0]->data.integer == 3);
    expr_free(e5); expr_free(res5);

    // Sqrt[-9] -> 3*I
    Expr* e6 = parse_expression("Sqrt[-9]");
    Expr* res6 = evaluate(e6);
    ASSERT(res6 && res6->type == EXPR_FUNCTION && strcmp(res6->data.function.head->data.symbol, "Complex") == 0);
    ASSERT(res6->data.function.args[1]->type == EXPR_INTEGER && res6->data.function.args[1]->data.integer == 3);
    expr_free(e6); expr_free(res6);
}

void test_apply() {
    symtab_clear_symbol("myApplyFunc");
    // Plus @@ {1, 2, 3} -> 6
    Expr* e1 = parse_expression("Plus @@ {1, 2, 3}");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_INTEGER && res1->data.integer == 6);
    expr_free(e1); expr_free(res1);

    // myApplyFunc @@ {a, b} -> myApplyFunc[a, b]
    Expr* e2 = parse_expression("myApplyFunc @@ {a, b}");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && strcmp(res2->data.function.head->data.symbol, "myApplyFunc") == 0);
    ASSERT(res2->data.function.arg_count == 2);
    expr_free(e2); expr_free(res2);

    // myApplyFunc @@@ {{a, b}, {c, d}} -> {myApplyFunc[a, b], myApplyFunc[c, d]}
    Expr* e3 = parse_expression("myApplyFunc @@@ {{a, b}, {c, d}}");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && strcmp(res3->data.function.head->data.symbol, "List") == 0);
    ASSERT(res3->data.function.arg_count == 2);
    ASSERT(res3->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res3->data.function.args[0]->data.function.head->data.symbol, "myApplyFunc") == 0);
    expr_free(e3); expr_free(res3);
}

void test_map() {
    symtab_clear_symbol("f");
    // f /@ {a, b, c} -> {f[a], f[b], f[c]}
    Expr* e1 = parse_expression("f /@ {a, b, c}");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol, "List") == 0);
    ASSERT(res1->data.function.arg_count == 3);
    ASSERT(res1->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res1->data.function.args[0]->data.function.head->data.symbol, "f") == 0);
    expr_free(e1); expr_free(res1);

    // Map[f, {{a, b}, {c, d}}, {2}] -> {{f[a], f[b]}, {f[c], f[d]}}
    Expr* e2 = parse_expression("Map[f, {{a, b}, {c, d}}, {2}]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 2);
    ASSERT(res2->data.function.args[0]->data.function.args[0]->type == EXPR_FUNCTION);
    expr_free(e2); expr_free(res2);
}

void test_map_all() {
    symtab_clear_symbol("f");
    // f //@ {a, {b}} -> f[{f[a], f[{f[b]}]}]
    Expr* e1 = parse_expression("f //@ {a, {b}}");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol, "f") == 0);
    expr_free(e1); expr_free(res1);
}

void test_select() {
    // Select elements that are even:
    // Select[{1,2,4,7,6,2},EvenQ] -> {2,4,6,2}
    Expr* e1 = parse_expression("Select[{1,2,4,7,6,2},EvenQ]");
    Expr* res1 = evaluate(e1);
    char* s1 = expr_to_string_fullform(res1);
    ASSERT_STR_EQ(s1, "List[2, 4, 6, 2]");
    free(s1);
    expr_free(e1); expr_free(res1);

    // Use a pure function to test each element:
    // Select[{1,2,4,7,6,2},#>2&] -> {4,7,6}
    Expr* e2 = parse_expression("Select[{1,2,4,7,6,2},#>2&]");
    Expr* res2 = evaluate(e2);
    char* s2 = expr_to_string_fullform(res2);
    ASSERT_STR_EQ(s2, "List[4, 7, 6]");
    free(s2);
    expr_free(e2); expr_free(res2);

    // Return only the first expression selected:
    // Select[{1,2,4,7,6,2},#>2&,1] -> {4}
    Expr* e3 = parse_expression("Select[{1,2,4,7,6,2},#>2&,1]");
    Expr* res3 = evaluate(e3);
    char* s3 = expr_to_string_fullform(res3);
    ASSERT_STR_EQ(s3, "List[4]");
    free(s3);
    expr_free(e3); expr_free(res3);
}

void test_matchq() {
    // MatchQ[1, _] -> True
    Expr* e1 = parse_expression("MatchQ[1, _]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "True") == 0);
    expr_free(e1); expr_free(res1);

    // MatchQ[f[x], f[_]] -> True
    Expr* e2 = parse_expression("MatchQ[f[x], f[_]]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol, "True") == 0);
    expr_free(e2); expr_free(res2);

    // MatchQ[f[x], g[_]] -> False
    Expr* e3 = parse_expression("MatchQ[f[x], g[_]]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol, "False") == 0);
    expr_free(e3); expr_free(res3);
}

void test_compoundexpression() {
    Expr* e1 = parse_expression("a; b; c;");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol, "Null") == 0);
    expr_free(res1); expr_free(e1);

    Expr* e2 = parse_expression("x = 11; y = 12");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_INTEGER && res2->data.integer == 12);
    expr_free(res2); expr_free(e2);
}

int main() {
    symtab_init();
    core_init();
    TEST(test_hold);
    TEST(test_orderless);
    TEST(test_flat);
    TEST(test_flat_and_orderless);
    TEST(test_set_and_ownvalues);
    TEST(test_setdelayed_and_downvalues);
    TEST(test_length);
    TEST(test_dimensions);
    TEST(test_part);
    TEST(test_builtin_head);
    TEST(test_first_last_most_rest);
    TEST(test_append_prepend);
    TEST(test_append_prepend_to);
    TEST(test_attributes);
    TEST(test_listable);
    TEST(test_oneidentity);
    TEST(test_list_assignment);
    TEST(test_plus);
    TEST(test_divide);
    TEST(test_rational_arithmetic);
    TEST(test_power);
    TEST(test_sqrt);
    TEST(test_apply);
    TEST(test_map);
    TEST(test_map_all);
    TEST(test_select);
    TEST(test_matchq);
    TEST(test_compoundexpression);
    printf("All tests passed!\n");
    return 0;
}