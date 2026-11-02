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
    ASSERT(strcmp(res->data.function.head->data.symbol.name, "Hold") == 0);
    
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
    ASSERT(strcmp(res->data.function.args[0]->data.symbol.name, "a") == 0);
    
    ASSERT(res->data.function.args[1]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[1]->data.symbol.name, "b") == 0);
    
    ASSERT(res->data.function.args[2]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[2]->data.symbol.name, "c") == 0);
    
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
    ASSERT(strcmp(res->data.function.args[0]->data.symbol.name, "a") == 0);
    
    ASSERT(res->data.function.args[1]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[1]->data.symbol.name, "b") == 0);
    
    ASSERT(res->data.function.args[2]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[2]->data.symbol.name, "c") == 0);
    
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
    ASSERT(strcmp(res1->data.function.head->data.symbol.name, "f") == 0);
    ASSERT(res1->data.function.args[0]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res1->data.function.args[0]->data.symbol.name, "x") == 0);
    expr_free(e1);
    expr_free(res1);

    Expr* e2 = parse_expression("f[f[x], f[y]]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 2);
    ASSERT(strcmp(res2->data.function.head->data.symbol.name, "f") == 0);
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
    ASSERT(strcmp(res->data.function.args[0]->data.symbol.name, "a") == 0);
    
    ASSERT(res->data.function.args[1]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[1]->data.symbol.name, "b") == 0);
    
    ASSERT(res->data.function.args[2]->type == EXPR_SYMBOL);
    ASSERT(strcmp(res->data.function.args[2]->data.symbol.name, "c") == 0);
    
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
    ASSERT(res1 != NULL && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "Null") == 0);
    
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
    ASSERT(strcmp(res1->data.function.head->data.symbol.name, "List") == 0);
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
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "b") == 0);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("{a, b, c}[[-1]]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol.name, "c") == 0);
    expr_free(e2); expr_free(res2);

    Expr* e3 = parse_expression("{a, b, c}[[0]]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol.name, "List") == 0);
    expr_free(e3); expr_free(res3);
}

void test_builtin_head() {
    Expr* e1 = parse_expression("Head[g[x]]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "g") == 0);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Head[123]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol.name, "Integer") == 0);
    expr_free(e2); expr_free(res2);

    /* Atomic heads. Use fresh symbols since earlier tests bind x, y, etc. */
    assert_eval_eq("Head[3.14]", "Real", 0);
#ifdef USE_MPFR
    /* Regression: an arbitrary-precision real (EXPR_MPFR) reports head Real,
     * same as a machine real. N[Pi, 24] previously left Head unevaluated
     * because expr_head had no EXPR_MPFR case. */
    assert_eval_eq("Head[N[Pi, 24]]", "Real", 0);
    assert_eval_eq("Head[N[3, 50]]", "Real", 0);
    assert_eval_eq("Head[N[Pi, 24], headWrap]", "headWrap[Real]", 0);
#endif
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
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "a") == 0);
    expr_free(e1); expr_free(res1);

    Expr* e2 = parse_expression("Last[{a, b, c}]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol.name, "c") == 0);
    expr_free(e2); expr_free(res2);

    Expr* e3 = parse_expression("Most[{a, b, c}]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && res3->data.function.arg_count == 2);
    ASSERT(strcmp(res3->data.function.args[0]->data.symbol.name, "a") == 0);
    ASSERT(strcmp(res3->data.function.args[1]->data.symbol.name, "b") == 0);
    expr_free(e3); expr_free(res3);

    Expr* e4 = parse_expression("Rest[{a, b, c}]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && res4->data.function.arg_count == 2);
    ASSERT(strcmp(res4->data.function.args[0]->data.symbol.name, "b") == 0);
    ASSERT(strcmp(res4->data.function.args[1]->data.symbol.name, "c") == 0);
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
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol.name, "List") == 0);
    // Should contain Flat, Listable, NumericFunction, OneIdentity, Orderless
    bool found_flat = false;
    for (size_t i = 0; i < res1->data.function.arg_count; i++) {
        if (strcmp(res1->data.function.args[i]->data.symbol.name, "Flat") == 0) found_flat = true;
    }
    ASSERT(found_flat);
    expr_free(e1); expr_free(res1);
}

void test_listable() {
    // Plus[{1, 2}, {3, 4}] -> {4, 6}
    Expr* e1 = parse_expression("Plus[{1, 2}, {3, 4}]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(res1->data.function.arg_count == 2);
    ASSERT(res1->data.function.args[0]->type == EXPR_INTEGER && res1->data.function.args[0]->data.integer == 4);
    ASSERT(res1->data.function.args[1]->type == EXPR_INTEGER && res1->data.function.args[1]->data.integer == 6);
    expr_free(e1); expr_free(res1);

    // {a,b,c} + 1 -> {Plus[1, a], Plus[1, b], Plus[1, c]}
    Expr* e2 = parse_expression("{a, b, c} + 1");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && res2->data.function.arg_count == 3);
    // Elements are sorted because Plus is Orderless: Plus[1, a]
    ASSERT(res2->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res2->data.function.args[0]->data.function.head->data.symbol.name, "Plus") == 0);
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
    ASSERT(res5->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res5->data.function.args[0]->data.function.head->data.symbol.name, "Times") == 0);
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
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "unbound") == 0);
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
    ASSERT(strcmp(res3->data.function.head->data.symbol.name, "Times") == 0);
    ASSERT(res3->data.function.args[0]->data.integer == 2);
    expr_free(e3); expr_free(res3);

    // Plus[2*x^2, 3*x^2] -> 5*x^2
    // Note: Power is symbolic for now
    Expr* e4 = parse_expression("2*x^2 + 3*x^2");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION);
    ASSERT(strcmp(res4->data.function.head->data.symbol.name, "Times") == 0);
    ASSERT(res4->data.function.args[0]->data.integer == 5);
    expr_free(e4); expr_free(res4);

    // Plus[1, 2*x, x^2, 4*x] -> Plus[1, 6*x, x^2]
    Expr* e6 = parse_expression("1 + 2*x + x^2 + 4*x");
    Expr* res6 = evaluate(e6);
    ASSERT(res6 && res6->type == EXPR_FUNCTION);
    ASSERT(strcmp(res6->data.function.head->data.symbol.name, "Plus") == 0);
    ASSERT(res6->data.function.arg_count == 3);
    // Check canonical order: 1, 6*x, x^2
    ASSERT(res6->data.function.args[0]->type == EXPR_INTEGER && res6->data.function.args[0]->data.integer == 1);
    ASSERT(res6->data.function.args[1]->type == EXPR_FUNCTION && strcmp(res6->data.function.args[1]->data.function.head->data.symbol.name, "Times") == 0);
    ASSERT(res6->data.function.args[2]->type == EXPR_FUNCTION && strcmp(res6->data.function.args[2]->data.function.head->data.symbol.name, "Power") == 0);
    expr_free(e6); expr_free(res6);

    // Plus[1, x, -1] -> x
    Expr* e5 = parse_expression("1 + x - 1");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_SYMBOL && strcmp(res5->data.symbol.name, "x") == 0);
    expr_free(e5); expr_free(res5);
}

void test_divide() {
    // Divide[1, 2] -> Rational[1, 2]
    Expr* e1 = parse_expression("Divide[1, 2]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol.name, "Rational") == 0);
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
    ASSERT(res4 && res4->type == EXPR_FUNCTION && strcmp(res4->data.function.head->data.symbol.name, "Rational") == 0);
    ASSERT(res4->data.function.args[0]->data.integer == 21);
    ASSERT(res4->data.function.args[1]->data.integer == 11);
    expr_free(e4); expr_free(res4);

    // Divide[x, y] -> Times[x, Power[y, -1]]
    Expr* e5 = parse_expression("Divide[x, y]");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && strcmp(res5->data.function.head->data.symbol.name, "Times") == 0);
    ASSERT(res5->data.function.arg_count == 2);
    ASSERT(strcmp(res5->data.function.args[1]->data.function.head->data.symbol.name, "Power") == 0);
    expr_free(e5); expr_free(res5);
}

void test_rational_arithmetic() {
    // 1 + 1/2 -> 3/2
    Expr* e1 = parse_expression("1 + 1/2");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol.name, "Rational") == 0);
    ASSERT(res1->data.function.args[0]->data.integer == 3);
    ASSERT(res1->data.function.args[1]->data.integer == 2);
    expr_free(e1); expr_free(res1);

    // 3/2 * 1/2 -> 3/4
    Expr* e2 = parse_expression("3/2 * 1/2");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_FUNCTION && strcmp(res2->data.function.head->data.symbol.name, "Rational") == 0);
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
    ASSERT(res4 && res4->type == EXPR_FUNCTION && strcmp(res4->data.function.head->data.symbol.name, "Times") == 0);
    // Elements should be Power[a, 2] and Power[b, 2]
    expr_free(e4); expr_free(res4);

    // (a^2)^3 -> a^6
    Expr* e5 = parse_expression("(a^2)^3");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && strcmp(res5->data.function.head->data.symbol.name, "Power") == 0);
    ASSERT(res5->data.function.args[1]->type == EXPR_INTEGER && res5->data.function.args[1]->data.integer == 6);
    expr_free(e5); expr_free(res5);

    // (-1)^(1/2) -> Complex[0, 1]
    Expr* e6 = parse_expression("(-1)^(1/2)");
    Expr* res6 = evaluate(e6);
    ASSERT(res6 && res6->type == EXPR_FUNCTION && strcmp(res6->data.function.head->data.symbol.name, "Complex") == 0);
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
    ASSERT(res2 && res2->type == EXPR_FUNCTION && strcmp(res2->data.function.head->data.symbol.name, "Power") == 0);
    expr_free(e2); expr_free(res2);

    // Sqrt[{1, 4, 9}] -> {1, 2, 3}
    Expr* e3 = parse_expression("Sqrt[{1, 4, 9}]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && strcmp(res3->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(res3->data.function.arg_count == 3);
    ASSERT(res3->data.function.args[0]->data.integer == 1);
    ASSERT(res3->data.function.args[1]->data.integer == 2);
    ASSERT(res3->data.function.args[2]->data.integer == 3);
    expr_free(e3); expr_free(res3);

    // Sqrt[x^2] -> Power[Power[x, 2], 1/2]
    Expr* e4 = parse_expression("Sqrt[x^2]");
    Expr* res4 = evaluate(e4);
    ASSERT(res4 && res4->type == EXPR_FUNCTION && strcmp(res4->data.function.head->data.symbol.name, "Power") == 0);
    // Should NOT automatically simplify to x
    expr_free(e4); expr_free(res4);

    // Sqrt[45] -> 3*Sqrt[5]
    Expr* e5 = parse_expression("Sqrt[45]");
    Expr* res5 = evaluate(e5);
    ASSERT(res5 && res5->type == EXPR_FUNCTION && strcmp(res5->data.function.head->data.symbol.name, "Times") == 0);
    ASSERT(res5->data.function.arg_count == 2);
    ASSERT(res5->data.function.args[0]->type == EXPR_INTEGER && res5->data.function.args[0]->data.integer == 3);
    expr_free(e5); expr_free(res5);

    // Sqrt[-9] -> 3*I
    Expr* e6 = parse_expression("Sqrt[-9]");
    Expr* res6 = evaluate(e6);
    ASSERT(res6 && res6->type == EXPR_FUNCTION && strcmp(res6->data.function.head->data.symbol.name, "Complex") == 0);
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
    ASSERT(res2 && res2->type == EXPR_FUNCTION && strcmp(res2->data.function.head->data.symbol.name, "myApplyFunc") == 0);
    ASSERT(res2->data.function.arg_count == 2);
    expr_free(e2); expr_free(res2);

    // myApplyFunc @@@ {{a, b}, {c, d}} -> {myApplyFunc[a, b], myApplyFunc[c, d]}
    Expr* e3 = parse_expression("myApplyFunc @@@ {{a, b}, {c, d}}");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_FUNCTION && strcmp(res3->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(res3->data.function.arg_count == 2);
    ASSERT(res3->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res3->data.function.args[0]->data.function.head->data.symbol.name, "myApplyFunc") == 0);
    expr_free(e3); expr_free(res3);
}

void test_map() {
    symtab_clear_symbol("f");
    // f /@ {a, b, c} -> {f[a], f[b], f[c]}
    Expr* e1 = parse_expression("f /@ {a, b, c}");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol.name, "List") == 0);
    ASSERT(res1->data.function.arg_count == 3);
    ASSERT(res1->data.function.args[0]->type == EXPR_FUNCTION && strcmp(res1->data.function.args[0]->data.function.head->data.symbol.name, "f") == 0);
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
    ASSERT(res1 && res1->type == EXPR_FUNCTION && strcmp(res1->data.function.head->data.symbol.name, "f") == 0);
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
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "True") == 0);
    expr_free(e1); expr_free(res1);

    // MatchQ[f[x], f[_]] -> True
    Expr* e2 = parse_expression("MatchQ[f[x], f[_]]");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_SYMBOL && strcmp(res2->data.symbol.name, "True") == 0);
    expr_free(e2); expr_free(res2);

    // MatchQ[f[x], g[_]] -> False
    Expr* e3 = parse_expression("MatchQ[f[x], g[_]]");
    Expr* res3 = evaluate(e3);
    ASSERT(res3 && res3->type == EXPR_SYMBOL && strcmp(res3->data.symbol.name, "False") == 0);
    expr_free(e3); expr_free(res3);
}

void test_compoundexpression() {
    Expr* e1 = parse_expression("a; b; c;");
    Expr* res1 = evaluate(e1);
    ASSERT(res1 && res1->type == EXPR_SYMBOL && strcmp(res1->data.symbol.name, "Null") == 0);
    expr_free(res1); expr_free(e1);

    Expr* e2 = parse_expression("x = 11; y = 12");
    Expr* res2 = evaluate(e2);
    ASSERT(res2 && res2->type == EXPR_INTEGER && res2->data.integer == 12);
    expr_free(res2); expr_free(e2);
}

/* An assignment evaluates its LHS head: h = g; h[x_] := ... defines g, and so
 * does a scoped local holding a computed symbol -- a package's cache builds
 * ch = ToExpression["cacheFunc" <> name] inside a Block and defines
 * ch[a_, b_, c_] := ...; before, the rule landed on the local and every cached
 * call stayed unevaluated. */
void test_assignment_evaluates_lhs_head() {
    assert_eval_eq("hh = gg; hh[x_] := x^2; hh[1] = 5; {gg[3], gg[1], Length[DownValues[gg]], DownValues[hh]}",
                   "{9, 5, 2, {}}", 0);
    assert_eval_eq("hh[1] =.; {gg[1], Length[DownValues[gg]]}", "{1, 1}", 0);
    assert_eval_eq("Block[{k = q2}, k[x_] := x^3]; q2[2]", "8", 0);
    assert_eval_eq("f1[a_] := Block[{ch}, ch = ToExpression[\"cf\" <> ToString[a]]; "
                   "ch[u_, v_] := {u, v, a}; ch[1, 2]]; {f1[zz], cfzz[3, 4]}",
                   "{{1, 2, zz}, {3, 4, zz}}", 0);
    assert_eval_eq("Module[{loc}, loc = gg2; loc[1] = 9]; gg2[1]", "9", 0);
    /* a head without an OwnValue is untouched */
    assert_eval_eq("ff3[x_] := x + 1; ff3[1]", "2", 0);
}

/* Semantics that FeynCalc's TID / FCApart / DiracTrace / Collect2 paths
 * rely on, each a kernel divergence found by running the package (v0.586):
 * rule order (a literal-headed call outranks a blank), Orderless DownValue
 * dispatch, a guard on a delayed rule's RHS (ReplaceAll, Cases, ReplaceList,
 * scoped), Unevaluated replacements, Delete/Extract with {} positions, nested
 * Part assignment, Together/ExpandAll leaving held parts alone, an aliased
 * option key, GroebnerBasis`DistributedTermsList and valued symbols inside an
 * assignment's LHS patterns. */
void test_feyncalc_semantics() {
    /* earlier groups leave x = 11, y = 12 and friends: start from clean symbols */
    assert_eval_eq("ClearAll[x, y, a, b, c, d, e, f, g, m, p, q, M]", "Null", 0);
    assert_eval_eq("c2[e_, {f___}, o:OptionsPattern[]] := c2 @@ Prepend[{f, o}, e]; c2[e_, f_, o:OptionsPattern[]] := {1, f}; c2[e_, f_, g__, o:OptionsPattern[]] := {2, f, g}; {c2[x, {F, G}], c2[x, {G}]}",
                   "List[List[2, F, G], List[1, G]]", 1);
    assert_eval_eq("SetAttributes[pr, Orderless]; pr[n_ x_mo, y_] := n pr[x, y]; {pr[mo[1], lm mo[1]], pr[lm mo[1], lm mo[1]]}",
                   "List[Times[lm, pr[mo[1], mo[1]]], Times[Power[lm, 2], pr[mo[1], mo[1]]]]", 1);
    assert_eval_eq("{f[2] /. f[x_] :> (x /; x > 1), f[0] /. f[x_] :> (x /; x > 1), f[2] /. f[x_] -> (x /; x > 1)}",
                   "List[2, f[0], Condition[2, Greater[2, 1]]]", 1);
    assert_eval_eq("Cases[{1, 2, 3}, x_ :> (x^2 /; x > 1)]",
                   "List[4, 9]", 1);
    assert_eval_eq("ReplaceList[a + b + c, x_ + y_ :> ({x, y} /; x === b)]",
                   "List[List[b, Plus[a, c]]]", 1);
    assert_eval_eq("{f[2] /. f[x_] :> Module[{t = x}, t /; t > 1], f[0] /. f[x_] :> Module[{t = x}, t /; t > 1]}",
                   "List[2, f[0]]", 1);
    assert_eval_eq("{f[2], g[3]} /. {f[x_] :> (x /; x > 5), f[x_] :> -x}",
                   "List[-2, g[3]]", 1);
    assert_eval_eq("(2 M[p] - M[q]) /. (n3_. M[x_] + n4_. M[y_]) :> (M[n3 x + n4 y] /; NumberQ[n3] && NumberQ[n4])",
                   "M[Plus[Times[2, p], Times[-1, q]]]", 1);
    assert_eval_eq("{d^2 f /. d -> Unevaluated[Sequence[]], d^2 /. d :> Unevaluated[Sequence[]], Hold[d] /. d :> Unevaluated[1 + 1], Hold[d, e] /. d :> Unevaluated[Sequence[a, b]]}",
                   "List[Times[2, f], 2, Hold[Plus[1, 1]], Hold[a, b, e]]", 1);
    assert_eval_eq("{Delete[{a, b}, {}], Delete[f[a, b], {}], Extract[{a, b}, {}], Extract[{a, b}, {{}}]}",
                   "List[List[a, b], f[a, b], List[], List[List[a, b]]]", 1);
    assert_eval_eq("vs = {{a, b}, {c}, {-1, 1, 1}}; vs[[3]][[1]] = vs[[3]][[1]] + 1; m = {{1, 2}, {3, 4}}; m[[2]][[1]] += 9; {vs, m}",
                   "List[List[List[a, b], List[c], List[0, 1, 1]], List[List[1, 2], List[12, 4]]]", 1);
    assert_eval_eq("ic[1] = a + b; {ExpandAll[f[HoldForm[ic[1]]]], Together[f[HoldForm[ic[1]]]]}",
                   "List[f[HoldForm[ic[1]]], f[HoldForm[ic[1]]]]", 1);
    assert_eval_eq("Options[fo3] = {al2 -> True}; al2 = real2; fo3[OptionsPattern[]] := OptionValue[al2]; fo3[]",
                   "True", 1);
    assert_eval_eq("{GroebnerBasis`DistributedTermsList[a x^2 + 3 x y + 5 y^3, {x, y}], GroebnerBasis`DistributedTermsList[y^2 + x, {x, y}, MonomialOrder -> DegreeReverseLexicographic], GroebnerBasis`DistributedTermsList[0, {x, y}]}",
                   "List[List[List[List[List[2, 0], a], List[List[1, 1], 3], List[List[0, 3], 5]], List[x, y]], List[List[List[List[0, 2], 1], List[List[1, 0], 1]], List[x, y]], List[List[], List[x, y]]]", 1);
    assert_eval_eq("PDx = PropD; fs[l_List][a__PDx] := {l, a}; fs[{1}][PropD[1], PropD[2]]",
                   "List[List[1], PropD[1], PropD[2]]", 1);
}

/* The pattern-free parts of an assignment's LHS pattern are evaluated
 * (lhs_fold_pattern_free): h[y_ - 4] is stored as h[-4 + y_], not as the
 * parse's Plus[y_, Times[-1, 4]], which matched nothing.  FeynCalc's BMHV Pair
 * rules on dim_Symbol - 4 depend on it.  Pattern names, guards and HoldPattern
 * bodies stay as written. */
void test_lhs_pattern_constants_folded() {
    assert_eval_eq("ClearAll[lf1, lf2, lf3, lf4, lf5, lf6, lfpr, li, mo, x, y, D4]", "Null", 0);
    assert_eval_eq("lf1[y_ - 4] := ok; {lf1[x - 4], lf1[x + 4]}", "List[ok, lf1[Plus[4, x]]]", 1);
    assert_eval_eq("lf2[y_, z_ - 4] := {y, z}; lf2[D4, D4 - 4]", "List[D4, D4]", 1);
    assert_eval_eq("lf3[x_/2] := x; lf3[y/2]", "y", 0);
    assert_eval_eq("lf4[x_ /; x > 2 - 1] := big; {lf4[3], lf4[0]}", "List[big, lf4[0]]", 1);
    assert_eval_eq("lf5[HoldPattern[y_ - 4]] := held; DownValues[lf5]",
                   "List[RuleDelayed[HoldPattern[lf5[HoldPattern[Plus[Pattern[y, Blank[]], -4]]]], held]]", 1);   /* a - 4 parses to Plus[a, -4] (v0.6139) */
    assert_eval_eq("SetAttributes[lfpr, Orderless]; "
                   "lfpr[(a : li | mo)[x_, dim_Symbol], (b : li | mo)[y_, dim_Symbol - 4]] := lfpr[a[x, dim - 4], b[y, dim - 4]]; "
                   "{lfpr[li[mu, D4], mo[p, D4 - 4]], lfpr[mo[p, D4 - 4], li[mu, D4]]}",
                   "List[lfpr[li[mu, Plus[-4, D4]], mo[p, Plus[-4, D4]]], lfpr[li[mu, Plus[-4, D4]], mo[p, Plus[-4, D4]]]]", 1);
}

/* Kernel semantics found by FeynCalc's unit tests (v0.591), each checked
 * against the reference: the sign of a -1 or -1/n coefficient distributes
 * over a sum; HoldPattern[f[...]] := ... is a rule for f; a protected tag is
 * Set::write; N uses N-values (x /: N[x] = v, N[x] = v, h /: N[h[x_]] := ...)
 * inside an expression and numericalizes their values. */
void test_feyncalc_semantics_v591() {
    assert_eval_eq("ClearAll[x, y, z, w, g, k, q, m, h, z2, z4]", "Null", 0);
    assert_eval_eq("{-(-1 + x), -(x + y) z, Times[-1/2, x + y], -(x + y)/2, -(2 x + 3)}",
                   "List[Plus[1, Times[-1, x]], Times[Plus[Times[-1, x], Times[-1, y]], z], "
                   "Times[Rational[1, 2], Plus[Times[-1, x], Times[-1, y]]], "
                   "Times[Rational[1, 2], Plus[Times[-1, x], Times[-1, y]]], Plus[-3, Times[-2, x]]]", 1);
    assert_eval_eq("{Times[-1, x + y, z], Times[-3/2, x + y], Times[-2, x + y], -(x + y)^2, Times[-I, x + y]}",
                   "List[Times[-1, Plus[x, y], z], Times[Rational[-3, 2], Plus[x, y]], Times[-2, Plus[x, y]], "
                   "Times[-1, Power[Plus[x, y], 2]], Times[Complex[0, -1], Plus[x, y]]]", 1);
    assert_eval_eq("HoldPattern[g[__, _]] := False; {g[1, 2], g[1], DownValues[g], DownValues[HoldPattern]}",
                   "List[False, g[1], List[RuleDelayed[HoldPattern[HoldPattern[g[BlankSequence[], Blank[]]]], False]], List[]]", 1);
    assert_eval_eq("HoldPattern[k[x_]] = 3; HoldPattern[q] = 5; HoldPattern[m[x_] /; x > 1] := big; {k[1], q, m[2], m[0]}",
                   "{3, 5, big, m[0]}", 0);
    assert_eval_eq("h /: N[h[x_]] := 2 x; z2 /: N[z2] = 3/2; {N[h[3]], N[1 + h[3]], N[2 + z2], N[z2]}",
                   "{6.0, 7.0, 3.5, 1.5}", 0);
    assert_eval_eq("N[z4] = 5; {N[z4], N[1 + z4], z4}", "{5.0, 6.0, z4}", 0);
    /* an N-value without a precision argument serves machine N only */
    assert_eval_eq("FreeQ[N[z2, 20], z2]", "False", 0);

    /* Indexed iterator variables bind v[k] as a DownValue, restored after. */
    assert_eval_eq("ClearAll[ll, kk, rr]; Table[ll[1], {ll[1], 0, 2}]", "{0, 1, 2}", 0);
    assert_eval_eq("rr = {}; Do[rr = Append[rr, {ll[1], ll[2]}], {ll[1], 0, 1}, {ll[2], 0, ll[1]}]; rr",
                   "{{0, 0}, {1, 0}, {1, 1}}", 0);
    assert_eval_eq("{Sum[kk[1]^2, {kk[1], 1, 3}], Product[kk[2], {kk[2], 1, 4}], Sum[f[kk[1]], {kk[1], 0}]}",
                   "{14, 24, 0}", 0);
    assert_eval_eq("ll[1] = 99; {Table[ll[1], {ll[1], 2}], ll[1]}", "{{1, 2}, 99}", 0);
    /* A head symbol with a value inside an LHS pattern is evaluated. */
    assert_eval_eq("ClearAll[hh9, gg9, f9, DD9, ds9]; hh9 = gg9; f9[hh9[a_]] := {a}; f9[gg9[1]]", "{1}", 0);
    assert_eval_eq("DD9 = Dot; ds9[DD9[a_, b_]] := ds9[a, b]; ds9[x . y]", "ds9[x, y]", 0);
    /* OptionsPattern under a Condition; a named OptionsPattern binds the
     * arguments as written. */
    assert_eval_eq("ClearAll[iso9, f8]; iso9[ex_, (o:OptionsPattern[]) /; {o} =!= {}] := {ex, {o}}; {iso9[x, aa -> 3], iso9[x]}",
                   "{{x, {aa -> 3}}, iso9[x]}", 0);
    assert_eval_eq("f8[x_, o:OptionsPattern[]] := {o}; f8[1, {a -> 1}, b -> 2]", "{{a -> 1}, b -> 2}", 0);
    assert_eval_eq("MatchQ[f[1], f[x_, (o:OptionsPattern[]) /; True]]", "True", 0);
}

/* Normal forms the sign distribution leans on (v0.591), checked against the
 * reference: Root keeps a positive leading coefficient; a negative numeric
 * factor under a radical comes out as its magnitude; FactorList reports the
 * -1 of a factor whose sign the product absorbed; series coefficients and
 * FunctionExpand pairs that cancel only after expansion are recognised. */
void test_sign_normal_forms_v591() {
    assert_eval_eq("ClearAll[x, y, s, a, b]", "Null", 0);
    assert_eval_eq("Coefficient[First[First[Root[-#^3 - 2 # + 1 &, 1]]] /. Slot[1] -> y, y, 3]", "1", 0);
    assert_eval_eq("{Sqrt[-4 (1 - x^2)], Sqrt[-4 a], (-8 a)^(1/3)}",
                   "{2 Sqrt[-1 + x^2], 2 Sqrt[-a], 2 (-a)^(1/3)}", 0);
    assert_eval_eq("Solve[1 - x^2 + y^2 == 0, y]", "{{y -> -Sqrt[-1 + x^2]}, {y -> Sqrt[-1 + x^2]}}", 0);
    assert_eval_eq("{FactorList[-x - 1], FactorList[1 - x], FactorList[x - y], FactorList[y - x], FactorList[-x/2 - 1/2]}",
                   "{{{-1, 1}, {1 + x, 1}}, {{-1, 1}, {-1 + x, 1}}, {{1, 1}, {x - y, 1}}, {{-1, 1}, {x - y, 1}}, {{-1, 1}, {2, -1}, {1 + x, 1}}}", 0);
    assert_eval_eq("Series[Beta[s, 1], {s, 0, 1}]", "1/s + O[s]^2", 0);
    assert_eval_eq("FunctionExpand[PolyLog[2, E^(I Pi/3)] + PolyLog[2, E^(-I Pi/3)]]", "1/18 Pi^2", 0);
}

/* Substitution splices a Sequence only where a pattern variable's value goes
 * (match.c rb_rec); a literal Sequence in the replacement is the evaluator's
 * to remove, and Unevaluated / HoldComplete keep it. Checked against the
 * reference. */
void test_sequence_substitution_v592() {
    assert_eval_eq("ClearAll[hz, h5, q2, hc, g, a, b, x, y]; hz[n_] := Table[If[i > 5, {i}, Unevaluated[Sequence[]]], {i, 1, n}]; hz[3]",
                   "{}", 0);
    assert_eval_eq("h5[n_] := HoldComplete[Sequence[], n]; h5[3]", "HoldComplete[Sequence[], 3]", 0);
    assert_eval_eq("q2[] := Hold[Unevaluated[Sequence[]]]; q2[]", "Hold[Unevaluated[Sequence[]]]", 0);
    assert_eval_eq("hc[x__] := HoldComplete[x]; hc[1, 2]", "HoldComplete[1, 2]", 0);
    assert_eval_eq("g[a, b] /. g[x_, y_] :> g[x, If[True, Unevaluated[Sequence[]], 1], y]", "g[a, b]", 0);
}

/* Coefficient as the reference computes it, and nested option lists. */
/* Numeric bases under one symbolic exponent merge, as in the reference.
 * Comparisons by === so the checks do not depend on printed factor order. */
void test_numeric_base_power_merge_v594() {
    assert_eval_eq("ClearAll[a, b, x, y, f]", "Null", 0);
    assert_eval_eq("{2^a 3^a === 6^a, 3^-a 5^a === (5/3)^a, 2^a 3^a 5^-a === (6/5)^a, "
                   "3^-a 5^-a === 15^(-a), 2^(a b) 3^(a b) === 6^(a b)}",
                   "{True, True, True, True, True}", 0);
    /* One base that is not a positive real may join; two may not. */
    assert_eval_eq("{(-2)^a 3^a === (-6)^a, (-2)^-a 3^a === (-2/3)^(-a), "
                   "Head[(-2)^a (-3)^a], 2^a I^a === (2 I)^a}",
                   "{True, True, Times, True}", 0);
    /* Numeric-quantity bases merge; symbolic ones and unequal exponents do not. */
    assert_eval_eq("{Sinh[1]^a Cosh[1]^a === (Cosh[1] Sinh[1])^a, "
                   "Log[2]^a Log[3]^a === (Log[2] Log[3])^a}", "{True, True}", 0);
    assert_eval_eq("{Head[x^a y^a], Head[2^a f[1]^a], Head[2^a 3^(2 a)], Head[2^a 3^b]}",
                   "{Times, Times, Times, Times}", 0);
    assert_eval_eq("Simplify[Sinh[1]^a (1 + Cosh[1])^(-a)]", "Tanh[1/2]^a", 0);
}

/* A power of a product splits off its positive numeric factors inside a
 * product, as in the reference (every expected FullForm is the reference's);
 * a lone power of a product, and a numeric base with nothing to combine
 * with, stay. FeynCalc's FCFeynmanParametrize needs Pi^(D/2) (2 Pi)^(-D) as
 * 2^(-D) Pi^(-D/2). */
void test_power_of_product_split_o3() {
    assert_eval_eq("ClearAll[a, b, c, x, y, D]", "Null", 0);
    assert_eval_eq("y (2 x)^b", "Times[Power[2, b], Power[x, b], y]", 1);
    assert_eval_eq("x^a (2 x)^b", "Times[Power[2, b], Power[x, Plus[a, b]]]", 1);
    assert_eval_eq("x^a (-2 x)^b", "Times[Power[2, b], Power[Times[-1, x], b], Power[x, a]]", 1);
    assert_eval_eq("x^a (2 x y)^b", "Times[Power[2, b], Power[x, a], Power[Times[x, y], b]]", 1);
    assert_eval_eq("x^a (x/2)^b", "Times[Power[2, Times[-1, b]], Power[x, Plus[a, b]]]", 1);
    assert_eval_eq("y (3/2 x)^b", "Times[Power[Rational[3, 2], b], Power[x, b], y]", 1);
    assert_eval_eq("y (2 I x)^b", "Times[Power[2, b], Power[Times[Complex[0, 1], x], b], y]", 1);
    assert_eval_eq("y (2 Sqrt[2] x)^b", "Times[Power[2, Times[Rational[3, 2], b]], Power[x, b], y]", 1);
    assert_eval_eq("(2 x)^b (3 x)^b", "Times[Power[6, b], Power[x, Times[2, b]]]", 1);
    assert_eval_eq("(2 x)^b (2 x)^c", "Times[Power[2, Plus[b, c]], Power[x, Plus[b, c]]]", 1);
    assert_eval_eq("x^a (Pi x)^b", "Times[Power[Pi, b], Power[x, Plus[a, b]]]", 1);
    /* a numeric base splits only to combine with another factor */
    assert_eval_eq("Pi^a (2 Pi)^b", "Times[Power[2, b], Power[Pi, Plus[a, b]]]", 1);
    assert_eval_eq("2^a (2 Pi)^b", "Times[Power[2, Plus[a, b]], Power[Pi, b]]", 1);
    assert_eval_eq("Pi^a (-2 Pi)^b", "Times[Power[-2, b], Power[Pi, Plus[a, b]]]", 1);
    assert_eval_eq("Pi^(D/2) (2 Pi)^(-D)", "Times[Power[2, Times[-1, D]], Power[Pi, Times[Rational[-1, 2], D]]]", 1);
    assert_eval_eq("(-2 Pi)^b y", "Times[Power[Times[-2, Pi], b], y]", 1);
    assert_eval_eq("{(2 x)^b, (2 Pi)^b}", "List[Power[Times[2, x], b], Power[Times[2, Pi], b]]", 1);
    assert_eval_eq("(x^a (I x)^b)", "Times[Power[Times[Complex[0, 1], x], b], Power[x, a]]", 1);
}

/* Numeric powers whose bases share factors are written over one coprime
 * basis when an exponent is symbolic (expected values are the reference's);
 * bases that share nothing stay. */
void test_numeric_power_coprime_basis_o4() {
    assert_eval_eq("ClearAll[a, b, x, y]", "Null", 0);
    assert_eval_eq("{(3/4)^a (5/2)^b 2^(2a+b), 2^a 4^b, 2^a 8^a, 4^a 6^b, (1/4)^a 2^(2a), 12^a 2^-a}",
                   "List[Times[Power[3, a], Power[5, b]], Power[2, Plus[a, Times[2, b]]], Power[16, a], "
                   "Times[Power[2, Plus[Times[2, a], b]], Power[3, b]], 1, Power[6, a]]", 1);
    assert_eval_eq("{(2/3)^a (9/4)^b, 4^a 2^(-2a) 3, 2^(1/2+a) Sqrt[2], 4^a/2, 8 2^a 4^b, 4 2^e, 12 2^n, 2^e/8}",
                   "List[Power[Rational[2, 3], Plus[a, Times[-2, b]]], 3, Power[2, Plus[1, a]], "
                   "Power[2, Plus[-1, Times[2, a]]], Power[2, Plus[3, a, Times[2, b]]], Power[2, Plus[2, e]], "
                   "Times[3, Power[2, Plus[2, n]]], Power[2, Plus[-3, e]]]", 1);
    /* gcds, no factoring: huge prime-power bases collapse */
    assert_eval_eq("1/379037881^7 (p/379037881^4)^(-1 - 3 e) (r/379037881^3)^(-1 + 4 e) - p^(-1 - 3 e) r^(-1 + 4 e)", "0", 1);
    assert_eval_eq("{3 4^a, 2^a 3^b, 2^x (2*10000000019)^y}",
                   "List[Times[3, Power[4, a]], Times[Power[2, a], Power[3, b]], "
                   "Times[Power[2, Plus[x, y]], Power[10000000019, y]]]", 1);
    /* only a symbolic exponent calls for the basis: 17 Sqrt[17] as 17^(3/2)
     * evaluates straight back, an endless loop through Times (reclim) */
    assert_eval_eq("{17 Sqrt[17], 9 Sqrt[3] 3^a, -153 Sqrt[17] 23^(-2 + 4 Sqrt[17]) 10^(2 - 4 Sqrt[17])}",
                   "List[Times[17, Power[17, Rational[1, 2]]], Power[3, Plus[Rational[5, 2], a]], "
                   "Times[-153, Power[Rational[10, 23], Plus[2, Times[-4, Power[17, Rational[1, 2]]]]], "
                   "Power[17, Rational[1, 2]]]]", 1);
    /* opposite exponents merge for any non-rational exponent, numeric ones
     * included; a sum keeps the sign whose LAST term is negative */
    assert_eval_eq("{2^Sqrt[2] 3^-Sqrt[2], 2^Pi 3^Pi, 2^I 3^-I, 2^(1+a) 3^(-1-a), 2^(a-1) 3^(1-a), "
                   "2^(-a+b) 3^(a-b), 2^(1+Sqrt[3]) 3^(-1-Sqrt[3]), 2^(1-Sqrt[2]) 3^(Sqrt[2]-1)}",
                   "{(2/3)^Sqrt[2], 6^Pi, (2/3)^I, (3/2)^(-1 - a), (3/2)^(1 - a), (3/2)^(a - b), "
                   "(3/2)^(-1 - Sqrt[3]), (2/3)^(1 - Sqrt[2])}", 0);
    /* complex coefficients: in a sum the last term reads as negative as the
     * printer writes it (-I t, (-1 + I) t); a product flips only for a
     * negative imaginary coefficient; a complex number keeps im > 0.  The
     * first case looped through Simplify's RadicalCanon (a stack overflow)
     * while the merged sign disagreed with the split one's. */
    assert_eval_eq("{4^(5/8 - I t) 5^(-5/8 + I t), 2^(-I t) 3^(I t), 2^((-1 + I) t) 3^((1 - I) t), "
                   "2^(1 + (-1 + I) t) 3^(-1 + (1 - I) t), 2^(-I) 3^I, 2^(1 - I) 3^(-1 + I), "
                   "2^(2.5 t) 3^(-2.5 t)}",
                   "{(4/5)^(5/8 - I t), (3/2)^(I t), (2/3)^((-1 + I) t), (2/3)^(1 + (-1 + I) t), "
                   "(3/2)^I, (3/2)^(-1 + I), 2^(2.5 t) 3^(-2.5 t)}", 0);
}

/* A condition does not promote a rule past an earlier, at-least-as-specific
 * one; between rules that differ only in the condition the guarded one goes
 * first. Expected values are the reference's (FeynCalc's fromGFAD rules). */
void test_downvalue_order_conditions() {
    assert_eval_eq("ClearAll[fH, fC, g1, g2, x, y]", "Null", 0);
    assert_eval_eq("fH[p[a_,a_]] := \"spec\"; fH[p[a1_,b1_]] := \"gen\" /; True; fH[p[x,x]]", "\"spec\"", 0);
    assert_eval_eq("fC[p[a_,a_] + q[a_,b_]] := \"spec\" /; True; "
                   "fC[p[a1_,a1_] + q[a2_,b2_]] := \"gen\" /; True /; True; fC[p[x,x] + q[x,y]]", "\"spec\"", 0);
    assert_eval_eq("g1[x_] := 1; g1[y_] /; y > 0 := 2; {g1[3], g1[-3]}", "{2, 1}", 0);
    assert_eval_eq("g2[y_] /; y > 0 := 2; g2[x_] := 1; {g2[3], g2[-3]}", "{2, 1}", 0);
}

/* OptionsPattern[] in any position takes a prefix of option arguments, and
 * SyntacticNegativeQ reads a sum by its first term and a complex by its
 * imaginary part. Expected values are the reference's (FeynCalc's
 * SUNSimplify: SUND is Orderless, so SUND[args__, OptionsPattern[]] sorts to
 * SUND[OptionsPattern[], args__]). */
void test_optionspattern_position_and_syntactic_negative() {
    assert_eval_eq("ClearAll[s, t, u, f, g, a, b, o, x, y, z]", "Null", 0);
    assert_eval_eq("SetAttributes[s, Orderless]; {MatchQ[s[x], s[a__, OptionsPattern[]]], "
                   "MatchQ[t[x], t[OptionsPattern[], a_]], s[x, y, z] /. s[args__, OptionsPattern[]] :> g[args]}",
                   "{True, True, g[x, y, z]}", 0);
    assert_eval_eq("{t[a -> 1, x] /. t[o:OptionsPattern[], b_] :> {{o}, b}, MatchQ[f[x], f[OptionsPattern[]]], "
                   "MatchQ[u[x, a -> 1, y], u[_, OptionsPattern[], _]], MatchQ[t[a -> 1], t[OptionsPattern[], b_]]}",
                   "{{{a -> 1}, x}, False, True, True}", 0);
    assert_eval_eq("Internal`SyntacticNegativeQ /@ {-x - y, -1 + x, x - y, -2 x + y, -I, Complex[-1, 2], Complex[1, -2], (-x)^2}",
                   "{True, True, False, True, True, False, True, False}", 0);
}

void test_coefficient_v592() {
    assert_eval_eq("ClearAll[a, b, c, d, l, x, y, ep, x1, x2, x3, ff, gg]", "Null", 0);
    assert_eval_eq("{Coefficient[a/x, x, -1], Coefficient[a/x + b x, x, -1], Coefficient[a/x^2 + 3, x, -2]}",
                   "{a, a, a}", 0);
    assert_eval_eq("{Coefficient[x1/ep + x2 + ep x3, ep, 0], Coefficient[x/(1 + x), x], Coefficient[x/(1 + x), x, 0]}",
                   "{x2, 0, 0}", 0);
    assert_eval_eq("{Coefficient[x^3 + 1, x^2, 0], Coefficient[x Sin[x] + x, x], Coefficient[a x y + b x, x y]}",
                   "{1, 1 + Sin[x], a}", 0);
    assert_eval_eq("Coefficient[(a l^2 + b l^2 + c)/d, l, 2] === (a + b)/d", "True", 0);
    assert_eval_eq("Coefficient[(a + b) x/c + d x, x] === (a + b)/c + d", "True", 0);
    assert_eval_eq("{OptionValue[ff, {{a -> 1}}, a], OptionValue[ff, {{a -> 1}, b -> 2}, b], OptionValue[ff, {{{a -> 5}}}, a]}",
                   "{1, 2, 5}", 0);
    assert_eval_eq("Options[ff] = {a -> 0}; gg[x_, OptionsPattern[ff]] := {x, OptionValue[a]}; {gg[1, {{a -> 5}}], gg[3]}",
                   "{{1, 5}, {3, 0}}", 0);
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
    TEST(test_assignment_evaluates_lhs_head);
    TEST(test_feyncalc_semantics);
    TEST(test_lhs_pattern_constants_folded);
    TEST(test_feyncalc_semantics_v591);
    TEST(test_sign_normal_forms_v591);
    TEST(test_sequence_substitution_v592);
    TEST(test_coefficient_v592);
    TEST(test_numeric_base_power_merge_v594);
    TEST(test_power_of_product_split_o3);
    TEST(test_numeric_power_coprime_basis_o4);
    TEST(test_downvalue_order_conditions);
    TEST(test_optionspattern_position_and_syntactic_negative);
    printf("All tests passed!\n");
    return 0;
}