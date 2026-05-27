#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "arithmetic.h"
#include "parse.h"
#include <string.h>

void test_sameq_basic() {
    Expr* x1 = expr_new_symbol("x");
    Expr* x2 = expr_new_symbol("x");
    Expr* args[] = {x1, x2};
    Expr* sameq = expr_new_function(expr_new_symbol("SameQ"), args, 2);
    
    Expr* res = evaluate(sameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(sameq);
    expr_free(res);
}

void test_sameq_different() {
    Expr* x = expr_new_symbol("x");
    Expr* y = expr_new_symbol("y");
    Expr* args[] = {x, y};
    Expr* sameq = expr_new_function(expr_new_symbol("SameQ"), args, 2);
    
    Expr* res = evaluate(sameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    
    expr_free(sameq);
    expr_free(res);
}

void test_sameq_multiple() {
    Expr* x1 = expr_new_symbol("x");
    Expr* x2 = expr_new_symbol("x");
    Expr* x3 = expr_new_symbol("x");
    Expr* args[] = {x1, x2, x3};
    Expr* sameq = expr_new_function(expr_new_symbol("SameQ"), args, 3);
    
    Expr* res = evaluate(sameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(sameq);
    expr_free(res);
}

void test_sameq_multiple_false() {
    Expr* x1 = expr_new_symbol("x");
    Expr* x2 = expr_new_symbol("y");
    Expr* x3 = expr_new_symbol("x");
    Expr* args[] = {x1, x2, x3};
    Expr* sameq = expr_new_function(expr_new_symbol("SameQ"), args, 3);
    
    Expr* res = evaluate(sameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    
    expr_free(sameq);
    expr_free(res);
}

void test_unsameq_basic() {
    Expr* x = expr_new_symbol("x");
    Expr* y = expr_new_symbol("y");
    Expr* args[] = {x, y};
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), args, 2);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_false() {
    Expr* x1 = expr_new_symbol("x");
    Expr* x2 = expr_new_symbol("x");
    Expr* args[] = {x1, x2};
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), args, 2);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_multiple() {
    Expr* x = expr_new_symbol("x");
    Expr* y = expr_new_symbol("y");
    Expr* z = expr_new_symbol("z");
    Expr* args[] = {x, y, z};
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), args, 3);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_multiple_false() {
    Expr* x1 = expr_new_symbol("x");
    Expr* y = expr_new_symbol("y");
    Expr* x2 = expr_new_symbol("x");
    Expr* args[] = {x1, y, x2};
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), args, 3);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_zero_arg() {
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), NULL, 0);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_one_arg() {
    Expr* x = expr_new_symbol("x");
    Expr* args[] = {x};
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), args, 1);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_equal_identical() {
    Expr* x1 = expr_new_symbol("x");
    Expr* x2 = expr_new_symbol("x");
    Expr* args[] = {x1, x2};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 2);
    
    Expr* res = evaluate(equal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(equal);
    expr_free(res);
}

void test_equal_numeric() {
    Expr* two_int = expr_new_integer(2);
    Expr* two_real = expr_new_real(2.0);
    Expr* args[] = {two_int, two_real};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 2);
    
    Expr* res = evaluate(equal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(equal);
    expr_free(res);
}

void test_equal_unevaluated() {
    Expr* x = expr_new_symbol("x");
    Expr* y = expr_new_symbol("y");
    Expr* args[] = {x, y};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 2);
    
    Expr* res = evaluate(equal);
    // Returns original Equal[x, y]
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol, "Equal");
    
    expr_free(equal);
    expr_free(res);
}

void test_equal_numeric_false() {
    Expr* one = expr_new_integer(1);
    Expr* two = expr_new_integer(2);
    Expr* args[] = {one, two};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 2);
    
    Expr* res = evaluate(equal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    
    expr_free(equal);
    expr_free(res);
}

void test_equal_tolerance() {
    Expr* a = expr_new_real(1.0);
    Expr* b = expr_new_real(1.0 + 1e-15);
    Expr* args[] = {a, b};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 2);
    
    Expr* res = evaluate(equal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    
    expr_free(equal);
    expr_free(res);
}

void test_equal_mixed_false() {
    Expr* x = expr_new_symbol("x");
    Expr* one = expr_new_integer(1);
    Expr* two = expr_new_integer(2);
    Expr* args[] = {x, one, two};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 3);
    
    Expr* res = evaluate(equal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    
    expr_free(equal);
    expr_free(res);
}

void test_equal_mixed_unevaluated() {
    Expr* one = expr_new_integer(1);
    Expr* one_again = expr_new_integer(1);
    Expr* x = expr_new_symbol("x");
    Expr* args[] = {one, one_again, x};
    Expr* equal = expr_new_function(expr_new_symbol("Equal"), args, 3);
    
    Expr* res = evaluate(equal);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol, "Equal");
    
    expr_free(equal);
    expr_free(res);
}

void test_less_basic() {
    Expr* one = expr_new_integer(1);
    Expr* two = expr_new_integer(2);
    Expr* args[] = {one, two};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    expr_free(less); expr_free(res);
}

void test_less_false() {
    Expr* two = expr_new_integer(2);
    Expr* one = expr_new_integer(1);
    Expr* args[] = {two, one};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    expr_free(less); expr_free(res);
}

void test_less_rational_real() {
    Expr* half = make_rational(1, 2);
    Expr* zero_point_six = expr_new_real(0.6);
    Expr* args[] = {half, zero_point_six};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    expr_free(less); expr_free(res);
}

void test_less_exact_rational() {
    Expr* one_third = make_rational(1, 3);
    Expr* one_half = make_rational(1, 2);
    Expr* args[] = {one_third, one_half};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    expr_free(less); expr_free(res);
}

void test_lessequal_mixed() {
    Expr* a = make_rational(3, 4);
    Expr* b = expr_new_real(0.75);
    Expr* c = expr_new_integer(1);
    Expr* args[] = {a, b, c};
    Expr* le = expr_new_function(expr_new_symbol("LessEqual"), args, 3);
    Expr* res = evaluate(le);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    expr_free(le); expr_free(res);
}

void test_greater_unevaluated() {
    Expr* x = expr_new_symbol("x");
    Expr* one = expr_new_integer(1);
    Expr* args[] = {x, one};
    Expr* gt = expr_new_function(expr_new_symbol("Greater"), args, 2);
    Expr* res = evaluate(gt);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol, "Greater");
    expr_free(gt); expr_free(res);
}

/* Regression: BigInt-vs-BigInt comparisons used to coerce both sides to
 * double, which collapsed neighbouring values above 2^53 to a single
 * representation and returned False for Less[10^20, 10^20 + 1]. */
static void check_compare_yields(const char* expr_text, const char* expected) {
    Expr* in = parse_expression(expr_text);
    ASSERT(in != NULL);
    Expr* res = evaluate(in);
    ASSERT(res != NULL);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, expected);
    expr_free(in);
    expr_free(res);
}

void test_compare_bigint_adjacent() {
    check_compare_yields("Less[10^20, 10^20 + 1]",          "True");
    check_compare_yields("Less[10^30, 10^30 + 1]",          "True");
    check_compare_yields("Less[10^50 - 1, 10^50]",          "True");
    check_compare_yields("Greater[10^30, 10^30 - 1]",       "True");
    check_compare_yields("LessEqual[10^20 + 1, 10^20]",     "False");
    check_compare_yields("GreaterEqual[10^30 - 1, 10^30]",  "False");
    check_compare_yields("Equal[10^30, 10^30]",             "True");
    check_compare_yields("Equal[10^30, 10^30 + 1]",         "False");
    check_compare_yields("Unequal[10^30, 10^30 + 1]",       "True");
    /* Mixed Integer/BigInt: 5 < 10^30 must hold exactly. */
    check_compare_yields("Less[5, 10^30]",                  "True");
    check_compare_yields("Greater[10^30, 5]",               "True");
    /* Sanity: small-integer path still works. */
    check_compare_yields("Less[1, 2]",                      "True");
    check_compare_yields("Less[2, 1]",                      "False");
    check_compare_yields("LessEqual[1.5, 2.5]",             "True");
    check_compare_yields("Less[3/4, 4/5]",                  "True");
}

void test_unequal_basic() {
    Expr* one = expr_new_integer(1);
    Expr* two = expr_new_integer(2);
    Expr* args[] = {one, two};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 2);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "True");
    expr_free(unequal); expr_free(res);
}

void test_unequal_false() {
    Expr* one = expr_new_integer(1);
    Expr* one_point_zero = expr_new_real(1.0);
    Expr* args[] = {one, one_point_zero};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 2);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    expr_free(unequal); expr_free(res);
}

void test_unequal_multiple() {
    Expr* a = expr_new_integer(1);
    Expr* b = expr_new_integer(2);
    Expr* c = expr_new_integer(1);
    Expr* args[] = {a, b, c};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 3);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol, "False");
    expr_free(unequal); expr_free(res);
}

void test_unequal_unevaluated() {
    Expr* one = expr_new_integer(1);
    Expr* x = expr_new_symbol("x");
    Expr* args[] = {one, x};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 2);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol, "Unequal");
    expr_free(unequal); expr_free(res);
}

void test_not_basic() {
    Expr* t = expr_new_symbol("True");
    Expr* args_t[] = {t};
    Expr* not_t = expr_new_function(expr_new_symbol("Not"), args_t, 1);
    Expr* res_t = evaluate(not_t);
    ASSERT(res_t->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_t->data.symbol, "False");
    expr_free(not_t); expr_free(res_t);

    Expr* f = expr_new_symbol("False");
    Expr* args_f[] = {f};
    Expr* not_f = expr_new_function(expr_new_symbol("Not"), args_f, 1);
    Expr* res_f = evaluate(not_f);
    ASSERT(res_f->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_f->data.symbol, "True");
    expr_free(not_f); expr_free(res_f);
}

void test_and_basic() {
    Expr* and_empty = expr_new_function(expr_new_symbol("And"), NULL, 0);
    Expr* res_empty = evaluate(and_empty);
    ASSERT(res_empty->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_empty->data.symbol, "True");
    expr_free(and_empty); expr_free(res_empty);

    Expr* t = expr_new_symbol("True");
    Expr* f = expr_new_symbol("False");
    Expr* args1[] = {expr_copy(t), expr_copy(f)};
    Expr* and1 = expr_new_function(expr_new_symbol("And"), args1, 2);
    Expr* res1 = evaluate(and1);
    ASSERT(res1->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res1->data.symbol, "False");
    expr_free(and1); expr_free(res1);

    Expr* args2[] = {expr_copy(t), expr_copy(t)};
    Expr* and2 = expr_new_function(expr_new_symbol("And"), args2, 2);
    Expr* res2 = evaluate(and2);
    ASSERT(res2->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res2->data.symbol, "True");
    expr_free(and2); expr_free(res2);
    
    expr_free(t); expr_free(f);
}

void test_or_basic() {
    Expr* or_empty = expr_new_function(expr_new_symbol("Or"), NULL, 0);
    Expr* res_empty = evaluate(or_empty);
    ASSERT(res_empty->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_empty->data.symbol, "False");
    expr_free(or_empty); expr_free(res_empty);

    Expr* t = expr_new_symbol("True");
    Expr* f = expr_new_symbol("False");
    Expr* args1[] = {expr_copy(f), expr_copy(t)};
    Expr* or1 = expr_new_function(expr_new_symbol("Or"), args1, 2);
    Expr* res1 = evaluate(or1);
    ASSERT(res1->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res1->data.symbol, "True");
    expr_free(or1); expr_free(res1);

    Expr* args2[] = {expr_copy(f), expr_copy(f)};
    Expr* or2 = expr_new_function(expr_new_symbol("Or"), args2, 2);
    Expr* res2 = evaluate(or2);
    ASSERT(res2->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res2->data.symbol, "False");
    expr_free(or2); expr_free(res2);
    
    expr_free(t); expr_free(f);
}

int main() {
    symtab_init();
    core_init();
    
    TEST(test_sameq_basic);
    TEST(test_sameq_different);
    TEST(test_sameq_multiple);
    TEST(test_sameq_multiple_false);
    TEST(test_unsameq_basic);
    TEST(test_unsameq_false);
    TEST(test_unsameq_multiple);
    TEST(test_unsameq_multiple_false);
    TEST(test_unsameq_zero_arg);
    TEST(test_unsameq_one_arg);
    TEST(test_equal_identical);
    TEST(test_equal_numeric);
    TEST(test_equal_unevaluated);
    TEST(test_equal_numeric_false);
    TEST(test_equal_tolerance);
    TEST(test_equal_mixed_false);
    TEST(test_equal_mixed_unevaluated);
    TEST(test_unequal_basic);
    TEST(test_unequal_false);
    TEST(test_unequal_multiple);
    TEST(test_unequal_unevaluated);
    TEST(test_less_basic);
    TEST(test_less_false);
    TEST(test_less_rational_real);
    TEST(test_less_exact_rational);
    TEST(test_lessequal_mixed);
    TEST(test_greater_unevaluated);
    TEST(test_compare_bigint_adjacent);
    TEST(test_not_basic);
    TEST(test_and_basic);
    TEST(test_or_basic);
    
    printf("All comparisons tests passed!\n");
    return 0;
}
