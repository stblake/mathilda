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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_zero_arg() {
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), NULL, 0);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
    expr_free(unsameq);
    expr_free(res);
}

void test_unsameq_one_arg() {
    Expr* x = expr_new_symbol("x");
    Expr* args[] = {x};
    Expr* unsameq = expr_new_function(expr_new_symbol("UnsameQ"), args, 1);
    
    Expr* res = evaluate(unsameq);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "Equal");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    
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
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "Equal");
    
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    expr_free(less); expr_free(res);
}

void test_less_false() {
    Expr* two = expr_new_integer(2);
    Expr* one = expr_new_integer(1);
    Expr* args[] = {two, one};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    expr_free(less); expr_free(res);
}

void test_less_rational_real() {
    Expr* half = make_rational(1, 2);
    Expr* zero_point_six = expr_new_real(0.6);
    Expr* args[] = {half, zero_point_six};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    expr_free(less); expr_free(res);
}

void test_less_exact_rational() {
    Expr* one_third = make_rational(1, 3);
    Expr* one_half = make_rational(1, 2);
    Expr* args[] = {one_third, one_half};
    Expr* less = expr_new_function(expr_new_symbol("Less"), args, 2);
    Expr* res = evaluate(less);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "True");
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
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    expr_free(le); expr_free(res);
}

void test_greater_unevaluated() {
    Expr* x = expr_new_symbol("x");
    Expr* one = expr_new_integer(1);
    Expr* args[] = {x, one};
    Expr* gt = expr_new_function(expr_new_symbol("Greater"), args, 2);
    Expr* res = evaluate(gt);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "Greater");
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
    ASSERT_STR_EQ(res->data.symbol.name, expected);
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

/* Exact rationals beyond the double range compare exactly: a BigInt
 * denominator past ~1e308 made the value 0 at machine precision, and the
 * exact zero test Equal consulted then confirmed 2^2000/3^1300 == 0 (the
 * value is ~6.34e-19).  PossibleZeroQ, which shares the test, likewise. */
void test_compare_rational_beyond_double_range() {
    check_compare_yields("Equal[2^2000/3^1300, 0]",          "False");
    check_compare_yields("Equal[0, 2^2000/3^1300]",          "False");
    check_compare_yields("Unequal[2^2000/3^1300, 0]",        "True");
    check_compare_yields("Equal[1/3^700, 0]",                "False");
    check_compare_yields("Equal[-1/3^700, 0]",               "False");
    check_compare_yields("Equal[1/3^700, 1/3^701]",          "False");
    check_compare_yields("Equal[1/3^700, 3/3^701]",          "True");
    check_compare_yields("Greater[2^2000/3^1300, 0]",        "True");
    check_compare_yields("Less[2^2000/3^1300, 10^-18]",      "True");
    check_compare_yields("PossibleZeroQ[2^2000/3^1300]",     "False");
    check_compare_yields("PossibleZeroQ[1/3^700 - 3/3^701]", "True");
}

/* IEEE 754 / ISO 60559: a NaN is unordered with every value, itself
 * included.  Indeterminate is Mathilda's NaN, so Equal / Less / LessEqual /
 * Greater / GreaterEqual all yield False against it and Unequal yields True.
 * SameQ and UnsameQ are exempt — they compare structure, not value. */
void test_indeterminate_ieee_unordered() {
    /* The headline case: Equal must not take expr_eq's identity shortcut. */
    check_compare_yields("Equal[Indeterminate, Indeterminate]",        "False");
    check_compare_yields("Unequal[Indeterminate, Indeterminate]",      "True");

    /* Unordered against ordinary values, in either operand position. */
    check_compare_yields("Equal[Indeterminate, 1]",                    "False");
    check_compare_yields("Equal[1, Indeterminate]",                    "False");
    check_compare_yields("Equal[Indeterminate, x]",                    "False");
    check_compare_yields("Unequal[Indeterminate, 1]",                  "True");
    check_compare_yields("Unequal[Indeterminate, x]",                  "True");

    /* Every ordering predicate is False, including the reflexive ones. */
    check_compare_yields("Less[Indeterminate, Indeterminate]",         "False");
    check_compare_yields("LessEqual[Indeterminate, Indeterminate]",    "False");
    check_compare_yields("Greater[Indeterminate, Indeterminate]",      "False");
    check_compare_yields("GreaterEqual[Indeterminate, Indeterminate]", "False");
    check_compare_yields("Less[Indeterminate, 1]",                     "False");
    check_compare_yields("Greater[1, Indeterminate]",                  "False");
    check_compare_yields("LessEqual[Indeterminate, 1]",                "False");
    check_compare_yields("GreaterEqual[Indeterminate, 1]",             "False");

    /* Derived Indeterminates (0/0, Infinity - Infinity) behave identically. */
    check_compare_yields("Equal[0/0, 0/0]",                            "False");
    check_compare_yields("Equal[Infinity - Infinity, Indeterminate]",  "False");

    /* Chained comparisons route through Inequality: one NaN kills the chain. */
    check_compare_yields("1 < Indeterminate < 3",                      "False");
    check_compare_yields("1 < 2 < Indeterminate",                      "False");
    check_compare_yields("Indeterminate == Indeterminate == Indeterminate",
                                                                       "False");
    check_compare_yields("1 == 1 == Indeterminate",                    "False");

    /* Variadic Unequal is pairwise: a NaN pair is unequal, but an equal
     * non-NaN pair still decides the whole call False. */
    check_compare_yields("Unequal[Indeterminate, Indeterminate, 1]",   "True");
    check_compare_yields("Unequal[1, 1, Indeterminate]",               "False");

    /* Structural identity is a different question and keeps its answer. */
    check_compare_yields("SameQ[Indeterminate, Indeterminate]",        "True");
    check_compare_yields("UnsameQ[Indeterminate, Indeterminate]",      "False");

    /* Nothing to compare in the 0- and 1-argument forms; True by convention. */
    check_compare_yields("Equal[Indeterminate]",                       "True");
    check_compare_yields("Unequal[Indeterminate]",                     "True");

    /* Ordinary comparisons are untouched. */
    check_compare_yields("Equal[1, 1]",                                "True");
    check_compare_yields("Equal[x, x]",                                "True");
    check_compare_yields("Unequal[1, 2]",                              "True");
    check_compare_yields("1 < 2 < 3",                                  "True");
}

void test_unequal_basic() {
    Expr* one = expr_new_integer(1);
    Expr* two = expr_new_integer(2);
    Expr* args[] = {one, two};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 2);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "True");
    expr_free(unequal); expr_free(res);
}

void test_unequal_false() {
    Expr* one = expr_new_integer(1);
    Expr* one_point_zero = expr_new_real(1.0);
    Expr* args[] = {one, one_point_zero};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 2);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res->data.symbol.name, "False");
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
    ASSERT_STR_EQ(res->data.symbol.name, "False");
    expr_free(unequal); expr_free(res);
}

void test_unequal_unevaluated() {
    Expr* one = expr_new_integer(1);
    Expr* x = expr_new_symbol("x");
    Expr* args[] = {one, x};
    Expr* unequal = expr_new_function(expr_new_symbol("Unequal"), args, 2);
    Expr* res = evaluate(unequal);
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(res->data.function.head->data.symbol.name, "Unequal");
    expr_free(unequal); expr_free(res);
}

/* Parse src, evaluate, assert the result is the bare symbol `want`. */
static void expect_sym(const char* src, const char* want) {
    Expr* e = parse_expression(src);
    Expr* r = evaluate(e);
    ASSERT(r && r->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(r->data.symbol.name, want);
    expr_free(e); expr_free(r);
}
/* Parse src, evaluate, assert the result stays a call with head `head`. */
static void expect_head(const char* src, const char* head) {
    Expr* e = parse_expression(src);
    Expr* r = evaluate(e);
    ASSERT(r && r->type == EXPR_FUNCTION);
    ASSERT_STR_EQ(r->data.function.head->data.symbol.name, head);
    expr_free(e); expr_free(r);
}

/* Equal/Unequal decide closed-form numeric constants that carry no real
 * ordering (I, 2 I Pi, complex radicals) via the exact zero-test -- matching
 * Mathematica -- while free symbols stay unevaluated. */
void test_equal_unequal_numeric_constants() {
    expect_sym("I == 0", "False");
    expect_sym("I != 0", "True");
    expect_sym("2 I Pi != 0", "True");
    expect_sym("0 != 2 I Pi", "True");
    expect_sym("(1 + I) == (1 - I)", "False");
    expect_sym("1/2 (-1 - I Sqrt[3]) != 1", "True");
    expect_sym("(-1 + 1/8 (-1 - I Sqrt[3])^3)/(-1 + 1/2 (-1 - I Sqrt[3])) == 0", "True");
    /* free symbols must NOT be decided */
    expect_head("x == 0", "Equal");
    expect_head("x != 0", "Unequal");
    /* existing exact identities and machine reals unaffected */
    expect_sym("Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2] == 0", "True");
    expect_sym("3 != 5", "True");
}

void test_not_basic() {
    Expr* t = expr_new_symbol("True");
    Expr* args_t[] = {t};
    Expr* not_t = expr_new_function(expr_new_symbol("Not"), args_t, 1);
    Expr* res_t = evaluate(not_t);
    ASSERT(res_t->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_t->data.symbol.name, "False");
    expr_free(not_t); expr_free(res_t);

    Expr* f = expr_new_symbol("False");
    Expr* args_f[] = {f};
    Expr* not_f = expr_new_function(expr_new_symbol("Not"), args_f, 1);
    Expr* res_f = evaluate(not_f);
    ASSERT(res_f->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_f->data.symbol.name, "True");
    expr_free(not_f); expr_free(res_f);
}

void test_and_basic() {
    Expr* and_empty = expr_new_function(expr_new_symbol("And"), NULL, 0);
    Expr* res_empty = evaluate(and_empty);
    ASSERT(res_empty->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_empty->data.symbol.name, "True");
    expr_free(and_empty); expr_free(res_empty);

    Expr* t = expr_new_symbol("True");
    Expr* f = expr_new_symbol("False");
    Expr* args1[] = {expr_copy(t), expr_copy(f)};
    Expr* and1 = expr_new_function(expr_new_symbol("And"), args1, 2);
    Expr* res1 = evaluate(and1);
    ASSERT(res1->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res1->data.symbol.name, "False");
    expr_free(and1); expr_free(res1);

    Expr* args2[] = {expr_copy(t), expr_copy(t)};
    Expr* and2 = expr_new_function(expr_new_symbol("And"), args2, 2);
    Expr* res2 = evaluate(and2);
    ASSERT(res2->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res2->data.symbol.name, "True");
    expr_free(and2); expr_free(res2);
    
    expr_free(t); expr_free(f);
}

void test_or_basic() {
    Expr* or_empty = expr_new_function(expr_new_symbol("Or"), NULL, 0);
    Expr* res_empty = evaluate(or_empty);
    ASSERT(res_empty->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res_empty->data.symbol.name, "False");
    expr_free(or_empty); expr_free(res_empty);

    Expr* t = expr_new_symbol("True");
    Expr* f = expr_new_symbol("False");
    Expr* args1[] = {expr_copy(f), expr_copy(t)};
    Expr* or1 = expr_new_function(expr_new_symbol("Or"), args1, 2);
    Expr* res1 = evaluate(or1);
    ASSERT(res1->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res1->data.symbol.name, "True");
    expr_free(or1); expr_free(res1);

    Expr* args2[] = {expr_copy(f), expr_copy(f)};
    Expr* or2 = expr_new_function(expr_new_symbol("Or"), args2, 2);
    Expr* res2 = evaluate(or2);
    ASSERT(res2->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(res2->data.symbol.name, "False");
    expr_free(or2); expr_free(res2);
    
    expr_free(t); expr_free(f);
}

/* `a === b === c` is SameQ[a, b, c], not SameQ[SameQ[a, b], c]: the nested
 * form compared the Boolean True with c, so 1 === 1 === 1 answered False.
 * The same held for UnsameQ and Unequal ("all distinct"). */
void test_variadic_comparison_chains() {
    assert_eval_eq("1 === 1 === 1", "True", 0);
    assert_eval_eq("1 === 1 === 2", "False", 0);
    assert_eval_eq("x === x === x", "True", 0);
    assert_eval_eq("1 =!= 2 =!= 3", "True", 0);
    assert_eval_eq("1 =!= 2 =!= 1", "False", 0);
    assert_eval_eq("1 != 2 != 3", "True", 0);
    assert_eval_eq("1 != 2 != 1", "False", 0);
    /* A parenthesised operand is a Boolean, compared as such. */
    assert_eval_eq("(1 === 1) === 1", "False", 0);
    assert_eval_eq("(1 === 1) === True", "True", 0);
    /* Unequal mixed into an ordering chain joins its Inequality. */
    assert_eval_eq("1 < 2 != 3", "True", 0);
    assert_eval_eq("1 < 2 != 2", "False", 0);
    assert_eval_eq("1 != 2 < 3", "True", 0);
    assert_eval_eq("Inequality[1, Unequal, 2, Less, 3]", "True", 0);
    assert_eval_eq("Inequality[1, Less, x, Unequal, 3]", "Inequality[1, Less, x, Unequal, 3]", 1);
}

/* Equal on Lists applies Equal elementwise, recursively, so the machine
 * tolerance that makes 0.1 + 0.2 == 0.3 True holds inside lists too. */
void test_equal_lists_tolerance() {
    assert_eval_eq("0.1 + 0.2 == 0.3", "True", 0);
    assert_eval_eq("{0.1 + 0.2} == {0.3}", "True", 0);
    assert_eval_eq("{0.1 + 0.2, 1.} == {0.3, 1.}", "True", 0);
    assert_eval_eq("{{0.1 + 0.2}, {1}} == {{0.3}, {1.}}", "True", 0);
    assert_eval_eq("{1, 2} == {1, 3}", "False", 0);
    assert_eval_eq("{1, 2} == {1, 2, 3}", "False", 0);
    assert_eval_eq("{x, 1} == {y, 2}", "False", 0);
    assert_eval_eq("{x} == {y}", "Equal[List[x], List[y]]", 1);
    assert_eval_eq("{x, 1} == {x, 1.}", "True", 0);
    assert_eval_eq("{Sqrt[2]^2} == {2}", "True", 0);
    assert_eval_eq("{1, 2} != {1, 3}", "True", 0);
    assert_eval_eq("{0.1 + 0.2} != {0.3}", "False", 0);
    assert_eval_eq("{x} != {y}", "Unequal[List[x], List[y]]", 1);
    /* ...and inside a chain, which decides each == pair with Equal itself. */
    assert_eval_eq("{0.1 + 0.2} == {0.3} == {0.3}", "True", 0);
}

/* Approximate numbers are Equal within the uncertainty of the less accurate
 * one, ignoring its last 7 bits: |x - y| <= 2^(7 - a), a = the least of the
 * operands' accuracies in bits (precision bits - log2|x|; exact operands
 * have none to spend). N[x, p] carries ceil(p log2 10) bits, so at 1 the
 * window is 2^(7 - ceil(p log2 10)), between 68 10^-p (p = 50) and
 * 101 10^-p (p = 30): an offset of 50 10^-p is inside it at every p
 * below, 200 10^-p outside. An arbitrary-precision operand used to be read
 * as a double and compared at the machine tolerance (relative 2^-46), so
 * N[1 - 10^-20, 30] == 1 was True and N[1 - 10^-20, 30] < 1 False. */
static void test_equal_approximate_tolerance(void) {
    /* the cases that were wrong, and the two that were already right */
    assert_eval_eq("d = N[1 - 1192/10^30, 30]; {d < 1, d == 1, d <= 1, d != 1}",
                   "List[True, False, True, True]", 1);
    assert_eval_eq("e = N[1 - 10^-20, 30]; {e < 1, e == 1, e > 1, e != 1}",
                   "List[True, False, False, True]", 1);
    assert_eval_eq("f = N[1 - 10^-29, 30]; {f < 1, f == 1, f <= 1, f >= 1, f != 1}",
                   "List[False, True, True, True, False]", 1);
    assert_eval_eq("{1.0 - 10.^-15 < 1, 1.0 - 10.^-15 == 1}", "List[False, True]", 1);
    /* both sides of the window, at 20, 30, 50 and 100 digits: inside it Less
     * is False and LessEqual True, outside it the sign decides */
    assert_eval_eq("Union[Table[{N[1 - 50/10^p, p] == 1, N[1 + 50/10^p, p] == 1, "
                   "N[1 - 50/10^p, p] < 1, N[1 - 50/10^p, p] <= 1, N[1 + 50/10^p, p] > 1, "
                   "N[1 + 50/10^p, p] >= 1}, {p, {20, 30, 50, 100}}]]",
                   "List[List[True, True, False, True, False, True]]", 1);
    assert_eval_eq("Union[Table[{N[1 - 200/10^p, p] == 1, N[1 + 200/10^p, p] == 1, "
                   "N[1 - 200/10^p, p] < 1, N[1 + 200/10^p, p] > 1, N[1 - 200/10^p, p] != 1}, "
                   "{p, {20, 30, 50, 100}}]]",
                   "List[List[False, False, True, True, True]]", 1);
    /* negative numbers */
    assert_eval_eq("{N[-1 - 50/10^30, 30] == -1, N[-1 - 50/10^30, 30] < -1, "
                   "N[-1 - 200/10^30, 30] == -1, N[-1 - 200/10^30, 30] < -1, N[-1 + 200/10^30, 30] > -1}",
                   "List[True, False, False, True, True]", 1);
    /* mixed precisions: the lower one governs */
    assert_eval_eq("{N[1 - 50/10^30, 50] == N[1, 30], N[1 - 50/10^30, 30] == N[1, 50], "
                   "N[1 - 200/10^30, 50] == N[1, 30], N[1 - 50/10^30, 50] == N[1, 50], "
                   "N[1 - 50/10^30, 50] < N[1, 30], N[1 - 50/10^30, 50] < N[1, 50]}",
                   "List[True, True, False, False, False, True]", 1);
    /* against a machine real the machine rule (relative 2^-46 = 1.42 10^-14) governs */
    assert_eval_eq("{N[1 + 10^-14, 30] == 1., N[1 + 2 10^-14, 30] == 1., "
                   "N[1 + 10^-14, 30] > 1., N[1 + 2 10^-14, 30] > 1.}",
                   "List[True, False, False, True]", 1);
    /* an approximate zero compares by its accuracy (~30 digits here) alone */
    assert_eval_eq("z = N[1, 30] - N[1, 30]; {z == 0, z == 10^-31, z == 10^-29, "
                   "z < 10^-29, z < 10^-31, z == -10^-31}",
                   "List[True, True, False, True, False, True]", 1);
    /* an exact cancellation knows the zero as well as the operands knew their
     * values: N[1/34, 30] - N[1/34, 30] is ~2^-105 uncertain (accuracy above
     * 31 digits), N[10^10, 30] - N[10^10, 30] ~2^-67 (below 20 digits) */
    assert_eval_eq("a34 = N[1/34, 30]; b10 = N[10^10, 30]; {a34 - a34 < 10^-31, "
                   "a34 - N[1/34, 40] < 10^-31, b10 - b10 < 10^-21, b10 - b10 < 10^-19, "
                   "31 < Accuracy[a34 - a34] < 32, 19 < Accuracy[b10 - b10] < 20}",
                   "List[True, True, False, True, True, True]", 1);
    /* a small nonzero number is not an approximate zero: relative tolerance */
    assert_eval_eq("{N[10^-40, 30] == 0, N[10^-40, 30] > 0, "
                   "N[10^-40, 30] == N[10^-40 (1 + 50/10^30), 30], "
                   "N[10^-40, 30] == N[10^-40 (1 + 200/10^30), 30]}",
                   "List[False, True, True, False]", 1);
    /* complex numbers: part by part (Equal only) */
    assert_eval_eq("{N[1 + I + 50/10^30, 30] == 1 + I, N[1 + I + 200/10^30, 30] == 1 + I, "
                   "N[1 + I (1 + 50/10^30), 30] == 1 + I, N[1 + I (1 + 200/10^30), 30] == 1 + I, "
                   "N[1 + I (1 + 200/10^30), 30] != 1 + I, N[2 I, 30] == 2 I}",
                   "List[True, False, True, False, True, True]", 1);
    /* exact closed forms compare at the approximate operand's precision */
    assert_eval_eq("{N[Pi, 30] == Pi, N[Pi + 10^-26, 30] == Pi, N[Pi + 10^-26, 30] > Pi, "
                   "N[Sqrt[2], 40] == Sqrt[2], N[Sqrt[2], 40] < Sqrt[2], N[E, 50] >= E}",
                   "List[True, False, True, True, False, True]", 1);
    /* the tolerance is not transitive: every pair of a chain must be equal */
    assert_eval_eq("a = N[1, 30]; b = N[1 - 90/10^30, 30]; c = N[1 - 180/10^30, 30]; "
                   "{a == b, b == c, a == c, a == b == c, Equal[a, b, c], a >= b >= c, a > b, Unequal[a, b, c]}",
                   "List[True, True, False, False, False, True, False, False]", 1);
    assert_eval_eq("{1. == 1. + 10.^-14 == 1. + 2 10.^-14, 0.1 + 0.2 == 0.3}",
                   "List[False, True]", 1);
    /* beyond the double range the same relative rule holds */
    assert_eval_eq("{N[10^3000, 30] == N[10^3000 (1 + 10^-29), 30], "
                   "N[10^3000, 30] == N[10^3000 (1 + 10^-25), 30], "
                   "N[10^3000, 30] < N[10^3000 (1 + 10^-29), 30], "
                   "N[10^3000, 30] < N[10^3000 (1 + 10^-25), 30]}",
                   "List[True, False, False, True]", 1);
}

/* The real infinities order against reals and each other (reference values).
 * FeynCalc's MaxIterations -> Infinity loop guard counter < Infinity stayed
 * unevaluated, so NestWhileList stopped after one element. */
static void test_infinity_ordering(void) {
    assert_eval_eq("{0 < Infinity, 1/2 < Infinity, 0. < Infinity, Pi < Infinity, Infinity > 0, -Infinity < 0}",
                   "List[True, True, True, True, True, True]", 1);
    assert_eval_eq("{Infinity < Infinity, Infinity <= Infinity, -Infinity < Infinity, 1 < -Infinity, Infinity > Infinity}",
                   "List[False, True, True, False, False]", 1);
    assert_eval_eq("{1 < 2 < Infinity, 1 < Infinity < 2, Infinity >= -Infinity, -Infinity <= -Infinity, N[Infinity] < 1}",
                   "List[True, False, True, True, False]", 1);
    assert_eval_eq("{DirectedInfinity[1] > 5, DirectedInfinity[-1] < -5}", "List[True, True]", 1);
    /* not on the real line, or not numeric: undecided */
    assert_eval_eq("Quiet[{x < Infinity, I < Infinity, ComplexInfinity < 1, DirectedInfinity[I] < 1}]",
                   "List[Less[x, Infinity], Less[Complex[0, 1], Infinity], "
                   "Less[ComplexInfinity, 1], Less[DirectedInfinity[Complex[0, 1]], 1]]", 1);
    assert_eval_eq("cinf = 0; Length[NestWhileList[(cinf++; If[cinf > 3, {}, {#}]) &, {{}, p}, (# =!= {}) && (cinf < Infinity) &]]",
                   "5", 1);
}

/* Exact quantities closer than a double resolves, or outside the double
 * range, compare by their exact values.  Expected values from first
 * principles: Pi/3^700 ~ 3.3 10^-334 > 0; Pi > E; Exp[1 + 10^-20] - E =
 * E (e^10^-20 - 1) > 0; Pi + 10^-20 > Pi. */
static void test_compare_exact_quantities_precision_and_range(void) {
    assert_eval_eq("{Pi/3^700 == 0, Pi/3^700 != 0, Pi/3^700 > 0, Pi/3^700 < 0, -Pi/3^700 < 0}",
                   "List[False, True, True, False, True]", 1);
    assert_eval_eq("{I Pi/3^700 == 0, E/10^400 == 0, Pi 10^500 == 0, Zeta[3] 10^1000 == 0}",
                   "List[False, False, False, False]", 1);
    assert_eval_eq("{Pi/3^700 < E/3^700, Pi 10^500 > E 10^500, Pi/3^700 > 1/3^701}",
                   "List[False, True, True]", 1);
    assert_eval_eq("{Pi + 10^-20 == Pi, Exp[1 + 10^-20] == E, Exp[1 + 10^-20] != E, "
                   "Sqrt[2] + 10^-30 == Sqrt[2], Exp[1 + 10^-20] > E, Pi + 10^-20 > Pi}",
                   "List[False, False, True, False, True, True]", 1);
    assert_eval_eq("Exp[-10^4] - Exp[-10^4 - 10^-20] == 0", "False", 1);
    /* equal exact values stay equal */
    assert_eval_eq("{Sqrt[5 + 2 Sqrt[6]] == Sqrt[3] + Sqrt[2], Sin[Pi/7]^2 + Cos[Pi/7]^2 == 1, "
                   "Sqrt[5 + 2 Sqrt[6]] <= Sqrt[3] + Sqrt[2], Sqrt[5 + 2 Sqrt[6]] < Sqrt[3] + Sqrt[2]}",
                   "List[True, True, True, False]", 1);
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
    TEST(test_equal_unequal_numeric_constants);
    TEST(test_less_basic);
    TEST(test_less_false);
    TEST(test_less_rational_real);
    TEST(test_less_exact_rational);
    TEST(test_lessequal_mixed);
    TEST(test_greater_unevaluated);
    TEST(test_compare_bigint_adjacent);
    TEST(test_compare_rational_beyond_double_range);
    TEST(test_compare_exact_quantities_precision_and_range);
    TEST(test_indeterminate_ieee_unordered);
    TEST(test_variadic_comparison_chains);
    TEST(test_equal_lists_tolerance);
    TEST(test_equal_approximate_tolerance);
    TEST(test_not_basic);
    TEST(test_and_basic);
    TEST(test_or_basic);
    
    printf("All comparisons tests passed!\n");
    return 0;
}
