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
#include <math.h>

static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "FixedPointList %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* Compare each element of a List result element-wise, allowing tolerance on reals. */
static void run_real_list(const char* input, const double* expected, size_t n, double tol) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    ASSERT_MSG(r->type == EXPR_FUNCTION &&
               r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "List") == 0,
               "FixedPointList %s: result not a List", input);
    ASSERT_MSG(r->data.function.arg_count == n,
               "FixedPointList %s: expected length %zu, got %zu",
               input, n, r->data.function.arg_count);
    for (size_t i = 0; i < n; i++) {
        Expr* el = r->data.function.args[i];
        double got = 0.0;
        if (el->type == EXPR_REAL) got = el->data.real;
        else if (el->type == EXPR_INTEGER) got = (double)el->data.integer;
        else {
            char* s = expr_to_string(el);
            ASSERT_MSG(0, "FixedPointList %s: element %zu not numeric: %s", input, i, s);
            free(s);
        }
        ASSERT_MSG(fabs(got - expected[i]) <= tol,
                   "FixedPointList %s: element %zu expected ~%g, got %g",
                   input, i, expected[i], got);
    }
    expr_free(e);
    expr_free(r);
}

/* ---------- Spec examples ---------- */

static void test_fpl_halving_to_two(void) {
    /* FixedPointList[1+Floor[#/2]&,1000]
       -> {1000, 501, 251, 126, 64, 33, 17, 9, 5, 3, 2, 2} */
    run_full("FixedPointList[1 + Floor[#/2] &, 1000]",
             "List[1000, 501, 251, 126, 64, 33, 17, 9, 5, 3, 2, 2]");
}

static void test_fpl_gcd_via_rule(void) {
    /* FixedPointList[#/.{a_,b_}/;b!=0->{b,Mod[a,b]}&,{28,21}]
       -> {{28,21},{21,7},{7,0},{7,0}} */
    run_full("FixedPointList[# /. {a_, b_} /; b != 0 -> {b, Mod[a, b]} &, {28, 21}]",
             "List[List[28, 21], List[21, 7], List[7, 0], List[7, 0]]");
}

static void test_fpl_max_steps(void) {
    /* FixedPointList[1+Floor[#/2]&,1000,5]
       -> {1000, 501, 251, 126, 64, 33}  -- last two need NOT be equal */
    run_full("FixedPointList[1 + Floor[#/2] &, 1000, 5]",
             "List[1000, 501, 251, 126, 64, 33]");
}

static void test_fpl_newton_sqrt2(void) {
    /* FixedPointList[(# + 2/#)/2 &, 1.0]
       -> {1., 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421} */
    double expected[] = {1.0, 1.5, 1.41666666666666, 1.41421568627451,
                         1.41421356237469, 1.41421356237310, 1.41421356237310};
    run_real_list("FixedPointList[(# + 2/#)/2 &, 1.0]", expected, 7, 1e-10);
}

/* ---------- Additional behavior ---------- */

static void test_fpl_already_fixed(void) {
    /* If f[expr] == expr immediately, list is {expr, expr} */
    run_full("FixedPointList[# &, 5]", "List[5, 5]");
}

static void test_fpl_max_zero(void) {
    /* n = 0 means no applications; list is just {expr} */
    run_full("FixedPointList[1 + Floor[#/2] &, 1000, 0]", "List[1000]");
}

static void test_fpl_max_one(void) {
    /* n = 1 means at most one application of f */
    run_full("FixedPointList[1 + Floor[#/2] &, 1000, 1]", "List[1000, 501]");
}

static void test_fpl_max_infinity(void) {
    /* Infinity as max bound is allowed and equivalent to omitting it */
    run_full("FixedPointList[1 + Floor[#/2] &, 1000, Infinity]",
             "List[1000, 501, 251, 126, 64, 33, 17, 9, 5, 3, 2, 2]");
}

static void test_fpl_symbolic_no_fixed_point_with_bound(void) {
    /* Symbolic input with a bound: f stays unevaluated, no fixed point reached */
    run_full("FixedPointList[f, x, 3]",
             "List[x, f[x], f[f[x]], f[f[f[x]]]]");
}

static void test_fpl_floor_to_zero(void) {
    /* Successive halving converges to 0 */
    run_full("FixedPointList[Floor[#/2] &, 8]",
             "List[8, 4, 2, 1, 0, 0]");
}

static void test_fpl_equivalent_to_nestwhilelist(void) {
    /* FixedPointList[f, expr] == NestWhileList[f, expr, UnsameQ, 2] */
    run_full("FixedPointList[Floor[#/2] &, 10]",
             "List[10, 5, 2, 1, 0, 0]");
    run_full("NestWhileList[Floor[#/2] &, 10, UnsameQ, 2]",
             "List[10, 5, 2, 1, 0, 0]");
}

static void test_fpl_sametest_tolerance(void) {
    /* SameTest -> custom predicate: stop when consecutive iterates differ by < 0.01 */
    /* Newton iteration for sqrt(2). With tolerance, stop earlier. */
    Expr* e = parse_expression("FixedPointList[(# + 2/#)/2 &, 1.0, SameTest -> (Abs[#1 - #2] < 0.01 &)]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           r->data.function.head->type == EXPR_SYMBOL &&
           strcmp(r->data.function.head->data.symbol.name, "List") == 0);
    /* {1., 1.5, 1.41667, 1.41422} -- stop at step 3 since |1.41422 - 1.41667| < 0.01 */
    ASSERT_MSG(r->data.function.arg_count == 4,
               "expected 4 elements, got %zu", r->data.function.arg_count);
    expr_free(e);
    expr_free(r);
}

static void test_fpl_sametest_with_max(void) {
    /* SameTest combined with explicit n limit */
    run_full("FixedPointList[(# + 2/#)/2 &, 1.0, 2, SameTest -> (Abs[#1 - #2] < 0.0000001 &)]",
             "List[1.0, 1.5, 1.41667]");
}

static void test_fpl_sametest_always_true(void) {
    /* SameTest -> True& stops after first application (last two trivially "same") */
    run_full("FixedPointList[# + 1 &, 0, SameTest -> (True &)]", "List[0, 1]");
}

static void test_fpl_sametest_never_true(void) {
    /* SameTest -> False& with bounded n iterates until n hit */
    run_full("FixedPointList[# + 1 &, 0, 4, SameTest -> (False &)]",
             "List[0, 1, 2, 3, 4]");
}

static void test_fpl_collatz_one_loop(void) {
    /* Collatz sequence eventually hits cycle 1 -> 4 -> 2 -> 1, with x=1 the sequence
       1 -> 4 -> 2 -> 1 -> 4 -> ... never has consecutive equal values, so we need a bound.
       Use a simple monotone case instead: f = If[# > 0, # - 1, 0]. */
    run_full("FixedPointList[If[# > 0, # - 1, 0] &, 3]",
             "List[3, 2, 1, 0, 0]");
}

/* ---------- Unevaluated cases ---------- */

static void test_fpl_wrong_arg_count(void) {
    run_full("FixedPointList[f]", "FixedPointList[f]");
    run_full("FixedPointList[]", "FixedPointList[]");
}

static void test_fpl_bad_max(void) {
    run_full("FixedPointList[f, x, -1]", "FixedPointList[f, x, -1]");
    run_full("FixedPointList[f, x, y]", "FixedPointList[f, x, y]");
}

static void test_fpl_bad_option(void) {
    /* Unknown option-shaped argument */
    run_full("FixedPointList[f, x, Method -> Newton]",
             "FixedPointList[f, x, Rule[Method, Newton]]");
}

static void test_fpl_duplicate_sametest(void) {
    /* Two SameTest options is an error */
    run_full("FixedPointList[f, x, SameTest -> Equal, SameTest -> SameQ]",
             "FixedPointList[f, x, Rule[SameTest, Equal], Rule[SameTest, SameQ]]");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_fpl_halving_to_two);
    TEST(test_fpl_gcd_via_rule);
    TEST(test_fpl_max_steps);
    TEST(test_fpl_newton_sqrt2);

    /* Additional behavior */
    TEST(test_fpl_already_fixed);
    TEST(test_fpl_max_zero);
    TEST(test_fpl_max_one);
    TEST(test_fpl_max_infinity);
    TEST(test_fpl_symbolic_no_fixed_point_with_bound);
    TEST(test_fpl_floor_to_zero);
    TEST(test_fpl_equivalent_to_nestwhilelist);
    TEST(test_fpl_sametest_tolerance);
    TEST(test_fpl_sametest_with_max);
    TEST(test_fpl_sametest_always_true);
    TEST(test_fpl_sametest_never_true);
    TEST(test_fpl_collatz_one_loop);

    /* Unevaluated cases */
    TEST(test_fpl_wrong_arg_count);
    TEST(test_fpl_bad_max);
    TEST(test_fpl_bad_option);
    TEST(test_fpl_duplicate_sametest);

    printf("All FixedPointList tests passed!\n");
    return 0;
}
