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
    ASSERT_MSG(strcmp(s, expected) == 0, "NestWhileList %s: expected %s, got %s", input, expected, s);
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
               "NestWhileList %s: result not a List", input);
    ASSERT_MSG(r->data.function.arg_count == n,
               "NestWhileList %s: expected length %zu, got %zu",
               input, n, r->data.function.arg_count);
    for (size_t i = 0; i < n; i++) {
        Expr* el = r->data.function.args[i];
        double got = 0.0;
        if (el->type == EXPR_REAL) got = el->data.real;
        else if (el->type == EXPR_INTEGER) got = (double)el->data.integer;
        else {
            char* s = expr_to_string(el);
            ASSERT_MSG(0, "NestWhileList %s: element %zu not numeric: %s", input, i, s);
            free(s);
        }
        ASSERT_MSG(fabs(got - expected[i]) <= tol,
                   "NestWhileList %s: element %zu expected ~%g, got %g",
                   input, i, expected[i], got);
    }
    expr_free(e);
    expr_free(r);
}

/* ---------- Spec examples ---------- */

static void test_nwl_divide_by_two_until_odd(void) {
    /* NestWhileList[#/2 &, 123456, EvenQ] -> {123456,61728,30864,15432,7716,3858,1929} */
    run_full("NestWhileList[#/2 &, 123456, EvenQ]",
             "List[123456, 61728, 30864, 15432, 7716, 3858, 1929]");
}

static void test_nwl_log_until_nonpositive(void) {
    /* NestWhileList[Log, 100., # > 0 &] -> {100., 4.60517, 1.52718, 0.423423, -0.859384} */
    double expected[] = {100.0, 4.60517, 1.52718, 0.423423, -0.859384};
    run_real_list("NestWhileList[Log, 100., # > 0 &]", expected, 5, 1e-4);
}

static void test_nwl_max_iterations(void) {
    /* NestWhileList[Floor[#/2] &, 20, UnsameQ, 2, 4] -> {20,10,5,2,1} */
    run_full("NestWhileList[Floor[#/2] &, 20, UnsameQ, 2, 4]",
             "List[20, 10, 5, 2, 1]");
}

static void test_nwl_infinity_max(void) {
    /* NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity] -> {20,10,5,2,1} */
    run_full("NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity]",
             "List[20, 10, 5, 2, 1]");
}

static void test_nwl_extra_step(void) {
    /* NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1] -> {20,10,5,2,1,0} */
    run_full("NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]",
             "List[20, 10, 5, 2, 1, 0]");
}

static void test_nwl_drop_last(void) {
    /* NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1] -> {20,10,5,2} */
    run_full("NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]",
             "List[20, 10, 5, 2]");
}

static void test_nwl_next_prime(void) {
    /* NestWhileList[# + 1 &, 899, !PrimeQ[#] &] -> {899,...,907} */
    run_full("NestWhileList[# + 1 &, 899, !PrimeQ[#] &]",
             "List[899, 900, 901, 902, 903, 904, 905, 906, 907]");
}

static void test_nwl_mult_order(void) {
    /* NestWhileList[Mod[2 #, 19] &, 2, # != 1 &] -> sequence of powers of 2 mod 19 */
    run_full("NestWhileList[Mod[2 #, 19] &, 2, # != 1 &]",
             "List[2, 4, 8, 16, 13, 7, 14, 9, 18, 17, 15, 11, 3, 6, 12, 5, 10, 1]");
}

static void test_nwl_orbit_all(void) {
    /* NestWhileList[Mod[5 #, 7] &, 4, Unequal, All] -> {4,6,2,3,1,5,4} */
    run_full("NestWhileList[Mod[5 #, 7] &, 4, Unequal, All]",
             "List[4, 6, 2, 3, 1, 5, 4]");
}

static void test_nwl_collatz_all(void) {
    /* NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All]
       -> {400,200,...,1,2} (includes repeated 2) */
    run_full("NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All]",
             "List[400, 200, 100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1, 2]");
}

static void test_nwl_collatz_drop_repeat(void) {
    /* NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All, Infinity, -1]
       -> drops trailing repeat */
    run_full("NestWhileList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 400, Unequal, All, Infinity, -1]",
             "List[400, 200, 100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1]");
}

/* ---------- Additional behavior ---------- */

static void test_nwl_test_false_initially(void) {
    /* If test is never True, only expr is in the list */
    run_full("NestWhileList[f, x, False &]", "List[x]");
    run_full("NestWhileList[# + 1 &, 10, # < 5 &]", "List[10]");
}

static void test_nwl_symbolic_test(void) {
    /* Unevaluating test yields a symbolic expression, not True -> stop immediately */
    run_full("NestWhileList[f, x, g]", "List[x]");
}

static void test_nwl_max_zero(void) {
    /* max=0 means no applications; the list is just {expr} */
    run_full("NestWhileList[f, x, TrueQ, 1, 0]", "List[x]");
    run_full("NestWhileList[# + 1 &, 10, True &, 1, 0]", "List[10]");
}

static void test_nwl_fixedpointlist_equivalence(void) {
    /* NestWhileList[f, expr, UnsameQ, 2] == FixedPointList */
    run_full("NestWhileList[Floor[#/2] &, 8, UnsameQ, 2]",
             "List[8, 4, 2, 1, 0, 0]");
}

static void test_nwl_range_equivalence_all(void) {
    /* All is equivalent to {1, Infinity} */
    run_full("NestWhileList[Mod[5 #, 7] &, 4, Unequal, {1, Infinity}]",
             "List[4, 6, 2, 3, 1, 5, 4]");
}

static void test_nwl_range_min_delay(void) {
    /* {mmin, mmax} delays first test until mmin results exist */
    run_full("NestWhileList[# + 1 &, 0, False &, {3, 3}]",
             "List[0, 1, 2]");
}

static void test_nwl_extra_two(void) {
    /* Append 2 extra applications after loop ends */
    run_full("NestWhileList[# + 1 &, 0, # < 3 &, 1, Infinity, 2]",
             "List[0, 1, 2, 3, 4, 5]");
}

static void test_nwl_drop_two(void) {
    /* Drop 2 from end */
    run_full("NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -2]",
             "List[20, 10, 5]");
}

static void test_nwl_drop_all(void) {
    /* Drop more than length -> empty list */
    run_full("NestWhileList[# + 1 &, 0, # < 2 &, 1, Infinity, -10]",
             "List[]");
}

static void test_nwl_all_detects_any_repeat(void) {
    /* UnsameQ + All stops when any prior value repeats */
    run_full("NestWhileList[Mod[# + 3, 7] &, 0, UnsameQ, All]",
             "List[0, 3, 6, 2, 5, 1, 4, 0]");
}

static void test_nwl_pair_test(void) {
    /* 2-arg test compares adjacent pairs */
    run_full("NestWhileList[Floor[#/2] &, 10, UnsameQ, 2]",
             "List[10, 5, 2, 1, 0, 0]");
}

/* ---------- Unevaluated cases ---------- */

static void test_nwl_wrong_arg_count(void) {
    run_full("NestWhileList[f, x]", "NestWhileList[f, x]");
    run_full("NestWhileList[f, x, g, 1, 2, 3, 4]", "NestWhileList[f, x, g, 1, 2, 3, 4]");
}

static void test_nwl_bad_m(void) {
    run_full("NestWhileList[f, x, g, 0]", "NestWhileList[f, x, g, 0]");
    run_full("NestWhileList[f, x, g, -2]", "NestWhileList[f, x, g, -2]");
    run_full("NestWhileList[f, x, g, y]", "NestWhileList[f, x, g, y]");
}

static void test_nwl_bad_max(void) {
    run_full("NestWhileList[f, x, g, 1, -1]", "NestWhileList[f, x, g, 1, -1]");
    run_full("NestWhileList[f, x, g, 1, y]", "NestWhileList[f, x, g, 1, y]");
}

static void test_nwl_bad_n(void) {
    run_full("NestWhileList[f, x, g, 1, 2, y]", "NestWhileList[f, x, g, 1, 2, y]");
}

static void test_nwl_bad_range(void) {
    run_full("NestWhileList[f, x, g, {3, 2}]", "NestWhileList[f, x, g, List[3, 2]]");
    run_full("NestWhileList[f, x, g, {1, 2, 3}]", "NestWhileList[f, x, g, List[1, 2, 3]]");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_nwl_divide_by_two_until_odd);
    TEST(test_nwl_log_until_nonpositive);
    TEST(test_nwl_max_iterations);
    TEST(test_nwl_infinity_max);
    TEST(test_nwl_extra_step);
    TEST(test_nwl_drop_last);
    TEST(test_nwl_next_prime);
    TEST(test_nwl_mult_order);
    TEST(test_nwl_orbit_all);
    TEST(test_nwl_collatz_all);
    TEST(test_nwl_collatz_drop_repeat);

    /* Additional behavior */
    TEST(test_nwl_test_false_initially);
    TEST(test_nwl_symbolic_test);
    TEST(test_nwl_max_zero);
    TEST(test_nwl_fixedpointlist_equivalence);
    TEST(test_nwl_range_equivalence_all);
    TEST(test_nwl_range_min_delay);
    TEST(test_nwl_extra_two);
    TEST(test_nwl_drop_two);
    TEST(test_nwl_drop_all);
    TEST(test_nwl_all_detects_any_repeat);
    TEST(test_nwl_pair_test);

    /* Unevaluated cases */
    TEST(test_nwl_wrong_arg_count);
    TEST(test_nwl_bad_m);
    TEST(test_nwl_bad_max);
    TEST(test_nwl_bad_n);
    TEST(test_nwl_bad_range);

    printf("All NestWhileList tests passed!\n");
    return 0;
}
