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
    ASSERT_MSG(strcmp(s, expected) == 0, "FixedPoint %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_real(const char* input, double expected, double tol) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    double got = 0.0;
    if (r->type == EXPR_REAL) got = r->data.real;
    else if (r->type == EXPR_INTEGER) got = (double)r->data.integer;
    else {
        char* s = expr_to_string(r);
        ASSERT_MSG(0, "FixedPoint %s: result not numeric: %s", input, s);
        free(s);
    }
    ASSERT_MSG(fabs(got - expected) <= tol,
               "FixedPoint %s: expected ~%g, got %g", input, expected, got);
    expr_free(e);
    expr_free(r);
}

/* ---------- Spec examples ---------- */

static void test_fp_newton_sqrt2(void) {
    /* FixedPoint[(# + 2/#)/2 &, 1.] -> 1.41421... */
    run_real("FixedPoint[(# + 2/#)/2 &, 1.0]", 1.41421356237310, 1e-10);
}

static void test_fp_halving_to_two(void) {
    /* FixedPoint[1+Floor[#/2]&,1000] -> 2 */
    run_full("FixedPoint[1 + Floor[#/2] &, 1000]", "2");
}

static void test_fp_gcd_via_rule(void) {
    /* FixedPoint[#/.{a_,b_}/;b!=0->{b,Mod[a,b]}&,{28,21}] -> {7,0} */
    run_full("FixedPoint[# /. {a_, b_} /; b != 0 -> {b, Mod[a, b]} &, {28, 21}]",
             "List[7, 0]");
}

static void test_fp_max_steps_bounded(void) {
    /* FixedPoint[1+Floor[#/2]&,10^6,10] -> 978 */
    run_full("FixedPoint[1 + Floor[#/2] &, 10^6, 10]", "978");
}

static void test_fp_unbounded_to_two(void) {
    /* FixedPoint[1+Floor[#/2]&,10^6] -> 2 */
    run_full("FixedPoint[1 + Floor[#/2] &, 10^6]", "2");
}

static void test_fp_function_form(void) {
    /* FixedPoint[Function[x,(x+2/x)/2],1.0] -> 1.41421... */
    run_real("FixedPoint[Function[x, (x + 2/x)/2], 1.0]", 1.41421356237310, 1e-10);
}

/* ---------- Additional behavior ---------- */

static void test_fp_already_fixed(void) {
    /* If f[expr] == expr immediately, result is expr */
    run_full("FixedPoint[# &, 5]", "5");
}

static void test_fp_max_zero(void) {
    /* n = 0 means no applications; result is expr */
    run_full("FixedPoint[1 + Floor[#/2] &, 1000, 0]", "1000");
}

static void test_fp_max_one(void) {
    /* n = 1 means at most one application */
    run_full("FixedPoint[1 + Floor[#/2] &, 1000, 1]", "501");
}

static void test_fp_max_infinity(void) {
    /* Infinity max equivalent to omitting it */
    run_full("FixedPoint[1 + Floor[#/2] &, 1000, Infinity]", "2");
}

static void test_fp_symbolic_with_bound(void) {
    /* Symbolic input never converges; bound returns f^n[x] */
    run_full("FixedPoint[f, x, 3]", "f[f[f[x]]]");
}

static void test_fp_floor_to_zero(void) {
    /* Successive halving converges to 0 */
    run_full("FixedPoint[Floor[#/2] &, 8]", "0");
}

static void test_fp_equivalent_to_last_fpl(void) {
    /* FixedPoint[f, expr] == Last[FixedPointList[f, expr]] */
    run_full("FixedPoint[Floor[#/2] &, 10]", "0");
    run_full("Last[FixedPointList[Floor[#/2] &, 10]]", "0");
}

static void test_fp_sametest_tolerance(void) {
    /* Stop when |#1 - #2| < 0.01 — Newton sqrt(2) stops at ~1.41422 */
    Expr* e = parse_expression("FixedPoint[(# + 2/#)/2 &, 1.0, SameTest -> (Abs[#1 - #2] < 0.01 &)]");
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_REAL);
    ASSERT_MSG(fabs(r->data.real - 1.41421568627451) <= 1e-8,
               "expected ~1.41422, got %g", r->data.real);
    expr_free(e);
    expr_free(r);
}

static void test_fp_sametest_with_max(void) {
    /* SameTest combined with explicit n limit; converges in 2 steps */
    run_full("FixedPoint[(# + 2/#)/2 &, 1.0, 2, SameTest -> (Abs[#1 - #2] < 0.0000001 &)]",
             "1.41667");
}

static void test_fp_sametest_always_true(void) {
    /* SameTest -> True& stops after first application */
    run_full("FixedPoint[# + 1 &, 0, SameTest -> (True &)]", "1");
}

static void test_fp_sametest_never_true(void) {
    /* SameTest -> False& with bounded n iterates until n hit */
    run_full("FixedPoint[# + 1 &, 0, 4, SameTest -> (False &)]", "4");
}

static void test_fp_monotone_decrement(void) {
    /* f[#] = If[# > 0, # - 1, 0] converges to 0 */
    run_full("FixedPoint[If[# > 0, # - 1, 0] &, 3]", "0");
}

/* ---------- Throw exit ---------- */

static void test_fp_throw_exits(void) {
    /* Throw[#] inside f exits FixedPoint early, propagating the Throw. With no
     * enclosing Catch the uncaught Throw returns to top level wrapped in Hold
     * (matching Wolfram's Throw::nocatch behaviour: Hold[Throw[6]]). */
    run_full("FixedPoint[If[# > 5, Throw[#], # + 1] &, 0]", "Hold[Throw[6]]");
}

/* ---------- Unevaluated cases ---------- */

static void test_fp_wrong_arg_count(void) {
    run_full("FixedPoint[f]", "FixedPoint[f]");
    run_full("FixedPoint[]", "FixedPoint[]");
}

static void test_fp_bad_max(void) {
    run_full("FixedPoint[f, x, -1]", "FixedPoint[f, x, -1]");
    run_full("FixedPoint[f, x, y]", "FixedPoint[f, x, y]");
}

static void test_fp_bad_option(void) {
    run_full("FixedPoint[f, x, Method -> Newton]",
             "FixedPoint[f, x, Rule[Method, Newton]]");
}

static void test_fp_duplicate_sametest(void) {
    run_full("FixedPoint[f, x, SameTest -> Equal, SameTest -> SameQ]",
             "FixedPoint[f, x, Rule[SameTest, Equal], Rule[SameTest, SameQ]]");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_fp_newton_sqrt2);
    TEST(test_fp_halving_to_two);
    TEST(test_fp_gcd_via_rule);
    TEST(test_fp_max_steps_bounded);
    TEST(test_fp_unbounded_to_two);
    TEST(test_fp_function_form);

    /* Additional behavior */
    TEST(test_fp_already_fixed);
    TEST(test_fp_max_zero);
    TEST(test_fp_max_one);
    TEST(test_fp_max_infinity);
    TEST(test_fp_symbolic_with_bound);
    TEST(test_fp_floor_to_zero);
    TEST(test_fp_equivalent_to_last_fpl);
    TEST(test_fp_sametest_tolerance);
    TEST(test_fp_sametest_with_max);
    TEST(test_fp_sametest_always_true);
    TEST(test_fp_sametest_never_true);
    TEST(test_fp_monotone_decrement);

    /* Throw exit */
    TEST(test_fp_throw_exits);

    /* Unevaluated cases */
    TEST(test_fp_wrong_arg_count);
    TEST(test_fp_bad_max);
    TEST(test_fp_bad_option);
    TEST(test_fp_duplicate_sametest);

    printf("All FixedPoint tests passed!\n");
    return 0;
}
