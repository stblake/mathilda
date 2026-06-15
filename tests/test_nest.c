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
    ASSERT_MSG(strcmp(s, expected) == 0, "Nest %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "Nest %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* Evaluate and compare the result to an expected double within a tolerance. */
static void run_real_approx(const char* input, double expected, double tol) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    double got = 0.0;
    if (r->type == EXPR_REAL) {
        got = r->data.real;
    } else if (r->type == EXPR_INTEGER) {
        got = (double)r->data.integer;
    } else {
        char* s = expr_to_string(r);
        ASSERT_MSG(0, "Nest %s: expected real, got %s", input, s);
        free(s);
    }
    ASSERT_MSG(fabs(got - expected) <= tol,
               "Nest %s: expected ~%g, got %g", input, expected, got);
    expr_free(e);
    expr_free(r);
}

/* ---------- Examples from the spec ---------- */

static void test_nest_basic_symbolic(void) {
    /* Nest[f, x, 3] -> f[f[f[x]]] */
    run_full("Nest[f, x, 3]", "f[f[f[x]]]");
}

static void test_nest_pure_function_number(void) {
    /* Nest[(1+#)^2 &, 1, 3] -> 676 */
    run_full("Nest[(1 + #)^2 &, 1, 3]", "676");
}

static void test_nest_pure_function_symbolic(void) {
    /* Nest[(1+#)^2 &, x, 5] -> (1+(1+(1+(1+(1+x)^2)^2)^2)^2)^2 */
    run_infix("Nest[(1 + #)^2 &, x, 5]",
              "(1 + (1 + (1 + (1 + (1 + x)^2)^2)^2)^2)^2");
}

static void test_nest_sqrt_real(void) {
    /* Nest[Sqrt, 100.0, 4] -> ~1.33352 */
    run_real_approx("Nest[Sqrt, 100.0, 4]", 1.333521432, 1e-5);
}

static void test_nest_continued_fraction(void) {
    /* Nest[1/(1+#) &, x, 5] -> 1/(1+1/(1+1/(1+1/(1+1/(1+x))))) */
    run_infix("Nest[1/(1 + #) &, x, 5]",
              "1/(1 + 1/(1 + 1/(1 + 1/(1 + 1/(1 + x)))))");
}

static void test_nest_power_tower(void) {
    /* Nest[x^# &, x, 6] -> x^x^x^x^x^x^x */
    run_infix("Nest[x^# &, x, 6]", "x^x^x^x^x^x^x");
}

static void test_nest_compound_growth(void) {
    /* Nest[#(1 + 0.05) &, 1000, 10] -> ~1628.89 */
    run_real_approx("Nest[#(1 + 0.05) &, 1000, 10]", 1628.894627, 1e-3);
}

static void test_nest_newton_sqrt2(void) {
    /* Nest[(# + 2/#)/2 &, 1.0, 5] -> ~1.41421 */
    run_real_approx("Nest[(# + 2/#)/2 &, 1.0, 5]", 1.414213562, 1e-6);
}

static void test_nest_fibonacci_matrix(void) {
    /* Nest[{{1,1},{1,0}}.# &, {0,1}, 10] -> {55, 34} */
    run_full("Nest[{{1,1},{1,0}} . # &, {0,1}, 10]",
             "List[55, 34]");
}

/* ---------- Edge cases ---------- */

static void test_nest_zero(void) {
    /* Nest[f, x, 0] -> x */
    run_full("Nest[f, x, 0]", "x");
}

static void test_nest_zero_pure(void) {
    /* Pure function, n=0 returns expr unchanged */
    run_full("Nest[(# + 1) &, 7, 0]", "7");
}

static void test_nest_one(void) {
    /* Nest[f, x, 1] -> f[x] */
    run_full("Nest[f, x, 1]", "f[x]");
}

static void test_nest_evaluates_each_step(void) {
    /* Each step should evaluate. With numeric f, result reduces. */
    run_full("Nest[# + 1 &, 0, 10]", "10");
}

static void test_nest_with_symbol_head(void) {
    /* Nest with a named function that has a DownValue */
    run_full("Nest[Plus[#, 1] &, 5, 4]", "9");
}

static void test_nest_list_building(void) {
    /* Nest can build lists by appending */
    run_full("Nest[Append[#, x] &, {}, 3]",
             "List[x, x, x]");
}

static void test_nest_sqrt_symbolic(void) {
    /* Sqrt[a] reduces to Power[a, 1/2] internally, and nested rational powers
     * with |exponent| < 1 compose on the principal branch for any base, so
     * Sqrt[Sqrt[Sqrt[a]]] folds to a^(1/8). */
    run_full("Nest[Sqrt, a, 3]",
             "Power[a, Rational[1, 8]]");
}

static void test_nest_negative_n_unevaluated(void) {
    /* Negative n should leave Nest unevaluated */
    run_full("Nest[f, x, -1]", "Nest[f, x, -1]");
}

static void test_nest_non_integer_n_unevaluated(void) {
    /* Non-integer n should leave Nest unevaluated */
    run_full("Nest[f, x, n]", "Nest[f, x, n]");
}

static void test_nest_wrong_arg_count(void) {
    /* 2 args - unevaluated */
    run_full("Nest[f, x]", "Nest[f, x]");
}

static void test_nest_double_pure(void) {
    /* Nest[2# &, 1, 8] -> 256 */
    run_full("Nest[2 # &, 1, 8]", "256");
}

static void test_nest_factorial_like(void) {
    /* Simple chained multiplication via Nest */
    run_full("Nest[# * 2 &, 3, 4]", "48");
}

static void test_nest_list_rotation(void) {
    /* Rotating a list three times should yield the expected rotation */
    run_full("Nest[RotateLeft, {a, b, c, d}, 2]",
             "List[c, d, a, b]");
}

static void test_nest_head_changing(void) {
    /* Apply a function that changes the head */
    run_full("Nest[f, g[x], 2]", "f[f[g[x]]]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_nest_basic_symbolic);
    TEST(test_nest_pure_function_number);
    TEST(test_nest_pure_function_symbolic);
    TEST(test_nest_sqrt_real);
    TEST(test_nest_continued_fraction);
    TEST(test_nest_power_tower);
    TEST(test_nest_compound_growth);
    TEST(test_nest_newton_sqrt2);
    TEST(test_nest_fibonacci_matrix);

    TEST(test_nest_zero);
    TEST(test_nest_zero_pure);
    TEST(test_nest_one);
    TEST(test_nest_evaluates_each_step);
    TEST(test_nest_with_symbol_head);
    TEST(test_nest_list_building);
    TEST(test_nest_sqrt_symbolic);
    TEST(test_nest_negative_n_unevaluated);
    TEST(test_nest_non_integer_n_unevaluated);
    TEST(test_nest_wrong_arg_count);
    TEST(test_nest_double_pure);
    TEST(test_nest_factorial_like);
    TEST(test_nest_list_rotation);
    TEST(test_nest_head_changing);

    printf("All Nest tests passed!\n");
    return 0;
}
