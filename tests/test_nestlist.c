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
    ASSERT_MSG(strcmp(s, expected) == 0, "NestList %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "NestList %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* Check that the result is a real-valued list of the expected length whose
 * elements approximate `expected` within tolerance. */
static void run_real_list_approx(const char* input, const double* expected, size_t n, double tol) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    ASSERT_MSG(r->type == EXPR_FUNCTION, "NestList %s: expected List, got non-function", input);
    ASSERT_MSG(r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "List") == 0,
               "NestList %s: expected List head", input);
    ASSERT_MSG(r->data.function.arg_count == n,
               "NestList %s: expected length %zu, got %zu", input, n, r->data.function.arg_count);
    for (size_t i = 0; i < n; i++) {
        Expr* ei = r->data.function.args[i];
        double got = 0.0;
        if (ei->type == EXPR_REAL) got = ei->data.real;
        else if (ei->type == EXPR_INTEGER) got = (double)ei->data.integer;
        else {
            char* s = expr_to_string(ei);
            ASSERT_MSG(0, "NestList %s: element %zu not numeric: %s", input, i, s);
            free(s);
        }
        ASSERT_MSG(fabs(got - expected[i]) <= tol,
                   "NestList %s: element %zu expected ~%g, got %g", input, i, expected[i], got);
    }
    expr_free(e);
    expr_free(r);
}

/* ---------- Examples from the spec ---------- */

static void test_nestlist_basic_symbolic(void) {
    /* NestList[f, x, 4] -> {x, f[x], f[f[x]], f[f[f[x]]], f[f[f[f[x]]]]} */
    run_full("NestList[f, x, 4]",
             "List[x, f[x], f[f[x]], f[f[f[x]]], f[f[f[f[x]]]]]");
}

static void test_nestlist_cos_real(void) {
    /* NestList[Cos, 1.0, 10] produces Cos iterated values */
    double expected[11] = {
        1.0,
        0.5403023058681398,
        0.8575532158463934,
        0.6542897904977791,
        0.7934803587425656,
        0.7013687736227565,
        0.7639596829006542,
        0.7221024250267077,
        0.7504177617637605,
        0.7314040424225099,
        0.7442373549005569
    };
    run_real_list_approx("NestList[Cos, 1.0, 10]", expected, 11, 1e-5);
}

static void test_nestlist_pure_function_symbolic(void) {
    /* NestList[(1+#)^2 &, x, 3] -> {x, (1+x)^2, (1+(1+x)^2)^2, (1+(1+(1+x)^2)^2)^2} */
    run_infix("NestList[(1 + #)^2 &, x, 3]",
              "{x, (1 + x)^2, (1 + (1 + x)^2)^2, (1 + (1 + (1 + x)^2)^2)^2}");
}

static void test_nestlist_sqrt_real(void) {
    /* NestList[Sqrt, 100.0, 4] -> {100, 10, 3.16228, 1.77828, 1.33352} */
    double expected[5] = {100.0, 10.0, 3.1622776601683795, 1.7782794100389228, 1.3335214321633239};
    run_real_list_approx("NestList[Sqrt, 100.0, 4]", expected, 5, 1e-5);
}

static void test_nestlist_powers_of_two(void) {
    /* NestList[2# &, 1, 10] -> {1, 2, 4, 8, ..., 1024} */
    run_full("NestList[2 # &, 1, 10]",
             "List[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]");
}

static void test_nestlist_successive_integers(void) {
    /* NestList[# + 1 &, 0, 10] -> {0, 1, ..., 10} */
    run_full("NestList[# + 1 &, 0, 10]",
             "List[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]");
}

static void test_nestlist_successive_squaring(void) {
    /* NestList[#^2 &, 2, 6] -> {2, 4, 16, 256, 65536, 4294967296, 18446744073709551616} */
    run_full("NestList[#^2 &, 2, 6]",
             "List[2, 4, 16, 256, 65536, 4294967296, 18446744073709551616]");
}

static void test_nestlist_compound_growth(void) {
    /* NestList[#(1 + 0.05) &, 1000, 10] -> {1000, 1050, ..., ~1628.89} */
    double expected[11] = {
        1000.0, 1050.0, 1102.5, 1157.625, 1215.50625,
        1276.2815625, 1340.0956406, 1407.1004226, 1477.4554437,
        1551.328216, 1628.894627
    };
    run_real_list_approx("NestList[#(1 + 0.05) &, 1000, 10]", expected, 11, 1e-3);
}

static void test_nestlist_newton_sqrt2(void) {
    /* NestList[(# + 2/#)/2 &, 1.0, 5] converges to sqrt(2) */
    double expected[6] = {
        1.0,
        1.5,
        1.4166666666666667,
        1.4142156862745099,
        1.4142135623746899,
        1.4142135623730951
    };
    run_real_list_approx("NestList[(# + 2/#)/2 &, 1.0, 5]", expected, 6, 1e-6);
}

static void test_nestlist_collatz(void) {
    /* Iterates of the 3n+1 problem starting at 100, 20 steps */
    run_full("NestList[If[EvenQ[#], #/2, (3 # + 1)/2] &, 100, 20]",
             "List[100, 50, 25, 38, 19, 29, 44, 22, 11, 17, 26, 13, 20, 10, 5, 8, 4, 2, 1, 2, 1]");
}

static void test_nestlist_lcg(void) {
    /* Linear congruential pseudorandom generator */
    run_full("NestList[Mod[59 #, 101] &, 1, 15]",
             "List[1, 59, 47, 46, 88, 41, 96, 8, 68, 73, 65, 98, 25, 61, 64, 39]");
}

/* ---------- Edge cases ---------- */

static void test_nestlist_zero(void) {
    /* NestList[f, x, 0] -> {x} */
    run_full("NestList[f, x, 0]", "List[x]");
}

static void test_nestlist_zero_pure(void) {
    /* NestList pure function, n=0 returns a singleton list */
    run_full("NestList[(# + 1) &, 7, 0]", "List[7]");
}

static void test_nestlist_one(void) {
    /* NestList[f, x, 1] -> {x, f[x]} */
    run_full("NestList[f, x, 1]", "List[x, f[x]]");
}

static void test_nestlist_length(void) {
    /* Length[NestList[f, x, 7]] should be 8 */
    run_full("Length[NestList[f, x, 7]]", "8");
}

static void test_nestlist_evaluates_each_step(void) {
    /* Each step should evaluate to a concrete value */
    run_full("NestList[# + 1 &, 0, 5]", "List[0, 1, 2, 3, 4, 5]");
}

static void test_nestlist_with_symbol_head(void) {
    /* NestList with a pure function wrapping Plus */
    run_full("NestList[Plus[#, 1] &, 5, 3]", "List[5, 6, 7, 8]");
}

static void test_nestlist_list_building(void) {
    /* NestList accumulating via Append */
    run_full("NestList[Append[#, a] &, {}, 3]",
             "List[List[], List[a], List[a, a], List[a, a, a]]");
}

static void test_nestlist_negative_n_unevaluated(void) {
    /* Negative n should leave NestList unevaluated */
    run_full("NestList[f, x, -1]", "NestList[f, x, -1]");
}

static void test_nestlist_non_integer_n_unevaluated(void) {
    /* Non-integer n should leave NestList unevaluated */
    run_full("NestList[f, x, n]", "NestList[f, x, n]");
}

static void test_nestlist_wrong_arg_count(void) {
    /* 2 args - unevaluated */
    run_full("NestList[f, x]", "NestList[f, x]");
}

static void test_nestlist_head_changing(void) {
    /* Apply a function that changes the head each step */
    run_full("NestList[f, g[x], 2]", "List[g[x], f[g[x]], f[f[g[x]]]]");
}

static void test_nestlist_first_last(void) {
    /* First should be expr, Last should match Nest */
    run_full("First[NestList[# + 1 &, 10, 5]]", "10");
    run_full("Last[NestList[# + 1 &, 10, 5]]", "15");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_nestlist_basic_symbolic);
    TEST(test_nestlist_cos_real);
    TEST(test_nestlist_pure_function_symbolic);
    TEST(test_nestlist_sqrt_real);
    TEST(test_nestlist_powers_of_two);
    TEST(test_nestlist_successive_integers);
    TEST(test_nestlist_successive_squaring);
    TEST(test_nestlist_compound_growth);
    TEST(test_nestlist_newton_sqrt2);
    TEST(test_nestlist_collatz);
    TEST(test_nestlist_lcg);

    TEST(test_nestlist_zero);
    TEST(test_nestlist_zero_pure);
    TEST(test_nestlist_one);
    TEST(test_nestlist_length);
    TEST(test_nestlist_evaluates_each_step);
    TEST(test_nestlist_with_symbol_head);
    TEST(test_nestlist_list_building);
    TEST(test_nestlist_negative_n_unevaluated);
    TEST(test_nestlist_non_integer_n_unevaluated);
    TEST(test_nestlist_wrong_arg_count);
    TEST(test_nestlist_head_changing);
    TEST(test_nestlist_first_last);

    printf("All NestList tests passed!\n");
    return 0;
}
