#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

void test_mean() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"Mean[{1, 2, 3, 4}]", "5/2"},
        {"Mean[{1.2, 2.8}]", "2.0"},
        {"Mean[{a, b, c, d}]", "1/4 (a + b + c + d)"},
        {"Mean[{{a, u}, {b, v}, {c, w}}]", "{1/3 (a + b + c), 1/3 (u + v + w)}"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Mean test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_variance() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"Variance[{1, 2, 3}]", "1"},
        {"Variance[{1, 2, 3, 4}]", "5/3"},
        {"Variance[{{5.2, 7}, {5.3, 8}, {5.4, 9}}]", "{0.01, 1}"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("Variance test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}

void test_standard_deviation() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"StandardDeviation[{1, 2, 3}]", "1"},
        {"StandardDeviation[{{5.2, 7}, {5.3, 8}, {5.4, 9}}]", "{0.1, 1}"}
    };

    for (int i = 0; i < (int)(sizeof(tests) / sizeof(tests[0])); i++) {
        Expr* e = parse_expression(tests[i].input);
        Expr* res = evaluate(e);
        char* res_str = expr_to_string(res);
        if (strcmp(res_str, tests[i].expected) != 0) {
            printf("StandardDeviation test failed: %s expected %s, got %s\n", tests[i].input, tests[i].expected, res_str);
            ASSERT(0);
        }
        free(res_str);
        expr_free(e);
        expr_free(res);
    }
}


void test_median() {
    assert_eval_eq("Median[{1,2,3,4,5,6,7}]", "4", 0);
    assert_eval_eq("Median[{1,2,3,4,5,6,7,8}]", "9/2", 0);
    assert_eval_eq("Median[{1,2,3,4}]", "5/2", 0);
    assert_eval_eq("Median[{Pi,E,2}]", "E", 0);
    assert_eval_eq("Median[{1.,2.,3.,4.}]", "2.5", 0);
    assert_eval_eq("Median[{{1,11,3},{4,6,7}}]", "{5/2, 17/2, 5}", 0);
    assert_eval_eq("Median[{{{3,7},{2,1}},{{5,19},{12,4}}}]", "{{4, 13}, {7, 5/2}}", 0);
    assert_eval_eq("Median[{a,b,c}]", "Median[{a, b, c}]", 0);
}


void test_rootmeansquare() {
    assert_eval_eq("RootMeanSquare[{a,b,c,d}]", "1/2 Sqrt[a^2 + b^2 + c^2 + d^2]", 0);
    assert_eval_eq("RootMeanSquare[{{1,2},{5,10},{5,2},{4,8}}]", "{1/2 Sqrt[67], Sqrt[43]}", 0);
    assert_eval_eq("RootMeanSquare[{1,2,3,4}]", "Sqrt[15/2]", 0);
    assert_eval_eq("RootMeanSquare[{Pi,E,2}]", "Sqrt[1/3 (4 + E^2 + Pi^2)]", 0);
    assert_eval_eq("RootMeanSquare[{1.,2.,3.,4.}]", "2.73861", 0);
}

void test_quartiles() {
    assert_eval_eq("Quartiles[{1,3,4,2,5,6}]", "{2, 7/2, 5}", 0);
    assert_eval_eq("Quartiles[{1,2,3,4}]", "{3/2, 5/2, 7/2}", 0);
    assert_eval_eq("Quartiles[{1.,2.,3.,4.}]", "{1.5, 2.5, 3.5}", 0);
    assert_eval_eq("Quartiles[{-1,5,10,4,25,2,1}]", "{5/4, 4, 35/4}", 0);
    assert_eval_eq("Quartiles[{-1,5,10,4,25,2,1},{{0,0},{1,0}}]", "{1, 4, 10}", 0);
    assert_eval_eq("Quartiles[{{1,11,3},{4,6,7}}]", "{{1, 5/2, 4}, {6, 17/2, 11}, {3, 5, 7}}", 0);
    assert_eval_eq("Quartiles[{{{3,7},{2,1}},{{5,19},{12,4}}}]", "{{{3, 4, 5}, {7, 13, 19}}, {{2, 7, 12}, {1, 5/2, 4}}}", 0);
    assert_eval_eq("Quartiles[{a,b,c}]", "Quartiles[{a, b, c}]", 0);
}

void test_moving_average() {
    /* Pairwise simple moving average: symbolic data, factored output. */
    assert_eval_eq("MovingAverage[{a, b, c, d, e}, 2]",
                   "{1/2 (a + b), 1/2 (b + c), 1/2 (c + d), 1/2 (d + e)}", 0);

    /* Weighted moving average: spec-style distributed output. */
    assert_eval_eq("MovingAverage[{a, b, c, d, e}, {1, 2}]",
                   "{1/3 a + 2/3 b, 1/3 b + 2/3 c, 1/3 c + 2/3 d, 1/3 d + 2/3 e}", 0);

    /* Integer data yields exact rationals. */
    assert_eval_eq("MovingAverage[{1, 5, 7, 3, 6, 2}, 3]",
                   "{13/3, 5, 16/3, 11/3}", 0);

    /* Approximate (real) data yields approximate output. */
    assert_eval_eq("MovingAverage[{1.2, 5.2, 3.4, 4.5, 2.3, 4.5}, 3]",
                   "{3.26667, 4.36667, 3.4, 3.76667}", 0);

    /* Window length 1 is the identity. */
    assert_eval_eq("MovingAverage[{1, 2, 3, 4, 5}, 1]", "{1, 2, 3, 4, 5}", 0);

    /* Window length equal to list length gives a single-element list. */
    assert_eval_eq("MovingAverage[{1, 2, 3, 4, 5}, 5]", "{3}", 0);
    assert_eval_eq("MovingAverage[{a, b, c}, 3]", "{1/3 (a + b + c)}", 0);

    /* Output length is Length[list] - r + 1. */
    assert_eval_eq("Length[MovingAverage[Range[10], 4]]", "7", 0);
    assert_eval_eq("MovingAverage[Range[10], 4]",
                   "{5/2, 7/2, 9/2, 11/2, 13/2, 15/2, 17/2}", 0);

    /* Equal weights match the unweighted average value (distributed form). */
    assert_eval_eq("MovingAverage[{1, 2, 3, 4}, {1, 1, 1}]", "{2, 3}", 0);
    assert_eval_eq("MovingAverage[{a, b, c, d}, {1, 1, 1}]",
                   "{1/3 a + 1/3 b + 1/3 c, 1/3 b + 1/3 c + 1/3 d}", 0);

    /* Mixed rationals stay exact. */
    assert_eval_eq("MovingAverage[{1/2, 1/3, 1/6, 5/6}, 2]", "{5/12, 1/4, 1/2}", 0);

    /* Real-valued weights produce approximate output. */
    assert_eval_eq("MovingAverage[{1, 2, 3, 4}, {0.5, 0.5}]", "{1.5, 2.5, 3.5}", 0);

    /* Bignum support: large windows over arbitrary-precision integers. */
    assert_eval_eq("MovingAverage[{2^100, 2^101, 2^102, 2^103}, 2]",
                   "{1901475900342344102245054808064, "
                   "3802951800684688204490109616128, "
                   "7605903601369376408980219232256}", 0);
    assert_eval_eq("MovingAverage[{2^200, 2^200, 2^200, 2^200}, 3]",
                   "{1606938044258990275541962092341162602522202993782792835301376, "
                   "1606938044258990275541962092341162602522202993782792835301376}", 0);

    /* Edge cases: stay unevaluated when r is out of range or shape is wrong. */
    assert_eval_eq("MovingAverage[{1, 2, 3, 4, 5}, 6]",
                   "MovingAverage[{1, 2, 3, 4, 5}, 6]", 0);
    assert_eval_eq("MovingAverage[{1, 2, 3}, 0]",
                   "MovingAverage[{1, 2, 3}, 0]", 0);
    assert_eval_eq("MovingAverage[{1, 2, 3}, -1]",
                   "MovingAverage[{1, 2, 3}, -1]", 0);
    assert_eval_eq("MovingAverage[{1, 2, 3}, 2.5]",
                   "MovingAverage[{1, 2, 3}, 2.5]", 0);
    assert_eval_eq("MovingAverage[{}, 1]",
                   "MovingAverage[{}, 1]", 0);
    assert_eval_eq("MovingAverage[x, 2]",
                   "MovingAverage[x, 2]", 0);
    assert_eval_eq("MovingAverage[{1, 2, 3}, {}]",
                   "MovingAverage[{1, 2, 3}, {}]", 0);

    /* MovingAverage is Protected. */
    assert_eval_eq("MemberQ[Attributes[MovingAverage], Protected]", "True", 0);
}

void test_moving_median() {
    /* Vector cases — odd window. */
    assert_eval_eq("MovingMedian[{1,2,5,6,1,4,3},3]", "{2, 5, 5, 4, 3}", 0);

    /* Even window yields exact rationals. */
    assert_eval_eq("MovingMedian[{1,2,3,4},2]", "{3/2, 5/2, 7/2}", 0);

    /* Matrix input: column-wise medians within each row-window. */
    assert_eval_eq("MovingMedian[{{1,2},{5,3},{1,4},{3,2},{5,5}},2]",
                   "{{3, 5/2}, {3, 7/2}, {2, 3}, {4, 7/2}}", 0);

    /* Machine-precision (real) data preserves approximate output. */
    assert_eval_eq("MovingMedian[N[{1,5,7,3,6,2}],3]",
                   "{5.0, 5.0, 6.0, 3.0}", 0);
    assert_eval_eq("MovingMedian[{1.,5.,7.,3.,6.,2.},3]",
                   "{5.0, 5.0, 6.0, 3.0}", 0);

    /* Pi and E are accepted (NumericQ-real). r==Length yields a one-element list. */
    assert_eval_eq("MovingMedian[{Pi, E, 2}, 3]", "{E}", 0);

    /* Window length 1 is the identity (each element is its own median). */
    assert_eval_eq("MovingMedian[{1, 2, 3, 4, 5}, 1]", "{1, 2, 3, 4, 5}", 0);

    /* Window length equal to list length gives a single median. */
    assert_eval_eq("MovingMedian[{1, 2, 3, 4, 5}, 5]", "{3}", 0);

    /* Mixed exact rationals stay exact. */
    assert_eval_eq("MovingMedian[{1/2, 1/3, 1/6, 5/6}, 3]",
                   "{1/3, 1/3}", 0);

    /* Output length is Length[list] - r + 1. */
    assert_eval_eq("Length[MovingMedian[Range[10], 4]]", "7", 0);
    assert_eval_eq("MovingMedian[Range[10], 4]",
                   "{5/2, 7/2, 9/2, 11/2, 13/2, 15/2, 17/2}", 0);

    /* Bignum support. */
    assert_eval_eq("MovingMedian[{2^100, 2^101, 2^102, 2^103}, 2]",
                   "{1901475900342344102245054808064, "
                   "3802951800684688204490109616128, "
                   "7605903601369376408980219232256}", 0);
    assert_eval_eq("MovingMedian[{2^200, 2^201, 2^202}, 3]",
                   "{3213876088517980551083924184682325205044405987565585670602752}", 0);

    /* Edge cases: stay unevaluated when r is out of range or input shape is wrong. */
    assert_eval_eq("MovingMedian[{1, 2, 3, 4, 5}, 6]",
                   "MovingMedian[{1, 2, 3, 4, 5}, 6]", 0);
    assert_eval_eq("MovingMedian[{1, 2, 3}, 0]",
                   "MovingMedian[{1, 2, 3}, 0]", 0);
    assert_eval_eq("MovingMedian[{1, 2, 3}, -1]",
                   "MovingMedian[{1, 2, 3}, -1]", 0);
    assert_eval_eq("MovingMedian[{1, 2, 3}, 2.5]",
                   "MovingMedian[{1, 2, 3}, 2.5]", 0);
    assert_eval_eq("MovingMedian[{}, 1]",
                   "MovingMedian[{}, 1]", 0);
    assert_eval_eq("MovingMedian[x, 2]",
                   "MovingMedian[x, 2]", 0);

    /* Symbolic data triggers MovingMedian::arg1 and stays unevaluated. */
    assert_eval_eq("MovingMedian[{a, b, c}, 2]",
                   "MovingMedian[{a, b, c}, 2]", 0);

    /* MovingMedian is Protected. */
    assert_eval_eq("MemberQ[Attributes[MovingMedian], Protected]", "True", 0);
}

void test_exponential_moving_average() {
    /* Mathematica reference: exact-rational output for integer data + rational alpha. */
    assert_eval_eq("ExponentialMovingAverage[Range[10], 1/3]",
        "{1, 4/3, 17/9, 70/27, 275/81, 1036/243, 3773/729, 13378/2187, 46439/6561, 158488/19683}", 0);

    /* Machine-precision (real) data with rational alpha — fast double path. */
    assert_eval_eq("ExponentialMovingAverage[N[{1,5,7,3,6,2}], 1/2]",
        "{1.0, 3.0, 5.0, 4.0, 5.0, 3.5}", 0);
    assert_eval_eq("ExponentialMovingAverage[{1.,5.,7.,3.,6.,2.}, 1/2]",
        "{1.0, 3.0, 5.0, 4.0, 5.0, 3.5}", 0);
    /* Real alpha. */
    assert_eval_eq("ExponentialMovingAverage[{1, 2.5, 3, 4.5}, 0.5]",
        "{1.0, 1.75, 2.375, 3.4375}", 0);
    /* Mixed Integer/Real data with rational alpha promotes to Real. */
    assert_eval_eq("ExponentialMovingAverage[{1, 2.5, 3, 4.5}, 1/2]",
        "{1.0, 1.75, 2.375, 3.4375}", 0);

    /* alpha = 0: every output equals the first element (constant series). */
    assert_eval_eq("ExponentialMovingAverage[{a, b, c, d}, 0]", "{a, a, a, a}", 0);
    assert_eval_eq("ExponentialMovingAverage[Range[5], 0]", "{1, 1, 1, 1, 1}", 0);

    /* alpha = 1: output equals the input (no smoothing). */
    assert_eval_eq("ExponentialMovingAverage[{a, b, c, d}, 1]", "{a, b, c, d}", 0);
    assert_eval_eq("ExponentialMovingAverage[Range[5], 1]", "{1, 2, 3, 4, 5}", 0);

    /* Single-element list returns itself unchanged for any alpha. */
    assert_eval_eq("ExponentialMovingAverage[{5}, 1/2]", "{5}", 0);
    assert_eval_eq("ExponentialMovingAverage[{a}, x]", "{a}", 0);
    assert_eval_eq("ExponentialMovingAverage[{3.5}, 0.7]", "{3.5}", 0);

    /* Two-element with rational alpha. y2 = 1 + (1/4)(2-1) = 5/4. */
    assert_eval_eq("ExponentialMovingAverage[{1, 2}, 1/4]", "{1, 5/4}", 0);

    /* Six-element exact rational walk-through (Mathematica equivalent). */
    assert_eval_eq("ExponentialMovingAverage[{1, 5, 7, 3, 6, 2}, 1/2]",
        "{1, 3, 5, 4, 5, 7/2}", 0);

    /* Even integer input where every step stays integral. */
    assert_eval_eq("ExponentialMovingAverage[{2, 4, 6, 8, 10}, 1/2]",
        "{2, 3, 9/2, 25/4, 65/8}", 0);

    /* Symbolic alpha — recurrence kept (no Distribute). */
    assert_eval_eq("ExponentialMovingAverage[{a, b}, x]",
        "{a, a + x (-a + b)}", 0);
    assert_eval_eq("ExponentialMovingAverage[{a, b}, 1/2]",
        "{a, a + 1/2 (-a + b)}", 0);

    /* Output length is preserved. */
    assert_eval_eq("Length[ExponentialMovingAverage[Range[20], 1/3]]", "20", 0);
    assert_eval_eq("Length[ExponentialMovingAverage[Range[7], 0.3]]", "7", 0);
    assert_eval_eq("Length[ExponentialMovingAverage[{a, b, c, d, e, f}, x]]", "6", 0);

    /* Bignum support: alpha = 1 gives the input verbatim, including large GMP integers. */
    assert_eval_eq("ExponentialMovingAverage[{2^100, 2^101}, 1]",
        "{1267650600228229401496703205376, 2535301200456458802993406410752}", 0);

    /* Bignum + rational alpha: y2 = 2^99 + 2^199 stays exact. */
    assert_eval_eq("ExponentialMovingAverage[{2^100, 2^200}, 1/2]",
        "{1267650600228229401496703205376, "
        "803469022129495137770981046171215126561215611592144769253376}", 0);

    /* Edge cases: stay unevaluated when shape or arity is wrong. */
    assert_eval_eq("ExponentialMovingAverage[{}, 1/2]",
        "ExponentialMovingAverage[{}, 1/2]", 0);
    assert_eval_eq("ExponentialMovingAverage[x, 1/2]",
        "ExponentialMovingAverage[x, 1/2]", 0);
    assert_eval_eq("ExponentialMovingAverage[{1, 2, 3}]",
        "ExponentialMovingAverage[{1, 2, 3}]", 0);
    assert_eval_eq("ExponentialMovingAverage[{1, 2, 3}, 1/2, 7]",
        "ExponentialMovingAverage[{1, 2, 3}, 1/2, 7]", 0);

    /* ExponentialMovingAverage is Protected. */
    assert_eval_eq("MemberQ[Attributes[ExponentialMovingAverage], Protected]", "True", 0);
}

int main() {
    symtab_init();
    core_init();
    TEST(test_quartiles);
    TEST(test_median);

    TEST(test_mean);
    TEST(test_rootmeansquare);
    TEST(test_variance);
    TEST(test_standard_deviation);
    TEST(test_moving_average);
    TEST(test_moving_median);
    TEST(test_exponential_moving_average);

    printf("All stats tests passed!\n");
    return 0;
}
