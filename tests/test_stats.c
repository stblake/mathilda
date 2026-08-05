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
        "{a, a + (-a + b) x}", 0);
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

void test_moment() {
    /* --- scalar order on an exact vector: exact input -> exact output --- */
    assert_eval_eq("Moment[{1,2,3,4},2]", "15/2", 0);       /* (1+4+9+16)/4 */
    assert_eval_eq("Moment[{1,2,3,4},1]", "5/2", 0);        /* r=1 is the Mean */
    assert_eval_eq("Moment[{1,2,3,4},0]", "1", 0);          /* r=0 -> 1 */
    assert_eval_eq("Moment[{1,2,3,4},3]", "25", 0);         /* (1+8+27+64)/4 */
    assert_eval_eq("Moment[Range[10],2]", "77/2", 0);       /* 385/10 */

    /* --- approximate input -> approximate output --- */
    assert_eval_eq("Moment[{1.,2.,3.,4.},2]", "7.5", 0);
    /* arbitrary precision is preserved (routes through the MPFR symbolic path):
     * Moment[.,2] of {1,2,3} is 14/3 = 4.6666..., which shows its digits (unlike
     * an exact-in-binary value like 7.5, which prints short even at 30 digits). */
    assert_eval_startswith("Moment[N[{1,2,3},30],2]", "4.6666666666666666666666666");

    /* --- symbolic data --- */
    assert_eval_eq("Simplify[Moment[{a,b,c},2]]", "1/3 (a^2 + b^2 + c^2)", 0);
    /* r=1 is the mean, symbolically */
    assert_eval_eq("Moment[{Pi,E,2},1]", "1/3 (2 + E + Pi)", 0);

    /* --- matrix: columnwise moments --- */
    assert_eval_eq("Moment[{{1,2},{3,4},{5,6}},3]", "{51, 96}", 0);
    assert_eval_eq("Simplify[Moment[{{a11,a12},{a21,a22}},2]]",
                   "{1/2 (a11^2 + a21^2), 1/2 (a12^2 + a22^2)}", 0);

    /* --- rank-3 array: columnwise at the first level --- */
    assert_eval_eq("Moment[{{{1,2},{3,4}},{{5,6},{7,8}}},2]", "{{13, 20}, {29, 40}}", 0);

    /* --- multivariate (vector) order: sum the r_j-th power in the j-th column --- */
    assert_eval_eq("Moment[{{1,2},{3,5}},{1,1}]", "17/2", 0);    /* (1*2 + 3*5)/2 */
    assert_eval_eq("Simplify[Moment[{{a,b},{c,d}},{1,2}]]", "1/2 (a b^2 + c d^2)", 0);

    /* --- Association works on the values --- */
    assert_eval_eq("Moment[<|1->{1,2},2->{3,4},3->{5,6}|>,2]", "{35/3, 56/3}", 0);

    /* --- NDArray / packed surfaces --- */
    assert_eval_eq("Moment[NDArray[{1.,2.,3.,4.},DataType->\"float64\"],2]", "7.5", 0);
    /* a visible int64 buffer degrades to the exact Rational (like Variance) */
    assert_eval_eq("Moment[NDArray[{1,2,3,4},DataType->\"int64\"],2]", "15/2", 0);
    /* packed and unpacked agree on the same data */
    assert_eval_eq("Moment[NDArray[{1.,2.,4.,8.,16.},DataType->\"float64\"],2] "
                   "== Moment[{1.,2.,4.,8.,16.},2]", "True", 0);

    /* --- Compile[]: lowerable, and the compiled kernel matches the interpreter --- */
    assert_eval_eq("(CompileDiagnostics[{{v,_Real,1},{k,_Integer}}, Moment[v,k]] "
                   "/. List -> Association)[\"Compiled\"]", "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1}}, Moment[v,2]][{1.,2.,4.,8.,16.}] "
                   "== Moment[{1.,2.,4.,8.,16.},2]", "True", 0);

    /* --- attributes: NHoldAll and Protected --- */
    assert_eval_eq("MemberQ[Attributes[Moment], NHoldAll]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[Moment], Protected]", "True", 0);

    /* --- stays unevaluated when it should --- */
    assert_eval_eq("Moment[x,2]", "Moment[x, 2]", 0);
    assert_eval_eq("Moment[{1,2,3}]", "Moment[{1, 2, 3}]", 0);
}

/* Regression: the exact int64 rational fast paths in mean.c and variance.c
 * multiplied without checking for overflow, so a long list whose mean has a
 * large irreducible denominator returned a SILENTLY WRONG answer.
 *
 * Found by benchmarks/19-statistics, where Mathematica and NumPy agreed with
 * each other and Mathilda did not: Variance came back 0.0176 where the answer is
 * 84815.8, and some overflowed sums of squares were NEGATIVE, which a variance
 * cannot be.
 *
 * The trap is that friendly data hides it completely: a list whose mean reduces
 * to a small denominator (999/2, say) never overflows at any length. These cases
 * use Range[n]/d with d coprime to the sum so the denominator survives gcd
 * reduction and the accumulator actually grows.
 *
 * Each case pins the fast path against the definition computed through the
 * evaluator (exact via GMP) rather than against a written-down constant that
 * could be wrong in the same direction as the code. */
void test_exact_overflow_regression() {
    /* Mean of rationals with a surviving denominator, at a length that overflowed. */
    assert_eval_eq("Mean[Range[10000]/10007] == Total[Range[10000]/10007]/10000",
                   "True", 0);
    /* Variance squares an n-scaled numerator, so it overflows sooner than Mean. */
    assert_eval_eq("Variance[Range[4000]/4001] == "
                   "Total[(Range[4000]/4001 - Mean[Range[4000]/4001])^2]/3999",
                   "True", 0);
    /* A sum of squares over a positive count is never negative. */
    assert_eval_eq("Variance[Range[10000]/10007] > 0", "True", 0);
    /* The centred-and-powered list is what CentralMoment feeds back into Mean,
     * which is how Skewness/Kurtosis inherited the bug. */
    assert_eval_eq("CentralMoment[Range[4000]/4001, 2] == "
                   "Total[(Range[4000]/4001 - Mean[Range[4000]/4001])^2]/4000",
                   "True", 0);
    /* Exactness is preserved, not traded for safety: the fallback is GMP, so the
     * result is still an exact Rational rather than a Real. */
    assert_eval_eq("Head[Mean[Range[10000]/10007]]", "Rational", 0);
    assert_eval_eq("Head[Variance[Range[4000]/4001]]", "Rational", 0);
    /* Small inputs still take the int64 fast path, unchanged. */
    assert_eval_eq("Mean[{1/2, 1/3, 1/6}]", "1/3", 0);
    assert_eval_eq("Variance[{1, 2, 3, 4}]", "5/3", 0);
}

void test_central_moment() {
    /* --- scalar order on an exact vector: exact input -> exact output --- */
    assert_eval_eq("CentralMoment[{1,2,3,4},2]", "5/4", 0);
    assert_eval_eq("CentralMoment[{1,2,3,4},4]", "41/16", 0);
    assert_eval_eq("CentralMoment[Range[10],2]", "33/4", 0);
    assert_eval_eq("CentralMoment[{2,4,4,4,5,5,7,9},2]", "4", 0);
    assert_eval_eq("CentralMoment[{1,2,3,4},0]", "1", 0);   /* r=0 -> 1 */
    assert_eval_eq("CentralMoment[{1,2,3,4},1]", "0", 0);   /* first central moment is 0 */
    assert_eval_eq("CentralMoment[{1,2,3,4},3]", "0", 0);   /* odd moment of symmetric data */

    /* --- approximate input -> approximate output --- */
    assert_eval_eq("CentralMoment[{1.,2.,3.,4.},2]", "1.25", 0);
    /* arbitrary precision is preserved (routes through the MPFR symbolic path) */
    assert_eval_startswith("CentralMoment[N[{1,2,3},30],2]", "0.66666666666666666666666666666");

    /* --- symbolic data --- */
    assert_eval_eq("Simplify[CentralMoment[{a,b},2]]", "1/4 (a - b)^2", 0);
    assert_eval_eq("Together[CentralMoment[{Pi,E,2},2]]",
                   "1/9 (8 - 4 E + 2 E^2 - 4 Pi - 2 E Pi + 2 Pi^2)", 0);

    /* --- matrix: columnwise moments --- */
    assert_eval_eq("CentralMoment[{{1,2},{3,4},{5,6}},2]", "{8/3, 8/3}", 0);
    assert_eval_eq("Simplify[CentralMoment[{{a11,a12},{a21,a22}},4]]",
                   "{1/16 (a11 - a21)^4, 1/16 (a12 - a22)^4}", 0);

    /* --- rank-3 array: columnwise at the first level --- */
    assert_eval_eq("CentralMoment[{{{1,2},{3,4}},{{5,6},{7,8}}},2]", "{{4, 4}, {4, 4}}", 0);

    /* --- multivariate (vector) order --- */
    assert_eval_eq("CentralMoment[{{1,2},{3,5}},{1,1}]", "3/2", 0);
    assert_eval_eq("Simplify[CentralMoment[{{a,b},{c,d}},{2,2}]]", "1/16 (a - c)^2 (b - d)^2", 0);

    /* --- Association works on the values --- */
    assert_eval_eq("CentralMoment[<|1->{1,2},2->{3,4},3->{5,6}|>,2]", "{8/3, 8/3}", 0);

    /* --- NDArray / packed surfaces --- */
    assert_eval_eq("CentralMoment[NDArray[{1.,2.,3.,4.},DataType->\"float64\"],2]", "1.25", 0);
    /* a visible int64 buffer degrades to the exact Rational (like Variance) */
    assert_eval_eq("CentralMoment[NDArray[{1,2,3,4},DataType->\"int64\"],2]", "5/4", 0);
    /* packed and unpacked agree on the same data */
    assert_eval_eq("CentralMoment[NDArray[{1.,2.,4.,8.,16.},DataType->\"float64\"],2] "
                   "== CentralMoment[{1.,2.,4.,8.,16.},2]", "True", 0);

    /* --- Compile[]: lowerable, and the compiled kernel matches the interpreter --- */
    assert_eval_eq("(CompileDiagnostics[{{v,_Real,1},{k,_Integer}}, CentralMoment[v,k]] "
                   "/. List -> Association)[\"Compiled\"]", "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1}}, CentralMoment[v,2]][{1.,2.,4.,8.,16.}] "
                   "== CentralMoment[{1.,2.,4.,8.,16.},2]", "True", 0);

    /* --- stays unevaluated when it should --- */
    assert_eval_eq("CentralMoment[x,2]", "CentralMoment[x, 2]", 0);
    assert_eval_eq("CentralMoment[{1,2,3}]", "CentralMoment[{1, 2, 3}]", 0);
}

void test_skewness() {
    /* Skewness = CentralMoment[.,3] / CentralMoment[.,2]^(3/2) */
    assert_eval_eq("Skewness[{1,2,3,4,5}]", "0", 0);              /* symmetric */
    assert_eval_eq("Skewness[{1,2,3,10}]", "18/25 Sqrt[2]", 0);  /* exact radical */
    assert_eval_eq("Skewness[{1.,2.,3.,10.}]", "1.01823", 0);    /* approximate */
    assert_eval_eq("Simplify[Skewness[{a,b}]]", "0", 0);         /* two-point is symmetric */
    /* matrix: columnwise */
    assert_eval_eq("Skewness[{{1,2},{2,4},{3,6},{4,8},{5,10}}]", "{0, 0}", 0);
    /* equivalent to the CentralMoment ratio */
    assert_eval_eq("Skewness[{1,2,3,10}] == CentralMoment[{1,2,3,10},3]/CentralMoment[{1,2,3,10},2]^(3/2)",
                   "True", 0);
    /* Association works on the values */
    assert_eval_eq("Skewness[<|\"a\"->{1,2},\"b\"->{2,4},\"c\"->{3,6},\"d\"->{4,8},\"e\"->{5,10}|>]",
                   "{0, 0}", 0);
    /* NDArray: real fast path, and an int64 buffer degrades to the exact radical */
    assert_eval_eq("Skewness[NDArray[{1.,2.,3.,10.},DataType->\"float64\"]]", "1.01823", 0);
    assert_eval_eq("Skewness[NDArray[{1,2,3,10},DataType->\"int64\"]]", "18/25 Sqrt[2]", 0);
    /* Compile: lowerable, and the compiled kernel matches the interpreter */
    assert_eval_eq("(CompileDiagnostics[{{v,_Real,1}}, Skewness[v]] /. List -> Association)[\"Compiled\"]",
                   "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1}}, Skewness[v]][{1.,2.,4.,8.,16.,3.,7.}] "
                   "== Skewness[{1.,2.,4.,8.,16.,3.,7.}]", "True", 0);
    /* unevaluated when it should be */
    assert_eval_eq("Skewness[x]", "Skewness[x]", 0);
    assert_eval_eq("Skewness[{}]", "Skewness[{}]", 0);
}

void test_kurtosis() {
    /* Kurtosis = CentralMoment[.,4] / CentralMoment[.,2]^2 (Pearson, not excess) */
    assert_eval_eq("Kurtosis[{1,2,3,4,5}]", "17/10", 0);
    assert_eval_eq("Kurtosis[{1,2,4,8}]", "25141/13225", 0);     /* exact Rational */
    assert_eval_eq("Kurtosis[{1.,2.,3.,4.,5.}]", "1.7", 0);
    assert_eval_eq("Simplify[Kurtosis[{a,b}]]", "1", 0);         /* two-point is flat */
    assert_eval_eq("Kurtosis[{{1,2},{2,4},{3,6},{4,8},{5,10}}]", "{17/10, 17/10}", 0);
    assert_eval_eq("Kurtosis[{1,2,4,8}] == CentralMoment[{1,2,4,8},4]/CentralMoment[{1,2,4,8},2]^2",
                   "True", 0);
    assert_eval_eq("Kurtosis[NDArray[{1.,2.,3.,4.,5.},DataType->\"float64\"]]", "1.7", 0);
    assert_eval_eq("Kurtosis[NDArray[{1,2,3,4,5},DataType->\"int64\"]]", "17/10", 0);
    assert_eval_eq("(CompileDiagnostics[{{v,_Real,1}}, Kurtosis[v]] /. List -> Association)[\"Compiled\"]",
                   "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1}}, Kurtosis[v]][{1.,2.,4.,8.,16.,3.,7.}] "
                   "== Kurtosis[{1.,2.,4.,8.,16.,3.,7.}]", "True", 0);
    assert_eval_eq("Kurtosis[x]", "Kurtosis[x]", 0);
}

void test_covariance() {
    /* Exact input yields exact output. */
    assert_eval_eq("Covariance[{1,3/2},{2,11}]", "9/4", 0);
    /* Approximate input yields approximate output. */
    assert_eval_eq("Covariance[{1.5,3,5,10},{2,1.25,15,8}]", "11.2604", 0);
    /* Symbolic (real transcendentals) reduces to the closed form. */
    assert_eval_eq("Simplify[Covariance[{1,Pi},{E,2}]]", "-1/2 (2 - E) (1 - Pi)", 0);
    /* Complex vectors: the conjugate lands on the SECOND argument. */
    assert_eval_eq("Covariance[{2+I,3-2I,5+4I},{I,1+2I,10-5I}]", "-7/3 + 56/3*I", 0);
    /* Auto-covariance matrix, exact. */
    assert_eval_eq("Covariance[{{1,2},{3,4},{5,7}}]", "{{4, 5}, {5, 19/3}}", 0);
    /* Cross-covariance matrix is p x q. */
    assert_eval_eq("Covariance[{{a,b},{c,d}},{{x},{y}}] // Dimensions", "{2, 1}", 0);
    /* NDArray / packed fast paths agree with the List path. */
    assert_eval_eq("Covariance[NDArray[{1.,2.,3.,4.}],NDArray[{2.,1.,4.,3.}]] "
                   "== Covariance[{1.,2.,3.,4.},{2.,1.,4.,3.}]", "True", 0);
    assert_eval_eq("Chop[Normal[Covariance[NDArray[{{1.,5.},{3.,9.},{5.,6.}}]]] "
                   "- Covariance[{{1.,5.},{3.,9.},{5.,6.}}]]", "{{0, 0}, {0, 0}}", 0);
    /* An integer NDArray covariance is exact (a Rational), so it degrades. */
    assert_eval_eq("Covariance[NDArray[{1,2,3,4},DataType->\"int64\"],"
                   "NDArray[{2,1,4,3},DataType->\"int64\"]]", "1", 0);
    /* Packed >= 250-element vectors take the buffer path. */
    assert_eval_eq("Covariance[Range[1.,300.],Range[300.,1.,-1.]] "
                   "== Covariance[Table[N[i],{i,300}],Table[N[301-i],{i,300}]]", "True", 0);
    /* Compile[] lowering, standalone and inside a larger body (the subset cliff). */
    assert_eval_eq("(CompileDiagnostics[{{v,_Real,1},{w,_Real,1}}, Covariance[v,w]] "
                   "/. List -> Association)[\"Compiled\"]", "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1},{w,_Real,1}}, Covariance[v,w]]"
                   "[{1.5,3.,5.,10.},{2.,1.25,15.,8.}] "
                   "== Covariance[{1.5,3.,5.,10.},{2.,1.25,15.,8.}]", "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1},{w,_Real,1}}, Covariance[v,w]+1.0]"
                   "[{1.5,3.,5.,10.},{2.,1.25,15.,8.}]", "12.2604", 0);
    /* Errors / unevaluated forms. */
    assert_eval_eq("Covariance[]", "Covariance[]", 0);                    /* argb -> unevaluated */
    assert_eval_eq("Covariance[{1,2,3}]", "Covariance[{1, 2, 3}]", 0);    /* single vector */
    assert_eval_eq("Covariance[x,y]", "Covariance[x, y]", 0);             /* non-list args */
}

void test_correlation() {
    /* Exact input yields exact output. */
    assert_eval_eq("Correlation[{5,3/4,1},{2,1/2,1}]", "2 Sqrt[3/13]", 0);
    /* Approximate input yields approximate output. */
    assert_eval_eq("Correlation[{1.5,3,5,10},{2,1.25,15,8}]", "0.475976", 0);
    /* Symbolic (real transcendentals) closed form. */
    assert_eval_eq("Simplify[Correlation[{1,Pi,2},{2,2,1}]]",
                   "(1/2 (-3 + Pi))/Sqrt[3 - 3 Pi + Pi^2]", 0);
    /* Complex vectors. */
    assert_eval_eq("Correlation[{2+I,3-2I,5+4I},{I,1+2I,10-5I}]", "(-7/2 + 28*I)/Sqrt[1139]", 0);
    /* Auto-correlation diagonal is exactly 1 for exact/symbolic data ... */
    assert_eval_eq("Correlation[{{a,b},{c,d}}][[1,1]]", "1", 0);
    assert_eval_eq("Correlation[{{a,b},{c,d}}][[2,2]]", "1", 0);
    /* ... and the real 1. for real data, keeping the matrix uniform and symmetric. */
    assert_eval_eq("Correlation[{{1.,2.},{3.,4.},{5.,7.}}][[1,1]]", "1.0", 0);
    assert_eval_eq("SymmetricMatrixQ[Correlation[{{1.,5.,2.},{3.,9.,1.},"
                   "{5.,6.,8.},{2.,3.,4.},{7.,1.,0.}}]]", "True", 0);
    /* A correlation matrix is the covariance scaled by the standard deviations. */
    assert_eval_eq("Module[{d={{1.,5.,2.},{3.,9.,1.},{5.,6.,8.},{2.,3.,4.},{7.,1.,0.}},s}, "
                   "s=DiagonalMatrix[1/StandardDeviation[d]]; "
                   "Chop[Correlation[d] - s.Covariance[d].s]]",
                   "{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}", 0);
    /* NDArray auto-correlation matches the List path. */
    assert_eval_eq("Chop[Normal[Correlation[NDArray[{{1.,5.,2.},{3.,9.,1.},"
                   "{5.,6.,8.},{2.,3.,4.},{7.,1.,0.}}]]] - Correlation[{{1.,5.,2.},{3.,9.,1.},"
                   "{5.,6.,8.},{2.,3.,4.},{7.,1.,0.}}]]", "{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}", 0);
    /* NDArray / Compile vector correlation. */
    assert_eval_eq("Correlation[NDArray[{1.,2.,3.,4.}],NDArray[{2.,1.,4.,3.}]] "
                   "== Correlation[{1.,2.,3.,4.},{2.,1.,4.,3.}]", "True", 0);
    assert_eval_eq("(CompileDiagnostics[{{v,_Real,1},{w,_Real,1}}, Correlation[v,w]] "
                   "/. List -> Association)[\"Compiled\"]", "True", 0);
    assert_eval_eq("Compile[{{v,_Real,1},{w,_Real,1}}, Correlation[v,w]]"
                   "[{1.5,3.,5.,10.},{2.,1.25,15.,8.}] "
                   "== Correlation[{1.5,3.,5.,10.},{2.,1.25,15.,8.}]", "True", 0);
    /* Errors / unevaluated forms. */
    assert_eval_eq("Correlation[]", "Correlation[]", 0);
    assert_eval_eq("Correlation[{1,2,3}]", "Correlation[{1, 2, 3}]", 0);
}

int main() {
    symtab_init();
    core_init();
    TEST(test_quartiles);
    TEST(test_median);

    TEST(test_mean);
    TEST(test_rootmeansquare);
    TEST(test_variance);
    TEST(test_moment);
    TEST(test_exact_overflow_regression);
    TEST(test_central_moment);
    TEST(test_skewness);
    TEST(test_kurtosis);
    TEST(test_standard_deviation);
    TEST(test_covariance);
    TEST(test_correlation);
    TEST(test_moving_average);
    TEST(test_moving_median);
    TEST(test_exponential_moving_average);

    printf("All stats tests passed!\n");
    return 0;
}
