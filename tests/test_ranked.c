/* Unit tests for RankedMin / RankedMax — the n-th smallest / largest element.
 *
 * RankedMin[list, n]  = n-th smallest (n<0: n-th largest)
 * RankedMax[list, n]  = n-th largest  (n<0: n-th smallest)
 * Identity under test throughout: RankedMax[list, k] == RankedMin[list, -k].
 *
 * Results are compared in FullForm so an Integer answer (12) is distinguished
 * from a Real (12.) and the exact element (E, Pi, 1/4) is verified verbatim. */
#include "sort.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e); expr_free(res);
}

/* Every worked example from the Wolfram documentation for RankedMin. */
void test_rankedmin_documented() {
    run_test("RankedMin[{12, 13, 11}, 2]", "12");            /* 2nd smallest of three */
    run_test("RankedMin[{Pi, Sqrt[2], E, 3}, 3]", "3");      /* 3rd smallest of four */
    run_test("RankedMin[{12.6555, 3.5265, 25.65}, 2]", "12.6555");
    run_test("RankedMin[{12, 13, 48.5, Pi}, 4]", "48.5");    /* 4th smallest = largest */
    run_test("RankedMin[{2.5, E, 12, 15, 485}, -2]", "15");  /* 2nd largest */
    run_test("RankedMin[{2.5, E, 12, 15, 485}, -4]", "E");   /* 4th largest */
    run_test("RankedMin[{2.5, E, 12, 15, 485}, -5]", "2.5"); /* 5th largest = smallest */
    run_test("RankedMin[{12, 52, E, Pi}, -4]", "E");         /* smallest of four */
}

/* Every worked example from the Wolfram documentation for RankedMax. */
void test_rankedmax_documented() {
    run_test("RankedMax[{12, 13, 11}, 2]", "12");            /* 2nd largest of three */
    run_test("RankedMax[{Pi, Sqrt[2], E, 3}, 3]", "E");      /* 3rd largest of four */
    run_test("RankedMax[{12.6555, 33.5265, 25.65}, 2]", "25.65");
    run_test("RankedMax[{12, 13, 48.5, Pi}, 4]", "Pi");      /* 4th largest = smallest */
    run_test("RankedMax[{2.5, E, 12, 15, 485}, -2]", "E");   /* 2nd smallest */
    run_test("RankedMax[{2.5, E, 12, 15, 485}, -4]", "15");  /* 4th smallest */
    run_test("RankedMax[{2.5, E, 12, 15, 485}, -5]", "485"); /* 5th smallest = largest */
}

/* Values at +-Infinity rank as +-infinity, not symbolically. */
void test_ranked_infinity() {
    run_test("RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]", "Infinity");
    run_test("RankedMax[{Infinity, 5, Infinity, -Infinity}, -1]", "Times[-1, Infinity]");
    run_test("RankedMin[{Infinity, 5, Infinity, -Infinity}, 1]", "Times[-1, Infinity]");
    run_test("RankedMin[{Infinity, 5, -Infinity}, 2]", "5");
    run_test("RankedMin[{-Infinity, Infinity}, -1]", "Infinity");
}

/* RankedMin[l,1] == Min[l], RankedMin[l,-1] == Max[l], and the general
 * RankedMax[l,k] == RankedMin[l,-k] identity. */
void test_ranked_minmax_identities() {
    run_test("RankedMin[{3, 1, 2}, 1] === Min[{3, 1, 2}]", "True");
    run_test("RankedMin[{3, 1, 2}, -1] === Max[{3, 1, 2}]", "True");
    run_test("RankedMax[{3, 1, 2}, 1] === Max[{3, 1, 2}]", "True");
    run_test("RankedMax[{3, 1, 2}, -1] === Min[{3, 1, 2}]", "True");
    run_test("RankedMax[{5, 2, 8, 1, 9, 3}, 2] === RankedMin[{5, 2, 8, 1, 9, 3}, -2]", "True");
    run_test("RankedMax[{5, 2, 8, 1, 9, 3}, -3] === RankedMin[{5, 2, 8, 1, 9, 3}, 3]", "True");
}

/* Exact ordering — arbitrary-precision integers and rationals go through the
 * exact (expr_compare) path, so two bignums past 2^53 are ordered by value. */
void test_ranked_exact() {
    run_test("RankedMin[{1/2, 1/3, 1/4, 1/5}, 2]", "Rational[1, 4]");
    run_test("RankedMax[{1/2, 1/3, 1/4, 1/5}, 1]", "Rational[1, 2]");
    run_test("RankedMin[{1000000000000000000000000000000, "
             "1000000000000000000000000000001, "
             "999999999999999999999999999999}, 2]",
             "1000000000000000000000000000000");
    run_test("RankedMax[{1000000000000000000000000000000, "
             "1000000000000000000000000000001, "
             "999999999999999999999999999999}, 1]",
             "1000000000000000000000000000001");
    /* Full sweep of ranks on a fixed list, both heads. */
    run_test("RankedMin[{40, 10, 30, 20, 50}, 1]", "10");
    run_test("RankedMin[{40, 10, 30, 20, 50}, 3]", "30");
    run_test("RankedMin[{40, 10, 30, 20, 50}, 5]", "50");
    run_test("RankedMax[{40, 10, 30, 20, 50}, 1]", "50");
    run_test("RankedMax[{40, 10, 30, 20, 50}, 3]", "30");
    run_test("RankedMax[{40, 10, 30, 20, 50}, 5]", "10");
}

/* Duplicates: a repeated value is a valid answer at consecutive ranks. */
void test_ranked_duplicates() {
    run_test("RankedMin[{3, 1, 4, 1, 5}, 1]", "1");
    run_test("RankedMin[{3, 1, 4, 1, 5}, 2]", "1");
    run_test("RankedMin[{3, 1, 4, 1, 5}, 3]", "3");
    run_test("RankedMin[{2, 2, 2}, 2]", "2");
}

/* Symbolic real constants order by numeric value; a mixed exact/symbolic list
 * still resolves. */
void test_ranked_symbolic_reals() {
    run_test("RankedMin[{Pi, E, Sqrt[2]}, 1]", "Power[2, Rational[1, 2]]");  /* Sqrt[2] */
    run_test("RankedMax[{Pi, E, Sqrt[2]}, 1]", "Pi");
    run_test("RankedMin[{Pi + E, Pi, E}, 2]", "Pi");   /* E < Pi < Pi+E */
    run_test("RankedMin[{10, Pi, 2}, 2]", "Pi");       /* 2 < Pi < 10 */
}

/* Cases that must stay unevaluated: bad rank, empty, non-list, non-real. */
void test_ranked_unevaluated() {
    run_test("RankedMin[{1, 2, 3}, 5]", "RankedMin[List[1, 2, 3], 5]");   /* out of range */
    run_test("RankedMin[{1, 2, 3}, -5]", "RankedMin[List[1, 2, 3], -5]");
    run_test("RankedMin[{1, 2, 3}, 0]", "RankedMin[List[1, 2, 3], 0]");   /* zero rank */
    run_test("RankedMin[{1, 2, 3}, 1.5]", "RankedMin[List[1, 2, 3], 1.5]"); /* non-integer */
    run_test("RankedMin[{}, 1]", "RankedMin[List[], 1]");                 /* empty */
    run_test("RankedMin[{x, 1, 2}, 1]", "RankedMin[List[x, 1, 2], 1]");   /* free symbol */
    run_test("RankedMin[{I, 1, 2}, 1]", "RankedMin[List[Complex[0, 1], 1, 2], 1]"); /* non-real */
    run_test("RankedMax[{1, 2, 3}, 9]", "RankedMax[List[1, 2, 3], 9]");
    run_test("RankedMin[5, 1]", "RankedMin[5, 1]");                       /* non-list */
    /* Wrong arg count -> error / unevaluated (both args required). */
    run_test("RankedMin[{1, 2, 3}]", "RankedMin[List[1, 2, 3]]");
}

/* NDArray argument takes the buffer fast path and keeps the dtype:
 * int64 selects an Integer, float64 a Real. */
void test_ranked_ndarray() {
    run_test("Head[RankedMin[NDArray[{3, 1, 4, 1, 5}, DataType -> \"int64\"], 2]]", "Integer");
    run_test("RankedMin[NDArray[{3, 1, 4, 1, 5}, DataType -> \"int64\"], 2]", "1");
    run_test("RankedMax[NDArray[{3, 1, 4, 1, 5}, DataType -> \"int64\"], 2]", "4");
    run_test("RankedMin[NDArray[{3, 1, 4, 1, 5}, DataType -> \"int64\"], -1]", "5");
    run_test("RankedMin[NDArray[{3., 1., 4., 1., 5.}], 2]", "1.0");
    run_test("RankedMax[NDArray[{3., 1., 4., 1., 5.}], 1]", "5.0");
    /* A packed integer list (Range) is on the buffer path and stays exact. */
    run_test("Head[RankedMin[Range[100], 5]]", "Integer");
    run_test("RankedMin[Range[100], 5]", "5");
    run_test("RankedMax[Range[100], 5]", "96");
    /* Rank-2 array degrades to the List path, which has no real elements. */
    run_test("RankedMin[{{1, 2}, {3, 4}}, 1]", "RankedMin[List[List[1, 2], List[3, 4]], 1]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_rankedmin_documented);
    TEST(test_rankedmax_documented);
    TEST(test_ranked_infinity);
    TEST(test_ranked_minmax_identities);
    TEST(test_ranked_exact);
    TEST(test_ranked_duplicates);
    TEST(test_ranked_symbolic_reals);
    TEST(test_ranked_unevaluated);
    TEST(test_ranked_ndarray);

    printf("All ranked tests passed!\n");
    return 0;
}
