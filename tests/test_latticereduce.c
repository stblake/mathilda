/* test_latticereduce.c -- unit tests for LatticeReduce (LLL reduction).
 *
 * Coverage:
 *   - Exact reduced-basis output reproducing the documented examples.
 *   - Lattice invariants that hold regardless of the exact reduced basis:
 *       * Abs[Det] preserved (square lattices),
 *       * right null-space relations preserved (rectangular relation
 *         lattices, integer and Gaussian),
 *   - arbitrary-precision (10^20) integer-relation finding,
 *   - rational and Gaussian-rational entries,
 *   - single-vector / trivial cases,
 *   - argument and entry validation diagnostics (argx / matrix / latm).
 */
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, res_str);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* Exact reduced bases (these reproduce the documented Wolfram output). */
static void test_exact_bases(void) {
    run_test("LatticeReduce[{{1, 0, 0, 1345}, {0, 1, 0, 35}, {0, 0, 1, 154}}]",
             "List[List[0, 9, -2, 7], List[1, 1, -9, -6], List[1, -3, -8, 8]]");
    run_test("LatticeReduce[{{12, 2}, {13, 4}}]",
             "List[List[1, 2], List[9, -4]]");
    run_test("LatticeReduce[{{1, 0, 0, -1}, {0, 1, 0, -2}, {0, 0, 1, -3}}]",
             "List[List[1, 0, 0, -1], List[-1, 1, 0, -1], List[-1, -1, 1, 0]]");
    /* Single generator: already reduced. */
    run_test("LatticeReduce[{{3, 4}}]", "List[List[3, 4]]");
}

/* Abs[Det] (the lattice covolume) is preserved by reduction. */
static void test_determinant_invariant(void) {
    run_test("Abs[Det[{{12, 2}, {13, 4}}]] == "
             "Abs[Det[LatticeReduce[{{12, 2}, {13, 4}}]]]", "True");
    run_test("Abs[Det[{{8, 4, 7}, {8, 3, 6}, {7, 6, 5}}]] == "
             "Abs[Det[LatticeReduce[{{8, 4, 7}, {8, 3, 6}, {7, 6, 5}}]]]", "True");
    run_test("Abs[Det[{{1, 2, 3}, {4, 9, 5}, {7, 2, 11}}]] == "
             "Abs[Det[LatticeReduce[{{1, 2, 3}, {4, 9, 5}, {7, 2, 11}}]]]", "True");
    /* Rational lattice covolume preserved. */
    run_test("Abs[Det[{{1/2, 0}, {0, 1/3}}]] == "
             "Abs[Det[LatticeReduce[{{1/2, 0}, {0, 1/3}}]]]", "True");
}

/* Linear relations in the right null space survive reduction. */
static void test_relations_preserved(void) {
    /* a . {1,2,3,1} == {0,0,0};  the reduced basis must give the same. */
    run_test("Dot[LatticeReduce[{{1, 0, 0, -1}, {0, 1, 0, -2}, {0, 0, 1, -3}}], "
             "{1, 2, 3, 1}]", "List[0, 0, 0]");
    /* Gaussian relation lattice over Z[i]. */
    run_test("Dot[LatticeReduce[{{1, 0, Complex[-2, -3]}, {0, 1, Complex[-1, 1]}}], "
             "{Complex[2, 3], Complex[1, -1], 1}]", "List[0, 0]");
}

/* Arbitrary-precision (bignum) integer-relation finding: the leading
 * reduced vector recovers Pi/4 + ArcTan[1/239] - 4 ArcTan[1/5] == 0. */
static void test_bignum_relation(void) {
    run_test("First[LatticeReduce[{"
             "{1, 0, 0, 0, -100000000000000000000}, "
             "{0, 1, 0, 0, -78539816339744830962}, "
             "{0, 0, 1, 0, -19739555984988075837}, "
             "{0, 0, 0, 1, -418407600207472386}}]]",
             "List[0, 1, -4, 1, 0]");
}

/* Gaussian / rational structured output. */
static void test_gaussian_rational(void) {
    run_test("LatticeReduce[{{1/2, 0}, {0, 1/3}}]",
             "List[List[0, Rational[1, 3]], List[Rational[1, 2], 0]]");
}

/* Validation diagnostics: the call is left unevaluated. */
static void test_errors(void) {
    run_test("LatticeReduce[]", "LatticeReduce[]");
    run_test("LatticeReduce[{}]", "LatticeReduce[List[]]");
    run_test("LatticeReduce[{{1, 2}, {3, 4.5}}]",
             "LatticeReduce[List[List[1, 2], List[3, 4.5]]]");
    /* Non-rectangular (jagged). */
    run_test("LatticeReduce[{{1, 2}, {3}}]",
             "LatticeReduce[List[List[1, 2], List[3]]]");
    /* Symbolic entry. */
    run_test("LatticeReduce[{{1, x}, {0, 1}}]",
             "LatticeReduce[List[List[1, x], List[0, 1]]]");
    /* Linearly dependent generating set (two vectors in 1-D). */
    run_test("LatticeReduce[{{6}, {10}}]",
             "LatticeReduce[List[List[6], List[10]]]");
}

int main(void) {
    symtab_init();
    core_init();

    test_exact_bases();
    test_determinant_invariant();
    test_relations_preserved();
    test_bignum_relation();
    test_gaussian_rational();
    test_errors();

    printf("\nAll LatticeReduce tests passed.\n");
    return 0;
}
