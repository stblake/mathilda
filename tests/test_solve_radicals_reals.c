/* test_solve_radicals_reals.c
 *
 * Pinned-FullForm regression for Solve dropping EXTRANEOUS and NON-REAL
 * Root[] solutions.  Two fixes are exercised:
 *
 *   (A) src/solverad.c -- the radical-equation specialist now back-substitutes
 *       each Root[] candidate (N[] residual) instead of accepting it verbatim,
 *       so the branches of x^(1/lcm) that do not satisfy the original Sqrt /
 *       x^(p/q) are rejected in BOTH the Reals and the default Complexes domain.
 *
 *   (B) src/solve.c -- Solve's shared post-dispatch funnel drops any solution
 *       whose bound value is a provably non-real number when the domain is
 *       Reals / Integers / Rationals, covering irreducible polynomials and
 *       polynomial systems that emit complex Root[] objects/tuples.
 *
 * Regression origin: Solve[Sqrt[x] + 3 x^(1/3) == 5, x, Reals] returned three
 * Root[] objects (one real + two complex extraneous); Mathematica returns only
 * the single valid one.  Companion form-invariant coverage lives in the
 * "O-*" rows of tests/solve_corpus.m.
 */

#include <stdio.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("FAIL: failed to parse: %s\n", input);
        ASSERT(0);
        return;
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        free(res_str);
        expr_free(res);
        expr_free(e);
        ASSERT(0);
        return;
    }
    printf("PASS: %s -> %s\n", input, res_str);
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* The reported bug: the cleared resultant is the irreducible cubic
 * x^3 + 6 x^2 + 8625 x - 15625; only its real root (index 1) satisfies the
 * original equation.  Over Reals AND over the default domain, the two complex
 * roots must be gone. */
#define CUBIC_ROOT(k) \
    "Root[Function[Plus[-15625, Times[8625, Slot[1]], " \
    "Times[6, Power[Slot[1], 2]], Power[Slot[1], 3]]], " #k "]"

static void test_reported_cbrt_sqrt_reals(void) {
    run_test("Solve[Sqrt[x] + 3 x^(1/3) == 5, x, Reals]",
             "List[List[Rule[x, " CUBIC_ROOT(1) "]]]");
}

static void test_reported_cbrt_sqrt_default_domain(void) {
    /* Default (Complexes): the extraneous roots fail by NOT SATISFYING the
     * equation (residual O(1)), so Fix A alone must still trim to one. */
    run_test("Solve[Sqrt[x] + 3 x^(1/3) == 5, x]",
             "List[List[Rule[x, " CUBIC_ROOT(1) "]]]");
}

static void test_more_radical_cubics_reals(void) {
    run_test("Solve[Sqrt[x] + x^(1/3) == 3, x, Reals]",
             "List[List[Rule[x, Root[Function[Plus[-729, Times[297, Slot[1]], "
             "Times[-10, Power[Slot[1], 2]], Power[Slot[1], 3]]], 1]]]]");
    run_test("Solve[Sqrt[x] - x^(1/3) == 1, x, Reals]",
             "List[List[Rule[x, Root[Function[Plus[-1, Slot[1], "
             "Times[-10, Power[Slot[1], 2]], Power[Slot[1], 3]]], 1]]]]");
}

/* Irreducible polynomials over Reals: complex Root objects dropped (Fix B). */
static void test_poly_reals_single_real_root(void) {
    run_test("Solve[x^3 + 6 x^2 + 8625 x - 15625 == 0, x, Reals]",
             "List[List[Rule[x, " CUBIC_ROOT(1) "]]]");
    run_test("Solve[x^5 - x - 1 == 0, x, Reals]",
             "List[List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 5]]], 1]]]]");
}

static void test_poly_reals_two_real_roots(void) {
    /* x^6 - x - 1 has exactly two real roots (indices 1 and 2). */
    run_test("Solve[x^6 - x - 1 == 0, x, Reals]",
             "List["
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 6]]], 1]]], "
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 6]]], 2]]]]");
}

/* No regression over the default (Complexes) domain: all roots retained. */
static void test_no_regression_complexes(void) {
    run_test("Solve[x^5 - x - 1 == 0, x]",
             "List["
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 5]]], 1]]], "
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 5]]], 2]]], "
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 5]]], 3]]], "
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 5]]], 4]]], "
             "List[Rule[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 5]]], 5]]]]");
}

/* Reality filter is conservative: real radical values survive, symbolic /
 * parametric values (Sqrt[a], which is not provably non-real) survive, and a
 * genuinely-empty real set stays empty. */
static void test_reals_filter_conservative(void) {
    run_test("Solve[x^2 - 2 == 0, x, Reals]",
             "List[List[Rule[x, Times[-1, Power[2, Rational[1, 2]]]]], "
             "List[Rule[x, Power[2, Rational[1, 2]]]]]");
    run_test("Solve[x^2 == a, x, Reals]",
             "List[List[Rule[x, Times[-1, Power[a, Rational[1, 2]]]]], "
             "List[Rule[x, Power[a, Rational[1, 2]]]]]");
    run_test("Solve[x^2 + 1 == 0, x, Reals]", "List[]");
}

/* The explicit-number radical path (already correct before the fix) stays
 * correct: extraneous concrete root dropped. */
static void test_explicit_radical_extraneous(void) {
    run_test("Solve[Sqrt[x] == x - 2, x, Reals]", "List[List[Rule[x, 4]]]");
}

int main(void) {
    symtab_init();
    core_init();
    printf("Running solve radicals/reals tests...\n");
    TEST(test_reported_cbrt_sqrt_reals);
    TEST(test_reported_cbrt_sqrt_default_domain);
    TEST(test_more_radical_cubics_reals);
    TEST(test_poly_reals_single_real_root);
    TEST(test_poly_reals_two_real_roots);
    TEST(test_no_regression_complexes);
    TEST(test_reals_filter_conservative);
    TEST(test_explicit_radical_extraneous);
    printf("All solve radicals/reals tests passed!\n");
    return 0;
}
