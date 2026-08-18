/*
 * test_solve_integers.c
 *
 * Unit tests for the integer-domain (Diophantine) solver
 * (src/solveint.c): Solve[eqns && constraints, vars, Integers].
 *
 * Outputs are compared against FullForm strings so the canonical
 * List[List[Rule[var, val]], ...] form (ascending by value tuple) is
 * asserted exactly.  Every case here is one the bounded engine solves
 * completely; the deferred families (Pell, Thue, ...) are asserted to
 * stay unevaluated rather than return a wrong answer.
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
    if (!e) { printf("FAIL: failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        free(got); expr_free(res); expr_free(e); ASSERT(0); return;
    }
    printf("PASS: %s\n", input);
    free(got); expr_free(res); expr_free(e);
}

/* The motivating case: positivity bounds x,y; the leaf x is solved by a
 * perfect-square test. */
static void test_motivating(void) {
    run_test("Solve[x^2 + 2 y^3 == 3681 && x > 0 && y > 0, {x, y}, Integers]",
        "List[List[Rule[x, 15], Rule[y, 12]], "
             "List[Rule[x, 41], Rule[y, 10]], "
             "List[Rule[x, 57], Rule[y, 6]]]");
}

/* Two quadratics; positivity of the first bounds all three variables and the
 * system has no integer solution. */
static void test_two_quadratics_empty(void) {
    run_test("Solve[x^2 + y^2 + z^2 == 1000 && x^2 + 2 y^2 - z^2 == 200 "
             "&& x > 0 && y > 0 && z > 0, {x, y, z}, Integers]",
        "List[]");
}

/* Explicit box + perfect-square / -power leaf. */
static void test_boxed_power(void) {
    run_test("Solve[x^2 + y^3 == z^5 && 0 < x < 1000 && 0 < y < 1000 "
             "&& 0 < z < 100, {x, y, z}, Integers]",
        "List[List[Rule[x, 104], Rule[y, 28], Rule[z, 8]], "
             "List[Rule[x, 654], Rule[y, 127], Rule[z, 19]]]");
}

/* Markov triples up to 1000: the quadratic-in-the-last-variable leaf. */
static void test_markov(void) {
    run_test("Length[Solve[x^2 + y^2 + z^2 == 3 x y z && 0 < x <= y <= z <= 1000, "
             "{x, y, z}, Integers]]", "13");
    run_test("Solve[x^2 + y^2 + z^2 == 3 x y z && 0 < x <= y <= z <= 5, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 1], Rule[z, 1]], "
             "List[Rule[x, 1], Rule[y, 1], Rule[z, 2]], "
             "List[Rule[x, 1], Rule[y, 2], Rule[z, 5]]]");
}

/* Taxicab: the smallest is 1729 = 1^3+12^3 = 9^3+10^3. */
static void test_taxicab_1729(void) {
    run_test("Solve[x^3 + y^3 == z^3 + w^3 && 0 < x < y && 0 < z < w && x != z "
             "&& x^3 + y^3 < 1730, {x, y, z, w}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 12], Rule[z, 9], Rule[w, 10]], "
             "List[Rule[x, 9], Rule[y, 10], Rule[z, 1], Rule[w, 12]]]");
}

/* Sum of four fifth powers equal to a fifth power, with d < 100 (which
 * excludes the (27,84,110,133;144) solution), solved via meet-in-the-middle. */
static void test_euler_sum_of_powers_empty(void) {
    run_test("Solve[a^5 + b^5 + c^5 + d^5 == e^5 && 0 < a <= b <= c <= d < 100 "
             "&& e < 150, {a, b, c, d, e}, Integers]", "List[]");
}

/* Integer-restriction: a non-integer rational / irrational solution is
 * dropped; an actual integer root is kept. */
static void test_integer_restriction(void) {
    run_test("Solve[2 x == 3, x, Integers]", "List[]");
    run_test("Solve[x^2 == 2, x, Integers]", "List[]");
    run_test("Solve[x^2 == 9, x, Integers]",
        "List[List[Rule[x, -3]], List[Rule[x, 3]]]");
    run_test("Solve[x^3 == 8, x, Integers]", "List[List[Rule[x, 2]]]");
}

/* Phase 2: Pythagorean-perimeter via linear elimination + bilinear divisor
 * factoring.  With z > 0 this is the classic "triangles of perimeter p". */
static void test_pythagorean_perimeter(void) {
    run_test("Solve[x^2 + y^2 == z^2 && x + y + z == 12 && 0 < x < y && z > 0, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 4], Rule[z, 5]]]");
    run_test("Solve[x^2 + y^2 == z^2 && x + y + z == 3000 && 0 < x < y && z > 0, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, 500], Rule[y, 1200], Rule[z, 1300]], "
             "List[Rule[x, 600], Rule[y, 1125], Rule[z, 1275]], "
             "List[Rule[x, 750], Rule[y, 1000], Rule[z, 1250]]]");
}

/* Phase 2: Egyptian fractions via the unit-fraction recursion. */
static void test_egyptian_fractions(void) {
    run_test("Solve[1 == 1/x + 1/y + 1/z && 0 < x <= y <= z, {x, y, z}, Integers]",
        "List[List[Rule[x, 2], Rule[y, 3], Rule[z, 6]], "
             "List[Rule[x, 2], Rule[y, 4], Rule[z, 4]], "
             "List[Rule[x, 3], Rule[y, 3], Rule[z, 3]]]");
    run_test("Length[Solve[4/2027 == 1/x + 1/y + 1/z && 0 < x <= y <= z, "
             "{x, y, z}, Integers]]", "73");
}

/* Phase 2d: separable odd-power sums via the divisor method (s = x+y | m).
 * These boxes are too large for the leaf search, so the divisor path runs. */
static void test_powersum_divisor(void) {
    /* Sum of two cubes: 1729 = 1^3+12^3 = 9^3+10^3 (Ramanujan). */
    run_test("Solve[x^3 + y^3 == 1729 && 0 < x <= y && x < 10^5, {x, y}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 12]], List[Rule[x, 9], Rule[y, 10]]]");
    /* Sum of three cubes == 3: the known small solutions 1,1,1 and 4,4,-5. */
    run_test("Solve[x^3 + y^3 + z^3 == 3 && Abs[x] < 8000 && Abs[y] < 8000 "
             "&& Abs[z] < 8000, {x, y, z}, Integers]",
        "List[List[Rule[x, -5], Rule[y, 4], Rule[z, 4]], "
             "List[Rule[x, 1], Rule[y, 1], Rule[z, 1]], "
             "List[Rule[x, 4], Rule[y, -5], Rule[z, 4]], "
             "List[Rule[x, 4], Rule[y, 4], Rule[z, -5]]]");
    /* No small representation of 42 as a sum of three cubes. */
    run_test("Solve[x^3 + y^3 + z^3 == 42 && Abs[x] < 8000 && Abs[y] < 8000 "
             "&& Abs[z] < 8000, {x, y, z}, Integers]", "List[]");
    /* Generality: two fifth powers, 1267 = 3^5 + 4^5. */
    run_test("Solve[x^5 + y^5 == 1267 && 0 < x <= y && x < 10^5, {x, y}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 4]]]");
}

/* Phase 3: Pell equations x^2 - D y^2 == +/-1 via continued fractions. */
static void test_pell(void) {
    /* The classic large fundamental solution of x^2 - 61 y^2 == 1. */
    run_test("Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0 && x < 10^10, {x, y}, Integers]",
        "List[List[Rule[x, 1766319049], Rule[y, 226153980]]]");
    /* The bounded orbit for D = 2. */
    run_test("Solve[x^2 - 2 y^2 == 1 && x > 0 && y > 0 && x < 100, {x, y}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 2]], "
             "List[Rule[x, 17], Rule[y, 12]], "
             "List[Rule[x, 99], Rule[y, 70]]]");
    /* Negative Pell: solvable for D = 2, unsolvable for D = 3. */
    run_test("Solve[x^2 - 2 y^2 == -1 && x > 0 && y > 0 && x < 100, {x, y}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 1]], "
             "List[Rule[x, 7], Rule[y, 5]], "
             "List[Rule[x, 41], Rule[y, 29]]]");
    run_test("Solve[x^2 - 3 y^2 == -1 && x > 0 && y > 0 && x < 100, {x, y}, Integers]",
        "List[]");
}

/* Phase 2b: an unconstrained linear equation gives the parametric family via
 * the gcd staircase; an unsolvable one gives {}. */
static void test_linear_parametric(void) {
    run_test("Solve[x + y == 10, {x, y}, Integers]",
        "List[List[Rule[x, C[1]], Rule[y, Plus[10, Times[-1, C[1]]]]]]");
    run_test("Solve[x + 2 y == 5, {x, y}, Integers]",
        "List[List[Rule[x, Plus[5, Times[2, C[1]]]], Rule[y, Times[-1, C[1]]]]]");
    /* gcd(2,2)=2 does not divide 3 -> no solution. */
    run_test("Solve[2 x + 2 y == 3, {x, y}, Integers]", "List[]");
    /* Bounded but unsolvable: gcd(314159265, 271828182, 161803398) = 3 does not
     * divide 1, so no solution regardless of the box. */
    run_test("Solve[314159265 x + 271828182 y + 161803398 z == 1 "
             "&& Abs[x] < 10^6 && Abs[y] < 10^6 && Abs[z] < 10^6, {x, y, z}, Integers]",
        "List[]");
    /* Bounded, solvable, large box -> LLL lattice enumeration.  The 2-variable
     * lattice is a single arithmetic progression; this box holds 2000 points. */
    run_test("Length[Solve[1000003 x + 999983 y == 7 && Abs[x] < 10^9 && Abs[y] < 10^9, "
             "{x, y}, Integers]]", "2000");
    /* Three variables -> two parameters. */
    run_test("Solve[x + y + z == 5, {x, y, z}, Integers]",
        "List[List[Rule[x, C[1]], "
             "Rule[y, Plus[C[2], Times[-1, C[1]]]], "
             "Rule[z, Plus[5, Times[-1, C[2]]]]]]");
}

/* Deferred families must stay unevaluated (never a wrong answer). */
static void test_deferred_unevaluated(void) {
    /* Mordell y^2 = x^3 + k: elliptic-integral-points phase, no bound. */
    run_test("Solve[y^2 == x^3 - 10000 && x > 0 && y > 0, {x, y}, Integers]",
        "Solve[And[Equal[Power[y, 2], Plus[-10000, Power[x, 3]]], "
             "Greater[x, 0], Greater[y, 0]], List[x, y], Integers]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_motivating);
    TEST(test_two_quadratics_empty);
    TEST(test_boxed_power);
    TEST(test_markov);
    TEST(test_taxicab_1729);
    TEST(test_euler_sum_of_powers_empty);
    TEST(test_pythagorean_perimeter);
    TEST(test_egyptian_fractions);
    TEST(test_powersum_divisor);
    TEST(test_pell);
    TEST(test_linear_parametric);
    TEST(test_integer_restriction);
    TEST(test_deferred_unevaluated);

    printf("All solve-integers tests passed!\n");
    return 0;
}
