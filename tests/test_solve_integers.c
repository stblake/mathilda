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

/* Phase 4: exponential Diophantine (variable exponents), handled before the
 * MPoly stage.  The Catalan shape x^a - y^b == +/-1 uses Mihailescu's theorem
 * (the unique solution is 3^2 - 2^3 = 1); other shapes over a finite box are
 * enumerated exactly. */
static void test_exponential(void) {
    /* Catalan / Mihailescu: the only solution with bases, exponents >= 2. */
    run_test("Solve[x^a - y^b == 1 && 1 < x < 100 && 1 < y < 100 && a > 1 && b > 1, "
             "{x, y, a, b}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 2], Rule[a, 2], Rule[b, 3]]]");
    run_test("Solve[x^a - y^b == -1 && 1 < x < 100 && 1 < y < 100 && a > 1 && b > 1, "
             "{x, y, a, b}, Integers]",
        "List[List[Rule[x, 2], Rule[y, 3], Rule[a, 3], Rule[b, 2]]]");
    /* The unique Catalan solution excluded by the box -> no solution. */
    run_test("Solve[x^a - y^b == 1 && 1 < x < 3 && 1 < y < 100 && a > 1 && b > 1, "
             "{x, y, a, b}, Integers]", "List[]");
    /* Bounded exponent search: 2^2 - 3^3 = -23. */
    run_test("Solve[2^a - 3^b == -23 && 0 < a < 10 && 0 < b < 10, {a, b}, Integers]",
        "List[List[Rule[a, 2], Rule[b, 3]]]");
    run_test("Solve[2^a - 3^b == 100 && 0 < a < 10 && 0 < b < 10, {a, b}, Integers]",
        "List[]");
}

/* Phase 4: elliptic / hyperelliptic curves over a finite box (y^m == f(x)) are
 * solved by the ordinary bounded search -- enumerate x, test that f(x) is a
 * perfect power.  Included here as regression coverage for the families. */
static void test_superelliptic_bounded(void) {
    /* Mordell y^2 = x^3 - 2: Fermat's unique integral point (3, +/-5). */
    run_test("Solve[y^2 == x^3 - 2 && 0 < x < 1000 && y > 0, {x, y}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 5]]]");
    /* Mordell y^2 = x^3 + 1: the five integral points. */
    run_test("Solve[y^2 == x^3 + 1 && -100 < x < 100, {x, y}, Integers]",
        "List[List[Rule[x, -1], Rule[y, 0]], "
             "List[Rule[x, 0], Rule[y, -1]], List[Rule[x, 0], Rule[y, 1]], "
             "List[Rule[x, 2], Rule[y, -3]], List[Rule[x, 2], Rule[y, 3]]]");
    /* y^2 = x^3 - 10000 has the integral point (25, 75). */
    run_test("Solve[y^2 == x^3 - 10000 && 0 < x < 100000 && y > 0, {x, y}, Integers]",
        "List[List[Rule[x, 25], Rule[y, 75]]]");
    /* Hyperelliptic y^2 = x^4 - 4x^3 + 5x^2 - 2x = x(x-1)^2(x-2): no y > 0. */
    run_test("Solve[x^4 - 4 x^3 + 5 x^2 - 2 x - y^2 == 0 && 0 < x < 100000 && y > 0, "
             "{x, y}, Integers]", "List[]");
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

/* Multi-leaf staged elimination: the Euler brick has three "determined"
 * variables (a,b,c), each an exact square root per (x,y,z) -- only x<y<z is
 * enumerated.  The smallest brick is (44,117,240; 125,244,267). */
static void test_euler_brick(void) {
    run_test("Solve[x^2+y^2==a^2 && x^2+z^2==b^2 && y^2+z^2==c^2 "
             "&& 0<x<y<z<250 && a>0 && b>0 && c>0, {x,y,z,a,b,c}, Integers]",
        "List[List[Rule[x, 44], Rule[y, 117], Rule[z, 240], "
             "Rule[a, 125], Rule[b, 244], Rule[c, 267]]]");
}

/* Ordering-reduced estimate + int64 fast leaf: a single quadratic-in-w leaf
 * over the ordered triple x<=y<=z (raw box is otherwise declined). */
static void test_quadratic_form(void) {
    run_test("Solve[2 (x^2+y^2+z^2+w^2) == (x+y+z+w)^2 && 0<x<=y<=z<=w<20, "
             "{x,y,z,w}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 1], Rule[z, 4], Rule[w, 12]], "
             "List[Rule[x, 2], Rule[y, 2], Rule[z, 3], Rule[w, 15]]]");
    /* The four-variable Markov-Hurwitz equation  Sum x_i^2 == Prod x_i. */
    run_test("Solve[x1^2+x2^2+x3^2+x4^2 == x1 x2 x3 x4 && 0<x1<=x2<=x3<=x4<=100, "
             "{x1,x2,x3,x4}, Integers]",
        "List[List[Rule[x1, 2], Rule[x2, 2], Rule[x3, 2], Rule[x4, 2]], "
             "List[Rule[x1, 2], Rule[x2, 2], Rule[x3, 2], Rule[x4, 6]], "
             "List[Rule[x1, 2], Rule[x2, 2], Rule[x3, 6], Rule[x4, 22]], "
             "List[Rule[x1, 2], Rule[x2, 2], Rule[x3, 22], Rule[x4, 82]]]");
}

/* Non-polynomial bounded power-leaf: Brocard's problem  n! + 1 == m^2.  The
 * only known solutions (Brown numbers) are n = 4, 5, 7. */
static void test_brocard(void) {
    run_test("Solve[Factorial[n] + 1 == m^2 && n > 0 && m > 0 && n < 100, "
             "{n, m}, Integers]",
        "List[List[Rule[n, 4], Rule[m, 5]], "
             "List[Rule[n, 5], Rule[m, 11]], "
             "List[Rule[n, 7], Rule[m, 71]]]");
}

/* Binary-quadratic conic  Y^2 == A X^2 + B X + C  with A a perfect square:
 * complete the square to a difference of squares and factor the constant. */
static void test_conic(void) {
    /* Euler's prime-generating polynomial  n^2 + n + 41 == y^2  (D = 163). */
    run_test("Solve[n^2 + n + 41 == y^2 && n > 0 && y > 0, {n, y}, Integers]",
        "List[List[Rule[n, 40], Rule[y, 41]]]");
    /* Difference of two squares. */
    run_test("Solve[x^2 - y^2 == 15 && x > 0 && y > 0, {x, y}, Integers]",
        "List[List[Rule[x, 4], Rule[y, 1]], List[Rule[x, 8], Rule[y, 7]]]");
}

/* Unbounded positive Pell -> parametric fundamental-unit family.  The family is
 * a ConditionalExpression in C[1] >= 1; substituting C[1] = k recovers the k-th
 * solution (checked here for k = 1, 2 on x^2 - 2 y^2 == 1). */
static void test_pell_parametric(void) {
    run_test("Simplify[Solve[x^2 - 2 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers] "
             "/. C[1] -> 1]", "List[List[Rule[x, 3], Rule[y, 2]]]");
    run_test("Simplify[Solve[x^2 - 2 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers] "
             "/. C[1] -> 2]", "List[List[Rule[x, 17], Rule[y, 12]]]");
}

/* Homogeneous linear SYSTEM with positivity -> parametric ray (the primitive
 * kernel vector times C[1] >= 1), or {} when the kernel is mixed-sign (no
 * positive solution).  The kernel is the generalised cross product. */
static void test_linear_system_ray(void) {
    /* Kernel (6,4,3): substituting C[1]=2 gives (12,8,6). */
    run_test("Solve[2 x == 3 y && 4 z == 3 y && x>0 && y>0 && z>0, {x,y,z}, Integers] "
             "/. C[1] -> 2", "List[List[Rule[x, 12], Rule[y, 8], Rule[z, 6]]]");
    /* This rational system forces x = -13 z / 35, so there is NO positive
     * solution -- correctly resolved to {} (was previously unevaluated). */
    run_test("Solve[w == 5/6 x + y && x == 9/20 y + z && y == 13/42 z + w "
             "&& w > 0 && x > 0 && y > 0 && z > 0, {w, x, y, z}, Integers]", "List[]");
}

/* Fixed-base exponential  P^m - Q^n == +/-1  (constant bases, variable
 * exponents), sound & complete via Mihailescu + the exponent-1 cases. */
static void test_fixed_base_exponential(void) {
    /* 3^m - 2^n == 1: only (1,1) and (2,3), even though m,n are unbounded. */
    run_test("Solve[3^m - 2^n == 1 && m > 0 && n > 0, {m, n}, Integers]",
        "List[List[Rule[m, 1], Rule[n, 1]], List[Rule[m, 2], Rule[n, 3]]]");
    /* 2^n - 3^m == 1: only 2^2 - 3 = 1. */
    run_test("Solve[2^n - 3^m == 1 && m > 0 && n > 0, {m, n}, Integers]",
        "List[List[Rule[m, 1], Rule[n, 2]]]");
}

/* Prouhet-Tarry-Escott: two triples with equal power sums for degrees 1,2,3
 * are the same multiset (Newton's identities), so the strict orderings force
 * a=d, and the disequation a!=d makes the system empty. */
static void test_prouhet_tarry_escott(void) {
    run_test("Solve[a + b + c == d + e + f && a^2 + b^2 + c^2 == d^2 + e^2 + f^2 "
             "&& a^3 + b^3 + c^3 == d^3 + e^3 + f^3 && 0 < a < b < c && 0 < d < e < f "
             "&& a != d, {a, b, c, d, e, f}, Integers]", "List[]");
}

/* Deferred families must stay unevaluated (never a wrong answer). */
static void test_deferred_unevaluated(void) {
    /* Mordell y^2 = x^3 + k: elliptic-integral-points phase, no bound. */
    run_test("Solve[y^2 == x^3 - 10000 && x > 0 && y > 0, {x, y}, Integers]",
        "Solve[And[Equal[Power[y, 2], Plus[-10000, Power[x, 3]]], "
             "Greater[x, 0], Greater[y, 0]], List[x, y], Integers]");
    /* Correctness guard: an UNBOUNDED nonlinear curve that solveint cannot
     * finitely bound must stay unevaluated, NOT collapse to {} via the
     * parametric dispatch's non-integer closed form.  (y^2 == x^3 - 2 is now
     * solved by the Z[sqrt k] Mordell path, so the guard is exercised here by
     * y^2 == x^3 + 1, whose k = +1 is outside the sound PID cases, and by the
     * genuinely-parametric y == x^2.) */
    run_test("Solve[y^2 == x^3 + 1, {x, y}, Integers]",
        "Solve[Equal[Power[y, 2], Plus[1, Power[x, 3]]], List[x, y], Integers]");
    run_test("Solve[y == x^2, {x, y}, Integers]",
        "Solve[Equal[y, Power[x, 2]], List[x, y], Integers]");
}

/* Unbounded Mordell y^2 == x^3 + k, solved COMPLETELY via factorisation in
 * Z[sqrt k] for every imaginary k = 2,3 (mod 4) with |k| squarefree and 3 not
 * dividing the class number.  Verified against brute force for all engaged
 * k in [-150,-2].  (k = +3 etc. -- the real-quadratic case -- are declined.) */
static void test_mordell(void) {
    run_test("Solve[y^2 == x^3 - 2, {x, y}, Integers]",
        "List[List[Rule[x, 3], Rule[y, -5]], List[Rule[x, 3], Rule[y, 5]]]");
    run_test("Solve[y^2 == x^3 - 2 && y > 0, {x, y}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 5]]]");
    run_test("Solve[y^2 == x^3 - 1, {x, y}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 0]]]");
    /* Generalised: k = -13 (class number 2, 3 does not divide it). */
    run_test("Solve[y^2 == x^3 - 13, {x, y}, Integers]",
        "List[List[Rule[x, 17], Rule[y, -70]], List[Rule[x, 17], Rule[y, 70]]]");
    /* k = -5: the descent proves there is NO integer point (not merely unfound). */
    run_test("Solve[y^2 == x^3 - 5, {x, y}, Integers]", "List[]");
    /* k = -7 (half-integer ring) and k = -4 (not squarefree) stay unevaluated. */
    run_test("Solve[y^2 == x^3 - 7, {x, y}, Integers]",
        "Solve[Equal[Power[y, 2], Plus[-7, Power[x, 3]]], List[x, y], Integers]");
}

/* Correctness (P0) + HNF (P0b): a general linear SYSTEM over the Integers.
 *
 * The underdetermined case used to return a silent wrong `{}` -- the
 * Complexes-oriented linear-system dispatch expressed the pivot variables as a
 * RATIONAL parametric family (y -> (4 + 2 x)/5), then rejected each RHS as "not
 * a concrete integer", emptying a set with infinitely many integer points.  The
 * HNF path (`si_solve_linear_system_hnf`) now returns the COMPLETE integer
 * family: a particular solution + the kernel lattice as C[k].  The particular
 * representative is whatever the deterministic HNF pivots produce; both cases
 * below are verified against their equations at several C[1] values. */
static void test_linear_system_integers_hnf(void) {
    /* Underdetermined: 2 equations, 3 unknowns -> a one-parameter family. */
    run_test("Solve[{x + 2 y + 3 z == 10, x - y + z == 2}, {x, y, z}, Integers]",
        "List[List[Rule[x, Plus[18, Times[5, C[1]]]], "
             "Rule[y, Plus[8, Times[2, C[1]]]], "
             "Rule[z, Plus[-8, Times[-3, C[1]]]]]]");
    /* Confirm the family satisfies both equations for a specific parameter. */
    run_test("{x + 2 y + 3 z, x - y + z} /. "
             "(Solve[{x + 2 y + 3 z == 10, x - y + z == 2}, {x, y, z}, "
             "Integers] /. C[1] -> 3)[[1]]",
        "List[10, 2]");
    /* A second underdetermined system. */
    run_test("Solve[{x + y + z == 6, x + 2 y + 4 z == 10}, {x, y, z}, Integers]",
        "List[List[Rule[x, Plus[14, Times[2, C[1]]]], "
             "Rule[y, Plus[-14, Times[-3, C[1]]]], "
             "Rule[z, Plus[6, C[1]]]]]");
    /* Determined, integer solution -> read off. */
    run_test("Solve[{x + y == 5, x - y == 1}, {x, y}, Integers]",
        "List[List[Rule[x, 3], Rule[y, 2]]]");
    /* Determined, non-integer solution (3/2, 1/2) -> {} is a real proof. */
    run_test("Solve[{x + y == 2, x - y == 1}, {x, y}, Integers]", "List[]");
    /* Inconsistent -> {} is a real proof. */
    run_test("Solve[{x + y == 1, x + y == 2}, {x, y}, Integers]", "List[]");
    /* Underdetermined but integer-inconsistent (2 x + 2 y == 3 forces an odd
     * sum of evens): HNF's divisibility test proves no integer solution. */
    run_test("Solve[{2 x + 2 y == 3, x - y == 0}, {x, y}, Integers]", "List[]");
}

/* Factorable binary quadratic (Runge's simplest case): a single 2-variable
 * degree-2 equation whose quadratic part has a CROSS term and a perfect-square
 * discriminant delta = B^2 - 4 A C factors into two rational linear forms, so
 * the hyperbola has finitely many integer points found by a divisor
 * enumeration (`si_solve_factorable_conic`).  Exhaustive, so `{}` is a PROOF. */
static void test_factorable_conic(void) {
    /* x^2 + x y - 2 y^2 == (x - y)(x + 2 y); == 4 has six integer points. */
    run_test("Solve[x^2 + x y - 2 y^2 == 4, {x, y}, Integers]",
        "List[List[Rule[x, -3], Rule[y, 1]], List[Rule[x, -2], Rule[y, -1]], "
             "List[Rule[x, -2], Rule[y, 0]], List[Rule[x, 2], Rule[y, 0]], "
             "List[Rule[x, 2], Rule[y, 1]], List[Rule[x, 3], Rule[y, -1]]]");
    /* Same form == 15: a mod-3 obstruction means NO integer solution -- the
     * divisor enumeration proves it (not a decline). */
    run_test("Solve[(x - y) (x + 2 y) == 15, {x, y}, Integers]", "List[]");
    run_test("Solve[x^2 + x y - 2 y^2 == 15, {x, y}, Integers]", "List[]");
    /* Non-unit leading coefficients: 2 x^2 + 3 x y - 2 y^2 = (2 x - y)(x + 2 y),
     * delta = 25.  (si_solve_conic could not: it needs a unit Y^2 coefficient.) */
    run_test("Solve[2 x^2 + 3 x y - 2 y^2 == 7, {x, y}, Integers]",
        "List[List[Rule[x, -3], Rule[y, 1]], List[Rule[x, 3], Rule[y, -1]]]");
    /* Constraints filter the exhaustive set to the positive point. */
    run_test("Solve[x^2 + x y - 2 y^2 == 4 && x > 0 && y > 0, {x, y}, Integers]",
        "List[List[Rule[x, 2], Rule[y, 1]]]");
    /* Non-square discriminant (delta = 5) is a genuine Pell-type conic, not
     * factorable: unbounded -> left unevaluated, never a wrong {}. */
    run_test("Solve[x^2 + 3 x y + y^2 == 11, {x, y}, Integers]",
        "Solve[Equal[Plus[Power[x, 2], Times[3, Times[x, y]], Power[y, 2]], 11], "
             "List[x, y], Integers]");
}

/* Generalised Pell  x^2 - D y^2 == N  (N != +1) with x > 0 && y > 0, unbounded:
 * one parametric family per solution class (Nagell fundamentals + fundamental-
 * unit orbit, `si_solve_genpell_parametric`).  Validated exhaustively against a
 * brute-force positive-orthant enumeration over ~30 (D,N) pairs during
 * development; the property assertions here pin the class count, the class
 * fundamentals (C[1] -> 0), that every parameter value satisfies the equation,
 * and that an unsolvable equation is a PROOF (empty), never a wrong family. */
static void test_generalized_pell(void) {
    /* x^2 - 2 y^2 == 7: two classes, fundamentals (3,1) and (5,3). */
    run_test("Length[Solve[x^2 - 2 y^2 == 7 && x > 0 && y > 0, {x, y}, Integers]]",
        "2");
    run_test("Simplify[{x, y} /. (Solve[x^2 - 2 y^2 == 7 && x > 0 && y > 0, "
             "{x, y}, Integers][[1]] /. C[1] -> 0)]", "List[3, 1]");
    run_test("Simplify[{x, y} /. (Solve[x^2 - 2 y^2 == 7 && x > 0 && y > 0, "
             "{x, y}, Integers][[2]] /. C[1] -> 0)]", "List[5, 3]");
    /* The next orbit member of the first class is (13, 9). */
    run_test("Simplify[{x, y} /. (Solve[x^2 - 2 y^2 == 7 && x > 0 && y > 0, "
             "{x, y}, Integers][[1]] /. C[1] -> 1)]", "List[13, 9]");
    /* Every parameter value of every class satisfies x^2 - 2 y^2 == 7. */
    run_test("Simplify[(x^2 - 2 y^2) /. (Solve[x^2 - 2 y^2 == 7 && x > 0 && "
             "y > 0, {x, y}, Integers][[1]] /. C[1] -> 3)]", "7");
    run_test("Simplify[(x^2 - 2 y^2) /. (Solve[x^2 - 2 y^2 == 7 && x > 0 && "
             "y > 0, {x, y}, Integers][[2]] /. C[1] -> 5)]", "7");
    /* x^2 - 5 y^2 == 4 has three classes. */
    run_test("Length[Solve[x^2 - 5 y^2 == 4 && x > 0 && y > 0, {x, y}, Integers]]",
        "3");
    /* Provably unsolvable (a mod-8 obstruction) -> {}, not a decline. */
    run_test("Solve[x^2 - 2 y^2 == 5 && x > 0 && y > 0, {x, y}, Integers]",
        "List[]");
    /* Negative Pell N = -1: solvable over D = 13, fundamental (18, 5). */
    run_test("Length[Solve[x^2 - 13 y^2 == -1 && x > 0 && y > 0, {x, y}, Integers]]",
        "1");
    run_test("Simplify[{x, y} /. (Solve[x^2 - 13 y^2 == -1 && x > 0 && y > 0, "
             "{x, y}, Integers][[1]] /. C[1] -> 0)]", "List[18, 5]");
    run_test("Simplify[(x^2 - 13 y^2) /. (Solve[x^2 - 13 y^2 == -1 && x > 0 && "
             "y > 0, {x, y}, Integers][[1]] /. C[1] -> 2)]", "-1");
    /* Negative Pell unsolvable over D = 3 (even CF period) -> {} proof. */
    run_test("Solve[x^2 - 3 y^2 == -1 && x > 0 && y > 0, {x, y}, Integers]",
        "List[]");
    /* Without the positive-orthant constraints the unbounded family is declined
     * (left unevaluated), matching the N = +1 parametric convention. */
    run_test("Head[Solve[x^2 - 2 y^2 == 7, {x, y}, Integers]]", "Solve");
}

/* Two fourth powers a^4+b^4 as a sum in two distinct ways.  The smallest is
 * 635318657 = 59^4+158^4 = 133^4+134^4 (Euler), so the < 10^8 box is a TRUE
 * negative -- an empty set that must not be mistaken for a false negative.
 * The < 10^9 box reaches the Euler pair, guarding the search's completeness. */
static void test_two_fourth_powers(void) {
    run_test("Solve[x^4 + y^4 == z^4 + w^4 && 0 < x < y && 0 < z < w && x != z "
             "&& x^4 + y^4 < 10^8, {x, y, z, w}, Integers]", "List[]");
    run_test("Solve[x^4 + y^4 == z^4 + w^4 && 0 < x < y && 0 < z < w && x != z "
             "&& x^4 + y^4 < 10^9, {x, y, z, w}, Integers]",
        "List[List[Rule[x, 59], Rule[y, 158], Rule[z, 133], Rule[w, 134]], "
             "List[Rule[x, 133], Rule[y, 134], Rule[z, 59], Rule[w, 158]]]");
}

/* Bounded superelliptic / mixed-power boxes solved by the ordinary leaf search
 * (enumerate the bounded variables, solve the leaf by a perfect-power test). */
static void test_bounded_mixed_power(void) {
    /* Superelliptic y^3 = x^5 - x + 1 over |x|,|y| < 1000: three points. */
    run_test("Solve[y^3 == x^5 - x + 1 && Abs[x] < 1000 && Abs[y] < 1000, "
             "{x, y}, Integers]",
        "List[List[Rule[x, -1], Rule[y, 1]], "
             "List[Rule[x, 0], Rule[y, 1]], List[Rule[x, 1], Rule[y, 1]]]");
    /* x^2 + y^3 == z^7 over a box; z bounded by the equation magnitude. */
    run_test("Solve[x^2 + y^3 == z^7 && 1 < x < 1000 && 1 < y < 1000 && z > 1, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, 8], Rule[y, 4], Rule[z, 2]], "
             "List[Rule[x, 250], Rule[y, 25], Rule[z, 5]], "
             "List[Rule[x, 729], Rule[y, 162], Rule[z, 9]], "
             "List[Rule[x, 832], Rule[y, 112], Rule[z, 8]]]");
    /* 239 needs nine positive cubes (Waring g(3)=9 witness): two representations. */
    run_test("Solve[x1^3 + x2^3 + x3^3 + x4^3 + x5^3 + x6^3 + x7^3 + x8^3 + x9^3 "
             "== 239 && 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= x6 <= x7 <= x8 <= x9, "
             "{x1, x2, x3, x4, x5, x6, x7, x8, x9}, Integers]",
        "List[List[Rule[x1, 1], Rule[x2, 1], Rule[x3, 1], Rule[x4, 3], Rule[x5, 3], "
             "Rule[x6, 3], Rule[x7, 3], Rule[x8, 4], Rule[x9, 4]], "
             "List[Rule[x1, 1], Rule[x2, 2], Rule[x3, 2], Rule[x4, 2], Rule[x5, 2], "
             "Rule[x6, 3], Rule[x7, 3], Rule[x8, 3], Rule[x9, 5]]]");
}

/* Unconstrained sum of two squares.  A variable that appears only with even
 * exponents is sign-symmetric: derive_even_only_bounds bounds |x|,|y| and the
 * search covers both signs, so the full 12-element set is returned rather than
 * the solver declining and Solve fabricating {} over the Integers. */
static void test_sum_of_two_squares(void) {
    run_test("Solve[x^2 + y^2 == 25, {x, y}, Integers]",
        "List[List[Rule[x, -5], Rule[y, 0]], "
             "List[Rule[x, -4], Rule[y, -3]], List[Rule[x, -4], Rule[y, 3]], "
             "List[Rule[x, -3], Rule[y, -4]], List[Rule[x, -3], Rule[y, 4]], "
             "List[Rule[x, 0], Rule[y, -5]], List[Rule[x, 0], Rule[y, 5]], "
             "List[Rule[x, 3], Rule[y, -4]], List[Rule[x, 3], Rule[y, 4]], "
             "List[Rule[x, 4], Rule[y, -3]], List[Rule[x, 4], Rule[y, 3]], "
             "List[Rule[x, 5], Rule[y, 0]]]");
    /* The degenerate sum-of-even-powers == 0: only the origin (a positive-
     * definite term pinned to <= 0 forces its variable to 0). */
    run_test("Solve[x^2 + y^2 == 0, {x, y}, Integers]",
        "List[List[Rule[x, 0], Rule[y, 0]]]");
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
    TEST(test_exponential);
    TEST(test_superelliptic_bounded);
    TEST(test_pell);
    TEST(test_linear_parametric);
    TEST(test_integer_restriction);
    TEST(test_sum_of_two_squares);
    TEST(test_linear_system_integers_hnf);
    TEST(test_factorable_conic);
    TEST(test_generalized_pell);
    TEST(test_two_fourth_powers);
    TEST(test_bounded_mixed_power);
    TEST(test_euler_brick);
    TEST(test_quadratic_form);
    TEST(test_brocard);
    TEST(test_conic);
    TEST(test_pell_parametric);
    TEST(test_linear_system_ray);
    TEST(test_fixed_base_exponential);
    TEST(test_prouhet_tarry_escott);
    TEST(test_mordell);
    TEST(test_deferred_unevaluated);

    printf("All solve-integers tests passed!\n");
    return 0;
}
