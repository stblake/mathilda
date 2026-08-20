/*
 * test_solve_integers.c
 *
 * Unit tests for the integer-domain (Diophantine) solver
 * (src/solve/): Solve[eqns && constraints, vars, Integers].
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
#include <math.h>

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

/* Abs-value ordering constraints.  |x| < |y| < |z| < B is the natural way to
 * ask for the ordered representatives of the sum-of-three-cubes solution set
 * (one per permutation orbit) instead of every permutation; the magnitude chain
 * bounds every variable and filters the result to the ordered subset. */
static void test_abs_ordering(void) {
    /* Full ordering: the six ordered representatives of x^3+y^3+z^3 == 63. */
    run_test("Solve[x^3 + y^3 + z^3 == 63 && Abs[x] < Abs[y] < Abs[z] < 10000, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, -4178], Rule[y, -6034], Rule[z, 6639]], "
             "List[Rule[x, -38], Rule[y, -58], Rule[z, 63]], "
             "List[Rule[x, -37], Rule[y, -63], Rule[z, 67]], "
             "List[Rule[x, -4], Rule[y, -6], Rule[z, 7]], "
             "List[Rule[x, 0], Rule[y, -1], Rule[z, 4]], "
             "List[Rule[x, 102], Rule[y, 146], Rule[z, -161]]]");
    /* Partial ordering |x| < |y| with z free: the |x|<|y| subset (z arbitrary). */
    run_test("Solve[x^3 + y^3 + z^3 == 63 && Abs[x] < Abs[y] < 10000 && Abs[z] < 10000, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, -6034], Rule[y, 6639], Rule[z, -4178]], "
             "List[Rule[x, -4178], Rule[y, -6034], Rule[z, 6639]], "
             "List[Rule[x, -4178], Rule[y, 6639], Rule[z, -6034]], "
             "List[Rule[x, -63], Rule[y, 67], Rule[z, -37]], "
             "List[Rule[x, -58], Rule[y, 63], Rule[z, -38]], "
             "List[Rule[x, -38], Rule[y, -58], Rule[z, 63]], "
             "List[Rule[x, -38], Rule[y, 63], Rule[z, -58]], "
             "List[Rule[x, -37], Rule[y, -63], Rule[z, 67]], "
             "List[Rule[x, -37], Rule[y, 67], Rule[z, -63]], "
             "List[Rule[x, -6], Rule[y, 7], Rule[z, -4]], "
             "List[Rule[x, -4], Rule[y, -6], Rule[z, 7]], "
             "List[Rule[x, -4], Rule[y, 7], Rule[z, -6]], "
             "List[Rule[x, -1], Rule[y, 4], Rule[z, 0]], "
             "List[Rule[x, 0], Rule[y, -1], Rule[z, 4]], "
             "List[Rule[x, 0], Rule[y, 4], Rule[z, -1]], "
             "List[Rule[x, 102], Rule[y, -161], Rule[z, 146]], "
             "List[Rule[x, 102], Rule[y, 146], Rule[z, -161]], "
             "List[Rule[x, 146], Rule[y, -161], Rule[z, 102]]]");
    /* Small box below the Booker gate: exercises the abs-ordering leaf search
     * (effective_bounds pruning + abs-aware enumeration order). */
    run_test("Solve[x^3 + y^3 + z^3 == 63 && Abs[x] < Abs[y] < Abs[z] < 10, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, -4], Rule[y, -6], Rule[z, 7]], "
             "List[Rule[x, 0], Rule[y, -1], Rule[z, 4]]]");
    /* Non-strict <= keeps ties: (1,1,1) and (4,4,-5) for k == 3. */
    run_test("Solve[x^3 + y^3 + z^3 == 3 && Abs[x] <= Abs[y] <= Abs[z] < 50, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, 1], Rule[y, 1], Rule[z, 1]], "
             "List[Rule[x, 4], Rule[y, 4], Rule[z, -5]]]");
    /* Strict ordering rejects those ties -> empty (1=1=1, 4=4<5 both fail). */
    run_test("Solve[x^3 + y^3 + z^3 == 3 && Abs[x] < Abs[y] < Abs[z] < 50, "
             "{x, y, z}, Integers]", "List[]");
    /* Ordering as a separate chain from the bound, box just above the gate. */
    run_test("Solve[x^3 + y^3 + z^3 == 63 && Abs[x] < Abs[y] < Abs[z] && Abs[z] < 200, "
             "{x, y, z}, Integers]",
        "List[List[Rule[x, -38], Rule[y, -58], Rule[z, 63]], "
             "List[Rule[x, -37], Rule[y, -63], Rule[z, 67]], "
             "List[Rule[x, -4], Rule[y, -6], Rule[z, 7]], "
             "List[Rule[x, 0], Rule[y, -1], Rule[z, 4]], "
             "List[Rule[x, 102], Rule[y, 146], Rule[z, -161]]]");
    /* No small representation of 42 -> empty even with the ordering. */
    run_test("Solve[x^3 + y^3 + z^3 == 42 && Abs[x] < Abs[y] < Abs[z] < 8000, "
             "{x, y, z}, Integers]", "List[]");
}

/* Integer cube root of t (t may be negative); writes it to *zout and returns 1
 * iff t is a perfect cube. */
static int test_int_cube_root(long long t, long long* zout) {
    long long s = t < 0 ? -1 : 1;
    long long a = t < 0 ? -t : t;
    long long r = (long long)cbrtl((long double)a);
    for (long long z = (r > 2 ? r - 2 : 0); z <= r + 2; z++)
        if (z * z * z == a) { *zout = s * z; return 1; }
    return 0;
}

/* Independent cross-check of the Booker cube-root-mod-d engine: brute-force the
 * box [-(B-1), B-1]^3, format the solution set exactly as build_result does
 * (ascending (x,y,z) tuples), and assert Solve returns it.  B is chosen just
 * above the Booker engagement gate (2B > 200) so the engine runs, yet small
 * enough to brute-force here -- proving completeness, not merely locking the
 * current output. */
static void check_three_cubes_brute(long long k, long long B) {
    long long lo = -(B - 1), hi = B - 1;
    size_t cap = 256, n = 0;
    long long (*sol)[3] = (long long(*)[3])malloc(cap * sizeof(*sol));
    for (long long x = lo; x <= hi; x++)
        for (long long y = lo; y <= hi; y++) {
            long long t = k - x * x * x - y * y * y, z;
            if (!test_int_cube_root(t, &z)) continue;
            if (z < lo || z > hi) continue;
            if (n == cap) { cap *= 2; sol = (long long(*)[3])realloc(sol, cap * sizeof(*sol)); }
            sol[n][0] = x; sol[n][1] = y; sol[n][2] = z; n++;
        }
    /* build expected FullForm string */
    size_t cap_s = n * 72 + 16;
    char* expected = (char*)malloc(cap_s);
    size_t p = 0;
    p += (size_t)snprintf(expected + p, cap_s - p, "List[");
    for (size_t i = 0; i < n; i++)
        p += (size_t)snprintf(expected + p, cap_s - p,
            "%sList[Rule[x, %lld], Rule[y, %lld], Rule[z, %lld]]",
            i ? ", " : "", sol[i][0], sol[i][1], sol[i][2]);
    p += (size_t)snprintf(expected + p, cap_s - p, "]");
    char input[192];
    snprintf(input, sizeof input,
        "Solve[x^3 + y^3 + z^3 == %lld && Abs[x] < %lld && Abs[y] < %lld "
        "&& Abs[z] < %lld, {x, y, z}, Integers]", k, B, B, B);
    run_test(input, expected);
    free(sol); free(expected);
}

/* Run the brute cross-check on several (k, box) just above the Booker gate. */
static void test_booker_brute_crosscheck(void) {
    long long ks[] = {2, 3, 29, -13, 9};
    long long Bs[] = {105, 130};
    for (size_t i = 0; i < sizeof(ks) / sizeof(ks[0]); i++)
        for (size_t j = 0; j < sizeof(Bs) / sizeof(Bs[0]); j++)
            check_three_cubes_brute(ks[i], Bs[j]);
}

/* Fermat's Last Theorem short-circuit: x^n + y^n == z^n (n >= 3, all strictly
 * positive) has no solutions, decided instantly -- both for a bounded box and
 * for the unbounded ("indefinite") form, and independent of which side carries
 * the isolated cube.  The guards below confirm it does NOT over-fire. */
static void test_fermat_last_theorem(void) {
    /* Bounded: the whole 10^4 box is dismissed without a search. */
    run_test("Solve[x^3 + y^3 == z^3 && z > x > y > 0 && x < 10000 && y < 10000 "
             "&& z < 10000, {x, y, z}, Integers]", "List[]");
    /* Indefinite (no upper bound): still decided, not left unevaluated. */
    run_test("Solve[x^3 + y^3 == z^3 && x > 0 && y > 0 && z > 0, "
             "{x, y, z}, Integers]", "List[]");
    /* Isolated cube on the left. */
    run_test("Solve[z^3 == x^3 + y^3 && x > 0 && y > 0 && z > 0, "
             "{x, y, z}, Integers]", "List[]");
    /* General exponent n >= 3. */
    run_test("Solve[x^5 + y^5 == z^5 && x > 0 && y > 0 && z > 0, "
             "{x, y, z}, Integers]", "List[]");

    /* Guards -- must NOT be swallowed by the FLT short-circuit. */
    /* n = 2 has Pythagorean triples. */
    run_test("Length[Solve[x^2 + y^2 == z^2 && 0 < x < y && z < 20 && y < 20, "
             "{x, y, z}, Integers]]", "11");
    /* k != 0 (x^3 + y^3 == z^3 + 1) is a different equation with solutions. */
    run_test("Length[Solve[x^3 + y^3 == z^3 + 1 && 0 < x <= 2 && 0 < y <= 2 "
             "&& 0 < z <= 2, {x, y, z}, Integers]]", "3");
    /* A box admitting 0 / negatives keeps the (0, a, a) family. */
    run_test("Length[Solve[x^3 + y^3 == z^3 && -2 <= x <= 2 && -2 <= y <= 2 "
             "&& -2 <= z <= 2, {x, y, z}, Integers]]", "13");
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

/* Definite binary quadratic (ellipse): a 2-variable degree-2 equation with a
 * negative discriminant B^2 - 4 A C < 0 is a compact conic with finitely many
 * integer points -- but a ROTATED one (cross term) is not bounded by the
 * interval bounder.  `si_solve_elliptic_bqf` treats it as a quadratic in x per
 * fixed y over the finite y-interval where an x is real, exhaustively, so `{}`
 * is a proof.  All counts below were checked against brute force over |.|<=200. */
static void test_elliptic_bqf(void) {
    /* Rotated ellipse x^2 + x y + y^2 == 3: six points. */
    run_test("Solve[x^2 + x y + y^2 == 3, {x, y}, Integers]",
        "List[List[Rule[x, -2], Rule[y, 1]], List[Rule[x, -1], Rule[y, -1]], "
             "List[Rule[x, -1], Rule[y, 2]], List[Rule[x, 1], Rule[y, -2]], "
             "List[Rule[x, 1], Rule[y, 1]], List[Rule[x, 2], Rule[y, -1]]]");
    run_test("Length[Solve[x^2 + x y + y^2 == 7, {x, y}, Integers]]", "12");
    run_test("Length[Solve[x^2 + x y + y^2 == 91, {x, y}, Integers]]", "24");
    /* Non-unit rotated ellipse. */
    run_test("Solve[3 x^2 + 2 x y + 3 y^2 == 24, {x, y}, Integers]",
        "List[List[Rule[x, -3], Rule[y, 1]], List[Rule[x, -1], Rule[y, 3]], "
             "List[Rule[x, 1], Rule[y, -3]], List[Rule[x, 3], Rule[y, -1]]]");
    /* Negative-definite form is normalised to positive-definite. */
    run_test("Length[Solve[-x^2 - x y - y^2 == -7, {x, y}, Integers]]", "12");
    /* Definite form with linear terms. */
    run_test("Solve[x^2 + x y + y^2 - 3 x == 7, {x, y}, Integers]", "List[]");
    /* Provably no representation (exhaustive over the ellipse) -> {} proof. */
    run_test("Solve[x^2 + x y + y^2 == 2, {x, y}, Integers]", "List[]");
    run_test("Solve[5 x^2 + 6 x y + 5 y^2 == 8, {x, y}, Integers]", "List[]");
    run_test("Solve[7 x^2 - 5 x y + 7 y^2 == 39, {x, y}, Integers]", "List[]");
    /* A constraint filters the exhaustive set. */
    run_test("Length[Solve[x^2 + x y + y^2 == 7 && x > 0, {x, y}, Integers]]", "6");
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

/* Ordering-aware, 128-bit meet-in-the-middle for a separable power sum.
 * Lander-Parkin  x^5+y^5+z^5+w^5 == r^5  over a box whose UNORDERED iterate side
 * exceeds SI_MAX_NODES (so the plain int64 mitm_solve declines) but whose ordered
 * combinations fit: si_solve_separable_mitm returns the complete set -- the
 * (27,84,110,133;144) counterexample and its 2x, 3x, 4x multiples below 700. */
static void test_power_sum_mitm_quintic(void) {
    run_test("Solve[x^5 + y^5 + z^5 + w^5 == r^5 && 0 < x < y < z < w < r < 700, "
             "{x, y, z, w, r}, Integers]",
        "List[List[Rule[x, 27], Rule[y, 84], Rule[z, 110], Rule[w, 133], Rule[r, 144]], "
             "List[Rule[x, 54], Rule[y, 168], Rule[z, 220], Rule[w, 266], Rule[r, 288]], "
             "List[Rule[x, 81], Rule[y, 252], Rule[z, 330], Rule[w, 399], Rule[r, 432]], "
             "List[Rule[x, 108], Rule[y, 336], Rule[z, 440], Rule[w, 532], Rule[r, 576]]]");
}

/* Frye's biquadrate search  x^4+y^4+z^4 == w^4:  the minimal counterexample to
 * Euler's conjecture, 95800^4+217519^4+414560^4 = 422481^4.  A narrow window
 * around w brackets it so the full modular-sieve + decompose machinery runs in
 * well under a second (the whole 0<w<10^6 search is a ~10-minute single-core
 * computation -- Frye needed a Connection Machine for it in 1988). */
static void test_frye_biquadrate(void) {
    run_test("Solve[x^4 + y^4 + z^4 == w^4 && 0 < x < y < z < w && 422400 < w < 422500, "
             "{x, y, z, w}, Integers]",
        "List[List[Rule[x, 95800], Rule[y, 217519], Rule[z, 414560], Rule[w, 422481]]]");
}

/* Modular-sieved leaf search for a large NON-separable box: the ordering-pruned
 * box (0<x<y<z<w<600, a cross-term equation) exceeds SI_MAX_NODES so the ordinary
 * leaf search declines; the residue sieve makes it exhaustible.  Self-validating:
 * the result restricted to the smaller box the ORDINARY engine solves exactly
 * must match that engine -- a completeness proof independent of any magic count. */
static void test_box_modsieve_nonseparable(void) {
    run_test("With[{a = Solve[x^2 + y^2 + z^2 + x y == w^2 && 0 < x < y < z < w < 500, "
             "{x, y, z, w}, Integers], "
             "b = Solve[x^2 + y^2 + z^2 + x y == w^2 && 0 < x < y < z < w < 600, "
             "{x, y, z, w}, Integers]}, "
             "Sort[a] === Sort[Select[b, (w /. #) < 500 &]]]", "True");
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
    TEST(test_abs_ordering);
    TEST(test_booker_brute_crosscheck);
    TEST(test_fermat_last_theorem);
    TEST(test_exponential);
    TEST(test_superelliptic_bounded);
    TEST(test_pell);
    TEST(test_linear_parametric);
    TEST(test_integer_restriction);
    TEST(test_sum_of_two_squares);
    TEST(test_linear_system_integers_hnf);
    TEST(test_factorable_conic);
    TEST(test_generalized_pell);
    TEST(test_elliptic_bqf);
    TEST(test_two_fourth_powers);
    TEST(test_power_sum_mitm_quintic);
    TEST(test_frye_biquadrate);
    TEST(test_box_modsieve_nonseparable);
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
