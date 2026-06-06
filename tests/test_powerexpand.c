#include "expand_power.h"
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
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* Default (Assumptions -> Automatic) transforms. */
static void test_automatic(void) {
    run_test("PowerExpand[Sqrt[x y]]",
             "Times[Power[x, Rational[1, 2]], Power[y, Rational[1, 2]]]");
    run_test("PowerExpand[(x y)^a]", "Times[Power[x, a], Power[y, a]]");
    run_test("PowerExpand[Sqrt[z^2]]", "z");
    run_test("PowerExpand[(a^x)^y]", "Power[a, Times[x, y]]");
    /* (E^x)^y now stays nested through the core evaluator (it no longer
     * auto-folds to E^(x y)), so PowerExpand performs the merge itself. */
    run_test("PowerExpand[(E^x)^y]", "Power[E, Times[x, y]]");
    run_test("PowerExpand[Log[z^a]]", "Times[a, Log[z]]");
    run_test("PowerExpand[Log[1/z]]", "Times[-1, Log[z]]");
    run_test("PowerExpand[Log[x y]]", "Plus[Log[x], Log[y]]");
    run_test("PowerExpand[ArcTan[Tan[x]]]", "x");
    run_test("PowerExpand[ArcSin[Sin[x]]]", "x");
    /* Nested: Log of a power of a product collapses without expanding the
     * power first. */
    run_test("PowerExpand[Log[(a b)^c]]",
             "Times[Plus[Log[a], Log[b]], c]");
    /* A negative numeric coefficient under a root is folded into a Plus factor
     * so the root stays real: Sqrt[-4 Dt[u]^2 (-1+u)] -> 2 Dt[u] Sqrt[1-u]
     * rather than 2 I Dt[u] Sqrt[-1+u]. */
    run_test("PowerExpand[Sqrt[-4 Dt[u]^2 (-1 + u)]]",
             "Times[2, Dt[u], Power[Plus[1, Times[-1, u]], Rational[1, 2]]]");
    /* No Plus factor to absorb the sign -> genuinely imaginary, left as is. */
    run_test("PowerExpand[Sqrt[-4 Dt[u]^2]]",
             "Times[Complex[0, 2], Dt[u]]");
    /* Integer exponent: no folding (the sign cannot turn imaginary). */
    run_test("PowerExpand[(-2 (1 + x))^2]",
             "Times[4, Power[Plus[1, x], 2]]");
    /* No-op on a plain sum. */
    run_test("PowerExpand[x + y]", "Plus[x, y]");
    run_test("PowerExpand[Sqrt[x] + Sqrt[y]]",
             "Plus[Power[x, Rational[1, 2]], Power[y, Rational[1, 2]]]");
}

/* Variable-restricted form PowerExpand[expr, {x1, ...}]. */
static void test_restriction(void) {
    run_test("PowerExpand[Sqrt[a b] + Sqrt[c d], {a, b}]",
             "Plus[Times[Power[a, Rational[1, 2]], Power[b, Rational[1, 2]]], "
             "Power[Times[c, d], Rational[1, 2]]]");
    /* A single variable as the second argument. */
    run_test("PowerExpand[Sqrt[x y], y]",
             "Times[Power[x, Rational[1, 2]], Power[y, Rational[1, 2]]]");
    /* No listed variable present -> no expansion. */
    run_test("PowerExpand[Sqrt[c d], {a, b}]",
             "Power[Times[c, d], Rational[1, 2]]");
}

/* Threading over lists / equations / inequalities / logic. */
static void test_threading(void) {
    run_test("PowerExpand[{Sqrt[x y], Log[a b]}]",
             "List[Times[Power[x, Rational[1, 2]], Power[y, Rational[1, 2]]], "
             "Plus[Log[a], Log[b]]]");
    run_test("PowerExpand[Sqrt[x y] == Sqrt[z]]",
             "Equal[Times[Power[x, Rational[1, 2]], Power[y, Rational[1, 2]]], "
             "Power[z, Rational[1, 2]]]");
    run_test("PowerExpand[1 < Sqrt[x y] < 2]",
             "Inequality[1, Less, Times[Power[x, Rational[1, 2]], "
             "Power[y, Rational[1, 2]]], Less, 2]");
    run_test("PowerExpand[Log[a b] && Log[c d]]",
             "And[Plus[Log[a], Log[b]], Plus[Log[c], Log[d]]]");
}

/* Assumptions -> True: universally-correct branch-correction formulas. */
static void test_true(void) {
    run_test("PowerExpand[Sqrt[x y], Assumptions -> True]",
             "Times[Power[x, Rational[1, 2]], Power[y, Rational[1, 2]], "
             "Power[E, Times[Complex[0, 1], Pi, Floor[Plus[Rational[1, 2], "
             "Times[Rational[-1, 2], Times[Power[Pi, -1], "
             "Plus[Arg[x], Arg[y]]]]]]]]]");
    run_test("PowerExpand[Log[x y], Assumptions -> True]",
             "Plus[Log[x], Log[y], Times[Complex[0, 2], Times[Pi, "
             "Floor[Plus[Rational[1, 2], Times[Rational[-1, 2], "
             "Times[Power[Pi, -1], Plus[Arg[x], Arg[y]]]]]]]]]");
    run_test("PowerExpand[Arg[x y], Assumptions -> True]",
             "Plus[Arg[x], Arg[y], Times[2, Times[Pi, "
             "Floor[Plus[Rational[1, 2], Times[Rational[-1, 2], "
             "Times[Power[Pi, -1], Plus[Arg[x], Arg[y]]]]]]]]]");
}

/* Assumptions -> assum: refinement of the correction terms. */
static void test_assumptions(void) {
    run_test("PowerExpand[Sqrt[z^2], Assumptions -> z < 0]", "Times[-1, z]");
    run_test("PowerExpand[(z^p)^(1/p), Assumptions -> 0 < p < 1]", "z");
    run_test("PowerExpand[Log[(a b)^c], "
             "Assumptions -> 3 < a < 5 && -2 < b < -1 && 7 < c < 9]",
             "Plus[Times[Complex[0, -8], Pi], Times[Plus[Log[a], Log[b]], c]]");
    /* Positive real: correction collapses to the naive expansion. */
    run_test("PowerExpand[Sqrt[z^2], Assumptions -> z > 0]", "z");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_automatic);
    TEST(test_restriction);
    TEST(test_threading);
    TEST(test_true);
    TEST(test_assumptions);

    printf("All PowerExpand tests passed!\n");
    return 0;
}
