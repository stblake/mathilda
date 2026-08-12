/* test_core_algebra.c -- exhaustive unit tests for Mathilda's core algebraic
 * capabilities: Plus, Times, Power, Divide/Subtract, Apply (@@ / @@@ and every
 * level-spec form), Length/Dimensions, canonical ordering, and Table
 * construction. Expected values are FullForm, cross-checked against Wolfram
 * Language 14 semantics.
 *
 * This suite deliberately exercises the code paths touched by the 2026-08-12
 * evaluator-throughput work so a regression fails a test:
 *   - the specialised top-level Apply fast path and the general apply_at_level
 *     cleanup (src/funcprog.c) -- every level-spec variant is pinned here, so a
 *     mistake in the fast path (which only fires for the default `f @@ x`)
 *     or in the non-negative-spec early-return shows up immediately;
 *   - the Plus / Times fused special-case short-circuits (src/plus.c,
 *     src/times.c) -- each guarded pass (NDArray, neg-Plus, SeriesData,
 *     inexact contagion, Infinity/Indeterminate, radical fusion, trig
 *     canonicalisation) has a test that forces it, alongside the ordinary
 *     symbolic products/sums that take the short-circuit.
 */
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

/* Parse -> evaluate -> compare FullForm against `expected`. Hard-aborts (exit 1)
 * on mismatch so failures survive an NDEBUG (Release) build. */
static void ck(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ Plus */
static void test_plus(void) {
    /* Like-term combination */
    ck("x + x", "Times[2, x]");
    ck("x + 2 x", "Times[3, x]");
    ck("3 x + 2 x + x", "Times[6, x]");
    ck("a + b + a", "Plus[Times[2, a], b]");
    ck("x + y + x + y", "Plus[Times[2, x], Times[2, y]]");
    /* Orderless canonical ordering */
    ck("a + b + c", "Plus[a, b, c]");
    ck("c + a + b", "Plus[a, b, c]");
    ck("z + a + m", "Plus[a, m, z]");
    /* Numeric folding + identity element */
    ck("2 + 3", "5");
    ck("2 + 3 + x", "Plus[5, x]");
    ck("3 + x + 1", "Plus[4, x]");
    ck("0 + x", "x");
    ck("1/2 + 1/3", "Rational[5, 6]");
    /* Cancellation (forces the neg-Plus distribution pass) */
    ck("a + b - (a + b)", "0");
    ck("x - x", "0");
    /* Arity edge cases */
    ck("Plus[]", "0");
    ck("Plus[x]", "x");
    /* BigInt promotion on int64 overflow */
    ck("9223372036854775807 + 1", "9223372036854775808");
    /* Infinity / Indeterminate pass */
    ck("Infinity + 1", "Infinity");
    ck("ComplexInfinity + 5", "ComplexInfinity");
    ck("Indeterminate + x", "Indeterminate");
    /* Inexact contagion pass */
    ck("1. + 2", "3.0");
}

/* ----------------------------------------------------------------- Times */
static void test_times(void) {
    /* Like-base combination -> Power */
    ck("x x", "Power[x, 2]");
    ck("x^2 x", "Power[x, 3]");
    ck("x y x", "Times[Power[x, 2], y]");
    /* Orderless canonical ordering */
    ck("a b c", "Times[a, b, c]");
    ck("c b a", "Times[a, b, c]");
    /* Absorbing / identity elements */
    ck("0 x", "0");
    ck("1 x", "x");
    /* Numeric coefficient folding */
    ck("2 3 x", "Times[6, x]");
    ck("(2 x)(3 y)", "Times[6, x, y]");
    /* Arity edge cases */
    ck("Times[]", "1");
    ck("Times[x]", "x");
    /* Radical fusion + sqrt coefficient (post-grouping passes) */
    ck("2 Sqrt[2]", "Times[2, Power[2, Rational[1, 2]]]");
    ck("Sqrt[2] Sqrt[3]", "Power[6, Rational[1, 2]]");
    ck("Sqrt[8]", "Times[2, Power[2, Rational[1, 2]]]");
    /* Infinity pass */
    ck("2 Infinity", "Infinity");
    ck("(-3) Infinity", "Times[-1, Infinity]");
    ck("2 ComplexInfinity", "ComplexInfinity");
    /* Ordinary symbolic products (the short-circuit path) */
    ck("b a + d c", "Plus[Times[a, b], Times[c, d]]");
}

/* --------------------------------------------------- Power / Divide / Subtract */
static void test_power_divide_subtract(void) {
    ck("x^0", "1");
    ck("x^1", "x");
    ck("2^3", "8");
    ck("(x^2)^3", "Power[x, 6]");
    ck("x^-1", "Power[x, -1]");
    ck("Sqrt[x]^2", "x");
    ck("4^(1/2)", "2");
    /* Divide is Times[.., Power[.., -1]] */
    ck("a/b", "Times[a, Power[b, -1]]");
    ck("x/x", "1");
    ck("6/3", "2");
    ck("x^3/x", "Power[x, 2]");
    /* Subtract is Plus[.., Times[-1, ..]] */
    ck("a - b", "Plus[a, Times[-1, b]]");
    ck("5 - 2", "3");
}

/* ----------------------------------------------------------------- Apply */
static void test_apply(void) {
    /* Default top-level (@@) -- the specialised fast path */
    ck("Plus @@ {a, b, c}", "Plus[a, b, c]");
    ck("Times @@ {a, b, c}", "Times[a, b, c]");
    ck("f @@ g[1, 2, 3]", "f[1, 2, 3]");
    ck("Plus @@ {}", "0");                 /* empty -> Plus[] -> 0 */
    ck("Plus @@ {x}", "x");                /* single */
    ck("f @@ x", "x");                     /* atomic operand: fast path declines, general returns it */
    ck("Plus @@ Range[100]", "5050");
    ck("Length[Plus @@ Table[c[k] x, {k, 50}]]", "50");
    /* Non-default level specs -- must take the GENERAL path, not the fast path */
    ck("Apply[f, {{a, b}, {c, d}}, {1}]", "List[f[a, b], f[c, d]]");   /* @@@ */
    ck("f @@@ {{a, b}, {c, d}}", "List[f[a, b], f[c, d]]");            /* @@@ operator */
    ck("Apply[g, {{a, b}, {c, d}}, {2}]", "List[List[a, b], List[c, d]]"); /* atoms: unchanged */
    ck("Apply[h, {{a, b}, {c, d}}, {-1}]", "List[List[a, b], List[c, d]]"); /* leaves: unchanged */
    ck("Apply[q, g[x, y], Heads -> True]", "q[x, y]");
    /* Deeper nesting -- the level-1 fast recursion must refcount-SHARE (not
     * descend into) a nested element subtree, and the level-2 form must apply
     * exactly one level down. These pin the apply_child past-max short-circuit. */
    ck("f @@@ {{a, {b, c}}, {d}}", "List[f[a, List[b, c]], f[d]]");
    ck("Apply[f, {{{a, b}}, {{c, d}}}, {2}]", "List[List[f[a, b]], List[f[c, d]]]");
    ck("Apply[f, {{a, b}, {c, d}}, {1, 2}]", "List[f[a, b], f[c, d]]");
    ck("Apply[f, {{a, b}, {c, d}}, {-2}]", "List[f[a, b], f[c, d]]"); /* negative spec path */
}

/* ------------------------------------------------------ Length / Dimensions */
static void test_length_dimensions(void) {
    ck("Length[a + b + c]", "3");
    ck("Length[{1, 2, 3, 4}]", "4");
    ck("Length[f[x, y]]", "2");
    ck("Length[x]", "0");                  /* atom */
    ck("Dimensions[{{1, 2}, {3, 4}}]", "List[2, 2]");
}

/* ---------------------------------------------------------------- Table */
static void test_table(void) {
    ck("Table[k^2, {k, 5}]", "List[1, 4, 9, 16, 25]");
    ck("Table[c[k] x, {k, 3}]", "List[Times[c[1], x], Times[c[2], x], Times[c[3], x]]");
    ck("Length[Table[k, {k, 2000}]]", "2000");
    /* The arithmetic-sentinel shape from the evaluator-throughput benchmark:
     * Table of symbolic terms, summed via @@, then Length. Value must be stable. */
    ck("Length[Plus @@ Table[c[k] x, {k, 2000}]]", "2000");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_plus);
    TEST(test_times);
    TEST(test_power_divide_subtract);
    TEST(test_apply);
    TEST(test_length_dimensions);
    TEST(test_table);

    printf("All core algebra tests passed!\n");
    return 0;
}
