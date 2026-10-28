/* FreeQ and MemberQ: an exhaustive table of reference behaviour.
 *
 * Every expected value below was produced by the reference system (14.0) on
 * the identical input, not written by hand.  The two heads differ in exactly
 * the way the table probes:
 *   - FreeQ asks whether `form` occurs anywhere a rule could rewrite it, so a
 *     Flat call form is also found as a SUBSET of a longer call on the same
 *     head (FreeQ[a + b + c, a + b] is False; issue #115).  Default level
 *     {0, Infinity}, Heads -> True.
 *   - MemberQ is a whole-element test (MatchQ semantics) at level 1 by
 *     default, Heads -> False.
 * Both are checked on atoms, patterns, level specs (positive, negative,
 * Infinity in either slot), heads, Rational/Complex atoms, Hold, associations,
 * packed arrays, operator forms, and calls wider than 512 arguments (the
 * matcher's old per-argument depth limit; issue #116). */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { const char* input; const char* expected; } Case;

static void run_cases(const Case* cs, size_t n) {
    for (size_t i = 0; i < n; i++) assert_eval_eq(cs[i].input, cs[i].expected, 0);
}

/* flatf is Flat, flatof Flat and Orderless; ordf has neither. */
static void setup_heads(void) {
    assert_eval_eq("SetAttributes[flatf, Flat]; SetAttributes[flatof, {Flat, Orderless}]; "
                   "{Attributes[flatf], Attributes[flatof], Attributes[ordf]}",
                   "{{Flat}, {Flat, Orderless}, {}}", 0);
}
static void test_freeq_atoms(void) {
    static const Case cs[] = {
        { "FreeQ[x, x]", "False" },
        { "FreeQ[x, y]", "True" },
        { "FreeQ[1, 1]", "False" },
        { "FreeQ[1.5, 1.5]", "False" },
        { "FreeQ[1.5, 1]", "True" },
        { "FreeQ[\"abc\", \"abc\"]", "False" },
        { "FreeQ[\"abc\", \"b\"]", "True" },
        { "FreeQ[{}, x]", "True" },
        { "FreeQ[{}, List]", "False" },
        { "FreeQ[{}, List, Heads -> False]", "True" },
        { "FreeQ[f[], f]", "False" },
        { "FreeQ[f[], f, Heads -> False]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_patterns(void) {
    static const Case cs[] = {
        { "FreeQ[{a, 1, b}, _Integer]", "False" },
        { "FreeQ[{a, 1.5, b}, _Integer]", "True" },
        { "FreeQ[{a, 1.5, b}, _Real]", "False" },
        { "FreeQ[{a, \"s\"}, _String]", "False" },
        { "FreeQ[{1, 2, 3}, x_ /; x > 2]", "False" },
        { "FreeQ[{1, 2, 3}, x_ /; x > 3]", "True" },
        { "FreeQ[{1, 3, 5}, _?EvenQ]", "True" },
        { "FreeQ[{1, 3, 6}, _?EvenQ]", "False" },
        { "FreeQ[{a, b}, c | d]", "True" },
        { "FreeQ[{a, b}, c | b]", "False" },
        { "FreeQ[{a, b}, Except[a | b | List]]", "False" },
        { "FreeQ[{a, b}, Except[a | List]]", "False" },
        { "FreeQ[f[g[1, 2]], g[_, _]]", "False" },
        { "FreeQ[f[g[1, 2]], g[_]]", "True" },
        { "FreeQ[f[g[1, 1]], g[x_, x_]]", "False" },
        { "FreeQ[f[g[1, 2]], g[x_, x_]]", "True" },
        { "FreeQ[f[g[1, 2, 3]], g[__]]", "False" },
        { "FreeQ[f[g[]], g[__]]", "True" },
        { "FreeQ[f[g[]], g[___]]", "False" },
        { "FreeQ[x^2 + y, x^_]", "False" },
        { "FreeQ[x^2 + y, z^_]", "True" },
        { "FreeQ[Sin[x] + Cos[y], _Sin]", "False" },
        { "FreeQ[Sin[x] + Cos[y], _Tan]", "True" },
        { "FreeQ[f[x], HoldPattern[f[_]]]", "False" },
        { "FreeQ[f[x_], x]", "False" },
        { "FreeQ[f[x_], Verbatim[x_]]", "False" },
        { "FreeQ[f[y], Verbatim[_]]", "True" },
        { "FreeQ[f[_], Verbatim[_]]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_levels(void) {
    static const Case cs[] = {
        { "FreeQ[{1, 2, 4, 1, 0}, 0, {2}]", "True" },
        { "FreeQ[{f[0], 2}, 0, {1}]", "True" },
        { "FreeQ[{f[0], 2}, 0, {2}]", "False" },
        { "FreeQ[{f[0], 2}, 0, 1]", "True" },
        { "FreeQ[{f[0], 2}, 0, 2]", "False" },
        { "FreeQ[{f[0], 2}, 0, Infinity]", "False" },
        { "FreeQ[{f[0], 2}, 0, {0, 1}]", "True" },
        { "FreeQ[{f[0], 2}, 0, {2, Infinity}]", "False" },
        { "FreeQ[{f[0], 2}, 0, {3, Infinity}]", "True" },
        { "FreeQ[{f[0], 2}, 0, {-1}]", "False" },
        { "FreeQ[{f[0], 2}, f[0], {-1}]", "True" },
        { "FreeQ[{f[0], 2}, f[0], {-2}]", "False" },
        { "FreeQ[{f[0], 2}, f[0], {-2, -1}]", "False" },
        { "FreeQ[{f[0], 2}, {f[0], 2}, {0}]", "False" },
        { "FreeQ[{f[0], 2}, {f[0], 2}, {1}]", "True" },
        { "FreeQ[{f[0], 2}, {f[0], 2}, {0, Infinity}]", "False" },
        { "FreeQ[{f[0], 2}, 0, {0, Infinity}]", "False" },
        { "FreeQ[g[h[k[1]]], 1, {3}]", "False" },
        { "FreeQ[g[h[k[1]]], 1, {4}]", "True" },
        { "FreeQ[g[h[k[1]]], 1, {2, 3}]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_heads(void) {
    static const Case cs[] = {
        { "FreeQ[f[x], f]", "False" },
        { "FreeQ[f[x], f, Heads -> False]", "True" },
        { "FreeQ[f[x], f, Heads -> True]", "False" },
        { "FreeQ[f[g][x], g]", "False" },
        { "FreeQ[f[g][x], g, Heads -> False]", "True" },
        { "FreeQ[f[g][x], f[g]]", "False" },
        { "FreeQ[f[g[a]], a, Heads -> False]", "False" },
        { "FreeQ[f[g[a]], g, Heads -> False]", "True" },
        { "FreeQ[f[g[a]], g, {1}]", "True" },
        { "FreeQ[f[g[a]], g, {1}, Heads -> False]", "True" },
        { "FreeQ[f[g[a]], f, {1}]", "False" },
        { "FreeQ[f[g[a]], g[a], {2, Infinity}]", "True" },
        { "FreeQ[f[g[a]], a, {2, Infinity}]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_numbers(void) {
    static const Case cs[] = {
        { "FreeQ[1/2, 2]", "True" },
        { "FreeQ[f[1/2], 2]", "True" },
        { "FreeQ[f[1/2], 1/2]", "False" },
        { "FreeQ[f[1/2], Rational]", "False" },
        { "FreeQ[f[1/2], Rational, Heads -> False]", "True" },
        { "FreeQ[1/2, Rational, Heads -> False]", "True" },
        { "FreeQ[I x, Complex]", "False" },
        { "FreeQ[3 + 4 I, 4]", "True" },
        { "FreeQ[f[3 + 4 I], _Complex]", "False" },
        { "FreeQ[f[2/3], _Rational]", "False" },
        { "FreeQ[{1, 2.}, _Integer]", "False" },
        { "FreeQ[{1., 2.}, _Integer]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_flat_subsets(void) {
    static const Case cs[] = {
        { "FreeQ[a f[1] f[2], x_f y_f]", "False" },
        { "FreeQ[a f[1] f[2], x_f y_f, {0}]", "False" },
        { "FreeQ[a f[1] f[2], x_f y_f, {0, Infinity}]", "False" },
        { "FreeQ[a f[1] f[2], x_f y_f, {1}]", "True" },
        { "FreeQ[a f[1] f[2], x_f y_f z_f]", "True" },
        { "FreeQ[a f[1] g[2], x_f y_f]", "True" },
        { "FreeQ[a + b + c, a + b]", "False" },
        { "FreeQ[a + b + c, a + d]", "True" },
        { "FreeQ[a + b + c, a + c]", "False" },
        { "FreeQ[a + b + c, a + b, {1}]", "True" },
        { "FreeQ[g[a + b + c], a + b]", "False" },
        { "FreeQ[g[a + b + c], a + b, {1}]", "False" },
        { "FreeQ[g[a + b + c], a + b, {2}]", "True" },
        { "FreeQ[a b c, a c]", "False" },
        { "FreeQ[2 a b, 2 b]", "False" },
        { "FreeQ[x^2 + 2 x y + y^2, _ x y]", "False" },
        { "FreeQ[a + b + c, _Symbol + _Symbol]", "True" },
        { "FreeQ[a + b + 1, _Integer + _Symbol]", "False" },
        { "FreeQ[a + b + c, (x_ + y_) /; x === c]", "False" },
        { "FreeQ[a + b + c, (x_ + y_) /; x === d]", "True" },
        { "FreeQ[d ** a ** b ** c, a ** b]", "False" },
        { "FreeQ[d ** a ** b ** c, b ** a]", "True" },
        { "FreeQ[d ** a ** b ** c, a ** c]", "True" },
        { "FreeQ[d ** a ** b ** c, b ** c]", "False" },
        { "FreeQ[d . a . b . c, a . b]", "False" },
        { "FreeQ[d . a . b . c, d . b]", "True" },
        { "FreeQ[a && b && c, a && c]", "True" },
        { "FreeQ[a && b && c, b && c]", "False" },
        { "FreeQ[a || b || c, a || b]", "False" },
        { "FreeQ[a || b || c, a || c]", "True" },
        { "FreeQ[flatf[a, b, c], flatf[a, b]]", "False" },
        { "FreeQ[flatf[a, b, c], flatf[a, c]]", "True" },
        { "FreeQ[flatf[a, b, c], flatf[b, c]]", "False" },
        { "FreeQ[h[flatf[a, b, c]], flatf[b, c]]", "False" },
        { "FreeQ[h[flatf[a, b, c]], flatf[b, c], {1}]", "False" },
        { "FreeQ[flatof[a, b, c], flatof[a, c]]", "False" },
        { "FreeQ[flatof[a, b, c], flatof[c, a]]", "False" },
        { "FreeQ[flatof[a, b, c], flatof[a, d]]", "True" },
        { "FreeQ[ordf[a, b, c], ordf[a, b]]", "True" },
        { "FreeQ[g[a, b, c], g[a, b]]", "True" },
        { "FreeQ[{a, b, c}, {a, b}]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_holds_and_associations(void) {
    static const Case cs[] = {
        { "FreeQ[Hold[x], x]", "False" },
        { "FreeQ[Hold[1 + 1], 2]", "True" },
        { "FreeQ[Hold[1 + 1], HoldPattern[1 + 1]]", "False" },
        { "FreeQ[x -> y, y]", "False" },
        { "FreeQ[x -> y, Rule]", "False" },
        { "FreeQ[<|a -> 1|>, a]", "True" },
        { "FreeQ[<|a -> 1|>, 1]", "False" },
        { "FreeQ[<|a -> 1|>, 2]", "True" },
        { "FreeQ[{<|a -> 1|>}, a]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_packed_and_wide(void) {
    static const Case cs[] = {
        { "FreeQ[Range[10], 5]", "False" },
        { "FreeQ[Range[10], 11]", "True" },
        { "FreeQ[Range[10], _Real]", "True" },
        { "FreeQ[Range[10], _Integer]", "False" },
        { "FreeQ[N[Range[5]], 3.]", "False" },
        { "FreeQ[N[Range[5]], 3]", "True" },
        { "FreeQ[N[Range[5]], _Integer]", "True" },
        { "FreeQ[Range[10000], 10001]", "True" },
        { "FreeQ[Range[10000], 10000]", "False" },
        { "FreeQ[Table[{i, i^2}, {i, 100}], 10000]", "False" },
        { "FreeQ[Table[{i, i^2}, {i, 100}], 10001]", "True" },
        { "FreeQ[Sum[f[i], {i, 1000}], f[1000]]", "False" },
        { "FreeQ[Sum[f[i], {i, 1000}], f[1001]]", "True" },
        { "FreeQ[Sum[f[i], {i, 1000}], Sum[f[i], {i, 1000}]]", "False" },
        { "FreeQ[Sum[f[i], {i, 1000}], f[1] + f[1000]]", "False" },
        { "FreeQ[Product[f[i], {i, 600}], Product[f[i], {i, 600}]]", "False" },
        { "FreeQ[g @@ Range[2000], g @@ Range[2000]]", "False" },
        { "FreeQ[k[g @@ Range[2000]], g @@ Table[_, 2000]]", "False" },
        { "FreeQ[k[g @@ Range[2000]], g @@ Table[_, 1999]]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_freeq_operator_form(void) {
    static const Case cs[] = {
        { "FreeQ[x][f[x]]", "False" },
        { "FreeQ[y][f[x]]", "True" },
        { "FreeQ[_Integer][{a, 1}]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_basic(void) {
    static const Case cs[] = {
        { "MemberQ[{1, 2, 3}, 2]", "True" },
        { "MemberQ[{1, 2, 3}, 4]", "False" },
        { "MemberQ[{}, x]", "False" },
        { "MemberQ[{1, 2, 3}, _Integer]", "True" },
        { "MemberQ[{1, 2, 3}, _Real]", "False" },
        { "MemberQ[{a, b}, a | c]", "True" },
        { "MemberQ[{a, b}, c | d]", "False" },
        { "MemberQ[{a, b}, Except[a]]", "True" },
        { "MemberQ[{a, a}, Except[a]]", "False" },
        { "MemberQ[{1, 5, 9}, x_ /; x > 8]", "True" },
        { "MemberQ[{1, 5, 9}, x_ /; x > 9]", "False" },
        { "MemberQ[{1, 2}, _?(# > 1 &)]", "True" },
        { "MemberQ[{1, 2}, _?(# > 5 &)]", "False" },
        { "MemberQ[{f[1], g[2]}, f[_]]", "True" },
        { "MemberQ[{f[1], g[2]}, h[_]]", "False" },
        { "MemberQ[{f[1, 1]}, f[x_, x_]]", "True" },
        { "MemberQ[{f[1, 2]}, f[x_, x_]]", "False" },
        { "MemberQ[{1.0}, 1]", "False" },
        { "MemberQ[{1}, 1.0]", "False" },
        { "MemberQ[{\"a\", \"b\"}, \"a\"]", "True" },
        { "MemberQ[{\"a\", \"b\"}, \"c\"]", "False" },
        { "MemberQ[{\"a\", \"b\"}, _String]", "True" },
        { "MemberQ[{x_}, x]", "False" },
        { "MemberQ[{x_}, Verbatim[x_]]", "True" },
        { "MemberQ[{x}, Verbatim[x_]]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_non_list_heads_and_atoms(void) {
    static const Case cs[] = {
        { "MemberQ[f[a, b], a]", "True" },
        { "MemberQ[f[a, b], c]", "False" },
        { "MemberQ[a + b, a]", "True" },
        { "MemberQ[x^2, 2]", "True" },
        { "MemberQ[x, x]", "False" },
        { "MemberQ[5, 5]", "False" },
        { "MemberQ[\"abc\", \"abc\"]", "False" },
        { "MemberQ[1/2, 2]", "False" },
        { "MemberQ[1/2, 1]", "False" },
        { "MemberQ[3 + 4 I, 4]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_levels(void) {
    static const Case cs[] = {
        { "MemberQ[{{1}}, 1]", "False" },
        { "MemberQ[{{1}}, 1, 2]", "True" },
        { "MemberQ[{{1}}, 1, {2}]", "True" },
        { "MemberQ[{{1}}, 1, {1}]", "False" },
        { "MemberQ[{{1}}, 1, Infinity]", "True" },
        { "MemberQ[{{1}}, {{1}}, {0}]", "True" },
        { "MemberQ[{{1}}, {1}, {-2}]", "True" },
        { "MemberQ[{{1}}, 1, {-1}]", "True" },
        { "MemberQ[{{{1}}}, 1, {2}]", "False" },
        { "MemberQ[{{{1}}}, 1, {3}]", "True" },
        { "MemberQ[{{{1}}}, 1, {2, 3}]", "True" },
        { "MemberQ[{{{1}}}, {1}, {2, Infinity}]", "True" },
        { "MemberQ[{{1/2}}, 2, Infinity]", "False" },
        { "MemberQ[{{1/2}}, 1/2, Infinity]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_heads(void) {
    static const Case cs[] = {
        { "MemberQ[f[x], f]", "False" },
        { "MemberQ[f[x], f, Heads -> True]", "True" },
        { "MemberQ[{1}, List, Heads -> True]", "True" },
        { "MemberQ[{1}, List]", "False" },
        { "MemberQ[{1/2}, Rational, Heads -> True]", "True" },
        { "MemberQ[{1/2}, _Rational]", "True" },
        { "MemberQ[{I}, Complex, Heads -> True]", "True" },
        { "MemberQ[{I}, _Complex]", "True" },
        { "MemberQ[{f[g][x]}, g, Infinity]", "False" },
        { "MemberQ[{f[g][x]}, g, Infinity, Heads -> True]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_never_subset_matches_flat_heads(void) {
    static const Case cs[] = {
        { "MemberQ[{a + b + c}, a + b]", "False" },
        { "MemberQ[a + b + c, a + b]", "False" },
        { "MemberQ[{a f[1] f[2]}, x_f y_f]", "False" },
        { "MemberQ[{a f[1] f[2]}, x_f y_f, Infinity]", "False" },
        { "MemberQ[{a f[1] f[2]}, x_f y_f a]", "True" },
        { "MemberQ[{d ** a ** b}, a ** b]", "False" },
        { "MemberQ[{flatf[a, b, c]}, flatf[a, b]]", "False" },
        { "MemberQ[{flatf[a, b, c]}, flatf[a, b, c]]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_holds_and_associations(void) {
    static const Case cs[] = {
        { "MemberQ[Hold[1 + 1], 1 + 1]", "False" },
        { "MemberQ[Hold[1 + 1], HoldPattern[1 + 1]]", "True" },
        { "MemberQ[Hold[x], x]", "True" },
        { "MemberQ[<|a -> 1, b -> 2|>, 2]", "True" },
        { "MemberQ[<|a -> 1, b -> 2|>, a]", "False" },
        { "MemberQ[<|a -> 1, b -> 2|>, _Rule]", "False" },
        { "MemberQ[{a -> 1}, _Rule]", "True" },
        { "MemberQ[{a -> 1}, a, Infinity]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_packed_and_wide(void) {
    static const Case cs[] = {
        { "MemberQ[Range[10], 5]", "True" },
        { "MemberQ[Range[10], 5.]", "False" },
        { "MemberQ[Range[10], 11]", "False" },
        { "MemberQ[Range[10], _Real]", "False" },
        { "MemberQ[N[Range[5]], 3]", "False" },
        { "MemberQ[N[Range[5]], 3.]", "True" },
        { "MemberQ[N[Range[5]], _Real]", "True" },
        { "MemberQ[{{1, 2}, {3, 4}}, {3, 4}]", "True" },
        { "MemberQ[{{1, 2}, {3, 4}}, 4]", "False" },
        { "MemberQ[{{1, 2}, {3, 4}}, 4, 2]", "True" },
        { "MemberQ[Range[100000], 99999]", "True" },
        { "MemberQ[Range[100000], 100001]", "False" },
        { "MemberQ[Table[f[i], {i, 2000}], f[2000]]", "True" },
        { "MemberQ[Table[f[i], {i, 2000}], f[2001]]", "False" },
        { "MemberQ[{Sum[f[i], {i, 600}]}, Sum[f[i], {i, 600}]]", "True" },
        { "MemberQ[{g @@ Range[2000]}, g @@ Range[2000]]", "True" },
        { "MemberQ[{g @@ Range[2000]}, g @@ Table[_Integer, 2000]]", "True" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

static void test_memberq_operator_form(void) {
    static const Case cs[] = {
        { "MemberQ[2][{1, 2}]", "True" },
        { "MemberQ[3][{1, 2}]", "False" },
    };
    run_cases(cs, sizeof cs / sizeof cs[0]);
}

int main() {
    symtab_init();
    core_init();
    TEST(setup_heads);
    TEST(test_freeq_atoms);
    TEST(test_freeq_patterns);
    TEST(test_freeq_levels);
    TEST(test_freeq_heads);
    TEST(test_freeq_numbers);
    TEST(test_freeq_flat_subsets);
    TEST(test_freeq_holds_and_associations);
    TEST(test_freeq_packed_and_wide);
    TEST(test_freeq_operator_form);
    TEST(test_memberq_basic);
    TEST(test_memberq_non_list_heads_and_atoms);
    TEST(test_memberq_levels);
    TEST(test_memberq_heads);
    TEST(test_memberq_never_subset_matches_flat_heads);
    TEST(test_memberq_holds_and_associations);
    TEST(test_memberq_packed_and_wide);
    TEST(test_memberq_operator_form);
    printf("All FreeQ/MemberQ tests passed!\n");
    return 0;
}
