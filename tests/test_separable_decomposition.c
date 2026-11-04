/* test_separable_decomposition.c
 *
 * SeparablePolynomialDecomposition[p, S, T]
 * (src/internal/simp/SeparablePolynomialDecomposition.m) and its Simplify stage
 * (try_separable_decomposition, src/simp/simp_search.c).
 *
 * The function writes p as Sum_j B_j[S] Q_j[T] with the fewest terms (the rank
 * of the (S, T) coefficient flattening), each factor Factor'd.  Covered here:
 *   - the docstring example and pinned output forms;
 *   - the degenerate inputs (0, constants, monomials, p in S or T alone,
 *     variables absent from p, symbolic / Pi / Gaussian coefficients);
 *   - declines: non-polynomial, empty or overlapping blocks, non-list blocks
 *     answer $Failed; a wrong arity stays unevaluated;
 *   - inexact input: Rationalize[p, 0] then decompose (an exact answer), or
 *     p itself unchanged when inexact numbers survive Rationalize;
 *   - properties over seeded random inputs: the answer is identically p, and
 *     its rank never exceeds the rank p was built with;
 *   - minimal rank against an independent MatrixRank of the flattening;
 *   - Simplify: the stage fires after Factor (strictly cheaper sums only, so
 *     a Factor tie, and any rank-1 answer, stays Factor's), never changes the value, never changes the
 *     precision of an inexact input, leaves univariate input alone, and
 *     works under assumptions;
 *   - symbol hygiene: the package's locals do not capture user symbols named
 *     like them (r, cs, c, M, Q, ...).
 * Timing gates live in bench_separable_decomposition.c (nightly tier).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* `e` evaluates to True. */
static void check(const char* e) { assert_eval_eq(e, "True", 0); }

/* lhs and rhs evaluate to the same expression. */
static void same(const char* lhs, const char* rhs) {
    char buf[8192];
    snprintf(buf, sizeof buf, "(%s) === (%s)", lhs, rhs);
    check(buf);
}

/* The InputForm of `e` is exactly `expected`. */
static void inputform(const char* e, const char* expected) {
    char buf[8192], quoted[8192];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);
    assert_eval_eq(buf, quoted, 0);
}

/* SeparablePolynomialDecomposition[p, S, T] is identically p. */
static void exact(const char* p, const char* S, const char* T) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "With[{pp = %s}, Expand[SeparablePolynomialDecomposition[pp, %s, %s] - pp] === 0]",
             p, S, T);
    check(buf);
}

/* Shared definitions for the tests below. */
static void setup(void) {
    assert_eval_eq(
        "spdP0 = Expand[(3 a b - b^2)^3 (x^2 - y + 1)^2 - (a - b)^2 (x - y)]; Null",
        "Null", 0);
    /* Random polynomial: nt terms in vars, each exponent in [0, deg]. */
    assert_eval_eq(
        "spdRP[vars_, deg_, nt_] := Total@Table[RandomInteger[{-9, 9}] "
        "Times @@ (vars^RandomInteger[{0, deg}, Length[vars]]), {k, nt}]; Null",
        "Null", 0);
    /* Rank of the (S, T) flattening, computed independently of the package:
     * rows = T-exponent vectors, columns = S-exponent vectors. */
    assert_eval_eq(
        "spdRank[p_, S_, T_] := Module[{cr = CoefficientRules[p, Join[S, T]], "
        "n = Length[S], rows, cols}, "
        "rows = Union[Drop[First[#], n] & /@ cr]; cols = Union[Take[First[#], n] & /@ cr]; "
        "MatrixRank[Table[Total[Cases[cr, (e_ -> v_) /; Drop[e, n] === rr && "
        "Take[e, n] === cc :> v]], {rr, rows}, {cc, cols}]]]; Null",
        "Null", 0);
    /* Number of separable products in an answer: B Q terms, counting a term
     * whose T-part is a constant once (it is a polynomial in S alone, and
     * Plus flattens it into the top level). */
    assert_eval_eq(
        "spdTerms[e_, T_] := Module[{ts = If[Head[e] === Plus, List @@ e, {e}]}, "
        "Length[Select[ts, !FreeQ[#, Alternatives @@ T] &]] + "
        "If[AnyTrue[ts, FreeQ[#, Alternatives @@ T] &], 1, 0]]; Null",
        "Null", 0);
}

static void test_docstring_example(void) {
    inputform("SeparablePolynomialDecomposition[spdP0, {a, b}, {x, y}]",
              "-((a - b)^2*(x - y)) + (3*a - b)^3*b^3*(1 + x^2 - y)^2");
    exact("spdP0", "{a, b}", "{x, y}");
    check("LeafCount[SeparablePolynomialDecomposition[spdP0, {a, b}, {x, y}]] == 38");
    check("LeafCount[spdP0] == 251");
    /* the swapped split finds the same rank-2 form */
    exact("spdP0", "{x, y}", "{a, b}");
    check("spdTerms[SeparablePolynomialDecomposition[spdP0, {x, y}, {a, b}], {a, b}] == 2");
}

static void test_pinned_forms(void) {
    inputform("SeparablePolynomialDecomposition[Expand[(a + b)^2 (x - y)^3], {a, b}, {x, y}]",
              "(a + b)^2*(x - y)^3");
    inputform("SeparablePolynomialDecomposition[Expand[(a + b)^2 (x - y)^3], {x, y}, {a, b}]",
              "(a + b)^2*(x - y)^3");
    inputform("SeparablePolynomialDecomposition[Expand[(a/2 + 1/3) (x - 3/4)], {a}, {x}]",
              "((-2 - 3*a)*(3 - 4*x))/24");
    /* rank 2 with basis {1, x} over T */
    inputform("SeparablePolynomialDecomposition[Expand[(a + 2) (x + 1) + (a - 3) (x - 1)], {a}, {x}]",
              "5 + (-1 + 2*a)*x");
    inputform("SeparablePolynomialDecomposition[Expand[(a + I) (x + 2)], {a}, {x}]",
              "(I + a)*(2 + x)");
    inputform("SeparablePolynomialDecomposition[Expand[Pi (a + 1) (x + 2)], {a}, {x}]",
              "(1 + a)*Pi*(2 + x)");
    /* a symbol outside S and T is a constant coefficient */
    inputform("SeparablePolynomialDecomposition[Expand[k (a + 1) (x + 2)], {a}, {x}]",
              "(1 + a)*k*(2 + x)");
}

static void test_degenerate_inputs(void) {
    same("SeparablePolynomialDecomposition[0, {a}, {x}]", "0");
    same("SeparablePolynomialDecomposition[7, {a}, {x}]", "7");
    same("SeparablePolynomialDecomposition[-3/4, {a, b}, {x, y}]", "-3/4");
    same("SeparablePolynomialDecomposition[3 a^2 x, {a}, {x}]", "3 a^2 x");
    same("SeparablePolynomialDecomposition[a^2 + 2 a + 1, {a}, {x}]", "(1 + a)^2");
    same("SeparablePolynomialDecomposition[x^2 - 1, {a}, {x}]", "(-1 + x) (1 + x)");
    /* a block variable absent from p */
    same("SeparablePolynomialDecomposition[Expand[(a + 1) (x + 1)], {a, zz}, {x}]",
         "(1 + a) (1 + x)");
    same("SeparablePolynomialDecomposition[Expand[(a + 1) (x + 1)], {a}, {x, zz}]",
         "(1 + a) (1 + x)");
    /* p independent of every block variable: a constant */
    same("SeparablePolynomialDecomposition[k^2 + 1, {a}, {x}]", "1 + k^2");
    /* an unexpanded input is accepted (CoefficientRules expands) */
    exact("(a + b)^3 (x + y)^2 + (a - b) (x - y)^2", "{a, b}", "{x, y}");
}

static void test_declines(void) {
    same("SeparablePolynomialDecomposition[Sqrt[a] x, {a}, {x}]", "$Failed");
    same("SeparablePolynomialDecomposition[Sin[x] a, {a}, {x}]", "$Failed");
    same("SeparablePolynomialDecomposition[1/(a + x), {a}, {x}]", "$Failed");
    same("SeparablePolynomialDecomposition[a^(1/2) + x, {a}, {x}]", "$Failed");
    same("SeparablePolynomialDecomposition[a x, {}, {a, x}]", "$Failed");
    same("SeparablePolynomialDecomposition[a x, {a, x}, {}]", "$Failed");
    same("SeparablePolynomialDecomposition[a x, {a, x}, {x}]", "$Failed");   /* overlap */
    same("SeparablePolynomialDecomposition[a x, a, x]", "$Failed");         /* not lists */
    /* wrong arity stays unevaluated */
    inputform("SeparablePolynomialDecomposition[a x, {a}]",
              "SeparablePolynomialDecomposition[a*x, {a}]");
    inputform("SeparablePolynomialDecomposition[]", "SeparablePolynomialDecomposition[]");
}

static void test_inexact_inputs(void) {
    /* machine reals: rationalised, decomposed, exact answer */
    inputform("SeparablePolynomialDecomposition[Expand[(a + 1.5) (x + 2)], {a}, {x}]",
              "((3 + 2*a)*(2 + x))/2");
    inputform("SeparablePolynomialDecomposition[Expand[(a + 0.25) (x - 0.5)], {a}, {x}]",
              "((-1 - 4*a)*(1 - 2*x))/8");
    /* arbitrary precision */
    inputform("SeparablePolynomialDecomposition[Expand[(a + N[1/3, 30]) (x + 2)], {a}, {x}]",
              "((1 + 3*a)*(2 + x))/3");
    /* inexact complex */
    inputform("SeparablePolynomialDecomposition[Expand[(a + (2.0 + 0.5 I)) (x + 2)], {a}, {x}]",
              "((4 + I + 2*a)*(2 + x))/2");
    /* the decomposition of the rationalised form is p exactly */
    check("With[{e = Expand[(a + 0.25 b) (x - 0.5 y) + (a - b)^2 (x + y)]}, "
          "Expand[SeparablePolynomialDecomposition[e, {a, b}, {x, y}] - Rationalize[e, 0]] === 0]");
    check("FreeQ[SeparablePolynomialDecomposition[Expand[(a + 0.25 b) (x - 0.5 y)], {a, b}, {x, y}], "
          "_?InexactNumberQ]");
    /* inexact numbers that survive Rationalize[p, 0]: p comes back unchanged */
    check("With[{e = 1.*^300 + 1.*^-300 a x}, SeparablePolynomialDecomposition[e, {a}, {x}] === e]");
    /* rationalised but not a polynomial: the ordinary decline */
    same("SeparablePolynomialDecomposition[1.5 Sqrt[a] x, {a}, {x}]", "$Failed");
    same("SeparablePolynomialDecomposition[0., {a}, {x}]", "0");
}

static void test_random_exactness_and_rank(void) {
    /* 1-3 variables per block, rank 1-4, 25 seeded draws per shape: the answer
     * is identically p and has at most the rank p was built with. */
    check("SeedRandom[42]; And @@ Flatten@Table["
          "Module[{SS = {a, b, c}[[;; ns]], TT = {x, y, z}[[;; nt]], r, p, res}, "
          "r = RandomInteger[{1, 4}]; "
          "p = Expand[Total@Table[spdRP[SS, 3, 3] spdRP[TT, 3, 3], {j, r}]]; "
          "res = SeparablePolynomialDecomposition[p, SS, TT]; "
          "Expand[res - p] === 0 && spdTerms[res, TT] <= r], "
          "{ns, 1, 3}, {nt, 1, 3}, {rep, 25}]");
    /* larger: degree 6, rank 8 */
    check("SeedRandom[7]; And @@ Table["
          "Module[{p = Expand[Total@Table[spdRP[{a, b}, 6, 6] spdRP[{x, y}, 6, 6], {j, 8}]], res}, "
          "res = SeparablePolynomialDecomposition[p, {a, b}, {x, y}]; "
          "Expand[res - p] === 0 && spdTerms[res, {x, y}] <= 8], {rep, 10}]");
}

static void test_minimal_rank(void) {
    /* the number of separable terms equals the rank of the flattening */
    check("SeedRandom[3]; And @@ Table["
          "Module[{p = Expand[Total@Table[spdRP[{a, b}, 3, 4] spdRP[{x, y}, 3, 4], {j, r}]], res}, "
          "res = SeparablePolynomialDecomposition[p, {a, b}, {x, y}]; "
          "spdTerms[res, {x, y}] == spdRank[p, {a, b}, {x, y}]], {r, 1, 5}, {rep, 4}] // Flatten // "
          "Apply[And]");
    check("spdRank[spdP0, {a, b}, {x, y}] == 2");
    /* full rank: no compression is possible, still exact */
    check("With[{p = Expand[(1 + a + b + x + y)^3]}, "
          "spdRank[p, {a, b}, {x, y}] == 4 && "
          "Expand[SeparablePolynomialDecomposition[p, {a, b}, {x, y}] - p] === 0]");
    /* x y + 1 is rank 2 and already minimal */
    same("SeparablePolynomialDecomposition[x y + 1, {x}, {y}]", "1 + x y");
}

static void test_simplify_stage(void) {
    /* the docstring example: Expand / Factor give 251 leaves, the stage 38 */
    check("With[{s = Simplify[spdP0]}, Expand[s - spdP0] === 0 && LeafCount[s] <= 38]");
    inputform("Simplify[spdP0]", "-((a - b)^2*(x - y)) + (3*a - b)^3*b^3*(1 + x^2 - y)^2");
    check("With[{e = Expand[(1 + a + b)^6 (1 + x + y)^6 - (a - b)^3 (x + y)^4]}, "
          "With[{s = Simplify[e]}, Expand[s - e] === 0 && LeafCount[s] < LeafCount[e]/10]]");
    check("With[{e = Expand[(a + b)^2 (x + y)^2 + (a - b) (x - y)]}, "
          "With[{s = Simplify[e]}, Expand[s - e] === 0 && LeafCount[s] < 25]]");
    /* Factor still wins its own cases (ties stay with Factor) */
    same("Simplify[Expand[(a + b) (x + y)]]", "(a + b) (x + y)");
    same("Simplify[Expand[(x - 1)^3 (y + 2)^2]]", "(-1 + x)^3 (2 + y)^2");
    same("Simplify[Expand[(1 + a + b + c + x + y)^4]]", "(1 + a + b + c + x + y)^4");
    /* a rank-1 decomposition is a factorisation: never offered, so Factor's
     * sign convention stands (2 (4 + 3 a) (1 - 3 x) scores 13 < 14) */
    same("Simplify[(-2 + 6 x) (-4 - 3 a)]", "-2 (4 + 3 a) (-1 + 3 x)");
    /* the expanded input already answers this through Factor (unchanged by
     * the stage) */
    same("Simplify[Expand[(-2 + 6 x) (-4 - 3 a)]]", "2 (4 + 3 a) (1 - 3 x)");
    /* below the gate: too few terms, one variable */
    same("Simplify[x y + x + y + 1]", "1 + x + y + x y");
    same("Simplify[a x + b y]", "a x + b y");
    same("Simplify[Expand[(1 + x)^4 + (1 + x)^2]]", "(1 + x)^2 (2 + 2 x + x^2)");
    /* inexact input keeps its precision (Simplify rationalises internally and
     * numericises its answer) and its value */
    check("With[{e = Expand[(a + 1.5) (x + 2) + (a - b)^2 (x - y)]}, "
          "With[{s = Simplify[e]}, Precision[s] === MachinePrecision && "
          "Chop[Expand[s - e]] === 0 && LeafCount[s] < LeafCount[e]]]");
    /* never changes the value: random rank-2 sums in 4 variables */
    check("SeedRandom[9]; And @@ Table["
          "Module[{e = Expand[spdRP[{a, b}, 2, 3] spdRP[{x, y}, 2, 3] + "
          "spdRP[{a, b}, 2, 3] spdRP[{x, y}, 2, 3]], s}, "
          "s = Simplify[e]; Expand[s - e] === 0 && LeafCount[s] <= LeafCount[e]], {rep, 20}]");
    /* a polynomial under assumptions.  (Not Simplify[Sin[spdP0]]: that runs
     * > 60 s with or without this stage, a separate Simplify cost.) */
    check("With[{s = Simplify[spdP0, a > 0]}, Expand[s - spdP0] === 0 && LeafCount[s] <= 38]");
}

static void test_symbol_hygiene(void) {
    /* user symbols named like the package's locals are not captured */
    check("Block[{r = 5, cs = 6, M = 7, Q = 8, B = 9, R = 10, piv = 11, mons = 12}, "
          "Expand[SeparablePolynomialDecomposition[spdP0, {a, b}, {x, y}] - spdP0] === 0]");
    /* block variables named like the package's arguments and locals */
    exact("Expand[(p + S)^2 (T + c)^3 + (p - S) (T - c)]", "{p, S}", "{T, c}");
    exact("Expand[(r + M)^2 (Q + B)^3 + (r - M) (Q - B)]", "{r, M}", "{Q, B}");
}

static void test_symbol_metadata(void) {
    same("Attributes[SeparablePolynomialDecomposition]", "{Protected}");
    check("StringQ[SeparablePolynomialDecomposition::usage] && "
          "StringContainsQ[SeparablePolynomialDecomposition::usage, \"B[S] Q[T]\"]");
    same("Context[SeparablePolynomialDecomposition]", "\"SeparablePolynomialDecomposition`\"");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    TEST(test_docstring_example);
    TEST(test_pinned_forms);
    TEST(test_degenerate_inputs);
    TEST(test_declines);
    TEST(test_inexact_inputs);
    TEST(test_random_exactness_and_rank);
    TEST(test_minimal_rank);
    TEST(test_simplify_stage);
    TEST(test_symbol_hygiene);
    TEST(test_symbol_metadata);
    printf("All SeparablePolynomialDecomposition tests passed.\n");
    return 0;
}
