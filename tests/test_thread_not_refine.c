/* test_thread_not_refine.c
 *
 * Algebraic transformations thread over a Rule, relations and logical
 * connectives (FeynCalc's Factor2 runs Together/Numerator/Factor over a list
 * of solution rules, as FCMatchSolve does); Together pulls the integer content
 * out of a denominator-free sum; Not flips a two-argument relation; Refine
 * does not falsify a symbolic equation (Refine[m == a] was False); Reduce over
 * Complexes treats a variable pinned by a linear ordering atom as real.
 * Expected values computed by the reference kernel (InputForm strings).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void eq_if(const char* e, const char* expected) {
    char buf[1024], quoted[1024];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);   /* a String prints quoted */
    assert_eval_eq(buf, quoted, 0);
}

static void test_rule_threading(void) {
    eq_if("Factor[x -> a^2 - 1]", "x -> (-1 + a)*(1 + a)");
    eq_if("Factor[x -> y -> a^2 - 1]", "x -> y -> (-1 + a)*(1 + a)");
    eq_if("Factor[x -> a^2 + 1, Extension -> I]", "x -> (-I + a)*(I + a)");
    eq_if("Factor[Xor[a^2 - 1 == 0, b]]", "Xor[b, (-1 + a)*(1 + a) == 0]");
    eq_if("Expand[x -> (a + b)^2, a]", "x -> a^2 + 2*a*b + b^2");
    eq_if("FactorTerms[x -> 2 a + 4 b x, x]", "x -> 2*(a + 2*b*x)");
    eq_if("Apart[x -> 1/(a (a + 1)), a]", "x -> a^(-1) - (1 + a)^(-1)");
    eq_if("Cancel[m -> (a^2 - 1)/(a - 1)]", "m -> 1 + a");
    eq_if("ExpandNumerator[x -> (a + 1)^2/b]", "x -> (1 + 2*a + a^2)/b");
    assert_eval_eq("Apart[Inequality[0, Less, 1/(a (a + 1)), Less, 1], a] === "
                   "Inequality[0, Less, 1/a - 1/(1 + a), Less, 1]", "True", 0);
    eq_if("Together[{m -> (3 a + 3 a z)/(8 n^2)}]", "{m -> (3*(a + a*z))/(8*n^2)}");
    eq_if("Factor[Together[{m -> (3 a + 3 a z)/(8 n^2)}]]", "{m -> (3*a*(1 + z))/(8*n^2)}");
    /* RuleDelayed is not threaded */
    eq_if("Factor[m :> a^2 - 1]", "m :> a^2 - 1");   /* held as parsed, Plus[a^2, -1] */
}

static void test_together_content(void) {
    eq_if("Together[2 a + 2 b]", "2*(a + b)");
    eq_if("Together[-2 a - 2 b]", "-2*(a + b)");
    eq_if("Together[-a - b]", "-a - b");
    eq_if("Together[2 b - 2 a]", "-2*(a - b)");
    eq_if("Together[2 - 2 x]", "-2*(-1 + x)");
    eq_if("Together[2 x^2 + 4 x + 2]", "2*(1 + 2*x + x^2)");
    eq_if("Together[(a + b) f[c] - a f[c]]", "b*f[c]");
    eq_if("Together[x (1 + x) + 2 x]", "3*x + x^2");
    eq_if("Together[4 a + 6 b + 2 Pi]", "2*(2*a + 3*b + Pi)");
}

static void test_not(void) {
    eq_if("{Not[x == 0], Not[x != 0], Not[x < 1], Not[x >= 1], Not[x > 1], Not[x <= 1]}",
          "{x != 0, x == 0, x >= 1, x < 1, x <= 1, x > 1}");
    assert_eval_eq("Not[Not[a]]", "a", 0);
    assert_eval_eq("Not[Element[x, Reals]]", "NotElement[x, Reals]", 1);
    assert_eval_eq("Not[a && b]", "Not[And[a, b]]", 1);
}

static void test_refine_equations(void) {
    eq_if("{Refine[m == Sqrt[a^2]], Refine[m == a], Refine[x == y + 1], Refine[m == Abs[a]], "
          "Refine[m == Sqrt[a^2], a > 0], Refine[x == -1, x > 0], Refine[x != -1, x > 0], "
          "Refine[2 == 3], Refine[Sqrt[2] == 1], Refine[m != a], Refine[x^2 == -1, x > 0], "
          "Refine[x == y, x > 0 && y < 0], Refine[x == 1, x == 1], Refine[x == 2, x == 1]}",
          "{m == Sqrt[a^2], m == a, x == 1 + y, m == Abs[a], m == a, False, True, False, False, "
          "m != a, False, False, True, False}");
}

static void test_reduce_complexes_pinned_real(void) {
    eq_if("Reduce[x > 0 && y < 0 && x == y, {x, y}, Complexes]", "False");
    eq_if("Reduce[x > 0 && y < 0 && x != y, {x, y}, Complexes]", "x > 0 && y < 0");
    eq_if("Reduce[x > 0 && x^2 == 4, x, Complexes]", "x == 2");
    eq_if("Reduce[Abs[x] < 1 && x^2 == 1/4, x, Complexes]", "x == -1/2 || x == 1/2");
}

static void test_det_rational_entries(void) {
    /* a sum in a denominator: the cofactor form is kept (no blow-up); the
     * expanded polynomial otherwise */
    eq_if("Det[{{a/(x + 1), b}, {c, d/(x + 2)}}]", "-(b*c) + (a*d)/((1 + x)*(2 + x))");
    eq_if("Det[{{a, b, c}, {d, e, f}, {g, h, 1/k}}]",
          "-(c*e*g) + b*f*g + c*d*h - a*f*h - (b*d)/k + (a*e)/k");
}

static void test_literal_symbol_in_pattern(void) {
    /* a literal b is literal even with a pattern variable b bound
     * (FeynCalc's Uncontract-ID48 writes FV[a, b_] FV[q, b_] + FV[b, c_] ...) */
    assert_eval_eq("MatchQ[f[a,i] f[q,i] + f[b,j] f[q,j], f[a, b_] f[q, b_] + f[b, c_] f[q, c_]]", "True", 0);
    assert_eval_eq("{MatchQ[{i, b}, {b_, b}], MatchQ[{i, i}, {b_, b}]}", "List[True, False]", 1);
    assert_eval_eq("Replace[{i, b}, {b_, b} :> {b}]", "List[i]", 1);
}

static void test_freeq_levelspec(void) {
    assert_eval_eq("{FreeQ[{f[0], 2}, 0, {0, 1}], FreeQ[{f[0], 2}, 0, {2}], FreeQ[{f[0], 2}, 0, {-1}], "
                   "FreeQ[{f[g[0]], 2}, g[0], {-2}], FreeQ[{f[g[0]], 2}, f[_], {1, -2}], "
                   "FreeQ[{f[g[0]], 2}, g[_], {1, -3}], FreeQ[{{x}}, x, 1], FreeQ[{{x}}, {x}, {-2}]}",
                   "List[True, False, False, False, False, True, True, False]", 1);
}

static void test_minus_literal_parse(void) {
    assert_eval_eq("Hold[a - 1, a - 2.5, a - 100000000000000000000, a - 2^3, a - 0, 1 - 1]",
                   "Hold[Plus[a, -1], Plus[a, -2.5], Plus[a, -100000000000000000000], "
                   "Plus[a, Times[-1, Power[2, 3]]], Plus[a, 0], Plus[1, -1]]", 1);
    eq_if("Hold[Plus[2, Times[-1, 1]]]", "Hold[2 - 1]");
    assert_eval_eq("{1 - 1, 2^3 - 2, 5 - 2.5}", "List[0, 6, 2.5]", 1);
}

static void test_pzq_pure_function_kept(void) {
    /* zero off the real axis, non-zero on 0 < w < 1: the pure function inside
     * RootSum is not an independent atom */
    assert_eval_eq("PossibleZeroQ[D[w^(1/4) - ArcTan[w^(1/4)]/2 + Log[1 - w^(1/4)]/4 - Log[1 + w^(1/4)]/4 + "
                   "RootSum[1/256 + #1^4 &, Log[1 + 64 w^(1/4) #1^3] #1 &], w] - "
                   "(1 + w^2)/(8 w^2) (1 + (1 + w^2)/(-1 + w^2))^(1/4)/((-1 + w^2)/(2 w))^(3/4), "
                   "Assumptions -> Element[w, Complexes]]", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_rule_threading);
    TEST(test_together_content);
    TEST(test_not);
    TEST(test_refine_equations);
    TEST(test_reduce_complexes_pinned_real);
    TEST(test_det_rational_entries);
    TEST(test_literal_symbol_in_pattern);
    TEST(test_freeq_levelspec);
    TEST(test_minus_literal_parse);
    TEST(test_pzq_pure_function_kept);
    printf("All thread/Not/Refine tests passed.\n");
    return 0;
}
