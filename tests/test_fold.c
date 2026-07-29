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

static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "Fold %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "Fold %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* An NDArray is atomic, so Fold's element walk used to look straight past it
 * and Fold[f, x, NDArray[...]] came back UNEVALUATED while the identical call on
 * a plain List folded fine.  The contract is that the two agree, so each case
 * below asserts the packed form against the List form's own answer.
 *
 * This is also what makes Fold compilable: `Compile[{{v,_Real,1}}, ...]` may be
 * handed a genuine NDArray, and a compiled path that answered where the
 * interpreter declined would be a divergence, not a speedup. */
static void test_fold_ndarray_matches_list(void) {
    run_infix("Fold[Plus, 0., NDArray[{1., 2., 3.}]]", "6.0");
    run_infix("Fold[Times, NDArray[{2., 3., 4.}]]", "24.0");
    /* Rank 2 folds over ROWS, exactly as the nested-List form does. */
    run_infix("Fold[Times, 1., NDArray[{{1., 2.}, {3., 4.}}]]", "{3.0, 8.0}");
    /* Seedless over a singleton is the element itself; no application happens. */
    run_infix("Fold[Plus, NDArray[{7.}]]", "7.0");
    /* A packed argument must not be mutated by the fold. */
    run_infix("ndf = NDArray[{1., 2., 3.}]; Fold[Plus, 0., ndf]; Normal[ndf]",
              "{1.0, 2.0, 3.0}");
    /* FoldList gives the history back packed, matching Map[f, NDArray[...]]. */
    run_infix("FoldList[Plus, 0., NDArray[{1., 2.}]]", "NDArray[{0.0, 1.0, 3.0}]");
    run_infix("Normal[FoldList[Plus, 0., NDArray[{1., 2.}]]] === FoldList[Plus, 0., {1., 2.}]",
              "True");
}

/* ---------- Spec examples ---------- */

static void test_fold_basic_symbolic(void) {
    /* Fold[f, x, {a, b, c, d}] -> f[f[f[f[x, a], b], c], d] */
    run_full("Fold[f, x, {a, b, c, d}]",
             "f[f[f[f[x, a], b], c], d]");
}

static void test_fold_list_head(void) {
    /* Fold[List, x, {a, b, c, d}] -> {{{{x, a}, b}, c}, d} */
    run_full("Fold[List, x, {a, b, c, d}]",
             "List[List[List[List[x, a], b], c], d]");
}

static void test_fold_times_one(void) {
    /* Fold[Times, 1, {a, b, c, d}] -> a b c d */
    run_full("Fold[Times, 1, {a, b, c, d}]",
             "Times[a, b, c, d]");
}

static void test_fold_no_seed(void) {
    /* Fold[f, {a, b, c, d}] -> f[f[f[a, b], c], d] */
    run_full("Fold[f, {a, b, c, d}]",
             "f[f[f[a, b], c], d]");
}

static void test_fold_pure_pair(void) {
    /* Fold[{2 #1, 3 #2}&, x, {a, b, c, d}]
       -> {{{{16 x, 24 a}, 12 b}, 6 c}, 3 d} */
    run_full("Fold[{2 #1, 3 #2} &, x, {a, b, c, d}]",
             "List[List[List[List[Times[16, x], Times[24, a]], Times[12, b]], Times[6, c]], Times[3, d]]");
}

static void test_fold_nonlist_head(void) {
    /* Head need not be List: Fold[f, x, p[a, b, c, d]] -> f[f[f[f[x, a], b], c], d] */
    run_full("Fold[f, x, p[a, b, c, d]]",
             "f[f[f[f[x, a], b], c], d]");
}

static void test_fold_horner_numeric(void) {
    /* Horner evaluation at x=2 of polynomial with coefficients {1,0,1,1}:
       2^3 + 2 + 1 = 11 */
    run_full("Fold[2 #1 + #2 &, 0, {1, 0, 1, 1}]", "11");
}

static void test_fold_digits_to_number(void) {
    /* Fold[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}] -> 451678 */
    run_full("Fold[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]", "451678");
}

static void test_fold_empty_list_with_seed(void) {
    /* Fold[f, x, {}] -> x (function never applied) */
    run_full("Fold[f, x, {}]", "x");
}

static void test_fold_singleton_no_seed(void) {
    /* Fold[f, {a}] -> a (function never applied) */
    run_full("Fold[f, {a}]", "a");
}

static void test_fold_empty_no_seed_unevaluated(void) {
    /* Fold[f, {}] stays unevaluated */
    run_full("Fold[f, {}]", "Fold[f, List[]]");
}

static void test_fold_last_equals_foldlist(void) {
    /* Fold[f, x, {a, b, c}] == Last[FoldList[f, x, {a, b, c}]] */
    run_full("Fold[f, x, {a, b, c}]",
             "f[f[f[x, a], b], c]");
    run_full("Last[FoldList[f, x, {a, b, c}]]",
             "f[f[f[x, a], b], c]");
}

static void test_fold_like_nest(void) {
    /* Fold[f[#1] &, x, Range[5]] -> f[f[f[f[f[x]]]]]
       (second argument ignored by pure function) */
    run_full("Fold[f[#1] &, x, Range[5]]",
             "f[f[f[f[f[x]]]]]");
    /* Should match Nest */
    run_full("Nest[f[#1] &, x, 5]",
             "f[f[f[f[f[x]]]]]");
}

static void test_fold_factorial_via_times(void) {
    /* Fold[Times, 1, Range[5]] -> 120 */
    run_full("Fold[Times, 1, Range[5]]", "120");
}

static void test_fold_plus_sum(void) {
    /* Fold[Plus, 0, Range[10]] -> 55 */
    run_full("Fold[Plus, 0, Range[10]]", "55");
}

static void test_fold_continued_fraction_infix(void) {
    /* Fold[1/(#2+#1)&, x, Reverse[{a,b,c,d}]]
       -> 1/(a+1/(b+1/(c+1/(d+x)))) */
    run_infix("Fold[1/(#2 + #1) &, x, Reverse[{a, b, c, d}]]",
              "1/(a + 1/(b + 1/(c + 1/(d + x))))");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_fold_wrong_arg_count(void) {
    run_full("Fold[f]", "Fold[f]");
    run_full("Fold[]", "Fold[]");
    run_full("Fold[f, x, {a}, extra]", "Fold[f, x, List[a], extra]");
}

static void test_fold_list_is_symbol(void) {
    /* Third arg is a bare symbol, not a compound expression -> unevaluated */
    run_full("Fold[f, x, y]", "Fold[f, x, y]");
}

static void test_fold_attributes_protected(void) {
    run_full("MemberQ[Attributes[Fold], Protected]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_fold_basic_symbolic);
    TEST(test_fold_list_head);
    TEST(test_fold_times_one);
    TEST(test_fold_no_seed);
    TEST(test_fold_pure_pair);
    TEST(test_fold_nonlist_head);
    TEST(test_fold_horner_numeric);
    TEST(test_fold_digits_to_number);
    TEST(test_fold_empty_list_with_seed);
    TEST(test_fold_singleton_no_seed);
    TEST(test_fold_empty_no_seed_unevaluated);
    TEST(test_fold_last_equals_foldlist);
    TEST(test_fold_like_nest);
    TEST(test_fold_factorial_via_times);
    TEST(test_fold_plus_sum);
    TEST(test_fold_continued_fraction_infix);
    TEST(test_fold_ndarray_matches_list);

    /* Unevaluated cases */
    TEST(test_fold_wrong_arg_count);
    TEST(test_fold_list_is_symbol);
    TEST(test_fold_attributes_protected);

    printf("All Fold tests passed!\n");
    return 0;
}
