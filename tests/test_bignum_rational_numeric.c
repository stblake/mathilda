/* Mathilda — regression tests for numeric builtins on BIGNUM RATIONALS.
 *
 * A rational whose numerator or denominator overflows int64 is stored as a
 * Rational[X, Y] with a BigInt component (e.g. Rational[1, 10^25]).  A whole
 * family of numeric heads recognised rationals only through is_rational(),
 * whose int64 out-parameters cannot represent such a value, so they silently
 * declined — leaving Abs / Sign / Re / Im / ReIm / Arg / Numerator /
 * Denominator unevaluated and Max / Min / Sort / Median misordering (they
 * read a bignum rational as 0).  The shared helpers expr_numeric_sign,
 * expr_compare, is_real_numeric and extract_num_den were the root causes.
 *
 * These tests pin the fixed behaviour.  They compare the STANDARD printed
 * form (not FullForm) so the expected strings read naturally, and use
 * ASSERT_MSG (test_utils.h) which aborts unconditionally — unlike the libc
 * assert() in assert_eval_eq, which -DNDEBUG turns into a no-op under a
 * CMake Release build.
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

static void chk(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e); expr_free(res);
}

/* A small negative rational (~ -1.24e-15): int64 numerator, bignum
 * denominator — the exact value from the original bug report. */
#define RS "-29043852094387/23409578234095283745029348750"
/* Both components bignum, does not reduce (~0.333). */
#define RR "100000000000000000000000000000/300000000000000000000000000001"
/* Bignum numerator, int64 denominator (~8.06e14). */
#define RP "23409578234095283745029348750/29043852094387"

static void test_abs(void) {
    chk("Abs[" RS "]", "29043852094387/23409578234095283745029348750");
    chk("Abs[-(" RR ")]", RR);
    chk("Abs[" RP "]", RP);
    chk("Abs[-(" RP ")]", RP);
    /* int64 rationals unchanged */
    chk("Abs[-3/4]", "3/4");
}

static void test_sign(void) {
    chk("Sign[" RS "]", "-1");
    chk("Sign[" RR "]", "1");
    chk("Sign[-(" RR ")]", "-1");
    chk("Sign[" RP "]", "1");
    /* underflow guard: |value| far below double resolution still signs right */
    chk("Sign[-1/10^400]", "-1");
    chk("Sign[-3/4]", "-1");
}

static void test_re_im_reim(void) {
    chk("Re[" RS "]", RS);
    chk("Im[" RS "]", "0");
    chk("ReIm[" RS "]", "{" RS ", 0}");
    /* plain bigint integer, which was also declining */
    chk("Re[10^30]", "1000000000000000000000000000000");
    chk("Im[10^30]", "0");
    chk("ReIm[10^30]", "{1000000000000000000000000000000, 0}");
    /* int64 path unchanged */
    chk("Re[5/7]", "5/7");
    chk("Im[5/7]", "0");
}

static void test_arg(void) {
    chk("Arg[" RS "]", "Pi");   /* negative */
    chk("Arg[" RR "]", "0");    /* positive */
    chk("Arg[10^30]", "0");     /* positive bigint */
    chk("Arg[-10^30]", "Pi");   /* negative bigint */
    chk("Arg[-1/10^400]", "Pi");/* tiny negative rational */
    /* int64 / real paths unchanged */
    chk("Arg[-3]", "Pi");
    chk("Arg[3]", "0");
    chk("Arg[0]", "0");
}

static void test_conjugate(void) {
    chk("Conjugate[" RS "]", RS);
    chk("Conjugate[10^30]", "1000000000000000000000000000000");
}

static void test_numerator_denominator(void) {
    chk("Numerator[" RS "]", "-29043852094387");
    chk("Denominator[" RS "]", "23409578234095283745029348750");
    chk("Numerator[" RR "]", "100000000000000000000000000000");
    chk("Denominator[" RR "]", "300000000000000000000000000001");
    /* int64 path unchanged */
    chk("Numerator[6/8]", "3");
    chk("Denominator[6/8]", "4");
}

static void test_min_max(void) {
    chk("Max[" RS ", " RR ", " RP ", 0]", RP);
    chk("Min[" RS ", " RR ", " RP ", 0]", RS);
    /* int64 path unchanged */
    chk("Max[1/2, 3/4, 1/3]", "3/4");
    chk("Min[1/2, 3/4, 1/3]", "1/3");
}

static void test_sort_median_ordering(void) {
    /* Canonical order: -10^30 < RS(~ -1.24e-15) < 0 < RR(~0.333) < RP(~8e14) < 10^30 */
    chk("Sort[{" RP ", " RR ", " RS ", 0, 10^30, -10^30}]",
        "{-1000000000000000000000000000000, " RS ", 0, " RR ", " RP
        ", 1000000000000000000000000000000}");
    /* Median of {RP, RS, 0} is the middle value 0 (was RP before the fix). */
    chk("Median[{" RP ", " RS ", 0}]", "0");
    /* Ordering returns index permutation that sorts ascending. */
    chk("Ordering[{" RP ", " RR ", " RS "}]", "{3, 2, 1}");
    /* int64 path unchanged */
    chk("Sort[{3/2, -1, 0, 5, 1/3}]", "{-1, 0, 1/3, 3/2, 5}");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_abs);
    TEST(test_sign);
    TEST(test_re_im_reim);
    TEST(test_arg);
    TEST(test_conjugate);
    TEST(test_numerator_denominator);
    TEST(test_min_max);
    TEST(test_sort_median_ordering);

    printf("All bignum-rational numeric tests passed!\n");
    return 0;
}
