/*
 * test_limit_assumptions.c -- Targeted Assumptions support for Limit.
 *
 * Bugs 6 & 7 -- Limit[x^n, n -> ±Infinity, Assumptions -> Abs[x] R c]
 *
 *   The standard Limit machinery has no assumption context, so picocas
 *   used to return E^DirectedInfinity[Log[x]] (an unevaluated form) for
 *   Limit[x^n, n -> Infinity] regardless of any Abs[x] constraint.
 *
 *   builtin_limit_impl now intercepts the specific (but useful) shape
 *      Limit[Power[base, lim_var], lim_var -> ±Infinity,
 *            Assumptions -> Abs[base] R c]
 *   and dispatches based on (R, c vs 1):
 *      Abs[base] < 1   -> 0   (x -> ComplexInfinity for -Infinity)
 *      Abs[base] > 1   -> ComplexInfinity   (-> 0 for -Infinity)
 *      Abs[base] == 1  -> Indeterminate
 *   Boundary cases (Abs[base] <=/>= 1) require the strict inequality
 *   on c to reduce; otherwise the caller falls through to the standard
 *   (assumption-free) limit machinery.
 *
 *   Strict scope: only Limit[Power[B, var], var -> ±Inf, Assumptions ->
 *   Abs[expr] R c] with `expr` matching `B` structurally. Everything
 *   else passes through unchanged.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_bug6_abs_lt_one(void) {
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> Abs[x] < 1]",
                   "0", 0);
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> Abs[x] < 1/2]",
                   "0", 0);
    /* Swapped inequality (c on left). */
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> 1 > Abs[x]]",
                   "0", 0);
}

static void test_bug7_abs_gt_one(void) {
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> Abs[x] > 1]",
                   "ComplexInfinity", 0);
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> Abs[x] > 2]",
                   "ComplexInfinity", 0);
    /* Swapped inequality. */
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> 1 < Abs[x]]",
                   "ComplexInfinity", 0);
}

static void test_neg_infinity_inversion(void) {
    /* For var -> -Infinity, x^var has the inverse behaviour:
     * |x| < 1 blows up, |x| > 1 vanishes. */
    assert_eval_eq("Limit[x^n, n -> -Infinity, Assumptions -> Abs[x] < 1]",
                   "ComplexInfinity", 0);
    assert_eval_eq("Limit[x^n, n -> -Infinity, Assumptions -> Abs[x] > 1]",
                   "0", 0);
}

static void test_boundary_abs_eq_one(void) {
    /* Abs[x] == 1 -> Indeterminate (limit value depends on x's argument). */
    assert_eval_eq("Limit[x^n, n -> Infinity, Assumptions -> Abs[x] == 1]",
                   "Indeterminate", 0);
}

static void test_passthrough_no_assumption(void) {
    /* Without an Abs constraint, the standard machinery still produces
     * the unevaluated DirectedInfinity-of-Log form. Pin to confirm
     * the assumption-dispatch only fires when applicable. */
    assert_eval_eq("Limit[x^n, n -> Infinity]",
                   "E^DirectedInfinity[Log[x]]", 0);
    /* Concrete |x|<1 base still works via the standard machinery. */
    assert_eval_eq("Limit[(1/2)^n, n -> Infinity]", "0", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bug6_abs_lt_one);
    TEST(test_bug7_abs_gt_one);
    TEST(test_neg_infinity_inversion);
    TEST(test_boundary_abs_eq_one);
    TEST(test_passthrough_no_assumption);

    printf("All limit_assumptions tests passed!\n");
    return 0;
}
