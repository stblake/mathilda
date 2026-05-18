/* Tests for Chop.
 *
 * Cover: machine reals at, below, and above the default 10^-10 tolerance;
 * positive and negative values; threading through lists, Plus, Times, and
 * arbitrary heads; Complex handling for both "machine complex" (both
 * components Real) and mixed Integer/Real complex shapes; the explicit
 * delta argument, including Integer and Rational deltas; pass-through for
 * exact numbers and symbolic input; structural pass-through for nested
 * subexpressions; the Protected attribute; and basic memory hygiene
 * (constructing many Chop expressions in a tight loop).
 *
 * The headline cases mirror the spec examples in the Chop docstring
 * (see docs/spec/changelog/2026-05.md).
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  Default tolerance (10^-10)
 * ---------------------------------------------------------------------- */

static void test_chop_below_default_tolerance(void) {
    /* Any Real below 10^-10 chops to the exact integer 0. */
    assert_eval_eq("Chop[1.0e-12]", "0", 0);
    assert_eval_eq("Chop[-1.0e-12]", "0", 0);
    assert_eval_eq("Chop[1.0e-15]", "0", 0);
    assert_eval_eq("Chop[1.22461e-16]", "0", 0);
}

static void test_chop_above_default_tolerance(void) {
    /* Reals above tolerance survive unchanged. */
    assert_eval_eq("Chop[1.5]", "1.5", 0);
    assert_eval_eq("Chop[-3.5]", "-3.5", 0);
    /* 10^-9 is just above the 10^-10 default tolerance. */
    assert_eval_eq("Chop[1.0e-9]", "1e-09", 0);
}

static void test_chop_at_boundary(void) {
    /* fabs(x) < delta is strict: x == delta does NOT chop. */
    assert_eval_eq("Chop[1.0e-10]", "1e-10", 0);
}

static void test_chop_zero(void) {
    /* 0.0 chops to exact integer 0 (|0| < delta). */
    assert_eval_eq("Chop[0.0]", "0", 0);
    /* Integer 0 stays as integer 0 (Chop is a no-op on exact numbers). */
    assert_eval_eq("Chop[0]", "0", 0);
}

/* ------------------------------------------------------------------------
 *  Explicit delta argument
 * ---------------------------------------------------------------------- */

static void test_chop_with_explicit_real_delta(void) {
    /* Larger delta widens the chop window. */
    assert_eval_eq("Chop[1.0e-5, 1.0e-4]", "0", 0);
    /* 1.5 is above 0.1 -> survives. */
    assert_eval_eq("Chop[1.5, 0.1]", "1.5", 0);
    /* 0.05 is below 0.1 -> chops. */
    assert_eval_eq("Chop[0.05, 0.1]", "0", 0);
}

static void test_chop_with_integer_delta(void) {
    /* Integer delta is accepted. Anything with |x| < 1 chops. */
    assert_eval_eq("Chop[0.5, 1]", "0", 0);
    assert_eval_eq("Chop[2.5, 1]", "2.5", 0);
}

static void test_chop_with_rational_delta(void) {
    /* Rational delta is accepted. */
    assert_eval_eq("Chop[0.4, 1/2]", "0", 0);
    assert_eval_eq("Chop[0.6, 1/2]", "0.6", 0);
}

static void test_chop_with_negative_delta(void) {
    /* Magnitude of delta is what matters; -1e-5 is treated as 1e-5. */
    assert_eval_eq("Chop[1.0e-7, -1.0e-5]", "0", 0);
}

/* ------------------------------------------------------------------------
 *  Sign
 * ---------------------------------------------------------------------- */

static void test_chop_negative_real(void) {
    assert_eval_eq("Chop[-1.0e-12]", "0", 0);
    assert_eval_eq("Chop[-1.5]", "-1.5", 0);
}

/* ------------------------------------------------------------------------
 *  Complex handling -- the trickiest part
 * ---------------------------------------------------------------------- */

static void test_chop_drops_small_imaginary_part(void) {
    /* Spec example: small imaginary part on a machine complex dissolves
     * the Complex wrapper, leaving just the machine real. */
    assert_eval_eq("Chop[2.0 + 1.0e-12 I]", "2.0", 0);
    assert_eval_eq("Chop[-3.5 + 1.0e-15 I]", "-3.5", 0);
}

static void test_chop_keeps_machine_zero_real_part(void) {
    /* Spec example: small real part on a machine complex becomes 0.0
     * (a machine zero), NOT exact 0, so the result remains a machine
     * complex. Prints as "0.0 + 2.0*I" in Mathilda's standard form. */
    assert_eval_eq("Chop[1.0e-12 + 2.0 I]", "0.0 + 2.0*I", 0);
}

static void test_chop_both_parts_small(void) {
    /* Both components below tolerance -> exact integer 0. */
    assert_eval_eq("Chop[1.0e-12 + 1.0e-12 I]", "0", 0);
}

static void test_chop_both_parts_above(void) {
    /* Neither component chops -> Complex unchanged. */
    assert_eval_eq("Chop[2.0 + 3.0 I]", "2.0 + 3.0*I", 0);
}

static void test_chop_complex_with_exact_real_part(void) {
    /* Complex[1, 1.e-12]: integer real part, tiny real imaginary part.
     * Not "machine complex" (re is Integer), so general recursion kicks
     * in: chop(1)=1, chop(1.e-12)=0 (exact integer), and builtin_complex
     * collapses Complex[1, 0] -> 1. */
    assert_eval_eq("Chop[Complex[1, 1.0e-12]]", "1", 0);
}

static void test_chop_complex_with_exact_imag_part(void) {
    /* Complex[1.e-12, 1]: tiny real real-part, integer imaginary part.
     * Not machine complex (im is Integer). chop(1.e-12)=0, chop(1)=1,
     * leaving Complex[0, 1] which is the canonical I. */
    assert_eval_eq("Chop[Complex[1.0e-12, 1]]", "I", 0);
}

static void test_chop_exact_complex_passthrough(void) {
    /* Pure exact-number Complex: nothing approximate to chop. */
    assert_eval_eq("Chop[Complex[1, 2]]", "1 + 2*I", 0);
    assert_eval_eq("Chop[Complex[1/2, 1/3]]", "1/2 + 1/3*I", 0);
}

/* ------------------------------------------------------------------------
 *  Threading through structure
 * ---------------------------------------------------------------------- */

static void test_chop_inside_list(void) {
    /* The spec's flagship list-chop example. */
    assert_eval_eq(
        "Chop[{-1.0 + 1.22461e-16 I, 1.0 - 2.44921e-16 I, "
        "-1.0 + 3.67382e-16 I, 1.0 - 4.89843e-16 I}]",
        "{-1.0, 1.0, -1.0, 1.0}", 0);
}

static void test_chop_inside_nested_list(void) {
    assert_eval_eq("Chop[{{1.0e-12, 1.5}, {2.0, 1.0e-15}}]",
                   "{{0, 1.5}, {2.0, 0}}", 0);
}

static void test_chop_inside_plus(void) {
    /* Plus[x, 1.e-12] chops to Plus[x, 0] which Plus auto-simplifies. */
    assert_eval_eq("Chop[x + 1.0e-12]", "x", 0);
}

static void test_chop_inside_times(void) {
    /* Times[x, 1.e-12] chops to Times[x, 0] -> 0. */
    assert_eval_eq("Chop[x * 1.0e-12]", "0", 0);
}

static void test_chop_inside_generic_head(void) {
    /* Arbitrary head f[...]: Chop just walks the args. */
    assert_eval_eq("Chop[f[1.0e-12, 1.5, x]]", "f[0, 1.5, x]", 0);
}

static void test_chop_listable_threading(void) {
    /* With a delta argument, the 2-arg form also walks the structure
     * via manual recursion, applying the same tolerance to every element. */
    assert_eval_eq("Chop[{1.0e-5, 1.0e-3, 1.5}, 1.0e-4]",
                   "{0, 0.001, 1.5}", 0);
}

/* ------------------------------------------------------------------------
 *  Exact-number pass-through
 * ---------------------------------------------------------------------- */

static void test_chop_exact_integers(void) {
    assert_eval_eq("Chop[5]", "5", 0);
    assert_eval_eq("Chop[-7]", "-7", 0);
    /* Integer 1 with a giant delta survives because Chop only touches
     * APPROXIMATE numbers. */
    assert_eval_eq("Chop[1, 100]", "1", 0);
}

static void test_chop_exact_rationals(void) {
    assert_eval_eq("Chop[1/2]", "1/2", 0);
    /* Even with a huge delta, exact rationals don't chop. */
    assert_eval_eq("Chop[1/1000000, 1]", "1/1000000", 0);
}

static void test_chop_bigint_passthrough(void) {
    /* A BigInt larger than int64 max passes through unchanged. */
    assert_eval_eq("Chop[10^30]",
                   "1000000000000000000000000000000", 0);
}

static void test_chop_symbolic_passthrough(void) {
    assert_eval_eq("Chop[x]", "x", 0);
    assert_eval_eq("Chop[Sin[x]]", "Sin[x]", 0);
    assert_eval_eq("Chop[Pi]", "Pi", 0);
}

/* ------------------------------------------------------------------------
 *  Rationalize-style equivalence check from the spec discussion.
 * ---------------------------------------------------------------------- */

static void test_chop_rationalize_difference_default(void) {
    /* x = N[Pi], r = Rationalize[x, 10^-12]; x - r is below 10^-10. */
    assert_eval_eq(
        "Chop[N[Pi] - Rationalize[N[Pi], 10^-12]] === 0",
        "True", 0);
}

static void test_chop_rationalize_difference_tight(void) {
    /* Tighten tolerance to 10^-14: the difference no longer chops. */
    assert_eval_eq(
        "Chop[N[Pi] - Rationalize[N[Pi], 10^-12], 10^-14] === 0",
        "False", 0);
}

/* ------------------------------------------------------------------------
 *  Argument-count / argument-shape edge cases
 * ---------------------------------------------------------------------- */

static void test_chop_no_args_unevaluated(void) {
    assert_eval_eq("Chop[]", "Chop[]", 0);
}

static void test_chop_too_many_args_unevaluated(void) {
    /* 3+ args: stays unevaluated. */
    assert_eval_eq("Chop[1.0, 1.0, 1.0]", "Chop[1.0, 1.0, 1.0]", 0);
}

static void test_chop_symbolic_delta_unevaluated(void) {
    /* Non-numeric delta: builtin can't coerce, stays unevaluated. */
    assert_eval_eq("Chop[1.5, x]", "Chop[1.5, x]", 0);
}

/* ------------------------------------------------------------------------
 *  Protected attribute
 * ---------------------------------------------------------------------- */

static void test_chop_is_protected(void) {
    /* The Chop symbol must carry Protected so users can't override it. */
    Expr* parsed = parse_expression("Attributes[Chop]");
    Expr* evaluated = evaluate(parsed);
    char* str = expr_to_string(evaluated);
    ASSERT_MSG(strstr(str, "Protected") != NULL,
               "expected Protected in attributes, got: %s", str);
    free(str);
    expr_free(parsed);
    expr_free(evaluated);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene
 * ---------------------------------------------------------------------- */

static void test_chop_memory_loop(void) {
    /* Repeatedly construct, evaluate, and free Chop expressions of every
     * shape exercised above.  Valgrind --leak-check=full should report
     * zero Mathilda-attributable losses across this loop. */
    const char* cases[] = {
        "Chop[1.0e-12]",
        "Chop[1.5]",
        "Chop[-1.0e-15]",
        "Chop[2.0 + 1.0e-12 I]",
        "Chop[1.0e-12 + 2.0 I]",
        "Chop[Complex[1, 1.0e-12]]",
        "Chop[{1.0e-12, 1.5, 2.0e-15}]",
        "Chop[f[1.0e-13, g[x, 1.0e-12]]]",
        "Chop[1.0e-5, 1.0e-4]",
        "Chop[1/2]",
        "Chop[Sin[x]]",
        NULL
    };
    for (int rep = 0; rep < 20; rep++) {
        for (int i = 0; cases[i]; i++) {
            Expr* p = parse_expression(cases[i]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
    }
}

/* ------------------------------------------------------------------------
 *  Main
 * ---------------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    /* Default tolerance */
    TEST(test_chop_below_default_tolerance);
    TEST(test_chop_above_default_tolerance);
    TEST(test_chop_at_boundary);
    TEST(test_chop_zero);

    /* Explicit delta */
    TEST(test_chop_with_explicit_real_delta);
    TEST(test_chop_with_integer_delta);
    TEST(test_chop_with_rational_delta);
    TEST(test_chop_with_negative_delta);

    /* Sign */
    TEST(test_chop_negative_real);

    /* Complex */
    TEST(test_chop_drops_small_imaginary_part);
    TEST(test_chop_keeps_machine_zero_real_part);
    TEST(test_chop_both_parts_small);
    TEST(test_chop_both_parts_above);
    TEST(test_chop_complex_with_exact_real_part);
    TEST(test_chop_complex_with_exact_imag_part);
    TEST(test_chop_exact_complex_passthrough);

    /* Structure */
    TEST(test_chop_inside_list);
    TEST(test_chop_inside_nested_list);
    TEST(test_chop_inside_plus);
    TEST(test_chop_inside_times);
    TEST(test_chop_inside_generic_head);
    TEST(test_chop_listable_threading);

    /* Exact passthrough */
    TEST(test_chop_exact_integers);
    TEST(test_chop_exact_rationals);
    TEST(test_chop_bigint_passthrough);
    TEST(test_chop_symbolic_passthrough);

    /* Rationalize cross-check from the spec */
    TEST(test_chop_rationalize_difference_default);
    TEST(test_chop_rationalize_difference_tight);

    /* Argument-shape */
    TEST(test_chop_no_args_unevaluated);
    TEST(test_chop_too_many_args_unevaluated);
    TEST(test_chop_symbolic_delta_unevaluated);

    /* Protected */
    TEST(test_chop_is_protected);

    /* Memory hygiene */
    TEST(test_chop_memory_loop);

    printf("All chop_tests passed.\n");
    return 0;
}
