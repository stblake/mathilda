/* Tests for Clip.
 *
 * Cover: the three documented call shapes (Clip[x], Clip[x, {min, max}],
 * Clip[x, {min, max}, {vmin, vmax}]); machine-precision and exact-number
 * inputs; symbolic constants (Pi, E) numericalized via N for the
 * comparison while returning the original symbolic form; Infinity and
 * -Infinity; pass-through for fully symbolic arguments; the complex
 * (non-real) rejection path with Clip::ncompl; explicit first-argument
 * threading over a List; nested lists; argument-shape errors; the
 * (NumericFunction, Protected) attributes; and a memory-hygiene loop.
 *
 * The headline cases mirror the Clip docstring (see
 * docs/spec/changelog/2026-05.md).
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
 *  Default interval [-1, 1]
 * ---------------------------------------------------------------------- */

static void test_clip_above_one(void) {
    assert_eval_eq("Clip[8.5]", "1", 0);
    assert_eval_eq("Clip[7.5]", "1", 0);
    assert_eval_eq("Clip[1.0001]", "1", 0);
}

static void test_clip_below_minus_one(void) {
    assert_eval_eq("Clip[-8.5]", "-1", 0);
    assert_eval_eq("Clip[-2]", "-1", 0);
}

static void test_clip_inside_default(void) {
    /* 0 is exactly in range -> returned unchanged. */
    assert_eval_eq("Clip[0]", "0", 0);
    assert_eval_eq("Clip[0.5]", "0.5", 0);
    assert_eval_eq("Clip[1]", "1", 0);
    assert_eval_eq("Clip[-1]", "-1", 0);
}

/* ------------------------------------------------------------------------
 *  Explicit {min, max}
 * ---------------------------------------------------------------------- */

static void test_clip_with_interval(void) {
    /* Spec example: Clip[-5/2, {-2, 2}] -> -2 because -5/2 = -2.5 < -2. */
    assert_eval_eq("Clip[-5/2, {-2, 2}]", "-2", 0);
    assert_eval_eq("Clip[3, {1, 5}]", "3", 0);
    assert_eval_eq("Clip[10, {1, 5}]", "5", 0);
    assert_eval_eq("Clip[-3, {1, 5}]", "1", 0);
}

static void test_clip_with_real_interval(void) {
    assert_eval_eq("Clip[0.25, {0.1, 0.5}]", "0.25", 0);
    assert_eval_eq("Clip[0.05, {0.1, 0.5}]", "0.1", 0);
    assert_eval_eq("Clip[0.7,  {0.1, 0.5}]", "0.5", 0);
}

static void test_clip_with_rational_interval(void) {
    /* All exact -- the rational comparison path. */
    assert_eval_eq("Clip[1/4, {1/8, 1/2}]", "1/4", 0);
    assert_eval_eq("Clip[1/16, {1/8, 1/2}]", "1/8", 0);
    assert_eval_eq("Clip[1, {1/8, 1/2}]", "1/2", 0);
}

/* ------------------------------------------------------------------------
 *  Explicit {min, max} with replacement {vmin, vmax}
 * ---------------------------------------------------------------------- */

static void test_clip_with_replacement_inside(void) {
    /* Spec example: Pi in [-9, 7] -> returns Pi unchanged (replacements
     * are only used outside the interval). */
    assert_eval_eq("Clip[Pi, {-9, 7}, {11, 28}]", "Pi", 0);
    /* 3 is inside [1, 5] -> 3 returned unchanged even though replacements
     * exist. */
    assert_eval_eq("Clip[3, {1, 5}, {0, 100}]", "3", 0);
}

static void test_clip_with_replacement_below(void) {
    assert_eval_eq("Clip[-3, {1, 5}, {0, 100}]", "0", 0);
    /* Symbolic replacement values are passed through verbatim. */
    assert_eval_eq("Clip[-3, {1, 5}, {a, b}]", "a", 0);
}

static void test_clip_with_replacement_above(void) {
    assert_eval_eq("Clip[6, {1, 5}, {0, 100}]", "100", 0);
    assert_eval_eq("Clip[6, {1, 5}, {a, b}]", "b", 0);
}

/* ------------------------------------------------------------------------
 *  Symbolic constants: numericalized only for the decision
 * ---------------------------------------------------------------------- */

static void test_clip_pi_default(void) {
    /* Pi ~ 3.14, outside [-1, 1], clamp up. */
    assert_eval_eq("Clip[Pi]", "1", 0);
}

static void test_clip_neg_pi_default(void) {
    assert_eval_eq("Clip[-Pi]", "-1", 0);
}

static void test_clip_pi_in_wider_range(void) {
    /* In [0, 10]: Pi survives as Pi (not 3.14...). */
    assert_eval_eq("Clip[Pi, {0, 10}]", "Pi", 0);
}

static void test_clip_e_default(void) {
    /* E ~ 2.71, outside [-1, 1], clamp up to 1. */
    assert_eval_eq("Clip[E]", "1", 0);
    /* In [0, 5], E survives. */
    assert_eval_eq("Clip[E, {0, 5}]", "E", 0);
}

/* ------------------------------------------------------------------------
 *  Infinity handling
 * ---------------------------------------------------------------------- */

static void test_clip_positive_infinity(void) {
    assert_eval_eq("Clip[Infinity]", "1", 0);
    assert_eval_eq("Clip[Infinity, {-10, 10}]", "10", 0);
    /* With explicit replacement: vmax is returned. */
    assert_eval_eq("Clip[Infinity, {-10, 10}, {a, b}]", "b", 0);
}

static void test_clip_negative_infinity(void) {
    assert_eval_eq("Clip[-Infinity]", "-1", 0);
    assert_eval_eq("Clip[-Infinity, {-10, 10}]", "-10", 0);
    assert_eval_eq("Clip[-Infinity, {-10, 10}, {a, b}]", "a", 0);
}

/* ------------------------------------------------------------------------
 *  Threading over a List in the first argument
 * ---------------------------------------------------------------------- */

static void test_clip_threads_first_arg(void) {
    /* Spec example. */
    assert_eval_eq("Clip[{-2, 0, 2}]", "{-1, 0, 1}", 0);
}

static void test_clip_threads_first_arg_with_interval(void) {
    /* The {1, 5} bounds list must NOT be threaded over. */
    assert_eval_eq("Clip[{0, 3, 10}, {1, 5}]", "{1, 3, 5}", 0);
}

static void test_clip_threads_first_arg_with_replacement(void) {
    assert_eval_eq("Clip[{-3, 3, 9}, {1, 5}, {0, 100}]",
                   "{0, 3, 100}", 0);
}

static void test_clip_threads_nested_lists(void) {
    /* The threading recurses through nested Lists because each Clip on
     * the inner List threads its first arg again. */
    assert_eval_eq("Clip[{{-2, 0}, {2, -10}}, {-1, 1}]",
                   "{{-1, 0}, {1, -1}}", 0);
}

/* ------------------------------------------------------------------------
 *  Complex (non-real) rejection
 * ---------------------------------------------------------------------- */

static void test_clip_rejects_complex(void) {
    /* Clip stays unevaluated and emits Clip::ncompl on stderr. */
    assert_eval_eq("Clip[2 - 3 I]", "Clip[2 - 3*I]", 0);
}

static void test_clip_complex_real_imag_workaround(void) {
    /* Spec example: Clip[Re[z]] + Clip[Im[z]] I -- chain works because
     * Re / Im strip the Complex first. */
    assert_eval_eq("Clip[Re[2 - 3 I]] + Clip[Im[2 - 3 I]] I", "1 - I", 0);
}

/* ------------------------------------------------------------------------
 *  Fully symbolic arguments -- stay unevaluated
 * ---------------------------------------------------------------------- */

static void test_clip_symbolic_passthrough(void) {
    assert_eval_eq("Clip[a]", "Clip[a]", 0);
    assert_eval_eq("Clip[a, {b, c}]", "Clip[a, {b, c}]", 0);
}

/* ------------------------------------------------------------------------
 *  Arbitrary precision via N
 * ---------------------------------------------------------------------- */

static void test_clip_arbitrary_precision_n(void) {
    /* From the spec discussion: 1/11 is below 1/7, so Clip returns 1/7,
     * then N[1/7, 50] gives 50+ digits.  We check the leading digits
     * because Mathilda's MPFR layer may print one trailing guard digit
     * over Mathematica's reference output. */
    assert_eval_startswith("N[Clip[1/11, {1/7, 5}], 50]",
                           "0.1428571428571428571428571428571428571428571428571");
}

static void test_clip_mpfr_bound(void) {
    /* An MPFR upper bound: 1/8 > 1/11, so the clip returns the upper
     * bound (1/11 at 100-digit precision). We assert only the prefix
     * because exact low-order digits depend on rounding. */
    assert_eval_startswith("Clip[1/8, {1/17, 1/11`100}]",
                           "0.0909090909090909090909");
}

/* ------------------------------------------------------------------------
 *  Argument-shape edge cases
 * ---------------------------------------------------------------------- */

static void test_clip_no_args_unevaluated(void) {
    assert_eval_eq("Clip[]", "Clip[]", 0);
}

static void test_clip_too_many_args_unevaluated(void) {
    assert_eval_eq("Clip[1, {0, 5}, {a, b}, foo]",
                   "Clip[1, {0, 5}, {a, b}, foo]", 0);
}

static void test_clip_interval_not_a_list(void) {
    /* Second arg must be a 2-element List. Anything else: unevaluated. */
    assert_eval_eq("Clip[1, foo]", "Clip[1, foo]", 0);
    assert_eval_eq("Clip[1, {1}]", "Clip[1, {1}]", 0);
    assert_eval_eq("Clip[1, {1, 2, 3}]", "Clip[1, {1, 2, 3}]", 0);
}

static void test_clip_replacement_not_a_list(void) {
    assert_eval_eq("Clip[1, {0, 5}, foo]", "Clip[1, {0, 5}, foo]", 0);
    assert_eval_eq("Clip[1, {0, 5}, {a}]", "Clip[1, {0, 5}, {a}]", 0);
}

/* ------------------------------------------------------------------------
 *  Attributes
 * ---------------------------------------------------------------------- */

static void test_clip_attributes(void) {
    Expr* parsed = parse_expression("Attributes[Clip]");
    Expr* evaluated = evaluate(parsed);
    char* str = expr_to_string(evaluated);
    ASSERT_MSG(strstr(str, "Protected") != NULL,
               "expected Protected in attributes, got: %s", str);
    ASSERT_MSG(strstr(str, "NumericFunction") != NULL,
               "expected NumericFunction in attributes, got: %s", str);
    free(str);
    expr_free(parsed);
    expr_free(evaluated);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene
 * ---------------------------------------------------------------------- */

static void test_clip_memory_loop(void) {
    /* Exercise every shape exercised above repeatedly so valgrind has
     * many chances to catch a leak in any path. The complex-rejection
     * branch and the threading branch are deliberately included. */
    const char* cases[] = {
        "Clip[8.5]",
        "Clip[-2]",
        "Clip[0.5]",
        "Clip[-5/2, {-2, 2}]",
        "Clip[Pi, {-9, 7}, {11, 28}]",
        "Clip[Infinity]",
        "Clip[-Infinity, {-10, 10}, {a, b}]",
        "Clip[{-2, 0, 2}]",
        "Clip[{0, 3, 10}, {1, 5}]",
        "Clip[{-3, 3, 9}, {1, 5}, {0, 100}]",
        "Clip[a]",
        "Clip[2 - 3 I]",
        "Clip[1/8, {1/17, 1/11}]",
        "N[Clip[1/11, {1/7, 5}], 50]",
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

    /* Default interval */
    TEST(test_clip_above_one);
    TEST(test_clip_below_minus_one);
    TEST(test_clip_inside_default);

    /* Explicit interval */
    TEST(test_clip_with_interval);
    TEST(test_clip_with_real_interval);
    TEST(test_clip_with_rational_interval);

    /* Replacement values */
    TEST(test_clip_with_replacement_inside);
    TEST(test_clip_with_replacement_below);
    TEST(test_clip_with_replacement_above);

    /* Symbolic constants */
    TEST(test_clip_pi_default);
    TEST(test_clip_neg_pi_default);
    TEST(test_clip_pi_in_wider_range);
    TEST(test_clip_e_default);

    /* Infinity */
    TEST(test_clip_positive_infinity);
    TEST(test_clip_negative_infinity);

    /* Threading */
    TEST(test_clip_threads_first_arg);
    TEST(test_clip_threads_first_arg_with_interval);
    TEST(test_clip_threads_first_arg_with_replacement);
    TEST(test_clip_threads_nested_lists);

    /* Complex */
    TEST(test_clip_rejects_complex);
    TEST(test_clip_complex_real_imag_workaround);

    /* Symbolic pass-through */
    TEST(test_clip_symbolic_passthrough);

    /* Arbitrary precision */
    TEST(test_clip_arbitrary_precision_n);
    TEST(test_clip_mpfr_bound);

    /* Argument shape */
    TEST(test_clip_no_args_unevaluated);
    TEST(test_clip_too_many_args_unevaluated);
    TEST(test_clip_interval_not_a_list);
    TEST(test_clip_replacement_not_a_list);

    /* Attributes */
    TEST(test_clip_attributes);

    /* Memory */
    TEST(test_clip_memory_loop);

    printf("All clip_tests passed.\n");
    return 0;
}
