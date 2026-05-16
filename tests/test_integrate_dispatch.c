/* test_integrate_dispatch.c
 *
 * Smoke tests for the three-stage Integrate cascade and the Method
 * option added 2026-05-15.  Covers:
 *  - Cascade routing: rational integrand goes through BronsteinRational, an
 *    elementary one through RischNorman, an explicit CRCTable call
 *    survives without infinite-looping on the formerly-divergent
 *    inputs (Formula 49 family).
 *  - Method option: strict passthrough; unknown method bubbles back
 *    with Integrate::method.
 *  - Termination: pathological inputs that would have looped under
 *    the pre-2026-05-15 CRC table all return within a small budget.
 *
 * What we deliberately DO NOT test here:
 *  - Numerical correctness of every CRC formula.  Most rules don't
 *    fire today because Mathilda's matcher does not fully support
 *    /;-guarded multi-arg patterns (a separate work item); the
 *    cascade is in place for when that lands.
 *  - End-to-end Risch correctness (covered by intrischnorman_tests).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static double elapsed_seconds(const char* input) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    if (result) expr_free(result);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    return (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
}

static void test_cascade_rational(void) {
    /* Default cascade — first stage closes. */
    assert_eval_eq("Integrate[x^2, x]", "1/3 x^3", 0);
    assert_eval_eq("Integrate[1/(x^2 + 1), x]", "ArcTan[x]", 0);
}

static void test_cascade_risch(void) {
    /* Default cascade — Risch closes once Rational gives up. */
    Expr* parsed = parse_expression("Integrate[Sin[x], x]");
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    /* Any non-NULL, non-Integrate-headed result is success — the
     * specific form differs across Risch-Norman builds. */
    ASSERT(result != NULL);
    ASSERT(!(result->type == EXPR_FUNCTION
             && result->data.function.head
             && result->data.function.head->type == EXPR_SYMBOL
             && strcmp(result->data.function.head->data.symbol, "Integrate") == 0));
    expr_free(result);
}

static void test_method_strict_rational(void) {
    /* Method -> "BronsteinRational" closes the polynomial case. */
    assert_eval_eq("Integrate[x^3, x, Method -> \"BronsteinRational\"]", "1/4 x^4", 0);

    /* Method -> "BronsteinRational" on a non-rational integrand bubbles back
     * with no Risch / CRC fallback. */
    Expr* parsed = parse_expression("Integrate[Sin[x], x, Method -> \"BronsteinRational\"]");
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL
        && result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol, "Integrate") == 0);
    expr_free(result);
}

static void test_method_invalid(void) {
    /* Method -> "Bogus" emits Integrate::method (to stderr) and bubbles
     * back unevaluated with no crash. */
    Expr* parsed = parse_expression("Integrate[x, y, Method -> \"Bogus\"]");
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol, "Integrate") == 0);
    expr_free(result);
}

static void test_crctable_termination(void) {
    /* These integrands match the LHS of Formula 49 (line 119 in the
     * CRC .m file) with non-positive integer exponents — the
     * pre-2026-05-15 rule would have spun forever.  With the
     * IntegerQ + bound guard added in this change, each call must
     * return promptly. */
    double t1 = elapsed_seconds(
        "Integrate[1/(x^2 - 1)^(-3), x, Method -> \"CRCTable\"]");
    double t2 = elapsed_seconds(
        "Integrate[1/(x^2 - 1)^(-5), x, Method -> \"CRCTable\"]");
    double t3 = elapsed_seconds(
        "Integrate[1/(x^2 - 1)^(-7), x, Method -> \"CRCTable\"]");

    /* A correct termination guard returns within a few hundred
     * milliseconds even at REPL warm-up.  Allow 5s wall-clock as a
     * very generous ceiling — Formula 49's old form would lock the
     * REPL indefinitely. */
    ASSERT_MSG(t1 < 5.0, "Formula 49 (n=-3) took %.2fs (>5s) — likely diverging.", t1);
    ASSERT_MSG(t2 < 5.0, "Formula 49 (n=-5) took %.2fs (>5s) — likely diverging.", t2);
    ASSERT_MSG(t3 < 5.0, "Formula 49 (n=-7) took %.2fs (>5s) — likely diverging.", t3);
}

void test_integrate_dispatch(void) {
    symtab_init();
    core_init();

    TEST(test_cascade_rational);
    TEST(test_cascade_risch);
    TEST(test_method_strict_rational);
    TEST(test_method_invalid);
    TEST(test_crctable_termination);

    printf("All Integrate dispatch tests passed!\n");
}

int main(void) {
    test_integrate_dispatch();
    return 0;
}
