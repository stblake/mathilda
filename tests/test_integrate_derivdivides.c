/* test_integrate_derivdivides.c
 *
 * Tests for the derivative-divides substitution method
 * (Integrate`DerivativeDivides, Method -> "DerivativeDivides"; added
 * 2026-06-06).  Correctness is asserted by the universal predicate
 * Simplify[D[Integrate[f, x], x] - f] === 0 rather than by fixed output
 * strings, so the tests survive surface-form changes -- the exception is a
 * numerical branch-check that pins down which inverse-function branch was
 * selected.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Assert that the evaluated `input` is an unevaluated call with head
 * `head` (used for strict no-fallback / no-match cases). */
static void assert_head_unevaluated(const char* input, const char* head) {
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* Automatic cascade uses the quiet, branch-correct direct quotient strategy.
 * Integrate[Sin[x] Sqrt[1 - Cos[x]], x] is the headline case the rational /
 * Risch / table stages all miss. */
static void test_automatic_direct(void) {
    assert_eval_eq(
        "Simplify[D[Integrate[Sin[x] Sqrt[1 - Cos[x]], x], x]"
        " - Sin[x] Sqrt[1 - Cos[x]]]", "0", 0);
    /* The direct method picks the +2/3 (1-Cos[x])^(3/2) branch; a wrong
     * branch would differentiate to -Sin[x] Sqrt[1-Cos[x]].  Pin the branch
     * numerically at an interior point. */
    assert_eval_eq(
        "Abs[N[(D[Integrate[Sin[x] Sqrt[1 - Cos[x]], x], x]"
        " - Sin[x] Sqrt[1 - Cos[x]]) /. x -> 1.3]] < 0.000001", "True", 0);

    /* Other clean compositions handled by the direct strategy. */
    assert_eval_eq(
        "Simplify[D[Integrate[2 x Exp[x^2], x], x] - 2 x Exp[x^2]]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[1/(x Log[x]), x], x] - 1/(x Log[x])]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate[2 x Cos[x^2], x], x] - 2 x Cos[x^2]]", "0", 0);
}

/* The explicit method additionally runs the Eliminate/Solve branch search,
 * which closes integrands the direct strategy cannot because Cancel
 * canonicalises 1/Cos[x] -> Sec[x], breaking the syntactic u-substitution.
 * Here u = Cos[x] is a literal subexpression and the Eliminate path recovers
 * the clean reduced integrand -1/u^k. */
static void test_explicit_eliminate(void) {
    assert_eval_eq(
        "Simplify[D[Integrate[Sin[x]/Cos[x]^2, x, Method -> \"DerivativeDivides\"], x]"
        " - Sin[x]/Cos[x]^2]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`DerivativeDivides[Sin[x]/Cos[x]^4, x], x]"
        " - Sin[x]/Cos[x]^4]", "0", 0);
}

/* Radical substitution u = Sqrt[Tan[x]].  The Eliminate/Solve search closes
 * Integrate[Sqrt[Tan[x]], x] (reducing to the rational integral 2u^2/(1+u^4)),
 * which the direct strategy cannot.  This is now reachable from BOTH the
 * Automatic cascade and the explicit method.  Correctness is verified by the
 * differentiation predicate; the simplification of the radical-trig derivative
 * back to Sqrt[Tan[x]] exercises the trigrat quadratic-radical generators. */
static void test_radical_substitution(void) {
    /* Automatic cascade closes it. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[Tan[x]], x], x] - Sqrt[Tan[x]]]", "0", 0);
    /* Explicit method and package head close it too. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[Tan[x]], x, Method -> \"DerivativeDivides\"], x]"
        " - Sqrt[Tan[x]]]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`DerivativeDivides[Sqrt[Tan[x]], x], x]"
        " - Sqrt[Tan[x]]]", "0", 0);
    /* Pin the branch numerically at an interior point. */
    assert_eval_eq(
        "Abs[N[(D[Integrate[Sqrt[Tan[x]], x], x] - Sqrt[Tan[x]]) /. x -> 0.7]]"
        " < 0.000001", "True", 0);
}

/* Radical substitution u = Sqrt[Cot[x]].  Mirror of the Sqrt[Tan[x]] case, but
 * the radicand Cot = Cos/Sin is rational with the odd generator Sin in its
 * denominator: the trigrat verification gate previously evaluated Power[0,-1]
 * (1/0) during the conjugate rationalisation and returned Indeterminate, so the
 * differentiation gate failed and Integrate returned unevaluated.  With the
 * inverse-odd-generator clearing it closes end to end. */
static void test_radical_substitution_cot(void) {
    /* Automatic cascade closes it. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[Cot[x]], x], x] - Sqrt[Cot[x]]]", "0", 0);
    /* Explicit method and package head close it too. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[Cot[x]], x, Method -> \"DerivativeDivides\"], x]"
        " - Sqrt[Cot[x]]]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`DerivativeDivides[Sqrt[Cot[x]], x], x]"
        " - Sqrt[Cot[x]]]", "0", 0);
    /* Pin the branch numerically at an interior point. */
    assert_eval_eq(
        "Abs[N[(D[Integrate[Sqrt[Cot[x]], x], x] - Sqrt[Cot[x]]) /. x -> 0.7]]"
        " < 0.000001", "True", 0);
}

/* Method plumbing: both the option string and the explicit package head
 * route to the routine and close the headline case. */
static void test_method_plumbing(void) {
    assert_eval_eq(
        "Simplify[D[Integrate[Sin[x] Sqrt[1 - Cos[x]], x,"
        " Method -> \"DerivativeDivides\"], x] - Sin[x] Sqrt[1 - Cos[x]]]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`DerivativeDivides[Sin[x] Sqrt[1 - Cos[x]], x], x]"
        " - Sin[x] Sqrt[1 - Cos[x]]]", "0", 0);
}

/* Strict: a non-substitutable / non-elementary integrand returns unevaluated
 * with no fallback when the method is forced. */
static void test_strict_no_fallback(void) {
    /* Exp[x^2] is non-elementary: direct fails (a lone x survives), the
     * Eliminate path bails on the transcendental, so it bubbles back. */
    assert_head_unevaluated(
        "Integrate[Exp[x^2], x, Method -> \"DerivativeDivides\"]", "Integrate");
    /* The explicit head with no usable kernel (Sin[x] has only itself and x)
     * stays unevaluated. */
    assert_head_unevaluated(
        "Integrate`DerivativeDivides[Sin[x], x]", "Integrate`DerivativeDivides");
    /* Free of x: nothing to substitute. */
    assert_head_unevaluated(
        "Integrate`DerivativeDivides[5, x]", "Integrate`DerivativeDivides");
}

void test_integrate_derivdivides(void) {
    symtab_init();
    core_init();

    TEST(test_automatic_direct);
    TEST(test_explicit_eliminate);
    TEST(test_radical_substitution);
    TEST(test_radical_substitution_cot);
    TEST(test_method_plumbing);
    TEST(test_strict_no_fallback);

    printf("All Integrate DerivativeDivides tests passed!\n");
}

int main(void) {
    test_integrate_derivdivides();
    return 0;
}
