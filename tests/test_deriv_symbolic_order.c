/*
 * test_deriv_symbolic_order.c -- Symbolic-order D + FactorialPower.
 *
 * Bug 8 -- D[x^n, {x, k}] = FactorialPower[n, k] x^(n-k)
 * Bug 9 -- D[x^n, {n, k}] = x^n Log[x]^k
 *
 *   The integer-only parse_var_spec previously fell back to treating
 *   the spec list itself as the variable when the order was symbolic,
 *   yielding 0. We now route symbolic orders through
 *   compute_deriv_symbolic_order, which handles:
 *     - Constant input: D[c, {x, k}] = 0 for c free of x.
 *     - Plus distribution: D[a + b, {x, k}] = D[a, {x, k}] + D[b, {x, k}].
 *     - Constant-factor pull-out: D[c f, {x, k}] = c D[f, {x, k}]
 *       when c is var-free.
 *     - D[Power[x, n], {x, k}] -> FactorialPower[n, k] * Power[x, n - k]
 *     - D[Power[b, x], {x, k}] -> Power[b, x] * Power[Log[b], k]
 *
 *   Forms outside that set fall back to leaving D[f, {var, k}]
 *   unevaluated rather than returning 0 silently.
 *
 * FactorialPower is also added as a first-class builtin: integer pairs
 * compute the product directly (GMP), symbolic n with concrete k <= 32
 * expands to an explicit product, and any other shape stays
 * unevaluated.
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

/* --- FactorialPower builtin ----------------------------------------- */
static void test_factorialpower_concrete(void) {
    /* FactorialPower[5, 3] = 5*4*3 = 60. */
    assert_eval_eq("FactorialPower[5, 3]", "60", 0);
    /* FactorialPower[10, 4] = 10*9*8*7 = 5040. */
    assert_eval_eq("FactorialPower[10, 4]", "5040", 0);
    /* k = 0 always 1. */
    assert_eval_eq("FactorialPower[n, 0]", "1", 0);
    assert_eval_eq("FactorialPower[5, 0]", "1", 0);
    /* Symbolic n with concrete k: expand to explicit product. */
    assert_eval_eq("FactorialPower[n, 3]", "n (-2 + n) (-1 + n)", 0);
    assert_eval_eq("FactorialPower[n, 1]", "n", 0);
    /* Both symbolic: stay unevaluated. */
    assert_eval_eq("FactorialPower[n, k]", "FactorialPower[n, k]", 0);
}

/* --- Bug 8: D[x^n, {x, k}] ------------------------------------------ */
static void test_bug8_dx_xn(void) {
    /* The Plus inside the exponent gets canonical "(-k+n)" ordering. */
    assert_eval_eq("D[x^n, {x, k}]",
                   "FactorialPower[n, k] x^(-k + n)", 0);
    /* Concrete n: FactorialPower[3, k] stays unevaluated. */
    assert_eval_eq("D[x^3, {x, k}]",
                   "FactorialPower[3, k] x^(3 - k)", 0);
    /* Constant pull-out: D[a x^n, {x, k}] = a x^(n-k) FactorialPower[n,k]. */
    assert_eval_eq("D[a*x^n, {x, k}]",
                   "a FactorialPower[n, k] x^(-k + n)", 0);
}

/* --- Bug 9: D[x^n, {n, k}] ------------------------------------------ */
static void test_bug9_dn_xn(void) {
    assert_eval_eq("D[x^n, {n, k}]", "Log[x]^k x^n", 0);
    /* Concrete base. */
    assert_eval_eq("D[2^n, {n, k}]", "Log[2]^k 2^n", 0);
}

/* --- Sin / Cos / Sinh / Cosh of a linear argument ------------------- */
static void test_trig_symbolic_order(void) {
    /* Basic D[Cos[x], {x, n}] = Cos[n Pi/2 + x] (canonical Plus order). */
    assert_eval_eq("D[Cos[x], {x, n}]",
                   "Cos[1/2 Pi n + x]", 0);
    assert_eval_eq("D[Sin[x], {x, n}]",
                   "Sin[1/2 Pi n + x]", 0);

    /* Linear chain rule: D[F[a x + b], {x, n}] = a^n F[n Pi/2 + a x + b]. */
    assert_eval_eq("D[Cos[3*x + 5], {x, n}]",
                   "3^n Cos[5 + 1/2 Pi n + 3 x]", 0);
    assert_eval_eq("D[Sin[a*x + b], {x, n}]",
                   "a^n Sin[b + 1/2 Pi n + a x]", 0);

    /* Hyperbolic: (-I)^n Cos[n Pi/2 - I x] and I (-I)^n Sin[n Pi/2 - I x].
     * The printer keeps the parentheses around -I to preserve the
     * (-I)^n vs -(I^n) distinction. */
    assert_eval_eq("D[Cosh[x], {x, n}]",
                   "(-I)^n Cos[1/2 Pi n - I x]", 0);
    assert_eval_eq("D[Sinh[x], {x, n}]",
                   "I (-I)^n Sin[1/2 Pi n - I x]", 0);

    /* Non-linear argument -> stays unevaluated. */
    assert_eval_eq("D[Cos[x^2], {x, n}]", "D[Cos[x^2], {x, n}]", 0);
    assert_eval_eq("D[Sin[Log[x]], {x, n}]", "D[Sin[Log[x]], {x, n}]", 0);
}

/* --- Algorithmic generality (Plus / Times pull-out) ----------------- */
static void test_symbolic_order_distributes(void) {
    /* Constant: 0 for any k >= 1. */
    assert_eval_eq("D[5, {x, k}]", "0", 0);
    assert_eval_eq("D[Pi, {x, k}]", "0", 0);
    /* Plus: distributes additively. */
    assert_eval_eq("D[a*x^n + b, {x, k}]",
                   "a FactorialPower[n, k] x^(-k + n)", 0);
    /* Multiple Power-in-x summands. */
    assert_eval_eq("D[x^n + x^m, {x, k}]",
                   "FactorialPower[m, k] x^(-k + m)"
                   " + FactorialPower[n, k] x^(-k + n)", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_factorialpower_concrete);
    TEST(test_bug8_dx_xn);
    TEST(test_bug9_dn_xn);
    TEST(test_trig_symbolic_order);
    TEST(test_symbolic_order_distributes);

    printf("All deriv_symbolic_order tests passed!\n");
    return 0;
}
