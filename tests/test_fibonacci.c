/* Tests for Fibonacci[n] (Fibonacci numbers) and Fibonacci[n, x]
 * (Fibonacci polynomials): exact integer orders, negative orders, the
 * polynomial recurrence, Listable threading, numeric (machine + MPFR /
 * complex) evaluation via the generalized closed form, and the symbolic
 * derivative rules. */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Exact Fibonacci numbers --------------------------------------- */

static void test_small_numbers(void) {
    assert_eval_eq("Fibonacci[0]", "0", 0);
    assert_eval_eq("Fibonacci[1]", "1", 0);
    assert_eval_eq("Fibonacci[2]", "1", 0);
    assert_eval_eq("Fibonacci[8]", "21", 0);
    assert_eval_eq("Fibonacci[10]", "55", 0);
    assert_eval_eq("Table[Fibonacci[n], {n, 20}]",
        "{1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765}", 0);
}

static void test_bigint(void) {
    assert_eval_eq("Fibonacci[100]", "354224848179261915075", 0);
    assert_eval_eq("Fibonacci[200]",
        "280571172992510140037611932413038677189525", 0);
}

static void test_negative_order(void) {
    /* F_{-n} = (-1)^{n+1} F_n */
    assert_eval_eq("Fibonacci[-1]", "1", 0);
    assert_eval_eq("Fibonacci[-2]", "-1", 0);
    assert_eval_eq("Fibonacci[-6]", "-8", 0);
    assert_eval_eq("Fibonacci[-7]", "13", 0);
}

static void test_listable(void) {
    assert_eval_eq("Fibonacci[{1, 2, 3, 4, 5}]", "{1, 1, 2, 3, 5}", 0);
}

/* ---- Fibonacci polynomials ----------------------------------------- */

static void test_polynomials(void) {
    assert_eval_eq("Fibonacci[0, x]", "0", 0);
    assert_eval_eq("Fibonacci[1, x]", "1", 0);
    assert_eval_eq("Fibonacci[2, x]", "x", 0);
    assert_eval_eq("Fibonacci[3, x]", "1 + x^2", 0);
    assert_eval_eq("Fibonacci[7, x]", "1 + 6 x^2 + 5 x^4 + x^6", 0);
    /* Values at fixed points and at zero. */
    assert_eval_eq("Fibonacci[1, 0]", "1", 0);
    assert_eval_eq("Fibonacci[0, 0]", "0", 0);
    assert_eval_eq("Fibonacci[3, 4]", "17", 0);
    assert_eval_eq("Table[Fibonacci[3, x], {x, 1, 5}]", "{2, 5, 10, 17, 26}", 0);
    /* Negative order: F_{-7}(x) = F_7(x). */
    assert_eval_eq("Fibonacci[-7, x]", "1 + 6 x^2 + 5 x^4 + x^6", 0);
    /* Exact complex argument stays exact (integer order). */
    assert_eval_eq("Fibonacci[5, 8 - I]", "3903 - 2064*I", 0);
}

/* ---- Non-integer (rational) order, exact closed form --------------- */

static void test_rational_order(void) {
    /* Exact non-integer order with x = 0 collapses to (1 - Cos[Pi n])/2. */
    assert_eval_eq("Fibonacci[1/2, 0]", "1/2", 0);
    assert_eval_eq("Fibonacci[1/3, 0]", "1/4", 0);
    assert_eval_eq("Fibonacci[2/3, 0]", "3/4", 0);
    assert_eval_eq("Fibonacci[3/2, 0]", "1/2", 0);
    /* Symbolic x keeps the rational-order call unevaluated. */
    assert_eval_eq("Fibonacci[1/2, x]", "Fibonacci[1/2, x]", 0);
}

/* ---- Symbolic (unevaluated) ---------------------------------------- */

static void test_symbolic(void) {
    assert_eval_eq("Fibonacci[n]", "Fibonacci[n]", 0);
    assert_eval_eq("Fibonacci[n, x]", "Fibonacci[n, x]", 0);
    /* Exact non-integer order without N stays symbolic. */
    assert_eval_eq("Fibonacci[15/17]", "Fibonacci[15/17]", 0);
}

/* ---- Numeric evaluation -------------------------------------------- */

static void test_numeric_machine(void) {
    /* Generalized Fibonacci polynomial at real order. */
    assert_eval_startswith("Fibonacci[5.8, 3]", "283.48");
    /* Rational order with an inexact x numericalizes the closed form. */
    assert_eval_startswith("Fibonacci[1/2, 3.2]", "0.49483");
    /* Numeric one-argument order: Cos[Pi/2] = 0, so this is
     * GoldenRatio^(1/2)/Sqrt[5] = 0.568864... */
    assert_eval_startswith("Fibonacci[0.5]", "0.5688");
}

#ifdef USE_MPFR
static void test_numeric_mpfr(void) {
    /* Arbitrary precision tracks the request. */
    assert_eval_startswith("N[Fibonacci[15/17], 50]", "0.9565199139243112250858226342769229864860");
    assert_eval_startswith("N[Fibonacci[15/17], 35]", "0.95651991392431122508582263427692298");
}
#endif

/* ---- Derivatives --------------------------------------------------- */

static void test_derivatives(void) {
    assert_eval_eq("D[Fibonacci[n, x], x]",
        "(2 n Fibonacci[-1 + n, x] + (-1 + n) x Fibonacci[n, x])/(4 + x^2)", 0);
    assert_eval_eq("D[Fibonacci[n], n]",
        "(GoldenRatio^(-n) (Log[GoldenRatio] GoldenRatio^(2 n) + Log[GoldenRatio] Cos[Pi n] + Pi Sin[Pi n]))/Sqrt[5]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_small_numbers);
    TEST(test_bigint);
    TEST(test_negative_order);
    TEST(test_listable);
    TEST(test_polynomials);
    TEST(test_rational_order);
    TEST(test_symbolic);
    TEST(test_numeric_machine);
#ifdef USE_MPFR
    TEST(test_numeric_mpfr);
#endif
    TEST(test_derivatives);

    printf("All Fibonacci tests passed.\n");
    return 0;
}
