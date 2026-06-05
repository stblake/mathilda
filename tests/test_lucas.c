/* Tests for LucasL[n] (Lucas numbers) and LucasL[n, x] (Lucas
 * polynomials): exact integer orders, negative orders, the polynomial
 * recurrence, Listable threading, numeric (machine + MPFR / complex)
 * evaluation via the generalized closed form, and the symbolic derivative
 * rules. */

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

/* ---- Exact Lucas numbers ------------------------------------------- */

static void test_small_numbers(void) {
    assert_eval_eq("LucasL[0]", "2", 0);
    assert_eval_eq("LucasL[1]", "1", 0);
    assert_eval_eq("LucasL[2]", "3", 0);
    assert_eval_eq("LucasL[7]", "29", 0);
    assert_eval_eq("LucasL[10]", "123", 0);
    assert_eval_eq("Table[LucasL[n], {n, 20}]",
        "{1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521, 843, 1364, "
        "2207, 3571, 5778, 9349, 15127}", 0);
}

static void test_bigint(void) {
    assert_eval_eq("LucasL[100]", "792070839848372253127", 0);
    /* L_143; also the value of LucasL[143, 1]. */
    assert_eval_eq("LucasL[143]", "767772505664398093937756525279", 0);
    assert_eval_eq("LucasL[143, 1]", "767772505664398093937756525279", 0);
}

static void test_negative_order(void) {
    /* L_{-n} = (-1)^n L_n. */
    assert_eval_eq("LucasL[-1]", "-1", 0);
    assert_eval_eq("LucasL[-2]", "3", 0);
    assert_eval_eq("LucasL[-6]", "18", 0);
    assert_eval_eq("LucasL[-7]", "-29", 0);
    assert_eval_eq("LucasL[-11]", "-199", 0);
}

static void test_listable(void) {
    assert_eval_eq("LucasL[{1, 2, 3, 4, 5}]", "{1, 3, 4, 7, 11}", 0);
}

/* ---- Lucas polynomials --------------------------------------------- */

static void test_polynomials(void) {
    assert_eval_eq("LucasL[0, x]", "2", 0);
    assert_eval_eq("LucasL[1, x]", "x", 0);
    assert_eval_eq("LucasL[2, x]", "2 + x^2", 0);
    assert_eval_eq("LucasL[3, x]", "3 x + x^3", 0);
    assert_eval_eq("LucasL[7, x]", "7 x + 14 x^3 + 7 x^5 + x^7", 0);
    assert_eval_eq("Table[LucasL[n, x], {n, 5}]",
        "{x, 2 + x^2, 3 x + x^3, 2 + 4 x^2 + x^4, 5 x + 5 x^3 + x^5}", 0);
    /* Values at zero and at fixed points. */
    assert_eval_eq("LucasL[1, 0]", "0", 0);
    assert_eval_eq("LucasL[0, 0]", "2", 0);
    assert_eval_eq("Table[LucasL[10, x], {x, 1, 5}]",
        "{123, 6726, 154451, 1860498, 14250627}", 0);
    /* Negative order: L_{-n}(x) = (-1)^n L_n(x). */
    assert_eval_eq("LucasL[-7, x]", "-(7 x + 14 x^3 + 7 x^5 + x^7)", 0);
    assert_eval_eq("LucasL[-6, x]", "2 + 9 x^2 + 6 x^4 + x^6", 0);
    /* Exact complex argument stays exact (integer order). */
    assert_eval_eq("LucasL[3, 2 + I]", "8 + 14*I", 0);
    assert_eval_eq("LucasL[5, 8 - I]", "30168 - 20801*I", 0);
}

/* ---- Symbolic (unevaluated) ---------------------------------------- */

static void test_symbolic(void) {
    assert_eval_eq("LucasL[n]", "LucasL[n]", 0);
    assert_eval_eq("LucasL[n, x]", "LucasL[n, x]", 0);
    /* Exact non-integer order without N stays symbolic. */
    assert_eval_eq("LucasL[5/2]", "LucasL[5/2]", 0);
    /* Wrong argument counts stay unevaluated (no message), as for Fibonacci. */
    assert_eval_eq("LucasL[]", "LucasL[]", 0);
    assert_eval_eq("LucasL[5, 2, 1]", "LucasL[5, 2, 1]", 0);
}

/* ---- Numeric evaluation -------------------------------------------- */

static void test_numeric_machine(void) {
    /* L_{-11} via the closed form at a machine-real order. */
    assert_eval_startswith("LucasL[-11.]", "-199.");
    /* Generalized Lucas polynomial at real order. */
    assert_eval_startswith("LucasL[5.8, 3]", "1022.1");
    assert_eval_startswith("LucasL[2.5]", "3.330");
}

static void test_numeric_complex(void) {
    assert_eval_startswith("N[LucasL[1 + I/2]]", "0.0653384 + 0.755095*I");
}

#ifdef USE_MPFR
static void test_numeric_mpfr(void) {
    /* Arbitrary precision tracks the request. */
    assert_eval_startswith("N[LucasL[11/3], 50]",
        "5.923962652961955410135697862194012");
    assert_eval_startswith("N[LucasL[11/3], 35]",
        "5.9239626529619554101356978621940127");
}
#endif

/* ---- Derivatives --------------------------------------------------- */

static void test_derivatives(void) {
    assert_eval_eq("D[LucasL[n, x], x]",
        "(n (x LucasL[n, x] + 2 LucasL[-1 + n, x]))/(4 + x^2)", 0);
    assert_eval_eq("D[LucasL[n], n]",
        "GoldenRatio^(-n) (Log[GoldenRatio] GoldenRatio^(2 n) - "
        "Log[GoldenRatio] Cos[Pi n] - Pi Sin[Pi n])", 0);
    /* Cross-check the x-derivative against the closed polynomial derivative
     * at fixed points: d/dx L_3(x) = 3 + 3 x^2 -> 15 at x = 2, and
     * d/dx L_5(x) = 5 + 15 x^2 + 5 x^4 -> 25 at x = 1. (A Simplify[... - ...]
     * cancellation check would be cleaner, but the Simplify->0 path has a
     * pre-existing 48-byte leak unrelated to LucasL, so we substitute.) */
    assert_eval_eq("D[LucasL[3, x], x] /. x -> 2", "15", 0);
    assert_eval_eq("D[LucasL[5, x], x] /. x -> 1", "25", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_small_numbers);
    TEST(test_bigint);
    TEST(test_negative_order);
    TEST(test_listable);
    TEST(test_polynomials);
    TEST(test_symbolic);
    TEST(test_numeric_machine);
    TEST(test_numeric_complex);
#ifdef USE_MPFR
    TEST(test_numeric_mpfr);
#endif
    TEST(test_derivatives);

    printf("All LucasL tests passed.\n");
    return 0;
}
