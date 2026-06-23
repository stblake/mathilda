/* Tests for LegendreP: Legendre polynomials P_n(x), the associated Legendre
 * functions P_n^m(x), and the Legendre functions of types 1/2/3.
 *
 * Covers exact integer polynomials (including the P_10 headline case and the
 * negative-order reflection P_{-1-n} = P_n), the special value P_n(1) = 1,
 * Listable threading, machine and arbitrary-precision (MPFR) numerics through
 * the Gauss 2F1 series (real and complex), precision tracking, the associated
 * Legendre derivative form, the type-2/3 regularized-2F1 forms, symbolic
 * fall-through for exact non-integer / symbolic / out-of-range arguments,
 * argument-count diagnostics, and attributes. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- numeric helpers ------------------------------------------------- */

static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a Real result", input);
    double v = r->data.real;
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_real(input);
    ASSERT_MSG(fabs(v - expected) <= tol,
               "%s: expected %.12g, got %.12g", input, expected, v);
}

/* ---- exact integer polynomials -------------------------------------- */

void test_legendre_polynomials() {
    assert_eval_eq("Table[LegendreP[n, x], {n, 0, 5}]",
                   "{1, x, -1/2 + 3/2 x^2, -3/2 x + 5/2 x^3, "
                   "3/8 - 15/4 x^2 + 35/8 x^4, "
                   "15/8 x - 35/4 x^3 + 63/8 x^5}", 0);
    assert_eval_eq("LegendreP[0, x]", "1", 0);
    assert_eval_eq("LegendreP[1, x]", "x", 0);
    assert_eval_eq("LegendreP[2, x]", "-1/2 + 3/2 x^2", 0);
    /* The tenth Legendre polynomial (headline example). */
    assert_eval_eq("LegendreP[10, x]",
                   "-63/256 + 3465/256 x^2 - 15015/128 x^4 + 45045/128 x^6 "
                   "- 109395/256 x^8 + 46189/256 x^10", 0);
}

void test_legendre_negative_order() {
    /* P_{-1-n} = P_n. */
    assert_eval_eq("LegendreP[-1, x]", "1", 0);
    assert_eval_eq("LegendreP[-3, x]", "-1/2 + 3/2 x^2", 0);   /* = P_2 */
    assert_eval_eq("LegendreP[-6, x]",                          /* = P_5 */
                   "15/8 x - 35/4 x^3 + 63/8 x^5", 0);
}

void test_legendre_special_values() {
    /* P_n(1) = 1 for any order, including symbolic. */
    assert_eval_eq("LegendreP[2, 1]", "1", 0);
    assert_eval_eq("LegendreP[100, 1]", "1", 0);
    assert_eval_eq("LegendreP[n, 1]", "1", 0);
    assert_eval_eq("LegendreP[3/2, 1]", "1", 0);
    /* P_n(-1) = (-1)^n for integer n. */
    assert_eval_eq("LegendreP[4, -1]", "1", 0);
    assert_eval_eq("LegendreP[5, -1]", "-1", 0);
    /* Exact values at exact arguments. */
    assert_eval_eq("LegendreP[2, 2]", "11/2", 0);
    assert_eval_eq("LegendreP[2, 1/3]", "-1/3", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_legendre_listable() {
    assert_eval_eq("LegendreP[{1, 2, 3}, x]",
                   "{x, -1/2 + 3/2 x^2, -3/2 x + 5/2 x^3}", 0);
    assert_eval_eq("LegendreP[2, {0, 1}]", "{-1/2, 1}", 0);
    assert_eval_eq("LegendreP[{}, x]", "{}", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_legendre_machine() {
    /* Non-integer order, real argument (Gauss 2F1 series). */
    assert_close("LegendreP[2.5, 2]", 9.583124034019361, 1e-9);
    assert_close("LegendreP[1.5, 2]", 3.243939666040805, 1e-9);
    /* Integer order, inexact argument routes through the polynomial. */
    assert_close("LegendreP[4, 0.3]", 0.0729375, 1e-12);
    /* The numeric series reproduces the exact polynomial value. */
    assert_close("LegendreP[2., 0.5]", -0.125, 1e-12);
    assert_close("LegendreP[3., 0.5]", -0.4375, 1e-12);
}

void test_legendre_complex() {
    /* LegendreP[3/2 + I, 1.5 - I] = 5.20466 + 0.299479 I. */
    assert_close("Re[LegendreP[3/2 + I, 1.5 - I]]", 5.204663, 1e-5);
    assert_close("Im[LegendreP[3/2 + I, 1.5 - I]]", 0.299479, 1e-5);
}

/* ---- arbitrary-precision numerics ----------------------------------- */

void test_legendre_highprec() {
    /* N[LegendreP[3/2, 2], 50]. */
    assert_eval_startswith("N[LegendreP[3/2, 2], 50]",
                           "3.243939666040804915450228792970455767207515411017");
    /* Precision tracks the input precision. */
    assert_eval_startswith("LegendreP[3/2, 2`30]",
                           "3.24393966604080491545022879297");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_legendre_symbolic() {
    /* Exact non-integer order with exact argument stays symbolic. */
    assert_eval_eq("LegendreP[3/2, 2]", "LegendreP[3/2, 2]", 0);
    assert_eval_eq("LegendreP[n, x]", "LegendreP[n, x]", 0);
    assert_eval_eq("LegendreP[1/2, y]", "LegendreP[1/2, y]", 0);
}

/* ---- associated Legendre (type 1) ----------------------------------- */

void test_legendre_associated() {
    /* P_2^1(x) = -3 x Sqrt[1 - x^2]. */
    assert_eval_eq("LegendreP[2, 1, x]", "-3 x Sqrt[1 - x^2]", 0);
    /* P_2^2(2) = 3 (1 - 4) = -9. */
    assert_eval_eq("LegendreP[2, 2, 2]", "-9", 0);
    /* P_3^3(x) = -15 (1 - x^2)^(3/2). */
    assert_eval_eq("LegendreP[3, 3, x]", "-15 (1 - x^2)^(3/2)", 0);
    /* P_4^2(x) = (1 - x^2)(-15/2 + 105/2 x^2). */
    assert_eval_eq("LegendreP[4, 2, x]",
                   "(1 - x^2) (-15/2 + 105/2 x^2)", 0);
    /* m > n gives 0. */
    assert_eval_eq("LegendreP[2, 3, x]", "0", 0);
    /* m == 0 reduces to the ordinary polynomial. */
    assert_eval_eq("LegendreP[3, 0, x]", "-3/2 x + 5/2 x^3", 0);
}

/* ---- Legendre function types ---------------------------------------- */

void test_legendre_types() {
    /* a == 1 is the default and equals the three-argument form. */
    assert_eval_eq("LegendreP[2, 1, 1, z]", "-3 z Sqrt[1 - z^2]", 0);
    /* Types 2 and 3 give the regularized-2F1 surface forms. They are
     * mathematically equal to type 1 on the cut but differ in radical
     * grouping; check the exact deterministic output. */
    assert_eval_eq("LegendreP[2, 1, 2, z]",
                   "(Sqrt[1 + z] (-3 z + 3 z^2))/Sqrt[1 - z]", 0);
    assert_eval_eq("LegendreP[2, 1, 3, z]",
                   "(Sqrt[1 + z] (-3 z + 3 z^2))/Sqrt[-1 + z]", 0);
    /* Type 2 at a numeric interior point matches type 1 (same value). */
    assert_close("LegendreP[2, 1, 2, 0.4]", eval_real("LegendreP[2, 1, 1, 0.4]"),
                 1e-12);
}

/* ---- deferred cases stay symbolic ----------------------------------- */

void test_legendre_deferred() {
    /* Non-integer order/degree associated form is deferred. */
    assert_eval_eq("LegendreP[3/2, 1/2, x]", "LegendreP[3/2, 1/2, x]", 0);
    /* Negative degree is deferred. */
    assert_eval_eq("LegendreP[2, -1, x]", "LegendreP[2, -1, x]", 0);
    /* Out-of-range type is deferred. */
    assert_eval_eq("LegendreP[2, 1, 4, z]", "LegendreP[2, 1, 4, z]", 0);
    /* Numeric series outside its disk of convergence is deferred. */
    assert_eval_eq("LegendreP[1.5, 5.0]", "LegendreP[1.5, 5.0]", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_legendre_argcount() {
    /* Wrong arg counts stay unevaluated (an argb message goes to stderr). */
    assert_eval_eq("LegendreP[]", "LegendreP[]", 0);
    assert_eval_eq("LegendreP[3]", "LegendreP[3]", 0);
    assert_eval_eq("LegendreP[1, 2, 3, 4, 5]", "LegendreP[1, 2, 3, 4, 5]", 0);
}

void test_legendre_attributes() {
    SymbolDef* d = symtab_get_def("LegendreP");
    ASSERT_MSG(d != NULL, "LegendreP not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "LegendreP not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "LegendreP not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "LegendreP not Protected");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_legendre_polynomials);
    TEST(test_legendre_negative_order);
    TEST(test_legendre_special_values);
    TEST(test_legendre_listable);
    TEST(test_legendre_machine);
    TEST(test_legendre_complex);
    TEST(test_legendre_highprec);
    TEST(test_legendre_symbolic);
    TEST(test_legendre_associated);
    TEST(test_legendre_types);
    TEST(test_legendre_deferred);
    TEST(test_legendre_argcount);
    TEST(test_legendre_attributes);

    printf("All LegendreP tests passed.\n");
    return 0;
}
