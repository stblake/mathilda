/* Tests for the Pochhammer symbol (rising factorial) Pochhammer[a, n].
 *
 * Covers exact integer expansion (numeric and symbolic a), negative n
 * (reciprocal product), the a = 0 and Infinity short-circuits, n = 0,
 * exact half-integer reduction via the Gamma ratio, machine real & complex
 * numerics, arbitrary-precision (MPFR) precision tracking, symbolic
 * fall-through, Listable threading, attributes and arity diagnostics. */

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

/* ---- numeric helpers (mirrors test_gamma.c) ------------------------- */

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

static void assert_complex_close(const char* input, double er, double ei, double tol) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_FUNCTION &&
               r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
               r->data.function.arg_count == 2,
               "%s: expected Complex[..], got something else", input);
    Expr* re = r->data.function.args[0];
    Expr* im = r->data.function.args[1];
    ASSERT(re->type == EXPR_REAL && im->type == EXPR_REAL);
    ASSERT_MSG(fabs(re->data.real - er) <= tol && fabs(im->data.real - ei) <= tol,
               "%s: expected %.6g %+.6g I, got %.6g %+.6g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact integer n: numeric a ------------------------------------- */

void test_poch_integer_numeric() {
    assert_eval_eq("Pochhammer[10, 6]", "3603600", 0);   /* 10*11*...*15 */
    assert_eval_eq("Pochhammer[5, 3]", "210", 0);         /* 5*6*7 */
    assert_eval_eq("Pochhammer[1, 5]", "120", 0);         /* = 5! */
    assert_eval_eq("Pochhammer[2, 4]", "120", 0);         /* 2*3*4*5 */
    /* BigInt territory: Pochhammer[1, 25] = 25!. */
    assert_eval_eq("Pochhammer[1, 25]", "15511210043330985984000000", 0);
    /* Rational a stays exact. */
    assert_eval_eq("Pochhammer[1/2, 3]", "15/8", 0);      /* (1/2)(3/2)(5/2) */
}

/* ---- exact integer n: symbolic a (polynomial product) --------------- */

void test_poch_integer_symbolic() {
    assert_eval_eq("Pochhammer[n, 5]", "n (1 + n) (2 + n) (3 + n) (4 + n)", 0);
    assert_eval_eq("Pochhammer[x, 4]", "x (1 + x) (2 + x) (3 + x)", 0);
    assert_eval_eq("Pochhammer[x, 1]", "x", 0);
    assert_eval_eq("Table[Pochhammer[x, 2], {x, 1, 5}]", "{2, 6, 12, 20, 30}", 0);
}

/* ---- n = 0 ----------------------------------------------------------- */

void test_poch_zero_order() {
    assert_eval_eq("Pochhammer[a, 0]", "1", 0);
    assert_eval_eq("Pochhammer[5, 0]", "1", 0);
    assert_eval_eq("Pochhammer[0, 0]", "1", 0);
    assert_eval_eq("Pochhammer[Infinity, 0]", "1", 0);
}

/* ---- negative n: reciprocal product --------------------------------- */

void test_poch_negative_order() {
    assert_eval_eq("Pochhammer[n, -5]",
                   "1/((-5 + n) (-4 + n) (-3 + n) (-2 + n) (-1 + n))", 0);
    assert_eval_eq("Pochhammer[10, -3]", "1/504", 0);   /* 1/(9*8*7) */
    assert_eval_eq("Pochhammer[5, -2]", "1/12", 0);     /* 1/(4*3) */
    assert_eval_eq("Pochhammer[0, -3]", "-1/6", 0);     /* 1/((-1)(-2)(-3)) */
}

/* ---- a = 0 and Infinity short-circuits ------------------------------ */

void test_poch_zero_base() {
    assert_eval_eq("Pochhammer[0, 1]", "0", 0);
    assert_eval_eq("Pochhammer[0, 5]", "0", 0);
    /* Beyond the product cap, the a = 0 short-circuit still gives 0. */
    assert_eval_eq("Pochhammer[0, 1285]", "0", 0);
}

void test_poch_infinity() {
    assert_eval_eq("Pochhammer[Infinity, 2]", "Infinity", 0);
    assert_eval_eq("Pochhammer[Infinity, 1]", "Infinity", 0);
}

/* ---- exact half-integer reduction via the Gamma ratio --------------- */

void test_poch_half_integer() {
    assert_eval_eq("Pochhammer[3/2, 1/2]", "2/Sqrt[Pi]", 0);
    assert_eval_eq("Pochhammer[1/2, 1/2]", "1/Sqrt[Pi]", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_poch_machine_real() {
    /* Pochhammer[2.4, 8.5] = Gamma[10.9]/Gamma[2.4] = 2310224.6732407078. */
    assert_close("Pochhammer[2.4, 8.5]", 2310224.6732407078, 1e-2);
    /* Integer n with a machine-real base uses the product path. */
    assert_close("Pochhammer[2.5, 3]", 2.5 * 3.5 * 4.5, 1e-9);
}

void test_poch_machine_complex() {
    /* Pochhammer[2.+5 I, 8 I] = 2.13868e-6 - 1.42187e-5 I. */
    assert_complex_close("Pochhammer[2.+5 I, 8 I]",
                         2.13868e-06, -1.42187e-05, 1e-9);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_poch_arbitrary_precision() {
    /* Exact rational product, numericalised to 50 digits. */
    assert_eval_startswith("N[Pochhammer[1/3, 7], 50]",
                           "505.971650663008687700045724737");
    /* Precision of the output tracks the precision of the MPFR input. */
    assert_eval_startswith("Pochhammer[1.011111111111000000000000000, 8]",
                           "41552.2758490877803808");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_poch_symbolic() {
    assert_eval_eq("Pochhammer[a, n]", "Pochhammer[a, n]", 0);
    assert_eval_eq("Pochhammer[a, b]", "Pochhammer[a, b]", 0);
    /* Non-integer order with symbolic base: no closed form. */
    assert_eval_eq("Pochhammer[a, 1/2]", "Pochhammer[a, 1/2]", 0);
    /* Exact rationals whose Gamma ratio does not reduce stay symbolic. */
    assert_eval_eq("Pochhammer[1/2, 1/3]", "Pochhammer[1/2, 1/3]", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_poch_listable() {
    assert_eval_eq("Pochhammer[{2, 3, 5, 7, 11}, 3]",
                   "{24, 60, 210, 504, 1716}", 0);
    assert_eval_eq("Pochhammer[{1, 2, 3}, 2]", "{2, 6, 12}", 0);
}

void test_poch_attributes() {
    SymbolDef* d = symtab_get_def("Pochhammer");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Pochhammer must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "Pochhammer must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Pochhammer must be Protected");
}

/* ---- arity diagnostics (stays unevaluated) -------------------------- */

void test_poch_arity() {
    /* Wrong argument counts emit Pochhammer::argrx and stay unevaluated. */
    assert_eval_eq("Pochhammer[]", "Pochhammer[]", 0);
    assert_eval_eq("Pochhammer[1, 2, 3]", "Pochhammer[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_poch_integer_numeric);
    TEST(test_poch_integer_symbolic);
    TEST(test_poch_zero_order);
    TEST(test_poch_negative_order);
    TEST(test_poch_zero_base);
    TEST(test_poch_infinity);
    TEST(test_poch_half_integer);
    TEST(test_poch_machine_real);
    TEST(test_poch_machine_complex);
    TEST(test_poch_arbitrary_precision);
    TEST(test_poch_symbolic);
    TEST(test_poch_listable);
    TEST(test_poch_attributes);
    TEST(test_poch_arity);

    printf("All Pochhammer tests passed.\n");
    return 0;
}
