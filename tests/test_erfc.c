/* Tests for the complementary error function Erfc[z].
 *
 * Covers exact special values (1, 0, 2, the -+I directed infinities),
 * symbolic fall-through (erfc is not odd, so Erfc[-x] stays put), machine
 * real & complex numerics, arbitrary-precision (MPFR) reals and complex,
 * precision tracking, the erfc(z) = 1 - erf(z) relation, derivatives, Series,
 * Listable threading, attributes, and arity errors. Reference values
 * cross-checked against mpmath / the running binary. */

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

/* ---- numeric helpers ------------------------------------------------ */

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
               strcmp(r->data.function.head->data.symbol, "Complex") == 0 &&
               r->data.function.arg_count == 2,
               "%s: expected Complex[..], got something else", input);
    Expr* re = r->data.function.args[0];
    Expr* im = r->data.function.args[1];
    ASSERT(re->type == EXPR_REAL && im->type == EXPR_REAL);
    ASSERT_MSG(fabs(re->data.real - er) <= tol && fabs(im->data.real - ei) <= tol,
               "%s: expected %.8g %+.8g I, got %.8g %+.8g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact special values ------------------------------------------- */

void test_erfc_exact() {
    assert_eval_eq("Erfc[0]", "1", 0);
    assert_eval_eq("Erfc[Infinity]", "0", 0);
    assert_eval_eq("Erfc[-Infinity]", "2", 0);
    assert_eval_eq("Erfc[ComplexInfinity]", "ComplexInfinity", 0);
    assert_eval_eq("Erfc[Indeterminate]", "Indeterminate", 0);
    /* Directed imaginary infinities are negated relative to Erf:
     * erfc(z) = 1 - erf(z), so Erfc[I Infinity] = DirectedInfinity[-I]. */
    assert_eval_eq("Erfc[I Infinity]", "DirectedInfinity[Complex[0, -1]]", 1);
    assert_eval_eq("Erfc[-I Infinity]", "DirectedInfinity[Complex[0, 1]]", 1);
    /* Listable over the standard infinity pair. */
    assert_eval_eq("Erfc[{Infinity, -Infinity}]", "{0, 2}", 0);
}

/* ---- symbolic fall-through (erfc is NOT odd) ------------------------ */

void test_erfc_symbolic() {
    assert_eval_eq("Erfc[x]", "Erfc[x]", 0);
    /* No odd-symmetry rewrite: Erfc[-x] stays as Erfc[-x]. */
    assert_eval_eq("Erfc[-x]", "Erfc[-x]", 0);
    assert_eval_eq("Erfc[a + b]", "Erfc[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_erfc_machine_real() {
    assert_close("Erfc[0.95]", erfc(0.95), 1e-12);
    assert_close("Erfc[1.5]", erfc(1.5), 1e-12);
    assert_close("Erfc[-1.0]", erfc(-1.0), 1e-12);   /* negatives via libm */
    assert_close("Erfc[3.0]", erfc(3.0), 1e-12);
    assert_close("Erfc[0.0]", 1.0, 0.0);
    /* Complementary relation at machine precision: erfc(x) + erf(x) = 1. */
    assert_close("Erfc[-2.5] + Erf[-2.5]", 1.0, 1e-12);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_erfc_arbitrary_precision() {
    /* N[Erfc[3/2], 50] (mpmath reference, prefix kept short of the last
     * rounding digit). */
    assert_eval_startswith("N[Erfc[3/2], 50]",
        "0.03389485352468927293302373835405214131858952074");
    /* N[Erfc[1], 35] = 1 - erf(1). */
    assert_eval_startswith("N[Erfc[1], 35]",
        "0.15729920705028513065877936491739");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("Erfc[1.500000000000000000000000`28]",
                           "0.0338948535246892729330237383");
}

/* ---- machine complex ------------------------------------------------ */

void test_erfc_machine_complex() {
    /* erfc(z) = 1 - erf(z). Erf[1.5 - I] = 1.0784 + 0.0279637 I, so
     * Erfc[1.5 - I] = -0.0783992 - 0.0279637 I. */
    assert_complex_close("Erfc[1.5 - I]", -0.0783992074989335, -0.0279637112386558, 1e-10);
    /* Erf[I] = 1.6504257587975428 I, so Erfc[I] = 1 - 1.6504257587975428 I. */
    assert_complex_close("Erfc[1.0 I]", 1.0, -1.6504257587975428, 1e-10);
    /* Conjugate symmetry: Erfc[conj z] = conj Erfc[z]. */
    assert_complex_close("Erfc[0.5 - 1.0 I]", -0.2048475583142180, 1.0244008816084459, 1e-10);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_erfc_arbitrary_complex() {
    /* N[Erfc[1/2 + I], 30] = 1 - Erf[1/2 + I] (mpmath reference). */
    assert_eval_startswith("N[Erfc[1/2 + I], 30]",
        "-0.20484755831421800270211268209");
}

/* ---- derivatives & series ------------------------------------------- */

void test_erfc_derivatives() {
    /* D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2) (note the minus sign vs Erf). */
    assert_eval_eq("D[Erfc[x], x]", "-(2 E^(-x^2))/Sqrt[Pi]", 0);
    /* Chain rule. */
    assert_eval_eq("D[Erfc[x^2], x]", "-(4 x E^(-x^4))/Sqrt[Pi]", 0);
    /* Series at the origin (generic D-based Taylor fallback); constant term
     * is the Erfc[0] = 1 exact value. */
    assert_eval_startswith("Series[Erfc[x], {x, 0, 3}]", "1 + -2/Sqrt[Pi] x");
}

/* ---- Listable & attributes ------------------------------------------ */

void test_erfc_listable() {
    assert_eval_eq("Erfc[{0}]", "{1}", 0);
    assert_eval_eq("Erfc[{}]", "{}", 0);
    /* Numeric threading. */
    Expr* e = parse_expression("Erfc[{0.5, 1.0, 1.5}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol, "List") == 0 &&
           r->data.function.arg_count == 3);
    double exp0[3] = { erfc(0.5), erfc(1.0), erfc(1.5) };
    for (int i = 0; i < 3; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= 1e-12,
                   "Erfc list element %d: expected %.12g got %.12g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_erfc_attributes() {
    SymbolDef* d = symtab_get_def("Erfc");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Erfc must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "Erfc must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Erfc must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_erfc_arity() {
    /* Wrong arity stays unevaluated (argx diagnostic to stderr). */
    assert_eval_eq("Erfc[]", "Erfc[]", 0);
    assert_eval_eq("Erfc[1, 2, 3]", "Erfc[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_erfc_exact);
    TEST(test_erfc_symbolic);
    TEST(test_erfc_machine_real);
    TEST(test_erfc_arbitrary_precision);
    TEST(test_erfc_machine_complex);
    TEST(test_erfc_arbitrary_complex);
    TEST(test_erfc_derivatives);
    TEST(test_erfc_listable);
    TEST(test_erfc_attributes);
    TEST(test_erfc_arity);

    printf("All Erfc tests passed.\n");
    return 0;
}
