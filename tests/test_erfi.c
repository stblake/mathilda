/* Tests for the imaginary error function Erfi[z].
 *
 * Covers exact special values (0, +-Infinity, and the FINITE imaginary-axis
 * limits +-I), symbolic odd symmetry (Erfi[-x] = -Erfi[x]), machine real &
 * complex numerics, arbitrary-precision (MPFR) reals and complex, precision
 * tracking, derivatives, Series, Listable threading, attributes, and arity
 * errors. Reference values cross-checked against the running binary and the
 * task specification (erfi(z) = -I Erf[I z]). */

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
               "%s: expected %.15g, got %.15g", input, expected, v);
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
               "%s: expected %.10g %+.10g I, got %.10g %+.10g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact special values ------------------------------------------- */

void test_erfi_exact() {
    assert_eval_eq("Erfi[0]", "0", 0);
    assert_eval_eq("Erfi[Infinity]", "Infinity", 0);
    assert_eval_eq("Erfi[-Infinity]", "-Infinity", 0);
    assert_eval_eq("Erfi[ComplexInfinity]", "ComplexInfinity", 0);
    assert_eval_eq("Erfi[Indeterminate]", "Indeterminate", 0);
    /* The imaginary-axis limits are FINITE (unlike Erf): erfi(i y) -> i,
     * erfi(-i y) -> -i as y -> +Infinity. */
    assert_eval_eq("Erfi[I Infinity]", "I", 0);
    assert_eval_eq("Erfi[-I Infinity]", "-I", 0);
    /* Listable over the full special list from the spec. */
    assert_eval_eq("Erfi[{-Infinity, Infinity, I Infinity, -I Infinity}]",
                   "{-Infinity, Infinity, I, -I}", 0);
}

/* ---- symbolic odd symmetry ------------------------------------------ */

void test_erfi_symbolic() {
    assert_eval_eq("Erfi[x]", "Erfi[x]", 0);
    assert_eval_eq("Erfi[-x]", "-Erfi[x]", 0);
    assert_eval_eq("Erfi[-2 x]", "-Erfi[2 x]", 0);
    assert_eval_eq("Erfi[a + b]", "Erfi[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_erfi_machine_real() {
    assert_close("Erfi[0.0]", 0.0, 0.0);
    assert_close("Erfi[0.5]", 0.6149520946965110, 1e-12);
    assert_close("Erfi[1.0]", 1.6504257587975429, 1e-12);
    assert_close("Erfi[1.5]", 4.584733257284427, 1e-11);
    assert_close("Erfi[2.5]", 130.39575501324693, 1e-9);
    /* Odd in x: numeric negatives evaluate directly. */
    assert_close("Erfi[-0.5]", -0.6149520946965110, 1e-12);
    assert_close("Erfi[-2.5]", -130.39575501324693, 1e-9);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_erfi_arbitrary_precision() {
    /* N[Erfi[1/2], 50] (reference kept short of the last rounding digit). */
    assert_eval_startswith("N[Erfi[1/2], 50]",
        "0.6149520946965109808396811856236413930513456178954");
    /* N[Erfi[5/2], 35]. */
    assert_eval_startswith("N[Erfi[5/2], 35]",
        "130.3957550132469268137315308322899");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("Erfi[0.5000000000000000000000000000000000]",
                           "0.614952094696510980839681185623641");
    /* High-precision input stays MPFR (prompt timing case, abbreviated). */
    assert_eval_startswith("Erfi[1/2`60]",
                           "0.6149520946965109808396811856236413930513");
}

/* ---- machine complex ------------------------------------------------ */

void test_erfi_machine_complex() {
    /* Erfi[3/2 - I] = -0.70136046... - 1.84683301... I. */
    assert_complex_close("Erfi[1.5 - I]",
                         -0.7013604642514806, -1.8468330146085419, 1e-9);
    /* Conjugate symmetry: Erfi[conj z] = conj Erfi[z]. */
    assert_complex_close("Erfi[1.5 + I]",
                         -0.7013604642514806, 1.8468330146085419, 1e-9);
    /* erfi(1/2 + I). */
    assert_complex_close("Erfi[0.5 + 1.0 I]",
                         0.18797346722338331, 0.95070972831895717, 1e-10);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_erfi_arbitrary_complex() {
    /* N[Erfi[1/2 + I], 30]. */
    assert_eval_startswith("N[Re[Erfi[1/2 + I]], 30]",
        "0.18797346722338331362826381007");
    assert_eval_startswith("N[Im[Erfi[1/2 + I]], 30]",
        "0.95070972831895717380461182637");
}

/* ---- derivatives & series ------------------------------------------- */

void test_erfi_derivatives() {
    /* D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2) (positive exponent vs Erf). */
    assert_eval_eq("D[Erfi[x], x]", "(2 E^x^2)/Sqrt[Pi]", 0);
    /* Chain rule. */
    assert_eval_eq("D[Erfi[x^2], x]", "(4 x E^x^4)/Sqrt[Pi]", 0);
    /* Series at the origin via the generic D-based Taylor fallback. */
    assert_eval_startswith("Series[Erfi[x], {x, 0, 3}]", "2/Sqrt[Pi] x");
}

/* ---- Listable & attributes ------------------------------------------ */

void test_erfi_listable() {
    assert_eval_eq("Erfi[{0}]", "{0}", 0);
    assert_eval_eq("Erfi[{}]", "{}", 0);
    /* Numeric threading. */
    Expr* e = parse_expression("Erfi[{0.5, 1.5, 2.5}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol.name, "List") == 0 &&
           r->data.function.arg_count == 3);
    double exp0[3] = { 0.6149520946965110, 4.584733257284427, 130.39575501324693 };
    double tol[3]  = { 1e-12, 1e-11, 1e-9 };
    for (int i = 0; i < 3; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= tol[i],
                   "Erfi list element %d: expected %.15g got %.15g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_erfi_attributes() {
    SymbolDef* d = symtab_get_def("Erfi");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Erfi must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "Erfi must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Erfi must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_erfi_arity() {
    /* Wrong arity stays unevaluated (argx diagnostic to stderr). */
    assert_eval_eq("Erfi[]", "Erfi[]", 0);
    assert_eval_eq("Erfi[1, 2, 3]", "Erfi[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_erfi_exact);
    TEST(test_erfi_symbolic);
    TEST(test_erfi_machine_real);
    TEST(test_erfi_arbitrary_precision);
    TEST(test_erfi_machine_complex);
    TEST(test_erfi_arbitrary_complex);
    TEST(test_erfi_derivatives);
    TEST(test_erfi_listable);
    TEST(test_erfi_attributes);
    TEST(test_erfi_arity);

    printf("All Erfi tests passed.\n");
    return 0;
}
