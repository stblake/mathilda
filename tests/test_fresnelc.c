/* Tests for the Fresnel integral FresnelC[z] = Int_0^z Cos[Pi t^2/2] dt.
 *
 * Covers exact special values (0, +-Infinity -> +-1/2, +-I Infinity -> +-I/2,
 * ComplexInfinity/Indeterminate), machine real (convergent series and, for
 * large |x|, the asymptotic expansion), arbitrary-precision (MPFR) reals with
 * precision tracking, machine & arbitrary complex, odd symmetry, derivatives,
 * Taylor and asymptotic Series, the symbolic-n SeriesCoefficient closed form,
 * Listable threading, attributes and arity errors. Reference values
 * cross-checked against mpmath. */

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
               "%s: expected %.12g %+.12g I, got %.12g %+.12g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact special values ------------------------------------------- */

void test_fc_exact() {
    assert_eval_eq("FresnelC[0]", "0", 0);
    assert_eval_eq("FresnelC[Infinity]", "1/2", 0);
    assert_eval_eq("FresnelC[-Infinity]", "-1/2", 0);
    assert_eval_eq("FresnelC[I Infinity]", "1/2*I", 0);
    assert_eval_eq("FresnelC[-I Infinity]", "-1/2*I", 0);
    assert_eval_eq("FresnelC[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("FresnelC[Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("FresnelC[{-Infinity, Infinity, -I Infinity, I Infinity}]",
                   "{-1/2, 1/2, -1/2*I, 1/2*I}", 0);
}

/* ---- symbolic passthrough & odd symmetry ---------------------------- */

void test_fc_symbolic() {
    assert_eval_eq("FresnelC[x]", "FresnelC[x]", 0);
    assert_eval_eq("FresnelC[2]", "FresnelC[2]", 0);
    assert_eval_eq("FresnelC[a + b]", "FresnelC[a + b]", 0);
    assert_eval_eq("FresnelC[-x]", "-FresnelC[x]", 0);
    assert_eval_eq("FresnelC[-2 y]", "-FresnelC[2 y]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_fc_machine_real() {
    assert_close("FresnelC[0.5]", 0.49234422587144639, 1e-12);
    assert_close("FresnelC[1.0]", 0.77989340037682283, 1e-12);
    assert_close("FresnelC[1.8]", 0.33363292722155710, 1e-12);
    assert_close("FresnelC[2.0]", 0.48825340607534075, 1e-12);
    assert_close("FresnelC[1.5]", 0.44526117603982154, 1e-12);
    assert_close("FresnelC[2.5]", 0.45741300964177705, 1e-12);
    assert_close("FresnelC[3.5]", 0.53257243502800085, 1e-12);
    /* Odd: negative arguments. */
    assert_close("FresnelC[-1.8]", -0.33363292722155710, 1e-12);
    assert_close("FresnelC[-2.0]", -0.48825340607534075, 1e-12);
    /* Large |x|: exercises the asymptotic path (approaches 1/2). */
    assert_close("FresnelC[50.0]", 0.49999918943072797, 1e-10);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_fc_arbitrary_precision() {
    /* N[FresnelC[2], 50] -- the task-spec reference value. */
    assert_eval_startswith("N[FresnelC[2], 50]",
        "0.4882534060753407545002235033572610376883671545092");
    assert_eval_startswith("N[FresnelC[1], 40]",
        "0.7798934003768228294742064136526901366306");
    assert_eval_startswith("N[FresnelC[1/2], 40]",
        "0.49234422587144639287884366515668163");
    /* Large argument at high precision (asymptotic path). */
    assert_eval_startswith("N[FresnelC[50], 30]",
        "0.499999189430727967955810163");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("FresnelC[2.0000000000000000000000000000000]",
                           "0.4882534060753407545002235033");
}

/* ---- machine complex ------------------------------------------------ */

void test_fc_machine_complex() {
    assert_complex_close("FresnelC[2.5 + I]", 116.648061381, -105.228735671, 1e-6);
    /* On the imaginary axis, C(I y) = I C(y). */
    assert_complex_close("FresnelC[2.0 I]", 0.0, 0.48825340607534075, 1e-10);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_fc_arbitrary_complex() {
    /* Real part printed first; validates the high-precision complex path. */
    assert_eval_startswith("N[FresnelC[1/2 + I/3], 30]",
        "0.51910069121417857776167347263");
}

/* ---- derivatives ---------------------------------------------------- */

void test_fc_derivatives() {
    assert_eval_eq("D[FresnelC[x], x]", "Cos[1/2 Pi x^2]", 0);
    /* Chain rule. */
    assert_eval_eq("D[FresnelC[x^2], x]", "2 x Cos[1/2 Pi x^4]", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_fc_series_at_zero() {
    assert_eval_eq("Series[FresnelC[x], {x, 0, 9}]",
        "x + -1/40 Pi^2 x^5 + 1/3456 Pi^4 x^9 + O[x]^10", 0);
    assert_eval_eq("Series[FresnelC[x], {x, 0, 20}]",
        "x + -1/40 Pi^2 x^5 + 1/3456 Pi^4 x^9 + -1/599040 Pi^6 x^13 + "
        "1/175472640 Pi^8 x^17 + O[x]^21", 0);
    /* Composed argument. */
    assert_eval_eq("Series[FresnelC[x^2], {x, 0, 10}]",
        "x^2 + -1/40 Pi^2 x^10 + O[x]^11", 0);
}

void test_fc_series_at_infinity() {
    assert_eval_eq("Normal[Series[FresnelC[x], {x, Infinity, 3}]]",
        "1/2 - Cos[1/2 Pi x^2]/(Pi^2 x^3) + Sin[1/2 Pi x^2]/(Pi x)", 0);
}

void test_fc_series_coefficient() {
    /* Symbolic-n general term (Piecewise closed form). */
    assert_eval_eq("SeriesCoefficient[FresnelC[x], {x, 0, n}]",
        "Piecewise[{{(-1)^(1/4 (n - 1)) 2^(1 - n) Pi^(1/2 n)/"
        "(n Gamma[1/4 (1 + n)] Gamma[1/4 (3 + n)]), Mod[-1 + n, 4] == 0 && n >= 1}}, 0]", 0);
    /* Concrete integer indices reduce to the explicit coefficients. */
    assert_eval_eq("SeriesCoefficient[FresnelC[x], {x, 0, 1}]", "1", 0);
    assert_eval_eq("SeriesCoefficient[FresnelC[x], {x, 0, 5}]", "-1/40 Pi^2", 0);
    assert_eval_eq("SeriesCoefficient[FresnelC[x], {x, 0, 9}]", "1/3456 Pi^4", 0);
    assert_eval_eq("SeriesCoefficient[FresnelC[x], {x, 0, 4}]", "0", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_fc_listable() {
    assert_eval_eq("FresnelC[{}]", "{}", 0);
    assert_eval_eq("FresnelC[{0}]", "{0}", 0);
}

void test_fc_listable_numeric() {
    assert_close("FresnelC[{1.5, 2.5, 3.5}][[1]]", 0.44526117603982154, 1e-10);
    assert_close("FresnelC[{1.5, 2.5, 3.5}][[2]]", 0.45741300964177705, 1e-10);
    assert_close("FresnelC[{1.5, 2.5, 3.5}][[3]]", 0.53257243502800085, 1e-10);
}

void test_fc_attributes() {
    SymbolDef* d = symtab_get_def("FresnelC");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "FresnelC must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "FresnelC must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "FresnelC must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_fc_arity() {
    assert_eval_eq("FresnelC[]", "FresnelC[]", 0);
    assert_eval_eq("FresnelC[1, 2, 3]", "FresnelC[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_fc_exact);
    TEST(test_fc_symbolic);
    TEST(test_fc_machine_real);
    TEST(test_fc_arbitrary_precision);
    TEST(test_fc_machine_complex);
    TEST(test_fc_arbitrary_complex);
    TEST(test_fc_derivatives);
    TEST(test_fc_series_at_zero);
    TEST(test_fc_series_at_infinity);
    TEST(test_fc_series_coefficient);
    TEST(test_fc_listable);
    TEST(test_fc_listable_numeric);
    TEST(test_fc_attributes);
    TEST(test_fc_arity);

    printf("All FresnelC tests passed.\n");
    return 0;
}
