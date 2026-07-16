/* Tests for the Fresnel integral FresnelS[z] = Int_0^z Sin[Pi t^2/2] dt.
 *
 * Covers exact special values (0, +-Infinity -> +-1/2, +-I Infinity -> -+I/2,
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

void test_fs_exact() {
    assert_eval_eq("FresnelS[0]", "0", 0);
    assert_eval_eq("FresnelS[Infinity]", "1/2", 0);
    assert_eval_eq("FresnelS[-Infinity]", "-1/2", 0);
    /* Note the sign flip vs FresnelC on the imaginary directions. */
    assert_eval_eq("FresnelS[I Infinity]", "-1/2*I", 0);
    assert_eval_eq("FresnelS[-I Infinity]", "1/2*I", 0);
    assert_eval_eq("FresnelS[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("FresnelS[Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("FresnelS[{-Infinity, Infinity, -I Infinity, I Infinity}]",
                   "{-1/2, 1/2, 1/2*I, -1/2*I}", 0);
}

/* ---- symbolic passthrough & odd symmetry ---------------------------- */

void test_fs_symbolic() {
    assert_eval_eq("FresnelS[x]", "FresnelS[x]", 0);
    assert_eval_eq("FresnelS[2]", "FresnelS[2]", 0);
    assert_eval_eq("FresnelS[a + b]", "FresnelS[a + b]", 0);
    assert_eval_eq("FresnelS[-x]", "-FresnelS[x]", 0);
    assert_eval_eq("FresnelS[-2 y]", "-FresnelS[2 y]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_fs_machine_real() {
    assert_close("FresnelS[0.5]", 0.064732432859999278, 1e-12);
    assert_close("FresnelS[1.0]", 0.43825914739035477, 1e-12);
    assert_close("FresnelS[1.8]", 0.45093876926758310, 1e-12);
    assert_close("FresnelS[2.0]", 0.34341567836369824, 1e-12);
    assert_close("FresnelS[1.5]", 0.69750496008209301, 1e-12);
    assert_close("FresnelS[2.5]", 0.61918175581959294, 1e-12);
    assert_close("FresnelS[3.5]", 0.41524801197243752, 1e-12);
    /* Odd: negative arguments. */
    assert_close("FresnelS[-1.8]", -0.45093876926758310, 1e-12);
    assert_close("FresnelS[-2.0]", -0.34341567836369824, 1e-12);
    /* Large |x|: exercises the asymptotic path (approaches 1/2). */
    assert_close("FresnelS[50.0]", 0.49363380258593874, 1e-10);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_fs_arbitrary_precision() {
    /* N[FresnelS[2], 50] -- the task-spec reference value. */
    assert_eval_startswith("N[FresnelS[2], 50]",
        "0.3434156783636982421953008159580684568865418122025");
    assert_eval_startswith("N[FresnelS[1], 40]",
        "0.43825914739035476607675669662515263749");
    /* FresnelS[3.4`45] -- the task-spec reference value. */
    assert_eval_startswith("FresnelS[3.4`45]",
        "0.42964946444392686440314980678606499145656");
    /* Large argument at high precision (asymptotic path). */
    assert_eval_startswith("N[FresnelS[50], 30]",
        "0.493633802585938741453268239798");
}

/* ---- machine complex ------------------------------------------------ */

void test_fs_machine_complex() {
    assert_complex_close("FresnelS[2.5 + I]", 105.728734983, 116.148016849, 1e-6);
    /* On the imaginary axis, S(I y) = -I S(y). */
    assert_complex_close("FresnelS[2.0 I]", 0.0, -0.34341567836369824, 1e-10);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_fs_arbitrary_complex() {
    /* FresnelS[1/2 + I/3] has a negative real part (printed first). */
    assert_eval_startswith("N[FresnelS[1/2 + I/3], 30]",
        "-0.020323606444949477286081733574");
}

/* ---- derivatives ---------------------------------------------------- */

void test_fs_derivatives() {
    assert_eval_eq("D[FresnelS[x], x]", "Sin[1/2 Pi x^2]", 0);
    /* Chain rule. */
    assert_eval_eq("D[FresnelS[x^2], x]", "2 x Sin[1/2 Pi x^4]", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_fs_series_at_zero() {
    assert_eval_eq("Series[FresnelS[x], {x, 0, 11}]",
        "1/6 Pi x^3 + -1/336 Pi^3 x^7 + 1/42240 Pi^5 x^11 + O[x]^12", 0);
    assert_eval_eq("Series[FresnelS[x], {x, 0, 20}]",
        "1/6 Pi x^3 + -1/336 Pi^3 x^7 + 1/42240 Pi^5 x^11 + -1/9676800 Pi^7 x^15 + "
        "1/3530096640 Pi^9 x^19 + O[x]^21", 0);
}

void test_fs_series_at_infinity() {
    assert_eval_eq("Normal[Series[FresnelS[x], {x, Infinity, 3}]]",
        "1/2 - Cos[1/2 Pi x^2]/(Pi x) - Sin[1/2 Pi x^2]/(Pi^2 x^3)", 0);
}

void test_fs_series_coefficient() {
    /* Symbolic-n general term (Piecewise closed form). */
    assert_eval_eq("SeriesCoefficient[FresnelS[x], {x, 0, n}]",
        "Piecewise[{{(-1)^(1/4 (n - 3)) 2^(1 - n) Pi^(1/2 n)/"
        "(n Gamma[1/4 (1 + n)] Gamma[1/4 (3 + n)]), Mod[-3 + n, 4] == 0 && n >= 3}}, 0]", 0);
    /* Concrete integer indices reduce to the explicit coefficients. */
    assert_eval_eq("SeriesCoefficient[FresnelS[x], {x, 0, 3}]", "1/6 Pi", 0);
    assert_eval_eq("SeriesCoefficient[FresnelS[x], {x, 0, 7}]", "-1/336 Pi^3", 0);
    assert_eval_eq("SeriesCoefficient[FresnelS[x], {x, 0, 1}]", "0", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_fs_listable() {
    assert_eval_eq("FresnelS[{}]", "{}", 0);
    assert_eval_eq("FresnelS[{0}]", "{0}", 0);
}

void test_fs_listable_numeric() {
    assert_close("FresnelS[{1.5, 2.5, 3.5}][[1]]", 0.69750496008209301, 1e-10);
    assert_close("FresnelS[{1.5, 2.5, 3.5}][[2]]", 0.61918175581959294, 1e-10);
    assert_close("FresnelS[{1.5, 2.5, 3.5}][[3]]", 0.41524801197243752, 1e-10);
}

void test_fs_attributes() {
    SymbolDef* d = symtab_get_def("FresnelS");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "FresnelS must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "FresnelS must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "FresnelS must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_fs_arity() {
    assert_eval_eq("FresnelS[]", "FresnelS[]", 0);
    assert_eval_eq("FresnelS[1, 2]", "FresnelS[1, 2]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_fs_exact);
    TEST(test_fs_symbolic);
    TEST(test_fs_machine_real);
    TEST(test_fs_arbitrary_precision);
    TEST(test_fs_machine_complex);
    TEST(test_fs_arbitrary_complex);
    TEST(test_fs_derivatives);
    TEST(test_fs_series_at_zero);
    TEST(test_fs_series_at_infinity);
    TEST(test_fs_series_coefficient);
    TEST(test_fs_listable);
    TEST(test_fs_listable_numeric);
    TEST(test_fs_attributes);
    TEST(test_fs_arity);

    printf("All FresnelS tests passed.\n");
    return 0;
}
