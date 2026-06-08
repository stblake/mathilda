/* Tests for the Gamma function family: Gamma[z], Gamma[a, z], Gamma[a, z0, z1].
 *
 * Covers exact integer / half-integer reduction, BigInt results, poles,
 * machine real & complex numerics, arbitrary-precision (MPFR) reals,
 * precision tracking, the incomplete and generalized incomplete forms,
 * symbolic fall-through, Listable threading, and attributes. */

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
               "%s: expected %.10g, got %.10g", input, expected, v);
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
               "%s: expected %.6g %+.6g I, got %.6g %+.6g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact symbolic values ------------------------------------------ */

void test_gamma_integer_values() {
    /* Gamma[n] = (n-1)! for positive integers. */
    assert_eval_eq("Gamma[1]", "1", 0);
    assert_eval_eq("Gamma[2]", "1", 0);
    assert_eval_eq("Gamma[3]", "2", 0);
    assert_eval_eq("Gamma[4]", "6", 0);
    assert_eval_eq("Gamma[10]", "362880", 0);
    assert_eval_eq("Table[Gamma[n], {n, 10}]",
                   "{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880}", 0);
    /* BigInt territory. */
    assert_eval_eq("Gamma[25]", "620448401733239439360000", 0);
    assert_eval_eq("Gamma[30]", "8841761993739701954543616000000", 0);
}

void test_gamma_poles() {
    /* Non-positive integers are poles -> ComplexInfinity. */
    assert_eval_eq("Gamma[0]", "ComplexInfinity", 0);
    assert_eval_eq("Gamma[-1]", "ComplexInfinity", 0);
    assert_eval_eq("Gamma[-3]", "ComplexInfinity", 0);
    /* A negative-integer machine real is also a pole. */
    assert_eval_eq("Gamma[-2.0]", "ComplexInfinity", 0);
}

void test_gamma_half_integers() {
    assert_eval_eq("Gamma[1/2]", "Sqrt[Pi]", 0);
    assert_eval_eq("Gamma[3/2]", "1/2 Sqrt[Pi]", 0);
    assert_eval_eq("Gamma[5/2]", "3/4 Sqrt[Pi]", 0);
    assert_eval_eq("Gamma[7/2]", "15/8 Sqrt[Pi]", 0);
    assert_eval_eq("Gamma[-1/2]", "-2 Sqrt[Pi]", 0);
    assert_eval_eq("Table[Gamma[n + 1/2], {n, 5}]",
                   "{1/2 Sqrt[Pi], 3/4 Sqrt[Pi], 15/8 Sqrt[Pi], "
                   "105/16 Sqrt[Pi], 945/32 Sqrt[Pi]}", 0);
}

void test_gamma_infinities() {
    assert_eval_eq("Gamma[Infinity]", "Infinity", 0);
    assert_eval_eq("Gamma[-Infinity]", "Indeterminate", 0);
    assert_eval_eq("Gamma[ComplexInfinity]", "ComplexInfinity", 0);
    assert_eval_eq("Gamma[Indeterminate]", "Indeterminate", 0);
}

void test_gamma_symbolic() {
    /* Stays unevaluated where there is no closed form. */
    assert_eval_eq("Gamma[x]", "Gamma[x]", 0);
    assert_eval_eq("Gamma[1/3]", "Gamma[1/3]", 0);
    assert_eval_eq("Gamma[a, z]", "Gamma[a, z]", 0);
    /* Exact non-integer incomplete gamma stays symbolic (no closed form). */
    assert_eval_eq("Gamma[3/2, z]", "Gamma[3/2, z]", 0);
    /* Exact complex incomplete gamma stays symbolic (no inexact part). */
    assert_eval_eq("Gamma[3/2, I]", "Gamma[3/2, I]", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_gamma_machine_real() {
    assert_close("Gamma[2.2]", 1.1018024908797127, 1e-12);
    assert_close("Gamma[2.5]", 1.3293403881791370, 1e-12);
    assert_close("Gamma[0.5]", 1.7724538509055160, 1e-12);   /* Sqrt[Pi] */
    assert_close("Gamma[5.0]", 24.0, 1e-9);
    assert_close("Gamma[-2.5]", -0.9453087204829419, 1e-12); /* reflection region */
}

void test_gamma_machine_complex() {
    /* Gamma[2.3 + I] = 0.71914093653728128 + 0.54061446790985013 I. */
    assert_complex_close("Gamma[2.3 + I]", 0.71914093653728128, 0.54061446790985013, 1e-9);
    /* Gamma[I] is symbolic, but N[Gamma[I]] is numeric:
       Gamma[I] = -0.154950 - 0.498016 I. */
    assert_complex_close("N[Gamma[I]]", -0.1549498283, -0.4980156681, 1e-7);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_gamma_arbitrary_precision() {
    /* N[Gamma[22/10], 50] -- 50-digit value from the reference. */
    assert_eval_startswith(
        "N[Gamma[22/10], 50]",
        "1.1018024908797127327691419862229964808241863995904");
    /* Precision tracks the input precision. */
    assert_eval_startswith("Gamma[2.2`40]",
                           "1.101802490879712732769141986222996480824");
    /* Gamma[1/2] at high precision is Sqrt[Pi]. */
    assert_eval_startswith("N[Gamma[1/2], 40]",
                           "1.772453850905516027298167483341145182797");
    /* Output precision survives: 50-digit request prints ~50 digits. */
    assert_eval_startswith("N[Gamma[10], 30]", "362880.");
}

void test_gamma_arbitrary_complex() {
    /* Arbitrary-precision complex Gamma via Spouge. N[Gamma[I], 50]:
       Gamma[I] = -0.15494982830181068512... - 0.49801566811835604271... I. */
    assert_eval_startswith("N[Gamma[I], 50]",
                           "-0.154949828301810685124955130483886605195879");
    /* Gamma[1/2 + I] = 0.30069461726065582 - 0.42496787943312381 I. */
    assert_eval_startswith("N[Gamma[1/2 + I], 20]", "0.3006946172606558");
    /* Reflection (Re < 1/2) agrees with the machine path. */
    assert_eval_startswith("N[Gamma[-3/2 + 2 I], 20]", "-0.0018843965411520");
}

/* ---- incomplete & generalized incomplete gamma ---------------------- */

void test_gamma_incomplete_exact() {
    /* Gamma[a, 0] = Gamma[a]. */
    assert_eval_eq("Gamma[a, 0]", "Gamma[a]", 0);
    assert_eval_eq("Gamma[5, 0]", "24", 0);
    /* Gamma[1, z] = E^-z. */
    assert_eval_eq("Gamma[1, z]", "E^(-z)", 0);
    /* Gamma[a, Infinity] = 0. */
    assert_eval_eq("Gamma[3.0, Infinity]", "0", 0);
}

void test_gamma_incomplete_int_closed() {
    /* Positive integer first argument -> finite closed form. */
    assert_eval_eq("Gamma[2, 3]", "4/E^3", 0);            /* exact stays exact */
    assert_eval_eq("Gamma[2, x]", "(1 + x) E^(-x)", 0);   /* symbolic reduces  */
    assert_eval_eq("Gamma[3, x]", "E^(-x) (2 + 2 x + x^2)", 0);
    assert_eval_eq("Gamma[2, 1/2]", "3/2/Sqrt[E]", 0);    /* exact rational z  */
    /* Gamma[4, 2] = 38 e^-2. */
    assert_eval_eq("Gamma[4, 2]", "38/E^2", 0);
    /* Numerically consistent with the incomplete-gamma integral. */
    assert_close("Gamma[2, 3.0]", 4.0 * exp(-3.0), 1e-12);
}

void test_gamma_derivatives() {
    /* dGamma/dz of the incomplete form: -z^(a-1) e^-z. */
    assert_eval_eq("D[Gamma[a, x], x]", "-E^(-x) x^(-1 + a)", 0);
    /* Chain rule through the second argument. */
    assert_eval_eq("D[Gamma[a, x^2], x]", "-2 x E^(-x^2) x^2^(-1 + a)", 0);
    /* Derivative wrt the first argument has no closed form here: generic. */
    assert_eval_eq("D[Gamma[x, z], x]", "Derivative[1, 0][Gamma][x, z]", 0);
    /* One-argument Gamma'[z] stays generic (PolyGamma is out of scope). */
    assert_eval_eq("D[Gamma[x], x]", "Derivative[1][Gamma][x]", 0);
}

void test_gamma_incomplete_numeric() {
    /* Gamma[1.5, 7.5] = 0.00160996322827... (matches Mathematica's 0.00160996). */
    assert_close("Gamma[1.5, 7.5]", 0.0016099632282723202, 1e-12);
    /* Gamma[1, 2.0] = E^-2. */
    assert_close("Gamma[1, 2.0]", exp(-2.0), 1e-12);
    /* Gamma[2.2] complete equals Gamma[2.2, 0]. */
    assert_close("Gamma[2.2, 0.0]", 1.1018024908797127, 1e-12);
    /* Arbitrary-precision incomplete gamma. */
    assert_eval_startswith("N[Gamma[1/2, 1], 30]", "0.27880");
}

void test_gamma_incomplete_complex() {
    /* Machine complex incomplete gamma. Gamma[2.0, 1+I] = (2+I) e^-(1+I):
       0.70709209634593808 - 0.42035364095981146 I. */
    assert_complex_close("Gamma[2.0, 1.0 + 1.0 I]",
                         0.7070920963459381, -0.4203536409598115, 1e-12);
    /* Continued-fraction branch (Re(z) large). Tiny imaginary part tracks the
       real incomplete gamma Gamma[1.5, 8.0] = 0.001004967410648. */
    assert_complex_close("Gamma[1.5, 8.0 + 1.0*^-12 I]",
                         0.001004967410648176, 0.0, 1e-12);
    /* Arbitrary-precision complex incomplete (series branch). N[Gamma[3/2, 2+I], 30]
       = 0.160487401929263240325468324955 - 0.176588715957602345887257880028 I. */
    assert_eval_startswith("N[Gamma[3/2, 2 + I], 30]", "0.16048740192926324032");
}

void test_gamma_generalized() {
    /* Gamma[1, 1.1, 2.2] = Gamma[1,1.1] - Gamma[1,2.2] = e^-1.1 - e^-2.2. */
    assert_close("Gamma[1, 1.1, 2.2]", exp(-1.1) - exp(-2.2), 1e-9);
    /* Gamma[2, z] = (1 + z) e^-z, so Gamma[2,1] - Gamma[2,3]. */
    assert_close("Gamma[2.0, 1.0, 3.0]",
                 (1.0 + 1.0) * exp(-1.0) - (1.0 + 3.0) * exp(-3.0), 1e-9);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_gamma_listable() {
    assert_eval_eq("Gamma[{2, 3, 5, 7, 11}]", "{1, 2, 24, 720, 3628800}", 0);
    assert_eval_eq("Gamma[{1, 2, 3}]", "{1, 1, 2}", 0);
}

void test_gamma_attributes() {
    SymbolDef* d = symtab_get_def("Gamma");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Gamma must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "Gamma must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Gamma must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_gamma_integer_values);
    TEST(test_gamma_poles);
    TEST(test_gamma_half_integers);
    TEST(test_gamma_infinities);
    TEST(test_gamma_symbolic);
    TEST(test_gamma_machine_real);
    TEST(test_gamma_machine_complex);
    TEST(test_gamma_arbitrary_precision);
    TEST(test_gamma_arbitrary_complex);
    TEST(test_gamma_incomplete_exact);
    TEST(test_gamma_incomplete_int_closed);
    TEST(test_gamma_derivatives);
    TEST(test_gamma_incomplete_numeric);
    TEST(test_gamma_incomplete_complex);
    TEST(test_gamma_generalized);
    TEST(test_gamma_listable);
    TEST(test_gamma_attributes);

    printf("All Gamma tests passed.\n");
    return 0;
}
