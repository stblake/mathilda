/* Tests for the cosine integral CosIntegral[z] = Ci(z) = -Int_z^Inf Cos[t]/t dt.
 *
 * Covers exact special values (0 -> -Infinity, Infinity -> 0, -Infinity -> I Pi,
 * +-I Infinity -> Infinity, ComplexInfinity/Indeterminate), machine real
 * (convergent series and, for large |x|, the asymptotic expansion), the
 * negative-real branch cut (Ci(-x) = Ci(x) + I Pi from above), arbitrary-
 * precision (MPFR) reals with precision tracking, machine & arbitrary complex,
 * the pure-imaginary axis (Ci(I y) = Chi(y) + I Pi/2), derivatives, Taylor and
 * asymptotic Series, Listable threading, attributes and arity errors. Reference
 * values cross-checked against the running binary and the Ei relation
 * Ci(z) = (1/2)(Ei(I z) + Ei(-I z)). */

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

void test_ci_exact() {
    assert_eval_eq("CosIntegral[0]", "-Infinity", 0);
    assert_eval_eq("CosIntegral[Infinity]", "0", 0);
    assert_eval_eq("CosIntegral[-Infinity]", "I Pi", 0);
    assert_eval_eq("CosIntegral[I Infinity]", "Infinity", 0);
    assert_eval_eq("CosIntegral[-I Infinity]", "Infinity", 0);
    assert_eval_eq("CosIntegral[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("CosIntegral[Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("CosIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]",
                   "{I Pi, 0, Infinity, Infinity}", 0);
}

/* ---- symbolic passthrough (no odd symmetry: Ci is not odd) ---------- */

void test_ci_symbolic() {
    assert_eval_eq("CosIntegral[x]", "CosIntegral[x]", 0);
    assert_eval_eq("CosIntegral[2]", "CosIntegral[2]", 0);
    assert_eval_eq("CosIntegral[a + b]", "CosIntegral[a + b]", 0);
    /* No fold: CosIntegral[-x] stays symbolic (unlike SinIntegral). */
    assert_eval_eq("CosIntegral[-x]", "CosIntegral[-x]", 0);
}

/* ---- machine real (positive) ---------------------------------------- */

void test_ci_machine_real() {
    assert_close("CosIntegral[0.5]", -0.17778407880661290, 1e-12);
    assert_close("CosIntegral[1.0]",  0.33740392290096813, 1e-12);
    assert_close("CosIntegral[2.0]",  0.42298082877486500, 1e-12);
    assert_close("CosIntegral[2.5]",  0.28587119636538350, 1e-12);
    assert_close("CosIntegral[2.8]",  0.18648838964317577, 1e-12);
    assert_close("CosIntegral[1.5]",  0.47035631719539989, 1e-12);
    assert_close("CosIntegral[3.5]", -0.03212854851248111, 1e-12);
    /* Large |x|: exercises the asymptotic path (decays toward 0). */
    assert_close("CosIntegral[50.0]",  -0.0056283863241163054, 1e-10);
    assert_close("CosIntegral[100.0]", -0.0051488251426104921, 1e-10);
}

/* ---- negative real: branch cut, Ci(-x) = Ci(x) + I Pi (from above) -- */

void test_ci_negative_real() {
    assert_complex_close("CosIntegral[-2.0]", 0.42298082877486500, M_PI, 1e-12);
    assert_complex_close("CosIntegral[-1.0]", 0.33740392290096813, M_PI, 1e-12);
    assert_complex_close("CosIntegral[-50.0]", -0.0056283863241163054, M_PI, 1e-10);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_ci_arbitrary_precision() {
    /* N[CosIntegral[2], 50] -- the task-spec reference value. */
    assert_eval_startswith("N[CosIntegral[2], 50]",
        "0.42298082877486499569856515319825589413573775630");
    assert_eval_startswith("N[CosIntegral[1], 40]",
        "0.33740392290096813466264620388915076999");
    assert_eval_startswith("N[CosIntegral[1/2], 40]",
        "-0.17778407880661290133581027107056907809");
    /* Large argument at high precision (asymptotic path). */
    assert_eval_startswith("N[CosIntegral[50], 30]",
        "-0.0056283863241163054401858954984");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("CosIntegral[2.0000000000000000000000]",
                           "0.42298082877486499569856");
    assert_eval_startswith("CosIntegral[2.8`30]",
                           "0.18648838964317576774819020");
}

/* ---- machine complex ------------------------------------------------ */

void test_ci_machine_complex() {
    assert_complex_close("CosIntegral[2.5 + I]",
                         0.33146500990509870, -0.38823676733822632, 1e-10);
    assert_complex_close("CosIntegral[2.5 + 1.5 I]",
                         0.36255769861805063, -0.70927876348396112, 1e-6);
}

/* ---- pure imaginary: Ci(I y) = Chi(y) + I Pi/2 ---------------------- */

void test_ci_pure_imaginary() {
    /* Real part = Chi(3) = 4.960392094765609760...; imaginary part = Pi/2. */
    assert_complex_close("CosIntegral[3.0 I]", 4.9603920947656098, M_PI_2, 1e-9);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_ci_arbitrary_complex() {
    assert_eval_startswith("N[CosIntegral[2.5 + I], 30]",
        "0.331465009905098695330849523088");
    /* Pure imaginary at high precision: imaginary part is exactly Pi/2. */
    assert_eval_startswith("N[CosIntegral[3 I], 30]",
        "4.96039209476560976029791763669");
}

/* ---- derivatives ---------------------------------------------------- */

void test_ci_derivatives() {
    assert_eval_eq("D[CosIntegral[x], x]", "Cos[x]/x", 0);
    /* Chain rule: d/dx Ci(x^2) = 2 x Cos[x^2]/x^2 = 2 Cos[x^2]/x. */
    assert_eval_eq("D[CosIntegral[x^2], x]", "(2 Cos[x^2])/x", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_ci_series_at_zero() {
    assert_eval_eq("Series[CosIntegral[x], {x, 0, 6}]",
        "EulerGamma + Log[x] - 1/4 x^2 + 1/96 x^4 - 1/4320 x^6 + O[x]^7", 0);
}

void test_ci_series_at_infinity() {
    assert_eval_eq("Normal[Series[CosIntegral[x], {x, Infinity, 3}]]",
                   "-Cos[x]/x^2 + Sin[x] (1/x - 2/x^3)", 0);
}

void test_ci_series_at_point() {
    /* Generic-point Taylor via the Cos[x]/x derivative rule. */
    assert_eval_eq("Series[CosIntegral[x], {x, Pi, 2}]",
        "CosIntegral[Pi] + -1/Pi (x - Pi) + 1/2/Pi^2 (x - Pi)^2 + O[x - Pi]^3", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_ci_listable() {
    assert_eval_eq("CosIntegral[{}]", "{}", 0);
    assert_eval_eq("CosIntegral[{0}]", "{-Infinity}", 0);
}

void test_ci_listable_numeric() {
    assert_close("CosIntegral[{1.5, 2.5, 3.5}][[1]]",  0.47035631719539989, 1e-10);
    assert_close("CosIntegral[{1.5, 2.5, 3.5}][[2]]",  0.28587119636538350, 1e-10);
    assert_close("CosIntegral[{1.5, 2.5, 3.5}][[3]]", -0.03212854851248111, 1e-10);
}

void test_ci_attributes() {
    SymbolDef* d = symtab_get_def("CosIntegral");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "CosIntegral must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "CosIntegral must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "CosIntegral must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_ci_arity() {
    assert_eval_eq("CosIntegral[]", "CosIntegral[]", 0);
    assert_eval_eq("CosIntegral[1, 2, 3]", "CosIntegral[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ci_exact);
    TEST(test_ci_symbolic);
    TEST(test_ci_machine_real);
    TEST(test_ci_negative_real);
    TEST(test_ci_arbitrary_precision);
    TEST(test_ci_machine_complex);
    TEST(test_ci_pure_imaginary);
    TEST(test_ci_arbitrary_complex);
    TEST(test_ci_derivatives);
    TEST(test_ci_series_at_zero);
    TEST(test_ci_series_at_infinity);
    TEST(test_ci_series_at_point);
    TEST(test_ci_listable);
    TEST(test_ci_listable_numeric);
    TEST(test_ci_attributes);
    TEST(test_ci_arity);

    printf("All CosIntegral tests passed.\n");
    return 0;
}
