/* Tests for the hyperbolic sine integral
 * SinhIntegral[z] = Shi(z) = Int_0^z Sinh[t]/t dt.
 *
 * Covers exact special values (0, +-Infinity, +-I Infinity,
 * ComplexInfinity/Indeterminate), machine real (convergent series and, for
 * large |x|, the asymptotic expansion, including the double-overflow -> MPFR
 * fallback), arbitrary-precision (MPFR) reals with precision tracking, machine
 * & arbitrary complex, the imaginary-axis Stokes value, odd symmetry,
 * derivatives, Taylor and asymptotic Series, Integrate emission, Listable
 * threading, attributes and arity errors. Reference values cross-checked
 * against the running binary and the exact identity Shi(z) = -i Si(i z). */

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

void test_shi_exact() {
    assert_eval_eq("SinhIntegral[0]", "0", 0);
    assert_eval_eq("SinhIntegral[Infinity]", "Infinity", 0);
    assert_eval_eq("SinhIntegral[-Infinity]", "-Infinity", 0);
    assert_eval_eq("SinhIntegral[I Infinity]", "(1/2*I) Pi", 0);
    assert_eval_eq("SinhIntegral[-I Infinity]", "(-1/2*I) Pi", 0);
    assert_eval_eq("SinhIntegral[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("SinhIntegral[Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("SinhIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]",
                   "{-Infinity, Infinity, (-1/2*I) Pi, (1/2*I) Pi}", 0);
}

/* ---- symbolic passthrough & odd symmetry ---------------------------- */

void test_shi_symbolic() {
    assert_eval_eq("SinhIntegral[x]", "SinhIntegral[x]", 0);
    assert_eval_eq("SinhIntegral[2]", "SinhIntegral[2]", 0);
    assert_eval_eq("SinhIntegral[a + b]", "SinhIntegral[a + b]", 0);
    /* Odd symmetry: SinhIntegral[-x] -> -SinhIntegral[x]. */
    assert_eval_eq("SinhIntegral[-x]", "-SinhIntegral[x]", 0);
    assert_eval_eq("SinhIntegral[-2 y]", "-SinhIntegral[2 y]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_shi_machine_real() {
    assert_close("SinhIntegral[0.5]", 0.506996749819667247, 1e-12);
    assert_close("SinhIntegral[1.0]", 1.05725087537572859, 1e-12);
    assert_close("SinhIntegral[2.0]", 2.50156743335497556, 1e-12);
    assert_close("SinhIntegral[2.8]", 4.34807650812718993, 1e-11);
    assert_close("SinhIntegral[1.5]", 1.7006525157682153, 1e-12);
    assert_close("SinhIntegral[3.5]", 6.96616206750494182, 1e-11);
    /* Odd: negative arguments. */
    assert_close("SinhIntegral[-1.0]", -1.05725087537572859, 1e-12);
    assert_close("SinhIntegral[-2.0]", -2.50156743335497556, 1e-12);
    /* Large |x|: exercises the asymptotic path (grows like e^x/(2x)). */
    assert_close("SinhIntegral[20.0]", 12807826.332028294, 1e-1);
}

/* ---- overflow -> MPFR real ------------------------------------------ */

void test_shi_overflow() {
    /* A machine real whose result overflows a C double (~1.5e434288) must come
     * back as a finite (MPFR) real, not Infinity. */
    assert_eval_startswith("SinhIntegral[10.^6]", "1.5166092150117");
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_shi_arbitrary_precision() {
    /* N[SinhIntegral[2], 50] -- the task-spec reference value. */
    assert_eval_startswith("N[SinhIntegral[2], 50]",
        "2.5015674333549756414733724827275423989162728736915");
    assert_eval_startswith("N[SinhIntegral[1], 40]",
        "1.057250875375728514571842354895877959024");
    /* Large argument at high precision (asymptotic path). */
    assert_eval_startswith("N[SinhIntegral[20], 30]",
        "12807826.3320282944594181868553");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("SinhIntegral[2.0000000000000000000000]",
                           "2.501567433354975641473");
}

/* ---- machine complex ------------------------------------------------ */

void test_shi_machine_complex() {
    assert_complex_close("SinhIntegral[2.5 + I]",
                         2.846494698401862, 2.221765034323259, 1e-10);
    /* Shi(3 I) = I Si(3). */
    assert_complex_close("SinhIntegral[3.0 I]", 0.0, 1.848652527999468, 1e-9);
    /* Large |z| off-axis: asymptotic path, upper-left quadrant. */
    assert_complex_close("SinhIntegral[-80.0 + 60.0 I]",
                         2.63289163945237e32, 9.31083310328351e31, 1e18);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_shi_arbitrary_complex() {
    assert_eval_startswith("N[SinhIntegral[2.5 + I], 30]",
        "2.846494698401862422088015591726");
    /* Imaginary-axis Stokes constant: Shi(i y) -> i Pi/2, residual is i Si-tail.
     * With an explicit-precision argument the ~1e-61 residual survives. */
    assert_eval_startswith("N[SinhIntegral[I*10^60`120] - I Pi/2, 20]",
        "0.0 + 5.5718294824856670897");
}

/* ---- derivatives ---------------------------------------------------- */

void test_shi_derivatives() {
    assert_eval_eq("D[SinhIntegral[x], x]", "Sinh[x]/x", 0);
    /* Chain rule. */
    assert_eval_eq("D[SinhIntegral[x^2], x]", "(2 Sinh[x^2])/x", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_shi_series_at_zero() {
    assert_eval_eq("Series[SinhIntegral[x], {x, 0, 7}]",
        "x + 1/18 x^3 + 1/600 x^5 + 1/35280 x^7 + O[x]^8", 0);
    /* Composed argument. */
    assert_eval_eq("Series[SinhIntegral[x^2], {x, 0, 7}]",
        "x^2 + 1/18 x^6 + O[x]^8", 0);
}

void test_shi_series_at_infinity() {
    assert_eval_eq("Normal[Series[SinhIntegral[x], {x, Infinity, 3}]]",
                   "(-1/2*I) Pi + Sinh[x]/x^2 + Cosh[x] (1/x + 2/x^3)", 0);
}

/* ---- Integrate emission --------------------------------------------- */

void test_shi_integrate() {
    assert_eval_eq("Integrate[Sinh[x]/x, x]", "SinhIntegral[x]", 0);
    assert_eval_eq("Integrate[Sinh[2 x]/x, x]", "SinhIntegral[2 x]", 0);
    /* Diff-back certifies the emitted antiderivative. */
    assert_eval_eq("Simplify[D[Integrate[Sinh[a x]/x, x], x] - Sinh[a x]/x]", "0", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_shi_listable() {
    assert_eval_eq("SinhIntegral[{}]", "{}", 0);
    assert_eval_eq("SinhIntegral[{0}]", "{0}", 0);
}

void test_shi_listable_numeric() {
    assert_close("SinhIntegral[{1.2, 1.5, 1.8}][[1]]", 1.30025036102205704, 1e-10);
    assert_close("SinhIntegral[{1.2, 1.5, 1.8}][[2]]", 1.7006525157682153, 1e-10);
    assert_close("SinhIntegral[{1.2, 1.5, 1.8}][[3]]", 2.1572903434259012, 1e-10);
}

void test_shi_attributes() {
    SymbolDef* d = symtab_get_def("SinhIntegral");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "SinhIntegral must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "SinhIntegral must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "SinhIntegral must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_shi_arity() {
    assert_eval_eq("SinhIntegral[]", "SinhIntegral[]", 0);
    assert_eval_eq("SinhIntegral[1, 2, 3]", "SinhIntegral[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_shi_exact);
    TEST(test_shi_symbolic);
    TEST(test_shi_machine_real);
    TEST(test_shi_overflow);
    TEST(test_shi_arbitrary_precision);
    TEST(test_shi_machine_complex);
    TEST(test_shi_arbitrary_complex);
    TEST(test_shi_derivatives);
    TEST(test_shi_series_at_zero);
    TEST(test_shi_series_at_infinity);
    TEST(test_shi_integrate);
    TEST(test_shi_listable);
    TEST(test_shi_listable_numeric);
    TEST(test_shi_attributes);
    TEST(test_shi_arity);

    printf("All SinhIntegral tests passed.\n");
    return 0;
}
