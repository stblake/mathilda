/* Tests for the hyperbolic cosine integral
 * CoshIntegral[z] = Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t]-1)/t dt.
 *
 * Covers exact special values (0 -> -Infinity, Infinity, +-I Infinity,
 * ComplexInfinity/Indeterminate), machine real (convergent series and, for
 * large |x|, the asymptotic expansion, including the double-overflow -> MPFR
 * fallback), the negative-real from-above branch value Complex[Chi(|x|), Pi],
 * arbitrary-precision (MPFR) reals with precision tracking, machine & arbitrary
 * complex, the imaginary-axis Stokes value, derivatives, log-series and
 * asymptotic Series, Integrate emission, Listable threading, attributes and
 * arity errors. Reference values cross-checked against the running binary and
 * the exact identity Chi(z) = Ci(i z) - i Pi/2. */

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

void test_chi_exact() {
    assert_eval_eq("CoshIntegral[0]", "-Infinity", 0);
    assert_eval_eq("CoshIntegral[Infinity]", "Infinity", 0);
    assert_eval_eq("CoshIntegral[-Infinity]", "Infinity", 0);
    assert_eval_eq("CoshIntegral[I Infinity]", "(1/2*I) Pi", 0);
    assert_eval_eq("CoshIntegral[-I Infinity]", "(-1/2*I) Pi", 0);
    assert_eval_eq("CoshIntegral[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("CoshIntegral[Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("CoshIntegral[{Infinity, -I Infinity, I Infinity}]",
                   "{Infinity, (-1/2*I) Pi, (1/2*I) Pi}", 0);
}

/* ---- symbolic passthrough ------------------------------------------- */

void test_chi_symbolic() {
    assert_eval_eq("CoshIntegral[x]", "CoshIntegral[x]", 0);
    assert_eval_eq("CoshIntegral[2]", "CoshIntegral[2]", 0);
    assert_eval_eq("CoshIntegral[a + b]", "CoshIntegral[a + b]", 0);
    /* Not odd, not even: no symmetry fold. */
    assert_eval_eq("CoshIntegral[-x]", "CoshIntegral[-x]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_chi_machine_real() {
    assert_close("CoshIntegral[0.5]", -0.0527768449564936168, 1e-12);
    assert_close("CoshIntegral[1.0]", 0.837866940980208197, 1e-12);
    assert_close("CoshIntegral[2.0]", 2.45266692264691466, 1e-12);
    assert_close("CoshIntegral[2.8]", 4.33122121568197382, 1e-11);
    assert_close("CoshIntegral[1.5]", 1.60063293336158252, 1e-12);
    assert_close("CoshIntegral[3.5]", 6.95919192764739325, 1e-11);
    /* Large |x|: exercises the asymptotic path (grows like e^x/(2x)). */
    assert_close("CoshIntegral[20.0]", 12807826.332028294, 1e-1);
}

/* ---- negative real: from-above branch value Complex[Chi(|x|), Pi] --- */

void test_chi_negative_real() {
    assert_complex_close("CoshIntegral[-2.0]", 2.45266692264691466, M_PI, 1e-11);
    assert_complex_close("CoshIntegral[-1.0]", 0.837866940980208197, M_PI, 1e-11);
}

/* ---- overflow -> MPFR real ------------------------------------------ */

void test_chi_overflow() {
    /* A machine real whose result overflows a C double (~1.5e434288) must come
     * back as a finite (MPFR) real, not Infinity. */
    assert_eval_startswith("CoshIntegral[10.^6]", "1.5166092150117");
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_chi_arbitrary_precision() {
    /* N[CoshIntegral[2], 50] -- the task-spec reference value. */
    assert_eval_startswith("N[CoshIntegral[2], 50]",
        "2.4526669226469145219061326474994928766017806887285");
    assert_eval_startswith("N[CoshIntegral[1], 40]",
        "0.83786694098020824089467857943575630999303");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("CoshIntegral[2.000000000000000000000000000000000000]",
                           "2.4526669226469145219061326474994928");
}

/* ---- machine complex ------------------------------------------------ */

void test_chi_machine_complex() {
    assert_complex_close("CoshIntegral[3.5 + I]",
                         5.366675935843618, 4.284227027938926, 1e-9);
    /* Chi(3 I) = Ci(3) + I Pi/2. */
    assert_complex_close("CoshIntegral[3.0 I]",
                         0.119629786008000, 1.5707963267948966, 1e-9);
    /* Large |z| off-axis: asymptotic path, upper-left quadrant. */
    assert_complex_close("CoshIntegral[-80.0 + 60.0 I]",
                         -2.63289163945237e32, -9.31083310328351e31, 1e18);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_chi_arbitrary_complex() {
    /* Exact input (7/2 = 3.5 exactly) so N[..., 30] genuinely exercises the
     * MPFR complex path: a machine 3.5 + I is MachinePrecision and N cannot
     * manufacture 30 digits from it (matches WL). */
    assert_eval_startswith("N[CoshIntegral[7/2 + I], 30]",
        "5.36667593584361750026705329022");
    /* Imaginary-axis Stokes constant: Chi(i y) -> i Pi/2, residual is Ci(y). */
    assert_eval_startswith("N[CoshIntegral[I*10^60`120] - I Pi/2, 20]",
        "8.3038976521934266466");
}

/* ---- derivatives ---------------------------------------------------- */

void test_chi_derivatives() {
    assert_eval_eq("D[CoshIntegral[x], x]", "Cosh[x]/x", 0);
    /* Chain rule. */
    assert_eval_eq("D[CoshIntegral[x^2], x]", "(2 Cosh[x^2])/x", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_chi_series_at_zero() {
    assert_eval_eq("Series[CoshIntegral[x], {x, 0, 6}]",
        "EulerGamma + Log[x] + 1/4 x^2 + 1/96 x^4 + 1/4320 x^6 + O[x]^7", 0);
    assert_eval_eq("Series[CoshIntegral[x], {x, 0, 5}]",
        "EulerGamma + Log[x] + 1/4 x^2 + 1/96 x^4 + O[x]^6", 0);
}

void test_chi_series_at_infinity() {
    assert_eval_eq("Normal[Series[CoshIntegral[x], {x, Infinity, 3}]]",
                   "(-1/2*I) Pi + Cosh[x]/x^2 + Sinh[x] (1/x + 2/x^3)", 0);
}

/* ---- Integrate emission --------------------------------------------- */

void test_chi_integrate() {
    assert_eval_eq("Integrate[Cosh[x]/x, x]", "CoshIntegral[x]", 0);
    assert_eval_eq("Integrate[Cosh[3 x]/x, x]", "CoshIntegral[3 x]", 0);
    /* Diff-back certifies the emitted antiderivative. */
    assert_eval_eq("Simplify[D[Integrate[Cosh[a x]/x, x], x] - Cosh[a x]/x]", "0", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_chi_listable() {
    assert_eval_eq("CoshIntegral[{}]", "{}", 0);
}

void test_chi_listable_numeric() {
    assert_close("CoshIntegral[{1.2, 1.5, 1.8}][[1]]", 1.14184192417059438, 1e-10);
    assert_close("CoshIntegral[{1.2, 1.5, 1.8}][[2]]", 1.60063293336158252, 1e-10);
    assert_close("CoshIntegral[{1.2, 1.5, 1.8}][[3]]", 2.09257721406203245, 1e-10);
}

void test_chi_attributes() {
    SymbolDef* d = symtab_get_def("CoshIntegral");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "CoshIntegral must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "CoshIntegral must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "CoshIntegral must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_chi_arity() {
    assert_eval_eq("CoshIntegral[]", "CoshIntegral[]", 0);
    assert_eval_eq("CoshIntegral[1, 2, 3]", "CoshIntegral[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_chi_exact);
    TEST(test_chi_symbolic);
    TEST(test_chi_machine_real);
    TEST(test_chi_negative_real);
    TEST(test_chi_overflow);
    TEST(test_chi_arbitrary_precision);
    TEST(test_chi_machine_complex);
    TEST(test_chi_arbitrary_complex);
    TEST(test_chi_derivatives);
    TEST(test_chi_series_at_zero);
    TEST(test_chi_series_at_infinity);
    TEST(test_chi_integrate);
    TEST(test_chi_listable);
    TEST(test_chi_listable_numeric);
    TEST(test_chi_attributes);
    TEST(test_chi_arity);

    printf("All CoshIntegral tests passed.\n");
    return 0;
}
