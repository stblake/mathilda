/* Tests for the sine integral SinIntegral[z] = Si(z) = Int_0^z Sin[t]/t dt.
 *
 * Covers exact special values (0, +-Infinity -> +-Pi/2, +-I Infinity,
 * ComplexInfinity/Indeterminate), machine real (convergent series and, for
 * large |x|, the asymptotic expansion), arbitrary-precision (MPFR) reals with
 * precision tracking, machine & arbitrary complex, odd symmetry, derivatives,
 * Taylor and asymptotic Series, Listable threading, attributes and arity
 * errors. Reference values cross-checked against the running binary and mpmath. */

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

void test_si_exact() {
    assert_eval_eq("SinIntegral[0]", "0", 0);
    assert_eval_eq("SinIntegral[Infinity]", "1/2 Pi", 0);
    assert_eval_eq("SinIntegral[-Infinity]", "-1/2 Pi", 0);
    assert_eval_eq("SinIntegral[I Infinity]", "I Infinity", 0);
    assert_eval_eq("SinIntegral[-I Infinity]", "-I Infinity", 0);
    assert_eval_eq("SinIntegral[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("SinIntegral[Indeterminate]", "Indeterminate", 0);
    assert_eval_eq("SinIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]",
                   "{-1/2 Pi, 1/2 Pi, -I Infinity, I Infinity}", 0);
}

/* ---- symbolic passthrough & odd symmetry ---------------------------- */

void test_si_symbolic() {
    assert_eval_eq("SinIntegral[x]", "SinIntegral[x]", 0);
    assert_eval_eq("SinIntegral[2]", "SinIntegral[2]", 0);
    assert_eval_eq("SinIntegral[a + b]", "SinIntegral[a + b]", 0);
    /* Odd symmetry: SinIntegral[-x] -> -SinIntegral[x]. */
    assert_eval_eq("SinIntegral[-x]", "-SinIntegral[x]", 0);
    assert_eval_eq("SinIntegral[-2 y]", "-SinIntegral[2 y]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_si_machine_real() {
    assert_close("SinIntegral[0.5]", 0.4931074180430667, 1e-12);
    assert_close("SinIntegral[1.0]", 0.9460830703671830, 1e-12);
    assert_close("SinIntegral[2.0]", 1.6054129768026949, 1e-12);
    assert_close("SinIntegral[2.8]", 1.8320965890813223, 1e-12);
    assert_close("SinIntegral[1.5]", 1.3246835311721197, 1e-12);
    assert_close("SinIntegral[3.5]", 1.8331253986659970, 1e-12);
    /* Odd: negative arguments. */
    assert_close("SinIntegral[-1.0]", -0.9460830703671830, 1e-12);
    assert_close("SinIntegral[-2.0]", -1.6054129768026949, 1e-12);
    /* Large |x|: exercises the asymptotic path (approaches Pi/2). */
    assert_close("SinIntegral[50.0]", 1.5516170724859359, 1e-10);
    assert_close("SinIntegral[100.0]", 1.5622254668890563, 1e-10);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_si_arbitrary_precision() {
    /* N[SinIntegral[2], 50] -- the task-spec reference value. */
    assert_eval_startswith("N[SinIntegral[2], 50]",
        "1.6054129768026948485767201481985889408485834223285");
    assert_eval_startswith("N[SinIntegral[1], 40]",
        "0.94608307036718301494135331382317965781");
    assert_eval_startswith("N[SinIntegral[1/2], 40]",
        "0.49310741804306668916162670757276465364");
    /* Large argument at high precision (asymptotic path). */
    assert_eval_startswith("N[SinIntegral[50], 30]",
        "1.55161707248593589472798559485");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("SinIntegral[2.0000000000000000000000]",
                           "1.6054129768026948485767");
}

/* ---- machine complex ------------------------------------------------ */

void test_si_machine_complex() {
    assert_complex_close("SinIntegral[2.5 + I]", 1.9954872853802745, 0.2229952966090301, 1e-10);
    /* Si(3 I) = I Shi(3). */
    assert_complex_close("SinIntegral[3.0 I]", 0.0, 4.973440475859807, 1e-9);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_si_arbitrary_complex() {
    assert_eval_startswith("N[SinIntegral[2.5 + I], 30]",
        "1.99548728538027453360825802519");
}

/* ---- derivatives ---------------------------------------------------- */

void test_si_derivatives() {
    assert_eval_eq("D[SinIntegral[x], x]", "Sinc[x]", 0);
    /* Chain rule. */
    assert_eval_eq("D[SinIntegral[x^2], x]", "2 x Sinc[x^2]", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_si_series_at_zero() {
    assert_eval_eq("Series[SinIntegral[x], {x, 0, 7}]",
        "x - 1/18 x^3 + 1/600 x^5 - 1/35280 x^7 + O[x]^8", 0);
    /* Composed argument. */
    assert_eval_eq("Series[SinIntegral[x^2], {x, 0, 7}]",
        "x^2 - 1/18 x^6 + O[x]^8", 0);
}

void test_si_series_at_infinity() {
    /* Normal form matches Pi/2 - Sin[x]/x^2 + Cos[x] (-1/x + 2/x^3). */
    assert_eval_eq("Normal[Series[SinIntegral[x], {x, Infinity, 3}]]",
                   "1/2 Pi - Sin[x]/x^2 + Cos[x] (-1/x + 2/x^3)", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_si_listable() {
    assert_eval_eq("SinIntegral[{}]", "{}", 0);
    assert_eval_eq("SinIntegral[{0}]", "{0}", 0);
}

void test_si_listable_numeric() {
    assert_close("SinIntegral[{1.5, 2.5, 3.5}][[1]]", 1.3246835311721197, 1e-10);
    assert_close("SinIntegral[{1.5, 2.5, 3.5}][[2]]", 1.7785201734438266, 1e-10);
    assert_close("SinIntegral[{1.5, 2.5, 3.5}][[3]]", 1.8331253986659970, 1e-10);
}

void test_si_attributes() {
    SymbolDef* d = symtab_get_def("SinIntegral");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "SinIntegral must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "SinIntegral must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "SinIntegral must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_si_arity() {
    assert_eval_eq("SinIntegral[]", "SinIntegral[]", 0);
    assert_eval_eq("SinIntegral[1, 2, 3]", "SinIntegral[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_si_exact);
    TEST(test_si_symbolic);
    TEST(test_si_machine_real);
    TEST(test_si_arbitrary_precision);
    TEST(test_si_machine_complex);
    TEST(test_si_arbitrary_complex);
    TEST(test_si_derivatives);
    TEST(test_si_series_at_zero);
    TEST(test_si_series_at_infinity);
    TEST(test_si_listable);
    TEST(test_si_listable_numeric);
    TEST(test_si_attributes);
    TEST(test_si_arity);

    printf("All SinIntegral tests passed.\n");
    return 0;
}
