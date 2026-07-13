/* Tests for the Euler beta function and its incomplete / generalized forms.
 *
 *   Beta[a, b]         -- Euler beta B(a,b) = Gamma(a)Gamma(b)/Gamma(a+b)
 *   Beta[z, a, b]      -- incomplete beta  Int_0^z t^(a-1)(1-t)^(b-1) dt
 *   Beta[z0, z1, a, b] -- generalized incomplete = Beta[z1,a,b] - Beta[z0,a,b]
 *
 * Covers exact rational reductions (integer and half-integer arguments, the
 * positive-integer Pochhammer collapse), the surviving / cancelling pole
 * structure on the integer lattice, machine and arbitrary-precision (MPFR)
 * real and complex numerics, the incomplete and generalized forms (special
 * values, exact closed forms, numeric evaluation), Listable threading,
 * symbolic fall-through, all three derivative families, attributes, and the
 * arity diagnostic. */

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

/* ---- numeric helpers (mirrors test_gamma.c / test_pochhammer.c) ----- */

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

/* ---- exact 2-arg: integer arguments --------------------------------- */

void test_beta_integer() {
    assert_eval_eq("Beta[2, 3]", "1/12", 0);
    assert_eval_eq("Beta[5, 4]", "1/280", 0);
    assert_eval_eq("Beta[1, 3]", "1/3", 0);
    assert_eval_eq("Beta[5, 1]", "1/5", 0);
    /* Symmetry B(a,b) = B(b,a). */
    assert_eval_eq("Beta[4, 5]", "1/280", 0);
}

/* ---- exact 2-arg: half-integer arguments ---------------------------- */

void test_beta_half_integer() {
    assert_eval_eq("Beta[5/2, 7/2]", "3/256 Pi", 0);
    assert_eval_eq("Beta[1/2, 1/2]", "Pi", 0);          /* Gamma(1/2)^2/Gamma(1) */
}

/* ---- positive-integer Pochhammer collapse to a rational ------------- */

void test_beta_pochhammer_rational() {
    assert_eval_eq("Beta[3, 1/3]", "27/14", 0);
    assert_eval_eq("Beta[1/3, 3]", "27/14", 0);          /* symmetric */
    assert_eval_eq("Beta[7/2, 2]", "4/63", 0);
}

/* ---- Gamma-form (neither argument a positive integer) --------------- */

void test_beta_gamma_form() {
    assert_eval_eq("Beta[1/3, 1/3]", "Gamma[1/3]^2/Gamma[2/3]", 0);
    assert_eval_eq("Beta[2/3, 4/3]", "Gamma[2/3] Gamma[4/3]", 0);
}

/* ---- poles: a surviving gamma pole gives ComplexInfinity ------------ */

void test_beta_poles() {
    assert_eval_eq("Beta[0, b]", "ComplexInfinity", 0);
    assert_eval_eq("Beta[a, 0]", "ComplexInfinity", 0);
    assert_eval_eq("Beta[0, n]", "ComplexInfinity", 0);
    assert_eval_eq("Beta[-2, 5]", "ComplexInfinity", 0);
    assert_eval_eq("Beta[Infinity, 0]", "ComplexInfinity", 0);
    assert_eval_eq("Beta[0, Infinity]", "ComplexInfinity", 0);
}

/* ---- cancelling poles: finite limit of the gamma ratio -------------- */

void test_beta_pole_cancellation() {
    assert_eval_eq("Beta[-2, 1]", "-1/2", 0);
    assert_eval_eq("Beta[2, -5]", "1/20", 0);            /* symmetric form */
}

/* ---- machine-precision real numerics -------------------------------- */

void test_beta_machine_real() {
    assert_close("Beta[2.3, 3.2]", 0.0540298, 1e-6);
    assert_close("Beta[1.1, 2.1]", 0.410722, 1e-6);
    assert_close("Beta[2.0, 3.0]", 1.0 / 12.0, 1e-9);
}

/* ---- machine-precision complex numerics ----------------------------- */

void test_beta_machine_complex() {
    assert_complex_close("Beta[2.5 + I, 1 - I]", 0.0831078, 0.142164, 1e-5);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_beta_arbitrary_precision() {
    assert_eval_startswith("N[Beta[22/10, 33/10], 50]",
                           "0.0564856913732825668070517540044914293695377770152");
}

/* ---- Listable threading --------------------------------------------- */

void test_beta_listable() {
    assert_eval_eq("Beta[2, {2, 3, 4, 5}]", "{1/6, 1/12, 1/20, 1/30}", 0);
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_beta_symbolic() {
    assert_eval_eq("Beta[a, b]", "Beta[a, b]", 0);
}

/* ---- incomplete beta: special values & exact closed forms ----------- */

void test_beta_incomplete_special() {
    assert_eval_eq("Beta[0, a, b]", "0", 0);
    assert_eval_eq("Beta[1, a, b]", "Beta[a, b]", 0);
    assert_eval_eq("Beta[z, a, 1]", "z^a/a", 0);
    /* Positive-integer b terminates the 2F1 to an exact rational. */
    assert_eval_eq("Beta[1/2, 2, 3]", "11/192", 0);
    assert_eval_eq("Beta[1/2, 3, 4]", "7/640", 0);
    /* Fully symbolic incomplete stays unevaluated. */
    assert_eval_eq("Beta[z, a, b]", "Beta[z, a, b]", 0);
}

void test_beta_incomplete_numeric() {
    assert_close("Beta[0.5, 2.2, 3.3]", 0.0392297, 1e-6);
    assert_close("N[Beta[1/2, 3, 4]]", 7.0 / 640.0, 1e-9);
}

/* ---- generalized (4-arg) incomplete beta ---------------------------- */

void test_beta_generalized() {
    /* Beta[z0,z1,a,b] = Beta[z1,a,b] - Beta[z0,a,b]; check the numeric value. */
    assert_close("Beta[0.2, 0.6, 2, 3]", 0.0533333, 1e-6);
    /* Symbolic four-argument form stays unevaluated. */
    assert_eval_eq("Beta[z0, z1, a, b]", "Beta[z0, z1, a, b]", 0);
}

/* ---- derivatives ----------------------------------------------------- */

void test_beta_derivatives_two_arg() {
    assert_eval_eq("D[Beta[a, b], a]",
                   "Beta[a, b] (PolyGamma[0, a] - PolyGamma[0, a + b])", 0);
    assert_eval_eq("D[Beta[a, b], b]",
                   "Beta[a, b] (PolyGamma[0, b] - PolyGamma[0, a + b])", 0);
    assert_eval_eq("D[Beta[1/2, x], x]",
                   "Beta[1/2, x] (PolyGamma[0, x] - PolyGamma[0, 1/2 + x])", 0);
    /* Second derivative composes through the product rule and PolyGamma. */
    assert_eval_eq("D[Beta[a, b], {b, 2}]",
                   "Beta[a, b] (PolyGamma[1, b] - PolyGamma[1, a + b]) + "
                   "Beta[a, b] (PolyGamma[0, b] - PolyGamma[0, a + b])^2", 0);
}

void test_beta_derivatives_incomplete() {
    /* d/dz Beta[z,a,b] = z^(a-1) (1-z)^(b-1). */
    assert_eval_eq("D[Beta[z, a, b], z]", "z^(-1 + a) (1 - z)^(-1 + b)", 0);
    /* d/da, d/db are the generic inert partials. */
    assert_eval_eq("D[Beta[z, a, b], a]", "Derivative[0, 1, 0][Beta][z, a, b]", 0);
    /* generalized: d/dz1 and d/dz0 are the integrand at the endpoints. */
    assert_eval_eq("D[Beta[z0, z1, a, b], z1]",
                   "z1^(-1 + a) (1 - z1)^(-1 + b)", 0);
    assert_eval_eq("D[Beta[z0, z1, a, b], z0]",
                   "-z0^(-1 + a) (1 - z0)^(-1 + b)", 0);
}

/* ---- attributes ------------------------------------------------------ */

void test_beta_attributes() {
    SymbolDef* d = symtab_get_def("Beta");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Beta must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "Beta must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Beta must be Protected");
    ASSERT_MSG((d->attributes & ATTR_READPROTECTED) != 0,
               "Beta must be ReadProtected");
}

/* ---- arity diagnostics (stays unevaluated) -------------------------- */

void test_beta_arity() {
    assert_eval_eq("Beta[]", "Beta[]", 0);
    assert_eval_eq("Beta[1]", "Beta[1]", 0);
    assert_eval_eq("Beta[1, 2, 3, 4, 5]", "Beta[1, 2, 3, 4, 5]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_beta_integer);
    TEST(test_beta_half_integer);
    TEST(test_beta_pochhammer_rational);
    TEST(test_beta_gamma_form);
    TEST(test_beta_poles);
    TEST(test_beta_pole_cancellation);
    TEST(test_beta_machine_real);
    TEST(test_beta_machine_complex);
    TEST(test_beta_arbitrary_precision);
    TEST(test_beta_listable);
    TEST(test_beta_symbolic);
    TEST(test_beta_incomplete_special);
    TEST(test_beta_incomplete_numeric);
    TEST(test_beta_generalized);
    TEST(test_beta_derivatives_two_arg);
    TEST(test_beta_derivatives_incomplete);
    TEST(test_beta_attributes);
    TEST(test_beta_arity);

    printf("All Beta tests passed.\n");
    return 0;
}
