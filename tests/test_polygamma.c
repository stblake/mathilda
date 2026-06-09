/* Tests for the PolyGamma function family: PolyGamma[z] and PolyGamma[n, z].
 *
 * Covers canonicalisation to the two-argument form, argument-count errors,
 * special values (poles, infinities), exact closed forms at positive integers
 * (digamma -> rational - EulerGamma; odd polygamma -> rational + rational*Pi^k;
 * even polygamma stays symbolic), the LogGamma reduction at order -1, machine and
 * arbitrary-precision (MPFR) real numerics with precision tracking, complex
 * numerics, derivative rules, Listable threading, and attributes. */

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

/* ---- canonicalisation & argument count ------------------------------ */

void test_polygamma_canonical() {
    /* The one-argument form rewrites to PolyGamma[0, z]. */
    assert_eval_eq("PolyGamma[z]", "PolyGamma[0, z]", 0);
    /* ... and then evaluates further for special arguments. */
    assert_eval_eq("PolyGamma[5]", "25/12 - EulerGamma", 0);
}

void test_polygamma_argt() {
    /* Wrong argument counts stay unevaluated (a ::argt message is printed). */
    assert_eval_eq("PolyGamma[]", "PolyGamma[]", 0);
    assert_eval_eq("PolyGamma[1, 3, 8, 9]", "PolyGamma[1, 3, 8, 9]", 0);
}

/* ---- special values ------------------------------------------------- */

void test_polygamma_poles() {
    /* psi^(n) has poles at every non-positive integer. */
    assert_eval_eq("PolyGamma[0]", "ComplexInfinity", 0);
    assert_eval_eq("PolyGamma[-2]", "ComplexInfinity", 0);
    assert_eval_eq("PolyGamma[0, 0]", "ComplexInfinity", 0);
    assert_eval_eq("PolyGamma[1, -3]", "ComplexInfinity", 0);
    assert_eval_eq("PolyGamma[2, -10]", "ComplexInfinity", 0);
}

void test_polygamma_infinities() {
    assert_eval_eq("PolyGamma[Infinity]", "Infinity", 0);
    assert_eval_eq("PolyGamma[0, Infinity]", "Infinity", 0);
    assert_eval_eq("PolyGamma[1, Infinity]", "0", 0);
    assert_eval_eq("PolyGamma[2, Infinity]", "0", 0);
    assert_eval_eq("PolyGamma[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("PolyGamma[3, Indeterminate]", "Indeterminate", 0);
}

/* ---- exact closed forms -------------------------------------------- */

void test_polygamma_digamma_integers() {
    /* psi(m) = H_{m-1} - EulerGamma. */
    assert_eval_eq("PolyGamma[1]", "-EulerGamma", 0);
    assert_eval_eq("PolyGamma[2]", "1 - EulerGamma", 0);
    assert_eval_eq("PolyGamma[0, 3]", "3/2 - EulerGamma", 0);
    assert_eval_eq("PolyGamma[4]", "11/6 - EulerGamma", 0);
    assert_eval_eq("PolyGamma[{1, 2, 3, 4, 5}]",
                   "{-EulerGamma, 1 - EulerGamma, 3/2 - EulerGamma, "
                   "11/6 - EulerGamma, 25/12 - EulerGamma}", 0);
}

void test_polygamma_odd_order_exact() {
    /* psi^(n)(m) for odd n closes in terms of Pi^(n+1). */
    assert_eval_eq("PolyGamma[1, 1]", "1/6 Pi^2", 0);     /* zeta(2) = Pi^2/6  */
    assert_eval_eq("PolyGamma[3, 1]", "1/15 Pi^4", 0);    /* 6 zeta(4) = Pi^4/15 */
    assert_eval_eq("PolyGamma[1, 2]", "-1 + 1/6 Pi^2", 0);
    /* PolyGamma[3, 5] = 6(-22369/20736 + Pi^4/90) = -22369/3456 + Pi^4/15. */
    assert_eval_eq("PolyGamma[3, 5]", "-22369/3456 + 1/15 Pi^4", 0);
}

void test_polygamma_even_order_symbolic() {
    /* Even order at an integer involves zeta(odd): no closed form, stays put. */
    assert_eval_eq("PolyGamma[2, 1]", "PolyGamma[2, 1]", 0);
    assert_eval_eq("PolyGamma[4, 1]", "PolyGamma[4, 1]", 0);
    assert_eval_eq("Table[PolyGamma[n, 1], {n, 0, 4}]",
                   "{-EulerGamma, 1/6 Pi^2, PolyGamma[2, 1], 1/15 Pi^4, "
                   "PolyGamma[4, 1]}", 0);
}

void test_polygamma_loggamma() {
    /* PolyGamma[-1, z] = LogGamma[z], which now evaluates: symbolic stays
     * LogGamma[z], exact integers reduce through LogGamma to Log[(n-1)!]. */
    assert_eval_eq("PolyGamma[-1, z]", "LogGamma[z]", 0);
    assert_eval_eq("PolyGamma[-1, 5]", "Log[24]", 0);
    /* Orders <= -2 stay symbolic. */
    assert_eval_eq("PolyGamma[-2, z]", "PolyGamma[-2, z]", 0);
}

void test_polygamma_symbolic() {
    /* Symbolic / non-integer-order arguments stay unevaluated. */
    assert_eval_eq("PolyGamma[0, x]", "PolyGamma[0, x]", 0);
    assert_eval_eq("PolyGamma[n, x]", "PolyGamma[n, x]", 0);
    /* Exact non-integer rational argument has no closed form here. */
    assert_eval_eq("PolyGamma[0, 5/2]", "PolyGamma[0, 5/2]", 0);
}

/* ---- machine-precision numerics ------------------------------------ */

void test_polygamma_machine_real() {
    assert_close("PolyGamma[0, 2.5]", 0.7031566406452432, 1e-12);   /* digamma(2.5) */
    assert_close("PolyGamma[100.5]", 4.6051743525818454, 1e-9);     /* digamma(100.5) */
    assert_close("PolyGamma[1, 3.5]", 0.3303577561002349, 1e-12);   /* trigamma(3.5) */
    assert_close("PolyGamma[1, 2.0]", 0.6449340668482264, 1e-12);   /* Pi^2/6 - 1 */
    assert_close("PolyGamma[2, 1.0]", -2.4041138063191885, 1e-9);   /* -2 zeta(3) */
}

/* ---- arbitrary-precision numerics ---------------------------------- */

void test_polygamma_mpfr_real() {
    /* N[PolyGamma[22/10], 50] -- digamma(2.2) to 50 digits. */
    assert_eval_startswith("N[PolyGamma[22/10], 50]",
                           "0.5442934367411450377861253708833812285077450591266");
    /* The two-argument numeric path at high precision: trigamma(5/2). */
    assert_eval_startswith("N[PolyGamma[1, 5/2], 40]",
                           "0.49035775610023486497280105549363");
}

/* ---- complex numerics ---------------------------------------------- */

void test_polygamma_complex() {
    /* digamma(2.5 + 3 I). */
    assert_complex_close("PolyGamma[2.5 + 3.0 I]", 1.281269, 0.979805, 1e-4);
    /* digamma(1 + I) = 0.0946503 + 1.076674 I. */
    assert_complex_close("PolyGamma[0, 1.0 + 1.0 I]", 0.0946503, 1.076674, 1e-5);
    /* trigamma(1 + I) = 0.463000 - 0.794234 I. */
    assert_complex_close("PolyGamma[1, 1.0 + 1.0 I]", 0.463000, -0.794234, 1e-5);
}

/* ---- derivatives ---------------------------------------------------- */

void test_polygamma_derivatives() {
    assert_eval_eq("D[PolyGamma[x], x]", "PolyGamma[1, x]", 0);
    assert_eval_eq("D[PolyGamma[n, x], x]", "PolyGamma[1 + n, x]", 0);
    assert_eval_eq("D[PolyGamma[2, x^2], x]", "2 x PolyGamma[3, x^2]", 0);
    /* Gamma now differentiates through PolyGamma. */
    assert_eval_eq("D[Gamma[x], x]", "Gamma[x] PolyGamma[0, x]", 0);
    assert_eval_eq("D[Gamma[z], {z, 2}]",
                   "Gamma[z] PolyGamma[1, z] + Gamma[z] PolyGamma[0, z]^2", 0);
}

/* ---- listability & attributes -------------------------------------- */

void test_polygamma_listable() {
    assert_eval_eq("PolyGamma[1, {1, 2}]", "{1/6 Pi^2, -1 + 1/6 Pi^2}", 0);
}

void test_polygamma_attributes() {
    SymbolDef* d = symtab_get_def("PolyGamma");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "PolyGamma must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "PolyGamma must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "PolyGamma must be Protected");

    SymbolDef* lg = symtab_get_def("LogGamma");
    ASSERT(lg != NULL);
    ASSERT_MSG((lg->attributes & ATTR_PROTECTED) != 0, "LogGamma must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_polygamma_canonical);
    TEST(test_polygamma_argt);
    TEST(test_polygamma_poles);
    TEST(test_polygamma_infinities);
    TEST(test_polygamma_digamma_integers);
    TEST(test_polygamma_odd_order_exact);
    TEST(test_polygamma_even_order_symbolic);
    TEST(test_polygamma_loggamma);
    TEST(test_polygamma_symbolic);
    TEST(test_polygamma_machine_real);
    TEST(test_polygamma_mpfr_real);
    TEST(test_polygamma_complex);
    TEST(test_polygamma_derivatives);
    TEST(test_polygamma_listable);
    TEST(test_polygamma_attributes);

    printf("All PolyGamma tests passed.\n");
    return 0;
}
