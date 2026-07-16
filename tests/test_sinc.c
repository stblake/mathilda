/* Tests for the cardinal sine Sinc[z] = Sin[z]/z, with Sinc[0] = 1.
 *
 * Covers exact special values (0 -> 1, +-Infinity -> 0, ComplexInfinity ->
 * Indeterminate), machine real, arbitrary-precision (MPFR) reals with precision
 * tracking, machine & arbitrary complex, symbolic passthrough, derivatives,
 * Taylor Series, Listable threading, attributes and arity errors. Reference
 * values cross-checked against the running binary. */

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

void test_sinc_exact() {
    assert_eval_eq("Sinc[0]", "1", 0);
    assert_eval_eq("Sinc[Infinity]", "0", 0);
    assert_eval_eq("Sinc[-Infinity]", "0", 0);
    assert_eval_eq("Sinc[ComplexInfinity]", "Indeterminate", 0);
}

/* ---- symbolic passthrough ------------------------------------------- */

void test_sinc_symbolic() {
    assert_eval_eq("Sinc[x]", "Sinc[x]", 0);
    /* Sinc is even, but WL does not auto-fold Sinc[-x]. */
    assert_eval_eq("Sinc[-x]", "Sinc[-x]", 0);
    assert_eval_eq("Sinc[a + b]", "Sinc[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_sinc_machine_real() {
    assert_close("Sinc[1.0]", 0.8414709848078965, 1e-12);   /* Sin[1]/1 */
    assert_close("Sinc[2.0]", 0.4546487134128409, 1e-12);   /* Sin[2]/2 */
    assert_close("Sinc[3.0]", 0.0470400026866224, 1e-12);   /* Sin[3]/3 */
    /* Continuity at the origin from the machine path. */
    assert_close("Sinc[0.0]", 1.0, 1e-15);
    assert_close("Sinc[1.0*10^-8]", 1.0, 1e-14);
    /* Even: Sinc[-x] = Sinc[x] numerically. */
    assert_close("Sinc[-2.0]", 0.4546487134128409, 1e-12);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_sinc_arbitrary_precision() {
    assert_eval_startswith("N[Sinc[2], 45]",
        "0.454648713412840847698009932955872421351127485");
    assert_eval_startswith("N[Sinc[1], 40]",
        "0.84147098480789650665250232163029899962");
    /* Precision tracking. */
    assert_eval_startswith("Sinc[2.0000000000000000000000]",
                           "0.45464871341284084769");
}

/* ---- complex -------------------------------------------------------- */

void test_sinc_machine_complex() {
    assert_complex_close("Sinc[1.0 + I]", 0.9667107481003567, -0.3317468333156206, 1e-10);
}

void test_sinc_arbitrary_complex() {
    assert_eval_startswith("N[Sinc[1 + I], 30]",
        "0.96671074810035670154056228439");
}

/* ---- derivatives ---------------------------------------------------- */

void test_sinc_derivatives() {
    assert_eval_eq("D[Sinc[x], x]", "Cos[x]/x - Sin[x]/x^2", 0);
}

/* ---- Series --------------------------------------------------------- */

void test_sinc_series() {
    assert_eval_eq("Series[Sinc[x], {x, 0, 6}]",
        "1 - 1/6 x^2 + 1/120 x^4 - 1/5040 x^6 + O[x]^7", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_sinc_listable() {
    assert_eval_eq("Sinc[{}]", "{}", 0);
    assert_eval_eq("Sinc[{0}]", "{1}", 0);
}

void test_sinc_attributes() {
    SymbolDef* d = symtab_get_def("Sinc");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Sinc must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "Sinc must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Sinc must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_sinc_arity() {
    assert_eval_eq("Sinc[]", "Sinc[]", 0);
    assert_eval_eq("Sinc[1, 2, 3]", "Sinc[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_sinc_exact);
    TEST(test_sinc_symbolic);
    TEST(test_sinc_machine_real);
    TEST(test_sinc_arbitrary_precision);
    TEST(test_sinc_machine_complex);
    TEST(test_sinc_arbitrary_complex);
    TEST(test_sinc_derivatives);
    TEST(test_sinc_series);
    TEST(test_sinc_listable);
    TEST(test_sinc_attributes);
    TEST(test_sinc_arity);

    printf("All Sinc tests passed.\n");
    return 0;
}
