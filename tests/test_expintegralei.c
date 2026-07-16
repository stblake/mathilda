/* Tests for the exponential integral ExpIntegralEi[z].
 *
 * Covers exact special values (0 -> -Infinity, +-Infinity, +-I Infinity ->
 * +-I Pi, ComplexInfinity/Indeterminate), machine real (positive via mpfr_eint,
 * negative via the on-cut real series), arbitrary-precision (MPFR) reals with
 * precision tracking, machine & arbitrary complex, branch behaviour across the
 * negative real axis, derivatives, Listable threading, attributes and arity
 * errors. Reference values cross-checked against the running binary and the
 * task specification. */

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

void test_ei_exact() {
    assert_eval_eq("ExpIntegralEi[0]", "-Infinity", 0);
    assert_eval_eq("ExpIntegralEi[Infinity]", "Infinity", 0);
    assert_eval_eq("ExpIntegralEi[-Infinity]", "0", 0);
    assert_eval_eq("ExpIntegralEi[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("ExpIntegralEi[Indeterminate]", "Indeterminate", 0);
    /* On the imaginary axis Ei tends to +-I Pi. */
    assert_eval_eq("ExpIntegralEi[I Infinity]", "I Pi", 0);
    assert_eval_eq("ExpIntegralEi[-I Infinity]", "-I Pi", 0);
    /* Listable over the full special list from the spec. */
    assert_eval_eq("ExpIntegralEi[{Infinity, -Infinity, I Infinity, -I Infinity}]",
                   "{Infinity, 0, I Pi, -I Pi}", 0);
}

/* ---- symbolic ------------------------------------------------------- */

void test_ei_symbolic() {
    assert_eval_eq("ExpIntegralEi[x]", "ExpIntegralEi[x]", 0);
    /* Exact non-zero arguments stay symbolic (no automatic exact value). */
    assert_eval_eq("ExpIntegralEi[2]", "ExpIntegralEi[2]", 0);
    assert_eval_eq("ExpIntegralEi[1/2]", "ExpIntegralEi[1/2]", 0);
    assert_eval_eq("ExpIntegralEi[a + b]", "ExpIntegralEi[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_ei_machine_real() {
    assert_close("ExpIntegralEi[0.5]", 0.45421990486317358, 1e-12);
    assert_close("ExpIntegralEi[1.0]", 1.8951178163559368, 1e-12);
    assert_close("ExpIntegralEi[1.2]", 2.4420922851926516, 1e-11);
    assert_close("ExpIntegralEi[1.5]", 3.3012854491297978, 1e-11);
    assert_close("ExpIntegralEi[1.8]", 4.2498675574879335, 1e-11);
    assert_close("ExpIntegralEi[2.0]", 4.9542343560018902, 1e-10);
    assert_close("ExpIntegralEi[3.0]", 9.9338325706254166, 1e-9);
    /* Negative arguments take the on-cut real series and stay real. */
    assert_close("ExpIntegralEi[-1.0]", -0.21938393439552027, 1e-12);
    assert_close("ExpIntegralEi[-2.0]", -0.048900510708061119, 1e-13);
    assert_close("ExpIntegralEi[-5.0]", -0.0011482955912753257, 1e-15);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_ei_arbitrary_precision() {
    /* N[ExpIntegralEi[2], 50] (reference kept short of the last digit). */
    assert_eval_startswith("N[ExpIntegralEi[2], 50]",
        "4.954234356001890163379505130227035275518053562420");
    /* N[ExpIntegralEi[1], 40]. */
    assert_eval_startswith("N[ExpIntegralEi[1], 40]",
        "1.89511781635593675546652093433");
    /* Negative argument at high precision (on-cut series). */
    assert_eval_startswith("N[ExpIntegralEi[-2], 35]",
        "-0.04890051070806111956723983");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("ExpIntegralEi[2.0000000000000000000000]",
                           "4.9542343560018901633");
}

/* ---- machine complex ------------------------------------------------ */

void test_ei_machine_complex() {
    /* ExpIntegralEi[2 + I] = 4.06998... + 3.40094... I. */
    assert_complex_close("ExpIntegralEi[2. + I]",
                         4.0699809478939277, 3.4009439698001216, 1e-9);
    /* Conjugate symmetry: Ei[conj z] = conj Ei[z]. */
    assert_complex_close("ExpIntegralEi[2. - I]",
                         4.0699809478939277, -3.4009439698001216, 1e-9);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_ei_arbitrary_complex() {
    /* N[Re/Im[ExpIntegralEi[2 + I]], 30]. */
    assert_eval_startswith("N[Re[ExpIntegralEi[2 + I]], 30]",
        "4.0699809478939277422876902552");
    assert_eval_startswith("N[Im[ExpIntegralEi[2 + I]], 30]",
        "3.4009439698001216216304046260");
}

/* ---- branch cut behaviour ------------------------------------------- */

void test_ei_branch() {
    /* Approaching the cut from above (+I0) gives +I Pi; from below, -I Pi. The
     * shared real part is the principal value Ei(-1) = -0.219384. */
    assert_complex_close("ExpIntegralEi[-1. + 10^-10 I]",
                         -0.21938393439552027, 3.14159265, 1e-4);
    assert_complex_close("ExpIntegralEi[-1. - 10^-10 I]",
                         -0.21938393439552027, -3.14159265, 1e-4);
}

/* ---- large |z|: asymptotic expansion -------------------------------- */

void test_ei_large_magnitude() {
    /* For large |z| the convergent series is infeasible (and its guard-bit
     * arithmetic used to overflow, aborting in mpfr_init2); the asymptotic
     * series takes over. At 8 digits z = -50 I routes through the asymptotic
     * path and must agree with the 25-digit convergent reference:
     *   Ei(-50 I) = Ci(50) - I (Pi/2 + Si(50)). */
    assert_eval_startswith("N[Re[ExpIntegralEi[-50 I]], 8]", "-0.0056283");
    assert_eval_startswith("N[Im[ExpIntegralEi[-50 I]], 8]", "-3.122413");
    /* +50 I is the conjugate (+I Pi branch). */
    assert_eval_startswith("N[Im[ExpIntegralEi[50 I]], 8]", "3.122413");

    /* Re z << 0: e^z underflows, the result collapses to the branch constant
     * -I Pi (Im z < 0) with a vanishing real part. */
    assert_eval_startswith("N[Im[ExpIntegralEi[-80 - 2 I]], 12]",
                           "-3.14159265358");

    /* Regression: the reported abort. N numericises -10^60 I, driving |z| to
     * 1e60; this must now yield a finite numeric result, not crash. */
    {
        Expr* e = parse_expression("N[ExpIntegralEi[-10^60 I] + I Pi, 20]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        expr_free(e);
        ASSERT_MSG(r != NULL && (r->type == EXPR_REAL ||
                   (r->type == EXPR_FUNCTION &&
                    strcmp(r->data.function.head->data.symbol.name, "Complex") == 0)),
                   "N[ExpIntegralEi[-10^60 I] + I Pi, 20] must be finite numeric");
        expr_free(r);
    }
    /* Real large-negative argument: Ei(-10^60) ~ 0 (e^x underflow), no crash. */
    assert_eval_startswith("N[ExpIntegralEi[-10^60], 20]", "-0.");
}

/* ---- derivatives ---------------------------------------------------- */

void test_ei_derivatives() {
    /* D[ExpIntegralEi[x], x] = E^x / x. */
    assert_eval_eq("D[ExpIntegralEi[x], x]", "E^x/x", 0);
    /* Chain rule: D[ExpIntegralEi[x^2], x] = 2 x E^(x^2)/x^2 = 2 E^(x^2)/x. */
    assert_eval_eq("D[ExpIntegralEi[x^2], x]", "(2 E^x^2)/x", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_ei_listable() {
    assert_eval_eq("ExpIntegralEi[{}]", "{}", 0);
    assert_eval_eq("ExpIntegralEi[{0}]", "{-Infinity}", 0);
    /* Numeric threading (spec list). */
    Expr* e = parse_expression("ExpIntegralEi[{1.2, 1.5, 1.8}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol.name, "List") == 0 &&
           r->data.function.arg_count == 3);
    double exp0[3] = { 2.4420922851926516, 3.3012854491297978, 4.2498675574879335 };
    for (int i = 0; i < 3; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= 1e-9,
                   "ExpIntegralEi list element %d: expected %.15g got %.15g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_ei_attributes() {
    SymbolDef* d = symtab_get_def("ExpIntegralEi");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "ExpIntegralEi must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "ExpIntegralEi must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "ExpIntegralEi must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_ei_arity() {
    /* Wrong arity stays unevaluated (argx diagnostic to stderr). */
    assert_eval_eq("ExpIntegralEi[]", "ExpIntegralEi[]", 0);
    assert_eval_eq("ExpIntegralEi[1, 2, 3]", "ExpIntegralEi[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ei_exact);
    TEST(test_ei_symbolic);
    TEST(test_ei_machine_real);
    TEST(test_ei_arbitrary_precision);
    TEST(test_ei_machine_complex);
    TEST(test_ei_arbitrary_complex);
    TEST(test_ei_branch);
    TEST(test_ei_large_magnitude);
    TEST(test_ei_derivatives);
    TEST(test_ei_listable);
    TEST(test_ei_attributes);
    TEST(test_ei_arity);

    printf("All ExpIntegralEi tests passed.\n");
    return 0;
}
