/* Tests for the inverse error function: InverseErf[s] and InverseErf[z0, s].
 *
 * Covers exact special values (0, +-1 -> 0, +-Infinity), symbolic
 * fall-through and odd symmetry, machine-real and arbitrary-precision (MPFR)
 * numerics with precision tracking, the two-argument generalized form,
 * derivatives and Series, out-of-domain handling, Listable threading, and
 * attributes. Numerical references cross-checked against the running binary
 * and the Wolfram Language documentation. */

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
               "%s: expected %.12g, got %.12g", input, expected, v);
}

/* ---- exact special values ------------------------------------------- */

void test_inverf_exact() {
    assert_eval_eq("InverseErf[0]", "0", 0);
    assert_eval_eq("InverseErf[1]", "Infinity", 0);
    assert_eval_eq("InverseErf[-1]", "-Infinity", 0);
    /* Real boundaries collapse to the same infinities. */
    assert_eval_eq("InverseErf[1.0]", "Infinity", 0);
    assert_eval_eq("InverseErf[-1.0]", "-Infinity", 0);
    /* Indeterminate passes through. */
    assert_eval_eq("InverseErf[Indeterminate]", "Indeterminate", 0);
}

/* ---- symbolic fall-through & odd symmetry --------------------------- */

void test_inverf_symbolic() {
    assert_eval_eq("InverseErf[x]", "InverseErf[x]", 0);
    /* Odd: InverseErf[-x] = -InverseErf[x], InverseErf[-2 x] = -InverseErf[2 x]. */
    assert_eval_eq("InverseErf[-x]", "-InverseErf[x]", 0);
    assert_eval_eq("InverseErf[-2 x]", "-InverseErf[2 x]", 0);
    /* Exact rationals stay symbolic (only N[] makes them numeric). */
    assert_eval_eq("InverseErf[1/2]", "InverseErf[1/2]", 0);
    /* Out-of-domain integers/reals stay symbolic (no complex evaluation). */
    assert_eval_eq("InverseErf[2]", "InverseErf[2]", 0);
    assert_eval_eq("InverseErf[1.3]", "InverseErf[1.3]", 0);
    assert_eval_eq("InverseErf[-2.5]", "InverseErf[-2.5]", 0);
    /* ComplexInfinity / unrelated symbols stay put. */
    assert_eval_eq("InverseErf[a + b]", "InverseErf[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_inverf_machine_real() {
    /* erf(InverseErf[s]) == s round-trips at machine precision. */
    assert_close("InverseErf[0.6]", 0.5951160814499948, 1e-12);
    assert_close("InverseErf[0.33]", 0.3013321461337058, 1e-12);
    assert_close("InverseErf[-0.6]", -0.5951160814499948, 1e-12);
    assert_close("InverseErf[0.0]", 0.0, 0.0);
    /* Documentation list values InverseErf[1/{2.,3.,4.,5.}]. */
    assert_close("InverseErf[0.5]", 0.4769362762044699, 1e-10);
    assert_close("InverseErf[1./3.]", 0.3045701941739856, 1e-10);
    assert_close("InverseErf[0.25]", 0.2253120550121781, 1e-10);
    assert_close("InverseErf[0.2]", 0.1791434546212918, 1e-10);
    /* Round-trip check near the edge of the domain. */
    assert_close("InverseErf[0.99]", 1.8213863677184496, 1e-9);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_inverf_arbitrary_precision() {
    /* N[InverseErf[33/100], 50] reference (digits short of the rounding edge). */
    assert_eval_startswith("N[InverseErf[33/100], 50]",
        "0.3013321461337058261285027181583947739658242828285");
    /* High-precision literal input, precision tracks the input. */
    assert_eval_startswith("InverseErf[0.6`100]",
        "0.595116081449994850019300360168108253439616886279848409714612569202207183989666902521700511902573");
    /* Precision-tracking on a 24-significant-digit input. */
    assert_eval_startswith("InverseErf[0.330000000000000000000000]",
        "0.301332146133705826128503");
    /* MPFR coverage case (numeric-builtins-cover-MPFR lesson). */
    assert_eval_startswith("N[InverseErf[1/2], 35]",
        "0.4769362762044698733814183536431305");
}

/* ---- two-argument generalized form ---------------------------------- */

void test_inverf_two_arg() {
    /* InverseErf[z0, s] = InverseErf[s + Erf[z0]]. */
    assert_close("InverseErf[0.4, 0.2]", 0.6317759030550062, 1e-9);
    /* InverseErf[0, 1.3] = InverseErf[1.3] (out of domain, stays symbolic). */
    assert_eval_eq("InverseErf[0, 1.3]", "InverseErf[1.3]", 0);
    /* Fully symbolic two-arg stays in two-arg form (Erf[z0] does not reduce). */
    assert_eval_eq("InverseErf[a, b]", "InverseErf[a, b]", 0);
}

/* ---- derivatives & series ------------------------------------------- */

void test_inverf_derivatives() {
    assert_eval_eq("D[InverseErf[x], x]",
                   "1/2 Sqrt[Pi] E^InverseErf[x]^2", 0);
    /* Second derivative (chain rule on the InverseErf-containing exponent). */
    assert_eval_eq("D[InverseErf[x], {x, 2}]",
                   "1/2 Pi InverseErf[x] E^(2 InverseErf[x]^2)", 0);
    /* Series at the origin via the generic D-based Taylor fallback. */
    assert_eval_startswith("Series[InverseErf[x], {x, 0, 3}]",
                           "1/2 Sqrt[Pi] x");
    assert_eval_startswith("Series[InverseErf[x], {x, 0, 8}]",
                           "1/2 Sqrt[Pi] x + 1/24 Pi^(3/2) x^3 + 7/960 Pi^(5/2) x^5 + "
                           "127/80640 Pi^(7/2) x^7");
}

/* ---- Listable & attributes ------------------------------------------ */

void test_inverf_listable() {
    assert_eval_eq("InverseErf[{0}]", "{0}", 0);
    assert_eval_eq("InverseErf[{}]", "{}", 0);
    assert_eval_eq("InverseErf[{0, -1, 1}]", "{0, -Infinity, Infinity}", 0);
    /* Numeric threading. */
    Expr* e = parse_expression("InverseErf[{0.2, 0.5, 0.99}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol.name, "List") == 0 &&
           r->data.function.arg_count == 3);
    double exp0[3] = { 0.1791434546212918, 0.4769362762044699, 1.8213863677184496 };
    for (int i = 0; i < 3; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= 1e-9,
                   "InverseErf list element %d: expected %.12g got %.12g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_inverf_attributes() {
    SymbolDef* d = symtab_get_def("InverseErf");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "InverseErf must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "InverseErf must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "InverseErf must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_inverf_arity() {
    /* Wrong arity stays unevaluated (argt diagnostic to stderr). */
    assert_eval_eq("InverseErf[]", "InverseErf[]", 0);
    assert_eval_eq("InverseErf[1, 2, 3, 4, 5]", "InverseErf[1, 2, 3, 4, 5]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_inverf_exact);
    TEST(test_inverf_symbolic);
    TEST(test_inverf_machine_real);
    TEST(test_inverf_arbitrary_precision);
    TEST(test_inverf_two_arg);
    TEST(test_inverf_derivatives);
    TEST(test_inverf_listable);
    TEST(test_inverf_attributes);
    TEST(test_inverf_arity);

    printf("All InverseErf tests passed.\n");
    return 0;
}
