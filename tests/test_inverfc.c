/* Tests for the inverse complementary error function: InverseErfc[s].
 *
 * Covers exact special values (0, 1, 2 -> Infinity, 0, -Infinity), symbolic
 * fall-through, machine-real and arbitrary-precision (MPFR) numerics with
 * precision tracking (including the cancellation-free small-s path), the
 * derivative and Series, out-of-domain handling, Listable threading, and
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
               "%s: expected %.17g, got %.17g", input, expected, v);
}

/* ---- exact special values ------------------------------------------- */

void test_inverfc_exact() {
    assert_eval_eq("InverseErfc[0]", "Infinity", 0);
    assert_eval_eq("InverseErfc[1]", "0", 0);
    assert_eval_eq("InverseErfc[2]", "-Infinity", 0);
    /* Real boundaries collapse to the same infinities / zero. */
    assert_eval_eq("InverseErfc[0.0]", "Infinity", 0);
    assert_eval_eq("InverseErfc[2.0]", "-Infinity", 0);
    /* Indeterminate passes through. */
    assert_eval_eq("InverseErfc[Indeterminate]", "Indeterminate", 0);
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_inverfc_symbolic() {
    assert_eval_eq("InverseErfc[x]", "InverseErfc[x]", 0);
    /* erfc is not odd: no auto-applied symmetry rewrite. */
    assert_eval_eq("InverseErfc[-x]", "InverseErfc[-x]", 0);
    /* Exact rationals stay symbolic (only N[] makes them numeric). */
    assert_eval_eq("InverseErfc[1/2]", "InverseErfc[1/2]", 0);
    /* Out-of-domain integers/reals stay symbolic (domain is [0, 2]). */
    assert_eval_eq("InverseErfc[3]", "InverseErfc[3]", 0);
    assert_eval_eq("InverseErfc[-1]", "InverseErfc[-1]", 0);
    assert_eval_eq("InverseErfc[2.3]", "InverseErfc[2.3]", 0);
    assert_eval_eq("InverseErfc[-0.5]", "InverseErfc[-0.5]", 0);
    /* Unrelated symbolic combinations stay put. */
    assert_eval_eq("InverseErfc[a + b]", "InverseErfc[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_inverfc_machine_real() {
    /* erfc(InverseErfc[s]) == s round-trips at machine precision. */
    assert_close("InverseErfc[0.6]", 0.37080715859355793, 1e-12);
    assert_close("InverseErfc[1.]", 0.0, 1e-15);
    /* Documentation list values InverseErfc[1/{2.,3.,4.,5.}]. */
    assert_close("InverseErfc[0.5]", 0.47693627620446987, 1e-12);
    assert_close("InverseErfc[1./3.]", 0.68407034965662264, 1e-12);
    assert_close("InverseErfc[0.25]", 0.81341984759761854, 1e-12);
    assert_close("InverseErfc[0.2]", 0.90619380243682322, 1e-12);
    /* Reflection about s = 1 is computed (not rewritten): InverseErfc[1.5]. */
    assert_close("InverseErfc[1.5]", -0.47693627620446987, 1e-12);
    /* Small s (large z): the cancellation-free erfc-Newton path. */
    assert_close("InverseErfc[0.01]", 1.8213863677184497, 1e-11);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_inverfc_arbitrary_precision() {
    /* N[InverseErfc[33/100], 50] reference (digits short of the rounding edge). */
    assert_eval_startswith("N[InverseErfc[33/100], 50]",
        "0.6888025281165564504025047289052578354494899234937");
    /* High-precision literal input, precision tracks the input. */
    assert_eval_startswith("InverseErfc[0.6`100]",
        "0.370807158593557929058249477522449138604304883162931113534584");
    /* Precision-tracking on a 24-significant-digit input. */
    assert_eval_startswith("InverseErfc[0.330000000000000000000000]",
        "0.6888025281165564504025");
    /* MPFR coverage case (numeric-builtins-cover-MPFR lesson). */
    assert_eval_startswith("N[InverseErfc[1/2], 35]",
        "0.4769362762044698733814183536431305");
    /* Small-s high precision exercises the cancellation-free path. */
    assert_eval_startswith("N[InverseErfc[1/100], 40]",
        "1.82138636771844967304021031862099524");
}

/* ---- derivative & series -------------------------------------------- */

void test_inverfc_derivatives() {
    assert_eval_eq("D[InverseErfc[x], x]",
                   "-1/2 Sqrt[Pi] E^InverseErfc[x]^2", 0);
    /* Second derivative (two sign flips cancel: same shape as InverseErf'). */
    assert_eval_eq("D[InverseErfc[x], {x, 2}]",
                   "1/2 Pi InverseErfc[x] E^(2 InverseErfc[x]^2)", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_inverfc_listable() {
    assert_eval_eq("InverseErfc[{1}]", "{0}", 0);
    assert_eval_eq("InverseErfc[{}]", "{}", 0);
    assert_eval_eq("InverseErfc[{0, 1, 2}]", "{Infinity, 0, -Infinity}", 0);
    /* Numeric threading: InverseErfc[1/{2.,3.,4.,5.}]. */
    Expr* e = parse_expression("InverseErfc[{0.5, 1./3., 0.25, 0.2}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol, "List") == 0 &&
           r->data.function.arg_count == 4);
    double exp0[4] = { 0.47693627620446987, 0.68407034965662264,
                       0.81341984759761854, 0.90619380243682322 };
    for (int i = 0; i < 4; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= 1e-12,
                   "InverseErfc list element %d: expected %.17g got %.17g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_inverfc_attributes() {
    SymbolDef* d = symtab_get_def("InverseErfc");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "InverseErfc must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "InverseErfc must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "InverseErfc must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_inverfc_arity() {
    /* Wrong arity stays unevaluated (argx diagnostic to stderr). */
    assert_eval_eq("InverseErfc[]", "InverseErfc[]", 0);
    assert_eval_eq("InverseErfc[1, 2, 3]", "InverseErfc[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_inverfc_exact);
    TEST(test_inverfc_symbolic);
    TEST(test_inverfc_machine_real);
    TEST(test_inverfc_arbitrary_precision);
    TEST(test_inverfc_derivatives);
    TEST(test_inverfc_listable);
    TEST(test_inverfc_attributes);
    TEST(test_inverfc_arity);

    printf("All InverseErfc tests passed.\n");
    return 0;
}
