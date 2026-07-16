/* Tests for the logarithmic integral LogIntegral[z].
 *
 * Covers exact special values (0 -> 0, 1 -> -Infinity, Infinity -> Infinity,
 * ComplexInfinity/Indeterminate), machine real (z > 1 and 0 < z < 1, plus the
 * complex result for z < 0 via the branch cut), arbitrary-precision (MPFR)
 * reals with precision tracking, machine & arbitrary complex, derivatives
 * (including the chain rule), Listable threading, attributes and arity errors.
 *
 * Implementation rests on li(z) = Ei(Log z); reference values cross-checked
 * against the running binary and the task specification. */

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

void test_li_exact() {
    assert_eval_eq("LogIntegral[0]", "0", 0);
    assert_eval_eq("LogIntegral[1]", "-Infinity", 0);
    assert_eval_eq("LogIntegral[Infinity]", "Infinity", 0);
    assert_eval_eq("LogIntegral[ComplexInfinity]", "Indeterminate", 0);
    assert_eval_eq("LogIntegral[Indeterminate]", "Indeterminate", 0);
    /* Listable over the special list. */
    assert_eval_eq("LogIntegral[{0, 1, Infinity}]", "{0, -Infinity, Infinity}", 0);
}

/* ---- symbolic ------------------------------------------------------- */

void test_li_symbolic() {
    assert_eval_eq("LogIntegral[x]", "LogIntegral[x]", 0);
    /* Exact non-special arguments stay symbolic (no automatic exact value). */
    assert_eval_eq("LogIntegral[2]", "LogIntegral[2]", 0);
    assert_eval_eq("LogIntegral[1/2]", "LogIntegral[1/2]", 0);
    assert_eval_eq("LogIntegral[a + b]", "LogIntegral[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_li_machine_real() {
    /* z > 1 routes through mpfr_eint(Log z); 0 < z < 1 through the on-cut
     * real series with Log z < 0. */
    assert_close("LogIntegral[20.]", 9.9052999776328203, 1e-9);
    assert_close("LogIntegral[2.]", 1.0451637801174927, 1e-12);
    assert_close("LogIntegral[1.2]", -0.93378729266725769, 1e-11);
    assert_close("LogIntegral[1.5]", 0.12506498631529636, 1e-11);
    assert_close("LogIntegral[1.8]", 0.73263703111392142, 1e-11);
    /* 0 < z < 1: li(1/2). */
    assert_close("LogIntegral[0.5]", -0.37867104306108798, 1e-12);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_li_arbitrary_precision() {
    /* N[LogIntegral[2], 50] (reference kept short of the last digit). */
    assert_eval_startswith("N[LogIntegral[2], 50]",
        "1.045163780117492784844588889194613136522615578151");
    /* N[LogIntegral[1/2], 40] on the on-cut branch. */
    assert_eval_startswith("N[LogIntegral[1/2], 40]",
        "-0.378671043061087976727207184636560980");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("LogIntegral[2.000000000000000000000]",
                           "1.0451637801174927848");
}

/* ---- machine complex ------------------------------------------------ */

void test_li_machine_complex() {
    /* LogIntegral[2 + I] = 1.411259... + 1.224707... I. */
    assert_complex_close("LogIntegral[2. + I]",
                         1.4112590420178010, 1.2247069384103027, 1e-9);
    /* Conjugate symmetry off the cut: li[conj z] = conj li[z]. */
    assert_complex_close("LogIntegral[2. - I]",
                         1.4112590420178010, -1.2247069384103027, 1e-9);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_li_arbitrary_complex() {
    assert_eval_startswith("N[Re[LogIntegral[2 + I]], 30]",
        "1.41125904201780100568439320706");
    assert_eval_startswith("N[Im[LogIntegral[2 + I]], 30]",
        "1.22470693841030271349717531822");
}

/* ---- derivatives ---------------------------------------------------- */

void test_li_derivatives() {
    /* D[LogIntegral[x], x] = 1/Log[x]. */
    assert_eval_eq("D[LogIntegral[x], x]", "1/Log[x]", 0);
    /* Chain rule: D[LogIntegral[x^2], x] = (2 x)/Log[x^2]. */
    assert_eval_eq("D[LogIntegral[x^2], x]", "(2 x)/Log[x^2]", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_li_listable() {
    assert_eval_eq("LogIntegral[{}]", "{}", 0);
    assert_eval_eq("LogIntegral[{0, 1}]", "{0, -Infinity}", 0);
    /* Numeric threading (spec list). */
    Expr* e = parse_expression("LogIntegral[{1.2, 1.5, 1.8}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol.name, "List") == 0 &&
           r->data.function.arg_count == 3);
    double exp0[3] = { -0.93378729266725769, 0.12506498631529636, 0.73263703111392142 };
    for (int i = 0; i < 3; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= 1e-9,
                   "LogIntegral list element %d: expected %.15g got %.15g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_li_attributes() {
    SymbolDef* d = symtab_get_def("LogIntegral");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "LogIntegral must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "LogIntegral must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "LogIntegral must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_li_arity() {
    /* Wrong arity stays unevaluated (argx diagnostic to stderr). */
    assert_eval_eq("LogIntegral[]", "LogIntegral[]", 0);
    assert_eval_eq("LogIntegral[1, 2, 3]", "LogIntegral[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_li_exact);
    TEST(test_li_symbolic);
    TEST(test_li_machine_real);
    TEST(test_li_arbitrary_precision);
    TEST(test_li_machine_complex);
    TEST(test_li_arbitrary_complex);
    TEST(test_li_derivatives);
    TEST(test_li_listable);
    TEST(test_li_attributes);
    TEST(test_li_arity);

    printf("All LogIntegral tests passed.\n");
    return 0;
}
