/* Tests for ProductLog[z] / ProductLog[k, z], the Lambert W function.
 *
 * Covers: exact special values and symbolic passthrough; machine- and
 * arbitrary-precision numerics (real and complex) verified both against the
 * documented reference values and by the self-consistency identity
 * W(z) e^W(z) == z; non-principal branches; list threading; the derivative
 * rule (one- and two-argument forms); Series at 0, at the branch point -1/E,
 * and at Infinity; SeriesCoefficient (numeric and the symbolic general term);
 * argument-count diagnostics; attributes.
 *
 * Numeric correctness is checked inside the language (via comparison/Abs
 * predicates that reduce to True) rather than by parsing printed reals. */

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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- numeric helpers ------------------------------------------------- */

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

/* ---- exact values and passthrough ----------------------------------- */

void test_productlog_exact() {
    assert_eval_eq("ProductLog[0]", "0", 0);
    assert_eval_eq("ProductLog[0, 0]", "0", 0);
    /* W_k(0) = -Infinity for k != 0. */
    assert_eval_eq("ProductLog[1, 0]", "-Infinity", 0);
    assert_eval_eq("ProductLog[-1, 0]", "-Infinity", 0);
    /* Closed-form fixed points. */
    assert_eval_eq("ProductLog[E]", "1", 0);
    assert_eval_eq("ProductLog[-1/E]", "-1", 0);
    /* Branches 0 and -1 coincide at the branch point. */
    assert_eval_eq("ProductLog[-1, -1/E]", "-1", 0);
    assert_eval_eq("ProductLog[Infinity]", "Infinity", 0);
    assert_eval_eq("ProductLog[ComplexInfinity]", "Infinity", 0);
    /* ProductLog[-Pi/2] = I Pi/2. */
    assert_eval_eq("Re[ProductLog[-Pi/2]]", "0", 0);
    assert_eval_eq("Im[ProductLog[-Pi/2]]", "1/2 Pi", 0);
}

void test_productlog_symbolic() {
    /* Exact non-special arguments stay symbolic (no auto-numericize). */
    assert_eval_eq("ProductLog[1/3]", "ProductLog[1/3]", 0);
    assert_eval_eq("ProductLog[2]", "ProductLog[2]", 0);
    assert_eval_eq("ProductLog[x]", "ProductLog[x]", 0);
    /* Non-integer branch index stays symbolic. */
    assert_eval_eq("ProductLog[x, 2.5]", "ProductLog[x, 2.5]", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_productlog_machine() {
    /* Documented reference values. */
    assert_close("ProductLog[1.0]", 0.56714329040978387, 1e-9);
    assert_close("ProductLog[2.5]", 0.95858635672870690, 1e-9);
    /* Inverse identity W(z) e^W(z) = z (real and complex). */
    assert_eval_eq("Abs[ProductLog[3.] Exp[ProductLog[3.]] - 3.] < 10^-12",
                   "True", 0);
    assert_eval_eq(
        "Abs[ProductLog[1 + 3.5 I] Exp[ProductLog[1 + 3.5 I]] - (1 + 3.5 I)] < 10^-12",
        "True", 0);
    /* Negative real < -1/e gives a genuinely complex principal value. */
    assert_eval_startswith("ProductLog[-1.5]", "-0.0327837");
}

void test_productlog_branches() {
    /* W_k(z) e^{W_k(z)} = z holds on every branch. */
    assert_eval_eq(
        "With[{w = ProductLog[-2, 2.3]}, Abs[w Exp[w] - 2.3] < 10^-10]", "True", 0);
    assert_eval_eq(
        "With[{w = ProductLog[2, 2.3]}, Abs[w Exp[w] - 2.3] < 10^-10]", "True", 0);
    assert_eval_eq(
        "With[{w = ProductLog[-1, 0.3]}, Abs[w Exp[w] - 0.3] < 10^-10]", "True", 0);
    /* Ordering by imaginary part: k = -2 is the most negative imaginary part. */
    assert_eval_startswith("Table[ProductLog[k, 2.3], {k, -2, 2}]", "{-1.56175");
}

void test_productlog_listable() {
    assert_eval_startswith("ProductLog[{1.5, 3.75, 5.5, 7.25}]", "{0.725861");
    /* Listable threads, leaving non-numeric entries symbolic. */
    assert_eval_eq("ProductLog[{0, E, -1/E}]", "{0, 1, -1}", 0);
}

/* ---- arbitrary precision -------------------------------------------- */

void test_productlog_highprec() {
    /* N[ProductLog[1/3], 100] -- check a generous leading prefix. */
    assert_eval_startswith(
        "N[ProductLog[1/3], 100]",
        "0.2576276530497367042829162016260977909096");
    /* Self-consistency at 100 digits on a rational argument. */
    assert_eval_eq(
        "With[{w = ProductLog[N[7/3, 100]]}, Abs[w Exp[w] - 7/3] < 10^-90]",
        "True", 0);
    /* Complex, high precision. */
    assert_eval_eq(
        "With[{z = N[1 + 3 I, 80], w = ProductLog[N[1 + 3 I, 80]]}, "
        "Abs[w Exp[w] - z] < 10^-70]", "True", 0);
}

/* ---- derivative ------------------------------------------------------ */

void test_productlog_derivative() {
    assert_eval_eq("D[ProductLog[z], z]",
                   "ProductLog[z]/(z (1 + ProductLog[z]))", 0);
    assert_eval_eq("D[ProductLog[2, z], z]",
                   "ProductLog[2, z]/(z (1 + ProductLog[2, z]))", 0);
    /* Chain rule through a scaled argument (the factor of 3 cancels). */
    assert_eval_eq("D[ProductLog[3 z], z]",
                   "ProductLog[3 z]/(z (1 + ProductLog[3 z]))", 0);
}

/* ---- series ---------------------------------------------------------- */

void test_productlog_series() {
    assert_eval_eq("Series[ProductLog[x], {x, 0, 5}]",
                   "x - x^2 + 3/2 x^3 - 8/3 x^4 + 125/24 x^5 + O[x]^6", 0);
    /* Branch point -1/E: Puiseux series in Sqrt[x + 1/E]. */
    assert_eval_startswith("Series[ProductLog[x], {x, -1/E, 2}]",
                           "-1 + Sqrt[2 E] Sqrt[x - -1/E]");
    /* Infinity: nested-logarithm asymptotic expansion (the x^0 coefficient). */
    assert_eval_eq("Normal[Series[ProductLog[x], {x, Infinity, 0}]]",
                   "Log[x] - Log[Log[x]] + Log[Log[x]]/Log[x] - "
                   "Log[Log[x]]/Log[x]^2 + 1/2 Log[Log[x]]^2/Log[x]^2", 0);
}

void test_productlog_seriescoefficient() {
    assert_eval_eq("SeriesCoefficient[ProductLog[x], {x, 0, 1}]", "1", 0);
    assert_eval_eq("SeriesCoefficient[ProductLog[x], {x, 0, 3}]", "3/2", 0);
    assert_eval_eq("SeriesCoefficient[ProductLog[x], {x, 0, 5}]", "125/24", 0);
    assert_eval_eq("SeriesCoefficient[ProductLog[x], {x, 0, 4}]", "-8/3", 0);
    /* Symbolic general term: (-n)^(n-1)/n! for n >= 1, else 0. */
    assert_eval_eq("SeriesCoefficient[ProductLog[x], {x, 0, n}]",
                   "Piecewise[{{(-n)^(-1 + n)/Factorial[n], n >= 1}}, 0]", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_productlog_argcount() {
    /* Wrong arg counts stay unevaluated (an argt message goes to stderr). */
    assert_eval_eq("ProductLog[]", "ProductLog[]", 0);
    assert_eval_eq("ProductLog[1, 2, 3]", "ProductLog[1, 2, 3]", 0);
}

void test_productlog_attributes() {
    SymbolDef* d = symtab_get_def("ProductLog");
    ASSERT_MSG(d != NULL, "ProductLog not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "ProductLog not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "ProductLog not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "ProductLog not Protected");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_productlog_exact);
    TEST(test_productlog_symbolic);
    TEST(test_productlog_machine);
    TEST(test_productlog_branches);
    TEST(test_productlog_listable);
    TEST(test_productlog_highprec);
    TEST(test_productlog_derivative);
    TEST(test_productlog_series);
    TEST(test_productlog_seriescoefficient);
    TEST(test_productlog_argcount);
    TEST(test_productlog_attributes);

    printf("All ProductLog tests passed.\n");
    return 0;
}
