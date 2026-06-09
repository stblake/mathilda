/* Tests for LogGamma[z], the log-gamma function log(Gamma(z)).
 *
 * Covers the exact closed forms (integers -> Log[(n-1)!], positive and negative
 * half-integers, with the negative-axis branch term), the Infinity poles, the
 * symbolic infinities, Listable threading (including the empty list), machine
 * and arbitrary-precision (MPFR) real and complex numerics, the continuous
 * (large-imaginary-part) branch that distinguishes LogGamma from Log[Gamma],
 * symbolic fall-through, the derivative rule, the argument-count diagnostic,
 * and attributes. */

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

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
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
               "%s: expected %.15g, got %.15g", input, expected, v);
}

/* ---- exact closed forms: integers ----------------------------------- */

void test_loggamma_integers() {
    /* LogGamma[n] = Log[(n-1)!]. */
    assert_eval_eq("LogGamma[1]", "0", 0);            /* Log[0!] = Log[1] = 0 */
    assert_eval_eq("LogGamma[2]", "0", 0);            /* Log[1!] = 0 */
    assert_eval_eq("LogGamma[5]", "Log[24]", 0);      /* Log[4!] */
    assert_eval_eq("LogGamma[6]", "Log[120]", 0);     /* Log[5!] */
    /* Non-positive integers are poles -> Infinity (not ComplexInfinity). */
    assert_eval_eq("LogGamma[0]",  "Infinity", 0);
    assert_eval_eq("LogGamma[-1]", "Infinity", 0);
    assert_eval_eq("LogGamma[-5]", "Infinity", 0);
}

/* ---- exact closed forms: half-integers ------------------------------ */

void test_loggamma_half_integers() {
    /* Positive half-integers: Log of the exact Sqrt[Pi] form. */
    assert_eval_eq("LogGamma[1/2]", "Log[Sqrt[Pi]]", 0);
    assert_eval_eq("LogGamma[3/2]", "Log[1/2 Sqrt[Pi]]", 0);
    assert_eval_eq("LogGamma[7/2]", "Log[15/8 Sqrt[Pi]]", 0);
    /* Negative half-integers carry the branch term -Ceiling[-z] Pi I. */
    assert_eval_eq("LogGamma[-1/2]", "-I Pi + Log[2 Sqrt[Pi]]", 0);
    assert_eval_eq("LogGamma[-3/2]", "-(2*I) Pi + Log[4/3 Sqrt[Pi]]", 0);
    assert_eval_eq("LogGamma[-7/2]", "-(4*I) Pi + Log[16/105 Sqrt[Pi]]", 0);
}

/* ---- symbolic special values ---------------------------------------- */

void test_loggamma_special() {
    assert_eval_eq("LogGamma[Infinity]",        "Infinity", 0);
    assert_eval_eq("LogGamma[-Infinity]",       "Indeterminate", 0);
    assert_eval_eq("LogGamma[I Infinity]",      "ComplexInfinity", 0);
    assert_eval_eq("LogGamma[ComplexInfinity]", "ComplexInfinity", 0);
    assert_eval_eq("LogGamma[Indeterminate]",   "Indeterminate", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_loggamma_listable() {
    assert_eval_eq("LogGamma[{2, 3, 4, 5, 6}]",
                   "{0, Log[2], Log[6], Log[24], Log[120]}", 0);
    /* Singular points all diverge to Infinity. */
    assert_eval_eq("LogGamma[{0, -1, -2, -3}]",
                   "{Infinity, Infinity, Infinity, Infinity}", 0);
    /* Threading over the empty list yields the empty list. */
    assert_eval_eq("LogGamma[{}]", "{}", 0);
}

/* ---- numeric: machine precision ------------------------------------- */

void test_loggamma_machine_real() {
    assert_close("LogGamma[2.5]",   0.28468287047291916, 1e-13);
    assert_close("LogGamma[0.5]",   0.572364942924700082, 1e-13);
    assert_close("LogGamma[10.5]",  13.9406252194037636, 1e-12);
    assert_close("LogGamma[100.5]", 361.43554046777763, 1e-10);
    /* Large positive argument stays finite where Gamma overflows. */
    assert_close("LogGamma[200.0]", 857.93366982585735, 1e-9);
}

void test_loggamma_machine_real_negative() {
    /* z < 0 non-integer: complex, Re = log|Gamma|, Im = -Pi Ceiling[-z].
     * Cross-checks the exact half-integer value LogGamma[-3/2]. */
    assert_close("Re[LogGamma[-1.5]]", 0.860047015376481011, 1e-12);
    assert_close("Im[LogGamma[-1.5]]", -2.0 * M_PI, 1e-12);
    assert_close("Im[LogGamma[-0.5]]", -1.0 * M_PI, 1e-12);
    assert_close("Im[LogGamma[-2.5]]", -3.0 * M_PI, 1e-12);
}

void test_loggamma_machine_complex() {
    assert_close("Re[LogGamma[2.5 + 3 I]]", -1.47095461034884169, 1e-11);
    assert_close("Im[LogGamma[2.5 + 3 I]]",  2.82261563826079945, 1e-11);
    /* Continuous branch: Im grows past Pi where Log[Gamma] would wrap. */
    assert_close("Im[LogGamma[10.0 + 8.0 I]]", 18.805440811251583, 1e-9);
    /* Reflection (Re < 1/2). */
    assert_close("Re[LogGamma[-0.5 + 2.0 I]]", -2.9461153555214209, 1e-10);
    assert_close("Im[LogGamma[-0.5 + 2.0 I]]", -2.4083119718987954, 1e-10);
}

/* ---- numeric: arbitrary precision (MPFR) ---------------------------- */

void test_loggamma_arbitrary_precision_real() {
    /* The prompt's 70-digit value of LogGamma[2.5]. */
    assert_eval_startswith(
        "LogGamma[2.5`70]",
        "0.284682870472919159632494669682701924320137695559894");
    /* N[LogGamma[22/10], 50] -- a safe common prefix (low digits unpinned). */
    assert_eval_startswith(
        "N[LogGamma[22/10], 50]",
        "0.0969474667906387764920151858546291862377217");
    /* High-precision value at a half-integer matches the 2.5 value. */
    assert_eval_startswith(
        "N[LogGamma[5/2], 40]",
        "0.284682870472919159632494669682701924320");
}

void test_loggamma_arbitrary_precision_complex() {
    /* Stirling MPFR complex path; verified self-consistent across precision
     * and against the machine result. */
    assert_eval_startswith(
        "Re[LogGamma[2.5`50 + 3`50 I]]",
        "-1.4709546103488416913054992949780750150");
    assert_eval_startswith(
        "Im[LogGamma[2.5`50 + 3`50 I]]",
        "2.8226156382607994500252655473187103700");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_loggamma_symbolic() {
    assert_eval_eq("LogGamma[x]",      "LogGamma[x]", 0);
    assert_eval_eq("LogGamma[1 + x]",  "LogGamma[1 + x]", 0);
    /* Exact non-(half-)integer rationals stay symbolic. */
    assert_eval_eq("LogGamma[1/3]",    "LogGamma[1/3]", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_loggamma_derivative() {
    /* D[LogGamma[x], x] = PolyGamma[0, x]. */
    assert_eval_eq("D[LogGamma[x], x]", "PolyGamma[0, x]", 0);
    /* Higher derivatives raise the PolyGamma order. */
    assert_eval_eq("D[LogGamma[x], {x, 2}]", "PolyGamma[1, x]", 0);
    assert_eval_eq("D[LogGamma[x], {x, 3}]", "PolyGamma[2, x]", 0);
    /* Chain rule. */
    assert_eval_eq("D[LogGamma[x^2], x]", "2 x PolyGamma[0, x^2]", 0);
}

/* ---- argument-count diagnostic -------------------------------------- */

void test_loggamma_argcount() {
    /* Wrong arities stay unevaluated (an argx message is printed to stderr). */
    assert_eval_eq("LogGamma[]",          "LogGamma[]", 0);
    assert_eval_eq("LogGamma[1, 2, 3, 4]","LogGamma[1, 2, 3, 4]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_loggamma_attributes() {
    SymbolDef* d = symtab_get_def("LogGamma");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "LogGamma must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "LogGamma must be a NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "LogGamma must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_loggamma_integers);
    TEST(test_loggamma_half_integers);
    TEST(test_loggamma_special);
    TEST(test_loggamma_listable);
    TEST(test_loggamma_machine_real);
    TEST(test_loggamma_machine_real_negative);
    TEST(test_loggamma_machine_complex);
    TEST(test_loggamma_arbitrary_precision_real);
    TEST(test_loggamma_arbitrary_precision_complex);
    TEST(test_loggamma_symbolic);
    TEST(test_loggamma_derivative);
    TEST(test_loggamma_argcount);
    TEST(test_loggamma_attributes);

    printf("All LogGamma tests passed.\n");
    return 0;
}
