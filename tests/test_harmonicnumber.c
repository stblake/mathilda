/* Tests for HarmonicNumber[n] and HarmonicNumber[n, r].
 *
 * Covers: exact rationals at non-negative integer n (order 1 and general r),
 * the explicit symbolic sum for symbolic order, HarmonicNumber[Infinity, r] =
 * Zeta[r], the Faulhaber polynomial at non-positive integer order, machine and
 * arbitrary-precision (MPFR) numerics with precision tracking, complex order,
 * the symbolic-passthrough cases, Listable threading, argument-count errors,
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

/* ---- numeric helpers (mirrors test_polygamma.c) --------------------- */

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

/* ---- exact integer n, order 1 --------------------------------------- */

void test_hn_integer_order1() {
    assert_eval_eq("HarmonicNumber[0]", "0", 0);
    assert_eval_eq("HarmonicNumber[1]", "1", 0);
    assert_eval_eq("HarmonicNumber[2]", "3/2", 0);
    assert_eval_eq("HarmonicNumber[3]", "11/6", 0);
    assert_eval_eq("HarmonicNumber[10]", "7381/2520", 0);
    /* The canonical first-ten table. */
    assert_eval_eq("Table[HarmonicNumber[n], {n, 10}]",
        "{1, 3/2, 11/6, 25/12, 137/60, 49/20, 363/140, 761/280, 7129/2520, "
        "7381/2520}", 0);
}

/* ---- exact integer n, general order r ------------------------------- */

void test_hn_integer_orderr() {
    /* 1 + 1/4 + 1/9 + 1/16 + 1/25 = 5269/3600 */
    assert_eval_eq("HarmonicNumber[5, 2]", "5269/3600", 0);
    assert_eval_eq("HarmonicNumber[3, 3]", "251/216", 0);
    assert_eval_eq("HarmonicNumber[0, 7]", "0", 0);
    /* symbolic order over an integer n expands to an explicit sum */
    assert_eval_eq("HarmonicNumber[4, r]", "1 + 2^(-r) + 3^(-r) + 4^(-r)", 0);
}

/* ---- HarmonicNumber[Infinity, r] = Zeta[r] -------------------------- */

void test_hn_infinity() {
    assert_eval_eq("HarmonicNumber[Infinity, 2]", "1/6 Pi^2", 0);
    assert_eval_eq("HarmonicNumber[Infinity, 4]", "1/90 Pi^4", 0);
    assert_eval_eq("HarmonicNumber[Infinity, 3]", "Zeta[3]", 0);
    assert_eval_eq("HarmonicNumber[Infinity]", "ComplexInfinity", 0);
}

/* ---- non-positive integer order: Faulhaber polynomial --------------- */

void test_hn_faulhaber() {
    assert_eval_eq("HarmonicNumber[n, 0]", "n", 0);
    assert_eval_eq("HarmonicNumber[n, -1]", "1/2 n + 1/2 n^2", 0);
    assert_eval_eq("HarmonicNumber[z, -4]",
        "-1/30 z + 1/3 z^3 + 1/2 z^4 + 1/5 z^5", 0);
    /* polynomial agrees with the direct sum at a positive integer */
    assert_eval_eq("HarmonicNumber[5, -2]", "55", 0);   /* 1+4+9+16+25 */
}

/* ---- machine-precision numerics ------------------------------------- */

void test_hn_machine() {
    assert_close("HarmonicNumber[.8, 3]", 0.940124, 1e-5);
    assert_close("HarmonicNumber[0.33]", 0.441524, 1e-5);
    /* inexact order forces numericization of a Constant n (E) */
    assert_close("HarmonicNumber[E, 1.]", 1.750019, 1e-5);
    assert_close("HarmonicNumber[3.0]", 11.0 / 6.0, 1e-9);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_hn_highprec() {
    /* precision tracks the input precision (>16 sig figs) */
    assert_eval_startswith("HarmonicNumber[0.33000000000000000000]",
                           "0.441523646937363528");
    /* N[HarmonicNumber[1/17, 5], 50] */
    assert_eval_startswith("N[HarmonicNumber[1/17, 5], 50]",
                           "0.253276152061187075210346261187542283134331408859");
    /* N[HarmonicNumber[13/7, 5], 100] */
    assert_eval_startswith("N[HarmonicNumber[13/7, 5], 100]",
                           "1.02983927571663964113052863402434138431036195317294");
}

/* ---- complex order -------------------------------------------------- */

void test_hn_complex() {
    assert_complex_close("N[HarmonicNumber[27, 5 - I]]", 1.02598, 0.0251513, 1e-4);
}

/* ---- symbolic passthrough ------------------------------------------- */

void test_hn_symbolic() {
    assert_eval_eq("HarmonicNumber[n]", "HarmonicNumber[n]", 0);
    assert_eval_eq("HarmonicNumber[n, 2]", "HarmonicNumber[n, 2]", 0);
    /* a generic free symbol must not be numericized by an inexact order */
    assert_eval_eq("HarmonicNumber[x, 2.5]", "HarmonicNumber[x, 2.5]", 0);
    /* exact non-integer n stays symbolic until N */
    assert_eval_eq("HarmonicNumber[1/17, 5]", "HarmonicNumber[1/17, 5]", 0);
    /* negative integer n is not a finite sum -- stays symbolic */
    assert_eval_eq("HarmonicNumber[-2, 3]", "HarmonicNumber[-2, 3]", 0);
}

/* ---- derivatives ---------------------------------------------------- */

void test_hn_derivatives() {
    /* d/dn H_n = Zeta[2] - HarmonicNumber[n, 2] = Pi^2/6 - HarmonicNumber[n,2]. */
    assert_eval_eq("D[HarmonicNumber[n], n]",
                   "1/6 Pi^2 - HarmonicNumber[n, 2]", 0);
    /* d/dn H_n^(r) = r (Zeta[r+1] - HarmonicNumber[n, r+1]). */
    assert_eval_eq("D[HarmonicNumber[n, r], n]",
                   "r (Zeta[1 + r] - HarmonicNumber[n, 1 + r])", 0);
    /* d/dr has no elementary closed form: generic partial. */
    assert_eval_eq("D[HarmonicNumber[n, r], r]",
                   "Derivative[0, 1][HarmonicNumber][n, r]", 0);
    /* chain rule through the index */
    assert_eval_eq("D[HarmonicNumber[n^2], n]",
                   "2 n (1/6 Pi^2 - HarmonicNumber[n^2, 2])", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_hn_listable() {
    assert_eval_eq("HarmonicNumber[{1, 2, 3, 4}]", "{1, 3/2, 11/6, 25/12}", 0);
    assert_eval_eq("HarmonicNumber[{2, 3}, 2]", "{5/4, 49/36}", 0);
}

/* ---- argument-count errors ------------------------------------------ */

void test_hn_argt() {
    /* 0 or >2 args: emit HarmonicNumber::argt and stay unevaluated */
    assert_eval_eq("HarmonicNumber[]", "HarmonicNumber[]", 0);
    assert_eval_eq("HarmonicNumber[1, 2, 3]", "HarmonicNumber[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_hn_attributes() {
    SymbolDef* d = symtab_get_def("HarmonicNumber");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0,
               "HarmonicNumber must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "HarmonicNumber must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0,
               "HarmonicNumber must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_hn_integer_order1);
    TEST(test_hn_integer_orderr);
    TEST(test_hn_infinity);
    TEST(test_hn_faulhaber);
    TEST(test_hn_machine);
    TEST(test_hn_highprec);
    TEST(test_hn_complex);
    TEST(test_hn_symbolic);
    TEST(test_hn_derivatives);
    TEST(test_hn_listable);
    TEST(test_hn_argt);
    TEST(test_hn_attributes);

    printf("All HarmonicNumber tests passed.\n");
    return 0;
}
