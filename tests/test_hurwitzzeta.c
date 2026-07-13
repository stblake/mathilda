/* Tests for the Hurwitz zeta function HurwitzZeta[s, a].
 *
 * Covers the exact reductions (a = 1 -> Zeta[s], a = 1/2 -> (2^s-1) Zeta[s],
 * positive integer a -> Zeta[s] minus a finite power sum), the poles at
 * non-positive integer a, the Bernoulli-polynomial values for non-positive
 * integer s, machine real & complex numerics, arbitrary-precision (MPFR) reals
 * and complexes, precision tracking, the principal-branch behaviour that
 * distinguishes HurwitzZeta from Zeta on the negative axis, the derivative
 * rules, symbolic fall-through, Listable threading, attributes and argument
 * counting. */

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

/* ---- exact reductions ----------------------------------------------- */

void test_hz_reduce_to_zeta() {
    /* HurwitzZeta[s, 1] = Zeta[s]. */
    assert_eval_eq("HurwitzZeta[s, 1]", "Zeta[s]", 0);
    assert_eval_eq("HurwitzZeta[3, 1]", "Zeta[3]", 0);
    assert_eval_eq("HurwitzZeta[2, 1]", "1/6 Pi^2", 0);
    assert_eval_eq("HurwitzZeta[0, 1]", "-1/2", 0);
    assert_eval_eq("HurwitzZeta[-1, 1]", "-1/12", 0);
    assert_eval_eq("HurwitzZeta[-2, 1]", "0", 0);
    /* Table[HurwitzZeta[s,1],{s,-2,2}] from the reference. */
    assert_eval_eq("Table[HurwitzZeta[s, 1], {s, -2, 2}]",
                   "{0, -1/12, -1/2, ComplexInfinity, 1/6 Pi^2}", 0);
}

void test_hz_half() {
    /* HurwitzZeta[s, 1/2] = (2^s - 1) Zeta[s]. */
    assert_eval_eq("HurwitzZeta[s, 1/2]", "(-1 + 2^s) Zeta[s]", 0);
    /* HurwitzZeta[2, 1/2] = 3 Zeta[2] = Pi^2/2. */
    assert_eval_eq("HurwitzZeta[2, 1/2]", "1/2 Pi^2", 0);
    /* Negative even integer s: Zeta[s] = 0 -> 0. */
    assert_eval_eq("HurwitzZeta[-2, 1/2]", "0", 0);
}

void test_hz_positive_integer_a() {
    /* zeta(s, m) = zeta(s) - Sum_{k=1}^{m-1} k^-s. */
    assert_eval_eq("HurwitzZeta[3, 2]", "-1 + Zeta[3]", 0);
    assert_eval_eq("HurwitzZeta[2, 2]", "-1 + 1/6 Pi^2", 0);
    assert_eval_eq("HurwitzZeta[2, 3]", "-5/4 + 1/6 Pi^2", 0);
    assert_eval_eq("HurwitzZeta[2, 5]", "-205/144 + 1/6 Pi^2", 0);
    assert_eval_eq("HurwitzZeta[4, 5]", "-22369/20736 + 1/90 Pi^4", 0);
    /* The pole carries through. */
    assert_eval_eq("HurwitzZeta[1, 2]", "ComplexInfinity", 0);
}

/* ---- poles and Bernoulli-polynomial values at non-positive integer a -- */

void test_hz_nonpositive_integer_a() {
    /* Positive integer s at a non-positive integer a -> pole. */
    assert_eval_eq("HurwitzZeta[2, 0]", "ComplexInfinity", 0);
    assert_eval_eq("HurwitzZeta[2, -2]", "ComplexInfinity", 0);
    assert_eval_eq("HurwitzZeta[3, -1]", "ComplexInfinity", 0);
    /* s = 1 is a pole everywhere. */
    assert_eval_eq("HurwitzZeta[1, 5]", "ComplexInfinity", 0);
    assert_eval_eq("HurwitzZeta[1, a]", "ComplexInfinity", 0);
    /* Non-positive integer s -> finite Bernoulli-polynomial value. */
    assert_eval_eq("HurwitzZeta[0, 0]", "1/2", 0);    /* 1/2 - a at a = 0 */
    assert_eval_eq("HurwitzZeta[0, 3]", "-5/2", 0);   /* 1/2 - 3 */
    assert_eval_eq("HurwitzZeta[-1, 0]", "-1/12", 0);
    assert_eval_eq("HurwitzZeta[-2, 0]", "0", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_hz_machine_real() {
    assert_close("HurwitzZeta[3, 0.2]", 125.739018057217966531, 1e-9);
    assert_close("HurwitzZeta[7., 5]", 1.84948548452090944698e-05, 1e-12);
    assert_close("HurwitzZeta[0.51, 0.87]", -1.3201550236949583755, 1e-9);
    /* Real Hurwitz at a positive integer a, inexact s. */
    assert_close("HurwitzZeta[3., 2]", 0.20205690315959429, 1e-12);
}

void test_hz_negative_axis_branch() {
    /* HurwitzZeta uses the principal branch (k+a)^-s, unlike Zeta. For an
     * integer s = 3 and a = -7/2, every summand is a real signed cube. */
    assert_close("HurwitzZeta[3, -3.5]", 0.0307784106605138472745, 1e-9);
}

void test_hz_machine_complex() {
    /* HurwitzZeta[2.3, 8 + I] = 0.0544701 - 0.00944852 I. */
    assert_complex_close("HurwitzZeta[2.3, 8 + I]",
                         0.0544700827321306593516, -0.00944851533647024666501, 1e-6);
    /* Complex order and complex argument. */
    assert_complex_close("HurwitzZeta[-1.5 + I, 2.5 - I]",
                         0.0184868091938106521569, 1.67553378430176592096, 1e-6);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_hz_arbitrary_precision() {
    /* N[HurwitzZeta[1/3, 8/7], 50] from the reference. */
    assert_eval_startswith("N[HurwitzZeta[1/3, 8/7], 50]",
                           "-1.1389367444490991746548674334535727810961919460755");
    /* HurwitzZeta[3, 2.1`60] from the reference. */
    assert_eval_startswith("HurwitzZeta[3, 2.1`60]",
                           "0.17941388827042539592296477445767280409529613750741");
    /* Precision tracks the lower-precision (machine) argument. */
    assert_eval_startswith("HurwitzZeta[2.3000000000000000000000000, 48]",
                           "0.00508546861589645115101714");
    /* Arbitrary-precision Hurwitz at integer a: equals Zeta[3] - 1. */
    assert_eval_startswith("N[HurwitzZeta[3, 2], 50]",
                           "0.20205690315959428539973816151144999076498629234");
}

void test_hz_arbitrary_complex() {
    assert_eval_startswith("N[HurwitzZeta[23/10, 8 + I], 30]", "0.05447008273213065935");
}

/* ---- derivatives ---------------------------------------------------- */

void test_hz_derivatives() {
    /* d/da HurwitzZeta[s, a] = -s HurwitzZeta[1+s, a]. */
    assert_eval_eq("D[HurwitzZeta[s, a], a]", "-s HurwitzZeta[1 + s, a]", 0);
    /* d/ds has no elementary closed form: generic partial. */
    assert_eval_eq("D[HurwitzZeta[s, a], s]",
                   "Derivative[1, 0][HurwitzZeta][s, a]", 0);
    /* Chain rule through the second argument. */
    assert_eval_eq("D[HurwitzZeta[s, x^2], x]",
                   "-2 s x HurwitzZeta[1 + s, x^2]", 0);
    /* Higher a-derivatives match the rising-factorial pattern. */
    assert_eval_eq("Table[D[HurwitzZeta[s, a], {a, k}], {k, 1, 4}]",
                   "{-s HurwitzZeta[1 + s, a], s (1 + s) HurwitzZeta[2 + s, a], "
                   "-s (1 + s) (2 + s) HurwitzZeta[3 + s, a], "
                   "s (1 + s) (2 + s) (3 + s) HurwitzZeta[4 + s, a]}", 0);
    /* Symbolic-order derivative in a:
     * D[HurwitzZeta[s,a],{a,n}] = (-1)^n Pochhammer[s,n] HurwitzZeta[n+s,a]. */
    assert_eval_eq("D[HurwitzZeta[s, a], {a, n}]",
                   "(-1)^n Pochhammer[s, n] HurwitzZeta[n + s, a]", 0);
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_hz_symbolic() {
    assert_eval_eq("HurwitzZeta[s, a]", "HurwitzZeta[s, a]", 0);
    assert_eval_eq("HurwitzZeta[3, a]", "HurwitzZeta[3, a]", 0);
    /* Symbolic s at a non-positive integer a: indeterminate -> stays. */
    assert_eval_eq("HurwitzZeta[s, 0]", "HurwitzZeta[s, 0]", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_hz_listable() {
    assert_eval_eq("HurwitzZeta[{3, 2}, 1]", "{Zeta[3], 1/6 Pi^2}", 0);
    assert_eval_eq("HurwitzZeta[2, {2, 5}]", "{-1 + 1/6 Pi^2, -205/144 + 1/6 Pi^2}", 0);
}

void test_hz_attributes() {
    SymbolDef* d = symtab_get_def("HurwitzZeta");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "HurwitzZeta must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "HurwitzZeta must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "HurwitzZeta must be Protected");
}

void test_hz_argcount() {
    /* Wrong argument counts emit HurwitzZeta::argrx and stay unevaluated. */
    assert_eval_eq("HurwitzZeta[]", "HurwitzZeta[]", 0);
    assert_eval_eq("HurwitzZeta[3]", "HurwitzZeta[3]", 0);
    assert_eval_eq("HurwitzZeta[1, 2, 3]", "HurwitzZeta[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_hz_reduce_to_zeta);
    TEST(test_hz_half);
    TEST(test_hz_positive_integer_a);
    TEST(test_hz_nonpositive_integer_a);
    TEST(test_hz_machine_real);
    TEST(test_hz_negative_axis_branch);
    TEST(test_hz_machine_complex);
    TEST(test_hz_arbitrary_precision);
    TEST(test_hz_arbitrary_complex);
    TEST(test_hz_derivatives);
    TEST(test_hz_symbolic);
    TEST(test_hz_listable);
    TEST(test_hz_attributes);
    TEST(test_hz_argcount);

    printf("All HurwitzZeta tests passed.\n");
    return 0;
}
