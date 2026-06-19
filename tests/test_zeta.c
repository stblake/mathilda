/* Tests for the zeta function family: Zeta[s] (Riemann) and Zeta[s, a]
 * (Hurwitz).
 *
 * Covers exact integer reduction (rational * Pi^(2n), negative-integer
 * rationals, the pole at 1), exact Hurwitz at integer a, machine real &
 * complex numerics, arbitrary-precision (MPFR) reals and complexes, precision
 * tracking, derivative rules, the Series expansions about s = 1 and s = 0,
 * symbolic fall-through, Listable threading, and attributes. */

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

/* ---- exact Riemann values ------------------------------------------- */

void test_zeta_exact_even() {
    /* Even positive integers: rational multiples of Pi^(2n). */
    assert_eval_eq("Zeta[2]", "1/6 Pi^2", 0);
    assert_eval_eq("Zeta[4]", "1/90 Pi^4", 0);
    assert_eval_eq("Zeta[6]", "1/945 Pi^6", 0);
    assert_eval_eq("Zeta[8]", "1/9450 Pi^8", 0);
}

void test_zeta_exact_nonpositive() {
    assert_eval_eq("Zeta[0]", "-1/2", 0);
    assert_eval_eq("Zeta[-1]", "-1/12", 0);
    assert_eval_eq("Zeta[-3]", "1/120", 0);
    assert_eval_eq("Zeta[-5]", "-1/252", 0);
    /* Negative even integers are trivial zeros. */
    assert_eval_eq("Zeta[-2]", "0", 0);
    assert_eval_eq("Zeta[-4]", "0", 0);
    /* The pole. */
    assert_eval_eq("Zeta[1]", "ComplexInfinity", 0);
}

void test_zeta_odd_symbolic() {
    /* Odd positive integers have no closed form -> stay symbolic. */
    assert_eval_eq("Zeta[3]", "Zeta[3]", 0);
    assert_eval_eq("Zeta[5]", "Zeta[5]", 0);
    assert_eval_eq("Zeta[7]", "Zeta[7]", 0);
}

void test_zeta_limits_symbolic() {
    assert_eval_eq("Zeta[Infinity]", "1", 0);
    assert_eval_eq("Zeta[s]", "Zeta[s]", 0);
    assert_eval_eq("Zeta[1/2]", "Zeta[1/2]", 0);
    assert_eval_eq("Zeta[s, a]", "Zeta[s, a]", 0);
}

/* ---- exact Hurwitz values ------------------------------------------- */

void test_zeta_hurwitz_exact() {
    /* zeta(s, 1) = zeta(s). */
    assert_eval_eq("Zeta[s, 1]", "Zeta[s]", 0);
    assert_eval_eq("Zeta[2, 1]", "1/6 Pi^2", 0);
    /* Positive integer a: zeta(s, m) = zeta(s) - Sum_{k=1}^{m-1} k^-s. */
    assert_eval_eq("Zeta[3, 2]", "-1 + Zeta[3]", 0);
    assert_eval_eq("Zeta[2, 2]", "-1 + 1/6 Pi^2", 0);
    assert_eval_eq("Zeta[2, 3]", "-5/4 + 1/6 Pi^2", 0);
    assert_eval_eq("Zeta[4, 5]", "-22369/20736 + 1/90 Pi^4", 0);
    /* The pole carries through: zeta(1, 2) = zeta(1) - 1 = ComplexInfinity. */
    assert_eval_eq("Zeta[1, 2]", "ComplexInfinity", 0);
    /* a = 1/2: Zeta[s, 1/2] = (2^s - 1) Zeta[s]. */
    assert_eval_eq("Zeta[s, 1/2]", "(-1 + 2^s) Zeta[s]", 0);
    assert_eval_eq("Zeta[2, 1/2]", "1/2 Pi^2", 0);   /* 3 Zeta[2] */
    assert_eval_eq("Zeta[3, 1/2]", "7 Zeta[3]", 0);
    assert_eval_eq("Zeta[-2, 1/2]", "0", 0);          /* Zeta[-2] = 0 */
}

/* ---- machine-precision numerics ------------------------------------- */

void test_zeta_machine_real() {
    assert_close("Zeta[2.]", 1.6449340668482264, 1e-12);   /* Pi^2/6 */
    assert_close("Zeta[4.]", 1.0823232337111382, 1e-12);   /* Pi^4/90 */
    assert_close("Zeta[3.]", 1.2020569031595943, 1e-12);   /* Apery */
    assert_close("Zeta[0.5]", -1.4603545088095868, 1e-12);
    assert_close("Zeta[0.]", -0.5, 1e-12);
    assert_close("Zeta[-1.]", -1.0 / 12.0, 1e-10);
    /* Hurwitz real: zeta(-1, 5) = -B_2(5)/2 = -(25-5+1/6)/2. */
    assert_close("Zeta[-1., 5]", -10.083333333333334, 1e-7);
    /* Hurwitz real: zeta(3, 2) = zeta(3) - 1. */
    assert_close("Zeta[3., 2]", 0.20205690315959429, 1e-12);
    /* Negative a uses the symmetric ((a+k)^2)^(-s/2) convention (|a+k|^-s on the
     * real axis): the left-half-plane head terms ADD rather than cancel, so this
     * differs from the principal-branch HurwitzZeta value. */
    assert_close("Zeta[3, -1.5]", 16.710694618413456, 1e-6);
    assert_close("Zeta[3, -3.5]", 16.798018233573806, 1e-6);
    assert_close("Zeta[2.5, -0.5]", 11.903964884061192, 1e-6);
    /* HurwitzZeta keeps the principal branch (negative terms cancel pairwise). */
    assert_close("HurwitzZeta[3, -1.5]", 0.1181020258, 1e-6);
}

void test_zeta_machine_complex() {
    /* N[Zeta[1/2 + 1/2 I]] = -0.459303 - 0.961254 I. */
    assert_complex_close("N[Zeta[1/2 + 1/2 I]]", -0.459303, -0.961254, 1e-5);
    /* Complex Hurwitz: Zeta[-1.5 + I, 2.5 - I] = 0.0184868 + 1.67553 I. */
    assert_complex_close("Zeta[-1.5 + I, 2.5 - I]", 0.0184868, 1.67553, 1e-4);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_zeta_arbitrary_precision() {
    /* N[Zeta[3], 50] -- Apery's constant. */
    assert_eval_startswith("N[Zeta[3], 50]",
                           "1.2020569031595942853997381615114499907649862923405");
    /* N[Zeta[5/4], 50] -- from the reference. */
    assert_eval_startswith("N[Zeta[5/4], 50]",
                           "4.5951118258429433806853780396946256522810297806048");
    /* N[Zeta[2], 50] -- Pi^2/6. */
    assert_eval_startswith("N[Zeta[2], 50]",
                           "1.6449340668482264364724151666460251892189499012067");
    /* Arbitrary-precision Hurwitz: N[Zeta[3, 2], 50] = zeta(3) - 1. */
    assert_eval_startswith("N[Zeta[3, 2], 50]",
                           "0.20205690315959428539973816151144999076498629234");
    /* Precision tracks the input precision. */
    assert_eval_startswith("Zeta[2`40]", "1.6449340668482264364724151666460251892");
}

void test_zeta_arbitrary_complex() {
    /* Arbitrary-precision complex zeta. N[Zeta[1/2 + I/2], 30] -- real part;
     * agrees with the machine path and the reference -0.459303 - 0.961254 I. */
    assert_eval_startswith("N[Zeta[1/2 + 1/2 I], 30]", "-0.45930289");
}

/* ---- derivatives ---------------------------------------------------- */

void test_zeta_derivatives() {
    /* d/da Zeta[s, a] = -s Zeta[1+s, a]. */
    assert_eval_eq("D[Zeta[s, a], a]", "-s Zeta[1 + s, a]", 0);
    /* d/ds has no elementary closed form: generic partial. */
    assert_eval_eq("D[Zeta[s, a], s]", "Derivative[1, 0][Zeta][s, a]", 0);
    /* One-argument derivative: generic. */
    assert_eval_eq("D[Zeta[x], x]", "Derivative[1][Zeta][x]", 0);
    /* Chain rule through the second argument. */
    assert_eval_eq("D[Zeta[s, x^2], x]", "-2 s x Zeta[1 + s, x^2]", 0);
    /* Higher a-derivatives match the rising-factorial pattern. */
    assert_eval_eq("Table[D[Zeta[s, a], {a, k}], {k, 1, 3}]",
                   "{-s Zeta[1 + s, a], s (1 + s) Zeta[2 + s, a], "
                   "-s (1 + s) (2 + s) Zeta[3 + s, a]}", 0);
}

/* ---- Series expansions ---------------------------------------------- */

void test_zeta_series() {
    /* Laurent expansion about s = 1 introduces the Stieltjes constants. */
    assert_eval_eq("Series[Zeta[x], {x, 1, 2}] // Normal",
                   "EulerGamma + 1/(-1 + x) - StieltjesGamma[1] (-1 + x) + "
                   "1/2 StieltjesGamma[2] (-1 + x)^2", 0);
    /* Taylor expansion about s = 0: zeta(0) = -1/2, zeta'(0) = -1/2 Log[2 Pi]. */
    assert_eval_eq("Series[Zeta[x], {x, 0, 1}] // Normal",
                   "-1/2 - 1/2 Log[2 Pi] x", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_zeta_listable() {
    assert_eval_eq("Zeta[{2, 4, 6}]", "{1/6 Pi^2, 1/90 Pi^4, 1/945 Pi^6}", 0);
    assert_eval_eq("Zeta[{0, -1, -3}]", "{-1/2, -1/12, 1/120}", 0);
}

void test_zeta_attributes() {
    SymbolDef* d = symtab_get_def("Zeta");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Zeta must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "Zeta must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Zeta must be Protected");
}

void test_zeta_argcount() {
    /* Wrong argument counts emit Zeta::argt and stay unevaluated. */
    assert_eval_eq("Zeta[]", "Zeta[]", 0);
    assert_eval_eq("Zeta[1, 2, 3]", "Zeta[1, 2, 3]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_zeta_exact_even);
    TEST(test_zeta_exact_nonpositive);
    TEST(test_zeta_odd_symbolic);
    TEST(test_zeta_limits_symbolic);
    TEST(test_zeta_hurwitz_exact);
    TEST(test_zeta_machine_real);
    TEST(test_zeta_machine_complex);
    TEST(test_zeta_arbitrary_precision);
    TEST(test_zeta_arbitrary_complex);
    TEST(test_zeta_derivatives);
    TEST(test_zeta_series);
    TEST(test_zeta_listable);
    TEST(test_zeta_attributes);
    TEST(test_zeta_argcount);

    printf("All Zeta tests passed.\n");
    return 0;
}
