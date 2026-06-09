/* Tests for BernoulliB: Bernoulli numbers BernoulliB[n] and Bernoulli
 * polynomials BernoulliB[n, x].
 *
 * Covers exact rational numbers (including odd-index zeros and a big-integer
 * case), the polynomial expansion and its B_n(0) = B_n boundary, Listable
 * threading, machine and arbitrary-precision (MPFR) numerics, symbolic
 * fall-through for non-integer / negative / symbolic arguments, argument-count
 * diagnostics, and attributes. */

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

/* ---- numeric helper -------------------------------------------------- */

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

/* ---- exact Bernoulli numbers ---------------------------------------- */

void test_bern_numbers_small() {
    /* The first eleven Bernoulli numbers (B_0 .. B_10). */
    assert_eval_eq("Table[BernoulliB[k], {k, 0, 10}]",
                   "{1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66}", 0);
    /* Individual values. */
    assert_eval_eq("BernoulliB[0]",  "1",       0);
    assert_eval_eq("BernoulliB[1]",  "-1/2",    0);
    assert_eval_eq("BernoulliB[2]",  "1/6",     0);
    assert_eval_eq("BernoulliB[12]", "-691/2730", 0);
    assert_eval_eq("BernoulliB[20]", "-174611/330", 0);
}

void test_bern_numbers_odd_zero() {
    /* All odd indices above 1 vanish exactly. */
    assert_eval_eq("BernoulliB[3]",  "0", 0);
    assert_eval_eq("BernoulliB[5]",  "0", 0);
    assert_eval_eq("BernoulliB[7]",  "0", 0);
    assert_eval_eq("BernoulliB[99]", "0", 0);
}

void test_bern_numbers_bigint() {
    /* Large index exercises the GMP big-rational path. */
    assert_eval_eq("BernoulliB[30]", "8615841276005/14322", 0);
}

/* ---- Bernoulli polynomials ------------------------------------------ */

void test_bern_polynomials() {
    /* B_n(z) for n = 0 .. 5. */
    assert_eval_eq("Table[BernoulliB[n, z], {n, 0, 5}]",
                   "{1, -1/2 + z, 1/6 - z + z^2, 1/2 z - 3/2 z^2 + z^3, "
                   "-1/30 + z^2 - 2 z^3 + z^4, -1/6 z + 5/3 z^3 - 5/2 z^4 + z^5}",
                   0);
    assert_eval_eq("BernoulliB[0, z]", "1", 0);
    assert_eval_eq("BernoulliB[2, z]", "1/6 - z + z^2", 0);
    assert_eval_eq("BernoulliB[6, x]",
                   "1/42 - 1/2 x^2 + 5/2 x^4 - 3 x^5 + x^6", 0);
}

void test_bern_polynomial_at_zero() {
    /* B_n(0) = B_n. */
    assert_eval_eq("BernoulliB[4, 0]", "-1/30", 0);
    assert_eval_eq("BernoulliB[5, 0]", "0",     0);
    assert_eval_eq("BernoulliB[2, 0]", "1/6",   0);
    /* The polynomial reduces to the plain number for every n. */
    assert_eval_eq("Table[BernoulliB[n, 0] - BernoulliB[n], {n, 0, 8}]",
                   "{0, 0, 0, 0, 0, 0, 0, 0, 0}", 0);
}

void test_bern_polynomial_identity() {
    /* Expand[B_n(x)] matches the closed form (sanity on the coefficients). */
    assert_eval_eq("Expand[BernoulliB[3, x] - (x^3 - 3/2 x^2 + 1/2 x)]", "0", 0);
    /* Difference identity B_n(x+1) - B_n(x) = n x^(n-1). */
    assert_eval_eq("Expand[BernoulliB[4, x + 1] - BernoulliB[4, x]]",
                   "4 x^3", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_bern_listable() {
    assert_eval_eq("BernoulliB[{2, 4, 6}]", "{1/6, -1/30, 1/42}", 0);
    /* Two-argument form threads over a list of orders. */
    assert_eval_eq("BernoulliB[{0, 1, 2}, z]",
                   "{1, -1/2 + z, 1/6 - z + z^2}", 0);
}

/* ---- numeric evaluation --------------------------------------------- */

void test_bern_machine_numeric() {
    /* Inexact integer-valued order -> machine real. */
    assert_close("BernoulliB[2.0]", 1.0 / 6.0, 1e-12);
    assert_close("BernoulliB[4.0]", -1.0 / 30.0, 1e-12);
    /* Inexact polynomial argument evaluates numerically. */
    assert_close("BernoulliB[2, 0.5]", 0.25 - 0.5 + 1.0 / 6.0, 1e-12);
}

void test_bern_arbitrary_precision() {
    /* BernoulliB can be evaluated to arbitrary numerical precision:
     * 691/2730 = 0.25311355311355... (repeating). */
    assert_eval_startswith("N[BernoulliB[12], 30]", "-0.25311355311355311355311355");
    /* Inexact MPFR-precision order also numericalises. */
    assert_eval_startswith("N[BernoulliB[2], 40]", "0.16666666666666666666");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_bern_symbolic() {
    /* Symbolic, negative-integer and non-integer arguments stay unevaluated. */
    assert_eval_eq("BernoulliB[n]",    "BernoulliB[n]",    0);
    assert_eval_eq("BernoulliB[-1]",   "BernoulliB[-1]",   0);
    assert_eval_eq("BernoulliB[3/2]",  "BernoulliB[3/2]",  0);
    assert_eval_eq("BernoulliB[n, x]", "BernoulliB[n, x]", 0);
}

void test_bern_empty_list() {
    /* Listable threading over an empty list yields an empty list. */
    assert_eval_eq("BernoulliB[{}]", "{}", 0);
    assert_eval_eq("BernoulliB[{}, z]", "{}", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_bern_argcount() {
    /* Wrong argument counts emit BernoulliB::argt and stay unevaluated. */
    assert_eval_eq("BernoulliB[]",        "BernoulliB[]",        0);
    assert_eval_eq("BernoulliB[4, 3, 2, 1]", "BernoulliB[4, 3, 2, 1]", 0);
}

void test_bern_attributes() {
    SymbolDef* d = symtab_get_def("BernoulliB");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "BernoulliB must be Listable");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "BernoulliB must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_bern_numbers_small);
    TEST(test_bern_numbers_odd_zero);
    TEST(test_bern_numbers_bigint);
    TEST(test_bern_polynomials);
    TEST(test_bern_polynomial_at_zero);
    TEST(test_bern_polynomial_identity);
    TEST(test_bern_listable);
    TEST(test_bern_machine_numeric);
    TEST(test_bern_arbitrary_precision);
    TEST(test_bern_symbolic);
    TEST(test_bern_empty_list);
    TEST(test_bern_argcount);
    TEST(test_bern_attributes);

    printf("All BernoulliB tests passed.\n");
    return 0;
}
