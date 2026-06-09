/* Tests for EulerE: Euler numbers EulerE[n] and Euler polynomials
 * EulerE[n, x].
 *
 * Covers exact integer numbers (including odd-index zeros and a big-integer
 * case), the polynomial expansion, structural identities, the EulerE[n, 1/2]
 * symbolic fold, Listable threading (including empty lists), machine and
 * arbitrary-precision (MPFR) numerics, symbolic fall-through for non-integer /
 * negative / symbolic arguments, argument-count diagnostics, and attributes. */

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

/* ---- exact Euler numbers -------------------------------------------- */

void test_euler_numbers_small() {
    /* The first eleven Euler numbers (E_0 .. E_10). */
    assert_eval_eq("Table[EulerE[k], {k, 0, 10}]",
                   "{1, 0, -1, 0, 5, 0, -61, 0, 1385, 0, -50521}", 0);
    /* Individual values. */
    assert_eval_eq("EulerE[0]",  "1",      0);
    assert_eval_eq("EulerE[2]",  "-1",     0);
    assert_eval_eq("EulerE[4]",  "5",      0);
    assert_eval_eq("EulerE[6]",  "-61",    0);
    assert_eval_eq("EulerE[8]",  "1385",   0);
    assert_eval_eq("EulerE[10]", "-50521", 0);
}

void test_euler_numbers_odd_zero() {
    /* All odd indices vanish exactly. */
    assert_eval_eq("EulerE[1]",  "0", 0);
    assert_eval_eq("EulerE[3]",  "0", 0);
    assert_eval_eq("EulerE[5]",  "0", 0);
    assert_eval_eq("EulerE[99]", "0", 0);
}

void test_euler_numbers_bigint() {
    /* Large index exercises the GMP big-integer path. */
    assert_eval_eq("EulerE[20]", "370371188237525", 0);
}

/* ---- Euler polynomials ---------------------------------------------- */

void test_euler_polynomials() {
    /* E_n(z) for n = 0 .. 5. */
    assert_eval_eq("Table[EulerE[n, z], {n, 0, 5}]",
                   "{1, -1/2 + z, -z + z^2, 1/4 - 3/2 z^2 + z^3, "
                   "z - 2 z^3 + z^4, -1/2 + 5/2 z^2 - 5/2 z^4 + z^5}",
                   0);
    assert_eval_eq("EulerE[0, z]", "1", 0);
    assert_eval_eq("EulerE[2, z]", "-z + z^2", 0);
    assert_eval_eq("EulerE[6, x]",
                   "-3 x + 5 x^3 - 3 x^5 + x^6", 0);
}

void test_euler_polynomial_number_relation() {
    /* E_n = 2^n E_n(1/2). */
    assert_eval_eq("Table[2^n EulerE[n, 1/2] - EulerE[n], {n, 0, 8}]",
                   "{0, 0, 0, 0, 0, 0, 0, 0, 0}", 0);
    /* The two-argument form at 1/2 gives the exact rational. */
    assert_eval_eq("EulerE[2, 1/2]", "-1/4", 0);
    assert_eval_eq("EulerE[4, 1/2]", "5/16", 0);
}

void test_euler_polynomial_identity() {
    /* Difference identity E_n(x+1) + E_n(x) = 2 x^n. */
    assert_eval_eq("Expand[EulerE[2, x + 1] + EulerE[2, x]]", "2 x^2", 0);
    assert_eval_eq("Expand[EulerE[5, x + 1] + EulerE[5, x]]", "2 x^5", 0);
    /* Derivative identity E_n'(x) = n E_{n-1}(x). */
    assert_eval_eq("Expand[D[EulerE[4, x], x] - 4 EulerE[3, x]]", "0", 0);
}

void test_euler_half_symbolic() {
    /* EulerE[n, 1/2] folds to 2^-n EulerE[n] for symbolic n. */
    assert_eval_eq("EulerE[n, 1/2]", "EulerE[n] 2^(-n)", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_euler_listable() {
    assert_eval_eq("EulerE[{2, 4, 6}]", "{-1, 5, -61}", 0);
    /* Two-argument form threads over a list of orders, broadcasting x. */
    assert_eval_eq("EulerE[{2, 4, 6}, x]",
                   "{-x + x^2, x - 2 x^3 + x^4, -3 x + 5 x^3 - 3 x^5 + x^6}", 0);
}

void test_euler_empty_list() {
    /* Listable threading over an empty list yields an empty list. */
    assert_eval_eq("EulerE[{}]", "{}", 0);
    assert_eval_eq("EulerE[{}, x]", "{}", 0);
}

/* ---- numeric evaluation --------------------------------------------- */

void test_euler_machine_numeric() {
    /* Inexact integer-valued order -> machine real. */
    assert_close("EulerE[2.0]", -1.0, 1e-12);
    assert_close("EulerE[6.0]", -61.0, 1e-12);
    /* Inexact polynomial argument evaluates numerically: E_2(0.5) = -1/4. */
    assert_close("EulerE[2, 0.5]", -0.25, 1e-12);
}

void test_euler_arbitrary_precision() {
    /* Polynomial at an exact rational, numericalised to 30 digits:
     * E_6(1/3) = -602/729 = -0.8257887517146776406035665294924... */
    assert_eval_startswith("N[EulerE[6, 1/3], 30]",
                           "-0.82578875171467764060356652949");
    /* MPFR-precision integer order numericalises the exact integer. */
    assert_eval_startswith("N[EulerE[10, N[1, 40]], 40]", "0.");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_euler_symbolic() {
    /* Symbolic, negative-integer and non-integer arguments stay unevaluated. */
    assert_eval_eq("EulerE[n]",    "EulerE[n]",    0);
    assert_eval_eq("EulerE[-1]",   "EulerE[-1]",   0);
    assert_eval_eq("EulerE[3/2]",  "EulerE[3/2]",  0);
    assert_eval_eq("EulerE[n, x]", "EulerE[n, x]", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_euler_argcount() {
    /* Wrong argument counts emit EulerE::argt and stay unevaluated. */
    assert_eval_eq("EulerE[]",            "EulerE[]",            0);
    assert_eval_eq("EulerE[1, 2, 3, x]",  "EulerE[1, 2, 3, x]",  0);
}

void test_euler_attributes() {
    SymbolDef* d = symtab_get_def("EulerE");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "EulerE must be Listable");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "EulerE must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_euler_numbers_small);
    TEST(test_euler_numbers_odd_zero);
    TEST(test_euler_numbers_bigint);
    TEST(test_euler_polynomials);
    TEST(test_euler_polynomial_number_relation);
    TEST(test_euler_polynomial_identity);
    TEST(test_euler_half_symbolic);
    TEST(test_euler_listable);
    TEST(test_euler_empty_list);
    TEST(test_euler_machine_numeric);
    TEST(test_euler_arbitrary_precision);
    TEST(test_euler_symbolic);
    TEST(test_euler_argcount);
    TEST(test_euler_attributes);

    printf("All EulerE tests passed.\n");
    return 0;
}
