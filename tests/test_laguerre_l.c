/* Tests for LaguerreL: the (generalized) Laguerre polynomial.
 *
 *   LaguerreL[n, x]      L_n(x)
 *   LaguerreL[n, a, x]   L_n^(a)(x)   (generalized; L_n[n,x] == L_n[n,0,x])
 *
 * Covers the exact integer polynomials (the L_5 / L_7 headline cases and the
 * exact rational table), the generalized symbolic-parameter polynomial, special
 * and zero values, Listable threading, machine and arbitrary-precision (MPFR)
 * numerics (real and complex, 2- and 3-argument, including the non-integer /
 * complex order via the confluent-hypergeometric form), the stable-recurrence
 * path that avoids the polynomial cancellation for a large integer order, the
 * D[] derivative rules, cross-surface identity (scalar == List == Compile),
 * symbolic fall-through, argument-count diagnostics, and attributes. */

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

/* ---- exact integer polynomials -------------------------------------- */

void test_laguerre_l_polynomials() {
    assert_eval_eq("LaguerreL[0, x]", "1", 0);
    assert_eval_eq("LaguerreL[1, x]", "1 - x", 0);
    assert_eval_eq("LaguerreL[2, x]", "1 - 2 x + 1/2 x^2", 0);
    /* Headline examples from the Wolfram reference. */
    assert_eval_eq("LaguerreL[5, x]",
                   "1 - 5 x + 5 x^2 - 5/3 x^3 + 5/24 x^4 - 1/120 x^5", 0);
    assert_eval_eq("LaguerreL[7, x]",
                   "1 - 7 x + 21/2 x^2 - 35/6 x^3 + 35/24 x^4 - 7/40 x^5 "
                   "+ 7/720 x^6 - 1/5040 x^7", 0);
}

/* ---- generalized (3-argument) polynomial ---------------------------- */

void test_laguerre_l_generalized() {
    /* Symbolic parameter a -> polynomial in a and x. */
    assert_eval_eq("LaguerreL[2, a, x]",
                   "1 + 3/2 a + 1/2 a^2 - 2 x - a x + 1/2 x^2", 0);
    /* The n/a grid from the reference. */
    assert_eval_eq("Table[LaguerreL[n, a, x], {n, 0, 3}, {a, 0, 1}]",
                   "{{1, 1}, {1 - x, 2 - x}, "
                   "{1 - 2 x + 1/2 x^2, 3 - 3 x + 1/2 x^2}, "
                   "{1 - 3 x + 3/2 x^2 - 1/6 x^3, 4 - 6 x + 2 x^2 - 1/6 x^3}}", 0);
    /* LaguerreL[n, x] == LaguerreL[n, 0, x]. */
    assert_eval_eq("LaguerreL[4, 0, x] == LaguerreL[4, x]", "True", 0);
}

/* ---- exact values at exact arguments -------------------------------- */

void test_laguerre_l_exact_values() {
    assert_eval_eq("LaguerreL[2, 3]", "-1/2", 0);        /* 1 - 6 + 9/2 */
    assert_eval_eq("LaguerreL[2, 3, 1]", "11/2", 0);     /* (a+1)(a+2)/2 - (a+2)x + x^2/2 */
    /* The exact rational table from the reference. */
    assert_eval_eq("Table[LaguerreL[10, x], {x, 1, 5}]",
                   "{168919/403200, -4381/14175, -31361/44800, "
                   "931/675, 254927/145152}", 0);
}

void test_laguerre_l_zero_values() {
    assert_eval_eq("LaguerreL[0, 0]", "1", 0);
    assert_eval_eq("LaguerreL[0, 0, 0]", "1", 0);
    assert_eval_eq("LaguerreL[n, 0]", "1", 0);           /* L_n(0) = 1 */
    assert_eval_eq("LaguerreL[5, 0]", "1", 0);
    assert_eval_eq("LaguerreL[2, 3, 0]", "10", 0);       /* Binomial[5, 2] */
}

/* ---- Listable threading --------------------------------------------- */

void test_laguerre_l_listable() {
    assert_eval_eq("LaguerreL[3, {1, 2, 3}]", "{-2/3, -1/3, 1}", 0);
    assert_eval_eq("LaguerreL[{0, 1, 2}, x]", "{1, 1 - x, 1 - 2 x + 1/2 x^2}", 0);
    assert_eval_eq("LaguerreL[2, {}]", "{}", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_laguerre_l_machine() {
    /* Inexact (integer-valued) order, exact argument. */
    assert_close("LaguerreL[6., 5]", -2.0902777777777777, 1e-4);
    /* Generalized, all machine arguments. */
    assert_close("LaguerreL[3, .51, .87]", -0.396576, 1e-5);
    /* An inexact integer-valued order yields a clean real (no spurious Im). */
    assert_close("LaguerreL[2., 3.]", -0.5, 1e-9);
    /* Cancellation case: the stable recurrence gives ~1.45553, NOT the
     * ~1.45577 that substituting into the explicit polynomial produces. */
    assert_close("LaguerreL[100, 0, 5.0]", 1.4555271625328812, 1e-6);
}

void test_laguerre_l_complex() {
    /* Complex (non-integer) order, 2-argument. */
    assert_close("Re[LaguerreL[6. - I, 5]]", -2.76039, 1e-4);
    assert_close("Im[LaguerreL[6. - I, 5]]", -2.46853, 1e-4);
    /* Complex argument, generalized. */
    assert_close("Re[LaguerreL[2.5, 1.3, 8 + I]]", 8.50391, 1e-4);
    assert_close("Im[LaguerreL[2.5, 1.3, 8 + I]]", 0.576867, 1e-4);
}

/* ---- arbitrary-precision numerics ----------------------------------- */

void test_laguerre_l_highprec() {
    /* Generalized, exact rational arguments, via N. */
    assert_eval_startswith(
        "N[LaguerreL[2/5, 1/3, 8/7], 50]",
        "0.68821856392835285644553764124515707164455717484011");
    /* Precision tracks the input precision. */
    assert_eval_startswith("LaguerreL[2.3000000000000000000000000, 3, 3]",
                           "-1.11292980686726510314");
    /* The cancellation case at 20-digit precision: stable recurrence. */
    assert_eval_startswith("LaguerreL[100, 0, N[5, 20]]",
                           "1.4555271625328812606");
}

/* ---- cross-surface identity ----------------------------------------- */

void test_laguerre_l_surfaces() {
    /* Scalar == element of the Listable (packed) result. */
    assert_eval_eq("LaguerreL[7, 2.5] == First[LaguerreL[7, {2.5}]]", "True", 0);
    /* Scalar == Compile[] output, bit-for-bit. */
    assert_eval_eq(
        "LaguerreL[7, 2.5] == Compile[{{x, _Real}}, LaguerreL[7, x]][2.5]",
        "True", 0);
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_laguerre_l_symbolic() {
    assert_eval_eq("LaguerreL[n, x]", "LaguerreL[n, x]", 0);
    assert_eval_eq("LaguerreL[n, a, x]", "LaguerreL[n, a, x]", 0);
    /* Exact non-integer order with a symbolic argument stays symbolic. */
    assert_eval_eq("LaguerreL[1/2, x]", "LaguerreL[1/2, x]", 0);
}

/* ---- derivatives ---------------------------------------------------- */

void test_laguerre_l_derivative() {
    assert_eval_eq("D[LaguerreL[n, x], x]", "-LaguerreL[-1 + n, 1, x]", 0);
    assert_eval_eq("D[LaguerreL[n, a, x], x]", "-LaguerreL[-1 + n, 1 + a, x]", 0);
    /* The derivative of an explicit polynomial is the polynomial's derivative. */
    assert_eval_eq("D[LaguerreL[2, x], x]", "-2 + x", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_laguerre_l_argcount() {
    /* Wrong arg counts stay unevaluated (an argt message goes to the funnel). */
    assert_eval_eq("LaguerreL[]", "LaguerreL[]", 0);
    assert_eval_eq("LaguerreL[3]", "LaguerreL[3]", 0);
    assert_eval_eq("LaguerreL[1, 2, 3, 4]", "LaguerreL[1, 2, 3, 4]", 0);
}

void test_laguerre_l_attributes() {
    SymbolDef* d = symtab_get_def("LaguerreL");
    ASSERT_MSG(d != NULL, "LaguerreL not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "LaguerreL not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "LaguerreL not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "LaguerreL not Protected");
}

/* The explicit polynomial agrees with the sum form
 *   L_n^(a)(x) = Sum_k (-1)^k Binomial[n + a, n - k] x^k / k!
 * (with Binomial[n + a, n - k] = Pochhammer[k + a + 1, n - k]/(n - k)! for a
 * rational a) on both construction paths: the O(n) coefficient ratio (a = 0, a rational,
 * a a negative integer below -n) and the recurrence kept for a negative
 * integer in [-n, -1]; and D of the degree-1000 expansion is exactly
 * -L_999^(1) at a rational point. */
static void test_laguerre_l_explicit_closed_form(void) {
    const char* cases[] = {
        "Expand[LaguerreL[300, x] - Sum[(-1)^k Binomial[300, k] x^k/k!, {k, 0, 300}]] === 0",
        /* Binomial[n + a, n - k] = Pochhammer[k + a + 1, n - k]/(n - k)! */
        "Expand[LaguerreL[40, 7/3, x] - Sum[(-1)^k Pochhammer[k + 7/3 + 1, 40 - k]/(40 - k)! x^k/k!, {k, 0, 40}]] === 0",
        "Expand[LaguerreL[30, -45, x] - Sum[(-1)^k Binomial[-15, 30 - k] x^k/k!, {k, 0, 30}]] === 0",
        "Expand[LaguerreL[30, -12, x] - Sum[(-1)^k Binomial[18, 30 - k] x^k/k!, {k, 0, 30}]] === 0",
        "(D[LaguerreL[1000, x], x] /. x -> 7/2) === -LaguerreL[999, 1, 7/2]",
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        Expr* e = parse_expression(cases[i]);
        Expr* r = evaluate(e);
        char* got = expr_to_string_fullform(r);
        if (strcmp(got, "True") != 0)
            fprintf(stdout, "FAIL: %s\n  got: %s\n", cases[i], got);
        ASSERT_STR_EQ(got, "True");
        free(got);
        expr_free(e);
        expr_free(r);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_laguerre_l_polynomials);
    TEST(test_laguerre_l_generalized);
    TEST(test_laguerre_l_exact_values);
    TEST(test_laguerre_l_zero_values);
    TEST(test_laguerre_l_listable);
    TEST(test_laguerre_l_machine);
    TEST(test_laguerre_l_complex);
    TEST(test_laguerre_l_highprec);
    TEST(test_laguerre_l_surfaces);
    TEST(test_laguerre_l_symbolic);
    TEST(test_laguerre_l_derivative);
    TEST(test_laguerre_l_explicit_closed_form);
    TEST(test_laguerre_l_argcount);
    TEST(test_laguerre_l_attributes);

    printf("All LaguerreL tests passed.\n");
    return 0;
}
