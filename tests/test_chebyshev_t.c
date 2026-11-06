/* Tests for ChebyshevT: the Chebyshev polynomial of the first kind T_n(x).
 *
 * Covers exact integer polynomials (including the T_10 / T_7 headline cases and
 * the even-order reflection T_{-n} = T_n), exact values at exact arguments, the
 * special values T_n(0) = Cos[n Pi/2] and T_n(1) = 1, Listable threading,
 * machine and arbitrary-precision (MPFR) numerics (real and complex, including
 * the non-integer/complex order via the 2F1 and ArcCos closed forms and the
 * branch cut), the stable-recurrence path that avoids polynomial cancellation
 * for large integer order, the D[] derivative rule (n U_{n-1}), symbolic
 * fall-through, argument-count diagnostics, and attributes. */

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

void test_chebyshev_t_polynomials() {
    assert_eval_eq("Table[ChebyshevT[n, x], {n, 0, 4}]",
                   "{1, x, -1 + 2 x^2, -3 x + 4 x^3, 1 - 8 x^2 + 8 x^4}", 0);
    assert_eval_eq("ChebyshevT[0, x]", "1", 0);
    assert_eval_eq("ChebyshevT[1, x]", "x", 0);
    assert_eval_eq("ChebyshevT[2, x]", "-1 + 2 x^2", 0);
    /* Headline example from the Wolfram reference. */
    assert_eval_eq("ChebyshevT[10, x]",
                   "-1 + 50 x^2 - 400 x^4 + 1120 x^6 - 1280 x^8 + 512 x^10", 0);
    assert_eval_eq("ChebyshevT[7, x]",
                   "-7 x + 56 x^3 - 112 x^5 + 64 x^7", 0);
}

void test_chebyshev_t_negative_order() {
    /* T is even in its order index: T_{-n} = T_n (no sign, no zero case). */
    assert_eval_eq("ChebyshevT[-1, x]", "x", 0);             /* T_1 */
    assert_eval_eq("ChebyshevT[-2, x]", "-1 + 2 x^2", 0);    /* T_2 */
    assert_eval_eq("ChebyshevT[-3, x]", "-3 x + 4 x^3", 0);  /* T_3 */
    assert_eval_eq("ChebyshevT[-5, x]", "5 x - 20 x^3 + 16 x^5", 0); /* T_5 */
    assert_close("ChebyshevT[-1, 2.5]", 2.5, 1e-12);         /* T_1(2.5) = 2.5 */
}

void test_chebyshev_t_exact_values() {
    assert_eval_eq("ChebyshevT[2, 3]", "17", 0);
    assert_eval_eq("ChebyshevT[2, 5]", "49", 0);
    assert_eval_eq("ChebyshevT[0, 0]", "1", 0);
    assert_eval_eq("ChebyshevT[2, 1/3]", "-7/9", 0);
    assert_eval_eq("ChebyshevT[3, 1/3]", "-23/27", 0);
    /* T_0 is constant for any argument. */
    assert_eval_eq("ChebyshevT[0, Infinity]", "1", 0);
}

void test_chebyshev_t_special_values() {
    /* T_n(0) = Cos[n Pi/2] for symbolic order. */
    assert_eval_eq("ChebyshevT[n, 0]", "Cos[1/2 n Pi]", 0);
    /* T_n(1) = 1 for symbolic order. */
    assert_eval_eq("ChebyshevT[n, 1]", "1", 0);
    /* Table at Pi/2 (integer order, symbolic argument). */
    assert_eval_eq("Table[ChebyshevT[n, Pi/2], {n, 0, 3}]",
                   "{1, 1/2 Pi, -1 + 1/2 Pi^2, -3/2 Pi + 1/2 Pi^3}", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_chebyshev_t_listable() {
    assert_eval_eq("ChebyshevT[2, {1, 2, 3}]", "{1, 7, 17}", 0);
    assert_eval_eq("ChebyshevT[{0, 1, 2, 3}, x]",
                   "{1, x, -1 + 2 x^2, -3 x + 4 x^3}", 0);
    assert_eval_eq("ChebyshevT[2, {}]", "{}", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_chebyshev_t_machine() {
    /* Integer order, inexact argument. */
    assert_close("ChebyshevT[12, 2.5]", 73180801.0, 1.0);
    /* Cancellation case: the stable recurrence gives ~0.980336, NOT 0. The
     * explicit degree-60 polynomial substituted at an inexact x loses all
     * significance here. */
    assert_close("ChebyshevT[60, 0.6666666666666666]", 0.98033630042969088, 1e-4);
    /* T_n(0) at inexact order evaluates through Cos[n Pi/2]. */
    assert_close("ChebyshevT[0.3, 0]", 0.89100652418836786, 1e-12);
    /* An inexact integer-valued order yields a clean real, not a complex with a
     * spurious tiny imaginary part. */
    assert_close("ChebyshevT[2., 3.]", 17.0, 1e-9);
    assert_close("ChebyshevT[12., 2.5]", 73180801.0, 1.0);
}

void test_chebyshev_t_complex() {
    /* Non-integer complex order at x = 3 (outside the 2F1 disk; the ArcCos form
     * answers): ChebyshevT[2-4I, 3] = 12.23044... - 11.78701... I. */
    assert_close("Re[ChebyshevT[2 - 4 I, 3.0]]", 12.230447270331584, 1e-6);
    assert_close("Im[ChebyshevT[2 - 4 I, 3.0]]", -11.787016192022614, 1e-6);
    /* Branch cut (-Infinity, -1): the ArcCos form takes Mathematica's branch,
     * giving an essentially imaginary value (Re ~ 0, Im ~ -4.85263). */
    assert_close("Re[ChebyshevT[-1.5, -2.4]]", 0.0, 1e-10);
    assert_close("Im[ChebyshevT[-1.5, -2.4]]", -4.8526281538976382, 1e-6);
}

/* ---- arbitrary-precision numerics ----------------------------------- */

void test_chebyshev_t_highprec() {
    /* Non-integer order, exact rational arguments, via N. */
    assert_eval_startswith(
        "N[ChebyshevT[1/7, 1/3], 50]",
        "0.984577955912817889705914449258059859431113237102");
    /* Complex non-integer order at MPFR precision, outside the 2F1 disk
     * (the ArcCos form answers with the Mathematica branch). */
    assert_eval_startswith("ChebyshevT[2 - 4 I, 3`40]",
                           "12.230447270331584157732321598448518916");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_chebyshev_t_symbolic() {
    assert_eval_eq("ChebyshevT[n, x]", "ChebyshevT[n, x]", 0);
    /* Exact non-integer order with exact/symbolic argument stays symbolic
     * (half-integer closed forms are deferred). */
    assert_eval_eq("ChebyshevT[1/2, x]", "ChebyshevT[1/2, x]", 0);
    assert_eval_eq("ChebyshevT[3/2, 2]", "ChebyshevT[3/2, 2]", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_chebyshev_t_derivative() {
    /* d/dx T_n(x) = n U_{n-1}(x). */
    assert_eval_eq("D[ChebyshevT[n, x], x]", "n ChebyshevU[-1 + n, x]", 0);
    /* The derivative of an explicit polynomial is the polynomial's derivative:
     * d/dx (2 x^2 - 1) = 4 x. */
    assert_eval_eq("D[ChebyshevT[2, x], x]", "4 x", 0);
    /* d/dx T_5(x) = 5 U_4(x) = 5 - 60 x^2 + 80 x^4. */
    assert_eval_eq("D[ChebyshevT[5, x], x]", "5 - 60 x^2 + 80 x^4", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_chebyshev_t_argcount() {
    /* Wrong arg counts stay unevaluated (an argrx/argx message goes to stderr). */
    assert_eval_eq("ChebyshevT[]", "ChebyshevT[]", 0);
    assert_eval_eq("ChebyshevT[3]", "ChebyshevT[3]", 0);
    assert_eval_eq("ChebyshevT[1, 2, 3]", "ChebyshevT[1, 2, 3]", 0);
}

void test_chebyshev_t_attributes() {
    SymbolDef* d = symtab_get_def("ChebyshevT");
    ASSERT_MSG(d != NULL, "ChebyshevT not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "ChebyshevT not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "ChebyshevT not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "ChebyshevT not Protected");
}

/* High degree at a complex arbitrary-precision argument (the closed form
 * (w^n + w^-n)/2 by binary powering): against Cos[n ArcCos[z]] at 60 digits,
 * off the real axis, next to it inside and outside [-1, 1], on the imaginary
 * axis and near the endpoints +-1, odd and even degree. */
void test_chebyshev_t_complex_high_degree() {
    const char* zs[] = { "3/5 Exp[I/3]", "3/10 + I/10^10", "2 + 3 I", "-1/2 + I/100",
                         "-7/10 - I/3", "I", "1 + I/10^6", "-1 - I/10^8" };
    const char* ns[] = { "65", "2000", "2001" };
    char buf[512];
    for (size_t i = 0; i < sizeof(zs) / sizeof(zs[0]); i++)
        for (size_t k = 0; k < sizeof(ns) / sizeof(ns[0]); k++) {
            snprintf(buf, sizeof buf,
                     "With[{v = ChebyshevT[%s, N[%s, 30]], r = N[Cos[%s ArcCos[%s]], 60]}, "
                     "Abs[v - r] <= 10^-28 Abs[r] && Precision[v] >= 29]",
                     ns[k], zs[i], ns[k], zs[i]);
            assert_eval_eq(buf, "True", 0);
        }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_chebyshev_t_complex_high_degree);

    TEST(test_chebyshev_t_polynomials);
    TEST(test_chebyshev_t_negative_order);
    TEST(test_chebyshev_t_exact_values);
    TEST(test_chebyshev_t_special_values);
    TEST(test_chebyshev_t_listable);
    TEST(test_chebyshev_t_machine);
    TEST(test_chebyshev_t_complex);
    TEST(test_chebyshev_t_highprec);
    TEST(test_chebyshev_t_symbolic);
    TEST(test_chebyshev_t_derivative);
    TEST(test_chebyshev_t_argcount);
    TEST(test_chebyshev_t_attributes);

    printf("All ChebyshevT tests passed.\n");
    return 0;
}
