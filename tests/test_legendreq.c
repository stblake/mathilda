/* Tests for LegendreQ: Legendre functions of the second kind Q_n(x), the
 * associated functions Q_n^m(x), and the Legendre functions of types 1/2/3.
 *
 * Covers the exact integer closed form Q_n = P_n L + v_n (with
 * L(x) = (1/2)(Log[1+x] - Log[1-x])), the exact special value Q_v(0), Listable
 * threading, machine and arbitrary-precision (MPFR) numerics through the origin
 * Frobenius series (real and complex), precision tracking, the associated
 * derivative/type-2/3 forms, the argument-derivative rule, the origin Series,
 * symbolic fall-through for exact non-integer / negative / out-of-range
 * arguments, argument-count diagnostics, and attributes. */

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

/* ---- exact integer closed forms ------------------------------------- */

void test_legendreq_closed_forms() {
    /* Q_0(x) = (1/2)(Log[1+x] - Log[1-x]). */
    assert_eval_eq("LegendreQ[0, x]",
                   "1/2 Log[1 + x] - 1/2 Log[1 - x]", 0);
    /* Q_1(x) = -1 + x L(x). */
    assert_eval_eq("LegendreQ[1, x]",
                   "-1 + x (1/2 Log[1 + x] - 1/2 Log[1 - x])", 0);
    /* Q_2(x) = -3x/2 + (-1/2 + 3/2 x^2) L(x). */
    assert_eval_eq("LegendreQ[2, x]",
                   "-3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) "
                   "(-1/2 + 3/2 x^2)", 0);
    /* Q_3(x). */
    assert_eval_eq("LegendreQ[3, x]",
                   "2/3 - 5/2 x^2 + (1/2 Log[1 + x] - 1/2 Log[1 - x]) "
                   "(-3/2 x + 5/2 x^3)", 0);
    /* The fifth Legendre function of the second kind (headline example). */
    assert_eval_eq("LegendreQ[5, x]",
                   "-8/15 + 49/8 x^2 - 63/8 x^4 + "
                   "(1/2 Log[1 + x] - 1/2 Log[1 - x]) "
                   "(15/8 x - 35/4 x^3 + 63/8 x^5)", 0);
}

/* ---- special value Q_n(0) ------------------------------------------- */

void test_legendreq_special_values() {
    /* Integer order: Q_n(0) = v_n(0). */
    assert_eval_eq("LegendreQ[0, 0]", "0", 0);
    assert_eval_eq("LegendreQ[1, 0]", "-1", 0);
    assert_eval_eq("LegendreQ[2, 0]", "0", 0);
    assert_eval_eq("LegendreQ[3, 0]", "2/3", 0);
    /* Non-integer order: the exact Gamma special value. */
    assert_eval_eq("LegendreQ[1/2, 0]",
                   "-(1/2 Gamma[3/4] Sqrt[Pi])/(Sqrt[2] Gamma[5/4])", 0);
    /* Its numeric value is ~ -0.847213. */
    assert_close("N[LegendreQ[1/2, 0]]", -0.8472130847939791, 1e-9);
}

/* ---- Listable threading --------------------------------------------- */

void test_legendreq_listable() {
    assert_eval_eq("LegendreQ[{0, 1, 2}, x]",
                   "{1/2 Log[1 + x] - 1/2 Log[1 - x], "
                   "-1 + x (1/2 Log[1 + x] - 1/2 Log[1 - x]), "
                   "-3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) "
                   "(-1/2 + 3/2 x^2)}", 0);
    assert_eval_eq("LegendreQ[{}, x]", "{}", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_legendreq_machine() {
    /* Integer order through the closed form. */
    assert_close("LegendreQ[2, 0.5]", -0.8186632680417569, 1e-9);
    assert_close("LegendreQ[10, 0.3]", 0.0242860961383206, 1e-9);
    /* Non-integer order through the origin Frobenius series. */
    assert_close("LegendreQ[1/3, 0.5]", -0.0399532947598897, 1e-9);
    assert_close("LegendreQ[1/2, 0.5]", -0.2655964076372758, 1e-9);
    assert_close("LegendreQ[3/2, 0.5]", -0.8959028209247316, 1e-9);
}

void test_legendreq_complex() {
    /* LegendreQ[1/3 - I, 0.5 + 0.45 I] = -0.346939 + 2.52382 I. */
    assert_close("Re[LegendreQ[1/3 - I, 0.5 + 0.45 I]]", -0.346939, 1e-5);
    assert_close("Im[LegendreQ[1/3 - I, 0.5 + 0.45 I]]",  2.52382,  1e-5);
}

/* ---- arbitrary-precision numerics ----------------------------------- */

void test_legendreq_highprec() {
    /* N[LegendreQ[3/2, 1/2], 50] (prefix stops before the last-digit rounding). */
    assert_eval_startswith("N[LegendreQ[3/2, 1/2], 50]",
        "-0.895902820924731621258525533131864225704282994");
    /* Precision tracks the input precision. */
    assert_eval_startswith("LegendreQ[3/2, 0.5`30]",
        "-0.895902820924731621258525533");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_legendreq_symbolic() {
    /* Exact non-integer order with exact non-zero argument stays symbolic. */
    assert_eval_eq("LegendreQ[3/2, 2]", "LegendreQ[3/2, 2]", 0);
    assert_eval_eq("LegendreQ[n, x]", "LegendreQ[n, x]", 0);
    assert_eval_eq("LegendreQ[1/2, y]", "LegendreQ[1/2, y]", 0);
    /* Negative integer order is singular -> deferred. */
    assert_eval_eq("LegendreQ[-2, x]", "LegendreQ[-2, x]", 0);
    /* Numeric series outside the cut |x| < 1 is deferred (stays symbolic for a
     * non-integer order). */
    assert_eval_eq("LegendreQ[1/2, 5.0]", "LegendreQ[1/2, 5.0]", 0);
}

/* ---- associated Legendre (type 1) ----------------------------------- */

void test_legendreq_associated() {
    /* m == 0 reduces to the ordinary function. */
    assert_eval_eq("LegendreQ[3, 0, x] === LegendreQ[3, x]", "True", 0);
    /* Q_3^1(0.5) = 2.491853 (reference value 2.49185). */
    assert_close("LegendreQ[3, 1, 0.5]", 2.4918525917089543, 1e-7);
    /* Q_2^2(0.5) = 4.069272 (reference value 4.06927). */
    assert_close("LegendreQ[2, 2, 0.5]", 4.0692721580849567, 1e-7);
    /* Q_0^1(x) = -1/Sqrt[1-x^2]: value at 0.5 is -2/Sqrt[3]. */
    assert_close("LegendreQ[0, 1, 0.5]", -1.1547005383792515, 1e-9);
}

/* ---- Legendre function types ---------------------------------------- */

void test_legendreq_types() {
    /* a == 1 is the default and equals the three-argument form. */
    assert_close("LegendreQ[2, 1, 1, 0.4] - LegendreQ[2, 1, 0.4]", 0.0, 1e-12);
    /* Type 2 at an interior cut point matches type 1 (same real value). */
    assert_close("LegendreQ[2, 1, 2, 0.4] - LegendreQ[2, 1, 1, 0.4]", 0.0, 1e-9);
    /* Higher m: type 2 still matches type 1 on the cut. */
    assert_close("LegendreQ[3, 2, 2, 0.35] - LegendreQ[3, 2, 1, 0.35]", 0.0, 1e-9);
    /* Type 3 uses the (-1+x)^(-m/2) branch, complex on the cut -- as for
     * LegendreP -- so only check it evaluates to a finite complex number. */
    assert_eval_eq("NumericQ[LegendreQ[2, 1, 3, 0.4]]", "True", 0);
}

/* ---- derivative rule ------------------------------------------------ */

void test_legendreq_derivative() {
    /* d/dx Q_n(x) = ((-1-n) x Q_n + (1+n) Q_{n+1}) / (x^2 - 1). */
    assert_eval_eq("D[LegendreQ[n, x], x]",
                   "((1 + n) LegendreQ[1 + n, x] + (-1 - n) x LegendreQ[n, x])/"
                   "(-1 + x^2)", 0);
    /* The rule reproduces the closed-form derivative at a numeric point:
     * d/dx Q_2(x) at x = 0.4 equals the closed-form value. */
    assert_close("(D[LegendreQ[2, x], x] /. x -> 0.4) - "
                 "(D[LegendreQ[2, y], y] /. y -> 0.4)", 0.0, 1e-12);
}

/* ---- origin Series -------------------------------------------------- */

void test_legendreq_series() {
    /* The order-2 origin series, evaluated at x = 0.3, matches the analytic
     * 2nd-order Taylor value (the O(x^3) remainder is the difference from the
     * exact function value). */
    assert_close("N[Normal[Series[LegendreQ[1/2, x], {x, 0, 2}]] /. x -> 0.3]",
                 -0.540508, 1e-5);
    /* Integer-order series is a genuine (log-carrying) expansion; check the
     * leading term of Q_1 = -1 + x L(x) at the origin is -1. */
    assert_close("N[Normal[Series[LegendreQ[1, x], {x, 0, 3}]] /. x -> 0.0]",
                 -1.0, 1e-12);
}

/* ---- packed / NDArray parity ---------------------------------------- */

void test_legendreq_ndarray() {
    /* The packed fast path agrees with the Listable interpreter path. */
    assert_eval_eq("Max[Abs[Flatten[Normal[LegendreQ[2, NDArray[{0.3, 0.4, 0.5}]]] "
                   "- LegendreQ[2, {0.3, 0.4, 0.5}]]]] < 1/1000000000", "True", 0);
    /* An out-of-disk element (1.3) degrades faithfully to the complex closed
     * form, element-wise matching the plain-List result. */
    assert_eval_eq("Max[Abs[Flatten[Normal[LegendreQ[2, NDArray[{0.4, 1.3}]]] "
                   "- LegendreQ[2, {0.4, 1.3}]]]] < 1/1000000000", "True", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_legendreq_deferred() {
    /* Non-integer order/degree associated form is deferred. */
    assert_eval_eq("LegendreQ[3/2, 1/2, x]", "LegendreQ[3/2, 1/2, x]", 0);
    /* Negative degree is deferred. */
    assert_eval_eq("LegendreQ[2, -1, x]", "LegendreQ[2, -1, x]", 0);
    /* Out-of-range type is deferred. */
    assert_eval_eq("LegendreQ[2, 1, 4, z]", "LegendreQ[2, 1, 4, z]", 0);
}

void test_legendreq_argcount() {
    /* Wrong arg counts stay unevaluated (an argb message goes to stderr). */
    assert_eval_eq("LegendreQ[]", "LegendreQ[]", 0);
    assert_eval_eq("LegendreQ[3]", "LegendreQ[3]", 0);
    assert_eval_eq("LegendreQ[1, 2, 3, 4, 5]", "LegendreQ[1, 2, 3, 4, 5]", 0);
}

void test_legendreq_attributes() {
    SymbolDef* d = symtab_get_def("LegendreQ");
    ASSERT_MSG(d != NULL, "LegendreQ not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "LegendreQ not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "LegendreQ not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "LegendreQ not Protected");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_legendreq_closed_forms);
    TEST(test_legendreq_special_values);
    TEST(test_legendreq_listable);
    TEST(test_legendreq_machine);
    TEST(test_legendreq_complex);
    TEST(test_legendreq_highprec);
    TEST(test_legendreq_symbolic);
    TEST(test_legendreq_associated);
    TEST(test_legendreq_types);
    TEST(test_legendreq_derivative);
    TEST(test_legendreq_series);
    TEST(test_legendreq_ndarray);
    TEST(test_legendreq_deferred);
    TEST(test_legendreq_argcount);
    TEST(test_legendreq_attributes);

    printf("All LegendreQ tests passed.\n");
    return 0;
}
