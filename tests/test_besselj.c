/* Tests for BesselJ[n, z], the Bessel function of the first kind J_n(z).
 *
 * Covers: exact values at the origin and symbolic/exact passthrough; machine-
 * and arbitrary-precision real numerics across the power-series (small |z|) and
 * asymptotic (large |z|) regimes; an independent mpfr_jn oracle for integer
 * order; a closed-form oracle for half-integer order (J_{1/2}, J_{-1/2}); a
 * recurrence-consistency oracle for generic non-integer order; complex order
 * and argument; precision tracking; list threading; the derivative rule;
 * Series at 0 and Infinity; SeriesCoefficient; stress (deep precision, large
 * order/argument, near a zero of J0); argument-count diagnostics; attributes.
 *
 * NB: the symbolic half-integer -> elementary rewrites and negative-integer
 * reflection live in src/internal/bessel.m, which is loaded by the REPL but NOT
 * by this test harness (no init.m). Those are verified separately in the REPL;
 * here the corresponding calls stay symbolic / use the numeric C path. */

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

#ifdef USE_MPFR
#include <mpfr.h>
#endif

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

void test_besselj_exact() {
    /* J_0(0) = 1, J_n(0) = 0 for integer n != 0 (either sign). */
    assert_eval_eq("BesselJ[0, 0]", "1", 0);
    assert_eval_eq("BesselJ[3, 0]", "0", 0);
    assert_eval_eq("BesselJ[5, 0]", "0", 0);
    assert_eval_eq("BesselJ[-2, 0]", "0", 0);
    /* Non-integer order at the origin: J_nu(0) = 0 for Re(nu) > 0, and
     * ComplexInfinity for non-integer Re(nu) < 0 ((z/2)^nu -> inf). */
    assert_eval_eq("BesselJ[1/2, 0]", "0", 0);
    assert_eval_eq("BesselJ[3/2, 0]", "0", 0);
    assert_eval_eq("BesselJ[-1/2, 0]", "ComplexInfinity", 0);
    assert_eval_eq("BesselJ[-3/2, 0]", "ComplexInfinity", 0);
    /* Exact non-zero / symbolic arguments stay symbolic (no auto-numericize). */
    assert_eval_eq("BesselJ[0, 2]", "BesselJ[0, 2]", 0);
    assert_eval_eq("BesselJ[n, x]", "BesselJ[n, x]", 0);
    /* Symbolic order at the origin stays symbolic (could be 0, 1, or inf). */
    assert_eval_eq("BesselJ[a, 0]", "BesselJ[a, 0]", 0);
    /* Half-integer order with exact numeric argument stays unevaluated. */
    assert_eval_eq("BesselJ[11/2, 1]", "BesselJ[11/2, 1]", 0);
}

/* ---- argument parity: J_n(-z) = (-1)^n J_n(z), integer n only ------- */

void test_besselj_argument_parity() {
    /* Even / odd integer order folds the negated argument. */
    assert_eval_eq("BesselJ[0, -z]", "BesselJ[0, z]", 0);
    assert_eval_eq("BesselJ[2, -z]", "BesselJ[2, z]", 0);
    assert_eval_eq("BesselJ[1, -z]", "-BesselJ[1, z]", 0);
    assert_eval_eq("BesselJ[3, -z]", "-BesselJ[3, z]", 0);
    /* Negated coefficient inside a product. */
    assert_eval_eq("BesselJ[1, -2 x]", "-BesselJ[1, 2 x]", 0);
    /* Exact negative integer argument folds (stays symbolic). */
    assert_eval_eq("BesselJ[1, -2]", "-BesselJ[1, 2]", 0);
    /* Non-integer / symbolic order does NOT fold (parity needs integer n). */
    assert_eval_eq("BesselJ[n, -z]", "BesselJ[n, -z]", 0);
    /* Positive-looking argument is untouched (no spurious fold / loop). */
    assert_eval_eq("BesselJ[1, z]", "BesselJ[1, z]", 0);
}

/* ---- machine-precision real ----------------------------------------- */

void test_besselj_machine_real() {
    assert_close("BesselJ[0, 5.2]", -0.11029043979098654, 1e-7);
    assert_close("BesselJ[0, 4.0]", -0.3971498098638474, 1e-7);
    assert_close("BesselJ[1, 1.0]",  0.4400505857449335, 1e-7);
    assert_close("BesselJ[2, 3.0]",  0.4860912605858911, 1e-7);
    assert_close("BesselJ[0, -5.2]", -0.11029043979098654, 1e-7); /* J0 even */
    assert_close("BesselJ[1, -1.0]", -0.4400505857449335, 1e-7);  /* J1 odd  */
    /* Half-integer order, inexact argument -> general numeric path (NOT the
     * catastrophically-cancelling elementary form). */
    assert_close("BesselJ[35/2, 1.]", 3.551525807680e-21, 1e-30);
    /* Non-integer (real) order. */
    assert_close("BesselJ[0.5, 2.0]", 0.5130161365618272, 1e-7);
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_besselj_arbitrary_real() {
    /* N[BesselJ[0,4],50] matches the reference to all 50 digits. */
    assert_eval_startswith("N[BesselJ[0, 4], 50]",
        "-0.39714980986384737228659076845169804197561868528939");
    /* Precision tracks the input. */
    assert_eval_startswith("BesselJ[0, 4.000000000000000000000000]",
        "-0.3971498098638473722865");
    assert_eval_startswith("BesselJ[0, 5.2`30]",
        "-0.11029043979098653962103299647");
}

/* ---- mpfr_jn oracle for integer order ------------------------------- */
/* Independently compute J_n(x) with MPFR's mpfr_jn over a grid spanning the
 * power-series (|x| small) and asymptotic (|x| large) regimes, for several
 * integer orders and both signs of x. This validates the full pipeline
 * (parse -> evaluate -> result construction -> precision) end to end. */
void test_besselj_mpfr_oracle() {
#ifdef USE_MPFR
    mpfr_t ax, ref;
    mpfr_init2(ax, 160);
    mpfr_init2(ref, 160);
    for (long order = 0; order <= 4; order++) {
        for (int t = -40; t <= 40; t++) {
            double x = t * 0.6;
            char in[64];
            snprintf(in, sizeof(in), "BesselJ[%ld, %.4f]", order, x);
            double got = eval_real(in);
            mpfr_set_d(ax, x, MPFR_RNDN);
            mpfr_jn(ref, order, ax, MPFR_RNDN);
            double want = mpfr_get_d(ref, MPFR_RNDN);
            ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                       "BesselJ[%ld, %.4f]: builtin %.14g vs mpfr_jn %.14g",
                       order, x, got, want);
        }
    }
    mpfr_clear(ax);
    mpfr_clear(ref);
#endif
}

/* ---- closed-form oracle for half-integer order ---------------------- */
/* J_{1/2}(x) = sqrt(2/(pi x)) sin x,  J_{-1/2}(x) = sqrt(2/(pi x)) cos x.
 * Validates the general (non-integer order) core against an independent
 * formula that mpfr_jn cannot provide. */
void test_besselj_halfinteger_oracle() {
    for (int t = 1; t <= 30; t++) {
        double x = t * 0.5;
        char in[64];
        snprintf(in, sizeof(in), "BesselJ[1/2, %.4f]", x);
        double got = eval_real(in);
        double want = sqrt(2.0 / (M_PI * x)) * sin(x);
        ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                   "BesselJ[1/2, %.4f]: builtin %.14g vs closed form %.14g",
                   x, got, want);

        snprintf(in, sizeof(in), "BesselJ[-1/2, %.4f]", x);
        got = eval_real(in);
        want = sqrt(2.0 / (M_PI * x)) * cos(x);
        ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                   "BesselJ[-1/2, %.4f]: builtin %.14g vs closed form %.14g",
                   x, got, want);
    }
}

/* ---- recurrence-consistency oracle for generic order ---------------- */
/* The true Bessel function satisfies J_{v-1}(x) + J_{v+1}(x) = (2v/x) J_v(x).
 * Checking this for a generic non-integer order (v = 1/3) tests the general
 * core's internal consistency independent of any tabulated value. */
void test_besselj_recurrence() {
    for (int t = 1; t <= 20; t++) {
        double x = t * 0.7;
        char lhs[96], rhs[96];
        snprintf(lhs, sizeof(lhs), "BesselJ[-2/3, %.4f] + BesselJ[4/3, %.4f]", x, x);
        snprintf(rhs, sizeof(rhs), "(2/3) BesselJ[1/3, %.4f] / %.4f", x, x);
        double l = eval_real(lhs);
        double r = eval_real(rhs);
        ASSERT_MSG(fabs(l - r) <= 1e-9 * (1.0 + fabs(r)),
                   "recurrence v=1/3 at x=%.4f: %.14g vs %.14g", x, l, r);
    }
}

/* ---- complex order and argument ------------------------------------- */

void test_besselj_complex() {
    /* Spec value: BesselJ[7/3 + I, 4.5 - I] = 1.18908 + 0.715653 I. */
    assert_close("Re[BesselJ[7/3 + I, 4.5 - I]]", 1.18908, 1e-4);
    assert_close("Im[BesselJ[7/3 + I, 4.5 - I]]", 0.715653, 1e-4);
    /* Integer order, complex argument (general core). */
    assert_close("Re[BesselJ[2, 4.5 - 1.5 I]]", 0.353675, 1e-4);
    assert_close("Im[BesselJ[2, 4.5 - 1.5 I]]", 0.647659, 1e-4);
    /* Complex consistency: J0 of a real argument given as Complex with 0 im
     * matches the real evaluation. */
    assert_close("Re[BesselJ[0, 3.0 + 0.0 I]]", eval_real("BesselJ[0, 3.0]"), 1e-12);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_besselj_listable() {
    assert_eval_eq("BesselJ[1, {a, b}]", "{BesselJ[1, a], BesselJ[1, b]}", 0);
    assert_eval_eq("BesselJ[1, {}]", "{}", 0);
    assert_close("BesselJ[1, {0.5, 1.0, 1.5}][[1]]", 0.2422684577, 1e-7);
    assert_close("BesselJ[1, {0.5, 1.0, 1.5}][[2]]", 0.4400505857, 1e-7);
    assert_close("BesselJ[1, {0.5, 1.0, 1.5}][[3]]", 0.5579365079, 1e-7);
    /* Threading also over the order. */
    assert_eval_eq("BesselJ[{0, 1}, x]", "{BesselJ[0, x], BesselJ[1, x]}", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_besselj_derivative() {
    assert_eval_eq("D[BesselJ[n, x], x]",
                   "1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])", 0);
    /* Chain rule through the argument. */
    assert_eval_eq("D[BesselJ[n, x^2], x]",
                   "x (BesselJ[-1 + n, x^2] - BesselJ[1 + n, x^2])", 0);
    /* No x-dependence -> 0. */
    assert_eval_eq("D[BesselJ[n, y], x]", "0", 0);
}

/* ---- Series at 0 ---------------------------------------------------- */

void test_besselj_series_zero() {
    assert_eval_eq("Series[BesselJ[0, x], {x, 0, 10}]",
        "1 - 1/4 x^2 + 1/64 x^4 - 1/2304 x^6 + 1/147456 x^8 - "
        "1/14745600 x^10 + O[x]^11", 0);
    assert_eval_eq("Series[BesselJ[2, x], {x, 0, 6}]",
        "1/8 x^2 - 1/96 x^4 + 1/3072 x^6 + O[x]^7", 0);
}

/* ---- Series at Infinity --------------------------------------------- */

void test_besselj_series_infinity() {
    assert_eval_eq("Series[BesselJ[0, x], {x, Infinity, 1}]",
        "Cos[1/4 Pi - x] (Sqrt[2/Pi] Sqrt[1/x] + O[1/x]^(5/2))", 0);
    assert_eval_eq("Series[BesselJ[0, x], {x, Infinity, 2}]",
        "Cos[1/4 Pi - x] (Sqrt[2/Pi] Sqrt[1/x] + O[1/x]^(5/2)) + "
        "(-1/8 Sqrt[2/Pi] (1/x)^(3/2) + O[1/x]^(7/2)) Sin[1/4 Pi - x]", 0);
}

/* ---- SeriesCoefficient ---------------------------------------------- */

void test_besselj_seriescoefficient() {
    assert_eval_eq("SeriesCoefficient[BesselJ[0, x], {x, 0, 0}]", "1", 0);
    assert_eval_eq("SeriesCoefficient[BesselJ[0, x], {x, 0, 3}]", "0", 0);
    assert_eval_eq("SeriesCoefficient[BesselJ[0, x], {x, 0, 4}]", "1/64", 0);
    assert_eval_eq("SeriesCoefficient[BesselJ[0, x], {x, 0, 6}]", "-1/2304", 0);
    assert_eval_eq("SeriesCoefficient[BesselJ[2, x], {x, 0, 2}]", "1/8", 0);
}

/* ---- stress --------------------------------------------------------- */

void test_besselj_stress() {
#ifdef USE_MPFR
    /* Deep precision: 200-bit J_0(4) must match mpfr_jn to ~200 bits. */
    {
        Expr* e = parse_expression("BesselJ[0, 4.0`70]");
        Expr* r = evaluate(e);
        expr_free(e);
        ASSERT_MSG(r->type == EXPR_MPFR, "BesselJ[0,4.0`70] not MPFR");
        mpfr_t ref, x4;
        mpfr_init2(ref, mpfr_get_prec(r->data.mpfr));
        mpfr_init2(x4, mpfr_get_prec(r->data.mpfr));
        mpfr_set_ui(x4, 4, MPFR_RNDN);
        mpfr_jn(ref, 0, x4, MPFR_RNDN);
        mpfr_sub(ref, ref, r->data.mpfr, MPFR_RNDN);
        mpfr_abs(ref, ref, MPFR_RNDN);
        double err = mpfr_get_d(ref, MPFR_RNDN);
        ASSERT_MSG(err < 1e-60, "BesselJ[0,4.0`70] error %.3e too large", err);
        mpfr_clear(ref); mpfr_clear(x4);
        expr_free(r);
    }
#endif
    /* Large order with small argument: tiny but finite and accurate. */
    assert_close("BesselJ[20, 1.0]", 3.873503008525e-25, 1e-33);
    /* Large argument (deep asymptotic regime). */
    assert_close("BesselJ[0, 100.0]", 0.0199858503042, 1e-9);
    assert_close("BesselJ[3, 80.0]", 0.0594743333305, 1e-9);
    /* Near the first zero of J0 (2.404825557695773): value ~ 0. */
    assert_close("BesselJ[0, 2.404825557695773]", 0.0, 1e-12);
}

/* ---- argument count ------------------------------------------------- */

void test_besselj_argcount() {
    assert_eval_eq("BesselJ[]", "BesselJ[]", 0);
    assert_eval_eq("BesselJ[1, 2, 3]", "BesselJ[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_besselj_attributes() {
    SymbolDef* d = symtab_get_def("BesselJ");
    ASSERT_MSG(d != NULL, "BesselJ not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "BesselJ not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "BesselJ not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "BesselJ not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_besselj_exact);
    TEST(test_besselj_argument_parity);
    TEST(test_besselj_machine_real);
    TEST(test_besselj_arbitrary_real);
    TEST(test_besselj_mpfr_oracle);
    TEST(test_besselj_halfinteger_oracle);
    TEST(test_besselj_recurrence);
    TEST(test_besselj_complex);
    TEST(test_besselj_listable);
    TEST(test_besselj_derivative);
    TEST(test_besselj_series_zero);
    TEST(test_besselj_series_infinity);
    TEST(test_besselj_seriescoefficient);
    TEST(test_besselj_stress);
    TEST(test_besselj_argcount);
    TEST(test_besselj_attributes);

    printf("All BesselJ tests passed.\n");
    return 0;
}
