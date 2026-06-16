/* Tests for BesselY[n, z], the Bessel function of the second kind Y_n(z).
 *
 * Covers: exact singular values at the origin (-Infinity / ComplexInfinity) and
 * symbolic/exact passthrough; machine- and arbitrary-precision real numerics
 * across the logarithmic-series (small |z|, integer order), connection-formula
 * (small |z|, non-integer order) and asymptotic (large |z|) regimes; an
 * independent mpfr_yn oracle for integer order on z > 0; a closed-form oracle
 * for half-integer order (Y_{1/2}, Y_{-1/2}); a recurrence-consistency oracle
 * for generic non-integer order; the genuinely-complex value on the negative
 * real branch cut; complex order and argument; precision tracking; list
 * threading; the derivative rule; Series at 0 (logarithmic) and at Infinity;
 * stress (deep precision, large argument); argument-count diagnostics;
 * attributes.
 *
 * NB: the symbolic half-integer -> elementary rewrites and the negative-integer
 * reflection live in src/internal/bessel.m, which is loaded by the REPL but NOT
 * by this harness (no init.m). Half-integer calls here therefore use the
 * inexact numeric C path; symbolic reflection is verified separately in the
 * REPL. */

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

void test_bessely_exact() {
    /* Y_0(0) = -Infinity, Y_n(0) = ComplexInfinity for integer n != 0. */
    assert_eval_eq("BesselY[0, 0]", "-Infinity", 0);
    assert_eval_eq("BesselY[3, 0]", "ComplexInfinity", 0);
    assert_eval_eq("BesselY[-2, 0]", "ComplexInfinity", 0);
    /* Exact non-zero / symbolic arguments stay symbolic (no auto-numericize). */
    assert_eval_eq("BesselY[0, 2]", "BesselY[0, 2]", 0);
    assert_eval_eq("BesselY[n, x]", "BesselY[n, x]", 0);
    /* Half-integer order with exact numeric argument stays unevaluated. */
    assert_eval_eq("BesselY[11/2, 1]", "BesselY[11/2, 1]", 0);
    /* Non-integer order at the origin stays symbolic. */
    assert_eval_eq("BesselY[1/2, 0]", "BesselY[1/2, 0]", 0);
}

/* ---- no argument-parity fold (branch cut on the negative axis) ------- */

void test_bessely_no_parity() {
    /* Unlike BesselJ, Y_n(-z) does NOT fold (branch cut), so symbolic negated
     * arguments stay put for every order. */
    assert_eval_eq("BesselY[0, -z]", "BesselY[0, -z]", 0);
    assert_eval_eq("BesselY[2, -z]", "BesselY[2, -z]", 0);
    assert_eval_eq("BesselY[n, -z]", "BesselY[n, -z]", 0);
}

/* ---- machine-precision real ----------------------------------------- */

void test_bessely_machine_real() {
    assert_close("BesselY[0, 1.0]",  0.0882569642156769583, 1e-7);
    assert_close("BesselY[1, 1.0]", -0.781212821300288719,  1e-7);
    assert_close("BesselY[0, 2.5]",  0.4980703596152319,    1e-7);
    assert_close("BesselY[1, 2.0]", -0.107032431540937547,  1e-7);
    assert_close("BesselY[2, 3.0]", -0.16040039348492373,   1e-7);
    assert_close("BesselY[0, 5.2]", -0.331250934819884435,  1e-7);
    /* Half-integer order, inexact argument -> general numeric path. */
    assert_close("BesselY[1/2, 3.0]", 0.456048820794633179, 1e-7);
    /* Non-integer (real) order via the connection formula.
     * Y_{1/2}(2) = -sqrt(1/Pi) cos(2). */
    assert_close("BesselY[0.5, 2.0]", 0.234785710406, 1e-7);
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_bessely_arbitrary_real() {
    /* N[BesselY[0,1],50] matches the reference (last digit unpinned by N). */
    assert_eval_startswith("N[BesselY[0, 1], 50]",
        "0.0882569642156769579829267660235151628278175230");
    /* Precision tracks the input. */
    assert_eval_startswith("BesselY[0, 1.0000000000000000000000000000000000000000]",
        "0.08825696421567695798292676602351516282");
    assert_eval_startswith("BesselY[0, 2.5`20]",
        "0.4980703596152318878");
}

/* ---- mpfr_yn oracle for integer order (z > 0) ----------------------- */
/* Independently compute Y_n(x) with MPFR's mpfr_yn over a positive-x grid
 * spanning the logarithmic-series (|x| small) and asymptotic (|x| large)
 * regimes, for several integer orders. mpfr_yn requires x > 0 (Y is complex on
 * the negative real branch cut). Validates the full pipeline end to end. */
void test_bessely_mpfr_oracle() {
#ifdef USE_MPFR
    mpfr_t ax, ref;
    mpfr_init2(ax, 160);
    mpfr_init2(ref, 160);
    for (long order = 0; order <= 4; order++) {
        for (int t = 1; t <= 40; t++) {
            double x = t * 0.6;
            char in[64];
            snprintf(in, sizeof(in), "BesselY[%ld, %.4f]", order, x);
            double got = eval_real(in);
            mpfr_set_d(ax, x, MPFR_RNDN);
            mpfr_yn(ref, order, ax, MPFR_RNDN);
            double want = mpfr_get_d(ref, MPFR_RNDN);
            ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                       "BesselY[%ld, %.4f]: builtin %.14g vs mpfr_yn %.14g",
                       order, x, got, want);
        }
    }
    mpfr_clear(ax);
    mpfr_clear(ref);
#endif
}

/* ---- closed-form oracle for half-integer order ---------------------- */
/* Y_{1/2}(x) = -sqrt(2/(pi x)) cos x,  Y_{-1/2}(x) = sqrt(2/(pi x)) sin x.
 * Validates the general (non-integer order) connection-formula core against an
 * independent formula that mpfr_yn cannot provide. */
void test_bessely_halfinteger_oracle() {
    for (int t = 1; t <= 30; t++) {
        double x = t * 0.5;
        char in[64];
        snprintf(in, sizeof(in), "BesselY[1/2, %.4f]", x);
        double got = eval_real(in);
        double want = -sqrt(2.0 / (M_PI * x)) * cos(x);
        ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                   "BesselY[1/2, %.4f]: builtin %.14g vs closed form %.14g",
                   x, got, want);

        snprintf(in, sizeof(in), "BesselY[-1/2, %.4f]", x);
        got = eval_real(in);
        want = sqrt(2.0 / (M_PI * x)) * sin(x);
        ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                   "BesselY[-1/2, %.4f]: builtin %.14g vs closed form %.14g",
                   x, got, want);
    }
}

/* ---- recurrence-consistency oracle for generic order ---------------- */
/* The true Bessel function satisfies Y_{v-1}(x) + Y_{v+1}(x) = (2v/x) Y_v(x).
 * Checking this for a generic non-integer order (v = 1/3) tests the connection-
 * formula core's internal consistency independent of any tabulated value. */
void test_bessely_recurrence() {
    for (int t = 1; t <= 20; t++) {
        double x = t * 0.7;
        char lhs[96], rhs[96];
        snprintf(lhs, sizeof(lhs), "BesselY[-2/3, %.4f] + BesselY[4/3, %.4f]", x, x);
        snprintf(rhs, sizeof(rhs), "(2/3) BesselY[1/3, %.4f] / %.4f", x, x);
        double l = eval_real(lhs);
        double r = eval_real(rhs);
        ASSERT_MSG(fabs(l - r) <= 1e-9 * (1.0 + fabs(r)),
                   "recurrence v=1/3 at x=%.4f: %.14g vs %.14g", x, l, r);
    }
}

/* ---- branch cut: genuinely complex on the negative real axis -------- */

void test_bessely_branch_cut() {
    /* Y_0(-1) = Y_0(1) + 2 I J_0(1) (DLMF 10.11.5) = 0.0882570 + 1.5303954 I,
     * with 2 J_0(1) = 1.5303953731. */
    assert_close("Re[BesselY[0, -1.0]]", 0.0882569642156769583, 1e-6);
    assert_close("Im[BesselY[0, -1.0]]", 1.5303953731, 1e-6);
}

/* ---- complex order and argument ------------------------------------- */

void test_bessely_complex() {
    /* Spec value: BesselY[0.5 I, 3 - I] = 1.04686 + 0.884784 I. */
    assert_close("Re[BesselY[0.5 I, 3 - I]]", 1.04686, 1e-4);
    assert_close("Im[BesselY[0.5 I, 3 - I]]", 0.884784, 1e-4);
    /* Integer order, complex argument (general core). */
    assert_close("Re[BesselY[0, 3.0 + 0.0 I]]", eval_real("BesselY[0, 3.0]"), 1e-12);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_bessely_listable() {
    assert_eval_eq("BesselY[1, {a, b}]", "{BesselY[1, a], BesselY[1, b]}", 0);
    assert_eval_eq("BesselY[1, {}]", "{}", 0);
    assert_close("BesselY[0, {1.0, 2.0, 3.0}][[1]]",  0.0882569642, 1e-7);
    assert_close("BesselY[0, {1.0, 2.0, 3.0}][[2]]",  0.5103756726, 1e-7);
    assert_close("BesselY[0, {1.0, 2.0, 3.0}][[3]]",  0.3768500100, 1e-7);
    /* Threading also over the order. */
    assert_eval_eq("BesselY[{0, 1}, x]", "{BesselY[0, x], BesselY[1, x]}", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_bessely_derivative() {
    assert_eval_eq("D[BesselY[n, x], x]",
                   "1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])", 0);
    /* Chain rule through the argument. */
    assert_eval_eq("D[BesselY[n, x^2], x]",
                   "x (BesselY[-1 + n, x^2] - BesselY[1 + n, x^2])", 0);
    /* No x-dependence -> 0. */
    assert_eval_eq("D[BesselY[n, y], x]", "0", 0);
}

/* ---- Series at 0 (logarithmic) -------------------------------------- */
/* The integer-order expansion carries Log[x] and EulerGamma. Rather than pin
 * the (messy) printed form, verify the truncated series reproduces the numeric
 * value: Normal[Series[Y_p,{x,0,m}]] /. x->x0  ==  Y_p(x0) for small x0. */
void test_bessely_series_zero() {
    assert_close(
        "(Normal[Series[BesselY[0, x], {x, 0, 8}]] /. x -> 0.35) - BesselY[0, 0.35]",
        0.0, 1e-7);
    assert_close(
        "(Normal[Series[BesselY[1, x], {x, 0, 8}]] /. x -> 0.35) - BesselY[1, 0.35]",
        0.0, 1e-7);
    assert_close(
        "(Normal[Series[BesselY[2, x], {x, 0, 9}]] /. x -> 0.45) - BesselY[2, 0.45]",
        0.0, 1e-7);
    assert_close(
        "(Normal[Series[BesselY[3, x], {x, 0, 10}]] /. x -> 0.55) - BesselY[3, 0.55]",
        0.0, 1e-7);
}

/* ---- Series at Infinity --------------------------------------------- */

void test_bessely_series_infinity() {
    /* Leading-order trig-prefactored asymptotic (DLMF 10.17.4). */
    assert_eval_eq("Series[BesselY[0, x], {x, Infinity, 1}]",
        "-(Sqrt[2/Pi] Sqrt[1/x] + O[1/x]^(5/2)) Sin[1/4 Pi - x]", 0);
    /* Higher-order asymptotic reproduces the numeric value at large x. */
    assert_close(
        "(Normal[Series[BesselY[0, x], {x, Infinity, 8}]] /. x -> 50.0) - BesselY[0, 50.0]",
        0.0, 1e-10);
    assert_close(
        "(Normal[Series[BesselY[2, x], {x, Infinity, 8}]] /. x -> 30.0) - BesselY[2, 30.0]",
        0.0, 1e-9);
}

/* ---- stress --------------------------------------------------------- */

void test_bessely_stress() {
#ifdef USE_MPFR
    /* Deep precision: Y_0(4) at 70 digits must match mpfr_yn. */
    {
        Expr* e = parse_expression("BesselY[0, 4.0`70]");
        Expr* r = evaluate(e);
        expr_free(e);
        ASSERT_MSG(r->type == EXPR_MPFR, "BesselY[0,4.0`70] not MPFR");
        mpfr_t ref, x4;
        mpfr_init2(ref, mpfr_get_prec(r->data.mpfr));
        mpfr_init2(x4, mpfr_get_prec(r->data.mpfr));
        mpfr_set_ui(x4, 4, MPFR_RNDN);
        mpfr_yn(ref, 0, x4, MPFR_RNDN);
        mpfr_sub(ref, ref, r->data.mpfr, MPFR_RNDN);
        mpfr_abs(ref, ref, MPFR_RNDN);
        double err = mpfr_get_d(ref, MPFR_RNDN);
        ASSERT_MSG(err < 1e-60, "BesselY[0,4.0`70] error %.3e too large", err);
        mpfr_clear(ref); mpfr_clear(x4);
        expr_free(r);
    }
#endif
    /* Large argument (deep asymptotic regime). */
    assert_close("BesselY[0, 100.0]", -0.0772443133, 1e-9);
    assert_close("BesselY[3, 80.0]", -0.0665281519379, 1e-8);
    /* Near the first zero of Y0 (0.8935769663): value ~ 0. */
    assert_close("BesselY[0, 0.8935769662791675]", 0.0, 1e-12);
}

/* ---- argument count ------------------------------------------------- */

void test_bessely_argcount() {
    assert_eval_eq("BesselY[]", "BesselY[]", 0);
    assert_eval_eq("BesselY[1, 2, 3]", "BesselY[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_bessely_attributes() {
    SymbolDef* d = symtab_get_def("BesselY");
    ASSERT_MSG(d != NULL, "BesselY not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "BesselY not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "BesselY not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "BesselY not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_bessely_exact);
    TEST(test_bessely_no_parity);
    TEST(test_bessely_machine_real);
    TEST(test_bessely_arbitrary_real);
    TEST(test_bessely_mpfr_oracle);
    TEST(test_bessely_halfinteger_oracle);
    TEST(test_bessely_recurrence);
    TEST(test_bessely_branch_cut);
    TEST(test_bessely_complex);
    TEST(test_bessely_listable);
    TEST(test_bessely_derivative);
    TEST(test_bessely_series_zero);
    TEST(test_bessely_series_infinity);
    TEST(test_bessely_stress);
    TEST(test_bessely_argcount);
    TEST(test_bessely_attributes);

    printf("All BesselY tests passed.\n");
    return 0;
}
