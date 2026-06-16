/* Tests for BesselI[n, z], the modified Bessel function of the first kind I_n(z).
 *
 * MPFR has no modified-Bessel routine, so there is no mpfr_in oracle. Instead
 * correctness is pinned by: reference constants from the Wolfram Language docs;
 * a closed-form oracle for half-integer order (I_{1/2}, I_{-1/2}, I_{3/2}) that
 * spans both the small-|z| power series and the large-|z| asymptotic regimes,
 * for real AND complex argument; a recurrence-consistency oracle (DLMF 10.29.1)
 * for integer and non-integer order across both regimes; the spec complex value;
 * precision tracking; list threading; the derivative rule; Series at 0 and at
 * Infinity (the two-exponential form); special values at the origin;
 * argument-count diagnostics; attributes.
 *
 * NB: the symbolic half-integer -> elementary rewrites and the I_{-n}=I_n
 * reflection live in src/internal/bessel.m, which is loaded by the REPL but NOT
 * by this harness (no init.m). Half-integer calls here therefore use the
 * inexact numeric C path; symbolic reflection/rewrites are verified separately
 * in the REPL. */

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

/* Evaluate to a double, accepting Integer / Real / MPFR results. */
static double eval_mag(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    double v;
    if (r->type == EXPR_REAL)          v = r->data.real;
    else if (r->type == EXPR_INTEGER)  v = (double)r->data.integer;
#ifdef USE_MPFR
    else if (r->type == EXPR_MPFR)     v = mpfr_get_d(r->data.mpfr, MPFR_RNDN);
#endif
    else { ASSERT_MSG(0, "%s: expected a numeric result", input); v = 0.0; }
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_mag(input);
    ASSERT_MSG(fabs(v - expected) <= tol,
               "%s: expected %.12g, got %.12g", input, expected, v);
}

/* Relative-error check, for oracle values spanning a wide dynamic range. */
static void assert_rel(const char* input, double want, double rel) {
    double got = eval_mag(input);
    double err = fabs(got - want) / (fabs(want) > 0 ? fabs(want) : 1.0);
    ASSERT_MSG(err <= rel, "%s: builtin %.14g vs oracle %.14g (rel %.3e)",
               input, got, want, err);
}

/* ---- special values at the origin ----------------------------------- */

void test_besseli_origin() {
    /* I_0(0) = 1, I_n(0) = 0 for integer n != 0. */
    assert_eval_eq("BesselI[0, 0]", "1", 0);
    assert_eval_eq("BesselI[1, 0]", "0", 0);
    assert_eval_eq("BesselI[3, 0]", "0", 0);
    assert_eval_eq("BesselI[-2, 0]", "0", 0);
    /* Exact non-zero / symbolic arguments stay symbolic. */
    assert_eval_eq("BesselI[0, 2]", "BesselI[0, 2]", 0);
    assert_eval_eq("BesselI[n, x]", "BesselI[n, x]", 0);
    /* Half-integer order with exact numeric argument stays unevaluated. */
    assert_eval_eq("BesselI[11/2, 1]", "BesselI[11/2, 1]", 0);
}

/* ---- argument parity: I_n(-z) = (-1)^n I_n(z), integer n only ------- */

void test_besseli_argument_parity() {
    assert_eval_eq("BesselI[0, -z]", "BesselI[0, z]", 0);
    assert_eval_eq("BesselI[2, -z]", "BesselI[2, z]", 0);
    assert_eval_eq("BesselI[1, -z]", "-BesselI[1, z]", 0);
    assert_eval_eq("BesselI[3, -z]", "-BesselI[3, z]", 0);
    assert_eval_eq("BesselI[1, -2 x]", "-BesselI[1, 2 x]", 0);
    assert_eval_eq("BesselI[1, -2]", "-BesselI[1, 2]", 0);
    /* Non-integer / symbolic order does NOT fold. */
    assert_eval_eq("BesselI[n, -z]", "BesselI[n, -z]", 0);
    assert_eval_eq("BesselI[1, z]", "BesselI[1, z]", 0);
}

/* ---- machine-precision real (reference constants) ------------------- */

void test_besseli_machine_real() {
    assert_close("BesselI[0, 1.0]", 1.2660658777520084, 1e-12);
    assert_close("BesselI[1, 1.0]", 0.5651591039924851, 1e-12);
    assert_close("BesselI[0, 2.0]", 2.2795853023360673, 1e-12);
    assert_close("BesselI[1, 2.0]", 1.5906368546373291, 1e-12);
    assert_close("BesselI[0, 4.0]", 11.301921952136330, 1e-10);
    assert_close("BesselI[2, 3.0]", 2.2452124409299512, 1e-12);
    /* Small argument. */
    assert_close("BesselI[0, 0.5]", 1.0634833707413235, 1e-12);
    /* Even/odd symmetry in the argument (integer order). */
    assert_close("BesselI[0, -2.0]", 2.2795853023360673, 1e-12);  /* I0 even */
    assert_close("BesselI[1, -1.0]", -0.5651591039924851, 1e-12); /* I1 odd  */
    /* Large argument -> asymptotic regime. */
    assert_rel("BesselI[0, 20.0]", 43558282.55955353, 1e-10);
    assert_rel("BesselI[0, 50.0]", 2.9325537838493363e20, 1e-9);
    /* Non-integer (real) order; spec: I_{1/2}(2) = 2.0462368630890552. */
    assert_close("BesselI[0.5, 2.0]", 2.0462368630890552, 1e-12);
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_besseli_arbitrary_real() {
    /* N[BesselI[0,4],65] matches the WL reference to all printed digits. */
    assert_eval_startswith("N[BesselI[0, 4], 65]",
        "11.301921952136330496356270183217102497412616594435337706006496193");
    /* Precision tracks the input (37-digit argument). */
    assert_eval_startswith("BesselI[0, 4.000000000000000000000000000000000000]",
        "11.3019219521363304963562701832171025");
    /* Half-integer at high precision; spec: I_{1/2}(2) to 40 digits. */
    assert_eval_startswith("N[BesselI[1/2, 2], 40]",
        "2.046236863089055036605183612020732319267");
}

/* ---- closed-form oracle: half-integer order ------------------------- */
/* I_{1/2}(x)  = Sqrt[2/(pi x)] Sinh[x],
 * I_{-1/2}(x) = Sqrt[2/(pi x)] Cosh[x],
 * I_{3/2}(x)  = Sqrt[2/(pi x)] (Cosh[x] - Sinh[x]/x).
 * The grid spans the power-series regime (small |x|) and the asymptotic regime
 * (|x| > ~24 at machine precision), validating both numeric paths against a
 * formula no library provides. */
void test_besseli_halfinteger_oracle() {
    for (int t = 1; t <= 18; t++) {
        double x = t * 1.5;                 /* 1.5 .. 27 */
        char in[80];
        double pref = sqrt(2.0 / (M_PI * x));

        snprintf(in, sizeof(in), "BesselI[1/2, %.6f]", x);
        assert_rel(in, pref * sinh(x), 1e-9);

        snprintf(in, sizeof(in), "BesselI[-1/2, %.6f]", x);
        assert_rel(in, pref * cosh(x), 1e-9);

        snprintf(in, sizeof(in), "BesselI[3/2, %.6f]", x);
        assert_rel(in, pref * (cosh(x) - sinh(x) / x), 1e-9);
    }
}

/* ---- recurrence-consistency oracle ---------------------------------- */
/* DLMF 10.29.1: I_{nu-1}(x) - I_{nu+1}(x) = (2 nu / x) I_nu(x). Checked for an
 * integer order (nu=1, the integer power-series / asymptotic path) and a generic
 * non-integer order (nu=1/3), across a grid spanning both the small- and
 * large-|x| regimes. */
void test_besseli_recurrence() {
    for (int t = 1; t <= 22; t++) {
        double x = t * 1.4;                 /* 1.4 .. 30.8 */
        char lhs[112], rhs[112];

        /* nu = 1 (integer): I_0 - I_2 = (2/x) I_1. */
        snprintf(lhs, sizeof(lhs), "BesselI[0, %.6f] - BesselI[2, %.6f]", x, x);
        snprintf(rhs, sizeof(rhs), "(2/%.6f) BesselI[1, %.6f]", x, x);
        double l = eval_mag(lhs), r = eval_mag(rhs);
        ASSERT_MSG(fabs(l - r) <= 1e-9 * (1.0 + fabs(r)),
                   "I recurrence nu=1 at x=%.4f: %.14g vs %.14g", x, l, r);

        /* nu = 1/3 (non-integer): I_{-2/3} - I_{4/3} = (2/3) I_{1/3} / x. */
        snprintf(lhs, sizeof(lhs), "BesselI[-2/3, %.6f] - BesselI[4/3, %.6f]", x, x);
        snprintf(rhs, sizeof(rhs), "(2/3) BesselI[1/3, %.6f] / %.6f", x, x);
        l = eval_mag(lhs); r = eval_mag(rhs);
        ASSERT_MSG(fabs(l - r) <= 1e-9 * (1.0 + fabs(r)),
                   "I recurrence nu=1/3 at x=%.4f: %.14g vs %.14g", x, l, r);
    }
}

/* ---- complex order and argument ------------------------------------- */

void test_besseli_complex() {
    /* Spec value: BesselI[3 + I, 1.5 - I] = -0.25665 + 0.0492771 I. */
    assert_close("Re[BesselI[3 + I, 1.5 - I]]", -0.25665, 1e-5);
    assert_close("Im[BesselI[3 + I, 1.5 - I]]",  0.0492771, 1e-5);
    /* Real argument given as Complex with 0 imaginary part matches the real
     * evaluation. */
    assert_close("Re[BesselI[0, 4.0 + 0.0 I]]", eval_mag("BesselI[0, 4.0]"), 1e-12);
    /* Half-integer closed form I_{1/2}(z) = Sqrt[2/(pi z)] Sinh[z] for a complex
     * argument: an implementation-independent oracle. Small |z| exercises the
     * power-series core, large |z| the asymptotic core. Relative comparison. */
    assert_close("Abs[BesselI[1/2, 2.0 + 1.0 I] - Sqrt[2/(Pi (2.0 + 1.0 I))] "
                 "Sinh[2.0 + 1.0 I]] / Abs[Sqrt[2/(Pi (2.0 + 1.0 I))] "
                 "Sinh[2.0 + 1.0 I]]", 0.0, 1e-12);
    assert_close("Abs[BesselI[1/2, 25.0 + 5.0 I] - Sqrt[2/(Pi (25.0 + 5.0 I))] "
                 "Sinh[25.0 + 5.0 I]] / Abs[Sqrt[2/(Pi (25.0 + 5.0 I))] "
                 "Sinh[25.0 + 5.0 I]]", 0.0, 1e-10);
    /* Recurrence consistency on a large complex |z| (the asymptotic path):
     * I_0(z) - I_2(z) = (2/z) I_1(z). */
    assert_close("Abs[(BesselI[0, 20.0 + 8.0 I] - BesselI[2, 20.0 + 8.0 I]) - "
                 "(2/(20.0 + 8.0 I)) BesselI[1, 20.0 + 8.0 I]] / "
                 "Abs[(2/(20.0 + 8.0 I)) BesselI[1, 20.0 + 8.0 I]]", 0.0, 1e-10);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_besseli_listable() {
    assert_eval_eq("BesselI[1, {a, b}]", "{BesselI[1, a], BesselI[1, b]}", 0);
    assert_eval_eq("BesselI[1, {}]", "{}", 0);
    /* Spec: BesselI[{0,1,2}, 1.] = {1.26607, 0.565159, 0.135748}. */
    assert_close("BesselI[{0, 1, 2}, 1.0][[1]]", 1.2660658777520084, 1e-9);
    assert_close("BesselI[{0, 1, 2}, 1.0][[2]]", 0.5651591039924851, 1e-9);
    assert_close("BesselI[{0, 1, 2}, 1.0][[3]]", 0.1357476697670383, 1e-9);
    /* Matrix threading. */
    assert_eval_eq("BesselI[{0, 1}, x]", "{BesselI[0, x], BesselI[1, x]}", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_besseli_derivative() {
    /* Like BesselK the two-term sum carries a '+', but the coefficient is +1/2
     * (vs BesselK's -1/2). */
    assert_eval_eq("D[BesselI[n, x], x]",
                   "1/2 (BesselI[-1 + n, x] + BesselI[1 + n, x])", 0);
    /* Chain rule through the argument. */
    assert_eval_eq("D[BesselI[n, x^2], x]",
                   "x (BesselI[-1 + n, x^2] + BesselI[1 + n, x^2])", 0);
    /* No x-dependence -> 0. */
    assert_eval_eq("D[BesselI[n, y], x]", "0", 0);
}

/* ---- Series at 0 (Taylor) ------------------------------------------- */

void test_besseli_series_zero() {
    /* Falls out of the generic Taylor path (derivative rule + I_n(0) values). */
    assert_eval_eq("Series[BesselI[0, x], {x, 0, 6}]",
        "1 + 1/4 x^2 + 1/64 x^4 + 1/2304 x^6 + O[x]^7", 0);
}

/* ---- Series at Infinity (two-exponential, DLMF 10.40.5) ------------- */

void test_besseli_series_infinity() {
    /* Reproduces Mathematica's Normal form: a dominant E^x series with
     * (-1)^k a_k sqrt(1/(2pi)) coefficients plus a subdominant I E^{-x} series
     * with a_k sqrt(1/(2pi)) coefficients (n = 0 here, so E^{I n Pi} = 1). */
    assert_eval_eq("Series[BesselI[0, x], {x, Infinity, 2}]",
        "E^x (Sqrt[1/x]/Sqrt[2 Pi] + 1/8/Sqrt[2 Pi] (1/x)^(3/2) + "
        "O[1/x]^(5/2)) + I E^(-x) (Sqrt[1/x]/Sqrt[2 Pi] + "
        "-1/8/Sqrt[2 Pi] (1/x)^(3/2) + O[1/x]^(5/2))", 0);
}

/* ---- stress --------------------------------------------------------- */

void test_besseli_stress() {
#ifdef USE_MPFR
    /* Deep precision: I_0(4) - I_2(4) = (1/2) I_1(4) (DLMF 10.29.1, nu=1) must
     * hold to ~60 digits -- an internal-consistency oracle, no library needed. */
    {
        Expr* e = parse_expression(
            "Abs[BesselI[0, 4`80] - BesselI[2, 4`80] - (1/2) BesselI[1, 4`80]]");
        Expr* r = evaluate(e);
        expr_free(e);
        double err;
        if (r->type == EXPR_MPFR)      err = fabs(mpfr_get_d(r->data.mpfr, MPFR_RNDN));
        else if (r->type == EXPR_REAL) err = fabs(r->data.real);
        else if (r->type == EXPR_INTEGER) err = fabs((double)r->data.integer);
        else { ASSERT_MSG(0, "deep-precision recurrence: non-numeric result"); err = 1; }
        ASSERT_MSG(err < 1e-60, "deep-precision I recurrence error %.3e", err);
        expr_free(r);
    }
    /* A 150-bit evaluation completes and tracks precision. */
    {
        Expr* e = parse_expression("BesselI[0, 4`150]");
        Expr* r = evaluate(e);
        expr_free(e);
        ASSERT_MSG(r->type == EXPR_MPFR, "BesselI[0,4`150] not MPFR");
        assert_eval_startswith("BesselI[0, 4`150]",
            "11.30192195213633049635627018321710249741261659443533770600649619");
        expr_free(r);
    }
    /* Deep precision in the asymptotic regime stays real and accurate. */
    assert_eval_startswith("BesselI[0, 30`60]",
        "781672297823.977");
#endif
    /* Large order, small argument: tiny but accurate. */
    assert_rel("BesselI[10, 1.0]", 2.7529480398368736e-10, 1e-8);
    /* Large argument (deep asymptotic regime). */
    assert_rel("BesselI[0, 50.0]", 2.9325537838493363e20, 1e-9);
    assert_rel("BesselI[3, 40.0]", 1.3291455664733660e16, 1e-9);
    /* Large negative argument: integer order stays real (I_n even/odd). */
    assert_rel("BesselI[0, -50.0]", 2.9325537838493363e20, 1e-9);
    assert_rel("BesselI[1, -50.0]", -2.9030785901035568e20, 1e-9);
}

/* ---- argument count ------------------------------------------------- */

void test_besseli_argcount() {
    assert_eval_eq("BesselI[]", "BesselI[]", 0);
    assert_eval_eq("BesselI[1, 2, 3]", "BesselI[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_besseli_attributes() {
    SymbolDef* d = symtab_get_def("BesselI");
    ASSERT_MSG(d != NULL, "BesselI not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "BesselI not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "BesselI not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "BesselI not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_besseli_origin);
    TEST(test_besseli_argument_parity);
    TEST(test_besseli_machine_real);
    TEST(test_besseli_arbitrary_real);
    TEST(test_besseli_halfinteger_oracle);
    TEST(test_besseli_recurrence);
    TEST(test_besseli_complex);
    TEST(test_besseli_listable);
    TEST(test_besseli_derivative);
    TEST(test_besseli_series_zero);
    TEST(test_besseli_series_infinity);
    TEST(test_besseli_stress);
    TEST(test_besseli_argcount);
    TEST(test_besseli_attributes);

    printf("All BesselI tests passed.\n");
    return 0;
}
