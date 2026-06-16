/* Tests for BesselK[n, z], the modified Bessel function of the second kind.
 *
 * MPFR has no modified-Bessel routine, so there is no mpfr_kn oracle. Instead
 * correctness is pinned by: reference constants from the Wolfram Language docs;
 * a closed-form oracle for half-integer order (K_{1/2}, K_{3/2}); a recurrence-
 * consistency oracle that exercises both the small-|z| series and the large-|z|
 * asymptotic regimes for integer and non-integer order; complex order/argument;
 * precision tracking; list threading; the derivative rule; Series at 0 (the
 * logarithmic expansion) and at Infinity; special values at the origin;
 * argument-count diagnostics; attributes.
 *
 * NB: the symbolic half-integer -> elementary rewrites and the K_{-n}=K_n
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

void test_besselk_origin() {
    /* K_0(0) = Infinity (real, directed); K_n(0) = ComplexInfinity. */
    assert_eval_eq("BesselK[0, 0]", "Infinity", 0);
    assert_eval_eq("BesselK[1, 0]", "ComplexInfinity", 0);
    assert_eval_eq("BesselK[2, 0]", "ComplexInfinity", 0);
    /* Non-integer and symbolic orders also diverge -> ComplexInfinity. */
    assert_eval_eq("BesselK[1/2, 0]", "ComplexInfinity", 0);
    assert_eval_eq("BesselK[-1/2, 0]", "ComplexInfinity", 0);
    assert_eval_eq("BesselK[a, 0]", "ComplexInfinity", 0);
    /* Pure-imaginary order oscillates boundedly -> Indeterminate. */
    assert_eval_eq("BesselK[I, 0]", "Indeterminate", 0);
    /* Exact non-zero / symbolic arguments stay symbolic. */
    assert_eval_eq("BesselK[0, 2]", "BesselK[0, 2]", 0);
    assert_eval_eq("BesselK[n, x]", "BesselK[n, x]", 0);
    /* Half-integer order with exact numeric argument stays unevaluated. */
    assert_eval_eq("BesselK[11/2, 1]", "BesselK[11/2, 1]", 0);
}

/* ---- no argument-parity fold: K_n has a branch cut on z < 0 --------- */

void test_besselk_no_argument_parity() {
    /* Unlike BesselJ/BesselI, K_n(-z) carries no clean parity (branch cut
     * along the negative real axis), so it must stay unevaluated. */
    assert_eval_eq("BesselK[0, -z]", "BesselK[0, -z]", 0);
    assert_eval_eq("BesselK[1, -z]", "BesselK[1, -z]", 0);
    assert_eval_eq("BesselK[2, -z]", "BesselK[2, -z]", 0);
}

/* ---- machine-precision real (reference constants) ------------------- */

void test_besselk_machine_real() {
    assert_close("BesselK[0, 0.53]", 0.87656038041648568, 1e-7);
    assert_close("BesselK[0, 4.0]", 0.011159676085853024, 1e-9);
    assert_close("BesselK[1.0, 4.0]", 0.012483498887268431, 1e-9);
    assert_close("BesselK[2, 1.0]", 1.6248388986351774, 1e-7);
    assert_close("BesselK[3, 1.0]", 7.1012628247379448, 1e-6);
    /* Large argument -> asymptotic regime. */
    assert_close("BesselK[0, 30.0]", 2.132477496463056e-14, 1e-22);
    /* Small argument -> logarithmic series dominates. */
    assert_close("BesselK[0, 0.001]", 7.0236888005623813, 1e-6);
    /* Non-integer (real) order. */
    assert_close("BesselK[0.5, 2.0]", 0.11993777196806145, 1e-9);
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_besselk_arbitrary_real() {
    /* N[BesselK[0,4],50] matches the WL reference to all 50 digits. */
    assert_eval_startswith("N[BesselK[0, 4], 50]",
        "0.011159676085853024269745195979833489225009023888474");
    /* Precision tracks the input (37-digit argument). */
    assert_eval_startswith("BesselK[1, 4.000000000000000000000000000000000000]",
        "0.01248349888726843147038417998080606848");
    /* Half-integer at high precision via the connection formula. */
    assert_eval_startswith("N[BesselK[1/2, 2], 40]",
        "0.119937771968061447368036501636793516");
}

/* ---- closed-form oracle: half-integer order ------------------------- */
/* K_{1/2}(x) = Sqrt[pi/(2x)] e^{-x},  K_{3/2}(x) = Sqrt[pi/(2x)] e^{-x}(1+1/x).
 * The grid spans the connection-formula regime (small |x|) and the asymptotic
 * regime (|x| > ~24 at machine precision), validating both numeric paths
 * against a formula no library provides. */
void test_besselk_halfinteger_oracle() {
    for (int t = 1; t <= 18; t++) {
        double x = t * 1.5;                 /* 1.5 .. 27 */
        char in[80];
        double pref = sqrt(M_PI / (2.0 * x)) * exp(-x);

        snprintf(in, sizeof(in), "BesselK[1/2, %.6f]", x);
        assert_rel(in, pref, 1e-9);

        snprintf(in, sizeof(in), "BesselK[3/2, %.6f]", x);
        assert_rel(in, pref * (1.0 + 1.0 / x), 1e-9);
    }
}

/* ---- recurrence-consistency oracle ---------------------------------- */
/* DLMF 10.29.1: K_{nu+1}(x) - K_{nu-1}(x) = (2 nu / x) K_nu(x). Checked for an
 * integer order (nu=1, exercising the logarithmic series / asymptotic integer
 * path) and a generic non-integer order (nu=1/3, the connection formula),
 * across a grid spanning both the small- and large-|x| regimes. */
void test_besselk_recurrence() {
    for (int t = 1; t <= 22; t++) {
        double x = t * 1.4;                 /* 1.4 .. 30.8 */
        char lhs[112], rhs[112];

        /* nu = 1 (integer). */
        snprintf(lhs, sizeof(lhs), "BesselK[2, %.6f] - BesselK[0, %.6f]", x, x);
        snprintf(rhs, sizeof(rhs), "(2/%.6f) BesselK[1, %.6f]", x, x);
        double l = eval_mag(lhs), r = eval_mag(rhs);
        ASSERT_MSG(fabs(l - r) <= 1e-9 * (1.0 + fabs(r)),
                   "K recurrence nu=1 at x=%.4f: %.14g vs %.14g", x, l, r);

        /* nu = 1/3 (non-integer). */
        snprintf(lhs, sizeof(lhs), "BesselK[4/3, %.6f] - BesselK[-2/3, %.6f]", x, x);
        snprintf(rhs, sizeof(rhs), "(2/3) BesselK[1/3, %.6f] / %.6f", x, x);
        l = eval_mag(lhs); r = eval_mag(rhs);
        ASSERT_MSG(fabs(l - r) <= 1e-9 * (1.0 + fabs(r)),
                   "K recurrence nu=1/3 at x=%.4f: %.14g vs %.14g", x, l, r);
    }
}

/* ---- complex order and argument ------------------------------------- */

void test_besselk_complex() {
    /* Spec value: BesselK[1 + I, 3.0 - 2 I] = -0.0225108 + 0.0169607 I. */
    assert_close("Re[BesselK[1 + I, 3.0 - 2 I]]", -0.0225108, 1e-5);
    assert_close("Im[BesselK[1 + I, 3.0 - 2 I]]",  0.0169607, 1e-5);
    /* Real argument given as Complex with 0 imaginary part matches the real
     * evaluation. */
    assert_close("Re[BesselK[0, 4.0 + 0.0 I]]", eval_mag("BesselK[0, 4.0]"), 1e-12);
    /* Half-integer closed form K_{1/2}(z) = Sqrt[pi/(2z)] e^{-z} for a complex
     * argument off the branch cut (all sqrt branches agree there): a clean,
     * convention-independent oracle for the complex connection-formula core.
     * Compared part-by-part (Re/Im each give a clean real). */
    assert_close("Re[BesselK[1/2, 2.0 + 1.0 I]] - "
                 "Re[Sqrt[Pi/(2 (2.0 + 1.0 I))] Exp[-(2.0 + 1.0 I)]]", 0.0, 1e-12);
    assert_close("Im[BesselK[1/2, 2.0 + 1.0 I]] - "
                 "Im[Sqrt[Pi/(2 (2.0 + 1.0 I))] Exp[-(2.0 + 1.0 I)]]", 0.0, 1e-12);
    assert_close("Re[BesselK[3/2, 1.0 + 2.0 I]] - Re[Sqrt[Pi/(2 (1.0 + 2.0 I))] "
                 "Exp[-(1.0 + 2.0 I)] (1 + 1/(1.0 + 2.0 I))]", 0.0, 1e-12);
    assert_close("Im[BesselK[3/2, 1.0 + 2.0 I]] - Im[Sqrt[Pi/(2 (1.0 + 2.0 I))] "
                 "Exp[-(1.0 + 2.0 I)] (1 + 1/(1.0 + 2.0 I))]", 0.0, 1e-12);
    /* Negative real argument lies on the branch cut: the principal value is
     * complex. Its modulus is branch-independent: |K_{1/2}(-2)| =
     * Sqrt[pi/4] e^2 = 6.54759..., a clean check that does not depend on the
     * (convention-dependent) sign of the imaginary part. */
    assert_close("Abs[BesselK[1/2, -2.0]]", sqrt(M_PI / 4.0) * exp(2.0), 1e-4);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_besselk_listable() {
    assert_eval_eq("BesselK[1, {a, b}]", "{BesselK[1, a], BesselK[1, b]}", 0);
    assert_eval_eq("BesselK[1, {}]", "{}", 0);
    /* Spec: BesselK[{1,2,3},1.0] = {0.601907, 1.62484, 7.10126}. */
    assert_close("BesselK[{1, 2, 3}, 1.0][[1]]", 0.6019072302, 1e-7);
    assert_close("BesselK[{1, 2, 3}, 1.0][[2]]", 1.6248388986, 1e-7);
    assert_close("BesselK[{1, 2, 3}, 1.0][[3]]", 7.1012628247, 1e-6);
    /* Matrix threading. */
    assert_eval_eq("BesselK[{0, 1}, x]", "{BesselK[0, x], BesselK[1, x]}", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_besselk_derivative() {
    /* Sign differs from BesselJ: -(K_{n-1} + K_{n+1})/2. */
    assert_eval_eq("D[BesselK[n, x], x]",
                   "-1/2 (BesselK[-1 + n, x] + BesselK[1 + n, x])", 0);
    /* Chain rule through the argument. */
    assert_eval_eq("D[BesselK[n, x^2], x]",
                   "-x (BesselK[-1 + n, x^2] + BesselK[1 + n, x^2])", 0);
    /* No x-dependence -> 0. */
    assert_eval_eq("D[BesselK[n, y], x]", "0", 0);
}

/* ---- Series at 0 (logarithmic) -------------------------------------- */

void test_besselk_series_zero() {
    assert_eval_eq("Series[BesselK[0, x], {x, 0, 3}]",
        "1/2 (-2 EulerGamma + 2 Log[2] - 2 Log[x]) + "
        "1/8 (2 - 2 EulerGamma + 2 Log[2] - 2 Log[x]) x^2 + O[x]^4", 0);
    /* Integer order >= 1 carries a pole and a logarithmic part. */
    assert_eval_eq("Series[BesselK[1, x], {x, 0, 1}]",
        "1/x + -1/4 (1 - 2 EulerGamma + 2 Log[2] - 2 Log[x]) x + O[x]^2", 0);
    assert_eval_eq("Series[BesselK[2, x], {x, 0, 0}]",
        "2/x^2 - 1/2 + O[x]^1", 0);
}

/* ---- Series at Infinity --------------------------------------------- */

void test_besselk_series_infinity() {
    /* Structural comparison: the coefficients -1/8, 9/128, -75/1024 are the
     * DLMF 10.40.2 a_k * Sqrt[Pi/2]. (The Exp[-x] * SeriesData product keeps
     * its O-term under Normal, so compare the raw Series.) */
    assert_eval_eq("Series[BesselK[0, x], {x, Infinity, 4}]",
        "E^(-x) (Sqrt[1/2 Pi] Sqrt[1/x] + -1/8 Sqrt[1/2 Pi] (1/x)^(3/2) + "
        "9/128 Sqrt[1/2 Pi] (1/x)^(5/2) + -75/1024 Sqrt[1/2 Pi] (1/x)^(7/2) + "
        "O[1/x]^(9/2))", 0);
    /* Limit at infinity is 0 (falls out of the asymptotic series). */
    assert_eval_eq("Limit[BesselK[n, x], x -> Infinity]", "0", 0);
}

/* ---- stress --------------------------------------------------------- */

void test_besselk_stress() {
#ifdef USE_MPFR
    /* Deep precision: the recurrence K_2(4) - K_0(4) = (1/2) K_1(4) must hold to
     * ~60 digits, validating that all three integer-order log-series agree at
     * high precision (an internal-consistency oracle, no library needed). */
    {
        Expr* e = parse_expression(
            "Abs[BesselK[2, 4`80] - BesselK[0, 4`80] - (1/2) BesselK[1, 4`80]]");
        Expr* r = evaluate(e);
        expr_free(e);
        double err;
        if (r->type == EXPR_MPFR)      err = fabs(mpfr_get_d(r->data.mpfr, MPFR_RNDN));
        else if (r->type == EXPR_REAL) err = fabs(r->data.real);
        else if (r->type == EXPR_INTEGER) err = fabs((double)r->data.integer);
        else { ASSERT_MSG(0, "deep-precision recurrence: non-numeric result"); err = 1; }
        ASSERT_MSG(err < 1e-60, "deep-precision K recurrence error %.3e", err);
        expr_free(r);
    }
    /* A 500-bit evaluation completes and tracks precision. */
    {
        Expr* e = parse_expression("BesselK[0, 4`150]");
        Expr* r = evaluate(e);
        expr_free(e);
        ASSERT_MSG(r->type == EXPR_MPFR, "BesselK[0,4`150] not MPFR");
        assert_eval_startswith("BesselK[0, 4`150]",
            "0.0111596760858530242697451959798334892250090238884743405382552615");
        expr_free(r);
    }
#endif
    /* Large order, small argument: huge but finite and accurate. */
    assert_rel("BesselK[10, 1.0]", 1.807132899010e8, 1e-9);
    /* Large argument (deep asymptotic regime). */
    assert_rel("BesselK[0, 50.0]", 3.410167749789e-23, 1e-9);
    assert_rel("BesselK[3, 40.0]", 9.378903724645e-19, 1e-9);
}

/* ---- argument count ------------------------------------------------- */

void test_besselk_argcount() {
    assert_eval_eq("BesselK[]", "BesselK[]", 0);
    assert_eval_eq("BesselK[1, 2, 3]", "BesselK[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_besselk_attributes() {
    SymbolDef* d = symtab_get_def("BesselK");
    ASSERT_MSG(d != NULL, "BesselK not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "BesselK not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "BesselK not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "BesselK not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_besselk_origin);
    TEST(test_besselk_no_argument_parity);
    TEST(test_besselk_machine_real);
    TEST(test_besselk_arbitrary_real);
    TEST(test_besselk_halfinteger_oracle);
    TEST(test_besselk_recurrence);
    TEST(test_besselk_complex);
    TEST(test_besselk_listable);
    TEST(test_besselk_derivative);
    TEST(test_besselk_series_zero);
    TEST(test_besselk_series_infinity);
    TEST(test_besselk_stress);
    TEST(test_besselk_argcount);
    TEST(test_besselk_attributes);

    printf("All BesselK tests passed.\n");
    return 0;
}
