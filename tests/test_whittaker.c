/* Tests for WhittakerM[k, m, z] and WhittakerW[k, m, z], the Whittaker functions
 *   M_{k,m}(z) = e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z],
 *   W_{k,m}(z) = e^(-z/2) z^(m+1/2) HypergeometricU[m-k+1/2, 2m+1, z].
 *
 * Covers: the z = 0 special value and symbolic/exact passthrough; machine- and
 * arbitrary-precision (MPFR) real numerics; complex parameters and argument;
 * list threading; the NDArray / packed-array fast path (agreement with the
 * scalar path); Compile[]; the derivative rule; Series at the origin (concrete
 * and symbolic m, integer and half-integer leading power) and the terminating
 * case at Infinity; SeriesCoefficient; argument-count diagnostics; attributes.
 */

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

/* ---- special values and passthrough --------------------------------- */

void test_whittaker_special() {
    /* M_{k,m}(0) = 0 for Re(m) > -1/2 (leading power z^(m+1/2)). */
    assert_eval_eq("WhittakerM[1, 1, 0]", "0", 0);
    assert_eval_eq("WhittakerM[0, 0, 0]", "0", 0);
    assert_eval_eq("WhittakerM[2, 0.5, 0]", "0", 0);
    /* Symbolic / all-exact arguments stay symbolic (no expansion). */
    assert_eval_eq("WhittakerM[2, 3, x]", "WhittakerM[2, 3, x]", 0);
    assert_eval_eq("WhittakerM[2, 3, 4]", "WhittakerM[2, 3, 4]", 0);
    /* Symbolic m with z == 0 is left alone (sign of Re(m) unknown). */
    assert_eval_eq("WhittakerM[k, m, 0]", "WhittakerM[k, m, 0]", 0);
}

/* ---- machine real numerics ------------------------------------------ */

void test_whittaker_machine_real() {
    assert_close("WhittakerM[2, 3, 4.5]", 68.9670344349435886, 1e-6);
    assert_close("WhittakerM[2, 3, 1.7]", 4.07202, 1e-4);
    /* Rational parameters with an inexact argument. */
    assert_close("WhittakerM[2, 1/3, 0.5]", 0.286592, 1e-5);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_whittaker_arbitrary() {
    /* N[WhittakerM[2,1/3,1/7],100] -- first ~60 digits are rounding-stable. */
    assert_eval_startswith("N[WhittakerM[2, 1/3, 1/7], 100]",
        "0.16564977190667829031036482807198710437606805187130125676351652");
    /* Precision tracks the input. */
    assert_eval_startswith("WhittakerM[2, 1/3, 1/7`40]", "0.1656497719066782903103648");
}

/* ---- complex parameters and argument -------------------------------- */

void test_whittaker_complex() {
    /* N[WhittakerM[4/5, 2 - I, 2]] = 3.27958 - 3.25602 I. */
    assert_close("Re[N[WhittakerM[4/5, 2 - I, 2]]]", 3.27958, 1e-4);
    assert_close("Im[N[WhittakerM[4/5, 2 - I, 2]]]", -3.25602, 1e-4);
    /* WhittakerM[2, 3 - I, 4.5] = -18.2886 - 74.6279 I. */
    assert_close("Re[WhittakerM[2, 3 - I, 4.5]]", -18.2886002968, 1e-3);
    assert_close("Im[WhittakerM[2, 3 - I, 4.5]]", -74.6278856690, 1e-3);
    /* A negative real argument lands off the real axis (z^(m+1/2) branch cut). */
    assert_close("Re[WhittakerM[2, 3, -1.5]]", 0.0, 1e-10);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_whittaker_listable() {
    assert_eval_eq("WhittakerM[k, m, {a, b}]",
                   "{WhittakerM[k, m, a], WhittakerM[k, m, b]}", 0);
    assert_eval_eq("WhittakerM[k, m, {}]", "{}", 0);
    assert_close("WhittakerM[2, 3, {1.5, 2.5, 3.5}][[1]]", 2.76149, 1e-4);
    assert_close("WhittakerM[2, 3, {1.5, 2.5, 3.5}][[2]]", 13.0096, 1e-3);
    assert_close("WhittakerM[2, 3, {1.5, 2.5, 3.5}][[3]]", 34.2394, 1e-3);
}

/* ---- NDArray / packed surface agreement ----------------------------- */

void test_whittaker_ndarray() {
    /* Visible NDArray: element-wise, agreeing with the scalar path. */
    assert_close("WhittakerM[2, 3, NDArray[{1.5, 2.5, 3.5}]][[2]]", 13.0096, 1e-3);
    /* A large numeric list auto-packs; the packed fast path shares the machine
     * helper with the scalar path, so an element must match bit-closely. */
    assert_close("WhittakerM[2, 3, N[Range[60]/13]][[7]]",
                 eval_real("WhittakerM[2, 3, N[7/13]]"), 1e-12);
    assert_close("WhittakerM[2, 3, N[Range[60]/13]][[40]]",
                 eval_real("WhittakerM[2, 3, N[40/13]]"), 1e-12);
}

/* ---- Compile[] ------------------------------------------------------ */

void test_whittaker_compile() {
    /* Scalar form lowers (OP_KERNN via the N-ary machine kernel). */
    assert_eval_eq(
        "CompileDiagnostics[{{x, _Real}, {y, _Real}, {z, _Real}}, "
        "WhittakerM[x, y, z]][[1]]", "\"Compiled\" -> True", 0);
    /* The compiled function runs and agrees with the scalar path. */
    assert_close(
        "Compile[{{x, _Real}, {y, _Real}, {z, _Real}}, WhittakerM[x, y, z]][2., 3., 4.5]",
        68.9670344349435886, 1e-6);
}

/* ---- derivative ----------------------------------------------------- */

void test_whittaker_derivative() {
    assert_eval_eq("D[WhittakerM[k, m, z], z]",
        "(1/2 - k/z) WhittakerM[k, m, "
        "z] + ((1/2 + k + m) WhittakerM[1 + k, m, z])/z", 0);
    /* No elementary form for d/dk, d/dm -> inert Derivative. */
    assert_eval_eq("D[WhittakerM[k, m, z], k]",
                   "Derivative[1, 0, 0][WhittakerM][k, m, z]", 0);
    assert_eval_eq("D[WhittakerM[k, m, z], m]",
                   "Derivative[0, 1, 0][WhittakerM][k, m, z]", 0);
    /* No z-dependence -> 0. */
    assert_eval_eq("D[WhittakerM[k, m, y], z]", "0", 0);
}

/* ---- Series at the origin ------------------------------------------- */

void test_whittaker_series_zero() {
    /* Concrete m = 1/2: integer leading power. */
    assert_eval_eq("Series[WhittakerM[2, 1/2, x], {x, 0, 6}]",
        "x - x^2 + 3/8 x^3 - 1/12 x^4 + 5/384 x^5 - 1/640 x^6 + O[x]^7", 0);
    /* Concrete m = 1: half-integer (Puiseux) leading power. */
    assert_eval_eq("Series[WhittakerM[3, 1, x], {x, 0, 5}]",
        "x^(3/2) - x^(5/2) + 13/32 x^(7/2) - 47/480 x^(9/2) + O[x]^(11/2)", 0);
    /* Symbolic m: the z^(1/2+m) prefactor form. */
    assert_eval_startswith("Series[WhittakerM[k, m, x], {x, 0, 4}]",
        "x^(1/2 + m) (1 + ");
}

/* ---- Series at Infinity (terminating case only) --------------------- */

void test_whittaker_series_infinity() {
    /* a = m-k+1/2 = -1 is a non-positive integer -> elementary closed form. */
    assert_eval_eq("Series[WhittakerM[2, 1/2, x], {x, Infinity, 5}]",
                   "((1 - 1/2 x) x)/E^(1/2 x)", 0);
    /* a = 0 -> M_{1/2,0}(x) = Sqrt[x] E^(-x/2). */
    assert_eval_eq("Series[WhittakerM[1/2, 0, x], {x, Infinity, 3}]",
                   "Sqrt[x]/E^(1/2 x)", 0);
    /* Non-terminating / symbolic: no asymptotic machinery -> stays symbolic. */
    assert_eval_eq("Series[WhittakerM[1/3, 1/2, x], {x, Infinity, 3}]",
                   "Series[WhittakerM[1/3, 1/2, x], {x, Infinity, 3}]", 0);
}

/* ---- SeriesCoefficient ---------------------------------------------- */

void test_whittaker_seriescoefficient() {
    assert_eval_eq("SeriesCoefficient[WhittakerM[2, 1/2, x], {x, 0, 1}]", "1", 0);
    assert_eval_eq("SeriesCoefficient[WhittakerM[2, 1/2, x], {x, 0, 2}]", "-1", 0);
    assert_eval_eq("SeriesCoefficient[WhittakerM[2, 1/2, x], {x, 0, 3}]", "3/8", 0);
}

/* ---- argument count ------------------------------------------------- */

void test_whittaker_argcount() {
    assert_eval_eq("WhittakerM[]", "WhittakerM[]", 0);
    assert_eval_eq("WhittakerM[1, 2]", "WhittakerM[1, 2]", 0);
    assert_eval_eq("WhittakerM[1, 2, 3, 4]", "WhittakerM[1, 2, 3, 4]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_whittaker_attributes() {
    SymbolDef* d = symtab_get_def("WhittakerM");
    ASSERT_MSG(d != NULL, "WhittakerM not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "WhittakerM not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "WhittakerM not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "WhittakerM not Protected");
}

/* ===================== WhittakerW =================================== */
/* W_{k,m}(z) = e^(-z/2) z^(m+1/2) HypergeometricU[m-k+1/2, 2m+1, z].
 * Numeric reference values from mpmath (mp.whitw), an independent source. */

void test_whittaker_w_special() {
    /* Values at z == 0 by the leading exponent 1/2 - |Re m|. */
    assert_eval_eq("WhittakerW[0, 0, 0]", "0", 0);              /* |Re m| < 1/2 */
    assert_eval_eq("WhittakerW[2, 0.5, 0]", "0", 0);           /* m=1/2 -> 1/Gamma[-1] = 0 */
    assert_eval_eq("WhittakerW[0, 1, 0]", "ComplexInfinity", 0);/* |Re m| > 1/2 */
    assert_eval_eq("WhittakerW[k, 1/2, 0]", "1/Gamma[1 - k]", 0);
    /* Symbolic / all-exact arguments stay symbolic. */
    assert_eval_eq("WhittakerW[2, 3, x]", "WhittakerW[2, 3, x]", 0);
    assert_eval_eq("WhittakerW[2, 3, 4]", "WhittakerW[2, 3, 4]", 0);
}

void test_whittaker_w_machine_real() {
    assert_close("WhittakerW[2, 3, 4.5]", 11.306294801893110, 1e-6);
    assert_close("WhittakerW[6, 4, 1.7]", 1740.546241809134, 1e-3);
    /* non-integer b = 2m+1 (Kummer connection regime), machine precision. */
    assert_close("WhittakerW[0.3, 0.4, 0.55]", 0.717440169793, 1e-9);
}

void test_whittaker_w_arbitrary() {
    /* N[W[3/5,1/3,1/7],30]; mpmath to 30 digits, first ~28 rounding-stable. */
    assert_eval_startswith("N[WhittakerW[3/5, 1/3, 1/7], 30]",
        "0.3970969447006262740475282170");
    /* A 19-digit literal raises the working precision above machine; the output
     * tracks it (non-integer b, mpmath-confirmed). */
    assert_eval_startswith("WhittakerW[2, 1.2222222222222222222, 1]",
        "-0.242289411116497");
}

void test_whittaker_w_complex() {
    /* N[W[1/5, 2 - I, 2]] = 0.314759 - 1.27345 I. */
    assert_close("Re[N[WhittakerW[1/5, 2 - I, 2]]]", 0.3147585139, 1e-5);
    assert_close("Im[N[WhittakerW[1/5, 2 - I, 2]]]", -1.273446047, 1e-5);
    /* Complex k and MPFR z (45-digit), mpmath-confirmed. */
    assert_eval_startswith("Re[WhittakerW[2 - 5 I, 3, 4.5`45]]",
        "-283.2166670975702708743866");
    assert_eval_startswith("Im[WhittakerW[2 - 5 I, 3, 4.5`45]]",
        "-399.7625181886174290903785");
}

void test_whittaker_w_listable() {
    assert_eval_eq("WhittakerW[k, m, {a, b}]",
                   "{WhittakerW[k, m, a], WhittakerW[k, m, b]}", 0);
    assert_eval_eq("WhittakerW[k, m, {}]", "{}", 0);
    assert_close("WhittakerW[2, 3, {1.5, 2.5, 3.5}][[1]]", 86.38286431, 1e-3);
    assert_close("WhittakerW[2, 3, {1.5, 2.5, 3.5}][[3]]", 17.99204277, 1e-4);
}

void test_whittaker_w_ndarray() {
    assert_close("WhittakerW[2, 3, NDArray[{1.5, 2.5, 3.5}]][[2]]", 32.87672353, 1e-4);
    /* a large numeric list auto-packs; the packed fast path shares the machine
     * helper, so an element matches the scalar path bit-closely. */
    assert_close("WhittakerW[2, 3, N[Range[60]/13]][[7]]",
                 eval_real("WhittakerW[2, 3, N[7/13]]"), 1e-12);
    assert_close("WhittakerW[2, 3, N[Range[60]/13]][[40]]",
                 eval_real("WhittakerW[2, 3, N[40/13]]"), 1e-12);
}

void test_whittaker_w_compile() {
    assert_eval_eq(
        "CompileDiagnostics[{{x, _Real}, {y, _Real}, {z, _Real}}, "
        "WhittakerW[x, y, z]][[1]]", "\"Compiled\" -> True", 0);
    assert_close(
        "Compile[{{x, _Real}, {y, _Real}, {z, _Real}}, WhittakerW[x, y, z]][2., 3., 4.5]",
        11.306294801893110, 1e-6);
}

/* Regression: W_{k,m} with m - k + 1/2 = -n is elementary (U is the polynomial
 * U(-n, 2m+1, z)), and returned 0.0 at machine precision because the U kernel
 * did. W_{2,1/2}(z) = (z - 2) z e^(-z/2) ... at z = 1 that is -e^(-1/2).
 * References: mpmath 1.4 whitw. */
void test_whittaker_w_terminating() {
    assert_close("WhittakerW[2, 1/2, 1.]", -0.606530659712633, 1e-14);
    assert_close("WhittakerW[2, 1/2, 1.5]", -0.354274914555761, 1e-14);
    assert_close("WhittakerW[3, 1/2, 1.]", 0.606530659712633, 1e-14);
    assert_close("WhittakerW[3, 1/2, 1.5]", -0.531412371833642, 1e-14);
    assert_close("N[WhittakerW[2, 1/2, 1]]", -0.606530659712633, 1e-14);
    assert_eval_eq(
        "With[{zs = 0.37 N[Range[12]]}, "
        "Table[WhittakerW[2, 1/2, z], {z, zs}] === WhittakerW[2, 1/2, zs] === "
        "Normal[WhittakerW[2, 1/2, NDArray[zs]]] === "
        "Compile[{{v, _Real, 1}}, WhittakerW[2, 1/2, v]][zs]]", "True", 0);
    assert_eval_eq(
        "CompileDiagnostics[{{v, _Real, 1}}, WhittakerM[2, 1/2, v]][[1]]",
        "\"Compiled\" -> True", 0);
}

void test_whittaker_w_derivative() {
    assert_eval_eq("D[WhittakerW[k, m, z], z]",
        "(1/2 - k/z) WhittakerW[k, m, z] - WhittakerW[1 + k, m, z]/z", 0);
    assert_eval_eq("D[WhittakerW[k, m, z], k]",
                   "Derivative[1, 0, 0][WhittakerW][k, m, z]", 0);
    assert_eval_eq("D[WhittakerW[k, m, z], m]",
                   "Derivative[0, 1, 0][WhittakerW][k, m, z]", 0);
    assert_eval_eq("D[WhittakerW[k, m, y], z]", "0", 0);
}

void test_whittaker_w_series() {
    /* Terminating case (a = m-k+1/2 = -1): origin Puiseux collapses to a Taylor. */
    assert_eval_eq("Series[WhittakerW[2, 1/2, x], {x, 0, 6}]",
        "-2 x + 2 x^2 - 3/4 x^3 + 1/6 x^4 - 5/192 x^5 + 1/320 x^6 + O[x]^7", 0);
    /* Terminating case at Infinity -> the exact elementary closed form. */
    assert_eval_eq("Series[WhittakerW[2, 1/2, x], {x, Infinity, 5}]",
                   "-((2 - x) x)/E^(1/2 x)", 0);
    /* Non-terminating / symbolic stays symbolic (deferred). */
    assert_eval_eq("Series[WhittakerW[1/3, 1/2, x], {x, Infinity, 3}]",
                   "Series[WhittakerW[1/3, 1/2, x], {x, Infinity, 3}]", 0);
}

void test_whittaker_w_argcount() {
    assert_eval_eq("WhittakerW[]", "WhittakerW[]", 0);
    assert_eval_eq("WhittakerW[1, 2]", "WhittakerW[1, 2]", 0);
    assert_eval_eq("WhittakerW[1, 2, 3, 4]", "WhittakerW[1, 2, 3, 4]", 0);
}

void test_whittaker_w_attributes() {
    SymbolDef* d = symtab_get_def("WhittakerW");
    ASSERT_MSG(d != NULL, "WhittakerW not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "WhittakerW not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "WhittakerW not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "WhittakerW not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_whittaker_special);
    TEST(test_whittaker_machine_real);
    TEST(test_whittaker_arbitrary);
    TEST(test_whittaker_complex);
    TEST(test_whittaker_listable);
    TEST(test_whittaker_ndarray);
    TEST(test_whittaker_compile);
    TEST(test_whittaker_derivative);
    TEST(test_whittaker_series_zero);
    TEST(test_whittaker_series_infinity);
    TEST(test_whittaker_seriescoefficient);
    TEST(test_whittaker_argcount);
    TEST(test_whittaker_attributes);

    TEST(test_whittaker_w_special);
    TEST(test_whittaker_w_machine_real);
    TEST(test_whittaker_w_arbitrary);
    TEST(test_whittaker_w_complex);
    TEST(test_whittaker_w_listable);
    TEST(test_whittaker_w_ndarray);
    TEST(test_whittaker_w_compile);
    TEST(test_whittaker_w_terminating);
    TEST(test_whittaker_w_derivative);
    TEST(test_whittaker_w_series);
    TEST(test_whittaker_w_argcount);
    TEST(test_whittaker_w_attributes);

    printf("All WhittakerM / WhittakerW tests passed.\n");
    return 0;
}
