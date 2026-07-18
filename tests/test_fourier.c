/* Tests for Fourier / InverseFourier: the discrete Fourier transform.
 *
 * Covers machine-precision 1-D and multidimensional transforms (golden vectors
 * from the Wolfram Language reference), the roundoff real-collapse, the
 * Fourier <-> InverseFourier round trip, the FourierParameters {a,b} option
 * conventions, arbitrary-precision (MPFR radix-2 + Bluestein) transforms, the
 * exact symbolic transform, the position-selection form, argument diagnostics,
 * and attributes. */

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

/* Evaluate a scalar expression to a double (Integer / Real / MPFR). */
static double num_of(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    double v;
    switch (r->type) {
        case EXPR_REAL:    v = r->data.real;            break;
        case EXPR_INTEGER: v = (double)r->data.integer; break;
#ifdef USE_MPFR
        case EXPR_MPFR:    v = mpfr_get_d(r->data.mpfr, MPFR_RNDN); break;
#endif
        default:
            ASSERT_MSG(0, "%s: expected a real number, got head/type %d",
                       input, (int)r->type);
            v = NAN;
    }
    expr_free(r);
    return v;
}

/* Check that the scalar expression `elem` equals exp_re + exp_im I (tol). */
static void chk(const char* elem, double exp_re, double exp_im, double tol) {
    char buf[512];
    snprintf(buf, sizeof(buf), "Re[%s]", elem);
    double re = num_of(buf);
    snprintf(buf, sizeof(buf), "Im[%s]", elem);
    double im = num_of(buf);
    ASSERT_MSG(fabs(re - exp_re) <= tol,
               "%s: Re expected %.10g, got %.10g", elem, exp_re, re);
    ASSERT_MSG(fabs(im - exp_im) <= tol,
               "%s: Im expected %.10g, got %.10g", elem, exp_im, im);
}

/* ---- machine-precision 1-D (WMA golden vectors) --------------------- */

void test_machine_1d() {
    /* Fourier[{1,1,2,2,1,1,0,0}] */
    chk("Fourier[{1,1,2,2,1,1,0,0}][[1]]",  2.8284271247, 0.0,          1e-6);
    chk("Fourier[{1,1,2,2,1,1,0,0}][[2]]", -0.5,          1.2071067812, 1e-6);
    chk("Fourier[{1,1,2,2,1,1,0,0}][[4]]",  0.5,         -0.2071067812, 1e-6);
    chk("Fourier[{1,1,2,2,1,1,0,0}][[8]]", -0.5,         -1.2071067812, 1e-6);

    /* Fourier[{1,2,3,4,5,6}] — conjugate symmetry. */
    chk("Fourier[{1,2,3,4,5,6}][[1]]", 8.5732141, 0.0, 1e-5);
    chk("Fourier[{1,2,3,4,5,6}][[2]]", -1.2247449, -2.1213203, 1e-5);
    chk("Fourier[{1,2,3,4,5,6}][[6]]", -1.2247449,  2.1213203, 1e-5);
}

void test_power_spectrum() {
    /* Abs[Fourier[{1,2,3,4,5,6}]]^2 = {73.5, 6, 2, 1.5, 2, 6}. */
    ASSERT_MSG(fabs(num_of("(Abs[Fourier[{1,2,3,4,5,6}]]^2)[[1]]") - 73.5) < 1e-4, "power[1]");
    ASSERT_MSG(fabs(num_of("(Abs[Fourier[{1,2,3,4,5,6}]]^2)[[4]]") - 1.5) < 1e-4, "power[4]");
    ASSERT_MSG(fabs(num_of("(Abs[Fourier[{1,2,3,4,5,6}]]^2)[[3]]") - 2.0) < 1e-4, "power[3]");
}

/* ---- roundoff real collapse ----------------------------------------- */

void test_real_collapse() {
    /* Fourier[{1,0,1,0,1,0}] is numerically real -> a real list. */
    assert_eval_eq("Head[Fourier[{1,0,1,0,1,0}]]", "List", 0);
    chk("Fourier[{1,0,1,0,1,0}][[1]]", 1.2247449, 0.0, 1e-6);
    chk("Fourier[{1,0,1,0,1,0}][[4]]", 1.2247449, 0.0, 1e-6);
    chk("Fourier[{1,0,1,0,1,0}][[2]]", 0.0, 0.0, 1e-6);
    /* A real coefficient prints as a bare Real (imaginary part collapsed). */
    assert_eval_eq("Head[Fourier[{1,0,1,0,1,0}][[1]]]", "Real", 0);
}

/* ---- inverse round trip --------------------------------------------- */

void test_inverse_roundtrip() {
    chk("InverseFourier[Fourier[{3,1,4,1,5,9}]][[1]]", 3.0, 0.0, 1e-6);
    chk("InverseFourier[Fourier[{3,1,4,1,5,9}]][[3]]", 4.0, 0.0, 1e-6);
    chk("InverseFourier[Fourier[{3,1,4,1,5,9}]][[6]]", 9.0, 0.0, 1e-6);
    /* Non-power-of-two length exercises a mixed-radix FFTW plan. */
    chk("InverseFourier[Fourier[{2,7,1,8,2,8,1,8,2,8,4}]][[5]]", 2.0, 0.0, 1e-6);
}

/* ---- multidimensional ----------------------------------------------- */

void test_multidim() {
    /* 2-D DC term = sum / Sqrt[6]. */
    chk("Fourier[{{1,2,3},{4,5,6}}][[1,1]]", 8.5732141, 0.0, 1e-5);
    /* 2-D round trip. */
    chk("InverseFourier[Fourier[{{1,2,3},{4,5,6}}]][[2,3]]", 6.0, 0.0, 1e-6);
    /* 3-D doc example: nonzero diagonal, DC-ish corner = 2/Sqrt[24]. */
    ASSERT_MSG(fabs(num_of(
        "Module[{x=ConstantArray[0,{2,3,4}]}, x[[1,1,1]]=1; x[[2,2,2]]=1; "
        "Re[Fourier[x][[1,1,1]]]]") - 0.40824829) < 1e-6, "3D corner");
}

/* ---- FourierParameters conventions ---------------------------------- */

void test_fourier_parameters() {
    /* {1,1}: no 1/Sqrt[n] normalisation -> DC is the raw sum = 4. */
    chk("Fourier[{1,0,1,0,0,1,0,0,0,1},FourierParameters->{1,1}][[1]]", 4.0, 0.0, 1e-5);
    chk("Fourier[{1,0,1,0,0,1,0,0,0,1},FourierParameters->{1,1}][[2]]",
        1.1180340, 0.3632712, 1e-5);
    /* {-1,1}: 1/n normalisation -> DC = 4/10 = 0.4. */
    chk("Fourier[{1,0,1,0,0,1,0,0,0,1},FourierParameters->{-1,1}][[1]]", 0.4, 0.0, 1e-5);
    /* b=2 is non-invertible: the transform matrix has rank 4. */
    assert_eval_eq(
        "MatrixRank[Map[Fourier[#,FourierParameters->{0,2}]&, IdentityMatrix[8]]]",
        "4", 0);
}

/* ---- arbitrary precision (MPFR FFT) --------------------------------- */

#ifdef USE_MPFR
void test_arbitrary_precision() {
    /* 24-digit golden values (WMA), radix-2 length 7 -> Bluestein. */
    chk("Fourier[N[{1,0,0,1,0,0,1},24]][[1]]", 1.13389341902768168164355, 0.0, 1e-18);
    chk("Fourier[N[{1,0,0,1,0,0,1},24]][[2]]",
        0.27308724404093304366422, -0.13151188545021087883096, 1e-18);
    /* The DC-symmetric coefficient's imaginary roundoff is chopped: chk above
     * asserts Im == 0 to 1e-18, i.e. the sub-precision noise was removed. */
    /* Power-of-two length exercises the radix-2 kernel directly. */
    chk("Fourier[N[{1,2,3,4},30]][[1]]", 5.0, 0.0, 1e-24);
    /* Non-power-of-two length 5 exercises Bluestein; DC = 15/Sqrt[5]. */
    chk("Fourier[N[{1,2,3,4,5},30]][[1]]", 6.70820393249936908923, 0.0, 1e-18);
    /* Arbitrary-precision round trip. */
    chk("InverseFourier[Fourier[N[{1,2,3,4,5,6,7},25]]][[4]]", 4.0, 0.0, 1e-18);
}
#endif

/* ---- symbolic transform --------------------------------------------- */

void test_symbolic() {
    assert_eval_eq("Fourier[{a,b}]", "{(a + b)/Sqrt[2], (a - b)/Sqrt[2]}", 0);
    assert_eval_eq("Fourier[{a,b,c,d}]",
        "{1/2 (a + b + c + d), 1/2 (a + I b - c - I d), "
        "1/2 (a - b + c - d), 1/2 (a - I b - c + I d)}", 0);
    /* Numeric substitution check (independent of print formatting). */
    chk("N[Fourier[{a,b}][[1]] /. {a->1, b->3}]", 2.8284271, 0.0, 1e-6);
    /* 1/2 (a + I b - c - I d) at {1,2,3,4} = 1/2 (-2 - 2 I) = -1 - I. */
    chk("N[Fourier[{a,b,c,d}][[2]] /. {a->1,b->2,c->3,d->4}]", -1.0, -1.0, 1e-9);
}

/* ---- position-selection form ---------------------------------------- */

void test_position_form() {
    /* Extracts the requested coefficients of the full transform. */
    assert_eval_eq("Length[Fourier[{1,1,2,2,1,1,0,0}, {{1},{2},{3}}]]", "3", 0);
    chk("Fourier[{1,1,2,2,1,1,0,0}, {{1},{2},{3}}][[1]]", 2.8284271, 0.0, 1e-6);
    chk("Fourier[{1,1,2,2,1,1,0,0}, {{1},{2},{3}}][[2]]", -0.5, 1.2071068, 1e-6);
}

/* ---- diagnostics + attributes --------------------------------------- */

void test_argcount() {
    /* Zero args stays unevaluated (an argt message is emitted to stderr). */
    assert_eval_eq("Fourier[]", "Fourier[]", 0);
    assert_eval_eq("InverseFourier[]", "InverseFourier[]", 0);
    /* Non-integer b on numeric data stays unevaluated. */
    assert_eval_eq("Fourier[{1,2,3}, FourierParameters->{0,1/2}]",
                   "Fourier[{1, 2, 3}, FourierParameters -> {0, 1/2}]", 0);
    /* Symbolic (non-numeric) list with a symbol stays exact, not machine. */
    assert_eval_eq("Head[Fourier[{a,b}]]", "List", 0);
}

/* ---- NDArray fast path ---------------------------------------------- */

void test_ndarray_fastpath() {
    /* An NDArray argument returns an NDArray (Head is NDArray, not List). */
    assert_eval_eq("Head[Fourier[NDArray[{1,1,2,2,1,1,0,0}]]]", "NDArray", 0);
    assert_eval_eq("NDArrayQ[Fourier[NDArray[{1.,2.,3.,4.,5.,6.}]]]", "True", 0);
    /* Same values as the List path. */
    chk("Fourier[NDArray[{1,1,2,2,1,1,0,0}]][[2]]", -0.5, 1.2071068, 1e-6);
    chk("Fourier[NDArray[{{1,2,3},{4,5,6}}]][[1,1]]", 8.5732141, 0.0, 1e-5);
    /* Round trip through NDArray. */
    chk("InverseFourier[Fourier[NDArray[{3,1,4,1,5,9}]]][[3]]", 4.0, 0.0, 1e-6);
    /* Real-collapse yields a real-dtype NDArray. */
    assert_eval_eq("DataType[Fourier[NDArray[{1,0,1,0,1,0}]]]", "\"float64\"", 0);
    assert_eval_eq("DataType[Fourier[NDArray[{1.,2.,3.,4.}]]]", "\"complex64\"", 0);
    /* Position form on an NDArray returns the selected coefficients as a list. */
    assert_eval_eq("Length[Fourier[NDArray[{1,1,2,2,1,1,0,0}], {{1},{2}}]]", "2", 0);
}

void test_attributes() {
    SymbolDef* f = symtab_get_def("Fourier");
    SymbolDef* g = symtab_get_def("InverseFourier");
    ASSERT(f != NULL && g != NULL);
    ASSERT_MSG((f->attributes & ATTR_PROTECTED) != 0, "Fourier must be Protected");
    ASSERT_MSG((g->attributes & ATTR_PROTECTED) != 0, "InverseFourier must be Protected");
    /* FourierParameters default option is registered. */
    assert_eval_eq("FourierParameters /. Options[Fourier]", "{0, 1}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_machine_1d);
    TEST(test_power_spectrum);
    TEST(test_real_collapse);
    TEST(test_inverse_roundtrip);
    TEST(test_multidim);
    TEST(test_fourier_parameters);
#ifdef USE_MPFR
    TEST(test_arbitrary_precision);
#endif
    TEST(test_symbolic);
    TEST(test_position_form);
    TEST(test_ndarray_fastpath);
    TEST(test_argcount);
    TEST(test_attributes);

    printf("All Fourier tests passed.\n");
    return 0;
}
