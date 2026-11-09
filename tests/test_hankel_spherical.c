/* Tests for HankelH1, HankelH2, SphericalBesselJ, SphericalBesselY,
 * SphericalHankelH1, SphericalHankelH2.
 *
 * Reference values: mpmath (hankel1, hankel2, besselj, bessely) at 50 digits,
 * with the reference system's spherical prefactor Sqrt[Pi]/Sqrt[2 z] (principal
 * Sqrt[2 z]) and parity for an integer order on the negative axis.
 *
 * Covers: z = 0 values, unevaluated exact calls; 60 mpmath points (every
 * function: real positive / negative z, integer and fractional order, complex z
 * in both half planes -- where J + I Y would cancel -- and large order);
 * arbitrary precision; identities (Wronskians, H1 + H2 = 2 J, the elementary
 * j_0, y_0); D; packed / NDArray buffers and Compile[]; messages; attributes.
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

/* {code (HS_* order), n, Re z, Im z, Re value, Im value} */
static const double HS_REF[][6] = {
    { 0, 2.5, 3.0, 0.0, 4.1271003220971599344e-1, -3.6904073007379789735e-1 },
    { 0, 0.0, 0.3, 0.0, 9.7762624653829608757e-1, -8.0727357780451946575e-1 },
    { 0, 7.0, 12.0, 0.0, -1.702538041272080471e-1, 1.8952069552168660076e-1 },
    { 0, 30.0, 5.0, 0.0, 2.6711772782507988106e-21, -4.0285684185540875716e18 },
    { 0, 1.5, -4.0, 0.0, -3.6711203246093415489e-1, -1.8528594835426895264e-1 },
    { 0, 3.0, -2.0, 0.0, 1.289432494744020511e-1, 1.1277837768404277861 },
    { 0, 0.3, 2.0, 5.0, 2.2635110337252285675e-3, 3.269889034197171841e-4 },
    { 0, 0.3, 2.0, -5.0, 1.2266343551032776778e1, 5.031965102900921971e1 },
    { 0, 4.0, -1.5, 0.5, 5.2753645006129399896, -2.361118067241569721 },
    { 0, 12.0, 80.0, 0.0, 3.6310262577614403695e-4, -8.9713411475517826835e-2 },
    { 1, 2.5, 3.0, 0.0, 4.1271003220971599344e-1, 3.6904073007379789735e-1 },
    { 1, 0.0, 0.3, 0.0, 9.7762624653829608757e-1, 8.0727357780451946575e-1 },
    { 1, 7.0, 12.0, 0.0, -1.702538041272080471e-1, -1.8952069552168660076e-1 },
    { 1, 30.0, 5.0, 0.0, 2.6711772782507988106e-21, 4.0285684185540875716e18 },
    { 1, 1.5, -4.0, 0.0, 3.6711203246093415489e-1, -1.8528594835426895264e-1 },
    { 1, 3.0, -2.0, 0.0, -3.868297484232061533e-1, -1.1277837768404277861 },
    { 1, 0.3, 2.0, 5.0, 1.2266343551032776778e1, -5.031965102900921971e1 },
    { 1, 0.3, 2.0, -5.0, 2.2635110337252285675e-3, -3.269889034197171841e-4 },
    { 1, 4.0, -1.5, 0.5, -5.2649901260851071238, 2.3335623298667686914 },
    { 1, 12.0, 80.0, 0.0, 3.6310262577614403695e-4, 8.9713411475517826835e-2 },
    { 2, 2.5, 3.0, 0.0, 2.2363817356109351996e-1, 0.0 },
    { 2, 0.0, 0.3, 0.0, 9.8506735553779858368e-1, 0.0 },
    { 2, 7.0, 12.0, 0.0, -2.4838740722100813262e-2, 0.0 },
    { 2, 30.0, 5.0, 0.0, 4.2827302172992125126e-22, 0.0 },
    { 2, 1.5, -4.0, 0.0, 0.0, -2.2818347649544163888e-1 },
    { 2, 3.0, -2.0, 0.0, -6.0722097662874828461e-2, 0.0 },
    { 2, 0.3, 2.0, 5.0, 5.2489998927987535891, -1.2209385178121296948e1 },
    { 2, 0.3, 2.0, -5.0, 5.2489998927987535891, 1.2209385178121296948e1 },
    { 2, 4.0, -1.5, 0.5, 2.0869218710664477998e-3, -5.6652323627594836809e-3 },
    { 2, 12.0, 80.0, 0.0, -8.1456778050019470477e-3, 0.0 },
    { 3, 2.5, 3.0, 0.0, -3.8968938904543784199e-1, 0.0 },
    { 3, 0.0, 0.3, 0.0, -3.1844549637520200655, 0.0 },
    { 3, 7.0, 12.0, 0.0, 9.0579445655870130622e-2, 0.0 },
    { 3, 30.0, 5.0, 0.0, -7.7607175697584787819e18, 0.0 },
    { 3, 1.5, -4.0, 0.0, 4.5636695299088327776e-1, -1.3529751370697261114e-1 },
    { 3, 3.0, -2.0, 0.0, -1.4843665574430799239, 0.0 },
    { 3, 0.3, 2.0, 5.0, 1.2208172330142952062e1, 5.2485571448726599003 },
    { 3, 0.3, 2.0, -5.0, 1.2208172330142952062e1, -5.2485571448726599003 },
    { 3, 4.0, -1.5, 0.5, 9.389005157654214987e-1, 1.2236154214412702734e1 },
    { 3, 12.0, 80.0, 0.0, -9.5831769927638470494e-3, 0.0 },
    { 4, 2.5, 3.0, 0.0, 2.2363817356109351996e-1, -3.8968938904543784199e-1 },
    { 4, 0.0, 0.3, 0.0, 9.8506735553779858368e-1, -3.1844549637520200655 },
    { 4, 7.0, 12.0, 0.0, -2.4838740722100813262e-2, 9.0579445655870130622e-2 },
    { 4, 30.0, 5.0, 0.0, 4.2827302172992125126e-22, -7.7607175697584787819e18 },
    { 4, 1.5, -4.0, 0.0, 1.3529751370697261114e-1, 2.2818347649544163888e-1 },
    { 4, 3.0, -2.0, 0.0, -6.0722097662874828461e-2, -1.4843665574430799239 },
    { 4, 0.3, 2.0, 5.0, 4.427479260936887451e-4, -1.2128479783448852906e-3 },
    { 4, 0.3, 2.0, -5.0, 1.0497557037671413489e1, 2.441755750826424901e1 },
    { 4, 4.0, -1.5, 0.5, -1.2234067292541636286e1, 9.3323528340266201502e-1 },
    { 4, 12.0, 80.0, 0.0, -8.1456778050019470477e-3, -9.5831769927638470494e-3 },
    { 5, 2.5, 3.0, 0.0, 2.2363817356109351996e-1, 3.8968938904543784199e-1 },
    { 5, 0.0, 0.3, 0.0, 9.8506735553779858368e-1, 3.1844549637520200655 },
    { 5, 7.0, 12.0, 0.0, -2.4838740722100813262e-2, -9.0579445655870130622e-2 },
    { 5, 30.0, 5.0, 0.0, 4.2827302172992125126e-22, 7.7607175697584787819e18 },
    { 5, 1.5, -4.0, 0.0, -1.3529751370697261114e-1, -6.8455042948632491665e-1 },
    { 5, 3.0, -2.0, 0.0, -6.0722097662874828461e-2, 1.4843665574430799239 },
    { 5, 0.3, 2.0, 5.0, 1.0497557037671413489e1, -2.441755750826424901e1 },
    { 5, 0.3, 2.0, -5.0, 4.427479260936887451e-4, 1.2128479783448852906e-3 },
    { 5, 4.0, -1.5, 0.5, 1.2238241136283769181e1, -9.4456574812818098238e-1 },
    { 5, 12.0, 80.0, 0.0, -8.1456778050019470477e-3, 9.5831769927638470494e-3 },
};
static const char* const HS_MP[2] = { "0.412710032209715993437496795941862714987261116`45", "-0.369040730073797897345226869879039177007875281`45" };
static const char* const HS_MPJ = "0.000000000224839496827855160966814693354793173200932058`45";

static const char* const NAMES[6] = {
    "HankelH1", "HankelH2", "SphericalBesselJ", "SphericalBesselY",
    "SphericalHankelH1", "SphericalHankelH2" };

static Expr* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static void eval_num(const char* input, double* re, double* im) {
    Expr* r = eval_str(input);
    if (r->type == EXPR_REAL) { *re = r->data.real; *im = 0.0; expr_free(r); return; }
    ASSERT_MSG(r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
               r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
               r->data.function.args[0]->type == EXPR_REAL &&
               r->data.function.args[1]->type == EXPR_REAL,
               "%s: expected a machine number", input);
    *re = r->data.function.args[0]->data.real;
    *im = r->data.function.args[1]->data.real;
    expr_free(r);
}

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

static const char* R(double v, char* buf) {
    char t[40];
    snprintf(t, sizeof t, "%.17g", v);
    char* e = strchr(t, 'e');
    if (e) { *e = 0; snprintf(buf, 64, "%s%s*^%s", t, strchr(t, '.') ? "" : ".", e + 1); }
    else snprintf(buf, 64, "%s%s", t, strchr(t, '.') ? "" : ".");
    return buf;
}

/* 128 ulp: the shared machine Bessel kernels carry up to ~60 ulp near zeros */
#define HS_TOL 2.9e-14

void test_hs_exact() {
    assert_eval_eq("HankelH1[0, 0]", "ComplexInfinity", 0);
    assert_eval_eq("HankelH2[1, 0]", "ComplexInfinity", 0);
    assert_eval_eq("SphericalBesselJ[0, 0]", "1", 0);
    assert_eval_eq("SphericalBesselJ[3, 0]", "0", 0);
    assert_eval_eq("SphericalBesselY[0, 0]", "ComplexInfinity", 0);
    assert_eval_eq("SphericalBesselJ[n, 0]", "SphericalBesselJ[n, 0]", 0);
    assert_eval_eq("HankelH1[1/2, z]", "HankelH1[1/2, z]", 0);
    assert_eval_eq("SphericalBesselJ[2, z]", "SphericalBesselJ[2, z]", 0);
    assert_true("NumericQ[SphericalBesselJ[2, 3]]");
    assert_true("N[SphericalBesselJ[0, 2]] == Sin[2.]/2.");
}

void test_hs_machine() {
    char buf[400], b0[64], b1[64], b2[64];
    for (size_t i = 0; i < sizeof(HS_REF) / sizeof(HS_REF[0]); i++) {
        const double* c = HS_REF[i];
        if (c[3] == 0.0)
            snprintf(buf, sizeof buf, "%s[%s, %s]", NAMES[(int)c[0]], R(c[1], b0), R(c[2], b1));
        else
            snprintf(buf, sizeof buf, "%s[%s, Complex[%s, %s]]", NAMES[(int)c[0]],
                     R(c[1], b0), R(c[2], b1), R(c[3], b2));
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(hypot(re - c[4], im - c[5]) <= HS_TOL * hypot(c[4], c[5]),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", buf, c[4], c[5], re, im);
    }
    /* real j, y for real input; complex Hankel values */
    assert_true("MachineNumberQ[SphericalBesselJ[3, 2.5]] && Head[HankelH1[1.5, 2.]] === Complex");
    assert_true("Head[SphericalBesselJ[2.5, -2.]] === Complex");
    /* a complex order: arb inside the kernel */
    assert_true("Abs[HankelH1[1.5 + 0.5 I, 2.] - (BesselJ[1.5 + 0.5 I, 2.] + I BesselY[1.5 + 0.5 I, 2.])] < 10^-13");
}

void test_hs_arbitrary() {
    char buf[400];
    snprintf(buf, sizeof buf, "With[{v = N[HankelH1[5/2, 3], 40], r = %s + %s I}, Abs[v - r] < 10^-38 Abs[r] && Precision[v] >= 39]", HS_MP[0], HS_MP[1]);
    assert_true(buf);
    snprintf(buf, sizeof buf, "With[{v = N[SphericalBesselJ[7, 1/3], 40], r = %s}, Abs[v - r] < 10^-38 Abs[r]]", HS_MPJ);
    assert_true(buf);
    /* the K route at 40 digits: H1 decays far up the imaginary axis */
    assert_true("Abs[N[HankelH1[2, 30 I], 40] - N[2/(Pi I^3) BesselK[2, 30], 45]] < 10^-50");
}

void test_hs_identities() {
    assert_true("Abs[HankelH1[1.3, 2.7] + HankelH2[1.3, 2.7] - 2 BesselJ[1.3, 2.7]] < 10^-15");
    /* Wronskian J_nu H1'_nu - J'_nu H1_nu = 2 I/(Pi z) */
    assert_true("With[{nu = 1/3`40, z = 7/5`40}, Abs[BesselJ[nu, z] (HankelH1[nu - 1, z] - HankelH1[nu + 1, z])/2 - (BesselJ[nu - 1, z] - BesselJ[nu + 1, z])/2 HankelH1[nu, z] - 2 I/(Pi z)] < 10^-37]");
    /* j_n y_{n-1} - j_{n-1} y_n = 1/z^2 */
    assert_true("With[{z = 3.7}, Abs[SphericalBesselJ[4, z] SphericalBesselY[3, z] - SphericalBesselJ[3, z] SphericalBesselY[4, z] - 1/z^2] < 10^-15]");
    assert_true("Abs[SphericalBesselY[0, 2.2] + Cos[2.2]/2.2] < 10^-16");
    assert_true("Abs[SphericalHankelH1[0, 1.7] - (-I) Exp[1.7 I]/1.7] < 10^-15");
    /* the recurrence f_{n-1} + f_{n+1} = (2n + 1)/z f_n, complex z */
    assert_true("With[{z = 2. + 1.5 I}, Abs[SphericalBesselJ[2, z] + SphericalBesselJ[4, z] - 7/z SphericalBesselJ[3, z]] < 10^-14]");
}

void test_hs_derivative() {
    assert_true("D[HankelH1[n, z], z] === (HankelH1[-1 + n, z] - HankelH1[1 + n, z])/2");
    assert_true("D[SphericalBesselJ[n, z], z] === -SphericalBesselJ[n, z]/(2 z) + (SphericalBesselJ[-1 + n, z] - SphericalBesselJ[1 + n, z])/2");
    assert_true("D[HankelH2[n, z], n] === Derivative[1, 0][HankelH2][n, z]");
    assert_true("Abs[(D[SphericalBesselY[3, z], z] /. z -> 2.3) - (SphericalBesselY[3, 2.3 + 10.^-6] - SphericalBesselY[3, 2.3 - 10.^-6])/(2 10.^-6)] < 10^-7");
}

void test_hs_buffers() {
    for (int k = 0; k < 6; k++) {
        char buf[400];
        snprintf(buf, sizeof buf,
                 "With[{v = Range[0.05, 40., 0.05]}, PackedArrayQ[%s[3, v]] && %s[3, v] === Table[%s[3, x], {x, FromPackedArray[v]}]]",
                 NAMES[k], NAMES[k], NAMES[k]);
        assert_true(buf);
    }
    assert_true("With[{v = Range[-10., 10., 0.5] + 0.3 I}, Normal[HankelH1[1.5, NDArray[v]]] === HankelH1[1.5, FromPackedArray[v]]]");
    assert_true("With[{v = Range[-5., -0.5, 0.5]}, SphericalBesselJ[2.5, v] === Table[SphericalBesselJ[2.5, x], {x, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[0., 20., 1.]}, SphericalBesselJ[v, 3.3] === Table[SphericalBesselJ[n, 3.3], {n, FromPackedArray[v]}]]");
}

void test_hs_compile() {
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, SphericalBesselJ[3, x]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, SphericalBesselY[2, v]][[1]]", "\"Compiled\" -> True", 0);
    assert_true("With[{f = Compile[{{v, _Real, 1}}, SphericalBesselJ[5, v]], v = Range[0.1, 30., 0.01]}, f[v] === SphericalBesselJ[5, v]]");
}

void test_hs_messages_attributes() {
    test_mute_messages();
    assert_eval_eq("Check[HankelH1[1.], \"failed\"]", "\"failed\"", 0);
    for (int k = 0; k < 6; k++) {
        char buf[200];
        snprintf(buf, sizeof buf, "Attributes[%s]", NAMES[k]);
        assert_eval_eq(buf, "{Listable, NumericFunction, Protected}", 0);
        SymbolDef* d = symtab_get_def(NAMES[k]);
        ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, NAMES[k]));
    }
}

int main() {
    symtab_init();
    core_init();
    TEST(test_hs_exact);
    TEST(test_hs_machine);
    TEST(test_hs_arbitrary);
    TEST(test_hs_identities);
    TEST(test_hs_derivative);
    TEST(test_hs_buffers);
    TEST(test_hs_compile);
    TEST(test_hs_messages_attributes);
    printf("All Hankel / spherical Bessel tests passed.\n");
    return 0;
}
