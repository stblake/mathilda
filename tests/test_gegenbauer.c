/* Tests for GegenbauerC[n, lambda, z] and GegenbauerC[n, z].
 *
 * Reference values: mpmath (mp.gegenbauer, mp.chebyt) at 50 digits; the
 * conventions (explicit-sum polynomial for integer n and every lambda, C_n^0 = 0,
 * the value above the cut z < -1, 2 ChebyshevT[n, z]/n) checked against the
 * reference system's documented outputs.
 *
 * Covers: exact rules and reductions (ChebyshevU, LegendreP, polynomial
 * expansion with symbolic lambda, lambda = -1/2 and -1), the two-argument
 * rewrite and its n = 0 message; machine values over every kernel branch
 * (double / double-double recurrence, the 2F1 series, the ODE stepping towards
 * -1, the general 2F1, arb), complex z and lambda, the cut z < -1, arbitrary
 * precision; identities (recurrence, D_z, orthogonality); D; Listable; packed /
 * NDArray buffers and Compile[] at both arities; messages; attributes.
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

static const double GEG_REAL[][4] = {
    { 3.0, 1.5, 0.3, -1.7774999999999999692 },
    { 10.0, 1.5, 0.3, 2.7275282300707032437 },
    { 20.0, 2.5, -0.7, 5.661096816728266444 },
    { 200.0, 0.75, 0.31, 2.418888721813442588e-1 },
    { 1468.0, 2.654701702662302, 0.6806962410453357, 2.7028132307283246616e3 },
    { 5.0, 0.25, 2.5, 2.070648193359375e2 },
    { 7.0, -0.3, -1.7, 4.7238438228298981107e1 },
    { 2.5, 1.5, 0.3, -1.5821725940398244591 },
    { 2.5, 1.5, 3.0, 1.7339404479791919371e2 },
    { 7.3, 1.5, 0.3, -2.0609300337876860513 },
    { 7.3, 1.5, -0.9, -6.0461884950767739242 },
    { 7.3, 1.5, -0.999, 2.6529264472737926016e2 },
    { 11.89, 5.857, -0.875, -1.6206837434359436911e3 },
    { -3.457, 1.889, 0.4, -7.2559467767113986814e-1 },
    { 0.5, 0.3, 0.9, 6.9569916691170834012e-1 },
    { 30.25, 0.6, 0.05, -9.1752475274779864212e-2 },
    { 2.5, -1.5, 0.3, 6.6665149162736388749e-1 },  /* the lambda -> -3/2 limit; mp.gegenbauer AT -3/2 is off from the 17th digit */
    { 3.5, 4.0, -0.5, -1.3125e1 },
};
/* {n, l, x, y, Re, Im}; y = 0 with x < -1 is the reference's value above the cut */
static const double GEG_CPLX[][6] = {
    { 2.5, 1.5, -3.0, 0.0, 3.1844013048926315061e-4, 1.7339404479791919371e2 },
    { 6.902632092380433, 0.9760664124948026, -2.651760304605201, 0.0, -7.2047803717424509076e4, 2.2752843796543976806e4 },
    { 4.0, 1.5, 0.3, -2.0, 6.4978143750000000642e2, 4.0099499999999998579e2 },
    { 2.5, 1.5, 0.3, -2.0, -6.5258622385229243581e1, 3.0236970855275929593e1 },
    { 12.0, 0.7, -1.1, 0.6, 2.4376760848226777943e2, -1.1630469161652278636e4 },
};
static const double GEG_TWO[][3] = {
    { 3.0, 0.4, -6.2933333333333334932e-1 },
    { 12.0, 0.9237249082406276, 8.0770848523714779922e-4 },
    { 2.5, 0.3, -7.9977596863121612496e-1 },
    { 2.5, 3.0, 3.2809754647055805132e1 },
    { 26.78, -0.2465, 3.0803603240443503917e-3 },
};
static const char* const GEG_MP = "-2.06093003378768626626726048828223585320113231`45";
static const char* const GEG_MPC[2] = { "-74.0258808859650086393068224549032672103878704`45", "-57.9699348799910289273358832685383973852036725`45" };

static Expr* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* The numeric value of `input` as (re, im); im = 0 for a machine Real. */
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

/* A double as a machine-Real literal (a "." forced, exponent as *^). */
static const char* R(double v, char* buf) {
    char t[40];
    snprintf(t, sizeof t, "%.17g", v);
    char* e = strchr(t, 'e');
    if (e) {
        *e = 0;
        snprintf(buf, 64, "%s%s*^%s", t, strchr(t, '.') ? "" : ".", e + 1);
    } else {
        snprintf(buf, 64, "%s%s", t, strchr(t, '.') ? "" : ".");
    }
    return buf;
}

/* 64 ulp; the non-integer-degree branches near -1 are looser (see the docs) */
#define GEG_TOL 1.43e-14

void test_geg_exact() {
    assert_true("GegenbauerC[3, a, x] === -2 a (1 + a) x + 4/3 a (1 + a) (2 + a) x^3");
    assert_eval_eq("GegenbauerC[0, a, x]", "1", 0);
    assert_eval_eq("GegenbauerC[2, 0, x]", "0", 0);
    assert_eval_eq("GegenbauerC[n, 0, x]", "0", 0);
    assert_eval_eq("GegenbauerC[0, 0, x]", "0", 0);
    assert_true("GegenbauerC[n, 1, x] === ChebyshevU[n, x]");
    assert_true("GegenbauerC[n, 1/2, x] === LegendreP[n, x]");
    assert_true("GegenbauerC[2, -1/2, x] === 1/2 - x^2/2");
    assert_eval_eq("GegenbauerC[3, -1, x]", "0", 0);
    assert_true("GegenbauerC[5, -3, x] === -6 x");
    assert_eval_eq("GegenbauerC[4, 3/2, 1/3]", "-5/9", 0);
    /* the reference leaves these unevaluated */
    assert_eval_eq("GegenbauerC[-2, a, x]", "GegenbauerC[-2, a, x]", 0);
    assert_eval_eq("GegenbauerC[n, a, 1]", "GegenbauerC[n, a, 1]", 0);
    assert_eval_eq("GegenbauerC[5/2, 3/2, 1/3]", "GegenbauerC[5/2, 3/2, 1/3]", 0);
    assert_true("NumericQ[GegenbauerC[5/2, 3/2, 1/3]]");
    /* GegenbauerC[n, z] = 2 ChebyshevT[n, z]/n */
    assert_true("GegenbauerC[n, x] === 2 ChebyshevT[n, x]/n");
    assert_true("GegenbauerC[3, x] === 2/3 (-3 x + 4 x^3)");
    test_mute_messages();
    assert_eval_eq("GegenbauerC[0, x]", "ComplexInfinity", 0);
}

void test_geg_machine() {
    char buf[512], b0[64], b1[64], b2[64], b3[64];
    for (size_t i = 0; i < sizeof(GEG_REAL) / sizeof(GEG_REAL[0]); i++) {
        snprintf(buf, sizeof buf, "GegenbauerC[%s, %s, %s]",
                 R(GEG_REAL[i][0], b0), R(GEG_REAL[i][1], b1), R(GEG_REAL[i][2], b2));
        double re, im, e = GEG_REAL[i][3];
        eval_num(buf, &re, &im);
        double tol = (GEG_REAL[i][0] != floor(GEG_REAL[i][0])) ? 256 * GEG_TOL : GEG_TOL;
        ASSERT_MSG(im == 0.0 && fabs(re - e) <= tol * fabs(e),
                   "%s: expected %.17g, got %.17g + %.3g I", buf, e, re, im);
    }
    for (size_t i = 0; i < sizeof(GEG_CPLX) / sizeof(GEG_CPLX[0]); i++) {
        const double* c = GEG_CPLX[i];
        if (c[3] == 0.0)
            snprintf(buf, sizeof buf, "GegenbauerC[%s, %s, %s]", R(c[0], b0), R(c[1], b1), R(c[2], b2));
        else
            snprintf(buf, sizeof buf, "GegenbauerC[%s, %s, Complex[%s, %s]]",
                     R(c[0], b0), R(c[1], b1), R(c[2], b2), R(c[3], b3));
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(hypot(re - c[4], im - c[5]) <= 4 * GEG_TOL * hypot(c[4], c[5]),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", buf, c[4], c[5], re, im);
    }
    for (size_t i = 0; i < sizeof(GEG_TWO) / sizeof(GEG_TWO[0]); i++) {
        snprintf(buf, sizeof buf, "GegenbauerC[%s, %s]", R(GEG_TWO[i][0], b0), R(GEG_TWO[i][1], b1));
        double re, im, e = GEG_TWO[i][2];
        eval_num(buf, &re, &im);
        ASSERT_MSG(fabs(re - e) <= 32 * GEG_TOL * fabs(e), "%s: expected %.17g, got %.17g", buf, e, re);
    }
    /* complex lambda */
    assert_true("Abs[GegenbauerC[2.5, 1.5 + I, 0.3 - 2 I] - (-74.02588088596501 - 57.96993487999103 I)] < 10^-12");
    assert_true("GegenbauerC[1.5, 0, 0.3] == 0");
    assert_true("GegenbauerC[2.5, 1.5, 1.] == 7.875");
    assert_true("MachineNumberQ[GegenbauerC[10, 3/2, 0.3]]");
}

void test_geg_arbitrary() {
    char buf[512];
    snprintf(buf, sizeof buf,
             "With[{v = N[GegenbauerC[73/10, 3/2, 3/10], 40], r = %s}, Abs[v - r] < 10^-38 Abs[r] && Precision[v] >= 39]",
             GEG_MP);
    assert_true(buf);
    snprintf(buf, sizeof buf,
             "With[{v = N[GegenbauerC[5/2, 3/2 + I, 3/10 - 2 I], 40], r = %s + %s I}, Abs[v - r] < 10^-37 Abs[r]]",
             GEG_MPC[0], GEG_MPC[1]);
    assert_true(buf);
    /* lambda = -3/2: the limit the plain 2F1 form cannot give */
    assert_true("Abs[N[GegenbauerC[5/2, -3/2, 3/10], 30] - 0.666651491627363887492635739587`30] < 10^-28");
    /* an integer degree at 40 digits agrees with the exact polynomial */
    assert_true("Abs[GegenbauerC[9, 7/3`40, 3/10`40] - N[GegenbauerC[9, 7/3, 3/10], 45]] < 10^-36");
}

void test_geg_identities() {
    /* (n+1) C_{n+1} = 2 (n + l) z C_n - (n + 2l - 1) C_{n-1}, a non-integer n */
    assert_true("With[{n = 37/10`40, l = 6/5`40, z = 2/5`40}, Abs[(n + 1) GegenbauerC[n + 1, l, z] - 2 (n + l) z GegenbauerC[n, l, z] + (n + 2 l - 1) GegenbauerC[n - 1, l, z]] < 10^-35]");
    /* orthogonality of the polynomials with weight (1 - x^2)^(l - 1/2) */
    assert_true("Integrate[GegenbauerC[2, 3/2, x] GegenbauerC[4, 3/2, x] (1 - x^2), {x, -1, 1}] === 0");
    /* the ODE (1 - z^2) y'' - (2l + 1) z y' + n (n + 2l) y = 0 */
    assert_true("With[{y = GegenbauerC[37/10, 6/5, z]}, Abs[((1 - z^2) D[y, {z, 2}] - (12/5 + 1) z D[y, z] + 37/10 (37/10 + 12/5) y) /. z -> 2/5`30] < 10^-25]");
    /* parity of the polynomials */
    assert_true("Abs[GegenbauerC[7, 1.3, -0.4] + GegenbauerC[7, 1.3, 0.4]] < 10^-14");
}

void test_geg_derivative() {
    assert_true("D[GegenbauerC[n, a, x], x] === 2 a GegenbauerC[-1 + n, 1 + a, x]");
    assert_true("D[GegenbauerC[n, x], x] === 2 ChebyshevU[-1 + n, x]");
    assert_true("D[GegenbauerC[n, a, x], a] === Derivative[0, 1, 0][GegenbauerC][n, a, x]");
    assert_true("D[GegenbauerC[n, a, x], n] === Derivative[1, 0, 0][GegenbauerC][n, a, x]");
}

void test_geg_buffers() {
    assert_true("PackedArrayQ[GegenbauerC[20, 1.5, Range[-1., 1., 0.001]]]");
    assert_true("With[{v = Range[-1., 1., 0.001]}, GegenbauerC[20, 1.5, v] === Table[GegenbauerC[20, 1.5, x], {x, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[-0.999, 1., 0.003]}, GegenbauerC[7.3, 1.5, v] === Table[GegenbauerC[7.3, 1.5, x], {x, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[-3., 3., 0.01]}, GegenbauerC[7.3, v] === Table[GegenbauerC[7.3, x], {x, FromPackedArray[v]}]]");
    /* a real buffer whose values go complex below -1 */
    assert_true("With[{v = Range[-3., 1., 0.25]}, Normal[GegenbauerC[2.5, 1.5, NDArray[v]]] === Table[GegenbauerC[2.5, 1.5, x], {x, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[-2., 2., 0.25] + 0.5 I}, Normal[GegenbauerC[4.2, 0.7, NDArray[v]]] === GegenbauerC[4.2, 0.7, FromPackedArray[v]]]");
    /* a buffer in the lambda slot */
    assert_true("With[{v = Range[0.1, 3., 0.1]}, GegenbauerC[6, v, 0.4] === Table[GegenbauerC[6, l, 0.4], {l, FromPackedArray[v]}]]");
}

void test_geg_compile() {
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, GegenbauerC[3.5, 1.2, x]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, GegenbauerC[3.5, 1.2, v]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, GegenbauerC[3.5, x]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, GegenbauerC[12, v]][[1]]", "\"Compiled\" -> True", 0);
    assert_true("With[{f = Compile[{{v, _Real, 1}}, GegenbauerC[20, 1.5, v]], v = Range[-1., 1., 0.01]}, f[v] === GegenbauerC[20, 1.5, v]]");
    assert_true("With[{f = Compile[{{x, _Real}}, GegenbauerC[7.3, 1.5, x]]}, And @@ Table[f[x] === GegenbauerC[7.3, 1.5, x], {x, -0.99, 1., 0.07}]]");
}

void test_geg_messages_attributes() {
    test_mute_messages();
    assert_eval_eq("Check[GegenbauerC[1.], \"failed\"]", "\"failed\"", 0);
    assert_eval_eq("Attributes[GegenbauerC]", "{Listable, NumericFunction, Protected}", 0);
    SymbolDef* d = symtab_get_def("GegenbauerC");
    ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, "GegenbauerC[n, lambda, z]"));
}

int main() {
    symtab_init();
    core_init();
    TEST(test_geg_exact);
    TEST(test_geg_machine);
    TEST(test_geg_arbitrary);
    TEST(test_geg_identities);
    TEST(test_geg_derivative);
    TEST(test_geg_buffers);
    TEST(test_geg_compile);
    TEST(test_geg_messages_attributes);
    printf("All GegenbauerC tests passed.\n");
    return 0;
}
