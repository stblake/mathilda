/* Tests for the Appell functions AppellF1, AppellF2, AppellF3, AppellF4.
 *
 * Reference values: mpmath appellf1..f4 at 60-100 digits (two precisions
 * agreeing to 30 digits); where mpmath has no analytic continuation (F1 with
 * both variables large, the Euler-integral region) its tanh-sinh quadrature of
 * the Euler integral at 60 digits, with the cuts taken from below. The
 * reference system prints AppellF1[0.3, 0.2, 0.5, 1.7, 10., -20.] as
 * 0.66114752 - 0.09761100 I (and 0.661158271... at 20 digits); the integral is
 * 0.66115898412418623437 - 0.09761940125683091707 I.
 *
 * Covers: exact rules and reductions; machine values real and complex, inside
 * and outside the double-series region, on the cuts, complex parameters, the
 * terminating (polynomial) cases; arbitrary precision real and complex;
 * contiguous relations and the F1 partial differential equation; D; packed,
 * NDArray and Compile[] agreement bit for bit; arbitrary-precision Compile[];
 * attributes and docstrings.
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

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* {call, Re value, Im value} */
static const struct { const char* call; double re, im; } AP_REF[] = {
    { "AppellF1[0.3, 0.2, 0.5, 1.7, 0.4, 0.3]", 1.0474505437456943176, 0.0 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, (0.4+0.1 I), -0.3]", 0.99090911180385493225, 0.0043337333412494106599 },
    { "AppellF1[1.5, -0.7, 2.2, 3.1, -0.8, 0.6]", 3.3185459492170322457, 0.0 },
    { "AppellF1[2.5, 1.5, 0.5, 0.7, 0.9, -0.95]", 3878.7612232527318095, 0.0 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, 0.5, 3.]", 1.0740172138762757162, -0.40392246296843720157 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, -5., 0.3]", 0.9359944671307625191, 0.0 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, 2., -3.]", 0.91896675744421024135, -0.055584516908244188717 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, (0.5+0.1 I), (3.+0.2 I)]", 1.0585213556423869201, 0.38965447732111431923 },
    { "AppellF1[1.2, 0.8, -0.3, 2.6, (-3.+1. I), (0.5-0.7 I)]", 0.50171224271888282603, 0.1138352517372519338 },
    { "AppellF1[0.3, 0.2, 0.5, -0.4, 0.2, 0.3]", 0.74616458277497021017, 0.0 },
    { "AppellF1[-0.582, -1.48, -0.251, 1.19, (0.89+0.0614 I), (0.945-0.233 I)]", 1.7094838944507778956, 0.013497459521894690663 },
    { "AppellF1[(0.3+0.2 I), 0.2, (0.5-0.1 I), 1.7, 0.4, 0.3]", 1.0511510415801153497, 0.026384432335697409512 },
    /* the Euler-integral region (no expansion converges) */
    { "AppellF1[0.3, 0.2, 0.5, 1.7, 10., -20.]", 0.66115898412418623440, -0.097619401256830917067 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, 2., (0.5+2. I)]", 1.0668133908101024148, 0.10227101751446976931 },
    { "AppellF1[0.3, 0.2, 0.5, 1.7, (0.5+2. I), 2.]", 1.2029765586022208764, -0.19715074563106889933 },
    { "AppellF2[0.3, 0.2, 0.5, 1.5, 1.7, 0.2, 0.3]", 1.0397537493086985947, 0.0 },
    { "AppellF2[1.3, 2.2, 0.5, 1.5, 2.7, -0.4, 0.5]", 0.57167149531626013663, 0.0 },
    { "AppellF2[0.3, 0.2, 0.5, 1.5, 1.7, (0.3+0.2 I), -0.25]", 0.99116525812774875366, 0.008887497735718002287 },
    { "AppellF2[0.7, 1.2, -0.5, 2.5, 0.7, -0.6, 0.35]", 0.71765751555009288205, 0.0 },
    { "AppellF3[0.3, 0.2, 0.5, 0.4, 1.7, 0.2, 0.3]", 1.0349149500815721007, 0.0 },
    { "AppellF3[1.3, 2.2, 0.5, 1.5, 2.7, -0.6, 0.85]", 8.4879904057454685599, 0.0 },
    { "AppellF3[0.3, 0.2, 0.5, 0.4, 1.7, (0.5-0.4 I), -0.7]", 1.0144486368653151791, -0.051153993629556914244 },
    { "AppellF3[0.3, 0.2, 0.5, 0.4, 1.7, 0.5, -3.]", 0.96995969534235321053, 0.0 },
    { "AppellF4[0.3, 0.2, 1.5, 1.7, 0.1, 0.2]", 1.0126802468419116382, 0.0 },
    { "AppellF4[1.3, 2.2, 2.5, 0.6, -0.3, 0.15]", 1.1305484644404232388, 0.0 },
    { "AppellF4[0.3, 0.2, 1.5, 1.7, (0.05+0.1 I), -0.2]", 0.9949495883206366256, 0.0034532812511124839904 },
    { "AppellF4[0.3, 0.2, 1.5, 1.7, -2., -3.]", 0.91547885873323880928, 0.0 },
    /* a terminating outer series, the inner variable on its cut */
    { "AppellF1[-3, 0.2, 0.5, 1.7, 0.5, 7.]", -19.395513160219042572, 0.0 },
};

#define AP_TOL 1e-13

static void test_ap_exact(void) {
    assert_true("AppellF1[a, b1, b2, c, x, y] === AppellF1[a, b1, b2, c, x, y]");
    assert_true("AppellF1[a, b1, b2, c, 0, y] === Hypergeometric2F1[a, b2, c, y]");
    assert_true("AppellF1[a, b1, b2, c, x, 0] === Hypergeometric2F1[a, b1, c, x]");
    assert_true("AppellF1[a, 0, b2, c, x, y] === Hypergeometric2F1[a, b2, c, y]");
    assert_true("AppellF1[a, b1, b2, c, x, x] === Hypergeometric2F1[a, b1 + b2, c, x]");
    assert_true("AppellF1[a, b1, b2, a, x, y] === (1 - x)^-b1 (1 - y)^-b2");
    assert_true("AppellF1[a, b1, b2, b1 + b2, x, y] === Hypergeometric2F1[a, b1, b1 + b2, (x - y)/(1 - y)]/(1 - y)^a");
    assert_true("AppellF1[0, b1, b2, c, x, y] === 1");
    assert_true("AppellF1[a, b1, b2, c, 1, y] === Hypergeometric2F1[a, b1, c, 1] Hypergeometric2F1[a, b2, c - b1, y]");
    assert_true("Expand[AppellF1[-2, b1, b2, c, x, y] - (1 - (2 (b1 x + b2 y))/c + (b1 (1 + b1) x^2 + 2 b1 b2 x y + b2 (1 + b2) y^2)/(c (1 + c)))] === 0");
    assert_true("AppellF1[-1, 1, 2, 3, 1/2, 1/3] === 11/18");
    assert_true("AppellF1[1/3, 1/2, 1/5, 2/3, 1/2, 1/3] === AppellF1[1/3, 1/2, 1/5, 2/3, 1/2, 1/3]");
    assert_true("AppellF2[a, b1, b2, c1, c2, 0, y] === Hypergeometric2F1[a, b2, c2, y]");
    assert_true("AppellF2[a, b1, b2, c1, c2, x, 0] === Hypergeometric2F1[a, b1, c1, x]");
    assert_true("AppellF2[0, b1, b2, c1, c2, x, y] === 1");
    assert_true("AppellF3[a1, a2, b1, b2, c, 0, y] === Hypergeometric2F1[a2, b2, c, y]");
    assert_true("AppellF3[a1, a2, 0, b2, c, x, y] === Hypergeometric2F1[a2, b2, c, y]");
    assert_true("AppellF3[a1, a2, b1, b2, c, x, 0] === Hypergeometric2F1[a1, b1, c, x]");
    assert_true("AppellF4[a, b, c1, c2, 0, y] === Hypergeometric2F1[a, b, c2, y]");
    assert_true("AppellF4[a, b, c1, c2, x, 0] === Hypergeometric2F1[a, b, c1, x]");
    assert_true("AppellF4[0, b, c1, c2, x, y] === 1");
    /* a pole in c stays unevaluated */
    assert_true("Head[AppellF1[0.3, 0.2, 0.5, -1, 0.2, 0.3]] === AppellF1");
    assert_true("NumericQ[AppellF1[1/3, 1/2, 1/5, 2/3, 1/2, 1/3]]");
    test_mute_messages();
    assert_eval_eq("Check[AppellF1[1, 2], \"argrx\"]", "\"argrx\"", 0);
}

static void test_ap_machine(void) {
    for (size_t i = 0; i < sizeof(AP_REF) / sizeof(AP_REF[0]); i++) {
        Expr* e = parse_expression(AP_REF[i].call);
        Expr* r = evaluate(e);
        expr_free(e);
        double re = 0, im = 0;
        if (r->type == EXPR_REAL) re = r->data.real;
        else {
            ASSERT_MSG(r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
                       r->data.function.args[0]->type == EXPR_REAL, "%s: not a number", AP_REF[i].call);
            re = r->data.function.args[0]->data.real;
            im = r->data.function.args[1]->data.real;
        }
        expr_free(r);
        ASSERT_MSG(hypot(re - AP_REF[i].re, im - AP_REF[i].im) <= AP_TOL * hypot(AP_REF[i].re, AP_REF[i].im),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", AP_REF[i].call,
                   AP_REF[i].re, AP_REF[i].im, re, im);
    }
    /* real inside the real region, complex beyond a cut */
    assert_true("Head[AppellF1[0.3, 0.2, 0.5, 1.7, 0.4, 0.3]] === Real");
    assert_true("Head[AppellF1[0.3, 0.2, 0.5, 1.7, 0.4, 3.]] === Complex");
    assert_true("Head[AppellF4[0.3, 0.2, 1.5, 1.7, -2., -3.]] === Real");
    /* the machine reductions agree with the Gauss function */
    assert_true("AppellF1[0.3, 0.2, 0.5, 1.7, 0., 0.3] == Hypergeometric2F1[0.3, 0.5, 1.7, 0.3]");
    assert_true("Abs[AppellF1[0.3, 0.2, 0.5, 1.7, 0.35, 0.35] - Hypergeometric2F1[0.3, 0.7, 1.7, 0.35]] < 10^-15");
    assert_true("Abs[AppellF1[0.3, 0.2, 0.5, 1.7, 1., 0.3] - Hypergeometric2F1[0.3, 0.2, 1.7, 1.] Hypergeometric2F1[0.3, 0.5, 1.5, 0.3]] < 10^-14");
    /* the convergence boundary of F4, sqrt x + sqrt y = 0.995 */
    assert_true("Abs[AppellF4[0.3, 0.2, 1.5, 1.7, 0.2, 0.3] - 1.0248673980162395192] < 10^-13");
    /* continuity across a convergence-region boundary: |x| = 1 */
    assert_true("Abs[AppellF1[0.3, 0.2, 0.5, 1.7, -0.999999, 0.3] - AppellF1[0.3, 0.2, 0.5, 1.7, -1.000001, 0.3]] < 10^-6");
}

static void test_ap_arbitrary(void) {
    assert_true("Abs[N[AppellF1[3/10, 1/5, 1/2, 17/10, 2/5, 3/10], 40] - 1.047450543745694317642557434117474210829827497443462`50] < 10^-38");
    assert_true("Abs[N[AppellF2[3/10, 1/5, 1/2, 3/2, 17/10, 1/5, 3/10], 40] - 1.039753749308698594703881667299414690847218760156`50] < 10^-38");
    assert_true("Abs[N[AppellF3[3/10, 1/5, 1/2, 2/5, 17/10, 1/5, 3/10], 40] - 1.034914950081572100741246131442868304235057855667`50] < 10^-38");
    assert_true("Abs[N[AppellF4[3/10, 1/5, 3/2, 17/10, 1/10, 1/5], 40] - 1.012680246841911638166466047242357991344084093702`50] < 10^-38");
    assert_true("Abs[N[AppellF1[3/10 + I/5, 1/5, 1/2, 17/10, 2/5 + I/10, 3/10], 40] - "
                "(1.042767524522606417877581850932776449512064590755`50 + 0.037493282971595458054623047121074050154317584148`50 I)] < 10^-38");
    /* the Euler integral at 30 digits */
    assert_true("Abs[N[AppellF1[3/10, 1/5, 1/2, 17/10, 10, -20], 30] - "
                "(0.6611589841241862343962620927147725804485`40 - 0.09761940125683091706686041998122671523739`40 I)] < 10^-28");
    assert_true("Abs[N[AppellF1[3/10, 1/5, 1/2, 17/10, 2, 1/2 + 2 I], 30] - "
                "(1.06681339081010241483059791802454917615970026407`50 + 0.102271017514469769305556648952261798461677939426`50 I)] < 10^-18");
    assert_true("Precision[N[AppellF4[3/10, 1/5, 3/2, 17/10, 1/10, 1/5], 50]] >= 49");
    assert_true("MachineNumberQ[AppellF2[0.3`30, 0.2`30, 0.5, 1.5, 1.7, 0.2, 0.3]]");
}

static void test_ap_identities(void) {
    /* F1(b1 + 1) - F1 = a x/c F1(a + 1; b1 + 1, b2; c + 1) and the F2/F3/F4
     * analogues, at 30 digits */
    assert_true("With[{a = 3/10, b1 = 1/5, b2 = 1/2, c = 17/10, x = 2/5, y = 3/10}, "
                "Abs[N[AppellF1[a, b1 + 1, b2, c, x, y] - AppellF1[a, b1, b2, c, x, y] - a x/c AppellF1[a + 1, b1 + 1, b2, c + 1, x, y], 30]] < 10^-28]");
    assert_true("With[{a = 3/10, b1 = 1/5, b2 = 1/2, c1 = 3/2, c2 = 17/10, x = 1/5, y = 3/10}, "
                "Abs[N[AppellF2[a, b1 + 1, b2, c1, c2, x, y] - AppellF2[a, b1, b2, c1, c2, x, y] - a x/c1 AppellF2[a + 1, b1 + 1, b2, c1 + 1, c2, x, y], 30]] < 10^-28]");
    assert_true("With[{a1 = 3/10, a2 = 1/5, b1 = 1/2, b2 = 2/5, c = 17/10, x = 1/5, y = 3/10}, "
                "Abs[N[AppellF3[a1, a2, b1 + 1, b2, c, x, y] - AppellF3[a1, a2, b1, b2, c, x, y] - a1 x/c AppellF3[a1 + 1, a2, b1 + 1, b2, c + 1, x, y], 30]] < 10^-28]");
    /* the F1 system (DLMF 16.14.1): x(1-x) F_xx + y(1-x) F_xy + (c - (a+b1+1) x) F_x - b1 y F_y - a b1 F = 0 */
    assert_true("With[{F = AppellF1[3/10, 1/5, 1/2, 17/10, x, y]}, "
                "Abs[N[(x (1 - x) D[F, x, x] + y (1 - x) D[F, x, y] + (17/10 - (3/10 + 1/5 + 1) x) D[F, x] - 1/5 y D[F, y] - 3/10 1/5 F) /. {x -> 2/5, y -> 3/10}, 30]] < 10^-27]");
    /* the F4 system: x(1-x) F_xx - y^2 F_yy - 2 x y F_xy + (c1 - (a+b+1) x) F_x - (a+b+1) y F_y - a b F = 0 */
    assert_true("With[{F = AppellF4[3/10, 1/5, 3/2, 17/10, x, y]}, "
                "Abs[N[(x (1 - x) D[F, x, x] - y^2 D[F, y, y] - 2 x y D[F, x, y] + (3/2 - (3/10 + 1/5 + 1) x) D[F, x] - (3/10 + 1/5 + 1) y D[F, y] - 3/10 1/5 F) /. {x -> 1/10, y -> 1/5}, 30]] < 10^-27]");
    /* the Euler substitution T3 against the identity across the cut-free region */
    assert_true("With[{x = -3. + 1. I, y = 0.5 - 0.7 I}, Abs[AppellF1[1.2, 0.8, -0.3, 2.6, x, y] - "
                "(1 - x)^-0.8 (1 - y)^0.3 AppellF1[2.6 - 1.2, 0.8, -0.3, 2.6, x/(x - 1), y/(y - 1)]] < 10^-13]");
}

static void test_ap_derivative(void) {
    assert_true("D[AppellF1[a, b1, b2, c, x, y], x] === a b1/c AppellF1[1 + a, 1 + b1, b2, 1 + c, x, y]");
    assert_true("D[AppellF1[a, b1, b2, c, x, y], y] === a b2/c AppellF1[1 + a, b1, 1 + b2, 1 + c, x, y]");
    assert_true("D[AppellF2[a, b1, b2, c1, c2, x, y], y] === a b2/c2 AppellF2[1 + a, b1, 1 + b2, c1, 1 + c2, x, y]");
    assert_true("D[AppellF3[a1, a2, b1, b2, c, x, y], x] === a1 b1/c AppellF3[1 + a1, a2, 1 + b1, b2, 1 + c, x, y]");
    assert_true("D[AppellF4[a, b, c1, c2, x, y], y] === a b/c2 AppellF4[1 + a, 1 + b, c1, 1 + c2, x, y]");
    assert_true("D[AppellF1[a, b1, b2, c, x, y], a] === Derivative[1, 0, 0, 0, 0, 0][AppellF1][a, b1, b2, c, x, y]");
    assert_true("With[{h = 10^-20}, Abs[N[(D[AppellF1[3/10, 1/5, 1/2, 17/10, x, 3/10], x] /. x -> 2/5) - "
                "(AppellF1[3/10, 1/5, 1/2, 17/10, 2/5 + h, 3/10] - AppellF1[3/10, 1/5, 1/2, 17/10, 2/5 - h, 3/10])/(2 h), 30]] < 10^-28]");
}

static void test_ap_buffers_compile(void) {
    static const char* const BODY[4] = {
        "AppellF1[0.3, 0.2, 0.5, 1.7, #, 0.3]&",
        "AppellF2[0.3, 0.2, 0.5, 1.5, 1.7, #, 0.2]&",
        "AppellF3[0.3, 0.2, 0.5, 0.4, 1.7, #, 0.3]&",
        "AppellF4[0.3, 0.2, 1.5, 1.7, #, 0.1]&",
    };
    char buf[1024];
    for (int k = 0; k < 4; k++) {
        /* scalar == packed == NDArray == Compile[] scalar and rank-1, bit for bit */
        snprintf(buf, sizeof buf,
                 "With[{f = %s, w = Range[-0.4, 0.4, 0.004]}, "
                 "With[{p = f[w], t = Table[f[x], {x, FromPackedArray[w]}], "
                 "c1 = Compile[{{x, _Real}}, Evaluate[f[x]]], c2 = Compile[{{v, _Real, 1}}, Evaluate[f[v]]]}, "
                 "PackedArrayQ[p] && p === t && Normal[f[NDArray[w]]] === t && (c1 /@ w) === t && c2[w] === t]]",
                 BODY[k]);
        assert_true(buf);
        snprintf(buf, sizeof buf, "Lookup[CompileDiagnostics[{{x, _Real}}, Evaluate[(%s)[x]]], \"Compiled\"]", BODY[k]);
        assert_true(buf);
        snprintf(buf, sizeof buf, "Lookup[CompileDiagnostics[{{v, _Real, 1}}, Evaluate[(%s)[v]]], \"Compiled\"]", BODY[k]);
        assert_true(buf);
    }
    /* complex arguments through the buffer and Compile[] */
    assert_true("With[{w = Table[0.3 Exp[I t], {t, 0., 6., 0.1}]}, "
                "AppellF1[0.3, 0.2, 0.5, 1.7, w, 0.3] === Table[AppellF1[0.3, 0.2, 0.5, 1.7, z, 0.3], {z, w}] && "
                "Compile[{{z, _Complex, 1}}, AppellF1[0.3, 0.2, 0.5, 1.7, z, 0.3]][w] === Table[AppellF1[0.3, 0.2, 0.5, 1.7, z, 0.3], {z, w}]]");
    /* outside the double-series region and on a cut, the buffer still agrees */
    assert_true("With[{w = Range[1.5, 6., 0.25]}, AppellF1[0.3, 0.2, 0.5, 1.7, 0.5, w] === Table[AppellF1[0.3, 0.2, 0.5, 1.7, 0.5, y], {y, FromPackedArray[w]}]]");
    /* arbitrary-precision Compile[]: compiled == interpreted */
    assert_true("With[{cf = Compile[{{x, _Real}}, AppellF1[1/3, 2/5, 3/7, 7/4, x/3, 1/4], WorkingPrecision -> 30]}, "
                "Abs[cf[N[7/10, 30]] - N[AppellF1[1/3, 2/5, 3/7, 7/4, 7/30, 1/4], 30]] < 10^-28]");
    assert_true("Lookup[CompileDiagnostics[{{x, _Real}}, AppellF4[1/3, 2/5, 3/2, 7/4, x/10, 1/10], WorkingPrecision -> 30], \"Compiled\"]");
    assert_true("With[{cf = Compile[{{x, _Complex}}, AppellF2[1/3, 2/5, 3/7, 3/2, 7/4, x/4, 1/5], WorkingPrecision -> 30]}, "
                "Abs[cf[N[7/10 + I/5, 30]] - N[AppellF2[1/3, 2/5, 3/7, 3/2, 7/4, (7/10 + I/5)/4, 1/5], 30]] < 10^-28]");
}

static void test_ap_attributes(void) {
    const char* names[4] = { "AppellF1", "AppellF2", "AppellF3", "AppellF4" };
    char buf[128];
    for (int k = 0; k < 4; k++) {
        snprintf(buf, sizeof buf, "Attributes[%s]", names[k]);
        assert_eval_eq(buf, "{Listable, NumericFunction, Protected}", 0);
        SymbolDef* d = symtab_get_def(names[k]);
        ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, names[k]));
    }
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_ap_exact);
    TEST(test_ap_machine);
    TEST(test_ap_arbitrary);
    TEST(test_ap_identities);
    TEST(test_ap_derivative);
    TEST(test_ap_buffers_compile);
    TEST(test_ap_attributes);
    printf("All Appell tests passed.\n");
    return 0;
}
