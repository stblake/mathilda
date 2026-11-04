/* Tests for MathieuC[a, q, z] and MathieuCPrime[a, q, z], the even Mathieu
 * function (the even solution of y'' + (a - 2 q cos 2z) y = 0) and its
 * z-derivative, and for MathieuS[a, q, z] / MathieuSPrime[a, q, z], the odd
 * one (the test_mathieus_* functions at the end).
 *
 * Reference values:
 *   - real a, q in a stability band, and every characteristic value: the
 *     reference CAS (25 digits), where the normalization agrees;
 *   - characteristic orders also cross-checked against scipy.special.mathieu_cem;
 *   - complex nu (instability bands, complex a or q): Mathilda uses the analytic
 *     normalization sum c_k^2 = 1, which differs from the reference by a
 *     z-INDEPENDENT factor; tested as (a) ratio-to-reference constant in z,
 *     (b) the ODE itself, (c) high-precision self-consistency;
 *   - an independent mpmath ODE integration for y_e(5i) (a cancellation case
 *     the reference itself gets wrong: it returns ~1e45 for a value ~1e-3).
 *
 * Covers: exact rules (q = 0, a = 0, parity, MathieuCPrime[a, q, 0]), symbolic
 * passthrough; machine and arbitrary precision, precision tracking, 500-digit
 * agreement; complex parameters and arguments; instability bands (real input,
 * purely imaginary value); characteristic values (DLMF normalization, snapping);
 * large |z|; the ODE residual; Listable threading; packed / NDArray buffers
 * (bit-for-bit with the scalar path, complex output buffers); Compile[] scalar
 * and array (CompileDiagnostics) and the complex-value fallback; D (table of
 * higher derivatives, chain rule, inert a / q derivatives); Series at 0 and at a
 * generic point; argument-count message under Check / Quiet; attributes and
 * docstrings.
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

/* ---- helpers --------------------------------------------------------- */

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
    ASSERT_MSG(fabs(v - expected) <= tol * fmax(1.0, fabs(expected)),
               "%s: expected %.17g, got %.17g", input, expected, v);
}

/* An expression that must evaluate to True. */
static void assert_true(const char* input) {
    assert_eval_eq(input, "True", 0);
}

/* ---- exact rules and passthrough ------------------------------------ */

void test_mathieu_exact() {
    /* q = 0: y'' + a y = 0 */
    assert_eval_eq("MathieuC[a, 0, z]", "Cos[Sqrt[a] z]", 0);
    assert_eval_eq("MathieuC[0, 0, z]", "1/Sqrt[2]", 0);
    assert_eval_eq("MathieuC[4, 0, z]", "Cos[2 z]", 0);
    assert_eval_eq("MathieuC[1, 0, z]", "Cos[z]", 0);
    assert_eval_eq("MathieuCPrime[a, 0, z]", "-Sqrt[a] Sin[Sqrt[a] z]", 0);
    assert_eval_eq("MathieuCPrime[0, 0, z]", "0", 0);
    /* ce is even: ce'(0) = 0, ce(-z) = ce(z), ce'(-z) = -ce'(z) */
    assert_eval_eq("MathieuCPrime[a, q, 0]", "0", 0);
    assert_eval_eq("MathieuC[a, q, -x]", "MathieuC[a, q, x]", 0);
    assert_eval_eq("MathieuCPrime[a, q, -x]", "-MathieuCPrime[a, q, x]", 0);
    assert_eval_eq("MathieuC[2, 1, -3]", "MathieuC[2, 1, 3]", 0);
    /* symbolic / all exact: unevaluated */
    assert_eval_eq("MathieuC[a, q, z]", "MathieuC[a, q, z]", 0);
    assert_eval_eq("MathieuC[2, 1, 3]", "MathieuC[2, 1, 3]", 0);
    assert_eval_eq("MathieuC[a, q, 0]", "MathieuC[a, q, 0]", 0);
    assert_eval_eq("MathieuC[2, 1, Infinity]", "MathieuC[2, 1, Infinity]", 0);
    /* even under the machine path: the value at -z is the value at z */
    assert_true("MathieuC[2, 1, -3.2] == MathieuC[2, 1, 3.2]");
    assert_true("MathieuCPrime[2, 1, -3.2] == -MathieuCPrime[2, 1, 3.2]");
}

/* ---- machine precision, real stability band (reference CAS) ---------- */

void test_mathieu_machine_real() {
    assert_close("MathieuC[2, 1, 3.2]", -0.89953116344146128307606187831, 1e-14);
    assert_close("MathieuC[2, 1, 3.7]", -0.70674734130376727046760504036, 1e-14);
    assert_close("MathieuC[3, 2, 0.]", 0.81408412726292871460751436273, 1e-14);
    assert_close("MathieuC[3, 2, 1.]", 0.73104040549976430240580334022, 1e-14);
    assert_close("MathieuC[1/2, 1/10, 20.]", -0.21015731344957690966423009706, 1e-13);
    assert_close("MathieuC[7, -3, 4.]", -0.32944103159019371059190812903, 1e-13);
    assert_close("MathieuC[50, 1, 1.]", 0.74740289221627282298791426598, 1e-14);
    assert_close("MathieuC[1/3, 1/7, -2.]", 0.45060409403831822374755735645, 1e-14);
    assert_close("MathieuC[3, 1, 1000.]", 1.0449929979990154020196493699, 1e-12);
    assert_close("MathieuC[2, 1, 0.]", 1.18244270657309635137377646072, 1e-14);
    assert_close("MathieuC[0, 0, 1.]", 0.70710678118654752440084436210, 1e-15);
    /* derivatives */
    assert_close("MathieuCPrime[2, 1, 3.2]", 0.33593784114337854670788229111, 1e-14);
    assert_close("MathieuCPrime[2, 1, 3.7]", 0.50336559268757053224923503277, 1e-14);
    assert_close("MathieuCPrime[3, 2, 1.]", -0.99705871903966667681570557936, 1e-14);
    assert_close("MathieuCPrime[7, -3, 4.]", 1.2914553022080034536915350588, 1e-13);
    assert_close("MathieuCPrime[50, 1, 1.]", -4.7042576442342536332873572077, 1e-13);
    assert_close("MathieuCPrime[2, 1, 1.]", -1.1394334463069862, 1e-13);
    /* a machine-number a or q works the same */
    assert_close("MathieuC[2., 1., 3.2]", -0.89953116344146128307606187831, 1e-14);
}

/* ---- characteristic values: DLMF normalization ----------------------- */

void test_mathieu_characteristic() {
    /* a = MathieuCharacteristicA[r, q] to 25 digits; machine a within an ulp
     * of a characteristic value is treated as one (reference behaves alike).
     * The reference values agree with scipy.special.mathieu_cem. */
    assert_close("MathieuC[-0.45513860410741354823, 1, 0.7]", 0.60361682181328065331, 1e-12);
    assert_close("MathieuC[-5.80004602085150841910, 5, 0.7]", 0.30434158630740801209, 1e-12);
    assert_close("MathieuC[1.85910807251436347233, 1, 0.7]", 0.82210699099760542957, 1e-12);
    assert_close("MathieuC[1.85818754154775088515, 5, 0.7]", 0.76937012801699253492, 1e-12);
    assert_close("MathieuC[4.37130098273508566117, 1, 0.7]", 0.45384183378664784025, 1e-12);
    assert_close("MathieuC[7.44910973952917781248, 5, 0.7]", 0.90961602310841831751, 1e-12);
    assert_close("MathieuC[9.07836884720310199254, 1, 0.7]", -0.33364463542107748773, 1e-12);
    assert_close("MathieuC[11.5488320363434007962, 5, 0.7]", 0.39696206537488828887, 1e-12);
    /* fractional order: a generic (non-integer) nu, ordinary normalization */
    assert_close("MathieuC[-0.30728535063196866908, 1, 0.7]", 0.67394512136592127415, 1e-12);
    assert_close("MathieuC[2.53718008711990169560, 1, 0.7]", 0.88038489821408165598, 1e-12);
    assert_close("MathieuC[1.97333341366726213923, 5, 0.7]", 0.79484559942375932724, 1e-12);
    /* the same characteristic value to 25 digits, in arbitrary precision */
    assert_true("Abs[MathieuC[18591080725143634723299175128414233451/10^37, 1, 7/10`25] - "
                "0.82210699099760542957016932632`25] < 10^-22");
}

/* ---- arbitrary precision -------------------------------------------- */

void test_mathieu_arbitrary() {
    assert_true("Abs[MathieuC[2, 1, 3.2`20] - (-0.899531163441461283076061878314)] < 10^-18");
    assert_true("Abs[N[MathieuC[2, 1, 37/10], 50] - "
                "(-0.70674734130376727046760504036431438649802416523011)] < 10^-48");
    /* precision tracks the input */
    assert_true("Precision[MathieuC[2, 1, 3.70000000000000000000000000000000]] > 30");
    assert_true("Precision[MathieuC[2, 1, 3.70000000000000000000000000000000]] < 34");
    assert_true("Abs[MathieuC[2, 1, 3.70000000000000000000000000000000] - "
                "(-0.7067473413037672704676050403643)] < 10^-30");
    /* 500 digits: the leading 200 against the reference */
    assert_true("Abs[MathieuC[2, 1, 3.7`500] - (-0.7067473413037672704676050403643143864980241652"
                "300004545407141204755156384628466124716621355813669573565303938859822238693"
                "00823863545750777017939516940003311631959255149895416815293887588878820538"
                "381579423)] < 10^-198");
    assert_true("Precision[MathieuC[2, 1, 3.7`500]] > 499");
    /* derivative at high precision */
    assert_true("Abs[N[MathieuCPrime[2, 1, 37/10], 30] - "
                "0.503365592687570532249235032774] < 10^-28");
    /* N[] with exact arguments; NHoldFirst keeps an exact a */
    assert_true("Abs[N[MathieuC[2, 1, 3]] - (-0.966547602955374)] < 10^-13");
    /* machine input is insufficient at z = 10^15; exact input with N is right */
    assert_true("Abs[N[MathieuC[3, 1, 10^15], 20] - (-0.78255601989076339546)] < 10^-18");
}

/* ---- complex arguments and parameters -------------------------------- */

/* Mathilda's value / the reference value must be the same constant for every z
 * at a fixed (a, q): both are the even solution, normalized differently. */
static void assert_ratio_const(const char* a, const char* q,
                               const char* z1, const char* ref1,
                               const char* z2, const char* ref2) {
    char buf[1024];
    snprintf(buf, sizeof buf,
             "Abs[MathieuC[%s, %s, %s]/(%s) - MathieuC[%s, %s, %s]/(%s)] < 10^-12",
             a, q, z1, ref1, a, q, z2, ref2);
    assert_true(buf);
}

void test_mathieu_complex() {
    /* complex z, real stable (a, q): normalization agrees with the reference */
    assert_true("Abs[MathieuC[2, 1, 3.2`20 + 0.5`20 I] - "
                "(-0.882739739030372925 + 0.156794630450643144 I)] < 10^-17");
    assert_true("Abs[MathieuC[2, 1, 1. + 2. I] - (-59.288189971519977 + 101.69613259289704 I)] < 10^-11");
    assert_true("Abs[MathieuCPrime[2, 1, 1. + 2. I] - (360.36259697842831 + 742.83547154776456 I)] < 10^-10");
    /* complex nu: ratio to the reference constant in z */
    assert_ratio_const("2 + I", "1", "0.", "0.23739244271613707620 - 1.0454038022076066924 I",
                       "1.5", "-0.64760121136912588880 - 0.022381283529922257106 I");
    assert_ratio_const("-5", "3", "0.", "0.57865600636230085968 + 0.21927845822645301171 I",
                       "0.7", "2.6426200568596460022 + 1.0014060951156436788 I");
    assert_ratio_const("4", "1", "0.", "0.35562381895235170837 + 0.56984687060280206544 I",
                       "0.5", "0.26403008517229570134 + 0.42307829161629624878 I");
    assert_ratio_const("0", "1", "0.", "0.32854925576409032318 + 0.14645488407660939157 I",
                       "1.", "0.57883351132778129104 + 0.25802217875677340411 I");
    /* the analytic normalization at the request's complex point, and its
     * agreement between machine and arbitrary precision */
    assert_true("Abs[MathieuC[I, -2, 2.5 + I] - (40.7538068411686 - 38.4189542779184 I)] < 10^-10");
    assert_true("Abs[MathieuC[I, -2, N[5/2 + I, 30]] - MathieuC[I, -2, 2.5 + I]] < 10^-11");
    /* complex parameters */
    assert_true("Abs[MathieuC[2, 1 + I, 1.] - MathieuC[2, 1 + I, 1`30]] < 10^-13");
    assert_true("Abs[MathieuC[2, I, 1.] - (-0.189893127964041 + 0.851380803526212 I)] < 10^-12");
}

/* ---- instability bands: real input, real or purely imaginary value ---- */

void test_mathieu_unstable() {
    /* a = 4, q = 1 lies in the first instability band: the even solution times
     * a normalization constant whose square is real -> purely imaginary here. */
    assert_true("Re[MathieuC[4, 1, 0.5]] == 0");
    assert_true("Abs[Im[MathieuC[4, 1, 0.5]] - 0.632114032416222] < 10^-13");
    assert_true("Head[MathieuC[-5, 3, 0.7]] === Real");
    assert_true("Abs[MathieuC[-5, 3, 0.7] - 3.27661473462632] < 10^-12");
    /* purely imaginary z with real (a, q) is real-or-imaginary too */
    assert_true("Re[MathieuC[-1, 2, 5. I]] == 0");
    /* y_e(5i) by an independent mpmath ODE integration: MathieuC = N y_e with
     * N = MathieuC[a, q, 0]. The Floquet sum cancels ~e^40 here (the machine
     * path declines and the ball kernel answers). */
    assert_true("Abs[N[MathieuC[-1, 2, 5 I], 25]/N[MathieuC[-1, 2, 0], 25] - "
                "(-0.0013952853058765040128828929104)] < 10^-22");
    assert_true("Abs[MathieuC[-1, 2, 5. I]/MathieuC[-1, 2, 0.] - "
                "(-0.0013952853058765040)] < 10^-13");
}

/* ---- the ODE itself -------------------------------------------------- */

void test_mathieu_ode() {
    /* (ce')' = (2 q Cos[2 z] - a) ce, by a central difference at 40 digits,
     * for real, unstable, complex-parameter and complex-argument cases */
    const char* cases[][3] = {
        { "2", "1", "13/10" }, { "4", "1", "1/2" }, { "-5", "3", "7/10" },
        { "I", "-2", "5/2 + I" }, { "2 + I", "1", "3/2" }, { "25", "30", "1" },
        { "1/2", "-7", "2 - I/3" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
            "Block[{a = %s, q = %s, z = %s, h = 10^-15}, "
            "Abs[(MathieuCPrime[a, q, N[z + h, 40]] - MathieuCPrime[a, q, N[z - h, 40]])/(2 h) - "
            "(2 q Cos[2 N[z, 40]] - a) MathieuC[a, q, N[z, 40]]] < 10^-20 * "
            "(1 + Abs[MathieuC[a, q, N[z, 40]]])]", cases[i][0], cases[i][1], cases[i][2]);
        assert_true(buf);
        /* ce' really is the derivative of ce */
        snprintf(buf, sizeof buf,
            "Block[{a = %s, q = %s, z = %s, h = 10^-15}, "
            "Abs[(MathieuC[a, q, N[z + h, 40]] - MathieuC[a, q, N[z - h, 40]])/(2 h) - "
            "MathieuCPrime[a, q, N[z, 40]]] < 10^-20 * (1 + Abs[MathieuCPrime[a, q, N[z, 40]]])]",
            cases[i][0], cases[i][1], cases[i][2]);
        assert_true(buf);
    }
}

/* ---- Listable -------------------------------------------------------- */

void test_mathieu_listable() {
    assert_eval_eq("MathieuC[a, q, {z1, z2, z3}]",
                   "{MathieuC[a, q, z1], MathieuC[a, q, z2], MathieuC[a, q, z3]}", 0);
    assert_true("MathieuC[2, 1, {1., 2.}] == {MathieuC[2, 1, 1.], MathieuC[2, 1, 2.]}");
    assert_true("MathieuC[{1, 2}, 1, 1.] == {MathieuC[1, 1, 1.], MathieuC[2, 1, 1.]}");
    assert_true("MathieuCPrime[2, {1, 2}, 1.] == {MathieuCPrime[2, 1, 1.], MathieuCPrime[2, 2, 1.]}");
}

/* ---- packed / NDArray ------------------------------------------------ */

void test_mathieu_buffers() {
    /* packed list: stays packed, and every element bit-for-bit the scalar */
    assert_true("PackedArrayQ[MathieuC[2, 1, Range[0., 10., 0.01]]]");
    assert_true("With[{v = Range[0., 10., 0.01]}, "
                "MathieuC[2, 1, v] === Table[MathieuC[2, 1, x], {x, v}]]");
    assert_true("With[{v = Range[0., 10., 0.01]}, "
                "MathieuCPrime[3, -2, v] === Table[MathieuCPrime[3, -2, x], {x, v}]]");
    /* unpacked List == packed == NDArray, exactly */
    assert_true("With[{v = Range[0., 3., 0.25]}, "
                "Normal[MathieuC[7, 2, NDArray[v]]] === MathieuC[7, 2, FromPackedArray[v]]]");
    assert_eval_eq("Head[MathieuC[2, 1, NDArray[{1.5, 2.5, 3.5}]]]", "NDArray", 0);
    /* a complex-valued result widens the buffer (never drops the imaginary part) */
    assert_eval_eq("Head[MathieuC[4, 1, NDArray[{0.5, 1.}]]]", "NDArray", 0);
    assert_true("Normal[MathieuC[4, 1, NDArray[{0.5, 1.}]]] == {MathieuC[4, 1, 0.5], MathieuC[4, 1, 1.]}");
    assert_true("PackedArrayQ[MathieuC[4, 1, Range[0., 10., 0.01]]]");
    assert_true("With[{v = Range[0., 1., 0.01]}, MathieuC[4, 1, v][[51]] === MathieuC[4, 1, v[[51]]]]");
    /* complex buffer */
    assert_true("With[{v = Range[0., 3., 0.25] + 0.3 I}, "
                "Normal[MathieuC[2, 1, NDArray[v]]] === MathieuC[2, 1, FromPackedArray[v]]]");
    /* an int64 buffer holds exact values: unevaluated element-wise, as a List */
    assert_eval_eq("MathieuC[2, 1, NDArray[{1, 2}, DataType -> \"int64\"]]",
                   "{MathieuC[2, 1, 1], MathieuC[2, 1, 2]}", 0);
    /* a buffer in the a slot */
    assert_true("MathieuC[N[Range[1, 120]/10], 1, 0.7][[37]] == MathieuC[3.7, 1, 0.7]");
}

/* ---- Compile[] -------------------------------------------------------- */

void test_mathieu_compile() {
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, MathieuC[2., 1., x]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, MathieuC[2., 1., v]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, MathieuCPrime[2., 1., x]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{a, _Real}, {q, _Real}, {x, _Real}}, "
                   "MathieuC[a, q, x]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{x, _Complex}}, MathieuC[2., 1., x]][[1]]",
                   "\"Compiled\" -> True", 0);
    /* compiled == scalar, bit-for-bit */
    assert_true("Compile[{{x, _Real}}, MathieuC[2., 1., x]][3.2] === MathieuC[2., 1., 3.2]");
    assert_true("Compile[{{x, _Real}}, MathieuCPrime[2., 1., x]][3.2] === MathieuCPrime[2., 1., 3.2]");
    assert_true("With[{w = Range[0., 5., 0.5]}, "
                "Max[Abs[Compile[{{v, _Real, 1}}, MathieuC[2., 1., v]][w] - MathieuC[2., 1., w]]] == 0]");
    assert_true("Compile[{{x, _Complex}}, MathieuC[2., 1., x]][1. + 0.5 I] === MathieuC[2., 1., 1. + 0.5 I]");
    /* a real slot whose value is imaginary hands the call to the interpreter */
    assert_true("Compile[{{x, _Real}}, MathieuC[4., 1., x]][0.5] == MathieuC[4., 1., 0.5]");
    /* auto-compiled numeric loops agree */
    assert_true("Abs[Total[Table[MathieuC[2, 1, x], {x, 0., 1., 0.25}]] - "
                "Total[MathieuC[2, 1, {0., 0.25, 0.5, 0.75, 1.}]]] < 10^-14");
}

/* ---- derivatives ------------------------------------------------------ */

void test_mathieu_derivative() {
    assert_eval_eq("D[MathieuC[a, q, x], x]", "MathieuCPrime[a, q, x]", 0);
    assert_eval_eq("D[MathieuCPrime[a, q, x], x]",
                   "(-a + 2 q Cos[2 x]) MathieuC[a, q, x]", 0);
    assert_eval_eq("Table[D[MathieuC[a, q, x], {x, n}], {n, 1, 4}]",
        "{MathieuCPrime[a, q, x], (-a + 2 q Cos[2 x]) MathieuC[a, q, x], "
        "(-a + 2 q Cos[2 x]) MathieuCPrime[a, q, x] - 4 q MathieuC[a, q, x] Sin[2 x], "
        "-8 q Cos[2 x] MathieuC[a, q, x] + (-a + 2 q Cos[2 x])^2 MathieuC[a, q, x] - "
        "8 q MathieuCPrime[a, q, x] Sin[2 x]}", 0);
    assert_eval_eq("D[MathieuC[a, q, x^2], x]", "2 x MathieuCPrime[a, q, x^2]", 0);
    assert_eval_eq("D[MathieuC[a, q, x], a]", "Derivative[1, 0, 0][MathieuC][a, q, x]", 0);
    assert_eval_eq("D[MathieuC[a, q, x], q]", "Derivative[0, 1, 0][MathieuC][a, q, x]", 0);
    assert_eval_eq("D[MathieuC[a, q, y], x]", "0", 0);
}

/* ---- Series ----------------------------------------------------------- */

void test_mathieu_series() {
    assert_eval_eq("Series[MathieuC[3, 2, z], {z, 0, 4}]",
        "MathieuC[3, 2, 0] + 1/2 MathieuC[3, 2, 0] z^2 - 5/8 MathieuC[3, 2, 0] z^4 + O[z]^5", 0);
    assert_eval_eq("Series[MathieuC[3, 2, x], {x, 0, 8}]",
        "MathieuC[3, 2, 0] + 1/2 MathieuC[3, 2, 0] x^2 - 5/8 MathieuC[3, 2, 0] x^4 - "
        "47/720 MathieuC[3, 2, 0] x^6 + 473/4480 MathieuC[3, 2, 0] x^8 + O[x]^9", 0);
    assert_eval_eq("Series[MathieuC[a, q, x], {x, x0, 2}]",
        "MathieuC[a, q, x0] + MathieuCPrime[a, q, x0] (x - x0) + "
        "1/2 (-a + 2 q Cos[2 x0]) MathieuC[a, q, x0] (x - x0)^2 + O[x - x0]^3", 0);
    /* the Taylor coefficients agree with the function numerically */
    assert_true("Abs[N[(Normal[Series[MathieuC[3, 2, x], {x, 0, 12}]] /. x -> 1/10), 30] - "
                "MathieuC[3, 2, 0.1`30]] < 10^-14");
}

/* ---- argument count, messages ---------------------------------------- */

void test_mathieu_argcount() {
    assert_eval_eq("MathieuC[]", "MathieuC[]", 0);
    assert_eval_eq("MathieuC[1, 2]", "MathieuC[1, 2]", 0);
    assert_eval_eq("MathieuCPrime[1, 2, 3, 4]", "MathieuCPrime[1, 2, 3, 4]", 0);
    /* the diagnostic routes through the message funnel: Check sees it, Quiet
     * silences it */
    assert_eval_eq("Check[MathieuC[], \"failed\"]", "\"failed\"", 0);
    assert_eval_eq("Quiet[MathieuC[1, 2]]", "MathieuC[1, 2]", 0);
}

/* ---- attributes and docstrings --------------------------------------- */

void test_mathieu_attributes() {
    assert_eval_eq("Attributes[MathieuC]",
                   "{Listable, NHoldFirst, NumericFunction, Protected}", 0);
    assert_eval_eq("Attributes[MathieuCPrime]",
                   "{Listable, NHoldFirst, NumericFunction, Protected}", 0);
    SymbolDef* d = symtab_get_def("MathieuC");
    ASSERT_MSG(d && d->docstring && strstr(d->docstring, "MathieuC[a, q, z]"),
               "MathieuC docstring missing");
    d = symtab_get_def("MathieuCPrime");
    ASSERT_MSG(d && d->docstring && strstr(d->docstring, "MathieuCPrime[a, q, z]"),
               "MathieuCPrime docstring missing");
}

/* ==================================================================== */
/* MathieuS / MathieuSPrime                                              */
/* ==================================================================== */
/* Reference values: the reference CAS (25 digits) wherever nu is real -- a
 * stability band, and every characteristic value b_n -- where MathieuC + I
 * MathieuS is exactly the Floquet solution in both systems. For complex nu the
 * reference's MathieuS is NOT a solution of the Mathieu equation (its ODE
 * residual at MathieuS[I, -2, 3/2 + 3 I/10] is ~3), so those cases are tested
 * by the ODE, oddness, the Floquet property and self-consistency. */

/* ---- exact rules and passthrough ------------------------------------ */

void test_mathieus_exact() {
    assert_eval_eq("MathieuS[a, 0, z]", "Sin[Sqrt[a] z]", 0);
    assert_eval_eq("MathieuS[4, 0, z]", "Sin[2 z]", 0);
    assert_eval_eq("MathieuS[1, 0, z]", "Sin[z]", 0);
    assert_eval_eq("MathieuS[0, 0, z]", "0", 0);
    assert_eval_eq("MathieuSPrime[a, 0, z]", "Sqrt[a] Cos[Sqrt[a] z]", 0);
    assert_eval_eq("MathieuSPrime[0, 0, z]", "0", 0);
    /* se is odd: se(0) = 0, se' even */
    assert_eval_eq("MathieuS[a, q, 0]", "0", 0);
    assert_eval_eq("MathieuS[2, 1, 0]", "0", 0);
    assert_eval_eq("MathieuSPrime[a, q, 0]", "MathieuSPrime[a, q, 0]", 0);
    assert_eval_eq("MathieuS[a, q, -x]", "-MathieuS[a, q, x]", 0);
    assert_eval_eq("MathieuSPrime[a, q, -x]", "MathieuSPrime[a, q, x]", 0);
    assert_eval_eq("MathieuS[2, 1, -3]", "-MathieuS[2, 1, 3]", 0);
    /* exact input stays symbolic */
    assert_eval_eq("MathieuS[a, q, z]", "MathieuS[a, q, z]", 0);
    assert_eval_eq("MathieuS[2, 1, 3]", "MathieuS[2, 1, 3]", 0);
    assert_eval_eq("MathieuS[2, 1, Infinity]", "MathieuS[2, 1, Infinity]", 0);
    /* numeric zero and parity on machine numbers */
    assert_true("MathieuS[2., 1., 0] == 0");
    assert_true("MathieuS[2, 1, -3.2] == -MathieuS[2, 1, 3.2]");
    assert_true("MathieuSPrime[2, 1, -3.2] == MathieuSPrime[2, 1, 3.2]");
    /* q = 0 numerically: sin / i sinh */
    assert_true("Abs[MathieuS[2, 0, 0.7] - Sin[Sqrt[2.] 0.7]] < 10^-15");
    assert_true("Re[MathieuS[-1, 0, 0.5]] == 0 && Abs[Im[MathieuS[-1, 0, 0.5]] - Sinh[0.5]] < 10^-15");
}

/* ---- machine precision, real stability band -------------------------- */

void test_mathieus_machine_real() {
    assert_close("MathieuS[2, 1, 3.2]", -0.7680982888758991508884760160, 1e-14);
    assert_close("MathieuS[2, 1, 3.7]", -0.94825674959404574160708085731, 1e-14);
    assert_close("MathieuS[2, 1, 0.5]", 0.26361344360795635701423299718, 1e-14);
    assert_close("MathieuS[3, 2, 1.]", 0.99031016385316265653888270297, 1e-14);
    assert_close("MathieuS[50, 1, 1.]", 0.65778463974512344586713129198, 1e-14);
    assert_close("MathieuS[1/2, 1/10, 20.]", 1.04382654891620502058886342995, 1e-13);
    assert_close("MathieuS[7, -3, 4.]", -1.02653385375864498691545447292, 1e-14);
    assert_close("MathieuS[1/3, 1/7, -2.]", -0.96817201670654653061403812589, 1e-14);
    assert_close("MathieuS[3, 1, 1000.]", -0.17419150349928755730461857289, 1e-12);
    assert_close("MathieuSPrime[2, 1, 0.5]", 0.50226074781711933771901987424, 1e-14);
    assert_close("MathieuSPrime[2, 1, 0.]", 0.53361615142920196105596353989, 1e-14);
    assert_close("MathieuSPrime[3, 2, 1.]", -0.02381538574086996863946798565, 1e-14);
    assert_close("MathieuSPrime[7, -3, 4.]", -1.50850750640838419057975888693, 1e-14);
    assert_close("MathieuSPrime[50, 1, 1.]", 5.31772615182975184137352076959, 1e-14);
    /* small z: relative accuracy (the Taylor arm), not the Floquet sum's
     * absolute accuracy */
    assert_true("Abs[MathieuS[2, 1, 1.*^-10]/(1.*^-10 0.53361615142920196105596354) - 1] < 10^-14");
    assert_true("Abs[MathieuS[2, 1, 1.*^-3]/MathieuS[2, 1, N[1/1000, 30]] - 1] < 10^-14");
    assert_true("Abs[MathieuS[50, 7, 0.01]/MathieuS[50, 7, 0.01`30] - 1] < 10^-14");
    assert_true("Abs[MathieuCPrime[3, 2, 1.*^-8]/MathieuCPrime[3, 2, N[10^-8, 30]] - 1] < 10^-14");
    assert_true("Abs[MathieuS[2, 1, 1.*^-10 (1 + I)]/MathieuS[2, 1, N[10^-10 (1 + I), 30]] - 1] < 10^-14");
}

/* ---- characteristic values ------------------------------------------- */

void test_mathieus_characteristic() {
    /* a = b_n(q): se_n with the DLMF normalization (b_n to 25 digits; a
     * machine a within an ulp snaps) */
    assert_close("MathieuS[3.917024772998471, 1., 0.5]", 0.76379180880612364073161728954, 1e-12);
    assert_close("MathieuS[2.0994604454866654, 5., 0.7]", 0.75768033899141663541702755179, 1e-12);
    assert_close("MathieuS[9.047739259809375, 1., 0.7]", 0.94726754585766302921814623858, 1e-12);
    assert_close("MathieuS[-5.790080598637771, 5., 0.7]", 0.30060395258588347683590532579, 1e-12);
    assert_close("MathieuS[-0.11024881699209517, 1., 0.5]", 0.36986839755591103147077965296, 1e-12);
    assert_true("Abs[MathieuS[2.09946044548666536395565456881965689508`25, 5, 7/10] - "
                "0.7576803389914166354170275518] < 10^-20");
    assert_true("Abs[MathieuS[9.04773925980937498237494653468680316843`25, 1, 7/10] - "
                "0.9472675458576630292181462386] < 10^-20");
    assert_true("Abs[MathieuSPrime[3.91702477299847118670341688529379991597`25, 1, 7/10] - "
                "0.6425728763755913524357839453] < 10^-20");
    /* a = a_n(q): no periodic odd solution; the Floquet-normalized odd
     * function is identically 0 (its limit from both sides) */
    assert_true("MathieuS[4.371300982735086, 1., 0.7] == 0");
    assert_true("MathieuSPrime[4.371300982735086, 1., 0.7] == 0");
    assert_true("MathieuS[4.37130098273508566117112091819948782529`25, 1, 7/10] == 0");
    assert_true("MathieuS[-0.45513860410741354823263318752888586664`25, 1, 7/10] == 0");
    assert_true("MathieuS[7.449109739529178, 5., 0.7] == 0");
    /* ... and symmetrically MathieuC at a = b_n */
    assert_true("MathieuC[3.917024772998471, 1., 0.7] == 0");
    assert_true("MathieuCPrime[3.917024772998471, 1., 0.7] == 0");
    assert_true("MathieuC[2.0994604454866654, 5., 0.7] == 0");
    assert_true("MathieuC[3.91702477299847118670341688529379991597`25, 1, 7/10] == 0");
    /* near (not at) a characteristic value the function is small, not 0 */
    assert_true("0 < Abs[MathieuS[4.3713, 1., 0.7]] < 10^-2");
    /* a band edge at larger q: nu = 2 + 0.0287 i is ill-conditioned, still
     * evaluated (MathieuC and MathieuS), and satisfies the ODE */
    assert_true("NumericQ[MathieuS[2.1, 5., 0.7]] && NumericQ[MathieuC[2.1, 5., 0.7]]");
    assert_true("Abs[MathieuS[2.1, 5., 0.7] - MathieuS[21/10, 5, 0.7`30]] < 10^-10");
    assert_true("Abs[MathieuC[2.1, 5., 0.7] - MathieuC[21/10, 5, 0.7`30]] < 10^-10");
}

/* ---- arbitrary precision --------------------------------------------- */

void test_mathieus_arbitrary() {
    assert_true("Abs[MathieuS[2, 1, 3.2`30] - (-0.7680982888758991508884760160)] < 10^-27");
    assert_true("Abs[N[MathieuS[2, 1, 37/10], 50] - "
                "(-0.94825674959404574160708085731501659004752955180031)] < 10^-48");
    assert_true("Precision[MathieuS[2, 1, 3.70000000000000000000000000000000]] > 30");
    assert_true("Precision[MathieuS[2, 1, 3.70000000000000000000000000000000]] < 34");
    assert_true("Abs[MathieuS[2, 1, 3.70000000000000000000000000000000] - "
                "(-0.9482567495940457416070808573150)] < 10^-30");
    assert_true("Abs[MathieuS[2, 1, 3.7`500] - (-0.94825674959404574160708085731501659004752955180032200691800800254810079970944587640169008500185104156959040913427876221900381404174519439217391658600839509861080967031186875016112982506290366176135128704259476145136047113713144310)] < 10^-228");
    assert_true("Precision[MathieuS[2, 1, 3.7`500]] > 499");
    /* 500 vs 600 digits agree to 500 */
    assert_true("MathieuS[2, 1, N[37/10, 500]] == MathieuS[2, 1, N[37/10, 600]]");
    assert_true("Abs[MathieuS[2, 1, N[37/10, 500]] - N[MathieuS[2, 1, N[37/10, 600]], 500]] < 10^-490 || "
                "MathieuS[2, 1, N[37/10, 500]] - MathieuS[2, 1, N[37/10, 600]] == 0");
    assert_true("Abs[N[MathieuSPrime[2, 1, 16/5], 30] - (-0.414590689653679482133477536123)] < 10^-28");
    assert_true("Abs[N[MathieuS[2, 1, 3]] - (-0.685047939537226094)] < 10^-13");
    /* machine input is insufficient at z = 10^12 (nu z to 1e-4); N is right */
    assert_true("Abs[N[MathieuS[3, 3, 10^12], 20] - (-1.1132442186984920324)] < 10^-18");
    assert_true("Abs[N[MathieuS[3, 1, 10^15], 25] - (-0.6596260006395244322902593)] < 10^-23");
    /* the ball kernel and the machine helper agree */
    assert_true("Abs[MathieuS[7, -3, 4`30] - MathieuS[7, -3, 4.]] < 10^-14");
    assert_true("Abs[MathieuSPrime[50, 1, 1`30] - MathieuSPrime[50, 1, 1.]] < 10^-13");
}

/* ---- complex arguments and parameters -------------------------------- */

void test_mathieus_complex() {
    /* complex z, real stable (a, q): the normalization agrees */
    assert_true("Abs[MathieuS[2, 1, 3.2`25 + 0.5`25 I] - "
                "(-0.75195327203699548458570281570 - 0.21251777957127455557463952437 I)] < 10^-22");
    assert_true("Abs[MathieuS[2, 1, 1. + 2. I] - (-54.015604428266988 + 11.511632826296836 I)] < 10^-12");
    assert_true("Abs[MathieuSPrime[2, 1, 1. + 2. I] - (-144.15963036897887 + 359.52975139493334 I)] < 10^-11");
    /* imaginary axis, real (a, q): se and ce' are i times a real there, se'
     * and ce real (machine and arbitrary precision alike) */
    assert_true("Abs[N[MathieuS[2, 1, I/2], 20] - 0.26335055509155453511 I] < 10^-19");
    assert_true("Abs[N[MathieuSPrime[2, 1, I/2], 20] - 0.49859470344438407657] < 10^-19");
    assert_true("Re[MathieuS[2, 1, 0.5 I]] == 0 && Im[MathieuSPrime[2, 1, 0.5 I]] == 0");
    /* regression: MathieuCPrime on the imaginary axis was zeroed (the
     * ball kernel assumed every function is real there) */
    assert_true("Abs[N[MathieuCPrime[2, 1, I/2], 20] - 0.20524473480830092561 I] < 10^-19");
    assert_true("Abs[MathieuCPrime[2, 1, 0.5 I] - 0.20524473480830092561 I] < 10^-14");
    /* complex parameters: machine and arbitrary precision agree */
    assert_true("Abs[MathieuS[I, -2, 2.5 + I] - MathieuS[I, -2, N[5/2 + I, 30]]] < 10^-11");
    assert_true("Abs[MathieuS[I, -2, 2.5 + I] - (38.494107004 + 40.641965471 I)] < 10^-8");
    assert_true("Abs[MathieuS[2, 1 + I, 1.] - MathieuS[2, 1 + I, 1`30]] < 10^-13");
    assert_true("Abs[MathieuSPrime[2 + I, 1, 1.5] - MathieuSPrime[2 + I, 1, 3/2`30]] < 10^-13");
}

/* ---- instability bands, the Floquet property, the Wronskian ---------- */

void test_mathieus_unstable() {
    /* real (a, q) in an instability band: a real normalization constant or a
     * purely imaginary one, each value accordingly */
    assert_true("Head[MathieuS[4, 1, 0.5]] === Real");
    assert_true("Re[MathieuS[-5, 3, 0.7]] == 0");
    assert_true("Re[MathieuS[2, 3, 0.5]] == 0");
    /* the 30-digit kernel decides the same kind */
    assert_true("Abs[MathieuS[-5, 3, 0.7] - MathieuS[-5, 3, 7/10`30]] < 10^-13");
    assert_true("Abs[MathieuS[4, 1, 0.5] - MathieuS[4, 1, 1/2`30]] < 10^-13");
    /* C + I S is the Floquet solution: F(z + Pi)/F(z) is the same multiplier
     * at every z, stable or not, real or complex parameters */
    const char* cases[][2] = {
        { "2", "1" }, { "4", "1" }, { "-5", "3" }, { "I", "-2" }, { "2 + I", "1" },
        { "1/2", "-7" }, { "25", "30" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
            "Block[{F = Function[z, MathieuC[%s, %s, z] + I MathieuS[%s, %s, z]]}, "
            "Abs[F[N[7/10 + Pi, 40]]/F[N[7/10, 40]] - F[N[13/10 + I/5 + Pi, 40]]/F[N[13/10 + I/5, 40]]] "
            "< 10^-25 * Abs[F[N[7/10 + Pi, 40]]/F[N[7/10, 40]]]]",
            cases[i][0], cases[i][1], cases[i][0], cases[i][1]);
        assert_true(buf);
        /* Wronskian C S' - C' S = C(0) S'(0), constant in z */
        snprintf(buf, sizeof buf,
            "Block[{a = %s, q = %s, z = N[13/10 + I/5, 40]}, "
            "Abs[MathieuC[a, q, z] MathieuSPrime[a, q, z] - MathieuCPrime[a, q, z] MathieuS[a, q, z] - "
            "MathieuC[a, q, N[0, 40]] MathieuSPrime[a, q, N[0, 40]]] < 10^-25 * "
            "(1 + Abs[MathieuC[a, q, N[0, 40]] MathieuSPrime[a, q, N[0, 40]]])]",
            cases[i][0], cases[i][1]);
        assert_true(buf);
    }
    /* the same at machine precision, real stable band */
    assert_true("Abs[MathieuC[2, 1, 1.3] MathieuSPrime[2, 1, 1.3] - "
                "MathieuCPrime[2, 1, 1.3] MathieuS[2, 1, 1.3] - "
                "MathieuC[2, 1, 0.] MathieuSPrime[2, 1, 0.]] < 10^-14");
}

/* ---- the ODE itself -------------------------------------------------- */

void test_mathieus_ode() {
    const char* cases[][3] = {
        { "2", "1", "13/10" }, { "4", "1", "1/2" }, { "-5", "3", "7/10" },
        { "I", "-2", "5/2 + I" }, { "I", "-2", "3/2 + 3 I/10" }, { "2 + I", "1", "3/2" },
        { "25", "30", "1" }, { "1/2", "-7", "2 - I/3" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
            "Block[{a = %s, q = %s, z = %s, h = 10^-15}, "
            "Abs[(MathieuSPrime[a, q, N[z + h, 40]] - MathieuSPrime[a, q, N[z - h, 40]])/(2 h) - "
            "(2 q Cos[2 N[z, 40]] - a) MathieuS[a, q, N[z, 40]]] < 10^-20 * "
            "(1 + Abs[MathieuS[a, q, N[z, 40]]])]", cases[i][0], cases[i][1], cases[i][2]);
        assert_true(buf);
        snprintf(buf, sizeof buf,
            "Block[{a = %s, q = %s, z = %s, h = 10^-15}, "
            "Abs[(MathieuS[a, q, N[z + h, 40]] - MathieuS[a, q, N[z - h, 40]])/(2 h) - "
            "MathieuSPrime[a, q, N[z, 40]]] < 10^-20 * (1 + Abs[MathieuSPrime[a, q, N[z, 40]]])]",
            cases[i][0], cases[i][1], cases[i][2]);
        assert_true(buf);
        /* odd in z at a complex point */
        snprintf(buf, sizeof buf,
            "Block[{a = %s, q = %s, z = N[%s, 30]}, "
            "Abs[MathieuS[a, q, z] + MathieuS[a, q, -z]] < 10^-25 (1 + Abs[MathieuS[a, q, z]])]",
            cases[i][0], cases[i][1], cases[i][2]);
        assert_true(buf);
    }
}

/* ---- Listable, buffers, Compile -------------------------------------- */

void test_mathieus_buffers() {
    assert_eval_eq("MathieuS[a, q, {z1, z2, z3}]",
                   "{MathieuS[a, q, z1], MathieuS[a, q, z2], MathieuS[a, q, z3]}", 0);
    assert_true("MathieuS[2, 1, {1., 2.}] == {MathieuS[2, 1, 1.], MathieuS[2, 1, 2.]}");
    assert_true("MathieuS[{1, 2}, 1, 1.] == {MathieuS[1, 1, 1.], MathieuS[2, 1, 1.]}");
    assert_true("MathieuSPrime[2, {1, 2}, 1.] == {MathieuSPrime[2, 1, 1.], MathieuSPrime[2, 2, 1.]}");
    /* packed buffer: bit-for-bit with the scalar path, stays packed */
    assert_true("PackedArrayQ[MathieuS[2, 1, Range[0., 10., 0.01]]]");
    assert_true("With[{v = Range[-3., 3., 0.01]}, "
                "MathieuS[2, 1, v] === Table[MathieuS[2, 1, x], {x, v}]]");
    assert_true("With[{v = Range[-3., 3., 0.01]}, "
                "MathieuSPrime[3, -2, v] === Table[MathieuSPrime[3, -2, x], {x, v}]]");
    assert_true("With[{v = Range[0., 3., 0.25]}, "
                "Normal[MathieuS[7, 2, NDArray[v]]] === MathieuS[7, 2, FromPackedArray[v]]]");
    assert_eval_eq("Head[MathieuS[2, 1, NDArray[{1.5, 2.5, 3.5}]]]", "NDArray", 0);
    /* an instability band with an imaginary normalization: complex output */
    assert_true("PackedArrayQ[MathieuS[-5, 3, Range[0., 10., 0.01]]]");
    assert_true("With[{v = Range[0., 1., 0.01]}, MathieuS[-5, 3, v][[51]] === MathieuS[-5, 3, v[[51]]]]");
    assert_true("Normal[MathieuS[-5, 3, NDArray[{0.5, 1.}]]] == {MathieuS[-5, 3, 0.5], MathieuS[-5, 3, 1.]}");
    /* characteristic values over a buffer: b_n gives se_n, a_n gives 0 */
    assert_true("Max[Abs[MathieuS[4.371300982735086, 1., Range[0., 3., 0.1]]]] == 0");
    assert_true("With[{v = Range[0., 3., 0.1]}, "
                "MathieuS[3.917024772998471, 1., v] === Table[MathieuS[3.917024772998471, 1., x], {x, v}]]");
    /* int64 NDArray: delisted, exact elements stay symbolic */
    assert_eval_eq("MathieuS[2, 1, NDArray[{1, 2}, DataType -> \"int64\"]]",
                   "{MathieuS[2, 1, 1], MathieuS[2, 1, 2]}", 0);
    /* Compile: scalar and array lower, values bit-equal to the interpreter */
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, MathieuS[2., 1., x]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, MathieuS[2., 1., v]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, MathieuSPrime[2., 1., x]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, MathieuSPrime[2., 1., v]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_true("Compile[{{x, _Real}}, MathieuS[2., 1., x]][3.2] === MathieuS[2., 1., 3.2]");
    assert_true("Compile[{{x, _Real}}, MathieuSPrime[2., 1., x]][3.2] === MathieuSPrime[2., 1., 3.2]");
    assert_true("With[{w = Range[0., 5., 0.05]}, "
                "Max[Abs[Compile[{{v, _Real, 1}}, MathieuS[2., 1., v]][w] - MathieuS[2., 1., w]]] == 0]");
    assert_true("Compile[{{x, _Complex}}, MathieuS[2., 1., x]][1. + 0.5 I] === MathieuS[2., 1., 1. + 0.5 I]");
    /* a purely imaginary value from a real-typed Compile falls back intact */
    assert_true("Compile[{{x, _Real}}, MathieuS[-5., 3., x]][0.7] == MathieuS[-5., 3., 0.7]");
    /* auto-compiled Table agrees with the buffer path */
    assert_true("Abs[Total[Table[MathieuS[2, 1, x], {x, 0., 1., 0.25}]] - "
                "Total[MathieuS[2, 1, {0., 0.25, 0.5, 0.75, 1.}]]] < 10^-14");
}

/* ---- D and Series ---------------------------------------------------- */

void test_mathieus_derivative() {
    assert_eval_eq("D[MathieuS[a, q, x], x]", "MathieuSPrime[a, q, x]", 0);
    assert_eval_eq("D[MathieuSPrime[a, q, x], x]",
                   "(-a + 2 q Cos[2 x]) MathieuS[a, q, x]", 0);
    assert_eval_eq("Table[D[MathieuS[a, q, x], {x, n}], {n, 1, 4}]",
        "{MathieuSPrime[a, q, x], (-a + 2 q Cos[2 x]) MathieuS[a, q, x], "
        "(-a + 2 q Cos[2 x]) MathieuSPrime[a, q, x] - 4 q MathieuS[a, q, x] Sin[2 x], "
        "-8 q Cos[2 x] MathieuS[a, q, x] + (-a + 2 q Cos[2 x])^2 MathieuS[a, q, x] - "
        "8 q MathieuSPrime[a, q, x] Sin[2 x]}", 0);
    assert_eval_eq("D[MathieuS[a, q, x^2], x]", "2 x MathieuSPrime[a, q, x^2]", 0);
    assert_eval_eq("D[MathieuS[a, q, x], a]", "Derivative[1, 0, 0][MathieuS][a, q, x]", 0);
    assert_eval_eq("D[MathieuSPrime[a, q, x], q]", "Derivative[0, 1, 0][MathieuSPrime][a, q, x]", 0);
    assert_eval_eq("D[MathieuS[a, q, y], x]", "0", 0);
    assert_eval_eq("Series[MathieuS[3, 2, z], {z, 0, 4}]",
        "MathieuSPrime[3, 2, 0] z + 1/6 MathieuSPrime[3, 2, 0] z^3 + O[z]^5", 0);
    assert_eval_eq("Series[MathieuS[3, 2, x], {x, 0, 8}]",
        "MathieuSPrime[3, 2, 0] x + 1/6 MathieuSPrime[3, 2, 0] x^3 - "
        "47/120 MathieuSPrime[3, 2, 0] x^5 + 113/5040 MathieuSPrime[3, 2, 0] x^7 + O[x]^9", 0);
    assert_eval_eq("Series[MathieuS[a, q, x], {x, x0, 2}]",
        "MathieuS[a, q, x0] + MathieuSPrime[a, q, x0] (x - x0) + "
        "1/2 (-a + 2 q Cos[2 x0]) MathieuS[a, q, x0] (x - x0)^2 + O[x - x0]^3", 0);
    /* the truncated series agrees numerically */
    assert_true("Abs[N[(Normal[Series[MathieuS[3, 2, x], {x, 0, 13}]] /. x -> 1/10), 30] - "
                "MathieuS[3, 2, 0.1`30]] < 10^-14");
}

/* ---- argument count, attributes, docstrings -------------------------- */

void test_mathieus_argcount_attributes() {
    assert_eval_eq("MathieuS[]", "MathieuS[]", 0);
    assert_eval_eq("MathieuS[1, 2]", "MathieuS[1, 2]", 0);
    assert_eval_eq("MathieuSPrime[1, 2, 3, 4]", "MathieuSPrime[1, 2, 3, 4]", 0);
    assert_eval_eq("Check[MathieuS[], \"failed\"]", "\"failed\"", 0);
    assert_eval_eq("Quiet[MathieuS[1, 2]]", "MathieuS[1, 2]", 0);
    assert_eval_eq("Check[MathieuSPrime[1], \"failed\"]", "\"failed\"", 0);
    assert_eval_eq("Attributes[MathieuS]",
                   "{Listable, NHoldFirst, NumericFunction, Protected}", 0);
    assert_eval_eq("Attributes[MathieuSPrime]",
                   "{Listable, NHoldFirst, NumericFunction, Protected}", 0);
    SymbolDef* d = symtab_get_def("MathieuS");
    ASSERT_MSG(d && d->docstring && strstr(d->docstring, "MathieuS[a, q, z]"),
               "MathieuS docstring missing");
    d = symtab_get_def("MathieuSPrime");
    ASSERT_MSG(d && d->docstring && strstr(d->docstring, "MathieuSPrime[a, q, z]"),
               "MathieuSPrime docstring missing");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_mathieu_exact);
    TEST(test_mathieu_machine_real);
    TEST(test_mathieu_characteristic);
    TEST(test_mathieu_arbitrary);
    TEST(test_mathieu_complex);
    TEST(test_mathieu_unstable);
    TEST(test_mathieu_ode);
    TEST(test_mathieu_listable);
    TEST(test_mathieu_buffers);
    TEST(test_mathieu_compile);
    TEST(test_mathieu_derivative);
    TEST(test_mathieu_series);
    TEST(test_mathieu_argcount);
    TEST(test_mathieu_attributes);

    TEST(test_mathieus_exact);
    TEST(test_mathieus_machine_real);
    TEST(test_mathieus_characteristic);
    TEST(test_mathieus_arbitrary);
    TEST(test_mathieus_complex);
    TEST(test_mathieus_unstable);
    TEST(test_mathieus_ode);
    TEST(test_mathieus_buffers);
    TEST(test_mathieus_derivative);
    TEST(test_mathieus_argcount_attributes);

    printf("All MathieuC / MathieuCPrime / MathieuS / MathieuSPrime tests passed.\n");
    return 0;
}
