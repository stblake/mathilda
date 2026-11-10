/* Tests for StruveL, AngerJ, WeberE and the Kelvin family
 * (KelvinBer / KelvinBei / KelvinKer / KelvinKei).
 *
 * Covers, per head: exact reductions (each one also checked NUMERICALLY against
 * an independent definition, so no identity is taken on trust); machine values
 * against mpmath 1.3 (an independent implementation: mp.struvel, mp.angerj,
 * mp.webere, mp.ber/bei/ker/kei); 30-digit values against mpmath and against
 * the integral definitions (NIntegrate at WorkingPrecision 40); complex
 * arguments, including the ker/kei branch cut; D against a central difference
 * quotient; bit-for-bit agreement of scalar, List, packed, NDArray and Compile;
 * CompileDiagnostics at scalar and rank-1 shapes; attributes.
 *
 * This harness loads core_init() only (no init.m): every path is C.
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

/* Machine value within `ulps`-ish relative tolerance of an mpmath reference. */
static void mclose(const char* input, double expected, double rel) {
    double v = eval_real(input);
    double err = fabs(v - expected) / fabs(expected);
    ASSERT_MSG(err <= rel, "%s: expected %.17g, got %.17g (rel err %.3g)",
               input, expected, v, err);
}
#define MCLOSE(s, v) mclose(s, v, 1.0e-14)

/* An expression that must evaluate to True. */
static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* ---- exact reductions ------------------------------------------------ */

void test_exact_forms() {
    assert_eval_eq("StruveL[nu, x]", "StruveL[nu, x]", 0);
    assert_eval_eq("StruveL[1, 0]", "0", 0);
    assert_eval_eq("StruveL[nu, 0]", "StruveL[nu, 0]", 0);
    assert_eval_eq("AngerJ[2, z]", "BesselJ[2, z]", 0);
    assert_eval_eq("AngerJ[nu, 0]", "Sinc[nu Pi]", 0);
    assert_eval_eq("AngerJ[1/2, 0]", "2/Pi", 0);
    assert_eval_eq("AngerJ[1/3, z]", "AngerJ[1/3, z]", 0);
    assert_eval_eq("WeberE[0, z]", "-StruveH[0, z]", 0);
    assert_eval_eq("WeberE[1, z]", "2/Pi - StruveH[1, z]", 0);
    assert_eval_eq("WeberE[3, z]", "2/(3 Pi) + (2 z^2)/(15 Pi) - StruveH[3, z]", 0);
    assert_eval_eq("WeberE[-2, z]", "(2 z)/(3 Pi) - StruveH[2, z]", 0);
    assert_eval_eq("WeberE[-3, z]", "-2/(3 Pi) - (2 z^2)/(15 Pi) + StruveH[3, z]", 0);
    assert_eval_eq("WeberE[1/2, 0]", "2/Pi", 0);
    assert_eval_eq("WeberE[2, 0]", "0", 0);
    assert_eval_eq("KelvinBer[z]", "KelvinBer[0, z]", 0);
    assert_eval_eq("KelvinKei[z]", "KelvinKei[0, z]", 0);
    assert_eval_eq("KelvinBer[-1, z]", "-KelvinBer[1, z]", 0);
    assert_eval_eq("KelvinKer[-2, z]", "KelvinKer[2, z]", 0);
    assert_eval_eq("KelvinBei[-3, z]", "-KelvinBei[3, z]", 0);
    assert_eval_eq("KelvinBer[0, 0]", "1", 0);
    assert_eval_eq("KelvinBer[1, 0]", "0", 0);
    assert_eval_eq("KelvinBei[0, 0]", "0", 0);
    assert_eval_eq("KelvinKer[0, 0]", "Infinity", 0);
    assert_eval_eq("KelvinKer[0, 0.]", "Infinity", 0);
    assert_eval_eq("KelvinKei[0, 0]", "-1/4 Pi", 0);
    assert_eval_eq("KelvinKei[1, 0]", "ComplexInfinity", 0);
    assert_eval_eq("KelvinBer[nu, z]", "KelvinBer[nu, z]", 0);
}

/* Every reduction above, checked numerically at 30 digits against the
 * integral definition / series (independent of the rewrite). */
void test_exact_forms_numeric() {
    /* StruveL half-odd-integer closed forms vs the 30-digit arb series */
    assert_true("Abs[N[StruveL[1/2, 13/10] - (Sqrt[2/Pi] (Cosh[x] - 1)/Sqrt[x] /. x -> 13/10), 30]] < 10^-28");
    assert_true("Abs[N[(StruveL[5/2, x] /. x -> 7/3) - StruveL[5/2, N[7/3, 40]], 30]] < 10^-25");
    assert_true("Abs[N[(StruveL[-7/2, x] /. x -> 7/3) - StruveL[-7/2, N[7/3, 40]], 30]] < 10^-25");
    assert_true("Abs[N[(StruveL[9/2, x] /. x -> 21/2) - StruveL[9/2, N[21/2, 40]], 30]] < 10^-20");
    /* AngerJ / WeberE vs the integral definitions, (1/pi) Int_0^pi of
     * cos / sin(nu t - z sin t), evaluated by mpmath.quad at 40 digits (the
     * built-in NIntegrate does not reach 25 digits on these oscillatory
     * integrands, so the reference values are literal). */
    assert_true("Abs[N[AngerJ[1/3, 5/2], 30] - 0.28887833221102066738875408860393`32] < 10^-28");
    assert_true("Abs[N[(WeberE[3, z] /. z -> 5/2), 30] - 0.28962943014838320217657065029185`32] < 10^-28");
    assert_true("Abs[N[(WeberE[4, z] /. z -> 5/2), 30] - 0.1463804420499113426007520966077`32] < 10^-28");
    assert_true("Abs[N[(WeberE[-3, z] /. z -> 5/2), 30] + 0.28962943014838320217657065029185`32] < 10^-28");
    /* the z = 0 values, whose integrals are elementary */
    assert_true("Abs[N[WeberE[2/3, 0] - Integrate[Sin[2/3 t], {t, 0, Pi}]/Pi, 30]] < 10^-28");
    assert_true("Abs[N[AngerJ[2/3, 0] - Integrate[Cos[2/3 t], {t, 0, Pi}]/Pi, 30]] < 10^-28");
    /* Kelvin negative integer orders vs the arb kernel at that order */
    assert_true("Abs[N[KelvinKei[-3, 27/10], 30] + N[KelvinKei[3, 27/10], 30]] < 10^-28");
    assert_true("Abs[KelvinBer[-3, N[27/10, 30]] + KelvinBer[3, N[27/10, 30]]] < 10^-28");
    assert_true("Abs[KelvinKer[-2, N[27/10, 30]] - KelvinKer[2, N[27/10, 30]]] < 10^-28");
}

/* ---- machine values vs mpmath ---------------------------------------- */

void test_machine_values() {
    MCLOSE("StruveL[1, 2.0]", 1.1027597873677158);
    MCLOSE("StruveL[0, 5.5]", 42.57386983077723);
    MCLOSE("StruveL[2.5, 30.0]", 703124015502.88738);
    MCLOSE("StruveL[-1.7, 3.3]", 3.6976667593837008);
    MCLOSE("StruveL[1/3, 12.0]", 18857.236537636667);
    MCLOSE("StruveL[1, 45.0]", 2.0601334620815772e+18);
    MCLOSE("AngerJ[1/3, 2.0]", 0.5505802268417606);
    MCLOSE("AngerJ[2.5, 7.5]", -0.26749840947648009);
    MCLOSE("AngerJ[-1.3, 15.0]", -0.15200463240340436);
    MCLOSE("AngerJ[0.3, -4.0]", -0.41029924328865747);
    MCLOSE("AngerJ[4.5, 0.6]", 0.062217256819454954);
    /* A negative order at small z: the Laplace-integral route's one-panel
     * rule was 1e-3 off here (and 1e-2 at (-0.3, 0.375)) while claiming 2e-16,
     * so the regime choice took it over the correct series.  The quadrature
     * error is now measured (one panel against two). */
    MCLOSE("AngerJ[-5/6 + 0., 5/6 + 0.]", -0.24614064761204724);
    MCLOSE("AngerJ[-0.3, 0.375]", 0.72333831342879523);
    MCLOSE("AngerJ[5/6 + 0., -5/6 + 0.]", -0.24614064761204724);
    assert_true("Abs[AngerJ[-5/6 + 0., 5/6 + 0.] - N[AngerJ[-5/6, 5/6], 30]] < 10^-14");
    MCLOSE("WeberE[1/3, 2.0]", -0.56546195207908119);
    MCLOSE("WeberE[2.5, 7.5]", -0.047184647791081139);
    MCLOSE("WeberE[-1.3, 15.0]", 0.11004674332791673);
    MCLOSE("WeberE[0.3, -4.0]", -0.080388923293129782);
    MCLOSE("WeberE[4.5, 0.6]", 0.082647392517346914);
    MCLOSE("KelvinBer[0, 1.5]", 0.92107218354625576);
    MCLOSE("KelvinBer[1.5, 7.0]", 14.19075299948257);
    /* x = 25: the rotated argument costs ~x/sqrt2 ulp of the envelope, and
     * |bei_2(25)| = 16 |ber_2(25)| */
    mclose("KelvinBer[2, 25.0]", -218564.01686031946, 3.0e-14);
    MCLOSE("KelvinBer[-0.5, 3.0]", 1.1009689024492468);
    MCLOSE("KelvinBei[0, 1.5]", 0.55756006230308669);
    MCLOSE("KelvinBei[1.5, 7.0]", 13.009374741922093);
    mclose("KelvinBei[2, 25.0]", 3592748.1260080393, 3.0e-14);
    MCLOSE("KelvinBei[-0.5, 3.0]", 1.5596596376212121);
    MCLOSE("KelvinKer[0, 1.5]", 0.052934915487710441);
    MCLOSE("KelvinKer[1.5, 7.0]", 0.00023677610845558903);
    mclose("KelvinKer[2, 25.0]", -4.1492913923345197e-9, 1.0e-13);
    MCLOSE("KelvinKer[-0.5, 3.0]", -0.013633041571314303);
    MCLOSE("KelvinKei[0, 1.5]", -0.3313955623385585);
    MCLOSE("KelvinKei[1.5, 7.0]", -0.0037033144628833505);
    mclose("KelvinKei[2, 25.0]", -3.6956633045041027e-9, 1.0e-13);
    MCLOSE("KelvinKei[-0.5, 3.0]", -0.085662378535217098);
    /* order-0 one-argument form */
    MCLOSE("KelvinBer[1.5]", 0.92107218354625576);
}

/* A machine sweep against the 40-digit arb path, normalised by the local
 * envelope (the functions oscillate through zeros, where only absolute
 * accuracy is attainable). Measured max ~2e-15 (Anger/Weber), ~4.4e-15
 * (Kelvin, x = 31), ~1.2e-15 relative (StruveL). */
void test_machine_sweep() {
    assert_true("Max[Flatten[Table[Abs[StruveL[nu, x] - N[StruveL[SetPrecision[nu, 50], SetPrecision[x, 50]], 40]]"
                "/Abs[N[StruveL[SetPrecision[nu, 50], SetPrecision[x, 50]], 40]], "
                "{nu, {-1.3, 0., 1/3, 1., 2.5}}, {x, {0.1, 0.9, 3.3, 7.7, 14.1, 19.9, 31., 60.}}]]] < 4*10^-15");
    assert_true("Max[Flatten[Table[Module[{a = N[AngerJ[SetPrecision[nu, 50], SetPrecision[x, 50]], 40], "
                "e = N[WeberE[SetPrecision[nu, 50], SetPrecision[x, 50]], 40]}, "
                "Max[Abs[AngerJ[nu, x] - a], Abs[WeberE[nu, x] - e]]/Sqrt[a^2 + e^2]], "
                "{nu, {-2.3, -0.5, 1/3, 1.5, 3.7}}, {x, {0.2, 1.1, 3.9, 8.8, 12.5, 18.6, -6.1}}]]] < 1*10^-14");
    assert_true("Max[Flatten[Table[Module[{r = N[KelvinBer[SetPrecision[nu, 50], SetPrecision[x, 50]], 40], "
                "i = N[KelvinBei[SetPrecision[nu, 50], SetPrecision[x, 50]], 40], "
                "kr = N[KelvinKer[SetPrecision[nu, 50], SetPrecision[x, 50]], 40], "
                "ki = N[KelvinKei[SetPrecision[nu, 50], SetPrecision[x, 50]], 40]}, "
                "Max[Max[Abs[KelvinBer[nu, x] - r], Abs[KelvinBei[nu, x] - i]]/Sqrt[r^2 + i^2], "
                "Max[Abs[KelvinKer[nu, x] - kr], Abs[KelvinKei[nu, x] - ki]]/Sqrt[kr^2 + ki^2]]], "
                "{nu, {-1.5, 0., 1/3, 1., 3.7}}, {x, {0.2, 1.1, 3.9, 8.8, 12.5, 18.6, 30.}}]]] < 1*10^-14");
}

/* ---- arbitrary precision and complex arguments ------------------------ */

void test_arbitrary_precision() {
    assert_true("Abs[N[StruveL[1/3, 7], 30] - 166.92846676891390012070671440806`32] < 10^-26");
    assert_true("Abs[N[AngerJ[1/3, 9], 30] - 0.07438154908011458698938997808829`32] < 10^-28");
    assert_true("Abs[N[WeberE[-5/2, 11/2], 30] - (-0.12411568528451102198027127813851`32)] < 10^-28");
    assert_true("Abs[N[KelvinBer[3/2, 7], 30] - 14.190752999482569998540316260158`32] < 10^-27");
    assert_true("Abs[N[KelvinKei[1/3, 4], 30] - 0.020441696204205762158658028321315`32] < 10^-29");
    assert_true("Abs[N[KelvinKer[2, 1/10], 30] - 0.49902224753375222259164152959315`32] < 10^-29");
    assert_true("Precision[N[KelvinBei[1, 3], 40]] >= 39");
    /* the order-zero StruveL[1/2, 3 I] and complex Anger/Weber vs mpmath */
    assert_true("Abs[N[AngerJ[1/2, 2 + I], 30] - (0.77282566746626662303970476079018`32 - "
                "0.46327089048184728661479414495297`32 I)] < 10^-28");
    assert_true("Abs[N[WeberE[1/3, -3 + 2 I], 30] - (0.32554918308071338406669896627286`32 + "
                "1.595791180387950480834291103415`32 I)] < 10^-28");
    assert_true("Abs[N[StruveL[1/2, 3 I], 30] - (-0.64821022167129365541960110455975`32 + "
                "0.64821022167129365541960110455975`32 I)] < 10^-28");
    /* machine complex input -> machine complex output */
    assert_true("Abs[AngerJ[0.5, 2. + I] - (0.7728256674662668 - 0.46327089048184733 I)] < 10^-14");
    /* real argument on the negative axis */
    assert_true("Abs[StruveL[0.5, -2.] - (-1.558402036629881 I)] < 10^-14");
    assert_true("Abs[KelvinBer[1/2, -2.] - (-0.3072383294073605 I)] < 10^-14");
    assert_true("Abs[KelvinBer[1, -2.] - 0.9970776519264285] < 10^-14");
    /* ker/kei: cut along the negative real axis, Schwarz-reflection symmetric,
     * value ON the cut taken from above (arg z = Pi) */
    assert_true("Abs[N[KelvinKer[1/3, -2 + I/10], 25] - Conjugate[N[KelvinKer[1/3, -2 - I/10], 25]]] < 10^-23");
    assert_true("Abs[N[KelvinKer[1/3, -2 + I/10], 25] - (-0.08855442635366594614287619 - 3.40955090629349418803986840 I)] < 10^-22");
    assert_true("Abs[N[KelvinKer[1/3, -2 + 3 I], 25] - (-10.59069788127108089847565 - 4.89677861808647898173069620 I)] < 10^-21");
    assert_true("Abs[KelvinKer[0, -2.] - (-0.0416645139915095 - 2.3616425858660275 I)] < 10^-13");
    /* ker + i kei = e^(-nu Pi I/2) K_nu(z e^(Pi I/4)) for Re z > 0 (an independent route) */
    assert_true("Abs[N[KelvinKer[2/3, 3 + I], 30] + I N[KelvinKei[2/3, 3 + I], 30] - "
                "N[Exp[-2/3 Pi I/2] BesselK[2/3, (3 + I) Exp[Pi I/4]], 30]] < 10^-27");
    assert_true("Abs[N[KelvinBer[2/3, 3 + I], 30] + I N[KelvinBei[2/3, 3 + I], 30] - "
                "N[BesselJ[2/3, (3 + I) Exp[3 Pi I/4]], 30]] < 10^-27");
}

/* ---- derivatives ------------------------------------------------------ */

void test_derivatives() {
    assert_eval_eq("D[AngerJ[n, z], z]", "1/2 (AngerJ[-1 + n, z] - AngerJ[1 + n, z])", 0);
    assert_eval_eq("D[WeberE[n, z], z]", "1/2 (WeberE[-1 + n, z] - WeberE[1 + n, z])", 0);
    assert_eval_eq("D[AngerJ[n, z], n]", "Derivative[1, 0][AngerJ][n, z]", 0);
    assert_eval_eq("D[KelvinBer[n, y], x]", "0", 0);
    /* each rule vs a central difference of the 40-digit function */
    const char* heads[] = { "StruveL", "AngerJ", "WeberE", "KelvinBer", "KelvinBei",
                            "KelvinKer", "KelvinKei" };
    char buf[512];
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        snprintf(buf, sizeof buf,
                 "Module[{h = 10^-12, d = D[%s[7/5, x], x] /. x -> 23/10}, "
                 "Abs[N[d, 30] - N[(%s[7/5, 23/10 + h] - %s[7/5, 23/10 - h])/(2 h), 40]] < 10^-15]",
                 heads[i], heads[i], heads[i]);
        assert_true(buf);
    }
    assert_true("Module[{h = 10^-12, d = D[KelvinKei[x], x] /. x -> 23/10}, "
                "Abs[N[d, 30] - N[(KelvinKei[0, 23/10 + h] - KelvinKei[0, 23/10 - h])/(2 h), 40]] < 10^-15]");
}

/* ---- the surfaces agree bit for bit ----------------------------------- */

void test_surfaces() {
    const char* heads[] = { "StruveL", "AngerJ", "WeberE", "KelvinBer", "KelvinBei",
                            "KelvinKer", "KelvinKei" };
    char buf[1024];
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        const char* h = heads[i];
        snprintf(buf, sizeof buf,
                 "Module[{xs = {0.3, 1.7, 4.2, 9.9, 15.5, 33.}, sc, cf, cfa},"
                 " sc = %s[1/3, #] & /@ xs;"
                 " cf = Compile[{{x, _Real}}, %s[1/3, x]];"
                 " cfa = Compile[{{v, _Real, 1}}, %s[1/3, v]];"
                 " sc === %s[1/3, Developer`FromPackedArray[xs]] &&"
                 " sc === %s[1/3, Developer`ToPackedArray[xs]] &&"
                 " sc === Normal[%s[1/3, NDArray[xs]]] &&"
                 " sc === (cf /@ xs) && sc === cfa[xs]]",
                 h, h, h, h, h, h);
        assert_true(buf);
        snprintf(buf, sizeof buf,
                 "Lookup[CompileDiagnostics[{{x, _Real}}, %s[1, x]], \"Compiled\"] &&"
                 " Lookup[CompileDiagnostics[{{v, _Real, 1}}, %s[1, v]], \"Compiled\"]", h, h);
        assert_true(buf);
    }
    assert_true("Normal[KelvinBer[NDArray[{0.3, 4.2}]]] === {KelvinBer[0, 0.3], KelvinBer[0, 4.2]}");
    assert_true("Lookup[CompileDiagnostics[{{x, _Real}}, KelvinKer[x]], \"Compiled\"]");
    assert_true("Developer`PackedArrayQ[StruveL[1, Developer`ToPackedArray[N[Range[100]/7]]]]");
}

void test_attributes() {
    const char* heads[] = { "StruveL", "AngerJ", "WeberE", "KelvinBer", "KelvinBei",
                            "KelvinKer", "KelvinKei" };
    char buf[256];
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        snprintf(buf, sizeof buf, "Attributes[%s]", heads[i]);
        assert_eval_eq(buf, "{Listable, NumericFunction, Protected}", 0);
    }
    assert_eval_eq("StruveL[1, {0, x}]", "{0, StruveL[1, x]}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_exact_forms);
    TEST(test_exact_forms_numeric);
    TEST(test_machine_values);
    TEST(test_machine_sweep);
    TEST(test_arbitrary_precision);
    TEST(test_derivatives);
    TEST(test_surfaces);
    TEST(test_attributes);

    printf("All StruveL / AngerJ / WeberE / Kelvin tests passed.\n");
    return 0;
}
