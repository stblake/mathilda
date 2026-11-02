/* test_mixed_numericalize.c -- a numeric builtin given an inexact argument
 * numericalizes its exact numeric ones (v0.477): the evaluator's mixed-argument
 * rule (numeric_mixed_args_call), and the special-function numeric paths the
 * sweep found missing -- Surd and KroneckerDelta (new), Max / Min over exact
 * constants, multi-argument HeavisideTheta / DiracDelta, StieltjesGamma[n, a]
 * and N[StieltjesGamma[n]], ExpIntegralE and PolyGamma of non-integer order,
 * LegendreP / LegendreQ of general degree and order, three-argument
 * FactorialPower, QPochhammer of non-integer n, BesselJZero and Nielsen's
 * PolyLog[n, p, z] -- plus exact rational comparison with big parts.
 *
 * Every literal expectation below was verified against the reference
 * implementation, not recalled. Machine values are compared with a tolerance;
 * the exact ones as printed. */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("mixed numericalize test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* |lhs - rhs| < tol, written as one input so the reference value is a literal. */
static void check_near(const char* lhs, const char* rhs, const char* tol) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Abs[(%s) - (%s)] < %s", lhs, rhs, tol);
    check(buf, "True");
}

/* ---- The evaluator rule: an irrational exact argument beside an inexact one ---- */
static void test_mixed_rule(void) {
    check_near("Zeta[2., Sqrt[2]]", "1.0114791294395424", "1*^-14");
    check_near("BesselJ[2., Sqrt[2]]", "0.2108524773255237", "1*^-14");
    check_near("BesselJ[Sqrt[2], 2.]", "0.5120936939951674", "1*^-14");
    check_near("Gamma[2., Sqrt[2]]", "0.5869357175109379", "1*^-14");
    check_near("PolyLog[2., Sqrt[2]/3]", "0.5430889214608642", "1*^-14");
    check_near("Hypergeometric2F1[1., 2, 3, Sqrt[2]/3]", "1.4951457759878553", "1*^-14");
    check_near("ArcTan[2., Sqrt[2]]", "0.6154797086703874", "1*^-14");
    check_near("Beta[2., Pi]", "0.07685687917856682", "1*^-15");
    check_near("LegendreP[2., E/3]", "0.7315093498217748", "1*^-14");
    check_near("Zeta[2., I]", "-0.5369999033772362 - 0.7942335427593189 I", "1*^-14");
    check_near("Binomial[2., Sqrt[2]]", "1.7882762940810417", "1*^-14");
    /* the lowest precision governs: 30 digits beside an exact constant */
    check("Precision[Zeta[2.`30, Sqrt[2]]] > 29", "True");
    check_near("Zeta[2.`30, Sqrt[2]]", "1.0114791294395424054028105339511`30", "1*^-28");
    /* a symbolic argument, or a decline with a message, stays as it was */
    check("Zeta[2., x]", "Zeta[2.0, x]");
    check("Zeta[2., Sqrt[2] x]", "Zeta[2.0, Sqrt[2] x]");
    check("ProductLog[1., Sqrt[2]]", "ProductLog[1.0, Sqrt[2]]");
    /* a user NumericFunction head is not a builtin: unevaluated, as in the reference */
    check("SetAttributes[mixg, NumericFunction]; mixg[2., Sqrt[2]]", "mixg[2.0, Sqrt[2]]");
}

/* ---- Max / Min order exact constants with the numbers ---- */
static void test_max_min(void) {
    check("{Max[1, Pi], Max[Pi, E], Min[1/2, Sqrt[2]], Max[1., Pi], Min[1., Pi]}",
          "{Pi, Pi, 1/2, Pi, 1.0}");
    check("{Max[1., Pi, x], Max[{1., Pi}], Max[2., 2]}", "{Max[Pi, x], Pi, 2.0}");
}

/* ---- Surd ---- */
static void test_surd(void) {
    check("{Surd[8, 3], Surd[-8, 3], Surd[2, 3], Surd[16, 4], Surd[-32, 5], Surd[8, 1], Surd[8, -3], Surd[-27, -3]}",
          "{2, -2, 2^(1/3), 2, -2, 8, 1/2, -1/3}");
    check("{Surd[-2, 3], Surd[-8/27, 3], Surd[2/3, 2], Surd[-Pi, 3], Surd[Sqrt[2], 3], Surd[2^(1/3), 3], Surd[E^2, 2]}",
          "{-2^(1/3), -2/3, Sqrt[2/3], -Pi^(1/3), 2^(1/6), 2^(1/9), E}");
    check("{Surd[0, 3], Surd[x, 3], Surd[8, k], Surd[Infinity, 2], Surd[-Infinity, 3]}",
          "{0, Surd[x, 3], Surd[8, k], Infinity, -Infinity}");
    check("Quiet[{Surd[-4, 2], Surd[8, 0], Surd[I, 3], Surd[8, 3/2]}]",
          "{Indeterminate, Indeterminate, Surd[I, 3], Surd[8, 3/2]}");
    check("Check[Surd[-4, 2], m, Surd::noneg]", "m");
    check("Check[Surd[8, 0], m, Surd::indet]", "m");
    check("Check[Surd[I, 3], m, Surd::preal]", "m");
    check("Check[Surd[8, 3/2], m, Surd::int]", "m");
    check("Check[Surd[8], m, Surd::argr]", "m");
    check("Check[Surd[1, 2, 3], m, Surd::argrx]", "m");
    check_near("Surd[-2.5, 3]", "-1.3572088082974532", "1*^-15");
    check_near("Surd[2.5, -2]", "0.6324555320336759", "1*^-15");
    check_near("Surd[2., 7]", "1.1040895136738123", "1*^-15");
    check("Surd[-8., 3]", "-2.0");
    check("Surd[{8, 27}, 3]", "{2, 3}");
    check("Attributes[Surd]", "{Listable, NumericFunction, Protected}");
    /* scalar == packed == NDArray == Compile */
    check("Surd[N[{-8, 2, 27}], 3] == Normal[Surd[NDArray[{-8., 2., 27.}], 3]] == "
          "Compile[{{v, _Real, 1}}, Surd[v, 3]][{-8., 2., 27.}] == {-2., Surd[2., 3], 3.}",
          "True");
}

/* ---- KroneckerDelta ---- */
static void test_kroneckerdelta(void) {
    check("{KroneckerDelta[], KroneckerDelta[0], KroneckerDelta[2], KroneckerDelta[x, x], "
          "KroneckerDelta[1, 1, 1], KroneckerDelta[1, 1, 2], KroneckerDelta[1, 1., 1]}",
          "{1, 1, 0, 1, 1, 0, 1}");
    check("{KroneckerDelta[Pi, Pi], KroneckerDelta[Pi, 3], KroneckerDelta[x, 1, 2], "
          "KroneckerDelta[I, I], KroneckerDelta[I, 1], KroneckerDelta[0.], "
          "KroneckerDelta[1/2, 0.5], KroneckerDelta[1, 1 + 10^-20]}",
          "{1, 0, 0, 1, 0, 1, 1, 0}");
    check("{KroneckerDelta[x], KroneckerDelta[x, y], KroneckerDelta[1, x], KroneckerDelta[x, x, y]}",
          "{KroneckerDelta[x], KroneckerDelta[x, y], KroneckerDelta[1, x], KroneckerDelta[x, x, y]}");
    check("Attributes[KroneckerDelta]", "{NumericFunction, Orderless, Protected}");
}

/* ---- Exact rationals with big parts compare exactly ---- */
static void test_rational_compare(void) {
    check("{1 == 1 + 10^-20, 1 < 1 + 10^-20, 1/3 == 1/3 + 10^-30}", "{False, True, False}");
}

/* ---- HeavisideTheta / DiracDelta with several arguments ---- */
static void test_distributions(void) {
    /* The HeavisideTheta / DiracDelta rules live in src/internal/distributions.m,
     * which the REPL loads at startup and this harness does not. The suite
     * runs from tests/ or tests/build/, so search upward as test_load_init_m. */
    check("Get[SelectFirst[{\"src\", \"../src\", \"../../src\", \"../../../src\"}, "
          "FileExistsQ[# <> \"/internal/distributions.m\"] &] <> "
          "\"/internal/distributions.m\"]; 0", "0");
    check("{HeavisideTheta[1, 2], HeavisideTheta[1, -2, x], HeavisideTheta[1., Sqrt[2]], "
          "HeavisideTheta[1, x], DiracDelta[1, x], DiracDelta[0., Pi], DiracDelta[x, y]}",
          "{1, 0, 1, HeavisideTheta[x], 0, 0, DiracDelta[x, y]}");
}

/* ---- StieltjesGamma ---- */
static void test_stieltjes(void) {
    check_near("StieltjesGamma[1, Sqrt[2]/3 + 0.]", "-1.560212950974281", "1*^-14");
    check_near("StieltjesGamma[1, 2.]", "-0.07281584548367671", "1*^-15");
    check_near("StieltjesGamma[0, 0.5]", "1.9635100260214235", "1*^-14");
    check_near("StieltjesGamma[1, -0.5]", "0.0328346803149491 - 6.283185307179586 I", "1*^-13");
    check_near("N[StieltjesGamma[3]]", "0.002053834420303346", "1*^-17");
    check_near("N[StieltjesGamma[3], 30]", "0.0020538344203033458661600465427533842857`30", "1*^-31");
    check_near("N[StieltjesGamma[2, 1/3]]", "3.6191639096189636", "1*^-14");
    check("{StieltjesGamma[1, 1], StieltjesGamma[1, 0.], StieltjesGamma[1, -1.], StieltjesGamma[2, 1/3]}",
          "{StieltjesGamma[1], ComplexInfinity, ComplexInfinity, StieltjesGamma[2, 1/3]}");
    check("Check[StieltjesGamma[-1], m, StieltjesGamma::intnm]", "m");
    check("Attributes[StieltjesGamma]", "{Listable, NHoldFirst, Protected}");
}

/* ---- ExpIntegralE / PolyGamma of non-integer order ---- */
static void test_order(void) {
    check_near("ExpIntegralE[Pi/3, 0.5]", "0.5432571240019218", "1*^-14");
    check_near("ExpIntegralE[0.5, -2.]", "-4.728907785610418 - 1.2533141373155003 I", "1*^-13");
    check_near("ExpIntegralE[1.5, 2. + I]", "0.009807468767872408 - 0.0388767535626179 I", "1*^-15");
    check_near("N[ExpIntegralE[1/3, 2], 25]", "0.06024889885134620552374433379765667381`25", "1*^-26");
    check("{ExpIntegralE[0.5, 0.], ExpIntegralE[1.5, 0.]}", "{ComplexInfinity, 2.0}");
    check_near("PolyGamma[0.5, 2.]", "0.7552610839464955", "1*^-14");
    check_near("PolyGamma[1/3, 2.]", "0.6769989218804935", "1*^-14");
    check_near("PolyGamma[-0.5, 2.]", "0.05438983957687823", "1*^-15");
    check_near("PolyGamma[0.5, -1.5]", "-2.989029111003179 + 3.220842847032668 I", "1*^-13");
    check_near("PolyGamma[0.5, 2. + I]", "0.8107022601528389 + 0.030985718114095007 I", "1*^-14");
    check_near("N[PolyGamma[1/2, 3], 25]", "0.77864761243032796861665855276810765594`25", "1*^-25");
    check_near("PolyGamma[2.5, 0.5`30]", "-31.027597757423173", "1*^-12");
}

/* ---- LegendreP / LegendreQ of general degree and order ---- */
static void test_legendre(void) {
    check_near("LegendreP[Pi/3, E/2, Sqrt[2]/3 + 0.]", "-1.2929259104618267", "1*^-14");
    check_near("LegendreQ[Pi/3, E/2, Sqrt[2]/3 + 0.]", "0.28032603167145926", "1*^-14");
    check_near("LegendreP[Pi/3, E/2, 3, Sqrt[2]/3 + 0.]",
               "0.6913105860002337 + 1.0925873345538822 I", "1*^-13");
    check_near("LegendreQ[Pi/3, E/2, 3, Sqrt[2]/3 + 0.]",
               "0.03187246470814304 + 2.0499307701966925 I", "1*^-13");
    check_near("N[LegendreP[Pi/3, E/2, 2, Sqrt[2]/3], 20]",
               "-1.29292591046182668792704533053032647075`20", "1*^-18");
    check_near("LegendreP[Pi/3, E/2, 1.5]", "-0.2156376105157518 - 0.3408061830300018 I", "1*^-13");
    check_near("LegendreQ[Pi/3, E/2, 1.5]", "-0.880562352942063 + 0.8843367728358586 I", "1*^-13");
    check_near("LegendreQ[Pi/3, 1.5]", "0.19537179424877609 - 2.4195355407695165 I", "1*^-13");
    check_near("LegendreP[Pi/3, E/2, 2. + I]", "-1.600215433155711 + 0.525039888422923 I", "1*^-13");
}

/* ---- FactorialPower with a step ---- */
static void test_factorialpower(void) {
    check("{FactorialPower[5, 3, 2], FactorialPower[5, -2], FactorialPower[5, -2, 2], "
          "FactorialPower[5, 3, 1], FactorialPower[7/2, 2, 1/2], FactorialPower[5, 0, h]}",
          "{15, 1/42, 1/63, 60, 21/2, 1}");
    check("FactorialPower[5, 1/2, 2] == 15 Sqrt[Pi/2]/8", "True");
    check("FactorialPower[x, 3, h]", "FactorialPower[x, 3, h]");
    check_near("FactorialPower[5., 0.5, 2]", "2.3499640074665633", "1*^-14");
    check_near("FactorialPower[Sqrt[5], 0.47, 0.866]", "1.5310964069693207", "1*^-13");
    check("{FactorialPower[2.5, 3, 0.5], FactorialPower[-2.5, 1.5, 0.5]}", "{7.5, ComplexInfinity}");
}

/* ---- QPochhammer ---- */
static void test_qpochhammer(void) {
    check_near("QPochhammer[0.5, 0.3]", "0.39808220430187763", "1*^-15");
    check_near("QPochhammer[0.5`30, 0.3`30]", "0.3980822043018776635608083355611187102950527623308`30", "1*^-28");
    check_near("QPochhammer[0.5 + I, 0.3]", "-0.008136743048371628 - 0.9484904560912539 I", "1*^-14");
    check_near("QPochhammer[0.5, 0.3, 2.5]", "0.41248848035749924", "1*^-14");
    check_near("QPochhammer[0.5, 0.3, 2.5 + I]", "0.4028188271856958 - 0.013333857685628588 I", "1*^-14");
    check_near("N[QPochhammer[1/2, 1/3], 25]", "0.38266319667903302328895503353319130868`25", "1*^-24");
    check("{QPochhammer[1/2, 1/3, -2], QPochhammer[2., 0.5], QPochhammer[0.5, 1.5]}",
          "{4/7, 0.0, QPochhammer[0.5, 1.5]}");
    check_near("QPochhammer[0.5, 0.3, -2]", "0.3292682926829269", "1*^-15");
}

/* ---- BesselJZero ---- */
static void test_besseljzero(void) {
    check_near("BesselJZero[1.0471975511965976, 2]", "7.08345235356954", "1*^-13");
    check_near("N[BesselJZero[1, 2]]", "7.015586669815622", "1*^-13");
    check_near("N[BesselJZero[0, 1], 30]", "2.4048255576957727686216318793264541313541553727`30", "1*^-28");
    check_near("BesselJZero[2.5, 100.]", "317.29140298173223", "1*^-11");
    check_near("BesselJZero[-0.5, 1.]", "Pi/2", "1*^-14");
    check_near("BesselJZero[50., 1]", "57.116899160119175", "1*^-12");
    check_near("BesselJZero[50., 3]", "67.69740841076478", "1*^-12");
    check_near("BesselJZero[0., 1000]", "3140.8072952250786", "1*^-10");
    check_near("BesselJZero[0.3, 25]", "78.22667965763888", "1*^-12");
    check("BesselJZero[1, 2]", "BesselJZero[1, 2]");
    check("BesselJZero[N[Range[3]], 2] == Compile[{{v, _Real, 1}}, BesselJZero[v, 2]][N[Range[3]]]", "True");
}

/* ---- Nielsen's generalized polylogarithm ---- */
static void test_nielsen(void) {
    check_near("PolyLog[2, 3, 0.4714045207910317]", "0.003132304552819349", "1*^-16");
    check_near("PolyLog[2, 1, 0.5]", "0.5372131936080403", "1*^-15");
    check_near("PolyLog[2.5, 2, -0.7]", "0.03484991408142621", "1*^-15");
    check_near("PolyLog[2, 2, 0.5 + 0.3 I]", "0.0160343247574304 + 0.05148738174866205 I", "1*^-15");
    check_near("PolyLog[1, 2, 0.95]", "0.7698743359411672", "1*^-14");
    check_near("PolyLog[Pi/3, 3, Sqrt[2]/3 + 0.]", "0.010572480514274287", "1*^-15");
    check_near("N[PolyLog[2, 3, 1/3], 30]", "0.0009427527623780048874763556343505861683936`30", "1*^-32");
    check_near("PolyLog[2, 2, 0.5`25]", "0.04075823915930925192075921743360984034`25", "1*^-25");
    check("PolyLog[2, 0, 0.5]", "PolyLog[2, 0, 0.5]");
}

/* ---- Nielsen beyond the unit disk and for a non-integer p (v0.478) ---- */
static void test_nielsen_continued(void) {
    /* |z| > 1, integer n >= 1 and p: the principal-logarithm integral */
    check_near("PolyLog[2, 2, 2.]", "0.4554289599591583 - 1.3189461296972319 I", "1*^-14");
    check_near("PolyLog[2, 2, -3.]", "0.5052321953754769", "1*^-15");
    check_near("PolyLog[2, 2, 2. + I]", "-0.46179573693753576 + 0.8126908557655475 I", "1*^-14");
    check_near("PolyLog[2, 2, 0.5 + 3 I]", "-0.6322864658786631 - 0.43026879618230067 I", "1*^-14");
    check_near("PolyLog[2, 3, 5.]", "-1.4264256992301926 + 2.158171317830972 I", "1*^-13");
    check_near("PolyLog[1, 2, 2.]", "-1.3169446513992683 - 2.583856390024985 I", "1*^-14");
    check_near("PolyLog[3, 2, 2.]", "0.4642952438592966 - 0.36813325559494618 I", "1*^-14");
    check_near("PolyLog[2, 2, 100.]", "-25.953340869693107 + 31.083771542159424 I", "1*^-12");
    check_near("PolyLog[2, 2, -1.]", "0.0877856715686553", "1*^-15");
    /* z = 1, and a non-integer or complex p in the closed disk */
    check_near("PolyLog[2, 2, 1.]", "0.27058080842778455", "1*^-15");
    check_near("PolyLog[2, 2.5, 1.]", "0.15686179900608218", "1*^-15");
    check_near("PolyLog[2, 0.5, 0.5]", "3.2411607336981399", "1*^-14");
    check_near("PolyLog[1.5, 1.5, 0.5]", "0.17858105367126854", "1*^-15");
    check_near("PolyLog[2, 1.5, -0.5]", "-0.10483906257910724 I", "1*^-15");
    check_near("PolyLog[2, 1. + I, 0.5]", "-0.2918380219781763 - 0.1824400742846398 I", "1*^-14");
    /* 30-digit quadrature; the reference's machine value is 1.3e-6 off */
    check_near("PolyLog[3, 0.5, 0.95]", "8.8877201600060307", "1*^-13");
    /* integer n <= 0: (z d/dz)^-n (-Log[1 - z])^p / Gamma[p + 1], any p and z */
    check_near("PolyLog[-1, 2, 2.]", "2 Pi I", "1*^-14");
    check_near("PolyLog[0, 2, 3.]", "-4.694575693585578 + 2.177586090303602 I", "1*^-14");
    check_near("PolyLog[-1, 2, 0.5]", "Log[2.]", "1*^-15");
    check_near("PolyLog[-1, 0.5, 1.5]", "-0.7356391804312229 - 0.5910237088861859 I", "1*^-13");
    /* arbitrary precision */
    check_near("N[PolyLog[2, 2, 3], 30]",
               "-0.48051063523149439645020401492659433081`30 - 2.30253895026584497656797406915324175659`30 I",
               "1*^-28");
    check_near("N[PolyLog[3/2, 5/2, 9/10], 40]", "0.1846209219293389738431779980488161567043718`40", "1*^-38");
    check_near("N[PolyLog[3, 2, 1], 30]", "0.0965511599894437344656455314289427640320`30", "1*^-29");
    /* deliberate difference: the reference says Indeterminate at machine 1. */
    check_near("PolyLog[1, 2, 1.]", "Zeta[3]", "1*^-14");
    /* left unevaluated, as in the reference: non-integer n or p off the disk, p <= 0 */
    check("{PolyLog[2, 2.5, 2.], PolyLog[2.5, 2, 3.], PolyLog[2, -0.5, 0.5], PolyLog[2, 2. + I, 3.]}",
          "{PolyLog[2, 2.5, 2.0], PolyLog[2.5, 2, 3.0], PolyLog[2, -0.5, 0.5], PolyLog[2, 2.0 + 1.0*I, 3.0]}");
}

/* ---- n-ary kernels over arrays and in Compile (v0.479) ---- */
static void test_nary_surfaces(void) {
    check("nzs = ToPackedArray[N[Range[-400, 99]/100]]; nns = ToPackedArray[N[Range[1, 40]/8]]; 0", "0");
    /* Nielsen PolyLog: packed == boxed == NDArray == Compile, the array anywhere */
    check("PolyLog[2, 2, nzs] === Map[PolyLog[2, 2, #] &, FromPackedArray[nzs]]", "True");
    check("PackedArrayQ[PolyLog[2, 2, nzs]]", "True");
    check("Normal[PolyLog[2, 2, NDArray[nzs]]] === PolyLog[2, 2, nzs]", "True");
    check("Compile[{{w, _Real, 1}}, PolyLog[2, 2, w]][nzs] === PolyLog[2, 2, nzs]", "True");
    check("Compile[{{w, _Real, 1}}, PolyLog[w, 2, 0.5]][nns] === Map[PolyLog[#, 2, 0.5] &, FromPackedArray[nns]]", "True");
    check("Compile[{{w, _Real, 1}}, PolyLog[2.5, w, 0.7]][nns] === Map[PolyLog[2.5, #, 0.7] &, FromPackedArray[nns]]", "True");
    check("Table[Compile[{{x, _Real}}, PolyLog[2, 3, x]][x], {x, -3., 0.9, 0.1}] === Table[PolyLog[2, 3, x], {x, -3., 0.9, 0.1}]", "True");
    check("Lookup[CompileDiagnostics[{{x, _Real}}, PolyLog[2, 2, x]], \"Compiled\"]", "True");
    check("Lookup[CompileDiagnostics[{{v, _Real, 1}}, PolyLog[2, 2, v]], \"Compiled\"]", "True");
    /* z > 1 is complex: the compiled real body hands it to the interpreter */
    check_near("Compile[{{w, _Real, 1}}, PolyLog[2, 2, w]][{2.}][[1]]", "0.4554289599591583 - 1.3189461296972319 I", "1*^-14");
    check("PolyLog[2, 2, {1/2, 1/3}]", "{PolyLog[2, 2, 1/2], PolyLog[2, 2, 1/3]}");
    /* LerchPhi and 2F1 with the array in a leading position */
    check("Lookup[CompileDiagnostics[{{v, _Real, 1}}, LerchPhi[v, 2., 1.5]], \"Compiled\"]", "True");
    check("Compile[{{w, _Real, 1}}, LerchPhi[w, 2., 1.5]][nns/6] === Map[LerchPhi[#, 2., 1.5] &, FromPackedArray[nns/6]]", "True");
    check("Compile[{{w, _Real, 1}}, LerchPhi[0.5, w + 1, 1.5]][nns] === Map[LerchPhi[0.5, # + 1, 1.5] &, FromPackedArray[nns]]", "True");
    check("Compile[{{w, _Real, 1}}, Hypergeometric2F1[w, 1., 1.5, 0.3]][nns] === Map[Hypergeometric2F1[#, 1., 1.5, 0.3] &, FromPackedArray[nns]]", "True");
    check("Compile[{{w, _Real, 1}}, WhittakerM[w, 0.25, 1.3]][nns] === Map[WhittakerM[#, 0.25, 1.3] &, FromPackedArray[nns]]", "True");
    /* the packed map for n-ary Listable heads */
    check("LerchPhi[nns/6, 2., 1.5] === Map[LerchPhi[#, 2., 1.5] &, FromPackedArray[nns/6]] && PackedArrayQ[LerchPhi[nns/6, 2., 1.5]]", "True");
    check("Hypergeometric2F1[0.5, 1., 1.5, nns/6] === Map[Hypergeometric2F1[0.5, 1., 1.5, #] &, FromPackedArray[nns/6]]", "True");
    check("LerchPhi[Range[3], 2, 1/2] // Length", "3");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_mixed_rule);
    TEST(test_max_min);
    TEST(test_surd);
    TEST(test_kroneckerdelta);
    TEST(test_rational_compare);
    TEST(test_distributions);
    TEST(test_stieltjes);
    TEST(test_order);
    TEST(test_legendre);
    TEST(test_factorialpower);
    TEST(test_qpochhammer);
    TEST(test_besseljzero);
    TEST(test_nielsen);
    TEST(test_nielsen_continued);
    TEST(test_nary_surfaces);

    printf("All mixed numericalize tests passed!\n");
    symtab_clear();
    return 0;
}
