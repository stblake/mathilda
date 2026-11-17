/* Tests for the machine kernels of the exponential-integral family
 * (src/special_functions/sf_machine_expint.c): SinIntegral, CosIntegral,
 * SinhIntegral, CoshIntegral, ExpIntegralEi, LogIntegral and ExpIntegralE.
 *
 * Covers: real spot values across every regime boundary of the kernels
 * (x = 2, 4, 8, 16, 64 for Ei; 0.6 / 2 / 8 for E1; 2 / 4 / 8 for Si/Ci; the
 * zeros of Ci, Chi, Ei and li, where the tolerance is scaled by |x f'(x)|);
 * complex values in both the series and the continued-fraction regions and on
 * the cuts; the declines (poles, overflow, underflow) that fall through to the
 * MPFR path; machine-in => machine-out; exact input staying exact; and the
 * five-surface bit-match scalar == unpacked List == packed List == NDArray ==
 * Compile, real and complex.
 *
 * Reference values: mpmath 1.3 at 30 digits for the exact double argument,
 * computed for this file (not recalled). Tolerances are in units of 2^-52
 * relative to the reference.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "arithmetic.h"
#include "special_functions/sf_machine.h"
#include "special_functions/sf_machine_expint.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Expr* ev(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static void assert_real(const char* input, double expected, double ulps) {
    Expr* r = ev(input);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a machine Real", input);
    double v = r->data.real;
    expr_free(r);
    ASSERT_MSG(fabs(v - expected) <= ulps * 2.220446049250313e-16 * fabs(expected),
               "%s: expected %.17g, got %.17g", input, expected, v);
}

static void assert_cplx(const char* input, double er, double ei, double ulps) {
    Expr* r = ev(input);
    Expr *re, *im;
    ASSERT_MSG(is_complex(r, &re, &im), "%s: expected Complex", input);
    ASSERT_MSG(re->type == EXPR_REAL && im->type == EXPR_REAL,
               "%s: expected MACHINE Complex parts", input);
    double dr = re->data.real - er, di = im->data.real - ei;
    double err = sqrt(dr * dr + di * di), mag = sqrt(er * er + ei * ei);
    double gr = re->data.real, gi = im->data.real;
    expr_free(r);
    ASSERT_MSG(err <= ulps * 2.220446049250313e-16 * mag,
               "%s: expected %.17g + %.17g I, got %.17g + %.17g I", input, er, ei, gr, gi);
}

static void assert_true(const char* input) {
    Expr* r = ev(input);
    bool ok = r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    expr_free(r);
    ASSERT_MSG(ok, "%s: expected True", input);
}

/* ---- real kernels -------------------------------------------------------- */

void test_expint_real_values() {
    assert_real("SinIntegral[0.5`]", 0.49310741804306668916, 4);
    assert_real("SinIntegral[3.9`]", 1.7765013604478054544, 4);
    assert_real("SinIntegral[4.0001`]", 1.7581842183061578669, 4);
    assert_real("SinIntegral[7.5`]", 1.5106815309433858782, 4);
    assert_real("SinIntegral[12.25`]", 1.4961780369236550291, 4);
    assert_real("SinIntegral[100.5`]", 1.5608558784061177879, 4);
    assert_real("SinIntegral[-2.5`]", -1.7785201734438266421, 4);
    assert_real("SinIntegral[1`*^-8]", 1.0000000000000000154e-8, 4);
    assert_real("CosIntegral[0.5`]", -0.17778407880661290134, 4);
    assert_real("CosIntegral[0.6165`]", -7.2598997999420000701e-6, 449542);
    assert_real("CosIntegral[3.3841804225511862`]", 5.6568522015712362632e-17, 68640266294124254);
    assert_real("CosIntegral[3.9823293187015207`]", -0.13805824385045862682, 4);
    assert_real("CosIntegral[7.5`]", 0.11563320323793427044, 4);
    assert_real("CosIntegral[25.5`]", 0.012615719381297764738, 296);
    assert_real("SinhIntegral[0.5`]", 0.50699674981966719583, 4);
    assert_real("SinhIntegral[3.9`]", 9.1599907458621255826, 4);
    assert_real("SinhIntegral[4.5`]", 13.967885049347147622, 4);
    assert_real("SinhIntegral[15.5`]", 186920.12661790259701, 4);
    assert_real("SinhIntegral[40.25`]", 3852842884896847.7999, 157);
    assert_real("SinhIntegral[-6.5`]", -63.873711831546402314, 4);
    assert_real("CoshIntegral[0.5`]", -0.052776844956493615913, 86);
    assert_real("CoshIntegral[0.5238`]", -0.000049138674235027499647, 92827);
    assert_real("CoshIntegral[3.9`]", 9.1557236005809070108, 4);
    assert_real("CoshIntegral[4.5`]", 13.965811648592433007, 4);
    assert_real("CoshIntegral[15.5`]", 186920.12661789131485, 4);
    assert_real("CoshIntegral[40.25`]", 3852842884896847.7999, 157);
    assert_real("ExpIntegralEi[0.3725`]", -0.00002887418318874596456, 201060);
    assert_real("ExpIntegralEi[1.5`]", 3.301285449129797838, 4);
    assert_real("ExpIntegralEi[2.5`]", 7.0737658945786007119, 4);
    assert_real("ExpIntegralEi[7.0`]", 191.50474333550139595, 4);
    assert_real("ExpIntegralEi[11.0`]", 6071.406374098611508, 4);
    assert_real("ExpIntegralEi[15.0`]", 234955.85249076830358, 4);
    assert_real("ExpIntegralEi[20.0`]", 25615652.66405658882, 4);
    assert_real("ExpIntegralEi[50.0`]", 1.0585636897131690963e+20, 196);
    assert_real("ExpIntegralEi[100.0`]", 2.7155527448538798219e+41, 396);
    assert_real("ExpIntegralEi[-0.3`]", -0.90567665167584673985, 4);
    assert_real("ExpIntegralEi[-1.5`]", -0.1000195824066326519, 4);
    assert_real("ExpIntegralEi[-7.0`]", -0.00011548173161033821643, 4);
    assert_real("ExpIntegralEi[-50.0`]", -3.7832640295504590187e-24, 204);
    assert_real("LogIntegral[0.25`]", -0.11866205644712310531, 4);
    assert_real("LogIntegral[0.9`]", -1.7758006834235251473, 4);
    assert_real("LogIntegral[1.1`]", -1.6757728303191987916, 4);
    assert_real("LogIntegral[2.5`]", 1.6672946675063239514, 4);
    assert_real("LogIntegral[20.0`]", 9.9052999776328211363, 4);
    assert_real("LogIntegral[10000.0`]", 1246.1372158993884597, 4);
    assert_real("LogIntegral[1000000000000.0`]", 37607950280.80486549, 4);
    assert_real("LogIntegral[1`*^30]", 1.4692398897720447927e+28, 4);
    assert_real("ExpIntegralE[2, 0.5`]", 0.32664386232455301773, 3);
    assert_real("ExpIntegralE[2, 1.0`]", 0.14849550677592204792, 3);
    assert_real("ExpIntegralE[2, 1.5`]", 0.07310078653848085108, 3);
    assert_real("ExpIntegralE[2, 2.367`]", 0.023420762139753148454, 3);
    assert_real("ExpIntegralE[2, 5.0`]", 0.00099646904270883810998, 3);
    assert_real("ExpIntegralE[2, 12.0`]", 4.4291416372121703957e-7, 3);
    assert_real("ExpIntegralE[2, 30.0`]", 2.9296693677373697047e-15, 3);
    assert_real("ExpIntegralE[1, 0.5`]", 0.55977359477616081175, 3);
    assert_real("ExpIntegralE[3, 0.25`]", 0.32468412597814364417, 3);
    assert_real("ExpIntegralE[3, 4.5`]", 0.0015524386995614347647, 3);
    assert_real("ExpIntegralE[10, 2.0`]", 0.012092085136400295867, 3);
    assert_real("ExpIntegralE[0, 3.0`]", 0.01659568945595464766, 3);
}

/* ---- complex kernels ----------------------------------------------------- */

void test_expint_complex_values() {
    /* Complex[x, 0.] collapses to the real x, so Ei of it is real: */
    assert_real("ExpIntegralEi[Complex[-5.`, 0.`]]", -0.0011482955912753257973, 4);
    assert_cplx("SinIntegral[Complex[1.5`, 1.0`]]", 1.5324237219775529809, 0.6883092845662321108, 4);
    assert_cplx("SinIntegral[Complex[5.5`, -2.0`]]", 1.0393467931108759644, 0.35483003398150615217, 4);
    assert_cplx("SinIntegral[Complex[-3.0`, 2.5`]]", -3.197655178639338831, -0.49620226984975673707, 4);
    assert_cplx("SinIntegral[Complex[50.0`, 10.0`]]", -214.30061584511994363, -18.939440913649358639, 4);
    assert_cplx("SinIntegral[Complex[0.5`, 20.0`]]", 5838478.6837651588591, 11394582.667845824077, 4);
    assert_cplx("CosIntegral[Complex[1.5`, 1.0`]]", 0.78391765510921827554, -0.081219410315058576689, 4);
    assert_cplx("CosIntegral[Complex[5.5`, 2.0`]]", -0.37524430990617789329, 0.5249586401364557965, 4);
    assert_cplx("CosIntegral[Complex[-3.0`, 2.5`]]", -0.48138403664422270845, 4.7582666769606495817, 4);
    assert_cplx("CosIntegral[Complex[-3.0`, -2.5`]]", -0.48138403664422270845, -4.7582666769606495817, 4);
    assert_cplx("CosIntegral[Complex[-2.0`, 0.0`]]", 0.4229808287748649957, 3.1415926535897932385, 4);
    assert_cplx("CosIntegral[Complex[50.0`, 10.0`]]", -18.939441324476748194, 215.87141138636322563, 4);
    assert_cplx("SinhIntegral[Complex[1.5`, 1.0`]]", 1.4054562505640779872, 1.3242960495975921074, 4);
    assert_cplx("SinhIntegral[Complex[2.0`, 5.5`]]", -0.35483003398150615217, 1.0393467931108759644, 4);
    assert_cplx("SinhIntegral[Complex[-3.0`, 2.5`]]", -0.068633050189996595508, 4.505735874563952267, 4);
    assert_cplx("SinhIntegral[Complex[10.0`, 50.0`]]", -18.939440913649358639, -214.30061584511994363, 4);
    assert_cplx("CoshIntegral[Complex[1.5`, 1.0`]]", 1.394215259191653174, 1.412783171911092401, 4);
    assert_cplx("CoshIntegral[Complex[2.0`, 5.5`]]", -0.37524430990617789329, 1.0458376866584408227, 4);
    assert_cplx("CoshIntegral[Complex[-3.0`, 2.5`]]", 0.079364420690696925877, -1.3647139361426454949, 4);
    assert_cplx("CoshIntegral[Complex[-3.0`, -2.5`]]", 0.079364420690696925877, 1.3647139361426454949, 4);
    assert_cplx("CoshIntegral[Complex[-2.0`, 0.0`]]", 2.4526669226469145219, 3.1415926535897932385, 4);
    assert_cplx("CoshIntegral[Complex[0.0`, 7.0`]]", 0.076695278482184518383, 1.5707963267948966192, 4);
    assert_cplx("ExpIntegralEi[Complex[1.5`, 1.0`]]", 2.7996715097557311612, 2.7370792215086845083, 4);
    assert_cplx("ExpIntegralEi[Complex[2.5`, -2.9`]]", -0.70290429723413769613, -6.5392455620159693101, 4);
    assert_cplx("ExpIntegralEi[Complex[-5.0`, 1.0`]]", -0.00043944980567556961704, 3.1405506204723477589, 4);
    assert_cplx("ExpIntegralEi[Complex[50.0`, 3.0`]]", -1.0349279039845694413e+20, 21280339215336758931.0, 4);
    assert_cplx("LogIntegral[Complex[9.4`, 1.4`]]", 5.9220379168577594618, 0.62287923057486396251, 4);
    assert_cplx("LogIntegral[Complex[-2.0`, 0.0`]]", 0.035532275913560467787, 3.7351504825512546398, 4);
    assert_cplx("LogIntegral[Complex[0.5`, -0.5`]]", -0.013953541766945212027, -2.6296841153703874488, 4);
    assert_cplx("LogIntegral[Complex[1000.0`, 1000.0`]]", 186.27709339553471677, 141.41031081818222198, 4);
    assert_cplx("ExpIntegralE[2, Complex[1.0`, 1.0`]]", 0.019159950855072582158, -0.13051696506573467662, 4);
    assert_cplx("ExpIntegralE[2, Complex[-1.0`, 0.0`]]", 0.82316401210310847989, -3.1415926535897932385, 4);
    assert_cplx("ExpIntegralE[3, Complex[-5.0`, 0.5`]]", -59.96874569000559025, -13.233756264867259102, 4);
    assert_cplx("ExpIntegralE[2, Complex[0.3`, -0.4`]]", 0.32027644698263528721, 0.28411272283597152139, 4);
    assert_cplx("ExpIntegralE[1, Complex[-2.0`, 0.0`]]", -4.9542343560018901634, -3.1415926535897932385, 4);
    assert_cplx("ExpIntegralE[5, Complex[10.0`, 20.0`]]", -8.9172881821120087109e-7, -1.5781765989324117111e-6, 4);
}

/* ---- declines, exact input, special values ------------------------------- */

void test_expint_special() {
    double v, vi;
    /* poles and complex-valued real arguments decline in the REAL kernels */
    ASSERT(!sf_machine_ci(0.0, &v));
    ASSERT(!sf_machine_ci(-1.0, &v));
    ASSERT(!sf_machine_chi(-1.0, &v));
    ASSERT(!sf_machine_ei(0.0, &v));
    ASSERT(!sf_machine_li(1.0, &v));
    ASSERT(sf_machine_li(0.0, &v) && v == 0.0);   /* li(0) = 0 */
    ASSERT(!sf_machine_li(-2.0, &v));
    ASSERT(!sf_machine_e1(0.0, &v));
    ASSERT(!sf_machine_si(NAN, &v));
    /* overflow and underflow decline, so the MPFR path answers them */
    ASSERT(!sf_machine_ei(800.0, &v));
    ASSERT(!sf_machine_ei(-800.0, &v));
    ASSERT(!sf_machine_shi(800.0, &v));
    ASSERT(!sf_machine_ei_c(0.0, 0.0, &v, &vi));
    ASSERT(!sf_machine_li_c(1.0, 0.0, &v, &vi));
    ASSERT(!sf_machine_expn(2.5, 1.0, &v));        /* non-integer order */
    ASSERT(!sf_machine_expn(2.0, 0.0, &v));
    /* ... and the interpreter answers them as before */
    assert_eval_eq("CosIntegral[0.]", "-Infinity", 0);
    assert_eval_eq("ExpIntegralEi[0.]", "-Infinity", 0);
    assert_eval_eq("LogIntegral[0.]", "0.0", 0);           /* li(0) = 0 */
    assert_true("ExpIntegralEi[800.] > 10^344");
    assert_true("ExpIntegralEi[-800.] < 0 && ExpIntegralEi[-800.] > -1.*^-300");
    /* below DBL_MIN: kept at machine precision in MPFR (numeric-tower review
     * S5), E_2(800) ~ e^-800/800, where it used to underflow to 0. */
    assert_true("ExpIntegralE[2, 800.] == N[ExpIntegralE[2, 800], 20] && !MachineNumberQ[ExpIntegralE[2, 800.]]");
    /* an exact order with a machine argument is a machine computation (exact
     * and inexact arguments numericalize together): E_{5/2}(1) =
     * Int_1^Infinity E^-t t^(-5/2) dt = 0.126487819593254 (NIntegrate agrees) */
    assert_true("MachineNumberQ[ExpIntegralE[5/2, 1.]] && "
                "Abs[ExpIntegralE[5/2, 1.] - 0.126487819593254] < 1.*^-14");
    /* exact input stays exact / symbolic */
    assert_eval_eq("SinIntegral[0]", "0", 0);
    assert_eval_eq("LogIntegral[2]", "LogIntegral[2]", 0);
    assert_eval_eq("ExpIntegralEi[x]", "ExpIntegralEi[x]", 0);
    assert_eval_eq("ExpIntegralE[2, 0]", "1", 0);
    /* machine in => machine out */
    assert_eval_eq("Precision[SinIntegral[Complex[1.5, 0.5]]]", "MachinePrecision", 0);
    assert_eval_eq("Precision[LogIntegral[7.5]]", "MachinePrecision", 0);
    /* arbitrary precision is untouched by the machine path */
    assert_true("Precision[SinIntegral[1.5`30]] > 29");
    assert_true("Precision[ExpIntegralEi[Complex[1.5`30, 2]]] > 29");
    assert_true("Precision[LogIntegral[7.5`25]] > 24");
}

/* ---- five surfaces: scalar == List == packed == NDArray == Compile -------- */

void test_expint_surfaces() {
    const char* defs =
        "xs$ = {0.5, 1.9, 2.1, 3.38, 3.9, 4.1, 7.9, 8.1, 15.9, 16.1, 33.3, 70.5, 0.37}; "
        "zs$ = {Complex[0.3, 0.2], Complex[-0.5, 0.7], Complex[2.5, 1.5], "
        "Complex[-7., -3.], Complex[5.5, -2.9], Complex[0.9, 4.], Complex[60., 2.]};";
    Expr* r = ev(defs); expr_free(r);
    const char* heads[6] = { "SinIntegral", "CosIntegral", "SinhIntegral",
                             "CoshIntegral", "ExpIntegralEi", "LogIntegral" };
    char buf[1024];
    for (int h = 0; h < 6; h++) {
        const char* f = heads[h];
        snprintf(buf, sizeof buf, "s$ = Map[%s, xs$]; "
                 "cr$ = Compile[{{x, _Real}}, %s[x]]; cv$ = Compile[{{v, _Real, 1}}, %s[v]];"
                 "cz$ = Compile[{{z, _Complex}}, %s[z]]; czv$ = Compile[{{v, _Complex, 1}}, %s[v]];",
                 f, f, f, f, f);
        r = ev(buf); expr_free(r);
        snprintf(buf, sizeof buf, "s$ === %s[FromPackedArray[xs$]]", f);         assert_true(buf);
        snprintf(buf, sizeof buf, "s$ === %s[ToPackedArray[xs$]]", f);           assert_true(buf);
        snprintf(buf, sizeof buf, "PackedArrayQ[%s[ToPackedArray[xs$]]]", f);   assert_true(buf);
        snprintf(buf, sizeof buf, "s$ === Normal[%s[NDArray[xs$]]]", f);        assert_true(buf);
        snprintf(buf, sizeof buf, "Head[%s[NDArray[xs$]]] === NDArray", f);      assert_true(buf);
        assert_true("s$ === Map[cr$, xs$]");
        assert_true("s$ === Normal[cv$[xs$]]");
        snprintf(buf, sizeof buf, "c$ = Map[%s, zs$];", f); r = ev(buf); expr_free(r);
        snprintf(buf, sizeof buf, "c$ === %s[zs$]", f);                          assert_true(buf);
        snprintf(buf, sizeof buf, "c$ === Normal[%s[NDArray[zs$]]]", f);         assert_true(buf);
        snprintf(buf, sizeof buf, "Head[%s[NDArray[zs$]]] === NDArray", f);      assert_true(buf);
        assert_true("c$ === Map[cz$, zs$]");
        assert_true("c$ === Normal[czv$[zs$]]");
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{x, _Real}}, %s[x]][[1, 2]]", f);
        assert_true(buf);
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{v, _Real, 1}}, %s[v]][[1, 2]]", f);
        assert_true(buf);
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{z, _Complex}}, %s[z]][[1, 2]]", f);
        assert_true(buf);
    }
    /* ExpIntegralE[2, .] -- the binary kernel, order fixed */
    r = ev("s$ = Map[ExpIntegralE[2, #] &, xs$]; "
           "ce$ = Compile[{{x, _Real}}, ExpIntegralE[2, x]]; "
           "cev$ = Compile[{{v, _Real, 1}}, ExpIntegralE[2, v]]; "
           "c$ = Map[ExpIntegralE[2, #] &, zs$]; "
           "cez$ = Compile[{{z, _Complex}}, ExpIntegralE[2, z]];");
    expr_free(r);
    assert_true("s$ === ExpIntegralE[2, FromPackedArray[xs$]]");
    assert_true("s$ === ExpIntegralE[2, ToPackedArray[xs$]]");
    assert_true("PackedArrayQ[ExpIntegralE[2, ToPackedArray[xs$]]]");
    assert_true("s$ === Normal[ExpIntegralE[2, NDArray[xs$]]]");
    assert_true("Head[ExpIntegralE[2, NDArray[xs$]]] === NDArray");
    assert_true("s$ === Map[ce$, xs$]");
    assert_true("s$ === Normal[cev$[xs$]]");
    assert_true("c$ === ExpIntegralE[2, zs$]");
    assert_true("c$ === Normal[ExpIntegralE[2, NDArray[zs$]]]");
    assert_true("c$ === Map[cez$, zs$]");
    assert_true("CompileDiagnostics[{{x, _Real}}, ExpIntegralE[2, x]][[1, 2]]");
    assert_true("CompileDiagnostics[{{v, _Real, 1}}, ExpIntegralE[2, v]][[1, 2]]");
    assert_true("CompileDiagnostics[{{z, _Complex}}, ExpIntegralE[2, z]][[1, 2]]");
    /* a real buffer with an element whose value is complex (Ci of a negative
     * argument) falls back, and still agrees with the scalar answers */
    assert_true("Map[CosIntegral, {0.5, -2.}] === Normal[CosIntegral[NDArray[{0.5, -2.}]]]");
    assert_true("Map[LogIntegral, {0.5, -2.}] === LogIntegral[{0.5, -2.}]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_expint_real_values);
    TEST(test_expint_complex_values);
    TEST(test_expint_special);
    TEST(test_expint_surfaces);

    printf("All exponential-integral machine-kernel tests passed.\n");
    return 0;
}
