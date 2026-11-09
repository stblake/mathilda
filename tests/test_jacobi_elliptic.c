/* Tests for the Jacobi elliptic family: JacobiSN, JacobiCN, JacobiDN, the nine
 * ratios JacobiCD ... JacobiSD, JacobiAmplitude, InverseJacobiSN / CN / DN.
 *
 * Reference values: mpmath (mp.ellipfun, mp.ellipf) at 50 digits; the amplitude
 * winding and the inverse branches are the reference system's (checked against
 * its outputs: JacobiAmplitude[10., 1.7] = 11.7593, InverseJacobiDN[0.4, 0.3] =
 * 1.71389 + 1.56202 I).
 *
 * Covers: exact rules (u = 0, m = 0, m = 1, parity, the inverses at v = 0, 1 and
 * m = 0, 1), 72 real and 48 complex mpmath points over all twelve ratios
 * (0 <= m <= 1, m > 1, m < 0, m near 0 and 1, |u| up to 40), the amplitude,
 * the inverses (real and the complex branches), arbitrary precision real and
 * complex; identities (sn^2 + cn^2 = 1, dn^2 + m sn^2 = 1, quarter-period
 * shifts, addition); D; Listable; packed / NDArray buffers and Compile[]; the
 * pole; messages; attributes.
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

/* {code, u, m, value} for the twelve ratios, real u and m */
static const double JAC_REAL[][4] = {
    { 0, 0.7, 0.3, 6.3230477631086451725e-1 },
    { 0, -2.9, 0.99, -9.9611939632856453354e-1 },
    { 0, 37.5, 0.5, 4.0102052106056030566e-1 },
    { 0, 1.3, 3.7, (0.3960007593647335057 + 6.5548385003428898406e-55j) },
    { 0, -5.1, -2.4, (-0.60800790675833782007 + 0.0j) },
    { 0, 0.2, 0.001, 1.9866803445626249854e-1 },
    { 1, 0.7, 0.3, 7.747197363269297698e-1 },
    { 1, -2.9, 0.99, 8.8012205165057502427e-2 },
    { 1, 37.5, 0.5, 9.1606907036986938562e-1 },
    { 1, 1.3, 3.7, (0.91825018300164494369 - 1.5673404481596962996e-55j) },
    { 1, -5.1, -2.4, (0.79393097012230501814 + 0.0j) },
    { 1, 0.2, 0.001, 9.8006684062123298461e-1 },
    { 2, 0.7, 0.3, 9.3811363968143021572e-1 },
    { 2, -2.9, 0.99, 1.3292361255787483891e-1 },
    { 2, 37.5, 0.5, 9.5895321619157123205e-1 },
    { 2, 1.3, 3.7, (-0.64790321403389515698 + 8.8262538870419963033e-55j) },
    { 2, -5.1, -2.4, 1.3737600500937467016 },
    { 2, 0.2, 0.001, 9.9998026531131367387e-1 },
    { 3, 0.7, 0.3, 8.2582717440289361384e-1 },
    { 3, -2.9, 0.99, 6.6212619015855483036e-1 },
    { 3, 37.5, 0.5, 9.5528025236516352714e-1 },
    { 3, 1.3, 3.7, (-1.4172644356624638432 - 1.3995255208532746553e-54j) },
    { 3, -5.1, -2.4, (0.57792550458002211397 + 0.0j) },
    { 3, 0.2, 0.001, 9.8008618231692677903e-1 },
    { 4, 0.7, 0.3, 1.2252315107391324966 },
    { 4, -2.9, 0.99, -8.8355076198141972561e-2 },
    { 4, 37.5, 0.5, 2.2843446214352925274 },
    { 4, 1.3, 3.7, (2.3188091469185733416 - 3.0004796731339512023e-54j) },
    { 4, -5.1, -2.4, (-1.3057905354475352332 + 0.0j) },
    { 4, 0.2, 0.001, 4.9331883878732304011 },
    { 5, 0.7, 0.3, 1.2109071134926510144 },
    { 5, -2.9, 0.99, 1.510286128027252387 },
    { 5, 37.5, 0.5, 1.0468132231605495297 },
    { 5, 1.3, 3.7, (-0.70558462827198207049 + 8.9716284607270185441e-55j) },
    { 5, -5.1, -2.4, (1.7303268190710825258 + 0.0j) },
    { 5, 0.2, 0.001, 1.0203184352992273538 },
    { 6, 0.7, 0.3, 1.4836415520293629743 },
    { 6, -2.9, 0.99, -1.3344144592284468726e-1 },
    { 6, 37.5, 0.5, 2.3912821559741439111 },
    { 6, 1.3, 3.7, (-1.6361160899622134304 + 3.9477389078815407678e-54j) },
    { 6, -5.1, -2.4, (-2.2594443835740592709 + 0.0j) },
    { 6, 0.2, 0.001, 5.0334230569511323282 },
    { 7, 0.7, 0.3, 1.2907893695094951454 },
    { 7, -2.9, 0.99, 1.1362060502002042457e1 },
    { 7, 37.5, 0.5, 1.0916207438334785238 },
    { 7, 1.3, 3.7, (1.0890278254354658926 + 1.7627047498206131673e-55j) },
    { 7, -5.1, -2.4, (1.2595553488056399415 + 0.0j) },
    { 7, 0.2, 0.001, 1.020338571363287859 },
    { 8, 0.7, 0.3, 1.0659689377713189359 },
    { 8, -2.9, 0.99, 7.5231178325416093195 },
    { 8, 37.5, 0.5, 1.0428037396562928914 },
    { 8, 1.3, 3.7, (-1.5434404064365157558 - 1.9568510787250834334e-54j) },
    { 8, -5.1, -2.4, 7.2792916050496521731e-1 },
    { 8, 0.2, 0.001, 1.0000197350781519497 },
    { 9, 0.7, 0.3, 1.5815158092501310657 },
    { 9, -2.9, 0.99, -1.0038957214222896936 },
    { 9, 37.5, 0.5, 2.4936379748231999381 },
    { 9, 1.3, 3.7, (2.5252476828686016732 - 3.1207442018538809984e-54j) },
    { 9, -5.1, -2.4, (-1.644715453342723587 + 0.0j) },
    { 9, 0.2, 0.001, 5.0335223919485330841 },
    { 10, 0.7, 0.3, 8.1617228355214317202e-1 },
    { 10, -2.9, 0.99, -1.1317968848302901432e1 },
    { 10, 37.5, 0.5, 4.3776231949261798076e-1 },
    { 10, 1.3, 3.7, (0.43125584584176893555 + 7.1062126239428610043e-55j) },
    { 10, -5.1, -2.4, (-0.76581961107358519929 + 0.0j) },
    { 10, 0.2, 0.001, 2.0270865845265532465e-1 },
    { 11, 0.7, 0.3, 6.7401725075182367823e-1 },
    { 11, -2.9, 0.99, -7.4939235938600067216 },
    { 11, 37.5, 0.5, 4.1818569904086744944e-1 },
    { 11, 1.3, 3.7, (-0.61120357298307315488 - 1.0918994786371609154e-54j) },
    { 11, -5.1, -2.4, (-0.44258668514697801712 + 0.0j) },
    { 11, 0.2, 0.001, 1.986719551854487871e-1 },
};
/* {code, ur, ui, mr, mi, Re, Im} complex u and / or m */
static const double JAC_CPLX[][7] = {
    { 0, 1.1, 0.4, 0.3, 0.0, 9.2588329999237267532e-1, 1.753361665838469129e-1 },
    { 0, -0.6, 1.7, 2.2, 0.0, -9.6582717509243375888e-1, -3.8229216611833445982e-1 },
    { 0, 0.7, 0.0, 0.3, 0.2, 6.3231365105086478467e-1, -7.9566853313217092273e-3 },
    { 0, 2.5, -0.3, -1.5, 0.0, -3.301431704526265343e-2, 2.9752625174849620506e-1 },
    { 1, 1.1, 0.4, 0.3, 0.0, 5.2038726647201111286e-1, -3.1196156974643425044e-1 },
    { 1, -0.6, 1.7, 2.2, 0.0, -7.0070565389013937187e-1, 5.2693761041055919009e-1 },
    { 1, 0.7, 0.0, 0.3, 0.2, 7.7478056407758312336e-1, 6.493607332678385779e-3 },
    { 1, 2.5, -0.3, -1.5, 0.0, -1.0428425789964588129, -9.4190879835055530074e-3 },
    { 2, 1.1, 0.4, 0.3, 0.0, 8.6901420855995002703e-1, -5.6043098119309497013e-2 },
    { 2, -0.6, 1.7, 2.2, 0.0, -7.2480128872169516319e-1, 1.1207236672124452286 },
    { 2, 0.7, 0.0, 0.3, 0.2, 9.3794579061834001708e-1, -4.1011312813066853399e-2 },
    { 2, 2.5, -0.3, -1.5, 0.0, 9.3225634745401267669e-1, -1.5804600362345121685e-2 },
    { 3, 1.1, 0.4, 0.3, 0.0, 6.1939970027247696664e-1, -3.1903792692690662499e-1 },
    { 3, -0.6, 1.7, 2.2, 0.0, 6.1662144096608671255e-1, 2.264408669584959921e-1 },
    { 3, 0.7, 0.0, 0.3, 0.2, 8.2416143034707273404e-1, 4.2959358594213962301e-2 },
    { 3, 2.5, -0.3, -1.5, 0.0, -1.1181294843260084866, -2.9059257907567648158e-2 },
    { 4, 1.1, 0.4, 0.3, 0.0, 4.8098925229519199659e-1, -4.2801979624772220366e-1 },
    { 4, -0.6, 1.7, 2.2, 0.0, 4.4052818132925971001e-1, -7.199508369813103945e-1 },
    { 4, 0.7, 0.0, 0.3, 0.2, 1.224987317605688931, 2.5684161533973080894e-2 },
    { 4, 2.5, -0.3, -1.5, 0.0, 3.5292566469686499232e-1, 3.4658823318905956923 },
    { 5, 1.1, 0.4, 0.3, 0.0, 1.2759520951294646765, 6.5721231558406041176e-1 },
    { 5, -0.6, 1.7, 2.2, 0.0, 1.429026789589252019, -5.2477913293855678218e-1 },
    { 5, 0.7, 0.0, 0.3, 0.2, 1.2100667622437920739, -6.3074647815390999253e-2 },
    { 5, 2.5, -0.3, -1.5, 0.0, -8.937471311935900519e-1, 2.3227746655083097193e-2 },
    { 6, 1.1, 0.4, 0.3, 0.0, 8.9501912560858814579e-1, -2.3002069550719875299e-1 },
    { 6, -0.6, 1.7, 2.2, 0.0, 2.5171139669910356549e-1, -1.2600090302664651595 },
    { 6, 0.7, 0.0, 0.3, 0.2, 1.4839364566480127608, -4.6186093447937993491e-2 },
    { 6, 2.5, -0.3, -1.5, 0.0, -3.9593093708889735833e-1, -3.0894247232541144618 },
    { 7, 1.1, 0.4, 0.3, 0.0, 1.4136236206503350934, 8.4743857534883645436e-1 },
    { 7, -0.6, 1.7, 2.2, 0.0, -9.1160364515804722475e-1, -6.8553499426515264881e-1 },
    { 7, 0.7, 0.0, 0.3, 0.2, 1.2905973722439920935, -1.0816782129681615907e-2 },
    { 7, 2.5, -0.3, -1.5, 0.0, -9.5883927972621014851e-1, 8.6603594058014901733e-3 },
    { 8, 1.1, 0.4, 0.3, 0.0, 1.1459630932060951765, 7.3903650183213408986e-2 },
    { 8, -0.6, 1.7, 2.2, 0.0, -4.0688121483550353215e-1, -6.2913989572856977779e-1 },
    { 8, 0.7, 0.0, 0.3, 0.2, 1.0641252583331264691, 4.6528460683228841832e-2 },
    { 8, 2.5, -0.3, -1.5, 0.0, 1.072358133009334172, 1.8179754724983622927e-2 },
    { 9, 1.1, 0.4, 0.3, 0.0, 1.0426582546766957685, -1.9745004735858420216e-1 },
    { 9, -0.6, 1.7, 2.2, 0.0, -8.9513858879577282412e-1, 3.5431232306555740276e-1 },
    { 9, 0.7, 0.0, 0.3, 0.2, 1.5812432331136152227, 1.9897490457872045324e-2 },
    { 9, 2.5, -0.3, -1.5, 0.0, -3.6841477678722590463e-1, -3.3201676556259701658 },
    { 10, 1.1, 0.4, 0.3, 0.0, 1.1602638716179567435, 1.0324885713220236272 },
    { 10, -0.6, 1.7, 2.2, 0.0, 6.1837691549941693346e-1, 1.0106072590869595709 },
    { 10, 0.7, 0.0, 0.3, 0.2, 8.1597627074854732037e-1, -1.7108476181416846842e-2 },
    { 10, 2.5, -0.3, -1.5, 0.0, 2.9078739703529433893e-2, -2.8556577277731622695e-1 },
    { 11, 1.1, 0.4, 0.3, 0.0, 1.0480701076874480949, 2.6935493132243993894e-1 },
    { 11, -0.6, 1.7, 2.2, 0.0, 1.5246168079321400811e-1, 7.6318790920179692847e-1 },
    { 11, 0.7, 0.0, 0.3, 0.2, 6.7323113959267109987e-1, 2.0953671018720931778e-2 },
    { 11, 2.5, -0.3, -1.5, 0.0, -4.0812125670267487625e-2, 3.1845450365998681613e-1 },
};
/* JacobiAmplitude {u, m, am}: sin am = sn, cos am = cn, the reference winding */
static const double JAC_AMP[][3] = {
    { 10.0, 0.3, 9.1425592025177380561 },
    { 10.0, 1.7, 1.1759291218739125424e1 },
    { 10.0, -2.3, 1.3676418360342770553e1 },
    { -25.0, 0.9, -1.4989314803345621939e1 },
    { 0.4, 0.5, 3.9485266598872287541e-1 },
};
/* {code 0..2 (SN, CN, DN), v, m, value}, real */
static const double JAC_INV[][4] = {
    { 0, 0.4, 0.3, 4.1496082248971368823e-1 },
    { 0, -0.95, 0.9, -1.6832446717979713534 },
    { 1, 0.4, 0.3, 1.2276118096732333249 },
    { 1, -0.7, 0.6, 3.0515454975849816346 },
    { 2, 0.9, 0.3, 9.5662560068323307025e-1 },
    { 2, 0.99, 0.05, 6.8511453134551220248e-1 },
};
static const char* const JAC_MP_SN = "0.326054647440372270290495028548404963658354776`45";
static const char* const JAC_MP_CD[2] = { "1.50064801209331926327799501991472572222914596`45", "-0.490237683752492257533212380660580303307365905`45" };

static const char* const NAMES[12] = {
    "JacobiSN", "JacobiCN", "JacobiDN", "JacobiCD", "JacobiCS", "JacobiDC",
    "JacobiDS", "JacobiNC", "JacobiND", "JacobiNS", "JacobiSC", "JacobiSD" };

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

/* 64 ulp */
#define JAC_TOL 1.43e-14

void test_jac_exact() {
    assert_true("JacobiSN[-u, m] === -JacobiSN[u, m]");
    assert_true("JacobiCN[-u, m] === JacobiCN[u, m]");
    assert_true("JacobiSC[-u, m] === -JacobiSC[u, m]");
    assert_true("JacobiDN[-u, m] === JacobiDN[u, m]");
    assert_true("JacobiSN[u, 0] === Sin[u] && JacobiCN[u, 0] === Cos[u] && JacobiDN[u, 0] === 1");
    assert_true("JacobiCS[u, 0] === Cot[u] && JacobiSC[u, 0] === Tan[u] && JacobiNS[u, 0] === Csc[u]");
    assert_true("JacobiSN[u, 1] === Tanh[u] && JacobiCN[u, 1] === Sech[u] && JacobiDN[u, 1] === Sech[u]");
    assert_true("JacobiCS[u, 1] === Csch[u] && JacobiNC[u, 1] === Cosh[u] && JacobiSD[u, 1] === Sinh[u]");
    assert_true("JacobiSN[0, m] === 0 && JacobiCN[0, m] === 1 && JacobiDN[0, m] === 1");
    assert_eval_eq("JacobiNS[0, m]", "ComplexInfinity", 0);
    assert_eval_eq("JacobiNS[0., 0.5]", "ComplexInfinity", 0);
    assert_true("JacobiAmplitude[u, 0] === u && JacobiAmplitude[0, m] === 0");
    assert_true("JacobiAmplitude[u, 1] === 2 ArcTan[Tanh[u/2]]");
    assert_true("InverseJacobiSN[v, 0] === ArcSin[v] && InverseJacobiSN[v, 1] === ArcTanh[v]");
    assert_true("InverseJacobiSN[1, m] === EllipticK[m] && InverseJacobiCN[0, m] === EllipticK[m]");
    assert_true("InverseJacobiCN[v, 0] === ArcCos[v] && InverseJacobiDN[v, 1] === ArcSech[v]");
    assert_true("InverseJacobiCN[1, m] === 0 && InverseJacobiDN[1, m] === 0 && InverseJacobiSN[0, m] === 0");
    assert_eval_eq("JacobiSN[u, 1/2]", "JacobiSN[u, 1/2]", 0);
    assert_true("NumericQ[JacobiSN[1, 1/2]]");
}

void test_jac_machine() {
    char buf[512], b0[64], b1[64], b2[64], b3[64];
    for (size_t i = 0; i < sizeof(JAC_REAL) / sizeof(JAC_REAL[0]); i++) {
        const double* c = JAC_REAL[i];
        snprintf(buf, sizeof buf, "%s[%s, %s]", NAMES[(int)c[0]], R(c[1], b0), R(c[2], b1));
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(im == 0.0 && fabs(re - c[3]) <= JAC_TOL * fabs(c[3]),
                   "%s: expected %.17g, got %.17g", buf, c[3], re);
    }
    for (size_t i = 0; i < sizeof(JAC_CPLX) / sizeof(JAC_CPLX[0]); i++) {
        const double* c = JAC_CPLX[i];
        char mstr[160];
        if (c[4] == 0.0) snprintf(mstr, sizeof mstr, "%s", R(c[3], b2));
        else snprintf(mstr, sizeof mstr, "Complex[%s, %s]", R(c[3], b2), R(c[4], b3));
        snprintf(buf, sizeof buf, "%s[Complex[%s, %s], %s]", NAMES[(int)c[0]],
                 R(c[1], b0), R(c[2], b1), mstr);
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(hypot(re - c[5], im - c[6]) <= JAC_TOL * hypot(c[5], c[6]),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", buf, c[5], c[6], re, im);
    }
    for (size_t i = 0; i < sizeof(JAC_AMP) / sizeof(JAC_AMP[0]); i++) {
        snprintf(buf, sizeof buf, "JacobiAmplitude[%s, %s]", R(JAC_AMP[i][0], b0), R(JAC_AMP[i][1], b1));
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(fabs(re - JAC_AMP[i][2]) <= JAC_TOL * fabs(JAC_AMP[i][2]),
                   "%s: expected %.17g, got %.17g", buf, JAC_AMP[i][2], re);
    }
    static const char* const INV[3] = { "InverseJacobiSN", "InverseJacobiCN", "InverseJacobiDN" };
    for (size_t i = 0; i < sizeof(JAC_INV) / sizeof(JAC_INV[0]); i++) {
        const double* c = JAC_INV[i];
        snprintf(buf, sizeof buf, "%s[%s, %s]", INV[(int)c[0]], R(c[1], b0), R(c[2], b1));
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(fabs(re - c[3]) <= JAC_TOL * fabs(c[3]), "%s: expected %.17g, got %.17g", buf, c[3], re);
    }
    /* the complex branches of the inverses, and the complex amplitude */
    assert_true("Abs[InverseJacobiSN[1.4, 0.3] - (1.713889448178791 - 1.1163042482101442 I)] < 10^-14");
    assert_true("Abs[InverseJacobiCN[-1.4, 0.3] - (3.427778896357582 - 0.8332190166425625 I)] < 10^-14");
    assert_true("Abs[InverseJacobiDN[0.4, 0.3] - (1.7138894481787912 + 1.5620200866917282 I)] < 10^-14");
    assert_true("Abs[JacobiAmplitude[1. + I, 0.5] - (1.0777172498501961 + 0.794587121938991 I)] < 10^-14");
    assert_true("MachineNumberQ[InverseJacobiSN[1.4, 0.3]] || MachineNumberQ[Re[InverseJacobiSN[1.4, 0.3]]]");
}

void test_jac_arbitrary() {
    char buf[512];
    snprintf(buf, sizeof buf, "With[{v = N[JacobiSN[1/3, 1/5], 40], r = %s}, Abs[v - r] < 10^-38 Abs[r] && Precision[v] >= 39]", JAC_MP_SN);
    assert_true(buf);
    snprintf(buf, sizeof buf, "With[{v = N[JacobiCD[1/2 + 2 I, 1/3 + I/4], 40], r = %s + %s I}, Abs[v - r] < 10^-37 Abs[r]]", JAC_MP_CD[0], JAC_MP_CD[1]);
    assert_true(buf);
    assert_true("Abs[N[JacobiAmplitude[10, 17/10], 30] - 11.759291218739126`16] < 10^-14");
    assert_true("Abs[N[InverseJacobiDN[2/5, 3/10], 30] - N[EllipticK[3/10] + I EllipticF[ArcSin[Sqrt[(1 - 3/10 - 4/25)/((7/10) (21/25))]], 7/10], 35]] < 10^-28");
    assert_true("Abs[JacobiSN[N[EllipticK[1/3], 40], 1/3`40] - 1] < 10^-38");
}

void test_jac_identities() {
    assert_true("With[{u = 2.37, m = 0.61}, Abs[JacobiSN[u, m]^2 + JacobiCN[u, m]^2 - 1] < 10^-15 && Abs[JacobiDN[u, m]^2 + m JacobiSN[u, m]^2 - 1] < 10^-15]");
    assert_true("With[{u = 1.1 + 0.7 I, m = 0.4}, Abs[JacobiSN[u, m]^2 + JacobiCN[u, m]^2 - 1] < 10^-14]");
    assert_true("With[{u = 7/10`40, m = 1/3`40}, Abs[JacobiSN[u + EllipticK[m], m] - JacobiCD[u, m]] < 10^-37]");
    /* the addition formula for sn */
    assert_true("With[{a = 0.37, b = 1.21, m = 0.83}, Abs[JacobiSN[a + b, m] - (JacobiSN[a, m] JacobiCN[b, m] JacobiDN[b, m] + JacobiSN[b, m] JacobiCN[a, m] JacobiDN[a, m])/(1 - m JacobiSN[a, m]^2 JacobiSN[b, m]^2)] < 10^-15]");
    /* inverses round trip */
    assert_true("Abs[JacobiSN[InverseJacobiSN[0.63, 0.4], 0.4] - 0.63] < 10^-15");
    assert_true("Abs[JacobiDN[InverseJacobiDN[0.4, 0.3], 0.3] - 0.4] < 10^-14");
    /* the amplitude inverts EllipticF */
    assert_true("Abs[EllipticF[JacobiAmplitude[5.5, 0.7], 0.7] - 5.5] < 10^-14");
}

void test_jac_derivative() {
    assert_true("D[JacobiSN[u, m], u] === JacobiCN[u, m] JacobiDN[u, m]");
    assert_true("D[JacobiDN[u, m], u] === -m JacobiCN[u, m] JacobiSN[u, m]");
    assert_true("D[JacobiCD[u, m], u] === (-1 + m) JacobiND[u, m] JacobiSD[u, m]");
    assert_true("D[JacobiAmplitude[u, m], u] === JacobiDN[u, m]");
    assert_true("D[InverseJacobiSN[v, m], v] === JacobiCD[InverseJacobiSN[v, m], m]/(1 - v^2)");
    assert_true("D[JacobiSN[u, m], m] === Derivative[0, 1][JacobiSN][u, m]");
    assert_true("D[JacobiSN[x^2, m], x] === 2 x JacobiCN[x^2, m] JacobiDN[x^2, m]");
    /* every ratio's derivative agrees with a difference quotient */
    for (int k = 0; k < 12; k++) {
        char buf[400];
        snprintf(buf, sizeof buf,
                 "Abs[(D[%s[u, 0.6], u] /. u -> 0.77) - (%s[0.77 + 10.^-6, 0.6] - %s[0.77 - 10.^-6, 0.6])/(2 10.^-6)] < 10^-8 Max[1, Abs[%s[0.77, 0.6]]]",
                 NAMES[k], NAMES[k], NAMES[k], NAMES[k]);
        assert_true(buf);
    }
}

void test_jac_buffers() {
    for (int k = 0; k < 12; k++) {
        char buf[400];
        snprintf(buf, sizeof buf,
                 "With[{v = Range[-30., 30., 0.1]}, PackedArrayQ[%s[v, 0.37]] && %s[v, 0.37] === Table[%s[x, 0.37], {x, FromPackedArray[v]}]]",
                 NAMES[k], NAMES[k], NAMES[k]);
        assert_true(buf);
    }
    assert_true("With[{v = Range[-5., 5., 0.1]}, JacobiSN[2.2, v] === Table[JacobiSN[2.2, m], {m, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[-3., 3., 0.25] + 0.5 I}, Normal[JacobiCN[NDArray[v], 0.3]] === JacobiCN[FromPackedArray[v], 0.3]]");
    assert_true("With[{v = Range[-20., 20., 0.1]}, JacobiAmplitude[v, 0.8] === Table[JacobiAmplitude[x, 0.8], {x, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[-1., 1., 0.01]}, InverseJacobiSN[v, 0.6] === Table[InverseJacobiSN[x, 0.6], {x, FromPackedArray[v]}]]");
    assert_true("With[{v = Range[0., 1., 0.01]}, InverseJacobiDN[v, 0.6] === Table[InverseJacobiDN[x, 0.6], {x, FromPackedArray[v]}]]");
}

void test_jac_compile() {
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, JacobiSN[x, 0.3]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, JacobiCD[v, 0.3]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, JacobiAmplitude[v, 0.3]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, InverseJacobiCN[v, 0.3]][[1]]", "\"Compiled\" -> True", 0);
    assert_true("With[{f = Compile[{{v, _Real, 1}}, JacobiDN[v, 1.7]], v = Range[-10., 10., 0.01]}, f[v] === JacobiDN[v, 1.7]]");
    assert_true("With[{f = Compile[{{x, _Real}, {m, _Real}}, JacobiSN[x, m] + JacobiCN[x, m]]}, f[0.7, -2.3] === JacobiSN[0.7, -2.3] + JacobiCN[0.7, -2.3]]");
}

void test_jac_messages_attributes() {
    test_mute_messages();
    assert_eval_eq("Check[JacobiSN[1.], \"failed\"]", "\"failed\"", 0);
    assert_eval_eq("Check[InverseJacobiCN[1., 2., 3.], \"failed\"]", "\"failed\"", 0);
    for (int k = 0; k < 12; k++) {
        char buf[200];
        snprintf(buf, sizeof buf, "Attributes[%s]", NAMES[k]);
        assert_eval_eq(buf, "{Listable, NumericFunction, Protected}", 0);
        SymbolDef* d = symtab_get_def(NAMES[k]);
        ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, NAMES[k]));
    }
    assert_eval_eq("Attributes[JacobiAmplitude]", "{Listable, NumericFunction, Protected}", 0);
    assert_eval_eq("Attributes[InverseJacobiDN]", "{Listable, NumericFunction, Protected}", 0);
}

int main() {
    symtab_init();
    core_init();
    TEST(test_jac_exact);
    TEST(test_jac_machine);
    TEST(test_jac_arbitrary);
    TEST(test_jac_identities);
    TEST(test_jac_derivative);
    TEST(test_jac_buffers);
    TEST(test_jac_compile);
    TEST(test_jac_messages_attributes);
    printf("All Jacobi elliptic tests passed.\n");
    return 0;
}
