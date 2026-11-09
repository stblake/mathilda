/* Tests for the Weierstrass elliptic functions: WeierstrassP, WeierstrassPPrime,
 * WeierstrassZeta, WeierstrassSigma [u, {g2, g3}], WeierstrassHalfPeriods,
 * WeierstrassInvariants.
 *
 * Reference values: mpmath at 50 digits through the theta-function formulas on
 * the lattice of the normalised half-periods (checked against the defining
 * identities wp'^2 = 4 wp^3 - g2 wp - g3, zeta' = -wp, sigma'/sigma = zeta to
 * 1e-60). The half-period normalisation was recovered from the reference
 * system's outputs on 380 sampled invariants; the values here are mpmath's.
 *
 * Covers: exact rules ({0, 0}, a double root, the pole at 0, parity, the exact
 * square / hexagonal invariants, NumericQ); 72 machine points (real and complex
 * u and invariants, |u| up to 45 periods, invariants 1e-5 .. 1e4); the
 * half-periods (the boundary and rotation rules); 40 digits; identities; D
 * (in u and in g2, g3, checked against central differences); buffers, NDArray,
 * Compile[]; attributes.
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

enum { WP_P, WP_PPRIME, WP_ZETA, WP_SIGMA };
static const char* const NAMES[4] = { "WeierstrassP", "WeierstrassPPrime",
                                      "WeierstrassZeta", "WeierstrassSigma" };

/* {kind, Re u, Im u, Re g2, Im g2, Re g3, Im g3, Re value, Im value} and
 * {Re g2, Im g2, Re g3, Im g3, Re w1, Im w1, Re w3, Im w3} */
static const double WP_REF[][10] = {
    { 0, 0.3, 0, 1.0, 0.0, 2.0, 0.0, 11.116190356306905838, 0.0 },
    { 1, 0.3, 0, 1.0, 0.0, 2.0, 0.0, -74.036345855448816823, 0.0 },
    { 2, 0.3, 0, 1.0, 0.0, 2.0, 0.0, 3.3328485908173352575, 0.0 },
    { 3, 0.3, 0, 1.0, 0.0, 2.0, 0.0, 0.29998935416205232314, 0.0 },
    { 0, 0.5, 0.3, 1.0, 1.0, 2.0, 0.0, 1.3724978877263547162, -2.5653819747410041358 },
    { 1, 0.5, 0.3, 1.0, 1.0, 2.0, 0.0, 0.52567954015986757912, 10.210993564823919092 },
    { 2, 0.5, 0.3, 1.0, 1.0, 2.0, 0.0, 1.4749220190322919658, -0.88589088367227488936 },
    { 3, 0.5, 0.3, 1.0, 1.0, 2.0, 0.0, 0.50041729180978807675, 0.30016729626401864367 },
    { 0, 1.7, -0.4, -3.0, 0.0, 0.5, 0.0, 0.046231798581325253029, 0.08948877669477910831 },
    { 1, 1.7, -0.4, -3.0, 0.0, 0.5, 0.0, -0.2093955737180954398, -0.63968638856437169937 },
    { 2, 1.7, -0.4, -3.0, 0.0, 0.5, 0.0, 0.73314411757768765287, 0.053290782912307596556 },
    { 3, 1.7, -0.4, -3.0, 0.0, 0.5, 0.0, 1.787236280736836457, -0.54915007321555094173 },
    { 0, 0.2, 0.9, 0.0, 0.0, 1.0, 0.0, -1.0491724137104673227, -0.51803434316544150678 },
    { 1, 0.2, 0.9, 0.0, 0.0, 1.0, 0.0, 1.4888361807414934327, -2.1112856759912157332 },
    { 2, 0.2, 0.9, 0.0, 0.0, 1.0, 0.0, 0.2310707121745542554, -1.0610125880499654067 },
    { 3, 0.2, 0.9, 0.0, 0.0, 1.0, 0.0, 0.20067349621042442631, 0.90002704259812801394 },
    { 0, 0.8, 0, 3.0, 0.0, 0.0, 0.0, 1.6604848332914031954, 0.0 },
    { 1, 0.8, 0, 3.0, 0.0, 0.0, 0.0, -3.6512691303787081522, 0.0 },
    { 2, 0.8, 0, 3.0, 0.0, 0.0, 0.0, 1.2241739445661270916, 0.0 },
    { 3, 0.8, 0, 3.0, 0.0, 0.0, 0.0, 0.79589651222748092874, 0.0 },
    { 0, 6.83, 0.776, 0.3279, -0.4894, 2.3239, 0.5105, -0.2188270357811235332, 0.32923858992720308018 },
    { 1, 6.83, 0.776, 0.3279, -0.4894, 2.3239, 0.5105, -0.22777281298910560926, 1.4907802380469006378 },
    { 2, 6.83, 0.776, 0.3279, -0.4894, 2.3239, 0.5105, 3.4490416971153185005, -0.31829850457036280411 },
    { 3, 6.83, 0.776, 0.3279, -0.4894, 2.3239, 0.5105, -179069.36441720808354, 3240.3355212642417706 },
    { 0, -1.748, 2.589, 155.6, 0.0, 296.36, 0.0, 24.32424974281364355, -9.6226368477837354916 },
    { 1, -1.748, 2.589, 155.6, 0.0, 296.36, 0.0, -217.97796043091847942, 145.10569213373912054 },
    { 2, -1.748, 2.589, 155.6, 0.0, 296.36, 0.0, -1.7568668104719594207, -5.0774909917870508764 },
    { 3, -1.748, 2.589, 155.6, 0.0, 296.36, 0.0, 5519.941011416189371, 5761.9121194334781567 },
    { 0, 10.0, 7.0, 1.0, 0.0, 2.0, 0.0, 0.80649749204912725441, 0.79426899395262941578 },
    { 1, 10.0, 7.0, 1.0, 0.0, 2.0, 0.0, -0.63308574653538694645, -2.6859753035251995159 },
    { 2, 10.0, 7.0, 1.0, 0.0, 2.0, 0.0, 5.6775991912927805568, -2.9658309321629684875 },
    { 3, 10.0, 7.0, 1.0, 0.0, 2.0, 0.0, 8011616381011163.9975, -2868391473107916.0893 },
    { 0, 2.5, 0, -2.0, 0.0, -3.0, 0.0, -0.38002152239977027807, 0.0 },
    { 1, 2.5, 0, -2.0, 0.0, -3.0, 0.0, 1.4214188893756277901, 0.0 },
    { 2, 2.5, 0, -2.0, 0.0, -3.0, 0.0, 1.4803346299750317656, 0.0 },
    { 3, 2.5, 0, -2.0, 0.0, -3.0, 0.0, 5.2803704452803052881, 0.0 },
    { 0, 0.01, 0, 1.0, 0.0, 2.0, 0.0, 10000.00000500071387, 0.0 },
    { 1, 0.01, 0, 1.0, 0.0, 2.0, 0.0, -1999999.9989997141603, 0.0 },
    { 2, 0.01, 0, 1.0, 0.0, 2.0, 0.0, 99.999999983331902679, 0.0 },
    { 3, 0.01, 0, 1.0, 0.0, 2.0, 0.0, 0.009999999999583309732, 0.0 },
    { 0, 25.3, -11.2, 4.0, 0.0, -1.0, 0.0, 0.27605074915294329195, -0.037961618009459529839 },
    { 1, 25.3, -11.2, 4.0, 0.0, -1.0, 0.0, 0.21808006865927964819, 0.26905568065987724284 },
    { 2, 25.3, -11.2, 4.0, 0.0, -1.0, 0.0, 8.0537964587541913177, 6.0156531341713676955 },
    { 3, 25.3, -11.2, 4.0, 0.0, -1.0, 0.0, -3.0241561118190924854e+58, -1.6445827611709359778e+58 },
    { 0, 0.7, 0.2, 1.3, 0.0, 2.1, -0.5, 1.644972145726602405, -0.9615908164709637007 },
    { 1, 0.7, 0.2, 1.3, 0.0, 2.1, -0.5, -3.2902490233356438353, 3.9385124739071484801 },
    { 2, 0.7, 0.2, 1.3, 0.0, 2.1, -0.5, 1.3138921595989292296, -0.38647080313906398505 },
    { 3, 0.7, 0.2, 1.3, 0.0, 2.1, -0.5, 0.69984282599305690008, 0.19863401904010131809 },
    { 0, 3.1, 1.2, -1.0, 2.0, 0.5, 1.0, 0.60570201426192062821, 0.06747333326005469602 },
    { 1, 3.1, 1.2, -1.0, 2.0, 0.5, 1.0, 1.2738324206762691771, -0.72541258264145383698 },
    { 2, 3.1, 1.2, -1.0, 2.0, 0.5, 1.0, 1.0072194084345398127, -0.55150804974529466021 },
    { 3, 3.1, 1.2, -1.0, 2.0, 0.5, 1.0, 9.6664122687837713506, 4.6163882240350488173 },
    { 0, 0.4, -0.05, 1000.0, 0.0, 10000.0, 0.0, 26.751507955325249196, -14.103775770098969102 },
    { 1, 0.4, -0.05, 1000.0, 0.0, 10000.0, 0.0, 193.30315454795219153, -247.78069726427407644 },
    { 2, 0.4, -0.05, 1000.0, 0.0, 10000.0, 0.0, 0.63022360192853771661, 1.5655747748131568985 },
    { 3, 0.4, -0.05, 1000.0, 0.0, 10000.0, 0.0, 0.34987460587869651382, -0.0066369970990896631879 },
    { 0, 40.0, 3.0, 0.001, 0.0, 1e-05, 0.0, 0.0070703380636549266881, -0.045627155184231788095 },
    { 1, 40.0, 3.0, 0.001, 0.0, 1e-05, 0.0, 0.011178493607080922269, 0.017811423380851004767 },
    { 2, 40.0, 3.0, 0.001, 0.0, 1e-05, 0.0, 0.53143883084125148151, -0.14304599471855535224 },
    { 3, 40.0, 3.0, 0.001, 0.0, 1e-05, 0.0, -3058.5059102887877687, 12716.538037262926099 },
    { 0, 1.2, 0, 10.0, 0.0, 1.0, 0.0, 1.8117383650661150814, 0.0 },
    { 1, 1.2, 0, 10.0, 0.0, 1.0, 0.0, 2.1610152012009069598, 0.0 },
    { 2, 1.2, 0, 10.0, 0.0, 1.0, 0.0, 0.47663341469101963788, 0.0 },
    { 3, 1.2, 0, 10.0, 0.0, 1.0, 0.0, 1.0888363069903560199, 0.0 },
    { 0, 0.6, 0.6, 2.0, 0.0, 0.0, 0.0, 0.0, -1.3181232013606045993 },
    { 1, 0.6, 0.6, 2.0, 0.0, 0.0, 0.0, 2.4286758288496534507, 2.4286758288496534507 },
    { 2, 0.6, 0.6, 2.0, 0.0, 0.0, 0.0, 0.84762722927391568745, -0.84762722927391568745 },
    { 3, 0.6, 0.6, 2.0, 0.0, 0.0, 0.0, 0.60258799998810291709, 0.60258799998810291709 },
    { 0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.96438364632990627616, 0.0 },
    { 1, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, -2.1418789253489088842, 0.0 },
    { 2, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0071339483266371232, 0.0 },
    { 3, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0011904415019613417, 0.0 },
};
static const double WP_HP[][8] = {
    { 2.0, 0.0, 3.0, 0.0, 1.1972208897783685677, 0.0, 0.59861044488918428387, 1.1751406146397711734 },
    { 4.0, 0.0, -1.0, 0.0, 0.0, -1.2256946909933950304, 1.4967293231159798149, 0.0 },
    { 10.0, 0.0, -1.0, 0.0, 0.0, -1.0217853591644661428, 1.0674251022120500683, 0.0 },
    { -3.0, 0.0, 1.0, 0.0, 0.84723513917405365084, 1.0500353672047947276, -0.84723513917405365084, 1.0500353672047947276 },
    { 1.0, 2.0, 3.0, -1.0, 0.4934326673920575794, 1.0959930980483127444, -1.2415592276212725522, 0.017378106212533579683 },
    { -1.0, 2.0, 0.5, 1.0, 0.78677127280667013854, 1.1064106993734195093, -1.3554174052609612291, 0.45408917727515466145 },
    { 3.0, 0.0, 0.0, 0.0, 1.408792103676543125, 0.0, 0.0, 1.408792103676543125 },
    { 0.0, 0.0, 1.0, 0.0, 1.5299540370571928749, 0.0, 0.76497701852859643746, 1.324979062714087468 },
    { 0.0, 0.0, -1.0, 0.0, 1.324979062714087468, -0.76497701852859643746, 0.0, 1.5299540370571928749 },
    { -3.0, 0.0, 0.0, 0.0, 0.99616644979174535367, -0.99616644979174535367, 0.99616644979174535367, 0.99616644979174535367 },
    { 0.0, 0.0, 0.0, 1.0, 1.4778221174287648292, -0.3959812429218893245, -0.3959812429218893245, 1.4778221174287648292 },
    { 1.0, 1.0, 0.0, 0.0, 1.6675252249896208332, -0.33169139016798233377, 0.33169139016798233377, 1.6675252249896208332 },
    { 2.0, 0.0, -3.0, 0.0, 0.0, -1.1972208897783685677, 1.1751406146397711734, -0.59861044488918428387 },
    { 12.0, 0.0, -8.0000001, 0.0, 0.0, -0.90689968132986962532, 3.602537212144561362, -0.45344984066493481266 },
};

static Expr* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool num_parts(const Expr* r, double* re, double* im) {
    if (r->type == EXPR_REAL) { *re = r->data.real; *im = 0.0; return true; }
    if (r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
        r->data.function.head->type == EXPR_SYMBOL &&
        strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
        r->data.function.args[0]->type == EXPR_REAL &&
        r->data.function.args[1]->type == EXPR_REAL) {
        *re = r->data.function.args[0]->data.real;
        *im = r->data.function.args[1]->data.real;
        return true;
    }
    return false;
}

static void eval_num(const char* input, double* re, double* im) {
    Expr* r = eval_str(input);
    ASSERT_MSG(num_parts(r, re, im), "%s: expected a machine number", input);
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

static const char* C(double re, double im, char* buf) {
    char a[64], b[64];
    if (im == 0.0) snprintf(buf, 160, "%s", R(re, a));
    else snprintf(buf, 160, "Complex[%s, %s]", R(re, a), R(im, b));
    return buf;
}

#define WP_TOL 1.5e-14     /* 64 ulp; measured max 33 ulp over 400 random points */

void test_wp_exact() {
    assert_true("(WeierstrassP[u, {0, 0}]) === (u^(-2))");
    assert_true("(WeierstrassPPrime[u, {0, 0}]) === (-2/u^3)");
    assert_true("(WeierstrassZeta[u, {0, 0}]) === (1/u)");
    assert_true("(WeierstrassSigma[u, {0, 0}]) === (u)");
    /* a double root: e = 1 > 0 (hyperbolic) and e = -1/2 < 0 (trigonometric) */
    assert_true("(WeierstrassP[u, {12, -8}]) === (1 + 3 Csch[Sqrt[3] u]^2)");
    assert_true("(WeierstrassPPrime[u, {12, -8}]) === (-6 Sqrt[3] Coth[Sqrt[3] u] Csch[Sqrt[3] u]^2)");
    assert_true("(WeierstrassZeta[u, {12, -8}]) === (-u + Sqrt[3] Coth[Sqrt[3] u])");
    assert_true("(WeierstrassP[u, {3, 1}]) === (1 + 3/2 Cot[Sqrt[3/2] u]^2)");
    assert_true("(WeierstrassZeta[u, {3, 1}]) === (u/2 + Sqrt[3/2] Cot[Sqrt[3/2] u])");
    assert_true("WeierstrassSigma[u, {3, 1}] === E^(u^2/4) u Sinc[Sqrt[3/2] u]");
    /* parity */
    assert_true("(WeierstrassP[-u, {g2, g3}]) === (WeierstrassP[u, {g2, g3}])");
    assert_true("(WeierstrassPPrime[-u, {g2, g3}]) === (-WeierstrassPPrime[u, {g2, g3}])");
    assert_true("(WeierstrassSigma[-2 u, {g2, g3}]) === (-WeierstrassSigma[2 u, {g2, g3}])");
    assert_true("(WeierstrassP[u, {g2, g3}]) === (WeierstrassP[u, {g2, g3}])");
    assert_true("(WeierstrassP[1, {1, 2}]) === (WeierstrassP[1, {1, 2}])");
    /* the pole at 0 */
    test_mute_messages();
    assert_true("(WeierstrassP[0, {1, 2}]) === (ComplexInfinity)");
    assert_true("(WeierstrassZeta[0., {1., 2.}]) === (ComplexInfinity)");
    assert_true("(WeierstrassSigma[0, {g2, g3}]) === (0)");
    assert_true("(WeierstrassSigma[0., {1., 2.}]) === (0.)");
    assert_true("NumericQ[WeierstrassP[1, {1, 2}]]");
    assert_true("!NumericQ[WeierstrassP[u, {1, 2}]]");
    /* threading over a List u */
    assert_true("WeierstrassP[{0.1, 0.2}, {1., 2.}] === {WeierstrassP[0.1, {1., 2.}], WeierstrassP[0.2, {1., 2.}]}");
    /* invariants of the square and hexagonal lattices; degenerate half-periods */
    assert_true("(WeierstrassInvariants[{1, I}]) === ({Gamma[1/4]^8/(256 Pi^2), 0})");
    assert_true("(WeierstrassInvariants[{2, 2 I}]) === ({Gamma[1/4]^8/(4096 Pi^2), 0})");
    assert_true("(WeierstrassInvariants[{1, Exp[I Pi/3]}]) === ({0, Gamma[1/3]^18/(4096 Pi^6)})");
    assert_true("(WeierstrassInvariants[{w1, w3}]) === (WeierstrassInvariants[{w1, w3}])");
    assert_true("(WeierstrassHalfPeriods[{3, 1}]) === ({Pi/Sqrt[6], DirectedInfinity[I]})");
    assert_true("(WeierstrassHalfPeriods[{12, -8}]) === ({-I Pi/(2 Sqrt[3]), Infinity})");
    assert_true("(WeierstrassHalfPeriods[{0, 0}]) === ({Infinity, DirectedInfinity[I]})");
    assert_true("(WeierstrassHalfPeriods[{1, 2}]) === (WeierstrassHalfPeriods[{1, 2}])");
    assert_eval_eq("Check[WeierstrassInvariants[{1., 2.}], \"nind\"]", "\"nind\"", 0);
    assert_eval_eq("Check[WeierstrassP[0.5], \"argr\"]", "\"argr\"", 0);
}

void test_wp_machine() {
    char buf[600], us[160], g2s[160], g3s[160];
    for (size_t i = 0; i < sizeof(WP_REF) / sizeof(WP_REF[0]); i++) {
        const double* c = WP_REF[i];
        snprintf(buf, sizeof buf, "%s[%s, {%s, %s}]", NAMES[(int)c[0]], C(c[1], c[2], us),
                 C(c[3], c[4], g2s), C(c[5], c[6], g3s));
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(hypot(re - c[7], im - c[8]) <= WP_TOL * hypot(c[7], c[8]),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", buf, c[7], c[8], re, im);
    }
    assert_true("MachineNumberQ[WeierstrassP[0.3, {1., 2.}]]");
    assert_true("Head[WeierstrassP[0.3, {1., 2. + I}]] === Complex");
    assert_true("Head[WeierstrassP[0.3, {1, 2.}]] === Real");
}

void test_wp_halfperiods() {
    char buf[400], g2s[160], g3s[160];
    for (size_t i = 0; i < sizeof(WP_HP) / sizeof(WP_HP[0]); i++) {
        const double* c = WP_HP[i];
        snprintf(buf, sizeof buf, "WeierstrassHalfPeriods[{%s, %s}]", C(c[0], c[1], g2s), C(c[2], c[3], g3s));
        Expr* r = eval_str(buf);
        ASSERT_MSG(r->type == EXPR_FUNCTION && r->data.function.arg_count == 2, "%s: not a pair", buf);
        double a, b, d, e;
        ASSERT_MSG(num_parts(r->data.function.args[0], &a, &b) && num_parts(r->data.function.args[1], &d, &e),
                   "%s: not numeric", buf);
        /* the nearly degenerate row is ill-conditioned: Delta ~ 1e-6 g2^3 */
        double tol = (c[2] == -8.0000001) ? 1e-8 : 1e-14;
        ASSERT_MSG(hypot(a - c[4], b - c[5]) <= tol * hypot(c[4], c[5]) &&
                   hypot(d - c[6], e - c[7]) <= tol * hypot(c[6], c[7]),
                   "%s: expected {%.17g + %.17g I, %.17g + %.17g I}, got {%.17g + %.17g I, %.17g + %.17g I}",
                   buf, c[4], c[5], c[6], c[7], a, b, d, e);
        expr_free(r);
    }
    /* the invariants of the half-periods give the invariants back */
    assert_true("With[{g = {1.3, -0.7 + 2.1 I}}, Max[Abs[WeierstrassInvariants[WeierstrassHalfPeriods[g]] - g]] < 10^-13]");
    assert_true("With[{g = {-2.5, 4.}}, Max[Abs[WeierstrassInvariants[WeierstrassHalfPeriods[g]] - g]] < 10^-13]");
    assert_true("Max[Abs[WeierstrassInvariants[{1., 0.5 + 1. I}] - {4.540098542453831, 8.385379356790661}]] < 10^-13");
    /* the half-periods are where wp' vanishes and wp takes the roots */
    assert_true("With[{h = WeierstrassHalfPeriods[{1., 2.}]}, Abs[WeierstrassPPrime[h[[1]], {1., 2.}]] < 10^-12 && Abs[4 WeierstrassP[h[[2]], {1., 2.}]^3 - WeierstrassP[h[[2]], {1., 2.}] - 2] < 10^-12]");
    /* arbitrary precision keeps the convention */
    assert_true("With[{h = N[WeierstrassHalfPeriods[{1, 2}], 30]}, Precision[h[[1]]] >= 29 && Abs[h[[1]] - 1.308364762637097`16] < 10^-15 && Abs[h[[2]] - (0.6541823813185484`16 + 1.229368751632915`16 I)] < 10^-15]");
}

void test_wp_arbitrary() {
    assert_true("With[{v = N[WeierstrassP[3/10, {1, 2}], 40]}, Abs[v - 11.1161903563069050162033093046098488646174819`45] < 10^-38 && Precision[v] >= 39]");
    assert_true("With[{v = WeierstrassSigma[7/10 + I/5, {13/10, 21/10 - I/2}] // N[#, 40]&}, Abs[v - (0.699842825993056944224805335385970816537316635`45 + 0.19863401904010130665393498505367739002560047`45 I)] < 10^-38]");
    assert_true("With[{v = N[WeierstrassZeta[7/10 + I/5, {13/10, 21/10 - I/2}], 40]}, Abs[v - (1.31389215959892916747090973532479286865340374`45 - 0.386470803139063923745292931204094268139368602`45 I)] < 10^-38]");
    /* MPFR in, MPFR out; the machine value agrees */
    assert_true("Abs[WeierstrassPPrime[0.5`30, {1, 2}] - WeierstrassPPrime[0.5, {1., 2.}]] < 10^-13");
    assert_true("Precision[WeierstrassZeta[0.5`30, {1, 2}]] >= 29");
}

void test_wp_identities() {
    assert_true("With[{u = 0.37 + 0.21 I, g = {1.3, -0.7 + 2.1 I}}, Abs[WeierstrassPPrime[u, g]^2 - (4 WeierstrassP[u, g]^3 - g[[1]] WeierstrassP[u, g] - g[[2]])] < 10^-12 Abs[WeierstrassPPrime[u, g]]^2]");
    /* periodicity and quasi-periodicity */
    assert_true("With[{g = {1., 2.}, h = WeierstrassHalfPeriods[{1., 2.}]}, Abs[WeierstrassP[0.4 + 2 h[[1]] + 6 h[[2]], g] - WeierstrassP[0.4, g]] < 10^-12]");
    assert_true("With[{g = {1., 2.}, w = WeierstrassHalfPeriods[{1., 2.}][[1]]}, Abs[WeierstrassZeta[0.3 + 2 w, g] - WeierstrassZeta[0.3, g] - 2 WeierstrassZeta[w, g]] < 10^-12]");
    assert_true("With[{g = {1., 2.}, w = WeierstrassHalfPeriods[{1., 2.}][[1]], u = 0.3}, Abs[WeierstrassSigma[u + 2 w, g] + Exp[2 WeierstrassZeta[w, g] (u + w)] WeierstrassSigma[u, g]] < 10^-12 Abs[WeierstrassSigma[u + 2 w, g]]]");
    /* homogeneity: wp(l u; l^-4 g2, l^-6 g3) = l^-2 wp(u; g2, g3) */
    assert_true("With[{l = 1.7 - 0.4 I}, Abs[WeierstrassP[l 0.5, {l^-4 1., l^-6 2.}] - l^-2 WeierstrassP[0.5, {1., 2.}]] < 10^-13 Abs[WeierstrassP[0.5, {1., 2.}]]]");
}

void test_wp_derivative() {
    assert_true("(D[WeierstrassP[u, {g2, g3}], u]) === (WeierstrassPPrime[u, {g2, g3}])");
    assert_true("(D[WeierstrassPPrime[u, {g2, g3}], u]) === (-g2/2 + 6 WeierstrassP[u, {g2, g3}]^2)");
    assert_true("(D[WeierstrassZeta[u, {g2, g3}], u]) === (-WeierstrassP[u, {g2, g3}])");
    assert_true("(D[WeierstrassSigma[u, {g2, g3}], u]) === (WeierstrassSigma[u, {g2, g3}] WeierstrassZeta[u, {g2, g3}])");
    assert_true("(D[WeierstrassP[x^2, {g2, g3}], x]) === (2 x WeierstrassPPrime[x^2, {g2, g3}])");
    /* the invariant derivatives against central differences at 40 digits */
    const char* f[4] = { "WeierstrassP", "WeierstrassPPrime", "WeierstrassZeta", "WeierstrassSigma" };
    for (int k = 0; k < 4; k++) {
        char buf[700];
        snprintf(buf, sizeof buf,
            "With[{d2 = D[%s[u, {a, b}], a] /. {u -> 7/10 + I/5, a -> 13/10, b -> 21/10 - I/2}, "
            "d3 = D[%s[u, {a, b}], b] /. {u -> 7/10 + I/5, a -> 13/10, b -> 21/10 - I/2}, h = 10^-15}, "
            "Abs[N[d2, 40] - N[(%s[7/10 + I/5, {13/10 + h, 21/10 - I/2}] - %s[7/10 + I/5, {13/10 - h, 21/10 - I/2}])/(2 h), 40]] < 10^-25 && "
            "Abs[N[d3, 40] - N[(%s[7/10 + I/5, {13/10, 21/10 - I/2 + h}] - %s[7/10 + I/5, {13/10, 21/10 - I/2 - h}])/(2 h), 40]] < 10^-25]",
            f[k], f[k], f[k], f[k], f[k], f[k]);
        assert_true(buf);
    }
}

void test_wp_buffers_compile() {
    for (int k = 0; k < 4; k++) {
        char buf[500];
        snprintf(buf, sizeof buf, "With[{v = Range[-6., 6., 0.01]}, PackedArrayQ[%s[v, {1., 2.}]] && %s[v, {1., 2.}] === Table[%s[x, {1., 2.}], {x, FromPackedArray[v]}]]",
                 NAMES[k], NAMES[k], NAMES[k]);
        assert_true(buf);
        snprintf(buf, sizeof buf, "With[{v = Range[-3., 3., 0.01] + 0.4 I}, %s[v, {1.5, -2. + I}] === Table[%s[x, {1.5, -2. + I}], {x, FromPackedArray[v]}]]",
                 NAMES[k], NAMES[k]);
        assert_true(buf);
        snprintf(buf, sizeof buf, "With[{f = Compile[{{v, _Real, 1}}, %s[v, {1., 2.}]], v = Range[-5., 5., 0.01]}, f[v] === %s[v, {1., 2.}]]",
                 NAMES[k], NAMES[k]);
        assert_true(buf);
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{x, _Real}}, %s[x, {1., 2.}]][[1]]", NAMES[k]);
        assert_eval_eq(buf, "\"Compiled\" -> True", 0);
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{v, _Real, 1}}, %s[v, {1., 2.}]][[1]]", NAMES[k]);
        assert_eval_eq(buf, "\"Compiled\" -> True", 0);
    }
    assert_true("With[{v = Range[-3., 3., 0.01]}, Normal[WeierstrassP[NDArray[v], {1., 2.}]] === WeierstrassP[v, {1., 2.}]]");
    assert_true("With[{f = Compile[{{x, _Real}}, WeierstrassZeta[x, {1., 2. + I}]]}, f[0.3] === WeierstrassZeta[0.3, {1., 2. + I}]]");
    assert_true("With[{f = Compile[{{x, _Complex}}, WeierstrassSigma[x, {-1., 0.5}]]}, f[0.3 + 2. I] === WeierstrassSigma[0.3 + 2. I, {-1., 0.5}]]");
}

void test_wp_attributes() {
    for (int k = 0; k < 4; k++) {
        char buf[200];
        snprintf(buf, sizeof buf, "Attributes[%s]", NAMES[k]);
        assert_eval_eq(buf, "{NumericFunction, Protected}", 0);
        SymbolDef* d = symtab_get_def(NAMES[k]);
        ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, NAMES[k]));
    }
    assert_true("(Attributes[WeierstrassHalfPeriods]) === ({Protected})");
    assert_true("(Attributes[WeierstrassInvariants]) === ({Protected})");
}

int main() {
    symtab_init();
    core_init();
    TEST(test_wp_exact);
    TEST(test_wp_machine);
    TEST(test_wp_halfperiods);
    TEST(test_wp_arbitrary);
    TEST(test_wp_identities);
    TEST(test_wp_derivative);
    TEST(test_wp_buffers_compile);
    TEST(test_wp_attributes);
    printf("All Weierstrass tests passed.\n");
    return 0;
}
