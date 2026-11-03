/* Tests for the raw-double Bessel / Airy kernels (src/special_functions/
 * sf_machine_bessel.c) and their wiring into the builtins.
 *
 * Covers:
 *   - spot values of every kernel (real and complex, all regimes: Temme /
 *     tabulated-node / Hankel K, series / Miller I, the cut, negative and
 *     non-integer orders; Airy tables, asymptotic tails, complex sheets)
 *     against mpmath at 30 digits, within a few ulps;
 *   - regime boundaries (no jump across x = 1/2, 2, 8, 32 for K; -10, 0, 2,
 *     10 for Airy; |z| = 2, 10, 21 for the complex core);
 *   - declines (underflow / overflow / complex order) and that the builtin
 *     then still answers through the MPFR path;
 *   - the five surfaces agreeing bit for bit:
 *       scalar == unpacked List == packed List == NDArray == Compile;
 *   - the symbolic rules that moved from src/internal/bessel.m into C.
 *
 * Reference values: mpmath 1.3 (mp.besselj / bessely / besseli / besselk,
 * airyai / airybi, 30 digits), rounded to double.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "special_functions/sf_machine_bessel.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* |got - want| in ulps of |want| (complex: of the modulus). */
static double ulps2(double gr, double gi, double wr, double wi) {
    double m = hypot(wr, wi);
    double u = ldexp(1.0, ilogb(m) - 52);
    return hypot(gr - wr, gi - wi) / u;
}

typedef struct { int kind; double nu, zr, zi; double wr, wi; } BCase;

static const BCase BESSEL_CASES[] = {
    { SF_BESSEL_J, 2.0, 3.7, 0, 0.42832965620657587, 0.0 },
    { SF_BESSEL_J, 2.0, 3.0, 1.5, 0.8510031487684113, 0.05963392444703343 },
    { SF_BESSEL_J, 0.0, 2.0, -0.5, 0.21560017471888615, 0.2967860645031962 },
    { SF_BESSEL_J, 2.5, 3.7, 0, 0.45685188411295335, 0.0 },
    { SF_BESSEL_J, 0.5, -2.0, 0, 0.0, 0.5130161365618278 },
    { SF_BESSEL_J, -3.61, 5.6, 0, 0.2285242021053817, 0.0 },
    { SF_BESSEL_J, 1.0, 15.0, 4.0, 5.353057628142319, -1.2610742191900528 },
    { SF_BESSEL_J, 2.0, 25.0, 3.0, -1.1246659853021834, -1.1162748491737255 },
    { SF_BESSEL_Y, 1.0, 3.7, 0, 0.41667437268380747, 0.0 },
    { SF_BESSEL_Y, 1.0, 3.0, 1.5, 0.8100138711450333, 0.47424284380058346 },
    { SF_BESSEL_Y, 1.3, 3.7, 0, 0.37147155857107167, 0.0 },
    { SF_BESSEL_Y, 1.0, -3.0, 0, -0.3246744247918, -0.678117917051873 },
    { SF_BESSEL_Y, 0.0, 7.5, 0.2, 0.12036234253815103, 0.05213860263826909 },
    { SF_BESSEL_Y, 2.0, 30.0, -2.0, 0.4692259666096587, -0.2675851659654811 },
    { SF_BESSEL_I, 0.0, 3.7, 0, 8.738617524169397, 0.0 },
    { SF_BESSEL_I, 0.0, 1.0, 2.0, 0.18785372808246173, 0.6461694351539807 },
    { SF_BESSEL_I, 0.0, 0.3, 0, 1.022626879351597, 0.0 },
    { SF_BESSEL_I, 1.0, 40.0, 0, 1.4707396163259352e+16, 0.0 },
    { SF_BESSEL_I, 0.4, 3.7, 0, 8.507933616049828, 0.0 },
    { SF_BESSEL_I, 0.3333333333333333, -2.0, 0, 1.0793912906864316, 1.8695605567162463 },
    { SF_BESSEL_I, 2.0, 2.5, -2.5, -1.4966843189724506, -0.9438637778901421 },
    { SF_BESSEL_I, 0.0, 0.5, 9.0, -0.10539897677933792, 0.1274853001034855 },
    { SF_BESSEL_K, 0.0, 3.7, 0, 0.01563065992162666, 0.0 },
    { SF_BESSEL_K, 0.0, 1.0, 2.0, -0.242345104491872, -0.17626718909269973 },
    { SF_BESSEL_K, 0.0, 0.2, 0, 1.7527038555281458, 0.0 },
    { SF_BESSEL_K, 0.0, 40.0, 0, 8.392861100099567e-19, 0.0 },
    { SF_BESSEL_K, 0.7, 3.7, 0, 0.0165811853935738, 0.0 },
    { SF_BESSEL_K, 0.0, -2.0, 0, 0.11389387274953344, -7.161528439050256 },
    { SF_BESSEL_K, 3.0, 2.2, 1.1, -0.0987385701977924, -0.33103727060630894 },
    { SF_BESSEL_K, 1.0, 9.9, 0.5, 1.7897042103826625e-05, -1.0413437605313617e-05 },
    { SF_BESSEL_K, 0.0, 10.1, 3.0, -1.568904493322715e-05, -5.252370886162643e-09 },
    { SF_BESSEL_K, 0.0, 25.0, 1.0, 1.8129118262617387e-12, -2.9503149232990907e-12 },
};

/* Real Airy: x, Ai, Ai', Bi, Bi'. */
static const double AIRY_REAL[][5] = {
    { -12.5, -0.27627456138116024, -0.41933133041950515, 0.11703336725739277, -0.974516536167174 },
    { -7.3, 0.3357703705151473, -0.18009580448329365, 0.07087411376989647, 0.9099842704363246 },
    { -2.2, 0.09614537800766888, 0.6862448249090017, -0.45036098416820736, 0.09622918593856436 },
    { 0.0, 0.3550280538878172, -0.2588194037928068, 0.6149266274460007, 0.4482883573538264 },
    { 0.6, 0.20980006166637946, -0.21279325938915852, 0.9110633416949405, 0.5931444786342857 },
    { 1.7, 0.05432479273291947, -0.07737488952532504, 2.3194075069389246, 2.555849356900438 },
    { 2.0, 0.03492413042327438, -0.05309038443365363, 3.2980949999782148, 4.10068204993289 },
    { 2.9, 0.007886312304121233, -0.014042089387786427, 11.942554067750535, 19.097832882836588 },
    { 9.99, 1.1405176956374916e-10, -3.632831449485575e-10, 441574174.82229245, 1384401415.5170975 },
    { 10.0, 1.1047532552898686e-10, -3.5206336767389237e-10, 455641153.54822516, 1429236134.4828658 },
    { 13.1, 2.769285274543043e-15, -1.007530290857721e-14, 15879872826068.336, 57168309370941.07 },
};

/* Complex Airy: zr, zi, then (re, im) of Ai, Ai', Bi, Bi'. */
static const double AIRY_CPLX[][10] = {
    { 0.5, 0.3, 0.22634795458107734, -0.06800141109668117, -0.23013706202248152, 0.03652315800475668, 0.8354902757584077, 0.1583392688855427, 0.4944923355194312, 0.12238319472309553 },
    { 1.5, -0.8, 0.03687042771720862, 0.07099319041852437, -0.07005485418925331, -0.08796483595181803, 1.1029169395397445, -1.1470501312854517, 0.6226085081097154, -1.6507115892270638 },
    { -2.2, 1.0, 0.2928603069467843, 0.9415064636720185, 1.48895456356486, -0.6059380022287009, -1.0303315671937539, 0.24753594768323053, 0.6404785874036298, 1.3312562039848255 },
    { -3.0, 0.0001, -0.3788142993598725, 3.145837714208748e-05, 0.31458377582942687, 0.00011364428856665761, -0.19828962934927094, -6.75611225732832e-05, -0.6756112318279788, 5.9486888435116144e-05 },
    { 0.5, -2.5, -0.7330217556442689, 0.33455739531309947, 0.6238253430890084, -1.0972445357204907, 0.27967861068084815, 0.622903013979301, -1.2859256528307605, -0.6983732639787367 },
    { 2.5, 1.0, -0.0019120892713783827, -0.018032905765349206, -0.0018792086096351543, 0.03102762428412435, 0.5125437840170126, 5.334995579271627, -1.2050484049806622, 8.290971678096389 },
    { -1.0, 0.2, 0.5463063011122739, -0.0027654937674220655, -0.021218909606911092, -0.10779467882663936, 0.1061635423358261, 0.1191267506719781, 0.6021543912932872, -0.022526466310674852 },
    { -6.0, -4.0, 1055.0261137598354, -3527.2741500363663, 8209.045127708676, 5330.393722476991, -3527.274163465412, -1055.0261048801983, 5330.393734646752, -8209.045085778602 },
    { 8.0, 3.0, -7.058972727591045e-08, -7.31731810800639e-08, 1.6747620674528355e-07, 2.485534123454018e-07, -297016.05311856186, 445697.65842548985, -1083279.778960837, 1111317.904024112 },
};

/* ---- kernel spot values --------------------------------------------- */

void test_bessel_kernels(void) {
    size_t n = sizeof(BESSEL_CASES) / sizeof(BESSEL_CASES[0]);
    for (size_t i = 0; i < n; i++) {
        const BCase* c = &BESSEL_CASES[i];
        double rr, ri;
        bool is_real;
        bool ok = sf_bessel_scalar(c->kind, c->nu, c->zr, c->zi, &rr, &ri, &is_real);
        ASSERT_MSG(ok, "kind %d nu %g z %g%+gi declined", c->kind, c->nu, c->zr, c->zi);
        double e = ulps2(rr, ri, c->wr, c->wi);
        ASSERT_MSG(e <= 8.0, "kind %d nu %g z %g%+gi: got %.17g%+.17gi want %.17g%+.17gi (%.1f ulp)",
                   c->kind, c->nu, c->zr, c->zi, rr, ri, c->wr, c->wi, e);
    }
}

void test_airy_kernels(void) {
    for (size_t i = 0; i < sizeof(AIRY_REAL) / sizeof(AIRY_REAL[0]); i++) {
        for (int w = 0; w < 4; w++) {
            double rr, ri; bool is_real;
            ASSERT(sf_airy_scalar(w, AIRY_REAL[i][0], 0.0, &rr, &ri, &is_real));
            ASSERT(is_real && ri == 0.0);
            double e = ulps2(rr, 0.0, AIRY_REAL[i][1 + w], 0.0);
            ASSERT_MSG(e <= 4.0, "Airy[%d](%g): %.17g vs %.17g (%.1f ulp)", w,
                       AIRY_REAL[i][0], rr, AIRY_REAL[i][1 + w], e);
        }
    }
    for (size_t i = 0; i < sizeof(AIRY_CPLX) / sizeof(AIRY_CPLX[0]); i++) {
        for (int w = 0; w < 4; w++) {
            double rr, ri; bool is_real;
            ASSERT(sf_airy_scalar(w, AIRY_CPLX[i][0], AIRY_CPLX[i][1], &rr, &ri, &is_real));
            double e = ulps2(rr, ri, AIRY_CPLX[i][2 + 2 * w], AIRY_CPLX[i][3 + 2 * w]);
            ASSERT_MSG(e <= 16.0, "Airy[%d](%g%+gi): %.17g%+.17gi (%.1f ulp)", w,
                       AIRY_CPLX[i][0], AIRY_CPLX[i][1], rr, ri, e);
        }
    }
}

/* ---- regime boundaries: no jump across a switch point ----------------- */

static void no_jump_real(int kind, double nu, double x) {
    double a, b, ai, bi; bool r;
    double xl = nextafter(x, -INFINITY), xr = nextafter(x, INFINITY);
    ASSERT(sf_bessel_scalar(kind, nu, xl, 0.0, &a, &ai, &r));
    ASSERT(sf_bessel_scalar(kind, nu, xr, 0.0, &b, &bi, &r));
    /* two adjacent doubles: the function moves by ~ x f'/f ulps at most */
    ASSERT_MSG(fabs(a - b) <= 8.0 * ldexp(1.0, ilogb(fabs(a)) - 52) + fabs(a) * 4e-16 * fabs(x),
               "kind %d nu %g: jump at x = %g (%.17g vs %.17g)", kind, nu, x, a, b);
}

void test_regime_boundaries(void) {
    double kx[] = { 0.5, 1.0, 2.0, 8.0, 32.0 };
    for (int i = 0; i < 5; i++) { no_jump_real(SF_BESSEL_K, 0.0, kx[i]); no_jump_real(SF_BESSEL_K, 1.0, kx[i]); }
    double ix[] = { 2.0, 8.0, 32.0 };
    for (int i = 0; i < 3; i++) { no_jump_real(SF_BESSEL_I, 0.0, ix[i]); no_jump_real(SF_BESSEL_I, 1.0, ix[i]); }
    /* Airy: the step between two adjacent doubles must match the derivative
     * (Ai' for Ai, x Ai for Ai', Bi' for Bi, x Bi for Bi') to a few ulps --
     * at x = -10 that step is already ~100 ulps of Ai itself. */
    double ax[] = { -10.0, 0.0, 2.0, 10.0 };
    for (int i = 0; i < 4; i++) {
        double xl = nextafter(ax[i], -INFINITY), xr = nextafter(ax[i], INFINITY);
        for (int w = 0; w < 4; w++) {
            double a, b, t, f, d; bool r;
            ASSERT(sf_airy_scalar(w, xl, 0.0, &a, &t, &r));
            ASSERT(sf_airy_scalar(w, xr, 0.0, &b, &t, &r));
            if (w == 0 || w == 2) ASSERT(sf_airy_scalar(w + 1, ax[i], 0.0, &d, &t, &r));
            else { ASSERT(sf_airy_scalar(w - 1, ax[i], 0.0, &f, &t, &r)); d = ax[i] * f; }
            double slack = 8.0 * ldexp(1.0, ilogb(fmax(fabs(a), fabs(d))) - 52);
            ASSERT_MSG(fabs((b - a) - d * (xr - xl)) <= slack,
                       "Airy[%d] jumps at %g: %.17g vs %.17g", w, ax[i], a, b);
        }
    }
    /* complex core switch radii: the step across each switch must match the
     * derivative from the order recurrences (DLMF 10.6.2, 10.29.2):
     * J' = (J_{n-1} - J_{n+1})/2 (same for Y), I' = (I_{n-1} + I_{n+1})/2,
     * K' = -(K_{n-1} + K_{n+1})/2.  A step of 1e-9 makes the first-order
     * prediction exact to ~1e-18, far below a jump between two methods. */
    double rad[] = { 1.25, 2.0, 10.0, 21.0 };
    for (int i = 0; i < 4; i++) {
        for (int kind = 0; kind < 4; kind++) {
            for (int q = 0; q < 3; q++) {
                double th = 0.15 + 0.6 * q, re1, im1, re2, im2, am, bm, ap, bp, t; bool r;
                double h = 1e-9, r1 = rad[i] - h, r2 = rad[i] + h;
                double c = cos(th), sn = sin(th);
                ASSERT(sf_bessel_scalar(kind, 1.0, r1 * c, r1 * sn, &re1, &im1, &r));
                ASSERT(sf_bessel_scalar(kind, 1.0, r2 * c, r2 * sn, &re2, &im2, &r));
                ASSERT(sf_bessel_scalar(kind, 0.0, rad[i] * c, rad[i] * sn, &am, &bm, &r));
                ASSERT(sf_bessel_scalar(kind, 2.0, rad[i] * c, rad[i] * sn, &ap, &bp, &r));
                double sg = (kind == SF_BESSEL_I) ? 1.0 : (kind == SF_BESSEL_K) ? -1.0 : -1.0;
                double dr = (kind == SF_BESSEL_K) ? -0.5 * (am + ap) : 0.5 * (am + sg * ap);
                double di = (kind == SF_BESSEL_K) ? -0.5 * (bm + bp) : 0.5 * (bm + sg * bp);
                /* predicted step: f' * dz, dz = 2h e^{i th} */
                double pr = (dr * c - di * sn) * 2.0 * h, pi = (dr * sn + di * c) * 2.0 * h;
                double err = hypot((re2 - re1) - pr, (im2 - im1) - pi);
                (void)t;
                ASSERT_MSG(err <= 16.0 * DBL_EPSILON * hypot(re1, im1) + 1e-300,
                           "kind %d jumps at |z| = %g, arg %g (%g)", kind, rad[i], th, err);
            }
        }
    }
}

/* ---- declines, and the builtin still answering ------------------------ */

void test_declines(void) {
    double v, rr, ri; bool r;
    ASSERT(!sf_bessel_k_r(0.0, 800.0, &v));          /* underflows a double */
    ASSERT(!sf_bessel_i_r(0.0, 800.0, &v));          /* overflows */
    ASSERT(!sf_airy_bi_r(120.0, &v));                /* overflows */
    ASSERT(!sf_airy_ai_r(120.0, &v));                /* underflows */
    ASSERT(!sf_bessel_nd(SF_BESSEL_J, 1.0, 0.5, 2.0, 0.0, &rr, &ri));   /* complex order */
    ASSERT(!sf_bessel_nd(SF_BESSEL_J, 0.5, 0.0, -2.0, 0.0, &rr, &ri));  /* real in, complex out */
    ASSERT(!sf_bessel_scalar(SF_BESSEL_K, 0.0, 0.0, 0.0, &rr, &ri, &r)); /* pole */
    /* the builtins fall through to MPFR with an extended exponent */
    assert_eval_startswith("BesselK[0, 800.]", "1.62503083544");
    assert_eval_startswith("AiryBi[120.]", "6.7340932779644");
    assert_eval_startswith("BesselI[0, 800.]", "3.84608161541921");
    assert_eval_startswith("AiryAi[120.]", "2.15749809103");    /* was flushed to 0. */
    assert_eval_eq("BesselK[0, 0]", "Infinity", 0);
}

/* ---- the five surfaces agree bit for bit ------------------------------ */

static void five_surfaces(const char* f, const char* pts) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Module[{xs = FromPackedArray[%s], sc, li, pk, nd, cs, cv},"
        " sc = Map[%s, xs]; li = %s[xs]; pk = FromPackedArray[%s[ToPackedArray[xs]]];"
        " nd = Normal[%s[NDArray[xs]]];"
        " cs = Map[Compile[{{x, _Real}}, %s[x]], xs]; cv = Compile[{{v, _Real, 1}}, %s[v]][xs];"
        " {sc === li, sc === pk, sc === nd, sc === cs, sc === cv}]",
        pts, f, f, f, f, f, f);
    assert_eval_eq(buf, "{True, True, True, True, True}", 0);
}

static void five_surfaces_c(const char* f, const char* pts) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Module[{xs = %s, sc, li, nd, cs, cv},"
        " sc = Map[%s, xs]; li = %s[xs]; nd = Normal[%s[NDArray[xs]]];"
        " cs = Map[Compile[{{x, _Complex}}, %s[x]], xs]; cv = Compile[{{v, _Complex, 1}}, %s[v]][xs];"
        " {sc === li, sc === nd, sc === cs, sc === cv}]",
        pts, f, f, f, f, f);
    assert_eval_eq(buf, "{True, True, True, True}", 0);
}

void test_surfaces(void) {
    const char* xs = "{0.13, 0.77, 1.9, 2.3, 3.7, 6.1, 9.9, 15.2, 33.3, 45.5}";
    const char* ax = "{-12.5, -7.3, -2.2, -0.4, 0., 0.6, 1.7, 2.9, 4.4, 9.7, 13.1}";
    five_surfaces("BesselJ[2, #] &", xs);
    five_surfaces("BesselY[1, #] &", xs);
    five_surfaces("BesselI[0, #] &", xs);
    five_surfaces("BesselK[0, #] &", xs);
    five_surfaces("BesselK[3, #] &", xs);
    five_surfaces("BesselJ[2.5, #] &", xs);
    five_surfaces("BesselK[1/3, #] &", xs);
    five_surfaces("AiryAi", ax);
    five_surfaces("AiryBi", ax);
    five_surfaces("AiryAiPrime", ax);
    five_surfaces("AiryBiPrime", ax);
    const char* zs = "{0.5 + 0.3 I, 1.5 - 0.8 I, -2.2 + 1. I, 0.5 - 2.5 I, 2.5 + 1. I, 7.1 + 3.3 I, 12.5 - 4. I, 30. + 2. I}";
    five_surfaces_c("BesselJ[2, #] &", zs);
    five_surfaces_c("BesselY[1, #] &", zs);
    five_surfaces_c("BesselI[0, #] &", zs);
    five_surfaces_c("BesselK[0, #] &", zs);
    five_surfaces_c("AiryAi", zs);
    five_surfaces_c("AiryBiPrime", zs);
    /* Compile lowers all of them at scalar and rank-1 shapes, real and complex */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, BesselJ[2, x]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, BesselK[0, v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Complex}}, BesselY[1, x]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Complex, 1}}, AiryBi[v]]", "True", 0);
}

/* ---- the scalar builtin answers: Real for real, Complex otherwise ------- */

void test_builtin_shapes(void) {
    assert_eval_eq("Head[BesselJ[2, 3.7]]", "Real", 0);
    assert_eval_eq("Head[BesselJ[1/2, -2.]]", "Complex", 0);
    assert_eval_eq("Head[BesselK[0, -2.]]", "Complex", 0);
    assert_eval_eq("Head[BesselJ[2, 3. + 0.5 I]]", "Complex", 0);
    assert_eval_eq("Head[AiryAi[1. - 2. I]]", "Complex", 0);
    assert_eval_eq("Head[BesselY[1/3, 2. I]]", "Complex", 0);
    assert_eval_eq("Precision[BesselI[0, 2.5]]", "MachinePrecision", 0);
    /* arbitrary precision is untouched (MPFR / arb path) */
    assert_eval_startswith("BesselJ[2, 3.7`30]", "0.42832965620657589712");
    assert_eval_startswith("AiryAi[1.3`25]", "0.09347466577150270986722");
}

/* ---- symbolic rules (formerly src/internal/bessel.m) ------------------- */

void test_symbolic_rules(void) {
    assert_eval_eq("BesselJ[1/2, x]", "Sqrt[2/(Pi x)] Sin[x]", 0);
    assert_eval_eq("BesselJ[-1/2, x]", "Sqrt[2/(Pi x)] Cos[x]", 0);
    assert_eval_eq("BesselY[1/2, x]", "-Sqrt[2/(Pi x)] Cos[x]", 0);
    assert_eval_eq("BesselI[-1/2, x]", "Sqrt[2/(Pi x)] Cosh[x]", 0);
    assert_eval_eq("BesselK[1/2, x]", "Sqrt[Pi/(2 x)]/E^x", 0);
    assert_eval_eq("BesselK[-3/2, x]", "Sqrt[Pi/(2 x)]/E^x + Sqrt[Pi/(2 x)]/(E^x x)", 0);
    assert_eval_eq("BesselJ[5/2, x]",
        "-Sqrt[2/(Pi x)] Sin[x] + (3 (-Sqrt[2/(Pi x)] Cos[x] + (Sqrt[2/(Pi x)] Sin[x])/x))/x", 0);
    assert_eval_eq("BesselJ[-2, x]", "BesselJ[2, x]", 0);
    assert_eval_eq("BesselJ[-3, x]", "-BesselJ[3, x]", 0);
    assert_eval_eq("BesselY[-1, x]", "-BesselY[1, x]", 0);
    assert_eval_eq("BesselI[-2, x]", "BesselI[2, x]", 0);
    assert_eval_eq("BesselK[-2, x]", "BesselK[2, x]", 0);
    /* a numeric argument keeps the half-integer call (Mathematica) */
    assert_eval_eq("BesselJ[11/2, 1]", "BesselJ[11/2, 1]", 0);
    assert_eval_eq("BesselJ[-3, 2]", "-BesselJ[3, 2]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bessel_kernels);
    TEST(test_airy_kernels);
    TEST(test_regime_boundaries);
    TEST(test_declines);
    TEST(test_surfaces);
    TEST(test_builtin_shapes);
    TEST(test_symbolic_rules);

    printf("All sf_machine_bessel tests passed.\n");
    return 0;
}
