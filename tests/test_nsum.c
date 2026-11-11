/* Tests for NSum — numerical summation.
 *
 * Cover: Direct (small finite) sums; infinite sums via Automatic dispatch
 * (Euler–Maclaurin for monotone, Cohen–Villegas–Zagier for alternating, Wynn
 * epsilon fallback); the explicit Method settings; step di; complex summands;
 * long finite sums (direct and verified finite Euler-Maclaurin, never
 * extrapolated); multidimensional sums with
 * independent and dependent inner bounds; VerifyConvergence and divergence;
 * arbitrary precision via WorkingPrecision; options and unevaluated forms; the
 * HoldAll/Protected attributes; and memory hygiene.
 *
 * Expected values come from closed forms or Keiper 1992 "The N functions of
 * Mathematica".  Numerical results are compared *inside* the language
 * (N[Abs[result - expected]] < tol) rather than by parsing the printed form.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

/* |input - expected| < tol, evaluated at full internal precision; the L1
 * magnitude keeps a spurious tiny imaginary residual from leaving it complex. */
static bool close_to(const char* input, const char* expected, double tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "N[Abs[Re[(%s) - (%s)]] + Abs[Im[(%s) - (%s)]]] < %.17g",
             input, expected, input, expected, tol);
    char* s = eval_str(buf);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  close_to FAIL: |%s - (%s)| < %g  =>  %s\n",
                     input, expected, tol, s);
    free(s);
    return ok;
}

#define ASSERT_CLOSE(input, expected, tol)                                   \
    ASSERT_MSG(close_to((input), (expected), (tol)),                         \
               "%s ~= %s (tol %g)", (input), (expected), (tol))

/* ---------------------------------------------------------------------- */

static void test_direct_finite(void) {
    ASSERT_CLOSE("NSum[1/i^2,{i,1,5}]", "5269/3600", 1e-12);
    ASSERT_CLOSE("NSum[i,{i,1,100}]", "5050", 1e-9);
    ASSERT_CLOSE("NSum[1/i,{i,1,10}]", "7381/2520", 1e-12);
    ASSERT_CLOSE("NSum[i^2,{i,3,7}]", "135", 1e-9);
    /* implicit lower bound of 1 */
    ASSERT_CLOSE("NSum[1/i^2,{i,3}]", "49/36", 1e-12);
}

static void test_infinite_monotone(void) {
    ASSERT_CLOSE("NSum[1/i^2,{i,1,Infinity}]", "Pi^2/6", 1e-9);
    ASSERT_CLOSE("NSum[1/i^4,{i,1,Infinity}]", "Pi^4/90", 1e-9);
    ASSERT_CLOSE("NSum[1/n^(11/10),{n,1,Infinity}]", "Zeta[11/10]", 1e-7);
    ASSERT_CLOSE("NSum[1/(x^2+5),{x,0,Infinity}]",
                 "(Csch[Sqrt[5] Pi](5 Pi Cosh[Sqrt[5] Pi]+Sqrt[5] Sinh[Sqrt[5] Pi]))/(10 Sqrt[5])",
                 1e-8);
    /* geometric and exponential-tail series */
    ASSERT_CLOSE("NSum[1/2^i,{i,0,Infinity}]", "2", 1e-10);
    ASSERT_CLOSE("NSum[(-5)^i/i!,{i,0,Infinity},NSumTerms->25]", "Exp[-5]", 1e-12);
    ASSERT_CLOSE("NSum[1/Fibonacci[i],{i,Infinity}]", "3.359885666243178", 1e-6);
}

/* Integer-only summands (Prime, PartitionsP, …) have no continuation off the
 * integers, so Euler–Maclaurin — which samples a continuous tail integral and a
 * complex derivative contour — cannot be applied: it would sample Prime[16.16]
 * (spewing Prime::intpp) and then bail to a weakly-seeded Wynn.  NSum must detect
 * this and route straight to partial-sum extrapolation with enough head terms.
 * References: prime zeta P(2)=0.4522474200…, P(4)=0.0769931397…. */
static void test_integer_only_summand(void) {
    /* Convergent prime-zeta: slow, irregular tail — a few digits, not a wild miss
     * (the pre-fix Euler–Maclaurin fallback gave 0.4485). */
    ASSERT_CLOSE("NSum[1./Prime[n]^2,{n,1,Infinity}]", "0.4522474200410655", 2e-3);
    ASSERT_CLOSE("NSum[1/Prime[n]^4,{n,1,Infinity}]", "0.0769931397642468", 1e-6);
    /* Another integer-only head, fast-converging. */
    ASSERT_CLOSE("NSum[1./(n PartitionsP[n]),{n,1,Infinity}]", "1.4815", 1e-3);
    /* The continuous sibling must still take Euler–Maclaurin and be exact. */
    ASSERT_CLOSE("NSum[1/n^2,{n,1,Infinity}]", "Pi^2/6", 1e-9);
}

static void test_alternating(void) {
    ASSERT_CLOSE("NSum[(-1)^k/(2k+1),{k,0,Infinity}]", "Pi/4", 1e-9);
    ASSERT_CLOSE("NSum[(-1)^k/(k+1),{k,0,Infinity}]", "Log[2]", 1e-9);
    /* peaked alternating summand: the dedicated method is very precise */
    ASSERT_CLOSE("NSum[(-1)^x/(1+(x-12)^2),{x,0,Infinity},Method->\"AlternatingSigns\","
                 "WorkingPrecision->30]",
                 "0.27519385941395303956897156158", 1e-24);
}

static void test_methods_explicit(void) {
    ASSERT_CLOSE("NSum[1/i^2,{i,1,Infinity},Method->EulerMaclaurin]", "Pi^2/6", 1e-10);
    /* Wynn epsilon excels on geometric / alternating tails (not monotone ones). */
    ASSERT_CLOSE("NSum[1/2^i,{i,0,Infinity},Method->WynnEpsilon]", "2", 1e-9);
    ASSERT_CLOSE("NSum[(-1)^k/(2k+1),{k,0,Infinity},Method->WynnEpsilon]", "Pi/4", 1e-5);
    ASSERT_CLOSE("NSum[(-1)^k/(2k+1),{k,0,Infinity},Method->AlternatingSigns]", "Pi/4", 1e-12);
}

static void test_levin(void) {
    /* u-transform (string, bare symbol) on monotone and alternating series. */
    ASSERT_CLOSE("NSum[1/i^2,{i,1,Infinity},Method->\"Levin\"]", "Pi^2/6", 1e-5);
    ASSERT_CLOSE("NSum[1/i^4,{i,1,Infinity},Method->Levin]", "Pi^4/90", 1e-7);
    ASSERT_CLOSE("NSum[(-1)^k/(k+1),{k,0,Infinity},Method->\"LevinU\"]", "Log[2]", 1e-8);
    /* t-transform is near-exact on a pure geometric tail. */
    ASSERT_CLOSE("NSum[1/2^i,{i,0,Infinity},Method->\"LevinT\"]", "2", 1e-10);
    /* MPFR: reaches full precision on a smooth series. */
    ASSERT_CLOSE("NSum[1/i^2,{i,1,Infinity},Method->\"Levin\",WorkingPrecision->30]",
                 "Pi^2/6", 1e-25);
}

static void test_step_di(void) {
    ASSERT_CLOSE("NSum[1/2^i,{i,0,Infinity,2}]", "4/3", 1e-10);
    ASSERT_CLOSE("NSum[1/2^(2j),{j,0,Infinity}]", "4/3", 1e-10);
    ASSERT_CLOSE("NSum[1/i^2,{i,2,Infinity,2}]", "Pi^2/24", 1e-8);
}

static void test_complex_summand(void) {
    /* Keiper: -0.182175 - 0.136618 I */
    ASSERT_CLOSE("NSum[Log[x]/x^(2+2I),{x,1,Infinity}]",
                 "-0.18217548 - 0.13661758 I", 1e-5);
}

static void test_large_finite(void) {
    /* difference of two infinite tails */
    ASSERT_CLOSE("NSum[1/i^2,{i,100,1000000}]",
                 "NSum[1/i^2,{i,100,Infinity}] - NSum[1/i^2,{i,1000000,Infinity}]", 1e-9);
    /* Sum_{100}^inf 1/i^2 = PolyGamma[1,100] */
    ASSERT_CLOSE("NSum[1/i^2,{i,100,Infinity}]", "PolyGamma[1,100]", 1e-9);
    ASSERT_CLOSE("NSum[1/i^2,{i,100,1000000}]", "0.0100492", 1e-6);
}

/* Long finite ranges (v0.6513).  A finite sum has an exact answer and is never
 * extrapolated: direct compensated summation while affordable, then a finite
 * Euler-Maclaurin estimate with the boundary corrections at both ends that is
 * accepted only when two independent splits agree, then direct again.  Through
 * v0.6511 such a sum was the difference of two "infinite" tails whenever the
 * head terms looked non-increasing: for 1/k both tails diverge, and
 * HarmonicNumber[10^6] came back as 3.81593 with only NSum::accgl; a constant
 * summand failed with NSum::ncvg. */
static void test_long_finite(void) {
    /* the three reported cases */
    ASSERT_CLOSE("NSum[1/k,{k,1,10^6}]", "N[HarmonicNumber[10^6],30]", 1e-12);
    ASSERT_CLOSE("NSum[1,{k,1,10^5}]", "100000", 1e-8);
    ASSERT_CLOSE("NSum[1,{k,1,10^7}]", "10^7", 1e-6);
    /* short ranges */
    ASSERT_CLOSE("NSum[1/k,{k,1,1000}]", "N[PolyGamma[0,1001]+EulerGamma,30]", 1e-13);
    /* past the direct budget: the verified finite Euler-Maclaurin path */
    ASSERT_CLOSE("NSum[1/k,{k,1,10^7}]", "N[PolyGamma[0,10^7+1]+EulerGamma,30]", 1e-12);
    ASSERT_CLOSE("NSum[1/k,{k,1,10^8}]", "N[PolyGamma[0,10^8+1]+EulerGamma,30]", 1e-11);
    ASSERT_CLOSE("NSum[1,{k,1,10^12}]", "10^12", 1e-1);
    ASSERT_CLOSE("NSum[Log[k],{k,1,10^7}]", "LogGamma[10000001`30]", 1e-5);
    ASSERT_CLOSE("NSum[1/Sqrt[k],{k,1,10^7}]",
                 "N[Zeta[1/2]-HurwitzZeta[1/2,10^7+1],30]", 1e-9);
    /* polynomial growth over a long range */
    ASSERT_CLOSE("NSum[k^2,{k,10^5,10^6}]/10^17",
                 "(10^6 (10^6+1)(2 10^6+1) - (10^5-1) 10^5 (2 10^5-1))/6/10^17", 1e-13);
    ASSERT_CLOSE("NSum[k^2,{k,10^5,10^7}]/10^20",
                 "(10^7 (10^7+1)(2 10^7+1) - (10^5-1) 10^5 (2 10^5-1))/6/10^20", 1e-12);
}

/* Non-monotone, oscillating and aliased summands over long ranges: the
 * Euler-Maclaurin model does not fit, the smoothness probe or the two-split
 * check rejects it, and direct summation answers. */
static void test_long_finite_oscillating(void) {
    /* sum_{k=1}^{2m} (-1)^k/k = H(m) - H(2m) */
    ASSERT_CLOSE("NSum[(-1)^k/k,{k,1,10^6}]",
                 "N[PolyGamma[0,5 10^5+1]-PolyGamma[0,10^6+1],30]", 1e-13);
    ASSERT_CLOSE("NSum[(-1)^k/k,{k,1,10^7}]",
                 "N[PolyGamma[0,5 10^6+1]-PolyGamma[0,10^7+1],30]", 1e-13);
    /* 40-digit reference sums */
    ASSERT_CLOSE("NSum[Sin[k]/k,{k,1,10^6}]", "1.07079529444189234462687620342038446", 1e-13);
    ASSERT_CLOSE("NSum[Sin[k]/k,{k,1,10^7}]", "1.07079643085964681354032959984767176", 1e-13);
    /* Cos[2 Pi k] = 1 on the integers: an aliased continuous extension */
    ASSERT_CLOSE("NSum[Cos[2 Pi k]/k,{k,1,10^7}]", "N[PolyGamma[0,10^7+1]+EulerGamma,30]", 1e-12);
    /* a unit-width peak in the middle of the range: the sum over all integers
     * is Pi Sinh[2Pi]/(Cosh[2Pi]-Cos[2Pi/3]); the two tails cut off ~4e-7 */
    ASSERT_CLOSE("NSum[1/(1+(k-5*10^6-1/3)^2),{k,1,10^7}]",
                 "N[Pi Sinh[2 Pi]/(Cosh[2 Pi]-Cos[2 Pi/3]),20]", 1e-6);
    /* a broad peak (k e^(-k/10^5) peaks at 10^5): sum_{k=1}^n k r^k in closed form */
    ASSERT_CLOSE("NSum[k Exp[-k/10^5],{k,1,10^7}]/10^10",
                 "N[(Exp[-1/10^5] - (10^7+1) Exp[-(10^7+1)/10^5] + 10^7 Exp[-(10^7+2)/10^5])"
                 "/(1-Exp[-1/10^5])^2/10^10, 30]", 1e-13);
    /* complex summand */
    ASSERT_CLOSE("NSum[(1+I)/k,{k,1,10^7}]", "(1+I) N[PolyGamma[0,10^7+1]+EulerGamma,30]", 1e-11);
}

/* Non-unit and negative steps, negative ranges, the empty range. */
static void test_long_finite_steps(void) {
    ASSERT_CLOSE("NSum[1/k^2,{k,1,10^6,2}]",
                 "N[(PolyGamma[1,1/2]-PolyGamma[1,1/2+5 10^5])/4,30]", 1e-14);
    ASSERT_CLOSE("NSum[1/k,{k,1,3*10^7,3}]",
                 "N[(PolyGamma[0,1/3+10^7]-PolyGamma[0,1/3])/3,30]", 1e-12);
    ASSERT_CLOSE("NSum[1/k,{k,10^7,1,-1}]", "N[PolyGamma[0,10^7+1]+EulerGamma,30]", 1e-12);
    ASSERT_CLOSE("NSum[1/k^2,{k,-10^7,-1}]", "N[Zeta[2]-HurwitzZeta[2,10^7+1],30]", 1e-14);
    ASSERT_CLOSE("NSum[1/k,{k,5,1}]", "0", 1e-300);
}

/* Summands that do not compile (PrimePi, Prime): the interpreter serves
 * every term, the direct budget is smaller, and Euler-Maclaurin still applies
 * where the summand extends off the integers. */
static void test_long_finite_interpreted(void) {
    ASSERT_CLOSE("NSum[1/(k+PrimePi[3]),{k,1,10^6}]",
                 "N[PolyGamma[0,10^6+3]+EulerGamma-3/2,30]", 1e-12);
    ASSERT_CLOSE("NSum[1/Prime[k],{k,1,2000}]", "N[Total[1/Prime[Range[2000]]],30]", 1e-13);
    /* an integer-only summand past the hard cap has no verified method: the
     * call stays unevaluated with NSum::nsum, never an extrapolation */
    char* s = eval_str("Quiet[NSum[1/Prime[k],{k,1,10^9}]]");
    ASSERT_MSG(strncmp(s, "NSum[", 5) == 0, "unevaluated, got %s", s);
    free(s);
    s = eval_str("Check[NSum[1/Prime[k],{k,1,10^9}], failed]");
    ASSERT_MSG(strcmp(s, "failed") == 0, "Check sees NSum::nsum, got %s", s);
    free(s);
}

/* WorkingPrecision -> 30 over a long finite range: compiled MPFR direct sum. */
static void test_long_finite_working_precision(void) {
    ASSERT_CLOSE("NSum[1/k,{k,1,10^5},WorkingPrecision->30]",
                 "N[PolyGamma[0,10^5+1]+EulerGamma,40]", 1e-26);
    char* s = eval_str("Precision[NSum[1/k,{k,1,10^5},WorkingPrecision->30]]");
    ASSERT_MSG(strcmp(s, "30.") == 0 || strcmp(s, "30.0") == 0, "precision 30, got %s", s);
    free(s);
}

/* Timing bound: the long-range paths stay fast (generous for CI load; each
 * takes 2-30 ms on an M-series machine). */
static void test_long_finite_timing(void) {
    const char* cases[] = {
        "First[AbsoluteTiming[NSum[1/k,{k,1,10^8}]]] < 1.5",
        "First[AbsoluteTiming[NSum[1,{k,1,10^7}]]] < 1.5",
        "First[AbsoluteTiming[NSum[1/(k+PrimePi[3]),{k,1,10^6}]]] < 1.5",
        "First[AbsoluteTiming[NSum[1/k,{k,1,10^6}]]] < 1.5",
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char* s = eval_str(cases[i]);
        ASSERT_MSG(strcmp(s, "True") == 0, "%s => %s", cases[i], s);
        free(s);
    }
}

/* Randomised agreement with exact closed forms: fixed seed, summand families
 * with known finite sums, spans 10^3..10^7 (both the direct and the
 * Euler-Maclaurin paths).  Every case must agree to 1e-11 relative. */
static void test_long_finite_random(void) {
    const char* prog =
        "Module[{fam, res = {}},"
        " SeedRandom[6513];"
        " fam = {"
        "  {Function[{k, c}, 1/(k + c)], Function[{a, b, c}, PolyGamma[0, b + 1 + c] - PolyGamma[0, a + c]]},"
        "  {Function[{k, c}, 1/(k + c)^2], Function[{a, b, c}, PolyGamma[1, a + c] - PolyGamma[1, b + 1 + c]]},"
        "  {Function[{k, c}, (k + c)^(-3/2)], Function[{a, b, c}, HurwitzZeta[3/2, a + c] - HurwitzZeta[3/2, b + 1 + c]]},"
        "  {Function[{k, c}, Log[k + c]], Function[{a, b, c}, LogGamma[b + 1 + c] - LogGamma[a + c]]},"
        "  {Function[{k, c}, Sqrt[k + c]], Function[{a, b, c}, HurwitzZeta[-1/2, a + c] - HurwitzZeta[-1/2, b + 1 + c]]},"
        "  {Function[{k, c}, k^2 - c k], Function[{a, b, c}, (b (b + 1) (2 b + 1) - (a - 1) a (2 a - 1))/6 - c (b (b + 1) - (a - 1) a)/2]}};"
        " Do[Module[{f = RandomChoice[fam], c = RandomInteger[{1, 99}]/RandomInteger[{1, 9}],"
        "     a = RandomInteger[{1, 50}], b, v, x},"
        "   b = a + Round[10^RandomReal[{3, 7}]] - 1;"
        "   v = NSum[Evaluate[f[[1]][k, c]], {k, a, b}];"
        "   x = N[f[[2]][N[a, 40], N[b, 40], N[c, 40]], 30];"
        "   AppendTo[res, If[NumberQ[v], N[Abs[SetPrecision[v, 40] - x]/Abs[x]], Infinity]]],"
        "  {24}];"
        " Max[res]]";
    char buf[4096];
    snprintf(buf, sizeof buf, "(%s) < 10^-11", prog);
    char* s = eval_str(buf);
    if (strcmp(s, "True") != 0) {
        char* m = eval_str(prog);
        fprintf(stderr, "  random finite NSum: max rel err %s\n", m);
        free(m);
    }
    ASSERT_MSG(strcmp(s, "True") == 0, "random finite NSum agreement, got %s", s);
    free(s);
}

static void test_multidim(void) {
    /* Keiper multidimensional examples */
    ASSERT_CLOSE("NSum[(-1)^n (2/n)^k/k^2,{n,2,Infinity},{k,1,Infinity}]", "1.14434", 1e-4);
    ASSERT_CLOSE("NSum[(-1)^n (2/n)^k/k^2,{n,2,Infinity},{k,1,n}]", "0.770188", 1e-4);
}

/* Rapidly decaying summands are summed directly (ns_fast_decay): no
 * extrapolation, no phase split.  Each case was either wrong in method choice
 * or orders of magnitude slow before (values: mpmath at 70 digits). */
static void test_fast_decay(void) {
    /* Alternating with terms ~1e-14k: the adjacent-sign product underflowed,
     * so the series was not seen as alternating; the Euler-Maclaurin far probe
     * then evaluated (-1)^x at x = 2^(2^29) and never returned. */
    ASSERT_CLOSE("4/Pi NSum[(-1)^k/(2 k + 1) Sech[5 (2 k + 1) Pi], {k, 0, Infinity}, "
                 "WorkingPrecision -> 60]",
                 "3.837587979251226103407133186204839100793005594072509569030023`65*^-7", 1e-63);
    ASSERT_CLOSE("4/Pi NSum[(-1)^k/(2 k + 1) Sech[5 (2 k + 1) Pi], {k, 0, Infinity}]",
                 "3.837587979251226103407133186204839100793005594072509569030023`65*^-7", 1e-20);
    /* heat eigen-series (super-geometric decay) */
    ASSERT_CLOSE("4/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 (42/100)/4], {k, 0, Infinity}, "
                 "WorkingPrecision -> 60]",
                 "0.451659083985231656057607964274935125243134852153997739921580122`65", 1e-57);
    /* a summand peaking at k = 10 before it decays */
    ASSERT_CLOSE("NSum[10^k/k!, {k, 0, Infinity}, WorkingPrecision -> 40]/N[E^10, 50]", "1", 1e-38);
    /* a fast head hiding a slow tail: the route abandons when the decay slows */
    ASSERT_CLOSE("NSum[2^-k + 10^-30/k^2, {k, 1, Infinity}, WorkingPrecision -> 60]",
                 "1 + 10^-30 Pi^2/6", 1e-55);
    /* isolated zero terms (Sin[j Pi/2] vanishes at every even j) */
    ASSERT_CLOSE("NSum[Sin[j Pi/2] Exp[-j], {j, 0, Infinity}]", "E/(1 + E^2)", 1e-14);
    /* two-dimensional, the inner summand compiled with the outer index as a
     * parameter */
    ASSERT_CLOSE("NSum[(-1)^(j + k)/((2 j + 1) (2 k + 1)) Sech[3/2 Pi Sqrt[(2 j + 1)^2 + (2 k + 1)^2]], "
                 "{j, 0, Infinity}, {k, 0, Infinity}, WorkingPrecision -> 60]",
                 "0.002550782687931167381813295907469973865484160198991840612034191`65", 1e-58);
    ASSERT_CLOSE("16/Pi^2 NSum[(-1)^(j + k)/((2 j + 1) (2 k + 1)) (Exp[-Pi Sqrt[(2 j + 1)^2 + (2 k + 1)^2]/4]"
                 " - Exp[-7 Pi Sqrt[(2 j + 1)^2 + (2 k + 1)^2]/4])/(1 - Exp[-2 Pi Sqrt[(2 j + 1)^2 + (2 k + 1)^2]]), "
                 "{j, 0, Infinity}, {k, 0, Infinity}, WorkingPrecision -> 60]",
                 "0.4580868145806731115172877286282238150718731301039713896853592`65", 1e-55);
}

/* An integer-only summand (BesselJZero has no value off the integers) used to
 * go to Wynn's epsilon, which fails on partial sums already constant to the
 * working precision; a saturated sequence now returns its last partial sum. */
static void test_saturated_sequence(void) {
    ASSERT_CLOSE("1 - NSum[2/(BesselJZero[0, n] BesselJ[1, BesselJZero[0, n]]) Exp[-BesselJZero[0, n]^2 0.15], "
                 "{n, 1, Infinity}]",
                 "0.338165656302720936129493448598455049118926235927033732789354", 1e-13);
    ASSERT_CLOSE("1 - NSum[2/(BesselJZero[0, n] BesselJ[1, BesselJZero[0, n]]) Exp[-BesselJZero[0, n]^2 15/100], "
                 "{n, 1, Infinity}, WorkingPrecision -> 40]",
                 "0.338165656302720936129493448598455049118926235927033732789354`60", 1e-38);
    ASSERT_CLOSE("NSum[Exp[-k^2], {k, 0, Infinity}, Method -> \"WynnEpsilon\"]",
                 "1.386318602413326076515625275578929240671693022686123026915795", 1e-14);
    ASSERT_CLOSE("NSum[Exp[-k^2], {k, 0, Infinity}, Method -> \"WynnEpsilon\", WorkingPrecision -> 30]",
                 "1.386318602413326076515625275578929240671693022686123026915795`60", 1e-28);
    ASSERT_CLOSE("NSum[Exp[-k^2], {k, 0, Infinity}, Method -> \"Levin\"]",
                 "1.386318602413326076515625275578929240671693022686123026915795", 1e-14);
}

static void test_verify_convergence(void) {
    /* divergent: default verification returns ComplexInfinity */
    char* s = eval_str("NSum[2^i,{i,0,Infinity}]");
    ASSERT_MSG(strcmp(s, "ComplexInfinity") == 0, "2^i should be ComplexInfinity: %s", s);
    free(s);
    /* without verification the formal Shanks value -1 is returned */
    ASSERT_CLOSE("NSum[2^i,{i,0,Infinity},VerifyConvergence->False]", "-1", 1e-6);
}

static void test_working_precision(void) {
    ASSERT_CLOSE("NSum[1/n^(11/10),{n,1,Infinity},WorkingPrecision->40]", "Zeta[11/10]", 1e-25);
    ASSERT_CLOSE("NSum[1/n^(51/50),{n,1,Infinity},WorkingPrecision->40]", "Zeta[51/50]", 1e-25);
    ASSERT_CLOSE("NSum[1/i^2,{i,1,Infinity},WorkingPrecision->30]", "Pi^2/6", 1e-25);
    ASSERT_CLOSE("NSum[(-1)^k/(2k+1),{k,0,Infinity},WorkingPrecision->30,Method->AlternatingSigns]",
                 "Pi/4", 1e-26);
}

static void test_peaked_summand(void) {
    /* A summand peaking past the head-probe window (peak at k=20) must not be
     * mis-flagged divergent and must sum correctly.  Regression: the engine only
     * probed the first ~18 terms, so it returned ComplexInfinity (or ~0.085).
     * Reference: an independent direct partial sum (tail ~ 2e-5 at 50000). */
    ASSERT_CLOSE("NSum[1/(1+(k-20)^2),{k,0,Infinity}]",
                 "Total[Table[1.0/(1+(k-20)^2),{k,0,50000}]]", 5e-5);
    char* s = eval_str("NSum[1/(1+(k-20)^2),{k,0,Infinity}]");
    ASSERT_MSG(strcmp(s, "ComplexInfinity") != 0,
               "peaked convergent sum must not be flagged divergent: %s", s);
    free(s);
}

static void test_composite_precision(void) {
    /* Composite (Log-near-1) and rational summands now reach the working
     * precision: WP-scaled tail integral + numerical contour corrections +
     * guard digits against the float-sampling cancellation.  Regression: the
     * Log forms capped at ~8-12 digits and rationals at ~15. */
    ASSERT_CLOSE("NSum[Log[1+1/n^2],{n,1,Infinity},WorkingPrecision->30]",
                 "Log[Sinh[Pi]/Pi]", 1e-25);
    ASSERT_CLOSE("NSum[1/(n^2+1),{n,1,Infinity},WorkingPrecision->30]",
                 "(Pi Coth[Pi]-1)/2", 1e-25);
    ASSERT_CLOSE("NSum[Log[E^(-1/(2n)) (1+1/(2n))],{n,1,Infinity},WorkingPrecision->30]",
                 "Log[1/(E^(EulerGamma/2) Gamma[3/2])]", 1e-24);
}

/* Slowly converging algebraic series at MACHINE precision.  Terms ~ n^-p with
 * p just above 1 put almost all of the sum in the Euler–Maclaurin tail
 * integral, whose remainder ∫_X^∞ ~ X^(1-p)/(p-1) is still ~1e-4 at the
 * largest X a double holds; exp-sinh quadrature lost it while its levels still
 * agreed (the ISSAC '97 problem-3 sum came out 21.192745, 5e-5 relative).
 * Each case must reach 10 significant digits.  References are closed forms
 * (Zeta, Hurwitz zeta), the binomial–Zeta expansion
 *     Σ (n^2+1)^-s = 2^-s + Σ_k Binomial[-s,k] (Zeta[2s+2k] - 1),
 * or a published 20-digit value. */
static void test_slow_algebraic_machine(void) {
    ASSERT_CLOSE("NSum[1/n^1.05,{n,1,Infinity}]", "Zeta[21/20]", 21e-10);
    ASSERT_CLOSE("NSum[1/n^(11/10),{n,1,Infinity}]", "Zeta[11/10]", 11e-10);
    ASSERT_CLOSE("NSum[1/n^(51/50),{n,1,Infinity}]", "Zeta[51/50]", 51e-10);
    ASSERT_CLOSE("NSum[1/n^1.5,{n,1,Infinity}]", "Zeta[3/2]", 3e-10);
    ASSERT_CLOSE("NSum[1/(n+1/2)^(3/2),{n,0,Infinity}]", "Zeta[3/2,1/2]", 5e-10);
    ASSERT_CLOSE("NSum[1/n^(21/20)+1/n^(41/20),{n,1,Infinity}]",
                 "Zeta[21/20]+Zeta[41/20]", 22e-10);
    ASSERT_CLOSE("NSum[1/(n^2+1)^(11/20),{n,1,Infinity}]",
                 "2^(-11/20)+Sum[Binomial[-11/20,k](Zeta[11/10+2k]-1),{k,0,80}]", 1e-9);
    ASSERT_CLOSE("NSum[1/(n^2+1)^(3/4),{n,1,Infinity}]",
                 "2^(-3/4)+Sum[Binomial[-3/4,k](Zeta[3/2+2k]-1),{k,0,80}]", 3e-10);
    /* ISSAC '97 Systems Challenge, problem 3. */
    ASSERT_CLOSE("NSum[1/(n^Pi+n^2+n^Sqrt[2]+1)^(1/3),{n,1,Infinity}]",
                 "21.19324037771154094", 21e-10);
}

/* A slowly varying factor defeats any power-law tail: Σ_{n≥2} 1/(n Log[n]^2)
 * = 2.10974280123689197448 has a remainder 1/Log[X] beyond the double range.
 * NSum must then either be right to 10 digits or say that it is not (an
 * NSum::accgl message, observed through Check). */
static void test_slow_log_honest(void) {
    char* s = eval_str("Quiet[Check[If[Abs[NSum[1/(n Log[n]^2),{n,2,Infinity}]"
                       " - 2.10974280123689197448] < 3*10^-10, ok, bad], warned]]");
    ASSERT_MSG(strcmp(s, "ok") == 0 || strcmp(s, "warned") == 0,
               "1/(n Log[n]^2): silently inaccurate: %s", s);
    free(s);
}

/* Linear-phase summands, split by frequency (Abel–Plana per oscillatory
 * group).  Closed forms: sum sin(a n)/n = (Pi - a)/2 for 0 < a < 2 Pi (a = 7
 * aliases to 7 - 2 Pi), sum cos(n)/n = -Log[2 Sin[1/2]], the Clausen-type
 * sum cos(n t)/n^2 = Pi^2/6 - Pi t/2 + t^2/4, sum e^(i n)/n^2 = Li2(e^i),
 * sum 2^-n cos(n)/n = -Re Log[1 - e^i/2]. */
static void test_trig_oscillatory(void) {
    ASSERT_CLOSE("NSum[Sin[n]/n,{n,1,Infinity}]", "(Pi-1)/2", 1e-12);
    ASSERT_CLOSE("NSum[Sin[2 n]/n,{n,1,Infinity}]", "(Pi-2)/2", 1e-12);
    ASSERT_CLOSE("NSum[Sin[7 n]/n,{n,1,Infinity}]", "(3 Pi-7)/2", 1e-12);
    ASSERT_CLOSE("NSum[Cos[n]/n,{n,1,Infinity}]", "-Log[2 Sin[1/2]]", 1e-12);
    ASSERT_CLOSE("NSum[Cos[n]/n^2,{n,1,Infinity}]", "Pi^2/6 - Pi/2 + 1/4", 1e-12);
    ASSERT_CLOSE("NSum[Cos[3 n]/n^2,{n,1,Infinity}]", "Pi^2/6 - 3 Pi/2 + 9/4", 1e-12);
    ASSERT_CLOSE("NSum[Cos[Pi n/3]/n,{n,1,Infinity}]", "0", 1e-12);
    /* a non-oscillatory part (1/(2 n^2)) next to the e^(+-2 i n) groups */
    ASSERT_CLOSE("NSum[Sin[n]^2/n^2,{n,1,Infinity}]", "(Pi-1)/2", 1e-12);
    ASSERT_CLOSE("NSum[E^(I n)/n^2,{n,1,Infinity}]", "PolyLog[2, E^I]", 1e-12);
    ASSERT_CLOSE("NSum[2^-n Cos[n]/n,{n,1,Infinity}]", "-Re[Log[1 - E^I/2]]", 1e-12);
    /* 1/Sqrt[n] amplitude: Re Li_(1/2)(e^i), independently evaluated through
     * the Hurwitz-zeta reflection formula */
    ASSERT_CLOSE("NSum[Cos[n]/Sqrt[n],{n,1,Infinity}]",
                 "-0.1941089350921820497391492449", 1e-12);
    /* arbitrary precision */
    ASSERT_CLOSE("NSum[Sin[n]/n,{n,1,Infinity},WorkingPrecision->30]", "(Pi-1)/2", 1e-26);
    ASSERT_CLOSE("NSum[Cos[n]/Sqrt[n],{n,1,Infinity},WorkingPrecision->30]",
                 "-0.1941089350921820497391492449`30", 1e-25);
    ASSERT_CLOSE("NSum[Cos[n]/n^2,{n,1,Infinity},WorkingPrecision->30]",
                 "Pi^2/6 - Pi/2 + 1/4", 1e-26);
}

/* Dirichlet series of arithmetic functions: the summand is not a smooth
 * function of n, so no smooth tail model applies (the extrapolators returned
 * 2.6421 for zeta(2)^2 = 2.70581 and claimed 1e-7).  Closed forms are the
 * Dirichlet generating functions: sum d(n)/n^s = zeta(s)^2, sum sigma_a(n)/n^s
 * = zeta(s) zeta(s - a), sum phi(n)/n^s = zeta(s - 1)/zeta(s), sum mu(n)/n^s =
 * 1/zeta(s), sum lambda(n)/n^s = zeta(2 s)/zeta(s); a lower bound above 1
 * removes the head terms d(1)/1 + d(2)/2^2 = 3/2, phi(1) + phi(2)/8 = 9/8. */
static void test_dirichlet_series(void) {
    ASSERT_CLOSE("NSum[DivisorSigma[0, n]/n^2, {n, 1, Infinity}]", "Pi^4/36", 1e-13);
    ASSERT_CLOSE("NSum[DivisorSigma[0, n]/n^3, {n, 1, Infinity}]",
                 "1.44494079843363423391368507881", 1e-13);      /* zeta(3)^2 */
    ASSERT_CLOSE("NSum[DivisorSigma[1, n]/n^3, {n, 1, Infinity}]",
                 "Zeta[3] Pi^2/6", 1e-13);
    ASSERT_CLOSE("NSum[EulerPhi[n]/n^3, {n, 1, Infinity}]",
                 "1.368432777620205875736766", 1e-13);           /* zeta(2)/zeta(3) */
    ASSERT_CLOSE("NSum[MoebiusMu[n]/n^2, {n, 1, Infinity}]", "6/Pi^2", 1e-13);
    ASSERT_CLOSE("NSum[LiouvilleLambda[n]/n^2, {n, 1, Infinity}]", "(Pi^4/90)/(Pi^2/6)", 1e-13);
    ASSERT_CLOSE("NSum[DivisorSigma[0, n]/n^2, {n, 3, Infinity}]", "Pi^4/36 - 3/2", 1e-13);
    ASSERT_CLOSE("NSum[EulerPhi[n]/n^3, {n, 3, Infinity}]",
                 "1.368432777620205875736766 - 9/8", 1e-13);
    /* arbitrary precision: the closed form at the working precision */
    ASSERT_CLOSE("NSum[DivisorSigma[0, n]/n^2, {n, 1, Infinity}, WorkingPrecision -> 30]",
                 "Pi^4/36", 1e-27);
    ASSERT_CLOSE("NSum[EulerPhi[n]/n^3, {n, 1, Infinity}, WorkingPrecision -> 30]",
                 "1.36843277762020587573676585398`30", 1e-27);
    /* outside the half-plane of convergence the zeta ratio is not the sum:
     * Sum declines (EulerPhi[n]/n^(3/2) diverges, Zeta[1/2]/Zeta[3/2] < 0) */
    char* s = eval_str("Head[Sum[EulerPhi[n]/n^(3/2), {n, 1, Infinity}]]");
    ASSERT_MSG(strcmp(s, "Sum") == 0, "divergent phi series closed: %s", s);
    free(s);
    s = eval_str("Head[Sum[DivisorSigma[2, n]/n^3, {n, 1, Infinity}]]");
    ASSERT_MSG(strcmp(s, "Sum") == 0, "divergent sigma_2 series closed: %s", s);
    free(s);
    /* von Mangoldt: sum Lambda(n)/n^s = -zeta'(s)/zeta(s) (mpmath values); the
     * terms log p/p^s on prime powers only, the tail ~ N^(1-s), so no direct
     * or extrapolated sum gets past a few digits at s = 3/2 */
    ASSERT_CLOSE("NSum[MangoldtLambda[n]/n^2, {n, 1, Infinity}]",
                 "0.5699609930945328063998643600197300024035", 1e-13);
    ASSERT_CLOSE("NSum[MangoldtLambda[n]/n^(3/2), {n, 1, Infinity}, WorkingPrecision -> 35]",
                 "1.50523535578826791942204360303927497104`38", 1e-32);
    ASSERT_CLOSE("NSum[MangoldtLambda[n]/n^2, {n, 3, Infinity}]",
                 "0.5699609930945328063998643600197300024035 - Log[2]/4", 1e-13);
    s = eval_str("Head[Sum[MangoldtLambda[n]/n, {n, 1, Infinity}]]");
    ASSERT_MSG(strcmp(s, "Sum") == 0, "divergent Lambda(n)/n series closed: %s", s);
    free(s);
    /* the conditionally convergent boundary s = 1 of mu is kept: sum mu(n)/n = 0 */
    s = eval_str("Sum[MoebiusMu[n]/n, {n, 1, Infinity}]");
    ASSERT_MSG(strcmp(s, "0") == 0, "sum mu(n)/n: %s", s);
    free(s);
}

/* An integer-only summand with no closed form: the partial-sum extrapolation
 * must either be right or say it is not.  sum Omega(n)/n^3 = zeta(3) P_omega,
 * P_omega = sum_p 1/(p^3 - 1) = sum_k P(3k) = 0.194118169832633792 (prime zeta
 * values, computed independently), i.e. 0.2333410860760, and the
 * prime zeta P(2) = 0.452247420041065498506543. */
static void test_nonsmooth_honest(void) {
    char* s = eval_str("Quiet[Check[If[Abs[NSum[PrimeOmega[n]/n^3,{n,1,Infinity}]"
                       " - 0.233341086076023955302] < 10^-10, ok, bad], warned]]");
    ASSERT_MSG(strcmp(s, "ok") == 0 || strcmp(s, "warned") == 0,
               "PrimeOmega[n]/n^3: silently inaccurate: %s", s);
    free(s);
    s = eval_str("Quiet[Check[If[Abs[NSum[1/Prime[n]^2,{n,1,Infinity}]"
                 " - 0.452247420041065498506543] < 10^-10, ok, bad], warned]]");
    ASSERT_MSG(strcmp(s, "ok") == 0 || strcmp(s, "warned") == 0,
               "1/Prime[n]^2: silently inaccurate: %s", s);
    free(s);
}

/* Step-function and digit-length summands (nsum_blocks.c).  A summand defined
 * through Round[n] has values off the integers but is a staircase there, so
 * Euler–Maclaurin saw the jumps: 1/Round[n]^2 came back 1.62533 (Pi^2/6 =
 * 1.64493) and Sum 1/(n IntegerLength[n]^2) as 3.32572 (4.32597), with no
 * warning.  Now: the step is detected, smooth integer-only heads go to Levin u
 * from the best start, and heads smooth on blocks [b^(k-1), b^k) are summed
 * block by block (Euler–Maclaurin in log n per block, the block series by the
 * same routine).  Reference values: grouped by digit length, mpmath dps 90
 * (benchmarks/237, truth.tsv). */
static void test_step_and_block_summands(void) {
    char* s = eval_str("h1[n_?NumericQ] := 1/Round[n]^2; SetAttributes[h1, NHoldAll];"
                       " d1[n_?NumericQ] := 1/(Round[n] IntegerLength[Round[n]]^2);"
                       " SetAttributes[d1, NHoldAll];"
                       " d2[n_?NumericQ] := 1/(Round[n] BitLength[Round[n]]^2);"
                       " SetAttributes[d2, NHoldAll];"
                       " d9[n_?NumericQ] := (-1)^BitLength[Round[n]]/(Round[n] BitLength[Round[n]]);"
                       " SetAttributes[d9, NHoldAll];"
                       " g1[n_Integer] := 1/(n IntegerLength[n]^2); ok");
    free(s);
    /* smooth on the integers: Levin u, to machine precision, no warning */
    s = eval_str("Check[If[Abs[NSum[h1[n], {n, 1, Infinity}] - Pi^2/6] < 10^-13, ok, bad], warned]");
    ASSERT_MSG(strcmp(s, "ok") == 0, "1/Round[n]^2: %s", s);
    free(s);
    ASSERT_CLOSE("NSum[h1[n], {n, 1, Infinity}, WorkingPrecision -> 30]", "Pi^2/6", 1e-27);
    /* blocked: base 10, base 2, alternating by block, integer-only head */
    s = eval_str("Check[If[Abs[NSum[d1[n], {n, 1, Infinity}]"
                 " - 4.32597084672048487201066701504] < 10^-13, ok, bad], warned]");
    ASSERT_MSG(strcmp(s, "ok") == 0, "1/(n IntegerLength[n]^2): %s", s);
    free(s);
    ASSERT_CLOSE("NSum[d2[n], {n, 1, Infinity}]", "1.49245452678426373572652175665", 1e-12);
    ASSERT_CLOSE("NSum[d9[n], {n, 1, Infinity}]", "-0.73353427573011417741968544586", 1e-13);
    ASSERT_CLOSE("NSum[g1[n], {n, 1, Infinity}]", "4.32597084672048487201066701504", 1e-13);
    ASSERT_CLOSE("NSum[d1[n], {n, 1, Infinity}, WorkingPrecision -> 60]",
                 "4.3259708467204848720106670150410465454229605594715427602057521256", 1e-48);
    /* a sum from a later start takes the partial first block */
    ASSERT_CLOSE("NSum[d1[n], {n, 7, Infinity}] + Sum[1/n, {n, 1, 6}]",
                 "4.32597084672048487201066701504", 1e-12);
    free(eval_str("ClearAll[h1, d1, d2, d9, g1]"));
}

/* A self-similar summand nested past any reachable index (Putnam 2002 A6:
 * f(n) = n f(BitLength[n])): the block series is followed as far as the
 * index budget allows and the unbounded remainder is reported, not hidden —
 * the value carries NSum::accgl with an error that brackets the true sum. */
static void test_block_unbounded_honest(void) {
    char* s = eval_str("f177[n_Integer] := If[n <= 2, n, n f177[BitLength[n]]];"
                       " Quiet[Check[NSum[1/f177[n], {n, 1, Infinity}]; ok, warned]]");
    ASSERT_MSG(strcmp(s, "warned") == 0, "nested f(BitLength): no accgl: %s", s);
    free(s);
    free(eval_str("ClearAll[f177]"));
}

/* Bessel summands: J_nu(c n), Y_nu(c n) oscillate like e^(+-i c n), a phase the
 * symbolic split finds only through the Hankel decomposition.  Closed forms
 * from Poisson summation: sum_{n>=1} J_0(n x) = -1/2 + (1/x)(1 + 2 sum_{m >= 1,
 * 2 Pi m < x} (1 - (2 Pi m/x)^2)^(-1/2)), so 1/2 at x = 1, 0 at x = 2 and one
 * aliased frequency at x = 10; sum J_1(n x)/n = 1 - x/4 (0 < x < 2 Pi); and
 * sum Y_0(n) = 0.61698948158458515325 (an independent high-precision sum). */
static void test_bessel_oscillatory(void) {
    ASSERT_CLOSE("NSum[BesselJ[0, n], {n, 1, Infinity}]", "1/2", 1e-12);
    ASSERT_CLOSE("NSum[BesselJ[0, 2 n], {n, 1, Infinity}]", "0", 1e-12);
    ASSERT_CLOSE("NSum[BesselJ[0, 10 n], {n, 1, Infinity}]",
                 "-1/2 + (1/10) (1 + 2/Sqrt[1 - (2 Pi/10)^2])", 1e-12);
    ASSERT_CLOSE("NSum[BesselJ[1, n]/n, {n, 1, Infinity}]", "3/4", 1e-12);
    ASSERT_CLOSE("NSum[BesselJ[1, 3 n]/n, {n, 1, Infinity}]", "1/4", 1e-12);
    ASSERT_CLOSE("NSum[BesselY[0, n], {n, 1, Infinity}]", "0.61698948158458515325", 1e-12);
    ASSERT_CLOSE("NSum[BesselJ[0, n], {n, 1, Infinity}, WorkingPrecision -> 30]",
                 "1/2", 1e-26);
}

/* Alternating summands the numeric sign probe misses (a zero first term,
 * a magnitude that rises before it falls).  -eta'(s) closed forms:
 * sum (-1)^n Log[n]/n = EulerGamma Log[2] - Log[2]^2/2, and
 * sum (-1)^(n+1) Log[n]/n^2 = -(Log[2] Pi^2/12 + Zeta'(2)/2) with
 * Zeta'(2) = -0.937548254315843753702574094568. */
static void test_alternating_log(void) {
    ASSERT_CLOSE("NSum[(-1)^n Log[n]/n,{n,1,Infinity}]",
                 "EulerGamma Log[2] - Log[2]^2/2", 1e-12);
    ASSERT_CLOSE("NSum[Cos[Pi n] Log[n]/n,{n,1,Infinity}]",
                 "EulerGamma Log[2] - Log[2]^2/2", 1e-12);
    ASSERT_CLOSE("NSum[(-1)^(n+1) Log[n]/n^2,{n,1,Infinity}]",
                 "-(Log[2] Pi^2/12 - 0.93754825431584375370257409456786/2)", 1e-12);
    ASSERT_CLOSE("NSum[(-1)^n Log[n]/n,{n,1,Infinity},WorkingPrecision->30]",
                 "EulerGamma Log[2] - Log[2]^2/2", 1e-26);
    /* the eta values themselves stay exact through the split */
    ASSERT_CLOSE("NSum[(-1)^(n+1)/Sqrt[n],{n,1,Infinity},WorkingPrecision->30]",
                 "(1 - Sqrt[2]) Zeta[1/2]", 1e-26);
}

/* Logarithmically slow tails: the Euler–Maclaurin tail integral in
 * log-substituted variables.  Expected values: Euler–Maclaurin at base 10^4
 * with the tail integral in closed form (int_N^oo dx/(x Log[x]^k) =
 * 1/((k-1) Log[N]^(k-1))), the first terms summed directly — an independent
 * route (the h3 values are the ISSAC holdout's, from N = 10^4 and 2·10^4). */
static void test_slow_log_tails(void) {
    ASSERT_CLOSE("NSum[1/(n Log[n]^2),{n,2,Infinity}]", "2.109742801236891974479257198", 1e-12);
    ASSERT_CLOSE("NSum[1/(n Log[n] Log[Log[n]]^2),{n,3,Infinity}]",
                 "38.40676809282178631849374770", 1e-10);
    ASSERT_CLOSE("NSum[1/(n Log[n]^3),{n,2,Infinity}]",
                 "N[Total[Table[1/(n Log[n]^3), {n, 2, 9999}]] + 1/(2 Log[10^4]^2)"
                 " + (1/(n Log[n]^3) /. n -> 10^4)/2"
                 " - (1/12) (D[1/(n Log[n]^3), n] /. n -> 10^4)"
                 " + (1/720) (D[1/(n Log[n]^3), {n, 3}] /. n -> 10^4), 20]", 1e-12);
    ASSERT_CLOSE("NSum[1/(n Log[n]^2),{n,2,Infinity},WorkingPrecision->30]",
                 "2.109742801236891974479257198`30", 1e-24);
    ASSERT_CLOSE("NSum[1/(n Log[n] Log[Log[n]]^2),{n,3,Infinity},WorkingPrecision->30]",
                 "38.40676809282178631849374770`30", 1e-22);
    /* oscillatory, through the phase split's arbitrary-precision half-line
     * quadrature (nodes of each level reused by the next): (Pi - 1)/2 */
    ASSERT_CLOSE("NSum[Sin[n]/n,{n,1,Infinity},WorkingPrecision->30]",
                 "N[(Pi - 1)/2, 40]", 1e-26);
    ASSERT_CLOSE("NSum[Cos[n]/n^2,{n,1,Infinity},WorkingPrecision->30]",
                 "N[Pi^2/6 - Pi/2 + 1/4, 40]", 1e-26);
    /* 1/Log[n] - 1/Log[n+1] telescopes to 1/Log[2], but its continuous
     * extension cancels to nothing at the huge x the tail needs: no method can
     * certify it, so it must warn rather than be silently wrong. */
    char* s = eval_str("Quiet[Check[If[Abs[NSum[1/Log[n] - 1/Log[n+1],{n,2,Infinity}]"
                       " - 1/Log[2]] < 10^-10, ok, bad], warned]]");
    ASSERT_MSG(strcmp(s, "ok") == 0 || strcmp(s, "warned") == 0,
               "1/Log[n] - 1/Log[n+1]: silently inaccurate: %s", s);
    free(s);
}

static void test_unevaluated_forms(void) {
    char* s1 = eval_str("NSum[1/i^2,{i,1,n}]");         /* symbolic upper bound */
    ASSERT_MSG(strstr(s1, "NSum[") != NULL, "symbolic bound should stay: %s", s1);
    free(s1);
    char* s2 = eval_str("NSum[1/i^2]");                 /* too few args */
    ASSERT_MSG(strstr(s2, "NSum[") != NULL, "1-arg NSum should stay: %s", s2);
    free(s2);
    char* s3 = eval_str("NSum[1/i^2,i]");               /* not an iterator spec */
    ASSERT_MSG(strstr(s3, "NSum[") != NULL, "non-list spec should stay: %s", s3);
    free(s3);
}

static void test_attributes(void) {
    char* s1 = eval_str("MemberQ[Attributes[NSum], Protected]");
    ASSERT_MSG(strcmp(s1, "True") == 0, "NSum should be Protected: %s", s1);
    free(s1);
    char* s2 = eval_str("MemberQ[Attributes[NSum], HoldAll]");
    ASSERT_MSG(strcmp(s2, "True") == 0, "NSum should be HoldAll: %s", s2);
    free(s2);
}

static void test_memory_loop(void) {
    const char* inputs[] = {
        "NSum[1/i^2,{i,1,5}]",
        "NSum[1/i^2,{i,1,Infinity}]",
        "NSum[(-1)^k/(2k+1),{k,0,Infinity}]",
        "NSum[1/n^(11/10),{n,1,Infinity}]",
        "NSum[Log[x]/x^(2+2I),{x,1,Infinity}]",
        "NSum[2^i,{i,0,Infinity}]",
        "NSum[1/2^i,{i,0,Infinity,2}]",
        "NSum[(-1)^n (2/n)^k/k^2,{n,2,Infinity},{k,1,n}]",
        "NSum[1/i^2,{i,1,Infinity},WorkingPrecision->30]",
        "NSum[1./Prime[i]^2,{i,1,Infinity}]",   /* integer-only: continuity probe */
        "NSum[Sin[n]/n,{n,1,Infinity}]",        /* phase split, Abel–Plana */
        "NSum[(-1)^n Log[n]/n,{n,1,Infinity}]", /* phase split, CVZ */
        "NSum[1/(n Log[n]^2),{n,2,Infinity}]",  /* log-substituted tail */
        "NSum[1/(n IntegerLength[n]^2),{n,1,Infinity}]",  /* digit blocks */
    };
    /* 3 passes suffice for leak detection; the complex Log[x]/x^(2+2I) summand
     * is the slowest path (oscillatory Euler–Maclaurin), so keep the count low. */
    for (int i = 0; i < 3; i++)
        for (size_t k = 0; k < sizeof(inputs)/sizeof(inputs[0]); k++) {
            Expr* p = parse_expression(inputs[k]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
}

/* N of an unevaluated definite Sum / Product / Integrate falls back to NSum /
 * NProduct / NIntegrate (it used to come back as 2.0 Sum[...]).  Expected:
 * mpmath 2 nsum(besselj(5/8 + 2k + 1, 7/4), [0, inf]) = 0.892429323529416225329596,
 * Prod (1 + 1/k^3) = Cosh[Sqrt[3] Pi/2]/Pi.  A free parameter keeps it symbolic. */
static void test_n_of_unevaluated_sum(void) {
    ASSERT_CLOSE("N[2 Sum[BesselJ[5/8 + 2 k + 1, 7/4], {k, 0, Infinity}]]",
                 "0.892429323529416225329596", 1e-12);
    ASSERT_CLOSE("N[2 Sum[BesselJ[5/8 + 2 k + 1, 7/4], {k, 0, Infinity}], 20]",
                 "0.892429323529416225329596", 1e-18);
    ASSERT_CLOSE("N[Product[1 + 1/k^3, {k, 1, Infinity}]]", "Cosh[Sqrt[3] Pi/2]/Pi", 1e-9);
    ASSERT_CLOSE("N[Integrate[BesselJ[0, x] Cos[x/2], {x, 0, Infinity}]]", "2/Sqrt[3]", 1e-9);
    char* s = eval_str("Head[N[Sum[f[k], {k, 1, Infinity}]]]");
    ASSERT_MSG(strcmp(s, "Sum") == 0, "symbolic summand must stay a Sum, got %s", s);
    free(s);
}

/* Single-element iterators {c}, by the Sum convention (v0.6511): a real c is
 * Floor[c] terms, each evaluated once ({Pi} is 3); anything else (a free
 * symbol: a numeric sum has no indefinite form) is NSum::itraw. */
static void test_count_iterator(void) {
    ASSERT_CLOSE("NSum[1/2, {4}]", "2", 1e-12);
    ASSERT_CLOSE("NSum[1, {Pi}]", "3", 1e-12);
    ASSERT_CLOSE("NSum[1, {2.5}]", "2", 1e-12);
    ASSERT_CLOSE("NSum[1, {0}]", "0", 1e-12);
    /* the same accuracy as the named form {j, 1, 2} gives */
    ASSERT_CLOSE("NSum[1/k^2, {k, 1, Infinity}, {2}]", "Pi^2/3", 1e-6);
    char* s = eval_str("Block[{c = 0}, {NSum[c++, {4}], c}]");
    ASSERT_MSG(strcmp(s, "{6.0, 4}") == 0, "side effects: %s", s);
    free(s);
    s = eval_str("Check[NSum[1, {n}], \"msg\"]");
    ASSERT_MSG(strcmp(s, "\"msg\"") == 0, "NSum[1, {n}] -> %s", s);
    free(s);
    s = eval_str("Quiet[NSum[1, {\"s\"}]]");
    ASSERT_MSG(strncmp(s, "NSum[", 5) == 0, "NSum[1, {\"s\"}] -> %s", s);
    free(s);
}

/* Infinite sums whose range starts far out (the tail of NProduct's log-sum
 * from 10^6 is one): the Euler–Maclaurin base, tail integral and alternating
 * and geometric routes must hold at a large start index.  Closed forms:
 * PolyGamma[1, n] = Sum_{k>=n} 1/k^2, Zeta[s, n] = Sum_{k>=n} 1/k^s,
 * Sum_{k>=n} (-1)^k/k = (-1)^n (PolyGamma[0, (n+1)/2] - PolyGamma[0, n/2])/2,
 * Sum_{k>=1} r^k = r/(1 - r).  No call may emit a message. */
static void test_far_start_tails(void) {
    ASSERT_CLOSE("Check[NSum[1/k^2, {k, 10^6 + 1, Infinity}], \"msg\"]",
                 "N[PolyGamma[1, 10^6 + 1], 30]", 1e-21);
    ASSERT_CLOSE("Check[NSum[1/k^3, {k, 10^4, Infinity}], \"msg\"]",
                 "N[Zeta[3, 10^4], 30]", 1e-23);
    ASSERT_CLOSE("Check[NSum[1/k^1.5, {k, 10^5, Infinity}], \"msg\"]",
                 "N[Zeta[3/2, 10^5], 30]", 1e-16);
    ASSERT_CLOSE("Check[NSum[(-1)^k/k, {k, 10^6, Infinity}], \"msg\"]",
                 "N[(PolyGamma[0, (10^6 + 1)/2] - PolyGamma[0, 10^6/2])/2, 30]", 1e-21);
    ASSERT_CLOSE("Check[NSum[Exp[-k/10^5], {k, 1, Infinity}], \"msg\"]",
                 "N[1/(Exp[1/10^5] - 1), 30]", 1e-9);
    ASSERT_CLOSE("Check[NSum[1/k^2, {k, 10^6 + 1, Infinity}, WorkingPrecision -> 30], \"msg\"]",
                 "N[PolyGamma[1, 10^6 + 1], 40]", 1e-34);
    ASSERT_CLOSE("Check[NSum[1/k^3, {k, 10^4, Infinity}, WorkingPrecision -> 30], \"msg\"]",
                 "N[Zeta[3, 10^4], 40]", 1e-36);
    ASSERT_CLOSE("Check[NSum[Exp[-k/10^5], {k, 1, Infinity}, WorkingPrecision -> 30], \"msg\"]",
                 "N[1/(Exp[1/10^5] - 1), 40]", 1e-24);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_far_start_tails);

    TEST(test_direct_finite);
    TEST(test_count_iterator);
    TEST(test_infinite_monotone);
    TEST(test_integer_only_summand);
    TEST(test_alternating);
    TEST(test_methods_explicit);
    TEST(test_levin);
    TEST(test_step_di);
    TEST(test_complex_summand);
    TEST(test_large_finite);
    TEST(test_long_finite);
    TEST(test_long_finite_oscillating);
    TEST(test_long_finite_steps);
    TEST(test_long_finite_interpreted);
    TEST(test_long_finite_working_precision);
    TEST(test_long_finite_timing);
    TEST(test_long_finite_random);
    TEST(test_multidim);
    TEST(test_fast_decay);
    TEST(test_saturated_sequence);
    TEST(test_verify_convergence);
    TEST(test_peaked_summand);
    TEST(test_composite_precision);
    TEST(test_working_precision);
    TEST(test_slow_algebraic_machine);
    TEST(test_slow_log_honest);
    TEST(test_trig_oscillatory);
    TEST(test_dirichlet_series);
    TEST(test_nonsmooth_honest);
    TEST(test_step_and_block_summands);
    TEST(test_block_unbounded_honest);
    TEST(test_bessel_oscillatory);
    TEST(test_alternating_log);
    TEST(test_slow_log_tails);
    TEST(test_unevaluated_forms);
    TEST(test_attributes);
    TEST(test_memory_loop);
    TEST(test_n_of_unevaluated_sum);

    printf("All nsum_tests passed.\n");
    return 0;
}
