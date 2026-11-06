/* Tests for NSum — numerical summation.
 *
 * Cover: Direct (small finite) sums; infinite sums via Automatic dispatch
 * (Euler–Maclaurin for monotone, Cohen–Villegas–Zagier for alternating, Wynn
 * epsilon fallback); the explicit Method settings; step di; complex summands;
 * large finite sums (difference of tails); multidimensional sums with
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

static void test_multidim(void) {
    /* Keiper multidimensional examples */
    ASSERT_CLOSE("NSum[(-1)^n (2/n)^k/k^2,{n,2,Infinity},{k,1,Infinity}]", "1.14434", 1e-4);
    ASSERT_CLOSE("NSum[(-1)^n (2/n)^k/k^2,{n,2,Infinity},{k,1,n}]", "0.770188", 1e-4);
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

int main(void) {
    symtab_init();
    core_init();

    TEST(test_direct_finite);
    TEST(test_infinite_monotone);
    TEST(test_integer_only_summand);
    TEST(test_alternating);
    TEST(test_methods_explicit);
    TEST(test_levin);
    TEST(test_step_di);
    TEST(test_complex_summand);
    TEST(test_large_finite);
    TEST(test_multidim);
    TEST(test_verify_convergence);
    TEST(test_peaked_summand);
    TEST(test_composite_precision);
    TEST(test_working_precision);
    TEST(test_slow_algebraic_machine);
    TEST(test_slow_log_honest);
    TEST(test_trig_oscillatory);
    TEST(test_alternating_log);
    TEST(test_slow_log_tails);
    TEST(test_unevaluated_forms);
    TEST(test_attributes);
    TEST(test_memory_loop);

    printf("All nsum_tests passed.\n");
    return 0;
}
