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
    };
    for (int i = 0; i < 8; i++)
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
    TEST(test_alternating);
    TEST(test_methods_explicit);
    TEST(test_step_di);
    TEST(test_complex_summand);
    TEST(test_large_finite);
    TEST(test_multidim);
    TEST(test_verify_convergence);
    TEST(test_working_precision);
    TEST(test_unevaluated_forms);
    TEST(test_attributes);
    TEST(test_memory_loop);

    printf("All nsum_tests passed.\n");
    return 0;
}
