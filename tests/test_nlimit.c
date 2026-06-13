/* Tests for NLimit — numerical limit.
 *
 * Cover: the default EulerSum (Richardson/Romberg) method on finite, infinite,
 * complex-valued and complex-point limits; the SequenceLimit (Wynn epsilon)
 * method and WynnDegree scaling; the Direction option (one-sided and arbitrary
 * complex rays); Scale and Terms; arbitrary precision via WorkingPrecision; the
 * noise / cannot-recognise diagnostic; option/argument-shape edge cases; the
 * Protected attribute; stress batches; and memory hygiene.
 *
 * Expected values are taken from closed forms where available, otherwise from
 * Mathematica's NumericalCalculus`NLimit documentation. Numerical results are
 * compared *inside* the language (N[Abs[result - expected]] < tol) rather than
 * by parsing the printed form, which is rounded to ~6 significant figures.
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

/* True if input is within `tol` of expected, evaluated at full internal
 * precision. Uses the L1 magnitude |Re d| + |Im d| of the difference so a
 * spurious tiny imaginary residual does not leave the modulus complex. */
static bool close_to(const char* input, const char* expected, double tol) {
    char buf[2048];
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

/* ------------------------------------------------------------------------
 *  EulerSum (default) — elementary finite limits
 * ---------------------------------------------------------------------- */

static void test_finite_real_limits(void) {
    ASSERT_CLOSE("NLimit[Sin[x]/x,x->0]", "1", 1e-6);
    ASSERT_CLOSE("NLimit[Tan[x]/x,x->0]", "1", 1e-6);
    ASSERT_CLOSE("NLimit[Log[1+x]/x,x->0]", "1", 1e-5);
    ASSERT_CLOSE("NLimit[(Cos[x]-1)/x^2,x->0]", "-1/2", 1e-9);
    ASSERT_CLOSE("NLimit[(1-Cos[x])/x^2,x->0]", "1/2", 1e-9);
    ASSERT_CLOSE("NLimit[(E^x-1-x)/x^2,x->0]", "1/2", 1e-9);
    ASSERT_CLOSE("NLimit[(Sqrt[1+x]-1)/x,x->0]", "1/2", 1e-6);
}

/* ------------------------------------------------------------------------
 *  Infinite limit points (approach outward on a ray)
 * ---------------------------------------------------------------------- */

static void test_infinite_limits(void) {
    ASSERT_CLOSE("NLimit[(1+1/n)^n,n->Infinity]", "E", 1e-6);
    ASSERT_CLOSE("NLimit[n Sin[1/n],n->Infinity]", "1", 1e-8);
    ASSERT_CLOSE("NLimit[Tanh[x],x->Infinity]", "1", 1e-4);
    /* manifestly complex along a real ray */
    ASSERT_CLOSE("NLimit[(1+I/x)^x,x->Infinity]", "Exp[I]", 1e-5);
}

/* ------------------------------------------------------------------------
 *  Complex limit point and complex-valued result
 * ---------------------------------------------------------------------- */

static void test_complex_point(void) {
    /* doc: -8.17617e-7 - 1.5708 I; the real part is a roundoff residual. */
    ASSERT_CLOSE("NLimit[Tanh[Pi x]/(1+x^2),x->I]", "-Pi/2 I", 1e-4);
}

static void test_complex_result(void) {
    ASSERT_CLOSE("NLimit[(Exp[I x]-1)/x,x->0]", "I", 1e-8);
}

/* ------------------------------------------------------------------------
 *  SequenceLimit (Wynn epsilon) and WynnDegree
 * ---------------------------------------------------------------------- */

static void test_sequencelimit(void) {
    ASSERT_CLOSE("NLimit[(1+1/n)^n,n->Infinity,Method->SequenceLimit]", "E", 3e-3);
    ASSERT_CLOSE("NLimit[(10^x-1)/x,x->0,Method->SequenceLimit]", "Log[10]", 1e-2);
    ASSERT_CLOSE("NLimit[(10^x-1)/x,x->0,Terms->10,Method->SequenceLimit]",
                 "Log[10]", 1e-3);
    /* a poorly-behaved default-method case where SequenceLimit is the right
     * tool: Tan[z] approaching I Infinity. doc: 0.+1. I */
    ASSERT_CLOSE("NLimit[Tan[z],z->Infinity I,Method->SequenceLimit]", "I", 1e-4);
}

static void test_wynndegree_improves(void) {
    /* Higher WynnDegree must reduce the error on a smooth geometric tail. */
    ASSERT_CLOSE("NLimit[(10^x-1)/x,x->0,Terms->10,Method->SequenceLimit,"
                 "WynnDegree->4]", "Log[10]", 1e-6);
    char* s = eval_str(
        "N[Abs[NLimit[(10^x-1)/x,x->0,Terms->10,Method->SequenceLimit,WynnDegree->4]-Log[10]]"
        " < Abs[NLimit[(10^x-1)/x,x->0,Terms->10,Method->SequenceLimit,WynnDegree->1]-Log[10]]]");
    ASSERT_MSG(strcmp(s, "True") == 0, "WynnDegree 4 should beat 1: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Direction — one-sided and arbitrary complex rays
 * ---------------------------------------------------------------------- */

static void test_direction(void) {
    /* z + Conjugate[z]/z at 0: the value depends on the direction of approach. */
    ASSERT_CLOSE("NLimit[z+Conjugate[z]/z,z->0,Direction->1]", "1", 1e-9);
    ASSERT_CLOSE("NLimit[z+Conjugate[z]/z,z->0,Direction->-I]", "-1", 1e-9);
    ASSERT_CLOSE("NLimit[z+Conjugate[z]/z,z->0,Direction->-Exp[225 Degree I]]",
                 "-I", 1e-6);
    /* default (Automatic == -1) approaches from larger values. */
    ASSERT_CLOSE("NLimit[z+Conjugate[z]/z,z->0]", "1", 1e-9);
}

/* ------------------------------------------------------------------------
 *  Scale and Terms options
 * ---------------------------------------------------------------------- */

static void test_scale_terms(void) {
    ASSERT_CLOSE("NLimit[Sin[x]/x,x->0,Scale->1/4]", "1", 1e-6);
    ASSERT_CLOSE("NLimit[Sin[x]/x,x->0,Terms->12]", "1", 1e-9);
    /* Tanh approaches its limit exponentially => few terms suffice. */
    ASSERT_CLOSE("NLimit[Tanh[x],x->Infinity,Terms->5]", "1", 1e-2);
}

/* ------------------------------------------------------------------------
 *  Arbitrary precision (MPFR)
 * ---------------------------------------------------------------------- */

static void test_working_precision(void) {
    ASSERT_CLOSE("NLimit[Sin[x]/x,x->0,WorkingPrecision->35,Terms->20]", "1", 1e-25);
    ASSERT_CLOSE("NLimit[(2^x-1)/x,x->0,WorkingPrecision->30,Terms->14]",
                 "Log[2]", 1e-18);
    ASSERT_CLOSE("NLimit[(10^x-1)/x,x->0,WorkingPrecision->30,Terms->14,"
                 "Method->SequenceLimit,WynnDegree->4]", "Log[10]", 1e-12);
    /* complex result at high precision */
    ASSERT_CLOSE("NLimit[(Exp[I x]-1)/x,x->0,WorkingPrecision->30,Terms->14]",
                 "I", 1e-18);
}

/* ------------------------------------------------------------------------
 *  Noise / cannot-recognise diagnostic
 * ---------------------------------------------------------------------- */

static void test_noise(void) {
    /* power-law approach to infinity: NLimit cannot recognise a value. */
    char* s = eval_str("NLimit[1/x,x->0]");
    ASSERT_MSG(strstr(s, "NLimit[") != NULL, "1/x should stay unevaluated: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Edge cases / unevaluated forms
 * ---------------------------------------------------------------------- */

static void test_unevaluated_forms(void) {
    char* s1 = eval_str("NLimit[Sin[x]/x,x->a]");          /* non-numeric point */
    ASSERT_MSG(strstr(s1, "NLimit[") != NULL, "symbolic point should stay: %s", s1);
    free(s1);
    char* s2 = eval_str("NLimit[Sin[x]/x]");               /* too few args */
    ASSERT_MSG(strstr(s2, "NLimit[") != NULL, "1-arg NLimit should stay: %s", s2);
    free(s2);
    char* s3 = eval_str("NLimit[Sin[x]/x,{x,0}]");         /* list, not a rule */
    ASSERT_MSG(strstr(s3, "NLimit[") != NULL, "list spec should stay: %s", s3);
    free(s3);
    char* s4 = eval_str("NLimit[Sin[x]/x,3->0]");          /* non-symbol var */
    ASSERT_MSG(strstr(s4, "NLimit[") != NULL, "non-symbol var should stay: %s", s4);
    free(s4);
}

static void test_protected(void) {
    char* s = eval_str("MemberQ[Attributes[NLimit], Protected]");
    ASSERT_MSG(strcmp(s, "True") == 0, "NLimit should be Protected: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Stress: a broad battery of limits, each checked for correctness
 * ---------------------------------------------------------------------- */

static void test_stress_battery(void) {
    struct { const char* in; const char* exp; double tol; } cases[] = {
        { "NLimit[(Exp[x]-1)/x,x->0]",            "1",         1e-7 },
        { "NLimit[(Cos[x]-1)/x^2,x->0]",          "-1/2",      1e-9 },
        { "NLimit[Sin[3x]/Sin[5x],x->0]",         "3/5",       1e-3 },
        { "NLimit[(Tan[x]-Sin[x])/x^3,x->0]",     "1/2",       1e-4 },
        { "NLimit[x (Sqrt[1+1/x]-1),x->Infinity]","1/2",       1e-5 },
        { "NLimit[(1+2/n)^n,n->Infinity]",        "Exp[2]",    2e-4 },
        { "NLimit[x ArcTan[1/x],x->Infinity]",    "1",         1e-6 },
        { "NLimit[n (E^(1/n)-1),n->Infinity]",    "1",         1e-7 },
        { "NLimit[(Sinh[x])/x,x->0]",             "1",         1e-7 },
        { "NLimit[(1-1/n)^n,n->Infinity]",        "1/E",       1e-6 },
        { "NLimit[ArcTan[x],x->-Infinity]",       "-Pi/2",     1e-4 },
    };
    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++)
        ASSERT_CLOSE(cases[i].in, cases[i].exp, cases[i].tol);
}

/* Stress the extrapolators directly across a sweep of Terms / WynnDegree. */
static void test_stress_terms_sweep(void) {
    const char* tmpl_euler =
        "NLimit[(Cos[x]-1)/x^2,x->0,Terms->%d]";
    char buf[256];
    for (int t = 3; t <= 16; t++) {
        snprintf(buf, sizeof buf, tmpl_euler, t);
        ASSERT_CLOSE(buf, "-1/2", 1e-3);
    }
    for (int d = 1; d <= 4; d++) {
        snprintf(buf, sizeof buf,
                 "NLimit[(1+1/n)^n,n->Infinity,Terms->%d,Method->SequenceLimit,WynnDegree->%d]",
                 2 * (d + 1) + 2, d);
        ASSERT_CLOSE(buf, "E", 1e-2);
    }
}

static void test_memory_loop(void) {
    const char* inputs[] = {
        "NLimit[Sin[x]/x,x->0]",
        "NLimit[(1+1/n)^n,n->Infinity]",
        "NLimit[(10^x-1)/x,x->0,Method->SequenceLimit,Terms->10]",
        "NLimit[Tanh[Pi x]/(1+x^2),x->I]",
        "NLimit[(Exp[I x]-1)/x,x->0,WorkingPrecision->30]",
        "NLimit[1/x,x->0]",                    /* noise path */
        "NLimit[z+Conjugate[z]/z,z->0,Direction->-I]",
    };
    for (int i = 0; i < 20; i++) {
        for (size_t k = 0; k < sizeof(inputs)/sizeof(inputs[0]); k++) {
            Expr* p = parse_expression(inputs[k]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
    }
}

/* ------------------------------------------------------------------------
 *  Main
 * ---------------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_finite_real_limits);
    TEST(test_infinite_limits);
    TEST(test_complex_point);
    TEST(test_complex_result);

    TEST(test_sequencelimit);
    TEST(test_wynndegree_improves);

    TEST(test_direction);
    TEST(test_scale_terms);

    TEST(test_working_precision);

    TEST(test_noise);
    TEST(test_unevaluated_forms);
    TEST(test_protected);

    TEST(test_stress_battery);
    TEST(test_stress_terms_sweep);
    TEST(test_memory_loop);

    printf("All nlimit_tests passed.\n");
    return 0;
}
