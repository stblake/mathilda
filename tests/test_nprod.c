/* Tests for NProduct — numerical product.
 *
 * Cover: infinite products via Euler–Maclaurin (default) and explicit Wynn
 * epsilon; finite products (positive factors -> real, no spurious imaginary);
 * step di and the equivalent reindexing; multidimensional products with
 * independent and dependent inner bounds; complex infinite products;
 * AccuracyGoal/PrecisionGoal; NProductFactors; the Sin product identity;
 * arbitrary precision via WorkingPrecision; options and unevaluated forms; the
 * HoldAll/Protected attributes; and memory hygiene.
 *
 * Per Keiper 1992 ("The N functions of Mathematica") NProduct is evaluated as
 * Exp[NSum[Log[f], ...]].  Expected values come from closed forms or Keiper.
 * Numerical results are compared *inside* the language
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

static void test_infinite_euler_maclaurin(void) {
    /* Infinite, Euler-Maclaurin (default). Keiper: Product[1+1/i^2] = Sinh[Pi]/Pi. */
    ASSERT_CLOSE("NProduct[1+1/i^2,{i,1,Infinity}]", "Sinh[Pi]/Pi", 1e-7);
    ASSERT_CLOSE("NProduct[1+1/n^2,{n,1,Infinity}]", "Sinh[Pi]/Pi", 1e-7);
}

static void test_finite(void) {
    /* finite products: positive factors -> real, no spurious imaginary */
    ASSERT_CLOSE("NProduct[1+(-1)^i/i^2,{i,100,10000}]", "1.00005", 1e-4);
    ASSERT_CLOSE("NProduct[i,{i,1,5}]", "120", 1e-9);          /* exact 5! */
    ASSERT_CLOSE("NProduct[1+1/i,{i,1,9}]", "10", 1e-9);       /* telescopes to 10 */
}

static void test_step_di(void) {
    ASSERT_CLOSE("NProduct[1+1/2^i,{i,0,Infinity,2}]", "2.71182", 1e-4);
    ASSERT_CLOSE("NProduct[1+1/2^(2j),{j,0,Infinity}]", "2.71182", 1e-4);
}

static void test_method_wynn(void) {
    /* Wynn's epsilon excels on geometric tails (it is documented-weak on
     * monotone ones): the explicit method routes through and agrees with the
     * default to machine precision on q = 1/2. */
    ASSERT_CLOSE("NProduct[1+1/2^n,{n,1,Infinity},Method->\"WynnEpsilon\"]",
                 "NProduct[1+1/2^n,{n,1,Infinity}]", 1e-9);
}

static void test_multidim(void) {
    /* Multidimensional recursion (the outer body is an inner NProduct over the
     * remaining specs), exact on finite ranges -- including a dependent inner
     * bound {j,1,i} that sees the outer index.  NProduct[i+j,{i,1,3},{j,1,2}]
     * = Prod_i (i+1)(i+2) = 6*12*20; with {j,1,i}: 2 * (3*4) * (4*5*6).
     * (Infinite-*outer* multidimensional products are only approximate -- see
     * the changelog note -- because the outer NSum must extrapolate over a
     * numeric black-box inner product; not asserted here.) */
    ASSERT_CLOSE("NProduct[i+j,{i,1,3},{j,1,2}]", "1440", 1e-9);
    ASSERT_CLOSE("NProduct[i+j,{i,1,3},{j,1,i}]", "2880", 1e-9);
}

static void test_complex(void) {
    /* Keiper: 1.43706 + 1.07945 I. */
    ASSERT_CLOSE("NProduct[1+E^(I n 2/3)/n^2,{n,1,Infinity}]",
                 "1.43706 + 1.07945 I", 1e-4);
}

static void test_goals(void) {
    /* AccuracyGoal/PrecisionGoal -> 8 target ~8 digits; on the machine path the
     * result tracks the closed form to that accuracy. */
    ASSERT_CLOSE("NProduct[1+1/n^2,{n,1,Infinity},AccuracyGoal->8,PrecisionGoal->8]",
                 "Sinh[Pi]/Pi", 1e-7);
}

static void test_factors(void) {
    /* NProductFactors maps to NSum's NSumTerms: taking more leading factors
     * explicitly leaves a convergent result unchanged. */
    ASSERT_CLOSE("NProduct[1+1/n^2,{n,1,Infinity},NProductFactors->30]",
                 "Sinh[Pi]/Pi", 1e-7);
}

static void test_sin_identity(void) {
    ASSERT_CLOSE("Block[{x=1/2}, x NProduct[1-x^2/(Pi^2 k^2),{k,1,Infinity}]]",
                 "Sin[0.5]", 1e-5);
}

static void test_working_precision(void) {
    /* MPFR path: a finite product is exact at arbitrary working precision. */
    ASSERT_CLOSE("NProduct[i,{i,1,10},WorkingPrecision->30]", "3628800", 1e-15);
    /* Infinite product at 30 digits.  With NSum's WP-scaled Euler-Maclaurin tail
     * integral, numerical contour corrections and guard digits, the log-sum now
     * reaches the working precision, which Exp carries into the product. */
    ASSERT_CLOSE("NProduct[1+1/n^2,{n,1,Infinity},WorkingPrecision->30]",
                 "Sinh[Pi]/Pi", 1e-25);
}

static void test_unevaluated_forms(void) {
    char* s1 = eval_str("NProduct[1+1/i^2,{i,1,m}]");      /* symbolic upper bound */
    ASSERT_MSG(strstr(s1, "NProduct[") != NULL, "symbolic bound should stay: %s", s1);
    free(s1);
    char* s2 = eval_str("NProduct[1+1/i^2]");              /* too few args */
    ASSERT_MSG(strstr(s2, "NProduct[") != NULL, "1-arg NProduct should stay: %s", s2);
    free(s2);
    char* s3 = eval_str("NProduct[1+1/i^2,i]");            /* not an iterator spec */
    ASSERT_MSG(strstr(s3, "NProduct[") != NULL, "non-list spec should stay: %s", s3);
    free(s3);
    char* s4 = eval_str("NProduct[1+1/i^2,{i,1,Infinity},Bogus->1]"); /* unknown option */
    ASSERT_MSG(strstr(s4, "NProduct[") != NULL, "unknown option should stay: %s", s4);
    free(s4);
}

static void test_attributes(void) {
    char* s1 = eval_str("MemberQ[Attributes[NProduct], Protected]");
    ASSERT_MSG(strcmp(s1, "True") == 0, "NProduct should be Protected: %s", s1);
    free(s1);
    char* s2 = eval_str("MemberQ[Attributes[NProduct], HoldAll]");
    ASSERT_MSG(strcmp(s2, "True") == 0, "NProduct should be HoldAll: %s", s2);
    free(s2);
}

static void test_memory_loop(void) {
    /* A couple of passes over varied code paths suffices for leak detection.
     * The complex infinite product is intentionally omitted here -- it is the
     * slowest path (NSum's complex Euler-Maclaurin) and is already exercised by
     * test_complex; repeating it would push the suite past the harness alarm. */
    const char* inputs[] = {
        "NProduct[i,{i,1,5}]",
        "NProduct[1+1/i^2,{i,1,Infinity}]",
        "NProduct[1+1/2^i,{i,0,Infinity,2}]",
        "NProduct[i+j,{i,1,3},{j,1,i}]",
        "NProduct[1+1/n^2,{n,1,Infinity},WorkingPrecision->25]",
        "NProduct[1+1/i^2,{i,1,m}]",
    };
    for (int i = 0; i < 2; i++)
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

    TEST(test_infinite_euler_maclaurin);
    TEST(test_finite);
    TEST(test_step_di);
    TEST(test_method_wynn);
    TEST(test_multidim);
    TEST(test_complex);
    TEST(test_goals);
    TEST(test_factors);
    TEST(test_sin_identity);
    TEST(test_working_precision);
    TEST(test_unevaluated_forms);
    TEST(test_attributes);
    TEST(test_memory_loop);

    printf("All nprod_tests passed.\n");
    return 0;
}
