/* Tests for NSeries — numerical Taylor/Laurent series via circle sampling + DFT.
 *
 * Cover: Taylor expansion of Exp about real and complex centres, arbitrary
 * precision (WorkingPrecision), Laurent expansion about an essential
 * singularity (Sin[x + 1/x]), Radius-controlled annulus selection of a rational
 * function's Laurent series, the SeriesData result shape, option/argument-shape
 * edge cases, the Protected (non-Listable) attribute, and memory hygiene.
 *
 * Coefficients are checked numerically inside the language (not by string
 * match): the k-th coefficient of (x - x0)^e is Part[<series>, 3, e + n + 1].
 * eval_real / eval_imag wrap that in Re[...] / Im[...] and parse the value.
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

/* Evaluate `input` and return its printed form (caller frees). */
static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

/* |Re[coef_e] - expected| and |Im[coef_e] - expected|, where `expected` is a
 * Mathilda expression. The deviation is computed INSIDE the language so the
 * printer's 6-significant-figure rounding of machine reals does not limit the
 * comparison — only the (tiny) magnitude of the difference is parsed. */
static double coef_err(const char* part, const char* ns, int n, int e,
                       const char* expected) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Abs[%s[Part[%s,3,%d]] - (%s)]",
             part, ns, e + n + 1, expected);
    char* s = eval_str(buf);
    double v = strtod(s, NULL);
    free(s);
    return v;
}
static double coef_re_err(const char* ns, int n, int e, const char* expected) {
    return coef_err("Re", ns, n, e, expected);
}
static double coef_im_err(const char* ns, int n, int e, const char* expected) {
    return coef_err("Im", ns, n, e, expected);
}

#define ASSERT_COEF_RE(ns, n, e, expected, tol, what)                        \
    ASSERT_MSG(coef_re_err((ns), (n), (e), (expected)) < (double)(tol),      \
               "%s: Re deviation %.3g exceeds tol %g", (what),               \
               coef_re_err((ns), (n), (e), (expected)), (double)(tol))
#define ASSERT_COEF_IM(ns, n, e, expected, tol, what)                        \
    ASSERT_MSG(coef_im_err((ns), (n), (e), (expected)) < (double)(tol),      \
               "%s: Im deviation %.3g exceeds tol %g", (what),               \
               coef_im_err((ns), (n), (e), (expected)), (double)(tol))

/* ------------------------------------------------------------------------
 *  Taylor expansion of Exp about 0
 * ---------------------------------------------------------------------- */

static void test_exp_taylor_machine(void) {
    const char* ns = "NSeries[Exp[x],{x,0,5}]";
    /* a_k = 1/k! for k >= 0; negative-index coefficients are spurious ~0. */
    ASSERT_COEF_RE(ns, 5, 0, "1",     1e-10, "Exp a0");
    ASSERT_COEF_RE(ns, 5, 1, "1",     1e-10, "Exp a1");
    ASSERT_COEF_RE(ns, 5, 2, "1/2",   1e-10, "Exp a2");
    ASSERT_COEF_RE(ns, 5, 3, "1/6",   1e-10, "Exp a3");
    ASSERT_COEF_RE(ns, 5, 4, "1/24",  1e-10, "Exp a4");
    ASSERT_COEF_RE(ns, 5, 5, "1/120", 1e-10, "Exp a5");
    /* Spurious parts must be tiny: negative-index coeff and imaginary part. */
    ASSERT_COEF_RE(ns, 5, -1, "0", 1e-9, "Exp a_-1 spurious");
    ASSERT_COEF_IM(ns, 5, 2,  "0", 1e-9, "Exp a2 imag spurious");
}

static void test_exp_taylor_mpfr(void) {
    /* Extended precision drives the coefficients (and residuals) much tighter. */
    const char* ns = "NSeries[Exp[x],{x,0,5},WorkingPrecision->30]";
    ASSERT_COEF_RE(ns, 5, 3, "1/6",   1e-20, "Exp@30 a3");
    ASSERT_COEF_RE(ns, 5, 5, "1/120", 1e-20, "Exp@30 a5");
    ASSERT_COEF_IM(ns, 5, 3, "0",     1e-20, "Exp@30 a3 imag");
    ASSERT_COEF_RE(ns, 5, -1, "0",    1e-20, "Exp@30 a_-1");
}

static void test_exp_complex_center(void) {
    /* About x0 = I: a_k = Exp[I]/k! = (cos 1 + i sin 1)/k!. */
    const char* ns = "NSeries[Exp[x],{x,I,5}]";
    /* a_k = Exp[I]/k! = (cos 1 + i sin 1)/k!; compare against in-language refs. */
    ASSERT_COEF_RE(ns, 5, 0, "Cos[1]",   1e-9, "Exp@I a0 re");
    ASSERT_COEF_IM(ns, 5, 0, "Sin[1]",   1e-9, "Exp@I a0 im");
    ASSERT_COEF_RE(ns, 5, 2, "Cos[1]/2", 1e-9, "Exp@I a2 re");
    ASSERT_COEF_IM(ns, 5, 2, "Sin[1]/2", 1e-9, "Exp@I a2 im");
    /* The SeriesData base must be (x - I): check it prints with x0 = I. */
    char* str = eval_str("Part[NSeries[Exp[x],{x,I,5}],2]");
    ASSERT_STR_EQ(str, "I");
    free(str);
}

/* ------------------------------------------------------------------------
 *  Laurent expansion about an essential singularity
 * ---------------------------------------------------------------------- */

static void test_laurent_essential(void) {
    /* Sin[x + 1/x] about 0: odd-power coefficients nonzero, even ~ 0.
     * From the reference: a_1 = a_-1 ~ 0.576725, a_3 = a_-3 ~ -0.128943. */
    const char* ns = "NSeries[Sin[x+1/x],{x,0,10}]";
    ASSERT_COEF_RE(ns, 10,  1, "0.576725", 1e-5, "Sin a1");
    ASSERT_COEF_RE(ns, 10, -1, "0.576725", 1e-5, "Sin a_-1");
    ASSERT_COEF_RE(ns, 10,  3, "-0.128943", 1e-5, "Sin a3");
    ASSERT_COEF_RE(ns, 10, -3, "-0.128943", 1e-5, "Sin a_-3");
    /* Even powers are zero by symmetry. */
    ASSERT_COEF_RE(ns, 10, 2, "0", 1e-5, "Sin a2 (even)");
    ASSERT_COEF_RE(ns, 10, 0, "0", 1e-5, "Sin a0 (even)");
}

/* ------------------------------------------------------------------------
 *  Radius-controlled annulus selection of a rational Laurent series
 * ---------------------------------------------------------------------- */

static void test_radius_far_annulus(void) {
    /* 1/((1+x)(3+x)) with Radius->5 selects the |x|>=3 annulus, where the
     * Laurent series has the integer coefficients ... 121, -40, 13, -4, 1
     * for x^-6 .. x^-2. */
    const char* ns = "NSeries[1/((1+x)(3+x)),{x,0,10},Radius->5]";
    ASSERT_COEF_RE(ns, 10, -2,  "1",    1e-4, "rat a_-2");
    ASSERT_COEF_RE(ns, 10, -3,  "-4",   1e-4, "rat a_-3");
    ASSERT_COEF_RE(ns, 10, -4,  "13",   1e-4, "rat a_-4");
    ASSERT_COEF_RE(ns, 10, -10, "9841", 1e-1, "rat a_-10");
    /* Non-negative powers vanish in this annulus. */
    ASSERT_COEF_RE(ns, 10, 0, "0", 1e-4, "rat a0 far annulus");
}

static void test_radius_near_annulus(void) {
    /* Radius->2 selects 1<|x|<3: a0 = -1/6, a_-1 = 1/2. */
    const char* ns = "NSeries[1/((1+x)(3+x)),{x,0,10},Radius->2]";
    ASSERT_COEF_RE(ns, 10,  0, "-1/6", 1e-5, "rat near a0");
    ASSERT_COEF_RE(ns, 10, -1, "1/2",  1e-5, "rat near a_-1");
}

/* ------------------------------------------------------------------------
 *  SeriesData result shape
 * ---------------------------------------------------------------------- */

static void test_result_shape(void) {
    char* s;
    s = eval_str("Head[NSeries[Exp[x],{x,0,5}]]");
    ASSERT_STR_EQ(s, "SeriesData"); free(s);
    s = eval_str("Length[NSeries[Exp[x],{x,0,5}]]");
    ASSERT_STR_EQ(s, "6"); free(s);
    /* coefficient list length = 2n+1. */
    s = eval_str("Length[Part[NSeries[Exp[x],{x,0,5}],3]]");
    ASSERT_STR_EQ(s, "11"); free(s);
    /* nmin = -n, nmax = n+1, den = 1. */
    s = eval_str("Part[NSeries[Exp[x],{x,0,5}],4]");
    ASSERT_STR_EQ(s, "-5"); free(s);
    s = eval_str("Part[NSeries[Exp[x],{x,0,5}],5]");
    ASSERT_STR_EQ(s, "6"); free(s);
    s = eval_str("Part[NSeries[Exp[x],{x,0,5}],6]");
    ASSERT_STR_EQ(s, "1"); free(s);
}

/* ------------------------------------------------------------------------
 *  Argument-shape / option edge cases (return unevaluated)
 * ---------------------------------------------------------------------- */

static void test_unevaluated_forms(void) {
    char* s;
    /* Negative order: unevaluated. */
    s = eval_str("NSeries[Exp[x],{x,0,-1}]");
    ASSERT_MSG(strstr(s, "NSeries") != NULL, "expected unevaluated, got: %s", s);
    free(s);
    /* Spec without an order: unevaluated. */
    s = eval_str("NSeries[Exp[x],{x,0}]");
    ASSERT_MSG(strstr(s, "NSeries") != NULL, "expected unevaluated, got: %s", s);
    free(s);
    /* Second argument not a list: unevaluated. */
    s = eval_str("NSeries[Exp[x],x]");
    ASSERT_MSG(strstr(s, "NSeries") != NULL, "expected unevaluated, got: %s", s);
    free(s);
    /* Variable slot not a symbol: unevaluated. */
    s = eval_str("NSeries[Exp[x],{2,0,5}]");
    ASSERT_MSG(strstr(s, "NSeries") != NULL, "expected unevaluated, got: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Protected attribute (and NOT Listable)
 * ---------------------------------------------------------------------- */

static void test_protected(void) {
    char* s = eval_str("Attributes[NSeries]");
    ASSERT_MSG(strstr(s, "Protected") != NULL,
               "expected Protected in attributes, got: %s", s);
    ASSERT_MSG(strstr(s, "Listable") == NULL,
               "NSeries must not be Listable, got: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene — valgrind --leak-check=full should be clean
 * ---------------------------------------------------------------------- */

static void test_memory_loop(void) {
    const char* cases[] = {
        "NSeries[Exp[x],{x,0,5}]",
        "NSeries[Exp[x],{x,0,5},WorkingPrecision->30]",
        "NSeries[Exp[x],{x,I,5}]",
        "NSeries[Sin[x+1/x],{x,0,10}]",
        "NSeries[1/((1+x)(3+x)),{x,0,10},Radius->5]",
        "NSeries[1/((1+x)(3+x)),{x,0,10},Radius->2]",
        "NSeries[Exp[x],{x,0,-1}]",   /* unevaluated path */
        "NSeries[Exp[x],{x,0}]",      /* unevaluated path */
        NULL
    };
    for (int rep = 0; rep < 3; rep++) {
        for (int i = 0; cases[i]; i++) {
            Expr* p = parse_expression(cases[i]);
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

    TEST(test_exp_taylor_machine);
    TEST(test_exp_taylor_mpfr);
    TEST(test_exp_complex_center);

    TEST(test_laurent_essential);

    TEST(test_radius_far_annulus);
    TEST(test_radius_near_annulus);

    TEST(test_result_shape);

    TEST(test_unevaluated_forms);
    TEST(test_protected);
    TEST(test_memory_loop);

    printf("All nseries_tests passed.\n");
    return 0;
}
