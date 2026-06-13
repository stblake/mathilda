/* Tests for NResidue — numerical residue via contour integration.
 *
 * Cover: simple poles (at the origin and shifted, real and complex z0),
 * higher-order / rational poles, essential singularities (Exp[1/x],
 * Sin[1/x], Cos[1/x] where the residue is 0), radius-based pole isolation,
 * manual list-threading over arg 1, the Radius -> Automatic adaptive search,
 * arbitrary precision (WorkingPrecision), the Derivative-via-residue Zeta
 * identity, option/argument-shape edge cases, the Protected attribute, and
 * memory hygiene.
 *
 * Floating-point results are checked with a numeric tolerance (not exact
 * string match) via the eval_real() helper, which wraps the expression in
 * Re[...] and parses the printed value.
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

/* Evaluate Re[input] to a machine double. Robust to Real / Complex / MPFR
 * results; returns 0.0 if the result is non-numeric. */
static double eval_real(const char* input) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Re[%s]", input);
    char* s = eval_str(buf);
    double v = strtod(s, NULL);
    free(s);
    return v;
}

/* Evaluate Im[input] to a machine double. */
static double eval_imag(const char* input) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Im[%s]", input);
    char* s = eval_str(buf);
    double v = strtod(s, NULL);
    free(s);
    return v;
}

#define ASSERT_NEAR(input, expected, tol) do {                               \
    double _v = eval_real(input);                                            \
    ASSERT_MSG(fabs(_v - (double)(expected)) < (double)(tol),                \
               "%s\n  Re = %.17g, expected %.17g (tol %g)",                  \
               (input), _v, (double)(expected), (double)(tol));              \
} while (0)

/* ------------------------------------------------------------------------
 *  Simple poles
 * ---------------------------------------------------------------------- */

static void test_simple_pole_origin(void) {
    /* Res(1/z, 0) = 1. */
    ASSERT_NEAR("NResidue[1/x,{x,0}]", 1.0, 1e-9);
    /* Imaginary part is a spurious residual; must be tiny. */
    ASSERT_MSG(fabs(eval_imag("NResidue[1/x,{x,0}]")) < 1e-9,
               "imaginary part of Res(1/x,0) should be ~0");
}

static void test_simple_pole_shifted(void) {
    /* Res(1/(z-2), 2) = 1. */
    ASSERT_NEAR("NResidue[1/(x-2),{x,2}]", 1.0, 1e-9);
    /* Res(3/(z+5), -5) = 3. */
    ASSERT_NEAR("NResidue[3/(x+5),{x,-5}]", 3.0, 1e-9);
}

static void test_rational_two_poles(void) {
    /* 1/(z^2 - 2.7 z + 1.7) has poles at z=1 and z=1.7.  Res at z=1 is
     * 1/(1 - 1.7) = -10/7 ~ -1.4285714. The default radius 1/100 isolates
     * the pole at 1 (the other is 0.7 away). z0 = 1. (machine real). */
    ASSERT_NEAR("NResidue[1/(1.7-2.7z+z^2),{z,1.}]", -10.0 / 7.0, 1e-5);
}

static void test_complex_z0(void) {
    /* Res(1/(z - I), I) = 1, with a genuinely complex centre. */
    ASSERT_NEAR("NResidue[1/(x-I),{x,I}]", 1.0, 1e-9);
}

/* ------------------------------------------------------------------------
 *  Essential singularities (where symbolic Residue fails)
 * ---------------------------------------------------------------------- */

static void test_essential_sin(void) {
    /* Res(Sin[1/(10 z)], 0) = 1/10. */
    ASSERT_NEAR("NResidue[Sin[1/(10x)],{x,0}]", 0.1, 1e-6);
}

static void test_essential_exp_radius1(void) {
    /* Res(Exp[1/z], 0) = 1, needs a larger radius to converge. */
    ASSERT_NEAR("NResidue[Exp[1/x],{x,0},Radius->1]", 1.0, 1e-6);
}

static void test_essential_cos_is_zero(void) {
    /* Cos[1/z] has only even negative powers => residue 0. */
    ASSERT_NEAR("NResidue[Cos[1/x],{x,0},Radius->1]", 0.0, 1e-6);
}

/* ------------------------------------------------------------------------
 *  Radius-based pole isolation
 * ---------------------------------------------------------------------- */

static void test_radius_isolation(void) {
    /* 1/z + 1/(z + 0.005): two simple poles, each residue 1.  The default
     * radius 1/100 encloses BOTH (sum = 2); a radius of 0.001 isolates the
     * pole at the origin (residue 1). */
    ASSERT_NEAR("NResidue[1/x+1/(x+0.005),{x,0}]", 2.0, 1e-6);
    ASSERT_NEAR("NResidue[1/x+1/(x+0.005),{x,0},Radius->0.001]", 1.0, 1e-6);
}

/* ------------------------------------------------------------------------
 *  Adaptive radius (Radius -> Automatic)
 * ---------------------------------------------------------------------- */

static void test_radius_automatic(void) {
    /* Automatic finds a converging radius without manual tuning, even for
     * the essential singularity that fails at the default radius. */
    ASSERT_NEAR("NResidue[Exp[1/x],{x,0},Radius->Automatic]", 1.0, 1e-6);
    ASSERT_NEAR("NResidue[1/x,{x,0},Radius->Automatic]", 1.0, 1e-9);
}

/* ------------------------------------------------------------------------
 *  Manual list-threading over arg 1
 * ---------------------------------------------------------------------- */

static void test_list_threading(void) {
    /* NResidue threads element-wise over arg 1, keeping {x, 0} fixed. The
     * {z, z0} spec must NOT be split (the whole reason NResidue is not
     * ATTR_LISTABLE). */
    char* s = eval_str("Head[NResidue[{Exp[1/x],Sin[1/x],Cos[1/x]},{x,0},Radius->1]]");
    ASSERT_STR_EQ(s, "List");
    free(s);
    s = eval_str("Length[NResidue[{Exp[1/x],Sin[1/x],Cos[1/x]},{x,0},Radius->1]]");
    ASSERT_STR_EQ(s, "3");
    free(s);
    ASSERT_NEAR("Part[NResidue[{Exp[1/x],Sin[1/x],Cos[1/x]},{x,0},Radius->1],1]", 1.0, 1e-6);
    ASSERT_NEAR("Part[NResidue[{Exp[1/x],Sin[1/x],Cos[1/x]},{x,0},Radius->1],2]", 1.0, 1e-6);
    ASSERT_NEAR("Part[NResidue[{Exp[1/x],Sin[1/x],Cos[1/x]},{x,0},Radius->1],3]", 0.0, 1e-6);
}

/* ------------------------------------------------------------------------
 *  Arbitrary precision (WorkingPrecision)
 * ---------------------------------------------------------------------- */

static void test_mpfr_simple_pole(void) {
    /* The MPFR path must reach the simple-pole residue 1 at 30 digits. */
    ASSERT_NEAR("NResidue[1/(x-2),{x,2},WorkingPrecision->30]", 1.0, 1e-12);
}

static void test_mpfr_zeta_derivative(void) {
    /* Derivative[10][Zeta][0] = 10! * a_{-1} of Zeta[x]/x^11.  The known
     * value is approximately -3628799.99945676588...  (requires complex
     * MPFR Zeta on the contour). */
    ASSERT_NEAR("10! NResidue[Zeta[x]/x^11,{x,0},Radius->1/2,WorkingPrecision->30]",
                -3628799.99945676588, 1e-2);
}

/* ------------------------------------------------------------------------
 *  Argument-shape / option edge cases (return unevaluated)
 * ---------------------------------------------------------------------- */

static void test_unevaluated_forms(void) {
    /* One argument: unevaluated. */
    char* s = eval_str("NResidue[1/x]");
    ASSERT_MSG(strstr(s, "NResidue") != NULL, "expected unevaluated, got: %s", s);
    free(s);
    /* Second argument not a {z, z0} list: unevaluated. */
    s = eval_str("NResidue[1/x, x]");
    ASSERT_MSG(strstr(s, "NResidue") != NULL, "expected unevaluated, got: %s", s);
    free(s);
    /* Variable slot not a symbol: unevaluated. */
    s = eval_str("NResidue[1/x, {2, 0}]");
    ASSERT_MSG(strstr(s, "NResidue") != NULL, "expected unevaluated, got: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Protected attribute
 * ---------------------------------------------------------------------- */

static void test_protected(void) {
    char* s = eval_str("Attributes[NResidue]");
    ASSERT_MSG(strstr(s, "Protected") != NULL,
               "expected Protected in attributes, got: %s", s);
    /* Must NOT be Listable — generic threading would split {z, z0}. */
    ASSERT_MSG(strstr(s, "Listable") == NULL,
               "NResidue must not be Listable, got: %s", s);
    free(s);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene — valgrind --leak-check=full should be clean
 * ---------------------------------------------------------------------- */

static void test_memory_loop(void) {
    const char* cases[] = {
        "NResidue[1/x,{x,0}]",
        "NResidue[1/(x-2),{x,2}]",
        "NResidue[Sin[1/(10x)],{x,0}]",
        "NResidue[Exp[1/x],{x,0},Radius->1]",
        "NResidue[1/x+1/(x+0.005),{x,0},Radius->0.001]",
        "NResidue[{Exp[1/x],Sin[1/x],Cos[1/x]},{x,0},Radius->1]",
        "NResidue[1/x,{x,0},Radius->Automatic]",
        "NResidue[1/(x-2),{x,2},WorkingPrecision->30]",
        "NResidue[Sqrt[x],{x,0},Radius->1]",   /* branch-cut path */
        "NResidue[1/x]",                        /* unevaluated path */
        NULL
    };
    for (int rep = 0; rep < 5; rep++) {
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

    TEST(test_simple_pole_origin);
    TEST(test_simple_pole_shifted);
    TEST(test_rational_two_poles);
    TEST(test_complex_z0);

    TEST(test_essential_sin);
    TEST(test_essential_exp_radius1);
    TEST(test_essential_cos_is_zero);

    TEST(test_radius_isolation);
    TEST(test_radius_automatic);

    TEST(test_list_threading);

    TEST(test_mpfr_simple_pole);
    TEST(test_mpfr_zeta_derivative);

    TEST(test_unevaluated_forms);
    TEST(test_protected);
    TEST(test_memory_loop);

    printf("All nresidue_tests passed.\n");
    return 0;
}
