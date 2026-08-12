/* Tests for the numeric-domain fixes to Mod, Divide, QuotientRemainder,
 * Factorial, Factorial2, FactorialPower, BarnesG and Hyperfactorial.
 *
 * Two bug classes are covered:
 *   (A) silently-wrong machine coercion -- Mod / Divide read a Rational (or
 *       MPFR / Complex) sibling of a machine real as 0, so Mod[7/3, 0.5] was
 *       0.0 and Divide[2.5, 1/3] was ComplexInfinity;
 *   (B) missing non-integer / complex / arbitrary-precision continuations --
 *       Factorial (complex), Factorial2, FactorialPower (non-integer k),
 *       BarnesG and Hyperfactorial had no numeric path off the integers.
 *
 * Reference values are cross-checked against the equivalent Gamma / recurrence
 * constructions (e.g. Factorial[z] == Gamma[z+1], the Barnes recurrence
 * G(w+1) = Gamma(w) G(w), and Hyperfactorial[z] == Gamma[z+1]^z / BarnesG[z+1]).
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

/* ---- numeric helpers (mirrors test_fresnelc.c) ---------------------- */

static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a Real result", input);
    double v = r->data.real;
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_real(input);
    ASSERT_MSG(fabs(v - expected) <= tol,
               "%s: expected %.15g, got %.15g", input, expected, v);
}

static void assert_complex_close(const char* input, double er, double ei, double tol) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_FUNCTION &&
               r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
               r->data.function.arg_count == 2,
               "%s: expected Complex[..], got something else", input);
    Expr* re = r->data.function.args[0];
    Expr* im = r->data.function.args[1];
    ASSERT(re->type == EXPR_REAL && im->type == EXPR_REAL);
    ASSERT_MSG(fabs(re->data.real - er) <= tol && fabs(im->data.real - ei) <= tol,
               "%s: expected %.12g %+.12g I, got %.12g %+.12g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* =====================================================================
 * Group 1 -- silently-wrong coercion (Mod / Divide / QuotientRemainder)
 * ===================================================================== */

void test_mod_rational_sibling(void) {
    /* Machine real + exact Rational: the Rational must not be read as 0. */
    assert_close("Mod[7/3, 0.5]", 1.0 / 3.0, 1e-12);   /* was 0.0 */
    assert_close("Mod[2.5, 1/3]", 1.0 / 6.0, 1e-12);   /* was unevaluated */
    assert_close("Mod[2.5, 1/2]", 0.0, 1e-12);
    assert_close("Mod[1/3, 2.5]", 1.0 / 3.0, 1e-12);
    /* Three-argument Mod with an offset and a Rational modulus. */
    assert_close("Mod[8.5, 1/3, 1]", 7.0 / 6.0, 1e-12);
}

void test_mod_mpfr_sibling(void) {
    /* MPFR real + exact Rational stays at full precision (was 0.0). */
    assert_eval_startswith("Mod[N[7/3, 30], 1/2]", "0.33333333333333333333");
}

void test_mod_regressions(void) {
    assert_eval_eq("Mod[7, 3]", "1", 0);
    assert_eval_eq("Mod[7/3, 1/2]", "1/3", 0);      /* both exact -> rational path */
    assert_close("Mod[2.5, 0.5]", 0.0, 1e-12);      /* both machine */
    assert_eval_eq("Head[Mod[2.0, 1 + I]]", "Mod", 0);  /* complex modulus stays symbolic */
}

void test_divide_rational_mpfr_complex(void) {
    assert_close("Divide[2.5, 1/3]", 7.5, 1e-12);            /* was ComplexInfinity */
    assert_close("Divide[1/3, 2.5]", 1.0 / 3.0 / 2.5, 1e-12); /* was 0.0 */
    assert_complex_close("Divide[2.0, 1 + I]", 1.0, -1.0, 1e-12); /* was ComplexInfinity */
    assert_eval_startswith("Divide[N[1/4, 30], 2.0]", "0.125");   /* was 0.0 */
}

void test_divide_regressions(void) {
    assert_eval_eq("Divide[6, 3]", "2", 0);
    assert_close("Divide[2.0, 4.0]", 0.5, 1e-12);
    assert_close("2.5/(1/3)", 7.5, 1e-12);   /* operator form routes through Times/Power */
}

void test_quotientremainder(void) {
    /* Inherits the Mod fix: the remainder half read the Rational as 0. */
    assert_eval_eq("QuotientRemainder[7.5, 1/2][[1]]", "15", 0);
    assert_close("QuotientRemainder[7.5, 1/2][[2]]", 0.0, 1e-12);   /* was unevaluated */
    assert_eval_eq("QuotientRemainder[1/3, 7.5][[1]]", "0", 0);
    assert_close("QuotientRemainder[1/3, 7.5][[2]]", 1.0 / 3.0, 1e-12); /* was 0.0 */
}

/* =====================================================================
 * Group 2 -- Factorial complex
 * ===================================================================== */

void test_factorial_complex(void) {
    /* Factorial[z] == Gamma[z+1]; machine + arbitrary precision + symbolic. */
    assert_complex_close("N[Factorial[1/2 + I/3]]",
                         0.84184564078675363, 0.014450865220296312, 1e-9);
    assert_eval_startswith("N[Factorial[1/2 + I/3], 25]", "0.84184564078675362776113");
    assert_eval_eq("Head[Factorial[1/2 + I/3]]", "Factorial", 0);  /* exact -> symbolic */
    /* Real / integer regressions unchanged. */
    assert_eval_eq("Factorial[5]", "120", 0);
    assert_close("N[Factorial[2.5]]", 3.3233509704478426, 1e-9);
}

/* =====================================================================
 * Group 3 -- Factorial2 / FactorialPower non-integer continuations
 * ===================================================================== */

void test_factorial2_numeric(void) {
    /* Reproduces integer double factorials, extends to real/complex. */
    assert_close("N[Factorial2[5.0]]", 15.0, 1e-9);
    assert_close("N[Factorial2[4.0]]", 8.0, 1e-9);
    assert_close("N[Factorial2[3.5]]", 4.8323193861368527, 1e-9);
    assert_eval_startswith("N[Factorial2[7/2], 25]", "4.8323193861368526656583");
    assert_complex_close("N[Factorial2[1.0 + I]]", 0.250651, 0.100474, 1e-4);
    /* exact / integer stays as-is. */
    assert_eval_eq("Head[Factorial2[7/2]]", "Factorial2", 0);
    assert_eval_eq("Factorial2[6]", "48", 0);
    assert_eval_eq("Factorial2[7]", "105", 0);
}

void test_factorialpower_numeric(void) {
    /* Non-integer k: Gamma[n+1]/Gamma[n-k+1]. */
    assert_close("N[FactorialPower[3.5, 1.5]]", 5.8158641982837245, 1e-9);
    assert_eval_startswith("N[FactorialPower[7/2, 3/2], 25]", "5.8158641982837244645721");
    assert_eval_eq("Head[FactorialPower[7/2, 3/2]]", "FactorialPower", 0);  /* exact -> symbolic */
    /* Integer-k regressions (exact product path). */
    assert_eval_eq("FactorialPower[7, 3]", "210", 0);
    assert_eval_eq("FactorialPower[7/2, 2]", "35/4", 0);
    assert_close("N[FactorialPower[3.5, 2]]", 8.75, 1e-9);
}

/* =====================================================================
 * Group 3 -- BarnesG / Hyperfactorial continuations
 * ===================================================================== */

void test_barnesg_numeric(void) {
    /* Numeric at integer arguments must reproduce the exact superfactorials. */
    assert_close("N[BarnesG[5.0]]", 12.0, 1e-6);
    assert_close("N[BarnesG[6.0]]", 288.0, 1e-5);
    assert_close("N[BarnesG[7.0]]", 34560.0, 1e-3);
    assert_close("N[BarnesG[8.0]]", 24883200.0, 1.0);
    /* half-integer, machine + arbitrary precision. */
    assert_close("N[BarnesG[3.5]]", 1.2596482574951921, 1e-9);
    assert_eval_startswith("N[BarnesG[13/2], 30]", "2548.7457695684989897359061");
    assert_complex_close("N[BarnesG[2.5 + 1.0 I]]", 0.743798, -0.0953168, 1e-4);
    /* exact / integer stays as-is. */
    assert_eval_eq("Head[BarnesG[7/2]]", "BarnesG", 0);
    assert_eval_eq("BarnesG[6]", "288", 0);
}

void test_barnesg_recurrence(void) {
    /* G(w+1) = Gamma(w) G(w): a strong internal consistency check that the
     * asymptotic value and the Gamma recurrence agree. */
    assert_close("N[BarnesG[4.5]/BarnesG[3.5]] - Gamma[3.5]", 0.0, 1e-9);
    assert_close("N[BarnesG[8.3]/BarnesG[7.3]] - Gamma[7.3]", 0.0, 1e-6);
}

void test_hyperfactorial_numeric(void) {
    /* Hyperfactorial[z] = Gamma[z+1]^z / BarnesG[z+1]. */
    assert_close("N[Hyperfactorial[4.0]]", 27648.0, 1e-4);
    assert_close("N[Hyperfactorial[5.0]]", 86400000.0, 1.0);
    assert_close("N[Hyperfactorial[3.5]]", 1282.1220994534575, 1e-6);
    assert_eval_startswith("N[Hyperfactorial[7/2], 30]", "1282.1220994534574594154227");
    assert_complex_close("N[Hyperfactorial[2.5 + 1.0 I]]", -5.53952, -0.477309, 1e-4);
    /* Recurrence H(z)/H(z-1) = z^z. */
    assert_close("N[Hyperfactorial[4.3]/Hyperfactorial[3.3]] - 4.3^4.3", 0.0, 1e-6);
    /* exact / integer stays as-is. */
    assert_eval_eq("Head[Hyperfactorial[7/2]]", "Hyperfactorial", 0);
    assert_eval_eq("Hyperfactorial[4]", "27648", 0);
}

void test_binomial_complex_numeric(void) {
    /* Class-B continuation for Binomial: an inexact operand alongside a
     * complex operand must evaluate through the Gamma quotient, with the
     * exact Gaussian sibling carried along by the Times/Plus numeric
     * contagion.  Both were previously left unevaluated (or, before the
     * binomial_to_double fix, returned nan.0). */
    /* Exact Gaussian n, machine-real m: Binomial[1 + I, 5] = -1/12 - I/12. */
    assert_complex_close("Binomial[1 + I, 5.]", -1.0 / 12.0, -1.0 / 12.0, 1e-9);
    /* Inexact-complex n, exact-Gaussian m. */
    assert_complex_close("Binomial[2. + I, 7 - 3 I]",
                         -75.4683473822230435511, 106.815265970790547809, 1e-4);
    /* Machine-real n, exact-Gaussian m -- the mirror direction. */
    assert_complex_close("Binomial[5., 1 + I]", 3.89232, 6.48720, 1e-3);
    /* Existing inexact-complex path is unchanged. */
    assert_complex_close("N[Binomial[1/2 + I/3, 1/4]]",
                         1.08986784071993926042, 0.0929283046772024343122, 1e-9);

    /* Guards: no inexact operand => stays symbolic (Mathematica-faithful). */
    assert_eval_eq("Head[Binomial[1 + I, 2 + I]]", "Binomial", 0);  /* both exact complex */
    assert_eval_eq("Binomial[1 + I, 5]", "-1/12 - 1/12*I", 0);      /* exact int m: polynomial */
    assert_eval_eq("Head[Binomial[x, 5.]]", "Binomial", 0);         /* symbolic n, no complex */
    assert_eval_eq("Head[Binomial[7/3, 1/5]]", "Binomial", 0);      /* exact rationals */
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_mod_rational_sibling);
    TEST(test_mod_mpfr_sibling);
    TEST(test_mod_regressions);
    TEST(test_divide_rational_mpfr_complex);
    TEST(test_divide_regressions);
    TEST(test_quotientremainder);
    TEST(test_factorial_complex);
    TEST(test_factorial2_numeric);
    TEST(test_factorialpower_numeric);
    TEST(test_barnesg_numeric);
    TEST(test_barnesg_recurrence);
    TEST(test_hyperfactorial_numeric);
    TEST(test_binomial_complex_numeric);

    printf("All numeric-domain tests passed.\n");
    return 0;
}
