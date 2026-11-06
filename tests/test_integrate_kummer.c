/* test_integrate_kummer.c
 *
 * The KummerLogGamma definite-integration stage
 * (src/internal/KummerLogGamma.m, wired in src/calculus/integrate.c):
 * INT_0^1 LogGamma[x] T(x) dx for T a trigonometric polynomial of period 1,
 * from the moments of Kummer's Fourier series
 *   LogGamma[x] = Log[2 Pi]/2 + Sum_n Cos[2 Pi n x]/(2 n)
 *                 + (EulerGamma + Log[2 Pi n]) Sin[2 Pi n x]/(Pi n).
 * Expected values follow from that series by orthogonality
 * (INT Cos^2 = INT Sin^2 = 1/2 over a period) and are cross-checked against a
 * 30-digit NIntegrate, an independent quadrature of the integrand itself.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* Integrate[f, {x, 0, 1}] is exact, equals the closed form `ref` derived from
 * the series, and agrees with a 30-digit quadrature of f to 20 digits. */
static void check_kummer(const char* f, const char* ref) {
    char buf[1024];
    snprintf(buf, sizeof buf, "FreeQ[Integrate[%s, {x, 0, 1}], Integrate]", f);
    check_eq(buf, "True");
    snprintf(buf, sizeof buf,
             "PossibleZeroQ[Integrate[%s, {x, 0, 1}] - (%s)]", f, ref);
    check_eq(buf, "True");
    snprintf(buf, sizeof buf,
             "Abs[N[%s, 30] - NIntegrate[%s, {x, 0, 1}, WorkingPrecision -> 30]] < 10^-20",
             ref, f);
    check_eq(buf, "True");
}

static void test_first_moments(void) {
    /* the ISSAC '97 holdout case: the Sin coefficient of n = 1 halved */
    check_kummer("LogGamma[x] Sin[2 Pi x]", "(EulerGamma + Log[2 Pi])/(2 Pi)");
    check_kummer("LogGamma[x] Cos[2 Pi x]", "1/4");
    /* Raabe, in the LogGamma spelling the C Raabe stage does not read */
    check_kummer("LogGamma[x]", "Log[2 Pi]/2");
    check_kummer("Log[Gamma[x]] Sin[6 Pi x]", "(EulerGamma + Log[6 Pi])/(6 Pi)");
    check_kummer("LogGamma[x] Cos[10 Pi x]", "1/20");
}

static void test_trig_polynomials(void) {
    /* Sin^2 = (1 - Cos[4 Pi x])/2: Raabe/2 - (1/8)/2 */
    check_kummer("Log[Gamma[x]] Sin[2 Pi x]^2", "Log[2 Pi]/4 - 1/16");
    /* Sin[4 Pi x] Cos[2 Pi x] = (Sin[6 Pi x] + Sin[2 Pi x])/2 */
    check_kummer("LogGamma[x] Sin[4 Pi x] Cos[2 Pi x]",
                 "((EulerGamma + Log[6 Pi])/(6 Pi) + (EulerGamma + Log[2 Pi])/(2 Pi))/2");
    /* a term free of LogGamma is integrated by the cascade */
    check_kummer("x + LogGamma[x] (1 + Sin[2 Pi x])",
                 "1/2 + Log[2 Pi]/2 + (EulerGamma + Log[2 Pi])/(2 Pi)");
    /* symbolic coefficient */
    check_eq("PossibleZeroQ[Integrate[a LogGamma[x] Cos[4 Pi x], {x, 0, 1}] - a/8]", "True");
}

static void test_reflection(void) {
    /* x -> 1 - x: Sin[6 Pi (1 - x)] = -Sin[6 Pi x] */
    check_kummer("LogGamma[1 - x] Sin[6 Pi x]", "-(EulerGamma + Log[6 Pi])/(6 Pi)");
    check_kummer("LogGamma[1 - x] Cos[2 Pi x]", "1/4");
}

static void test_declines(void) {
    /* not a period-1 trigonometric weight: no Kummer moment */
    check_eq("KummerLogGamma`Run[LogGamma[x] Sin[x], x]", "$Failed");
    check_eq("KummerLogGamma`Run[LogGamma[x] x, x]", "$Failed");
    check_eq("KummerLogGamma`Run[LogGamma[x] Sin[Pi x], x]", "$Failed");
    /* degree two in LogGamma, another argument */
    check_eq("KummerLogGamma`Run[LogGamma[x]^2, x]", "$Failed");
    check_eq("KummerLogGamma`Run[LogGamma[2 x] Sin[2 Pi x], x]", "$Failed");
    /* another interval is not this stage's: the series holds on (0, 1) */
    check_eq("FreeQ[Integrate[LogGamma[x] Sin[2 Pi x], {x, 1, 2}], Integrate]", "False");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_first_moments);       /* loads the package through the stage */
    TEST(test_trig_polynomials);
    TEST(test_reflection);
    TEST(test_declines);

    printf("All Integrate KummerLogGamma tests passed!\n");
    return 0;
}
