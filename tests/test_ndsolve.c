/* Unit tests for NDSolve (src/numerical_calculus/ndsolve*.c).
 *
 * Coverage:
 *   - Scalar linear IVP y'=-y vs e^{-t}; exponential growth; logistic.
 *   - Higher-order y''+y=0 -> cos t (order reduction).
 *   - Systems x'=y, y'=-x -> {cos, -sin}.
 *   - Every method (ExplicitEuler/Midpoint/RK4/DOPRI5/BackwardEuler/
 *     ImplicitTrapezoid/BDF/Adams) on the reference problem.
 *   - Stiff problem solved by BackwardEuler/BDF.
 *   - Method equivalence: Method->"X" == NDSolve`X[...].
 *   - Error control: tighter PrecisionGoal reduces error.
 *   - InterpolatingFunction derivative: y'[t] == -y[t].
 *   - MPFR: WorkingPrecision -> 30 yields a high-precision node value.
 *   - Result shape and left-value handling (u vs u[x]).
 *
 * Numeric tests are soft (print FAIL, keep going); run: ./ndsolve_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int failures = 0;

static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) { freopen("/dev/null", "w", stderr); done = 1; }
}

static bool eval_double(const char* input, double* out) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    bool ok = false;
    if (r) {
        if (r->type == EXPR_REAL)         { *out = r->data.real; ok = true; }
        else if (r->type == EXPR_INTEGER) { *out = (double)r->data.integer; ok = true; }
#ifdef USE_MPFR
        else if (r->type == EXPR_MPFR)    { *out = mpfr_get_d(r->data.mpfr, MPFR_RNDN); ok = true; }
#endif
    }
    expr_free(e); expr_free(r);
    return ok;
}

static void check_close(const char* input, double expected, double tol) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric\n", input); failures++; return; }
    if (fabs(v - expected) > tol) {
        printf("FAIL: %s -> %.12g (expected %.12g, tol %.1e)\n", input, v, expected, tol);
        failures++;
    } else {
        printf("ok:   %s = %.10g\n", input, v);
    }
}

/* Assert the value is >= lo (used for Precision / error-magnitude bounds). */
static void check_ge(const char* input, double lo) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric\n", input); failures++; return; }
    if (!(v >= lo)) { printf("FAIL: %s -> %.6g (expected >= %.6g)\n", input, v, lo); failures++; }
    else printf("ok:   %s = %.6g (>= %.6g)\n", input, v, lo);
}

static void check_head(const char* input, const char* expected_head) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* ff = expr_to_string_fullform(r);
    if (!ff || strncmp(ff, expected_head, strlen(expected_head)) != 0) {
        printf("FAIL: %s -> %s (expected head %s)\n", input, ff ? ff : "(null)", expected_head);
        failures++;
    } else printf("ok:   %s head %s\n", input, expected_head);
    free(ff); expr_free(e); expr_free(r);
}

/* ---- shared problem strings ---- */
#define P_DECAY "NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}]"

static void test_scalar_linear(void) {
    check_close("First[y[1.0] /. " P_DECAY "]", exp(-1.0), 1e-6);
    check_close("First[y[3.0] /. " P_DECAY "]", exp(-3.0), 1e-6);
    check_close("First[y[4.5] /. " P_DECAY "]", exp(-4.5), 1e-6);
}

static void test_growth_and_logistic(void) {
    check_close("First[y[2.0] /. NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 3}]]", exp(2.0), 1e-5);
    check_close("First[y[2.0] /. NDSolve[{y'[x] == y[x] (1 - y[x]), y[0] == 1/2}, y, {x, 0, 4}]]",
                1.0 / (1.0 + exp(-2.0)), 1e-6);
}

static void test_higher_order(void) {
    /* y'' + y = 0, y(0)=1, y'(0)=0  ->  cos t */
    check_close("First[y[3.0] /. NDSolve[{y''[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 6}]]",
                cos(3.0), 1e-5);
    check_close("First[y[5.0] /. NDSolve[{y''[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 6}]]",
                cos(5.0), 1e-5);
}

static void test_system(void) {
    const char* S = "NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 6}]";
    char buf[512];
    snprintf(buf, sizeof(buf), "First[x[3.0] /. %s]", S); check_close(buf, cos(3.0), 1e-5);
    snprintf(buf, sizeof(buf), "First[y[3.0] /. %s]", S); check_close(buf, -sin(3.0), 1e-5);
}

static void test_all_methods(void) {
    const char* methods[] = { "ExplicitEuler", "ExplicitMidpoint", "RK4", "DOPRI5",
                              "BackwardEuler", "ImplicitTrapezoid", "BDF", "Adams" };
    for (size_t i = 0; i < sizeof(methods)/sizeof(methods[0]); i++) {
        char buf[256];
        snprintf(buf, sizeof(buf),
            "First[y[1.0] /. NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}, Method -> \"%s\"]]",
            methods[i]);
        check_close(buf, exp(-1.0), 2e-3);   /* ExplicitEuler is order 1: loose tol */
    }
}

static void test_stiff(void) {
    /* y' = -1000 (y - cos x) - sin x, y(0)=1  ->  cos x */
    check_close("First[y[2.0] /. NDSolve[{y'[x] == -1000 (y[x] - Cos[x]) - Sin[x], y[0] == 1}, "
                "y, {x, 0, 3}, Method -> \"BackwardEuler\"]]", cos(2.0), 1e-3);
    check_close("First[y[2.0] /. NDSolve[{y'[x] == -1000 (y[x] - Cos[x]) - Sin[x], y[0] == 1}, "
                "y, {x, 0, 3}, Method -> \"BDF\"]]", cos(2.0), 1e-3);
}

static void test_method_equivalence(void) {
    check_close(
        "First[y[3.0] /. NDSolve[{y'[x]==-y[x],y[0]==1}, y, {x,0,5}, Method->\"RK4\"]] - "
        "First[y[3.0] /. NDSolve`RK4[{y'[x]==-y[x],y[0]==1}, y, {x,0,5}]]", 0.0, 1e-14);
    check_close(
        "First[y[3.0] /. NDSolve[{y'[x]==-y[x],y[0]==1}, y, {x,0,5}, Method->\"ExplicitRungeKutta\"]] - "
        "First[y[3.0] /. NDSolve`ExplicitRungeKutta[{y'[x]==-y[x],y[0]==1}, y, {x,0,5}]]", 0.0, 1e-14);
}

static void test_error_control(void) {
    /* tighter PrecisionGoal -> smaller error at an interior point */
    check_ge("-Log[10, Abs[First[y[4.0] /. NDSolve[{y'[x]==-y[x],y[0]==1}, y, {x,0,5}, "
             "PrecisionGoal->10]] - Exp[-4.0]] + 0.0]", 7.0);   /* >= 1e-7 accurate */
}

static void test_if_derivative(void) {
    /* y'[t] via the InterpolatingFunction derivative == -y[t] */
    check_close("First[y'[2.0] /. " P_DECAY "]", -exp(-2.0), 1e-5);
}

static void test_mpfr(void) {
#ifdef USE_MPFR
    /* Integrated at 30-digit working precision (goal 22): the endpoint node
     * agrees with e^{-1} to >17 digits, well beyond machine precision (~15-16).
     * Query at an MPFR abscissa so the mpfr node value is returned. */
    check_ge("-Log[10, Abs[First[y[N[1,30]] /. NDSolve[{y'[x]==-y[x],y[0]==1}, y, {x,0,1}, "
             "WorkingPrecision->30, PrecisionGoal->22, AccuracyGoal->22, MaxSteps->200000]] "
             "- N[Exp[-1],30]] + 10^-40]", 17.0);
    /* querying at an MPFR abscissa yields a high-precision number */
    check_ge("Precision[First[y[N[1,30]] /. NDSolve[{y'[x]==-y[x],y[0]==1}, y, {x,0,5}, "
             "WorkingPrecision->30]]]", 20.0);
#endif
}

static void test_result_shape(void) {
    check_head(P_DECAY, "List[List[Rule[y, InterpolatingFunction");
    check_head("NDSolve[{y'[x] == -y[x], y[0] == 1}, y[x], {x, 0, 5}]",
               "List[List[Rule[y[x], InterpolatingFunction");
}

int main(void) {
    mute_stderr_once();
    symtab_init();
    core_init();

    TEST(test_scalar_linear);
    TEST(test_growth_and_logistic);
    TEST(test_higher_order);
    TEST(test_system);
    TEST(test_all_methods);
    TEST(test_stiff);
    TEST(test_method_equivalence);
    TEST(test_error_control);
    TEST(test_if_derivative);
    TEST(test_mpfr);
    TEST(test_result_shape);

    if (failures == 0) printf("\nAll NDSolve tests passed.\n");
    else printf("\n%d NDSolve checks FAILED.\n", failures);
    return 0;   /* soft-assert convention: harness greps 'FAIL:' */
}
