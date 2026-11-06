/* Unit tests for NDSolve's evaluated arguments, its arbitrary-precision
 * Taylor-series integrator and high-order dense output, InterpolatingFunction
 * property queries, and the FindMinimum forms they are combined with.
 *
 * Expected values are closed forms (1/(1-x), Exp[Sin[x]], ArcTan, Sec, Cosh,
 * Exp[Exp[x]], ...) or, for y'' = x^3 + y^3 + y'^3, y(0) = y'(0) = 0 (whose
 * y' blows up at x = r), an independent 60-digit computation of
 *     r    = 1.64437669033884072019384083323...
 *     y(r) = 0.93193876511028105335417439414...
 * (inverse function x(y) across the branch point, maximum by the root of
 * x'(y)), which the m1/m2 rows below must reproduce to 13 digits.
 *
 * Hard asserts: the process exits non-zero on the first failure. */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int failures = 0;

static Expr* run(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool eval_double(const char* input, double* out) {
    Expr* r = run(input);
    bool ok = false;
    if (r) {
        if (r->type == EXPR_REAL)         { *out = r->data.real; ok = true; }
        else if (r->type == EXPR_INTEGER) { *out = (double)r->data.integer; ok = true; }
#ifdef USE_MPFR
        else if (r->type == EXPR_MPFR)    { *out = mpfr_get_d(r->data.mpfr, MPFR_RNDN); ok = true; }
#endif
    }
    expr_free(r);
    return ok;
}

static void check_close(const char* input, double expected, double tol) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric\n", input); failures++; return; }
    if (!(fabs(v - expected) <= tol)) {
        printf("FAIL: %s -> %.15g (expected %.15g, tol %.1e)\n", input, v, expected, tol);
        failures++;
    } else printf("ok:   %s = %.12g\n", input, v);
}

/* value >= lo, e.g. a -Log10 error bound or a Precision. */
static void check_ge(const char* input, double lo) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric\n", input); failures++; return; }
    if (!(v >= lo)) { printf("FAIL: %s -> %.6g (expected >= %.6g)\n", input, v, lo); failures++; }
    else printf("ok:   %s = %.6g (>= %.6g)\n", input, v, lo);
}

static void check_eq(const char* input, const char* expected) {
    Expr* r = run(input);
    char* got = expr_to_string_fullform(r);
    if (!got || strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got ? got : "(null)");
        failures++;
    } else printf("ok:   %s == %s\n", input, expected);
    free(got);
    expr_free(r);
}

/* digits of agreement: -Log10 |a - b|, floored at 10^-80 */
#define DIG(a, b) "-Log[10, Abs[(" a ") - (" b ")] + 10^-80]"

/* ------------------------------------------------------------------ */
/* 1. NDSolve evaluates its arguments                                  */
/* ------------------------------------------------------------------ */
static void test_evaluated_arguments(void) {
    check_eq("Attributes[NDSolve]", "List[Protected]");
    /* the system held in a symbol (the HoldAll attribute made this nodeqs) */
    run("odeT = {y'[x] == z[x], z'[x] == -y[x], y[0] == 0, z[0] == 1}");
    check_close("y[1.] /. First[NDSolve[odeT, {y, z}, {x, 0, 2}]]", sin(1.0), 1e-7);
    /* an equation stored in a symbol inside the list, the range in another */
    run("eqT = y''[x] == -y[x]; rgT = {x, 0, 3}");
    check_close("y[2.] /. First[NDSolve[{eqT, y[0] == 1, y'[0] == 0}, y, rgT]]", cos(2.0), 1e-7);
    /* the dependent-variable list from a symbol */
    run("varsT = {y, z}");
    check_close("z[1.] /. First[NDSolve[odeT, varsT, {x, 0, 2}]]", cos(1.0), 1e-7);
    /* NDSolve`Method builtins too */
    check_close("y[1.] /. First[NDSolve`RK4[odeT, {y, z}, {x, 0, 2}]]", sin(1.0), 1e-5);
    /* a solution then used inside a Module (local symbols evaluate to the ODE) */
    check_close("Module[{e = {w'[s] == w[s], w[0] == 1}}, w[1.] /. First[NDSolve[e, w, {s, 0, 1}]]]",
                exp(1.0), 1e-6);
}

/* ------------------------------------------------------------------ */
/* 2. InterpolatingFunction properties                                 */
/* ------------------------------------------------------------------ */
static void test_if_properties(void) {
    run("ifT = Interpolation[{{0, 1.}, {1, 3.}, {2, 7.}, {4, 9.}}]");
    check_eq("ifT[\"Domain\"]", "List[List[0, 4]]");
    check_eq("ifT[\"Coordinates\"]", "List[List[0.0, 1.0, 2.0, 4.0]]");
    check_eq("ifT[\"Grid\"]", "List[List[0.0], List[1.0], List[2.0], List[4.0]]");
    check_eq("ifT[\"ValuesOnGrid\"]", "List[1.0, 3.0, 7.0, 9.0]");
    check_eq("ifT[\"InterpolationOrder\"]", "List[3]");
    check_eq("Length[ifT[\"Properties\"]]", "6");
    check_eq("MatchQ[ifT[\"NoSuchProperty\"], _InterpolatingFunction[_String]]", "True");
    /* exact data keeps exact coordinates */
    run("ifE = Interpolation[{{0, 1}, {1/2, 3}, {1, 7}}]");
    check_eq("ifE[\"Coordinates\"]", "List[List[0, Rational[1, 2], 1]]");
    check_eq("ifE[\"ValuesOnGrid\"]", "List[1, 3, 7]");
    check_eq("ifE[\"InterpolationOrder\"]", "List[2]");
    /* 2-D tensor grid */
    run("if2 = Interpolation[Flatten[Table[{{i, j}, i + 10 j}, {i, 0, 2}, {j, 0, 1}], 1]]");
    check_eq("if2[\"Domain\"]", "List[List[0, 2], List[0, 1]]");
    check_eq("if2[\"Coordinates\"]", "List[List[0, 1, 2], List[0, 1]]");
    check_eq("if2[\"ValuesOnGrid\"]", "List[List[0, 10], List[1, 11], List[2, 12]]");
    check_eq("Dimensions[if2[\"Grid\"]]", "List[3, 2, 2]");
    /* an NDSolve solution: its domain is the integration interval, and a
     * solution stopped at a singularity reports where it stopped */
    check_eq("(y /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 2}]])[\"Domain\"]",
             "List[List[0.0, 2.0]]");
    check_close("Last[First[(y /. First[NDSolve[{y'[x] == z[x], z'[x] == x^3 + y[x]^3 + z[x]^3, "
                "y[0] == 0, z[0] == 0}, {y, z}, {x, 0, 2}]])[\"Domain\"]]]",
                1.6443766903388407, 1e-4);
    /* y' = 1 + y^2, y(0) = 0: y = Tan[x] blows up at Pi/2 */
    check_close("Last[First[(y /. First[NDSolve[{y'[x] == 1 + y[x]^2, y[0] == 0}, y, {x, 0, 3}]])"
                "[\"Domain\"]]]", M_PI / 2, 1e-4);
}

/* ------------------------------------------------------------------ */
/* 3. MPFR Taylor integrator against closed forms                      */
/* ------------------------------------------------------------------ */
#define OPT40 "WorkingPrecision -> 40, PrecisionGoal -> 35, AccuracyGoal -> 35"
static void test_taylor_closed_forms(void) {
#ifdef USE_MPFR
    /* y' = y^2: 1/(1-x), pole at 1 */
    check_ge(DIG("y[9/10] /. First[NDSolve[{y'[x] == y[x]^2, y[0] == 1}, y, {x, 0, 9/10}, " OPT40 "]]", "10"), 33);
    /* the value at an exact point is a 40-digit number */
    check_ge("Precision[y[9/10] /. First[NDSolve[{y'[x] == y[x]^2, y[0] == 1}, y, {x, 0, 9/10}, " OPT40 "]]]", 39);
    /* a machine argument gives a machine number, computed from the 40-digit
     * table (so correctly rounded, not limited by a double Hermite basis) */
    check_eq("Precision[y[0.5] /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 1}, " OPT40 "]]]",
             "MachinePrecision");
    check_close("y[0.5] /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 1}, " OPT40 "]]",
                exp(-0.5), 2e-16);
    /* exp/sin vocabulary, between nodes (dense output) */
    check_ge(DIG("y[13/10] /. First[NDSolve[{y'[x] == Cos[x] y[x], y[0] == 1}, y, {x, 0, 3}, " OPT40 "]]", "Exp[Sin[13/10]]"), 34);
    check_ge(DIG("y[7/5] /. First[NDSolve[{y'[x] == 1/(1 + x^2), y[0] == 0}, y, {x, 0, 2}, " OPT40 "]]", "ArcTan[7/5]"), 34);
    check_ge(DIG("y[3/7] /. First[NDSolve[{y'[x] == Sqrt[y[x]], y[0] == 1}, y, {x, 0, 2}, " OPT40 "]]", "(3/14 + 1)^2"), 34);
    check_ge(DIG("y[5/7] /. First[NDSolve[{y'[x] == y[x] Log[y[x]], y[0] == E}, y, {x, 0, 1}, " OPT40 "]]", "Exp[Exp[5/7]]"), 33);
    check_ge(DIG("y[11/7] /. First[NDSolve[{y'[x] == Tanh[x] y[x], y[0] == 2}, y, {x, 0, 2}, " OPT40 "]]", "2 Cosh[11/7]"), 34);
    check_ge(DIG("y[4/3] /. First[NDSolve[{y'[x] == Tan[x] y[x], y[0] == 1}, y, {x, 0, 3/2}, " OPT40 "]]", "Sec[4/3]"), 33);
    check_ge(DIG("y[13/9] /. First[NDSolve[{y'[x] == Sinh[x], y[0] == 0}, y, {x, 0, 2}, " OPT40 "]]", "Cosh[13/9] - 1"), 34);
    /* second order, 50 digits, a long interval */
    check_ge(DIG("y[37/5] /. First[NDSolve[{y''[x] == -y[x], y[0] == 1, y'[0] == 0}, y, {x, 0, 10}, "
                 "WorkingPrecision -> 50, PrecisionGoal -> 40, AccuracyGoal -> 40]]", "Cos[37/5]"), 39);
    /* the derivative of the solution is as accurate (exact node derivatives) */
    check_ge(DIG("y'[37/5] /. First[NDSolve[{y''[x] == -y[x], y[0] == 1, y'[0] == 0}, y, {x, 0, 10}, "
                 "WorkingPrecision -> 50, PrecisionGoal -> 40, AccuracyGoal -> 40]]", "-Sin[37/5]"), 38);
    /* initial conditions at the working precision, not rounded to double:
     * y' = y, y(1/3) = 1/7 on [1/3, 1] -> e^(2/3)/7 */
    check_ge(DIG("y[1] /. First[NDSolve[{y'[x] == y[x], y[N[1/3, 50]] == N[1/7, 50]}, y, "
                 "{x, N[1/3, 50], 1}, WorkingPrecision -> 50, PrecisionGoal -> 40, AccuracyGoal -> 40]]",
                 "Exp[2/3]/7"), 39);
    /* backward from an interior initial point */
    check_ge(DIG("y[0] /. First[NDSolve[{y'[x] == -2 x y[x], y[1] == Exp[-1]}, y, {x, 0, 2}, " OPT40 "]]", "1"), 34);
    /* a named method is honoured (DOPRI5 at modest goals) */
    check_ge(DIG("y[1] /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 1}, WorkingPrecision -> 30, "
                 "PrecisionGoal -> 12, AccuracyGoal -> 12, Method -> \"ExplicitRungeKutta\"]]", "Exp[-1]"), 11);
    /* the higher-order Hermite table: InterpolationOrder 2K+1 > 3 */
    check_ge("First[(y /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 1}, " OPT40 "]])[\"InterpolationOrder\"]]", 9);
#endif
}

/* ------------------------------------------------------------------ */
/* 4. Machine dense output (exact node derivatives, degree 7)          */
/* ------------------------------------------------------------------ */
static void test_machine_dense_output(void) {
    /* between nodes the error is the node error, not the cubic's O(h^4):
     * max over a fine sample was 1.5e-7 (cubic), now ~1e-9 */
    check_ge("-Log[10, Max[Table[Abs[(y[q] /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}]]) "
             "- Exp[-q]], {q, 0., 5., 0.05}]]]", 8.5);
    check_ge("-Log[10, Max[Table[Abs[(y'[q] /. #) + Sin[q]], {q, 0., 20., 0.05}]]] &[First[NDSolve[{y''[x] == -y[x], "
             "y[0] == 1, y'[0] == 0}, y, {x, 0, 20}]]]", 7.0);
    check_eq("(y /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}]])[\"InterpolationOrder\"]", "List[7]");
    /* InterpolationOrder -> 3 keeps the cubic */
    check_eq("(y /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}, InterpolationOrder -> 3]])"
             "[\"InterpolationOrder\"]", "List[3]");
}

/* ------------------------------------------------------------------ */
/* 5. FindMinimum forms                                                */
/* ------------------------------------------------------------------ */
static void test_findmin_forms(void) {
    /* {x, a, b} with rational starts: two starting values of ONE variable
     * (Times[-1, Power[10, -1]] was read as an indexed variable) */
    check_eq("Length[Last[FindMinimum[(x - 1/30)^2 + 2, {x, -1/10, 1/10}]]]", "1");
    check_close("x /. Last[FindMinimum[(x - 1/30)^2 + 2, {x, -1/10, 1/10}]]", 1.0 / 30.0, 1e-7);
    check_close("x /. Last[FindMinimum[(x - Pi/7)^2, {x, Pi/9, Pi/5}]]", M_PI / 7, 1e-7);
    /* two starts that do not bracket the minimum: the search leaves [a, b] */
    check_close("x /. Last[FindMinimum[(x - 3)^2, {x, 1/2, 1}]]", 3.0, 1e-7);
    /* a numeric start is not a variable in a list of specs either */
    check_close("y /. Last[FindMinimum[(x - 1)^2 + (y + 2)^2, {{x, 1/3}, {y, -1/3}}]]", -2.0, 1e-6);
#ifdef USE_MPFR
    /* the minimum value carries the working precision */
    check_ge("Precision[First[FindMinimum[(t - 1)^2 + 1, {t, 0}, WorkingPrecision -> 30]]]", 29.9);
    check_ge("Precision[First[FindMinimum[(t - 1)^2 + 1/3, {t, -1/10, 1/10}, WorkingPrecision -> 30]]]", 29.9);
    check_ge(DIG("First[FindMinimum[(t - 1)^2 + 1/3, {t, -1/10, 1/10}, WorkingPrecision -> 30]]", "1/3"), 29);
#endif
}

/* ------------------------------------------------------------------ */
/* 6. The ODE singularity problem end to end                           */
/* ------------------------------------------------------------------ */
#define R_TRUE  "1.64437669033884072019384083323"
#define YR_TRUE "0.93193876511028105335417439414"
static void test_singularity_problem(void) {
#ifdef USE_MPFR
    run("odeS = {y'[x] == z[x], z'[x] == x^3 + y[x]^3 + z[x]^3, y[0] == 0, z[0] == 0}");
    run("m1S[wp_] := Module[{s1, y0, z0, s2, fm, o = {WorkingPrecision -> wp, PrecisionGoal -> 30, "
        "AccuracyGoal -> 20, MaxSteps -> 3000}}, "
        "s1 = NDSolve[odeS, {y, z}, {x, 0, 6/5}, StartingStepSize -> 10^-3, Sequence @@ o]; "
        "{y0, z0} = {y[6/5], z[6/5]} /. First[s1]; "
        "s2 = NDSolve[{-xx''[t] == (xx[t]^3 + t^3) xx'[t]^3 + 1, xx[y0] == 6/5, xx'[y0] == 1/z0}, "
        "xx, {t, y0, 1}, Sequence @@ o]; "
        "fm = FindMinimum[Evaluate[-xx[t] /. First[s2]], {t, 9/10, 1}, WorkingPrecision -> wp, "
        "PrecisionGoal -> 30, AccuracyGoal -> 30]; {-First[fm], t /. Last[fm]}]");
    run("m1R = m1S[40]");
    check_ge(DIG("m1R[[1]]", R_TRUE), 20);
    check_ge(DIG("m1R[[2]]", YR_TRUE), 13);
    /* polar coordinates about (0, 1) */
    run("polS = Module[{X, Y, yp, ypp}, X = rho[p] Cos[p]; Y = rho[p] Sin[p] + 1; "
        "yp = D[Y, p]/D[X, p]; ypp = D[yp, p]/D[X, p]; "
        "Factor[Numerator[Together[ypp - (X^3 + Y^3 + yp^3)]]] == 0]");
    run("m2S[wp_] := Module[{s, fm, p0}, s = NDSolve[{polS, rho[-Pi/2] == 1, rho'[-Pi/2] == 0}, rho, "
        "{p, -Pi/2, Pi/4}, WorkingPrecision -> wp, PrecisionGoal -> wp, AccuracyGoal -> 20, MaxSteps -> 5000]; "
        "fm = FindMinimum[Evaluate[-rho[p] Cos[p] /. First[s]], {p, -1/10, 1/10}, WorkingPrecision -> wp, "
        "AccuracyGoal -> 80]; p0 = p /. Last[fm]; {-First[fm], (rho[p0] Sin[p0] + 1) /. First[s]}]");
    run("m2R = m2S[30]");
    check_ge(DIG("m2R[[1]]", R_TRUE), 20);
    check_ge(DIG("m2R[[2]]", YR_TRUE), 13.5);
    /* machine precision: both methods to ~7 digits */
    run("m2M = Module[{s, fm, p0}, s = NDSolve[{polS, rho[-Pi/2] == 1, rho'[-Pi/2] == 0}, rho, "
        "{p, -Pi/2, Pi/4}, MaxSteps -> 5000]; fm = FindMinimum[Evaluate[-rho[p] Cos[p] /. First[s]], "
        "{p, -1/10, 1/10}]; p0 = p /. Last[fm]; {-First[fm], (rho[p0] Sin[p0] + 1) /. First[s]}]");
    check_ge(DIG("m2M[[1]]", R_TRUE), 8.5);
    check_ge(DIG("m2M[[2]]", YR_TRUE), 7);
#endif
}

/* ------------------------------------------------------------------ */
/* 7. ODE singularities at working precision: closed forms checked to   */
/*    the requested digits (Runge-Kutta start data, Abs kinks, blow-up, */
/*    a singular endpoint, extrapolation past the data)                 */
/* ------------------------------------------------------------------ */
static void test_wp_singular_family(void) {
#ifdef USE_MPFR
    /* The Runge-Kutta path at working precision reads its start, initial
     * values and end at that precision and uses the exact DOPRI5 tableau: an
     * exact y(0) = 1/3 used to be rounded to a double (16 digits at most). */
    check_ge(DIG("y[1] /. First[NDSolve[{y'[x] == -y[x], y[0] == 1/3}, y, {x, 0, 1}, WorkingPrecision -> 30, "
                 "PrecisionGoal -> 20, AccuracyGoal -> 20, Method -> \"ExplicitRungeKutta\"]]", "Exp[-1]/3"), 19);
    /* y' = y^2, y(1/7) = 7/13: y = 1/(2 - x) */
    check_ge(DIG("y[3/2] /. First[NDSolve[{y'[x] == y[x]^2, y[1/7] == 7/13}, y, {x, 1/7, 3/2}, WorkingPrecision -> 30, "
                 "PrecisionGoal -> 20, AccuracyGoal -> 20, Method -> \"ExplicitRungeKutta\"]]", "2"), 18);
    /* Lane-Emden n = 1 written with Abs: y = Sin[x]/x up to the kink at Pi,
     * and y = -Sinh[x - Pi]/x past it (|y| = -y there).  The step lands on
     * the zero of the Abs argument, so the series never straddles it. */
    run("leS[wp_, pg_] := NDSolve[{y''[x] + 2 y'[x]/x + Abs[y[x]] == 0, y[1/10] == Sin[1/10] 10, "
        "y'[1/10] == (Cos[1/10]/10 - Sin[1/10]) 100}, y, {x, 1/10, 7/2}, WorkingPrecision -> wp, "
        "PrecisionGoal -> pg, AccuracyGoal -> pg]");
    run("le30 = leS[30, 20]; le50 = leS[50, 40];");
    check_ge(DIG("x /. FindRoot[Evaluate[y[x] /. First[le30]], {x, 3}, WorkingPrecision -> 30, "
                 "AccuracyGoal -> 25, PrecisionGoal -> 25]", "Pi"), 19);
    check_ge(DIG("x /. FindRoot[Evaluate[y[x] /. First[le50]], {x, 3}, WorkingPrecision -> 50, "
                 "AccuracyGoal -> 45, PrecisionGoal -> 45]", "Pi"), 39);
    check_ge(DIG("y[3] /. First[le50]", "Sin[3]/3"), 39);
    check_ge(DIG("y[7/2] /. First[le30]", "-Sinh[7/2 - Pi]/(7/2)"), 19);
    check_ge(DIG("y[7/2] /. First[le50]", "-Sinh[7/2 - Pi]/(7/2)"), 39);
    check_ge(DIG("y'[33/10] /. First[le50]", "(D[-Sinh[t - Pi]/t, t] /. t -> 33/10)"), 38);
    /* Lane-Emden n = 3 (Abs[y]^3, a sign flip at the kink) and n = 3/2
     * (Abs[y]^(3/2): no series at the zero, crossed by Runge-Kutta): first
     * zeros xi_1 (published to 27 digits) from the series start at 10^-6 */
    run("leN[n_, xe_, g_] := Module[{x0 = 10^-6, s}, s = NDSolve[{y''[x] + 2 y'[x]/x + Abs[y[x]]^n == 0, "
        "y[x0] == 1 - x0^2/6 + n x0^4/120, y'[x0] == -x0/3 + n x0^3/30}, y, {x, x0, xe}, "
        "WorkingPrecision -> 30, PrecisionGoal -> 20, AccuracyGoal -> 20]; "
        "x /. FindRoot[Evaluate[y[x] /. First[s]], {x, g}, WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 25]]");
    check_ge(DIG("leN[3, 7, 69/10]", "6.896848619376960375454528187`30"), 19);
    check_ge(DIG("leN[3/2, 37/10, 365/100]", "3.653753736219122424609094280`30"), 19);
    /* Riccati blow-up y' = 1 + y^2 (y = Tan[x]): the naive solution stops
     * short of Pi/2; the reciprocal u = 1/y, u' = -(1 + u^2) from x = 1 is
     * regular through it and its zero is Pi/2 to the working precision. */
    run("ric[wp_] := Module[{s1, y1, s2}, s1 = NDSolve[{y'[x] == 1 + y[x]^2, y[0] == 0}, y, {x, 0, 1}, "
        "WorkingPrecision -> wp, PrecisionGoal -> wp - 10, AccuracyGoal -> wp - 10]; y1 = y[1] /. First[s1]; "
        "s2 = NDSolve[{u'[x] == -(1 + u[x]^2), u[1] == 1/y1}, u, {x, 1, 2}, WorkingPrecision -> wp, "
        "PrecisionGoal -> wp - 10, AccuracyGoal -> wp - 10]; "
        "x /. FindRoot[Evaluate[u[x] /. First[s2]], {x, 3/2}, WorkingPrecision -> wp, AccuracyGoal -> wp - 5, PrecisionGoal -> wp - 5]]");
    check_ge(DIG("ric[30]", "Pi/2"), 19);
    check_ge(DIG("ric[50]", "Pi/2"), 39);
    check_ge(DIG("Quiet[Last[First[(y /. First[NDSolve[{y'[x] == 1 + y[x]^2, y[0] == 0}, y, {x, 0, 2}, "
                 "WorkingPrecision -> 30, PrecisionGoal -> 20, AccuracyGoal -> 20]])[\"Domain\"]]]]", "Pi/2"), 5);
    /* Bessel J0 from exact data at x0 = 1/1000: y(10) = BesselJ[0, 10] */
    check_ge(DIG("y[10] /. First[NDSolve[{x y''[x] + y'[x] + x y[x] == 0, y[1/1000] == BesselJ[0, 1/1000], "
                 "y'[1/1000] == -BesselJ[1, 1/1000]}, y, {x, 1/1000, 10}, WorkingPrecision -> 50, "
                 "PrecisionGoal -> 40, AccuracyGoal -> 40]]", "BesselJ[0, 10]"), 39);
    /* Gauss hypergeometric equation integrated to its singular point x = 1
     * (y ~ A + B (1 - x)^(7/6)): the steps shrink geometrically, the solution
     * stops a few hundred ulps short of 1 at the returned precision (nodes no
     * longer collapse onto 1), and y[1] is extrapolated from the end node's
     * own series: 2F1(1/3, 1/2; 2; 1). */
    run("gS[wp_] := NDSolve[{x (1 - x) y''[x] + (2 - 11/6 x) y'[x] - 1/6 y[x] == 0, "
        "y[1/1000] == Hypergeometric2F1[1/3, 1/2, 2, 1/1000], "
        "y'[1/1000] == 1/12 Hypergeometric2F1[4/3, 3/2, 3, 1/1000]}, y, {x, 1/1000, 1}, "
        "WorkingPrecision -> wp, PrecisionGoal -> wp - 10, AccuracyGoal -> wp - 10]");
    run("g30 = Quiet[gS[30]]; g50 = Quiet[gS[50]];");
    check_ge(DIG("Quiet[y[1] /. First[g30]]", "Hypergeometric2F1[1/3, 1/2, 2, 1]"), 19);
    check_ge(DIG("Quiet[y[1] /. First[g50]]", "Hypergeometric2F1[1/3, 1/2, 2, 1]"), 39);
    check_ge(DIG("y[1/2] /. First[g50]", "Hypergeometric2F1[1/3, 1/2, 2, 1/2]"), 39);
    /* the domain ends short of 1, but within 10^-20 of it */
    check_ge("-Log[10, 1 - Last[First[(y /. First[g30])[\"Domain\"]]]]", 20);
    check_eq("Sign[1 - Last[First[(y /. First[g30])[\"Domain\"]]]]", "1");
    /* evaluating past the end says so (dmval, which Check sees) */
    check_eq("Check[y[1] /. First[g30], \"ext\"]", "\"ext\"");
    /* ...while a table of an entire function keeps the last interval's
     * interpolant a width or so out, which is far better there than the end
     * node's Taylor polynomial (~1e-31 against ~1e-14 for Exp at 11/10). */
    run("ipE = Interpolation[Table[{{x}, Sequence @@ Table[N[Exp[x], 40], {9}]}, {x, 0, 1, 1/10}]];");
    check_ge(DIG("Quiet[ipE[11/10]]", "Exp[11/10]"), 29);
    check_ge(DIG("Quiet[ipE[-1/20]]", "Exp[-1/20]"), 33);
#endif
}

int main(void) {
    test_mute_messages();
    symtab_init();
    core_init();

    TEST(test_evaluated_arguments);
    TEST(test_if_properties);
    TEST(test_taylor_closed_forms);
    TEST(test_machine_dense_output);
    TEST(test_findmin_forms);
    TEST(test_singularity_problem);
    TEST(test_wp_singular_family);

    if (failures == 0) printf("\nAll NDSolve Taylor / dense-output / property tests passed.\n");
    else printf("\n%d checks FAILED.\n", failures);
    return failures == 0 ? 0 : 1;
}
