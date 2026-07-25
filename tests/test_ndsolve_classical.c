/* NDSolve on classical ODE (and method-of-lines PDE) problems, compared to
 * exact solutions.
 *
 * NDSolve is an ODE initial-value solver; PDEs are covered by discretizing the
 * spatial operator (method of lines) into an ODE system and comparing to the
 * exact solution of the *semi-discrete* problem (an exact eigenmode), which
 * isolates the time-integration accuracy from spatial-discretization error.
 *
 * The default method is Dormand-Prince 5(4) with FSAL — the same algorithm as
 * MATLAB's ode45 — so accuracy and cost are directly comparable (this suite
 * checks accuracy against closed forms; the unit suite checks the ~6-evals/step
 * FSAL efficiency).
 *
 * Soft asserts: prints FAIL and keeps going. Run: ./ndsolve_classical_tests */
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

static const double PI = 3.14159265358979323846;
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

/* evaluate a side-effecting statement (e.g. an assignment) and discard */
static void run(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    expr_free(e); expr_free(r);
}

static void check_close_msg(const char* label, const char* input, double expected, double tol) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric [%s]\n", label, input); failures++; return; }
    if (fabs(v - expected) > tol) {
        printf("FAIL: %s -> %.12g (expected %.12g, |err|=%.2e, tol %.1e)\n",
               label, v, expected, fabs(v - expected), tol);
        failures++;
    } else {
        printf("ok:   %-42s = %.10g  (|err|=%.1e)\n", label, v, fabs(v - expected));
    }
}
#define CHECK(label, input, expected, tol) check_close_msg(label, input, expected, tol)

static void check_le(const char* label, const char* input, double hi) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric\n", label); failures++; return; }
    if (!(v <= hi)) { printf("FAIL: %s -> %.6g (expected <= %.6g)\n", label, v, hi); failures++; }
    else printf("ok:   %-42s = %.6g (<= %.6g)\n", label, v, hi);
}

/* ================= first-order ODEs ================= */
static void test_first_order(void) {
    /* exponential decay / growth: y' = a y  ->  y0 e^{a t} */
    CHECK("decay y'=-y", "First[y[3.0]/.NDSolve[{y'[x]==-y[x],y[0]==1},y,{x,0,5}]]", exp(-3.0), 1e-6);
    CHECK("growth y'=2y", "First[y[2.0]/.NDSolve[{y'[x]==2 y[x],y[0]==1},y,{x,0,3}]]", exp(4.0), 1e-4);

    /* Newton cooling: y' = -k(y-A) -> A + (y0-A)e^{-kt} */
    CHECK("Newton cooling",
          "First[y[4.0]/.NDSolve[{y'[x]==-(1/2)(y[x]-20),y[0]==100},y,{x,0,8}]]",
          20.0 + 80.0*exp(-2.0), 1e-4);

    /* linear with sinusoidal forcing: y'+y=sin t, y(0)=0 -> (sin t - cos t + e^{-t})/2 */
    CHECK("forced linear y'+y=sin",
          "First[y[3.0]/.NDSolve[{y'[x]==Sin[x]-y[x],y[0]==0},y,{x,0,6}]]",
          (sin(3.0) - cos(3.0) + exp(-3.0)) / 2.0, 1e-6);

    /* logistic: y'=y(1-y), y0=1/2 -> 1/(1+e^{-t}) */
    CHECK("logistic",
          "First[y[2.0]/.NDSolve[{y'[x]==y[x](1-y[x]),y[0]==1/2},y,{x,0,5}]]",
          1.0/(1.0+exp(-2.0)), 1e-6);

    /* y'=y^2, y0=1 -> 1/(1-t), sampled before the t=1 blow-up */
    CHECK("y'=y^2 (finite-time blowup)",
          "First[y[0.8]/.NDSolve[{y'[x]==y[x]^2,y[0]==1},y,{x,0,0.9}]]",
          1.0/(1.0-0.8), 1e-4);

    /* separable non-autonomous: y'=cos(t) y -> e^{sin t} */
    CHECK("separable y'=cos(t)y",
          "First[y[2.0]/.NDSolve[{y'[x]==Cos[x] y[x],y[0]==1},y,{x,0,4}]]",
          exp(sin(2.0)), 1e-5);

    /* Riccati y'=1+y^2 -> tan t (near the t=pi/2 singularity) */
    CHECK("Riccati y'=1+y^2 -> tan",
          "First[y[1.4]/.NDSolve[{y'[x]==1+y[x]^2,y[0]==0},y,{x,0,1.5}]]",
          tan(1.4), 1e-5);

    /* exact eqn: y'=2t(1+y), y0=0 -> e^{t^2}-1 */
    CHECK("exact y'=2t(1+y)",
          "First[y[1.3]/.NDSolve[{y'[x]==2 x (1+y[x]),y[0]==0},y,{x,0,1.5}]]",
          exp(1.3*1.3)-1.0, 1e-4);
}

/* ================= second-order linear ================= */
static void test_second_order(void) {
    /* SHO: y''+y=0 -> cos / sin */
    CHECK("SHO cos", "First[y[3.0]/.NDSolve[{y''[x]+y[x]==0,y[0]==1,y'[0]==0},y,{x,0,6}]]", cos(3.0), 1e-5);
    CHECK("SHO sin", "First[y[3.0]/.NDSolve[{y''[x]+y[x]==0,y[0]==0,y'[0]==1},y,{x,0,6}]]", sin(3.0), 1e-5);
    CHECK("SHO omega=2", "First[y[2.0]/.NDSolve[{y''[x]+4 y[x]==0,y[0]==1,y'[0]==0},y,{x,0,5}]]", cos(4.0), 1e-5);

    /* underdamped: y''+2y'+5y=0, y0=1,y'0=0 -> e^{-t}(cos2t + 1/2 sin2t) */
    CHECK("underdamped",
          "First[y[1.5]/.NDSolve[{y''[x]+2 y'[x]+5 y[x]==0,y[0]==1,y'[0]==0},y,{x,0,4}]]",
          exp(-1.5)*(cos(3.0)+0.5*sin(3.0)), 1e-6);

    /* critically damped: y''+2y'+y=0 -> (1+t)e^{-t} */
    CHECK("critically damped",
          "First[y[2.0]/.NDSolve[{y''[x]+2 y'[x]+y[x]==0,y[0]==1,y'[0]==0},y,{x,0,5}]]",
          (1.0+2.0)*exp(-2.0), 1e-6);

    /* overdamped: y''+5y'+6y=0 -> 3e^{-2t}-2e^{-3t} */
    CHECK("overdamped",
          "First[y[1.0]/.NDSolve[{y''[x]+5 y'[x]+6 y[x]==0,y[0]==1,y'[0]==0},y,{x,0,4}]]",
          3.0*exp(-2.0)-2.0*exp(-3.0), 1e-6);

    /* nonhomogeneous: y''+y=t, y0=0,y'0=0 -> t - sin t */
    CHECK("nonhomog y''+y=t",
          "First[y[3.0]/.NDSolve[{y''[x]+y[x]==x,y[0]==0,y'[0]==0},y,{x,0,6}]]",
          3.0 - sin(3.0), 1e-5);

    /* resonance: y''+y=cos t -> (t/2) sin t */
    CHECK("resonance y''+y=cos t",
          "First[y[4.0]/.NDSolve[{y''[x]+y[x]==Cos[x],y[0]==0,y'[0]==0},y,{x,0,6}]]",
          (4.0/2.0)*sin(4.0), 1e-4);
}

/* ================= higher order & variable coefficient ================= */
static void test_higher_and_varcoef(void) {
    /* third order y'''=y', ICs pick cosh: y(0)=1,y'(0)=0,y''(0)=1 -> cosh t */
    CHECK("3rd order -> cosh",
          "First[y[2.0]/.NDSolve[{y'''[x]==y'[x],y[0]==1,y'[0]==0,y''[0]==1},y,{x,0,3}]]",
          cosh(2.0), 1e-5);

    /* Cauchy-Euler: x^2 y'' + x y' - y = 0 on [1,3], y(1)=2,y'(1)=0 -> t + 1/t */
    CHECK("Cauchy-Euler -> t+1/t",
          "First[y[2.5]/.NDSolve[{x^2 y''[x]+x y'[x]-y[x]==0,y[1]==2,y'[1]==0},y,{x,1,3}]]",
          2.5 + 1.0/2.5, 1e-5);

    /* fourth order y''''=y, ICs pick cosh: cosh satisfies y''''=y */
    CHECK("4th order -> cosh",
          "First[y[1.5]/.NDSolve[{y''''[x]==y[x],y[0]==1,y'[0]==0,y''[0]==1,y'''[0]==0},y,{x,0,3}]]",
          cosh(1.5), 1e-5);
}

/* ================= systems ================= */
static void test_systems(void) {
    const char* HARM = "NDSolve[{x'[t]==y[t],y'[t]==-x[t],x[0]==1,y[0]==0},{x,y},{t,0,8}]";
    char buf[512];
    snprintf(buf, sizeof buf, "First[x[5.0]/.%s]", HARM); CHECK("harmonic x=cos", buf, cos(5.0), 1e-5);
    snprintf(buf, sizeof buf, "First[y[5.0]/.%s]", HARM); CHECK("harmonic y=-sin", buf, -sin(5.0), 1e-5);

    /* decoupled: x'=x, y'=-2y -> {e^t, e^{-2t}} */
    const char* DEC = "NDSolve[{x'[t]==x[t],y'[t]==-2 y[t],x[0]==1,y[0]==1},{x,y},{t,0,3}]";
    snprintf(buf, sizeof buf, "First[x[2.0]/.%s]", DEC); CHECK("decoupled x=e^t", buf, exp(2.0), 1e-4);
    snprintf(buf, sizeof buf, "First[y[2.0]/.%s]", DEC); CHECK("decoupled y=e^-2t", buf, exp(-4.0), 1e-6);

    /* coupled linear: x'=y, y'=-2x-3y -> x = 2e^{-t}-e^{-2t} (roots -1,-2) */
    CHECK("coupled linear",
          "First[x[2.0]/.NDSolve[{x'[t]==y[t],y'[t]==-2 x[t]-3 y[t],x[0]==1,y[0]==0},{x,y},{t,0,4}]]",
          2.0*exp(-2.0)-exp(-4.0), 1e-6);

    /* Lotka-Volterra: no closed form, but H = x - ln x + y - ln y is conserved.
     * Check H(4) == H(0) after integrating the nonlinear system. */
    run("lv = NDSolve[{x'[t]==x[t]-x[t] y[t],y'[t]==-y[t]+x[t] y[t],x[0]==1,y[0]==2},"
        "{x,y},{t,0,4}];");
    run("Hf[a_,b_]:=a-Log[a]+b-Log[b];");
    double H0 = 1.0 - log(1.0) + 2.0 - log(2.0);
    CHECK("Lotka-Volterra invariant H(t)",
          "Hf[First[x[4.0]/.lv], First[y[4.0]/.lv]]", H0, 1e-3);
}

/* ================= initial point not at a boundary / backward integration === */
static void test_backward_and_interior_ic(void) {
    /* y'=y, y(0)=1 -> e^x on [-2, 2]: query both sides of x0 */
    const char* B = "NDSolve[{y'[x]==y[x],y[0]==1},y,{x,-2,2}]";
    char buf[256];
    snprintf(buf, sizeof buf, "First[y[1.5]/.%s]", B);  CHECK("forward  e^{1.5}",  buf, exp(1.5), 1e-4);
    snprintf(buf, sizeof buf, "First[y[-1.5]/.%s]", B); CHECK("backward e^{-1.5}", buf, exp(-1.5), 1e-6);

    /* IC given in the interior: y(1)=E so y=e^x, range {x,0,2} */
    const char* M = "NDSolve[{y'[x]==y[x],y[1]==E},y,{x,0,2}]";
    snprintf(buf, sizeof buf, "First[y[0.0]/.%s]", M); CHECK("interior IC, y[0]=1", buf, 1.0, 1e-5);
    snprintf(buf, sizeof buf, "First[y[2.0]/.%s]", M); CHECK("interior IC, y[2]=e^2", buf, exp(2.0), 1e-4);
}

/* ================= method-of-lines PDEs ================= */

/* Build the heat-equation MOL ODE system on N interior points of [0,1] with
 * Dirichlet BCs and initial data sin(pi x); store into `hsol`.  The mode
 * sin(pi x_i) is an exact eigenvector of the discrete Laplacian, so the
 * semi-discrete solution is u_i(t) = e^{lambda t} sin(pi x_i) exactly, with
 * lambda = -(2/h^2)(1 - cos(pi h)).  Returns lambda. */
static double build_heat(char* buf, size_t bufsz, int N, double T) {
    double h = 1.0 / (N + 1), c = (double)(N + 1) * (N + 1);   /* c = 1/h^2 */
    int p = snprintf(buf, bufsz, "hsol = NDSolve[{");
    for (int i = 1; i <= N; i++) {
        char L[24], R[24];
        if (i == 1) strcpy(L, "0"); else snprintf(L, sizeof L, "u%d[t]", i - 1);
        if (i == N) strcpy(R, "0"); else snprintf(R, sizeof R, "u%d[t]", i + 1);
        p += snprintf(buf + p, bufsz - p, "u%d'[t]==%.17g*(%s-2*u%d[t]+%s), ", i, c, L, i, R);
    }
    for (int i = 1; i <= N; i++)
        p += snprintf(buf + p, bufsz - p, "u%d[0]==%.17g%s", i, sin(PI * i * h), (i == N ? "" : ", "));
    p += snprintf(buf + p, bufsz - p, "}, {");
    for (int i = 1; i <= N; i++) p += snprintf(buf + p, bufsz - p, "u%d%s", i, (i == N ? "" : ", "));
    p += snprintf(buf + p, bufsz - p, "}, {t, 0, %.17g}, Method->\"BDF\", MaxSteps->200000];", T);
    return -2.0 * c * (1.0 - cos(PI * h));
}

/* Wave equation u_tt = u_xx MOL; the mode gives u_i(t)=cos(omega t) sin(pi x_i),
 * omega = sqrt((2/h^2)(1-cos(pi h))).  Returns omega. */
static double build_wave(char* buf, size_t bufsz, int N, double T) {
    double h = 1.0 / (N + 1), c = (double)(N + 1) * (N + 1);
    int p = snprintf(buf, bufsz, "wsol = NDSolve[{");
    for (int i = 1; i <= N; i++) {
        char L[24], R[24];
        if (i == 1) strcpy(L, "0"); else snprintf(L, sizeof L, "u%d[t]", i - 1);
        if (i == N) strcpy(R, "0"); else snprintf(R, sizeof R, "u%d[t]", i + 1);
        p += snprintf(buf + p, bufsz - p, "u%d''[t]==%.17g*(%s-2*u%d[t]+%s), ", i, c, L, i, R);
    }
    for (int i = 1; i <= N; i++)
        p += snprintf(buf + p, bufsz - p, "u%d[0]==%.17g, u%d'[0]==0, ", i, sin(PI * i * h), i);
    /* strip trailing ", " and close */
    p -= 2;
    p += snprintf(buf + p, bufsz - p, "}, {");
    for (int i = 1; i <= N; i++) p += snprintf(buf + p, bufsz - p, "u%d%s", i, (i == N ? "" : ", "));
    p += snprintf(buf + p, bufsz - p, "}, {t, 0, %.17g}, MaxSteps->200000];", T);
    return sqrt(2.0 * c * (1.0 - cos(PI * h)));
}

static void test_pde_method_of_lines(void) {
    enum { N = 6 };
    double h = 1.0 / (N + 1);
    char buf[4096], q[64];

    /* Heat equation, BDF (stiff): compare each grid value to the exact
     * semi-discrete eigenmode e^{lambda T} sin(pi x_i). */
    double T = 0.05;
    double lambda = build_heat(buf, sizeof buf, N, T);
    run(buf);
    for (int i = 1; i <= N; i++) {
        double exact = exp(lambda * T) * sin(PI * i * h);
        snprintf(q, sizeof q, "First[u%d[%.4f]/.hsol]", i, T);
        char lbl[48]; snprintf(lbl, sizeof lbl, "heat MOL u%d(T)", i);
        CHECK(lbl, q, exact, 1e-4);
    }

    /* Wave equation, DOPRI5 (default): compare to cos(omega T) sin(pi x_i). */
    double Tw = 0.5;
    double omega = build_wave(buf, sizeof buf, N, Tw);
    run(buf);
    for (int i = 1; i <= N; i++) {
        double exact = cos(omega * Tw) * sin(PI * i * h);
        snprintf(q, sizeof q, "First[u%d[%.4f]/.wsol]", i, Tw);
        char lbl[48]; snprintf(lbl, sizeof lbl, "wave MOL u%d(T)", i);
        CHECK(lbl, q, exact, 1e-4);
    }
}

/* result head is exactly List (i.e. NDSolve solved, not $Failed/unevaluated) */
static bool result_is_list(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    bool ok = r && r->type == EXPR_FUNCTION
              && r->data.function.head->type == EXPR_SYMBOL
              && strcmp(r->data.function.head->data.symbol.name, "List") == 0
              && r->data.function.arg_count >= 1;
    expr_free(e); expr_free(r);
    return ok;
}
static void check_list(const char* label, const char* input) {
    if (result_is_list(input)) printf("ok:   %-42s = solved (List)\n", label);
    else { printf("FAIL: %s -> did not solve to a List [%s]\n", label, input); failures++; }
}

/* check a >= lo (used for step-count bounds) */
static void check_ge(const char* label, const char* input, double lo) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric\n", label); failures++; return; }
    if (!(v >= lo)) { printf("FAIL: %s -> %.6g (expected >= %.6g)\n", label, v, lo); failures++; }
    else printf("ok:   %-42s = %.6g (>= %.6g)\n", label, v, lo);
}

/* Exact solution of the discontinuous-corner heat problem
 *   u_t = u_xx,  u(0,x)=0,  u(t,0)=0,  u(t,1)=1
 * u(t,x) = x + (2/pi) sum_{n>=1} ((-1)^n/n) sin(n pi x) e^{-n^2 pi^2 t}. */
static double corner_heat_ref(double t, double x) {
    double s = x;
    for (int n = 1; n <= 500; n++) {
        double term = (((n & 1) ? -1.0 : 1.0) / n) * sin(n * PI * x)
                      * exp(-(double)n * n * PI * PI * t);
        s += (2.0 / PI) * term;
        if (fabs(term) < 1e-18) break;
    }
    return s;
}

/* ================= adaptive implicit (stiff) stepping ================= *
 *
 * BDF is now an adaptive variable-step method (orders 1-2) with predictor-
 * corrector local error control and Newton-failure step recovery.  These tests
 * exercise the stiff regimes where the old fixed-step BDF either lost accuracy
 * or (at incompatible IC/BC corners) diverged in Newton and returned nothing. */
static void test_adaptive_implicit(void) {
    /* Stiff linear scalar: y' = -1000(y - cos t) - sin t, y0=1 -> cos t.
     * The fast -1000 mode makes this stiff; adaptive BDF nails it. */
    CHECK("stiff BDF y'=-1000(y-cos)-sin",
          "First[y[1.0]/.NDSolve[{y'[t]==-1000 (y[t]-Cos[t])-Sin[t],y[0]==1},y,{t,0,1},"
          "Method->\"BDF\"]]", cos(1.0), 1e-7);

    /* Same via the StiffnessSwitching alias (resolves to BDF). */
    CHECK("stiff StiffnessSwitching -> cos",
          "First[y[1.0]/.NDSolve[{y'[t]==-1000 (y[t]-Cos[t])-Sin[t],y[0]==1},y,{t,0,1},"
          "Method->\"StiffnessSwitching\"]]", cos(1.0), 1e-7);

    /* Stiff manufactured variable-coefficient problem: y' = -50(y - t^2) + 2 t,
     * y0=0 -> y = t^2 exactly (transient e^{-50t} has zero coefficient).
     * Exercises the symbolic RHS sampler under stiff BDF. */
    CHECK("stiff BDF y'=-50(y-t^2)+2t -> t^2",
          "First[y[2.0]/.NDSolve[{y'[t]==-50 (y[t]-t^2)+2 t,y[0]==0},y,{t,0,2},"
          "Method->\"BDF\"]]", 4.0, 1e-5);

    /* Stiff linear system, well-separated eigenvalues -1000 and -1:
     *   x'=-1000 x, y'=-y  ->  x=e^{-1000 t}, y=e^{-t}.
     * The slow (accurately resolved) component is checked tightly; the fast
     * component has decayed below the absolute tolerance, so only that it is
     * ~0 is meaningful. */
    run("stf=NDSolve[{x'[t]==-1000 x[t],y'[t]==-y[t],x[0]==1,y[0]==1},{x,y},"
        "{t,0,0.02},Method->\"BDF\"];");
    CHECK("stiff system slow comp y=e^{-t}", "First[y[0.02]/.stf]", exp(-0.02), 1e-6);
    check_le("stiff system fast comp ~0",   "Abs[First[x[0.02]/.stf]]", 1e-6);

    /* ImplicitTrapezoid (order 2, one-step) on the stiff forced problem. */
    CHECK("ImplicitTrapezoid stiff -> cos",
          "First[y[1.0]/.NDSolve[{y'[t]==-50 (y[t]-Cos[t])+(-Sin[t]),y[0]==1},y,{t,0,1},"
          "Method->\"ImplicitTrapezoid\"]]", cos(1.0), 1e-4);

    /* BackwardEuler (order 1) on a gentle short stiff problem where its step
     * budget suffices: y'=-20(y-1), y0=0 -> 1 - e^{-20 t}. */
    CHECK("BackwardEuler gentle stiff",
          "First[y[0.4]/.NDSolve[{y'[t]==-20 (y[t]-1),y[0]==0},y,{t,0,0.5},"
          "Method->\"BackwardEuler\",PrecisionGoal->6]]", 1.0 - exp(-8.0), 1e-4);

    /* Adaptivity: a looser PrecisionGoal must take strictly fewer steps than a
     * tight one (and both must stay well within the MaxSteps budget). */
    run("nlo=Length[NDSolve[{y'[t]==-100 (y[t]-1),y[0]==0},y,{t,0,5},Method->\"BDF\","
        "PrecisionGoal->4][[1,1,2,2]]];");
    run("nhi=Length[NDSolve[{y'[t]==-100 (y[t]-1),y[0]==0},y,{t,0,5},Method->\"BDF\","
        "PrecisionGoal->9][[1,1,2,2]]];");
    check_ge("BDF adaptivity: tighter goal -> more steps", "N[nhi-nlo]", 1.0);
    check_le("BDF loose goal well under budget",           "N[nlo]", 9000.0);

    /* Tight vs loose accuracy under BDF. */
    check_le("BDF tight-goal error",
        "Abs[First[y[3.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,5},Method->\"BDF\","
        "PrecisionGoal->9]]-Exp[-3.0]]", 1e-6);
}

/* ================= higher-order (variable-order) BDF ================= *
 *
 * BDF is a variable-step VARIABLE-ORDER method (orders 1-5).  High order shows
 * up as high accuracy at few steps on smooth stiff problems, and as accuracy
 * that tracks the tolerance goals — an order-2-only method would need ~50x more
 * steps for the tight-tolerance cases below. */
static void test_bdf_high_order(void) {
    /* Accuracy tracks the goals: PG/AG 10 on stiff decay -> ~1e-9 or better. */
    CHECK("BDF hi-order accuracy y'=-y",
          "First[y[10.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,10},Method->\"BDF\","
          "PrecisionGoal->10,AccuracyGoal->10]]", exp(-10.0), 1e-8);

    /* Efficiency: tight tolerance is reached in few steps (variable order).  An
     * order-2 method needs O((1/rtol)^(1/3)) ~ 5000 steps for rtol 1e-10; the
     * variable-order method uses far fewer.  Bound generously but well below. */
    check_le("BDF hi-order step count (tight goal)",
        "N[Length[NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,10},Method->\"BDF\","
        "PrecisionGoal->10,AccuracyGoal->10][[1,1,2,2]]]]", 1500.0);

    /* Tighter goals give a strictly smaller error (order/step adapt to tol). */
    check_le("BDF loose-goal error",
        "Abs[First[y[8.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,10},Method->\"BDF\","
        "PrecisionGoal->4,AccuracyGoal->6]]-Exp[-8.0]]", 1e-4);
    check_le("BDF tight-goal error (few digits more)",
        "Abs[First[y[8.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,10},Method->\"BDF\","
        "PrecisionGoal->11,AccuracyGoal->13]]-Exp[-8.0]]", 1e-9);

    /* Smooth non-stiff ODE at high order: sin over several periods. */
    CHECK("BDF hi-order harmonic x=cos",
          "First[x[12.0]/.NDSolve[{x'[t]==y[t],y'[t]==-x[t],x[0]==1,y[0]==0},{x,y},"
          "{t,0,12},Method->\"BDF\",PrecisionGoal->10,AccuracyGoal->10]]", cos(12.0), 1e-6);

    /* Stiff linear system (eigenvalues -500, -1) at high order + tight goal. */
    run("hs=NDSolve[{x'[t]==-500 x[t]+499 y[t],y'[t]==-y[t],x[0]==2,y[0]==1},{x,y},"
        "{t,0,3},Method->\"BDF\",PrecisionGoal->9,AccuracyGoal->9];");
    /* y = e^{-t}; x = e^{-t} (particular) since x'=-500x+499 e^{-t}, x(0)=2:
     * homogeneous e^{-500t} (coeff 1) + particular e^{-t}: x=e^{-500t}+e^{-t}. */
    CHECK("BDF stiff-system slow part", "First[y[3.0]/.hs]", exp(-3.0), 1e-7);
    CHECK("BDF stiff-system x=e^{-500t}+e^{-t}", "First[x[3.0]/.hs]",
          exp(-500.0*3.0) + exp(-3.0), 1e-7);
}

/* ================= incompatible IC/BC corner PDE (the headline fix) ======= *
 *
 * u_t = u_xx with u(0,x)=0 but u(t,1)=1 has a discontinuous corner at (0,1).
 * Fixed-step BDF Newton diverged here (ndcf) and NDSolve returned nothing; the
 * adaptive BDF (Newton-failure step recovery) now solves it.  We compare the
 * method-of-lines solution to the exact PDE Fourier series (loose tolerance,
 * since the comparison includes spatial-discretization error). */
static void test_corner_pde(void) {
    const char* SOL =
        "csol=NDSolve[{D[u[t,x],t]==D[u[t,x],x,x],u[0,x]==0,u[t,0]==0,u[t,1]==1},"
        "u,{t,0,0.2},{x,0,1}];";
    run(SOL);
    /* headline: it solved at all (was $Failed/unevaluated before). */
    check_list("corner heat solves (was ndcf)", "csol");

    double T = 0.2;
    CHECK("corner heat u(T,0.5)", "First[u[0.2,0.5]/.csol]", corner_heat_ref(T, 0.5), 5e-3);
    CHECK("corner heat u(T,0.9)", "First[u[0.2,0.9]/.csol]", corner_heat_ref(T, 0.9), 5e-3);
    CHECK("corner heat u(T,0.25)", "First[u[0.2,0.25]/.csol]", corner_heat_ref(T, 0.25), 5e-3);

    /* Boundary values are honoured exactly. */
    CHECK("corner heat BC u(T,1)=1", "First[u[0.2,1.0]/.csol]", 1.0, 1e-6);
    CHECK("corner heat BC u(T,0)=0", "First[u[0.2,0.0]/.csol]", 0.0, 1e-6);
}

/* ================= adaptive Adams (non-stiff multistep) ================= */
static void test_adaptive_adams(void) {
    /* accuracy on decay */
    CHECK("Adams decay y'=-y",
          "First[y[3.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,5},Method->\"Adams\"]]",
          exp(-3.0), 1e-5);

    /* adaptive step responds to the goal: tighter goal is at least as accurate */
    check_le("Adams loose-goal error",
        "Abs[First[y[5.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,5},Method->\"Adams\","
        "PrecisionGoal->4]]-Exp[-5.0]]", 1e-4);
    check_le("Adams tight-goal error",
        "Abs[First[y[5.0]/.NDSolve[{y'[t]==-y[t],y[0]==1},y,{t,0,5},Method->\"Adams\","
        "PrecisionGoal->9]]-Exp[-5.0]]", 1e-5);

    /* nonlinear (logistic) under adaptive Adams */
    CHECK("Adams logistic",
          "First[y[2.0]/.NDSolve[{y'[t]==y[t](1-y[t]),y[0]==1/2},y,{t,0,5},Method->\"Adams\"]]",
          1.0/(1.0+exp(-2.0)), 1e-5);

    /* oscillatory system (harmonic): adaptive Adams tracks cos over many periods */
    CHECK("Adams harmonic x=cos",
          "First[x[5.0]/.NDSolve[{x'[t]==y[t],y'[t]==-x[t],x[0]==1,y[0]==0},{x,y},"
          "{t,0,8},Method->\"Adams\"]]", cos(5.0), 1e-3);
}

/* ================= accuracy/efficiency parity + precision ================= */
static void test_parity_and_precision(void) {
    /* ode45-parity: default (adaptive DOPRI5) accuracy on a full period of the
     * harmonic oscillator is well within machine tolerance. */
    CHECK("ode45-parity accuracy",
          "First[y[10.0]/.NDSolve[{y''[x]+y[x]==0,y[0]==1,y'[0]==0},y,{x,0,10}]]",
          cos(10.0), 1e-5);

    /* FSAL efficiency: ~6 function evals per accepted step (ode45 parity). */
    run("cc=0; es=NDSolve[{y'[x]==-y[x],y[0]==1},y,{x,0,10},EvaluationMonitor:>(cc=cc+1)];");
    check_le("evals per step (FSAL ~6)", "N[cc/Length[es[[1,1,2,2]]]]", 6.6);

    /* tighter PrecisionGoal reduces error */
    check_le("loose goal error",
        "Abs[First[y[5.0]/.NDSolve[{y'[x]==-y[x],y[0]==1},y,{x,0,5},PrecisionGoal->4]]-Exp[-5.0]]", 1e-3);
    check_le("tight goal error",
        "Abs[First[y[5.0]/.NDSolve[{y'[x]==-y[x],y[0]==1},y,{x,0,5},PrecisionGoal->11]]-Exp[-5.0]]", 1e-8);

#ifdef USE_MPFR
    /* arbitrary precision beyond machine: >17 correct digits at 30-digit WP */
    check_le("MPFR error (30-digit WP)",
        "Abs[First[y[N[1,30]]/.NDSolve[{y'[x]==-y[x],y[0]==1},y,{x,0,1},"
        "WorkingPrecision->30,PrecisionGoal->22,AccuracyGoal->22,MaxSteps->500000]]-N[Exp[-1],30]]",
        1e-17);
#endif
}

int main(void) {
    mute_stderr_once();
    symtab_init();
    core_init();

    TEST(test_first_order);
    TEST(test_second_order);
    TEST(test_higher_and_varcoef);
    TEST(test_systems);
    TEST(test_backward_and_interior_ic);
    TEST(test_pde_method_of_lines);
    TEST(test_adaptive_implicit);
    TEST(test_bdf_high_order);
    TEST(test_corner_pde);
    TEST(test_adaptive_adams);
    TEST(test_parity_and_precision);

    if (failures == 0) printf("\nAll classical NDSolve tests passed.\n");
    else printf("\n%d classical NDSolve checks FAILED.\n", failures);
    return 0;
}
