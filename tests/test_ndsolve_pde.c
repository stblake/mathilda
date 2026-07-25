/* NDSolve on partial differential equations by the method of lines (Phase 1).
 *
 * A PDE is discretized in space (2nd-order central differences) into a large
 * ODE system that the shared time integrator solves.  Two kinds of reference
 * are used:
 *
 *   (1) Semi-discrete eigenmode.  For linear constant-coefficient problems the
 *       initial profile sin(pi x) is an EXACT eigenvector of the discrete
 *       Laplacian, so the semi-discrete solution is U_i(t) = e^{lambda t}
 *       sin(pi x_i) (heat) or cos(omega t) sin(pi x_i) (wave) *exactly* — with
 *       lambda = -(2/h^2)(1-cos(pi h)), omega = sqrt(-lambda).  Comparing to this
 *       isolates the time-integration error from spatial-discretization error.
 *
 *   (2) Manufactured solution.  For nonlinear / forced problems we pick a target
 *       U(t,x), add the forcing that makes it exact, and check the solver
 *       converges to U at the O(h^2) rate of the spatial discretization.
 *
 * Soft asserts: prints FAIL and keeps going. Run: ./ndsolve_pde_tests */
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

static void run(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    expr_free(e); expr_free(r);
}

/* head symbol of the evaluated expression, copied into `out` (empty if none) */
static void eval_head(const char* input, char* out, size_t n) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    out[0] = '\0';
    if (r) {
        const char* nm = NULL;
        if (r->type == EXPR_SYMBOL) nm = r->data.symbol.name;
        else if (r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL)
            nm = r->data.function.head->data.symbol.name;
        if (nm) { strncpy(out, nm, n - 1); out[n - 1] = '\0'; }
    }
    expr_free(e); expr_free(r);
}

static void check_close(const char* label, const char* input, double expected, double tol) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric [%s]\n", label, input); failures++; return; }
    if (fabs(v - expected) > tol) {
        printf("FAIL: %s -> %.12g (expected %.12g, |err|=%.2e, tol %.1e)\n",
               label, v, expected, fabs(v - expected), tol);
        failures++;
    } else {
        printf("ok:   %-46s = %.10g  (|err|=%.1e)\n", label, v, fabs(v - expected));
    }
}
#define CHECK(label, input, expected, tol) check_close(label, input, expected, tol)

static void check_true(const char* label, bool cond, const char* detail) {
    if (cond) printf("ok:   %-46s [%s]\n", label, detail);
    else { printf("FAIL: %s [%s]\n", label, detail); failures++; }
}

/* semi-discrete discrete-Laplacian eigenvalue for the sin(pi x) mode */
static double disc_lambda(int nx) {
    double h = 1.0 / (nx - 1);
    return -2.0 * (1.0 - cos(PI * h)) / (h * h);
}

/* ============================================================= *
 *  1. Heat equation  u_t = u_xx  (Dirichlet, BDF)               *
 *     vs exact semi-discrete eigenmode e^{lambda t} sin(pi x)   *
 * ============================================================= */
static void test_heat_eigenmode(void) {
    const int nx = 11;                 /* 9 interior unknowns */
    const double T = 0.05, h = 1.0 / (nx - 1), lam = disc_lambda(nx);
    char buf[1024], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "hs = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, MaxSteps->3000];", T, nx);
    run(buf);
    for (int i = 1; i <= nx - 2; i++) {
        double xi = i * h, exact = exp(lam * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.6f, %.10f] /. hs]", T, xi);
        snprintf(lbl, sizeof lbl, "heat eigenmode u(T,x%d)", i);
        CHECK(lbl, q, exact, 1e-4);
    }
}

/* ============================================================= *
 *  2. Spatial convergence: as the grid refines, the numerical  *
 *     solution approaches the TRUE PDE mode e^{-pi^2 t} sin.    *
 * ============================================================= */
static void test_heat_spatial_convergence(void) {
    const double T = 0.03, xq = 0.5, exact = exp(-PI * PI * T) * sin(PI * xq);
    int grids[3] = { 11, 21, 41 };
    double err[3];
    for (int g = 0; g < 3; g++) {
        int nx = grids[g];
        char buf[1024], q[96];
        snprintf(buf, sizeof buf,
            "cs = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
            "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
            "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
            "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, MaxSteps->3000];", T, nx);
        run(buf);
        snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. cs]", T, xq);
        double v; err[g] = eval_double(q, &v) ? fabs(v - exact) : 1e9;
        printf("ok:   heat conv nx=%-3d err=%.3e\n", nx, err[g]);
    }
    /* 2nd-order: halving h (11->21->41) quarters the error, roughly. */
    check_true("heat spatial 2nd-order (nx 11->21)", err[1] < err[0] * 0.5,
               "error at least halves on refinement");
    check_true("heat spatial 2nd-order (nx 21->41)", err[2] < err[1] * 0.5,
               "error at least halves on refinement");
    check_true("heat finest grid near true PDE", err[2] < 5e-4, "|err| < 5e-4 at nx=41");
}

/* ============================================================= *
 *  3. Wave equation  u_tt = u_xx  (Dirichlet, default DOPRI5)   *
 *     vs exact semi-discrete cos(omega t) sin(pi x)            *
 * ============================================================= */
static void test_wave_eigenmode(void) {
    const int nx = 21;
    const double T = 0.5, h = 1.0 / (nx - 1);
    const double omega = sqrt(-disc_lambda(nx));
    char buf[1024], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "ws = NDSolve[{D[u[t,x],{t,2}]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "Derivative[1,0][u][0,x]==0, u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, "
        "{x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, MaxSteps->100000];", T, nx);
    run(buf);
    for (int i = 4; i <= nx - 2; i += 5) {
        double xi = i * h, exact = cos(omega * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.6f, %.10f] /. ws]", T, xi);
        snprintf(lbl, sizeof lbl, "wave eigenmode u(T,x%d)", i);
        CHECK(lbl, q, exact, 5e-4);
    }
}

/* ============================================================= *
 *  4. Reaction-diffusion  u_t = u_xx + a u  (Dirichlet, BDF)    *
 *     eigenmode e^{(lambda+a) t} sin(pi x)                      *
 * ============================================================= */
static void test_reaction_diffusion(void) {
    const int nx = 11;
    const double T = 0.08, a = 3.0, h = 1.0 / (nx - 1);
    const double rate = disc_lambda(nx) + a;
    char buf[1024], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "rs = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}]+%.1f u[t,x], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, Method->\"BDF\", MaxSteps->3000];",
        a, T, nx);
    run(buf);
    for (int i = 2; i <= nx - 2; i += 3) {
        double xi = i * h, exact = exp(rate * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.6f, %.10f] /. rs]", T, xi);
        snprintf(lbl, sizeof lbl, "reaction-diffusion u(T,x%d)", i);
        CHECK(lbl, q, exact, 2e-4);
    }
}

/* ============================================================= *
 *  5. Inhomogeneous Dirichlet BC -> linear steady state u=x    *
 *     Smooth compatible data u(0,x)=x+sin(pi x) relaxes to x:  *
 *     the discrete Laplacian is exact on the linear profile,   *
 *     and the sin transient decays, so u(T,x_i) -> x_i.        *
 * ============================================================= */
static void test_inhomogeneous_bc(void) {
    const int nx = 11;
    const double T = 3.0;             /* long enough to reach steady state */
    char buf[1024], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "is = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==x + Sin[Pi x], "
        "u[t,0]==0, u[t,1]==1}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, Method->\"BDF\", MaxSteps->3000];",
        T, nx);
    run(buf);
    double h = 1.0 / (nx - 1);
    for (int i = 2; i <= nx - 2; i += 3) {
        double xi = i * h;
        snprintf(q, sizeof q, "First[u[%.4f, %.10f] /. is]", T, xi);
        snprintf(lbl, sizeof lbl, "steady-state u(inf,x%d)=x", i);
        CHECK(lbl, q, xi, 5e-3);
    }
    /* the right boundary value must equal the Dirichlet datum exactly */
    snprintf(q, sizeof q, "First[u[%.4f, 1.0] /. is]", T);
    CHECK("right Dirichlet boundary held", q, 1.0, 1e-9);
}

/* ============================================================= *
 *  6. Time-dependent Dirichlet BC: boundary value is honoured  *
 * ============================================================= */
static void test_time_dependent_bc(void) {
    const int nx = 11;
    char buf[1024];
    snprintf(buf, sizeof buf,
        "ts = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==0, "
        "u[t,0]==Sin[t], u[t,1]==0}, u, {t,0,1.0}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, Method->\"BDF\", MaxSteps->3000];", nx);
    run(buf);
    CHECK("time-dependent left BC u(0.7,0)=Sin[0.7]",
          "First[u[0.7, 0.0] /. ts]", sin(0.7), 1e-6);
}

/* ============================================================= *
 *  7. Nonlinear manufactured solution (stress, symbolic sampler)*
 *     u_t = u_xx - u^3 + S,  U(t,x)=e^{-t} sin(pi x),           *
 *     S = (pi^2-1) e^{-t} sin(pi x) + (e^{-t} sin(pi x))^3      *
 *     -> O(h^2) spatial convergence to U.                       *
 * ============================================================= */
static void nonlin_solve(int nx, double T) {
    char buf[1400];
    snprintf(buf, sizeof buf,
        "ns = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}] - u[t,x]^3 "
        "+ (Pi^2-1) Exp[-t] Sin[Pi x] + (Exp[-t] Sin[Pi x])^3, "
        "u[0,x]==Sin[Pi x], u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, Method->\"BDF\", MaxSteps->2000];",
        T, nx);
    run(buf);
}

static void test_nonlinear_manufactured(void) {
    const double T = 0.3, xq = 0.5, exact = exp(-T) * sin(PI * xq);
    char q[96];
    nonlin_solve(11, T);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. ns]", T, xq);
    double v11; double e11 = eval_double(q, &v11) ? fabs(v11 - exact) : 1e9;
    CHECK("nonlinear manufactured u(T,0.5) nx=11", q, exact, 7e-3);

    nonlin_solve(21, T);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. ns]", T, xq);
    double v21; double e21 = eval_double(q, &v21) ? fabs(v21 - exact) : 1e9;
    printf("ok:   nonlinear conv nx=11 err=%.2e  nx=21 err=%.2e  ratio=%.2f\n",
           e11, e21, e21 > 0 ? e11 / e21 : 0.0);
    check_true("nonlinear O(h^2) convergence", e21 < e11 * 0.5,
               "refining the grid at least halves the error");
}

/* ============================================================= *
 *  8. Structure / return shape / method equivalence            *
 * ============================================================= */
static void test_structure_and_equivalence(void) {
    /* return shape: {{u -> InterpolatingFunction[...]}} */
    run("sh = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,0.02}, {x,0,1}, Method->\"BDF\", MaxSteps->2000];");
    char hd[64];
    eval_head("sh", hd, sizeof hd);
    check_true("PDE result is a List", strcmp(hd, "List") == 0, hd[0] ? hd : "(null)");
    eval_head("First[First[sh]]", hd, sizeof hd);        /* the Rule */
    check_true("solution rule present", strcmp(hd, "Rule") == 0, hd[0] ? hd : "(null)");
    eval_head("Last[First[First[sh]]]", hd, sizeof hd);  /* rhs: InterpolatingFunction */
    check_true("rhs is InterpolatingFunction",
               strcmp(hd, "InterpolatingFunction") == 0, hd[0] ? hd : "(null)");

    /* NDSolve`MethodOfLines[...] == NDSolve[..., Method->"MethodOfLines"] == auto */
    const char* eqn = "{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
                      "u[t,0]==0, u[t,1]==0}";
    char buf[1024];
    snprintf(buf, sizeof buf,
        "e1 = First[u[0.03,0.5] /. NDSolve[%s, u, {t,0,0.03}, {x,0,1}, "
        "Method->\"BDF\", MaxSteps->2000]];", eqn);
    run(buf);
    snprintf(buf, sizeof buf,
        "e2 = First[u[0.03,0.5] /. NDSolve`MethodOfLines[%s, u, {t,0,0.03}, {x,0,1}, "
        "Method->\"BDF\", MaxSteps->2000]];", eqn);
    run(buf);
    double d1, d2;
    bool ok = eval_double("e1", &d1) && eval_double("e2", &d2);
    check_true("NDSolve`MethodOfLines == NDSolve (Method->MoL)",
               ok && fabs(d1 - d2) < 1e-9, "backtick form matches option form");

    /* applied form NDSolve[eqns, u[t,x], ...] -> {{u[t,x] -> IF[t,x]}} */
    run("ap = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u[t,x], {t,0,0.02}, {x,0,1}, Method->\"BDF\", MaxSteps->2000];");
    eval_head("ap", hd, sizeof hd);
    check_true("applied-form result is a List", strcmp(hd, "List") == 0, hd[0] ? hd : "(null)");
}

/* ============================================================= *
 *  9. Interpolation off grid nodes + at interior time          *
 * ============================================================= */
static void test_offgrid_interpolation(void) {
    const int nx = 21;
    const double T = 0.04, h = 1.0 / (nx - 1), lam = disc_lambda(nx);
    char buf[1024];
    snprintf(buf, sizeof buf,
        "os = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, Method->\"BDF\", MaxSteps->3000];", T, nx);
    run(buf);
    /* midway between grid nodes: still ~e^{lam t} sin(pi x) to interpolation order */
    double xq = 2.5 * h;                               /* between nodes 2 and 3 */
    double exact = exp(lam * T) * sin(PI * xq);
    CHECK("off-grid x interpolation",
          "First[u[0.04, 0.125] /. os]", exact, 5e-4);
    /* interior time between accepted steps */
    double tq = 0.017, xn = 10 * h;
    double exact2 = exp(lam * tq) * sin(PI * xn);
    char q[96]; snprintf(q, sizeof q, "First[u[%.4f, %.10f] /. os]", tq, xn);
    CHECK("interior-time interpolation", q, exact2, 5e-4);
    (void)xq;
}

/* ============================================================= *
 *  10. Advection-diffusion  u_t = u_xx - c u_x  (exercises the  *
 *      first-order spatial stencil).  Manufactured exact        *
 *      U(t,x)=e^{-t} sin(pi x); S makes it exact. O(h^2).       *
 * ============================================================= */
static void adv_diff_solve(int nx, double T) {
    char buf[1400];
    snprintf(buf, sizeof buf,
        "as = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}] - 2 D[u[t,x],x] "
        "+ (Pi^2-1) Exp[-t] Sin[Pi x] + 2 Pi Exp[-t] Cos[Pi x], "
        "u[0,x]==Sin[Pi x], u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, Method->\"BDF\", MaxSteps->2000];",
        T, nx);
    run(buf);
}
static void test_advection_diffusion(void) {
    const double T = 0.3, xq = 0.5, exact = exp(-T) * sin(PI * xq);
    char q[96];
    adv_diff_solve(11, T);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. as]", T, xq);
    double v11; double e11 = eval_double(q, &v11) ? fabs(v11 - exact) : 1e9;
    CHECK("advection-diffusion u(T,0.5) nx=11", q, exact, 1e-2);
    adv_diff_solve(21, T);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. as]", T, xq);
    double v21; double e21 = eval_double(q, &v21) ? fabs(v21 - exact) : 1e9;
    printf("ok:   adv-diff conv nx=11 err=%.2e  nx=21 err=%.2e  ratio=%.2f\n",
           e11, e21, e21 > 0 ? e11 / e21 : 0.0);
    check_true("advection-diffusion first-deriv stencil O(h^2)", e21 < e11 * 0.6,
               "refining halves error (validates u_x stencil)");
}

/* ============================================================= *
 *  11. Stress: larger grids stay exact on the eigenmode.        *
 *      Wave (non-stiff, adaptive DOPRI5) on a fine grid.        *
 * ============================================================= */
static void test_large_grid_stress(void) {
    const int nx = 41;
    const double T = 0.25, h = 1.0 / (nx - 1);
    const double omega = sqrt(-disc_lambda(nx));
    char buf[1024], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "ls = NDSolve[{D[u[t,x],{t,2}]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "Derivative[1,0][u][0,x]==0, u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, "
        "{x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, MaxSteps->200000];", T, nx);
    run(buf);
    for (int i = 8; i <= nx - 2; i += 12) {
        double xi = i * h, exact = cos(omega * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.6f, %.10f] /. ls]", T, xi);
        snprintf(lbl, sizeof lbl, "large-grid wave u(T,x%d) nx=41", i);
        CHECK(lbl, q, exact, 1e-3);
    }
}

/* ============================================================= *
 *  12. DifferenceOrder: higher-order stencils converge faster.  *
 *      Heat vs true PDE e^{-pi^2 t} sin(pi x); error ~ h^q.      *
 * ============================================================= */
static double heat_order_err(int nx, int order, double T, double xq, double exact) {
    char buf[1024], q[96];
    snprintf(buf, sizeof buf,
        "os = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->%d}}, "
        "Method->\"BDF\", MaxSteps->3000];", T, nx, order);
    run(buf);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. os]", T, xq);
    double v;
    return eval_double(q, &v) ? fabs(v - exact) : 1e9;
}
static void test_difference_order(void) {
    const double T = 0.02, xq = 0.5, exact = exp(-PI * PI * T) * sin(PI * xq);
    double e2 = heat_order_err(21, 2, T, xq, exact);
    double e4 = heat_order_err(21, 4, T, xq, exact);
    double e6 = heat_order_err(21, 6, T, xq, exact);
    printf("ok:   DifferenceOrder err nx=21: o2=%.2e o4=%.2e o6=%.2e\n", e2, e4, e6);
    check_true("order 4 beats order 2", e4 < e2 * 0.2, "higher order -> smaller error");
    check_true("order 6 beats order 4", e6 < e4, "still smaller at order 6");
    /* 4th-order spatial convergence: halving h drops the error by ~2^4. */
    double c11 = heat_order_err(11, 4, T, xq, exact);
    double c21 = heat_order_err(21, 4, T, xq, exact);
    printf("ok:   order-4 conv nx=11 err=%.2e nx=21 err=%.2e ratio=%.1f\n",
           c11, c21, c21 > 0 ? c11 / c21 : 0.0);
    check_true("order-4 stencil converges >O(h^2)", c21 < c11 * 0.16,
               "ratio well above 4 (2nd-order); ~16 expected");
}

/* ============================================================= *
 *  13. Neumann BC: insulated heat  u_x(t,0)=u_x(t,1)=0.         *
 *      cos(pi x) satisfies homogeneous Neumann; exact solution  *
 *      e^{-pi^2 t} cos(pi x).  One-sided boundary elimination    *
 *      gives O(h^q) error -> check accuracy + convergence.       *
 * ============================================================= */
static double neumann_err(int nx, int order, double T, double xq, double exact) {
    char buf[1100], q[96];
    snprintf(buf, sizeof buf,
        "nm = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Cos[Pi x], "
        "Derivative[0,1][u][t,0]==0, Derivative[0,1][u][t,1]==0}, u, {t,0,%.4f}, "
        "{x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->%d}}, "
        "Method->\"BDF\", MaxSteps->3000];", T, nx, order);
    run(buf);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. nm]", T, xq);
    double v;
    return eval_double(q, &v) ? fabs(v - exact) : 1e9;
}
static void test_neumann_bc(void) {
    const double T = 0.03, xq = 0.25, exact = exp(-PI * PI * T) * cos(PI * xq);
    double e = neumann_err(41, 4, T, xq, exact);
    printf("ok:   Neumann heat u(T,0.25) err=%.2e (exact=%.6f)\n", e, exact);
    check_true("Neumann BC accurate", e < 5e-3, "|err| < 5e-3 at nx=41, order 4");
    double c21 = neumann_err(21, 4, T, xq, exact);
    double c41 = neumann_err(41, 4, T, xq, exact);
    printf("ok:   Neumann conv nx=21 err=%.2e nx=41 err=%.2e ratio=%.1f\n",
           c21, c41, c41 > 0 ? c21 / c41 : 0.0);
    check_true("Neumann converges on refinement", c41 < c21 * 0.6, "error shrinks");
}

/* ============================================================= *
 *  14. Robin BC: Dirichlet left, Robin right  u+u_x=2.          *
 *      Steady state is u=x (exact discrete fixed point); with    *
 *      IC u(0,x)=x the solution must stay x -> validates the     *
 *      Robin coefficient extraction + elimination.               *
 * ============================================================= */
static void test_robin_bc(void) {
    const int nx = 15;
    char buf[1100], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "rb = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==x, "
        "u[t,0]==0, u[t,1]+Derivative[0,1][u][t,1]==2}, u, {t,0,1.0}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "Method->\"BDF\", MaxSteps->3000];", nx);
    run(buf);
    double hh = 1.0 / (nx - 1);
    for (int i = 3; i <= nx - 2; i += 4) {
        double xi = i * hh;
        snprintf(q, sizeof q, "First[u[1.0, %.10f] /. rb]", xi);
        snprintf(lbl, sizeof lbl, "Robin steady u(1,x%d)=x", i);
        CHECK(lbl, q, xi, 2e-3);
    }
    /* the Robin edge value is held at the eliminated value (u(1)=1) */
    CHECK("Robin right value", "First[u[1.0, 1.0] /. rb]", 1.0, 5e-3);
}

/* ============================================================= *
 *  15. Periodic BC: heat with u(t,0)==u(t,1), IC sin(2 pi x).   *
 *      sin(2 pi x) is an exact eigenmode of the cyclic discrete  *
 *      Laplacian -> exact semi-discrete e^{lambda t} sin(2 pi x),*
 *      lambda = -(2/h^2)(1-cos(2 pi h)).                         *
 * ============================================================= */
static void test_periodic_bc(void) {
    const int nx = 21;
    const double T = 0.02, h = 1.0 / (nx - 1);
    const double lam = -2.0 * (1.0 - cos(2.0 * PI * h)) / (h * h);
    char buf[1100], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "pd = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[2 Pi x], "
        "u[t,0]==u[t,1]}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "Method->\"BDF\", MaxSteps->3000];", T, nx);
    run(buf);
    for (int i = 2; i <= nx - 2; i += 4) {
        double xi = i * h, exact = exp(lam * T) * sin(2.0 * PI * xi);
        snprintf(q, sizeof q, "First[u[%.4f, %.10f] /. pd]", T, xi);
        snprintf(lbl, sizeof lbl, "periodic heat u(T,x%d)", i);
        CHECK(lbl, q, exact, 1e-4);
    }
    /* periodicity: the two ends carry the same value */
    double v0, v1;
    bool ok = eval_double("First[u[0.01, 0.0] /. pd]", &v0) &&
              eval_double("First[u[0.01, 1.0] /. pd]", &v1);
    check_true("periodic endpoints match", ok && fabs(v0 - v1) < 1e-9, "u(t,0)==u(t,1)");
}

/* ============================================================= *
 *  16. Compiled operator: fast path == symbolic path.          *
 *      Same linear PDE with Compiled->True vs Compiled->False    *
 *      must agree; and the auto-selected method (no Method given) *
 *      solves stiff heat correctly.                              *
 * ============================================================= */
static double advdiff_compiled(int nx, int comp, double T, double xq) {
    char buf[1300], q[96];
    snprintf(buf, sizeof buf,
        "cp = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}] - D[u[t,x],x], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->4}}, "
        "Method->\"BDF\", Compiled->%s, MaxSteps->3000];", T, nx, comp ? "True" : "False");
    run(buf);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. cp]", T, xq);
    double v;
    return eval_double(q, &v) ? v : NAN;
}
static void test_compiled_operator(void) {
    double vt = advdiff_compiled(21, 1, 0.2, 0.5);   /* compiled operator */
    double vf = advdiff_compiled(21, 0, 0.2, 0.5);   /* symbolic sampler  */
    printf("ok:   compiled=%.12g  symbolic=%.12g  diff=%.2e\n", vt, vf, fabs(vt - vf));
    check_true("compiled operator == symbolic RHS", fabs(vt - vf) < 1e-9,
               "fast path agrees with symbolic path");

    /* auto-method: parabolic heat with NO time-integration method specified
     * (only the MoL controller) must auto-select BDF and stay accurate. */
    const int nx = 11;
    const double T = 0.05, h = 1.0 / (nx - 1);
    const double lam = disc_lambda(nx);
    char buf[1024], q[96];
    snprintf(buf, sizeof buf,
        "au = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "MaxSteps->5000];", T, nx);
    run(buf);
    snprintf(q, sizeof q, "First[u[%.4f, %.10f] /. au]", T, 5 * h);
    CHECK("auto-selected method (stiff heat)", q, exp(lam * T) * sin(PI * 5 * h), 1e-4);
}

/* ============================================================= *
 *  17. Efficiency/scale: a larger stiff grid solves quickly and *
 *      accurately with the compiled banded operator + BDF.      *
 *      (Without it this would be O(d^3) dense + symbolic Jacobian *
 *      per step.)                                               *
 * ============================================================= */
static void test_operator_scale(void) {
    const int nx = 81;
    const double T = 0.02, h = 1.0 / (nx - 1), lam = disc_lambda(nx);
    char buf[1024], q[96], lbl[64];
    snprintf(buf, sizeof buf,
        "sc = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "Method->\"BDF\", MaxSteps->4000];", T, nx);
    run(buf);
    for (int i = 20; i <= nx - 2; i += 20) {
        double xi = i * h, exact = exp(lam * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.4f, %.10f] /. sc]", T, xi);
        snprintf(lbl, sizeof lbl, "scale nx=81 u(T,x%d)", i);
        CHECK(lbl, q, exact, 1e-4);
    }
}

/* ============================================================= *
 *  18. Two spatial dimensions: heat u_t = u_xx + u_yy and wave  *
 *      u_tt = u_xx + u_yy on the unit square, Dirichlet 0.       *
 *      sin(pi x) sin(pi y) is an exact eigenmode of the discrete *
 *      2-D Laplacian -> exact semi-discrete solution.           *
 * ============================================================= */
static void test_pde_2d_heat(void) {
    const int n = 13;
    const double h = 1.0 / (n - 1);
    const double lam = 2.0 * (-2.0 * (1.0 - cos(PI * h)) / (h * h));  /* lx+ly */
    const double T = 0.015;
    char buf[1400], q[128], lbl[64];
    snprintf(buf, sizeof buf,
        "h2 = NDSolve[{D[u[t,x,y],t]==D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}], "
        "u[0,x,y]==Sin[Pi x] Sin[Pi y], u[t,0,y]==0, u[t,1,y]==0, "
        "u[t,x,0]==0, u[t,x,1]==0}, u, {t,0,%.4f}, {x,0,1}, {y,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "MaxSteps->4000];", T, n);
    run(buf);
    int pts[3] = { 3, 6, 9 };
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            int ix = pts[a], iy = pts[b];
            double xi = ix * h, yi = iy * h;
            double exact = exp(lam * T) * sin(PI * xi) * sin(PI * yi);
            snprintf(q, sizeof q, "First[u[%.4f, %.8f, %.8f] /. h2]", T, xi, yi);
            snprintf(lbl, sizeof lbl, "2D heat u(T,x%d,y%d)", ix, iy);
            CHECK(lbl, q, exact, 1e-4);
        }
}
static void test_pde_2d_wave(void) {
    const int n = 13;
    const double h = 1.0 / (n - 1);
    const double lam = 2.0 * (-2.0 * (1.0 - cos(PI * h)) / (h * h));
    const double omega = sqrt(-lam);
    const double T = 0.3;
    char buf[1500], q[128], lbl[64];
    snprintf(buf, sizeof buf,
        "w2 = NDSolve[{D[u[t,x,y],{t,2}]==D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}], "
        "u[0,x,y]==Sin[Pi x] Sin[Pi y], Derivative[1,0,0][u][0,x,y]==0, "
        "u[t,0,y]==0, u[t,1,y]==0, u[t,x,0]==0, u[t,x,1]==0}, u, {t,0,%.4f}, "
        "{x,0,1}, {y,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "MaxSteps->100000];", T, n);
    run(buf);
    int pts[2] = { 4, 8 };
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++) {
            int ix = pts[a], iy = pts[b];
            double xi = ix * h, yi = iy * h;
            double exact = cos(omega * T) * sin(PI * xi) * sin(PI * yi);
            snprintf(q, sizeof q, "First[u[%.4f, %.8f, %.8f] /. w2]", T, xi, yi);
            snprintf(lbl, sizeof lbl, "2D wave u(T,x%d,y%d)", ix, iy);
            CHECK(lbl, q, exact, 2e-3);
        }
}

/* ============================================================= *
 *  19. Viscous Burgers (nonlinear advection):                   *
 *      u_t = nu u_xx - u u_x + S, manufactured U=e^{-t}sin(pi x).*
 *      Exercises the u*u_x nonlinearity via the symbolic sampler.*
 * ============================================================= */
static double burgers_err(int nx, double T, double xq, double exact) {
    char buf[1500], q[96];
    snprintf(buf, sizeof buf,
        "bg = NDSolve[{D[u[t,x],t]==(1/10)D[u[t,x],{x,2}] - u[t,x] D[u[t,x],x] "
        "+ ((1/10)Pi^2-1)Exp[-t]Sin[Pi x] + Pi Exp[-2t]Sin[Pi x]Cos[Pi x], "
        "u[0,x]==Sin[Pi x], u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->4}}, "
        "Method->\"BDF\", MaxSteps->3000];", T, nx);
    run(buf);
    snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. bg]", T, xq);
    double v;
    return eval_double(q, &v) ? fabs(v - exact) : 1e9;
}
static void test_burgers(void) {
    const double T = 0.3, xq = 0.5, exact = exp(-T) * sin(PI * xq);
    double e21 = burgers_err(21, T, xq, exact);
    double e41 = burgers_err(41, T, xq, exact);
    printf("ok:   Burgers nx=21 err=%.2e  nx=41 err=%.2e\n", e21, e41);
    check_true("viscous Burgers accurate", e41 < 1e-4, "|err| < 1e-4 at nx=41");
    check_true("Burgers converges (u u_x nonlinearity)", e41 < e21, "refines");
}

/* ============================================================= *
 *  20. Arbitrary precision (MPFR) PDE — non-stiff wave.         *
 *      The output is an MPFR-valued InterpolatingFunction; node  *
 *      values carry the working precision.                      *
 * ============================================================= */
static void test_mpfr_pde(void) {
    char hd[64];
    run("wv = NDSolve[{D[u[t,x],{t,2}]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "Derivative[1,0][u][0,x]==0, u[t,0]==0, u[t,1]==0}, u, {t,0,1/50}, "
        "{x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->7,\"DifferenceOrder\"->2}}, "
        "WorkingPrecision->25, PrecisionGoal->8, MaxSteps->5000];");
    eval_head("wv", hd, sizeof hd);
    check_true("MPFR PDE returns a solution", strcmp(hd, "List") == 0, hd[0] ? hd : "(null)");
    /* correctness: the value at the initial slice reproduces the IC */
    CHECK("MPFR PDE u(0,1/3) = Sin[Pi/3]",
          "First[u[N[0,25], N[1/3,25]] /. wv]", sin(PI / 3.0), 1e-12);
    /* high precision: the node value matches the reference well beyond machine
     * epsilon (the subtraction is done in MPFR, then read back). */
    double dv;
    bool ok = eval_double("First[u[N[0,25], N[1/2,25]] /. wv] - N[1, 25]", &dv);
    check_true("MPFR node value beyond machine precision", ok && fabs(dv) < 1e-17,
               "|u(0,1/2) - 1| < 1e-17");
}

/* ============================================================= *
 *  21. Complex PDE — free-particle Schrödinger  i psi_t = -psi_xx.*
 *      Semi-discrete eigenmode psi = e^{i lam t} sin(pi x),      *
 *      lam = -(2/h^2)(1-cos(pi h)):  Re = cos(lam t) sin(pi x),  *
 *      Im = sin(lam t) sin(pi x), and |psi|^2 = sin^2(pi x)      *
 *      conserved at every node.                                  *
 * ============================================================= */
static void test_schrodinger(void) {
    const int nx = 13;
    const double T = 0.03, h = 1.0 / (nx - 1), lam = disc_lambda(nx);
    char buf[1200], q[128], lbl[64];
    snprintf(buf, sizeof buf,
        "sc = NDSolve[{I D[u[t,x],t]==-D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
        "u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "MaxSteps->100000];", T, nx);
    run(buf);
    char hd[64]; eval_head("sc", hd, sizeof hd);
    check_true("Schrödinger result is a List", strcmp(hd, "List") == 0, hd[0] ? hd : "(null)");
    for (int i = 3; i <= nx - 2; i += 3) {
        double xi = i * h, s = sin(PI * xi);
        snprintf(q, sizeof q, "Re[First[u[%.4f, %.10f] /. sc]]", T, xi);
        snprintf(lbl, sizeof lbl, "Schrodinger Re(T,x%d)", i);
        CHECK(lbl, q, cos(lam * T) * s, 1e-4);
        snprintf(q, sizeof q, "Im[First[u[%.4f, %.10f] /. sc]]", T, xi);
        snprintf(lbl, sizeof lbl, "Schrodinger Im(T,x%d)", i);
        CHECK(lbl, q, sin(lam * T) * s, 1e-4);
        snprintf(q, sizeof q, "Abs[First[u[%.4f, %.10f] /. sc]]^2", T, xi);
        snprintf(lbl, sizeof lbl, "Schrodinger |psi|^2(T,x%d)", i);
        CHECK(lbl, q, s * s, 1e-4);   /* norm conserved per node */
    }
}

/* ============================================================= *
 *  22. Schrödinger with a potential well  i psi_t=-psi_xx+V psi. *
 *      No closed form, but the total norm sum|psi|^2 is a        *
 *      conserved quantity (unitary evolution) — a strong stress. *
 * ============================================================= */
static void test_schrodinger_potential(void) {
    const int nx = 15;
    const double T = 0.02, h = 1.0 / (nx - 1);
    char buf[1300], q[128];
    snprintf(buf, sizeof buf,
        "sp = NDSolve[{I D[u[t,x],t]==-D[u[t,x],{x,2}] + 200 (x-1/2)^2 u[t,x], "
        "u[0,x]==Sin[Pi x], u[t,0]==0, u[t,1]==0}, u, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}, "
        "MaxSteps->200000];", T, nx);
    run(buf);
    double nT = 0.0, n0 = 0.0;
    for (int i = 1; i <= nx - 2; i++) {
        double xi = i * h, re, im;
        snprintf(q, sizeof q, "Re[First[u[%.4f, %.10f] /. sp]]", T, xi);
        if (!eval_double(q, &re)) { printf("FAIL: potential Re query\n"); failures++; return; }
        snprintf(q, sizeof q, "Im[First[u[%.4f, %.10f] /. sp]]", T, xi);
        if (!eval_double(q, &im)) { printf("FAIL: potential Im query\n"); failures++; return; }
        nT += re * re + im * im;
        double s = sin(PI * xi);
        n0 += s * s;
    }
    printf("ok:   Schrodinger+V total norm: t=0 %.8g  T %.8g  rel-drift %.2e\n",
           n0, nT, fabs(nT - n0) / n0);
    check_true("Schrödinger potential: norm conserved", fabs(nT - n0) / n0 < 2e-3,
               "sum|psi|^2 conserved under unitary evolution");
}

int main(void) {
    mute_stderr_once();
    core_init();

    test_heat_eigenmode();
    test_heat_spatial_convergence();
    test_wave_eigenmode();
    test_reaction_diffusion();
    test_inhomogeneous_bc();
    test_time_dependent_bc();
    test_nonlinear_manufactured();
    test_structure_and_equivalence();
    test_offgrid_interpolation();
    test_advection_diffusion();
    test_large_grid_stress();
    test_difference_order();
    test_neumann_bc();
    test_robin_bc();
    test_periodic_bc();
    test_compiled_operator();
    test_operator_scale();
    test_pde_2d_heat();
    test_pde_2d_wave();
    test_burgers();
    test_mpfr_pde();
    test_schrodinger();
    test_schrodinger_potential();

    if (failures == 0) printf("\nAll NDSolve PDE tests passed.\n");
    else printf("\n%d NDSolve PDE test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
