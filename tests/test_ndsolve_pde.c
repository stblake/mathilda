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
    /* The compiled operator (A*Y matvec, exact constant Jacobian) and the
     * symbolic sampler evaluate the same linear system but differ at the last
     * bit; the adaptive variable-order BDF amplifies that into slightly
     * different meshes, so the two solutions agree to the solution tolerance
     * (~1e-8 here) rather than bit-for-bit.  This still verifies the fast path
     * computes the same operator — a wrong A would diverge by O(1). */
    check_true("compiled operator == symbolic RHS", fabs(vt - vf) < 1e-6,
               "fast path agrees with symbolic path (to solution tolerance)");

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
 *  18b. 2-D insulated heat (homogeneous Neumann on all edges).  *
 *       IC cos(pi x) cos(pi y); the exact PDE solution decays as *
 *       e^{-2 pi^2 t} cos(pi x) cos(pi y).  Checks accuracy and  *
 *       grid convergence (spatial error shrinks as nx grows).    *
 * ============================================================= */
static double heat2d_neumann_err(int n, double T, double xq, double yq, double exact) {
    char buf[1500], q[128];
    snprintf(buf, sizeof buf,
        "hn = NDSolve[{D[u[t,x,y],t]==D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}], "
        "u[0,x,y]==Cos[Pi x] Cos[Pi y], "
        "Derivative[0,1,0][u][t,0,y]==0, Derivative[0,1,0][u][t,1,y]==0, "
        "Derivative[0,0,1][u][t,x,0]==0, Derivative[0,0,1][u][t,x,1]==0}, u, "
        "{t,0,%.4f}, {x,0,1}, {y,0,1}, Method->{\"MethodOfLines\","
        "\"SpatialDiscretization\"->{\"TensorProductGrid\",\"MinPoints\"->%d,"
        "\"DifferenceOrder\"->4}}, MaxSteps->4000];", T, n);
    run(buf);
    snprintf(q, sizeof q, "First[u[%.4f, %.8f, %.8f] /. hn]", T, xq, yq);
    double v;
    return eval_double(q, &v) ? fabs(v - exact) : 1e9;
}
static void test_pde_2d_neumann(void) {
    const double T = 0.03;
    double xq = 0.25, yq = 0.25;
    double exact = exp(-2.0 * PI * PI * T) * cos(PI * xq) * cos(PI * yq);
    double e15 = heat2d_neumann_err(15, T, xq, yq, exact);
    double e29 = heat2d_neumann_err(29, T, xq, yq, exact);
    printf("ok:   2D Neumann heat nx=15 err=%.2e nx=29 err=%.2e\n", e15, e29);
    check_true("2D Neumann accurate", e29 < 5e-4, "|err| < 5e-4 at nx=29");
    check_true("2D Neumann converges", e29 < e15, "finer grid smaller error");
}

/* ============================================================= *
 *  18c. 2-D Robin/Neumann steady state.  u = x + y is harmonic  *
 *       (u_xx=u_yy=0), so u_t=0 and it is stationary.  Two edges  *
 *       carry a genuine Robin condition a*u + b*u_n + r == 0     *
 *       (both coefficients nonzero), two carry Neumann.  The FD   *
 *       elimination is exact for a linear field, so interior,    *
 *       edge, and corner values all reproduce x+y to ~machine.   *
 * ============================================================= */
static void test_pde_2d_robin_steady(void) {
    const double T = 0.1;
    run("rb = NDSolve[{D[u[t,x,y],t]==D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}], "
        "u[0,x,y]==x+y, "
        "Derivative[0,1,0][u][t,0,y]==1, "
        "u[t,1,y]+Derivative[0,1,0][u][t,1,y]==2+y, "
        "Derivative[0,0,1][u][t,x,0]==1, "
        "u[t,x,1]+Derivative[0,0,1][u][t,x,1]==2+x}, u, "
        "{t,0,0.1}, {x,0,1}, {y,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->13,\"DifferenceOrder\"->4}}, "
        "MaxSteps->4000];");
    struct { double x, y; const char* tag; } probes[] = {
        { 0.3, 0.7, "interior" },   /* interior node               */
        { 0.0, 0.5, "x=0 Neumann" },/* left edge (Neumann)         */
        { 1.0, 0.5, "x=1 Robin"    },/* right edge (Robin)          */
        { 0.5, 1.0, "y=1 Robin"    },/* top edge (Robin)            */
        { 0.0, 0.0, "corner"       },/* corner (two edges meet)     */
    };
    char q[128], lbl[64];
    for (size_t i = 0; i < sizeof probes / sizeof probes[0]; i++) {
        snprintf(q, sizeof q, "First[u[%.4f, %.6f, %.6f] /. rb]", T, probes[i].x, probes[i].y);
        snprintf(lbl, sizeof lbl, "2D Robin steady (%s)", probes[i].tag);
        CHECK(lbl, q, probes[i].x + probes[i].y, 1e-8);
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

/* ============================================================= *
 *  23. Multi-mode heat: superposition of two eigenmodes vs the  *
 *      TRUE PDE.  IC sin(pi x) + 1/2 sin(3 pi x) evolves to      *
 *      e^{-pi^2 t} sin(pi x) + 1/2 e^{-9 pi^2 t} sin(3 pi x).    *
 *      The modes decay at rates differing 9x (multi-scale        *
 *      stiffness); a fine 4th-order grid + BDF resolves both.    *
 * ============================================================= */
static void test_multimode_heat(void) {
    const double T = 0.02;
    run("mm = NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], "
        "u[0,x]==Sin[Pi x]+(1/2)Sin[3 Pi x], u[t,0]==0, u[t,1]==0}, u, "
        "{t,0,0.02}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->81,\"DifferenceOrder\"->4}}, "
        "Method->\"BDF\", MaxSteps->5000];");
    double xs[3] = { 0.25, 0.5, 0.75 };
    char q[96], lbl[64];
    for (int i = 0; i < 3; i++) {
        double x = xs[i];
        double exact = exp(-PI * PI * T) * sin(PI * x)
                     + 0.5 * exp(-9 * PI * PI * T) * sin(3 * PI * x);
        snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. mm]", T, x);
        snprintf(lbl, sizeof lbl, "multimode heat u(T,%.2f)", x);
        CHECK(lbl, q, exact, 1e-3);
    }
}

/* ============================================================= *
 *  24. Variable-coefficient diffusion  u_t = (1+x^2) u_xx + S.   *
 *      Manufactured U = e^{-t} sin(pi x); S = U_t - (1+x^2)U_xx  *
 *      makes it exact.  Exercises a spatially varying diffusion  *
 *      coefficient through the symbolic sampler; order-4 stencil.*
 * ============================================================= */
static void test_variable_coefficient(void) {
    const double T = 0.3;
    run("vc = NDSolve[{D[u[t,x],t]==(1+x^2)D[u[t,x],{x,2}] "
        "- Exp[-t]Sin[Pi x] + (1+x^2)Pi^2 Exp[-t]Sin[Pi x], "
        "u[0,x]==Sin[Pi x], u[t,0]==0, u[t,1]==0}, u, {t,0,0.3}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->41,\"DifferenceOrder\"->4}}, "
        "Method->\"BDF\", MaxSteps->4000];");
    double xs[3] = { 0.25, 0.5, 0.75 };
    char q[96], lbl[64];
    for (int i = 0; i < 3; i++) {
        double x = xs[i], exact = exp(-T) * sin(PI * x);
        snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. vc]", T, x);
        snprintf(lbl, sizeof lbl, "var-coeff diffusion u(T,%.2f)", x);
        CHECK(lbl, q, exact, 2e-3);
    }
}

/* ============================================================= *
 *  25. Klein-Gordon  u_tt = u_xx - m^2 u  (m=2).  sin(pi x) is   *
 *      an exact discrete-Laplacian eigenvector and the mass term *
 *      is diagonal, so the semi-discrete solution is exactly     *
 *      cos(omega t) sin(pi x), omega = sqrt(m^2 + |lambda_disc|).*
 *      Isolates time-integration error at temporal order 2.      *
 * ============================================================= */
static void test_klein_gordon(void) {
    const int nx = 21;
    const double T = 0.5, h = 1.0 / (nx - 1), m = 2.0;
    const double omega = sqrt(m * m - disc_lambda(nx));
    run("kg = NDSolve[{D[u[t,x],{t,2}]==D[u[t,x],{x,2}] - 4 u[t,x], "
        "u[0,x]==Sin[Pi x], Derivative[1,0][u][0,x]==0, u[t,0]==0, u[t,1]==0}, "
        "u, {t,0,0.5}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->21,\"DifferenceOrder\"->2}}, MaxSteps->100000];");
    char q[96], lbl[64];
    for (int i = 4; i <= nx - 2; i += 5) {
        double xi = i * h, exact = cos(omega * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.6f, %.10f] /. kg]", T, xi);
        snprintf(lbl, sizeof lbl, "Klein-Gordon u(T,x%d)", i);
        CHECK(lbl, q, exact, 1e-4);
    }
}

/* ============================================================= *
 *  26. Damped (telegraph) wave  u_tt = u_xx - a u_t  (a=1).      *
 *      A first-order time derivative sits INSIDE a second-order- *
 *      in-time PDE: the solver must introduce v=u_t and feed it  *
 *      back into the RHS.  sin(pi x) stays an eigenmode; each     *
 *      node obeys  U'' + a U' + |lam| U = 0  exactly, so          *
 *      U(t) = e^{-a t/2}(cos(wd t) + a/(2 wd) sin(wd t)),         *
 *      wd = sqrt(|lam| - a^2/4).  Exact semi-discrete reference.  *
 * ============================================================= */
static void test_damped_wave(void) {
    const int nx = 21;
    const double T = 0.5, h = 1.0 / (nx - 1), a = 1.0;
    const double absl = -disc_lambda(nx), wd = sqrt(absl - a * a / 4.0);
    const double U = exp(-a * T / 2.0) * (cos(wd * T) + (a / (2.0 * wd)) * sin(wd * T));
    run("dw = NDSolve[{D[u[t,x],{t,2}]==D[u[t,x],{x,2}] - D[u[t,x],t], "
        "u[0,x]==Sin[Pi x], Derivative[1,0][u][0,x]==0, u[t,0]==0, u[t,1]==0}, "
        "u, {t,0,0.5}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->21,\"DifferenceOrder\"->2}}, MaxSteps->100000];");
    char q[96], lbl[64];
    for (int i = 4; i <= nx - 2; i += 5) {
        double xi = i * h, exact = U * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.6f, %.10f] /. dw]", T, xi);
        snprintf(lbl, sizeof lbl, "damped wave u(T,x%d)", i);
        CHECK(lbl, q, exact, 2e-4);
    }
}

/* ============================================================= *
 *  27. Nonlinear (porous-medium-type) diffusion                 *
 *      u_t = u u_xx + u_x^2 + S = (u u_x)_x + S.                 *
 *      The diffusion coefficient depends on the solution ->      *
 *      full nonlinear symbolic sampler + Newton per BDF step.    *
 *      Manufactured U = e^{-t} sin(pi x) + 2 (kept positive;     *
 *      inhomogeneous Dirichlet u=2 at both ends).                *
 * ============================================================= */
static void test_nonlinear_diffusion(void) {
    const double T = 0.1;
    run("nd = NDSolve[{D[u[t,x],t]==u[t,x] D[u[t,x],{x,2}] + D[u[t,x],x]^2 "
        "- Exp[-t]Sin[Pi x] + (Exp[-t]Sin[Pi x]+2)Pi^2 Exp[-t]Sin[Pi x] "
        "- Pi^2 Exp[-2t]Cos[Pi x]^2, u[0,x]==Sin[Pi x]+2, u[t,0]==2, u[t,1]==2}, "
        "u, {t,0,0.1}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->21,\"DifferenceOrder\"->4}}, "
        "Method->\"BDF\", MaxSteps->10000];");
    double xs[3] = { 0.25, 0.5, 0.75 };
    char q[96], lbl[64];
    for (int i = 0; i < 3; i++) {
        double x = xs[i], exact = exp(-T) * sin(PI * x) + 2.0;
        snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. nd]", T, x);
        snprintf(lbl, sizeof lbl, "nonlinear diffusion u(T,%.2f)", x);
        CHECK(lbl, q, exact, 2e-3);
    }
}

/* ============================================================= *
 *  28. Wave with insulated (Neumann) ends  u_tt = u_xx,         *
 *      u_x(t,0)=u_x(t,1)=0, IC cos(pi x), zero velocity.         *
 *      True PDE solution cos(pi t) cos(pi x).  A fine 4th-order  *
 *      grid tracks it and validates one-sided Neumann            *
 *      elimination in a second-order-in-time PDE.                *
 * ============================================================= */
static void test_wave_neumann(void) {
    const double T = 0.2;
    run("wn = NDSolve[{D[u[t,x],{t,2}]==D[u[t,x],{x,2}], u[0,x]==Cos[Pi x], "
        "Derivative[1,0][u][0,x]==0, Derivative[0,1][u][t,0]==0, "
        "Derivative[0,1][u][t,1]==0}, u, {t,0,0.2}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->41,\"DifferenceOrder\"->4}}, MaxSteps->200000];");
    double xs[3] = { 0.1, 0.25, 0.75 };
    char q[96], lbl[64];
    for (int i = 0; i < 3; i++) {
        double x = xs[i], exact = cos(PI * T) * cos(PI * x);
        snprintf(q, sizeof q, "First[u[%.4f, %.4f] /. wn]", T, x);
        snprintf(lbl, sizeof lbl, "wave Neumann u(T,%.2f)", x);
        CHECK(lbl, q, exact, 2e-3);
    }
}

/* ============================================================= *
 *  29. 2-D reaction-diffusion  u_t = u_xx + u_yy + c u  (c=3).   *
 *      sin(pi x) sin(pi y) is an exact eigenmode of the discrete *
 *      2-D Laplacian, so the semi-discrete solution is           *
 *      e^{(2 lam + c) t} sin sin exactly (lx = ly = lam).        *
 * ============================================================= */
static void test_2d_reaction_diffusion(void) {
    const int n = 13;
    const double h = 1.0 / (n - 1), T = 0.02, c = 3.0;
    const double rate = 2.0 * disc_lambda(n) + c;
    run("r2 = NDSolve[{D[u[t,x,y],t]==D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}]+3 u[t,x,y], "
        "u[0,x,y]==Sin[Pi x]Sin[Pi y], u[t,0,y]==0, u[t,1,y]==0, u[t,x,0]==0, u[t,x,1]==0}, "
        "u, {t,0,0.02}, {x,0,1}, {y,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->13,\"DifferenceOrder\"->2}}, MaxSteps->5000];");
    int pts[3] = { 3, 6, 9 };
    char q[128], lbl[64];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            int ix = pts[a], iy = pts[b];
            double xi = ix * h, yi = iy * h;
            double exact = exp(rate * T) * sin(PI * xi) * sin(PI * yi);
            snprintf(q, sizeof q, "First[u[%.4f, %.8f, %.8f] /. r2]", T, xi, yi);
            snprintf(lbl, sizeof lbl, "2D reaction u(T,x%d,y%d)", ix, iy);
            CHECK(lbl, q, exact, 1e-4);
        }
}

/* ============================================================= *
 *  30. 2-D anisotropic heat  u_t = 2 u_xx + u_yy.  Unequal       *
 *      diffusion coefficients: the sin sin eigenmode decays at   *
 *      rate 2 lx + ly = 3 lam on the square grid (lx = ly = lam).*
 *      Exact semi-discrete reference.                            *
 * ============================================================= */
static void test_2d_anisotropic(void) {
    const int n = 13;
    const double h = 1.0 / (n - 1), T = 0.02;
    const double rate = 3.0 * disc_lambda(n);
    run("a2 = NDSolve[{D[u[t,x,y],t]==2 D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}], "
        "u[0,x,y]==Sin[Pi x]Sin[Pi y], u[t,0,y]==0, u[t,1,y]==0, u[t,x,0]==0, u[t,x,1]==0}, "
        "u, {t,0,0.02}, {x,0,1}, {y,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->13,\"DifferenceOrder\"->2}}, MaxSteps->5000];");
    int pts[3] = { 3, 6, 9 };
    char q[128], lbl[64];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            int ix = pts[a], iy = pts[b];
            double xi = ix * h, yi = iy * h;
            double exact = exp(rate * T) * sin(PI * xi) * sin(PI * yi);
            snprintf(q, sizeof q, "First[u[%.4f, %.8f, %.8f] /. a2]", T, xi, yi);
            snprintf(lbl, sizeof lbl, "2D aniso u(T,x%d,y%d)", ix, iy);
            CHECK(lbl, q, exact, 1e-4);
        }
}

/* ============================================================= *
 *  31. Time-integrator equivalence.  The same mild heat problem *
 *      solved with an explicit RK4 stepper and with stiff BDF    *
 *      must both track the semi-discrete eigenmode and agree     *
 *      with each other -- verifies Method-> routes to genuinely  *
 *      different integrators computing the same field.           *
 * ============================================================= */
static void test_method_equivalence_rk_bdf(void) {
    const int nx = 11;
    const double T = 0.02, lam = disc_lambda(nx);
    const char* stub =
        "NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], u[t,0]==0, u[t,1]==0}, "
        "u, {t,0,0.02}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->11,\"DifferenceOrder\"->2}}, ";
    char buf[1200];
    snprintf(buf, sizeof buf, "mrk = %sMethod->\"RungeKutta\", MaxSteps->50000];", stub);
    run(buf);
    snprintf(buf, sizeof buf, "mbd = %sMethod->\"BDF\", MaxSteps->5000];", stub);
    run(buf);
    double vr, vb;
    bool ok = eval_double("First[u[0.02,0.5] /. mrk]", &vr) &&
              eval_double("First[u[0.02,0.5] /. mbd]", &vb);
    double exact = exp(lam * T) * sin(PI * 0.5);
    check_true("RK4 tracks semi-discrete", ok && fabs(vr - exact) < 1e-4, "explicit RK accurate");
    check_true("BDF tracks semi-discrete", ok && fabs(vb - exact) < 1e-4, "stiff BDF accurate");
    check_true("RK4 == BDF (integrator equivalence)", ok && fabs(vr - vb) < 1e-4,
               "two different steppers agree");
}

/* ============================================================= *
 *  32. Graceful degradation.  Cases outside the MoL scope must  *
 *      return the NDSolve expression UNEVALUATED (head NDSolve)  *
 *      -- not crash, not fabricate a solution.  Systems of PDEs, *
 *      > 2 spatial dims, temporal order > 2, and 2-D periodic    *
 *      BCs are all deferred and must be left symbolic.           *
 * ============================================================= */
/* ============================================================= *
 *  33. Coupled system: reaction-diffusion, manufactured        *
 *      u_t = u_xx + w + Su,  w_t = w_xx + u + Sv  with exact    *
 *      u = e^{-t} sin(pi x), w = e^{-2t} sin(pi x).  Each        *
 *      equation reads the OTHER function (genuine coupling).    *
 * ============================================================= */
static void test_system_reaction_diffusion(void) {
    const double T = 0.2;
    char buf[1400], q[128], lbl[64];
    snprintf(buf, sizeof buf,
        "rd = NDSolve[{"
        "D[u[t,x],t]==D[u[t,x],{x,2}]+w[t,x]+((Pi^2-1)Exp[-t]Sin[Pi x]-Exp[-2t]Sin[Pi x]), "
        "D[w[t,x],t]==D[w[t,x],{x,2}]+u[t,x]+((Pi^2-2)Exp[-2t]Sin[Pi x]-Exp[-t]Sin[Pi x]), "
        "u[0,x]==Sin[Pi x], w[0,x]==Sin[Pi x], u[t,0]==0,u[t,1]==0,w[t,0]==0,w[t,1]==0}, "
        "{u,w}, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->41,\"DifferenceOrder\"->2}}, MaxSteps->20000];", T);
    run(buf);
    const double xs[3] = { 0.25, 0.5, 0.75 };
    for (int i = 0; i < 3; i++) {
        double xi = xs[i];
        double ue = exp(-T) * sin(PI * xi), we = exp(-2 * T) * sin(PI * xi);
        snprintf(q, sizeof q, "First[u[%.4f,%.4f] /. rd]", T, xi);
        snprintf(lbl, sizeof lbl, "reac-diff system u(T,%.2f)", xi);
        CHECK(lbl, q, ue, 2e-3);
        snprintf(q, sizeof q, "First[w[%.4f,%.4f] /. rd]", T, xi);
        snprintf(lbl, sizeof lbl, "reac-diff system w(T,%.2f)", xi);
        CHECK(lbl, q, we, 2e-3);
    }
}

/* ============================================================= *
 *  34. Shallow-water system, linearized standing gravity wave.  *
 *      h_t+(h u)_x=0, u_t+u u_x+g h_x=0 (g=1), periodic.  Small  *
 *      A: h = 1 + A cos(2 pi x) => standing wave at c=sqrt(gH)=1;*
 *      at t=1/2 (half period) h = 1 - A cos(2 pi x), u ~ 0.      *
 * ============================================================= */
static void test_system_shallow_water(void) {
    const double A = 1e-3;
    char buf[1400], q[128], lbl[64];
    snprintf(buf, sizeof buf,
        "sw = NDSolve[{D[h[t,x],t]+D[h[t,x] u[t,x],x]==0, "
        "D[u[t,x],t]+u[t,x] D[u[t,x],x]+D[h[t,x],x]==0, "
        "h[0,x]==1+%.6f Cos[2 Pi x], u[0,x]==0, h[t,0]==h[t,1], u[t,0]==u[t,1]}, "
        "{h,u}, {t,0,0.5}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->81}}];", A);
    run(buf);
    const double xs[3] = { 0.0, 0.5, 0.25 };
    for (int i = 0; i < 3; i++) {
        double xi = xs[i], he = 1.0 - A * cos(2 * PI * xi);
        snprintf(q, sizeof q, "First[h[0.5,%.4f] /. sw]", xi);
        snprintf(lbl, sizeof lbl, "shallow-water h(1/2,%.2f)", xi);
        CHECK(lbl, q, he, 5e-4);
    }
    double uval;
    if (eval_double("First[u[0.5,0.3] /. sw]", &uval))
        check_true("shallow-water u ~ 0 at half period", fabs(uval) < 1e-3, "|u|<1e-3");
}

/* ============================================================= *
 *  35. Coupled 2nd-order (wave) system via normal modes.        *
 *      u_tt=u_xx-u+w, w_tt=w_xx-w+u; ICs u=sin(pi x),w=0,rest 0. *
 *      s=u+w solves the wave eqn, d=u-w solves u_tt=u_xx-2u, so  *
 *      u=(cos os t + cos od t)/2 sin(pi x), w=(cos os t - cos od *
 *      t)/2 sin(pi x), os=sqrt(-lam), od=sqrt(-lam+2) EXACT on   *
 *      the order-2 semi-discrete operator.                       *
 * ============================================================= */
static void test_system_coupled_wave(void) {
    const int nx = 21;
    const double T = 0.3, lam = disc_lambda(nx);
    const double os = sqrt(-lam), od = sqrt(-lam + 2.0);
    char buf[1400], q[128];
    snprintf(buf, sizeof buf,
        "cw = NDSolve[{D[u[t,x],t,t]==D[u[t,x],{x,2}]-u[t,x]+w[t,x], "
        "D[w[t,x],t,t]==D[w[t,x],{x,2}]-w[t,x]+u[t,x], "
        "u[0,x]==Sin[Pi x], w[0,x]==0, Derivative[1,0][u][0,x]==0, Derivative[1,0][w][0,x]==0, "
        "u[t,0]==0,u[t,1]==0,w[t,0]==0,w[t,1]==0}, {u,w}, {t,0,%.4f}, {x,0,1}, "
        "Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d,\"DifferenceOrder\"->2}}];", T, nx);
    run(buf);
    double ue = 0.5 * (cos(os * T) + cos(od * T)) * sin(PI * 0.5);
    double we = 0.5 * (cos(os * T) - cos(od * T)) * sin(PI * 0.5);
    snprintf(q, sizeof q, "First[u[%.4f,0.5] /. cw]", T);
    CHECK("coupled wave u(T,1/2)", q, ue, 1e-3);
    snprintf(q, sizeof q, "First[w[%.4f,0.5] /. cw]", T);
    CHECK("coupled wave w(T,1/2)", q, we, 1e-3);
}

/* ============================================================= *
 *  36. Scalar linear advection u_t + u_x = 0, periodic.         *
 *      Exact u(t,x)=sin(2 pi(x-t)); at t=1/4, u(x)=-cos(2 pi x). *
 *      Centered is near-exact; donor-cell upwind advects with    *
 *      correct phase but first-order amplitude diffusion that    *
 *      shrinks as the grid refines; the reversed wind stays      *
 *      stable (wrong upwind direction would blow up).            *
 * ============================================================= */
static void test_advection_upwind(void) {
    const char* base =
        "%s = NDSolve[{D[u[t,x],t]%sD[u[t,x],x]==0, u[0,x]==Sin[2 Pi x], u[t,0]==u[t,1]}, "
        "u, {t,0,0.25}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->%d%s}}];";
    char buf[900];
    /* centered: near-exact -1 */
    snprintf(buf, sizeof buf, base, "aC", "+", 81, "");
    run(buf);
    CHECK("advection centered u(1/4,0)", "First[u[0.25,0.0] /. aC]", -1.0, 1e-2);
    /* donor-cell upwind, coarse and fine */
    snprintf(buf, sizeof buf, base, "aU", "+", 81, ",\"DifferenceOrder\"->1");
    run(buf);
    snprintf(buf, sizeof buf, base, "aF", "+", 321, ",\"DifferenceOrder\"->1");
    run(buf);
    double uc, uf, ur;
    bool okc = eval_double("First[u[0.25,0.0] /. aU]", &uc);
    bool okf = eval_double("First[u[0.25,0.0] /. aF]", &uf);
    check_true("upwind advects with correct phase, diffused",
               okc && uc < -0.85 && uc > -0.999, "value in (-0.999,-0.85)");
    check_true("upwind amplitude converges on refinement",
               okc && okf && fabs(uf) > fabs(uc), "|fine| > |coarse|");
    /* reversed wind: donor-cell must stay stable (correct upwind direction) */
    snprintf(buf, sizeof buf, base, "aR", "-", 81, ",\"DifferenceOrder\"->1");
    run(buf);
    bool okr = eval_double("First[u[0.25,0.0] /. aR]", &ur);
    check_true("reversed-wind upwind stays stable",
               okr && ur > 0.5 && ur < 1.0, "value in (0.5,1.0)");
}

/* ============================================================= *
 *  37. Upwind monotonicity: advecting a sharp top-hat, centered *
 *      high-order stencils ring (Gibbs over/undershoot) while    *
 *      the upwind schemes stay within the data bounds [0,1].     *
 * ============================================================= */
static void test_upwind_monotonicity(void) {
    const char* base =
        "%s = NDSolve[{D[u[t,x],t]+D[u[t,x],x]==0, "
        "u[0,x]==(Tanh[(x-3/10)/(1/50)]-Tanh[(x-6/10)/(1/50)])/2, u[t,0]==u[t,1]}, "
        "u, {t,0,0.3}, {x,0,1}, Method->{\"MethodOfLines\",\"SpatialDiscretization\"->"
        "{\"TensorProductGrid\",\"MinPoints\"->81%s}}];";
    char buf[900];
    snprintf(buf, sizeof buf, base, "hC", "");
    run(buf);
    snprintf(buf, sizeof buf, base, "hU", ",\"Upwind\"->True");
    run(buf);
    double cmax, cmin, umax, umin;
    /* build min/max over a sampled grid for each solution */
    run("cMax = Max[Table[(u/.First[hC])[0.3,x], {x,0.01,0.99,0.01}]];"
        "cMin = Min[Table[(u/.First[hC])[0.3,x], {x,0.01,0.99,0.01}]];"
        "uMax = Max[Table[(u/.First[hU])[0.3,x], {x,0.01,0.99,0.01}]];"
        "uMin = Min[Table[(u/.First[hU])[0.3,x], {x,0.01,0.99,0.01}]];");
    bool o1 = eval_double("cMax", &cmax), o2 = eval_double("cMin", &cmin);
    bool o3 = eval_double("uMax", &umax), o4 = eval_double("uMin", &umin);
    check_true("centered rings on sharp front (Gibbs)",
               o1 && o2 && (cmax > 1.01 || cmin < -0.01), "over/undershoot present");
    check_true("upwind stays within data bounds [0,1]",
               o3 && o4 && umax < 1.005 && umin > -0.005, "no over/undershoot");
}

static void test_unsupported_guards(void) {
    char hd[64];
    /* well-posed systems now solve; malformed ones must stay unevaluated. */
    eval_head("NDSolve[{D[u[t,x],t]==D[u[t,x],{x,2}], D[w[t,x],t]==D[w[t,x],{x,2}], "
              "u[0,x]==Sin[Pi x], w[0,x]==Sin[Pi x], u[t,0]==0,u[t,1]==0,w[t,0]==0}, "
              "{u,w}, {t,0,0.1}, {x,0,1}]", hd, sizeof hd);
    check_true("system missing a BC stays unevaluated", strcmp(hd, "NDSolve") == 0, hd[0] ? hd : "(null)");

    eval_head("NDSolve[{D[u[t,x],t]+D[w[t,x],t]==D[u[t,x],{x,2}], D[w[t,x],t]==D[w[t,x],{x,2}], "
              "u[0,x]==Sin[Pi x], w[0,x]==Sin[Pi x], u[t,0]==0,u[t,1]==0,w[t,0]==0,w[t,1]==0}, "
              "{u,w}, {t,0,0.1}, {x,0,1}]", hd, sizeof hd);
    check_true("coupled mass matrix stays unevaluated", strcmp(hd, "NDSolve") == 0, hd[0] ? hd : "(null)");

    eval_head("NDSolve[{I D[u[t,x],t]==D[u[t,x],{x,2}], D[w[t,x],t]==D[w[t,x],{x,2}], "
              "u[0,x]==Sin[Pi x], w[0,x]==Sin[Pi x], u[t,0]==0,u[t,1]==0,w[t,0]==0,w[t,1]==0}, "
              "{u,w}, {t,0,0.1}, {x,0,1}]", hd, sizeof hd);
    check_true("complex-valued system stays unevaluated", strcmp(hd, "NDSolve") == 0, hd[0] ? hd : "(null)");

    eval_head("NDSolve[{D[u[t,x,y,z],t]==D[u[t,x,y,z],{x,2}], u[0,x,y,z]==1, "
              "u[t,0,y,z]==0,u[t,1,y,z]==0,u[t,x,0,z]==0,u[t,x,1,z]==0,"
              "u[t,x,y,0]==0,u[t,x,y,1]==0}, u, {t,0,0.1},{x,0,1},{y,0,1},{z,0,1}]",
              hd, sizeof hd);
    check_true(">2 spatial dims stays unevaluated", strcmp(hd, "NDSolve") == 0, hd[0] ? hd : "(null)");

    eval_head("NDSolve[{D[u[t,x],{t,3}]==D[u[t,x],{x,2}], u[0,x]==Sin[Pi x], "
              "u[t,0]==0,u[t,1]==0}, u, {t,0,0.1},{x,0,1}]", hd, sizeof hd);
    check_true("temporal order 3 stays unevaluated", strcmp(hd, "NDSolve") == 0, hd[0] ? hd : "(null)");

    eval_head("NDSolve[{D[u[t,x,y],t]==D[u[t,x,y],{x,2}]+D[u[t,x,y],{y,2}], "
              "u[0,x,y]==Sin[2 Pi x]Sin[Pi y], u[t,0,y]==u[t,1,y], "
              "u[t,x,0]==0, u[t,x,1]==0}, u, {t,0,0.01},{x,0,1},{y,0,1}]", hd, sizeof hd);
    check_true("2-D periodic BC stays unevaluated", strcmp(hd, "NDSolve") == 0, hd[0] ? hd : "(null)");
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
    test_pde_2d_neumann();
    test_pde_2d_robin_steady();
    test_burgers();
    test_mpfr_pde();
    test_schrodinger();
    test_schrodinger_potential();
    test_multimode_heat();
    test_variable_coefficient();
    test_klein_gordon();
    test_damped_wave();
    test_nonlinear_diffusion();
    test_wave_neumann();
    test_2d_reaction_diffusion();
    test_2d_anisotropic();
    test_method_equivalence_rk_bdf();
    test_system_reaction_diffusion();
    test_system_shallow_water();
    test_system_coupled_wave();
    test_advection_upwind();
    test_upwind_monotonicity();
    test_unsupported_guards();

    if (failures == 0) printf("\nAll NDSolve PDE tests passed.\n");
    else printf("\n%d NDSolve PDE test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
