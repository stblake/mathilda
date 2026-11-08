/* The finite element method on planar regions: NDEigenvalues, NDEigensystem
 * and NDSolve over {x, y} in a region (src/numerical_calculus/fem*.c).
 *
 * References are closed forms: Dirichlet eigenvalues of the square (pi^2
 * (m^2 + n^2)), the disk (squared Bessel zeros), the equilateral triangle
 * (16 pi^2/3 for side 1), the L-shape of three unit squares (lambda_1 =
 * 9.6397238440219..., Fox-Henrici-Moler; lambda_3 = 2 pi^2) and the annulus;
 * Poisson problems with polynomial solutions (reproduced by P2 elements to
 * rounding) and manufactured smooth ones. Also the mesher's invariants
 * (area, angles) and the sparse Cholesky against a direct solve.
 *
 * Soft asserts. Run: ./fem_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "fem_mesh.h"
#include "fem_chol.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static const double PI = 3.14159265358979323846;
static int failures = 0;

static bool eval_double(const char* input, double* out) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    bool ok = false;
    if (r) {
        if (r->type == EXPR_REAL)         { *out = r->data.real; ok = true; }
        else if (r->type == EXPR_INTEGER) { *out = (double)r->data.integer; ok = true; }
    }
    expr_free(e); expr_free(r);
    return ok;
}

static void run(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    expr_free(e); expr_free(r);
}

static void check_rel(const char* label, const char* input, double expected, double rtol) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric [%s]\n", label, input); failures++; return; }
    double err = fabs(v - expected) / (fabs(expected) > 1e-300 ? fabs(expected) : 1.0);
    if (!(err <= rtol)) {
        printf("FAIL: %s -> %.15g (expected %.15g, rel err %.2e, tol %.1e)\n", label, v, expected, err, rtol);
        failures++;
    } else {
        printf("ok:   %-52s = %.12g  (rel %.1e)\n", label, v, err);
    }
}

static void check_abs(const char* label, const char* input, double expected, double tol) {
    double v;
    if (!eval_double(input, &v)) { printf("FAIL: %s -> not numeric [%s]\n", label, input); failures++; return; }
    if (!(fabs(v - expected) <= tol)) {
        printf("FAIL: %s -> %.15g (expected %.15g, |err| %.2e, tol %.1e)\n", label, v, expected, fabs(v - expected), tol);
        failures++;
    } else {
        printf("ok:   %-52s = %.12g  (|err| %.1e)\n", label, v, fabs(v - expected));
    }
}

static void check_true(const char* label, bool cond, const char* detail) {
    if (cond) printf("ok:   %-52s [%s]\n", label, detail);
    else { printf("FAIL: %s [%s]\n", label, detail); failures++; }
}

#define EIG(reg, k) "NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, " \
                    "u[x, y], {x, y} \\[Element] " reg ", " k "]"

/* ---------------- mesher invariants ---------------- */
static void test_mesher(void) {
    FemGeom g;
    fem_geom_init(&g);
    double L[12] = { 0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2 };
    fem_geom_add_loop(&g, L, NULL, 6, false);
    FemMeshOpts o = { 1e-3, -1, 1000000, true };
    FemMesh m; const char* why;
    bool ok = fem_mesh_build(&g, &o, &m, &why);
    check_true("L-shape meshes", ok, ok ? "ok" : why);
    if (ok) {
        double area = 0, minang = 180, amax = 0;
        bool ccw = true;
        for (int64_t e = 0; e < m.ne; e++) {
            const int64_t* t = m.tri + 6 * e;
            double P[3][2];
            for (int k = 0; k < 3; k++) { P[k][0] = m.x[2 * t[k]]; P[k][1] = m.x[2 * t[k] + 1]; }
            double a = 0.5 * ((P[1][0] - P[0][0]) * (P[2][1] - P[0][1]) - (P[2][0] - P[0][0]) * (P[1][1] - P[0][1]));
            if (a <= 0) ccw = false;
            area += a;
            if (a > amax) amax = a;
            for (int k = 0; k < 3; k++) {
                double ux = P[(k + 1) % 3][0] - P[k][0], uy = P[(k + 1) % 3][1] - P[k][1];
                double vx = P[(k + 2) % 3][0] - P[k][0], vy = P[(k + 2) % 3][1] - P[k][1];
                double ang = acos((ux * vx + uy * vy) / sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy))) * 180 / PI;
                if (ang < minang) minang = ang;
            }
        }
        char buf[128];
        snprintf(buf, sizeof buf, "area %.15g", area);
        check_true("mesh covers the region exactly", fabs(area - 3.0) < 1e-12, buf);
        check_true("every element counter-clockwise", ccw, "orientation");
        snprintf(buf, sizeof buf, "min angle %.2f deg", minang);
        check_true("quality bound (min angle > 20 deg)", minang > 20.0, buf);
        snprintf(buf, sizeof buf, "max area %.3g", amax);
        check_true("MaxCellMeasure respected", amax <= 1e-3, buf);
        /* boundary edges add up to the perimeter */
        double per = 0;
        for (int64_t b = 0; b < m.nb; b++) {
            int64_t i = m.bedge[3 * b], j = m.bedge[3 * b + 1];
            per += hypot(m.x[2 * j] - m.x[2 * i], m.x[2 * j + 1] - m.x[2 * i + 1]);
        }
        snprintf(buf, sizeof buf, "perimeter %.15g", per);
        check_true("boundary edges trace the perimeter", fabs(per - 8.0) < 1e-12, buf);
        fem_mesh_free(&m);
    }
    fem_geom_free(&g);

    /* a disk: midside nodes of the boundary edges lie on the circle */
    fem_geom_init(&g);
    int ai = fem_geom_add_arc(&g, 0, 0, 1);
    double d[32]; int arcs[16];
    for (int k = 0; k < 16; k++) { d[2 * k] = cos(2 * PI * k / 16); d[2 * k + 1] = sin(2 * PI * k / 16); arcs[k] = ai; }
    fem_geom_add_loop(&g, d, arcs, 16, false);
    ok = fem_mesh_build(&g, &o, &m, &why);
    check_true("disk meshes", ok, ok ? "ok" : why);
    if (ok) {
        double worst = 0;
        for (int64_t b = 0; b < m.nb; b++)
            for (int k = 0; k < 3; k++) {
                int64_t i = m.bedge[3 * b + k];
                double r = hypot(m.x[2 * i], m.x[2 * i + 1]);
                if (fabs(r - 1) > worst) worst = fabs(r - 1);
            }
        char buf[64]; snprintf(buf, sizeof buf, "max |r - 1| = %.1e", worst);
        check_true("boundary nodes on the circle (isoparametric)", worst < 1e-14, buf);
        fem_mesh_free(&m);
    }
    fem_geom_free(&g);
}

/* ---------------- sparse Cholesky vs a dense solve ---------------- */
static void test_cholesky(void) {
    /* the 5-point Laplacian on an m x m grid, plus the identity */
    const int mm = 30, n = mm * mm;
    int64_t* rp = malloc(sizeof(int64_t) * (n + 1));
    int64_t* ci = malloc(sizeof(int64_t) * 5 * n);
    double* x = malloc(sizeof(double) * 5 * n);
    double* xy = malloc(sizeof(double) * 2 * n);
    int64_t nz = 0;
    for (int i = 0; i < mm; i++)
        for (int j = 0; j < mm; j++) {
            int r = i * mm + j;
            rp[r] = nz;
            xy[2 * r] = j; xy[2 * r + 1] = i;
            if (i > 0) { ci[nz] = r - mm; x[nz++] = -1; }
            if (j > 0) { ci[nz] = r - 1; x[nz++] = -1; }
            ci[nz] = r; x[nz++] = 5;
            if (j < mm - 1) { ci[nz] = r + 1; x[nz++] = -1; }
            if (i < mm - 1) { ci[nz] = r + mm; x[nz++] = -1; }
        }
    rp[n] = nz;
    const char* why;
    FemChol* F = fem_chol_factor(n, rp, ci, x, xy, &why);
    check_true("nested-dissection Cholesky factors", F != NULL, F ? "ok" : why);
    if (F) {
        double* b = malloc(sizeof(double) * n);
        double* s = malloc(sizeof(double) * n);
        for (int i = 0; i < n; i++) b[i] = sin(0.37 * i) + 1.0;
        fem_chol_solve(F, b, s);
        double res = 0;
        for (int r = 0; r < n; r++) {
            double acc = -b[r];
            for (int64_t p = rp[r]; p < rp[r + 1]; p++) acc += x[p] * s[ci[p]];
            if (fabs(acc) > res) res = fabs(acc);
        }
        char buf[64]; snprintf(buf, sizeof buf, "residual %.1e, nnz(L) %lld", res, (long long)fem_chol_nnz(F));
        check_true("Cholesky solve residual", res < 1e-12, buf);
        /* nested dissection keeps the factor far below the dense n^2/2 */
        check_true("fill stays O(n log n)", fem_chol_nnz(F) < (int64_t)n * 40, buf);
        free(b); free(s);
        fem_chol_free(F);
    }
    /* indefinite: refused */
    for (int r = 0; r < n; r++)
        for (int64_t p = rp[r]; p < rp[r + 1]; p++) if (ci[p] == r) x[p] = -5;
    F = fem_chol_factor(n, rp, ci, x, xy, &why);
    check_true("indefinite matrix refused", F == NULL, F ? "factored" : why);
    fem_chol_free(F);
    free(rp); free(ci); free(x); free(xy);
}

/* ---------------- eigenvalues ---------------- */
static void test_eigenvalues(void) {
    check_rel("square lambda1 = 2 pi^2", "First[" EIG("Rectangle[]", "1") "]", 2 * PI * PI, 2e-6);
    check_rel("square lambda4 = 8 pi^2", "Last[" EIG("Rectangle[]", "4") "]", 8 * PI * PI, 2e-5);
    check_rel("square lambda2 = lambda3 = 5 pi^2 (double)",
              "#[[3]] - #[[2]] &[" EIG("Rectangle[]", "3") "] + 5 Pi^2", 5 * PI * PI, 2e-6);
    check_rel("rectangle 2 x 1 lambda1 = 5 pi^2/4", "First[" EIG("Rectangle[{0, 0}, {2, 1}]", "1") "]",
              1.25 * PI * PI, 2e-6);
    check_rel("disk lambda1 = j01^2", "First[" EIG("Disk[]", "1") "]", 5.783185962946784, 2e-6);
    check_rel("disk lambda2 = j11^2 (double)", "Last[" EIG("Disk[]", "3") "]", 14.681970642123893, 2e-5);
    check_rel("half disk lambda1 = j11^2", "First[" EIG("Disk[{0, 0}, 1, {0, Pi}]", "1") "]",
              14.681970642123893, 2e-5);
    check_rel("annulus 1/2 < r < 1 lambda1", "First[" EIG("Annulus[{0, 0}, {1/2, 1}]", "1") "]",
              39.01328849900286, 2e-5);
    check_rel("equilateral triangle lambda1 = 16 pi^2/3",
              "First[" EIG("Polygon[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}]", "1") "]", 16 * PI * PI / 3, 2e-6);
    check_rel("right isosceles triangle lambda1 = 5 pi^2", "First[" EIG("Triangle[]", "1") "]", 5 * PI * PI, 2e-6);
    /* the L-shape: a re-entrant corner (graded mesh) */
    check_rel("L-shape lambda1 = 9.6397238440...",
              "First[" EIG("Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}]", "1") "]",
              9.6397238440219410527, 1e-5);
    check_rel("L-shape lambda3 = 2 pi^2",
              "Last[" EIG("Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}]", "3") "]",
              2 * PI * PI, 1e-5);
    /* a finer mesh converges (h^4 on a smooth problem) */
    check_rel("square, MaxCellMeasure 1/4000",
              "First[NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, "
              "u[x, y], {x, y} \\[Element] Rectangle[], 1, Method -> {\"PDEDiscretization\" -> "
              "{\"FiniteElement\", {\"MeshOptions\" -> {\"MaxCellMeasure\" -> 1/4000}}}}]]", 2 * PI * PI, 1e-7);
    /* natural (Neumann) conditions everywhere: 0, pi^2, pi^2 */
    check_abs("Neumann square lambda1 = 0",
              "First[NDEigenvalues[-Laplacian[u[x, y], {x, y}], u[x, y], {x, y} \\[Element] Rectangle[], 3]]", 0, 1e-8);
    check_rel("Neumann square lambda2 = pi^2",
              "#[[2]] &[NDEigenvalues[-Laplacian[u[x, y], {x, y}], u[x, y], {x, y} \\[Element] Rectangle[], 3]]",
              PI * PI, 2e-6);
    /* a potential: -Lap u + 10 u shifts every eigenvalue by 10 */
    check_rel("operator with a reaction term",
              "First[NDEigenvalues[{-Laplacian[u[x, y], {x, y}] + 10 u[x, y], DirichletCondition[u[x, y] == 0, True]}, "
              "u[x, y], {x, y} \\[Element] Rectangle[], 1]]", 2 * PI * PI + 10, 2e-6);
}

/* ---------------- eigenfunctions ---------------- */
static void test_eigensystem(void) {
    run("es = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, "
        "u[x, y], {x, y} \\[Element] Rectangle[], 2];");
    check_rel("NDEigensystem value", "es[[1, 1]]", 2 * PI * PI, 2e-6);
    /* unit L2 norm, positive mean: u1 = 2 sin(pi x) sin(pi y) */
    check_abs("eigenfunction normalized and signed", "es[[2, 1]] /. {x -> 0.5, y -> 0.5}", 2.0, 1e-5);
    check_abs("eigenfunction shape", "es[[2, 1]] /. {x -> 0.25, y -> 0.5}", 2.0 * sin(PI / 4), 3e-5);
    check_abs("eigenfunction vanishes on the boundary", "es[[2, 1]] /. {x -> 0.3, y -> 0}", 0.0, 1e-14);
    check_abs("L2 norm 1", "NIntegrate[(es[[2, 1]])^2, {x, 0, 1}, {y, 0, 1}]", 1.0, 1e-5);
    run("es2 = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, "
        "u, {x, y} \\[Element] Disk[], 1];");
    check_abs("u (not u[x, y]) gives the function object", "es2[[2, 1]][0, 0]^2 Pi BesselJ[1, BesselJZero[0, 1]]^2",
              1.0, 1e-4);
}

/* ---------------- NDSolve over a region ---------------- */
static void test_ndsolve(void) {
    check_abs("Poisson on the disk: u = (1 - r^2)/4",
              "NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 1, DirichletCondition[u[x, y] == 0, True]}, u, "
              "{x, y} \\[Element] Disk[]][0.3, 0.4]", (1 - 0.25) / 4, 1e-7);
    check_abs("harmonic x^2 - y^2 reproduced (Laplacian == 0 sign)",
              "NDSolveValue[{Laplacian[u[x, y], {x, y}] == 0, DirichletCondition[u[x, y] == x^2 - y^2, True]}, u, "
              "{x, y} \\[Element] Rectangle[]][0.3, 0.7]", 0.09 - 0.49, 1e-12);
    check_abs("natural condition: u = x (2 - x)",
              "NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 2, DirichletCondition[u[x, y] == 0, x == 0]}, u, "
              "{x, y} \\[Element] Rectangle[]][1, 0.2]", 1.0, 1e-12);
    check_abs("NeumannValue flux: u = x",
              "NDSolveValue[{-Laplacian[u[x, y], {x, y}] == NeumannValue[1, x == 1], "
              "DirichletCondition[u[x, y] == 0, x == 0]}, u, {x, y} \\[Element] Rectangle[]][0.7, 0.4]", 0.7, 1e-12);
    check_abs("Robin NeumannValue[2 - u]: u = x",
              "NDSolveValue[{-Laplacian[u[x, y], {x, y}] == NeumannValue[2 - u[x, y], x == 1], "
              "DirichletCondition[u[x, y] == 0, x == 0]}, u, {x, y} \\[Element] Rectangle[]][0.6, 0.4]", 0.6, 1e-12);
    run("s6 = NDSolveValue[{-Laplacian[u[x, y], {x, y}] + 3 u[x, y] == (2 Pi^2 + 3) Sin[Pi x] Sin[Pi y], "
        "DirichletCondition[u[x, y] == 0, True]}, u, {x, y} \\[Element] Rectangle[]];");
    check_abs("reaction + variable source", "s6[0.3, 0.6]", sin(0.3 * PI) * sin(0.6 * PI), 2e-5);
    check_abs("x-derivative of the solution", "Derivative[1, 0][s6][0.3, 0.6]",
              PI * cos(0.3 * PI) * sin(0.6 * PI), 2e-3);
    check_abs("advection term b . grad u (u = x y (1-x)(1-y) manufactured)",
              "NDSolveValue[{-Laplacian[u[x, y], {x, y}] + D[u[x, y], x] == "
              "2 y (1 - y) + 2 x (1 - x) + (1 - 2 x) y (1 - y), DirichletCondition[u[x, y] == 0, True]}, u, "
              "{x, y} \\[Element] Rectangle[]][0.3, 0.4]", 0.3 * 0.4 * 0.7 * 0.6, 5e-6);
    check_abs("applied form u[x, y]",
              "NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 1, DirichletCondition[u[x, y] == 0, True]}, u[x, y], "
              "{x, y} \\[Element] Disk[]] /. {x -> 0, y -> 0}", 0.25, 1e-7);
    /* the pentagon of SIAM Problem 14: mixed Dirichlet / natural corners graded */
    check_abs("pentagon, mixed corners (graded mesh)",
              "With[{s = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == Exp[-2.46 (x^2 + y^2)] + "
              "NeumannValue[0, y < N[-Cos[Pi/5]] + 10.^-9], DirichletCondition[u[x, y] == 0, y > N[-Cos[Pi/5]] + 10.^-9]}, "
              "u, {x, y} \\[Element] Polygon[Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]], "
              "Method -> {\"FiniteElement\", \"MeshOptions\" -> {\"MaxCellMeasure\" -> 10.^-4}}]}, "
              "NIntegrate[s[x, N[-Cos[Pi/5]]], {x, -Sin[Pi/5], Sin[Pi/5]}]]", 0.085256341, 3e-8);
}

/* ---------------- NIntegrate of mesh functions over their region ---------------- */
static void test_nintegrate(void) {
    /* Saint-Venant torsion: Int u over the equilateral triangle of side 1 = Sqrt[3]/320 */
    run("tri = Triangle[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}]; "
        "st = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 1, DirichletCondition[u[x, y] == 0, True]}, u, "
        "{x, y} \\[Element] tri];");
    check_rel("NIntegrate over the mesh's region (default mesh)",
              "NIntegrate[st[x, y], {x, y} \\[Element] tri]", sqrt(3.0) / 320, 1e-6);
    check_rel("NIntegrate of u^2 over the region",
              "NIntegrate[st[x, y]^2, {x, y} \\[Element] tri] / NIntegrate[st[x, y]^2, {x, y} \\[Element] tri]", 1.0, 1e-15);
    /* rectangle by ranges: torsion of the unit square = 0.0351442... */
    run("sq = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 1, DirichletCondition[u[x, y] == 0, True]}, u, "
        "{x, y} \\[Element] Rectangle[]];");
    check_rel("NIntegrate by ranges over the rectangle", "NIntegrate[sq[x, y], {x, 0, 1}, {y, 0, 1}]",
              0.035144253738653, 5e-6);
    /* eigenfunction L2 norm by the region form */
    run("ed = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], "
        "{x, y} \\[Element] Disk[], 1];");
    check_abs("eigenfunction norm over the disk", "NIntegrate[ed[[2, 1]]^2, {x, y} \\[Element] Disk[]]", 1.0, 1e-12);
    double v;
    check_true("another region stays unevaluated",
               !eval_double("NIntegrate[sq[x, y], {x, y} \\[Element] Disk[]]", &v), "unevaluated");
}

/* ---------------- ContourPlot over a region ---------------- */
static void test_contourplot(void) {
    check_abs("circle contour length (pieces on the level set)",
              "Total[ArcLength /@ Cases[Normal[ContourPlot[x^2 + y^2 == 1/4, {x, y} \\[Element] "
              "Rectangle[{-1, -1}, {1, 1}]]], _Line, Infinity]]", PI, 1e-5);
    check_abs("straight contour", "Total[ArcLength /@ Cases[Normal[ContourPlot[x + y == 1/2, {x, y} \\[Element] "
              "Rectangle[]]], _Line, Infinity]]", sqrt(2.0) / 2, 1e-12);
    /* SIAM 17: the nodal line of the third L-shape eigenfunction, length 2 */
    run("esL = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], "
        "{x, y} \\[Element] Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}], 3];");
    check_abs("L-shape third eigenfunction nodal length = 2",
              "Total[ArcLength /@ Cases[Normal[ContourPlot[Evaluate[esL[[2, 3]] == 0], {x, y} \\[Element] "
              "Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}]]], _Line, Infinity]]", 2.0, 1e-5);
    /* the ground state vanishes only on the boundary: no interior contour */
    check_abs("boundary zero set left out",
              "Length[Cases[Normal[ContourPlot[Evaluate[esL[[2, 1]] == 0], {x, y} \\[Element] "
              "Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}]]], _Line, Infinity]]", 0.0, 0.0);
    /* a disk mode: the nodal circle r = j01/j02 of J0(j02 r) */
    run("esD = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], "
        "{x, y} \\[Element] Disk[], 6];");
    check_rel("disk nodal circle", "Total[ArcLength /@ Cases[Normal[ContourPlot[Evaluate[esD[[2, 6]] == 0], "
              "{x, y} \\[Element] Disk[]]], _Line, Infinity]]", 2 * PI * 2.404825557695773 / 5.520078110286311, 2e-5);
}

/* ---------------- 1-D NIntegrate along a path through a mesh ---------------- */
static void test_nintegrate_path(void) {
    /* u = x (2 - x) exactly (natural condition at x = 1): Int_0^1 u(1, y) dy = 1 */
    run("snat = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 2, DirichletCondition[u[x, y] == 0, x == 0]}, u, "
        "{x, y} \\[Element] Rectangle[]];");
    check_abs("NIntegrate along a side (piecewise exact)", "NIntegrate[snat[1, y], {y, 0, 1}]", 1.0, 1e-12);
    check_abs("NIntegrate along a diagonal", "NIntegrate[snat[t, t], {t, 0, 1}]", 2.0 / 3.0, 1e-12);
    /* a circle of radius 1/2 in the disk: u = (1 - r^2)/4, Int_0^2Pi u dp = 2 Pi 3/16 */
    run("sdisk = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == 1, DirichletCondition[u[x, y] == 0, True]}, u, "
        "{x, y} \\[Element] Disk[]];");
    check_abs("NIntegrate around a circle", "NIntegrate[sdisk[Cos[p]/2, Sin[p]/2], {p, 0, 2 Pi}]",
              2 * PI * 3.0 / 16.0, 1e-7);
    /* WorkingPrecision does not send it to the adaptive rule (which chased the kinks for minutes) */
    check_abs("NIntegrate along a side, WorkingPrecision -> 30",
              "NIntegrate[snat[1, y], {y, 0, 1}, WorkingPrecision -> 30]", 1.0, 1e-12);
}

/* ---------------- messages / declines ---------------- */
static void test_declines(void) {
    double v;
    check_true("nonlinear equation declined",
               !eval_double("NDSolveValue[{-Laplacian[u[x, y], {x, y}] == u[x, y]^2, "
                            "DirichletCondition[u[x, y] == 0, True]}, u, {x, y} \\[Element] Disk[]][0, 0]", &v),
               "unevaluated");
    check_true("too many eigenvalues declined",
               !eval_double("First[NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, "
                            "u[x, y], {x, y} \\[Element] Rectangle[], 10^7]]", &v), "unevaluated");
    check_rel("WorkingPrecision is ignored (with a message), as the reference does",
              "First[" "NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, "
              "u[x, y], {x, y} \\[Element] Rectangle[], 1, WorkingPrecision -> 30]]", 2 * PI * PI, 2e-6);
}

int main(void) {
    test_mute_messages();
    core_init();
    test_mesher();
    test_cholesky();
    test_eigenvalues();
    test_eigensystem();
    test_ndsolve();
    test_nintegrate();
    test_contourplot();
    test_nintegrate_path();
    test_declines();
    if (failures == 0) printf("\nAll FEM tests passed.\n");
    else printf("\n%d FEM test(s) FAILED.\n", failures);
    return failures ? 1 : 0;
}
