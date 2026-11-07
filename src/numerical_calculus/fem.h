/* Mathilda -- the finite element method for planar elliptic problems.
 *
 *   NDEigenvalues[{op, DirichletCondition[...]}, u[x, y], {x, y} in reg, n]
 *   NDEigensystem[{op, DirichletCondition[...]}, u[x, y], {x, y} in reg, n]
 *   NDSolve[{op == f + NeumannValue[g, pred], DirichletCondition[u == h, pred]},
 *           u, {x, y} in reg]
 *
 * for the linear second-order operator
 *
 *   op u = -div(c grad u) + b . grad u + d u,
 *
 * c a constant 2x2 matrix, b and d functions of (x, y), on a region given
 * as Rectangle, Polygon (convex or not), Triangle, Disk (and its sectors) or
 * Annulus. The region is meshed by Delaunay refinement (fem_mesh.c), the
 * problem discretized with quadratic (P2) Lagrange triangles -- isoparametric
 * along curved boundaries -- and the sparse system or generalized eigenproblem
 * solved directly (CHOLMOD / UMFPACK) or by shift-and-invert Krylov iteration
 * (ARPACK or the native Krylov-Schur). The solution and the eigenfunctions are
 * InterpolatingFunction objects over the mesh (fem_interp_apply).
 */
#ifndef MATHILDA_FEM_H
#define MATHILDA_FEM_H

#include "expr.h"
#include <stdbool.h>

void fem_init(void);

/* NDSolve[eqns, u, {x, y} in reg, opts] when the third argument is a region
 * specification: the {{u -> InterpolatingFunction[...]}} answer, or NULL
 * (unevaluated, after a message). `res` is borrowed. */
Expr* fem_ndsolve(Expr* res);
bool fem_is_region_spec(const Expr* e);

/* InterpolatingFunction[domain, NDSolve`FEM`MeshFunction[...], ders...]:
 * is `table` such a mesh function, and evaluate the object at call_args
 * (two coordinates, numbers or equal-length packed vectors). */
bool fem_is_mesh_function(const Expr* table);
Expr* fem_interp_apply(Expr* ifun, Expr** call_args, size_t argc);

/* NIntegrate[f, {x, y} in reg] (or {x, a, b}, {y, c, d} over a rectangle)
 * where the evaluated integrand f involves mesh functions on one mesh over
 * that same region: the integral element by element. NULL otherwise. */
Expr* fem_nintegrate(Expr* body, Expr** specs, size_t nspec);

/* ContourPlot[f, {x, y} in reg, opts]: contours on a mesh of the region
 * (Graphics), or NULL. `res` is the held call. */
Expr* fem_contourplot(Expr* res);

#endif /* MATHILDA_FEM_H */
