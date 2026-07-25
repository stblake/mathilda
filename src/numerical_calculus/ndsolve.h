/* Mathilda — NDSolve: numeric solver for ODE initial-value problems.
 *
 * Registers the NDSolve builtin and one NDSolve`Method builtin per integrator.
 * The numerical machinery lives in ndsolve_common.c and the per-method modules
 * (ndsolve_euler.c, ndsolve_rk.c, ndsolve_implicit.c, ndsolve_adams.c). */
#ifndef NDSOLVE_H
#define NDSOLVE_H

#include "../expr.h"

/* NDSolve[eqns, u, {x, xmin, xmax}] — numeric ODE solver returning a list of
 * rules {{u -> InterpolatingFunction[...]}}. */
Expr* builtin_ndsolve(Expr* res);

/* Register NDSolve and the NDSolve`Method submethod builtins. */
void ndsolve_init(void);

#endif /* NDSOLVE_H */
