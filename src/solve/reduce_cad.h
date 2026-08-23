/*
 * reduce_cad.h
 *
 * Multivariate real-inequality engine for `Reduce` (REDUCE_PLAN.md, Phase 6):
 * a Cylindrical Algebraic Decomposition over the reals.
 *
 * This pass implements the 2-variable case (McCallum projection eliminating the
 * last listed variable first, so the output reads in the given variable order),
 * architected so the projection/lift core extends to n variables later.  The
 * method:
 *
 *   1. factor the atom polynomials into a distinct-irreducible squarefree basis;
 *   2. project (discriminant, leading coefficient, pairwise resultant of the
 *      basis, all factored) to eliminate the last variable;
 *   3. build the 1-D sign diagram of the projection over the first variable
 *      (the base), then LIFT each base cell -- substitute its sample and isolate
 *      the fibre's real roots -- to a full cylindrical decomposition;
 *   4. evaluate the input formula's truth at one sample per full-dimensional
 *      cell (exact real-algebraic sign oracle), and emit the union of the
 *      satisfying cells with the fibre bounds expressed symbolically as
 *      functions of the outer variable.
 *
 * Hard invariant: any undecidable sign/ordering (the qqbar oracle returning
 * "unknown", FLINT absent, a nullification risking McCallum unsoundness, or a
 * fibre whose roots cannot be cleanly isolated/ordered/emitted) makes the engine
 * return NULL, leaving Reduce unevaluated rather than emitting a wrong formula.
 */
#ifndef REDUCE_CAD_H
#define REDUCE_CAD_H

#include "expr.h"
#include "reduce_form.h"

/* Solve the DNF formula `F` over the reals in the variables `vars[0..nv-1]` by
 * cylindrical algebraic decomposition.  Returns a freshly-owned Expr (True /
 * False / a logical combination of relations describing the whole solution set),
 * or NULL to decline (a non-polynomial atom, an undecidable sign/ordering, an
 * un-emittable fibre, or nv outside the supported range) -- in which case Reduce
 * stays unevaluated. */
Expr* reduce_cad(const RForm* F, Expr** vars, int nv);

#endif /* REDUCE_CAD_H */
