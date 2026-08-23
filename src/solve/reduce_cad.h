/*
 * reduce_cad.h
 *
 * Multivariate real-inequality engine for `Reduce` (REDUCE_PLAN.md, Phase 6):
 * a Cylindrical Algebraic Decomposition over the reals.
 *
 * McCallum projection eliminates the last listed variable first, so the output
 * reads in the given variable order.  The method:
 *
 *   1. factor the atom polynomials into a distinct-irreducible squarefree basis;
 *   2. project (discriminant, leading coefficient, pairwise resultant, all
 *      factored) once per variable to build the projection stack, eliminating
 *      the variables from last to first;
 *   3. build the 1-D sign diagram of the base projection over the first
 *      variable, then LIFT recursively -- substitute each cell's sample, isolate
 *      the next variable's real roots, and recurse -- to a full cylindrical
 *      decomposition (the innermost dimension merged into a symbolic y-region);
 *   4. evaluate the input formula's truth at one sample per full-dimensional
 *      cell (exact real-algebraic sign oracle), and emit the union of the
 *      satisfying cells with the fibre bounds expressed symbolically as
 *      functions of the outer variables.
 *
 * Two engines share this file: the 2-variable driver (`reduce_cad`, with the
 * boundary-merge that closes a non-strict outer range, e.g. -1<=x<=1) and the
 * n-variable recursive driver (`reduce_cad_nvar`, Phase 6d) for nu>=3.  The
 * recursive engine emits a correct DNF of cells; STRICT inequalities read as a
 * clean nested form, while CLOSED regions currently emit a correct but verbose
 * union of cells (the n-D boundary merge that would close an outer range is a
 * later cosmetic pass).  v1 scope is the rational-fibre regime: a breakpoint at
 * any non-innermost level, given the rational assignment above it, must be
 * rational, else the engine declines (irrational real-algebraic-coefficient
 * fibre isolation is Phase 6b).
 *
 * Hard invariant: any undecidable sign/ordering (the qqbar oracle returning
 * "unknown", FLINT absent, a nullification risking McCallum unsoundness, an
 * irrational non-innermost breakpoint, or a fibre whose roots cannot be cleanly
 * isolated/ordered/emitted) makes the engine return NULL, leaving Reduce
 * unevaluated rather than emitting a wrong formula.
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
