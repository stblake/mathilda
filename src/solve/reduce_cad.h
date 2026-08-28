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
 * clean nested form, while CLOSED regions merge (Stage B) to a non-strict nested
 * conjunction where the sampling comparison can prove it.
 *
 * Phase 6b (landed): a breakpoint at a non-innermost level may be an irrational
 * algebraic number.  Such a section pins the outer variable to an algebraic
 * value, so the deeper fibre has real-algebraic-number coefficients; it is
 * isolated by rru_algebraic_fiber_roots (reduce_algfiber.c) via iterated-
 * resultant tower projection back to Q + an exact qqbar filter of the conjugate-
 * spurious roots.  The all-rational assignment keeps its unchanged fast path.
 *
 * Hard invariant: any undecidable sign/ordering (the qqbar oracle returning
 * "unknown", FLINT absent, a nullification risking McCallum unsoundness, a
 * fibre whose roots cannot be cleanly isolated/ordered/emitted, a transcendental
 * (non-algebraic) breakpoint, or the algebraic-fibre resource budget overrun)
 * makes the engine return NULL, leaving Reduce unevaluated rather than emitting
 * a wrong formula.
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

/* Phase 7, Case B -- single-free-variable quantifier elimination.  Eliminate the
 * bound variables boundvars[0..nbound-1] from the DNF `F` (a statement in freevar
 * and the bound vars) under the quantifier `quant` (0 = Exists, 1 = ForAll),
 * returning the quantifier-free description in `freevar`, or NULL to decline
 * (a non-polynomial atom, a non-rational free-variable breakpoint, an undecidable
 * sign, or any CAD decline).  freevar is the outermost CAD level so the bound
 * vars are projected out first; the free variable's cells then carry the per-cell
 * Exists/ForAll verdict, emitted by the shared 1-D sign diagram.  F and every
 * Expr* argument are BORROWED; the caller must have stripped nbound==0 already. */
Expr* reduce_cad_qe(const RForm* F, Expr* freevar,
                    Expr** boundvars, int nbound, int quant);

#endif /* REDUCE_CAD_H */
