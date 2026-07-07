#ifndef RESIDUE_H
#define RESIDUE_H

#include "expr.h"

/*
 * Residue[f, {z, z0}]
 *   Gives the residue of f at the isolated singularity z = z0, i.e. the
 *   coefficient of (z - z0)^-1 in the Laurent expansion of f. Computed by
 *   power-series expansion (Series[f, {z, z0, 0}]); a residue can therefore be
 *   found only where f admits a Laurent series at z0. Returns unevaluated at
 *   branch points (fractional-power / Puiseux expansions) and when no series
 *   can be produced. Attribute: Protected. See NResidue for a numerical
 *   alternative that also handles essential singularities.
 */
Expr* builtin_residue(Expr* res);

/*
 * residue_compute(f, z, z0)
 *   The residue of f at z = z0 as a freshly-allocated Expr*, i.e. the coefficient
 *   of (z - z0)^-1 in the Laurent expansion, computed by expanding about z0 with
 *   an EXPANDED denominator (so algebraic pole locations such as z0 = -2+Sqrt[3]
 *   are detected).  f, z, z0 are borrowed (not consumed); z must be a symbol.
 *   Returns NULL at a branch point (Puiseux expansion), when no series can be
 *   produced, or on malformed input.  This is the body behind builtin_residue,
 *   exported so callers such as the contour/residue integrator can sum residues
 *   without re-parsing through the Residue head.
 */
Expr* residue_compute(Expr* f, Expr* z, Expr* z0);

/* Register the Residue builtin and its attributes. Called from core_init(). */
void residue_init(void);

#endif /* RESIDUE_H */
