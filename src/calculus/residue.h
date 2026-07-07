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

/* Register the Residue builtin and its attributes. Called from core_init(). */
void residue_init(void);

#endif /* RESIDUE_H */
