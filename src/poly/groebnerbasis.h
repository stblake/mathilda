/* groebnerbasis.h
 *
 * Mathilda's `GroebnerBasis` builtin.  See groebner.h for the underlying
 * Buchberger machinery; this file owns only the user-facing entry point
 * and the symbol-table wiring.
 */
#ifndef GROEBNERBASIS_H
#define GROEBNERBASIS_H

#include "expr.h"

/* `GroebnerBasis[polys, vars]` or `GroebnerBasis[polys, mainVars, elimVars]`
 * with optional trailing Rule[] options.  Returns a fresh `List[...]` of
 * polynomial Expr* on success, or NULL when the input does not satisfy
 * the polynomial-over-Q contract (non-polynomial substructure, bad
 * arity, ...).  Takes ownership of `res` per the standard builtin
 * contract. */
Expr* builtin_groebner_basis(Expr* res);

/* Register the builtin, attach the Protected attribute, and install
 * the docstring placeholder. */
void groebner_init(void);

#endif /* GROEBNERBASIS_H */
