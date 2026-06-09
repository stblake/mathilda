/* Mathilda -- the PolyGamma function family.
 *
 *   PolyGamma[z]      digamma psi(z)            (canonicalises to PolyGamma[0, z])
 *   PolyGamma[n, z]   n-th polygamma psi^(n)(z) = d^n/dz^n psi(z)
 *
 * See polygamma.c for the layered evaluation strategy.
 */
#ifndef POLYGAMMA_H
#define POLYGAMMA_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the standard contract:
 * returns a new Expr* on success, or NULL to leave the call symbolic. */
Expr* builtin_polygamma(Expr* res);

/* Registers PolyGamma (and the inert LogGamma symbol) in the symbol table. */
void polygamma_init(void);

#endif /* POLYGAMMA_H */
