#ifndef MATEIGEN_H
#define MATEIGEN_H

#include <stddef.h>
#include <stdint.h>
#include "expr.h"

/* Eigenvalues / Eigenvectors of a square matrix.
 *
 * Current implementation (Phase 0 extraction from linalg.c, unchanged):
 *   - Compute characteristic polynomial det(lambda*I - m) via
 *     Faddeev-Leverrier-Souriau for the ordinary case, or Laplace
 *     expansion of det(m - lambda*a) for the generalised case.
 *   - Route the resulting univariate polynomial through the public
 *     Solve builtin so its rationalise -> solve -> numericalize
 *     pipeline handles approximate (Real / MPFR) input automatically.
 *
 * Forthcoming numeric Method dispatch (Direct / Arnoldi / Banded /
 * FEAST) lands in subsequent phases; see plan
 * lovely-roaming-diffie.md.
 */
Expr* builtin_eigenvalues(Expr* res);
Expr* builtin_eigenvectors(Expr* res);

/* Registers Eigenvalues and Eigenvectors builtins and their attributes.
 * Must be called from core_init() (after symtab_init()). */
void  mateigen_init(void);

#endif /* MATEIGEN_H */
