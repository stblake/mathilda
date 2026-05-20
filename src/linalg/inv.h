#ifndef MATINV_H
#define MATINV_H

#include "expr.h"

/* Matrix inversion and pseudo-inversion.
 *
 * `Inverse[m]`            -- inverse of a non-empty square matrix m.
 * `PseudoInverse[m]`      -- Moore-Penrose pseudoinverse of a rectangular
 *                            (or rank-deficient square) matrix m.
 * `PseudoInverse[m,
 *      Tolerance -> t]`   -- on inexact-valued m, singular values below
 *                            tolerance are dropped.  With Tolerance ->
 *                            Automatic the rationalisation precision of
 *                            the input is used.
 *
 * Both builtins share the same source module so the fraction-free
 * Gauss-Jordan helpers can be inlined.  Registration is done by
 * `matinv_init`, which is called from `core_init` in `core.c`.
 */

Expr* builtin_inverse(Expr* res);
Expr* builtin_pseudoinverse(Expr* res);
void  matinv_init(void);

#endif /* MATINV_H */
