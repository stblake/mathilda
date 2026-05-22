#ifndef LUDECOMP_H
#define LUDECOMP_H

#include "expr.h"

/* LUDecomposition[m]
 *     gives the LU decomposition of a square matrix m as a list
 *     {lu, p, c} where:
 *       lu  - combined Doolittle L/U matrix.  The strictly lower
 *             triangle of lu is L (with implicit unit diagonal);
 *             the upper triangle (including the diagonal) is U.
 *       p   - 1-indexed row-permutation vector such that
 *             m[[p]] == l . u, where
 *               l = LowerTriangularize[lu, -1] + IdentityMatrix[n]
 *               u = UpperTriangularize[lu]
 *       c   - L-infinity condition-number estimate for approximate
 *             numerical matrices; the exact Integer 0 for exact or
 *             symbolic inputs.
 *
 * Algorithm.  Doolittle's elimination with partial row pivoting.
 *
 *   - Exact / symbolic matrices: zero-only pivot search (advance to
 *     the next non-zero pivot if the natural choice is provably zero);
 *     output stays exact.
 *
 *   - Machine-precision matrices: LAPACK fast path
 *     (`dgetrf` / `zgetrf` for the factorisation, `dgecon` / `zgecon`
 *     + `dlange` / `zlange` for the condition number).
 *
 *   - Arbitrary-precision MPFR matrices: Householder-free Doolittle
 *     elimination over MPFR arrays at the input precision; condition
 *     number from the explicit inverse via back-substitution.
 *
 * On any fast-kernel failure we fall through to the symbolic core,
 * which itself handles inexact inputs through the standard
 * rationalise -> exact pipeline -> numericalise round-trip used by
 * QRDecomposition / PseudoInverse.
 */
Expr* builtin_ludecomposition(Expr* res);
void  ludecomp_init(void);

#endif /* LUDECOMP_H */
