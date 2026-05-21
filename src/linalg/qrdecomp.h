#ifndef QRDECOMP_H
#define QRDECOMP_H

#include "expr.h"

/* QRDecomposition[m]
 *     gives the QR decomposition of m as a list {q, r}, where q is
 *     a row-orthonormal (or row-unitary in the complex case) matrix
 *     and r is upper triangular.  The original matrix m is equal
 *     to ConjugateTranspose[q] . r.
 *
 * QRDecomposition[m, Pivoting -> True]
 *     gives a list {q, r, p} where p is a permutation matrix such
 *     that m . p == ConjugateTranspose[q] . r.  Pivoting is most
 *     useful for numerical matrices where it puts the diagonal of r
 *     in order of decreasing magnitude.
 *
 * Result shape:
 *     Let n = Length[m] (rows), p = Length[First[m]] (cols), and
 *     r_rank = MatrixRank[m] ("thin" QR factorisation).  Then q is
 *     r_rank x n and r is r_rank x p; both have MatrixRank[m] rows.
 *
 * Algorithm:
 *     Modified Gram-Schmidt on the columns of m.  For inexact
 *     matrices (Real or MPFR leaves) we rationalise to exact
 *     arithmetic at the input precision, run the exact pipeline,
 *     then numericalise the result back to the input precision -
 *     the same pattern used by PseudoInverse / Eigenvalues.  This
 *     gives us a single code path that handles symbolic, exact
 *     rational, complex, machine-precision Real, and arbitrary-
 *     precision MPFR inputs.
 */
Expr* builtin_qrdecomposition(Expr* res);
void  qrdecomp_init(void);

#endif /* QRDECOMP_H */
