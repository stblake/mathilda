/* rref_numeric.h -- tolerance-aware Gauss-Jordan on inexact matrices.
 *
 * One kernel behind RowReduce, NullSpace, LinearSolve and Inverse of a
 * machine-precision or arbitrary-precision (MPFR) matrix, real or complex,
 * List or packed.  The exact paths (fraction-free Bareiss, FLINT) are the
 * wrong tool for an inexact matrix: they decide "is this pivot zero" exactly,
 * so a pivot that is zero up to rounding is used as a divisor -- RowReduce of a
 * singular 30-digit matrix answered the identity (and Overflow[] entries once
 * the 10^-30 residue was divided through), and a 54 x 54 complex MPFR solve
 * ground through boxed symbolic arithmetic for tens of seconds.
 *
 * Here a candidate pivot counts as zero when its magnitude is at most
 *     tol = f max(rows, cols) eps max|a_ij|,   eps = 2^(1 - bits)
 * (bits = 53 for machine numbers, the input's minimum precision for MPFR).
 * f = 1 for a machine matrix, the reference's own rule (probed: diag(1, 1, e)
 * has rank 3 exactly when e > 3 eps); f = 64 at arbitrary precision, where
 * the reference tracks the error instead -- a few dozen rounding errors of the
 * largest entry, the residue an exactly singular matrix leaves.  Partial
 * pivoting (largest magnitude in the column) keeps the elimination stable.
 * A machine RowReduce that used a pivot below 10^6 eps max|a_ij| warns
 * RowReduce::luc, as the reference does.
 *
 * Result conventions follow the reference: RowReduce puts exact 1 at each
 * pivot, exact 0 below the pivots and in the zero rows; above a pivot the
 * eliminated entry is an inexact zero for a machine matrix and exact 0 for an
 * arbitrary-precision one.  Every other entry carries the input precision.
 *
 * All entry points borrow their arguments and return a fresh Expr*, or NULL
 * to decline (not a numeric matrix, wrong shape, no MPFR support). */
#ifndef MATHILDA_LINALG_RREF_NUMERIC_H
#define MATHILDA_LINALG_RREF_NUMERIC_H

#include "expr.h"
#include <stdbool.h>

/* Working precision of an inexact matrix (and optional right-hand side): 53
 * when any inexact leaf is machine precision (contagion), else the minimum
 * MPFR precision; 0 when the operands are exact or not all numeric. */
long rref_numeric_bits(const Expr* m, const Expr* b);

/* RowReduce[m]. */
Expr* rref_numeric_rowreduce(const Expr* m, long bits);

/* MatrixRank[m]: the number of pivots the elimination above accepts, or -1
 * to decline.  At arbitrary precision this is the rank to working precision:
 * N[HilbertMatrix[40], 100] has rank 40 (its condition number, ~1e58, is far
 * beyond 2^53 but far below 10^100). */
int rref_numeric_rank(const Expr* m, long bits);

/* NullSpace[m]: an orthonormal basis (the reference's normalisation), {} for
 * a matrix of full column rank. */
Expr* rref_numeric_nullspace(const Expr* m, long bits);

/* LinearSolve[m, b] (b a vector or a matrix): the solution with every free
 * variable 0.  NULL with *nosol set when the system is inconsistent. */
Expr* rref_numeric_linearsolve(const Expr* m, const Expr* b, long bits, bool* nosol);

/* Inverse[m] of a square matrix.  NULL with *singular set when it is. */
Expr* rref_numeric_inverse(const Expr* m, long bits, bool* singular);

#endif /* MATHILDA_LINALG_RREF_NUMERIC_H */
