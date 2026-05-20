#ifndef MATRANK_H
#define MATRANK_H

#include "expr.h"

/* MatrixRank[m]                                  -- gives the rank of m.
 * MatrixRank[m, Method -> "<name>"]              -- explicit RREF method.
 * MatrixRank[m, Tolerance -> t]                  -- numerical tolerance.
 *
 * The rank of a matrix is the number of linearly independent rows
 * (equivalently, of linearly independent columns).
 *
 * MatrixRank works on both numerical and symbolic matrices and on
 * square or rectangular matrices.
 *
 * Method names accepted (same grammar as NullSpace / RowReduce /
 * LinearSolve / Inverse via matsol_parse_method_option):
 *   "Automatic"                 -- alias for "DivisionFreeRowReduction"
 *   "DivisionFreeRowReduction"  -- Bareiss-like fraction-free RREF
 *   "OneStepRowReduction"       -- classical Gauss-Jordan with division
 *   "CofactorExpansion"         -- identity-if-invertible (falls back to
 *                                  DivisionFreeRowReduction on singular /
 *                                  rectangular m).
 *
 * Tolerance accepted forms:
 *   Tolerance -> Automatic       -- machine-precision default for inexact m
 *   Tolerance -> 0               -- structural zero only
 *   Tolerance -> <non-neg num>   -- treat |entry| <= t as zero.
 *
 * For exact (integer / rational / symbolic) input the Method option
 * routes through RowReduce and the rank is the number of non-zero
 * rows.  For numerical input (or when Tolerance is supplied) MatrixRank
 * runs a partial-pivot Gaussian elimination over `double _Complex` with
 * tolerance-aware pivot selection.
 */
Expr* builtin_matrixrank(Expr* res);
void  matrank_init(void);

#endif /* MATRANK_H */
