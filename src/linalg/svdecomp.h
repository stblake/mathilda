#ifndef SVDECOMP_H
#define SVDECOMP_H

#include "expr.h"

/* SingularValueDecomposition[m]
 *     gives the singular value decomposition of a numerical or
 *     symbolic matrix m as a list {u, sigma, v}, where sigma is a
 *     diagonal matrix and m == u . sigma . ConjugateTranspose[v].
 *     u and v have orthonormal columns.
 *
 * SingularValueDecomposition[{m, a}]
 *     gives the generalized singular value decomposition of m with
 *     respect to a as {{u, ua}, {sigma, sigma_a}, v} such that
 *     m == u . sigma . ConjugateTranspose[v] and
 *     a == ua . sigma_a . ConjugateTranspose[v].
 *
 * SingularValueDecomposition[m, k]
 *     gives the SVD associated with the k largest singular values
 *     (or |k| smallest when k is negative).
 *
 * SingularValueDecomposition[m, UpTo[k]]
 *     gives the SVD for as many of the k largest singular values as
 *     are available (up to MatrixRank[m]).
 *
 * Options:
 *
 *   Tolerance -> t        :  zero out singular values below t.  The
 *                            default is 0 for exact input and
 *                            max(n, p) * 2^(-prec) * sigma[0] for
 *                            approximate input.
 *
 *   TargetStructure -> "Dense"       :  u, sigma, v are all dense
 *                                        List-of-List matrices.  Default.
 *   TargetStructure -> "Structured"  :  sigma is returned as
 *                                        DiagonalMatrix[{sigma_1, ...}];
 *                                        u and v stay dense.
 *
 * Dispatch.  The implementation picks one of three kernels based on
 * the minimum leaf precision of the input:
 *
 *   - Exact / symbolic                -> svd_symbolic_dispatch
 *                                        (eigendecomposition of m^H m;
 *                                        2x2 closed form fast path)
 *   - Inexact, min_bits <= 53          -> svd_machine_dispatch
 *                                        (LAPACK dgesdd / zgesdd, or
 *                                        dggsvd3 / zggsvd3 for the
 *                                        generalized pair form)
 *   - Inexact, min_bits > 53           -> svd_mpfr_dispatch
 *                                        (one-sided Jacobi SVD at
 *                                        the input precision)
 *
 * Each fast kernel returns NULL on any soft failure (USE_LAPACK / USE_MPFR
 * off, non-numeric leaf, LAPACK info != 0, etc.) and the symbolic
 * dispatcher takes over -- which itself rationalises inexact input
 * through the standard pipeline used by QRDecomposition / PseudoInverse
 * / Eigenvalues.  The public contract is invariant under build mode.
 *
 * Memory contract: standard builtin -- this file does NOT call
 * expr_free(res); the evaluator frees the wrapper on a non-NULL return.
 */
Expr* builtin_singularvaluedecomposition(Expr* res);
void  svdecomp_init(void);

#endif /* SVDECOMP_H */
