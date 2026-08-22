/* pca.h -- column statistics, standardisation and principal components.
 *
 * Second module in src/ml/, and buffer-level for the same reason gmm.h is: every
 * consumer here has more than one caller in prospect. Column means and standard
 * deviations are wanted by Standardize, by PrincipalComponents, by a future
 * Classify's feature normalisation and by LearnDistribution's Multinormal; the
 * symmetric eigendecomposition is wanted by PCA, by MultidimensionalScaling and by
 * LatentSemanticAnalysis. Writing them against a row-major n x dim double buffer
 * rather than against an Expr keeps all of that reachable.
 *
 * Convention throughout: `x` is n observations by dim variables, row-major, so
 * x[i * dim + j] is observation i's variable j. That matches Wolfram's data
 * orientation for PrincipalComponents (rows are examples, columns are features).
 */
#ifndef ML_PCA_H
#define ML_PCA_H

#include <stddef.h>
#include <stdbool.h>

/* Column means into `mean` (dim entries). */
void ml_column_mean(const double* x, size_t n, size_t dim, double* mean);

/* Column SAMPLE standard deviations (divisor n - 1) into `sd`.
 *
 * n - 1 rather than n because that is what Wolfram's Standardize and
 * StandardDeviation use, and the two must agree or `Standardize` would disagree with
 * `(x - Mean[x]) / StandardDeviation[x]` written out by hand. A constant column gets
 * sd 0, which callers must treat as "do not rescale" rather than dividing. */
void ml_column_sd(const double* x, size_t n, size_t dim, const double* mean,
                  double* sd);

/* Standardise `x` in place: subtract the column mean, and when `rescale` also divide
 * by the column sample standard deviation.
 *
 * A constant column is left at exactly zero instead of becoming NaN. That is a real
 * decision, not a guard: zero variance carries no information, so the honest
 * standardised value is "no deviation from the mean", and propagating NaN would
 * poison every downstream reduction over the row. */
void ml_standardize(double* x, size_t n, size_t dim, bool rescale);

/* Symmetric eigendecomposition with eigenvalues in DESCENDING order.
 *
 * `a` is dim x dim row-major and is read, not modified. `eval` receives dim
 * eigenvalues descending; `evec` receives eigenvector j in ROW j, so projecting a
 * centred observation is a plain row-by-matrix-transpose contraction.
 *
 * Uses LAPACK's dsyev when it is linked and falls back to cyclic Jacobi otherwise --
 * no #ifdef at the call site, because the no-LAPACK build links a stub that returns
 * a nonzero status, so "did it work" is the same question either way.
 *
 * SIGNS ARE CANONICALISED: an eigenvector is only defined up to sign, and LAPACK and
 * Jacobi do not agree on which they return. Each row is flipped so that its
 * largest-magnitude component is positive, which makes the output reproducible
 * across builds and testable at all. */
bool ml_sym_eigen_desc(const double* a, size_t dim, double* eval, double* evec);

/* Principal components of `x`.
 *
 * `correlation` selects the correlation matrix over the covariance matrix, i.e.
 * standardises each column to unit variance first -- the right choice when the
 * variables have incommensurable units, since otherwise the component with the
 * largest raw scale dominates for no statistical reason.
 *
 * `out` (n x dim, may be NULL) receives the observations in principal-component
 * coordinates; `eval` (dim, may be NULL) the component variances descending;
 * `evec` (dim x dim, may be NULL) the loadings, component j in row j. */
bool ml_pca(const double* x, size_t n, size_t dim, bool correlation,
            double* out, double* eval, double* evec);

/* ------------------------------------------------------------------------- */
/* Dimensionality reduction                                                   */
/* ------------------------------------------------------------------------- */

/* The three reducers implemented here are one algorithm with three ways of forming
 * the symmetric matrix to decompose, which is why they share ml_sym_eigen_desc
 * rather than each carrying its own linear algebra:
 *
 *   PCA  -- centre the columns, decompose the covariance, project onto the leading
 *           eigenvectors.
 *   LSA  -- do NOT centre; decompose the Gram matrix X'X instead. That is a
 *           truncated SVD, and skipping the centring is the whole difference: a
 *           term-document matrix is sparse and non-negative, and centring it
 *           destroys both properties along with the meaning of a zero entry.
 *   MDS  -- classical (Torgerson) scaling: double-centre the squared distance
 *           matrix, decompose that, and scale each axis by the square root of its
 *           eigenvalue. Its eigenproblem is n x n rather than dim x dim, because it
 *           works from distances between OBSERVATIONS.
 */
typedef enum {
    ML_REDUCE_PCA = 0,
    ML_REDUCE_LSA,
    ML_REDUCE_MDS
} MlReduceMethod;

/* Reduce `x` (row-major n x dim) to `target` dimensions into `out` (n x target).
 *
 * `target` is clamped to dim for PCA and LSA, whose output axes are directions in
 * feature space and so cannot outnumber the features. MDS can in principle produce
 * up to n - 1 axes, but a non-positive eigenvalue means the requested dimension
 * carries no real structure, and those columns come back as zero rather than as the
 * square root of a negative number.
 *
 * Returns false on allocation failure, a degenerate decomposition, or an n above the
 * MDS ceiling (its matrix is n x n). */
bool ml_reduce(const double* x, size_t n, size_t dim, size_t target,
               MlReduceMethod method, double* out);

void ml_init(void);

#endif /* ML_PCA_H */
