/* gmm.h -- Gaussian mixture fitting by expectation-maximisation.
 *
 * FIRST MODULE IN src/ml/, and the reason it is here rather than being another
 * static function in find_clusters.c is that it is the first clustering kernel with
 * two genuine consumers: `FindClusters[..., Method -> "GaussianMixture"]` and the
 * `LearnDistribution` learner of the same name. find_clusters.c currently holds 55
 * static functions and exports exactly one symbol, so nothing in it can be reached
 * from a second builtin; this module is the point at which that stops being free.
 *
 * The API is therefore deliberately NOT expressed in terms of FindClusters' FcData.
 * It takes a row-major n x dim buffer of machine doubles and returns fitted
 * parameters plus a hard assignment, so a distribution learner that has no notion
 * of a spanning tree or a cluster count can call it unchanged.
 */
#ifndef ML_GMM_H
#define ML_GMM_H

#include <stddef.h>
#include <stdbool.h>

/* A fitted mixture. Owns its buffers; release with ml_gmm_free.
 *
 * `chol` holds the lower Cholesky factor of each component's covariance and
 * `logdet` its log-determinant, both derived from `cov` and cached because the
 * E-step needs them once per point per component per iteration. `diagonal[j]` is
 * true when component j fell back to a diagonal covariance because the full matrix
 * would not factorise. */
typedef struct {
    size_t  k;          /* component count */
    size_t  dim;
    double* w;          /* k       -- mixing weights, sum to 1 */
    double* mu;         /* k x dim -- component means, row-major */
    double* cov;        /* k x dim x dim -- covariances, row-major per component */
    double* chol;       /* k x dim x dim -- lower Cholesky factor of cov */
    double* logdet;     /* k       -- log|cov| */
    bool*   diagonal;   /* k       -- true if this component fell back to diagonal */
    double  loglik;     /* log-likelihood of the data under the fit */
} MlGmm;

void ml_gmm_free(MlGmm* g);

/* Fit `k` components to `pts` (row-major, n x dim) by EM.
 *
 * `var_floor` is a variance ridge added to every covariance diagonal, in squared
 * data units. It is load-bearing rather than defensive: a Gaussian mixture's
 * likelihood is unbounded above, since a component collapsing onto a single point
 * drives its density to infinity, so without a floor at the resolvable scale a
 * model-selection search buys arbitrarily many near-singular spikes and reports one
 * component per point. Pass the squared typical inter-point spacing.
 *
 * `assign` (n entries, may be NULL) receives each point's highest-responsibility
 * component, numbered 0..k-1 by component index.
 *
 * Returns NULL if allocation fails or the fit degenerates. A non-NULL result always
 * carries a finite `loglik`. */
MlGmm* ml_gmm_fit(const double* pts, size_t n, size_t dim, size_t k,
                  double var_floor, size_t* assign);

/* Free parameter count for BIC, full covariance:
 *   (k - 1) weights + k*dim means + k*dim*(dim+1)/2 covariance entries.
 * The 1-D case collapses to the familiar 3k - 1. Exposed because model selection
 * belongs to the caller -- FindClusters wants a cluster count, a distribution
 * learner wants a density, and they need not agree on the criterion. */
double ml_gmm_param_count(size_t k, size_t dim);

/* Log-density of one point under a fitted mixture. The distribution-learner entry
 * point: it is what makes the fitted object a density rather than a partition. */
double ml_gmm_logpdf(const MlGmm* g, const double* x);

#endif /* ML_GMM_H */
