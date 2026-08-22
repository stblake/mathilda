/* gmm.c -- Gaussian mixture fitting by expectation-maximisation, any dimension.
 *
 * See gmm.h for why this lives in src/ml/ rather than inside find_clusters.c.
 *
 * Three things here are not mechanical generalisations of the one-dimensional
 * kernel this replaces for vector input, and each is commented where it happens:
 *
 *   - the covariance floor becomes a RIDGE on the diagonal rather than a clamp,
 *     because "clamp the variance up to a minimum" has no direct matrix reading;
 *   - the E-step runs in LOG space through log-sum-exp, because a product of dim
 *     Gaussian factors underflows to zero for a point a few standard deviations out
 *     once dim is more than a handful, and the linear-space kernel's response to a
 *     zero total is to hand the point uniform responsibilities -- i.e. to silently
 *     stop distinguishing exactly the points that most need it;
 *   - a component whose covariance will not factorise falls back to DIAGONAL rather
 *     than failing the fit, which keeps a rank-deficient component (fewer points
 *     than dimensions, or collinear ones) from taking the whole model with it.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mlutil.h"   /* ml_chol, ml_mahalanobis -- shared with dist.c */
#include "gmm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ML_GMM_MAX_ITER 200

void ml_gmm_free(MlGmm* g) {
    if (!g) return;
    free(g->w); free(g->mu); free(g->cov); free(g->chol);
    free(g->logdet); free(g->diagonal);
    free(g);
}

double ml_gmm_param_count(size_t k, size_t dim) {
    double kk = (double)k, dd = (double)dim;
    return (kk - 1.0) + kk * dd + kk * dd * (dd + 1.0) / 2.0;
}

/* The Cholesky factorisation and Mahalanobis distance that used to live here moved to
 * src/ml/mlutil.c when LearnDistribution's Multinormal became their second real
 * consumer -- the same "extract at the second consumer" rule that promoted the list
 * builders, and that declined to extract find_clusters.c's metric layer. */

/* Refresh chol/logdet from cov for component j, falling back to a diagonal
 * covariance if the full matrix will not factorise. Always succeeds: the diagonal
 * carries the ridge, so it is positive definite by construction. */
static void ml_refresh(MlGmm* g, size_t j, double var_floor) {
    size_t dim = g->dim;
    double* cv = g->cov + j * dim * dim;
    double* l  = g->chol + j * dim * dim;
    if (ml_chol(cv, l, dim)) {
        g->diagonal[j] = false;
    } else {
        /* Keep the variances, drop the covariances. A component with fewer points
         * than dimensions, or with collinear points, has a singular scatter matrix;
         * that is a statement about the data, not an error, and a diagonal
         * component still models it. */
        for (size_t a = 0; a < dim; a++)
            for (size_t b = 0; b < dim; b++)
                if (a != b) cv[a * dim + b] = 0.0;
        for (size_t a = 0; a < dim; a++)
            if (!(cv[a * dim + a] > var_floor)) cv[a * dim + a] = var_floor;
        (void)ml_chol(cv, l, dim);
        g->diagonal[j] = true;
    }
    double ld = 0.0;
    for (size_t a = 0; a < dim; a++) ld += log(l[a * dim + a]);
    g->logdet[j] = 2.0 * ld;
}

double ml_gmm_logpdf(const MlGmm* g, const double* x) {
    size_t dim = g->dim;
    double* y = malloc(sizeof(double) * (dim ? dim : 1));
    if (!y) return -INFINITY;
    double best = -INFINITY;
    double* lp = malloc(sizeof(double) * (g->k ? g->k : 1));
    if (!lp) { free(y); return -INFINITY; }
    for (size_t j = 0; j < g->k; j++) {
        if (!(g->w[j] > 0.0)) { lp[j] = -INFINITY; continue; }
        double q = ml_mahalanobis(g->chol + j * dim * dim, g->mu + j * dim, x, dim, y);
        lp[j] = log(g->w[j]) - 0.5 * ((double)dim * log(2.0 * M_PI) + g->logdet[j] + q);
        if (lp[j] > best) best = lp[j];
    }
    double s = 0.0;
    if (best > -INFINITY)
        for (size_t j = 0; j < g->k; j++)
            if (lp[j] > -INFINITY) s += exp(lp[j] - best);
    free(lp); free(y);
    return (s > 0.0) ? best + log(s) : -INFINITY;
}

MlGmm* ml_gmm_fit(const double* pts, size_t n, size_t dim, size_t k,
                  double var_floor, size_t* assign) {
    if (!pts || n == 0 || dim == 0 || k == 0 || k > n) return NULL;
    if (!(var_floor > 0.0)) var_floor = 1e-300;

    MlGmm* g = calloc(1, sizeof(MlGmm));
    if (!g) return NULL;
    g->k = k; g->dim = dim;
    g->w        = malloc(sizeof(double) * k);
    g->mu       = malloc(sizeof(double) * k * dim);
    g->cov      = malloc(sizeof(double) * k * dim * dim);
    g->chol     = malloc(sizeof(double) * k * dim * dim);
    g->logdet   = malloc(sizeof(double) * k);
    g->diagonal = malloc(sizeof(bool) * k);
    double* lp  = malloc(sizeof(double) * n * k);
    double* y   = malloc(sizeof(double) * dim);
    if (!g->w || !g->mu || !g->cov || !g->chol || !g->logdet || !g->diagonal
        || !lp || !y) {
        free(lp); free(y); ml_gmm_free(g); return NULL;
    }

    /* Deterministic farthest-first initialisation, the same choice fc_lloyd_ndim
     * makes and for the same reasons: no RandomVariate exists in the tree yet, and a
     * mean seeded from the data set rather than from its first element does not
     * depend on the order the points were given in. Each covariance starts
     * isotropic at the floor. */
    {
        double* nd = malloc(sizeof(double) * n);
        if (!nd) { free(lp); free(y); ml_gmm_free(g); return NULL; }
        double* c0 = calloc(dim, sizeof(double));
        if (!c0) { free(nd); free(lp); free(y); ml_gmm_free(g); return NULL; }
        for (size_t i = 0; i < n; i++)
            for (size_t a = 0; a < dim; a++) c0[a] += pts[i * dim + a];
        for (size_t a = 0; a < dim; a++) c0[a] /= (double)n;
        size_t seed = 0; double bd = -1.0;
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t a = 0; a < dim; a++) {
                double t = pts[i * dim + a] - c0[a]; s += t * t;
            }
            if (bd < 0.0 || s < bd) { bd = s; seed = i; }
        }
        free(c0);
        memcpy(g->mu, pts + seed * dim, sizeof(double) * dim);
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t a = 0; a < dim; a++) {
                double t = pts[i * dim + a] - g->mu[a]; s += t * t;
            }
            nd[i] = s;
        }
        for (size_t j = 1; j < k; j++) {
            size_t pick = 0; double far = -1.0;
            for (size_t i = 0; i < n; i++) if (nd[i] > far) { far = nd[i]; pick = i; }
            memcpy(g->mu + j * dim, pts + pick * dim, sizeof(double) * dim);
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t a = 0; a < dim; a++) {
                    double t = pts[i * dim + a] - g->mu[j * dim + a]; s += t * t;
                }
                if (s < nd[i]) nd[i] = s;
            }
        }
        free(nd);
    }
    for (size_t j = 0; j < k; j++) {
        g->w[j] = 1.0 / (double)k;
        for (size_t a = 0; a < dim; a++)
            for (size_t b = 0; b < dim; b++)
                g->cov[j * dim * dim + a * dim + b] = (a == b) ? var_floor : 0.0;
        ml_refresh(g, j, var_floor);
    }

    double loglik = -INFINITY, prev = -INFINITY;
    for (int it = 0; it < ML_GMM_MAX_ITER; it++) {
        /* ---- E step, in log space ----
         * A component's density is a product of dim factors, so in linear space a
         * point a few standard deviations from every mean underflows the total to
         * zero. The 1-D kernel answers a zero total by assigning uniform
         * responsibilities, which is survivable on a line and actively wrong here:
         * the outlying points are the informative ones. log-sum-exp keeps them. */
        loglik = 0.0;
        for (size_t i = 0; i < n; i++) {
            double best = -INFINITY;
            for (size_t j = 0; j < k; j++) {
                if (!(g->w[j] > 0.0)) { lp[i * k + j] = -INFINITY; continue; }
                double q = ml_mahalanobis(g->chol + j * dim * dim, g->mu + j * dim,
                                          pts + i * dim, dim, y);
                lp[i * k + j] = log(g->w[j])
                              - 0.5 * ((double)dim * log(2.0 * M_PI)
                                       + g->logdet[j] + q);
                if (lp[i * k + j] > best) best = lp[i * k + j];
            }
            double s = 0.0;
            if (best > -INFINITY)
                for (size_t j = 0; j < k; j++)
                    if (lp[i * k + j] > -INFINITY) s += exp(lp[i * k + j] - best);
            if (!(s > 0.0)) {           /* every component vanished for this point */
                for (size_t j = 0; j < k; j++) lp[i * k + j] = -log((double)k);
                loglik += -1e300 / (double)n;   /* finite, and hopeless: EM will move */
                continue;
            }
            double lse = best + log(s);
            loglik += lse;
            for (size_t j = 0; j < k; j++) lp[i * k + j] -= lse;   /* log-resp */
        }

        /* ---- M step ---- */
        for (size_t j = 0; j < k; j++) {
            double sw = 0.0;
            for (size_t i = 0; i < n; i++) sw += exp(lp[i * k + j]);
            if (!(sw > 1e-12)) {        /* component died; retire it cleanly */
                g->w[j] = 0.0;
                continue;
            }
            double* m = g->mu + j * dim;
            for (size_t a = 0; a < dim; a++) m[a] = 0.0;
            for (size_t i = 0; i < n; i++) {
                double r = exp(lp[i * k + j]);
                for (size_t a = 0; a < dim; a++) m[a] += r * pts[i * dim + a];
            }
            for (size_t a = 0; a < dim; a++) m[a] /= sw;

            double* cv = g->cov + j * dim * dim;
            for (size_t a = 0; a < dim * dim; a++) cv[a] = 0.0;
            for (size_t i = 0; i < n; i++) {
                double r = exp(lp[i * k + j]);
                for (size_t a = 0; a < dim; a++) {
                    double da = pts[i * dim + a] - m[a];
                    for (size_t b = 0; b <= a; b++)
                        cv[a * dim + b] += r * da * (pts[i * dim + b] - m[b]);
                }
            }
            for (size_t a = 0; a < dim; a++)
                for (size_t b = 0; b <= a; b++) {
                    cv[a * dim + b] /= sw;
                    cv[b * dim + a] = cv[a * dim + b];   /* symmetrise */
                }
            /* THE FLOOR, as a ridge rather than a clamp.
             *
             * The 1-D kernel writes `if (var < floor) var = floor`, and a matrix has
             * no equally direct reading of that: the honest translation is "no
             * direction may have variance below the floor", which is a clamp on the
             * eigenvalues and costs an eigendecomposition per component per
             * iteration. Adding floor*I instead raises EVERY direction's variance by
             * exactly the floor, which enforces the same lower bound (it can only
             * overshoot, never undershoot) for one add, and is what makes the
             * factorisation below succeed on a component with fewer points than
             * dimensions. The two agree when the covariance is already isotropic. */
            for (size_t a = 0; a < dim; a++) cv[a * dim + a] += var_floor;

            g->w[j] = sw / (double)n;
            ml_refresh(g, j, var_floor);
        }

        if (it > 0 && fabs(loglik - prev) <= 1e-10 * (fabs(loglik) + 1.0)) break;
        prev = loglik;
    }

    g->loglik = loglik;

    if (assign) {
        for (size_t i = 0; i < n; i++) {
            size_t best = 0; double bp = -INFINITY;
            for (size_t j = 0; j < k; j++)
                if (lp[i * k + j] > bp) { bp = lp[i * k + j]; best = j; }
            assign[i] = best;
        }
    }

    free(lp); free(y);
    if (!(g->loglik > -INFINITY) || g->loglik != g->loglik) { ml_gmm_free(g); return NULL; }
    return g;
}
