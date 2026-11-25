/* dist.h -- distribution objects, RandomVariate, and PDF.
 *
 * The substrate family 4 needs. Two things live here that nothing in the tree had:
 *
 *   - A GAUSSIAN DEVIATE. There was none anywhere -- `RandomReal` is uniform and the
 *     xoshiro stream underneath it produces uniforms. Box-Muller on top of
 *     random_uniform_01 is the standard construction and needs no tables.
 *   - DISTRIBUTION OBJECTS as first-class expressions, so that
 *     `NormalDistribution[0, 1]` can be passed around, given to RandomVariate, and
 *     asked for a PDF.
 *
 * Distribution objects are NOT the trained-model representation from
 * src/ml/predict.h, and the difference is worth stating because the two look alike.
 * A distribution is SPECIFIED by its parameters -- a user writes
 * NormalDistribution[mu, sigma] directly -- so it is an ordinary symbolic expression
 * whose arguments mean what they say, and it must print in full because those
 * arguments are the information. A fitted model is DERIVED from data, its parameters
 * are an implementation detail, and it prints elided. So: same mechanism, opposite
 * convention on visibility, and no shared representation.
 *
 * Sampling draws from the user-visible xoshiro stream via random_uniform_01, which is
 * what makes RandomVariate reproducible under SeedRandom. A sampler with its own
 * generator would silently ignore SeedRandom while RandomReal honoured it.
 */
#ifndef ML_DIST_H
#define ML_DIST_H

#include <stddef.h>
#include <stdbool.h>
#include "expr.h"

/* One standard normal deviate, mean 0 variance 1, from the user-visible stream.
 *
 * Box-Muller, generating a pair and caching the second: the transform naturally
 * produces two independent deviates and throwing one away would double the cost of
 * every draw. The cache is per-process state, which means a draw can be served
 * without touching the stream -- that is fine for reproducibility, because the cache
 * is also reset by SeedRandom (see ml_dist_reset_cache) so a reseeded run replays
 * exactly. Without that reset, SeedRandom would leave a stale deviate in hand and the
 * first draw after reseeding would differ. */
double ml_normal_deviate(void);

/* Drop any cached Box-Muller deviate. Called when the stream is reseeded. */
void ml_dist_reset_cache(void);

void ml_dist_init(void);

/* Closed-form properties of the specified distributions NormalDistribution and
 * UniformDistribution (dist_formula.c): the evaluated formula, exact for exact
 * input, or NULL when `dist` is neither family, a numeric parameter set is not a
 * distribution, or (Quantile) the probability is not a number in [0, 1].
 * `x` is the point (PDF, CDF) or probability (Quantile) -- a List threads --
 * and is ignored (may be NULL) for the moments. Arguments are borrowed. The
 * stats heads Mean / Median / Variance / Quantile consult this first. */
typedef enum {
    ML_PROP_PDF, ML_PROP_CDF, ML_PROP_QUANTILE,
    ML_PROP_MEAN, ML_PROP_MEDIAN, ML_PROP_VARIANCE
} MlDistProp;
Expr* ml_dist_closed_form(MlDistProp prop, const Expr* dist, const Expr* x);

/* ml_dist_closed_form, else -- for a one-dimensional SmoothKernel learned
 * distribution -- its machine-number Mean, Variance, Median, CDF or Quantile
 * (dist.c). NULL when neither applies. The hook the stats heads call. */
Expr* ml_dist_property(MlDistProp prop, Expr* dist, const Expr* x);

#endif /* ML_DIST_H */
