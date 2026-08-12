/* skewness.c -- Skewness[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout.
 *
 * Skewness[data] -- the coefficient of skewness, a measure of asymmetry.
 * Equivalent to CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2). For a
 * matrix or array it is taken columnwise (the CentralMoment ratio threads). The
 * shared body lives in stats_common.c (stats_standardized_moment), which also
 * routes NDArray / packed inputs to the buffer kernel ndred_skewness. */

#include "stats.h"
#include "stats_common.h"

Expr* builtin_skewness(Expr* res) {
    return stats_standardized_moment(res, 3);
}
