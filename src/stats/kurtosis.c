/* kurtosis.c -- Kurtosis[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout.
 *
 * Kurtosis[data] -- the coefficient of kurtosis, a measure of peak/tail vs flank
 * concentration. Equivalent to CentralMoment[data, 4] / CentralMoment[data, 2]^2
 * (Pearson kurtosis, not the excess form). For a matrix or array it is taken
 * columnwise (the CentralMoment ratio threads). The shared body lives in
 * stats_common.c (stats_standardized_moment), which also routes NDArray / packed
 * inputs to the buffer kernel ndred_kurtosis. */

#include "stats.h"
#include "stats_common.h"

Expr* builtin_kurtosis(Expr* res) {
    return stats_standardized_moment(res, 4);
}
