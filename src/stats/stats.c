/* stats.c -- registration hub for the statistics builtins.
 *
 * Each builtin lives in its own translation unit under src/stats/ (one builtin
 * per file, filename = lowercased builtin name); shared helpers live in
 * stats_common.c. This file only registers the builtins and sets their
 * attributes, mirroring numbertheory.c in src/numbertheory/. */

#include "stats.h"
#include "symtab.h"
#include "attr.h"

void stats_init(void) {
    symtab_add_builtin("Mean", builtin_mean);
    symtab_add_builtin("RootMeanSquare", builtin_rootmeansquare);
    symtab_add_builtin("Median", builtin_median);
    symtab_get_def("Median")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Quartiles", builtin_quartiles);
    symtab_get_def("Quartiles")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Variance", builtin_variance);
    symtab_add_builtin("CentralMoment", builtin_central_moment);
    symtab_add_builtin("Skewness", builtin_skewness);
    symtab_add_builtin("Kurtosis", builtin_kurtosis);
    symtab_add_builtin("StandardDeviation", builtin_standard_deviation);
    symtab_add_builtin("Covariance", builtin_covariance);
    symtab_add_builtin("Correlation", builtin_correlation);
    symtab_add_builtin("MovingAverage", builtin_moving_average);
    symtab_add_builtin("MovingMedian", builtin_moving_median);
    symtab_add_builtin("ExponentialMovingAverage", builtin_exponential_moving_average);

    symtab_get_def("Mean")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RootMeanSquare")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Variance")->attributes |= ATTR_PROTECTED;
    symtab_get_def("CentralMoment")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Skewness")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Kurtosis")->attributes |= ATTR_PROTECTED;
    symtab_get_def("StandardDeviation")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Covariance")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Correlation")->attributes |= ATTR_PROTECTED;
    symtab_get_def("MovingAverage")->attributes |= ATTR_PROTECTED;
    symtab_get_def("MovingMedian")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ExponentialMovingAverage")->attributes |= ATTR_PROTECTED;
}
