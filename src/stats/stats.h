#ifndef STATS_H
#define STATS_H

#include "expr.h"

Expr* builtin_mean(Expr* res);
Expr* builtin_rootmeansquare(Expr* res);
Expr* builtin_median(Expr* res);
Expr* builtin_quartiles(Expr* res);
Expr* builtin_variance(Expr* res);
Expr* builtin_central_moment(Expr* res);
Expr* builtin_skewness(Expr* res);
Expr* builtin_kurtosis(Expr* res);
Expr* builtin_standard_deviation(Expr* res);
Expr* builtin_covariance(Expr* res);
Expr* builtin_correlation(Expr* res);
Expr* builtin_moving_average(Expr* res);
Expr* builtin_moving_median(Expr* res);
Expr* builtin_exponential_moving_average(Expr* res);

void stats_init(void);

#endif
