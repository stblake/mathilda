#ifndef MATHILDA_NDREDUCE_H
#define MATHILDA_NDREDUCE_H

/* ---------------------------------------------------------------------------
 * NDArray reduction & order-statistic fast paths.
 *
 * Total / Mean / Variance / StandardDeviation / RootMeanSquare / Median /
 * Quartiles / Max / Min / Accumulate operate on a List by walking the Expr
 * tree; an NDArray (a distinct EXPR_NDARRAY node, not a List) slips past those
 * List gates and the call is left unevaluated. These entry points close that
 * gap: each computes directly on the flat machine-precision buffer and returns
 * a scalar (Real / Complex[re,im]) or a lower-rank NDArray.
 *
 * Reductions collapse the LEADING axis (matching the List semantics: Total of a
 * matrix adds its rows; Mean/Variance/... of a matrix are columnwise), so a
 * rank-r input yields a rank-(r-1) result, and a vector yields a scalar.
 *
 * Faithful degrade: any case a fast path does not handle (an unsupported level
 * spec, a complex dtype where the op is real-only, rank too high, ...) returns
 * ndarray_delist_and_reeval(res) — the call rebuilt with every NDArray arg
 * replaced by its nested List and re-evaluated — so the result is always
 * identical to the equivalent List call, never wrong, only (in that case) not
 * accelerated. Consequently these never return NULL for an NDArray call.
 *
 * `res` is the whole builtin call expression; it is borrowed (never freed).
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>

/* True when `res` is a function call whose first argument is an NDArray — the
 * signal for a builtin to hand off to the matching ndred_* entry point. */
bool ndred_call_has_ndarray(const Expr* res);

Expr* ndred_total(Expr* res);       /* Total[a] / Total[a,n] / Total[a,Infinity] */
Expr* ndred_mean(Expr* res);        /* Mean[a] */
Expr* ndred_max(Expr* res);         /* Max[a] */
Expr* ndred_min(Expr* res);         /* Min[a] */
Expr* ndred_accumulate(Expr* res);  /* Accumulate[a] */
Expr* ndred_variance(Expr* res);    /* Variance[a] */
Expr* ndred_std(Expr* res);         /* StandardDeviation[a] */
Expr* ndred_rms(Expr* res);         /* RootMeanSquare[a] */
Expr* ndred_median(Expr* res);      /* Median[a] */
Expr* ndred_quartiles(Expr* res);   /* Quartiles[a] */

Expr* ndred_moving_average(Expr* res);  /* MovingAverage[a, r] */
Expr* ndred_moving_median(Expr* res);   /* MovingMedian[a, r] */
Expr* ndred_ema(Expr* res);             /* ExponentialMovingAverage[a, alpha] */

#endif /* MATHILDA_NDREDUCE_H */
