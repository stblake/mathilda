#ifndef STATS_COMMON_H
#define STATS_COMMON_H

/* Cross-cutting helpers shared between the per-builtin translation units of
 * src/stats/ (mirrors the numbertheory_internal.h pattern in
 * src/numbertheory/). Each helper is defined once in stats_common.c. */

#include "expr.h"
#include <stdbool.h>

/* True for a real/rational/complex numeric leaf. On success stores the real
 * approximation in *val (for real inputs) and whether the value is complex in
 * *out_complex; either out-pointer may be NULL. (Named out_complex, not
 * complex, because <complex.h> makes `complex` a reserved macro.) */
bool stats_is_numeric(Expr* e, double* val, bool* out_complex);

/* True when a machine-double statistic is a machine number under the
 * boundary policy (docs/design/numeric_tower_review.md): finite and not a
 * nonzero value below DBL_MIN. A fast path whose result fails this declines
 * to its symbolic path, which keeps the value at machine precision in MPFR:
 * Mean[{1.*^308, 1.*^308}] was Overflow[] where the answer is 1.*^308. */
bool stats_machine_ok(double v);

/* Apply the builtin named func_name column-wise to a matrix, i.e.
 * Map[func_name, Transpose[matrix]]. Returns NULL if the transpose is not a
 * function. */
Expr* stats_apply_columnwise(const char* func_name, Expr* matrix);

/* True iff e is NumericQ and free of the imaginary unit I (a real number). */
bool stats_is_real_numeric(Expr* e);

/* One parameterized quantile point over an already-sorted arg vector: Wolfram's
 * Quantile[list, q, {{a,b},{c,d}}] definition. h = a + (n+b)q, edge-clamped to
 * the first/last element; at INTEGER h the Floor/Ceiling neighbors coincide and
 * the result is sorted_args[h-1] regardless of c,d (Quartiles' own c=0 default
 * made the two forms agree, so the pre-extraction loop never exercised this
 * case; Quantile's c=1,d=0 default does on every integer (n+b)q+a). Otherwise
 * A[j] + (c + d(h-j))(A[j+1] - A[j]), j = Floor[h] clamped to [1, n-1].
 * Arithmetic stays in the evaluator (exact in, exact out); doubles are read
 * only for clamping and the integer test, mirroring the original loop.
 * Returns a fresh Expr*; Indeterminate when h does not evaluate numerically. */
Expr* stats_quantile_point(Expr** sorted_args, size_t n, Expr* q,
                           Expr* a, Expr* b, Expr* c, Expr* d);

/* Standardized central moment m_p / m_2^(p/2), the shared body of Skewness
 * (p=3) and Kurtosis (p=4). Builds CentralMoment[data, p] / CentralMoment[data,
 * 2]^(p/2), which threads columnwise over matrices/arrays and stays exact or
 * symbolic through the evaluator; an NDArray argument takes the buffer fast path
 * (ndred_skewness / ndred_kurtosis). Returns NULL — leaving the caller's head
 * unevaluated — when the data does not reduce (a bare symbol, an empty list). */
Expr* stats_standardized_moment(Expr* res, int p);

/* The vector cores over n elements el[0..n) (borrowed): a List's arguments or
 * a LinkedList's elements (LLElems), so both containers run one algorithm.
 * NULL when the statistic is undefined there (n too small). */
Expr* stats_mean_vec(Expr** el, size_t n);
Expr* stats_variance_vec(Expr** el, size_t n);
/* The exact / symbolic trees of Moment and CentralMoment over elements: Mean
 * of el^r, and Mean of (el - Mean[el])^r, one element at a time -- what the
 * Listable Power and Plus do over a List (and a matrix's columns: the List
 * code reduces a matrix through these trees, not the vector fast loops). */
Expr* stats_moment_elements(Expr** el, size_t n, Expr* order);
Expr* stats_cm_elements(Expr** el, size_t n, Expr* order);

/* head @@ seq, unevaluated: head[elements...] of a List or a LinkedList (read
 * from the slots), else Apply[head, seq]. Consumes seq. */
Expr* stats_apply_head(const char* head, Expr* seq);

/* The moving window el[i, i + r) of `data` (List or LinkedList) as a sequence
 * of the same kind: a List, or a LinkedList whose slots are copied. */
Expr* stats_moving_window(const Expr* data, Expr** el, size_t i, size_t r);

/* The statistics with a LinkedList arm (stats_ll.c). */
typedef enum {
    STATS_MEAN, STATS_VARIANCE, STATS_STD, STATS_RMS, STATS_MOMENT,
    STATS_CENTRAL_MOMENT, STATS_SKEWNESS, STATS_KURTOSIS, STATS_MEDIAN,
    STATS_QUARTILES, STATS_IQR, STATS_MEAN_DEVIATION, STATS_MEDIAN_DEVIATION,
    STATS_MOVING_AVERAGE, STATS_MOVING_MEDIAN, STATS_EMA
} StatsOp;

/* The LinkedList arm of statistic `op`, tried before the List code: a machine
 * LinkedList vector through the packed kernels on its chunks, or a matrix with
 * a LinkedList in it columnwise. True (answer, or NULL for unevaluated, in
 * *out) when handled; false otherwise -- a LinkedList vector then continues
 * into the builtin's element path, read through LLElems (linkedlist.h). */
bool stats_ll_arm(Expr* res, StatsOp op, Expr** out);

#endif
