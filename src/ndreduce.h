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

/* Full reduction of NDArray `a` (borrowed) to a scalar — every axis collapsed,
 * i.e. Total[a, Infinity]. Unlike the entry points above this takes the ARRAY
 * rather than the enclosing call, so a caller that already holds the operand
 * (the Compile VM's OP_V_TOTAL) needs no wrapper call node. Shares
 * ndred_total's summation, so the rounding is identical. Caller owns the
 * result; NULL if `a` is not an NDArray. */
Expr* ndred_total_all(const Expr* a);
Expr* ndred_mean(Expr* res);        /* Mean[a] */
Expr* ndred_max(Expr* res);         /* Max[a] */
Expr* ndred_min(Expr* res);         /* Min[a] */
Expr* ndred_accumulate(Expr* res);  /* Accumulate[a] */
Expr* ndred_variance(Expr* res);    /* Variance[a] */
Expr* ndred_std(Expr* res);         /* StandardDeviation[a] */
Expr* ndred_rms(Expr* res);         /* RootMeanSquare[a] */
Expr* ndred_median(Expr* res);      /* Median[a] */
Expr* ndred_quartiles(Expr* res);   /* Quartiles[a] */

/* Tally[a] — an IRREGULAR reduction: a scatter-add into a keyed table rather
 * than a fixed-stride sweep. Distinct enough from the reductions above to be
 * worth naming: the win is not vectorising a loop but never boxing the input at
 * all, since the generic path must materialise one Expr per element and then
 * hash and compare Exprs. Returns the {value, count} pairs in first-appearance
 * order, exactly as the List path does. */
Expr* ndred_tally(Expr* res);           /* Tally[a] */

/* ------------------------------------------------------------------- scans
 *
 * Fold / FoldList over an ASSOCIATIVE machine operator, straight on the buffer.
 *
 * A scan is the one shape an array language cannot vectorise away, so it is
 * where an evaluator is most exposed: `FoldList[Max, First[v], Rest[v]]` — the
 * running maximum, i.e. numpy's `np.maximum.accumulate` — cost 605 ms on a 10^6
 * float64 vector against 2.7 ms in NumPy, because every step applied `Max`
 * through the evaluator and every element was boxed on the way in and out.
 *
 * Only four operators matter in practice (running sum, product, max, min) and
 * all four are a two-line C loop. The general case is NOT this function's job:
 * a non-associative or user-written body stays with numloop's compiled VM, which
 * reads the same buffer.
 *
 * `seed` may be NULL for the 2-arg spelling FoldList[f, a], which seeds from
 * element 0 and scans the rest. `as_list` selects FoldList over Fold.
 *
 * DECLINES (returns NULL, caller keeps its existing path) for rank >= 2, a
 * complex dtype, a seed that is not a machine number, and — this is the
 * exactness gate — any combination whose List answer would NOT be uniform:
 * Max returns one of its ARGUMENTS, so an exact seed beside a Real buffer
 * (Fold[Max, 1, {0.5}] -> the exact 1) escapes the buffer's dtype. */
typedef enum { ND_SCAN_PLUS, ND_SCAN_TIMES, ND_SCAN_MAX, ND_SCAN_MIN } NDScanOp;

/* True, setting *op, when `f` names one of the four scan operators — either as
 * a bare symbol (Max) or as the equivalent pure function (Max[#1, #2] &). Both
 * spellings mean the same scan and used to differ by 130 ms. */
bool ndred_scan_op_for(const Expr* f, NDScanOp* op);

Expr* ndred_scan(const Expr* a, NDScanOp op, const Expr* seed, bool as_list);

Expr* ndred_moving_average(Expr* res);  /* MovingAverage[a, r] */
Expr* ndred_moving_median(Expr* res);   /* MovingMedian[a, r] */
Expr* ndred_ema(Expr* res);             /* ExponentialMovingAverage[a, alpha] */

#endif /* MATHILDA_NDREDUCE_H */
