#ifndef NUMERIC_ORDER_H
#define NUMERIC_ORDER_H

/* numeric_order -- order expressions by NUMERIC VALUE (the reference's
 * NumericalOrder), for the heads that rank data: Median, Quantile, Quartiles,
 * InterquartileRange, MedianDeviation, TakeLargest/TakeSmallest(+By),
 * RankedMin/RankedMax.
 *
 * Canonical order (expr_compare, what Sort uses) is NOT numeric order once an
 * exact non-rational number is involved: it puts every explicit number before
 * Sqrt[2], E, Pi, -Pi or -Infinity, so Median[{Sqrt[2], 1, 3}] was 3 and
 * TakeLargest[{-Infinity, 3}, 1] was {-Infinity}.
 *
 * The order:
 *   - real numbers (every real-valued numeric quantity, +/-Infinity included)
 *     come first, ordered by value; two of equal value TIE (0) -- the caller
 *     settles a tie (canonical order, or position);
 *   - everything else follows, in canonical order.
 *
 * The value comparison is CERTIFIED: two machine numbers (or integers that a
 * double holds exactly) compare in O(1); otherwise the cached machine
 * approximations decide only a clear separation of operands that cannot
 * cancel, and anything closer goes to the comparison heads' certified
 * comparator (exact GMP, then MPFR at increasing precision). A pair that
 * comparator cannot separate ties.
 *
 * Use: compute one NumOrdKey per element (numord_key, which may numericalize
 * a symbolic element ONCE), then compare with numord_cmp. */

#include "expr.h"
#include <stdbool.h>

enum { NUMORD_REAL = 0, NUMORD_OTHER = 1 };

typedef struct {
    double v;              /* the value: exact when `exact`, else approximate */
    unsigned char cls;     /* NUMORD_REAL or NUMORD_OTHER */
    unsigned char exact;   /* v IS the value (machine real, small int, +/-Inf) */
    unsigned char sep;     /* v is accurate to ~1 ulp: a clear gap is decisive */
} NumOrdKey;

/* The key of e (borrowed). O(1) for a number atom; a symbolic element is
 * numericalized once at machine precision. */
void numord_key(Expr* e, NumOrdKey* k);

/* Numeric order of a and b given their keys: <0, 0 (tie), >0. */
int numord_cmp(Expr* a, const NumOrdKey* ka, Expr* b, const NumOrdKey* kb);

/* numord_cmp with the keys computed here: for one-off comparisons. */
int numord_compare(Expr* a, Expr* b);

#endif
