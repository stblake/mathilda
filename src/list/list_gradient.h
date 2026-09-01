#ifndef LIST_GRADIENT_H
#define LIST_GRADIENT_H

#include "expr.h"

/* ListGradient[f]                     — numerical gradient of a sampled array,
 * ListGradient[f, spacing]              a from-first-principles port of
 * ListGradient[f, spacing, opts...]     numpy.gradient. Returns second-order
 *                                       central differences in the interior and
 * one-sided differences at the edges, per axis, for an array of any rank.
 *
 * Options: Method -> "Centered"|"Forward"|"Backward", DifferenceOrder -> p,
 * WindowLength -> Automatic|m, Axis -> All|a|{a,...}. See src/list/list_gradient.c
 * and docs/spec/builtins/lists.md. */
Expr* builtin_list_gradient(Expr* res);

/* Buffer fast path AND the Compile[] ND_FNS delegate (one function, two callers).
 * Reads a float64/float32 NDArray (packed or visible) and returns a same-shape
 * float array — a single NDArray when one axis is computed, else a List of them.
 * Returns NULL to decline (non-float dtype, symbolic/non-numeric spacing, an
 * unreadable option, or an axis shorter than 2), so the caller falls back to the
 * exact/symbolic List path via ndstruct_delist_repack. */
Expr* list_gradient_ndarray(Expr* res);

#endif /* LIST_GRADIENT_H */
