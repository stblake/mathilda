#ifndef MATHILDA_NDSTRUCT_H
#define MATHILDA_NDSTRUCT_H

/* ---------------------------------------------------------------------------
 * NDArray structural fast paths.
 *
 * Sort / Reverse / Transpose / Flatten / Take / Drop / Clip walk the List tree;
 * an NDArray (a distinct EXPR_NDARRAY node) slips past their List gates and is
 * returned unchanged. These entry points operate directly on the flat buffer
 * and return an NDArray of the appropriate shape (a copy — the input buffer is
 * never mutated).
 *
 * Same faithful-degrade contract as ndreduce.h: any case not handled (custom
 * comparator, complex where the op has no meaning, a level/axis spec beyond the
 * fast domain) returns ndarray_delist_and_reeval(res), so the result always
 * matches the equivalent List call. `res` is borrowed (never freed).
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>

/* True when `res`'s first argument is an NDArray. */
bool ndstruct_call_has_ndarray(const Expr* res);

Expr* ndstruct_sort(Expr* res);       /* Sort[a] */
Expr* ndstruct_reverse(Expr* res);    /* Reverse[a] */
Expr* ndstruct_transpose(Expr* res);  /* Transpose[a] (rank 2) */
Expr* ndstruct_flatten(Expr* res);    /* Flatten[a] */
Expr* ndstruct_take(Expr* res);       /* Take[a, spec] */
Expr* ndstruct_drop(Expr* res);       /* Drop[a, spec] */
Expr* ndstruct_clip(Expr* res);       /* Clip[a] / Clip[a, {min,max}] */

#endif /* MATHILDA_NDSTRUCT_H */
