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
Expr* ndstruct_ordering(Expr* res);   /* Ordering[a] / Ordering[a, seq] -> int64 positions */
Expr* ndstruct_reverse(Expr* res);    /* Reverse[a] */
Expr* ndstruct_transpose(Expr* res);  /* Transpose[a] (rank 2) */
Expr* ndstruct_diagonal(Expr* res);   /* Diagonal[a] / Diagonal[a, k] (rank 2 -> rank 1) */
Expr* ndstruct_flatten(Expr* res);    /* Flatten[a] */
Expr* ndstruct_take(Expr* res);       /* Take[a, spec] */
Expr* ndstruct_drop(Expr* res);       /* Drop[a, spec] */
Expr* ndstruct_clip(Expr* res);       /* Clip[a] / Clip[a, {min,max}] */

/* RotateLeft[a, n] / RotateRight[a, n], with n an Integer or a List of Integers
 * (one per leading axis). `left` selects the direction.
 *
 * Returns NULL -- not a degrade call -- when the spec is outside the fast domain
 * (a non-Integer entry, or more entries than the array has axes), so the caller
 * keeps its existing ndstruct_delist_repack fallback and the List path continues
 * to define those cases.
 *
 * A rotation permutes contiguous blocks, so this is memcpy work, not an element
 * walk. It matters because rotation is how array-style code names a shifted
 * neighbour: going through delist_repack cost 42.6 ms on a 512x512 float64
 * matrix (85x an elementwise add of the same array) and dominated every stencil
 * built out of it. */
Expr* ndstruct_rotate(Expr* res, bool left);

/* Join / Partition / Differences / Riffle / PadLeft / PadRight on the buffer.
 *
 * Same NULL-to-decline contract as ndstruct_rotate: the caller keeps its
 * ndstruct_delist_repack fallback, which still defines every form these do not
 * handle (level specs, offsets, higher-order differences, a cycling List
 * separator, rank >= 2 padding).
 *
 * A head that introduces a NEW element (Riffle's separator, Pad's fill) only
 * uses the buffer when that element is exactly representable at the buffer's
 * dtype with a matching head. PadLeft[{1.,2.,3.}, 5] is {0, 0, 1., 2., 3.} in
 * both Mathilda and Mathematica -- exact zeros beside Reals, which no uniform
 * buffer holds -- so it declines and the List path answers. */
/* First / Last / Most / Rest on the buffer.
 *
 * Each is a leading-axis slice — one row, or all-but-one — so each is a pointer
 * read or a single memcpy. They previously went through ndstruct_delist_repack,
 * which made First and Last asymptotically wrong: reading ONE element of a 10^6
 * float64 vector cost 123 ms, next to 0.88 ms for the identical Drop[v, 250].
 *
 * Degrades (delist_and_reeval) for the empty results the buffer layer has no
 * shape for: Rest/Most of a single row is {}. */
typedef enum { ND_FIRST, ND_LAST, ND_MOST, ND_REST } NDHeadTail;
Expr* ndstruct_head_tail(Expr* res, NDHeadTail which);

Expr* ndstruct_join(Expr* res);
Expr* ndstruct_partition(Expr* res);
Expr* ndstruct_differences(Expr* res);
Expr* ndstruct_riffle(Expr* res);
Expr* ndstruct_pad(Expr* res, bool left);

/* The ninth round's additions: heads that had no buffer path on EITHER surface,
 * so a visible NDArray came back unevaluated and a packed List was materialised
 * first. Each returns NULL for anything outside its fast domain, and the caller
 * degrades so the answer still matches the List path exactly.
 *
 *   ndstruct_ratios         Ratios[a]                inexact dtypes only:
 *                           ratios of integers are Rationals, which no buffer
 *                           holds
 *   ndstruct_append         Append[a, x] (prepend=false) / Prepend[a, x]
 *                           one allocation, two memcpys. Both were on pack.c's
 *                           "correct by omission" list, which cost a full
 *                           materialise to add a single element
 *   ndstruct_catenate       Catenate[arr] (rank >= 2: a reshape, since a
 *                           row-major buffer already stores the rows
 *                           consecutively) and Catenate[{a1, a2, ...}] (Join)
 *   ndstruct_take_extreme   TakeLargest[a, k] / TakeSmallest[a, k] in
 *                           O(n log k) through a bounded heap -- ahead of the
 *                           full sort NumPy's own idiom for this uses
 */
Expr* ndstruct_ratios(Expr* res);
Expr* ndstruct_append(Expr* res, bool prepend);
Expr* ndstruct_catenate(Expr* res);
Expr* ndstruct_take_extreme(Expr* res, bool largest);

/* Faithful-degrade primitive for a head that has NO native buffer walk.
 *
 * Re-runs `call` with every NDArray argument materialised into a nested List,
 * then repacks a rectangular machine-numeric result with `src`'s dtype. The
 * ANSWER therefore comes from the ordinary List implementation and is identical
 * to the List call by construction; only the packing of the result is new.
 *
 * This exists because an NDArray is ATOMIC, so a builtin that walks
 * `arg->data.function.args` looks straight past one and returns the call
 * UNEVALUATED — which is what Select, Join, First, Differences and a dozen
 * others used to do, while the identical List call worked fine.
 *
 * A scalar, symbolic, empty or ragged result is returned exactly as the List
 * implementation produced it (there is nothing to pack), matching what Map does
 * with a non-numeric result. `call` and `src` are borrowed. */
Expr* ndstruct_delist_repack(const Expr* call, const Expr* src);

#endif /* MATHILDA_NDSTRUCT_H */
