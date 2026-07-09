#ifndef MATHILDA_NDARRAY_H
#define MATHILDA_NDARRAY_H

/* ---------------------------------------------------------------------------
 * NDArray — a first-class, visibly-distinct dense machine-precision ndarray.
 *
 * Modelled on numpy's ndarray: an N-dimensional (rank 1..NDARRAY_MAX_RANK),
 * rectangular, C-order (row-major) dense array of machine-precision values.
 * Unlike Mathematica's packed arrays, which are an invisible internal
 * optimization a list may or may not have, `NDArray[...]` is always what it
 * says it is: `Head[NDArray[{{1,2},{3,4}}]]` is `NDArray`, never `List`, and
 * `ArrayQ`/`ListQ` never disagree with each other about which one a value is.
 * It is a real `EXPR_NDARRAY` node (see expr.h's `NDArrayData`), storing a
 * flat row-major `double*` buffer plus rank/dims directly in the `Expr` —
 * not a nested `List[List[...]]` tree — so builtins that recognize it can
 * operate on the buffer directly instead of paying the flatten/rebuild cost
 * every List-based tensor operation pays today (see src/linalg/dot.c).
 *
 * numpy correspondence: `Dimensions` is `shape`, `ArrayDepth` is `ndim`,
 * `Length` is `len` (= shape[0]). Rank is unbounded up to NDARRAY_MAX_RANK.
 *
 * Machine precision only. Any operation that would need to store a
 * non-double entry (a symbol, an exact/rational number, a BigInt, an MPFR
 * value) degrades to an ordinary nested List instead of forcing a lossy
 * conversion — see ndarray_from_nested_list.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* Maximum supported rank (number of axes). */
#define NDARRAY_MAX_RANK 64

/* True when `e` is an EXPR_NDARRAY node. */
bool is_ndarray(const Expr* e);

/* Total element count (product of dims). */
size_t ndarray_size(const Expr* a);

/* Probe a nested List for machine-precision-packable shape: rejects jagged
 * lists and any leaf that isn't EXPR_INTEGER/EXPR_REAL. On success returns
 * a newly built NDArray[...] (caller owns it); on failure (jagged, empty,
 * non-machine-precision leaf, or not a List at all) returns NULL and the
 * caller keeps the original list unevaluated. */
Expr* ndarray_from_nested_list(const Expr* list);

/* Inverse of ndarray_from_nested_list: rebuilds the equivalent nested
 * List[...] tree from an NDArray's flat buffer. Caller owns the result. */
Expr* ndarray_to_nested_list(const Expr* a);

/* NDArray[nested_list] constructor builtin: packs `res`'s sole argument, or
 * returns NULL (leave unevaluated) if it can't be packed. */
Expr* builtin_ndarray(Expr* res);

/* NDArrayQ[expr] builtin: True iff expr is an NDArray value, else False. */
Expr* builtin_ndarrayq(Expr* res);

/* Fast C-level Dot for two NDArray operands of rank 1 or 2: contracts the
 * trailing axis of `a` with the leading axis of `b` using raw double loops
 * (no symbolic Times/Plus). Returns a new EXPR_NDARRAY (rank 1 or 2), or a
 * bare EXPR_REAL for a vector.vector contraction (rank 0). Returns NULL if
 * either rank is > 2 (caller should fall back to the generic tensor path);
 * sets *shape_error and returns NULL if the contracted dimensions disagree. */
Expr* ndarray_dot2(const Expr* a, const Expr* b, bool* shape_error);

/* Fast C-level elementwise Plus (is_plus=true) or Times (is_plus=false) over
 * n NDArray operands: loops over the flat double buffers directly, no
 * symbolic Times/Plus per element. Returns a new EXPR_NDARRAY with the same
 * shape, or NULL if any operand isn't EXPR_NDARRAY or the shapes disagree
 * (caller falls back to the generic symbolic path, which treats mismatched
 * NDArray operands as opaque non-numeric terms). Broadcasting (numpy-style
 * shape compatibility) is intentionally not yet handled here. */
Expr* ndarray_elementwise(Expr** args, size_t n, bool is_plus);

/* When every operand in args[0..n) is an NDArray but they are not all the same
 * shape, print a one-line `NDArray::shape` warning naming the two offending
 * shapes (verb is the operation, e.g. "added"/"multiplied") and return true.
 * Returns false when the operands are not all NDArrays or the shapes are
 * uniform (nothing to warn about). Respects the arithmetic-warning mute so it
 * stays quiet inside internal computations. */
bool ndarray_warn_shape_mismatch(Expr** args, size_t n, const char* verb);

/* Register NDArray's builtins, attributes, and docstring. */
void ndarray_init(void);

#endif /* MATHILDA_NDARRAY_H */
