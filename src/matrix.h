#ifndef MATHILDA_MATRIX_H
#define MATHILDA_MATRIX_H

/* ---------------------------------------------------------------------------
 * Matrix — a first-class, visibly-distinct dense machine-precision ndarray.
 *
 * Unlike Mathematica's packed arrays, which are an invisible internal
 * optimization a list may or may not have, `Matrix[...]` is always what it
 * says it is: `Head[Matrix[{{1,2},{3,4}}]]` is `Matrix`, never `List`, and
 * `MatrixQ`/`ListQ` never disagree with each other about which one a value
 * is. It is a real `EXPR_MATRIX` node (see expr.h's `MatrixData`), storing a
 * flat row-major `double*` buffer plus rank/dims directly in the `Expr` —
 * not a nested `List[List[...]]` tree — so builtins that recognize it can
 * operate on the buffer directly instead of paying the flatten/rebuild cost
 * every List-based linalg operation pays today (see src/linalg/dot.c).
 *
 * Machine precision only. Any operation that would need to store a
 * non-double entry (a symbol, an exact/rational number, a BigInt, an MPFR
 * value) degrades to an ordinary nested List instead of forcing a lossy
 * conversion — see matrix_from_nested_list / matrix_degrade.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* True when `e` is an EXPR_MATRIX node. */
bool is_matrix(const Expr* e);

/* Total element count (product of dims). */
size_t matrix_size(const Expr* m);

/* Probe a nested List for machine-precision-packable shape: rejects jagged
 * lists and any leaf that isn't EXPR_INTEGER/EXPR_REAL. On success returns
 * a newly built Matrix[...] (caller owns it); on failure (jagged, empty,
 * non-machine-precision leaf, or not a List at all) returns NULL and the
 * caller keeps the original list unevaluated. */
Expr* matrix_from_nested_list(const Expr* list);

/* Inverse of matrix_from_nested_list: rebuilds the equivalent nested
 * List[...] tree from a Matrix's flat buffer. Caller owns the result. */
Expr* matrix_to_nested_list(const Expr* m);

/* Matrix[nested_list] constructor builtin: packs `res`'s sole argument, or
 * returns NULL (leave unevaluated) if it can't be packed. */
Expr* builtin_matrix(Expr* res);

/* Fast C-level Dot for two Matrix operands of rank 1 or 2: contracts the
 * trailing axis of `a` with the leading axis of `b` using raw double loops
 * (no symbolic Times/Plus). Returns a new EXPR_MATRIX (rank 1 or 2), or a
 * bare EXPR_REAL for a vector.vector contraction (rank 0). Returns NULL if
 * either rank is > 2 (caller should fall back to the generic tensor path);
 * sets *shape_error and returns NULL if the contracted dimensions disagree. */
Expr* matrix_dot2(const Expr* a, const Expr* b, bool* shape_error);

/* Fast C-level elementwise Plus (is_plus=true) or Times (is_plus=false) over
 * n Matrix operands: loops over the flat double buffers directly, no
 * symbolic Times/Plus per element. Returns a new EXPR_MATRIX with the same
 * shape, or NULL if any operand isn't EXPR_MATRIX or the shapes disagree
 * (caller falls back to the generic symbolic path, which treats mismatched
 * Matrix operands as opaque non-numeric terms). */
Expr* matrix_elementwise(Expr** args, size_t n, bool is_plus);

/* Register Matrix's builtins, attributes, and docstring. */
void matrix_init(void);

#endif /* MATHILDA_MATRIX_H */
