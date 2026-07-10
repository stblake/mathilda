/*
 * src/linalg/ndlinalg.h
 *
 * NDArray fast-path layer for the linear-algebra builtins.
 *
 * `NDArray[...]` (EXPR_NDARRAY) is Mathilda's dense, row-major, machine-
 * precision numeric container. Every linalg builtin routes its NDArray inputs
 * through this module so the routine body never has to know about NDArray:
 * each builtin adds a single guard at the top,
 *
 *     if (linalg_call_has_ndarray(res)) return ndla_<op>(res);
 *
 * and `ndla_<op>` either
 *   - runs a genuine fast path over the flat `double` buffer (loaded via
 *     numarray.c's `na_load_*`, computed with an in-house double kernel or a
 *     `mat_lapack_*` routine, and rebuilt as an NDArray so the numeric linalg
 *     surface stays a closed system), or
 *   - falls back to `linalg_delist_and_reeval(res)` — the universal safety net
 *     that converts every top-level NDArray argument to a nested `List` and
 *     re-evaluates the call, so the routine's ordinary numeric-List path runs.
 *
 * All fast paths honour the builtin ownership contract: they never free `res`
 * (the evaluator does), return a freshly-owned `Expr*` on success, or `NULL`
 * (leaving the call unevaluated) on a genuine error after emitting the routine's
 * diagnostic. Loaded `double` buffers are always freed; the `NULL`/fall-through
 * branches allocate nothing that leaks.
 */
#ifndef MATHILDA_LINALG_NDLINALG_H
#define MATHILDA_LINALG_NDLINALG_H

#include "expr.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* True when the call `res` (an EXPR_FUNCTION) has at least one direct argument
 * that is an EXPR_NDARRAY. This is the cheap predicate each builtin guards on. */
bool linalg_call_has_ndarray(const Expr* res);

/* Universal fallback: rebuild `res` with every top-level NDArray argument
 * replaced by its nested-List form and evaluate the result. Returns a freshly
 * owned Expr* (the evaluated call), or NULL if the rebuilt call cannot be
 * evaluated. Never frees `res`. Because the rebuilt call has no NDArray
 * arguments it cannot re-enter a `linalg_call_has_ndarray` guard, so there is
 * no recursion. */
Expr* linalg_delist_and_reeval(Expr* res);

/* Per-routine NDArray fast paths. Each takes the full call `res`, parses the
 * arguments it needs, and returns a fast-path result, a genuine-error NULL, or
 * `linalg_delist_and_reeval(res)` when the buffer path does not apply. */
Expr* ndla_det(Expr* res);          /* Det[m]                         */
Expr* ndla_inverse(Expr* res);      /* Inverse[m]                     */
Expr* ndla_linearsolve(Expr* res);  /* LinearSolve[m, b]              */
Expr* ndla_matrixrank(Expr* res);   /* MatrixRank[m]                  */
Expr* ndla_tr(Expr* res);           /* Tr[m] (2-D, default Plus)      */
Expr* ndla_norm(Expr* res);         /* Norm[v] / Norm[v, p] / Norm[m] */
Expr* ndla_normalize(Expr* res);    /* Normalize[v]                   */
Expr* ndla_cross(Expr* res);        /* Cross[u, v] (3-vectors)        */

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_NDLINALG_H */
