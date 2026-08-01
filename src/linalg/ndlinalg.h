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
/* Induced matrix norm of an NDArray, straight to LAPACK (dlange/zlange, or
 * gesdd for the spectral 2-norm). `pe` may be NULL for the default 2-norm.
 *
 * Unlike ndla_norm this NEVER defers via linalg_delist_and_reeval -- it returns
 * NULL on anything it does not handle. src/linalg/norm.c calls it to give a
 * plain machine-numeric matrix the same answer a packed one gets, and deferring
 * there would re-enter builtin_norm and recurse without bound. */
Expr* ndla_matrix_norm_direct(Expr* v, Expr* pe);
Expr* ndla_normalize(Expr* res);    /* Normalize[v]                   */
Expr* ndla_cross(Expr* res);        /* Cross[u, v] (3-vectors)        */

/* ------------------------------------------------------------------------
 * The SVD pair: PseudoInverse[m] and LeastSquares[m, b].
 *
 * These do NOT take the call and do NOT defer via linalg_delist_and_reeval,
 * for the reason plans/HPC_IMPROVEMENT_PLAN.md item 10.2 records: a fast path
 * that declines by re-evaluating the call re-enters the builtin that dispatched
 * to it. `NULL` means "I decline"; the caller then runs its own path with the
 * arguments it already has, and there is no way back in.
 *
 * Both are restricted to a REAL rank-2 float64 buffer, which is what makes them
 * safe rather than merely fast: the exact pipeline in inv.c is the right answer
 * for an exact matrix, and an integer buffer never arrives here because neither
 * head is on pack.c's INT64_OK list. Complex input declines.
 *
 * `tol_automatic` / `tol_value` carry the call's Tolerance option. Automatic is
 * the LAPACK convention max(m,n) * eps * sigma_max; an explicit value is taken
 * as a fraction of sigma_max, matching Mathematica's documented meaning.
 * ---------------------------------------------------------------------- */
Expr* ndla_pseudoinverse_direct(const Expr* a, bool tol_automatic, double tol_value);
Expr* ndla_leastsquares_direct(const Expr* a, const Expr* b,
                               bool tol_automatic, double tol_value);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_NDLINALG_H */
