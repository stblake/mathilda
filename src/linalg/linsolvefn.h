/* Mathilda -- LinearSolveFunction: LinearSolve[m] as a reusable solver.
 *
 *   f = LinearSolve[m]            (m dense or SparseArray; options allowed)
 *   f[b]  f[B]  f[b, "T"]         solve m.x == b, m.X == B, Transpose[m].x == b
 *   f[b, "C"]  f[b, "J"]          conjugate transpose / conjugate systems
 *   f["Properties"], f["MatrixDimensions"], f["ConditionNumber"], f["getL"],
 *   f["getU"], f["getPermutations"]
 *
 * The object is an ordinary expression,
 *
 *   LinearSolveFunction[{rows, cols}, {key, m, {opts...}}]
 *
 * (Length 2, as the reference's), so it prints, copies, compares and is
 * stored like any other expression -- no new node type. `key` names a
 * machine factorisation held in a small LRU cache (src/linalg/linsolvefn.c);
 * key 0 means "no factorisation" (exact, MPFR or symbolic entries, a
 * singular or rectangular matrix, the Krylov method), and every application
 * then solves LinearSolve[m, b, opts] afresh, which is exactly the reference's
 * behaviour for those cases. Because the object carries m itself, an evicted
 * factor is simply recomputed: the cache bounds memory without ever changing
 * an answer.
 */
#ifndef MATHILDA_LINSOLVEFN_H
#define MATHILDA_LINSOLVEFN_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* LinearSolve[m] / LinearSolve[m, opts...] -> LinearSolveFunction[...].
 * `call` is borrowed. Returns a new expression, or NULL to decline (m is not a
 * matrix). Messages: LinearSolve::sing1, ::sqmat1. */
Expr* linsolvefn_make(const Expr* call);

/* True when `head` is LinearSolveFunction[{r, c}, {key, m, {opts}}]. */
bool linsolvefn_is(const Expr* head);

/* Apply LinearSolveFunction `head` to `args` (borrowed). NULL leaves the call
 * unevaluated (after a message where the reference gives one). */
Expr* linsolvefn_apply(const Expr* head, Expr* const* args, size_t nargs);

/* Release every cached factorisation (tests and leak checks). */
void linsolvefn_cache_clear(void);

void linsolvefn_init(void);

#endif /* MATHILDA_LINSOLVEFN_H */
