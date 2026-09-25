#ifndef MATHILDA_ASSOC_OPS_H
#define MATHILDA_ASSOC_OPS_H

/* ---------------------------------------------------------------------------
 * assoc_ops — the second tier of the Association family (see assoc.c for the
 * data structure and the first-tier heads).
 *
 * This module owns:
 *   - new heads: KeyIntersection, KeyComplement, MissingQ, JoinAcross,
 *     ApplyTo, Discard, CountDistinct, CountDistinctBy, SubsetQ,
 *     AssociationComap, and the Splice construct;
 *   - extended forms of existing heads (Keys/Values with a function and over
 *     lists, KeyUnion with a fill function, AssociationMap over an association,
 *     PositionIndex of an association, Merge over rules, nested rule lists in
 *     Association, deep Normal, DeleteMissing with levels, Lookup over mixed
 *     lists, Transpose of associations).  These are installed as thin wrappers
 *     around the builtin already registered for the head, so the first-tier
 *     implementation keeps handling every form it handled before;
 *   - Listable threading over association values (the evaluator hook below).
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>

/* Listable threading over associations, called by the evaluator for a Listable
 * head none of whose arguments is a List.  If at least one argument is a
 * well-formed association, every association argument must have the same keys
 * in the same order; the result is <|k -> f[..., v_k, ...], ...|> with the
 * non-association arguments repeated.  Returns a fresh (unevaluated)
 * association, or NULL when no association argument is present, an
 * association is malformed, or the keys disagree (Association::incmp).
 * Borrows `e`. */
Expr* assoc_thread_listable(Expr* e);

/* Splice[list] / Splice[list, h] inside the arguments of `e`: replace each
 * Splice whose head spec matches e's head (List or Association by default)
 * by the elements of its list.  Called by the evaluator's Sequence-flattening
 * step only when a Splice argument was seen.  Rewrites `e` in place; returns
 * true iff the arguments changed. */
bool eval_splice_args(Expr* e);

/* Register the heads and wrappers above.  Runs after every other subsystem has
 * registered its builtins (it wraps Normal, DeleteMissing and Transpose too). */
void assoc_ops_init(void);

#endif /* MATHILDA_ASSOC_OPS_H */
