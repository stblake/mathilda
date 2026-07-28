
#ifndef PART_H
#define PART_H

#include "expr.h"

Expr* expr_part(Expr* expr, Expr** indices, size_t nindices);
Expr* expr_part_assign(Expr* lhs, Expr* rhs);
Expr* builtin_part(Expr* res);
Expr* builtin_extract(Expr* res);
Expr* expr_head(Expr* e);
Expr* builtin_head(Expr* res);
Expr* builtin_first(Expr* res);
Expr* builtin_last(Expr* res);
Expr* builtin_most(Expr* res);
Expr* builtin_rest(Expr* res);
Expr* expr_insert(Expr* expr, Expr* elem, Expr* pos);
Expr* builtin_insert(Expr* res);
Expr* expr_delete(Expr* expr, Expr* pos);
Expr* builtin_delete(Expr* res);

/* ---------------- Position-path walker (MapAt / ReplaceAt) ----------------
 *
 * The third face of the Part position vocabulary: `expr_part` reads a part,
 * `expr_part_assign` overwrites one, and this walker rewrites one *through a
 * callback*.  MapAt and ReplaceAt differ only in what they do at the addressed
 * leaf, so they share everything else.
 */

/* Applied at each addressed leaf.  `leaf` is borrowed; the callback returns a
 * NEW owned expression (or NULL to abort the whole walk). */
typedef Expr* (*PartLeafFn)(void* ctx, Expr* leaf);

/*
 * Walk `path` (plen position indices) into `expr`, replacing every addressed
 * part with fn(ctx, part).  A path element may be:
 *
 *   - integer k   (>0 from the start, <0 from the end, 0 targets the head)
 *   - All         (every child at this level)
 *   - Span[a, b] / Span[a, b, step]
 *   - Key[k], or a bare literal key, when the node is an Association
 *
 * Returns a new owned expression, or NULL when the position does not exist:
 * an out-of-range integer, an absent key, a malformed spec, or a path that
 * runs into an atom.  Callers surface NULL by leaving the call unevaluated,
 * matching Part[{a,b,c}, 5].  An All or Span that happens to select nothing is
 * a successful no-op, not a failure.
 */
Expr* expr_apply_at_path(Expr* expr, Expr** path, size_t plen,
                         PartLeafFn fn, void* ctx);

/*
 * Position-argument dispatch shared by MapAt and ReplaceAt.  `pos` is either
 * a single path (a bare index, or a List of indices) or a List of paths, in
 * the form Position returns.  Positions are applied left to right over the
 * running result, so a repeated position applies fn twice.
 *
 * The empty list {} is an empty list of *positions* and maps nothing; {{}} is
 * one position -- the empty path -- and so addresses the whole expression.
 */
Expr* expr_apply_at_positions(Expr* expr, Expr* pos, PartLeafFn fn, void* ctx);

#endif
