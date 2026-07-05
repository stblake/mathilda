#ifndef MATHILDA_ASSOC_H
#define MATHILDA_ASSOC_H

/* ---------------------------------------------------------------------------
 * Association — the Wolfram-Language <| key -> value, ... |> data structure.
 *
 * Representation.  An association is stored as an ordinary expression:
 *
 *     <|a -> 1, b -> 2|>   ==   Association[Rule[a, 1], Rule[b, 2]]
 *
 * i.e. an EXPR_FUNCTION node with head `Association` whose arguments are
 * two-argument `Rule[k, v]` nodes.  Keeping associations as plain expressions
 * means the whole "everything is an expression" toolchain — expr_copy,
 * expr_free, expr_eq, expr_hash, pattern matching, ReplaceAll, FullForm — all
 * keep working with no special cases.
 *
 * Keys are unique and insertion-ordered: the first occurrence fixes a key's
 * position, the last occurrence its value (matching Wolfram semantics).
 *
 * Performance.  The naive way to build/query an association is an O(n) linear
 * scan per key, giving O(n^2) construction and O(n*m) bulk lookups.  Instead
 * this module builds a transient open-addressing hash index (keyed by
 * expr_hash / expr_eq) so that construction, de-duplication, Merge, KeyDrop,
 * bulk Lookup, Counts and GroupBy are all amortised O(n) — see assoc.c.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* True when `e` is an Association[...] node (head == SYM_Association). */
bool is_association(const Expr* e);

/* Build a canonical Association from `count` Rule[k,v] nodes.  The rules are
 * copied (the caller keeps ownership of `rules`).  Duplicate keys collapse
 * with last-value-wins while preserving first-occurrence order.  O(count). */
Expr* assoc_from_rules(Expr** rules, size_t count);

/* Register every Association-family builtin, with attributes and docstrings. */
void assoc_init(void);

#endif /* MATHILDA_ASSOC_H */
