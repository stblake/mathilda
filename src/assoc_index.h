#ifndef MATHILDA_ASSOC_INDEX_H
#define MATHILDA_ASSOC_INDEX_H

#include <stddef.h>
#include <stdint.h>

/* ---------------------------------------------------------------------------
 * assoc_index — a persistent key->position hash index cached on a canonical
 * Association node.
 *
 * Mathematica's Associations are hash maps with amortised O(1) single-key
 * lookup.  Mathilda stores an association as Association[Rule[k,v], ...] (see
 * assoc.h); historically every single-key Lookup / KeyExistsQ / Part was an
 * O(n) linear scan because the transient KeyIndex built for bulk operations
 * was thrown away.  This module keeps an index alive: it hangs off the
 * association node (in the free slot of the EXPR_FUNCTION union arm, so
 * sizeof(Expr) is unchanged), shared via the node's refcount, and freed with
 * the node.
 *
 * The index is built LAZILY, by the first single-key read (assoc_lookup_value
 * in assoc.c).  Eager construction at canonicalisation does not survive:
 * evaluate()'s fixed-point step keeps the ORIGINAL association node and frees
 * the rebuilt one an eager index would attach to, so the index would be
 * discarded before any reader saw it.  The reader instead attaches to the node
 * that actually survives, and evaluate()'s in-loop timestamp short-circuit keeps
 * that node stable across a Do/Table/Sum loop, so the O(n) build happens once
 * and every later probe is O(1).
 *
 * Correctness rests on canonical associations being immutable in place — every
 * update rebuilds a fresh node (assoc_entry_with_value / assoc_from_rules) — so
 * a cached key->position map can never go stale.  The index is pure
 * acceleration metadata: expr_eq / expr_hash / expr_compare never read it, and
 * a physical copy (expr_unshare) resets it to NULL rather than aliasing it.
 * The lazy build mutates a `const`/shared node, which is safe in the
 * single-threaded interpreter (the same discipline as last_evaluated_at); the
 * parallel compiled evaluator must pre-build the index at its marshalling
 * boundary so a shared association never reaches the lazy path from a worker
 * thread. -------------------------------------------------------------------- */

struct Expr;                       /* opaque here; the .c includes expr.h */
typedef struct AssocIndex AssocIndex;

/* Build an index over `n` canonical entries — each a two-argument Rule whose
 * key is args[0].  Returns NULL for n == 0 or on allocation failure; callers
 * MUST treat a NULL index as "not indexed" and fall back to a linear scan,
 * which is always correct (just O(n)). */
AssocIndex* assoc_index_build(struct Expr* const* entries, size_t n);

/* Entry position of `key` among the entries, or -1 if absent.  `entries` MUST
 * be the same array the index was built over (the association's args). */
int64_t assoc_index_lookup(const AssocIndex* idx, struct Expr* const* entries,
                           const struct Expr* key);

void assoc_index_free(AssocIndex* idx);

#endif /* MATHILDA_ASSOC_INDEX_H */
