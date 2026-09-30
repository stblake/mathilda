#ifndef MATHILDA_ASSOC_INDEX_H
#define MATHILDA_ASSOC_INDEX_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* ---------------------------------------------------------------------------
 * assoc_index — a persistent key->position hash index cached on an
 * Association node.
 *
 * Mathematica's Associations are hash maps with amortised O(1) single-key
 * lookup and update.  Mathilda stores an association as
 * Association[Rule[k,v], ...] (see assoc.h); this module keeps a key->position
 * hash table alive on the node itself (in the free slot of the EXPR_FUNCTION
 * union arm, so sizeof(Expr) is unchanged).  It is freed with the node and is
 * never aliased: a physical copy (expr_unshare) starts with no index.
 *
 * Lifecycle.
 *   - BUILT lazily by the first single-key read (assoc_lookup_value), by
 *     builtin_association's canonicality check, or by the in-place writers
 *     (assoc_write_path & co. in assoc.c) the first time they touch a node.
 *   - A successful build is a CERTIFICATE: every entry is a two-argument
 *     Rule/RuleDelayed and every key is distinct.  assoc_index_build returns
 *     NULL otherwise, and every caller treats NULL as "not indexed" and falls
 *     back to a linear scan, which is always correct.
 *   - UPDATED incrementally by the in-place writers, which are the only code
 *     allowed to rewrite a live association's args[] (and only when the node is
 *     uniquely referenced, refcount == 1): assoc_index_insert after an append,
 *     assoc_index_remove after a deletion; an overwrite keeps the key and the
 *     position, so the index needs no change.  Anything else that rewrites an
 *     association's args in place must free the index first.
 *   - `n` records how many entries the table covers; a reader whose node's
 *     arg_count disagrees treats the index as stale and rebuilds it.
 *
 * The index also records the node's args[] pointer and its capacity
 * (`args_cap`), so appends can grow args geometrically (amortised O(1)) rather
 * than reallocating by one slot each time.  args_cap >= arg_count always; a
 * freshly built index records the array it was built over with capacity equal
 * to the entry count (the array was allocated exact).  assoc_index_matches
 * compares both the pointer and the count, so an index is never trusted for a
 * node whose args[] something else replaced.
 *
 * The lazy build writes into a node that may be shared (refcount > 1); that is
 * benign metadata, safe in the single-threaded interpreter (the same
 * discipline as last_evaluated_at).  The parallel compiled evaluator pre-builds
 * the index at its marshalling boundary so a worker thread never builds one.
 * -------------------------------------------------------------------------- */

struct Expr;                       /* opaque here; the .c includes expr.h */
typedef struct AssocIndex AssocIndex;

/* Build an index over `n` entries.  Returns NULL when any entry is not a
 * two-argument Rule/RuleDelayed, when two entries share a key, or on allocation
 * failure (a NULL index means "not indexed": scan instead).  n == 0 yields a
 * valid empty index (the in-place writers need somewhere to record args_cap). */
AssocIndex* assoc_index_build(struct Expr* const* entries, size_t n);

/* Entry position of `key` among the entries, or -1 if absent.  `entries` MUST
 * be the array the index describes (the association's args). */
int64_t assoc_index_lookup(const AssocIndex* idx, struct Expr* const* entries,
                           const struct Expr* key);

/* Record that entries[pos] (pos == the old entry count) was appended with a
 * key not already present.  Grows the table when needed.  Returns false on
 * allocation failure, in which case the caller must free the index. */
bool assoc_index_insert(AssocIndex* idx, struct Expr* const* entries, size_t pos);

/* Remove the entry that USED to be at `pos`, whose key hash is `key_hash`,
 * after the caller has already shifted entries[pos+1..] down by one (so
 * `entries` now holds the n-1 remaining entries).  O(table size), no hashing
 * beyond the probe chain. */
void assoc_index_remove(AssocIndex* idx, struct Expr* const* entries,
                        size_t pos, uint64_t key_hash);

/* True when `idx` still describes the node whose args array is `args` with
 * `n` entries: the index remembers the array it was built over, so a node
 * whose args[] was swapped or resized by code outside the in-place writers is
 * detected and re-indexed instead of being read through a stale table. */
bool assoc_index_matches(const AssocIndex* idx, struct Expr* const* args, size_t n);

size_t assoc_index_count(const AssocIndex* idx);          /* entries covered */
/* Capacity of the node's args[] array (>= its arg_count). */
size_t assoc_index_args_cap(const AssocIndex* idx);
/* Record that the node's args[] is now `args`, with room for `cap` slots. */
void   assoc_index_set_args(AssocIndex* idx, struct Expr* const* args, size_t cap);

void assoc_index_free(AssocIndex* idx);

#endif /* MATHILDA_ASSOC_INDEX_H */
