/* ---------------------------------------------------------------------------
 * assoc_index.c — persistent key->position index for canonical Associations.
 * See assoc_index.h for the design and the immutability contract.
 *
 * The table is a fixed-capacity open-addressing hash set of 1-based entry
 * indices (0 = empty slot), sized once for the entry count and never rehashed —
 * the same discipline as the transient KeyIndex in assoc.c, but heap-allocated
 * and owned by the association node.  Keys are BORROWED from the entries; the
 * index stores only positions.
 * -------------------------------------------------------------------------- */

#include "assoc_index.h"
#include "expr.h"

#include <stdlib.h>

struct AssocIndex {
    size_t* pos;   /* pos[slot] = (entry index) + 1;  0 = empty */
    size_t  cap;   /* power of two */
    size_t  mask;  /* cap - 1 */
    size_t  n;     /* entries indexed (for reference / debugging) */
};

/* Smallest power of two strictly greater than 2*n (load factor < 0.5), floor 8
 * — identical sizing to assoc.c's KeyIndex so probe behaviour matches. */
static size_t ai_capacity_for(size_t n) {
    size_t want = (n < 4 ? 4 : n) * 2 + 1;
    size_t cap = 8;
    while (cap < want) cap <<= 1;
    return cap;
}

/* Key of a canonical entry: entries are two-argument Rule[k, v] nodes. */
static const Expr* entry_key(const Expr* rule) {
    return rule->data.function.args[0];
}

AssocIndex* assoc_index_build(Expr* const* entries, size_t n) {
    if (n == 0) return NULL;
    /* Guard against a non-canonical association (e.g. a user-typed
     * Association[1, 2]): the index keys off each entry's args[0], so every
     * entry must be a two-argument node.  Bail to a NULL index (the caller then
     * falls back to a linear scan, which stays correct). */
    for (size_t i = 0; i < n; i++)
        if (entries[i]->type != EXPR_FUNCTION || entries[i]->data.function.arg_count != 2)
            return NULL;
    AssocIndex* idx = malloc(sizeof(*idx));
    if (!idx) return NULL;
    idx->cap  = ai_capacity_for(n);
    idx->mask = idx->cap - 1;
    idx->n    = n;
    idx->pos  = calloc(idx->cap, sizeof(size_t));
    if (!idx->pos) { free(idx); return NULL; }

    /* Canonical entries have unique keys, so no key comparison is needed here:
     * the probe only steps past hash collisions between DISTINCT keys, and the
     * first-inserted entry sits earliest in probe order (first-occurrence
     * semantics, matching assoc_scan) should a duplicate ever slip through. */
    for (size_t i = 0; i < n; i++) {
        size_t slot = (size_t)expr_hash(entry_key(entries[i])) & idx->mask;
        while (idx->pos[slot] != 0) slot = (slot + 1) & idx->mask;
        idx->pos[slot] = i + 1;
    }
    return idx;
}

int64_t assoc_index_lookup(const AssocIndex* idx, Expr* const* entries,
                           const Expr* key) {
    if (!idx) return -1;
    size_t slot = (size_t)expr_hash(key) & idx->mask;
    while (idx->pos[slot] != 0) {
        size_t e = idx->pos[slot] - 1;
        if (expr_eq(entry_key(entries[e]), key)) return (int64_t)e;
        slot = (slot + 1) & idx->mask;
    }
    return -1;
}

void assoc_index_free(AssocIndex* idx) {
    if (!idx) return;
    free(idx->pos);
    free(idx);
}
