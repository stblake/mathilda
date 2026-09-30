/* ---------------------------------------------------------------------------
 * assoc_index.c — persistent key->position index for Associations.
 * See assoc_index.h for the design and the lifecycle contract.
 *
 * The table is an open-addressing (linear-probe) hash set of 1-based entry
 * positions (0 = empty slot) stored as uint32_t, so a slot costs 4 bytes: at
 * the < 0.5 load factor that is 8-16 bytes per entry.  Keys are BORROWED from
 * the entries — the index stores only positions, and every key's structural
 * hash is memoised on the key node itself (hash_cache), so rehashing on growth
 * and the backward-shift on removal cost no subtree walks.
 * -------------------------------------------------------------------------- */

#include "assoc_index.h"
#include "expr.h"
#include "sym_names.h"

#include <stdlib.h>

struct AssocIndex {
    uint32_t* pos;      /* pos[slot] = (entry index) + 1;  0 = empty */
    size_t    cap;      /* power of two */
    size_t    mask;     /* cap - 1 */
    size_t    n;        /* entries covered */
    size_t    args_cap; /* capacity of the owning node's args[] array */
    Expr* const* args;  /* the args[] array this index describes (identity only) */
};

/* Positions are stored as uint32_t (+1), so an index covers < 2^32 - 1 entries;
 * a larger association is simply left unindexed (scan fallback). */
#define AI_MAX_ENTRIES ((size_t)UINT32_MAX - 1)

/* Smallest power of two strictly greater than 2*n (load factor < 0.5), floor 8
 * — the same sizing as assoc.c's transient KeyIndex. */
static size_t ai_capacity_for(size_t n) {
    size_t want = (n < 4 ? 4 : n) * 2 + 1;
    size_t cap = 8;
    while (cap < want) cap <<= 1;
    return cap;
}

/* A well-formed entry: a two-argument Rule or RuleDelayed. */
static bool ai_is_entry(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

static const Expr* entry_key(const Expr* rule) {
    return rule->data.function.args[0];
}

/* Place entry `i` (known absent) into `pos`. */
static void ai_place(uint32_t* pos, size_t mask, Expr* const* entries, size_t i) {
    size_t slot = (size_t)expr_hash(entry_key(entries[i])) & mask;
    while (pos[slot] != 0) slot = (slot + 1) & mask;
    pos[slot] = (uint32_t)(i + 1);
}

AssocIndex* assoc_index_build(Expr* const* entries, size_t n) {
    if (n > AI_MAX_ENTRIES) return NULL;
    for (size_t i = 0; i < n; i++)
        if (!ai_is_entry(entries[i])) return NULL;
    AssocIndex* idx = malloc(sizeof(*idx));
    if (!idx) return NULL;
    idx->cap      = ai_capacity_for(n);
    idx->mask     = idx->cap - 1;
    idx->n        = n;
    idx->args_cap = n;
    idx->args     = entries;
    idx->pos      = calloc(idx->cap, sizeof(uint32_t));
    if (!idx->pos) { free(idx); return NULL; }

    /* Insert with a duplicate check: a key equal to one already placed means the
     * node is not canonical, so no index (see the certificate rule in the
     * header).  expr_eq runs only on a full-hash match, i.e. almost only on a
     * genuine duplicate. */
    for (size_t i = 0; i < n; i++) {
        const Expr* k = entry_key(entries[i]);
        uint64_t h = expr_hash(k);
        size_t slot = (size_t)h & idx->mask;
        while (idx->pos[slot] != 0) {
            const Expr* other = entry_key(entries[idx->pos[slot] - 1]);
            if (expr_hash(other) == h && expr_eq(other, k)) {
                assoc_index_free(idx);
                return NULL;
            }
            slot = (slot + 1) & idx->mask;
        }
        idx->pos[slot] = (uint32_t)(i + 1);
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

bool assoc_index_insert(AssocIndex* idx, Expr* const* entries, size_t pos) {
    if (!idx || pos != idx->n || pos >= AI_MAX_ENTRIES) return false;
    size_t need = pos + 1;
    if (need * 2 >= idx->cap) {
        /* Grow: double the table and re-place every entry.  Keys' hashes are
         * memoised, so this is a pass of cached loads — amortised O(1) per
         * insert. */
        size_t ncap = idx->cap << 1;
        uint32_t* np = calloc(ncap, sizeof(uint32_t));
        if (!np) return false;
        for (size_t i = 0; i < pos; i++) ai_place(np, ncap - 1, entries, i);
        free(idx->pos);
        idx->pos  = np;
        idx->cap  = ncap;
        idx->mask = ncap - 1;
    }
    ai_place(idx->pos, idx->mask, entries, pos);
    idx->n = need;
    return true;
}

void assoc_index_remove(AssocIndex* idx, Expr* const* entries,
                        size_t pos, uint64_t key_hash) {
    if (!idx || pos >= idx->n) return;
    uint32_t target = (uint32_t)(pos + 1);

    /* 1. Find the removed entry's slot by position (no key compare needed). */
    size_t s = (size_t)key_hash & idx->mask;
    while (idx->pos[s] != target) {
        if (idx->pos[s] == 0) return;          /* not found: index was stale */
        s = (s + 1) & idx->mask;
    }
    idx->pos[s] = 0;

    /* 2. Every later entry moved down one position. */
    for (size_t i = 0; i < idx->cap; i++)
        if (idx->pos[i] > target) idx->pos[i]--;
    idx->n--;

    /* 3. Backward-shift deletion: close the hole so every remaining key stays
     * reachable from its home slot without tombstones.  `entries` already
     * reflects the shifted positions, so pos values index it directly. */
    size_t i = s, j = s;
    for (;;) {
        j = (j + 1) & idx->mask;
        if (idx->pos[j] == 0) break;
        size_t h = (size_t)expr_hash(entry_key(entries[idx->pos[j] - 1])) & idx->mask;
        bool stays = (i <= j) ? (i < h && h <= j) : (i < h || h <= j);
        if (stays) continue;
        idx->pos[i] = idx->pos[j];
        idx->pos[j] = 0;
        i = j;
    }
}

bool assoc_index_matches(const AssocIndex* idx, Expr* const* args, size_t n) {
    return idx && idx->n == n && idx->args == args;
}

size_t assoc_index_count(const AssocIndex* idx)    { return idx ? idx->n : 0; }
size_t assoc_index_args_cap(const AssocIndex* idx) { return idx ? idx->args_cap : 0; }
void   assoc_index_set_args(AssocIndex* idx, Expr* const* args, size_t cap) {
    if (!idx) return;
    idx->args = args;
    idx->args_cap = cap;
}

void assoc_index_free(AssocIndex* idx) {
    if (!idx) return;
    free(idx->pos);
    free(idx);
}
