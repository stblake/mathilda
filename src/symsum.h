/* symsum.h -- per-node binder summaries for the scoping walkers.
 *
 * Every substituting walk -- Module/With renaming (substitute_scoping), rule
 * application (replace_bindings), pure-function application -- used to visit
 * every node of the tree it was given, though on the slow integration cases
 * 97% of those nodes come back unchanged.  A summary lets a walker skip a
 * subtree that provably contains none of the names it substitutes.
 *
 * WHAT A SUMMARY IS.  A 64-bit word stored in a FUNCTION node's
 * `data.function.index` slot (the slot an Association uses for its AssocIndex
 * pointer; that pointer is malloc-aligned, so the low bit tells them apart):
 *
 *     bit  0       tag: 1 = summary (0 = AssocIndex pointer or empty)
 *     bits 1..16   epoch the summary was computed in
 *     bits 17..62  46 binder bits
 *     bit  63      SCOPE: the subtree contains a Module / With / Function head
 *
 * Bits are given only to BINDER names -- names some Module/With/Function
 * binds, or some rule binds as a pattern variable -- plus the SCOPE bit, which
 * the three scoping heads share.  Every other symbol contributes nothing, so a
 * large data value (a context Association, a matrix) summarises to almost
 * nothing and is skipped at its root.  A summary is a one-hash Bloom filter
 * over the subtree's binder symbols, heads included: a bit absent means the
 * name is absent, so a skip is exact; a collision only costs a descent.
 *
 * WHEN IT IS VALID.  A symbol summarised before it became a binder is missing
 * from that summary.  Registering a binder therefore advances the global
 * epoch -- invalidating every summary -- when the symbol was seen as a bitless
 * leaf by a summary of the CURRENT epoch (SymbolDef.ss_seen_epoch).  A name
 * registered before any summary sees it, the usual case (pattern variables are
 * registered when the rule is stored), costs nothing.  If the 16-bit epoch is
 * ever exhausted the fast paths turn off for the rest of the session.
 *
 * WHERE IT IS FILLED.  Only on FROZEN nodes: nodes reached through a node with
 * refcount >= 2 (immutable, with everything below it, by the copy-on-write
 * discipline) or belonging to a stored rule.  A fresh node gets none.  The
 * summary is computed bottom-up from the children's, during the walk that is
 * happening anyway, and only when every FUNCTION child already has a valid
 * one.  expr_invalidate_hash and expr_unshare clear it.
 *
 * SWITCHES.  MATHILDA_NO_SYMSUM=1 turns every consumer off (the A/B baseline
 * and the field fallback); MATHILDA_SYMSUM_DIFF=1 re-runs each fast walk
 * slowly and aborts on any difference; -DMATHILDA_SYMSUM_VERIFY recomputes
 * every summary a walker relies on; -DMTH_NO_SYMSUM compiles the consumers
 * out.  See docs/design/scoping_substitution.md. */
#ifndef SYMSUM_H
#define SYMSUM_H

#include <stdint.h>
#include <stdbool.h>
#include "expr.h"
#include "symtab.h"

/* The tagged slot needs 64-bit pointers.  Elsewhere the consumers compile out. */
#if !defined(MTH_NO_SYMSUM) && (!defined(UINTPTR_MAX) || UINTPTR_MAX < UINT64_C(0xFFFFFFFFFFFFFFFF))
#define MTH_NO_SYMSUM 1
#endif

#define SYMSUM_TAG          UINT64_C(1)
#define SYMSUM_EPOCH_SHIFT  1
#define SYMSUM_EPOCH_MAX    0xFFFFu
#define SYMSUM_BIT_BASE     16          /* ss_bit k (1..47) is word bit 16 + k */
#define SYMSUM_NBINDER      46          /* ss_bit 1..46: binder names */
#define SYMSUM_SCOPE_SSBIT  47          /* ss_bit 47 = bit 63: the scoping heads */
#define SYMSUM_SCOPE        (UINT64_C(1) << 63)
#define SYMSUM_BITS         (~UINT64_C(0) << 17)

/* 1 while the fast paths run: not disabled by MATHILDA_NO_SYMSUM, not killed
 * by epoch exhaustion. */
extern int      g_symsum_on;
extern int      g_symsum_diff;     /* MATHILDA_SYMSUM_DIFF=1 */
extern uint32_t g_symsum_epoch;    /* 1..SYMSUM_EPOCH_MAX */

void symsum_init(void);
/* Give `d` a binder bit if it has none (advancing the epoch when a current
 * summary has already seen it bitless).  Idempotent. */
void symsum_register_binder(SymbolDef* d);
/* Statistics for MATHILDA_DISPATCH_STATS-style reports and tests. */
uint64_t symsum_epoch_bumps(void);

/* The bit a symbol contributes to a MASK (no side effects): 0 for a non-binder. */
static inline uint64_t symsum_def_bit(const SymbolDef* d) {
    return (d && d->ss_bit) ? (UINT64_C(1) << (SYMSUM_BIT_BASE + d->ss_bit)) : 0;
}
/* The SymbolDef of a symbol node WITHOUT materialising a definition: the
 * node's cached def, else the interned node (not cached -- symtab_node_def's
 * cache must keep meaning "materialised"). */
static inline SymbolDef* symsum_node_def(const Expr* sym) {
    SymbolDef* d = sym->data.symbol.def;
    return d ? d : intern_symbol_node(sym->data.symbol.name);
}
static inline uint64_t symsum_sym_bit(const Expr* sym) {
    return symsum_def_bit(symsum_node_def(sym));
}

/* The bits a symbol LEAF contributes to a SUMMARY being filled; a bitless
 * leaf is remembered as seen in this epoch (see symsum_register_binder). */
static inline uint64_t symsum_leaf_bits(const Expr* sym) {
    SymbolDef* d = symsum_node_def(sym);
    if (!d) return 0;
    if (d->ss_bit) return UINT64_C(1) << (SYMSUM_BIT_BASE + d->ss_bit);
    d->ss_seen_epoch = (uint16_t)g_symsum_epoch;
    return 0;
}

/* Walk counters, compiled in only with -DMTH_SCOPE_STATS (report at exit):
 * nodes visited, subtrees skipped by a summary, summaries filled, per walker. */
#ifdef MTH_SCOPE_STATS
typedef struct {
    uint64_t ss_top, ss_nodes, ss_skip, ss_fill;
    uint64_t rb_top, rb_nodes, rb_skip, rb_fill, rb_renames;
    uint64_t dg_builds, dg_avoided;
    uint64_t pf_top, pf_nodes, pf_skip, pf_fill;
} SymsumStats;
extern SymsumStats g_scope_stats;
#define SCOPE_STAT(f) (g_scope_stats.f++)
#else
#define SCOPE_STAT(f) ((void)0)
#endif

#ifdef MATHILDA_SYMSUM_VERIFY
void symsum_verify_hit(const Expr* e, uint64_t bits);
#endif

/* A FUNCTION node's valid summary bits, or false. */
static inline bool symsum_get(const Expr* e, uint64_t* bits) {
    uint64_t w = (uint64_t)(uintptr_t)e->data.function.index;
    if (!(w & SYMSUM_TAG)) return false;
    if (((w >> SYMSUM_EPOCH_SHIFT) & SYMSUM_EPOCH_MAX) != g_symsum_epoch) return false;
    *bits = w & SYMSUM_BITS;
#ifdef MATHILDA_SYMSUM_VERIFY
    symsum_verify_hit(e, *bits);
#endif
    return true;
}

/* Store a summary on a FUNCTION node, unless the slot holds an AssocIndex. */
static inline void symsum_put(Expr* e, uint64_t bits) {
    uintptr_t w = (uintptr_t)e->data.function.index;
    if (w && !(w & 1)) return;
    uint64_t s = SYMSUM_TAG | ((uint64_t)g_symsum_epoch << SYMSUM_EPOCH_SHIFT) | (bits & SYMSUM_BITS);
    e->data.function.index = (struct AssocIndex*)(uintptr_t)s;
}

/* Fill `e`'s summary from its children's, if every FUNCTION child has a
 * valid one.  Call only on a frozen node the walk found unchanged. */
static inline void symsum_fill(Expr* e) {
    uint64_t acc = 0;
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i <= n; i++) {
        const Expr* c = i == 0 ? e->data.function.head : e->data.function.args[i - 1];
        if (!c) return;
        if (c->type == EXPR_SYMBOL) acc |= symsum_leaf_bits(c);
        else if (c->type == EXPR_FUNCTION) {
            uint64_t b;
            if (!symsum_get(c, &b)) return;
            acc |= b;
        }
    }
    symsum_put(e, acc);
}

#endif /* SYMSUM_H */
