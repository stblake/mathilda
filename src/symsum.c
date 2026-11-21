/* symsum.c -- binder registry, epoch and switches for the per-node binder
 * summaries (symsum.h; docs/design/scoping_substitution.md). */
#include "symsum.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef MTH_NO_SYMSUM
/* The summary shares a pointer slot and is told apart by bit 0. */
typedef char symsum_assert_ptr64[(sizeof(void*) == 8 && sizeof(uintptr_t) == 8) ? 1 : -1];
int g_symsum_on = 1;
#else
int g_symsum_on = 0;
#endif
int      g_symsum_diff = 0;
uint32_t g_symsum_epoch = 1;

static unsigned g_next_bit = 0;      /* round robin over the 46 binder bits */
static uint64_t g_epoch_bumps = 0;

uint64_t symsum_epoch_bumps(void) { return g_epoch_bumps; }

/* Invalidate every summary.  On exhaustion of the epoch field the fast paths
 * are switched off for the session: an old summary could otherwise alias a
 * new epoch and be trusted. */
static void symsum_bump(void) {
    g_epoch_bumps++;
    if (g_symsum_epoch >= SYMSUM_EPOCH_MAX) {
        if (g_symsum_on)
            fprintf(stderr, "symsum: epoch exhausted after %llu binder registrations; "
                            "scoping fast paths off for this session\n",
                    (unsigned long long)g_epoch_bumps);
        g_symsum_on = 0;
        return;
    }
    g_symsum_epoch++;
}

void symsum_register_binder(SymbolDef* d) {
    if (!d || d->ss_bit) return;
    d->ss_bit = (uint8_t)(g_next_bit % SYMSUM_NBINDER + 1);
    g_next_bit++;
    /* A summary of this epoch may have seen d as a bitless leaf. */
    if (d->ss_seen_epoch == (uint16_t)g_symsum_epoch) symsum_bump();
}

static void scope_head(const char* name) {
    SymbolDef* d = intern_symbol_node(name);
    if (!d) return;
    if (d->ss_bit != SYMSUM_SCOPE_SSBIT && d->ss_seen_epoch == (uint16_t)g_symsum_epoch)
        symsum_bump();
    d->ss_bit = SYMSUM_SCOPE_SSBIT;
}

#ifdef MTH_SCOPE_STATS
SymsumStats g_scope_stats;
static void symsum_report(void) {
    const SymsumStats* t = &g_scope_stats;
    fprintf(stderr, "\n==== scoping walks (MTH_SCOPE_STATS) ====\n"
        "substitute_scoping: walks %llu nodes %llu skipped %llu filled %llu\n"
        "replace_bindings:   walks %llu nodes %llu skipped %llu filled %llu renames %llu\n"
        "danger set:         built %llu avoided by summaries %llu\n"
        "pure functions:     walks %llu nodes %llu skipped %llu filled %llu\n"
        "epoch %u bumps %llu fast paths %s\n",
        (unsigned long long)t->ss_top, (unsigned long long)t->ss_nodes,
        (unsigned long long)t->ss_skip, (unsigned long long)t->ss_fill,
        (unsigned long long)t->rb_top, (unsigned long long)t->rb_nodes,
        (unsigned long long)t->rb_skip, (unsigned long long)t->rb_fill,
        (unsigned long long)t->rb_renames,
        (unsigned long long)t->dg_builds, (unsigned long long)t->dg_avoided,
        (unsigned long long)t->pf_top, (unsigned long long)t->pf_nodes,
        (unsigned long long)t->pf_skip, (unsigned long long)t->pf_fill,
        g_symsum_epoch, (unsigned long long)g_epoch_bumps, g_symsum_on ? "on" : "off");
}
#endif

void symsum_init(void) {
#ifdef MTH_SCOPE_STATS
    atexit(symsum_report);
#endif
    const char* e = getenv("MATHILDA_NO_SYMSUM");
    if (e && *e == '1') g_symsum_on = 0;
    e = getenv("MATHILDA_SYMSUM_DIFF");
    g_symsum_diff = (e && *e == '1');
    /* The heads expr_is_binding_scope recognises share the SCOPE bit, so a
     * rule application can skip a subtree only if it holds no scope a
     * capture-avoiding rename could touch (match.c, rb_rec). */
    scope_head(SYM_Module);
    scope_head(SYM_With);
    scope_head(SYM_Function);
}

#ifdef MATHILDA_SYMSUM_VERIFY
/* The exact bits of a subtree, ignoring every stored summary. */
static uint64_t symsum_exact(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_SYMBOL) return symsum_sym_bit(e);
    if (e->type != EXPR_FUNCTION) return 0;
    uint64_t acc = symsum_exact(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        acc |= symsum_exact(e->data.function.args[i]);
    return acc;
}

/* A stale summary -- one missing a bit its subtree now has -- is the one way
 * a skip could be wrong; fail loudly on it. */
void symsum_verify_hit(const Expr* e, uint64_t bits) {
    uint64_t x = symsum_exact(e);
    if (x & ~bits) {
        fprintf(stderr, "symsum VERIFY: stale summary on %p: stored %016llx exact %016llx\n",
                (const void*)e, (unsigned long long)bits, (unsigned long long)x);
        abort();
    }
}
#endif
