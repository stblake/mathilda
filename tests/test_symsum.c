/* test_symsum.c -- the binder-summary substrate (src/symsum.h): word layout,
 * the tag against a live AssocIndex, epoch validity, ss_seen gating,
 * invalidation (expr_invalidate_hash, expr_unshare), and epoch exhaustion. */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "symsum.h"
#include "assoc_index.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Summarise a tree bottom-up, as a walker does on frozen nodes. */
static void fill_all(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return;
    fill_all(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++) fill_all(e->data.function.args[i]);
    symsum_fill(e);
}

static SymbolDef* def_of(const char* name) { return intern_symbol_node(name); }

static void test_layout_and_fill(void) {
    ASSERT(sizeof(Expr) == 56);
    Expr* e = parse_expression("fq[uq, gq[vq, 3], \"s\"]");
    uint64_t b;
    ASSERT(!symsum_get(e, &b));                       /* fresh: no summary */
    symsum_register_binder(def_of("uq"));
    ASSERT(def_of("uq")->ss_bit >= 1 && def_of("uq")->ss_bit <= SYMSUM_NBINDER);
    fill_all(e);
    ASSERT(symsum_get(e, &b));
    ASSERT(b == symsum_def_bit(def_of("uq")));        /* only the binder contributes */
    Expr* g = e->data.function.args[1];
    ASSERT(symsum_get(g, &b) && b == 0);              /* gq[vq, 3]: no binder */
    ASSERT(def_of("vq")->ss_seen_epoch == (uint16_t)g_symsum_epoch);
    /* The word is odd and carries the epoch. */
    uint64_t w = (uint64_t)(uintptr_t)e->data.function.index;
    ASSERT((w & 1) && ((w >> 1) & 0xFFFF) == g_symsum_epoch);
    ASSERT(expr_fn_assoc_index(e) == NULL);
    expr_free(e);
}

static void test_epoch_gating(void) {
    Expr* e = parse_expression("hq[aq, bq]");
    fill_all(e);
    uint64_t b;
    ASSERT(symsum_get(e, &b) && b == 0);
    /* A name never seen by a summary: registering it costs no bump. */
    uint32_t ep = g_symsum_epoch;
    symsum_register_binder(def_of("neverseenq"));
    ASSERT(g_symsum_epoch == ep);
    ASSERT(symsum_get(e, &b));
    /* A name a current summary saw bitless: registering it invalidates. */
    symsum_register_binder(def_of("aq"));
    ASSERT(g_symsum_epoch == ep + 1);
    ASSERT(!symsum_get(e, &b));
    /* Refilled, it now carries aq's bit. */
    fill_all(e);
    ASSERT(symsum_get(e, &b) && b == symsum_def_bit(def_of("aq")));
    /* Idempotent registration. */
    uint8_t bit = def_of("aq")->ss_bit;
    symsum_register_binder(def_of("aq"));
    ASSERT(def_of("aq")->ss_bit == bit && g_symsum_epoch == ep + 1);
    /* The scoping heads share the SCOPE bit. */
    ASSERT(def_of("Module")->ss_bit == SYMSUM_SCOPE_SSBIT);
    ASSERT(symsum_def_bit(def_of("With")) == SYMSUM_SCOPE);
    ASSERT(symsum_def_bit(def_of("Function")) == SYMSUM_SCOPE);
    ASSERT(symsum_def_bit(def_of("Block")) == 0);
    expr_free(e);
}

static void test_invalidation(void) {
    Expr* e = parse_expression("kq[aq, {1, 2}]");
    fill_all(e);
    uint64_t b;
    ASSERT(symsum_get(e, &b));
    expr_invalidate_hash(e);
    ASSERT(!symsum_get(e, &b));
    fill_all(e);
    ASSERT(symsum_get(e, &b));
    /* Exclusive unshare: about to be mutated, so the summary goes. */
    ASSERT(e->refcount == 1);
    Expr* same = expr_unshare(e);
    ASSERT(same == e && !symsum_get(e, &b));
    /* Shared unshare: a fresh node with an empty slot; the original keeps its. */
    fill_all(e);
    Expr* extra = expr_copy(e);
    Expr* fresh = expr_unshare(e);
    ASSERT(fresh != e && fresh->data.function.index == NULL);
    ASSERT(symsum_get(extra, &b));
    expr_free(fresh);
    expr_free(extra);
}

static void test_assoc_slot(void) {
    /* An Association whose slot holds a summary still builds and uses its
     * index on a keyed read, and the summary never masks the pointer. */
    Expr* src = parse_expression("<|\"p\" -> aq, \"q\" -> 2|>");
    Expr* a = evaluate(src);
    expr_free(src);
    if (expr_fn_assoc_index(a)) {
        uint64_t b;
        symsum_fill(a);                                /* no-op over a live index */
        ASSERT(!symsum_get(a, &b));
    }
    expr_free(a);

    Expr* parsed = parse_expression("Association[Rule[\"p\", aq], Rule[\"q\", 2]]");
    fill_all(parsed);
    uint64_t b;
    ASSERT(symsum_get(parsed, &b) && b == symsum_def_bit(def_of("aq")));
    Expr* key = expr_new_string("q");
    extern Expr* assoc_lookup_value(const Expr* assoc, const Expr* key);
    Expr* v = assoc_lookup_value(parsed, key);
    ASSERT(v && v->type == EXPR_INTEGER && v->data.integer == 2);
    ASSERT(expr_fn_assoc_index(parsed) != NULL);       /* the index replaced the summary */
    ASSERT(!symsum_get(parsed, &b));
    v = assoc_lookup_value(parsed, key);
    ASSERT(v && v->data.integer == 2);
    expr_free(key);
    expr_free(parsed);                                 /* frees the index, not a summary */
}

/* Last: it switches the fast paths off for the rest of the process. */
static void test_epoch_exhaustion(void) {
    char buf[64];
    Expr* leaf = NULL;
    int n = 0;
    while (g_symsum_on && n < 70000) {
        snprintf(buf, sizeof buf, "wrapq%d", n++);
        leaf = expr_new_symbol(buf);
        (void)symsum_leaf_bits(leaf);                  /* seen in this epoch */
        symsum_register_binder(symsum_node_def(leaf)); /* -> bump */
        expr_free(leaf);
    }
    ASSERT(!g_symsum_on);
    ASSERT(g_symsum_epoch == SYMSUM_EPOCH_MAX);
}

int main(void) {
    symtab_init();
    core_init();
    if (!g_symsum_on) { printf("symsum disabled (MATHILDA_NO_SYMSUM / MTH_NO_SYMSUM): skipped\n"); return 0; }
    TEST(test_layout_and_fill);
    TEST(test_epoch_gating);
    TEST(test_invalidation);
    TEST(test_assoc_slot);
    TEST(test_epoch_exhaustion);
    printf("All symsum tests passed.\n");
    return 0;
}
