/*
 * linkedlist_seq.c -- one sequence reader over a List or a LinkedList.
 *
 * The sorting, set, grouping, counting and restructuring heads (Sort, Union,
 * Tally, Gather, Partition, Split, Riffle, PadLeft, ...) each have ONE
 * algorithm, written against element indices of an LLSeq:
 *
 *   - a List lends its argument array: llq_peek(q, i) is args[i], borrowed,
 *     and the List path costs what it always did;
 *   - a LinkedList is read through a snapshot (ll_copy -- a LinkedList, not a
 *     List: code the head runs per element may change the original), one
 *     element at a time on demand, straight from the slots: a boxed slot is
 *     lent as it is, a typed slot is boxed into a small ring of temporaries.
 *
 * No array of the elements is ever built -- an algorithm that needs random
 * access keeps INDICES (or computed keys) -- and answers are assembled from
 * indices (llq_emit / llq_range / llq_push: slot to slot, so a typed slot is
 * never boxed to be stored). The conversion counter therefore never moves.
 * Machine-number LinkedLists whose answer has a packed raw kernel take that
 * kernel on an ll_gather scratch buffer instead; see the individual heads.
 */
#include "linkedlist_seq.h"

#include <stdio.h>
#include <stdlib.h>

#include "sym_names.h"

bool llq_open(const Expr* coll, LLSeq* q) {
    q->a = NULL; q->snap = NULL; q->n = 0; q->head = NULL; q->ri = 0;
    for (int k = 0; k < LLQ_RING; k++) q->ring[k] = NULL;
    if (!coll) return false;
    if (coll->type == EXPR_FUNCTION) {
        q->a = coll->data.function.args;
        q->n = coll->data.function.arg_count;
        q->head = coll->data.function.head;
        return true;
    }
    const LinkedList* ll = expr_linkedlist(coll);
    if (!ll) return false;
    q->snap = ll_copy(ll);
    q->n = (size_t)ll_length(q->snap);
    return true;
}

void llq_close(LLSeq* q) {
    for (int k = 0; k < LLQ_RING; k++) {
        if (q->ring[k]) expr_free(q->ring[k]);
        q->ring[k] = NULL;
    }
    if (q->snap) ll_release(q->snap);
    q->snap = NULL;
    q->a = NULL;
    q->n = 0;
}

LLKind llq_kind(const LLSeq* q) { return q->snap ? ll_kind(q->snap) : LL_BOXED; }

Expr* llq_peek(LLSeq* q, size_t i) {
    if (!q->snap) return q->a[i];
    const Expr* b = ll_peek_boxed(q->snap, (int64_t)i);
    if (b) return (Expr*)b;
    unsigned k = q->ri++ % LLQ_RING;
    if (q->ring[k]) expr_free(q->ring[k]);
    q->ring[k] = ll_get(q->snap, (int64_t)i);
    return q->ring[k];
}

void llq_push(LinkedList* dst, LLSeq* q, size_t i) {
    if (q->snap) ll_push_back_from(dst, q->snap, (int64_t)i);
    else         ll_push_back(dst, expr_copy(q->a[i]));
}

Expr* llq_emit(LLSeq* q, const size_t* idx, size_t k) {
    if (q->snap) {
        LinkedList* out = ll_new();
        for (size_t t = 0; t < k; t++) ll_push_back_from(out, q->snap, (int64_t)(idx ? idx[t] : t));
        return expr_new_linkedlist(out);
    }
    Expr** items = malloc(sizeof(Expr*) * (k ? k : 1));
    if (!items) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    for (size_t t = 0; t < k; t++) items[t] = expr_copy(q->a[idx ? idx[t] : t]);
    Expr* r = expr_new_function(expr_copy((Expr*)q->head), items, k);
    free(items);
    return r;
}

Expr* llq_range(LLSeq* q, size_t lo, size_t hi) {
    if (q->snap) {
        LinkedList* out = ll_new();
        if (hi > lo) ll_append_range(out, q->snap, (int64_t)lo, (int64_t)hi);
        return expr_new_linkedlist(out);
    }
    size_t k = hi > lo ? hi - lo : 0;
    Expr** items = malloc(sizeof(Expr*) * (k ? k : 1));
    if (!items) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    for (size_t t = 0; t < k; t++) items[t] = expr_copy(q->a[lo + t]);
    Expr* r = expr_new_function(expr_copy((Expr*)q->head), items, k);
    free(items);
    return r;
}

Expr* llq_result(const LLSeq* q, Expr** items, size_t n) {
    if (!q->snap) return expr_new_function(expr_copy((Expr*)q->head), items, n);
    return ll_seq_of(true, items, n);
}

Expr* ll_seq_of(bool linked, Expr** items, size_t n) {
    if (!linked) return expr_new_function(expr_new_symbol(SYM_List), items, n);
    LinkedList* out = ll_new();
    for (size_t i = 0; i < n; i++) ll_push_back(out, items[i]);
    return expr_new_linkedlist(out);
}

/* ------------------------------------------------------------ index sets */

/* Two elements match when their hashes are equal AND expr_eq holds: the
 * hash is part of the identity, as for the List heads' hash tables (whose
 * buckets kept 0. and -0. apart -- equal under expr_eq, distinct hashes --
 * and as the packed kernels, which key on the bit pattern). */
void llis_init(LLIdxSet* s, size_t n) {
    size_t cap = 16;
    while (cap < 2 * n + 2) cap <<= 1;
    s->pos = calloc(cap, sizeof(size_t));
    s->hash = malloc(cap * sizeof(uint64_t));
    if (!s->pos || !s->hash) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    s->mask = cap - 1;
}

void llis_free(LLIdxSet* s) { free(s->pos); free(s->hash); s->pos = NULL; s->hash = NULL; }

Expr* llq_peek_fn(void* q, size_t i) { return llq_peek((LLSeq*)q, i); }

/* `e` is held for the probe: a ring box of the same sequence could otherwise
 * be recycled by the peeks of a long probe chain. */
size_t llis_find(const LLIdxSet* s, LLPeekFn peek, void* ctx, const Expr* e0) {
    Expr* e = expr_copy((Expr*)e0);
    uint64_t h = expr_hash(e);
    size_t slot = (size_t)h & s->mask, out = SIZE_MAX;
    while (s->pos[slot] != 0) {
        size_t j = s->pos[slot] - 1;
        if (s->hash[slot] == h && expr_eq(peek(ctx, j), e)) { out = j; break; }
        slot = (slot + 1) & s->mask;
    }
    expr_free(e);
    return out;
}

size_t llis_add(LLIdxSet* s, LLPeekFn peek, void* ctx, size_t i) {
    Expr* e = expr_copy(peek(ctx, i));
    uint64_t h = expr_hash(e);
    size_t slot = (size_t)h & s->mask, out = SIZE_MAX;
    while (s->pos[slot] != 0) {
        size_t j = s->pos[slot] - 1;
        if (s->hash[slot] == h && expr_eq(peek(ctx, j), e)) { out = j; break; }
        slot = (slot + 1) & s->mask;
    }
    if (out == SIZE_MAX) { s->pos[slot] = i + 1; s->hash[slot] = h; }
    expr_free(e);
    return out;
}
