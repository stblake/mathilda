/*
 * linkedlist_seq.c -- one element-array view over a List or a LinkedList.
 *
 * The sorting, set, grouping, counting and restructuring heads (Sort, Union,
 * Tally, Gather, Partition, Split, Riffle, PadLeft, Accumulate, Max, ...) each
 * have one algorithm over `(Expr* const* elems, size_t n)`. LLView hands that
 * algorithm the elements of either container:
 *
 *   - a List (any EXPR_FUNCTION) lends its own argument array, borrowed: the
 *     List path is exactly what it was, no copy, no allocation;
 *   - a LinkedList lends an array of its elements read from the slots: the
 *     boxed slots by reference count (a snapshot -- a callback the head runs
 *     may change the list), typed slots as fresh boxes.
 *
 * No List of the LinkedList's elements is built, so the conversion counter
 * (Internal`LinkedListConversionCount[]) does not move, and both containers
 * run the same code and cannot drift. ll_view_result closes the loop: the
 * List head of the source for a List, a new LinkedList for a LinkedList (the
 * result-container rule: a head that returns a List for a List returns a
 * LinkedList for a LinkedList, sub-sequences included).
 *
 * Machine-number LinkedLists whose answer is arithmetic over a buffer take
 * the raw-buffer kernels the packed List reaches instead (ll_gather into a
 * scratch buffer, the ND kernel, ll_from_buffer); see the individual heads.
 */
#include "linkedlist_seq.h"

#include <stdio.h>
#include <stdlib.h>

#include "sym_names.h"

bool ll_view(const Expr* coll, LLView* v) {
    v->a = NULL; v->n = 0; v->head = NULL; v->owned = false;
    if (!coll) return false;
    if (coll->type == EXPR_FUNCTION) {
        v->a = coll->data.function.args;
        v->n = coll->data.function.arg_count;
        v->head = coll->data.function.head;
        return true;
    }
    const LinkedList* ll = expr_linkedlist(coll);
    if (!ll) return false;
    int64_t n = ll_length(ll);
    v->n = (size_t)n;
    v->owned = true;
    if (n == 0) return true;
    /* An array of every element -- borrowed refs of a boxed list as much as
     * boxes of a typed one -- is a List in all but its head: counted, so the
     * audit names every head still taking this path. Read on demand instead
     * (ll_seq_view over a snapshot), or a raw kernel for a typed list. */
    if (ll_kind(ll) == LL_BOXED) ll_note_materialised();
    v->a = malloc(sizeof(Expr*) * (size_t)n);
    if (!v->a) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    if (ll_kind(ll) == LL_BOXED) {
        ll_gather(ll, 0, n, v->a);                     /* borrowed ... */
        for (int64_t i = 0; i < n; i++) v->a[i] = expr_copy(v->a[i]);   /* ... held */
    } else {
        /* Boxing every typed slot into one array IS a List of the elements
         * in all but its head: counted, so the audit names every head that
         * still takes this path (a typed list wants a raw kernel or a span
         * copy instead). */
        ll_note_materialised();
        for (int64_t i = 0; i < n; i++) v->a[i] = ll_get(ll, i);
    }
    return true;
}

void ll_view_done(LLView* v) {
    if (v->owned) {
        for (size_t i = 0; i < v->n; i++) expr_free(v->a[i]);
        free(v->a);
    }
    v->a = NULL; v->n = 0; v->owned = false;
}

Expr* ll_view_result(const LLView* v, Expr** items, size_t n) {
    if (v->head) return expr_new_function(expr_copy((Expr*)v->head), items, n);
    return ll_seq_of(true, items, n);
}

Expr* ll_seq_of(bool linked, Expr** items, size_t n) {
    if (!linked) return expr_new_function(expr_new_symbol(SYM_List), items, n);
    LinkedList* out = ll_new();
    for (size_t i = 0; i < n; i++) ll_push_back(out, items[i]);
    return expr_new_linkedlist(out);
}
