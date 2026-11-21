/*
 * linkedlist_seq.h -- one element-array view over a List or a LinkedList
 * (linkedlist_seq.c). See the file comment there.
 */
#ifndef MATHILDA_LINKEDLIST_SEQ_H
#define MATHILDA_LINKEDLIST_SEQ_H

#include "linkedlist.h"

/* ---- element-array views (linkedlist_seq.c) ----------------------------
 * One algorithm over (elems, n) serves a List and a LinkedList alike: a List
 * (any EXPR_FUNCTION) lends its argument array (borrowed, nothing copied); a
 * LinkedList lends its elements read from the slots (boxed slots held by
 * refcount -- a snapshot, so a callback may change the list -- typed slots as
 * fresh boxes). No List is built. `head` is the List's head, NULL for a
 * LinkedList. */
typedef struct {
    Expr**      a;
    size_t      n;
    const Expr* head;
    bool        owned;
} LLView;
/* False (v empty) when coll is neither an EXPR_FUNCTION nor a LinkedList. */
bool  ll_view(const Expr* coll, LLView* v);
void  ll_view_done(LLView* v);
static inline bool ll_view_linked(const LLView* v) { return v->head == NULL; }
/* The container of the view's kind holding items[0..n) (consumed; the array
 * itself is the caller's): head[items...] for a List, a new LinkedList for a
 * LinkedList. */
Expr* ll_view_result(const LLView* v, Expr** items, size_t n);
/* A LinkedList (linked) or a List of items[0..n) (consumed). */
Expr* ll_seq_of(bool linked, Expr** items, size_t n);


#endif /* MATHILDA_LINKEDLIST_SEQ_H */
