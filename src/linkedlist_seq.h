/*
 * linkedlist_seq.h -- one sequence reader over a List or a LinkedList
 * (linkedlist_seq.c). See the file comment there.
 */
#ifndef MATHILDA_LINKEDLIST_SEQ_H
#define MATHILDA_LINKEDLIST_SEQ_H

#include "linkedlist.h"

/* ---- sequence readers ----------------------------------------------------
 * The sorting, set, grouping and restructuring heads run ONE algorithm over a
 * List and a LinkedList, written against element INDICES: the algorithm reads
 * element i on demand (llq_peek), keeps indices (never an array of the
 * elements), and builds its answer from indices (llq_emit / llq_push), so a
 * LinkedList's elements are read one at a time straight from the slots and
 * its answer is copied slot to slot (a typed slot is never boxed to be
 * stored). A List lends its argument array, so the List path costs what it
 * did. A LinkedList is read through a snapshot (ll_copy, a LinkedList): code
 * the head runs per element may change the list. */
#define LLQ_RING 8
typedef struct {
    Expr* const* a;          /* a List's arguments (borrowed); NULL for a LinkedList */
    LinkedList*  snap;       /* a LinkedList: the snapshot read from (owned) */
    size_t       n;
    const Expr*  head;       /* a List's head; NULL for a LinkedList */
    Expr*        ring[LLQ_RING];   /* boxes of typed slots handed out by llq_peek */
    unsigned     ri;
} LLSeq;

/* False (q empty) when coll is neither an EXPR_FUNCTION nor a LinkedList. */
bool  llq_open(const Expr* coll, LLSeq* q);
void  llq_close(LLSeq* q);
static inline bool llq_linked(const LLSeq* q) { return q->snap != NULL; }
/* The storage kind of a LinkedList's snapshot; LL_BOXED for a List. */
LLKind llq_kind(const LLSeq* q);
/* Element i, BORROWED: a List's argument or a boxed slot of the snapshot
 * (valid until llq_close), or a box of a typed slot (valid for the next
 * LLQ_RING - 1 peeks). Never store it; store i. */
Expr* llq_peek(LLSeq* q, size_t i);
/* The container of q's kind holding elements idx[0..k) of q: head[...] of
 * copies for a List, a new LinkedList filled slot to slot for a LinkedList.
 * idx NULL means 0..k-1. */
Expr* llq_emit(LLSeq* q, const size_t* idx, size_t k);
/* Elements [lo, hi) of q as a container of q's kind (span-wise for a
 * LinkedList). */
Expr* llq_range(LLSeq* q, size_t lo, size_t hi);
/* Append element i of q to the LinkedList dst (slot to slot for a LinkedList,
 * a copy of the argument for a List). */
void  llq_push(LinkedList* dst, LLSeq* q, size_t i);
/* The container of q's kind holding computed items[0..n) (consumed; the array
 * is the caller's). */
Expr* llq_result(const LLSeq* q, Expr** items, size_t n);
/* A LinkedList (linked) or a List of items[0..n) (consumed). */
Expr* ll_seq_of(bool linked, Expr** items, size_t n);

/* ---- an index hash set over elements read on demand (SameQ) -------------
 * Holds INDICES; hashing and equality read the elements through `peek`
 * (borrowed element i of whatever indexed collection ctx names -- one LLSeq
 * via llq_peek_fn, or several operands numbered consecutively). llis_add: the
 * stored index whose element equals element i, else SIZE_MAX after inserting
 * i. llis_find: the stored index whose element equals `e`, or SIZE_MAX. */
typedef Expr* (*LLPeekFn)(void* ctx, size_t i);
Expr* llq_peek_fn(void* q, size_t i);         /* llq_peek as an LLPeekFn */
typedef struct { size_t* pos; uint64_t* hash; size_t mask; } LLIdxSet;
void   llis_init(LLIdxSet* s, size_t n);
void   llis_free(LLIdxSet* s);
size_t llis_add(LLIdxSet* s, LLPeekFn peek, void* ctx, size_t i);
size_t llis_find(const LLIdxSet* s, LLPeekFn peek, void* ctx, const Expr* e);

#endif /* MATHILDA_LINKEDLIST_SEQ_H */
