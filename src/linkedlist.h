/*
 * linkedlist.h -- LinkedList: a mutable, reference-semantics sequence object.
 *
 * Design: LINKED_LIST_DESIGN.md (repo root). Summary:
 *
 *   - An EXPR_LINKEDLIST atom wraps a reference-counted `LinkedList` payload.
 *     Assignment and argument passing SHARE the payload; mutators (AppendTo,
 *     PrependTo, Push, Pop, PopFirst, Part assignment) change it in place.
 *     expr_eq / expr_hash / expr_compare use the payload's stable id, so a
 *     mutation never invalidates a cached hash anywhere.
 *
 *   - Storage is a deque of two "runs" (front grows on prepend, back grows on
 *     append). Each run is an unrolled linked list of chunks whose capacities
 *     double: chunk k holds LL_B << k slots, located through a fixed directory,
 *     so slot j of a run is found with one count-leading-zeros (no pointer
 *     chase). Append/prepend/pop at either end and index/set are all O(1);
 *     elements never move on growth.
 *
 *   - While every element is a machine number of one kind the chunks hold raw
 *     int64_t / double / (re, im) doubles (LL_I64 / LL_REAL / LL_CPLX); the
 *     first nonconforming element converts the list to LL_BOXED (owned Expr*)
 *     once. Typed storage is what lets Normal[] hand back a packed buffer and
 *     Compile[] read/write raw slots through these same helpers.
 *
 * Ownership conventions (match the rest of the kernel):
 *   - Functions taking `Expr* e` documented "consumes" take ownership of e.
 *   - Functions returning `Expr*` return a new reference the caller owns.
 *   - Indices are 0-based and must be in range (callers validate); the public
 *     builtins translate Part's 1-based / negative specs before calling in.
 */
#ifndef MATHILDA_LINKEDLIST_H
#define MATHILDA_LINKEDLIST_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

typedef struct LinkedList LinkedList;

/* Element storage kind. LL_EMPTY only while the list has never held an
 * element (or was emptied and is retyped by the next push). */
typedef enum {
    LL_EMPTY = 0,
    LL_I64,     /* int64_t slots: machine Integers            */
    LL_REAL,    /* double slots: machine Reals                */
    LL_CPLX,    /* 2x double slots: Complex[Real, Real]       */
    LL_BOXED    /* Expr* slots, each owned                    */
} LLKind;

/* ---- lifecycle -------------------------------------------------------- */
LinkedList* ll_new(void);
LinkedList* ll_retain(LinkedList* ll);       /* returns ll                      */
void        ll_release(LinkedList* ll);      /* frees at refcount 0             */
uint64_t    ll_id(const LinkedList* ll);     /* stable identity (hash / order)  */

/* Shallow copy: an independent list with the same elements and kind (boxed
 * elements shared by Expr refcount). O(n). */
LinkedList* ll_copy(const LinkedList* ll);

/* ---- inspection ------------------------------------------------------- */
int64_t ll_length(const LinkedList* ll);
LLKind  ll_kind(const LinkedList* ll);
/* Bytes held by the payload (header + chunk buffers; not boxed elements). */
int64_t ll_payload_bytes(const LinkedList* ll);

/* The storage kind an element would have on its own. */
LLKind  ll_classify(const Expr* e);

/* ---- element access (0-based, i in [0, length)) ----------------------- */
Expr* ll_get(const LinkedList* ll, int64_t i);          /* new reference   */
/* Borrowed pointer to a boxed element (LL_BOXED only), else NULL. */
const Expr* ll_peek_boxed(const LinkedList* ll, int64_t i);
/* Typed reads for the compiled path; false when the list is not that kind. */
bool ll_get_i64 (const LinkedList* ll, int64_t i, int64_t* out);
bool ll_get_real(const LinkedList* ll, int64_t i, double* out);
bool ll_get_cplx(const LinkedList* ll, int64_t i, double* re, double* im);
void ll_set(LinkedList* ll, int64_t i, Expr* e);        /* consumes e      */

/* ---- growth and removal ---------------------------------------------- */
void  ll_push_back (LinkedList* ll, Expr* e);           /* consumes e      */
void  ll_push_front(LinkedList* ll, Expr* e);           /* consumes e      */
void  ll_push_back_i64 (LinkedList* ll, int64_t v);
void  ll_push_back_real(LinkedList* ll, double v);
void  ll_push_back_cplx(LinkedList* ll, double re, double im);
void  ll_push_front_i64 (LinkedList* ll, int64_t v);
void  ll_push_front_real(LinkedList* ll, double v);
void  ll_push_front_cplx(LinkedList* ll, double re, double im);
/* Remove and return the last / first element. NULL when empty. */
Expr* ll_pop_back (LinkedList* ll);
Expr* ll_pop_front(LinkedList* ll);

/* Append every element of a List, packed List or (any-rank) NDArray `src`
 * (borrowed). Rank-1 int64/float64/complex64 buffers are copied chunk-wise.
 * Returns false (list unchanged) when `src` is none of those. */
bool  ll_append_elements(LinkedList* ll, const Expr* src);
/* True when `src` is something ll_append_elements accepts. */
bool  ll_is_list_source(const Expr* src);

/* The elements as a List: a packed List when the storage is typed and the
 * result is large enough to pack, otherwise a plain List. New reference. */
Expr* ll_to_list(const LinkedList* ll);

/* Would storing `e` (borrowed) into `target` make `target` reachable from
 * itself -- a reference cycle no refcount could free? Mutators reject it. */
bool  ll_would_cycle(const LinkedList* target, const Expr* e);

/* ---- Expr-level helpers ----------------------------------------------- */
/* The payload of an EXPR_LINKEDLIST node, or NULL for anything else. */
LinkedList* expr_linkedlist(const Expr* e);

/* Registers LinkedList, Push, Pop and PopFirst. */
void linkedlist_init(void);

#endif /* MATHILDA_LINKEDLIST_H */
