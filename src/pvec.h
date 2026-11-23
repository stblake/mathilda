/*
 * pvec.h -- PVec: a persistent (immutable, structurally shared) vector.
 *
 * A 32-way bitmapped trie with a separate tail leaf (the Clojure / Scala
 * Vector layout).  Every operation that "changes" a vector returns a new
 * vector that shares all untouched nodes with the old one, so a snapshot is
 * free to keep and the old value is never disturbed.
 *
 *   Part get          O(log32 n)   (<= 7 levels for 10^10 elements)
 *   Part set          O(log32 n)   path copy; O(1) when the vector is held once
 *   Append / push     O(1) amortised (one leaf copy per 32 pushes)
 *   Most / pop        O(1) amortised
 *   Length            O(1)
 *
 * Ownership ("consume and return"): the updating functions take the caller's
 * reference to `v` and return a vector the caller owns.  When `v` is held only
 * by the caller (refcount 1) the update happens IN PLACE, down to each node it
 * touches that is itself held once, so a loop `v = pv_push(v, x)` is a plain
 * amortised-O(1) append.  To keep the old value, pv_retain it first:
 *
 *     PVec* w = pv_push(pv_retain(v), x);   // v is unchanged, w shares with v
 *
 * Storage is typed: while every element is a machine Integer the leaves hold
 * raw int64_t, while every element is a machine Real they hold raw double;
 * the first element of any other kind converts the vector to boxed (owned
 * Expr*) leaves once.  All slots are one 8-byte word, so a leaf is 32 words.
 *
 * Elements are 0-based and callers validate the range (the Part-level
 * builtins translate Mathematica's 1-based / negative specs first).
 */
#ifndef MATHILDA_PVEC_H
#define MATHILDA_PVEC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

typedef struct PVec PVec;

typedef enum {
    PV_EMPTY = 0,   /* never held an element */
    PV_I64,         /* int64_t slots */
    PV_REAL,        /* double slots */
    PV_BOXED        /* Expr* slots, each owned */
} PVKind;

#define PV_BITS  5
#define PV_WIDTH (1 << PV_BITS)   /* 32 */

/* ---- lifecycle ---------------------------------------------------------- */
PVec*   pv_new(void);
PVec*   pv_retain(PVec* v);                 /* returns v */
void    pv_release(PVec* v);                /* frees at refcount 0 */
bool    pv_exclusive(const PVec* v);        /* held by exactly one reference */

/* ---- inspection --------------------------------------------------------- */
int64_t pv_length(const PVec* v);
PVKind  pv_kind(const PVec* v);
uint64_t pv_id(const PVec* v);              /* stable identity of THIS value */
/* Do `a` and `b` share their whole structure (same root and tail)? */
bool    pv_same(const PVec* a, const PVec* b);
/* Count of trie/tail nodes reachable from v, and of those held by more than
 * one reference (a measure of structural sharing, for tests). */
void    pv_node_stats(const PVec* v, int64_t* nodes, int64_t* shared);

/* Value equality / hash / order: same length and equal elements (typed and
 * boxed storage of equal values compare equal). */
bool     pv_equal(const PVec* a, const PVec* b);
uint64_t pv_hash(const PVec* v);
int      pv_compare(const PVec* a, const PVec* b);

/* ---- access ------------------------------------------------------------- */
Expr*   pv_get(const PVec* v, int64_t i);           /* new reference */
/* Typed reads; false when the vector is not of that kind. */
bool    pv_get_i64 (const PVec* v, int64_t i, int64_t* out);
bool    pv_get_real(const PVec* v, int64_t i, double* out);
/* The leaf cursor every scan uses: elements [i, i+k) occupy k consecutive
 * slots starting at *p (kind-sized words).  Returns k >= 1 (i in range). */
int64_t pv_span_at(const PVec* v, int64_t i, const void** p);

/* ---- updates (consume v, return the new vector) ------------------------- */
PVec*   pv_set (PVec* v, int64_t i, Expr* e);       /* consumes e */
PVec*   pv_push(PVec* v, Expr* e);                  /* consumes e */
PVec*   pv_push_i64 (PVec* v, int64_t x);
PVec*   pv_push_real(PVec* v, double x);
/* Remove the last element; *out (if non-NULL) receives it (new reference).
 * An empty vector is returned unchanged with *out = NULL. */
PVec*   pv_pop (PVec* v, Expr** out);

/* ---- conversion --------------------------------------------------------- */
PVec*   pv_from_elems(Expr* const* elems, size_t n);   /* borrows elems */
Expr*   pv_to_list(const PVec* v);                      /* plain List, new ref */

#endif /* MATHILDA_PVEC_H */
