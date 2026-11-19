/*
 * linkedlist.c -- LinkedList storage and the LinkedList[] constructor.
 *
 * See linkedlist.h for the contract and LINKED_LIST_DESIGN.md for the design.
 *
 * Storage layout
 * --------------
 * The logical sequence is   reverse(front run live slots) ++ back run live slots.
 *
 *     front run  (PrependTo writes slot front.end++)
 *     back  run  (AppendTo  writes slot back.end++)
 *
 * Each run is a chain of chunks with doubling capacity: chunk k holds
 * LL_B << k slots and covers run slots [LL_B (2^k - 1), LL_B (2^(k+1) - 1)).
 * `dir[k]` points at chunk k (the directory makes the "linked list" of chunks
 * randomly addressable): slot j lives in chunk floor(log2(j / LL_B + 1)) at
 * offset j - LL_B (2^k - 1). A run's live slots are [start, end); removing at
 * the run's outer end decrements `end`, removing at its inner end (the other
 * end of the whole list, once the opposite run is exhausted) increments
 * `start`.
 *
 * A queue (push at the back, pop at the front) walks `start` and `end` up
 * forever, so a dead prefix accumulates below `start`. When it is at least as
 * large as the live part, the run is REBASED: live slots move down to [0, n).
 * That costs O(live) after at least `live` pops, i.e. amortised O(1), and
 * keeps memory O(live). Chunks wholly past the end (beyond one spare) are
 * freed as the run shrinks.
 *
 * Elements never move on growth (only on a rebase or a typed->boxed
 * conversion), which is what makes append O(1) without the copy a realloc'd
 * array would pay.
 */
#include "linkedlist.h"

#include "common.h"       /* head_is, builtin_arg_error */
#include "eval.h"         /* eval_clock_tick */
#include "memacct.h"      /* MTH_MEM_CHARGE_BLOCK / MTH_MEM_RELEASE */
#include "message.h"
#include "ndarray.h"
#include "pack.h"         /* ndbuild_open, pack_offer */
#include "print.h"        /* expr_to_string */
#include "sym_names.h"
#include "symtab.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Capacity of chunk 0. Small so a short list costs little; doubling makes the
 * chunk count logarithmic, so the base barely matters for large lists. */
#define LL_B 8
/* Run slot indices are < 2^62, so chunk indices are < 60. */
#define LL_MAXCHUNK 60

typedef struct {
    unsigned char* dir[LL_MAXCHUNK];   /* chunk k: LL_B << k slots, or NULL */
    int64_t        start, end;          /* live slots [start, end)          */
} LLRun;

struct LinkedList {
    int64_t  refcount;
    uint64_t id;
    LLKind   kind;
    LLRun    front, back;
    int64_t  chunk_bytes;               /* bytes in allocated chunks        */
    uint64_t visit;                     /* cycle-check epoch (ll_would_cycle) */
};

static uint64_t g_ll_next_id = 1;

/* ------------------------------------------------------------------------
 *  Geometry
 * ---------------------------------------------------------------------- */

static inline int ll_floor_log2(uint64_t q) {          /* q >= 1 */
#if defined(__GNUC__)
    return 63 - __builtin_clzll((unsigned long long)q);
#else
    int k = 0;
    while (q >>= 1) k++;
    return k;
#endif
}

static inline int     chunk_of(int64_t j)   { return ll_floor_log2((uint64_t)(j / LL_B) + 1); }
static inline int64_t chunk_base(int k)     { return (int64_t)LL_B * ((INT64_C(1) << k) - 1); }
static inline int64_t chunk_cap(int k)      { return (int64_t)LL_B << k; }

static inline size_t kind_size(LLKind k) {
    switch (k) {
        case LL_I64:   return sizeof(int64_t);
        case LL_REAL:  return sizeof(double);
        case LL_CPLX:  return 2 * sizeof(double);
        case LL_BOXED: return sizeof(Expr*);
        default:       return 0;
    }
}

static inline int64_t run_len(const LLRun* r) { return r->end - r->start; }

/* Address of run slot j; the chunk must exist. */
static inline unsigned char* slot_ptr(const LLRun* r, int64_t j, size_t esz) {
    int k = chunk_of(j);
    return r->dir[k] + (size_t)(j - chunk_base(k)) * esz;
}

/* Make sure the chunk holding slot j exists. */
static void run_ensure(LinkedList* ll, LLRun* r, int64_t j) {
    int k = chunk_of(j);
    if (r->dir[k]) return;
    size_t bytes = (size_t)chunk_cap(k) * kind_size(ll->kind);
    r->dir[k] = malloc(bytes);
    if (!r->dir[k]) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    ll->chunk_bytes += (int64_t)bytes;
    MTH_MEM_CHARGE_BLOCK(bytes);
}

static void run_free_chunk(LinkedList* ll, LLRun* r, int k) {
    if (!r->dir[k]) return;
    size_t bytes = (size_t)chunk_cap(k) * kind_size(ll->kind);
    free(r->dir[k]);
    r->dir[k] = NULL;
    ll->chunk_bytes -= (int64_t)bytes;
    MTH_MEM_RELEASE(bytes);
}

/* Free chunks wholly beyond the run's end, keeping one spare chunk after the
 * one the next write lands in (hysteresis against push/pop thrash). */
static void run_trim_tail(LinkedList* ll, LLRun* r) {
    int keep = chunk_of(r->end) + 1;
    for (int k = keep + 1; k < LL_MAXCHUNK; k++) run_free_chunk(ll, r, k);
}

/* Move the live slots of a run down to [0, n) once the dead prefix is at
 * least as large as the live part (see the file comment). */
static void run_maybe_rebase(LinkedList* ll, LLRun* r) {
    int64_t n = run_len(r);
    if (r->start < n + LL_B) return;
    size_t esz = kind_size(ll->kind);
    if (n == 0) {
        r->start = r->end = 0;
    } else {
        /* Ascending copy: the destination is always below the source. */
        for (int64_t j = 0; j < n; j++) {
            run_ensure(ll, r, j);
            memcpy(slot_ptr(r, j, esz), slot_ptr(r, r->start + j, esz), esz);
        }
        r->start = 0;
        r->end = n;
    }
    run_trim_tail(ll, r);
}

/* ------------------------------------------------------------------------
 *  Mutation bookkeeping
 * ---------------------------------------------------------------------- */

/* Every change to a list's contents advances the evaluation clock, so no
 * fixed point stamped before the change is reused after it. (Nodes holding a
 * LinkedList are never INERT or GROUND -- an unknown atom is not a leaf to
 * those tests -- so the clock is the only cache that could otherwise see a
 * stale list.) One increment; no epoch is touched. */
static inline void ll_mutated(void) { eval_clock_tick(); }

/* ------------------------------------------------------------------------
 *  Lifecycle
 * ---------------------------------------------------------------------- */

LinkedList* ll_new(void) {
    LinkedList* ll = calloc(1, sizeof *ll);
    if (!ll) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    ll->refcount = 1;
    ll->id = g_ll_next_id++;
    ll->kind = LL_EMPTY;
    MTH_MEM_CHARGE(sizeof *ll);
    return ll;
}

LinkedList* ll_retain(LinkedList* ll) {
    if (ll) ll->refcount++;
    return ll;
}

/* Release the boxed elements of a run (LL_BOXED only). */
static void run_release_boxed(LLRun* r) {
    for (int64_t j = r->start; j < r->end; j++)
        expr_free(*(Expr**)slot_ptr(r, j, sizeof(Expr*)));
}

static void ll_free_storage(LinkedList* ll) {
    if (ll->kind == LL_BOXED) {
        run_release_boxed(&ll->front);
        run_release_boxed(&ll->back);
    }
    for (int k = 0; k < LL_MAXCHUNK; k++) {
        run_free_chunk(ll, &ll->front, k);
        run_free_chunk(ll, &ll->back, k);
    }
    ll->front.start = ll->front.end = 0;
    ll->back.start = ll->back.end = 0;
}

/* Release without recursion through nesting. Freeing a boxed list frees its
 * elements, and an element may be the last reference to another LinkedList:
 * recursing there costs C stack per level of nesting, and a 10^5-deep chain
 * of lists overflowed it. Instead, a payload that dies while another is being
 * torn down is queued, and the outermost release drains the queue -- constant
 * stack for any depth. (A plain List's own nesting is expr_free's business.) */
static LinkedList** g_ll_dying;
static size_t       g_ll_dying_n, g_ll_dying_cap;
static int          g_ll_releasing;

void ll_release(LinkedList* ll) {
    if (!ll || --ll->refcount > 0) return;
    if (g_ll_dying_n == g_ll_dying_cap) {
        size_t nc = g_ll_dying_cap ? 2 * g_ll_dying_cap : 16;
        LinkedList** nv = realloc(g_ll_dying, nc * sizeof *nv);
        if (!nv) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
        g_ll_dying = nv;
        g_ll_dying_cap = nc;
    }
    g_ll_dying[g_ll_dying_n++] = ll;
    if (g_ll_releasing) return;          /* the outer drain will free it */
    g_ll_releasing = 1;
    while (g_ll_dying_n > 0) {
        LinkedList* x = g_ll_dying[--g_ll_dying_n];
        ll_free_storage(x);              /* may queue more */
        MTH_MEM_RELEASE(sizeof *x);
        free(x);
    }
    g_ll_releasing = 0;
}

uint64_t ll_id(const LinkedList* ll)          { return ll ? ll->id : 0; }
int64_t  ll_length(const LinkedList* ll)      { return run_len(&ll->front) + run_len(&ll->back); }
LLKind   ll_kind(const LinkedList* ll)        { return ll->kind; }
int64_t  ll_payload_bytes(const LinkedList* ll) {
    return (int64_t)sizeof *ll + ll->chunk_bytes;
}

LinkedList* expr_linkedlist(const Expr* e) {
    return (e && e->type == EXPR_LINKEDLIST) ? e->data.llist : NULL;
}

/* ------------------------------------------------------------------------
 *  Element classification and typed <-> boxed
 * ---------------------------------------------------------------------- */

LLKind ll_classify(const Expr* e) {
    if (!e) return LL_BOXED;
    if (e->type == EXPR_INTEGER) return LL_I64;
    if (e->type == EXPR_REAL)    return LL_REAL;
    /* Machine complex: Complex[Real, Real], the same rule as packing (an exact
     * Gaussian integer stays boxed: widening it would change its head). */
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        head_is(e, SYM_Complex) &&
        e->data.function.args[0]->type == EXPR_REAL &&
        e->data.function.args[1]->type == EXPR_REAL)
        return LL_CPLX;
    return LL_BOXED;
}

/* Write a typed element held in `e` (already classified as ll->kind) into a
 * slot; consumes e. */
static void slot_write_typed(unsigned char* p, LLKind k, Expr* e) {
    switch (k) {
        case LL_I64:  { int64_t v = e->data.integer; memcpy(p, &v, sizeof v); break; }
        case LL_REAL: { double v = e->data.real;    memcpy(p, &v, sizeof v); break; }
        case LL_CPLX: {
            double v[2] = { e->data.function.args[0]->data.real,
                            e->data.function.args[1]->data.real };
            memcpy(p, v, sizeof v);
            break;
        }
        default: break;
    }
    expr_free(e);
}

/* Materialise the element in a slot as a new Expr. */
static Expr* slot_read(const unsigned char* p, LLKind k) {
    switch (k) {
        case LL_I64:  return ndarray_buffer_element_to_expr(p, 0, NDT_INT64);
        case LL_REAL: return ndarray_buffer_element_to_expr(p, 0, NDT_FLOAT64);
        case LL_CPLX: return ndarray_buffer_element_to_expr(p, 0, NDT_COMPLEX64);
        case LL_BOXED: return expr_copy(*(Expr* const*)p);
        default: return NULL;
    }
}

/* Re-store every live slot of `ll` at kind `to` (LL_BOXED from a typed kind),
 * keeping the slot numbering so the runs' start/end stay valid. */
static void ll_convert(LinkedList* ll, LLKind to) {
    LLKind from = ll->kind;
    size_t fsz = kind_size(from), tsz = kind_size(to);
    LLRun* runs[2] = { &ll->front, &ll->back };
    for (int ri = 0; ri < 2; ri++) {
        LLRun* r = runs[ri];
        LLRun nr;
        memset(&nr, 0, sizeof nr);
        nr.start = r->start;
        nr.end = r->end;
        /* Allocate the new chunks under the new kind's element size. */
        ll->kind = to;
        for (int64_t j = r->start; j < r->end; j++) run_ensure(ll, &nr, j);
        ll->kind = from;
        for (int64_t j = r->start; j < r->end; j++) {
            Expr* e = slot_read(slot_ptr(r, j, fsz), from);
            *(Expr**)slot_ptr(&nr, j, tsz) = e;
        }
        for (int k = 0; k < LL_MAXCHUNK; k++) run_free_chunk(ll, r, k);
        *r = nr;
    }
    ll->kind = to;
}

/* Settle the list's kind before storing an element of class `c`. After this
 * either ll->kind == c (typed) or ll->kind == LL_BOXED. */
static void ll_admit(LinkedList* ll, LLKind c) {
    if (ll->kind == c) return;
    if (ll_length(ll) == 0) {
        /* Empty: retype for free (drop chunks sized for the old kind). */
        ll_free_storage(ll);
        ll->kind = c;
        return;
    }
    if (ll->kind != LL_BOXED) ll_convert(ll, LL_BOXED);
}

/* Store `e` (consumed) into an existing slot pointer of the current kind. */
static void slot_store(LinkedList* ll, unsigned char* p, Expr* e) {
    if (ll->kind == LL_BOXED) *(Expr**)p = e;
    else slot_write_typed(p, ll->kind, e);
}

/* ------------------------------------------------------------------------
 *  Logical index -> (run, slot)
 * ---------------------------------------------------------------------- */

static inline unsigned char* ll_slot(const LinkedList* ll, int64_t i) {
    int64_t F = run_len(&ll->front);
    size_t esz = kind_size(ll->kind);
    if (i < F) return slot_ptr(&ll->front, ll->front.end - 1 - i, esz);
    return slot_ptr(&ll->back, ll->back.start + (i - F), esz);
}

int64_t ll_span_at(const LinkedList* ll, int64_t i, const void** p, int* dir) {
    int64_t F = run_len(&ll->front);
    size_t esz = kind_size(ll->kind);
    if (i < F) {                    /* front run: logical order descends slots */
        int64_t j = ll->front.end - 1 - i;
        int k = chunk_of(j);
        int64_t lo = chunk_base(k) > ll->front.start ? chunk_base(k) : ll->front.start;
        *p = slot_ptr(&ll->front, j, esz);
        *dir = -1;
        return j - lo + 1;
    }
    int64_t j = ll->back.start + (i - F);
    int k = chunk_of(j);
    int64_t hi = chunk_base(k) + chunk_cap(k);
    if (hi > ll->back.end) hi = ll->back.end;
    *p = slot_ptr(&ll->back, j, esz);
    *dir = 1;
    return hi - j;
}

Expr* ll_get(const LinkedList* ll, int64_t i) {
    return slot_read(ll_slot(ll, i), ll->kind);
}

const Expr* ll_peek_boxed(const LinkedList* ll, int64_t i) {
    if (ll->kind != LL_BOXED) return NULL;
    return *(Expr* const*)ll_slot(ll, i);
}

bool ll_get_i64(const LinkedList* ll, int64_t i, int64_t* out) {
    if (ll->kind != LL_I64) return false;
    memcpy(out, ll_slot(ll, i), sizeof *out);
    return true;
}

bool ll_get_real(const LinkedList* ll, int64_t i, double* out) {
    if (ll->kind != LL_REAL) return false;
    memcpy(out, ll_slot(ll, i), sizeof *out);
    return true;
}

bool ll_get_cplx(const LinkedList* ll, int64_t i, double* re, double* im) {
    if (ll->kind != LL_CPLX) return false;
    double v[2];
    memcpy(v, ll_slot(ll, i), sizeof v);
    *re = v[0];
    *im = v[1];
    return true;
}

void ll_set(LinkedList* ll, int64_t i, Expr* e) {
    ll_admit(ll, ll_classify(e));
    unsigned char* p = ll_slot(ll, i);
    if (ll->kind == LL_BOXED) {
        Expr* old = *(Expr**)p;
        *(Expr**)p = e;
        expr_free(old);
    } else {
        slot_write_typed(p, ll->kind, e);
    }
    ll_mutated();
}

bool ll_set_i64(LinkedList* ll, int64_t i, int64_t v) {
    if (ll->kind != LL_I64) return false;
    memcpy(ll_slot(ll, i), &v, sizeof v);
    ll_mutated();
    return true;
}

bool ll_set_real(LinkedList* ll, int64_t i, double v) {
    if (ll->kind != LL_REAL) return false;
    memcpy(ll_slot(ll, i), &v, sizeof v);
    ll_mutated();
    return true;
}

/* ------------------------------------------------------------------------
 *  Growth
 * ---------------------------------------------------------------------- */

/* Store e at the outer end of run r. The element is written BEFORE `end` is
 * advanced, so an abort can never expose an unwritten slot. */
static void run_push(LinkedList* ll, LLRun* r, Expr* e) {
    run_ensure(ll, r, r->end);
    slot_store(ll, slot_ptr(r, r->end, kind_size(ll->kind)), e);
    r->end++;
}

void ll_push_back(LinkedList* ll, Expr* e) {
    ll_admit(ll, ll_classify(e));
    run_push(ll, &ll->back, e);
    ll_mutated();
}

void ll_push_front(LinkedList* ll, Expr* e) {
    ll_admit(ll, ll_classify(e));
    run_push(ll, &ll->front, e);
    ll_mutated();
}

/* Typed pushes: no Expr is built when the list already has that kind. */
static void run_push_raw(LinkedList* ll, LLRun* r, const void* v) {
    run_ensure(ll, r, r->end);
    memcpy(slot_ptr(r, r->end, kind_size(ll->kind)), v, kind_size(ll->kind));
    r->end++;
}

#define LL_TYPED_PUSH(NAME, RUN, KIND, CTYPE, MAKE)                          \
    void NAME(LinkedList* ll, CTYPE v) {                                     \
        if (ll->kind == KIND || ll_length(ll) == 0) { \
            ll_admit(ll, KIND);                                              \
            run_push_raw(ll, &ll->RUN, &v);                                  \
            ll_mutated();                                                    \
        } else {                                                             \
            ll_push_##RUN(ll, MAKE);                                         \
        }                                                                    \
    }

LL_TYPED_PUSH(ll_push_back_i64,   back,  LL_I64,  int64_t, expr_new_integer(v))
LL_TYPED_PUSH(ll_push_back_real,  back,  LL_REAL, double,  expr_new_real(v))
LL_TYPED_PUSH(ll_push_front_i64,  front, LL_I64,  int64_t, expr_new_integer(v))
LL_TYPED_PUSH(ll_push_front_real, front, LL_REAL, double,  expr_new_real(v))

void ll_push_back_from(LinkedList* dst, const LinkedList* src, int64_t i) {
    if (src->kind != LL_BOXED && (dst->kind == src->kind || ll_length(dst) == 0)) {
        ll_admit(dst, src->kind);
        run_push_raw(dst, &dst->back, ll_slot(src, i));
        ll_mutated();
    } else {
        ll_push_back(dst, ll_get(src, i));
    }
}

static Expr* make_machine_complex(double re, double im) {
    Expr* a[2] = { expr_new_real(re), expr_new_real(im) };
    return expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
}

void ll_push_back_cplx(LinkedList* ll, double re, double im) {
    if (ll->kind == LL_CPLX || ll_length(ll) == 0) {
        double v[2] = { re, im };
        ll_admit(ll, LL_CPLX);
        run_push_raw(ll, &ll->back, v);
        ll_mutated();
    } else {
        ll_push_back(ll, make_machine_complex(re, im));
    }
}

void ll_push_front_cplx(LinkedList* ll, double re, double im) {
    if (ll->kind == LL_CPLX || ll_length(ll) == 0) {
        double v[2] = { re, im };
        ll_admit(ll, LL_CPLX);
        run_push_raw(ll, &ll->front, v);
        ll_mutated();
    } else {
        ll_push_front(ll, make_machine_complex(re, im));
    }
}

/* ------------------------------------------------------------------------
 *  Removal
 * ---------------------------------------------------------------------- */

/* Take the element out of slot p: a boxed element changes owner (no copy, no
 * free), a typed one is materialised. */
static Expr* slot_take(const LinkedList* ll, unsigned char* p) {
    if (ll->kind == LL_BOXED) return *(Expr**)p;
    return slot_read(p, ll->kind);
}

Expr* ll_pop_back(LinkedList* ll) {
    size_t esz = kind_size(ll->kind);
    Expr* e;
    if (run_len(&ll->back) > 0) {
        ll->back.end--;
        e = slot_take(ll, slot_ptr(&ll->back, ll->back.end, esz));
        run_trim_tail(ll, &ll->back);
    } else if (run_len(&ll->front) > 0) {
        /* The front run's lowest live slot is the list's last element. */
        e = slot_take(ll, slot_ptr(&ll->front, ll->front.start, esz));
        ll->front.start++;
        run_maybe_rebase(ll, &ll->front);
    } else {
        return NULL;
    }
    ll_mutated();
    return e;
}

Expr* ll_pop_front(LinkedList* ll) {
    size_t esz = kind_size(ll->kind);
    Expr* e;
    if (run_len(&ll->front) > 0) {
        ll->front.end--;
        e = slot_take(ll, slot_ptr(&ll->front, ll->front.end, esz));
        run_trim_tail(ll, &ll->front);
    } else if (run_len(&ll->back) > 0) {
        e = slot_take(ll, slot_ptr(&ll->back, ll->back.start, esz));
        ll->back.start++;
        run_maybe_rebase(ll, &ll->back);
    } else {
        return NULL;
    }
    ll_mutated();
    return e;
}

/* ------------------------------------------------------------------------
 *  Bulk
 * ---------------------------------------------------------------------- */

/* Append n raw elements of the list's (typed) kind, chunk-wise memcpy. */
static void ll_push_back_buffer(LinkedList* ll, const unsigned char* src, int64_t n) {
    size_t esz = kind_size(ll->kind);
    LLRun* r = &ll->back;
    while (n > 0) {
        run_ensure(ll, r, r->end);
        int k = chunk_of(r->end);
        int64_t room = chunk_base(k) + chunk_cap(k) - r->end;
        int64_t m = n < room ? n : room;
        memcpy(slot_ptr(r, r->end, esz), src, (size_t)m * esz);
        r->end += m;
        src += (size_t)m * esz;
        n -= m;
    }
}

LinkedList* ll_copy(const LinkedList* ll) {
    LinkedList* c = ll_new();
    int64_t n = ll_length(ll);
    c->kind = ll->kind;            /* same shape, even when empty */
    if (n == 0) return c;
    size_t esz = kind_size(ll->kind);
    /* Front run, reversed into logical order. */
    for (int64_t j = ll->front.end - 1; j >= ll->front.start; j--) {
        const unsigned char* p = slot_ptr(&ll->front, j, esz);
        if (ll->kind == LL_BOXED) {
            Expr* e = expr_copy(*(Expr* const*)p);
            run_push_raw(c, &c->back, &e);
        } else {
            run_push_raw(c, &c->back, p);
        }
    }
    /* Back run: contiguous within each chunk. */
    for (int64_t j = ll->back.start; j < ll->back.end; ) {
        int k = chunk_of(j);
        int64_t stop = chunk_base(k) + chunk_cap(k);
        if (stop > ll->back.end) stop = ll->back.end;
        const unsigned char* p = slot_ptr(&ll->back, j, esz);
        if (ll->kind == LL_BOXED) {
            for (int64_t t = j; t < stop; t++) {
                Expr* e = expr_copy(*(Expr* const*)slot_ptr(&ll->back, t, esz));
                run_push_raw(c, &c->back, &e);
            }
        } else {
            ll_push_back_buffer(c, p, stop - j);
        }
        j = stop;
    }
    return c;
}

bool ll_is_list_source(const Expr* src) {
    if (!src) return false;
    if (src->type == EXPR_NDARRAY) return true;          /* packed or visible */
    return src->type == EXPR_FUNCTION && head_is(src, SYM_List);
}

bool ll_append_elements(LinkedList* ll, const Expr* src) {
    if (!ll_is_list_source(src)) return false;
    if (src->type == EXPR_NDARRAY) {
        const NDArrayData* a = &src->data.ndarray;
        int64_t n = a->rank >= 1 ? a->dims[0] : 0;
        LLKind k = a->dtype == NDT_INT64     ? LL_I64
                 : a->dtype == NDT_FLOAT64   ? LL_REAL
                 : a->dtype == NDT_COMPLEX64 ? LL_CPLX : LL_EMPTY;
        if (a->rank == 1 && k != LL_EMPTY) {
            if (n == 0) return true;
            ll_admit(ll, k);
            if (ll->kind == k) {
                ll_push_back_buffer(ll, (const unsigned char*)a->data, n);
                ll_mutated();
                return true;
            }
            /* Boxed list: fall through to element-wise. */
        }
        /* Any other dtype or rank: element-wise through the List view. Rows of
         * a rank >= 2 array are offered back to the packer, so a list of
         * packed rows stays packed. */
        Expr* nested = ndarray_to_nested_list(src);
        size_t m = nested->data.function.arg_count;
        for (size_t i = 0; i < m; i++) {
            Expr* el = expr_copy(nested->data.function.args[i]);
            if (a->rank >= 2) el = pack_offer(el);
            ll_push_back(ll, el);
        }
        expr_free(nested);
        return true;
    }
    size_t m = src->data.function.arg_count;
    for (size_t i = 0; i < m; i++)
        ll_push_back(ll, expr_copy(src->data.function.args[i]));
    return true;
}

Expr* ll_to_list(const LinkedList* ll) {
    int64_t n = ll_length(ll);
    if (n > 0 && ll->kind != LL_BOXED) {
        NDType dt = ll->kind == LL_I64 ? NDT_INT64
                  : ll->kind == LL_REAL ? NDT_FLOAT64 : NDT_COMPLEX64;
        void* buf = NULL;
        Expr* out = ndbuild_open(1, &n, dt, &buf);
        if (out) {
            size_t esz = kind_size(ll->kind);
            unsigned char* dst = buf;
            for (int64_t j = ll->front.end - 1; j >= ll->front.start; j--, dst += esz)
                memcpy(dst, slot_ptr(&ll->front, j, esz), esz);
            for (int64_t j = ll->back.start; j < ll->back.end; ) {
                int k = chunk_of(j);
                int64_t stop = chunk_base(k) + chunk_cap(k);
                if (stop > ll->back.end) stop = ll->back.end;
                size_t bytes = (size_t)(stop - j) * esz;
                memcpy(dst, slot_ptr(&ll->back, j, esz), bytes);
                dst += bytes;
                j = stop;
            }
            return out;
        }
        /* Too small to pack (or packing off): a plain List below. */
    }
    Expr** args = n ? malloc((size_t)n * sizeof(Expr*)) : NULL;
    if (n && !args) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
    for (int64_t i = 0; i < n; i++) args[i] = ll_get(ll, i);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)n);
    free(args);
    return out;
}

/* ------------------------------------------------------------------------
 *  The LinkedList[] constructor
 * ---------------------------------------------------------------------- */

/* LinkedList[]          -> a new empty list
 * LinkedList[list]      -> a new list holding the elements of list
 * LinkedList[ll]        -> an independent shallow copy of ll
 * anything else         -> LinkedList::list, unevaluated */
static Expr* builtin_linkedlist(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc > 1) return builtin_arg_error("LinkedList", argc, 0, 1);
    if (argc == 0) return expr_new_linkedlist(ll_new());
    const Expr* a = res->data.function.args[0];
    LinkedList* src = expr_linkedlist(a);
    if (src) return expr_new_linkedlist(ll_copy(src));
    if (!ll_is_list_source(a)) {
        char* s = expr_to_string(res);
        mth_message("LinkedList", "list", "List expected at position 1 in %s.", s ? s : "");
        free(s);
        return NULL;
    }
    LinkedList* ll = ll_new();
    ll_append_elements(ll, a);
    return expr_new_linkedlist(ll);
}

/* ------------------------------------------------------------------------
 *  Cycle check (mutators only)
 * ---------------------------------------------------------------------- */

/* Does storing `e` into `target` create a reference cycle -- is `target`
 * reachable from a LinkedList inside `e`? Each list is visited once per check
 * (its `visit` stamp is set to the check's epoch, O(1) per test), so the walk
 * is O(|e| + the boxed lists reachable from it). Typed lists hold only
 * numbers and are not entered. */
static uint64_t g_ll_visit_epoch = 0;

/* Iterative (an explicit worklist): a deep chain of nested lists must not
 * cost C stack -- a 10^6-deep chain overflowed the recursive form. */
typedef struct { const Expr** v; size_t n, cap; } LLWork;

static void work_push(LLWork* w, const Expr* e) {
    if (!e) return;
    if (w->n == w->cap) {
        size_t nc = w->cap ? 2 * w->cap : 64;
        const Expr** nv = realloc(w->v, nc * sizeof *nv);
        if (!nv) { fprintf(stderr, "LinkedList: out of memory\n"); abort(); }
        w->v = nv;
        w->cap = nc;
    }
    w->v[w->n++] = e;
}

static bool reaches(const Expr* root, const LinkedList* target, uint64_t epoch) {
    LLWork w = { NULL, 0, 0 };
    bool found = false;
    work_push(&w, root);
    while (w.n > 0 && !found) {
        const Expr* e = w.v[--w.n];
        if (e->type == EXPR_LINKEDLIST) {
            LinkedList* l = e->data.llist;
            if (l == target) { found = true; break; }
            if (l->kind != LL_BOXED || l->visit == epoch) continue;
            l->visit = epoch;
            int64_t n = ll_length(l);
            for (int64_t i = 0; i < n; i++) work_push(&w, ll_peek_boxed(l, i));
        } else if (e->type == EXPR_FUNCTION) {
            work_push(&w, e->data.function.head);
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                work_push(&w, e->data.function.args[i]);
        }
    }
    free(w.v);
    return found;
}

bool ll_would_cycle(const LinkedList* target, const Expr* e) {
    if (ll_classify(e) != LL_BOXED) return false;      /* a machine number */
    return reaches(e, target, ++g_ll_visit_epoch);
}

/* Emit Head::llcyc for storing `e`, and return NULL (call left unevaluated). */
void ll_cycle_message(const char* head, Expr* e) {
    char* s = expr_to_string(e);
    mth_message(head, "llcyc", "Storing %s would make the linked list contain itself.",
                s ? s : "");
    free(s);
}

static Expr* ll_cycle_error(const char* head, Expr* e) {
    char* s = expr_to_string(e);
    mth_message(head, "llcyc", "Storing %s would make the linked list contain itself.",
                s ? s : "");
    free(s);
    return NULL;
}

/* ByteCount of a list and every LinkedList nested in it, iteratively: nested
 * lists go on a worklist instead of the C stack (a 10^6-deep chain overflowed
 * the recursive count). Every other boxed element is measured by `elem`. */
int64_t ll_byte_count(const LinkedList* root, int64_t (*elem)(Expr*)) {
    LLWork w = { NULL, 0, 0 };
    int64_t total = ll_payload_bytes(root);
    const LinkedList* l = root;
    for (;;) {
        if (l->kind == LL_BOXED) {
            int64_t n = ll_length(l);
            for (int64_t i = 0; i < n; i++) {
                const Expr* e = ll_peek_boxed(l, i);
                if (e->type == EXPR_LINKEDLIST) {
                    total += (int64_t)sizeof(Expr) + ll_payload_bytes(e->data.llist);
                    work_push(&w, e);
                } else {
                    total += elem((Expr*)e);
                }
            }
        }
        if (w.n == 0) break;
        l = w.v[--w.n]->data.llist;
    }
    free(w.v);
    return total;
}

/* ------------------------------------------------------------------------
 *  Push / Pop / PopFirst
 * ---------------------------------------------------------------------- */

/* The value a stack head acts on. `target` is the held first argument: a
 * symbol (which must have a value) or any expression; it is evaluated.
 * Returns the owned value, or NULL after Head::rvalue. */
static Expr* stack_target_value(const char* head, Expr* target) {
    if (target->type == EXPR_SYMBOL && !symtab_get_own_values(target->data.symbol.name)) {
        mth_message(head, "rvalue",
                    "%s is not a variable with a value, so its value cannot be changed.",
                    target->data.symbol.name);
        return NULL;
    }
    return evaluate(target);
}

/* Head::normal for a target whose value is not a list; returns NULL. */
static Expr* stack_not_list(const char* head, Expr* res) {
    char* s = expr_to_string(res);
    mth_message(head, "normal", "Nonatomic expression expected at position 1 in %s.",
                s ? s : "");
    free(s);
    return NULL;
}

/* A value a stack head can treat as a plain list: a List or packed List. */
static bool stack_plain_list(const Expr* v) {
    return (v->type == EXPR_FUNCTION && head_is(v, SYM_List)) || is_packed_list(v);
}

/* Evaluate h[arg] (arg borrowed). */
static Expr* eval1(const char* h, const Expr* arg) {
    Expr* a[1] = { expr_copy((Expr*)arg) };
    Expr* call = expr_new_function(expr_new_symbol(h), a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Push[ll, e]: append e to the linked list ll in place and return ll.
 * Push[s, e] for a symbol (or Part expression) holding a List is AppendTo. */
static Expr* builtin_push(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("Push", argc, 2, 2);
    Expr* target = res->data.function.args[0];
    Expr* e = res->data.function.args[1];
    Expr* val = stack_target_value("Push", target);
    if (!val) return NULL;
    LinkedList* ll = expr_linkedlist(val);
    if (ll) {
        if (ll_would_cycle(ll, e)) { expr_free(val); return ll_cycle_error("Push", e); }
        ll_push_back(ll, expr_copy(e));
        return val;
    }
    if (!stack_plain_list(val) || target->type == EXPR_LINKEDLIST) {
        expr_free(val);
        return stack_not_list("Push", res);
    }
    expr_free(val);
    /* A plain list: exactly AppendTo (O(n), rebinding the variable). */
    Expr* a[2] = { expr_copy(target), expr_copy(e) };
    Expr* call = expr_new_function(expr_new_symbol("AppendTo"), a, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Pop[ll] / PopFirst[ll]: remove the last / first element in place and
 * return it. For a symbol holding a List: take the element and rebind the
 * symbol to Most / Rest of it. */
static Expr* pop_impl(Expr* res, bool first) {
    const char* head = first ? "PopFirst" : "Pop";
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error(head, argc, 1, 1);
    Expr* target = res->data.function.args[0];
    Expr* val = stack_target_value(head, target);
    if (!val) return NULL;
    LinkedList* ll = expr_linkedlist(val);
    bool empty = ll ? ll_length(ll) == 0
                    : stack_plain_list(val) &&
                      (val->type == EXPR_FUNCTION ? val->data.function.arg_count == 0
                                                  : val->data.ndarray.dims[0] == 0);
    if (empty) {
        char* s = expr_to_string(target);
        mth_message(head, "empty", "Cannot pop from an empty list %s.", s ? s : "");
        free(s);
        expr_free(val);
        return NULL;
    }
    if (ll) {
        Expr* out = first ? ll_pop_front(ll) : ll_pop_back(ll);
        expr_free(val);
        return out;
    }
    if (!stack_plain_list(val) || target->type == EXPR_LINKEDLIST) {
        expr_free(val);
        return stack_not_list(head, res);
    }
    Expr* out = eval1(first ? "First" : "Last", val);
    Expr* rest = eval1(first ? "Rest" : "Most", val);
    expr_free(val);
    Expr* sa[2] = { expr_copy(target), rest };
    Expr* set = expr_new_function(expr_new_symbol(SYM_Set), sa, 2);
    expr_free(evaluate(set));
    expr_free(set);
    return out;
}

static Expr* builtin_pop(Expr* res)      { return pop_impl(res, false); }
static Expr* builtin_popfirst(Expr* res) { return pop_impl(res, true); }

void linkedlist_init(void) {
    symtab_add_builtin(SYM_LinkedList, builtin_linkedlist);
    symtab_get_def(SYM_LinkedList)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LinkedList",
        "LinkedList[]\n\tCreates an empty mutable linked list. LinkedList[list] creates one\n"
        "\tholding the elements of list (a List, packed array or NDArray);\n"
        "\tLinkedList[ll] returns an independent shallow copy of the linked list ll.\n"
        "\tAssignment and argument passing share the same list rather than copying\n"
        "\tit: AppendTo, PrependTo, Push, Pop, PopFirst and Part assignment change it\n"
        "\tin place, while Part, Take, Drop, Append, Delete, Map, ... return a new\n"
        "\tLinkedList. Normal[ll] gives the elements as a List. Prints as\n"
        "\tLinkedList[<n>]; appending, prepending, popping and indexing are O(1).");

    symtab_add_builtin(SYM_Push, builtin_push);
    symtab_get_def(SYM_Push)->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("Push",
        "Push[ll, e]\n\tAppends e to the end of the linked list ll in place, in O(1), and\n"
        "\treturns ll. Push[s, e] for a variable s holding a List appends e and\n"
        "\tresets s (as AppendTo, O(n)). With Pop, a last-in first-out stack; with\n"
        "\tPopFirst, a first-in first-out queue. Storing a linked list into itself,\n"
        "\tdirectly or through other lists, gives Push::llcyc.");

    symtab_add_builtin(SYM_Pop, builtin_pop);
    symtab_get_def(SYM_Pop)->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("Pop",
        "Pop[ll]\n\tRemoves the last element of the linked list ll in place, in O(1), and\n"
        "\treturns it. Pop[s] for a variable s holding a List removes its last\n"
        "\telement and resets s (O(n)). An empty list gives Pop::empty and the call\n"
        "\tstays unevaluated.");

    symtab_add_builtin(SYM_PopFirst, builtin_popfirst);
    symtab_get_def(SYM_PopFirst)->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("PopFirst",
        "PopFirst[ll]\n\tRemoves the first element of the linked list ll in place, in O(1),\n"
        "\tand returns it; Push with PopFirst gives a first-in first-out queue.\n"
        "\tPopFirst[s] for a variable s holding a List removes its first element\n"
        "\tand resets s (O(n)). An empty list gives PopFirst::empty and the call\n"
        "\tstays unevaluated.");
}
