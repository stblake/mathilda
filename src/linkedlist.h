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

/* A new list of `n` raw slots of typed `kind` (I64, REAL or CPLX) copied from
 * `buf` chunk-wise: where a kernel that computed into a scratch buffer hands
 * its answer over. */
LinkedList* ll_from_buffer(LLKind kind, const void* buf, int64_t n);

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
/* The storage cursor every LinkedList algorithm walks: elements [i, i+k) of
 * the logical order occupy k consecutive slots of one chunk. Returns k >= 1
 * (i must be in range); *p is element i's slot and *dir is +1 when the next
 * logical element is at the next higher address, -1 when lower (the front
 * run, which grows downward). Slots are kind-sized: int64_t, double, a
 * (re, im) double pair, or Expr*. A full scan is O(n) with O(log n) cursor
 * calls -- chunk capacities double. */
int64_t ll_span_at(const LinkedList* ll, int64_t i, const void** p, int* dir);
/* Span cursor over elements [lo, hi): the loop every native LinkedList
 * algorithm runs, so none hand-rolls ll_span_at bookkeeping.
 *
 *     LLCursor c; const void* p; int dir; int64_t base, k;
 *     ll_cursor(&c, ll, 0, ll_length(ll));
 *     while ((k = ll_cursor_next(&c, &p, &dir, &base)) > 0)
 *         for (int64_t t = 0; t < k; t++)
 *             use(LL_SLOT(double, p, dir, t));       (element base + t)
 *
 * A complex slot is a (re, im) pair: LL_CSLOT(p, dir, t) points at re. */
typedef struct { const LinkedList* ll; int64_t i, end; } LLCursor;
static inline void ll_cursor(LLCursor* c, const LinkedList* ll, int64_t lo, int64_t hi) {
    c->ll = ll; c->i = lo; c->end = hi;
}
static inline int64_t ll_cursor_next(LLCursor* c, const void** p, int* dir, int64_t* base) {
    if (c->i >= c->end) return 0;
    int64_t k = ll_span_at(c->ll, c->i, p, dir);
    if (k > c->end - c->i) k = c->end - c->i;
    *base = c->i;
    c->i += k;
    return k;
}
#define LL_SLOT(T, p, dir, t) (((const T*)(p))[(dir) > 0 ? (t) : -(t)])
#define LL_CSLOT(p, dir, t)   ((const double*)(p) + 2 * ((dir) > 0 ? (t) : -(t)))

/* Bytes per slot of the list's storage kind (0 for an empty list). */
size_t ll_slot_size(const LinkedList* ll);
/* Elements [lo, lo + n) copied into dst in logical order as raw slots
 * (int64_t, double, (re, im) pairs, or BORROWED Expr* for a boxed list):
 * the scratch buffer of an algorithm that needs random access across the
 * whole sequence (Sort, Median, Fourier). Not a List: nothing is boxed and
 * no Expr is built. */
void ll_gather(const LinkedList* ll, int64_t lo, int64_t n, void* dst);

/* Borrowed pointer to a boxed element (LL_BOXED only), else NULL. */
const Expr* ll_peek_boxed(const LinkedList* ll, int64_t i);
/* Typed reads for the compiled path; false when the list is not that kind. */
bool ll_get_i64 (const LinkedList* ll, int64_t i, int64_t* out);
bool ll_get_real(const LinkedList* ll, int64_t i, double* out);
bool ll_get_cplx(const LinkedList* ll, int64_t i, double* re, double* im);
void ll_set(LinkedList* ll, int64_t i, Expr* e);        /* consumes e      */
/* Raw typed writes: false (nothing written) unless the list is that kind. */
bool ll_set_i64 (LinkedList* ll, int64_t i, int64_t v);
bool ll_set_real(LinkedList* ll, int64_t i, double v);

/* Element i for code that only reads it (a matcher, a level walker): the
 * slot itself, borrowed, for a boxed list; a fresh box of a typed slot
 * otherwise. Release with ll_elem_done. No List is built either way. */
static inline Expr* ll_elem_view(const LinkedList* ll, int64_t i, bool* owned) {
    const Expr* b = ll_peek_boxed(ll, i);
    if (b) { *owned = false; return (Expr*)b; }
    *owned = true;
    return ll_get(ll, i);
}
static inline void ll_elem_done(Expr* e, bool owned) { if (owned) expr_free(e); }

/* ---- growth and removal ---------------------------------------------- */
void  ll_push_back (LinkedList* ll, Expr* e);           /* consumes e      */
void  ll_push_front(LinkedList* ll, Expr* e);           /* consumes e      */
void  ll_push_back_i64 (LinkedList* ll, int64_t v);
void  ll_push_back_real(LinkedList* ll, double v);
void  ll_push_back_cplx(LinkedList* ll, double re, double im);
void  ll_push_front_i64 (LinkedList* ll, int64_t v);
void  ll_push_front_real(LinkedList* ll, double v);
void  ll_push_front_cplx(LinkedList* ll, double re, double im);
/* Append element i of src to dst; typed slots copy without boxing. */
void  ll_push_back_from(LinkedList* dst, const LinkedList* src, int64_t i);
/* Remove and return the last / first element. NULL when empty. */
Expr* ll_pop_back (LinkedList* ll);
Expr* ll_pop_front(LinkedList* ll);
/* Typed pop for the compiled path: the removed slot's bytes (int64_t, double
 * or a (re, im) pair) into `out`. False when empty or boxed. */
bool  ll_pop_typed(LinkedList* ll, bool front, void* out);
bool  ll_set_cplx(LinkedList* ll, int64_t i, double re, double im);
/* The elements as a packed rank-1 array, any length; NULL if boxed or empty. */
Expr* ll_to_ndarray(const LinkedList* ll);

/* ---- undo frames (compiled runs; see linkedlist.c, "Undo log") ---------
 * While a frame is open every change to a list created BEFORE the frame was
 * opened is logged. ll_undo_rollback reverts the innermost frame's changes;
 * ll_undo_commit closes it (an inner commit keeps the records for the outer
 * frame). Frames nest; every begin is paired with one commit or rollback. */
void     ll_undo_begin(void);
void     ll_undo_commit(void);
void     ll_undo_rollback(void);
uint64_t ll_next_id(void);

/* Elements [lo, hi) of src appended to dst, in order / in reverse order.
 * Typed spans copy without boxing. src and dst must differ. */
void  ll_append_range    (LinkedList* dst, const LinkedList* src, int64_t lo, int64_t hi);
void  ll_append_range_rev(LinkedList* dst, const LinkedList* src, int64_t lo, int64_t hi);
/* Size an EMPTY list to n slots of typed `kind` whose contents are NOT yet
 * written: the caller fills every slot through ll_span_at (the back run, so
 * the cursor's direction is +1 and casting its pointer to writable is sound
 * for a list the caller alone holds) before the list escapes. For kernels
 * that write their output straight into the chunks. */
void  ll_reserve_typed(LinkedList* ll, LLKind kind, int64_t n);

/* Append every element of a List, packed List or (any-rank) NDArray `src`
 * (borrowed). Rank-1 int64/float64/complex64 buffers are copied chunk-wise.
 * Returns false (list unchanged) when `src` is none of those. */
bool  ll_append_elements(LinkedList* ll, const Expr* src);
/* True when `src` is something ll_append_elements accepts. */
bool  ll_is_list_source(const Expr* src);

/* The elements as a List: a packed List when the storage is typed and the
 * result is large enough to pack, otherwise a plain List. New reference. */
Expr* ll_to_list(const LinkedList* ll);
/* Whole-list materialisations since startup (ll_to_list, ll_block). Only
 * Normal may materialise a LinkedList; every other head reads the chunks. */
uint64_t ll_materialise_count(void);

/* Would storing `e` (borrowed) into `target` make `target` reachable from
 * itself -- a reference cycle no refcount could free? Mutators reject it. */
bool  ll_would_cycle(const LinkedList* target, const Expr* e);

/* Head::llcyc: storing `e` would make a list contain itself. */
void  ll_cycle_message(const char* head, Expr* e);

/* Bytes held by `ll` and, recursively, by its elements -- iterative through
 * nested LinkedLists; `elem` measures every other boxed element. */
int64_t ll_byte_count(const LinkedList* ll, int64_t (*elem)(Expr*));

/* A walker that may rebuild elements (ReplaceAll, Map at a level, N, ...):
 * fn(elem, i, ctx) returns a NEW element, or NULL to keep elem as it is
 * (elem is borrowed, see ll_elem_view). The answer is a new LinkedList when
 * some element changed and NULL when none did, so the caller can hand back
 * the list itself. Reads a snapshot: fn may run code that changes `src`.
 * *unwind receives a Throw/Abort sentinel fn returned (the walk stops; the
 * caller hands it on). */
Expr* ll_rebuild(const LinkedList* src, Expr* (*fn)(Expr* elem, int64_t i, void* ctx),
                 void* ctx, Expr** unwind);

/* The sequences the threading heads (MapThread, Inner, Outer, ...) index into
 * alike: a LinkedList, a List, or an array along its leading axis.
 * ll_seq_length is -1 for anything else; ll_seq_element(a, i) is element i
 * (0-based, in range) as a new reference -- read from the slots of a
 * LinkedList, a row (sub-array) of a rank >= 2 array. */
int64_t ll_seq_length(const Expr* a);
Expr*   ll_seq_element(const Expr* a, int64_t i);

/* head[e1, ..., en] (unevaluated) from the elements of `ll`, read from the
 * slots; consumes `head`. The call Apply, Thread and their kin build from a
 * List's arguments -- here from the chunks, with no List in between. */
Expr* ll_call_with_elements(Expr* head, const LinkedList* ll);

/* ---- Expr-level helpers ----------------------------------------------- */
/* The payload of an EXPR_LINKEDLIST node, or NULL for anything else. */
LinkedList* expr_linkedlist(const Expr* e);

/* Total of the elements of a typed list, in place (linkedlist_ops.c): exact
 * int64 promoting to GMP on overflow, or the packed Total's pairwise float
 * tree to the bit. NULL for a boxed or empty list, or a float sum that left
 * the machine range (the caller's Plus then answers). */
Expr* ll_sum_machine(const LinkedList* ll);

/* ---- copying heads (linkedlist_struct.c, funcprog.c, replace.c) ------
 * Each takes the evaluated call `res` (borrowed) whose list argument is an
 * EXPR_LINKEDLIST, and returns a NEW LinkedList (the argument is untouched),
 * or NULL after a message (the call stays unevaluated). The head's own
 * builtin dispatches here before its List code. */
Expr* ll_take_drop(Expr* res, bool take);      /* Take / Drop[ll, seq, ...]    */
Expr* ll_most_rest(Expr* res, bool most);      /* Most / Rest[ll]              */
Expr* ll_reverse(Expr* res);                   /* Reverse[ll], Reverse[ll, lv] */
Expr* ll_append_prepend(Expr* res, bool prepend);
Expr* ll_join(Expr* res);                      /* any argument a LinkedList    */
Expr* ll_insert(Expr* res);                    /* Insert[ll, e, pos]           */
Expr* ll_delete(Expr* res);                    /* Delete[ll, pos]              */
Expr* ll_select(Expr* res);                    /* Select[ll, crit (, n)]       */
Expr* ll_map(Expr* res);                       /* funcprog.c                   */
Expr* ll_replace_part(Expr* res);              /* replace.c                    */
/* A Listable call with a LinkedList argument, threaded as over a List; the
 * answer is a LinkedList. NULL when no argument is a LinkedList; the call
 * itself back (a copy) when it settles after Thread::tdlen. */
Expr* ll_thread_listable(Expr* res);
/* The direct arm of ll_thread_listable (linkedlist_ops.c): a machine-real
 * list through a head's unary ND kernel, or Plus/Times of machine-real lists
 * and machine scalars, computed from the chunks into the chunks of the new
 * list, threaded. Bit-identical to the packed List; NULL declines (any
 * element the packed path would not answer at machine precision), and the
 * caller then streams. */
Expr* ll_listable_direct(Expr* res);
/* The whole of Listable threading once the operand lengths are known to agree
 * (`len` elements): ll_listable_direct, else the machine engine a packed List
 * would reach (ndarray.h, NDOperand) over the slots, else head[a1_i, ...]
 * evaluated element by element from snapshots of the list operands. */
Expr* ll_listable_threaded(Expr* res, int64_t len);
/* Map[f, ll] at level 1 / Select[ll, crit] through numloop's compiled function
 * or predicate, straight from the typed slots; NULL declines (the caller maps
 * element by element). */
Expr* ll_map_compiled(const Expr* f, const LinkedList* src);
Expr* ll_select_compiled(const Expr* crit, const LinkedList* src);

/* Registers LinkedList, Push, Pop and PopFirst. */
void linkedlist_init(void);

/* ---- element views and per-element threading (numeric / symbolic heads) --
 * For a List core written over (Expr* const* elems, size_t n), so a List and
 * a LinkedList run one algorithm. A List lends its argument array (borrowed);
 * a LinkedList's elements are read from the slots into a malloc'd array of
 * owned references (a refcount bump for a boxed slot, a fresh box for a typed
 * one) -- an element array, not a List: no List Expr is built. */
typedef struct { Expr** v; int64_t n; bool owned; bool ll; } LLElems;
/* False (out cleared) when `e` is neither a List nor a LinkedList. */
bool  ll_elems_of(const Expr* e, LLElems* out);
void  ll_elems_release(LLElems* x);
/* Element i of a List (its argument, borrowed) or a LinkedList (ll_elem_view);
 * release with ll_elem_done. NULL for anything else. */
Expr* ll_seq_view(const Expr* e, int64_t i, bool* owned);
/* True for a List or a LinkedList. */
bool  ll_is_seq(const Expr* e);
/* H[a1, ..., ll, ..., an] for a head whose List form maps over the elements of
 * argument k: H[a1, ..., e_i, ..., an] evaluated for each element of a
 * snapshot of the LinkedList at position k, the answers in a new LinkedList
 * (the list itself, a new reference, when every answer is its element). A
 * Throw/Abort sentinel is handed back as it came. NULL when argument k is not
 * a LinkedList. */
Expr* ll_map_arg(const Expr* res, size_t k);
/* A List answer as the LinkedList of its elements (consumes `list`): the
 * result side of a head that computed through a shared element core. A rank-1
 * packed List is copied chunk-wise; anything else is returned unchanged. */
Expr* ll_from_list_result(Expr* list);
/* True when `e` holds a LinkedList anywhere (one walk, no evaluation). */
bool  ll_expr_has_linkedlist(const Expr* e);
/* A copying tree walker (Rationalize, Chop, N, ...) at a LinkedList node:
 * fn(el, ctx) is the walk of one element -- a new reference, the element
 * itself (expr_copy) when unchanged, or NULL to abandon the walk. The answer
 * is a new LinkedList of the walked elements, the list itself (a new
 * reference) when no element changed, or NULL when fn abandoned. */
Expr* ll_walk(const Expr* e, Expr* (*fn)(const Expr* el, void* ctx), void* ctx);
/* ---- rectangular nests of Lists and LinkedLists (matrix-shaped heads) ----
 * A level may be a List, a LinkedList, or (below the top) an NDArray block.
 * ll_tensor_dims: the rank of the rectangular nest `e` (0 for a leaf), dims
 * in `dims` (at most `maxrank` levels are descended); -1 when ragged.
 * ll_tensor_leaves: its leaves in row-major order as new references (the
 * caller sized `out` from the dims). ll_tensor_levels: lvl[k] true when level
 * k (0 = outermost) is a LinkedList along the first-element spine.
 * ll_tensor_build: the nest of shape dims[0..rank) over `leaves` (consumed,
 * row-major), each level a LinkedList where lvl[k] says so, else a List.
 * ll_has_ll_level1: `e` is a LinkedList, or a List with a LinkedList element:
 * the cheap gate a List head's LinkedList arm tests before anything else. */
int   ll_tensor_dims(const Expr* e, int64_t* dims, int maxrank);
void  ll_tensor_leaves(const Expr* e, int rank, Expr** out, size_t* idx);
void  ll_tensor_levels(const Expr* e, int rank, bool* lvl);
Expr* ll_tensor_build(const int64_t* dims, int rank, Expr** leaves, const bool* lvl);
bool  ll_has_ll_level1(const Expr* e);
/* Would Normal[ll] hand back a packed List (typed storage, packing on, at
 * least pack_min_elements)? Where the List heads answer a packed and a plain
 * List differently, a LinkedList follows its own Normal. */
bool  ll_normal_is_packed(const LinkedList* ll);

#endif /* MATHILDA_LINKEDLIST_H */
