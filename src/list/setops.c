#include "eval.h"
#include "list_common.h"
#include "setops.h"
#include "message.h"    /* mth_message: Quiet/Check funnel */
#include "assoc.h"
#include "ndarray.h"    /* is_ndarray */
#include "ndarray_internal.h"  /* nd_sort_i64_asc — the packed set-op fast path */
#include "ndreduce.h"   /* ndred_tally — Tally's packed-buffer fast path */
#include "linkedlist_seq.h" /* LLSeq: LinkedList operands, read element by element */

/* ---------------------------------------------------------------------------
 *  Packed integer set operations.
 *
 *  Union / Intersection / Complement are written generically above: expr_copy
 *  every element, qsort through a function pointer calling expr_compare, then
 *  dedup with expr_eq. On a packed integer list that is one Expr allocation and
 *  two indirect calls per element -- Union of 10^6 integers cost 745 ms, against
 *  ~60 ms for np.unique on the same data -- and it is the shape that matters,
 *  because a set operation over integer LABELS is what a graph traversal and a
 *  k-mer count are made of. Breadth-first search is literally
 *      Complement[Union[Flatten[adj[[frontier]]]], visited]
 *  once per level.
 *
 *  Restricted to int64 on purpose. Two Reals that compare equal but print
 *  differently (0. and -0.) would let the buffer path keep a different
 *  representative than the List path keeps, and a set operation whose ANSWER
 *  depends on which of two equal elements survived is not a fast path, it is a
 *  second implementation. Integers have no such pair.
 * ------------------------------------------------------------------------- */

/* True when `e` is a rank-1 int64 buffer -- the only shape below.
 *
 * `is_ndarray`, NOT `is_packed_list`: a VISIBLE NDArray[...] holds the same
 * int64 buffer and is the same fast domain. Presentation follows the input:
 * every packed result is built with expr_new_ndarray_like(args[0], ...). */
static bool setop_i64(const Expr* e) {
    return is_ndarray(e) && e->data.ndarray.rank == 1 &&
           e->data.ndarray.dtype == NDT_INT64;
}

/* ---- raw int64 cores ------------------------------------------------------
 * Every packed int64 set operation runs on raw (values, count) spans, so a
 * packed List and a typed LinkedList (its slots gathered into one scratch
 * buffer) run the same kernel and cannot disagree. Each core writes an owned
 * buffer to *out (count *m) and returns false only to decline (range too wide
 * for a direct index, allocation failure). */
typedef struct { const int64_t* v; size_t n; } I64Span;

/* Sorted, deduplicated copy of a span. */
static bool raw_sorted_unique(const int64_t* src, size_t n, int64_t** out, size_t* n_out) {
    int64_t* v = malloc(sizeof(int64_t) * (n ? n : 1));
    if (!v) return false;
    memcpy(v, src, sizeof(int64_t) * n);
    nd_sort_i64_asc(v, n);
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        if (m == 0 || v[i] != v[m - 1]) v[m++] = v[i];
    *out = v; *n_out = m;
    return true;
}

/* ---- Direct-index membership -------------------------------------------
 *
 * When the values come from a bounded range -- labels, ids, indices -- set
 * membership is a byte lookup, and walking the range in order emits the unique
 * values ALREADY SORTED. The table is only built when the range is no wider
 * than the input it replaces. Declines otherwise; the sort-merge path has no
 * range restriction. */
#define SETOP_DIRECT_MAX_RANGE ((uint64_t)1 << 26)

static uint64_t setop_direct_range(const int64_t* v, size_t n, int64_t* mn_out) {
    if (n == 0) return 0;
    int64_t mn = v[0], mx = v[0];
    for (size_t i = 1; i < n; i++) {
        if (v[i] < mn) mn = v[i];
        if (v[i] > mx) mx = v[i];
    }
    /* Unsigned: mx - mn as int64 overflows on a full-width spread. */
    uint64_t span = (uint64_t)mx - (uint64_t)mn;
    if (span >= SETOP_DIRECT_MAX_RANGE || span + 1 > (uint64_t)n) return 0;
    *mn_out = mn;
    return span + 1;
}

/* Deduplicate by direct indexing: `sorted` walks the range (ascending --
 * Union), else the input (first appearance -- DeleteDuplicates). */
static bool raw_direct_unique(const int64_t* v, size_t n, bool sorted, int64_t** outp, size_t* mp) {
    int64_t mn = 0;
    uint64_t range = setop_direct_range(v, n, &mn);
    if (!range) return false;
    unsigned char* seen = calloc((size_t)range, 1);
    if (!seen) return false;
    int64_t* out = malloc(sizeof(int64_t) * n);
    if (!out) { free(seen); return false; }
    size_t m = 0;
    if (sorted) {
        for (size_t i = 0; i < n; i++) seen[(uint64_t)v[i] - (uint64_t)mn] = 1;
        for (uint64_t d = 0; d < range; d++)
            if (seen[d]) out[m++] = mn + (int64_t)d;
    } else {
        for (size_t i = 0; i < n; i++) {
            uint64_t d = (uint64_t)v[i] - (uint64_t)mn;
            if (!seen[d]) { seen[d] = 1; out[m++] = v[i]; }
        }
    }
    free(seen);
    *outp = out; *mp = m;
    return true;
}

/* The many-list forms by the same direct index: one int32 per value holds a
 * BITMASK of the lists containing it, so Union (mask != 0), Intersection
 * (mask == all) and Complement (mask == 1) are one ascending walk. At most 32
 * lists. */
static bool raw_direct_multi(const I64Span* sp, size_t nl, int mode, int64_t** outp, size_t* mp) {
    if (nl < 2 || nl > 32) return false;
    int64_t mn = 0, mx = 0;
    size_t total = 0;
    bool seeded = false;
    for (size_t i = 0; i < nl; i++) {
        total += sp[i].n;
        for (size_t j = 0; j < sp[i].n; j++) {
            int64_t x = sp[i].v[j];
            if (!seeded) { mn = mx = x; seeded = true; continue; }
            if (x < mn) mn = x;
            if (x > mx) mx = x;
        }
    }
    if (!seeded || total == 0) { *outp = NULL; *mp = 0; return true; }   /* every list empty */
    uint64_t span = (uint64_t)mx - (uint64_t)mn;
    if (span >= SETOP_DIRECT_MAX_RANGE || span + 1 > (uint64_t)total) return false;
    uint64_t range = span + 1;
    uint32_t* mask = calloc((size_t)range, sizeof(uint32_t));
    if (!mask) return false;
    for (size_t i = 0; i < nl; i++) {
        uint32_t bit = (uint32_t)1 << i;
        for (size_t j = 0; j < sp[i].n; j++)
            mask[(uint64_t)sp[i].v[j] - (uint64_t)mn] |= bit;
    }
    uint32_t all = (nl == 32) ? 0xFFFFFFFFu : (((uint32_t)1 << nl) - 1);
    int64_t* out = malloc(sizeof(int64_t) * (size_t)range);
    if (!out) { free(mask); return false; }
    size_t m = 0;
    for (uint64_t d = 0; d < range; d++) {
        uint32_t k = mask[d];
        bool keep = (mode == 0) ? (k != 0) : (mode == 1) ? (k == all) : (k == 1u);
        if (keep) out[m++] = mn + (int64_t)d;
    }
    free(mask);
    *outp = out; *mp = m;
    return true;
}

/* Wide-range DeleteDuplicates: an open-addressed int64 set, the key in the
 * slot so a probe is one dependent load. */
typedef struct { int64_t key; uint8_t used; } SetSlot;

static bool raw_hash_first_unique(const int64_t* v, size_t n, int64_t** outp, size_t* mp) {
    size_t cap = 1024;
    while (cap * 7 < n * 10) {              /* keep load factor under 0.7 */
        if (cap > ((size_t)1 << 40)) return false;
        cap *= 2;
    }
    size_t mask = cap - 1;
    SetSlot* tab = calloc(cap, sizeof(SetSlot));
    if (!tab) return false;
    int64_t* out = malloc(sizeof(int64_t) * (n ? n : 1));
    if (!out) { free(tab); return false; }
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        /* splitmix64 finalizer: consecutive integers would otherwise collide */
        uint64_t z = (uint64_t)v[i];
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        size_t h = (size_t)(z ^ (z >> 31)) & mask;
        for (;;) {
            if (!tab[h].used) { tab[h].used = 1; tab[h].key = v[i]; out[m++] = v[i]; break; }
            if (tab[h].key == v[i]) break;
            h = (h + 1) & mask;
        }
    }
    free(tab);
    *outp = out; *mp = m;
    return true;
}

/* DeleteDuplicates of a span, first-appearance order. */
static bool raw_delete_duplicates(const int64_t* v, size_t n, int64_t** outp, size_t* mp) {
    if (n == 0) { *outp = NULL; *mp = 0; return true; }
    return raw_direct_unique(v, n, false, outp, mp) || raw_hash_first_unique(v, n, outp, mp);
}

/* Union / Intersection / Complement (mode 0, 1, 2) over `nl` spans. */
static bool raw_setop(const I64Span* sp, size_t nl, int mode, int64_t** outp, size_t* mp) {
    if (nl == 0) return false;
    if (nl == 1 && mode == 0 && sp[0].n > 0 && raw_direct_unique(sp[0].v, sp[0].n, true, outp, mp))
        return true;
    if (nl > 1 && raw_direct_multi(sp, nl, mode, outp, mp)) return true;

    int64_t* acc; size_t accn;
    if (!raw_sorted_unique(sp[0].v, sp[0].n, &acc, &accn)) return false;
    for (size_t i = 1; i < nl; i++) {
        int64_t* b; size_t bn;
        if (!raw_sorted_unique(sp[i].v, sp[i].n, &b, &bn)) { free(acc); return false; }
        size_t p = 0, q = 0, m = 0;
        if (mode == 0) {                              /* Union: merge */
            size_t cap = accn + bn;
            int64_t* merged = malloc(sizeof(int64_t) * (cap ? cap : 1));
            if (!merged) { free(acc); free(b); return false; }
            while (p < accn || q < bn) {
                int64_t v;
                if (q >= bn)                 v = acc[p++];
                else if (p >= accn)          v = b[q++];
                else if (acc[p] < b[q])      v = acc[p++];
                else if (b[q] < acc[p])      v = b[q++];
                else                       { v = acc[p++]; q++; }
                if (m == 0 || v != merged[m - 1]) merged[m++] = v;
            }
            free(acc); free(b); acc = merged; accn = m;
        } else if (mode == 1) {                       /* Intersection */
            while (p < accn && q < bn) {
                if (acc[p] < b[q]) p++;
                else if (b[q] < acc[p]) q++;
                else { acc[m++] = acc[p++]; q++; }
            }
            free(b); accn = m;
        } else {                                      /* Complement: acc \ b */
            while (p < accn) {
                while (q < bn && b[q] < acc[p]) q++;
                if (q < bn && b[q] == acc[p]) p++;
                else acc[m++] = acc[p++];
            }
            free(b); accn = m;
        }
    }
    *outp = acc; *mp = accn;
    return true;
}

/* ---- packed wrappers ---- */

/* An owned result buffer as a packed list shaped like `src` ({} when empty).
 * Frees `v`; NULL on allocation failure, so callers can degrade. */
static Expr* setop_emit(const Expr* src, int64_t* v, size_t m) {
    int64_t dims[1] = { (int64_t)m };
    if (m == 0) { free(v); return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }
    int64_t* buf = realloc(v, sizeof(int64_t) * m);   /* shrink to fit */
    if (!buf) buf = v;
    return expr_new_ndarray_like(src, 1, dims, buf, NDT_INT64);
}

static Expr* setop_packed_delete_duplicates(const Expr* a) {
    if (!setop_i64(a)) return NULL;
    int64_t* out; size_t m;
    if (!raw_delete_duplicates((const int64_t*)a->data.ndarray.data, (size_t)a->data.ndarray.dims[0], &out, &m))
        return NULL;
    return setop_emit(a, out, m);
}

static Expr* setop_packed(Expr** args, size_t nl, int mode) {
    if (nl == 0) return NULL;
    for (size_t i = 0; i < nl; i++) if (!setop_i64(args[i])) return NULL;
    I64Span* sp = malloc(sizeof(I64Span) * nl);
    if (!sp) return NULL;
    for (size_t i = 0; i < nl; i++) {
        sp[i].v = (const int64_t*)args[i]->data.ndarray.data;
        sp[i].n = (size_t)args[i]->data.ndarray.dims[0];
    }
    int64_t* out; size_t m;
    bool ok = raw_setop(sp, nl, mode, &out, &m);
    free(sp);
    return ok ? setop_emit(args[0], out, m) : NULL;
}

/* True when any of args[0..n) is an ndarray -- the guard every generic set
 * operation needs, because the code below would read an NDArray as an atom. */
static bool setop_any_nd(Expr** args, size_t n) {
    for (size_t i = 0; i < n; i++) if (is_ndarray(args[i])) return true;
    return false;
}

/* ---- typed LinkedList wrappers ----
 * The int64 slots of LinkedList operands gathered into scratch spans (raw
 * buffers, not Lists), for the raw cores above. False unless every operand is
 * a LinkedList of machine integers (or empty). */
static bool setop_ll_spans(Expr** args, size_t nl, I64Span* sp) {
    for (size_t i = 0; i < nl; i++) {
        const LinkedList* ll = expr_linkedlist(args[i]);
        if (!ll || (ll_kind(ll) != LL_I64 && ll_length(ll) != 0)) return false;
    }
    for (size_t i = 0; i < nl; i++) {
        const LinkedList* ll = args[i]->data.llist;
        size_t n = (size_t)ll_length(ll);
        int64_t* b = malloc(sizeof(int64_t) * (n ? n : 1));
        if (!b) { while (i-- > 0) free((void*)sp[i].v); return false; }
        if (n) ll_gather(ll, 0, (int64_t)n, b);
        sp[i].v = b; sp[i].n = n;
    }
    return true;
}

static void setop_ll_spans_free(I64Span* sp, size_t nl) {
    for (size_t i = 0; i < nl; i++) free((void*)sp[i].v);
}

/* An owned int64 buffer as a new LinkedList (consumes v). */
static Expr* setop_ll_emit(int64_t* v, size_t m) {
    LinkedList* ll = m ? ll_from_buffer(LL_I64, v, (int64_t)m) : ll_new();
    free(v);
    return expr_new_linkedlist(ll);
}

/* Union / Intersection / Complement of LinkedLists of machine integers on the
 * packed kernels; NULL for any other operand shape. */
static Expr* setop_ll_typed(Expr** args, size_t nl, int mode) {
    if (nl == 0) return NULL;
    I64Span* sp = malloc(sizeof(I64Span) * nl);
    if (!sp) return NULL;
    if (!setop_ll_spans(args, nl, sp)) { free(sp); return NULL; }
    int64_t* out; size_t m;
    bool ok = raw_setop(sp, nl, mode, &out, &m);
    setop_ll_spans_free(sp, nl);
    free(sp);
    return ok ? setop_ll_emit(out, m) : NULL;
}

/* ---- the operands of a generic set operation -------------------------------
 * Lists of one head, or LinkedLists (each standing for the List of its
 * elements; the answer is a LinkedList when any operand is one), read element
 * by element. The operands' elements are numbered consecutively ("global"
 * indices) so one index set or one sort can span them all. */
typedef struct { LLSeq* q; size_t n; size_t* off; bool linked; const Expr* head; } SetOps;

/* 0 with every operand open, the 1-based position of an atom operand, or
 * minus the 1-based position of the first operand whose head differs from
 * the first one's (nothing left open either way). */
static int setops_open(Expr** args, size_t nl, SetOps* so) {
    so->q = malloc(sizeof(LLSeq) * (nl ? nl : 1));
    so->off = malloc(sizeof(size_t) * (nl + 1));
    so->n = 0; so->linked = false; so->head = NULL;
    int bad = 0;
    for (size_t i = 0; i < nl && !bad; i++) {
        if (!llq_open(args[i], &so->q[i])) { bad = (int)i + 1; break; }
        so->n = i + 1;
        if (llq_linked(&so->q[i])) so->linked = true;
        else if (!so->head) so->head = so->q[i].head;
        else if (!expr_eq((Expr*)so->head, (Expr*)so->q[i].head)) bad = -(int)(i + 1);
    }
    if (!bad && so->linked && so->head &&
        !(so->head->type == EXPR_SYMBOL && so->head->data.symbol.name == SYM_List)) {
        /* a LinkedList (a List) beside an operand of another head */
        bool first_linked = llq_linked(&so->q[0]);
        for (size_t i = 1; i < nl && !bad; i++)
            if (llq_linked(&so->q[i]) != first_linked) bad = -(int)(i + 1);
        if (!bad) bad = -2;
    }
    if (bad) {
        for (size_t i = 0; i < so->n; i++) llq_close(&so->q[i]);
        free(so->q); free(so->off);
        so->q = NULL; so->off = NULL;
        return bad;
    }
    so->off[0] = 0;
    for (size_t i = 0; i < nl; i++) so->off[i + 1] = so->off[i] + so->q[i].n;
    return 0;
}

static void setops_close(SetOps* so) {
    for (size_t i = 0; i < so->n; i++) llq_close(&so->q[i]);
    free(so->q); free(so->off);
}

static size_t setops_total(const SetOps* so) { return so->off[so->n]; }

/* The operand holding global index g, and g's index within it. */
static size_t setops_which(const SetOps* so, size_t g, size_t* local) {
    size_t lo = 0, hi = so->n;                  /* off[lo] <= g < off[hi] */
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (so->off[mid] <= g) lo = mid; else hi = mid;
    }
    *local = g - so->off[lo];
    return lo;
}

static Expr* setops_peek(void* ctx, size_t g) {
    SetOps* so = ctx;
    size_t local, k = setops_which(so, g, &local);
    return llq_peek(&so->q[k], local);
}

/* The answer holding the elements at global indices g[0..k): a LinkedList
 * copied slot to slot, else head[copies]. */
static Expr* setops_emit(SetOps* so, const size_t* g, size_t k) {
    if (so->linked) {
        LinkedList* out = ll_new();
        for (size_t t = 0; t < k; t++) {
            size_t local, w = setops_which(so, g[t], &local);
            llq_push(out, &so->q[w], local);
        }
        return expr_new_linkedlist(out);
    }
    Expr** items = malloc(sizeof(Expr*) * (k ? k : 1));
    for (size_t t = 0; t < k; t++) items[t] = expr_copy(setops_peek(so, g[t]));
    Expr* r = expr_new_function(expr_copy((Expr*)so->head), items, k);
    free(items);
    return r;
}

/* qsort over global indices by canonical order (file-static context, saved
 * and restored: expr_compare runs no user code, but stay reentrant). */
static SetOps* g_cmp_so = NULL;
static int setops_idx_cmp(const void* a, const void* b) {
    size_t ia = *(const size_t*)a, ib = *(const size_t*)b;
    Expr* ea = expr_copy(setops_peek(g_cmp_so, ia));
    int c = expr_compare(ea, setops_peek(g_cmp_so, ib));
    expr_free(ea);
    return c;
}
static void setops_sort(SetOps* so, size_t* idx, size_t k) {
    SetOps* saved = g_cmp_so;
    g_cmp_so = so;
    qsort(idx, k, sizeof(size_t), setops_idx_cmp);
    g_cmp_so = saved;
}

/* Evaluate test[a, b] and report whether it yields True (a and b borrowed;
 * the call is built from fresh copies). */
static bool same_test_equal(Expr* test, Expr* a, Expr* b) {
    Expr* call_args[2] = { expr_copy(a), expr_copy(b) };
    Expr* call = expr_new_function(expr_copy(test), call_args, 2);
    Expr* r = eval_callback(call);
    bool eq = (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    expr_free(r);
    expr_free(call);
    return eq;
}

/* test[element ga, element gb] of the operands. */
static bool setops_test(SetOps* so, Expr* test, size_t ga, size_t gb) {
    Expr* a = expr_copy(setops_peek(so, ga));
    Expr* b = expr_copy(setops_peek(so, gb));
    bool eq = same_test_equal(test, a, b);
    expr_free(a); expr_free(b);
    return eq;
}

/* Head::heads for operands 1 and k (k >= 2) of `res` whose heads differ:
 * "Heads List and f at positions 1 and 2 are expected to be the same."
 * (audits Union#3, Intersection#2, Complement#2, Join#3). */
static Expr* setops_heads_error(const char* name, Expr* res, size_t k) {
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[k - 1];
    Expr* sl = expr_new_symbol(SYM_List);
    char* ha = expr_to_string_outputform(a->type == EXPR_FUNCTION ? a->data.function.head : sl);
    char* hb = expr_to_string_outputform(b->type == EXPR_FUNCTION ? b->data.function.head : sl);
    expr_free(sl);
    mth_message(name, "heads", "Heads %s and %s at positions 1 and %zu are expected to be the same.",
                ha ? ha : "", hb ? hb : "", k);
    free(ha); free(hb);
    return NULL;
}

/* A trailing SameTest -> f (or SameTest :> f) option: *last_arg is where the
 * operands end (audit Complement#3: the delayed rule was not seen). */
static Expr* setops_same_test(Expr* res, size_t* last_arg) {
    Expr* same_test = NULL;
    *last_arg = res->data.function.arg_count;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            (arg->data.function.head->data.symbol.name == SYM_Rule ||
             arg->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
            arg->data.function.arg_count == 2 &&
            arg->data.function.args[0]->type == EXPR_SYMBOL &&
            arg->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            same_test = arg->data.function.args[1];
            if (i < *last_arg) *last_arg = i;
        }
    }
    /* SameTest -> Automatic means "use the default comparison" (audit
     * Intersection#1: it was taken as a test that never held). */
    if (same_test && same_test->type == EXPR_SYMBOL && same_test->data.symbol.name == SYM_Automatic)
        same_test = NULL;
    return same_test;
}

/* Union[l1, l2, ...]: the sorted distinct elements of all the li. Default
 * comparison: the hash set over every element (O(total)), then a sort of the
 * distinct ones only (issue #112: 946k elements reduced to 91). SameTest -> f:
 * sort everything, then drop an element f calls equal to ANY element already
 * kept -- not only the last one, since f need not respect the canonical order
 * (Union[Range[6], SameTest -> (Mod[#1 - #2, 3] == 0 &)] is {1, 2, 3}; audit
 * Union#1). */
Expr* builtin_union(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    size_t last_arg;
    Expr* same_test = setops_same_test(res, &last_arg);
    if (last_arg == 0) return NULL;

    /* Packed operands: take the buffer merge, or degrade. */
    if (setop_any_nd(res->data.function.args, last_arg)) {
        if (!same_test) {
            Expr* fast = setop_packed(res->data.function.args, last_arg, 0);
            if (fast) return fast;
        }
        return ndarray_delist_and_reeval(res);
    }
    if (!same_test) {               /* LinkedLists of machine integers: the same kernel */
        Expr* fast = setop_ll_typed(res->data.function.args, last_arg, 0);
        if (fast) return fast;
    }

    SetOps so;
    int bad = setops_open(res->data.function.args, last_arg, &so);
    if (bad > 0) return builtin_normal_error("Union", res, bad);
    if (bad < 0) return setops_heads_error("Union", res, (size_t)(-bad)); /* heads must match */

    size_t total = setops_total(&so);
    if (total == 0) {
        Expr* r = so.linked ? setops_emit(&so, NULL, 0) : expr_copy(res->data.function.args[0]);
        setops_close(&so);
        return r;
    }
    size_t* keep = malloc(sizeof(size_t) * total);
    size_t nk = 0;
    if (same_test == NULL) {
        LLIdxSet seen;
        llis_init(&seen, total);
        for (size_t g = 0; g < total; g++)
            if (llis_add(&seen, setops_peek, &so, g) == SIZE_MAX) keep[nk++] = g;
        llis_free(&seen);
        setops_sort(&so, keep, nk);
    } else {
        size_t* all = malloc(sizeof(size_t) * total);
        for (size_t g = 0; g < total; g++) all[g] = g;
        setops_sort(&so, all, total);
        for (size_t i = 0; i < total; i++) {
            bool dup = false;
            for (size_t k = nk; k > 0 && !dup; k--)      /* the latest kept first */
                dup = setops_test(&so, same_test, all[i], keep[k - 1]);
            if (!dup) keep[nk++] = all[i];
        }
        free(all);
    }
    Expr* result = setops_emit(&so, keep, nk);
    free(keep);
    setops_close(&so);
    return result;
}

/* Intersection[l1, l2, ...] (mode 1) and Complement[eall, e1, ...] (mode 2):
 * the sorted distinct elements of the first operand present in every later
 * one / absent from all of them, under the head of the first.
 *
 * Default comparison: canonical structural equality through hash sets,
 * O(total). SameTest -> f is O(n^2) with f[a,b]===True as the equivalence;
 * Intersection keeps the canonically-GREATEST member of each class as its
 * representative, Complement the SMALLEST (the reference's choices), so the
 * first operand is walked in ascending canonical order. */
static Expr* setops_filter(Expr* res, int mode, const char* name) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1)
        return builtin_arg_error(name, (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0,
                                 1, SIZE_MAX);
    size_t last_arg;
    Expr* same_test = setops_same_test(res, &last_arg);
    if (last_arg == 0) return NULL;   /* only an option, no list operands */

    if (setop_any_nd(res->data.function.args, last_arg)) {
        if (!same_test) {
            Expr* fast = setop_packed(res->data.function.args, last_arg, mode);
            if (fast) return fast;
        }
        return ndarray_delist_and_reeval(res);
    }
    if (!same_test) {
        Expr* fast = setop_ll_typed(res->data.function.args, last_arg, mode);
        if (fast) return fast;
    }

    /* Every operand must share the head of the first; an atom is ::normal. */
    SetOps so;
    int bad = setops_open(res->data.function.args, last_arg, &so);
    if (bad > 0) return builtin_normal_error(name, res, bad);
    if (bad < 0) return setops_heads_error(name, res, (size_t)(-bad));

    size_t n0 = so.q[0].n;              /* operand 0 holds global indices 0..n0-1 */
    size_t* cand = malloc(sizeof(size_t) * (n0 ? n0 : 1));
    size_t ncand = 0;

    if (same_test == NULL) {
        LLIdxSet seen;
        llis_init(&seen, n0);
        for (size_t j = 0; j < n0; j++)
            if (llis_add(&seen, llq_peek_fn, &so.q[0], j) == SIZE_MAX) cand[ncand++] = j;
        llis_free(&seen);
        for (size_t i = 1; i < last_arg && ncand > 0; i++) {
            LLSeq* qi = &so.q[i];
            LLIdxSet h;
            llis_init(&h, qi->n);
            for (size_t j = 0; j < qi->n; j++) llis_add(&h, llq_peek_fn, qi, j);
            size_t w = 0;
            for (size_t k = 0; k < ncand; k++) {
                bool present = llis_find(&h, llq_peek_fn, qi, llq_peek(&so.q[0], cand[k])) != SIZE_MAX;
                if (present == (mode == 1)) cand[w++] = cand[k];
            }
            ncand = w;
            llis_free(&h);
        }
    } else {
        size_t* sorted0 = malloc(sizeof(size_t) * (n0 ? n0 : 1));
        for (size_t j = 0; j < n0; j++) sorted0[j] = j;
        setops_sort(&so, sorted0, n0);
        for (size_t j = 0; j < n0; j++) {
            size_t match = SIZE_MAX;
            for (size_t k = 0; k < ncand; k++)
                if (setops_test(&so, same_test, sorted0[j], cand[k])) { match = k; break; }
            if (match == SIZE_MAX) cand[ncand++] = sorted0[j];
            else if (mode == 1) cand[match] = sorted0[j];   /* the greater one */
        }
        free(sorted0);
        for (size_t i = 1; i < last_arg && ncand > 0; i++) {
            size_t w = 0;
            for (size_t k = 0; k < ncand; k++) {
                bool present = false;
                for (size_t g = so.off[i]; g < so.off[i + 1] && !present; g++)
                    present = setops_test(&so, same_test, cand[k], g);
                if (present == (mode == 1)) cand[w++] = cand[k];
            }
            ncand = w;
        }
    }

    setops_sort(&so, cand, ncand);
    Expr* result = setops_emit(&so, cand, ncand);
    free(cand);
    setops_close(&so);
    return result;
}

Expr* builtin_intersection(Expr* res) { return setops_filter(res, 1, "Intersection"); }
Expr* builtin_complement(Expr* res)   { return setops_filter(res, 2, "Complement"); }

/* ---- Tally / Commonest / DeleteDuplicates --------------------------------- */

/* The distinct elements of q (first appearance) and their multiplicities:
 * *uniq holds element INDICES. With a test, element i joins the first class
 * whose representative r gives test[i, r] === True (O(n g)). */
static size_t tally_indices(LLSeq* q, Expr* test, size_t** uniq, int64_t** mult) {
    size_t n = q->n;
    *uniq = malloc(sizeof(size_t) * (n ? n : 1));
    *mult = malloc(sizeof(int64_t) * (n ? n : 1));
    size_t nu = 0;
    if (!test) {
        LLIdxSet seen;
        llis_init(&seen, n);
        size_t* cls = malloc(sizeof(size_t) * (n ? n : 1));   /* element index -> class */
        for (size_t i = 0; i < n; i++) {
            size_t j = llis_add(&seen, llq_peek_fn, q, i);
            if (j == SIZE_MAX) { cls[i] = nu; (*uniq)[nu] = i; (*mult)[nu] = 1; nu++; }
            else { cls[i] = cls[j]; (*mult)[cls[j]]++; }
        }
        free(cls);
        llis_free(&seen);
        return nu;
    }
    for (size_t i = 0; i < n; i++) {
        size_t found = SIZE_MAX;
        for (size_t u = 0; u < nu && found == SIZE_MAX; u++) {
            Expr* a = expr_copy(llq_peek(q, i));
            Expr* b = expr_copy(llq_peek(q, (*uniq)[u]));
            if (same_test_equal(test, a, b)) found = u;
            expr_free(a); expr_free(b);
        }
        if (found != SIZE_MAX) (*mult)[found]++;
        else { (*uniq)[nu] = i; (*mult)[nu] = 1; nu++; }
    }
    return nu;
}

/* Tally of a typed LinkedList on the packed Tally's own kernel (the gathered
 * slots as machine words): a LinkedList of {value, count} Lists. NULL when the
 * list is boxed / complex / empty or a word is not keyable (non-finite). */
static Expr* tally_ll_machine(const LinkedList* ll) {
    LLKind k = ll_kind(ll);
    int64_t n = ll_length(ll);
    if (n == 0 || (k != LL_I64 && k != LL_REAL)) return NULL;
    void* buf = malloc(sizeof(double) * (size_t)n);
    if (!buf) return NULL;
    ll_gather(ll, 0, n, buf);
    uint64_t* keys; int64_t* cnts; size_t nu;
    bool ok = ndred_tally_words_raw(buf, k == LL_I64 ? NDT_INT64 : NDT_FLOAT64, (size_t)n, &keys, &cnts, &nu);
    free(buf);
    if (!ok) return NULL;
    LinkedList* out = ll_new();
    for (size_t u = 0; u < nu; u++) {
        Expr* pair[2];
        if (k == LL_I64) pair[0] = expr_new_integer((int64_t)keys[u]);
        else { double x; memcpy(&x, &keys[u], sizeof x); pair[0] = expr_new_real(x); }
        pair[1] = expr_new_integer(cnts[u]);
        ll_push_back(out, expr_new_function(expr_new_symbol(SYM_List), pair, 2));
    }
    free(keys); free(cnts);
    return expr_new_linkedlist(out);
}

Expr* builtin_tally(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    /* Tally[assoc] (and Tally[assoc, test]) tallies the association's values. */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* A packed buffer hashes as machine words; a user test needs expressions. */
    if (is_ndarray(list))
        return test ? ndarray_delist_and_reeval(res) : ndred_tally(res);
    if (!test && list->type == EXPR_LINKEDLIST) {
        Expr* fast = tally_ll_machine(list->data.llist);
        if (fast) return fast;
    }

    /* A LinkedList tallies its elements read from its slots; the answer is a
     * LinkedList of {element, count} Lists. */
    LLSeq q;
    if (!llq_open(list, &q)) {
        /* Tally[5]: an atom has no elements to tally (audit Tally#1) */
        char* cs = expr_to_string(res);
        mth_message("Tally", "list", "List expected at position 1 in %s.", cs ? cs : "");
        free(cs);
        return NULL;
    }
    bool linked = llq_linked(&q);
    size_t* uniq; int64_t* mult;
    size_t nu = tally_indices(&q, test, &uniq, &mult);
    Expr** pairs = malloc(sizeof(Expr*) * (nu ? nu : 1));
    for (size_t u = 0; u < nu; u++) {
        Expr* pair[2] = { expr_copy(llq_peek(&q, uniq[u])), expr_new_integer(mult[u]) };
        pairs[u] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
    }
    Expr* result = ll_seq_of(linked, pairs, nu);
    free(pairs); free(uniq); free(mult);
    llq_close(&q);
    return result;
}

/* DeleteDuplicates of a typed LinkedList of machine integers on the packed
 * kernel; NULL otherwise. */
static Expr* dedup_ll_machine(const LinkedList* ll) {
    if (ll_kind(ll) != LL_I64) return NULL;
    size_t n = (size_t)ll_length(ll);
    int64_t* v = malloc(sizeof(int64_t) * (n ? n : 1));
    if (!v) return NULL;
    ll_gather(ll, 0, (int64_t)n, v);
    int64_t* out; size_t m;
    bool ok = raw_delete_duplicates(v, n, &out, &m);
    free(v);
    return ok ? setop_ll_emit(out, m) : NULL;
}

Expr* builtin_deleteduplicates(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    /* DeleteDuplicates[assoc] keeps the first entry for each distinct value. */
    if (test == NULL && is_association(list)) return assoc_delete_duplicate_values(list);
    /* ... and DeleteDuplicates[assoc, test] for each class of values under
     * the test (audit DeleteDuplicates#2: the test was ignored). */
    if (test != NULL && is_association(list)) {
        Expr* vals = assoc_values_list(list);
        LLSeq vq;
        if (!vals || !llq_open(vals, &vq)) { if (vals) expr_free(vals); return NULL; }
        size_t* uniq; int64_t* mult;
        size_t nu = tally_indices(&vq, test, &uniq, &mult);
        llq_close(&vq);
        Expr** ents = malloc(sizeof(Expr*) * (nu ? nu : 1));
        for (size_t u = 0; u < nu; u++) ents[u] = list->data.function.args[uniq[u]];
        Expr* out = assoc_from_rules(ents, nu);
        free(ents); free(uniq); free(mult);
        expr_free(vals);
        return out;
    }

    /* Packed int64 with the default test: dedup on the machine words. Every
     * other array shape must delist here, not fall through (it would be read
     * as an atom and returned unchanged). */
    if (is_ndarray(list)) {
        if (test == NULL) {
            Expr* fast = setop_packed_delete_duplicates(list);
            if (fast) return fast;
        }
        return ndarray_delist_and_reeval(res);
    }
    if (test == NULL && list->type == EXPR_LINKEDLIST) {
        Expr* fast = dedup_ll_machine(list->data.llist);
        if (fast) return fast;
    }

    LLSeq q;
    if (!llq_open(list, &q)) return builtin_normal_error("DeleteDuplicates", res, 1);
    if (q.n == 0 && !llq_linked(&q)) { llq_close(&q); return expr_copy(list); }
    size_t* uniq; int64_t* mult;
    size_t nu = tally_indices(&q, test, &uniq, &mult);
    Expr* result = llq_emit(&q, uniq, nu);
    free(uniq); free(mult);
    llq_close(&q);
    return result;
}

/* DeleteDuplicatesBy[expr, f] keeps the first element for each distinct value
 * of f[element], preserving order. Over an association f is applied to each
 * value and the surviving entries form an association. The distinct f-values
 * (keys, not elements) are compared directly. */
static Expr* dedupby_key_peek(void* ctx, size_t i) { return ((Expr**)ctx)[i]; }

Expr* builtin_deleteduplicatesby(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    bool assoc = is_association(coll);
    LLSeq q;
    if (!llq_open(coll, &q)) return NULL;
    size_t n = q.n;
    size_t* kept = malloc(sizeof(size_t) * (n ? n : 1));  /* surviving element indices */
    Expr** keys = malloc(sizeof(Expr*) * (n ? n : 1));    /* f[element i] (owned) */
    size_t nkept = 0;
    /* The f-values are deduplicated through the SameQ hash set, so the pass is
     * O(n) in the number of distinct keys, not O(n^2) (audit
     * DeleteDuplicatesBy#1: 10^5 distinct elements took over 10 s). */
    LLIdxSet seen;
    llis_init(&seen, n);
    for (size_t i = 0; i < n; i++) {
        Expr* elem = llq_peek(&q, i);
        Expr* val  = assoc ? elem->data.function.args[1] : elem;  /* value for assoc */
        Expr* fcall = expr_new_function(expr_copy(f), (Expr*[]){ expr_copy(val) }, 1);
        keys[i] = eval_callback(fcall);
        expr_free(fcall);
        if (llis_add(&seen, dedupby_key_peek, keys, i) == SIZE_MAX) kept[nkept++] = i;
    }
    llis_free(&seen);
    Expr* out = llq_emit(&q, kept, nkept);
    for (size_t j = 0; j < n; j++) expr_free(keys[j]);
    free(keys); free(kept);
    llq_close(&q);
    return out;
}

typedef struct {
    size_t elem;          /* the element's index */
    int64_t count;
    size_t first_index;
} CommonestItem;

static int compare_commonest_items_desc(const void* a, const void* b) {
    const CommonestItem* item_a = (const CommonestItem*)a;
    const CommonestItem* item_b = (const CommonestItem*)b;
    if (item_a->count != item_b->count) {
        return (item_b->count > item_a->count) ? 1 : -1;
    }
    return (item_a->first_index > item_b->first_index) ? 1 : -1;
}

static int compare_commonest_items_index(const void* a, const void* b) {
    const CommonestItem* item_a = (const CommonestItem*)a;
    const CommonestItem* item_b = (const CommonestItem*)b;
    if (item_a->first_index == item_b->first_index) return 0;
    return (item_a->first_index > item_b->first_index) ? 1 : -1;
}

Expr* builtin_commonest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* list = res->data.function.args[0];

    /* Commonest[assoc] (and Commonest[assoc, n]) uses the association's values. */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* A packed buffer counts as machine words (ndred_commonest shares Tally's
     * counting routine, so count ties break the same way). */
    if (is_ndarray(list)) return ndred_commonest(res);

    Expr* n_arg = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    int64_t n = 0;
    bool have_n = (n_arg != NULL);
    bool n_upto = false;
    if (n_arg) {
        if (n_arg->type == EXPR_INTEGER) {
            n = n_arg->data.integer;
        } else if (n_arg->type == EXPR_FUNCTION && n_arg->data.function.head->type == EXPR_SYMBOL &&
                   n_arg->data.function.head->data.symbol.name == SYM_UpTo && n_arg->data.function.arg_count == 1) {
            if (n_arg->data.function.args[0]->type == EXPR_INTEGER) {
                n = n_arg->data.function.args[0]->data.integer;
                n_upto = true;
            } else return NULL;
        } else return NULL;
    }

    LLSeq q;
    if (!llq_open(list, &q)) {
        /* Commonest[x]: Commonest::list, unevaluated, as Tally (audit Commonest#2). */
        if (!expr_is_plain_atom(list)) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        char* cs = expr_to_string(res);
        mth_message("Commonest", "list", "List expected at position 1 in %s.", cs ? cs : "");
        free(cs);
        return NULL;
    }
    bool linked = llq_linked(&q);
    if (q.n == 0) { llq_close(&q); return ll_seq_of(linked, NULL, 0); }

    /* The tally: a typed LinkedList's on the packed kernel (the same first-
     * appearance order, so ties break as on the packed List), else by index. */
    size_t* uniq = NULL; int64_t* mult = NULL; size_t unique_count;
    uint64_t* words = NULL;
    LLKind kind = llq_kind(&q);
    if (linked && (kind == LL_I64 || kind == LL_REAL)) {
        void* buf = malloc(sizeof(double) * q.n);
        ll_gather(q.snap, 0, (int64_t)q.n, buf);
        if (!ndred_tally_words_raw(buf, kind == LL_I64 ? NDT_INT64 : NDT_FLOAT64, q.n, &words, &mult,
                                   &unique_count))
            words = NULL;
        free(buf);
    }
    if (!words) unique_count = tally_indices(&q, NULL, &uniq, &mult);

    CommonestItem* items = malloc(sizeof(CommonestItem) * unique_count);
    for (size_t i = 0; i < unique_count; i++) {
        items[i].elem = uniq ? uniq[i] : i;      /* word i, for the machine tally */
        items[i].count = mult[i];
        items[i].first_index = i;
    }
    free(mult); free(uniq);

    qsort(items, unique_count, sizeof(CommonestItem), compare_commonest_items_desc);

    /* `have_n` rather than a sentinel value of n: Commonest[list, -1] is {}. */
    size_t target_n;
    if (!have_n) {
        int64_t max_count = items[0].count;
        target_n = 0;
        while (target_n < unique_count && items[target_n].count == max_count) target_n++;
    } else {
        if (n < 0) n = 0;
        if ((size_t)n > unique_count) {
            if (!n_upto) {
                printf("Commonest::dstlms: The requested number of elements %" PRId64 " is greater than the number of distinct elements %zu. Only %zu elements will be returned.\n", n, unique_count, unique_count);
            }
            target_n = unique_count;
        } else {
            target_n = (size_t)n;
        }
    }
    if (target_n > 0) qsort(items, target_n, sizeof(CommonestItem), compare_commonest_items_index);

    Expr* out;
    if (words) {
        LinkedList* ol = ll_new();
        for (size_t i = 0; i < target_n; i++) {
            uint64_t w = words[items[i].elem];
            if (kind == LL_I64) ll_push_back_i64(ol, (int64_t)w);
            else { double x; memcpy(&x, &w, sizeof x); ll_push_back_real(ol, x); }
        }
        out = expr_new_linkedlist(ol);
        free(words);
    } else if (linked) {
        size_t* idx = malloc(sizeof(size_t) * (target_n ? target_n : 1));
        for (size_t i = 0; i < target_n; i++) idx[i] = items[i].elem;
        out = llq_emit(&q, idx, target_n);
        free(idx);
    } else {
        Expr** result_args = malloc(sizeof(Expr*) * (target_n ? target_n : 1));
        for (size_t i = 0; i < target_n; i++) result_args[i] = expr_copy(llq_peek(&q, items[i].elem));
        out = expr_new_function(expr_new_symbol(SYM_List), result_args, target_n);
        free(result_args);
    }
    free(items);
    llq_close(&q);
    return out;
}

/* ---------------- IntersectingQ / DuplicateFreeQ ---------------- */

/* Length of a set-op fast-path result: a packed list, or the plain {} the
 * fast paths return for an empty answer. */
static size_t setop_result_len(const Expr* r) {
    return is_ndarray(r) ? (size_t)r->data.ndarray.dims[0] : r->data.function.arg_count;
}

/* "Nonatomic expression expected at position 1 in call." style messages. */
static void setop_call_message(const char* head, const char* tag, const char* what, Expr* call) {
    char* txt = expr_to_string_outputform(call);
    mth_message(head, tag, "%s in %s.", what, txt ? txt : head);
    free(txt);
}

/* IntersectingQ[a, b] / IntersectingQ[a, b, SameTest -> f]
 *
 * True when a and b share at least one element: SameQ by default (so 1 and
 * 1. do not), f[x, y] === True under SameTest. Both must be non-atomic with
 * the same head (a LinkedList stands for a List). The default test hashes b
 * and probes each element of a; two int64 buffers (packed, or a typed
 * LinkedList's slots) intersect on the machine words. */
Expr* builtin_intersectingq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return builtin_arg_error("IntersectingQ", argc, 2, 2);
    Expr* test = NULL;
    if (argc == 3) {
        Expr* o = res->data.function.args[2];
        if (o->type != EXPR_FUNCTION || o->data.function.arg_count != 2
            || o->data.function.head->type != EXPR_SYMBOL
            || (o->data.function.head->data.symbol.name != SYM_Rule
                && o->data.function.head->data.symbol.name != SYM_RuleDelayed)
            || o->data.function.args[0]->type != EXPR_SYMBOL
            || o->data.function.args[0]->data.symbol.name != SYM_SameTest)
            return NULL;
        test = o->data.function.args[1];
        if (test->type == EXPR_SYMBOL && test->data.symbol.name == SYM_Automatic) test = NULL;
    }
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    if (is_ndarray(a) || is_ndarray(b)) {
        if (!test && setop_i64(a) && setop_i64(b)) {
            Expr* pair[2] = { a, b };
            Expr* common = setop_packed(pair, 2, 1);
            if (common) {
                bool any = setop_result_len(common) > 0;
                expr_free(common);
                return expr_new_symbol(any ? SYM_True : SYM_False);
            }
        }
        return ndarray_delist_and_reeval(res);
    }
    if (!test) {
        Expr* common = setop_ll_typed(res->data.function.args, 2, 1);
        if (common) {
            bool any = ll_length(common->data.llist) > 0;
            expr_free(common);
            return expr_new_symbol(any ? SYM_True : SYM_False);
        }
    }
    for (size_t i = 0; i < 2; i++) {
        Expr* x = res->data.function.args[i];
        if (x->type != EXPR_FUNCTION && x->type != EXPR_LINKEDLIST) {
            char* txt = expr_to_string_outputform(x);
            mth_message("IntersectingQ", "lsa", "List or association expected instead of %s.",
                        txt ? txt : "");
            free(txt);
            return NULL;
        }
    }
    SetOps so;
    if (setops_open(res->data.function.args, 2, &so) != 0) {
        Expr* sl = expr_new_symbol(SYM_List);
        char* ha = expr_to_string_outputform(a->type == EXPR_FUNCTION ? a->data.function.head : sl);
        char* hb = expr_to_string_outputform(b->type == EXPR_FUNCTION ? b->data.function.head : sl);
        expr_free(sl);
        mth_message("IntersectingQ", "heads",
                    "Heads %s and %s at positions 1 and 2 are expected to be the same.",
                    ha ? ha : "", hb ? hb : "");
        free(ha); free(hb);
        return NULL;
    }
    size_t na = so.q[0].n, nb = so.q[1].n;
    bool found = false;
    if (test) {
        for (size_t i = 0; i < na && !found; i++)
            for (size_t j = 0; j < nb && !found; j++)
                found = setops_test(&so, test, i, so.off[1] + j);
    } else if (na && nb) {
        LLIdxSet h;
        llis_init(&h, nb);
        for (size_t j = 0; j < nb; j++) llis_add(&h, llq_peek_fn, &so.q[1], j);
        for (size_t i = 0; i < na && !found; i++)
            found = llis_find(&h, llq_peek_fn, &so.q[1], llq_peek(&so.q[0], i)) != SIZE_MAX;
        llis_free(&h);
    }
    setops_close(&so);
    return expr_new_symbol(found ? SYM_True : SYM_False);
}

/* DuplicateFreeQ[list] / DuplicateFreeQ[list, test]
 *
 * True when no two elements are the same: SameQ by default, test[x, y] ===
 * True otherwise. Any head is accepted; an atom gives DuplicateFreeQ::normal.
 * The default test is one hash pass that stops at the first repeat; an int64
 * buffer compares its deduplicated length on the machine words. */
Expr* builtin_duplicatefreeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error("DuplicateFreeQ", argc, 1, 2);
    Expr* list = res->data.function.args[0];
    Expr* test = argc == 2 ? res->data.function.args[1] : NULL;

    if (is_ndarray(list)) {
        if (!test && setop_i64(list)) {
            Expr* uniq = setop_packed_delete_duplicates(list);
            if (uniq) {
                bool free_q = setop_result_len(uniq) == (size_t)list->data.ndarray.dims[0];
                expr_free(uniq);
                return expr_new_symbol(free_q ? SYM_True : SYM_False);
            }
        }
        return ndarray_delist_and_reeval(res);
    }
    if (!test && list->type == EXPR_LINKEDLIST) {
        Expr* uniq = dedup_ll_machine(list->data.llist);
        if (uniq) {
            bool free_q = ll_length(uniq->data.llist) == ll_length(list->data.llist);
            expr_free(uniq);
            return expr_new_symbol(free_q ? SYM_True : SYM_False);
        }
    }
    LLSeq q;
    if (!llq_open(list, &q)) {
        setop_call_message("DuplicateFreeQ", "normal", "Nonatomic expression expected at position 1", res);
        return NULL;
    }
    size_t n = q.n;
    bool dup = false;
    if (test) {
        for (size_t i = 1; i < n && !dup; i++)
            for (size_t j = 0; j < i && !dup; j++) {
                Expr* x = expr_copy(llq_peek(&q, j));
                Expr* y = expr_copy(llq_peek(&q, i));
                dup = same_test_equal(test, x, y);
                expr_free(x); expr_free(y);
            }
    } else if (n > 1) {
        LLIdxSet h;
        llis_init(&h, n);
        for (size_t i = 0; i < n && !dup; i++) dup = llis_add(&h, llq_peek_fn, &q, i) != SIZE_MAX;
        llis_free(&h);
    }
    llq_close(&q);
    return expr_new_symbol(dup ? SYM_False : SYM_True);
}
