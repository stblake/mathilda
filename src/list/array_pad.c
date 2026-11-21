#include "list_common.h"
#include "array_pad.h"
#include "pad_schemes.h"
#include "options.h"
#include "ndarray.h"     /* is_ndarray, ndarray_to_nested_list */
#include "ndstruct.h"    /* ndstruct_arraypad */
#include "arithmetic.h"  /* is_infinity_sym */
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "linkedlist_seq.h" /* LLSeq: LinkedList levels pad into LinkedLists */

/* ArrayPad[array, m]            pad m elements of padding on every side/level.
 * ArrayPad[array, {m, n}]       m before, n after, on every dimension.
 * ArrayPad[array, {{m1,n1},...}] per-level amounts ({mi} means {mi,mi}).
 * ArrayPad[array, amounts, padding]  use the given padding scheme (default 0).
 *
 * A negative amount removes elements from that side.  Padding may be a constant,
 * a cyclic list of constants, or one of the named schemes handled by
 * pad_schemes.c ("Fixed"/"Periodic"/"Reflected"/"Reversed"/"ReversedNegation"/
 * "ReflectedDifferences"/"ReversedDifferences"/"Extrapolated", the last taking
 * InterpolationOrder). */

#define AP_MAX_RANK 64

/* Nesting depth of a nested-List array (0 for an atom). {{},{}} has depth 2. */
static bool ap_is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* A level of the array: a List, or a LinkedList (read through an LLSeq). */
static bool ap_row(const Expr* e) { return ap_is_list(e) || e->type == EXPR_LINKEDLIST; }

static size_t ap_depth(const Expr* e) {
    if (!ap_row(e)) return 0;
    LLSeq v;
    llq_open(e, &v);
    size_t best = 0;
    for (size_t i = 0; i < v.n; i++) {
        size_t d = ap_depth(llq_peek(&v, i));
        if (d > best) best = d;
    }
    llq_close(&v);
    return 1 + best;
}

/* Rectangular dimension along each level, following the first-child spine. */
static void ap_dims(const Expr* e, int64_t* dims, size_t num_levels, size_t level) {
    if (level >= num_levels) return;
    if (!ap_row(e)) { dims[level] = 0; return; }
    LLSeq v;
    llq_open(e, &v);
    dims[level] = (int64_t)v.n;
    if (v.n > 0)
        ap_dims(llq_peek(&v, 0), dims, num_levels, level + 1);
    else
        for (size_t k = level + 1; k < num_levels; k++) dims[k] = 0;
    llq_close(&v);
}

/* True for a List[...] expression. */

/* ---- amounts parsing --------------------------------------------------- *
 * Fill lo[]/hi[] for each of *num_levels levels.  Returns false on a bad spec. */
static bool ap_parse_amounts(const Expr* spec, size_t rank,
                             int64_t* lo, int64_t* hi, size_t* num_levels) {
    if (spec->type == EXPR_INTEGER) {
        for (size_t k = 0; k < rank; k++) { lo[k] = hi[k] = spec->data.integer; }
        *num_levels = rank;
        return true;
    }
    if (!ap_is_list(spec)) return false;
    size_t n = spec->data.function.arg_count;
    if (n == 0) return false;
    bool list_of_lists = ap_is_list(spec->data.function.args[0]);
    if (list_of_lists) {
        if (n > AP_MAX_RANK) return false;
        for (size_t k = 0; k < n; k++) {
            const Expr* in = spec->data.function.args[k];
            if (!ap_is_list(in)) return false;
            size_t ic = in->data.function.arg_count;
            const Expr* a = (ic >= 1) ? in->data.function.args[0] : NULL;
            const Expr* b = (ic == 2) ? in->data.function.args[1] : a;
            if (ic < 1 || ic > 2 || a->type != EXPR_INTEGER
                || b->type != EXPR_INTEGER) return false;
            lo[k] = a->data.integer; hi[k] = b->data.integer;
        }
        *num_levels = n;
        return true;
    }
    /* flat integer list: {m} -> {m,m}, {m,n} -> {m,n}, applied to every dim */
    if (n > 2) return false;
    const Expr* a = spec->data.function.args[0];
    const Expr* b = (n == 2) ? spec->data.function.args[1] : a;
    if (a->type != EXPR_INTEGER || b->type != EXPR_INTEGER) return false;
    for (size_t k = 0; k < rank; k++) { lo[k] = a->data.integer; hi[k] = b->data.integer; }
    *num_levels = rank;
    return true;
}

/* ---- constant / cyclic builder ----------------------------------------- *
 * Cyclically index the padding block by the accumulated coordinates (matches
 * PadLeft/PadRight's pr_pad_at). A scalar block ignores coords. */
static Expr* ap_pad_at(const Expr* block, const int64_t* coords, size_t n) {
    if (!block) return expr_new_integer(0);
    if (n == 0 || !ap_is_list(block) || block->data.function.arg_count == 0)
        return expr_copy((Expr*)block);
    int64_t len = (int64_t)block->data.function.arg_count;
    int64_t idx = ((coords[0] % len) + len) % len;
    return ap_pad_at(block->data.function.args[idx], coords + 1, n - 1);
}

/* Build the padded array for a constant/cyclic block. `node` is the original
 * sub-expression at this position (NULL for a pure-padding subtree). Output
 * width at each level is orig_dim[level]+lo+hi (rectangular target). */
static Expr* ap_build_const(const Expr* node, const int64_t* lo, const int64_t* hi,
                            const int64_t* orig_dim, size_t num_levels, size_t level,
                            const Expr* block, int64_t* coords, bool linked) {
    if (level == num_levels)
        return node ? expr_copy((Expr*)node) : ap_pad_at(block, coords, num_levels);

    int64_t N = orig_dim[level] + lo[level] + hi[level];
    if (N < 0) N = 0;
    /* A LinkedList level pads into a LinkedList; a pure-padding level takes
     * the kind of the level above it. */
    LLSeq v;
    llq_open(NULL, &v);
    bool row = node && ap_row(node);
    if (row) { llq_open(node, &v); linked = llq_linked(&v); }
    int64_t L = (int64_t)v.n;
    Expr* head = (row && !linked) ? expr_copy((Expr*)v.head) : expr_new_symbol(SYM_List);
    if (N == 0) {
        llq_close(&v);
        if (linked) { expr_free(head); return ll_seq_of(true, NULL, 0); }
        return expr_new_function(head, NULL, 0);
    }

    if (linked && row && level + 1 == num_levels) {
        /* The last level of a LinkedList: its slots copied straight across
         * (span-wise for the original run), the padding beside them. */
        LinkedList* ol = ll_new();
        int64_t i = 0;
        for (; i < N && i - lo[level] < 0; i++) {
            coords[level] = i - lo[level];
            ll_push_back(ol, ap_pad_at(block, coords, num_levels));
        }
        int64_t a = i - lo[level], b = N - lo[level] < L ? N - lo[level] : L;
        if (a < b) { ll_append_range(ol, v.snap, a, b); i += b - a; }
        for (; i < N; i++) {
            coords[level] = i - lo[level];
            ll_push_back(ol, ap_pad_at(block, coords, num_levels));
        }
        expr_free(head);
        llq_close(&v);
        return expr_new_linkedlist(ol);
    }

    Expr** out = malloc(sizeof(Expr*) * (size_t)N);
    for (int64_t i = 0; i < N; i++) {
        int64_t src = i - lo[level];               /* original element index */
        const Expr* child = (row && src >= 0 && src < L) ? llq_peek(&v, src) : NULL;
        coords[level] = src;
        out[i] = ap_build_const(child, lo, hi, orig_dim, num_levels, level + 1,
                                block, coords, linked);
    }
    Expr* r;
    if (linked) { expr_free(head); r = ll_seq_of(true, out, (size_t)N); }
    else        r = expr_new_function(head, out, (size_t)N);
    free(out);
    llq_close(&v);
    return r;
}

/* ---- value-dependent builder (outer-to-inner) -------------------------- *
 * At each level, extend the fiber of sub-nodes by the scheme, then recurse into
 * every resulting sub-node to pad the deeper levels. Returns NULL on failure
 * (e.g. a too-short axis for a difference scheme). */
/* The fiber [start, ...) of a level, read on demand (pad_scheme_extend_at). */
typedef struct { LLSeq* q; size_t start; } APFiber;
static Expr* ap_fiber_at(void* c, size_t i) {
    APFiber* f = c;
    return llq_peek(f->q, f->start + i);
}

static Expr* ap_pad_valdep(const Expr* node, const int64_t* lo, const int64_t* hi,
                           size_t num_levels, size_t level, PadScheme sc,
                           int64_t order) {
    if (level == num_levels || !ap_row(node)) return expr_copy((Expr*)node);

    LLSeq v;                      /* a LinkedList level: read on demand */
    llq_open(node, &v);
    int64_t L = (int64_t)v.n;
    int64_t l = lo[level], h = hi[level];
    int64_t start = l < 0 ? -l : 0;
    int64_t end = L - (h < 0 ? -h : 0);
    if (end < start) end = start;
    int64_t src_len = end - start;
    int64_t front = l > 0 ? l : 0;
    int64_t back = h > 0 ? h : 0;
    /* A value-dependent scheme cannot derive padding from an empty axis, so
     * such a dimension is left unpadded (matching Wolfram's ArrayPad). */
    if (src_len == 0) { front = 0; back = 0; }

    Expr** cur = NULL;     /* the padded fiber (owned), or NULL: the slice itself */
    size_t cur_len = (size_t)src_len;
    if (front > 0 || back > 0) {
        APFiber fib = { &v, (size_t)start };
        if (!pad_scheme_extend_at(ap_fiber_at, &fib, (size_t)src_len, front, back, NULL, sc, order,
                                  &cur, &cur_len)) {
            llq_close(&v);
            return NULL;
        }
    }

    Expr** out = malloc(sizeof(Expr*) * (cur_len == 0 ? 1 : cur_len));
    bool failed = false;
    size_t built = 0;
    for (size_t i = 0; i < cur_len; i++) {
        Expr* el = cur ? cur[i] : expr_copy(llq_peek(&v, (size_t)start + i));
        out[i] = ap_pad_valdep(el, lo, hi, num_levels, level + 1, sc, order);
        if (!cur) expr_free(el);
        if (!out[i]) { failed = true; break; }
        built++;
    }

    Expr* r = NULL;
    if (!failed)
        r = llq_result(&v, out, cur_len);
    else
        for (size_t i = 0; i < built; i++) expr_free(out[i]);
    free(out);
    if (cur) { for (size_t i = 0; i < cur_len; i++) expr_free(cur[i]); free(cur); }
    llq_close(&v);
    return r;
}

/* Emit `ArrayPad::mindimsize` when a value-difference scheme is asked for
 * positive padding on an axis of length < 2. */
static Expr* ap_mindimsize(const Expr* scheme, const Expr* amounts, const Expr* arr) {
    /* Wolfram prints the scheme name without quotes in this message. */
    char* ss = (scheme && scheme->type == EXPR_STRING)
               ? NULL : expr_to_string((Expr*)scheme);
    const char* sname = ss ? ss : scheme->data.string;
    char* as = expr_to_string((Expr*)amounts);
    char* xs = expr_to_string((Expr*)arr);
    mth_message("ArrayPad", "mindimsize", "With padding %s, the padding amount %s should specify positive padding only in dimensions of length at least 2 in array %s.", sname, as, xs);
    free(ss); free(as); free(xs);
    return NULL;
}

Expr* builtin_array_pad(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    /* Strip a trailing InterpolationOrder option. */
    const Expr* ord_opt = NULL;
    bool ord_given = false;
    const OptEntry ents[1] = { { "InterpolationOrder", &ord_opt, &ord_given } };
    size_t argc = res->data.function.arg_count;
    if (!options_extract(res, "ArrayPad", ents, 1, &argc)) return NULL;
    if (argc < 2 || argc > 3) return builtin_arg_error("ArrayPad", argc, 2, 3);

    int64_t order = 1;                     /* default linear extrapolation */
    if (ord_given && ord_opt) {
        if (is_infinity_sym((Expr*)ord_opt)) order = -1;
        else if (ord_opt->type == EXPR_INTEGER && ord_opt->data.integer >= 0)
            order = ord_opt->data.integer;
        else return NULL;                  /* bad InterpolationOrder value */
    }

    Expr* array = res->data.function.args[0];
    const Expr* amounts = res->data.function.args[1];
    const Expr* padding = (argc == 3) ? res->data.function.args[2] : NULL;

    /* NDArray input: rank-1 constant fast path, else materialize once. */
    Expr* owned_array = NULL;
    if (is_ndarray(array)) {
        Expr* fast = ndstruct_arraypad(res);
        if (fast) return fast;
        owned_array = ndarray_to_nested_list(array);
        if (!owned_array) return NULL;
        array = owned_array;
    }

    if (!ap_row(array)) { if (owned_array) expr_free(owned_array); return NULL; }

    size_t rank = ap_depth(array);
    if (rank == 0 || rank > AP_MAX_RANK) { if (owned_array) expr_free(owned_array); return NULL; }

    int64_t lo[AP_MAX_RANK], hi[AP_MAX_RANK];
    for (size_t k = 0; k < AP_MAX_RANK; k++) { lo[k] = hi[k] = 0; }
    size_t num_levels = 0;
    if (!ap_parse_amounts(amounts, rank, lo, hi, &num_levels)) {
        if (owned_array) expr_free(owned_array);
        return NULL;
    }

    PadScheme sc = pad_scheme_classify(padding);
    if (sc == PS_INVALID) { if (owned_array) expr_free(owned_array); return NULL; }

    Expr* result;
    if (!pad_scheme_value_dependent(sc)) {
        int64_t orig_dim[AP_MAX_RANK];
        for (size_t k = 0; k < AP_MAX_RANK; k++) orig_dim[k] = 0;
        ap_dims(array, orig_dim, num_levels, 0);
        int64_t coords[AP_MAX_RANK];
        Expr* block = padding ? (Expr*)padding : NULL;
        result = ap_build_const(array, lo, hi, orig_dim, num_levels, 0, block, coords, false);
    } else {
        /* mindimsize guard for the difference/reflected schemes. */
        if (pad_scheme_needs_two(sc)) {
            int64_t dm[AP_MAX_RANK];
            for (size_t k = 0; k < AP_MAX_RANK; k++) dm[k] = 0;
            ap_dims(array, dm, num_levels, 0);
            for (size_t k = 0; k < num_levels; k++)
                if ((lo[k] > 0 || hi[k] > 0) && dm[k] < 2) {
                    if (owned_array) expr_free(owned_array);
                    return ap_mindimsize(padding, amounts, array);
                }
        }
        result = ap_pad_valdep(array, lo, hi, num_levels, 0, sc, order);
    }

    if (owned_array) expr_free(owned_array);
    return result;   /* may be NULL -> leave unevaluated */
}
