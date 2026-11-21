#include "list_common.h"
#include "array_reshape.h"
#include "pad_schemes.h"
#include "ndarray.h"     /* is_ndarray */
#include "ndstruct.h"    /* ndstruct_arrayreshape */
#include "linkedlist_seq.h" /* LLSeq: LinkedList input levels */

/* ArrayReshape[list, dims]           arrange the flattened elements of list into
 *                                    a rectangular dims array, dropping extras.
 * ArrayReshape[list, dims, padding]  pad with the given scheme (default 0) when
 *                                    list has too few elements.
 *
 * The input is fully flattened (List levels only), so up to the shared length
 * Flatten[ArrayReshape[list, dims]] == Flatten[list].  A packed/NDArray first
 * argument takes the buffer fast path in ndstruct_arrayreshape. */

#define AR_MAX_RANK 64

/* Parse a dims spec (a non-negative integer, or a List of them) into d[] /
 * *count / *total.  Returns false (leave unevaluated) on a bad spec, an empty
 * List, a negative or non-integer entry, rank overflow, or a total that would
 * overflow int64. */
static bool ar_parse_dims(const Expr* spec, int64_t* d, size_t* count, int64_t* total) {
    size_t k;
    if (spec->type == EXPR_INTEGER) {
        if (spec->data.integer < 0) return false;
        d[0] = spec->data.integer; k = 1;
    } else if ((spec->type == EXPR_FUNCTION
                && spec->data.function.head->type == EXPR_SYMBOL
                && spec->data.function.head->data.symbol.name == SYM_List)
               || spec->type == EXPR_LINKEDLIST) {
        k = (size_t)ll_seq_length(spec);
        if (k == 0 || k > AR_MAX_RANK) return false;
        for (size_t i = 0; i < k; i++) {
            Expr* di = ll_seq_element(spec, (int64_t)i);
            bool ok = di->type == EXPR_INTEGER && di->data.integer >= 0;
            if (ok) d[i] = di->data.integer;
            expr_free(di);
            if (!ok) return false;
        }
    } else {
        return false;
    }
    int64_t t = 1;
    for (size_t i = 0; i < k; i++) {
        if (d[i] != 0 && t > INT64_MAX / d[i]) return false;   /* overflow */
        t *= d[i];
    }
    *count = k; *total = t;
    return true;
}

/* The leaves of the input, in row-major order, as REFERENCES: a List's leaf
 * is borrowed in place; a LinkedList level's leaf is (reader, index), read on
 * demand from that level's snapshot -- never collected into an element
 * array. */
typedef struct { Expr* e; size_t view; size_t idx; } ARLeaf;
typedef struct {
    ARLeaf* leaf; size_t n, cap;
    LLSeq** view; size_t nv, vcap;    /* the LinkedList levels' readers */
} ARLeaves;

static void ar_push(ARLeaves* L, Expr* e, size_t view, size_t idx) {
    if (L->n == L->cap) { L->cap = L->cap ? 2 * L->cap : 16; L->leaf = realloc(L->leaf, sizeof(ARLeaf) * L->cap); }
    L->leaf[L->n].e = e; L->leaf[L->n].view = view; L->leaf[L->n].idx = idx;
    L->n++;
}

static void ar_collect(Expr* e, ARLeaves* L) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) ar_collect(e->data.function.args[i], L);
        return;
    }
    if (e->type == EXPR_LINKEDLIST) {          /* a level, read from its slots */
        if (L->nv == L->vcap) {
            L->vcap = L->vcap ? 2 * L->vcap : 4;
            L->view = realloc(L->view, sizeof(LLSeq*) * L->vcap);
        }
        size_t k = L->nv++;
        LLSeq* v = malloc(sizeof(LLSeq));
        L->view[k] = v;
        llq_open(e, v);
        if (llq_kind(v) != LL_BOXED) {          /* machine numbers: all leaves */
            for (size_t i = 0; i < v->n; i++) ar_push(L, NULL, k, i);
            return;
        }
        for (size_t i = 0; i < v->n; i++) {
            Expr* x = llq_peek(v, i);
            if (x->type == EXPR_LINKEDLIST || head_is(x, SYM_List)) ar_collect(x, L);
            else ar_push(L, NULL, k, i);
        }
        return;
    }
    ar_push(L, e, 0, 0);
}

/* Leaf t, borrowed (a typed slot: a short-lived box). */
static Expr* ar_leaf_at(void* c, size_t t) {
    ARLeaves* L = c;
    ARLeaf* f = &L->leaf[t];
    return f->e ? f->e : llq_peek(L->view[f->view], f->idx);
}

static void ar_leaves_free(ARLeaves* L) {
    for (size_t k = 0; k < L->nv; k++) { llq_close(L->view[k]); free(L->view[k]); }
    free(L->view); free(L->leaf);
}

/* Build the nested rectangular array from leaves 0.. in row-major order
 * (*cursor): an original leaf from its reference (slot to slot into a
 * LinkedList), a padding leaf adopted from pad[] (indices >= L->n). */
static Expr* ar_build(ARLeaves* L, Expr** pad, size_t* cursor, const int64_t* d,
                      size_t dim_count, size_t level, bool linked) {
    if (level == dim_count) {
        size_t t = (*cursor)++;
        return t < L->n ? expr_copy(ar_leaf_at(L, t)) : pad[t];   /* adopt a pad */
    }
    int64_t n = d[level];
    if (linked && level + 1 == dim_count) {    /* a LinkedList row of leaves */
        LinkedList* row = ll_new();
        for (int64_t i = 0; i < n; i++) {
            size_t t = (*cursor)++;
            if (t >= L->n)              ll_push_back(row, pad[t]);
            else if (L->leaf[t].e)      ll_push_back(row, expr_copy(L->leaf[t].e));
            else                        llq_push(row, L->view[L->leaf[t].view], L->leaf[t].idx);
        }
        return expr_new_linkedlist(row);
    }
    if (n == 0) return ll_seq_of(linked, NULL, 0);
    Expr** kids = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t i = 0; i < n; i++)
        kids[i] = ar_build(L, pad, cursor, d, dim_count, level + 1, linked);
    Expr* r = ll_seq_of(linked, kids, (size_t)n);    /* LinkedLists for a LinkedList */
    free(kids);
    return r;
}

Expr* builtin_array_reshape(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3)
        return builtin_arg_error("ArrayReshape", argc, 2, 3);

    Expr* list = res->data.function.args[0];

    /* Packed / NDArray input: buffer fast path (declines to materialize). */
    if (is_ndarray(list)) return ndstruct_arrayreshape(res);

    if (list->type != EXPR_FUNCTION && list->type != EXPR_LINKEDLIST)
        return NULL;   /* not a list: leave alone */
    bool linked = list->type == EXPR_LINKEDLIST;

    int64_t d[AR_MAX_RANK];
    size_t dim_count;
    int64_t total;
    if (!ar_parse_dims(res->data.function.args[1], d, &dim_count, &total))
        return NULL;

    /* The input's leaves, as references. */
    ARLeaves L = { NULL, 0, 0, NULL, 0, 0 };
    if (linked) {
        ar_collect(list, &L);
    } else {
        for (size_t i = 0; i < list->data.function.arg_count; i++)
            ar_collect(list->data.function.args[i], &L);
    }

    /* Exactly `total` leaves: the first `total`, or all of them and a padded
     * tail (the scheme reads the leaves on demand; the tail is adopted). */
    Expr** pad = NULL;
    size_t pad_len = 0;
    if ((int64_t)L.n < total) {
        const Expr* padding = (argc == 3) ? res->data.function.args[2] : NULL;
        PadScheme sc = pad_scheme_classify(padding);
        if (!pad_scheme_extend_at(ar_leaf_at, &L, L.n, 0, total - (int64_t)L.n, padding, sc,
                                  /*order=*/1, &pad, &pad_len)) {
            ar_leaves_free(&L);
            return NULL;
        }
        for (size_t t = 0; t < L.n; t++) expr_free(pad[t]);   /* the originals are read in place */
    }

    size_t cursor = 0;
    Expr* out = ar_build(&L, pad, &cursor, d, dim_count, 0, linked);
    free(pad);
    ar_leaves_free(&L);
    return out;
}
