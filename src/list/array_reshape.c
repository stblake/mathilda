#include "list_common.h"
#include "array_reshape.h"
#include "pad_schemes.h"
#include "ndarray.h"     /* is_ndarray */
#include "ndstruct.h"    /* ndstruct_arrayreshape */
#include "linkedlist_seq.h" /* LLView: LinkedList input levels */

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

/* Collect the leaves of `e` (descending List heads only) as borrowed pointers
 * into a growing buffer. */
/* The LinkedList views a collection borrows its leaves from, kept open until
 * the leaves have been copied out. */
typedef struct { LLView* v; size_t n, cap; } ARViews;

static void ar_collect(Expr* e, Expr*** buf, size_t* n, size_t* cap, ARViews* hold) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            ar_collect(e->data.function.args[i], buf, n, cap, hold);
        return;
    }
    if (e->type == EXPR_LINKEDLIST) {          /* a level, read from its slots */
        if (hold->n == hold->cap) {
            hold->cap = hold->cap ? 2 * hold->cap : 4;
            hold->v = realloc(hold->v, sizeof(LLView) * hold->cap);
        }
        LLView* v = &hold->v[hold->n++];
        ll_view(e, v);
        for (size_t i = 0; i < v->n; i++) ar_collect(v->a[i], buf, n, cap, hold);
        return;
    }
    if (*n == *cap) { *cap *= 2; *buf = realloc(*buf, sizeof(Expr*) * (*cap)); }
    (*buf)[(*n)++] = e;
}

/* Build the nested rectangular array, adopting owned leaves from flat[] in
 * row-major order via *cursor. */
static Expr* ar_build(Expr** flat, size_t* cursor, const int64_t* d,
                      size_t dim_count, size_t level, bool linked) {
    if (level == dim_count) return flat[(*cursor)++];   /* adopt */
    int64_t n = d[level];
    if (n == 0) return ll_seq_of(linked, NULL, 0);
    Expr** kids = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t i = 0; i < n; i++)
        kids[i] = ar_build(flat, cursor, d, dim_count, level + 1, linked);
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

    /* Flatten the input to a borrowed leaf sequence. */
    size_t cap = 16, nleaves = 0;
    Expr** leaves = malloc(sizeof(Expr*) * cap);
    ARViews hold = { NULL, 0, 0 };
    if (linked) {
        ar_collect(list, &leaves, &nleaves, &cap, &hold);
    } else {
        for (size_t i = 0; i < list->data.function.arg_count; i++)
            ar_collect(list->data.function.args[i], &leaves, &nleaves, &cap, &hold);
    }

    /* Produce exactly `total` owned leaves: copy/truncate, or pad the tail. */
    Expr** flat;
    size_t flat_len;
    if ((int64_t)nleaves >= total) {
        flat = malloc(sizeof(Expr*) * (total == 0 ? 1 : (size_t)total));
        for (int64_t i = 0; i < total; i++) flat[i] = expr_copy(leaves[i]);
        flat_len = (size_t)total;
    } else {
        const Expr* padding = (argc == 3) ? res->data.function.args[2] : NULL;
        PadScheme sc = pad_scheme_classify(padding);
        if (!pad_scheme_extend((const Expr**)leaves, nleaves, 0,
                               total - (int64_t)nleaves, padding, sc, /*order=*/1,
                               &flat, &flat_len)) {
            free(leaves);
            for (size_t i = 0; i < hold.n; i++) ll_view_done(&hold.v[i]);
            free(hold.v);
            return NULL;
        }
    }
    free(leaves);
    for (size_t i = 0; i < hold.n; i++) ll_view_done(&hold.v[i]);
    free(hold.v);

    size_t cursor = 0;
    Expr* out = ar_build(flat, &cursor, d, dim_count, 0, linked);
    free(flat);
    return out;
}
