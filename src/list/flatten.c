#include "list_common.h"
#include "flatten.h"
#include "ndarray.h"
#include "ndstruct.h"
#include "message.h"
#include "linkedlist_seq.h" /* LLSeq: LinkedList levels flatten like List levels */

static void flatten_rec(Expr* e, const char* h, int64_t level, Expr*** results, size_t* count, size_t* cap) {
    /* A visible NDArray nested inside an ordinary List. `head_is` is false for
     * it -- an EXPR_NDARRAY is not an EXPR_FUNCTION headed List -- so it was
     * collected as one opaque element and `Flatten[{ndA, ndB}]` answered with
     * the list unchanged where `Flatten[{{1., 2.}, {3., 4.}}]` gives
     * `{1., 2., 3., 4.}`. It IS a list of values by every other measure
     * (`ArrayQ`, `Dimensions`, `Length` all say so), so descending is the
     * answer that agrees with the other surface.
     *
     * This is why the sweep read `Dot`/`Inverse`/`LinearSolve` on a 6x6 as
     * ND-UNSUPPORTED: all three answered correctly and it was the probe's own
     * checksum, `N[Total[Flatten[{r}]]]` over a List of 200 result arrays, that
     * could not reduce to a number.
     *
     * The PACKED form cannot reach here -- the no-nesting invariant means a
     * packed node never sits inside a plain List (pack.c) -- so this is a
     * visible-surface repair only, and Flatten[ndarray] at the top level still
     * takes ndstruct_flatten's reshape below. */
    if (level != 0 && is_ndarray(e) && h == SYM_List) {
        Expr* rows = ndarray_to_nested_list(e);
        for (size_t i = 0; i < rows->data.function.arg_count; i++)
            flatten_rec(rows->data.function.args[i], h,
                        level == -1 ? -1 : level - 1, results, count, cap);
        expr_free(rows);
        return;
    }
    if (level != 0 && e->type == EXPR_LINKEDLIST && h == SYM_List) {
        /* a LinkedList is a List level, read from its slots */
        LLSeq v;
        llq_open(e, &v);
        for (size_t i = 0; i < v.n; i++)
            flatten_rec(llq_peek(&v, i), h, level == -1 ? -1 : level - 1, results, count, cap);
        llq_close(&v);
        return;
    }
    if (level != 0 && head_is(e, intern_symbol(h))) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_rec(e->data.function.args[i], h, level == -1 ? -1 : level - 1, results, count, cap);
        }
    } else {
        if (*count == *cap) {
            *cap *= 2;
            *results = realloc(*results, sizeof(Expr*) * (*cap));
        }
        (*results)[(*count)++] = expr_copy(e);
    }
}

/* flatten_rec for a LinkedList answer: each kept element pushed into `out`
 * as it is reached -- a typed LinkedList level appended span-wise, a boxed one
 * read element by element from a snapshot -- so no element array is built. */
static void flatten_ll_rec(LinkedList* out, Expr* e, const char* h, int64_t level) {
    int64_t next = level == -1 ? -1 : level - 1;
    if (level != 0 && e->type == EXPR_LINKEDLIST && h == SYM_List) {
        const LinkedList* ll = e->data.llist;
        if (ll_kind(ll) != LL_BOXED) {             /* machine numbers: atoms */
            ll_append_range(out, ll, 0, ll_length(ll));
            return;
        }
        LLSeq v;
        llq_open(e, &v);
        for (size_t i = 0; i < v.n; i++) flatten_ll_rec(out, llq_peek(&v, i), h, next);
        llq_close(&v);
        return;
    }
    if (level != 0 && (head_is(e, intern_symbol(h)) || (is_ndarray(e) && h == SYM_List))) {
        size_t cap = 16, count = 0;
        Expr** tmp = malloc(sizeof(Expr*) * cap);
        flatten_rec(e, h, level, &tmp, &count, &cap);     /* a List level */
        for (size_t i = 0; i < count; i++) ll_push_back(out, tmp[i]);
        free(tmp);
        return;
    }
    ll_push_back(out, expr_copy(e));
}

bool flatten_level_is_infinity(const Expr* lv) {
    if (lv->type == EXPR_SYMBOL) return lv->data.symbol.name == SYM_Infinity;
    return head_is(lv, SYM_DirectedInfinity) && lv->data.function.arg_count == 1 &&
           lv->data.function.args[0]->type == EXPR_INTEGER &&
           lv->data.function.args[0]->data.integer == 1;
}

void flatten_flev_message(const Expr* call) {
    char* ls = expr_to_string_outputform(call->data.function.args[1]);
    char* cs = expr_to_string_outputform((Expr*)call);
    mth_message("Flatten", "flev",
                "The level argument %s in position 2 of %s should be a non-negative integer or "
                "Infinity giving the levels to flatten through or a list of lists of levels to "
                "flatten together.", ls ? ls : "", cs ? cs : "");
    free(ls);
    free(cs);
}

/* ---- Flatten[list, {{s11, s12, ...}, {s21, ...}, ...}] ------------------
 *
 * Level i of the result combines the levels s_i1, s_i2, ... of list; levels up
 * to the deepest one named that no group names are appended as groups of
 * their own, in order.  (audit Flatten#2: the form stayed unevaluated.)
 *
 * Every part at depth D = max s_ij is a leaf, found at a position
 * (p_1, ..., p_D).  Its key is that position read group by group --
 * (p_s11, p_s12, ..., p_s21, ...) -- and the result is the leaves sorted by
 * key, grouped on each group's part of it in turn: level i holds one part per
 * distinct value of group i's indices.  For a full array this is the
 * transpose / reshape; for a ragged one it is the reference's ragged
 * transpose, Flatten[{{a, b, c}, {d, e}, {f}}, {{2}, {1}}] being
 * {{a, d, f}, {b, e}, {c}}.  The new levels take the head of list.  A part
 * above depth D that is not a list (under that head) leaves the call
 * unevaluated. */

typedef struct {
    Expr**   leaf;   /* borrowed */
    int64_t* pos;    /* L x D positions, row-major */
    size_t   n, cap;
    int      D;
} FlatLeaves;

static bool fg_collect(FlatLeaves* fl, Expr* e, const char* h, int depth, int64_t* path) {
    if (depth == fl->D) {
        if (fl->n == fl->cap) {
            size_t nc = fl->cap ? 2 * fl->cap : 16;
            Expr** nl = realloc(fl->leaf, nc * sizeof(Expr*));
            if (!nl) return false;
            fl->leaf = nl;
            int64_t* np = realloc(fl->pos, nc * (size_t)fl->D * sizeof(int64_t));
            if (!np) return false;
            fl->pos = np;
            fl->cap = nc;
        }
        fl->leaf[fl->n] = e;
        memcpy(fl->pos + fl->n * (size_t)fl->D, path, (size_t)fl->D * sizeof(int64_t));
        fl->n++;
        return true;
    }
    if (!head_is(e, intern_symbol(h))) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        path[depth] = (int64_t)i;
        if (!fg_collect(fl, e->data.function.args[i], h, depth + 1, path)) return false;
    }
    return true;
}

/* Compare leaves a and b by key: their positions read in the order `ord`. */
static int fg_cmp(const FlatLeaves* fl, const int* ord, size_t a, size_t b) {
    const int64_t* pa = fl->pos + a * (size_t)fl->D;
    const int64_t* pb = fl->pos + b * (size_t)fl->D;
    for (int k = 0; k < fl->D; k++) {
        int64_t x = pa[ord[k]], y = pb[ord[k]];
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}

/* Bottom-up merge sort of idx[0..n) by key (stable; no qsort context in C99). */
static bool fg_sort(const FlatLeaves* fl, const int* ord, size_t* idx, size_t n) {
    size_t* tmp = malloc((n ? n : 1) * sizeof(size_t));
    if (!tmp) return false;
    for (size_t w = 1; w < n; w *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * w) {
            size_t mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
            size_t i = lo, j = mid, k = lo;
            while (i < mid && j < hi) tmp[k++] = fg_cmp(fl, ord, idx[j], idx[i]) < 0 ? idx[j++] : idx[i++];
            while (i < mid) tmp[k++] = idx[i++];
            while (j < hi) tmp[k++] = idx[j++];
        }
        memcpy(idx, tmp, n * sizeof(size_t));
    }
    free(tmp);
    return true;
}

/* The part for sorted leaves idx[a..b) at result level g (0-based): h[...] of
 * one entry per distinct value of group g's indices. */
static Expr* fg_build(const FlatLeaves* fl, const int* ord, const int* gstart, int G,
                      const size_t* idx, size_t a, size_t b, int g, const char* h) {
    size_t cap = 8, cnt = 0;
    Expr** items = malloc(cap * sizeof(Expr*));
    int k0 = gstart[g], k1 = gstart[g + 1];
    size_t i = a;
    while (i < b) {
        size_t j = i + 1;
        const int64_t* pi = fl->pos + idx[i] * (size_t)fl->D;
        while (j < b) {
            const int64_t* pj = fl->pos + idx[j] * (size_t)fl->D;
            int same = 1;
            for (int k = k0; k < k1 && same; k++) same = pi[ord[k]] == pj[ord[k]];
            if (!same) break;
            j++;
        }
        if (cnt == cap) { cap *= 2; items = realloc(items, cap * sizeof(Expr*)); }
        items[cnt++] = g + 1 == G ? expr_copy(fl->leaf[idx[i]])
                                  : fg_build(fl, ord, gstart, G, idx, i, j, g + 1, h);
        i = j;
    }
    Expr* r = expr_new_function(expr_new_symbol(h), items, cnt);
    free(items);
    return r;
}

/* 1: `spec` is a list of non-empty lists of positive machine integers, with
 * the key order and group starts filled (caller frees); 0: it is not. */
static int fg_parse(const Expr* spec, int** ord_out, int** gstart_out, int* G_out, int* D_out) {
    if (!head_is((Expr*)spec, SYM_List) || spec->data.function.arg_count == 0) return 0;
    size_t ng = spec->data.function.arg_count;
    int64_t D = 0;
    size_t total = 0;
    for (size_t g = 0; g < ng; g++) {
        const Expr* grp = spec->data.function.args[g];
        if (!head_is((Expr*)grp, SYM_List) || grp->data.function.arg_count == 0) return 0;
        for (size_t k = 0; k < grp->data.function.arg_count; k++) {
            const Expr* s = grp->data.function.args[k];
            if (s->type != EXPR_INTEGER || s->data.integer < 1 || s->data.integer > 4096) return 0;
            if (s->data.integer > D) D = s->data.integer;
            total++;
        }
    }
    char* seen = calloc((size_t)D + 1, 1);
    int* ord = malloc(((size_t)D + total) * sizeof(int));
    int* gstart = malloc(((size_t)D + ng + 1) * sizeof(int));
    if (!seen || !ord || !gstart) { free(seen); free(ord); free(gstart); return 0; }
    int no = 0, G = 0;
    for (size_t g = 0; g < ng; g++) {
        const Expr* grp = spec->data.function.args[g];
        gstart[G++] = no;
        for (size_t k = 0; k < grp->data.function.arg_count; k++) {
            int64_t s = grp->data.function.args[k]->data.integer;
            if (seen[s]) { free(seen); free(ord); free(gstart); return 0; }   /* repeated */
            seen[s] = 1;
            ord[no++] = (int)(s - 1);
        }
    }
    for (int64_t s = 1; s <= D; s++)
        if (!seen[s]) { gstart[G++] = no; ord[no++] = (int)(s - 1); }
    gstart[G] = no;
    free(seen);
    *ord_out = ord; *gstart_out = gstart; *G_out = G; *D_out = (int)D;
    return 1;
}

/* Flatten[list, groups] for a parsed `groups`; NULL leaves the call alone. */
static Expr* flatten_groups(Expr* list, const int* ord, const int* gstart, int G, int D) {
    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = list->data.function.head->data.symbol.name;
    FlatLeaves fl = { NULL, NULL, 0, 0, D };
    int64_t* path = malloc((size_t)D * sizeof(int64_t));
    Expr* out = NULL;
    if (path && fg_collect(&fl, list, h, 0, path)) {
        size_t* idx = malloc((fl.n ? fl.n : 1) * sizeof(size_t));
        if (idx) {
            for (size_t i = 0; i < fl.n; i++) idx[i] = i;
            if (fg_sort(&fl, ord, idx, fl.n))
                out = fg_build(&fl, ord, gstart, G, idx, 0, fl.n, 0, h);
            free(idx);
        }
    }
    free(path); free(fl.leaf); free(fl.pos);
    return out;
}

Expr* builtin_flatten(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    Expr* list = res->data.function.args[0];

    /* The level: a non-negative integer, Infinity (every level; a positive
     * BigInt is as deep), or a list of lists of levels to flatten together.
     * Anything else -- a negative integer, a symbol, a Real -- is
     * Flatten::flev and the call stays as it is (audit Flatten#4: a symbol
     * or Real was silent), checked before any fast path, so a buffer argument
     * gets the same answer. */
    int64_t n = -1; // -1 means infinity
    if (res->data.function.arg_count >= 2) {
        Expr* lv = res->data.function.args[1];
        if (flatten_level_is_infinity(lv) ||
            (lv->type == EXPR_BIGINT && mpz_sgn(lv->data.bigint) > 0)) {
            n = -1;
        } else if (head_is(lv, SYM_List)) {
            int *ord, *gstart, G, D;
            if (!fg_parse(lv, &ord, &gstart, &G, &D)) {
                if (is_ndarray(list)) return ndarray_delist_and_reeval(res);
                flatten_flev_message(res);
                return NULL;
            }
            Expr* r = NULL;
            if (is_ndarray(list)) {
                free(ord); free(gstart);
                return ndarray_delist_and_reeval(res);
            }
            if (res->data.function.arg_count == 2) r = flatten_groups(list, ord, gstart, G, D);
            free(ord); free(gstart);
            return r;
        } else {
            if (lv->type != EXPR_INTEGER) {
                if (is_ndarray(list)) return ndarray_delist_and_reeval(res);
                flatten_flev_message(res);
                return NULL;
            }
            n = lv->data.integer;
            if (n < 0) {
                /* A buffer is reported on its List form, once: declining with
                 * the buffer in place would let the evaluator's post-gate
                 * materialise it and ask (and report) a second time. */
                if (is_ndarray(list)) return ndarray_delist_and_reeval(res);
                flatten_flev_message(res);
                return NULL;
            }
        }
    }

    /* Flatten[ndarray]: reshape the flat buffer to rank 1 (see ndstruct.c). */
    if (is_ndarray(list)) return ndstruct_flatten(res);

    /* a Rational/Complex number is an atom: ::normal, unevaluated. (Other atoms
     * keep the identity result, which the reference also rejects with ::normal:
     * internal code relies on Flatten[atom] === atom, so that stays.) */
    if (expr_is_number_atom(list)) return builtin_normal_error("Flatten", res, 1);
    if (list->type != EXPR_FUNCTION && list->type != EXPR_LINKEDLIST) return expr_copy(list);

    /* The head flattened through is, by default, the head of expr itself:
     * Flatten[f[a, f[b, c], f[f[d]]]] is f[a, b, c, d]. (A LinkedList is a
     * List level.) A compound head (g[x][...]) has no symbol to compare by and
     * keeps the List default. */
    const char* h = SYM_List;
    if (list->type == EXPR_FUNCTION && list->data.function.head->type == EXPR_SYMBOL)
        h = list->data.function.head->data.symbol.name;
    if (res->data.function.arg_count == 3) {
        if (res->data.function.args[2]->type != EXPR_SYMBOL) return NULL;
        h = res->data.function.args[2]->data.symbol.name;
    }

    /* A typed LinkedList has no levels below its own: a copy of its slots. */
    if (list->type == EXPR_LINKEDLIST && ll_kind(list->data.llist) != LL_BOXED)
        return expr_new_linkedlist(ll_copy(list->data.llist));

    if (list->type == EXPR_LINKEDLIST) {          /* a LinkedList, pushed as reached */
        LinkedList* out = ll_new();
        LLSeq v;
        llq_open(list, &v);
        for (size_t i = 0; i < v.n; i++) flatten_ll_rec(out, llq_peek(&v, i), h, n);
        llq_close(&v);
        return expr_new_linkedlist(out);
    }

    size_t cap = 16;
    size_t count = 0;
    Expr** results = malloc(sizeof(Expr*) * cap);

    // Initial call: we flatten children of the head if they also have head h.
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        flatten_rec(list->data.function.args[i], h, n, &results, &count, &cap);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), results, count);
    free(results);
    return result;
}
