#include "list_common.h"
#include "ndarray.h"    /* is_ndarray */
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "partition.h"
#include "assoc.h"      /* is_association */

static Expr* partition_rec(Expr* list, Expr* n_spec, Expr* d_spec, size_t level_idx) {
    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    // Get n and d for this level
    int64_t n = -1;
    bool n_upto = false;
    Expr* n_e = (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol.name == SYM_List) ? 
                (level_idx < n_spec->data.function.arg_count ? n_spec->data.function.args[level_idx] : NULL) : 
                (level_idx == 0 ? n_spec : NULL);
    
    if (!n_e) return expr_copy(list);

    if (n_e->type == EXPR_INTEGER) {
        n = n_e->data.integer;
    } else if (n_e->type == EXPR_FUNCTION && n_e->data.function.head->data.symbol.name == SYM_UpTo && n_e->data.function.arg_count == 1) {
        if (n_e->data.function.args[0]->type == EXPR_INTEGER) {
            n = n_e->data.function.args[0]->data.integer;
            n_upto = true;
        }
    }
    if (n <= 0) return expr_copy(list);

    int64_t d = n;
    if (d_spec) {
        Expr* d_e = (d_spec->type == EXPR_FUNCTION && d_spec->data.function.head->data.symbol.name == SYM_List) ? 
                    (level_idx < d_spec->data.function.arg_count ? d_spec->data.function.args[level_idx] : NULL) : 
                    (level_idx == 0 ? d_spec : NULL);
        if (d_e && d_e->type == EXPR_INTEGER) {
            d = d_e->data.integer;
        }
    }
    if (d <= 0) return expr_copy(list);

    size_t len = list->data.function.arg_count;
    size_t num_sublists = 0;
    if (n_upto) {
        num_sublists = (len > 0) ? (len - 1) / d + 1 : 0;
    } else {
        if (len >= (size_t)n) {
            num_sublists = (len - (size_t)n) / (size_t)d + 1;
        }
    }

    Expr** sublists = malloc(sizeof(Expr*) * num_sublists);
    for (size_t i = 0; i < num_sublists; i++) {
        size_t start = i * (size_t)d;
        size_t end = start + (size_t)n;
        if (end > len) end = len;
        
        size_t sub_count = end - start;
        Expr** sub_args = malloc(sizeof(Expr*) * sub_count);
        for (size_t j = 0; j < sub_count; j++) {
            sub_args[j] = partition_rec(list->data.function.args[start + j], n_spec, d_spec, level_idx + 1);
        }
        sublists[i] = expr_new_function(expr_copy(list->data.function.head), sub_args, sub_count);
        free(sub_args);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), sublists, num_sublists);
    free(sublists);
    return result;
}

/* ---------------------------------------------------------------------------
 * Overhang forms: Partition[list, n, d, {kL, kR}] and
 * Partition[list, n, d, {kL, kR}, padding]  (a bare k means {k, k}).
 *
 * The first element of `list` sits at position kL of the first sublist and the
 * last element at or after position kR of the last one; a negative k counts
 * from the end of a sublist (k -> n + 1 + k). With 1-based list indices the
 * first sublist therefore starts at s0 = 2 - kL, sublist i at s0 + i d, and the
 * sublists kept are exactly those in which the last element L still sits at a
 * position >= kR: L - (s0 + i d) + 1 >= kR. An index outside 1..L is padding:
 *   - no padding argument: cyclic, from the list itself (index j -> list[j mod L]);
 *   - padding {x1, ..., xm}: cyclic over the xi on the same absolute index
 *     (j -> x[(j - 1) mod m + 1]), so Partition[{a,...,g}, 3, 3, 1, {x, y, z}]
 *     ends {g, y, z};
 *   - padding {}: no padding -- the element is dropped, so edge sublists are
 *     shorter;
 *   - any other padding expression x: x itself.
 *
 * Without this the four-argument cyclic form came back unevaluated, and
 * Function[{a, b}, ...] @@@ Partition[pts, 2, 1, 1] then applied the function to
 * the unevaluated Partition's own argument list -- a closed-polygon edge walk
 * silently became a computation over the wrong objects, which in a Nest of such
 * steps grew the expression geometrically.
 *
 * Only the one-level form (n and d integers, or one-element lists of them) is
 * handled; a multi-level spec or UpTo with an overhang returns NULL. */
static bool partition_int_spec(const Expr* e, int64_t* out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_List &&
        e->data.function.arg_count == 1)
        e = e->data.function.args[0];
    if (e->type != EXPR_INTEGER) return false;
    *out = e->data.integer;
    return true;
}

/* Floor modulus for a possibly negative index. */
static int64_t partition_mod(int64_t a, int64_t m) {
    int64_t r = a % m;
    return r < 0 ? r + m : r;
}

static Expr* partition_overhang(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr* list = res->data.function.args[0];
    int64_t n, d;
    if (!partition_int_spec(res->data.function.args[1], &n) || n <= 0) return NULL;
    if (!partition_int_spec(res->data.function.args[2], &d) || d <= 0) return NULL;

    /* Alignment {kL, kR} or a bare k. */
    Expr* ks = res->data.function.args[3];
    int64_t kl, kr;
    if (ks->type == EXPR_INTEGER) {
        kl = kr = ks->data.integer;
    } else if (ks->type == EXPR_FUNCTION && ks->data.function.head->type == EXPR_SYMBOL &&
               ks->data.function.head->data.symbol.name == SYM_List &&
               ks->data.function.arg_count == 2 &&
               ks->data.function.args[0]->type == EXPR_INTEGER &&
               ks->data.function.args[1]->type == EXPR_INTEGER) {
        kl = ks->data.function.args[0]->data.integer;
        kr = ks->data.function.args[1]->data.integer;
    } else {
        return NULL;
    }
    if (kl < 0) kl = n + 1 + kl;
    if (kr < 0) kr = n + 1 + kr;
    if (kl < 1 || kr < 1) return NULL;        /* k = 0 / beyond -n: not handled */

    /* Padding mode. */
    Expr* pad = (argc == 5) ? res->data.function.args[4] : NULL;
    Expr** pad_items = NULL;
    int64_t pad_m = 0;
    bool pad_cyclic_self = (pad == NULL);
    bool pad_none = false;
    if (pad && pad->type == EXPR_FUNCTION && pad->data.function.head->type == EXPR_SYMBOL &&
        pad->data.function.head->data.symbol.name == SYM_List) {
        pad_m = (int64_t)pad->data.function.arg_count;
        pad_items = pad->data.function.args;
        if (pad_m == 0) pad_none = true;
    }

    int64_t len = (int64_t)list->data.function.arg_count;
    Expr* head = list->data.function.head;
    if (len == 0) return expr_new_function(expr_copy(head), NULL, 0);

    /* Sublists i = 0 .. count-1 with s0 + i d <= len - kr + 1. */
    int64_t s0 = 2 - kl;
    int64_t span = len - kr + 1 - s0;
    int64_t count = span < 0 ? 0 : span / d + 1;
    if (count > (int64_t)100000000 / (n > 0 ? n : 1)) return NULL;  /* absurd size */

    Expr** subs = malloc(sizeof(Expr*) * (size_t)(count ? count : 1));
    Expr** items = malloc(sizeof(Expr*) * (size_t)n);
    if (!subs || !items) { free(subs); free(items); return NULL; }
    for (int64_t i = 0; i < count; i++) {
        int64_t start = s0 + i * d;
        size_t m = 0;
        for (int64_t j = start; j < start + n; j++) {
            Expr* src;
            if (j >= 1 && j <= len)       src = list->data.function.args[j - 1];
            else if (pad_none)            continue;
            else if (pad_cyclic_self)     src = list->data.function.args[partition_mod(j - 1, len)];
            else if (pad_items)           src = pad_items[partition_mod(j - 1, pad_m)];
            else                          src = pad;
            items[m++] = expr_copy(src);
        }
        subs[i] = expr_new_function(expr_copy(head), items, m);
    }
    free(items);
    Expr* out = expr_new_function(expr_copy(head), subs, (size_t)count);
    free(subs);
    return out;
}

Expr* builtin_partition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 5) return NULL;

    if (res->data.function.arg_count >= 4) {
        Expr* list = res->data.function.args[0];
        if (is_ndarray(list)) return ndstruct_delist_repack(res, list);
        if (list->type != EXPR_FUNCTION || is_association(list)) return NULL;
        return partition_overhang(res);
    }
    
    Expr* list = res->data.function.args[0];
    Expr* n_spec = res->data.function.args[1];
    Expr* d_spec = (res->data.function.arg_count == 3) ? res->data.function.args[2] : NULL;

    /* Native buffer path first; ndstruct_delist_repack below stays the fallback
     * for every form it declines. See ndstruct.h. */
    if (is_ndarray(list)) {
        Expr* fast = ndstruct_partition(res);
        if (fast) return fast;
        return ndstruct_delist_repack(res, list);
    }

    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    /* Mathematica 15 does not partition an association: Partition[<|...|>, 2]
     * stays unevaluated (it used to come back as the input, truncated). */
    if (is_association(list)) return NULL;

    return partition_rec(list, n_spec, d_spec, 0);
}
