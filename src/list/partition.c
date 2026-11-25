#include "list_common.h"
#include "message.h"   /* mth_message: Partition::ilsmp */
#include "common.h"     /* head_is */
#include "ndarray.h"    /* is_ndarray */
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "partition.h"
#include "assoc.h"      /* is_association */
#include "linkedlist_seq.h" /* LLSeq: a LinkedList partitions from its slots */

/* Level k of a Partition size/offset spec -- element k of a List or a
 * LinkedList of them, or a bare spec at level 0. New reference, or NULL. */
static Expr* partition_spec_at(Expr* spec, size_t k) {
    int64_t len = (spec->type == EXPR_LINKEDLIST || head_is(spec, SYM_List)) ? ll_seq_length(spec) : -1;
    if (len >= 0) return (int64_t)k < len ? ll_seq_element(spec, (int64_t)k) : NULL;
    return k == 0 ? expr_copy(spec) : NULL;
}

/* The deeper levels of one block, re-nested into the block grid.
 *
 * At a level with a deeper spec, the block of rows k = 0 .. n-1 is partitioned
 * row by row: rows[k] = partition_rec(row k) has `g` grid levels (one per
 * deeper spec level) above its block levels. The reference's block layout puts
 * every grid index before every block index -- Partition[m, {n1, n2}] is an
 * N1 x N2 grid of n1 x n2 blocks -- so the row index k must move from the top
 * down past the g grid levels: out[j1]...[jg] = {rows[0][j1..jg], ...,
 * rows[n-1][j1..jg]}. A grid level holds the shortest of the rows there (a
 * ragged input yields only the blocks every row reaches). Borrows rows. */
static Expr* partition_interleave(Expr** rows, size_t n, size_t g, bool linked) {
    if (g == 0) {
        Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t k = 0; k < n; k++) items[k] = expr_copy(rows[k]);
        Expr* r = ll_seq_of(linked, items, n);
        free(items);
        return r;
    }
    int64_t m = -1;
    for (size_t k = 0; k < n; k++) {
        int64_t len = ll_seq_length(rows[k]);
        if (len < 0) len = 0;                          /* an atom reaches no block */
        if (m < 0 || len < m) m = len;
    }
    if (m < 0) m = 0;
    Expr** out = malloc(sizeof(Expr*) * (m ? (size_t)m : 1));
    Expr** sub = malloc(sizeof(Expr*) * (n ? n : 1));
    for (int64_t j = 0; j < m; j++) {
        /* A List row is read in place (borrowed); anything else (a
         * LinkedList, an array) through ll_seq_element, a new reference. */
        for (size_t k = 0; k < n; k++)
            sub[k] = rows[k]->type == EXPR_FUNCTION ? rows[k]->data.function.args[j]
                                                    : ll_seq_element(rows[k], j);
        out[j] = partition_interleave(sub, n, g - 1, linked);
        for (size_t k = 0; k < n; k++)
            if (rows[k]->type != EXPR_FUNCTION) expr_free(sub[k]);
    }
    free(sub);
    Expr* r = ll_seq_of(linked, out, (size_t)m);
    free(out);
    return r;
}

/* How many spec levels lie below level_idx (each one a grid level of the
 * partitioned rows). */
static size_t partition_deeper_levels(Expr* n_spec, size_t level_idx) {
    size_t g = 0;
    for (;;) {
        Expr* e = partition_spec_at(n_spec, level_idx + 1 + g);
        if (!e) return g;
        expr_free(e);
        g++;
    }
}

/* The block of rows q[start .. start + count) partitioned at level `lvl` when
 * that level is the last of the spec and its n (and d) are plain positive
 * integers, built straight into the block layout: column j is the List of
 * row_k[[j d + 1 ;; j d + n]] over the rows k -- what partition_interleave
 * makes of the rows' own partitions, without building them. Only List rows
 * of a List argument; NULL (use the general path) otherwise. */
static Expr* partition_direct_last(LLSeq* q, size_t start, size_t count,
                                   Expr* n_spec, Expr* d_spec, size_t lvl) {
    if (q->snap || count == 0) return NULL;              /* List argument only */
    Expr* nx = partition_spec_at(n_spec, lvl + 1);
    if (nx) { expr_free(nx); return NULL; }              /* not the last level */
    Expr* ne = partition_spec_at(n_spec, lvl);
    if (!ne) return NULL;
    bool n_ok = ne->type == EXPR_INTEGER && ne->data.integer > 0;
    int64_t n = n_ok ? ne->data.integer : 0;
    expr_free(ne);
    if (!n_ok) return NULL;
    int64_t d = n;
    if (d_spec) {
        Expr* de = partition_spec_at(d_spec, lvl);
        if (de) {
            bool d_ok = de->type == EXPR_INTEGER && de->data.integer > 0;
            if (d_ok) d = de->data.integer;
            expr_free(de);
            if (!d_ok) return NULL;
        }
    }
    int64_t m = -1;                                       /* grid columns */
    for (size_t k = 0; k < count; k++) {
        Expr* row = llq_peek(q, start + k);
        if (!(row->type == EXPR_FUNCTION && head_is(row, SYM_List))) return NULL;
        int64_t len = (int64_t)row->data.function.arg_count;
        int64_t c = len >= n ? (len - n) / d + 1 : 0;
        if (m < 0 || c < m) m = c;
    }
    Expr** cols = malloc(sizeof(Expr*) * (m ? (size_t)m : 1));
    Expr** blk = malloc(sizeof(Expr*) * count);
    Expr** seg = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t j = 0; j < m; j++) {
        for (size_t k = 0; k < count; k++) {
            Expr* row = llq_peek(q, start + k);
            for (int64_t t = 0; t < n; t++)
                seg[t] = expr_copy(row->data.function.args[j * d + t]);
            blk[k] = expr_new_function(expr_new_symbol(SYM_List), seg, (size_t)n);
        }
        cols[j] = expr_new_function(expr_new_symbol(SYM_List), blk, count);
    }
    free(seg);
    free(blk);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), cols, (size_t)m);
    free(cols);
    return r;
}

static Expr* partition_rec(Expr* list, Expr* n_spec, Expr* d_spec, size_t level_idx) {
    if (list->type != EXPR_FUNCTION && list->type != EXPR_LINKEDLIST) return expr_copy(list);

    // Get n and d for this level
    int64_t n = -1;
    bool n_upto = false;
    Expr* n_e = partition_spec_at(n_spec, level_idx);
    if (!n_e) return expr_copy(list);

    if (n_e->type == EXPR_INTEGER) {
        n = n_e->data.integer;
    } else if (n_e->type == EXPR_FUNCTION && n_e->data.function.head->data.symbol.name == SYM_UpTo && n_e->data.function.arg_count == 1) {
        if (n_e->data.function.args[0]->type == EXPR_INTEGER) {
            n = n_e->data.function.args[0]->data.integer;
            n_upto = true;
        }
    }
    expr_free(n_e);
    if (n <= 0) return expr_copy(list);

    int64_t d = n;
    if (d_spec) {
        Expr* d_e = partition_spec_at(d_spec, level_idx);
        if (d_e && d_e->type == EXPR_INTEGER) {
            d = d_e->data.integer;
        }
        if (d_e) expr_free(d_e);
    }
    if (d <= 0) return expr_copy(list);

    /* A LinkedList (at any level) partitions into LinkedLists. */
    LLSeq lv;
    llq_open(list, &lv);
    size_t len = lv.n;
    size_t num_sublists = 0;
    if (n_upto) {
        num_sublists = (len > 0) ? (len - 1) / d + 1 : 0;
    } else {
        if (len >= (size_t)n) {
            num_sublists = (len - (size_t)n) / (size_t)d + 1;
        }
    }

    /* With no specification for the next level, the elements are copied as
     * they are: each sublist is a range of the source (span-wise for a
     * LinkedList, so a typed slot is never boxed). */
    Expr* deeper = partition_spec_at(n_spec, level_idx + 1);
    if (deeper) expr_free(deeper);

    Expr** sublists = malloc(sizeof(Expr*) * (num_sublists ? num_sublists : 1));
    for (size_t i = 0; i < num_sublists; i++) {
        size_t start = i * (size_t)d;
        size_t end = start + (size_t)n;
        if (end > len) end = len;
        if (!deeper) { sublists[i] = llq_range(&lv, start, end); continue; }
        
        size_t sub_count = end - start;
        /* The common case -- the next level is the last, with plain integer
         * n2 and d2, over List rows -- built directly: grid column j is
         * {row_k[[j d2 + 1 ;; j d2 + n2]] for each row k}, with no
         * per-row partition to interleave afterwards. */
        {
            Expr* blk = partition_direct_last(&lv, start, sub_count, n_spec, d_spec, level_idx + 1);
            if (blk) { sublists[i] = blk; continue; }
        }
        /* Each row of the block is partitioned at the deeper levels, then
         * the row index is moved below the grid (partition_interleave). */
        Expr** sub_args = malloc(sizeof(Expr*) * (sub_count ? sub_count : 1));
        for (size_t j = 0; j < sub_count; j++) {
            sub_args[j] = partition_rec(llq_peek(&lv, start + j), n_spec, d_spec, level_idx + 1);
        }
        /* A row that is not a list (a spec deeper than the array) has no
         * grid to interleave: such a block keeps its rows as they are. */
        bool rows_are_lists = true;
        for (size_t j = 0; j < sub_count && rows_are_lists; j++)
            rows_are_lists = ll_seq_length(sub_args[j]) >= 0;
        if (!rows_are_lists) {
            sublists[i] = llq_result(&lv, sub_args, sub_count);   /* consumes */
            free(sub_args);
            continue;
        }
        sublists[i] = partition_interleave(sub_args, sub_count,
                                           partition_deeper_levels(n_spec, level_idx),
                                           list->type == EXPR_LINKEDLIST);
        for (size_t j = 0; j < sub_count; j++) expr_free(sub_args[j]);
        free(sub_args);
    }

    Expr* result = llq_result(&lv, sublists, num_sublists);
    free(sublists);
    llq_close(&lv);
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

    LLSeq lv;
    llq_open(list, &lv);
    int64_t len = (int64_t)lv.n;
    if (len == 0) { Expr* r = llq_result(&lv, NULL, 0); llq_close(&lv); return r; }

    /* Sublists i = 0 .. count-1 with s0 + i d <= len - kr + 1. */
    int64_t s0 = 2 - kl;
    int64_t span = len - kr + 1 - s0;
    int64_t count = span < 0 ? 0 : span / d + 1;
    if (count > (int64_t)100000000 / (n > 0 ? n : 1)) { llq_close(&lv); return NULL; }  /* absurd size */

    Expr** subs = malloc(sizeof(Expr*) * (size_t)(count ? count : 1));
    Expr** items = malloc(sizeof(Expr*) * (size_t)n);
    if (!subs || !items) { free(subs); free(items); llq_close(&lv); return NULL; }
    for (int64_t i = 0; i < count; i++) {
        int64_t start = s0 + i * d;
        if (llq_linked(&lv)) {
            /* a LinkedList sublist: the slots copied across, padding beside */
            LinkedList* sub = ll_new();
            for (int64_t j = start; j < start + n; j++) {
                if (j >= 1 && j <= len)       llq_push(sub, &lv, (size_t)(j - 1));
                else if (pad_none)            continue;
                else if (pad_cyclic_self)     llq_push(sub, &lv, (size_t)partition_mod(j - 1, len));
                else if (pad_items)           ll_push_back(sub, expr_copy(pad_items[partition_mod(j - 1, pad_m)]));
                else                          ll_push_back(sub, expr_copy(pad));
            }
            subs[i] = expr_new_linkedlist(sub);
            continue;
        }
        size_t m = 0;
        for (int64_t j = start; j < start + n; j++) {
            Expr* src;
            if (j >= 1 && j <= len)       src = llq_peek(&lv, j - 1);
            else if (pad_none)            continue;
            else if (pad_cyclic_self)     src = llq_peek(&lv, partition_mod(j - 1, len));
            else if (pad_items)           src = pad_items[partition_mod(j - 1, pad_m)];
            else                          src = pad;
            items[m++] = expr_copy(src);
        }
        subs[i] = llq_result(&lv, items, m);
    }
    free(items);
    Expr* out = llq_result(&lv, subs, (size_t)count);
    free(subs);
    llq_close(&lv);
    return out;
}

/* One level of a Partition block size (position 2) or offset (position 3): a
 * positive machine integer; the block size may also be UpTo[k], k > 0. */
static bool partition_level_ok(const Expr* e, bool size) {
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (size && head_is((Expr*)e, SYM_UpTo) && e->data.function.arg_count == 1) {
        const Expr* k = e->data.function.args[0];
        return k->type == EXPR_INTEGER && k->data.integer > 0;
    }
    return false;
}

/* A block size or offset spec: one level, or a nonempty List of levels. */
static bool partition_spec_ok(const Expr* spec, bool size) {
    if (head_is((Expr*)spec, SYM_List)) {
        if (spec->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < spec->data.function.arg_count; i++)
            if (!partition_level_ok(spec->data.function.args[i], size)) return false;
        return true;
    }
    if (spec->type == EXPR_LINKEDLIST) return true;      /* read level by level */
    if (is_ndarray(spec)) {
        /* A packed spec ({2} or {3, 1} packed int64): the same rule per entry. */
        if (spec->data.ndarray.rank != 1 || spec->data.ndarray.dtype != NDT_INT64 ||
            spec->data.ndarray.dims[0] == 0) return false;
        const int64_t* d = (const int64_t*)spec->data.ndarray.data;
        for (int64_t i = 0; i < spec->data.ndarray.dims[0]; i++) if (d[i] <= 0) return false;
        return true;
    }
    return partition_level_ok(spec, size);
}

Expr* builtin_partition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 5) return NULL;

    /* A block size or offset that is not positive (0, -1, UpTo[0]), not an
     * integer (x, 1.5) or not machine-sized (10^30) has no partition: the call
     * stays unevaluated with Partition::ilsmp, never answered with the input
     * (audit Partition#2/#3). */
    for (size_t pos = 1; pos <= 2 && pos < res->data.function.arg_count; pos++) {
        if (!partition_spec_ok(res->data.function.args[pos], pos == 1)) {
            char* cs = expr_to_string(res);
            mth_message("Partition", "ilsmp", "Single or list of positive machine-size integers "
                        "expected at position %zu of %s.", pos + 1, cs ? cs : "");
            free(cs);
            return NULL;
        }
    }

    if (res->data.function.arg_count >= 4) {
        Expr* list = res->data.function.args[0];
        if (is_ndarray(list)) return ndstruct_delist_repack(res, list);
        if ((list->type != EXPR_FUNCTION && list->type != EXPR_LINKEDLIST) || is_association(list))
            return NULL;
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

    /* Partition[x, 2], Partition["abcd", 2]: Partition::normal (audit
     * Partition#4); any other payload atom is left as it was. */
    if (expr_is_plain_atom(list)) return builtin_normal_error("Partition", res, 1);
    if (list->type != EXPR_FUNCTION && list->type != EXPR_LINKEDLIST) return expr_copy(list);

    /* The reference does not partition an association: Partition[<|...|>, 2]
     * stays unevaluated (it used to come back as the input, truncated). */
    if (is_association(list)) return NULL;

    return partition_rec(list, n_spec, d_spec, 0);
}
