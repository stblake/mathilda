/* ---------------------------------------------------------------------------
 * assoc_packed.c — machine-buffer paths of the Association builtins.
 *
 * See assoc_packed.h for the contract and for why every nested output is boxed.
 * The shared substrate is nd_group_words (src/ndreduce.c): one pass that keys
 * each element by its machine word -- direct-indexed over a narrow int64 range,
 * hashed otherwise -- and hands back the distinct words in first-appearance
 * order plus every element's group id. That is PositionIndex, GroupBy, GatherBy,
 * CountsBy and the key de-duplication of AssociationThread / AssociationMap,
 * each a different read of the same three arrays.
 * -------------------------------------------------------------------------- */

#include "assoc_packed.h"
#include "assoc.h"
#include "ndarray.h"
#include "ndreduce.h"
#include "pack.h"
#include "sym_names.h"
#include "eval.h"
#include "compile/autocompile.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ---------------------------------------------------------------- helpers */

Expr* assoc_list_view(Expr* a, Expr** owned) {
    *owned = NULL;
    if (a && is_ndarray(a)) {
        *owned = ndarray_to_nested_list(a);
        if (*owned) return *owned;
    }
    return a;
}

/* A rank-1 buffer whose words nd_group_words can key. */
static bool keyable_vector(const Expr* a) {
    if (!a || !is_ndarray(a) || a->data.ndarray.rank != 1) return false;
    NDType dt = a->data.ndarray.dtype;
    return (dt == NDT_INT64 || dt == NDT_FLOAT64) && a->data.ndarray.dims[0] > 0;
}

/* The Expr a distinct machine word stands for -- the head ndarray_buffer_
 * element_to_expr would give that element, so a key built here is SameQ to the
 * one the List path builds. */
static Expr* word_to_expr(uint64_t w, NDType dt) {
    if (dt == NDT_INT64) return expr_new_integer((int64_t)w);
    if (dt == NDT_BOOL)  return expr_new_symbol(w ? SYM_True : SYM_False);
    double x;
    memcpy(&x, &w, sizeof x);
    return expr_new_real(x);
}

/* Mark a node this module built as ALREADY at its fixed point under the current
 * evaluation clock, so the evaluator's entry / in-loop short-circuit returns it
 * as-is instead of re-walking it. Only ever applied to trees of machine numbers
 * (and True/False) under List / Rule / Association with unique keys -- which is
 * what those constructors leave unchanged -- or to values that evaluate() itself
 * has just returned. Without it, re-evaluating a 10^6-element PositionIndex
 * result cost as much as building it. The stamp is benign metadata (see
 * last_evaluated_at in expr.h); any later symbol-table mutation bumps the clock
 * and simply makes the node evaluate normally again. */
static Expr* stamp_fixed(Expr* e) {
    if (e) e->last_evaluated_at = eval_clock_get();
    return e;
}

static Expr* make_rule2(Expr* k, Expr* v) {
    Expr* args[2] = { k, v };
    return expr_new_function(expr_new_symbol(SYM_Rule), args, 2);
}

static Expr* make_assoc(Expr** rules, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_Association), rules, n);
}

/* Element i of a buffer, boxed. */
static Expr* elem_at(const Expr* arr, size_t i) {
    return ndarray_buffer_element_to_expr(arr->data.ndarray.data, i, arr->data.ndarray.dtype);
}

static bool group_permutation(const NDWordGroups* g, size_t n, size_t** perm_out,
                              size_t** off_out);

/* The elements of each group, boxed: lists[u] is a List with g->cnts[u] slots
 * holding the elements of group u in input order. `arr` supplies the elements
 * (it need not be the buffer that was grouped). Caller owns every lists[u].
 * Elements are allocated group by group -- the order the result is later freed
 * in, which the allocator handles far better than a strided free. */
static Expr** build_group_lists(const Expr* arr, const NDWordGroups* g) {
    size_t nu = g->nuniq;
    size_t n = (size_t)arr->data.ndarray.dims[0];
    Expr** lists = malloc(sizeof(Expr*) * nu);
    size_t *perm = NULL, *off = NULL;
    if (!lists || !group_permutation(g, n, &perm, &off)) { free(lists); return NULL; }
    const void* data = arr->data.ndarray.data;
    NDType dt = arr->data.ndarray.dtype;
    for (size_t u = 0; u < nu; u++) {
        size_t cnt = off[u + 1] - off[u];
        const size_t* pos = perm + off[u];
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), NULL, cnt);
        Expr** slot = list->data.function.args;
        if (dt == NDT_INT64) {
            const int64_t* iv = (const int64_t*)data;
            for (size_t k = 0; k < cnt; k++) slot[k] = expr_new_integer(iv[pos[k]]);
        } else if (dt == NDT_FLOAT64) {
            const double* dv = (const double*)data;
            for (size_t k = 0; k < cnt; k++) slot[k] = expr_new_real(dv[pos[k]]);
        } else {
            for (size_t k = 0; k < cnt; k++) slot[k] = elem_at(arr, pos[k]);
        }
        lists[u] = stamp_fixed(list);
    }
    free(perm); free(off);
    return lists;
}

/* The input positions of every group, contiguously: a stable counting sort of
 * 0..n-1 by group id. Group u's positions are perm[off[u] .. off[u+1]). */
static bool group_permutation(const NDWordGroups* g, size_t n, size_t** perm_out,
                              size_t** off_out) {
    size_t nu = g->nuniq;
    size_t* off = malloc(sizeof(size_t) * (nu + 1));
    size_t* cur = malloc(sizeof(size_t) * (nu ? nu : 1));
    size_t* perm = malloc(sizeof(size_t) * n);
    if (!off || !cur || !perm) { free(off); free(cur); free(perm); return false; }
    off[0] = 0;
    for (size_t u = 0; u < nu; u++) { off[u + 1] = off[u] + (size_t)g->cnts[u]; cur[u] = off[u]; }
    for (size_t i = 0; i < n; i++) perm[cur[g->gid[i]]++] = i;
    free(cur);
    *perm_out = perm; *off_out = off;
    return true;
}

/* One group's elements as a top-level list for a reducer call: PACKED when it
 * clears the packing threshold (it is handed straight to the reducer, never
 * nested), else boxed. */
static Expr* build_group_arg(const Expr* arr, const size_t* pos, size_t cnt) {
    NDType dt = arr->data.ndarray.dtype;
    int64_t len = (int64_t)cnt;
    void* buf = NULL;
    Expr* out = ndbuild_open(1, &len, dt, &buf);   /* NULL below threshold / packing off */
    if (out) {
        size_t esz = ndt_elem_size(dt);
        const char* src = (const char*)arr->data.ndarray.data;
        for (size_t k = 0; k < cnt; k++) memcpy((char*)buf + esz * k, src + esz * pos[k], esz);
        return out;
    }
    out = expr_new_function(expr_new_symbol(SYM_List), NULL, cnt);
    for (size_t k = 0; k < cnt; k++) out->data.function.args[k] = elem_at(arr, pos[k]);
    return out;
}

/* ------------------------------------------------------------ PositionIndex */

Expr* assoc_packed_positionindex(const Expr* arr) {
    if (!keyable_vector(arr)) return NULL;
    size_t n = (size_t)arr->data.ndarray.dims[0];
    NDWordGroups g;
    if (!nd_group_words(arr->data.ndarray.data, arr->data.ndarray.dtype, n, &g)) return NULL;
    size_t nu = g.nuniq;
    Expr** rules = malloc(sizeof(Expr*) * nu);
    size_t *perm = NULL, *off = NULL;
    if (!rules || !group_permutation(&g, n, &perm, &off)) {
        free(rules); nd_word_groups_free(&g); return NULL;
    }
    /* Allocate the position Integers group by group, in the order the result
     * will later be freed: freeing 10^6 nodes in allocation order is markedly
     * cheaper for the allocator than freeing them strided across the heap,
     * which is what filling by input position produced. */
    for (size_t u = 0; u < nu; u++) {
        size_t cnt = off[u + 1] - off[u];
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), NULL, cnt);
        Expr** slot = list->data.function.args;
        const size_t* pos = perm + off[u];
        for (size_t k = 0; k < cnt; k++) slot[k] = expr_new_integer((int64_t)pos[k] + 1);
        rules[u] = stamp_fixed(make_rule2(word_to_expr(g.keys[u], arr->data.ndarray.dtype),
                                          stamp_fixed(list)));
    }
    Expr* assoc = stamp_fixed(make_assoc(rules, nu));
    free(rules); free(perm); free(off);
    nd_word_groups_free(&g);
    return assoc;
}

/* --------------------------------------------------------- GroupBy family */

/* Run the key function over the buffer and group by the resulting words.
 * On true the caller owns *keys (the key buffer) and *g. */
static bool group_by_keyfn(const Expr* arr, const Expr* f, Expr** keys, NDWordGroups* g) {
    *keys = NULL;
    if (!keyable_vector(arr) || is_ndarray(f)) return false;
    Expr* kb = autocompile_map_unary(f, arr);
    if (!kb) return false;
    if (!nd_group_words(kb->data.ndarray.data, kb->data.ndarray.dtype,
                        (size_t)kb->data.ndarray.dims[0], g)) {
        expr_free(kb);
        return false;
    }
    *keys = kb;
    return true;
}

Expr* assoc_packed_groupby(const Expr* arr, const Expr* f, const Expr* reducer) {
    Expr* kb; NDWordGroups g;
    if (reducer && is_ndarray(reducer)) return NULL;
    if (!group_by_keyfn(arr, f, &kb, &g)) return NULL;
    size_t nu = g.nuniq;
    NDType kdt = kb->data.ndarray.dtype;
    Expr** rules = malloc(sizeof(Expr*) * nu);
    if (!rules) { expr_free(kb); nd_word_groups_free(&g); return NULL; }

    /* Without a reducer every group sits inside a Rule, so all are boxed. With
     * one, a large group is handed to the reducer as a PACKED top-level
     * argument -- Total/Mean/Max/... then reduce the buffer -- and the reducer's
     * answer is materialised if it is itself a buffer (pack_eval_plain), so
     * nothing packed is left inside the association. */
    Expr** lists = NULL;
    size_t *perm = NULL, *off = NULL;
    bool ok = reducer ? group_permutation(&g, (size_t)arr->data.ndarray.dims[0], &perm, &off)
                      : (lists = build_group_lists(arr, &g)) != NULL;
    if (!ok) { free(rules); expr_free(kb); nd_word_groups_free(&g); return NULL; }
    for (size_t u = 0; u < nu; u++) {
        Expr* value;
        if (!reducer) {
            value = lists[u];
        } else {
            Expr* grp = build_group_arg(arr, perm + off[u], off[u + 1] - off[u]);
            Expr* call = expr_new_function(expr_copy((Expr*)reducer), &grp, 1);
            value = pack_eval_plain(call);
            expr_free(call);
        }
        rules[u] = make_rule2(word_to_expr(g.keys[u], kdt), value);
    }
    /* Groups of machine numbers, or reducer answers evaluate() just returned:
     * nothing here reduces further, but a reducer may have mutated the symbol
     * table (and so bumped the clock), in which case the stamps just lapse. */
    for (size_t u = 0; u < nu; u++) stamp_fixed(rules[u]);
    Expr* assoc = stamp_fixed(make_assoc(rules, nu));
    free(rules); free(lists); free(perm); free(off);
    expr_free(kb);
    nd_word_groups_free(&g);
    return assoc;
}

Expr* assoc_packed_gatherby(const Expr* arr, const Expr* f) {
    Expr* kb; NDWordGroups g;
    if (!group_by_keyfn(arr, f, &kb, &g)) return NULL;
    Expr** lists = build_group_lists(arr, &g);
    Expr* out = lists ? stamp_fixed(expr_new_function(expr_new_symbol(SYM_List), lists, g.nuniq))
                      : NULL;
    free(lists);
    expr_free(kb);
    nd_word_groups_free(&g);
    return out;
}

static Expr* counts_from_groups(const NDWordGroups* g, NDType kdt) {
    Expr** rules = malloc(sizeof(Expr*) * (g->nuniq ? g->nuniq : 1));
    if (!rules) return NULL;
    for (size_t u = 0; u < g->nuniq; u++)
        rules[u] = stamp_fixed(make_rule2(word_to_expr(g->keys[u], kdt),
                                          expr_new_integer(g->cnts[u])));
    Expr* out = stamp_fixed(make_assoc(rules, g->nuniq));
    free(rules);
    return out;
}

Expr* assoc_packed_countsby(const Expr* arr, const Expr* f) {
    Expr* kb; NDWordGroups g;
    if (!group_by_keyfn(arr, f, &kb, &g)) return NULL;
    Expr* out = counts_from_groups(&g, kb->data.ndarray.dtype);
    expr_free(kb);
    nd_word_groups_free(&g);
    return out;
}

Expr* assoc_packed_counts_words(const Expr* keys) {
    if (!keys || !is_ndarray(keys) || keys->data.ndarray.rank != 1) return NULL;
    NDWordGroups g;
    if (!nd_group_words(keys->data.ndarray.data, keys->data.ndarray.dtype,
                        (size_t)keys->data.ndarray.dims[0], &g)) return NULL;
    Expr* out = counts_from_groups(&g, keys->data.ndarray.dtype);
    nd_word_groups_free(&g);
    return out;
}

/* ------------------------------------------- AssociationThread / Map */

/* Rules for de-duplicated buffer keys: the first occurrence of a key fixes its
 * position, the LAST its value (assoc_from_rules' rule). `value_at(i)` boxes
 * the value of entry i. */
typedef Expr* (*ValueAt)(const void* ctx, size_t i);

static Expr* thread_by_words(const Expr* keys, ValueAt value_at, const void* ctx,
                             bool stamp) {
    size_t n = (size_t)keys->data.ndarray.dims[0];
    NDWordGroups g;
    if (!nd_group_words(keys->data.ndarray.data, keys->data.ndarray.dtype, n, &g)) return NULL;
    size_t nu = g.nuniq;
    size_t* last = malloc(sizeof(size_t) * nu);
    Expr** rules = malloc(sizeof(Expr*) * nu);
    if (!last || !rules) { free(last); free(rules); nd_word_groups_free(&g); return NULL; }
    for (size_t i = 0; i < n; i++) last[g.gid[i]] = i;
    for (size_t u = 0; u < nu; u++) {
        rules[u] = make_rule2(word_to_expr(g.keys[u], keys->data.ndarray.dtype),
                              value_at(ctx, last[u]));
        if (stamp) stamp_fixed(rules[u]);
    }
    Expr* assoc = make_assoc(rules, nu);
    if (stamp) stamp_fixed(assoc);
    free(last); free(rules);
    nd_word_groups_free(&g);
    return assoc;
}

static Expr* value_from_buffer(const void* ctx, size_t i) { return elem_at((const Expr*)ctx, i); }
static Expr* value_from_list(const void* ctx, size_t i) {
    return expr_copy(((const Expr*)ctx)->data.function.args[i]);
}

Expr* assoc_packed_thread(const Expr* keys, const Expr* vals) {
    if (!keyable_vector(keys)) return NULL;
    size_t n = (size_t)keys->data.ndarray.dims[0];
    if (is_ndarray(vals)) {
        if (vals->data.ndarray.rank != 1 || (size_t)vals->data.ndarray.dims[0] != n) return NULL;
        return thread_by_words(keys, value_from_buffer, vals, true);
    }
    if (vals->type != EXPR_FUNCTION || vals->data.function.head->type != EXPR_SYMBOL ||
        vals->data.function.head->data.symbol.name != SYM_List ||
        vals->data.function.arg_count != n)
        return NULL;
    return thread_by_words(keys, value_from_list, vals, false);
}

typedef struct { const Expr* f; const Expr* keys; } MapCtx;

/* f[k], left for the evaluator to reduce -- the List path's shape. */
static Expr* value_unevaluated(const void* ctx, size_t i) {
    const MapCtx* m = (const MapCtx*)ctx;
    Expr* k = elem_at(m->keys, i);
    return expr_new_function(expr_copy((Expr*)m->f), &k, 1);
}

Expr* assoc_packed_map(const Expr* f, const Expr* arr) {
    if (!keyable_vector(arr) || is_ndarray(f)) return NULL;
    /* All of f[k] in one compiled loop when the exactness gate admits f; the
     * values are then read off that buffer. Otherwise the rules carry f[k]
     * unevaluated, exactly as the List path builds them. */
    Expr* vb = autocompile_map_unary(f, arr);
    Expr* out;
    if (vb) {
        out = thread_by_words(arr, value_from_buffer, vb, true);
        expr_free(vb);
    } else {
        MapCtx m = { f, arr };
        out = thread_by_words(arr, value_unevaluated, &m, false);
    }
    return out;
}

/* ----------------------------------------------------------------- Lookup */

Expr* assoc_packed_lookup(const Expr* assoc, const Expr* keys, const Expr* deflt) {
    if (!keyable_vector(keys) || !is_association(assoc)) return NULL;
    size_t m = (size_t)keys->data.ndarray.dims[0];
    NDType dt = keys->data.ndarray.dtype;
    const void* kd = keys->data.ndarray.data;
    Expr** found = malloc(sizeof(Expr*) * m);        /* borrowed values, NULL = absent */
    if (!found) return NULL;
    bool all_int = true, all_real = true;
    for (size_t j = 0; j < m; j++) {
        Expr* v = (dt == NDT_INT64) ? assoc_lookup_value_i64(assoc, ((const int64_t*)kd)[j])
                                    : assoc_lookup_value_real(assoc, ((const double*)kd)[j]);
        if (!v) v = (Expr*)deflt;
        found[j] = v;
        if (!v || v->type != EXPR_INTEGER) all_int = false;
        if (!v || v->type != EXPR_REAL)    all_real = false;
    }
    Expr* out = NULL;
    if ((all_int || all_real) && pack_enabled() && m >= pack_min_elements()) {
        int64_t len = (int64_t)m;
        void* buf = NULL;
        out = ndbuild_open(1, &len, all_int ? NDT_INT64 : NDT_FLOAT64, &buf);
        if (out) {
            if (all_int) for (size_t j = 0; j < m; j++) ((int64_t*)buf)[j] = found[j]->data.integer;
            else         for (size_t j = 0; j < m; j++) ((double*)buf)[j] = found[j]->data.real;
        }
    }
    if (!out) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, m);
        for (size_t j = 0; j < m; j++) {
            if (found[j]) { out->data.function.args[j] = expr_copy(found[j]); continue; }
            Expr* margs[2] = { expr_new_string("KeyAbsent"), elem_at(keys, j) };
            out->data.function.args[j] = expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
        }
    }
    free(found);
    return out;
}

/* ---------------------------------------------------------- Keys / Values */

Expr* assoc_pack_column(const Expr* assoc, bool want_keys) {
    if (!pack_enabled()) return NULL;
    size_t n = assoc->data.function.arg_count;
    if (n == 0 || n < pack_min_elements()) return NULL;
    Expr** a = assoc->data.function.args;
    int col = want_keys ? 0 : 1;
    const Expr* first = a[0]->data.function.args[col];
    NDType dt;
    if (first->type == EXPR_INTEGER)      dt = NDT_INT64;
    else if (first->type == EXPR_REAL)    dt = NDT_FLOAT64;
    else if (first->type == EXPR_SYMBOL &&
             (first->data.symbol.name == SYM_True || first->data.symbol.name == SYM_False))
        dt = NDT_BOOL;
    else return NULL;

    /* Class check first, so a mixed column costs no allocation. */
    for (size_t i = 0; i < n; i++) {
        const Expr* x = a[i]->data.function.args[col];
        bool ok = (dt == NDT_INT64)   ? x->type == EXPR_INTEGER
                : (dt == NDT_FLOAT64) ? x->type == EXPR_REAL
                : (x->type == EXPR_SYMBOL &&
                   (x->data.symbol.name == SYM_True || x->data.symbol.name == SYM_False));
        if (!ok) return NULL;
    }
    int64_t len = (int64_t)n;
    void* buf = NULL;
    Expr* out = ndbuild_open(1, &len, dt, &buf);
    if (!out) return NULL;
    if (dt == NDT_INT64) {
        int64_t* iv = (int64_t*)buf;
        for (size_t i = 0; i < n; i++) iv[i] = a[i]->data.function.args[col]->data.integer;
    } else if (dt == NDT_FLOAT64) {
        double* dv = (double*)buf;
        for (size_t i = 0; i < n; i++) dv[i] = a[i]->data.function.args[col]->data.real;
    } else {
        uint8_t* bv = (uint8_t*)buf;
        for (size_t i = 0; i < n; i++)
            bv[i] = a[i]->data.function.args[col]->data.symbol.name == SYM_True;
    }
    return out;
}
