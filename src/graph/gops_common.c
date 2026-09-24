/* gops_common.c - shared internals of the graph "ops" stream (graph_ops.h).
 *
 *   gops_view / gops_view_own / gops_view_free
 *       a graph's parts as raw arrays plus the validated-graph memo's integer
 *       endpoints (eu/ev/edir), so every edit is an integer pass over the edges;
 *   gops_graph_new
 *       the single result builder: assembles the canonical Graph and seeds the
 *       memo from the endpoint arrays the operation already computed;
 *   gops_edge / gops_heads_free
 *       edge nodes sharing one head symbol node per call;
 *   gops_parse_edge
 *       DirectedEdge/UndirectedEdge and the ->/<-> sugar;
 *   gops_items
 *       the "v | {v..}" argument convention.
 *
 * Why seed: a builtin's result is evaluated again, and Graph[...] on an unknown
 * node validates it from scratch -- a vertex hash index, two endpoint lookups
 * per edge and a parallel-edge set. Seeding hands graph_memo_seed the endpoint
 * indices directly, so only the vertex index is built (one hash per vertex),
 * and the constructor's check on return is an O(1) memo hit.
 *
 * Memory (SPEC section 4): nothing here frees or mutates a caller's graph.
 */

#include "graph_ops.h"
#include "sym_names.h"
#include "match.h"
#include "eval.h"
#include "symtab.h"
#include <stdlib.h>
#include <string.h>

void* gops_calloc(size_t n, size_t size) { return calloc(n > 0 ? n : 1, size); }
void* gops_malloc(size_t n, size_t size) { return malloc((n > 0 ? n : 1) * size); }

Expr* gops_truth(int b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

Expr* gops_list_take(Expr** items, size_t n) {
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return l;
}

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

int gops_view(const Expr* g, GopsView* v) {
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return 0;
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    memset(v, 0, sizeof(*v));
    v->g = g;
    v->verts = verts->data.function.args;
    v->nv = verts->data.function.arg_count;
    v->edges = edges->data.function.args;
    v->ne = edges->data.function.arg_count;
    if (g->data.function.arg_count == 3)
        v->weights = g->data.function.args[2]->data.function.args[1]->data.function.args;
    v->eu = eu; v->ev = ev; v->edir = edir;
    v->ndir = (size_t)graph_directed_edge_count(g);
    return 1;
}

int gops_view_own(GopsView* v) {
    if (v->owned) return 1;
    size_t ne = v->ne;
    int* eu = gops_malloc(ne, sizeof(int));
    int* ev = gops_malloc(ne, sizeof(int));
    unsigned char* ed = gops_malloc(ne, 1);
    if (!eu || !ev || !ed) { free(eu); free(ev); free(ed); return 0; }
    if (ne) {
        memcpy(eu, v->eu, ne * sizeof(int));
        memcpy(ev, v->ev, ne * sizeof(int));
        memcpy(ed, v->edir, ne);
    }
    v->eu = eu; v->ev = ev; v->edir = ed;
    v->owned = 1;
    return 1;
}

void gops_view_free(GopsView* v) {
    if (v && v->owned) {
        free((void*)v->eu); free((void*)v->ev); free((void*)v->edir);
        v->owned = 0;
    }
}

Expr* gops_edge(GopsHeads* hs, int directed, Expr* u, Expr* v) {
    int k = directed ? 1 : 0;
    if (!hs->h[k]) hs->h[k] = expr_new_symbol(k ? SYM_DirectedEdge : SYM_UndirectedEdge);
    Expr* a[2] = { u, v };
    return expr_new_function(expr_copy(hs->h[k]), a, 2);
}

void gops_heads_free(GopsHeads* hs) {
    if (hs->h[0]) expr_free(hs->h[0]);
    if (hs->h[1]) expr_free(hs->h[1]);
    hs->h[0] = hs->h[1] = NULL;
}

int gops_parse_edge(const Expr* e, const Expr** u, const Expr** v) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return -1;
    int d;
    if (head_is(e, SYM_DirectedEdge) || head_is(e, SYM_Rule)) d = 1;
    else if (head_is(e, SYM_UndirectedEdge) || head_is(e, SYM_TwoWayRule)) d = 0;
    else return -1;
    *u = e->data.function.args[0];
    *v = e->data.function.args[1];
    return d;
}

Expr* const* gops_items(Expr* const* slot, size_t* n) {
    if (graph_is_list(*slot)) {
        *n = (*slot)->data.function.arg_count;
        return (*slot)->data.function.args;
    }
    *n = 1;
    return slot;
}

/* True iff none of the heads a graph result is built from carries user rules,
 * so re-evaluating the result would reproduce it exactly. */
static int graph_heads_inert(void) {
    const char* heads[6];
    heads[0] = SYM_Graph; heads[1] = SYM_List; heads[2] = SYM_DirectedEdge;
    heads[3] = SYM_UndirectedEdge; heads[4] = SYM_Rule; heads[5] = SYM_EdgeWeight;
    for (int i = 0; i < 6; i++) {
        SymbolDef* d = symtab_lookup(heads[i]);
        if (d && (d->down_values || d->own_values)) return 0;
    }
    return 1;
}

static void free_exprs(Expr** a, size_t n) {
    if (!a) return;
    for (size_t i = 0; i < n; i++) if (a[i]) expr_free(a[i]);
    free(a);
}

Expr* gops_graph_new(Expr** verts, size_t nv, Expr** edges, size_t ne,
                     Expr** weights, int* eu, int* ev, unsigned char* edir) {
    GraphVIdx* ix = graph_vidx_new(nv);
    if (!ix) goto fail;
    for (size_t i = 0; i < nv; i++) {
        if (!graph_vidx_put(ix, verts[i], (int)i)) goto fail;   /* repeated vertex */
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), edges, ne);
    free(verts); free(edges);
    Expr* g;
    if (weights) {
        Expr* wlist = expr_new_function(expr_new_symbol(SYM_List), weights, ne);
        free(weights);
        Expr* ra[2] = { expr_new_symbol(SYM_EdgeWeight), wlist };
        Expr* ga[3] = { vlist, elist, expr_new_function(expr_new_symbol(SYM_Rule), ra, 2) };
        g = expr_new_function(expr_new_symbol(SYM_Graph), ga, 3);
    } else {
        Expr* ga[2] = { vlist, elist };
        g = expr_new_function(expr_new_symbol(SYM_Graph), ga, 2);
    }
    /* graph_memo_seed owns ix/eu/ev/edir from here, success or not. */
    if (!graph_memo_seed(g, ix, eu, ev, edir)) { expr_free(g); return NULL; }
    /* Stamp the result as evaluated under the live clock. It is a canonical,
     * validated Graph over already-evaluated parts, so evaluating it again is
     * the identity -- but the evaluator would rebuild the node (a fresh
     * pointer) and Graph[] would re-validate that copy from scratch, O(V + E)
     * with hashing, before settling on it. The stamp makes the evaluator's
     * fixed-point check return THIS (memo-seeded) node at once. Skipped if a
     * user has attached rules to any of the heads involved. */
    if (graph_heads_inert()) g->last_evaluated_at = eval_clock_get();
    return g;

fail:
    graph_vidx_free(ix);
    free_exprs(verts, nv);
    free_exprs(edges, ne);
    free_exprs(weights, ne);
    free(eu); free(ev); free(edir);
    return NULL;
}

/* ---- Incidence CSR -------------------------------------------------------- */
int gops_inc_build(const GopsView* v, int mode, GopsInc* inc) {
    int want_eid = (mode != GOPS_INC_OUT_NBR);
    if (mode == GOPS_INC_OUT_NBR) mode = GOPS_INC_OUT;
    size_t n = v->nv, ne = v->ne;
    memset(inc, 0, sizeof(*inc));
    inc->n = (int)n;
    inc->start = gops_calloc(n + 1, sizeof(int));
    if (!inc->start) return 0;
    size_t slots = 0;
    for (size_t k = 0; k < ne; k++) {
        int d = v->edir[k];
        if (mode == GOPS_INC_ALL || !d) {
            inc->start[v->eu[k] + 1]++; inc->start[v->ev[k] + 1]++; slots += 2;
        } else if (mode == GOPS_INC_OUT) {
            inc->start[v->eu[k] + 1]++; slots++;
        } else {                                   /* GOPS_INC_IN */
            inc->start[v->ev[k] + 1]++; slots++;
        }
    }
    for (size_t i = 0; i < n; i++) inc->start[i + 1] += inc->start[i];
    inc->nbr = gops_malloc(slots, sizeof(int));
    if (want_eid) inc->eid = gops_malloc(slots, sizeof(int));
    if (!inc->nbr || (want_eid && !inc->eid)) { gops_inc_free(inc); return 0; }
    if (!want_eid) {                       /* neighbours only: half the writes */
        int* cur = inc->start + 1;
        for (size_t i = n; i > 0; i--) cur[i - 1] = inc->start[i - 1];
        int* nb = inc->nbr;
        for (size_t k = 0; k < ne; k++) {
            int a = v->eu[k], b = v->ev[k];
            nb[cur[a]++] = b;
            if (!v->edir[k]) nb[cur[b]++] = a;
        }
        inc->start[0] = 0;
        return 1;
    }
    /* Fill using start[i+1] as row i's cursor (rows stay in EdgeList order);
     * afterwards start[i+1] has advanced to row i's end = row i+1's start,
     * and the counts shifted one slot recover the exclusive prefix. */
    int* cur = inc->start + 1;
    for (size_t i = n; i > 0; i--) cur[i - 1] = inc->start[i - 1];
    int* nb = inc->nbr; int* ei = inc->eid;
    const int* eu = v->eu; const int* ev = v->ev; const unsigned char* ed = v->edir;
    for (size_t k = 0; k < ne; k++) {
        int a = eu[k], b = ev[k], d = ed[k];
        if (mode == GOPS_INC_ALL || !d) {
            nb[cur[a]] = b; ei[cur[a]++] = (int)k;
            nb[cur[b]] = a; ei[cur[b]++] = (int)k;
        } else if (mode == GOPS_INC_OUT) {
            nb[cur[a]] = b; ei[cur[a]++] = (int)k;
        } else {
            nb[cur[b]] = a; ei[cur[b]++] = (int)k;
        }
    }
    inc->start[0] = 0;
    return 1;
}

void gops_inc_free(GopsInc* inc) {
    free(inc->start); free(inc->nbr); free(inc->eid);
    memset(inc, 0, sizeof(*inc));
}

/* ---- Stable counting sort ------------------------------------------------- */
int gops_csort(int* perm, size_t m, const int* key, int nkeys) {
    if (m < 2) return 1;
    int* cnt = gops_calloc((size_t)nkeys + 1, sizeof(int));
    int* tmp = gops_malloc(m, sizeof(int));
    if (!cnt || !tmp) { free(cnt); free(tmp); return 0; }
    for (size_t i = 0; i < m; i++) cnt[key[perm[i]] + 1]++;
    for (int k = 0; k < nkeys; k++) cnt[k + 1] += cnt[k];
    for (size_t i = 0; i < m; i++) tmp[cnt[key[perm[i]]]++] = perm[i];
    memcpy(perm, tmp, m * sizeof(int));
    free(cnt); free(tmp);
    return 1;
}

/* ---- Patterns ------------------------------------------------------------- */
int gops_has_pattern(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL) {
        const char* s = h->data.symbol.name;
        if (s == SYM_Blank || s == SYM_BlankSequence || s == SYM_BlankNullSequence
            || s == SYM_Pattern || s == SYM_PatternTest || s == SYM_Condition
            || s == SYM_Alternatives || s == SYM_Except)
            return 1;
    } else if (gops_has_pattern(h)) {
        return 1;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (gops_has_pattern(e->data.function.args[i])) return 1;
    return 0;
}

int gops_matchq(Expr* e, Expr* patt) {
    MatchEnv* env = env_new();
    int ok = match(e, patt, env) ? 1 : 0;
    env_free(env);
    return ok;
}

/* ---- Infinity / integer arguments ----------------------------------------- */
int gops_is_infinity(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Infinity;
    return head_is(e, SYM_DirectedInfinity) && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer == 1;
}

/* ---- Integer edge-key set -------------------------------------------------- */
uint64_t gops_edge_key(int a, int b, int directed) {
    if (!directed && a > b) { int t = a; a = b; b = t; }
    return ((uint64_t)(uint32_t)a << 33) | ((uint64_t)(uint32_t)b << 1) | (uint64_t)(directed ? 1 : 0);
}

int gops_keyset_init(GopsKeySet* s, size_t n) {
    size_t cap = 16;
    while (cap < 2 * (n + 1)) cap <<= 1;
    s->key = malloc(cap * sizeof(uint64_t));
    s->val = malloc(cap * sizeof(int));
    if (!s->key || !s->val) { free(s->key); free(s->val); s->key = NULL; s->val = NULL; return 0; }
    memset(s->key, 0xFF, cap * sizeof(uint64_t));
    s->mask = cap - 1;
    return 1;
}

void gops_keyset_free(GopsKeySet* s) { free(s->key); free(s->val); s->key = NULL; s->val = NULL; }

static size_t keyset_slot(const GopsKeySet* s, uint64_t k) {
    size_t h = (size_t)((k * 0x9E3779B97F4A7C15ULL) >> 29) & s->mask;
    while (s->key[h] != UINT64_MAX && s->key[h] != k) h = (h + 1) & s->mask;
    return h;
}

int gops_keyset_put(GopsKeySet* s, uint64_t k, int val) {
    size_t h = keyset_slot(s, k);
    if (s->key[h] == UINT64_MAX) { s->key[h] = k; s->val[h] = val; }
    return s->val[h];
}

int gops_keyset_get(const GopsKeySet* s, uint64_t k) {
    size_t h = keyset_slot(s, k);
    return s->key[h] == UINT64_MAX ? -1 : s->val[h];
}

/* ---- Per-graph incidence cache ---------------------------------------------
 * One slot per mode, keyed on the graph node POINTER and holding a reference
 * to it -- sound for the same reason the validated-graph memo is: a referenced
 * node cannot be freed or recycled, and a shared node is immutable. So a
 * repeated FindPath / FindEulerianCycle / FindCycle on one graph skips the
 * O(V + E) CSR build, as Mathematica's Graph object caches its adjacency.
 * Keeps at most four graphs alive past their last user reference. */
static struct { Expr* g; GopsInc inc; } g_inc_cache[4];

const GopsInc* gops_inc_cached(const GopsView* v, int mode) {
    if (mode < 0 || mode > 3) return NULL;
    if (g_inc_cache[mode].g == v->g) return &g_inc_cache[mode].inc;
    GopsInc fresh;
    if (!gops_inc_build(v, mode, &fresh)) return NULL;
    if (g_inc_cache[mode].g) {
        expr_free(g_inc_cache[mode].g);
        gops_inc_free(&g_inc_cache[mode].inc);
    }
    g_inc_cache[mode].g = expr_copy((Expr*)v->g);
    g_inc_cache[mode].inc = fresh;
    return &g_inc_cache[mode].inc;
}
