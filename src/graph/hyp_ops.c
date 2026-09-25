/* hyp_ops.c - structural hypergraph operations.
 *
 *   HyperedgeSizes[h]                 arity (Length) of each hyperedge
 *   HypergraphRank[h] / Corank[h]     max / min arity (0 for no hyperedges)
 *   UniformHypergraphQ[h(, k)]        all arities equal (to k)
 *   HypergraphDual[h]                 vertices <-> hyperedges
 *   HypergraphCliqueExpansion[h]      2-section Graph
 *   HypergraphStarExpansion[h]        incidence bipartite Graph
 *   HypergraphToGraph[h]              ordered-hyperedge directed Graph (FR)
 *   HypergraphLineGraph[h(, s)]       s-line graph (Aksoy et al.)
 *   HypergraphConnectedComponents[h]  vertex components
 *   ConnectedHypergraphQ[h]           (FR name)
 *   HyperedgeConnectedComponents[h(, s)]  s-connected components of hyperedges
 *   HyperedgeDistance[h, i, j(, s)]   s-walk distance between hyperedges
 *   HypergraphDistance[h, u, v]       vertex distance (hyperedge hops)
 *   HypergraphVertexAdd/Delete, HypergraphEdgeAdd/Delete, Subhypergraph,
 *   HypergraphRestriction             edits and sub-hypergraphs
 *
 * Every set-theoretic operation reads hyperedge j as the set of its distinct
 * vertices (HypView.sv); only HyperedgeSizes/Rank/Corank/Uniform (arity) and
 * HypergraphToGraph (order) read the raw List. All algorithms are linear in the
 * total incidence count sum|e|, except the s-line graph family, which is
 * O(sum_v deg(v)^2) -- the number of hyperedge pairs meeting at a vertex, which
 * bounds the output of the 1-line graph anyway.
 *
 * Hyperedge indices in inputs and outputs are 1-based positions in EdgeList.
 *
 * Memory (SPEC section 4): results are freshly built; `res` is never consumed.
 */

#include "graph_hyper.h"
#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "sym_intern.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers -------------------------------------------------------------- */

static Expr* mk_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

static Expr* mk_hyp(Expr* vlist, Expr* elist) {
    Expr* a[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), a, 2);
}

static Expr* mk_graph(Expr* vlist, Expr* elist) {
    Expr* a[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), a, 2);
}

/* Positive machine integer argument, else -1. */
static long pos_int(const Expr* e) {
    if (!e || e->type != EXPR_INTEGER || e->data.integer < 1) return -1;
    return (long)e->data.integer;
}

/* Accept a Hypergraph, or (FR compatibility) a plain List of hyperedge Lists,
 * which is turned into Hypergraph[list] and evaluated. On success *owned is
 * the object to free afterwards (NULL when the argument itself was used). */
static const Expr* hyp_arg(const Expr* a, Expr** owned) {
    *owned = NULL;
    if (hypergraph_is_valid(a)) return a;
    if (!graph_is_list(a)) return NULL;
    Expr* arg = expr_copy((Expr*)a);
    Expr* call = expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), &arg, 1);
    Expr* h = evaluate(call);
    expr_free(call);
    if (!hypergraph_is_valid(h)) { expr_free(h); return NULL; }
    *owned = h;
    return h;
}

/* Growable open-addressing set of nonzero uint64 keys. */
typedef struct { uint64_t* k; size_t mask, count; } U64Set;

static int u64_init(U64Set* s, size_t hint) {
    size_t cap = 16;
    while (cap < hint * 2 + 2) cap <<= 1;
    s->k = calloc(cap, sizeof(uint64_t));
    s->mask = cap - 1; s->count = 0;
    return s->k != NULL;
}
static size_t u64_slot(const U64Set* s, uint64_t key) {
    size_t i = (size_t)((key * 0x9E3779B97F4A7C15ULL) >> 20) & s->mask;
    while (s->k[i] && s->k[i] != key) i = (i + 1) & s->mask;
    return i;
}
/* 1 if inserted, 0 if present, -1 on allocation failure. key must be != 0. */
static int u64_insert(U64Set* s, uint64_t key) {
    size_t i = u64_slot(s, key);
    if (s->k[i]) return 0;
    s->k[i] = key; s->count++;
    if (s->count * 2 > s->mask + 1) {
        U64Set b;
        if (!u64_init(&b, (s->mask + 1))) return -1;
        for (size_t t = 0; t <= s->mask; t++)
            if (s->k[t]) b.k[u64_slot(&b, s->k[t])] = s->k[t];
        b.count = s->count;
        free(s->k);
        *s = b;
    }
    return 1;
}

static uint64_t pair_key(int a, int b) {           /* never 0: +1 offsets */
    return ((uint64_t)(uint32_t)(a + 1) << 32) | (uint64_t)(uint32_t)(b + 1);
}

/* Growable Expr* vector. */
typedef struct { Expr** a; size_t n, cap; } EVec;
static int evec_push(EVec* v, Expr* e) {
    if (v->n == v->cap) {
        size_t nc = v->cap ? v->cap * 2 : 64;
        Expr** na = realloc(v->a, nc * sizeof(Expr*));
        if (!na) { expr_free(e); return 0; }
        v->a = na; v->cap = nc;
    }
    v->a[v->n++] = e;
    return 1;
}
static void evec_free(EVec* v) {
    for (size_t i = 0; i < v->n; i++) expr_free(v->a[i]);
    free(v->a);
}
static Expr* evec_to_list(EVec* v) {
    Expr* out = mk_list(v->a, v->n);
    free(v->a);
    v->a = NULL; v->n = v->cap = 0;
    return out;
}

static Expr* edge2(Expr* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_copy(head), args, 2);
}

static int cmp_int(const void* a, const void* b) {
    int x = *(const int*)a, y = *(const int*)b;
    return (x > y) - (x < y);
}

/* ---- sizes / rank / uniformity -------------------------------------------- */

Expr* builtin_hyperedge_sizes(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    int64_t* d = malloc((size_t)(v.m > 0 ? v.m : 1) * sizeof(int64_t));
    if (!d) return NULL;
    for (int j = 0; j < v.m; j++) d[j] = v.eoff[j + 1] - v.eoff[j];
    Expr* out = hyp_int_list(d, (size_t)v.m);
    free(d);
    return out;
}

static Expr* rank_impl(Expr* res, int want_max) {
    if (res->data.function.arg_count != 1) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    int best = 0;
    for (int j = 0; j < v.m; j++) {
        int k = v.eoff[j + 1] - v.eoff[j];
        if (j == 0 || (want_max ? k > best : k < best)) best = k;
    }
    return expr_new_integer(best);
}
Expr* builtin_hypergraph_rank(Expr* res)   { return rank_impl(res, 1); }
Expr* builtin_hypergraph_corank(Expr* res) { return rank_impl(res, 0); }

Expr* builtin_uniform_hypergraph_q(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0))
        return argc == 1 ? expr_new_symbol(SYM_False) : NULL;
    long k = -1;
    if (argc == 2) {
        const Expr* ke = res->data.function.args[1];
        if (!ke || ke->type != EXPR_INTEGER || ke->data.integer < 0) return NULL;
        k = (long)ke->data.integer;
    }
    int ok = 1;
    for (int j = 0; j < v.m && ok; j++) {
        long a = v.eoff[j + 1] - v.eoff[j];
        if (k < 0) k = a;
        if (a != k) ok = 0;
    }
    return expr_new_symbol(ok ? SYM_True : SYM_False);
}

/* ---- dual ----------------------------------------------------------------- */

Expr* builtin_hypergraph_dual(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 1)) return NULL;
    Expr** es = malloc((size_t)(v.n > 0 ? v.n : 1) * sizeof(Expr*));
    Expr** tmp = malloc((size_t)(v.m > 0 ? v.m : 1) * sizeof(Expr*));
    Expr** ints = malloc((size_t)(v.m > 0 ? v.m : 1) * sizeof(Expr*));
    if (!es || !tmp || !ints) { free(es); free(tmp); free(ints); return NULL; }
    for (int j = 0; j < v.m; j++) ints[j] = expr_new_integer(j + 1);
    for (int i = 0; i < v.n; i++) {
        int c = 0;
        for (int t = v.voff[i]; t < v.voff[i + 1]; t++) tmp[c++] = expr_copy(ints[v.ve[t]]);
        es[i] = mk_list(tmp, (size_t)c);
    }
    Expr* vl = mk_list(ints, (size_t)v.m);        /* takes the ints' own refs */
    Expr* out = mk_hyp(vl, mk_list(es, (size_t)v.n));
    free(es); free(tmp); free(ints);
    return out;
}

/* ---- graph views ---------------------------------------------------------- */

/* 2-section: u <-> w whenever u and w share a hyperedge. Edge order: first
 * co-occurrence, oriented as the pair first appears in its hyperedge. */
Expr* builtin_hypergraph_clique_expansion(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned;
    const Expr* h = hyp_arg(res->data.function.args[0], &owned);
    HypView v;
    if (!h || !hyp_view(h, &v, 0)) { expr_free(owned); return NULL; }
    U64Set seen;
    if (!u64_init(&seen, (size_t)v.soff[v.m])) { expr_free(owned); return NULL; }
    EVec ev = {0};
    Expr* head = expr_new_symbol(SYM_UndirectedEdge);
    Expr* const* vx = v.verts->data.function.args;
    int ok = 1;
    for (int j = 0; j < v.m && ok; j++) {
        for (int a = v.soff[j]; a < v.soff[j + 1] && ok; a++)
            for (int b = a + 1; b < v.soff[j + 1] && ok; b++) {
                int x = v.sv[a], y = v.sv[b];
                int r = u64_insert(&seen, x < y ? pair_key(x, y) : pair_key(y, x));
                if (r < 0) ok = 0;
                else if (r == 1)
                    ok = evec_push(&ev, edge2(head, expr_copy(vx[x]), expr_copy(vx[y])));
            }
    }
    expr_free(head);
    free(seen.k);
    if (!ok) { evec_free(&ev); expr_free(owned); return NULL; }
    Expr* out = mk_graph(expr_copy((Expr*)v.verts), evec_to_list(&ev));
    expr_free(owned);
    return out;
}

/* Incidence (star) expansion: the bipartite Graph on VertexList[h] together
 * with one node Hyperedge[j] per hyperedge, v <-> Hyperedge[j] iff v is in e_j.
 * Unevaluated if some vertex of h is itself a Hyperedge[j] node. */
Expr* builtin_hypergraph_star_expansion(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned;
    const Expr* h = hyp_arg(res->data.function.args[0], &owned);
    HypView v;
    if (!h || !hyp_view(h, &v, 0)) { expr_free(owned); return NULL; }
    size_t nv = (size_t)v.n + (size_t)v.m;
    Expr** vs = malloc((nv ? nv : 1) * sizeof(Expr*));
    size_t ne = (size_t)v.soff[v.m];
    Expr** es = malloc((ne ? ne : 1) * sizeof(Expr*));
    if (!vs || !es) { free(vs); free(es); expr_free(owned); return NULL; }
    Expr* const* vx = v.verts->data.function.args;
    for (int i = 0; i < v.n; i++) vs[i] = expr_copy(vx[i]);
    Expr* hhead = expr_new_symbol(intern_symbol("Hyperedge"));
    int clash = 0;
    for (int j = 0; j < v.m; j++) {
        Expr* a = expr_new_integer(j + 1);
        Expr* node = expr_new_function(expr_copy(hhead), &a, 1);
        if (graph_vidx_get(v.ix, node) >= 0) clash = 1;
        vs[v.n + j] = node;
    }
    expr_free(hhead);
    if (clash) {
        for (size_t i = 0; i < nv; i++) expr_free(vs[i]);
        free(vs); free(es); expr_free(owned);
        return NULL;
    }
    Expr* head = expr_new_symbol(SYM_UndirectedEdge);
    size_t c = 0;
    for (int j = 0; j < v.m; j++)
        for (int k = v.soff[j]; k < v.soff[j + 1]; k++)
            es[c++] = edge2(head, expr_copy(vx[v.sv[k]]), expr_copy(vs[v.n + j]));
    expr_free(head);
    Expr* out = mk_graph(mk_list(vs, nv), mk_list(es, c));
    free(vs); free(es);
    expr_free(owned);
    return out;
}

/* Function Repository HypergraphToGraph: ordered hyperedge {v1, ..., vk}
 * contributes v_a -> v_b for every a < b (by position). Mathilda graphs are
 * simple, so self-loops (from a repeated vertex) are dropped and parallel
 * copies merged; every vertex of h is kept, including those only in unary
 * hyperedges (the FR function drops those). */
Expr* builtin_hypergraph_to_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned;
    const Expr* h = hyp_arg(res->data.function.args[0], &owned);
    HypView v;
    if (!h || !hyp_view(h, &v, 0)) { expr_free(owned); return NULL; }
    U64Set seen;
    if (!u64_init(&seen, (size_t)v.eoff[v.m])) { expr_free(owned); return NULL; }
    EVec ev = {0};
    Expr* head = expr_new_symbol(SYM_DirectedEdge);
    Expr* const* vx = v.verts->data.function.args;
    int ok = 1;
    for (int j = 0; j < v.m && ok; j++)
        for (int a = v.eoff[j]; a < v.eoff[j + 1] && ok; a++)
            for (int b = a + 1; b < v.eoff[j + 1] && ok; b++) {
                int x = v.ev[a], y = v.ev[b];
                if (x == y) continue;
                int r = u64_insert(&seen, pair_key(x, y));
                if (r < 0) ok = 0;
                else if (r == 1)
                    ok = evec_push(&ev, edge2(head, expr_copy(vx[x]), expr_copy(vx[y])));
            }
    expr_free(head);
    free(seen.k);
    if (!ok) { evec_free(&ev); expr_free(owned); return NULL; }
    Expr* out = mk_graph(expr_copy((Expr*)v.verts), evec_to_list(&ev));
    expr_free(owned);
    return out;
}

/* ---- s-overlap machinery -------------------------------------------------- *
 * For hyperedge i, the hyperedges j sharing at least s vertices with it: walk
 * i's distinct vertices and, through the incidence lists, count each touched
 * j. `cnt` must be all-zero on entry and is left all-zero. Writes the
 * qualifying j (> i when after_only, else != i) to out[], returns how many. */
static int s_neighbours(const HypView* v, int i, int s, int after_only,
                        int* cnt, int* touched, int* out) {
    int nt = 0, no = 0;
    for (int k = v->soff[i]; k < v->soff[i + 1]; k++) {
        int x = v->sv[k];
        for (int t = v->voff[x]; t < v->voff[x + 1]; t++) {
            int j = v->ve[t];
            if (j == i || (after_only && j < i)) continue;
            if (cnt[j]++ == 0) touched[nt++] = j;
        }
    }
    for (int t = 0; t < nt; t++) {
        int j = touched[t];
        if (cnt[j] >= s) out[no++] = j;
        cnt[j] = 0;
    }
    return no;
}

static int s_arg(const Expr* res, size_t at, int* s) {
    *s = 1;
    if (res->data.function.arg_count <= at) return 1;
    long x = pos_int(res->data.function.args[at]);
    if (x < 1 || x > INT32_MAX) return 0;
    *s = (int)x;
    return 1;
}

/* HypergraphLineGraph[h, s]: Graph on 1..m, i <-> j iff |e_i ∩ e_j| >= s. */
Expr* builtin_hypergraph_line_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    int s;
    if (!s_arg(res, 1, &s)) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 1)) return NULL;
    int mm = v.m > 0 ? v.m : 1;
    int* cnt = calloc((size_t)mm, sizeof(int));
    int* touched = malloc((size_t)mm * sizeof(int));
    int* nb = malloc((size_t)mm * sizeof(int));
    Expr** ints = malloc((size_t)mm * sizeof(Expr*));
    if (!cnt || !touched || !nb || !ints) {
        free(cnt); free(touched); free(nb); free(ints); return NULL;
    }
    for (int j = 0; j < v.m; j++) ints[j] = expr_new_integer(j + 1);
    EVec ev = {0};
    Expr* head = expr_new_symbol(SYM_UndirectedEdge);
    int ok = 1;
    for (int i = 0; i < v.m && ok; i++) {
        if (v.soff[i + 1] - v.soff[i] < s) continue;
        int no = s_neighbours(&v, i, s, 1, cnt, touched, nb);
        if (no > 1) qsort(nb, (size_t)no, sizeof(int), cmp_int);
        for (int t = 0; t < no && ok; t++)
            ok = evec_push(&ev, edge2(head, expr_copy(ints[i]), expr_copy(ints[nb[t]])));
    }
    expr_free(head);
    free(cnt); free(touched); free(nb);
    if (!ok) {
        evec_free(&ev);
        for (int j = 0; j < v.m; j++) expr_free(ints[j]);
        free(ints);
        return NULL;
    }
    Expr* vl = mk_list(ints, (size_t)v.m);
    free(ints);
    return mk_graph(vl, evec_to_list(&ev));
}

/* ---- components ----------------------------------------------------------- */

static int uf_find(int* p, int x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}
static void uf_union(int* p, int a, int b) {
    a = uf_find(p, a); b = uf_find(p, b);
    if (a != b) { if (a < b) p[b] = a; else p[a] = b; }
}

/* Group items 0..n-1 with active[i] by root, components ordered by their
 * smallest member, members ascending. items(i) renders member i. */
static Expr* groups_to_list(int* parent, int n, const unsigned char* active,
                            Expr* (*item)(const void*, int), const void* ctx) {
    int nn = n > 0 ? n : 1;
    int* id = malloc((size_t)nn * sizeof(int));
    int* size = calloc((size_t)nn, sizeof(int));
    if (!id || !size) { free(id); free(size); return NULL; }
    for (int i = 0; i < n; i++) id[i] = -1;
    int k = 0;
    int* cid = malloc((size_t)nn * sizeof(int));
    if (!cid) { free(id); free(size); return NULL; }
    for (int i = 0; i < n; i++) {
        if (active && !active[i]) { cid[i] = -1; continue; }
        int r = uf_find(parent, i);
        if (id[r] < 0) id[r] = k++;
        cid[i] = id[r];
        size[cid[i]]++;
    }
    Expr*** mem = malloc((size_t)(k > 0 ? k : 1) * sizeof(Expr**));
    int* fill = calloc((size_t)(k > 0 ? k : 1), sizeof(int));
    if (!mem || !fill) { free(id); free(size); free(cid); free(mem); free(fill); return NULL; }
    for (int c = 0; c < k; c++) mem[c] = malloc((size_t)size[c] * sizeof(Expr*));
    for (int i = 0; i < n; i++)
        if (cid[i] >= 0) mem[cid[i]][fill[cid[i]]++] = item(ctx, i);
    Expr** comps = malloc((size_t)(k > 0 ? k : 1) * sizeof(Expr*));
    for (int c = 0; c < k; c++) { comps[c] = mk_list(mem[c], (size_t)size[c]); free(mem[c]); }
    Expr* out = mk_list(comps, (size_t)k);
    free(comps); free(mem); free(fill); free(id); free(size); free(cid);
    return out;
}

static Expr* vertex_item(const void* ctx, int i) {
    return expr_copy(((const Expr*)ctx)->data.function.args[i]);
}
static Expr* index_item(const void* ctx, int i) {
    (void)ctx;
    return expr_new_integer(i + 1);
}

/* Union-find parent array over the vertices of v, unioned along hyperedges. */
static int* vertex_uf(const HypView* v) {
    int* p = malloc((size_t)(v->n > 0 ? v->n : 1) * sizeof(int));
    if (!p) return NULL;
    for (int i = 0; i < v->n; i++) p[i] = i;
    for (int j = 0; j < v->m; j++)
        for (int k = v->soff[j] + 1; k < v->soff[j + 1]; k++)
            uf_union(p, v->sv[v->soff[j]], v->sv[k]);
    return p;
}

Expr* builtin_hypergraph_connected_components(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned;
    const Expr* h = hyp_arg(res->data.function.args[0], &owned);
    HypView v;
    if (!h || !hyp_view(h, &v, 0)) { expr_free(owned); return NULL; }
    int* p = vertex_uf(&v);
    Expr* out = p ? groups_to_list(p, v.n, NULL, vertex_item, v.verts) : NULL;
    free(p);
    expr_free(owned);
    return out;
}

/* ConnectedHypergraphQ[h]: at least one vertex, and a single component. A
 * plain List of hyperedges is accepted, as the FR function takes. */
Expr* builtin_connected_hypergraph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned;
    const Expr* h = hyp_arg(res->data.function.args[0], &owned);
    HypView v;
    if (!h || !hyp_view(h, &v, 0)) { expr_free(owned); return expr_new_symbol(SYM_False); }
    int* p = vertex_uf(&v);
    if (!p) { expr_free(owned); return NULL; }
    int roots = 0;
    for (int i = 0; i < v.n; i++) if (p[i] == i) roots++;
    free(p);
    expr_free(owned);
    return expr_new_symbol(roots == 1 ? SYM_True : SYM_False);
}

/* HyperedgeConnectedComponents[h, s]: the s-connected components (Aksoy et
 * al.): hyperedges of at least s distinct vertices, joined when they share at
 * least s vertices, closed transitively. Lists of 1-based hyperedge indices.
 * Hyperedges with fewer than s vertices belong to no s-walk and are omitted. */
Expr* builtin_hyperedge_connected_components(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    int s;
    if (!s_arg(res, 1, &s)) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 1)) return NULL;
    int mm = v.m > 0 ? v.m : 1;
    int* p = malloc((size_t)mm * sizeof(int));
    unsigned char* act = malloc((size_t)mm);
    if (!p || !act) { free(p); free(act); return NULL; }
    for (int j = 0; j < v.m; j++) {
        p[j] = j;
        act[j] = (unsigned char)(v.soff[j + 1] - v.soff[j] >= s);
    }
    if (s == 1) {
        for (int x = 0; x < v.n; x++)
            for (int t = v.voff[x] + 1; t < v.voff[x + 1]; t++)
                uf_union(p, v.ve[v.voff[x]], v.ve[t]);
    } else {
        int* cnt = calloc((size_t)mm, sizeof(int));
        int* touched = malloc((size_t)mm * sizeof(int));
        int* nb = malloc((size_t)mm * sizeof(int));
        if (!cnt || !touched || !nb) {
            free(cnt); free(touched); free(nb); free(p); free(act); return NULL;
        }
        for (int i = 0; i < v.m; i++) {
            if (!act[i]) continue;
            int no = s_neighbours(&v, i, s, 1, cnt, touched, nb);
            for (int t = 0; t < no; t++) uf_union(p, i, nb[t]);
        }
        free(cnt); free(touched); free(nb);
    }
    Expr* out = groups_to_list(p, v.m, act, index_item, NULL);
    free(p); free(act);
    return out;
}

/* ---- distances ------------------------------------------------------------ */

static Expr* dist_expr(int d) {
    return d < 0 ? expr_new_symbol(SYM_Infinity) : expr_new_integer(d);
}

/* Distances as a List: packed ints when all finite, else Infinity entries. */
static Expr* dist_list(const int* d, int n) {
    int all = 1;
    for (int i = 0; i < n; i++) if (d[i] < 0) { all = 0; break; }
    if (all) {
        int64_t* b = malloc((size_t)(n > 0 ? n : 1) * sizeof(int64_t));
        if (!b) return NULL;
        for (int i = 0; i < n; i++) b[i] = d[i];
        Expr* out = hyp_int_list(b, (size_t)n);
        free(b);
        return out;
    }
    Expr** it = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    if (!it) return NULL;
    for (int i = 0; i < n; i++) it[i] = dist_expr(d[i]);
    Expr* out = mk_list(it, (size_t)n);
    free(it);
    return out;
}

/* BFS over hyperedges along s-adjacency from src; stops early at target >= 0.
 * dist[] (length m) receives the distances, -1 = unreachable. */
static int edge_bfs(const HypView* v, int src, int target, int s, int* dist) {
    int m = v->m, mm = m > 0 ? m : 1;
    int* q = malloc((size_t)mm * sizeof(int));
    if (!q) return 0;
    for (int j = 0; j < m; j++) dist[j] = -1;
    dist[src] = 0;
    int qh = 0, qt = 0;
    q[qt++] = src;
    if (src == target) { free(q); return 1; }
    if (s == 1) {
        unsigned char* vexp = calloc((size_t)(v->n > 0 ? v->n : 1), 1);
        if (!vexp) { free(q); return 0; }
        while (qh < qt) {
            int e = q[qh++];
            for (int k = v->soff[e]; k < v->soff[e + 1]; k++) {
                int x = v->sv[k];
                if (vexp[x]) continue;
                vexp[x] = 1;
                for (int t = v->voff[x]; t < v->voff[x + 1]; t++) {
                    int f = v->ve[t];
                    if (dist[f] >= 0) continue;
                    dist[f] = dist[e] + 1;
                    if (f == target) { free(vexp); free(q); return 1; }
                    q[qt++] = f;
                }
            }
        }
        free(vexp);
    } else if (v->soff[src + 1] - v->soff[src] >= s) {
        int* cnt = calloc((size_t)mm, sizeof(int));
        int* touched = malloc((size_t)mm * sizeof(int));
        int* nb = malloc((size_t)mm * sizeof(int));
        if (!cnt || !touched || !nb) { free(cnt); free(touched); free(nb); free(q); return 0; }
        while (qh < qt) {
            int e = q[qh++];
            int no = s_neighbours(v, e, s, 0, cnt, touched, nb);
            for (int t = 0; t < no; t++) {
                int f = nb[t];
                if (dist[f] >= 0) continue;
                dist[f] = dist[e] + 1;
                q[qt++] = f;
                if (f == target) { qh = qt; break; }
            }
        }
        free(cnt); free(touched); free(nb);
    }
    free(q);
    return 1;
}

/* HyperedgeDistance[h, i, j(, s)], HyperedgeDistance[h, i(, All, s)]. */
Expr* builtin_hyperedge_distance(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 4) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 1)) return NULL;
    long i = pos_int(res->data.function.args[1]);
    if (i < 1 || i > v.m) return NULL;
    long j = -1;       /* -1: all */
    if (argc >= 3) {
        const Expr* je = res->data.function.args[2];
        if (je->type == EXPR_SYMBOL && je->data.symbol.name == SYM_All) j = -1;
        else { j = pos_int(je); if (j < 1 || j > v.m) return NULL; }
    }
    int s;
    if (!s_arg(res, 3, &s)) return NULL;
    int* dist = malloc((size_t)(v.m > 0 ? v.m : 1) * sizeof(int));
    if (!dist) return NULL;
    if (!edge_bfs(&v, (int)i - 1, j > 0 ? (int)j - 1 : -1, s, dist)) { free(dist); return NULL; }
    Expr* out = j > 0 ? dist_expr(dist[j - 1]) : dist_list(dist, v.m);
    free(dist);
    return out;
}

/* HypergraphDistance[h, u, v] / [h, u]: least number of hyperedges in a chain
 * u = x0, x1, ..., xk = v with consecutive vertices sharing a hyperedge (the
 * distance in the clique expansion), Infinity if none. */
Expr* builtin_hypergraph_distance(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 1)) return NULL;
    int src = graph_vidx_get(v.ix, res->data.function.args[1]);
    if (src < 0) return NULL;
    int tgt = -1;
    if (argc == 3) {
        const Expr* te = res->data.function.args[2];
        tgt = graph_vidx_get(v.ix, te);
        if (tgt < 0) return NULL;
    }
    int nn = v.n > 0 ? v.n : 1;
    int* dist = malloc((size_t)nn * sizeof(int));
    int* q = malloc((size_t)nn * sizeof(int));
    unsigned char* eexp = calloc((size_t)(v.m > 0 ? v.m : 1), 1);
    if (!dist || !q || !eexp) { free(dist); free(q); free(eexp); return NULL; }
    for (int x = 0; x < v.n; x++) dist[x] = -1;
    dist[src] = 0;
    int qh = 0, qt = 0, done = (src == tgt);
    q[qt++] = src;
    while (qh < qt && !done) {
        int x = q[qh++];
        for (int t = v.voff[x]; t < v.voff[x + 1] && !done; t++) {
            int f = v.ve[t];
            if (eexp[f]) continue;
            eexp[f] = 1;
            for (int k = v.soff[f]; k < v.soff[f + 1]; k++) {
                int y = v.sv[k];
                if (dist[y] >= 0) continue;
                dist[y] = dist[x] + 1;
                q[qt++] = y;
                if (y == tgt) { done = 1; break; }
            }
        }
    }
    Expr* out = tgt >= 0 ? dist_expr(dist[tgt]) : dist_list(dist, v.n);
    free(dist); free(q); free(eexp);
    return out;
}

/* ---- edits ---------------------------------------------------------------- */

/* The items named by an edit argument: a List whose elements are all Lists is
 * a list of hyperedges (for the Edge heads); for the Vertex heads any List is a
 * list of vertices. Otherwise the argument is one item. Borrowed pointers. */
static const Expr* const* edit_items(const Expr* a, int edges, size_t* n, const Expr** single) {
    if (graph_is_list(a)) {
        int all = 1;
        if (edges)
            for (size_t k = 0; k < a->data.function.arg_count && all; k++)
                if (!graph_is_list(a->data.function.args[k])) all = 0;
        if (all) { *n = a->data.function.arg_count; return (const Expr* const*)a->data.function.args; }
    }
    *single = a;
    *n = 1;
    return single;
}

Expr* builtin_hypergraph_vertex_add(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    const Expr* single; size_t k;
    const Expr* const* xs = edit_items(res->data.function.args[1], 0, &k, &single);
    GraphVIdx* extra = graph_vidx_new(k);
    Expr** vs = malloc(((size_t)v.n + k + 1) * sizeof(Expr*));
    if (!extra || !vs) { graph_vidx_free(extra); free(vs); return NULL; }
    size_t c = 0;
    for (int i = 0; i < v.n; i++) vs[c++] = expr_copy(v.verts->data.function.args[i]);
    for (size_t t = 0; t < k; t++)
        if (graph_vidx_get(v.ix, xs[t]) < 0 && graph_vidx_put(extra, xs[t], 0))
            vs[c++] = expr_copy((Expr*)xs[t]);
    graph_vidx_free(extra);
    Expr* out = mk_hyp(mk_list(vs, c), expr_copy((Expr*)v.edges));
    free(vs);
    return out;
}

/* Keep the vertices with keep_v[i] and the hyperedges with keep_e[j]. */
static Expr* rebuild(const HypView* v, const unsigned char* keep_v, const unsigned char* keep_e) {
    Expr** vs = malloc((size_t)(v->n > 0 ? v->n : 1) * sizeof(Expr*));
    Expr** es = malloc((size_t)(v->m > 0 ? v->m : 1) * sizeof(Expr*));
    if (!vs || !es) { free(vs); free(es); return NULL; }
    size_t c = 0, d = 0;
    for (int i = 0; i < v->n; i++) if (keep_v[i]) vs[c++] = expr_copy(v->verts->data.function.args[i]);
    for (int j = 0; j < v->m; j++) if (keep_e[j]) es[d++] = expr_copy(v->edges->data.function.args[j]);
    Expr* out = mk_hyp(mk_list(vs, c), mk_list(es, d));
    free(vs); free(es);
    return out;
}

/* HypergraphVertexDelete[h, v | {v1, ...}]: removes the vertices and every
 * hyperedge incident to one of them, as VertexDelete does for a Graph.
 * Unevaluated if a named vertex is not in h. */
Expr* builtin_hypergraph_vertex_delete(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    const Expr* single; size_t k;
    const Expr* const* xs = edit_items(res->data.function.args[1], 0, &k, &single);
    unsigned char* kv = malloc((size_t)(v.n > 0 ? v.n : 1));
    unsigned char* ke = malloc((size_t)(v.m > 0 ? v.m : 1));
    if (!kv || !ke) { free(kv); free(ke); return NULL; }
    memset(kv, 1, (size_t)v.n);
    for (size_t t = 0; t < k; t++) {
        int i = graph_vidx_get(v.ix, xs[t]);
        if (i < 0) { free(kv); free(ke); return NULL; }
        kv[i] = 0;
    }
    for (int j = 0; j < v.m; j++) {
        ke[j] = 1;
        for (int q = v.soff[j]; q < v.soff[j + 1]; q++) if (!kv[v.sv[q]]) { ke[j] = 0; break; }
    }
    Expr* out = rebuild(&v, kv, ke);
    free(kv); free(ke);
    return out;
}

/* HypergraphEdgeAdd[h, e | {e1, ...}]: appends the hyperedges (repeats
 * allowed), adding any new vertices in first-appearance order. */
Expr* builtin_hypergraph_edge_add(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    const Expr* single; size_t k;
    const Expr* const* xs = edit_items(res->data.function.args[1], 1, &k, &single);
    size_t extra_n = 0;
    for (size_t t = 0; t < k; t++) {
        if (!graph_is_list(xs[t])) return NULL;
        extra_n += xs[t]->data.function.arg_count;
    }
    GraphVIdx* extra = graph_vidx_new(extra_n);
    Expr** vs = malloc(((size_t)v.n + extra_n + 1) * sizeof(Expr*));
    Expr** es = malloc(((size_t)v.m + k + 1) * sizeof(Expr*));
    if (!extra || !vs || !es) { graph_vidx_free(extra); free(vs); free(es); return NULL; }
    size_t c = 0;
    for (int i = 0; i < v.n; i++) vs[c++] = expr_copy(v.verts->data.function.args[i]);
    for (size_t t = 0; t < k; t++)
        for (size_t q = 0; q < xs[t]->data.function.arg_count; q++) {
            Expr* x = xs[t]->data.function.args[q];
            if (graph_vidx_get(v.ix, x) < 0 && graph_vidx_put(extra, x, 0)) vs[c++] = expr_copy(x);
        }
    graph_vidx_free(extra);
    size_t d = 0;
    for (int j = 0; j < v.m; j++) es[d++] = expr_copy(v.edges->data.function.args[j]);
    for (size_t t = 0; t < k; t++) es[d++] = expr_copy((Expr*)xs[t]);
    Expr* out = mk_hyp(mk_list(vs, c), mk_list(es, d));
    free(vs); free(es);
    return out;
}

/* HypergraphEdgeDelete[h, e | {e1, ...}]: removes every hyperedge SameQ to a
 * named one (the List as written, so {2, 1} does not delete {1, 2}).
 * Unevaluated if a named hyperedge does not occur in h. */
Expr* builtin_hypergraph_edge_delete(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    const Expr* single; size_t k;
    const Expr* const* xs = edit_items(res->data.function.args[1], 1, &k, &single);
    GraphVIdx* del = graph_vidx_new(k);
    unsigned char* hit = calloc(k ? k : 1, 1);
    unsigned char* kv = malloc((size_t)(v.n > 0 ? v.n : 1));
    unsigned char* ke = malloc((size_t)(v.m > 0 ? v.m : 1));
    if (!del || !hit || !kv || !ke) { graph_vidx_free(del); free(hit); free(kv); free(ke); return NULL; }
    for (size_t t = 0; t < k; t++) graph_vidx_put(del, xs[t], (int)t);
    memset(kv, 1, (size_t)v.n);
    for (int j = 0; j < v.m; j++) {
        int t = graph_vidx_get(del, v.edges->data.function.args[j]);
        ke[j] = (unsigned char)(t < 0);
        if (t >= 0) hit[t] = 1;
    }
    int ok = 1;
    for (size_t t = 0; t < k; t++) {
        /* A repeated name maps to its first slot; that slot's hit covers it. */
        int first = graph_vidx_get(del, xs[t]);
        if (!hit[first]) ok = 0;
    }
    Expr* out = ok ? rebuild(&v, kv, ke) : NULL;
    graph_vidx_free(del); free(hit); free(kv); free(ke);
    return out;
}

/* Membership mask for a vertex List argument (non-vertices ignored). */
static unsigned char* vertex_mask(const HypView* v, const Expr* s) {
    if (!graph_is_list(s)) return NULL;
    unsigned char* in = calloc((size_t)(v->n > 0 ? v->n : 1), 1);
    if (!in) return NULL;
    for (size_t t = 0; t < s->data.function.arg_count; t++) {
        int i = graph_vidx_get(v->ix, s->data.function.args[t]);
        if (i >= 0) in[i] = 1;
    }
    return in;
}

/* Subhypergraph[h, {v1, ...}]: the vertices of h among the v_i (VertexList
 * order) and the hyperedges lying entirely inside them -- Subgraph's rule. */
Expr* builtin_subhypergraph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    unsigned char* in = vertex_mask(&v, res->data.function.args[1]);
    unsigned char* ke = malloc((size_t)(v.m > 0 ? v.m : 1));
    if (!in || !ke) { free(in); free(ke); return NULL; }
    for (int j = 0; j < v.m; j++) {
        ke[j] = 1;
        for (int q = v.soff[j]; q < v.soff[j + 1]; q++) if (!in[v.sv[q]]) { ke[j] = 0; break; }
    }
    Expr* out = rebuild(&v, in, ke);
    free(in); free(ke);
    return out;
}

/* HypergraphRestriction[h, {v1, ...}]: Berge's induced sub-hypergraph -- every
 * hyperedge intersected with the vertex set (order and repeats of the kept
 * vertices preserved), hyperedges that miss it entirely dropped. */
Expr* builtin_hypergraph_restriction(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    HypView v;
    if (!hyp_view(res->data.function.args[0], &v, 0)) return NULL;
    unsigned char* in = vertex_mask(&v, res->data.function.args[1]);
    Expr** vs = malloc((size_t)(v.n > 0 ? v.n : 1) * sizeof(Expr*));
    Expr** es = malloc((size_t)(v.m > 0 ? v.m : 1) * sizeof(Expr*));
    Expr** tmp = malloc((size_t)(v.eoff[v.m] > 0 ? v.eoff[v.m] : 1) * sizeof(Expr*));
    if (!in || !vs || !es || !tmp) { free(in); free(vs); free(es); free(tmp); return NULL; }
    size_t c = 0, d = 0;
    for (int i = 0; i < v.n; i++) if (in[i]) vs[c++] = expr_copy(v.verts->data.function.args[i]);
    Expr* lh = expr_new_symbol(SYM_List);
    for (int j = 0; j < v.m; j++) {
        const Expr* e = v.edges->data.function.args[j];
        int all = 1, t = 0;
        for (int q = v.eoff[j]; q < v.eoff[j + 1]; q++) {
            if (in[v.ev[q]]) tmp[t++] = expr_copy(e->data.function.args[q - v.eoff[j]]);
            else all = 0;
        }
        if (t == 0 && v.eoff[j + 1] > v.eoff[j]) continue;       /* missed S   */
        if (t == 0) continue;                                     /* empty edge */
        if (all) {
            for (int q = 0; q < t; q++) expr_free(tmp[q]);
            es[d++] = expr_copy((Expr*)e);                        /* share it   */
        } else {
            es[d++] = expr_new_function(expr_copy(lh), tmp, (size_t)t);
        }
    }
    expr_free(lh);
    Expr* out = mk_hyp(mk_list(vs, c), mk_list(es, d));
    free(in); free(vs); free(es); free(tmp);
    return out;
}
