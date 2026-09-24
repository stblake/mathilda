/* gops_edit.c - graph editing builtins.
 *
 *   VertexAdd[g, v | {v..}]         append vertices not already present (in
 *                                   order, repeats ignored)
 *   VertexDelete[g, v | {v..} | p]  remove vertices and their incident edges;
 *                                   every listed vertex must exist; a pattern
 *                                   p removes every vertex matching it
 *   EdgeAdd[g, e | {e..}]           append edges (->/<-> sugar accepted); new
 *                                   endpoints become vertices (appended)
 *   EdgeDelete[g, e | {e..} | p]    remove edges; each must exist (an undirected
 *                                   edge matches either orientation)
 *   Subgraph[g, {v..} | p]          subgraph induced by the vertices
 *   NeighborhoodGraph[g, v | {v..}, k]   subgraph induced by the vertices within
 *                                   distance k (default 1, Infinity allowed),
 *                                   edge direction ignored
 *   VertexReplace[g, rules]         rename vertices by Replace
 *   EdgeRules[g]                    edges as u -> v rules
 *   VertexIndex[g, v | {v..}]       1-based VertexList position(s)
 *   EdgeIndex[g, e | {e..}]         1-based EdgeList position(s)
 *   IndexGraph[g] / IndexGraph[g, r]   vertices renamed r, r+1, ... (r = 1)
 *
 * Every edit is an integer pass over the memo's endpoint arrays (gops_view):
 * O(V + E) plus one hash per argument item, with vertex/edge/weight nodes
 * shared into the result. Weighted graphs keep their EdgeWeight list aligned
 * with the surviving edges; EdgeAdd gives a new edge weight 1 (as Mathematica).
 *
 * Ordering follows Mathematica 15: VertexDelete/EdgeDelete keep the original
 * orders; Subgraph and NeighborhoodGraph list vertices in the given (resp.
 * centres-then-balls) order and order edges by (later endpoint, earlier
 * endpoint) positions in that list -- the lower-triangular adjacency order --
 * with the edge out of the later vertex first on a tie.
 *
 * Deviations (Mathilda graphs are simple): an edit that would create a
 * self-loop or a parallel edge -- EdgeAdd of an existing edge, VertexReplace
 * merging two adjacent vertices -- is left unevaluated rather than returning a
 * multigraph. NeighborhoodGraph also accepts k = Infinity (Mathematica leaves
 * that unevaluated).
 *
 * Memory (SPEC section 4): results are fresh; res is borrowed, never modified.
 */

#include "graph_ops.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

/* Copy of the view's endpoint arrays for a result that keeps every edge. */
static int copy_endpoints(const GopsView* v, int** eu, int** ev, unsigned char** ed) {
    *eu = gops_malloc(v->ne, sizeof(int));
    *ev = gops_malloc(v->ne, sizeof(int));
    *ed = gops_malloc(v->ne, 1);
    if (!*eu || !*ev || !*ed) { free(*eu); free(*ev); free(*ed); return 0; }
    if (v->ne) {
        memcpy(*eu, v->eu, v->ne * sizeof(int));
        memcpy(*ev, v->ev, v->ne * sizeof(int));
        memcpy(*ed, v->edir, v->ne);
    }
    return 1;
}

static int is_rule_sugar(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rule;
}

static Expr** copy_exprs(Expr* const* src, size_t n, size_t extra) {
    Expr** a = gops_calloc(n + extra, sizeof(Expr*));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(src[i]);
    return a;
}

/* ---- VertexAdd ------------------------------------------------------------ */
Expr* builtin_vertex_add(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    size_t ni;
    Expr* const* items = gops_items(&res->data.function.args[1], &ni);

    /* Which items are new, deduplicated among themselves. */
    Expr** add = gops_malloc(ni, sizeof(Expr*));
    GraphVIdx* seen = graph_vidx_new(ni);
    if (!add || !seen) { free(add); graph_vidx_free(seen); return NULL; }
    size_t na = 0;
    for (size_t i = 0; i < ni; i++) {
        if (graph_vertex_position(v.g, items[i]) >= 0) continue;
        if (graph_vidx_put(seen, items[i], (int)na)) add[na++] = items[i];
    }
    graph_vidx_free(seen);
    if (na == 0) { free(add); return expr_copy(res->data.function.args[0]); }

    /* gops_view's arrays are still valid: graph_vertex_position hit the memo. */
    Expr** verts = copy_exprs(v.verts, v.nv, na);
    Expr** edges = copy_exprs(v.edges, v.ne, 0);
    Expr** ws = v.weights ? copy_exprs(v.weights, v.ne, 0) : NULL;
    int *eu, *ev; unsigned char* ed;
    if (!verts || !edges || (v.weights && !ws) || !copy_endpoints(&v, &eu, &ev, &ed)) {
        /* cannot fail in practice; release what was taken */
        if (verts) { for (size_t i = 0; i < v.nv; i++) expr_free(verts[i]); free(verts); }
        if (edges) { for (size_t i = 0; i < v.ne; i++) expr_free(edges[i]); free(edges); }
        if (ws) { for (size_t i = 0; i < v.ne; i++) expr_free(ws[i]); free(ws); }
        free(add);
        return NULL;
    }
    for (size_t i = 0; i < na; i++) verts[v.nv + i] = expr_copy(add[i]);
    free(add);
    return gops_graph_new(verts, v.nv + na, edges, v.ne, ws, eu, ev, ed);
}

/* ---- Keep-mask result builders -------------------------------------------- */

/* The graph on the vertices with vkeep[i] (order kept) and the edges with
 * ekeep[k] (order kept; both endpoints must be kept). */
static Expr* build_masked(const GopsView* v, const unsigned char* vkeep,
                          const unsigned char* ekeep) {
    int* nid = gops_malloc(v->nv, sizeof(int));
    if (!nid) return NULL;
    size_t nv2 = 0, ne2 = 0;
    for (size_t i = 0; i < v->nv; i++) nid[i] = vkeep[i] ? (int)nv2++ : -1;
    for (size_t k = 0; k < v->ne; k++) ne2 += ekeep[k] ? 1 : 0;

    Expr** verts = gops_malloc(nv2, sizeof(Expr*));
    Expr** edges = gops_malloc(ne2, sizeof(Expr*));
    Expr** ws = v->weights ? gops_malloc(ne2, sizeof(Expr*)) : NULL;
    int* eu = gops_malloc(ne2, sizeof(int));
    int* ev = gops_malloc(ne2, sizeof(int));
    unsigned char* ed = gops_malloc(ne2, 1);
    if (!verts || !edges || (v->weights && !ws) || !eu || !ev || !ed) {
        free(nid); free(verts); free(edges); free(ws); free(eu); free(ev); free(ed);
        return NULL;
    }
    size_t j = 0;
    for (size_t i = 0; i < v->nv; i++) if (vkeep[i]) verts[j++] = expr_copy(v->verts[i]);
    j = 0;
    for (size_t k = 0; k < v->ne; k++) {
        if (!ekeep[k]) continue;
        edges[j] = expr_copy(v->edges[k]);
        if (ws) ws[j] = expr_copy(v->weights[k]);
        eu[j] = nid[v->eu[k]]; ev[j] = nid[v->ev[k]]; ed[j] = v->edir[k];
        j++;
    }
    free(nid);
    return gops_graph_new(verts, nv2, edges, ne2, ws, eu, ev, ed);
}

/* ---- VertexDelete --------------------------------------------------------- */
Expr* builtin_vertex_delete(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* arg = res->data.function.args[1];
    GopsView v;
    if (!graph_is_valid(g)) return NULL;

    size_t nv = g->data.function.args[0]->data.function.arg_count;
    unsigned char* vkeep = gops_malloc(nv, 1);
    if (!vkeep) return NULL;
    memset(vkeep, 1, nv > 0 ? nv : 1);

    int p = graph_vertex_position(g, arg);
    if (p >= 0) {
        vkeep[p] = 0;
    } else if (graph_is_list(arg)) {
        /* Per item, as Mathematica: a literal vertex, or a pattern matched
         * against every vertex -- VertexDelete[g, {1, _?(# > 3 &)}]. */
        Expr* const* vs = g->data.function.args[0]->data.function.args;
        for (size_t i = 0; i < arg->data.function.arg_count; i++) {
            Expr* item = arg->data.function.args[i];
            int q = graph_vertex_position(g, item);
            if (q >= 0) { vkeep[q] = 0; continue; }
            if (!gops_has_pattern(item)) { free(vkeep); return NULL; }   /* not a vertex */
            for (size_t k = 0; k < nv; k++)
                if (gops_matchq(vs[k], item)) vkeep[k] = 0;
        }
    } else if (gops_has_pattern(arg)) {
        Expr* const* vs = g->data.function.args[0]->data.function.args;
        for (size_t i = 0; i < nv; i++)
            if (gops_matchq(vs[i], arg)) vkeep[i] = 0;
    } else {
        free(vkeep);
        return NULL;
    }

    /* Pattern matching may evaluate arbitrary code, so take the view last. */
    if (!gops_view(g, &v)) { free(vkeep); return NULL; }
    unsigned char* ekeep = gops_malloc(v.ne, 1);
    if (!ekeep) { free(vkeep); return NULL; }
    for (size_t k = 0; k < v.ne; k++) ekeep[k] = vkeep[v.eu[k]] && vkeep[v.ev[k]];
    Expr* out = build_masked(&v, vkeep, ekeep);
    free(vkeep); free(ekeep);
    return out;
}

/* ---- EdgeAdd -------------------------------------------------------------- */
Expr* builtin_edge_add(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* arg = res->data.function.args[1];
    const Expr *u0, *v0;
    size_t ni;
    Expr* const* items;
    if (gops_parse_edge(arg, &u0, &v0) >= 0) { items = &res->data.function.args[1]; ni = 1; }
    else if (graph_is_list(arg)) items = gops_items(&res->data.function.args[1], &ni);
    else return NULL;
    for (size_t i = 0; i < ni; i++)
        if (gops_parse_edge(items[i], &u0, &v0) < 0) return NULL;

    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    size_t nv = v.nv, ne = v.ne;

    /* New vertices, in first-appearance order; `fresh` indexes them. */
    Expr** verts = copy_exprs(v.verts, nv, 2 * ni);
    Expr** edges = copy_exprs(v.edges, ne, ni);
    Expr** ws = v.weights ? copy_exprs(v.weights, ne, ni) : NULL;
    int* eu = gops_malloc(ne + ni, sizeof(int));
    int* ev = gops_malloc(ne + ni, sizeof(int));
    unsigned char* ed = gops_malloc(ne + ni, 1);
    GraphVIdx* fresh = graph_vidx_new(2 * ni);
    GopsHeads hs = { { NULL, NULL } };
    if (!verts || !edges || (v.weights && !ws) || !eu || !ev || !ed || !fresh) goto oom;
    if (ne) {
        memcpy(eu, v.eu, ne * sizeof(int));
        memcpy(ev, v.ev, ne * sizeof(int));
        memcpy(ed, v.edir, ne);
    }
    size_t nv2 = nv;
    for (size_t i = 0; i < ni; i++) {
        const Expr *a, *b;
        int d = gops_parse_edge(items[i], &a, &b);
        /* u -> v sugar takes the graph's kind: undirected in an undirected
         * (or edgeless) graph, as Mathematica. DirectedEdge stays directed. */
        if (v.ndir == 0 && is_rule_sugar(items[i])) d = 0;
        const Expr* ends[2] = { a, b };
        int idx[2];
        for (int t = 0; t < 2; t++) {
            int q = graph_vertex_position(v.g, ends[t]);
            if (q < 0) {
                q = graph_vidx_get(fresh, ends[t]);
                if (q < 0) {
                    q = (int)nv2;
                    verts[nv2++] = expr_copy((Expr*)ends[t]);
                    graph_vidx_put(fresh, verts[q], q);
                }
            }
            idx[t] = q;
        }
        /* Share an already-canonical edge node; rebuild sugar. */
        const char* hk = graph_edge_kind(items[i]);
        edges[ne + i] = hk ? expr_copy(items[i])
                           : gops_edge(&hs, d, expr_copy(verts[idx[0]]), expr_copy(verts[idx[1]]));
        if (ws) ws[ne + i] = expr_new_integer(1);
        eu[ne + i] = idx[0]; ev[ne + i] = idx[1]; ed[ne + i] = (unsigned char)d;
    }
    graph_vidx_free(fresh);
    gops_heads_free(&hs);
    /* Self-loops and parallel edges are rejected by the memo seed (NULL). */
    return gops_graph_new(verts, nv2, edges, ne + ni, ws, eu, ev, ed);

oom:
    graph_vidx_free(fresh);
    if (verts) { for (size_t i = 0; i < nv; i++) expr_free(verts[i]); free(verts); }
    if (edges) { for (size_t i = 0; i < ne; i++) expr_free(edges[i]); free(edges); }
    if (ws) { for (size_t i = 0; i < ne; i++) expr_free(ws[i]); free(ws); }
    free(eu); free(ev); free(ed);
    return NULL;
}

/* ---- Edge lookup ---------------------------------------------------------- *
 * The edges named in an argument go into a small GopsKeySet, probed once per
 * graph edge: the table holds only the ARGUMENT's edges, so a pass over a
 * 10^5-edge graph costs one probe per edge into a tiny, cache-resident table. */
/* Resolve each edge item to its EdgeList index, writing idx[i] (-1 if the item
 * is not an edge of g). Returns 0 if some item is not an edge expression at
 * all, or on OOM. */
static int resolve_edges(const GopsView* v, Expr* const* items, size_t ni, int* idx,
                         int rule_follows_graph) {
    GopsKeySet s;
    if (!gops_keyset_init(&s, ni)) return 0;
    int* slot_of = gops_malloc(ni, sizeof(int));
    if (!slot_of) { gops_keyset_free(&s); return 0; }
    for (size_t i = 0; i < ni; i++) {
        const Expr *a, *b;
        int d = gops_parse_edge(items[i], &a, &b);
        if (d < 0) { free(slot_of); gops_keyset_free(&s); return 0; }
        if (rule_follows_graph && v->ndir == 0 && is_rule_sugar(items[i])) d = 0;
        idx[i] = -1;
        slot_of[i] = -1;
        int ia = graph_vertex_position(v->g, a), ib = graph_vertex_position(v->g, b);
        if (ia < 0 || ib < 0) continue;
        slot_of[i] = gops_keyset_put(&s, gops_edge_key(ia, ib, d), (int)i);
    }
    /* One pass over the graph: the first item naming each key gets the index;
     * repeats of a key share it via slot_of. */
    int* found = gops_malloc(ni, sizeof(int));
    if (!found) { free(slot_of); gops_keyset_free(&s); return 0; }
    for (size_t i = 0; i < ni; i++) found[i] = -1;
    for (size_t k = 0; k < v->ne; k++) {
        int it = gops_keyset_get(&s, gops_edge_key(v->eu[k], v->ev[k], v->edir[k]));
        if (it >= 0) found[it] = (int)k;
    }
    for (size_t i = 0; i < ni; i++) idx[i] = slot_of[i] >= 0 ? found[slot_of[i]] : -1;
    free(found); free(slot_of); gops_keyset_free(&s);
    return 1;
}

/* ---- EdgeDelete ----------------------------------------------------------- */
Expr* builtin_edge_delete(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* arg = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;
    const Expr *a, *b;
    size_t ni = 0;
    Expr* const* items = NULL;
    if (graph_is_list(arg)) items = gops_items(&res->data.function.args[1], &ni);
    else if (gops_has_pattern(arg) || gops_parse_edge(arg, &a, &b) >= 0) {
        items = &res->data.function.args[1]; ni = 1;
    } else return NULL;

    /* Classify per item, as Mathematica does for EdgeDelete[g, {1<->2, _[3,4]}]:
     * a pattern is matched against every edge; anything else must be a literal
     * edge of g. (Gating on the whole list containing a pattern used to match
     * each edge against the List itself, silently deleting nothing.) */
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    unsigned char* ekeep = gops_malloc(ne, 1);
    Expr** lits = gops_malloc(ni, sizeof(Expr*));
    if (!ekeep || !lits) { free(ekeep); free(lits); return NULL; }
    memset(ekeep, 1, ne > 0 ? ne : 1);
    size_t nl = 0;
    Expr* const* es = g->data.function.args[1]->data.function.args;
    for (size_t i = 0; i < ni; i++) {
        if (gops_has_pattern(items[i])) {
            for (size_t k = 0; k < ne; k++)
                if (gops_matchq(es[k], items[i])) ekeep[k] = 0;
        } else {
            lits[nl++] = items[i];
        }
    }

    /* Pattern matching may evaluate arbitrary code, so take the view last. */
    GopsView v;
    if (!gops_view(g, &v)) { free(lits); free(ekeep); return NULL; }
    if (nl > 0) {
        int* idx = gops_malloc(nl, sizeof(int));
        if (!idx || !resolve_edges(&v, lits, nl, idx, 0)) {
            free(idx); free(lits); free(ekeep); return NULL;
        }
        for (size_t i = 0; i < nl; i++) {
            if (idx[i] < 0) { free(idx); free(lits); free(ekeep); return NULL; }   /* not an edge */
            ekeep[idx[i]] = 0;
        }
        free(idx);
        if (!gops_view(g, &v)) { free(lits); free(ekeep); return NULL; }
    }
    free(lits);
    unsigned char* vkeep = gops_malloc(v.nv, 1);
    if (!vkeep) { free(ekeep); return NULL; }
    memset(vkeep, 1, v.nv > 0 ? v.nv : 1);
    Expr* out = build_masked(&v, vkeep, ekeep);
    free(vkeep); free(ekeep);
    return out;
}

/* ---- Induced subgraph in a given vertex order ----------------------------- *
 * pos[i] >= 0 is vertex i's position in the result (k vertices in all). Edge
 * order is Mathematica's: for each result vertex u in turn, u's incident edges
 * to EARLIER result vertices, in u's incidence order -- directed out-edges and
 * undirected edges in EdgeList order, then directed in-edges in EdgeList order.
 * Each edge is emitted once, at its later endpoint. O(V + E). */
static Expr* build_induced(const GopsView* v, const int* pos, size_t k) {
    size_t nv = v->nv, ne2 = 0;
    int* start = gops_calloc(nv + 1, sizeof(int));
    int* nout = gops_calloc(nv, sizeof(int));
    int* inv = gops_malloc(k, sizeof(int));        /* result position -> vertex */
    if (!start || !nout || !inv) { free(start); free(nout); free(inv); return NULL; }
    for (size_t e = 0; e < v->ne; e++) {
        int a = v->eu[e], b = v->ev[e];
        if (pos[a] < 0 || pos[b] < 0) continue;
        ne2++;
        start[a + 1]++; start[b + 1]++;
        nout[a]++;
        if (!v->edir[e]) nout[b]++;
    }
    for (size_t i = 0; i < nv; i++) start[i + 1] += start[i];
    int* inc = gops_malloc(2 * ne2, sizeof(int));  /* edge ids */
    int* fo = gops_malloc(nv, sizeof(int));
    int* fi = gops_malloc(nv, sizeof(int));
    Expr** verts = gops_malloc(k, sizeof(Expr*));
    Expr** edges = gops_malloc(ne2, sizeof(Expr*));
    Expr** ws = v->weights ? gops_malloc(ne2, sizeof(Expr*)) : NULL;
    int* eu = gops_malloc(ne2, sizeof(int));
    int* ev = gops_malloc(ne2, sizeof(int));
    unsigned char* ed = gops_malloc(ne2, 1);
    Expr* out = NULL;
    if (!inc || !fo || !fi || !verts || !edges || (v->weights && !ws) || !eu || !ev || !ed)
        goto done;
    for (size_t i = 0; i < nv; i++) { fo[i] = start[i]; fi[i] = start[i] + nout[i]; }
    for (size_t e = 0; e < v->ne; e++) {
        int a = v->eu[e], b = v->ev[e];
        if (pos[a] < 0 || pos[b] < 0) continue;
        inc[fo[a]++] = (int)e;
        if (v->edir[e]) inc[fi[b]++] = (int)e;
        else inc[fo[b]++] = (int)e;
    }
    for (size_t i = 0; i < nv; i++)
        if (pos[i] >= 0) { inv[pos[i]] = (int)i; verts[pos[i]] = expr_copy(v->verts[i]); }
    size_t j = 0;
    for (size_t r = 0; r < k; r++) {
        int u = inv[r];
        for (int t = start[u]; t < start[u + 1]; t++) {
            int e = inc[t];
            int w = v->eu[e] == u ? v->ev[e] : v->eu[e];
            if (pos[w] >= (int)r) continue;            /* emitted at w instead */
            edges[j] = expr_copy(v->edges[e]);
            if (ws) ws[j] = expr_copy(v->weights[e]);
            eu[j] = pos[v->eu[e]]; ev[j] = pos[v->ev[e]]; ed[j] = v->edir[e];
            j++;
        }
    }
    out = gops_graph_new(verts, k, edges, ne2, ws, eu, ev, ed);
    verts = NULL; edges = NULL; ws = NULL; eu = NULL; ev = NULL; ed = NULL;
done:
    free(start); free(nout); free(inv); free(inc); free(fo); free(fi);
    free(verts); free(edges); free(ws); free(eu); free(ev); free(ed);
    return out;
}

/* ---- Subgraph ------------------------------------------------------------- */
Expr* builtin_subgraph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* arg = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;
    size_t nv = g->data.function.args[0]->data.function.arg_count;
    int* pos = gops_malloc(nv, sizeof(int));
    if (!pos) return NULL;
    for (size_t i = 0; i < nv; i++) pos[i] = -1;
    size_t k = 0;
    if (graph_is_list(arg)) {
        for (size_t i = 0; i < arg->data.function.arg_count; i++) {
            const Expr* it = arg->data.function.args[i];
            int q = graph_vertex_position(g, it);
            if (q < 0) {
                const Expr *a, *b;
                /* An edge list selects an edge-induced subgraph: not handled. */
                if (gops_parse_edge(it, &a, &b) >= 0) { free(pos); return NULL; }
                continue;                              /* non-vertex: ignored */
            }
            if (pos[q] < 0) pos[q] = (int)k++;
        }
    } else if (gops_has_pattern(arg)) {
        Expr* const* vs = g->data.function.args[0]->data.function.args;
        for (size_t i = 0; i < nv; i++)
            if (gops_matchq(vs[i], arg)) pos[i] = (int)k++;
    } else {
        int q = graph_vertex_position(g, arg);
        if (q < 0) { free(pos); return NULL; }
        pos[q] = (int)k++;
    }
    GopsView v;
    Expr* out = gops_view(g, &v) ? build_induced(&v, pos, k) : NULL;
    free(pos);
    return out;
}

/* ---- NeighborhoodGraph ---------------------------------------------------- */
static int cmp_int(const void* a, const void* b) {
    int x = *(const int*)a, y = *(const int*)b;
    return (x > y) - (x < y);
}

Expr* builtin_neighborhood_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* carg = res->data.function.args[1];
    long depth = 1;                               /* -1 = unbounded */
    if (argc == 3) {
        const Expr* ka = res->data.function.args[2];
        if (gops_is_infinity(ka)) depth = -1;
        else if (ka->type == EXPR_INTEGER && ka->data.integer >= 0) depth = (long)ka->data.integer;
        else return NULL;
    }
    GopsView v;
    if (!gops_view(g, &v)) return NULL;
    size_t nv = v.nv;

    /* Centres: a vertex, or a list of them; non-vertices are ignored. */
    size_t nc = 0;
    int* centres = NULL;
    int q = graph_vertex_position(g, carg);
    if (q >= 0) {
        centres = gops_malloc(1, sizeof(int));
        if (!centres) return NULL;
        centres[nc++] = q;
    } else {
        size_t ni = 0;
        Expr* const* items = graph_is_list(carg) ? carg->data.function.args : NULL;
        if (items) ni = carg->data.function.arg_count;
        centres = gops_malloc(ni, sizeof(int));
        if (!centres) return NULL;
        for (size_t i = 0; i < ni; i++) {
            int p = graph_vertex_position(g, items[i]);
            if (p >= 0) centres[nc++] = p;
        }
    }
    if (!gops_view(g, &v)) { free(centres); return NULL; }

    int* pos = gops_malloc(nv, sizeof(int));
    int* mark = gops_calloc(nv, sizeof(int));       /* BFS stamp per centre   */
    int* dist = gops_malloc(nv, sizeof(int));
    int* queue = gops_malloc(nv, sizeof(int));
    int* ball = gops_malloc(nv, sizeof(int));
    GopsInc inc;
    int have_inc = 0;
    Expr* out = NULL;
    if (!pos || !mark || !dist || !queue || !ball) goto done;
    for (size_t i = 0; i < nv; i++) pos[i] = -1;
    size_t k = 0;
    for (size_t c = 0; c < nc; c++)
        if (pos[centres[c]] < 0) pos[centres[c]] = (int)k++;
    if (depth != 0 && nc > 0) {
        if (!gops_inc_build(&v, GOPS_INC_ALL, &inc)) goto done;
        have_inc = 1;
        for (size_t c = 0; c < nc; c++) {
            int s = centres[c], stamp = (int)c + 1;
            size_t head = 0, tail = 0, nb = 0;
            queue[tail++] = s; mark[s] = stamp; dist[s] = 0;
            while (head < tail) {
                int x = queue[head++];
                if (depth >= 0 && dist[x] >= depth) continue;
                for (int j = inc.start[x]; j < inc.start[x + 1]; j++) {
                    int y = inc.nbr[j];
                    if (mark[y] == stamp) continue;
                    mark[y] = stamp; dist[y] = dist[x] + 1;
                    queue[tail++] = y;
                    if (pos[y] < 0) ball[nb++] = y;
                }
            }
            /* This centre's new vertices, in VertexList order. */
            if (nb * 8 > nv) {
                size_t w = 0;
                for (size_t i = 0; i < nv; i++)
                    if (mark[i] == stamp && pos[i] < 0) ball[w++] = (int)i;
                nb = w;
            } else {
                qsort(ball, nb, sizeof(int), cmp_int);
            }
            for (size_t i = 0; i < nb; i++) pos[ball[i]] = (int)k++;
        }
    }
    if (!gops_view(g, &v)) goto done;
    out = build_induced(&v, pos, k);
done:
    if (have_inc) gops_inc_free(&inc);
    free(centres); free(pos); free(mark); free(dist); free(queue); free(ball);
    return out;
}

/* ---- VertexReplace -------------------------------------------------------- */
static int is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
            || e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

Expr* builtin_vertex_replace(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* rules = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;
    if (graph_is_list(rules)) {
        for (size_t i = 0; i < rules->data.function.arg_count; i++)
            if (!is_rule(rules->data.function.args[i])) return NULL;
    } else if (!is_rule(rules)) {
        return NULL;
    }

    /* Replace[VertexList, rules, {1}] -- one evaluation, Replace's semantics
     * (first matching rule, patterns, RuleDelayed) for free. */
    Expr* rargs[3];
    rargs[0] = expr_copy(g->data.function.args[0]);
    rargs[1] = expr_copy(rules);
    Expr* one = expr_new_integer(1);
    rargs[2] = expr_new_function(expr_new_symbol(SYM_List), &one, 1);
    Expr* call = expr_new_function(expr_new_symbol("Replace"), rargs, 3);
    Expr* nl = evaluate(call);
    expr_free(call);

    GopsView v;
    if (!graph_is_list(nl) || !gops_view(g, &v)
        || nl->data.function.arg_count != v.nv) { expr_free(nl); return NULL; }
    size_t nv = v.nv, ne = v.ne;
    Expr* const* img = nl->data.function.args;

    int* map = gops_malloc(nv, sizeof(int));
    Expr** verts = gops_malloc(nv, sizeof(Expr*));
    Expr** edges = gops_malloc(ne, sizeof(Expr*));
    Expr** ws = v.weights ? gops_malloc(ne, sizeof(Expr*)) : NULL;
    int* eu = gops_malloc(ne, sizeof(int));
    int* ev = gops_malloc(ne, sizeof(int));
    unsigned char* ed = gops_malloc(ne, 1);
    GraphVIdx* ix = graph_vidx_new(nv);
    GopsHeads hs = { { NULL, NULL } };
    if (!map || !verts || !edges || (v.weights && !ws) || !eu || !ev || !ed || !ix) {
        free(map); free(verts); free(edges); free(ws); free(eu); free(ev); free(ed);
        graph_vidx_free(ix); expr_free(nl);
        return NULL;
    }
    size_t nv2 = 0;
    for (size_t i = 0; i < nv; i++) {
        int p = graph_vidx_get(ix, img[i]);
        if (p < 0) {
            p = (int)nv2;
            verts[nv2++] = expr_copy(img[i]);
            graph_vidx_put(ix, verts[p], p);
        }
        map[i] = p;
    }
    graph_vidx_free(ix);
    for (size_t k = 0; k < ne; k++) {
        int a = v.eu[k], b = v.ev[k];
        /* eu[k] is the edge's first argument, so orientation is preserved. */
        if (img[a] == v.verts[a] && img[b] == v.verts[b])
            edges[k] = expr_copy(v.edges[k]);
        else
            edges[k] = gops_edge(&hs, v.edir[k], expr_copy(verts[map[a]]),
                                 expr_copy(verts[map[b]]));
        if (ws) ws[k] = expr_copy(v.weights[k]);
        eu[k] = map[a]; ev[k] = map[b]; ed[k] = v.edir[k];
    }
    gops_heads_free(&hs);
    free(map);
    expr_free(nl);
    /* A merge that creates a self-loop or parallel edge: unevaluated. */
    return gops_graph_new(verts, nv2, edges, ne, ws, eu, ev, ed);
}

/* ---- EdgeRules ------------------------------------------------------------ */
Expr* builtin_edge_rules(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* el = g->data.function.args[1];
    size_t ne = el->data.function.arg_count;
    Expr** out = gops_malloc(ne, sizeof(Expr*));
    if (!out) return NULL;
    Expr* rh = expr_new_symbol(SYM_Rule);
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = el->data.function.args[k];
        Expr* a[2] = { expr_copy(e->data.function.args[0]), expr_copy(e->data.function.args[1]) };
        out[k] = expr_new_function(expr_copy(rh), a, 2);
    }
    expr_free(rh);
    return gops_list_take(out, ne);
}

/* ---- VertexIndex / EdgeIndex ---------------------------------------------- */
Expr* builtin_vertex_index(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    Expr* arg = res->data.function.args[1];
    int p = graph_vertex_position(g, arg);
    if (p >= 0) return expr_new_integer(p + 1);
    if (p == -2 || !graph_is_list(arg)) return NULL;
    size_t n = arg->data.function.arg_count;
    Expr** out = gops_malloc(n, sizeof(Expr*));
    if (!out) return NULL;
    for (size_t i = 0; i < n; i++) {
        int q = graph_vertex_position(g, arg->data.function.args[i]);
        if (q < 0) {
            for (size_t j = 0; j < i; j++) expr_free(out[j]);
            free(out);
            return NULL;
        }
        out[i] = expr_new_integer(q + 1);
    }
    return gops_list_take(out, n);
}

Expr* builtin_edge_index(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* arg = res->data.function.args[1];
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    const Expr *a, *b;
    int single = gops_parse_edge(arg, &a, &b) >= 0;
    if (!single && !graph_is_list(arg)) return NULL;
    size_t ni;
    Expr* const* items = single ? &res->data.function.args[1]
                                : gops_items(&res->data.function.args[1], &ni);
    if (single) ni = 1;
    int* idx = gops_malloc(ni, sizeof(int));
    if (!idx || !resolve_edges(&v, items, ni, idx, 1)) { free(idx); return NULL; }
    for (size_t i = 0; i < ni; i++) if (idx[i] < 0) { free(idx); return NULL; }
    Expr* out;
    if (single) {
        out = expr_new_integer(idx[0] + 1);
    } else {
        Expr** l = gops_malloc(ni, sizeof(Expr*));
        if (!l) { free(idx); return NULL; }
        for (size_t i = 0; i < ni; i++) l[i] = expr_new_integer(idx[i] + 1);
        out = gops_list_take(l, ni);
    }
    free(idx);
    return out;
}

/* ---- IndexGraph ----------------------------------------------------------- */
Expr* builtin_index_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    int64_t r = 1;
    if (argc == 2) {
        const Expr* ra = res->data.function.args[1];
        if (ra->type != EXPR_INTEGER) return NULL;
        r = ra->data.integer;
    }
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    size_t nv = v.nv, ne = v.ne;
    if (nv > 0 && r > INT64_MAX - (int64_t)nv) return NULL;
    Expr** verts = gops_malloc(nv, sizeof(Expr*));
    Expr** edges = gops_malloc(ne, sizeof(Expr*));
    Expr** ws = v.weights ? copy_exprs(v.weights, ne, 0) : NULL;
    int *eu, *ev; unsigned char* ed;
    if (!verts || !edges || (v.weights && !ws) || !copy_endpoints(&v, &eu, &ev, &ed)) {
        free(verts); free(edges);
        if (ws) { for (size_t i = 0; i < ne; i++) expr_free(ws[i]); free(ws); }
        return NULL;
    }
    GopsHeads hs = { { NULL, NULL } };
    for (size_t i = 0; i < nv; i++) verts[i] = expr_new_integer(r + (int64_t)i);
    for (size_t k = 0; k < ne; k++) {
        /* the edge's own orientation: eu[k] is its first argument */
        edges[k] = gops_edge(&hs, v.edir[k], expr_copy(verts[v.eu[k]]), expr_copy(verts[v.ev[k]]));
    }
    gops_heads_free(&hs);
    return gops_graph_new(verts, nv, edges, ne, ws, eu, ev, ed);
}
