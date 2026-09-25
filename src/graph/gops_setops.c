/* gops_setops.c - graph set operations.
 *
 *   GraphUnion[g1, g2, ...]          vertices: the union, in canonical order
 *                                    (as Union). Edges: the distinct edges of
 *                                    all graphs (an undirected edge equals its
 *                                    reversal). All undirected: first-appearance
 *                                    order, each oriented by canonical vertex
 *                                    order. All directed: first-appearance
 *                                    order. Mixed: canonical (Sort) order.
 *                                    GraphUnion[g] is g.
 *   GraphIntersection[g1, g2, ...]   vertices: the union (canonical order);
 *                                    edges: those of g1 present in every graph,
 *                                    in canonical order (g1's orientation).
 *   GraphDifference[g1, g2]          vertices: the union (canonical order);
 *                                    edges: those of g1 absent from g2, in
 *                                    canonical order.
 *   GraphDisjointUnion[g1, g2, ...]  vertices relabelled 1..n (g1's first, in
 *                                    VertexList order, then g2's, ...), edges
 *                                    translated in order.
 * These orders are Mathematica 15's. Weights are dropped, as in Mathematica.
 *
 * Algorithm. Each graph's vertices map to union ids -- by an O(V) elementwise
 * SameQ check when a vertex list equals the first graph's (the common case of
 * two graphs on one vertex set: no hashing at all), else through one hash
 * index over the union. The union is ordered with expr_compare, skipped when it
 * is already sorted (O(V) check) and done on machine integers when every vertex
 * is one. Edges then become integer keys over result positions; duplicates are
 * found in a GopsKeySet and canonical edge order is a stable counting sort on
 * (kind, first, second) positions -- O(V + E) apart from the vertex sort.
 *
 * Memory (SPEC section 4): results are fresh; res is borrowed.
 */

#include "graph_ops.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

/* ---- Vertex union --------------------------------------------------------- */
typedef struct {
    size_t   ng;
    int**    map;       /* map[t][i]: result position of graph t's vertex i  */
    Expr**   verts;     /* result vertices, canonical order (borrowed nodes) */
    size_t   nv;
} Union;

static void union_free(Union* u) {
    if (u->map) for (size_t t = 0; t < u->ng; t++) free(u->map[t]);
    free(u->map);
    free(u->verts);
    memset(u, 0, sizeof(*u));
}

static int cmp_expr_ptr(const void* a, const void* b) {
    return expr_compare(*(Expr* const*)a, *(Expr* const*)b);
}

typedef struct { int64_t v; int id; } IntKey;
static int cmp_intkey(const void* a, const void* b) {
    const IntKey* x = a; const IntKey* y = b;
    return (x->v > y->v) - (x->v < y->v);
}

/* Build the canonical vertex union of the graphs gs[0..ng-1] (all valid). */
static int union_build(Expr* const* gs, size_t ng, Union* u) {
    memset(u, 0, sizeof(*u));
    u->ng = ng;
    u->map = gops_calloc(ng, sizeof(int*));
    if (!u->map) return 0;
    size_t cap = 0;
    for (size_t t = 0; t < ng; t++) cap += gs[t]->data.function.args[0]->data.function.arg_count;
    Expr** uv = gops_malloc(cap, sizeof(Expr*));      /* union, first-seen order */
    if (!uv) { union_free(u); return 0; }
    size_t nu = 0;
    GraphVIdx* ux = NULL;                               /* built on demand */
    const Expr* v0 = gs[0]->data.function.args[0];
    size_t n0 = v0->data.function.arg_count;

    for (size_t t = 0; t < ng; t++) {
        const Expr* vl = gs[t]->data.function.args[0];
        size_t n = vl->data.function.arg_count;
        int* m = gops_malloc(n, sizeof(int));
        if (!m) goto fail;
        u->map[t] = m;
        if (t == 0) {
            for (size_t i = 0; i < n; i++) { uv[i] = vl->data.function.args[i]; m[i] = (int)i; }
            nu = n;
            continue;
        }
        /* Fast path: the same vertex list as graph 0 (pointer or SameQ). */
        size_t i = 0;
        if (n == n0 && !ux) {
            while (i < n && (vl->data.function.args[i] == v0->data.function.args[i]
                             || expr_eq(vl->data.function.args[i], v0->data.function.args[i]))) {
                m[i] = (int)i; i++;
            }
            if (i == n) continue;
        }
        if (!ux) {
            ux = graph_vidx_new(nu + n);
            if (!ux) goto fail;
            for (size_t j = 0; j < nu; j++) graph_vidx_put(ux, uv[j], (int)j);
        }
        for (; i < n; i++) {
            Expr* x = vl->data.function.args[i];
            int id = graph_vidx_get(ux, x);
            if (id < 0) { id = (int)nu; uv[nu++] = x; graph_vidx_put(ux, x, id); }
            m[i] = id;
        }
    }
    graph_vidx_free(ux);
    ux = NULL;

    /* Canonical order: rank[id] = result position. */
    int* rank = gops_malloc(nu, sizeof(int));
    if (!rank) goto fail;
    int sorted = 1;
    for (size_t i = 0; i + 1 < nu && sorted; i++)
        if (expr_compare(uv[i], uv[i + 1]) > 0) sorted = 0;
    if (sorted) {
        for (size_t i = 0; i < nu; i++) rank[i] = (int)i;
        u->verts = uv;
    } else {
        Expr** sv = gops_malloc(nu, sizeof(Expr*));
        if (!sv) { free(rank); goto fail; }
        int allint = 1;
        for (size_t i = 0; i < nu && allint; i++) allint = (uv[i]->type == EXPR_INTEGER);
        if (allint) {
            IntKey* ks = gops_malloc(nu, sizeof(IntKey));
            if (!ks) { free(sv); free(rank); goto fail; }
            for (size_t i = 0; i < nu; i++) { ks[i].v = uv[i]->data.integer; ks[i].id = (int)i; }
            qsort(ks, nu, sizeof(IntKey), cmp_intkey);
            for (size_t i = 0; i < nu; i++) { sv[i] = uv[ks[i].id]; rank[ks[i].id] = (int)i; }
            free(ks);
        } else {
            /* Sort pointers, then recover ids through a fresh index. */
            memcpy(sv, uv, nu * sizeof(Expr*));
            qsort(sv, nu, sizeof(Expr*), cmp_expr_ptr);
            GraphVIdx* sx = graph_vidx_new(nu);
            if (!sx) { free(sv); free(rank); goto fail; }
            for (size_t i = 0; i < nu; i++) graph_vidx_put(sx, sv[i], (int)i);
            for (size_t i = 0; i < nu; i++) rank[i] = graph_vidx_get(sx, uv[i]);
            graph_vidx_free(sx);
        }
        free(uv);
        uv = NULL;
        u->verts = sv;
        for (size_t t = 0; t < ng; t++) {
            size_t n = gs[t]->data.function.args[0]->data.function.arg_count;
            for (size_t i = 0; i < n; i++) u->map[t][i] = rank[u->map[t][i]];
        }
    }
    free(rank);
    u->nv = nu;
    return 1;

fail:
    graph_vidx_free(ux);
    free(uv);
    u->verts = NULL;
    union_free(u);
    return 0;
}

/* ---- Result assembly ------------------------------------------------------ */

/* A kept edge: result positions a -> b (a is the first argument), its kind,
 * and the source edge node it may share. */
typedef struct { int a, b; unsigned char d; const Expr* src; int share; } OutEdge;

static Expr* union_graph(const Union* u, OutEdge* oe, const int* order, size_t m) {
    Expr** verts = gops_malloc(u->nv, sizeof(Expr*));
    Expr** edges = gops_malloc(m, sizeof(Expr*));
    int* eu = gops_malloc(m, sizeof(int));
    int* ev = gops_malloc(m, sizeof(int));
    unsigned char* ed = gops_malloc(m, 1);
    if (!verts || !edges || !eu || !ev || !ed) {
        free(verts); free(edges); free(eu); free(ev); free(ed);
        return NULL;
    }
    for (size_t i = 0; i < u->nv; i++) verts[i] = expr_copy(u->verts[i]);
    GopsHeads hs = { { NULL, NULL } };
    for (size_t i = 0; i < m; i++) {
        const OutEdge* e = &oe[order ? order[i] : (int)i];
        edges[i] = e->share ? expr_copy((Expr*)e->src)
                            : gops_edge(&hs, e->d, expr_copy(verts[e->a]), expr_copy(verts[e->b]));
        eu[i] = e->a; ev[i] = e->b; ed[i] = e->d;
    }
    gops_heads_free(&hs);
    return gops_graph_new(verts, u->nv, edges, m, NULL, eu, ev, ed);
}

/* Canonical (Sort) order of kept edges: DirectedEdge before UndirectedEdge,
 * then by first, then second argument position. Returns a permutation. */
static int* canonical_order(const OutEdge* oe, size_t m, size_t nv) {
    int* perm = gops_malloc(m, sizeof(int));
    int* ka = gops_malloc(m, sizeof(int));
    int* kb = gops_malloc(m, sizeof(int));
    int* kd = gops_malloc(m, sizeof(int));
    if (!perm || !ka || !kb || !kd) { free(perm); free(ka); free(kb); free(kd); return NULL; }
    for (size_t i = 0; i < m; i++) {
        perm[i] = (int)i; ka[i] = oe[i].a; kb[i] = oe[i].b; kd[i] = oe[i].d ? 0 : 1;
    }
    int ok = gops_csort(perm, m, kb, (int)nv) && gops_csort(perm, m, ka, (int)nv)
          && gops_csort(perm, m, kd, 2);
    free(ka); free(kb); free(kd);
    if (!ok) { free(perm); return NULL; }
    return perm;
}

/* Validate every argument as a graph; the argument array or NULL. */
static Expr* const* graph_args(Expr* res, size_t min_args) {
    size_t n = res->data.function.arg_count;
    if (n < min_args) return NULL;
    for (size_t t = 0; t < n; t++)
        if (!graph_is_valid(res->data.function.args[t])) return NULL;
    return res->data.function.args;
}

/* ---- GraphUnion ----------------------------------------------------------- */
Expr* builtin_graph_union(Expr* res) {
    Expr* const* gs = graph_args(res, 1);
    if (!gs) return NULL;
    size_t ng = res->data.function.arg_count;
    if (ng == 1) return expr_copy(gs[0]);
    Union u;
    if (!union_build(gs, ng, &u)) return NULL;
    size_t total = 0;
    for (size_t t = 0; t < ng; t++) total += gs[t]->data.function.args[1]->data.function.arg_count;
    OutEdge* oe = gops_malloc(total, sizeof(OutEdge));
    GopsKeySet ks;
    if (!oe || !gops_keyset_init(&ks, total)) { free(oe); union_free(&u); return NULL; }
    size_t m = 0, ndir = 0;
    for (size_t t = 0; t < ng; t++) {
        GopsView v;
        if (!gops_view(gs[t], &v)) { m = (size_t)-1; break; }
        const int* mp = u.map[t];
        for (size_t k = 0; k < v.ne; k++) {
            int a = mp[v.eu[k]], b = mp[v.ev[k]], d = v.edir[k];
            if (gops_keyset_put(&ks, gops_edge_key(a, b, d), (int)m) != (int)m) continue;
            oe[m].a = a; oe[m].b = b; oe[m].d = (unsigned char)d;
            oe[m].src = v.edges[k]; oe[m].share = 1;
            ndir += (size_t)d;
            m++;
        }
    }
    gops_keyset_free(&ks);
    Expr* out = NULL;
    if (m != (size_t)-1) {
        int* order = NULL;
        if (ndir == 0) {
            /* all undirected: orient by canonical vertex order */
            for (size_t i = 0; i < m; i++)
                if (oe[i].a > oe[i].b) { int x = oe[i].a; oe[i].a = oe[i].b; oe[i].b = x; oe[i].share = 0; }
        } else if (ndir < m) {
            order = canonical_order(oe, m, u.nv);
            if (!order) { free(oe); union_free(&u); return NULL; }
        }
        out = union_graph(&u, oe, order, m);
        free(order);
    }
    free(oe);
    union_free(&u);
    return out;
}

/* ---- GraphIntersection / GraphDifference ---------------------------------- */

/* Edges of gs[0] kept when (hits == ng - 1) for intersection, (hits == 0) for
 * difference, where hits counts the other graphs containing the edge. */
static Expr* filter_first(Expr* const* gs, size_t ng, int intersect) {
    Union u;
    if (!union_build(gs, ng, &u)) return NULL;
    GopsView v0;
    if (!gops_view(gs[0], &v0)) { union_free(&u); return NULL; }
    size_t ne0 = v0.ne;
    int* hits = gops_calloc(ne0, sizeof(int));
    OutEdge* oe = gops_malloc(ne0, sizeof(OutEdge));
    GopsKeySet ks;
    if (!hits || !oe || !gops_keyset_init(&ks, ne0)) {
        free(hits); free(oe); union_free(&u); return NULL;
    }
    for (size_t k = 0; k < ne0; k++) {
        oe[k].a = u.map[0][v0.eu[k]]; oe[k].b = u.map[0][v0.ev[k]];
        oe[k].d = v0.edir[k]; oe[k].src = v0.edges[k]; oe[k].share = 1;
        gops_keyset_put(&ks, gops_edge_key(oe[k].a, oe[k].b, oe[k].d), (int)k);
    }
    int ok = 1;
    for (size_t t = 1; t < ng && ok; t++) {
        GopsView v;
        if (!gops_view(gs[t], &v)) { ok = 0; break; }
        const int* mp = u.map[t];
        for (size_t k = 0; k < v.ne; k++) {
            int e = gops_keyset_get(&ks, gops_edge_key(mp[v.eu[k]], mp[v.ev[k]], v.edir[k]));
            if (e >= 0) hits[e]++;
        }
    }
    gops_keyset_free(&ks);
    Expr* out = NULL;
    if (ok) {
        size_t m = 0;
        int want = intersect ? (int)(ng - 1) : 0;
        for (size_t k = 0; k < ne0; k++) if (hits[k] == want) oe[m++] = oe[k];
        int* order = canonical_order(oe, m, u.nv);
        if (order) out = union_graph(&u, oe, order, m);
        free(order);
    }
    free(hits); free(oe);
    union_free(&u);
    return out;
}

Expr* builtin_graph_intersection(Expr* res) {
    Expr* const* gs = graph_args(res, 1);
    if (!gs) return NULL;
    if (res->data.function.arg_count == 1) return expr_copy(gs[0]);
    return filter_first(gs, res->data.function.arg_count, 1);
}

Expr* builtin_graph_difference(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* const* gs = graph_args(res, 2);
    return gs ? filter_first(gs, 2, 0) : NULL;
}

/* ---- GraphDisjointUnion --------------------------------------------------- */
Expr* builtin_graph_disjoint_union(Expr* res) {
    Expr* const* gs = graph_args(res, 1);
    if (!gs) return NULL;
    size_t ng = res->data.function.arg_count;
    if (ng == 1) return expr_copy(gs[0]);        /* nothing to separate */
    size_t nv = 0, ne = 0;
    for (size_t t = 0; t < ng; t++) {
        nv += gs[t]->data.function.args[0]->data.function.arg_count;
        ne += gs[t]->data.function.args[1]->data.function.arg_count;
    }
    Expr** verts = gops_malloc(nv, sizeof(Expr*));
    Expr** edges = gops_malloc(ne, sizeof(Expr*));
    int* eu = gops_malloc(ne, sizeof(int));
    int* ev = gops_malloc(ne, sizeof(int));
    unsigned char* ed = gops_malloc(ne, 1);
    if (!verts || !edges || !eu || !ev || !ed) {
        free(verts); free(edges); free(eu); free(ev); free(ed);
        return NULL;
    }
    for (size_t i = 0; i < nv; i++) verts[i] = expr_new_integer((int64_t)i + 1);
    GopsHeads hs = { { NULL, NULL } };
    size_t off = 0, j = 0;
    for (size_t t = 0; t < ng; t++) {
        GopsView v;
        if (!gops_view(gs[t], &v)) {           /* cannot happen: validated above */
            for (size_t i = 0; i < nv; i++) expr_free(verts[i]);
            for (size_t i = 0; i < j; i++) expr_free(edges[i]);
            free(verts); free(edges); free(eu); free(ev); free(ed);
            gops_heads_free(&hs);
            return NULL;
        }
        for (size_t k = 0; k < v.ne; k++, j++) {
            int a = (int)off + v.eu[k], b = (int)off + v.ev[k];
            edges[j] = gops_edge(&hs, v.edir[k], expr_copy(verts[a]), expr_copy(verts[b]));
            eu[j] = a; ev[j] = b; ed[j] = v.edir[k];
        }
        off += v.nv;
    }
    gops_heads_free(&hs);
    return gops_graph_new(verts, nv, edges, ne, NULL, eu, ev, ed);
}
