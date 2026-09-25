/* gops_transform.c - whole-graph transforms.
 *
 *   GraphComplement[g]   the complement on the same vertices. Undirected g:
 *                        every non-adjacent pair i < j as i <-> j. Directed or
 *                        mixed g: every ordered pair (i, j), i != j, with no
 *                        edge usable from i to j (a directed i -> j or an
 *                        undirected i <-> j) as i -> j. Row-major order by
 *                        VertexList position; weights are dropped.
 *   ReverseGraph[g]      every directed edge reversed, undirected edges kept;
 *                        order and weights kept.
 *   UndirectedGraph[g]   direction dropped; u -> v and v -> u merge into one
 *                        edge whose weight is the sum of theirs. Edges are
 *                        oriented and ordered by VertexList position (upper
 *                        triangle, row-major), as Mathematica does. An
 *                        undirected g is returned unchanged.
 *   DirectedGraph[g]     each undirected u <-> v becomes u -> v, v -> u (in
 *                        place, weight duplicated); directed edges kept.
 *   DirectedGraph[g, "Acyclic"]
 *                        each undirected edge oriented from the earlier to the
 *                        later vertex in VertexList, which gives a DAG when g
 *                        is undirected; then edges are sorted by (tail, head)
 *                        position. For a mixed g the edges stay in place (and
 *                        directed edges are kept as they are), as Mathematica.
 *   LineGraph[g]         vertices 1..m, one per edge (EdgeList position).
 *                        Undirected: j <-> i for edges sharing an endpoint,
 *                        listed by j, then by the shared endpoint (j's first,
 *                        then second argument), then by i < j. Directed:
 *                        i -> j when the head of edge i is the tail of edge j,
 *                        in (i, j) order. Mixed graphs are left unevaluated.
 *                        Deviation: Mathematica numbers directed line-graph
 *                        vertices in a traversal order of its own; Mathilda
 *                        always uses EdgeList position (an isomorphic graph).
 *
 * All O(V + E) (GraphComplement and LineGraph are output-bound), integer
 * passes over the memo's endpoints, with nodes shared where unchanged and the
 * result seeded into the validated-graph memo (gops_graph_new).
 *
 * Memory (SPEC section 4): results are fresh; res is borrowed.
 */

#include "graph_ops.h"
#include "eval.h"
#include "core.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

static Expr** share_verts(const GopsView* v) {
    Expr** a = gops_malloc(v->nv, sizeof(Expr*));
    if (!a) return NULL;
    for (size_t i = 0; i < v->nv; i++) a[i] = expr_copy(v->verts[i]);
    return a;
}

static void free_exprs(Expr** a, size_t n) {
    if (!a) return;
    for (size_t i = 0; i < n; i++) if (a[i]) expr_free(a[i]);
    free(a);
}

/* Growable edge buffer for output-bound builders. */
typedef struct {
    Expr** e; int* u; int* v; unsigned char* d;
    size_t n, cap;
} EdgeBuf;

static int eb_reserve(EdgeBuf* b, size_t want) {
    if (want <= b->cap) return 1;
    size_t cap = b->cap ? b->cap : 64;
    while (cap < want) cap *= 2;
    Expr** e = realloc(b->e, cap * sizeof(Expr*));
    if (!e) return 0;
    b->e = e;
    int* u = realloc(b->u, cap * sizeof(int));
    if (!u) return 0;
    b->u = u;
    int* v = realloc(b->v, cap * sizeof(int));
    if (!v) return 0;
    b->v = v;
    unsigned char* d = realloc(b->d, cap);
    if (!d) return 0;
    b->d = d;
    b->cap = cap;
    return 1;
}

static void eb_free(EdgeBuf* b) {
    free_exprs(b->e, b->n);
    free(b->u); free(b->v); free(b->d);
    memset(b, 0, sizeof(*b));
}

/* Hand the buffer to gops_graph_new (which owns everything afterwards). */
static Expr* eb_graph(EdgeBuf* b, Expr** verts, size_t nv) {
    if (!b->e) {                       /* no edges: allocate the 1-slot arrays */
        b->e = gops_malloc(0, sizeof(Expr*));
        b->u = gops_malloc(0, sizeof(int));
        b->v = gops_malloc(0, sizeof(int));
        b->d = gops_malloc(0, 1);
        if (!b->e || !b->u || !b->v || !b->d) { eb_free(b); free_exprs(verts, nv); return NULL; }
    }
    Expr* g = gops_graph_new(verts, nv, b->e, b->n, NULL, b->u, b->v, b->d);
    memset(b, 0, sizeof(*b));
    return g;
}

/* ---- GraphComplement ------------------------------------------------------ */
Expr* builtin_graph_complement(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    size_t n = v.nv;
    int undirected = (v.ndir == 0);
    GopsInc inc;
    if (!gops_inc_build(&v, GOPS_INC_OUT, &inc)) return NULL;
    Expr** verts = share_verts(&v);
    int* stamp = gops_malloc(n, sizeof(int));
    EdgeBuf b;
    memset(&b, 0, sizeof(b));
    size_t total = undirected ? (n * (n - (n > 0)) / 2 - v.ne)
                              : (n * (n - (n > 0)) - (inc.start[n]));
    if (!verts || !stamp || !eb_reserve(&b, total + 1)) {
        gops_inc_free(&inc); free_exprs(verts, n); free(stamp); eb_free(&b);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) stamp[i] = -1;
    GopsHeads hs = { { NULL, NULL } };
    for (size_t i = 0; i < n; i++) {
        for (int j = inc.start[i]; j < inc.start[i + 1]; j++) stamp[inc.nbr[j]] = (int)i;
        stamp[i] = (int)i;
        for (size_t j = undirected ? i + 1 : 0; j < n; j++) {
            if (stamp[j] == (int)i) continue;
            b.e[b.n] = gops_edge(&hs, !undirected, expr_copy(verts[i]), expr_copy(verts[j]));
            b.u[b.n] = (int)i; b.v[b.n] = (int)j; b.d[b.n] = (unsigned char)!undirected;
            b.n++;
        }
        if ((i & 63) == 63) tc_check_deadline();
    }
    gops_heads_free(&hs);
    gops_inc_free(&inc);
    free(stamp);
    return eb_graph(&b, verts, n);
}

/* ---- ReverseGraph --------------------------------------------------------- */
Expr* builtin_reverse_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    if (v.ndir == 0) return expr_copy(res->data.function.args[0]);
    size_t ne = v.ne;
    Expr** verts = share_verts(&v);
    Expr** edges = gops_malloc(ne, sizeof(Expr*));
    Expr** ws = v.weights ? gops_malloc(ne, sizeof(Expr*)) : NULL;
    int* eu = gops_malloc(ne, sizeof(int));
    int* ev = gops_malloc(ne, sizeof(int));
    unsigned char* ed = gops_malloc(ne, 1);
    if (!verts || !edges || (v.weights && !ws) || !eu || !ev || !ed) {
        free_exprs(verts, v.nv); free(edges); free(ws); free(eu); free(ev); free(ed);
        return NULL;
    }
    GopsHeads hs = { { NULL, NULL } };
    for (size_t k = 0; k < ne; k++) {
        if (v.edir[k]) {
            const Expr* e = v.edges[k];
            edges[k] = gops_edge(&hs, 1, expr_copy(e->data.function.args[1]),
                                 expr_copy(e->data.function.args[0]));
            eu[k] = v.ev[k]; ev[k] = v.eu[k];
        } else {
            edges[k] = expr_copy(v.edges[k]);
            eu[k] = v.eu[k]; ev[k] = v.ev[k];
        }
        ed[k] = v.edir[k];
        if (ws) ws[k] = expr_copy(v.weights[k]);
    }
    gops_heads_free(&hs);
    return gops_graph_new(verts, v.nv, edges, ne, ws, eu, ev, ed);
}

/* ---- UndirectedGraph ------------------------------------------------------ */
Expr* builtin_undirected_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    if (v.ndir == 0) return expr_copy(res->data.function.args[0]);
    /* Summing merged weights evaluates Plus, which could run arbitrary code
     * and evict g from the memo: work on private endpoint copies. */
    if (v.weights && !gops_view_own(&v)) return NULL;
    size_t ne = v.ne, n = v.nv;
    int* perm = gops_malloc(ne, sizeof(int));
    int* lo = gops_malloc(ne, sizeof(int));
    int* hi = gops_malloc(ne, sizeof(int));
    if (!perm || !lo || !hi) { free(perm); free(lo); free(hi); gops_view_free(&v); return NULL; }
    for (size_t k = 0; k < ne; k++) {
        int a = v.eu[k], b = v.ev[k];
        perm[k] = (int)k;
        lo[k] = a < b ? a : b;
        hi[k] = a < b ? b : a;
    }
    /* Row-major (lo, hi); a stable sort keeps EdgeList order among merges. */
    if (!gops_csort(perm, ne, hi, (int)n) || !gops_csort(perm, ne, lo, (int)n)) {
        free(perm); free(lo); free(hi); gops_view_free(&v); return NULL;
    }
    Expr** verts = share_verts(&v);
    Expr** edges = gops_malloc(ne, sizeof(Expr*));
    Expr** ws = v.weights ? gops_malloc(ne, sizeof(Expr*)) : NULL;
    int* eu = gops_malloc(ne, sizeof(int));
    int* ev = gops_malloc(ne, sizeof(int));
    unsigned char* ed = gops_calloc(ne, 1);
    if (!verts || !edges || (v.weights && !ws) || !eu || !ev || !ed) {
        free_exprs(verts, n); free(edges); free(ws); free(eu); free(ev); free(ed);
        free(perm); free(lo); free(hi); gops_view_free(&v);
        return NULL;
    }
    GopsHeads hs = { { NULL, NULL } };
    size_t m = 0;
    for (size_t i = 0; i < ne; ) {
        int k = perm[i];
        size_t j = i + 1;
        while (j < ne && lo[perm[j]] == lo[k] && hi[perm[j]] == hi[k]) j++;
        /* Share the node when it is already the undirected lo <-> hi edge. */
        if (!v.edir[k] && v.eu[k] == lo[k])
            edges[m] = expr_copy(v.edges[k]);
        else
            edges[m] = gops_edge(&hs, 0, expr_copy(verts[lo[k]]), expr_copy(verts[hi[k]]));
        if (ws) {
            if (j - i == 1) {
                ws[m] = expr_copy(v.weights[k]);
            } else {                   /* merged pair: the weights' sum */
                size_t c = j - i;
                Expr** terms = gops_malloc(c, sizeof(Expr*));
                if (!terms) { ws[m] = expr_copy(v.weights[k]); }
                else {
                    for (size_t t = 0; t < c; t++) terms[t] = expr_copy(v.weights[perm[i + t]]);
                    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, c);
                    free(terms);
                    ws[m] = evaluate(sum);
                    expr_free(sum);
                }
            }
        }
        eu[m] = lo[k]; ev[m] = hi[k];
        m++;
        i = j;
    }
    gops_heads_free(&hs);
    free(perm); free(lo); free(hi);
    gops_view_free(&v);
    return gops_graph_new(verts, n, edges, m, ws, eu, ev, ed);
}

/* ---- DirectedGraph -------------------------------------------------------- */
static int is_string(const Expr* e, const char* s) {
    return e && e->type == EXPR_STRING && strcmp(e->data.string, s) == 0;
}

Expr* builtin_directed_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    int acyclic = 0;
    if (argc == 2) {
        if (!is_string(res->data.function.args[1], "Acyclic")) return NULL;
        acyclic = 1;
    }
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    size_t ne = v.ne, n = v.nv;
    if (v.ndir == ne) return expr_copy(res->data.function.args[0]);
    size_t m = acyclic ? ne : ne + (ne - v.ndir);

    Expr** verts = share_verts(&v);
    Expr** edges = gops_malloc(m, sizeof(Expr*));
    Expr** ws = v.weights ? gops_malloc(m, sizeof(Expr*)) : NULL;
    int* eu = gops_malloc(m, sizeof(int));
    int* ev = gops_malloc(m, sizeof(int));
    unsigned char* ed = gops_malloc(m, 1);
    int* perm = gops_malloc(ne, sizeof(int));
    if (!verts || !edges || (v.weights && !ws) || !eu || !ev || !ed || !perm) {
        free_exprs(verts, n); free(edges); free(ws); free(eu); free(ev); free(ed); free(perm);
        return NULL;
    }
    GopsHeads hs = { { NULL, NULL } };
    size_t j = 0;
    if (!acyclic) {
        for (size_t k = 0; k < ne; k++) {
            int a = v.eu[k], b = v.ev[k];
            if (v.edir[k]) {
                edges[j] = expr_copy(v.edges[k]);
                eu[j] = a; ev[j] = b; ed[j] = 1;
                if (ws) ws[j] = expr_copy(v.weights[k]);
                j++;
                continue;
            }
            edges[j] = gops_edge(&hs, 1, expr_copy(verts[a]), expr_copy(verts[b]));
            eu[j] = a; ev[j] = b; ed[j] = 1;
            if (ws) ws[j] = expr_copy(v.weights[k]);
            j++;
            edges[j] = gops_edge(&hs, 1, expr_copy(verts[b]), expr_copy(verts[a]));
            eu[j] = b; ev[j] = a; ed[j] = 1;
            if (ws) ws[j] = expr_copy(v.weights[k]);
            j++;
        }
    } else {
        /* Orient each undirected edge from the earlier vertex; an undirected
         * g is then sorted by (tail, head), a mixed one keeps its order. */
        int* tl = gops_malloc(ne, sizeof(int));
        int* hd = gops_malloc(ne, sizeof(int));
        if (!tl || !hd) {
            free(tl); free(hd); free_exprs(verts, n); free(edges); free(ws);
            free(eu); free(ev); free(ed); free(perm);
            return NULL;
        }
        for (size_t k = 0; k < ne; k++) {
            int a = v.eu[k], b = v.ev[k];
            if (!v.edir[k] && a > b) { int t = a; a = b; b = t; }
            tl[k] = a; hd[k] = b; perm[k] = (int)k;
        }
        if (v.ndir == 0 && (!gops_csort(perm, ne, hd, (int)n) || !gops_csort(perm, ne, tl, (int)n))) {
            free(tl); free(hd); free_exprs(verts, n); free(edges); free(ws);
            free(eu); free(ev); free(ed); free(perm);
            return NULL;
        }
        for (size_t i = 0; i < ne; i++) {
            int k = perm[i];
            if (v.edir[k]) edges[j] = expr_copy(v.edges[k]);
            else edges[j] = gops_edge(&hs, 1, expr_copy(verts[tl[k]]), expr_copy(verts[hd[k]]));
            eu[j] = tl[k]; ev[j] = hd[k]; ed[j] = 1;
            if (ws) ws[j] = expr_copy(v.weights[k]);
            j++;
        }
        free(tl); free(hd);
    }
    free(perm);
    gops_heads_free(&hs);
    return gops_graph_new(verts, n, edges, j, ws, eu, ev, ed);
}

/* ---- LineGraph ------------------------------------------------------------ */
Expr* builtin_line_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return NULL;
    size_t ne = v.ne;
    if (v.ndir != 0 && v.ndir != ne) return NULL;          /* mixed */
    int directed = (v.ndir == ne && ne > 0);
    GopsInc inc;
    if (!gops_inc_build(&v, directed ? GOPS_INC_OUT : GOPS_INC_ALL, &inc)) return NULL;
    Expr** verts = gops_malloc(ne, sizeof(Expr*));
    EdgeBuf b;
    memset(&b, 0, sizeof(b));
    if (!verts) { gops_inc_free(&inc); return NULL; }
    for (size_t k = 0; k < ne; k++) verts[k] = expr_new_integer((int64_t)k + 1);
    GopsHeads hs = { { NULL, NULL } };
    int ok = 1;
    for (size_t j = 0; j < ne && ok; j++) {
        if (directed) {
            /* j -> i for every edge i leaving the head of edge j */
            int h = v.ev[j];
            for (int t = inc.start[h]; t < inc.start[h + 1]; t++) {
                int i = inc.eid[t];
                if (!eb_reserve(&b, b.n + 1)) { ok = 0; break; }
                b.e[b.n] = gops_edge(&hs, 1, expr_copy(verts[j]), expr_copy(verts[i]));
                b.u[b.n] = (int)j; b.v[b.n] = i; b.d[b.n] = 1; b.n++;
            }
        } else {
            int ends[2] = { v.eu[j], v.ev[j] };
            for (int s = 0; s < 2 && ok; s++) {
                int x = ends[s];
                for (int t = inc.start[x]; t < inc.start[x + 1]; t++) {
                    int i = inc.eid[t];
                    if (i >= (int)j) break;          /* incidence is in edge order */
                    if (!eb_reserve(&b, b.n + 1)) { ok = 0; break; }
                    b.e[b.n] = gops_edge(&hs, 0, expr_copy(verts[j]), expr_copy(verts[i]));
                    b.u[b.n] = (int)j; b.v[b.n] = i; b.d[b.n] = 0; b.n++;
                }
            }
        }
        if ((j & 255) == 255) tc_check_deadline();
    }
    gops_heads_free(&hs);
    gops_inc_free(&inc);
    if (!ok) { eb_free(&b); free_exprs(verts, ne); return NULL; }
    return eb_graph(&b, verts, ne);
}
