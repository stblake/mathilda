/* galg_matching.c - maximum matchings and edge covers.
 *
 *   FindIndependentEdgeSet[g]   a maximum independent edge set (matching)
 *   FindEdgeCover[g]            a minimum edge cover ({} if g has an isolated
 *                               vertex, as in Mathematica)
 *
 * Edge direction is ignored throughout (Wolfram's semantics); the returned
 * edges are g's own edge expressions, in EdgeList order.
 *
 * Algorithm (galg_max_matching):
 *   1. Karp-Sipser greedy start: repeatedly match a degree-1 vertex to its
 *      only neighbour (always safe -- some maximum matching contains that
 *      edge), and when none is left match an arbitrary remaining edge. On
 *      sparse random graphs this alone is optimal or within a handful of
 *      edges of it.
 *   2. Exact completion:
 *        bipartite graphs -- Hopcroft-Karp (shortest augmenting paths in
 *          phases; O(E sqrt V));
 *        general graphs -- Edmonds' blossom algorithm, one alternating-tree
 *          search per exposed vertex, with blossom bases kept in an array and
 *          only the vertices the search touched reset afterwards (so a search
 *          costs O(size of its tree), not O(V)). A search that fails leaves a
 *          Hungarian tree, whose vertices can never lie on an augmenting path
 *          again (Edmonds' lemma), so they are retired for the rest of the run.
 *   Both are exact maximum-cardinality algorithms.
 *
 * Minimum edge cover (Gallai): a maximum matching plus one arbitrary incident
 * edge for every exposed vertex; its size is n - nu(g).
 *
 * Memory (SPEC section 4): fresh results; scratch freed on every path.
 */

#include "graph.h"
#include "graph_algos.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

/* ---- Karp-Sipser greedy ---------------------------------------------------- */

static long gm_karp_sipser(const GalgUG* u, int* mate) {
    int n = u->n;
    int* deg = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* q = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!deg || !q) { free(deg); free(q); return -1; }
    long size = 0;
    int qh = 0, qt = 0;
    for (int v = 0; v < n; v++) {
        mate[v] = -1;
        deg[v] = u->off[v + 1] - u->off[v];
        if (deg[v] == 1) q[qt++] = v;
    }
    int scan = 0;
    for (;;) {
        int v = -1;
        while (qh < qt) {
            int x = q[qh++];
            if (mate[x] < 0 && deg[x] == 1) { v = x; break; }
        }
        int w = -1;
        if (v < 0) {
            /* no degree-1 vertex: an arbitrary live vertex of degree >= 2 */
            while (scan < n && (mate[scan] >= 0 || deg[scan] == 0)) scan++;
            if (scan >= n) break;
            v = scan;
        }
        for (int j = u->off[v]; j < u->off[v + 1]; j++)
            if (mate[u->adj[j]] < 0) { w = u->adj[j]; break; }
        if (w < 0) { deg[v] = 0; continue; }
        mate[v] = w; mate[w] = v; size++;
        /* both leave: update neighbour degrees */
        for (int s = 0; s < 2; s++) {
            int x = s ? w : v;
            for (int j = u->off[x]; j < u->off[x + 1]; j++) {
                int y = u->adj[j];
                if (mate[y] >= 0) continue;
                /* degrees only fall, so a vertex reaches 1 at most once:
                 * at most n pushes in all, and q (length n) never overflows */
                if (--deg[y] == 1) q[qt++] = y;
            }
        }
    }
    free(deg); free(q);
    return size;
}

/* ---- Bipartiteness --------------------------------------------------------- */

/* side[v] in {0,1}; returns 1 if bipartite. */
static int gm_bipartite(const GalgUG* u, char* side, int* queue) {
    int n = u->n;
    memset(side, -1, (size_t)n);
    for (int s = 0; s < n; s++) {
        if (side[s] >= 0) continue;
        int h = 0, t = 0;
        side[s] = 0; queue[t++] = s;
        while (h < t) {
            int v = queue[h++];
            for (int j = u->off[v]; j < u->off[v + 1]; j++) {
                int w = u->adj[j];
                if (side[w] < 0) { side[w] = (char)(1 - side[v]); queue[t++] = w; }
                else if (side[w] == side[v]) return 0;
            }
        }
    }
    return 1;
}

/* ---- Hopcroft-Karp ---------------------------------------------------------- */

#define GM_INF 0x3fffffff

static long gm_hopcroft_karp(const GalgUG* u, const char* side, int* mate, long size) {
    int n = u->n;
    int* dist = malloc((size_t)n * sizeof(int));
    int* queue = malloc((size_t)n * sizeof(int));
    int* it = malloc((size_t)n * sizeof(int));
    int* stk = malloc((size_t)n * sizeof(int));
    if (!dist || !queue || !it || !stk) { free(dist); free(queue); free(it); free(stk); return -1; }
    long phases = 0;
    for (;;) {
        if ((++phases & 15) == 0) galg_poll();
        /* BFS from free left vertices over alternating paths */
        int h = 0, t = 0, found = GM_INF;
        for (int v = 0; v < n; v++) {
            if (side[v] == 0 && mate[v] < 0) { dist[v] = 0; queue[t++] = v; }
            else dist[v] = GM_INF;
        }
        while (h < t) {
            int v = queue[h++];
            if (dist[v] >= found) continue;
            for (int j = u->off[v]; j < u->off[v + 1]; j++) {
                int w = u->adj[j], x = mate[w];
                if (x < 0) { if (found == GM_INF) found = dist[v] + 1; }
                else if (dist[x] == GM_INF) { dist[x] = dist[v] + 1; queue[t++] = x; }
            }
        }
        if (found == GM_INF) break;
        /* DFS (iterative) along dist layers from each free left vertex */
        for (int v = 0; v < n; v++) it[v] = u->off[v];
        for (int r = 0; r < n; r++) {
            if (side[r] != 0 || mate[r] >= 0 || dist[r] != 0) continue;
            int top = 0;
            stk[top++] = r;
            while (top > 0) {
                int v = stk[top - 1];
                int advanced = 0;
                while (it[v] < u->off[v + 1]) {
                    int w = u->adj[it[v]], x = mate[w];
                    if (x < 0) {
                        if (dist[v] + 1 == found) {
                            /* augment along the stack */
                            int cur = w;
                            for (int i = top - 1; i >= 0; i--) {
                                int a = stk[i], nxt = mate[a];
                                mate[a] = cur; mate[cur] = a;
                                cur = nxt;
                            }
                            size++;
                            top = 0; advanced = 1;
                            break;
                        }
                    } else if (dist[x] == dist[v] + 1) {
                        it[v]++;
                        stk[top++] = x; advanced = 1;
                        break;
                    }
                    it[v]++;
                }
                if (!advanced) { dist[v] = GM_INF; top--; }   /* dead: prune */
            }
        }
    }
    free(dist); free(queue); free(it); free(stk);
    return size;
}

/* ---- Edmonds' blossom algorithm -------------------------------------------- *
 * Blossom bases live in a union-find forest (dsu): contracting a blossom only
 * walks the two tree paths up to the base and links their bases under it, so
 * one search costs O(E alpha(V)) over the part of the graph it reaches. */

typedef struct {
    const GalgUG* u;
    int* mate;
    int* dsu;       /* blossom base = gm_find(v)                             */
    int* par;       /* tree parent links (through blossoms as well)          */
    signed char* label;  /* -1 unreached, 0 outer (even), 1 inner (odd)       */
    char* dead;     /* retired by a failed search                            */
    int* queue; int qt;
    int* touched; int ntouched;
    int* mark; int stamp;
} GmEd;

static int gm_find(int* dsu, int x) {
    while (dsu[x] != x) { dsu[x] = dsu[dsu[x]]; x = dsu[x]; }
    return x;
}

static void gm_label(GmEd* E, int v, int l) {
    if (E->label[v] < 0) E->touched[E->ntouched++] = v;
    E->label[v] = (signed char)l;
}

/* Nearest common base of the outer bases a and b in the search tree. */
static int gm_lca(GmEd* E, int a, int b) {
    E->stamp++;
    for (;;) {
        if (a >= 0) {
            if (E->mark[a] == E->stamp) return a;
            E->mark[a] = E->stamp;
            a = E->mate[a] < 0 ? -1 : gm_find(E->dsu, E->par[E->mate[a]]);
        }
        int t = a; a = b; b = t;
    }
}

/* Walk from outer vertex v up to base a, pointing par links across the new
 * edge (v, w), relabelling inner vertices outer and linking bases under a. */
static void gm_blossom(GmEd* E, int v, int w, int a) {
    while (gm_find(E->dsu, v) != a) {
        E->par[v] = w;
        w = E->mate[v];
        if (E->label[w] == 1) { E->label[w] = 0; E->queue[E->qt++] = w; }
        int bv = gm_find(E->dsu, v), bw = gm_find(E->dsu, w);
        E->dsu[bv] = a; E->dsu[bw] = a;
        v = E->par[w];
    }
}

/* One search from exposed root r: returns the exposed vertex that ends an
 * augmenting path (recover it through par/mate), or -1. */
static int gm_search(GmEd* E, int r) {
    const GalgUG* u = E->u;
    int qh = 0;
    E->qt = 0;
    gm_label(E, r, 0);
    E->queue[E->qt++] = r;
    while (qh < E->qt) {
        int v = E->queue[qh++];
        for (int j = u->off[v]; j < u->off[v + 1]; j++) {
            int x = u->adj[j];
            if (E->dead[x]) continue;
            if (E->label[x] < 0) {
                gm_label(E, x, 1);
                E->par[x] = v;
                if (E->mate[x] < 0) return x;
                gm_label(E, E->mate[x], 0);
                E->queue[E->qt++] = E->mate[x];
            } else if (E->label[x] == 0) {
                int bv = gm_find(E->dsu, v), bx = gm_find(E->dsu, x);
                if (bv == bx) continue;
                int a = gm_lca(E, bv, bx);
                gm_blossom(E, x, v, a);
                gm_blossom(E, v, x, a);
            }
        }
    }
    return -1;
}

static long gm_edmonds(const GalgUG* u, int* mate, long size) {
    int n = u->n;
    GmEd E;
    memset(&E, 0, sizeof(E));
    E.u = u; E.mate = mate;
    size_t nn = (size_t)(n > 0 ? n : 1);
    E.dsu = malloc(nn * sizeof(int));
    E.par = malloc(nn * sizeof(int));
    E.label = malloc(nn);
    E.dead = calloc(nn, 1);
    E.queue = malloc(nn * sizeof(int));
    E.touched = malloc(nn * sizeof(int));
    E.mark = calloc(nn, sizeof(int));
    if (!E.dsu || !E.par || !E.label || !E.dead || !E.queue || !E.touched || !E.mark) {
        size = -1; goto out;
    }
    for (int v = 0; v < n; v++) { E.dsu[v] = v; E.par[v] = -1; E.label[v] = -1; }
    for (int r = 0; r < n; r++) {
        if (mate[r] >= 0 || E.dead[r]) continue;
        if ((r & 255) == 0) galg_poll();
        int x = gm_search(&E, r);
        if (x >= 0) {
            while (x >= 0) {                     /* flip the augmenting path */
                int pv = E.par[x], nx = mate[pv];
                mate[x] = pv; mate[pv] = x;
                x = nx;
            }
            size++;
        } else {
            /* Hungarian tree: none of its vertices can lie on an augmenting
             * path later (Edmonds), so retire them all */
            for (int i = 0; i < E.ntouched; i++) E.dead[E.touched[i]] = 1;
        }
        for (int i = 0; i < E.ntouched; i++) {  /* reset what was touched */
            int v = E.touched[i];
            E.dsu[v] = v; E.par[v] = -1; E.label[v] = -1;
        }
        E.ntouched = 0;
    }
out:
    free(E.dsu); free(E.par); free(E.label); free(E.dead); free(E.queue);
    free(E.touched); free(E.mark);
    return size;
}

long galg_max_matching(const GalgUG* u, int* mate) {
    int n = u->n;
    long size = gm_karp_sipser(u, mate);
    if (size < 0 || n == 0) return size;
    char* side = malloc((size_t)n);
    int* queue = malloc((size_t)n * sizeof(int));
    if (!side || !queue) { free(side); free(queue); return -1; }
    int bip = gm_bipartite(u, side, queue);
    free(queue);
    if (bip) size = gm_hopcroft_karp(u, side, mate, size);
    else size = gm_edmonds(u, mate, size);
    free(side);
    return size;
}

/* ---- Builtins -------------------------------------------------------------- */

/* Matched edges of g in EdgeList order: for each matched pair, the first edge
 * of g joining it. */
static Expr* gm_matched_edges(const Expr* g, const int* mate) {
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    int n = galg_nv(g);
    long m = galg_ne(g), cnt = 0;
    char* used = calloc((size_t)(n > 0 ? n : 1), 1);
    long* idx = malloc((size_t)(m > 0 ? m : 1) * sizeof(long));
    Expr* out = NULL;
    if (used && idx) {
        for (long k = 0; k < m; k++) {
            int a = eu[k], b = ev[k];
            if (mate[a] == b && !used[a]) { used[a] = used[b] = 1; idx[cnt++] = k; }
        }
        out = galg_edge_list(g, idx, cnt);
    }
    free(used); free(idx);
    return out;
}

Expr* builtin_find_independent_edge_set(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    int* mate = malloc((size_t)(u->n > 0 ? u->n : 1) * sizeof(int));
    Expr* out = NULL;
    if (mate && galg_max_matching(u, mate) >= 0) out = gm_matched_edges(g, mate);
    free(mate); galg_ug_free(u);
    return out;
}

Expr* builtin_find_edge_cover(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    int n = u->n;
    for (int v = 0; v < n; v++)
        if (u->off[v + 1] == u->off[v]) {           /* isolated: no cover */
            galg_ug_free(u);
            return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        }
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    long m = galg_ne(g), cnt = 0;
    int* mate = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    char* cov = calloc((size_t)(n > 0 ? n : 1), 1);
    char* take = calloc((size_t)(m > 0 ? m : 1), 1);
    long* idx = malloc((size_t)(m > 0 ? m : 1) * sizeof(long));
    Expr* out = NULL;
    if (mate && cov && take && idx && galg_max_matching(u, mate) >= 0) {
        for (long k = 0; k < m; k++) {              /* the matching */
            int a = eu[k], b = ev[k];
            if (mate[a] == b && !cov[a]) { cov[a] = cov[b] = 1; take[k] = 1; }
        }
        for (long k = 0; k < m; k++) {              /* one edge per exposed vertex */
            int a = eu[k], b = ev[k];
            if (take[k]) continue;
            if (!cov[a] || !cov[b]) { take[k] = 1; cov[a] = cov[b] = 1; }
        }
        for (long k = 0; k < m; k++) if (take[k]) idx[cnt++] = k;
        out = galg_edge_list(g, idx, cnt);
    }
    free(mate); free(cov); free(take); free(idx); galg_ug_free(u);
    return out;
}
