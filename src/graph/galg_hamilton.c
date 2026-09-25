/* galg_hamilton.c - Hamiltonian cycles and paths.
 *
 *   FindHamiltonianCycle[g]          {c}: a Hamiltonian cycle as a list of
 *                                    edges, or {} if none exists
 *   FindHamiltonianCycle[g, n] / All up to n / all Hamiltonian cycles
 *   FindHamiltonianPath[g]           a Hamiltonian path as a vertex list, or {}
 *   FindHamiltonianPath[g, s, t]     one from s to t, or {}
 *   HamiltonianGraphQ[g]             True iff g has a Hamiltonian cycle
 *
 * Wolfram conventions reproduced (verified against Mathematica 15): directed
 * edges are followed forwards, undirected ones either way; the one-vertex
 * graph is Hamiltonian with the empty cycle {{}} and has the path {}; K2 has
 * no Hamiltonian cycle (a directed 2-cycle is one); a found cycle starts at
 * the first vertex of VertexList[g] and each edge is written in traversal
 * order, with the head of the edge of g it uses. Each cycle (as an edge set)
 * is reported once. Edge weights play no role (Mathematica returns {} for
 * FindHamiltonianPath on some weighted graphs that have one -- a defect this
 * implementation does not share).
 *
 * Algorithm: search over EDGE decisions with constraint propagation, the
 * formulation that makes small sparse instances easy (a random 400-vertex
 * cubic graph takes about 1 ms). It is not near-linear on large meshes: the
 * per-node biconnectivity check is linear, so GridGraph[{n, n}] grows roughly
 * quadratically in the vertex count (100 x 100 about 1.5 s, 150 x 150 about
 * 7 s, far larger grids impractical without TimeConstrained):
 *   - every vertex needs exactly two chosen edges (directed: one chosen in-arc
 *     and one chosen out-arc): a vertex with two chosen edges excludes the
 *     rest, a vertex with exactly two remaining edges forces both, fewer is a
 *     contradiction;
 *   - chosen edges form vertex-disjoint path fragments whose end-to-end links
 *     are kept, so an edge that would close a fragment into a cycle shorter
 *     than n is refused, and the edge joining a fragment's two ends is
 *     excluded as soon as it would do that;
 *   - at every search node the remaining (non-excluded) graph must be
 *     biconnected (directed: strongly connected) -- Tarjan, linear time;
 *   - branching extends a fragment end with the fewest remaining options,
 *     each option either chosen or excluded for the later siblings.
 * Every step is recorded on a trail and undone exactly, so the search is
 * complete: {} / False are proofs. Paths reduce to cycles through an added
 * vertex joined to the allowed endpoints.
 *
 * Budget: HAM_MAX_NODES search nodes, polled against the TimeConstrained
 * deadline; an exhausted budget leaves the call unevaluated -- never a wrong
 * "no".
 *
 * Memory (SPEC section 4): fresh results; scratch freed on every path.
 */

#include "graph.h"
#include "graph_algos.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

#define HAM_MAX_NODES 20000000L

/* ======================================================================== *
 * The engine. Undirected mode: items are edges {eu, ev}; directed mode: items
 * are arcs eu -> ev. Incidence lists: undirected -- inc of v lists every edge
 * at v (in off/inc); directed -- out-arcs (off/inc) and in-arcs (ioff/iinc).
 * ======================================================================== */

typedef struct { int e, a, b, oa, ob; char kind; } HamTr;

typedef struct {
    int n; long m;
    int directed;
    const int *eu, *ev;
    int *off, *inc;               /* undirected: all; directed: out-arcs     */
    int *ioff, *iinc;             /* directed: in-arcs                       */
    signed char* st;              /* 0 undecided, 1 chosen, -1 excluded      */
    int *chos, *avail;            /* undirected; directed: out counts        */
    int *ichos, *iavail;          /* directed: in counts                     */
    int* oth;                     /* fragment end <-> other end              */
    int nchosen;
    HamTr* tr; long tlen, tcap;
    int* q; char* inq; int qn;
    /* biconnectivity / reachability scratch */
    int *disc, *low, *stk, *itp, *par;
    /* results */
    long want, found;
    int* cycles; long ccap;
    long nodes;
    int aborted, stop;
} Ham;

static void ham_free(Ham* H) {
    free(H->off); free(H->inc); free(H->ioff); free(H->iinc); free(H->st);
    free(H->chos); free(H->avail); free(H->ichos); free(H->iavail); free(H->oth);
    free(H->tr); free(H->q); free(H->inq); free(H->disc); free(H->low); free(H->stk);
    free(H->itp); free(H->par); free(H->cycles);
}

static int ham_trail(Ham* H, char kind, int e, int a, int b, int oa, int ob) {
    if (H->tlen == H->tcap) {
        long nc = H->tcap ? 2 * H->tcap : 1024;
        HamTr* t = realloc(H->tr, (size_t)nc * sizeof(HamTr));
        if (!t) { H->aborted = 1; return 0; }
        H->tr = t; H->tcap = nc;
    }
    HamTr* x = &H->tr[H->tlen++];
    x->kind = kind; x->e = e; x->a = a; x->b = b; x->oa = oa; x->ob = ob;
    return 1;
}

static void ham_enq(Ham* H, int v) {
    if (!H->inq[v]) { H->inq[v] = 1; H->q[H->qn++] = v; }
}

static void ham_undo(Ham* H, long mark) {
    while (H->tlen > mark) {
        HamTr* x = &H->tr[--H->tlen];
        int e = x->e, u = H->eu[e], v = H->ev[e];
        H->st[e] = 0;
        if (x->kind == 'x') {
            if (H->directed) { H->avail[u]++; H->iavail[v]++; }
            else { H->avail[u]++; H->avail[v]++; }
        } else {
            if (H->directed) { H->chos[u]--; H->ichos[v]--; }
            else { H->chos[u]--; H->chos[v]--; }
            H->nchosen--;
            H->oth[x->a] = x->oa; H->oth[x->b] = x->ob;
        }
    }
}

static int ham_exclude(Ham* H, int e) {
    if (H->st[e] == -1) return 1;
    if (H->st[e] == 1) return 0;
    int u = H->eu[e], v = H->ev[e];
    H->st[e] = -1;
    if (H->directed) { H->avail[u]--; H->iavail[v]--; }
    else { H->avail[u]--; H->avail[v]--; }
    ham_enq(H, u); ham_enq(H, v);
    return ham_trail(H, 'x', e, 0, 0, 0, 0);
}

/* The undecided item joining x to y (x -> y when directed), or -1. */
static int ham_find(Ham* H, int x, int y) {
    for (int j = H->off[x]; j < H->off[x + 1]; j++) {
        int e = H->inc[j];
        if (H->st[e] != 0) continue;
        int o = H->directed ? H->ev[e] : (H->eu[e] == x ? H->ev[e] : H->eu[e]);
        if (o == y) return e;
    }
    return -1;
}

static int ham_choose(Ham* H, int e) {
    if (H->st[e] == 1) return 1;
    if (H->st[e] == -1) return 0;
    int u = H->eu[e], v = H->ev[e];
    int a, b;
    if (H->directed) {
        if (H->chos[u] || H->ichos[v]) return 0;
        a = H->oth[u];                 /* start of u's fragment (u is its end) */
        b = H->oth[v];                 /* end of v's fragment (v is its start)  */
        if (a == v && H->nchosen != H->n - 1) return 0;   /* short cycle */
    } else {
        if (H->chos[u] >= 2 || H->chos[v] >= 2) return 0;
        a = H->oth[u]; b = H->oth[v];
        if (a == v && H->nchosen != H->n - 1) return 0;
    }
    int closing = (a == v);
    H->st[e] = 1;
    if (H->directed) { H->chos[u]++; H->ichos[v]++; }
    else { H->chos[u]++; H->chos[v]++; }
    H->nchosen++;
    if (!ham_trail(H, 'c', e, a, b, H->oth[a], H->oth[b])) return 0;
    ham_enq(H, u); ham_enq(H, v);
    if (!closing) {
        H->oth[a] = b; H->oth[b] = a;
        if (H->nchosen < H->n - 1) {
            /* the item closing the new fragment early is now useless */
            int f = H->directed ? ham_find(H, b, a) : ham_find(H, a, b);
            if (f >= 0 && !ham_exclude(H, f)) return 0;
        }
        ham_enq(H, a); ham_enq(H, b);
    }
    return 1;
}

static int ham_propagate(Ham* H) {
    int ok = 1;
    while (H->qn > 0 && ok) {
        int v = H->q[--H->qn];
        H->inq[v] = 0;
        if (!H->directed) {
            if (H->avail[v] < 2) { ok = 0; break; }
            if (H->chos[v] == 2 && H->avail[v] > 2) {
                for (int j = H->off[v]; j < H->off[v + 1] && ok; j++)
                    if (H->st[H->inc[j]] == 0) ok = ham_exclude(H, H->inc[j]);
            } else if (H->avail[v] == 2 && H->chos[v] < 2) {
                for (int j = H->off[v]; j < H->off[v + 1] && ok; j++)
                    if (H->st[H->inc[j]] == 0) ok = ham_choose(H, H->inc[j]);
            }
        } else {
            if (H->avail[v] < 1 || H->iavail[v] < 1) { ok = 0; break; }
            if (H->chos[v] == 1 && H->avail[v] > 1) {
                for (int j = H->off[v]; j < H->off[v + 1] && ok; j++)
                    if (H->st[H->inc[j]] == 0) ok = ham_exclude(H, H->inc[j]);
            } else if (H->avail[v] == 1 && H->chos[v] == 0) {
                for (int j = H->off[v]; j < H->off[v + 1] && ok; j++)
                    if (H->st[H->inc[j]] == 0) ok = ham_choose(H, H->inc[j]);
            }
            if (!ok) break;
            if (H->ichos[v] == 1 && H->iavail[v] > 1) {
                for (int j = H->ioff[v]; j < H->ioff[v + 1] && ok; j++)
                    if (H->st[H->iinc[j]] == 0) ok = ham_exclude(H, H->iinc[j]);
            } else if (H->iavail[v] == 1 && H->ichos[v] == 0) {
                for (int j = H->ioff[v]; j < H->ioff[v + 1] && ok; j++)
                    if (H->st[H->iinc[j]] == 0) ok = ham_choose(H, H->iinc[j]);
            }
        }
    }
    if (!ok) while (H->qn > 0) H->inq[H->q[--H->qn]] = 0;
    return ok;
}

/* Undirected: is the non-excluded graph biconnected (connected, no cut
 * vertex)? Iterative Tarjan. */
static int ham_biconnected(Ham* H) {
    int n = H->n;
    if (n <= 2) return 1;
    for (int v = 0; v < n; v++) H->disc[v] = -1;
    int time = 0, top = 0, rootkids = 0;
    H->disc[0] = H->low[0] = time++; H->par[0] = -1; H->itp[0] = H->off[0];
    H->stk[top++] = 0;
    while (top > 0) {
        int v = H->stk[top - 1];
        if (H->itp[v] < H->off[v + 1]) {
            int e = H->inc[H->itp[v]++];
            if (H->st[e] == -1) continue;
            int w = H->eu[e] == v ? H->ev[e] : H->eu[e];
            if (H->disc[w] < 0) {
                H->disc[w] = H->low[w] = time++;
                H->par[w] = v; H->itp[w] = H->off[w];
                H->stk[top++] = w;
                if (v == 0) rootkids++;
            } else if (w != H->par[v] && H->disc[w] < H->low[v]) {
                H->low[v] = H->disc[w];
            }
        } else {
            top--;
            int p = H->par[v];
            if (p >= 0) {
                if (H->low[v] < H->low[p]) H->low[p] = H->low[v];
                if (p != 0 && H->low[v] >= H->disc[p]) return 0;   /* cut vertex */
            }
        }
    }
    if (time != n) return 0;
    return rootkids <= 1;
}

/* Directed: strongly connected over non-excluded arcs? */
static int ham_strong(Ham* H) {
    int n = H->n;
    for (int dir = 0; dir < 2; dir++) {
        int* off = dir ? H->ioff : H->off;
        int* inc = dir ? H->iinc : H->inc;
        for (int v = 0; v < n; v++) H->disc[v] = 0;
        int h = 0, t = 0;
        H->disc[0] = 1; H->stk[t++] = 0;
        while (h < t) {
            int x = H->stk[h++];
            for (int j = off[x]; j < off[x + 1]; j++) {
                int e = inc[j];
                if (H->st[e] == -1) continue;
                int y = dir ? H->eu[e] : H->ev[e];
                if (!H->disc[y]) { H->disc[y] = 1; H->stk[t++] = y; }
            }
        }
        if (t != n) return 0;
    }
    return 1;
}

static int ham_record(Ham* H) {
    int n = H->n;
    if ((H->found + 1) * n > H->ccap) {
        long nc = H->ccap ? 2 * H->ccap : (long)n * 4;
        int* t = realloc(H->cycles, (size_t)nc * sizeof(int));
        if (!t) { H->aborted = 1; return 0; }
        H->cycles = t; H->ccap = nc;
    }
    int* c = H->cycles + H->found * n;
    int prev = -1, cur = 0;
    for (int i = 0; i < n; i++) {
        c[i] = cur;
        int next = -1;
        for (int j = H->off[cur]; j < H->off[cur + 1]; j++) {
            int e = H->inc[j];
            if (H->st[e] != 1) continue;
            int o = H->directed ? H->ev[e] : (H->eu[e] == cur ? H->ev[e] : H->eu[e]);
            if (H->directed || o != prev || n == 2) { next = o; break; }
        }
        prev = cur; cur = next;
        if (cur < 0) break;
    }
    H->found++;
    if (H->want > 0 && H->found >= H->want) H->stop = 1;
    return 1;
}

static void ham_rec(Ham* H) {
    if (H->stop || H->aborted) return;
    if ((++H->nodes & 1023) == 0) {
        galg_poll();
        if (H->nodes > HAM_MAX_NODES) { H->aborted = 1; return; }
    }
    if (!ham_propagate(H)) return;
    if (H->nchosen == H->n) { ham_record(H); return; }
    if (!(H->directed ? ham_strong(H) : ham_biconnected(H))) return;
    /* branch vertex: a fragment end with the fewest undecided options (the
     * outgoing side when directed); else a vertex with the fewest options */
    int n = H->n, bv = -1, bopt = 0x7fffffff, bend = 0;
    for (int v = 0; v < n; v++) {
        int need, opt, isend;
        if (H->directed) {
            need = 1 - H->chos[v];
            if (need <= 0) continue;
            opt = H->avail[v] - H->chos[v];
            isend = H->ichos[v] == 1;
        } else {
            need = 2 - H->chos[v];
            if (need <= 0) continue;
            opt = H->avail[v] - H->chos[v];
            isend = H->chos[v] == 1;
        }
        if (isend > bend || (isend == bend && opt < bopt)) { bend = isend; bopt = opt; bv = v; }
    }
    if (bv < 0) return;
    /* binary branch on its first undecided item: choose it | exclude it
     * (the exclusion branch re-selects, so propagation that completes a
     * cycle between siblings is never skipped) */
    int e = -1;
    for (int j = H->off[bv]; j < H->off[bv + 1]; j++)
        if (H->st[H->inc[j]] == 0) { e = H->inc[j]; break; }
    if (e < 0) return;
    long m2 = H->tlen;
    if (ham_choose(H, e)) ham_rec(H);
    else while (H->qn > 0) H->inq[H->q[--H->qn]] = 0;
    ham_undo(H, m2);
    if (H->stop || H->aborted) return;
    m2 = H->tlen;
    if (ham_exclude(H, e)) ham_rec(H);
    else while (H->qn > 0) H->inq[H->q[--H->qn]] = 0;
    ham_undo(H, m2);
}

/* Build an instance over n vertices from items (eu[k], ev[k]). Returns 0 on
 * allocation failure. */
static int ham_init(Ham* H, int n, long m, const int* eu, const int* ev, int directed) {
    memset(H, 0, sizeof(*H));
    H->n = n; H->m = m; H->directed = directed; H->eu = eu; H->ev = ev;
    size_t nn = (size_t)n + 1, mm = (size_t)(m > 0 ? m : 1);
    H->off = calloc(nn + 1, sizeof(int));
    H->inc = malloc(2 * mm * sizeof(int));
    H->st = calloc(mm, 1);
    H->chos = calloc(nn, sizeof(int));
    H->avail = calloc(nn, sizeof(int));
    H->oth = malloc(nn * sizeof(int));
    H->q = malloc(nn * sizeof(int));
    H->inq = calloc(nn, 1);
    H->disc = malloc(nn * sizeof(int));
    H->low = malloc(nn * sizeof(int));
    H->stk = malloc(nn * sizeof(int));
    H->itp = malloc(nn * sizeof(int));
    H->par = malloc(nn * sizeof(int));
    if (directed) {
        H->ioff = calloc(nn + 1, sizeof(int));
        H->iinc = malloc(mm * sizeof(int));
        H->ichos = calloc(nn, sizeof(int));
        H->iavail = calloc(nn, sizeof(int));
        if (!H->ioff || !H->iinc || !H->ichos || !H->iavail) return 0;
    }
    if (!H->off || !H->inc || !H->st || !H->chos || !H->avail || !H->oth || !H->q || !H->inq
        || !H->disc || !H->low || !H->stk || !H->itp || !H->par) return 0;
    for (long k = 0; k < m; k++) {
        H->off[eu[k] + 1]++;
        if (directed) H->ioff[ev[k] + 1]++; else H->off[ev[k] + 1]++;
    }
    for (int v = 0; v < n; v++) {
        H->off[v + 1] += H->off[v];
        if (directed) H->ioff[v + 1] += H->ioff[v];
    }
    int* fo = malloc(nn * sizeof(int));
    int* fi = malloc(nn * sizeof(int));
    if (!fo || !fi) { free(fo); free(fi); return 0; }
    for (int v = 0; v < n; v++) { fo[v] = H->off[v]; if (directed) fi[v] = H->ioff[v]; }
    for (long k = 0; k < m; k++) {
        H->inc[fo[eu[k]]++] = (int)k;
        if (directed) H->iinc[fi[ev[k]]++] = (int)k;
        else H->inc[fo[ev[k]]++] = (int)k;
    }
    free(fo); free(fi);
    for (int v = 0; v < n; v++) {
        H->oth[v] = v;
        H->avail[v] = H->off[v + 1] - H->off[v];
        if (directed) H->iavail[v] = H->ioff[v + 1] - H->ioff[v];
    }
    return 1;
}

/* Run: returns 1 when the search completed (H->found cycles recorded), 0 if
 * it was aborted. */
static int ham_run(Ham* H, long want) {
    H->want = want;
    for (int v = 0; v < H->n; v++) ham_enq(H, v);
    ham_rec(H);
    return !H->aborted;
}

/* ======================================================================== *
 * Graph -> instance
 * ======================================================================== */

/* Items of g as used by the engine: for an all-undirected graph, its edges
 * (deduplicated pairs); otherwise arcs (both directions for undirected
 * edges, deduplicated). With virt, vertex n is added, joined to every vertex
 * (s < 0) or to s and t (z -> s and t -> z). *directed_out tells the mode.
 * The item arrays are malloc'd into *pu, *pv. Returns the item count, -1 on
 * failure. */
static long ham_items(const Expr* g, int virt, int s, int t, int** pu, int** pv, int* directed_out) {
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    int n = galg_nv(g);
    long m = galg_ne(g);
    int directed = graph_directed_edge_count(g) > 0;
    long cap = 2 * m + 2L * n + 8;
    int* a = malloc((size_t)cap * sizeof(int));
    int* b = malloc((size_t)cap * sizeof(int));
    if (!a || !b) { free(a); free(b); return -1; }
    long k = 0;
    for (long e = 0; e < m; e++) {
        if (!directed) { a[k] = eu[e] < ev[e] ? eu[e] : ev[e]; b[k] = eu[e] < ev[e] ? ev[e] : eu[e]; k++; }
        else {
            a[k] = eu[e]; b[k] = ev[e]; k++;
            if (!edir[e]) { a[k] = ev[e]; b[k] = eu[e]; k++; }
        }
    }
    if (virt) {
        int z = n;
        if (s < 0) {
            for (int v = 0; v < n; v++) {
                if (directed) { a[k] = z; b[k] = v; k++; a[k] = v; b[k] = z; k++; }
                else { a[k] = v; b[k] = z; k++; }
            }
        } else if (directed) {
            a[k] = z; b[k] = s; k++; a[k] = t; b[k] = z; k++;
        } else {
            a[k] = s; b[k] = z; k++; a[k] = t; b[k] = z; k++;
        }
    }
    /* deduplicate (a, b) pairs: sort by a then b via counting on a, and a
     * stamp per b */
    int nv = n + (virt ? 1 : 0);
    int* cnt = calloc((size_t)nv + 1, sizeof(int));
    int* ra = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
    int* rb = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
    int* stamp = malloc((size_t)nv * sizeof(int));
    if (!cnt || !ra || !rb || !stamp) { free(a); free(b); free(cnt); free(ra); free(rb); free(stamp); return -1; }
    for (long i = 0; i < k; i++) cnt[a[i] + 1]++;
    for (int v = 0; v < nv; v++) cnt[v + 1] += cnt[v];
    for (long i = 0; i < k; i++) { int p = cnt[a[i]]++; ra[p] = a[i]; rb[p] = b[i]; }
    for (int v = 0; v < nv; v++) stamp[v] = -1;
    long w = 0, i = 0;
    while (i < k) {
        int x = ra[i];
        while (i < k && ra[i] == x) {
            if (stamp[rb[i]] != x) { stamp[rb[i]] = x; a[w] = x; b[w] = rb[i]; w++; }
            i++;
        }
    }
    free(cnt); free(ra); free(rb); free(stamp);
    *pu = a; *pv = b; *directed_out = directed;
    return w;
}

/* Rotate cycle c (length n) to begin at vertex r; for undirected cycles
 * orient it so the second vertex is smaller than the last. */
static void ham_normalize(int* c, int n, int r, int undirected, int* tmp) {
    int p = 0;
    while (p < n && c[p] != r) p++;
    for (int i = 0; i < n; i++) tmp[i] = c[(p + i) % n];
    if (undirected && n >= 3 && tmp[1] > tmp[n - 1])
        for (int i = 1; i < n; i++) c[i] = tmp[n - i];
    else memcpy(c + 1, tmp + 1, (size_t)(n - 1) * sizeof(int));
    c[0] = r;
}

/* The edge list of cycle c in g: each edge in traversal order, with the head
 * UndirectedEdge when g joins the pair by an undirected edge, else
 * DirectedEdge. */
static Expr* ham_cycle_edges(const Expr* g, const int* c, int n) {
    Expr** items = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    if (!items) return NULL;
    Expr* hd = expr_new_symbol(SYM_DirectedEdge);
    Expr* hu = expr_new_symbol(SYM_UndirectedEdge);
    const Expr* verts = g->data.function.args[0];
    for (int i = 0; i < n; i++) {
        int a = c[i], b = c[(i + 1) % n];
        int und = graph_has_edge(g, verts->data.function.args[a], verts->data.function.args[b], 0) == 1;
        items[i] = galg_make_edge(g, und ? hu : hd, a, b);
    }
    expr_free(hd); expr_free(hu);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}

static Expr* ham_empty_list(void) { return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }

/* Two vertices: a Hamiltonian cycle needs an arc each way from two distinct
 * edges (so a single undirected edge -- K2 -- does not count). */
static Expr* ham_two(const Expr* g) {
    const int *eu, *ev;
    const unsigned char* edir;
    graph_edge_indices(g, &eu, &ev, &edir);
    long m = galg_ne(g);                   /* at most 3 edges on 2 vertices */
    for (long f = 0; f < m; f++) {
        if (edir[f] && eu[f] != 0) continue;          /* f must run 0 -> 1 */
        for (long r = 0; r < m; r++) {
            if (r == f || (edir[r] && eu[r] != 1)) continue;   /* r: 1 -> 0 */
            Expr* hf = expr_new_symbol(edir[f] ? SYM_DirectedEdge : SYM_UndirectedEdge);
            Expr* hr = expr_new_symbol(edir[r] ? SYM_DirectedEdge : SYM_UndirectedEdge);
            Expr* two[2] = { galg_make_edge(g, hf, 0, 1), galg_make_edge(g, hr, 1, 0) };
            expr_free(hf); expr_free(hr);
            Expr* ce = expr_new_function(expr_new_symbol(SYM_List), two, 2);
            return expr_new_function(expr_new_symbol(SYM_List), &ce, 1);
        }
    }
    return ham_empty_list();
}

/* Cycles of g (want <= 0: all). Returns the Wolfram list, NULL if aborted. */
static Expr* ham_cycles(const Expr* g, long want) {
    int n = galg_nv(g);
    if (n == 0) return ham_empty_list();
    if (n == 1) {
        Expr* e = ham_empty_list();
        return expr_new_function(expr_new_symbol(SYM_List), &e, 1);
    }
    if (n == 2) return ham_two(g);
    int *a = NULL, *b = NULL, directed;
    long k = ham_items(g, 0, -1, -1, &a, &b, &directed);
    if (k < 0) return NULL;
    Ham H;
    Expr* out = NULL;
    if (ham_init(&H, n, k, a, b, directed) && ham_run(&H, want)) {
        int* tmp = malloc((size_t)n * sizeof(int));
        Expr** items = malloc((size_t)(H.found > 0 ? H.found : 1) * sizeof(Expr*));
        if (tmp && items) {
            for (long i = 0; i < H.found; i++) {
                int* c = H.cycles + i * n;
                ham_normalize(c, n, 0, !directed, tmp);
                items[i] = ham_cycle_edges(g, c, n);
            }
            out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)H.found);
        }
        free(tmp); free(items);
    }
    ham_free(&H); free(a); free(b);
    return out;
}

Expr* builtin_find_hamiltonian_cycle(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    long want = 1;
    if (argc == 2 && !galg_parse_count(res->data.function.args[1], &want)) return NULL;
    return ham_cycles(g, want);
}

Expr* builtin_hamiltonian_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return galg_truth(0);
    Expr* c = ham_cycles(g, 1);
    if (!c) return NULL;
    int yes = c->data.function.arg_count > 0;
    expr_free(c);
    return galg_truth(yes);
}

Expr* builtin_find_hamiltonian_path(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int n = galg_nv(g);
    int s = -1, t = -1;
    if (argc == 3) {
        s = galg_vertex_arg(g, res->data.function.args[1]);
        t = galg_vertex_arg(g, res->data.function.args[2]);
        if (s < 0 || t < 0) return NULL;
        if (s == t) return ham_empty_list();
    }
    if (n <= 1) return ham_empty_list();
    int *a = NULL, *b = NULL, directed;
    long k = ham_items(g, 1, s, t, &a, &b, &directed);
    if (k < 0) return NULL;
    Ham H;
    Expr* out = NULL;
    /* n + 1 >= 3 vertices: the engine handles it; with n == 2 undirected and
     * s, t fixed the virtual vertex gives the triangle z-s-t */
    if (ham_init(&H, n + 1, k, a, b, directed) && ham_run(&H, 1)) {
        if (H.found == 0) out = ham_empty_list();
        else {
            int* c = H.cycles;             /* starts at vertex 0; rotate to z */
            int* p = malloc((size_t)(n + 1) * sizeof(int));
            if (p) {
                int zi = 0;
                while (c[zi] != n) zi++;
                for (int i = 0; i < n; i++) p[i] = c[(zi + 1 + i) % (n + 1)];
                /* undirected with fixed ends: read s .. t */
                if (s >= 0 && p[0] != s)
                    for (int i = 0, j = n - 1; i < j; i++, j--) { int x = p[i]; p[i] = p[j]; p[j] = x; }
                out = galg_vertex_list(g, p, n, 0);
            }
            free(p);
        }
    }
    ham_free(&H); free(a); free(b);
    return out;
}
