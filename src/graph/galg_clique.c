/* galg_clique.c - cliques.
 *
 *   FindClique[g]                    {c}: c a maximum clique
 *   FindClique[g, k]                 a largest maximal clique of at most k
 *                                    vertices (k may be Infinity)
 *   FindClique[g, {k}]               a maximal clique of exactly k vertices
 *   FindClique[g, {kmin, kmax}]      a largest maximal clique in the range
 *   FindClique[g, spec, n] / All     up to n / all maximal cliques meeting
 *                                    spec, largest first
 *   FindKClique[g, k]                {c}: a largest k-clique (vertices
 *                                    pairwise within distance k)
 *
 * Semantics (verified against Mathematica 15): the results of the spec forms
 * are MAXIMAL cliques (not extendable), so FindClique[g, 2] is {} when every
 * maximal clique has three vertices. In a directed graph a clique needs the
 * edges both ways (an UndirectedEdge counts as both). Each clique is listed in
 * VertexList order; several cliques come largest first, and within one size
 * Mathematica's order is reproduced (the n lexicographically smallest, listed
 * in reverse). FindClique[g] with no spec returns the clique the branch and
 * bound proves maximum -- Mathematica does not specify which maximum clique.
 *
 * Algorithms:
 *   - Maximum clique: bitset branch and bound with greedy-colouring bounds
 *     (Tomita's MCQ/MCS family in San Segundo's BBMC bitset form). Vertices
 *     are numbered by reverse degeneracy order. Graphs above GC_GLOBAL_MAX
 *     vertices are decomposed by degeneracy: every clique lies in
 *     {v} + N+(v) for its earliest vertex v, and |N+(v)| <= degeneracy, so
 *     each vertex gets a small local bitset problem (skipped outright when
 *     its core number or |N+(v)| cannot beat the incumbent).
 *   - Maximal cliques within a size range: Bron-Kerbosch with Tomita pivoting
 *     on bitsets, pruned by |R| + |P| < kmin and |R| >= kmax (a maximal clique
 *     through R would be too big).
 *   - FindKClique: a maximum clique of the k-th power graph (BFS to depth k).
 *
 * Budgets: GC_MAX_NODES search nodes, polled against the TimeConstrained
 * deadline every 4096 nodes; an exhausted budget leaves the call unevaluated
 * (never a non-maximum answer).
 *
 * Memory (SPEC section 4): fresh results; scratch freed on every path.
 */

#include "graph.h"
#include "graph_algos.h"
#include "expr.h"
#include "sym_names.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define GC_GLOBAL_MAX 3000        /* one global bitset problem up to this n    */
#define GC_BITSET_MAX 8192        /* largest n for complement / power bitsets  */
#define GC_MAX_NODES 50000000L    /* branch-and-bound / enumeration node budget */
#define GC_MAX_RESULT_INTS 20000000L

typedef uint64_t GcW;

static int gc_popcount(GcW x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0;
    while (x) { x &= x - 1; c++; }
    return c;
#endif
}

static int gc_ctz(GcW x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int c = 0;
    while (!(x & 1)) { x >>= 1; c++; }
    return c;
#endif
}

/* ---- Bitset graphs ---------------------------------------------------------- */

typedef struct {
    int n, w;          /* vertices, words per row */
    GcW* rows;         /* n * w                    */
    int* label;        /* local vertex -> caller's vertex id */
} GcBits;

static void gc_bits_free(GcBits* b) {
    if (!b) return;
    free(b->rows); free(b->label); free(b);
}

static GcBits* gc_bits_new(int n) {
    GcBits* b = calloc(1, sizeof(GcBits));
    if (!b) return NULL;
    b->n = n; b->w = (n + 63) / 64;
    if (b->w == 0) b->w = 1;
    b->rows = calloc((size_t)n * (size_t)b->w + 1, sizeof(GcW));
    b->label = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!b->rows || !b->label) { gc_bits_free(b); return NULL; }
    return b;
}

#define GC_ROW(b, i) ((b)->rows + (size_t)(i) * (size_t)(b)->w)
#define GC_SET(r, j) ((r)[(j) >> 6] |= (GcW)1 << ((j) & 63))
#define GC_TEST(r, j) (((r)[(j) >> 6] >> ((j) & 63)) & 1)

/* Induced bitset graph on vs[0..k-1] (in that local order) of u; complement
 * != 0 builds the complement instead. pos is scratch of size u->n, all -1 on
 * entry and restored on exit. */
static GcBits* gc_bits_induced(const GalgUG* u, const int* vs, int k, int complement, int* pos) {
    GcBits* b = gc_bits_new(k);
    if (!b) return NULL;
    for (int i = 0; i < k; i++) { pos[vs[i]] = i; b->label[i] = vs[i]; }
    for (int i = 0; i < k; i++) {
        GcW* r = GC_ROW(b, i);
        int v = vs[i];
        for (int j = u->off[v]; j < u->off[v + 1]; j++) {
            int p = pos[u->adj[j]];
            if (p >= 0) GC_SET(r, p);
        }
        if (complement) {
            for (int q = 0; q < b->w; q++) r[q] = ~r[q];
            int tail = k & 63;
            if (tail) r[b->w - 1] &= ((GcW)1 << tail) - 1;
            if (k == 0) r[0] = 0;
            r[i >> 6] &= ~((GcW)1 << (i & 63));
        }
    }
    for (int i = 0; i < k; i++) pos[vs[i]] = -1;
    return b;
}

/* ---- Degeneracy order ------------------------------------------------------- */

/* order[0..n-1]: vertices in removal order of the min-degree peeling (so the
 * densest core comes last); core[v] = core number. O(n + m) bucket queue. */
static int gc_degeneracy(const GalgUG* u, int* order, int* core) {
    int n = u->n, maxd = 0;
    int* deg = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* bin = NULL; int* pos = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* vert = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!deg || !pos || !vert) { free(deg); free(pos); free(vert); return 0; }
    for (int v = 0; v < n; v++) { deg[v] = u->off[v + 1] - u->off[v]; if (deg[v] > maxd) maxd = deg[v]; }
    bin = calloc((size_t)maxd + 2, sizeof(int));
    if (!bin) { free(deg); free(pos); free(vert); return 0; }
    for (int v = 0; v < n; v++) bin[deg[v]]++;
    for (int d = 0, s = 0; d <= maxd; d++) { int c = bin[d]; bin[d] = s; s += c; }
    for (int v = 0; v < n; v++) { pos[v] = bin[deg[v]]; vert[pos[v]] = v; bin[deg[v]]++; }
    for (int d = maxd; d > 0; d--) bin[d] = bin[d - 1];
    bin[0] = 0;
    for (int i = 0; i < n; i++) {
        int v = vert[i];
        order[i] = v; core[v] = deg[v];
        for (int j = u->off[v]; j < u->off[v + 1]; j++) {
            int w = u->adj[j];
            if (deg[w] > deg[v]) {
                int dw = deg[w], pw = pos[w], pstart = bin[dw], x = vert[pstart];
                if (x != w) { pos[w] = pstart; vert[pw] = x; pos[x] = pw; vert[pstart] = w; }
                bin[dw]++; deg[w]--;
            }
        }
    }
    free(deg); free(bin); free(pos); free(vert);
    return 1;
}

/* ---- Maximum clique: BBMC ---------------------------------------------------- */

typedef struct {
    const GcBits* b;
    int best;           /* size of the incumbent */
    int* bestc;         /* incumbent (local ids) */
    int* cur;           /* current clique        */
    GcW* pstack;        /* (depth) * w           */
    /* per-depth scratch, allocated on first use (depth <= clique size, so
     * this stays O(omega * n) rather than O(n^2)) */
    GcW** qbuf;         /* 2w: colouring work sets   */
    int** ubuf;         /* n: vertices in colour order */
    int** cbuf;         /* n: their colours          */
    long nodes;
    int aborted;
} GcMc;

static void gc_mc_expand(GcMc* S, int depth, GcW* P) {
    const GcBits* b = S->b;
    int w = b->w, n = b->n;
    if (S->aborted) return;
    if ((++S->nodes & 4095) == 0) {
        galg_poll();
        if (S->nodes > GC_MAX_NODES) { S->aborted = 1; return; }
    }
    if (!S->qbuf[depth]) {
        S->qbuf[depth] = malloc(2 * (size_t)w * sizeof(GcW));
        S->ubuf[depth] = malloc((size_t)n * sizeof(int));
        S->cbuf[depth] = malloc((size_t)n * sizeof(int));
        if (!S->qbuf[depth] || !S->ubuf[depth] || !S->cbuf[depth]) { S->aborted = 1; return; }
    }
    GcW* Q = S->qbuf[depth];
    GcW* Qk = Q + w;
    int* U = S->ubuf[depth];
    int* col = S->cbuf[depth];
    memcpy(Q, P, (size_t)w * sizeof(GcW));
    int kmin = S->best - depth;   /* need depth + colour > best */
    if (kmin < 1) kmin = 1;
    int m = 0, k = 0;
    for (;;) {
        int any = 0;
        for (int q = 0; q < w; q++) if (Q[q]) { any = 1; break; }
        if (!any) break;
        k++;
        memcpy(Qk, Q, (size_t)w * sizeof(GcW));
        for (int q = 0; q < w; q++) {
            while (Qk[q]) {
                int v = q * 64 + gc_ctz(Qk[q]);
                Qk[q] &= Qk[q] - 1;
                Q[q] &= ~((GcW)1 << (v & 63));
                const GcW* r = GC_ROW(b, v);
                for (int t = q; t < w; t++) Qk[t] &= ~r[t];
                if (k >= kmin) { U[m] = v; col[m] = k; m++; }
            }
        }
    }
    GcW* NP = S->pstack + (size_t)(depth + 1) * (size_t)w;
    for (int i = m - 1; i >= 0; i--) {
        if (depth + col[i] <= S->best) return;
        int v = U[i];
        S->cur[depth] = v;
        const GcW* r = GC_ROW(b, v);
        int any = 0;
        for (int q = 0; q < w; q++) { NP[q] = P[q] & r[q]; any |= NP[q] != 0; }
        if (!any) {
            if (depth + 1 > S->best) {
                S->best = depth + 1;
                memcpy(S->bestc, S->cur, (size_t)(depth + 1) * sizeof(int));
            }
        } else {
            gc_mc_expand(S, depth + 1, NP);
            if (S->aborted) return;
        }
        P[v >> 6] &= ~((GcW)1 << (v & 63));
    }
}

/* Maximum clique of b larger than lb (lb = size already known). Writes it to
 * out (local ids) and returns its size, or lb if none is larger, or -1 when
 * the budget ran out. *nodes accumulates the node count. */
static int gc_maxclique(const GcBits* b, int lb, int* out, long* nodes) {
    int n = b->n, w = b->w;
    if (n == 0) return lb;
    GcMc S;
    memset(&S, 0, sizeof(S));
    S.b = b; S.best = lb;
    S.bestc = malloc((size_t)n * sizeof(int));
    S.cur = malloc((size_t)n * sizeof(int));
    S.pstack = malloc((size_t)(n + 2) * (size_t)w * sizeof(GcW));
    S.qbuf = calloc((size_t)n + 1, sizeof(GcW*));
    S.ubuf = calloc((size_t)n + 1, sizeof(int*));
    S.cbuf = calloc((size_t)n + 1, sizeof(int*));
    S.nodes = *nodes;
    int ret = -1;
    if (S.bestc && S.cur && S.pstack && S.qbuf && S.ubuf && S.cbuf) {
        GcW* P = S.pstack;
        memset(P, 0, (size_t)w * sizeof(GcW));
        for (int i = 0; i < n; i++) GC_SET(P, i);
        gc_mc_expand(&S, 0, P);
        if (!S.aborted) {
            ret = S.best;
            if (S.best > lb) memcpy(out, S.bestc, (size_t)S.best * sizeof(int));
        }
    }
    *nodes = S.nodes;
    for (int d = 0; d <= n; d++) {
        if (S.qbuf) free(S.qbuf[d]);
        if (S.ubuf) free(S.ubuf[d]);
        if (S.cbuf) free(S.cbuf[d]);
    }
    free(S.bestc); free(S.cur); free(S.pstack); free(S.qbuf); free(S.ubuf); free(S.cbuf);
    return ret;
}

/* Maximum clique of the simple graph u (vertex ids of u). Returns the size and
 * fills out (sorted ascending), or -1 on budget exhaustion / allocation
 * failure. */
int galg_max_clique(const GalgUG* u, int* out) {
    int n = u->n;
    if (n == 0) return 0;
    int* order = calloc((size_t)n, sizeof(int));
    int* core = malloc((size_t)n * sizeof(int));
    int* pos = malloc((size_t)n * sizeof(int));
    int* vs = malloc((size_t)n * sizeof(int));
    int* loc = malloc((size_t)n * sizeof(int));
    int best = -1;
    long nodes = 0;
    if (!order || !core || !pos || !vs || !loc || !gc_degeneracy(u, order, core)) goto out;
    for (int v = 0; v < n; v++) pos[v] = -1;
    best = 1;
    out[0] = order[n - 1];
    if (n <= GC_GLOBAL_MAX) {
        /* one problem, numbered by reverse degeneracy order */
        for (int i = 0; i < n; i++) vs[i] = order[n - 1 - i];
        /* only vertices with core >= best can be in a larger clique; all are
         * kept here (best = 1) */
        GcBits* b = gc_bits_induced(u, vs, n, 0, pos);
        if (!b) { best = -1; goto out; }
        int r = gc_maxclique(b, 0, loc, &nodes);
        if (r < 0) best = -1;
        else if (r > 0) { best = r; for (int i = 0; i < r; i++) out[i] = b->label[loc[i]]; }
        gc_bits_free(b);
    } else {
        int* rank = malloc((size_t)n * sizeof(int));
        if (!rank) { best = -1; goto out; }
        for (int i = 0; i < n; i++) rank[order[i]] = i;
        /* densest cores first, so the incumbent grows early */
        for (int i = n - 1; i >= 0 && best >= 0; i--) {
            int v = order[i];
            if (core[v] + 1 <= best) continue;
            int k = 0;
            for (int j = u->off[v]; j < u->off[v + 1]; j++) {
                int x = u->adj[j];
                if (rank[x] > i && core[x] + 1 > best) vs[k++] = x;
            }
            if (k + 1 <= best) continue;
            /* local numbering: by descending core (reverse degeneracy rank) */
            for (int a = 1; a < k; a++) {
                int x = vs[a], b2 = a - 1;
                while (b2 >= 0 && rank[vs[b2]] < rank[x]) { vs[b2 + 1] = vs[b2]; b2--; }
                vs[b2 + 1] = x;
            }
            GcBits* b = gc_bits_induced(u, vs, k, 0, pos);
            if (!b) { best = -1; break; }
            int r = gc_maxclique(b, best - 1, loc, &nodes);
            if (r < 0) best = -1;
            else if (r > best - 1) {
                best = r + 1;
                out[0] = v;
                for (int t = 0; t < r; t++) out[t + 1] = b->label[loc[t]];
            }
            gc_bits_free(b);
        }
        free(rank);
    }
    if (best > 0) {
        /* sort ascending (small k: insertion sort) */
        for (int a = 1; a < best; a++) {
            int x = out[a], b2 = a - 1;
            while (b2 >= 0 && out[b2] > x) { out[b2 + 1] = out[b2]; b2--; }
            out[b2 + 1] = x;
        }
    }
out:
    free(order); free(core); free(pos); free(vs); free(loc);
    return best;
}

/* Maximum independent set of the subgraph of u induced by vs[0..k-1]: the
 * maximum clique of its complement, by the same branch and bound. Fills out
 * with vertex ids of u; returns the size or -1. Used by galg_mis.c for dense
 * components. */
int galg_mis_dense(const GalgUG* u, const int* vs, int k, int* out) {
    int n = u->n;
    int* pos = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* loc = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
    int r = -1;
    if (pos && loc) {
        for (int v = 0; v < n; v++) pos[v] = -1;
        /* number by ascending degree in g = descending degree in the
         * complement, the order the colouring bound likes */
        int* ord = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
        if (ord) {
            memcpy(ord, vs, (size_t)k * sizeof(int));
            for (int a = 1; a < k; a++) {
                int x = ord[a], dx = u->off[x + 1] - u->off[x], b2 = a - 1;
                while (b2 >= 0 && u->off[ord[b2] + 1] - u->off[ord[b2]] > dx) { ord[b2 + 1] = ord[b2]; b2--; }
                ord[b2 + 1] = x;
            }
            GcBits* b = gc_bits_induced(u, ord, k, 1, pos);
            if (b) {
                long nodes = 0;
                r = gc_maxclique(b, 0, loc, &nodes);
                if (r > 0) for (int i = 0; i < r; i++) out[i] = b->label[loc[i]];
                gc_bits_free(b);
            }
            free(ord);
        }
    }
    free(pos); free(loc);
    return r;
}

/* ---- Maximal cliques in a size range: pivoted Bron-Kerbosch ----------------- */

typedef struct {
    const GcBits* b;
    int kmin, kmax;
    int* R;
    GcW* buf;          /* per depth: P and X (2w)            */
    long nodes;
    int aborted;
    /* results: concatenated sorted lists */
    int* res; long rlen, rcap;
    long* rstart; long rcount, rscap;
} GcBk;

static int gc_bk_push(GcBk* S, int size) {
    if (S->rlen + size > S->rcap) {
        long nc = S->rcap ? 2 * S->rcap : 1024;
        while (nc < S->rlen + size) nc *= 2;
        if (nc > GC_MAX_RESULT_INTS) return 0;
        int* t = realloc(S->res, (size_t)nc * sizeof(int));
        if (!t) return 0;
        S->res = t; S->rcap = nc;
    }
    if (S->rcount + 1 > S->rscap) {
        long nc = S->rscap ? 2 * S->rscap : 256;
        long* t = realloc(S->rstart, (size_t)(nc + 1) * sizeof(long));
        if (!t) return 0;
        S->rstart = t; S->rscap = nc;
    }
    /* labels are ascending in local order only if the caller numbered so;
     * sort the reported clique by caller id */
    int* dst = S->res + S->rlen;
    for (int i = 0; i < size; i++) {
        int x = S->b->label[S->R[i]], j = i - 1;
        while (j >= 0 && dst[j] > x) { dst[j + 1] = dst[j]; j--; }
        dst[j + 1] = x;
    }
    S->rstart[S->rcount++] = S->rlen;
    S->rlen += size;
    S->rstart[S->rcount] = S->rlen;
    return 1;
}

static void gc_bk(GcBk* S, int depth, GcW* P, GcW* X) {
    const GcBits* b = S->b;
    int w = b->w;
    if (S->aborted) return;
    if ((++S->nodes & 4095) == 0) {
        galg_poll();
        if (S->nodes > GC_MAX_NODES) { S->aborted = 1; return; }
    }
    int pc = 0, xany = 0;
    for (int q = 0; q < w; q++) { pc += gc_popcount(P[q]); xany |= X[q] != 0; }
    if (pc == 0) {
        if (!xany && depth >= S->kmin && depth <= S->kmax)
            if (!gc_bk_push(S, depth)) S->aborted = 1;
        return;
    }
    if (depth + pc < S->kmin) return;
    if (depth >= S->kmax) return;          /* any maximal extension is too big */
    /* pivot: u in P u X maximizing |P n N(u)| */
    int piv = -1, bestc = -1;
    for (int q = 0; q < w; q++) {
        GcW m = P[q] | X[q];
        while (m) {
            int u = q * 64 + gc_ctz(m);
            m &= m - 1;
            const GcW* r = GC_ROW(b, u);
            int c = 0;
            for (int t = 0; t < w; t++) c += gc_popcount(P[t] & r[t]);
            if (c > bestc) { bestc = c; piv = u; }
        }
    }
    GcW* NP = S->buf + (size_t)(depth + 1) * 2 * (size_t)w;
    GcW* NX = NP + w;
    const GcW* pr = GC_ROW(b, piv);
    for (int q = 0; q < w; q++) {
        GcW cand = P[q] & ~pr[q];
        while (cand) {
            int v = q * 64 + gc_ctz(cand);
            cand &= cand - 1;
            const GcW* r = GC_ROW(b, v);
            for (int t = 0; t < w; t++) { NP[t] = P[t] & r[t]; NX[t] = X[t] & r[t]; }
            S->R[depth] = v;
            gc_bk(S, depth + 1, NP, NX);
            if (S->aborted) return;
            P[q] &= ~((GcW)1 << (v & 63));
            X[q] |= (GcW)1 << (v & 63);
        }
    }
}

/* Lexicographic compare of two sorted int lists (shorter-prefix first). */
static int gc_lexcmp(const int* a, int la, const int* b, int lb) {
    int k = la < lb ? la : lb;
    for (int i = 0; i < k; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return la < lb ? -1 : la > lb;
}

static const GcBk* gc_sort_ctx;
static int gc_res_cmp(const void* pa, const void* pb) {
    long a = *(const long*)pa, b = *(const long*)pb;
    const GcBk* S = gc_sort_ctx;
    int la = (int)(S->rstart[a + 1] - S->rstart[a]), lb = (int)(S->rstart[b + 1] - S->rstart[b]);
    if (la != lb) return la > lb ? -1 : 1;                  /* larger first */
    return gc_lexcmp(S->res + S->rstart[a], la, S->res + S->rstart[b], lb);
}

/* The Wolfram answer for a spec query on the bitset graph b whose labels are
 * vertex positions of g: up to `count` (<= 0: all) maximal cliques with size
 * in [kmin, kmax], largest first; if only_largest, only those of the largest
 * size present are candidates. Returns NULL on budget exhaustion. */
static Expr* gc_spec_answer(const Expr* g, const GcBits* b, int kmin, int kmax, long count) {
    int w = b->w, n = b->n;
    GcBk S;
    memset(&S, 0, sizeof(S));
    S.b = b; S.kmin = kmin; S.kmax = kmax;
    S.R = malloc((size_t)(n + 1) * sizeof(int));
    S.buf = malloc((size_t)(n + 2) * 2 * (size_t)w * sizeof(GcW));
    Expr* out = NULL;
    long* idx = NULL;
    if (!S.R || !S.buf) goto done;
    {
        GcW* P = S.buf; GcW* X = P + w;
        memset(P, 0, 2 * (size_t)w * sizeof(GcW));
        for (int i = 0; i < n; i++) GC_SET(P, i);
        if (n > 0) gc_bk(&S, 0, P, X);
    }
    if (S.aborted) goto done;
    idx = malloc((size_t)(S.rcount > 0 ? S.rcount : 1) * sizeof(long));
    if (!idx) goto done;
    for (long i = 0; i < S.rcount; i++) idx[i] = i;
    gc_sort_ctx = &S;
    qsort(idx, (size_t)S.rcount, sizeof(long), gc_res_cmp);
    {
        long take = S.rcount;
        if (count > 0 && count < take) take = count;
        /* reverse each run of equal size (Mathematica's listing order) */
        long i = 0;
        while (i < take) {
            long j = i;
            long sz = S.rstart[idx[i] + 1] - S.rstart[idx[i]];
            while (j < take && S.rstart[idx[j] + 1] - S.rstart[idx[j]] == sz) j++;
            for (long a = i, c = j - 1; a < c; a++, c--) { long t = idx[a]; idx[a] = idx[c]; idx[c] = t; }
            i = j;
        }
        Expr** items = malloc((size_t)(take > 0 ? take : 1) * sizeof(Expr*));
        if (!items) goto done;
        for (long t = 0; t < take; t++) {
            long r = idx[t];
            items[t] = galg_vertex_list(g, S.res + S.rstart[r],
                                        (int)(S.rstart[r + 1] - S.rstart[r]), 0);
        }
        out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)take);
        free(items);
    }
done:
    free(S.R); free(S.buf); free(S.res); free(S.rstart); free(idx);
    return out;
}

/* ---- Argument parsing shared with FindIndependentVertexSet ------------------ */

/* Parse a size spec: k (Integer >= 0 or Infinity) -> [1, k]; {k} -> [k, k];
 * {kmin, kmax}. Returns 1 on success. */
static int gc_parse_int_or_inf(const Expr* e, long* v) {
    if (e->type == EXPR_INTEGER && e->data.integer >= 0) { *v = (long)e->data.integer; return 1; }
    if (galg_is_symbol(e, "Infinity")) { *v = 0x3fffffff; return 1; }
    return 0;
}

int galg_parse_size_spec(const Expr* spec, int* kmin, int* kmax) {
    long a, b;
    if (gc_parse_int_or_inf(spec, &a)) { *kmin = 0; *kmax = (int)a; return 1; }
    if (graph_is_list(spec)) {
        size_t c = spec->data.function.arg_count;
        if (c == 1 && gc_parse_int_or_inf(spec->data.function.args[0], &a)) {
            *kmin = *kmax = (int)a; return 1;
        }
        if (c == 2 && gc_parse_int_or_inf(spec->data.function.args[0], &a)
            && gc_parse_int_or_inf(spec->data.function.args[1], &b) && a <= b) {
            *kmin = (int)a; *kmax = (int)b; return 1;
        }
    }
    return 0;
}

/* count argument: positive Integer or All (-> 0). */
int galg_parse_count(const Expr* e, long* count) {
    if (e->type == EXPR_INTEGER && e->data.integer >= 1) { *count = (long)e->data.integer; return 1; }
    if (galg_is_symbol(e, "All") || galg_is_symbol(e, "Infinity")) { *count = 0; return 1; }
    return 0;
}

Expr* galg_clique_spec_query(const Expr* g, const GalgUG* u, int complement,
                             int kmin, int kmax, long count) {
    int n = u->n;
    if (n > GC_BITSET_MAX) return NULL;
    if (kmin < 1) kmin = 1;
    if (kmax < kmin) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    int* vs = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* pos = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    Expr* out = NULL;
    if (vs && pos) {
        for (int v = 0; v < n; v++) { vs[v] = v; pos[v] = -1; }
        GcBits* b = gc_bits_induced(u, vs, n, complement, pos);
        if (b) {
            if (count == 1) {
                /* the first answer is a largest one: try sizes downward, so
                 * the pruning bound kmin is as strong as possible */
                int hi = kmax < n ? kmax : n, failed = 0;
                for (int s = hi; s >= kmin && !out; s--) {
                    Expr* r = gc_spec_answer(g, b, s, s, 1);
                    if (!r) { failed = 1; break; }          /* budget */
                    if (r->data.function.arg_count > 0) out = r;
                    else expr_free(r);
                }
                if (!out && !failed) out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            } else {
                out = gc_spec_answer(g, b, kmin, kmax, count);
            }
            gc_bits_free(b);
        }
    }
    free(vs); free(pos);
    return out;
}

/* ---- The mutual-adjacency graph for directed cliques ------------------------ */

/* Simple graph whose edges are the pairs joined both ways in g: an undirected
 * edge, or directed edges in both directions. */
static GalgUG* gc_mutual_graph(const Expr* g) {
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return NULL;
    int n = galg_nv(g);
    long m = galg_ne(g);
    if (graph_directed_edge_count(g) == 0) return galg_ug_from_edges(n, m, eu, ev);
    /* arcs a->b; a pair is mutual iff both arcs exist */
    long na = 0;
    for (long k = 0; k < m; k++) na += edir[k] ? 1 : 2;
    int* off = calloc((size_t)n + 1, sizeof(int));
    int* adj = malloc((size_t)(na > 0 ? na : 1) * sizeof(int));
    int* fill = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* ku = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
    int* kv = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
    GalgUG* res = NULL;
    if (off && adj && fill && ku && kv) {
        for (long k = 0; k < m; k++) { off[eu[k] + 1]++; if (!edir[k]) off[ev[k] + 1]++; }
        for (int v = 0; v < n; v++) off[v + 1] += off[v];
        for (int v = 0; v < n; v++) fill[v] = off[v];
        for (long k = 0; k < m; k++) {
            adj[fill[eu[k]]++] = ev[k];
            if (!edir[k]) adj[fill[ev[k]]++] = eu[k];
        }
        /* mark-based mutual test: for each a, stamp its out-neighbours, then
         * each b in out(a) with a in out(b) forms a mutual pair (a < b) */
        long cnt = 0;
        int* stamp = fill;
        for (int v = 0; v < n; v++) stamp[v] = -1;
        for (int a = 0; a < n; a++) {
            for (int j = off[a]; j < off[a + 1]; j++) {
                int b = adj[j];
                if (b <= a) continue;
                int mutual = 0;
                for (int t = off[b]; t < off[b + 1]; t++) if (adj[t] == a) { mutual = 1; break; }
                if (mutual && stamp[b] != a) { stamp[b] = a; ku[cnt] = a; kv[cnt] = b; cnt++; }
            }
        }
        res = galg_ug_from_edges(n, cnt, ku, kv);
    }
    free(off); free(adj); free(fill); free(ku); free(kv);
    return res;
}

/* ---- Builtins --------------------------------------------------------------- */

Expr* builtin_find_clique(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int kmin = 1, kmax = 0x3fffffff;
    long count = 1;
    if (argc >= 2 && !galg_parse_size_spec(res->data.function.args[1], &kmin, &kmax)) return NULL;
    if (argc == 3 && !galg_parse_count(res->data.function.args[2], &count)) return NULL;
    GalgUG* u = gc_mutual_graph(g);
    if (!u) return NULL;
    Expr* out = NULL;
    if (u->n == 0) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else if (argc == 1) {
        int* c = malloc((size_t)u->n * sizeof(int));
        if (c) {
            int k = galg_max_clique(u, c);
            if (k > 0) {
                Expr* one = galg_vertex_list(g, c, k, 0);
                out = expr_new_function(expr_new_symbol(SYM_List), &one, 1);
            }
        }
        free(c);
    } else {
        out = galg_clique_spec_query(g, u, 0, kmin, kmax, count);
    }
    galg_ug_free(u);
    return out;
}

Expr* builtin_find_k_clique(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* ka = res->data.function.args[1];
    if (!graph_is_valid(g) || ka->type != EXPR_INTEGER || ka->data.integer < 1) return NULL;
    long K = (long)ka->data.integer;
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    int n = u->n;
    Expr* out = NULL;
    if (n == 0) { galg_ug_free(u); return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }
    if (n > GC_BITSET_MAX) { galg_ug_free(u); return NULL; }
    /* power graph by truncated BFS from every vertex */
    int* dist = malloc((size_t)n * sizeof(int));
    int* q = malloc((size_t)n * sizeof(int));
    long cap = 1024, cnt = 0;
    int* pu = malloc((size_t)cap * sizeof(int));
    int* pv = malloc((size_t)cap * sizeof(int));
    int ok = dist && q && pu && pv;
    for (int v = 0; v < n && ok; v++) dist[v] = -1;
    for (int s = 0; s < n && ok; s++) {
        if ((s & 63) == 0) galg_poll();
        int h = 0, t = 0;
        dist[s] = 0; q[t++] = s;
        while (h < t) {
            int x = q[h++];
            if (dist[x] >= K) continue;
            for (int j = u->off[x]; j < u->off[x + 1]; j++) {
                int y = u->adj[j];
                if (dist[y] < 0) { dist[y] = dist[x] + 1; q[t++] = y; }
            }
        }
        for (int i = 0; i < t; i++) {
            int y = q[i];
            if (y > s) {
                if (cnt == cap) {
                    cap *= 2;
                    int* a = realloc(pu, (size_t)cap * sizeof(int));
                    if (a) pu = a;
                    int* b = realloc(pv, (size_t)cap * sizeof(int));
                    if (b) pv = b;
                    if (!a || !b) { ok = 0; break; }
                }
                pu[cnt] = s; pv[cnt] = y; cnt++;
            }
            dist[y] = -1;
        }
    }
    if (ok) {
        GalgUG* pw = galg_ug_from_edges(n, cnt, pu, pv);
        int* c = malloc((size_t)n * sizeof(int));
        if (pw && c) {
            int k = galg_max_clique(pw, c);
            if (k > 0) {
                Expr* one = galg_vertex_list(g, c, k, 0);
                out = expr_new_function(expr_new_symbol(SYM_List), &one, 1);
            }
        }
        free(c); galg_ug_free(pw);
    }
    free(dist); free(q); free(pu); free(pv); galg_ug_free(u);
    return out;
}
