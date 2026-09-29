/* spanningtree.c - FindSpanningTree[g], FindSpanningTree[{g, v}].
 *
 * Mathematica's semantics, checked case by case against Mathematica 15:
 *
 *   undirected, unweighted   a BFS spanning forest;
 *   undirected, EdgeWeight   a MINIMUM spanning forest (Kruskal): ties between
 *                            equal weights are broken by the edge's vertex
 *                            positions, lowest first, which reproduces
 *                            Mathematica's choice on tied inputs;
 *   directed,   unweighted   a BFS branching (edges followed forwards), its
 *                            roots taken in decreasing DFS finishing time, so
 *                            the forest needs as few roots as possible;
 *   directed,   EdgeWeight   a MINIMUM-WEIGHT spanning branching (Chu-Liu /
 *                            Edmonds) among those with the fewest roots -- a
 *                            minimum spanning arborescence whenever one vertex
 *                            reaches all;
 *   mixed                    left unevaluated, as Mathematica does ("not
 *                            implemented").
 *
 * {g, v} restricts the answer to the tree grown from v: v's component
 * (undirected) or the vertices v reaches (directed), rooted at v; the result's
 * VertexList is just those vertices. A v that is not a vertex of g emits
 * FindSpanningTree::inv and stays unevaluated. Trailing options (Method -> ...)
 * are accepted and ignored: every method yields the same optimum.
 *
 * Weights are compared EXACTLY. Integer, Rational, Real and arbitrary-precision
 * weights become GMP rationals (a Real's exact binary value), so 1/3 and
 * 0.3333333333333333 are told apart and nothing is rounded; other numeric
 * expressions (Sqrt[2], Pi) are compared through N[w, 40]. A weight that is not
 * a real number (a symbol, a Complex) leaves the call unevaluated, as in
 * Mathematica.
 *
 * The result is Graph[verts, treeEdges(, EdgeWeight -> w)(, EdgeCapacity -> c)]
 * with each tree edge's own weight and capacity carried over. As in Mathematica,
 * an undirected tree edge is written lower VertexList position first and the
 * edges are sorted by the positions of their endpoints.
 *
 * Everything runs on the validated-graph memo's integer endpoint arrays: no
 * expression is hashed or compared after the weights are read.
 *
 * Memory (SPEC section 4): returns a freshly allocated Graph built from copies
 * of g's parts; res is never modified (the evaluator frees it on success and
 * keeps it on NULL).
 */

#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ---- Exact weights --------------------------------------------------------- */

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static int is_int_like(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

/* dst = the integer e, into an ALREADY initialized mpz (expr_to_mpz itself
 * initializes its target, so it cannot write into an mpq's parts). */
static void set_mpz(mpz_t dst, const Expr* e) {
    mpz_t tmp;
    expr_to_mpz(e, tmp);
    mpz_set(dst, tmp);
    mpz_clear(tmp);
}

/* The value of an already-numeric atom (or Rational) as an exact rational.
 * Returns 1 on success, 0 if w is not one of those shapes. */
static int atom_to_mpq(const Expr* w, mpq_t q) {
    if (is_int_like(w)) {
        set_mpz(mpq_numref(q), w);
        mpz_set_ui(mpq_denref(q), 1);
        return 1;
    }
    if (w->type == EXPR_REAL) {
        if (!isfinite(w->data.real)) return 0;
        mpq_set_d(q, w->data.real);                 /* exact binary value */
        return 1;
    }
#ifdef USE_MPFR
    if (w->type == EXPR_MPFR) {
        if (!mpfr_number_p(w->data.mpfr)) return 0;
        mpfr_get_q(q, w->data.mpfr);                /* exact binary value */
        return 1;
    }
#endif
    if (head_is(w, SYM_Rational) && w->data.function.arg_count == 2
        && is_int_like(w->data.function.args[0]) && is_int_like(w->data.function.args[1])) {
        set_mpz(mpq_numref(q), w->data.function.args[0]);
        set_mpz(mpq_denref(q), w->data.function.args[1]);
        if (mpz_sgn(mpq_denref(q)) == 0) return 0;
        mpq_canonicalize(q);
        return 1;
    }
    return 0;
}

/* Exact rational value of an edge weight; a non-atomic numeric expression is
 * numericized first. Returns 0 if w is not a real number. */
static int weight_to_mpq(const Expr* w, mpq_t q) {
    if (atom_to_mpq(w, q)) return 1;
    /* Symbols too: Pi and E are numeric; a free x numericizes to itself and
     * is rejected below. */
    if ((w->type != EXPR_FUNCTION && w->type != EXPR_SYMBOL) || head_is(w, SYM_Complex))
        return 0;
    Expr* nargs[2] = { expr_copy((Expr*)w), expr_new_integer(40) };
    Expr* nv = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), nargs, 2));
    int ok = nv && atom_to_mpq(nv, q);
    if (nv) expr_free(nv);
    return ok;
}

/* ---- Small utilities ------------------------------------------------------- */

/* Stable merge sort of ids[0..n-1] under cmp(ctx, a, b) (C99 qsort has no
 * context argument, and a file-static context would not be re-entrant). */
typedef int (*IdCmp)(const void* ctx, int a, int b);

static void merge_sort_ids(int* ids, int* tmp, size_t n, IdCmp cmp, const void* ctx) {
    for (size_t width = 1; width < n; width *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * width) {
            size_t mid = lo + width < n ? lo + width : n;
            size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
            size_t i = lo, j = mid, k = lo;
            while (i < mid && j < hi) tmp[k++] = (cmp(ctx, ids[j], ids[i]) < 0) ? ids[j++] : ids[i++];
            while (i < mid) tmp[k++] = ids[i++];
            while (j < hi)  tmp[k++] = ids[j++];
        }
        memcpy(ids, tmp, n * sizeof(int));
    }
}

/* Out-incidence CSR with edge ids: every edge k contributes u -> v, and an
 * undirected one also v -> u, in EdgeList order (the order GraphAdj uses). */
typedef struct { int* start; int* nbr; int* eid; } Inc;

static void inc_free(Inc* c) { free(c->start); free(c->nbr); free(c->eid); }

static int inc_build(Inc* c, int n, size_t ne, const int* eu, const int* ev,
                     const unsigned char* edir) {
    size_t slots = 0;
    for (size_t k = 0; k < ne; k++) slots += edir[k] ? 1 : 2;
    c->start = calloc((size_t)n + 1, sizeof(int));
    c->nbr = malloc((slots ? slots : 1) * sizeof(int));
    c->eid = malloc((slots ? slots : 1) * sizeof(int));
    int* fill = calloc((size_t)n + 1, sizeof(int));
    if (!c->start || !c->nbr || !c->eid || !fill) { free(fill); inc_free(c); return 0; }
    for (size_t k = 0; k < ne; k++) {
        c->start[eu[k] + 1]++;
        if (!edir[k]) c->start[ev[k] + 1]++;
    }
    for (int i = 0; i < n; i++) c->start[i + 1] += c->start[i];
    for (size_t k = 0; k < ne; k++) {
        int p = c->start[eu[k]] + fill[eu[k]]++;
        c->nbr[p] = ev[k]; c->eid[p] = (int)k;
        if (!edir[k]) {
            p = c->start[ev[k]] + fill[ev[k]]++;
            c->nbr[p] = eu[k]; c->eid[p] = (int)k;
        }
    }
    free(fill);
    return 1;
}

/* Mark in `mask` every vertex reachable from r along the incidence. */
static int reach_from(const Inc* c, int n, int r, char* mask) {
    int* q = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!q) return 0;
    int head = 0, tail = 0;
    mask[r] = 1; q[tail++] = r;
    while (head < tail) {
        int u = q[head++];
        for (int p = c->start[u]; p < c->start[u + 1]; p++)
            if (!mask[c->nbr[p]]) { mask[c->nbr[p]] = 1; q[tail++] = c->nbr[p]; }
    }
    free(q);
    return 1;
}

/* ---- Unweighted: BFS forest / branching ------------------------------------ */

/* Vertices in decreasing DFS finishing time (iterative DFS, starts in
 * VertexList order, neighbours in incidence order). The first vertex of that
 * order lies in a source component, so BFS started from the vertices in this
 * order needs the fewest roots. */
static int* dfs_finish_desc(const Inc* c, int n) {
    size_t nn = (size_t)(n > 0 ? n : 1);
    int* order = malloc(nn * sizeof(int));
    int* stk = malloc(nn * sizeof(int));
    int* cur = malloc(nn * sizeof(int));
    char* seen = calloc(nn, 1);
    if (!order || !stk || !cur || !seen) { free(order); free(stk); free(cur); free(seen); return NULL; }
    int done = n;
    for (int s = 0; s < n; s++) {
        if (seen[s]) continue;
        int sp = 0; stk[sp++] = s; seen[s] = 1; cur[s] = c->start[s];
        while (sp > 0) {
            int u = stk[sp - 1];
            if (cur[u] < c->start[u + 1]) {
                int w = c->nbr[cur[u]++];
                if (!seen[w]) { seen[w] = 1; cur[w] = c->start[w]; stk[sp++] = w; }
            } else {
                order[--done] = u;       /* filled back to front: decreasing finish */
                sp--;
            }
        }
    }
    free(stk); free(cur); free(seen);
    return order;
}

/* BFS from each start (in order) not yet reached, recording the edge that
 * first reaches each vertex. Appends tree edge ids to *pick. */
static int bfs_forest(const Inc* c, int n, const int* starts, int nstarts,
                      int* pick, int* npick) {
    size_t nn = (size_t)(n > 0 ? n : 1);
    char* seen = calloc(nn, 1);
    int* q = malloc(nn * sizeof(int));
    if (!seen || !q) { free(seen); free(q); return 0; }
    for (int si = 0; si < nstarts; si++) {
        int s = starts[si];
        if (seen[s]) continue;
        int head = 0, tail = 0;
        seen[s] = 1; q[tail++] = s;
        while (head < tail) {
            int u = q[head++];
            for (int p = c->start[u]; p < c->start[u + 1]; p++) {
                int w = c->nbr[p];
                if (seen[w]) continue;
                seen[w] = 1; q[tail++] = w;
                pick[(*npick)++] = c->eid[p];
            }
        }
    }
    free(seen); free(q);
    return 1;
}

/* ---- Undirected, weighted: Kruskal ----------------------------------------- */

typedef struct { const mpq_t* w; const int* eu; const int* ev; } KrCtx;

static int kr_cmp(const void* vctx, int a, int b) {
    const KrCtx* k = vctx;
    int c = mpq_cmp(k->w[a], k->w[b]);
    if (c) return c;
    int alo = k->eu[a] < k->ev[a] ? k->eu[a] : k->ev[a];
    int ahi = k->eu[a] < k->ev[a] ? k->ev[a] : k->eu[a];
    int blo = k->eu[b] < k->ev[b] ? k->eu[b] : k->ev[b];
    int bhi = k->eu[b] < k->ev[b] ? k->ev[b] : k->eu[b];
    if (alo != blo) return alo < blo ? -1 : 1;
    if (ahi != bhi) return ahi < bhi ? -1 : 1;
    return 0;
}

static int uf_find(int* p, int x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

/* Minimum spanning forest over the edges whose endpoints are in `mask` (all
 * when NULL). */
static int kruskal(int n, size_t ne, const int* eu, const int* ev, const mpq_t* w,
                   const char* mask, int* pick, int* npick) {
    size_t nb = ne ? ne : 1;
    int* ids = malloc(nb * sizeof(int));
    int* tmp = malloc(nb * sizeof(int));
    int* parent = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* size = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!ids || !tmp || !parent || !size) { free(ids); free(tmp); free(parent); free(size); return 0; }
    size_t m = 0;
    for (size_t k = 0; k < ne; k++)
        if (!mask || mask[eu[k]]) ids[m++] = (int)k;   /* mask is a whole component */
    KrCtx ctx = { w, eu, ev };
    merge_sort_ids(ids, tmp, m, kr_cmp, &ctx);
    for (int i = 0; i < n; i++) { parent[i] = i; size[i] = 1; }
    for (size_t i = 0; i < m; i++) {
        int k = ids[i];
        int a = uf_find(parent, eu[k]), b = uf_find(parent, ev[k]);
        if (a == b) continue;
        if (size[a] < size[b]) { int t = a; a = b; b = t; }
        parent[b] = a; size[a] += size[b];
        pick[(*npick)++] = k;
    }
    free(ids); free(tmp); free(parent); free(size);
    return 1;
}

/* ---- Directed, weighted: Chu-Liu / Edmonds --------------------------------- *
 * The O(E log V) formulation (Tarjan; Gabow et al.) with a rollback union-find
 * for the contracted cycles and mergeable heaps of each super-vertex's incoming
 * edges carrying a lazy additive delta -- after KACTL's DirectedMST -- here with
 * exact GMP rational keys and LEFTIST heaps, whose merge recursion is bounded
 * by O(log E) where a skew heap's is not. Ties between equal keys prefer the
 * higher tail position, which is the choice Mathematica makes. */

typedef struct HNode {
    mpq_t key, delta;
    int eid, rank;
    struct HNode *l, *r;
} HNode;

typedef struct {
    const int* src;      /* tail of edge id (super-root edges: the super root) */
    const int* dst;      /* head of edge id */
} DmCtx;

static void hprop(HNode* a) {
    if (mpq_sgn(a->delta) == 0) return;
    mpq_add(a->key, a->key, a->delta);
    if (a->l) mpq_add(a->l->delta, a->l->delta, a->delta);
    if (a->r) mpq_add(a->r->delta, a->r->delta, a->delta);
    mpq_set_ui(a->delta, 0, 1);
}

/* a before b: smaller key, then the higher tail. Both already propagated. */
static int hbefore(const DmCtx* c, const HNode* a, const HNode* b) {
    int s = mpq_cmp(a->key, b->key);
    if (s) return s < 0;
    return c->src[a->eid] > c->src[b->eid];
}

static HNode* hmerge(const DmCtx* c, HNode* a, HNode* b) {
    if (!a) return b;
    if (!b) return a;
    hprop(a); hprop(b);
    if (hbefore(c, b, a)) { HNode* t = a; a = b; b = t; }
    a->r = hmerge(c, a->r, b);
    if (!a->l || a->l->rank < a->r->rank) { HNode* t = a->l; a->l = a->r; a->r = t; }
    a->rank = (a->r ? a->r->rank : 0) + 1;
    return a;
}

static void hpop(const DmCtx* c, HNode** a) {
    hprop(*a);
    *a = hmerge(c, (*a)->l, (*a)->r);
}

/* Rollback union-find: union by size, no path compression, history stack. */
typedef struct { int* e; int* st_a; int* st_v; int sp; } RUF;

static int ruf_find(const RUF* u, int x) { while (u->e[x] >= 0) x = u->e[x]; return x; }

static int ruf_join(RUF* u, int a, int b) {
    a = ruf_find(u, a); b = ruf_find(u, b);
    if (a == b) return 0;
    if (u->e[a] > u->e[b]) { int t = a; a = b; b = t; }
    u->st_a[u->sp] = a; u->st_v[u->sp++] = u->e[a];
    u->st_a[u->sp] = b; u->st_v[u->sp++] = u->e[b];
    u->e[a] += u->e[b]; u->e[b] = a;
    return 1;
}

static void ruf_rollback(RUF* u, int t) {
    while (u->sp > t) { u->sp--; u->e[u->st_a[u->sp]] = u->st_v[u->sp]; }
}

/* Minimum branching over N vertices rooted at `root` whose edges are ids
 * 0..m-1 (src/dst/w). skip[v] != 0 marks a vertex outside the problem. Fills
 * in_edge[v] with the chosen incoming edge id (-1 for the root and skipped
 * vertices). Returns 1, 0 on allocation failure, -1 if some vertex is
 * unreachable (cannot happen for the callers below). */
static int dmst(int N, int root, int m, const int* src, const int* dst, const mpq_t* w,
                const char* skip, int* in_edge) {
    size_t nn = (size_t)N, mm = (size_t)(m > 0 ? m : 1);
    DmCtx ctx = { src, dst };
    HNode* pool = malloc(mm * sizeof(HNode));
    HNode** heap = calloc(nn, sizeof(HNode*));
    int* seen = malloc(nn * sizeof(int));
    int* path = malloc(nn * sizeof(int));
    int* Q = malloc(nn * sizeof(int));
    RUF uf;
    uf.e = malloc(nn * sizeof(int));
    uf.st_a = malloc(2 * nn * sizeof(int));
    uf.st_v = malloc(2 * nn * sizeof(int));
    uf.sp = 0;
    /* cycle records: super-vertex, rollback time, and its edges in cbuf */
    int* cy_u = malloc(nn * sizeof(int));
    int* cy_t = malloc(nn * sizeof(int));
    int* cy_lo = malloc(nn * sizeof(int));
    int* cy_hi = malloc(nn * sizeof(int));
    size_t ccap = mm + nn, clen = 0;
    int* cbuf = malloc(ccap * sizeof(int));
    int ncy = 0, status = 0;
    int built = 0;
    if (!pool || !heap || !seen || !path || !Q || !uf.e || !uf.st_a || !uf.st_v
        || !cy_u || !cy_t || !cy_lo || !cy_hi || !cbuf)
        goto done;

    for (int k = 0; k < m; k++) {
        HNode* h = &pool[k];
        mpq_init(h->key); mpq_init(h->delta);
        mpq_set(h->key, w[k]);
        h->eid = k; h->rank = 1; h->l = h->r = NULL;
        built++;
        if (!skip[dst[k]] && !skip[src[k]] && dst[k] != root)
            heap[dst[k]] = hmerge(&ctx, heap[dst[k]], h);
    }
    for (int v = 0; v < N; v++) { uf.e[v] = -1; seen[v] = skip[v] ? v : -1; in_edge[v] = -1; }
    seen[root] = root;

    for (int s = 0; s < N; s++) {
        int u = s, qi = 0;
        while (seen[u] < 0) {
            if (!heap[u]) { status = -1; goto done; }
            hprop(heap[u]);
            HNode* top = heap[u];
            int e = top->eid;
            /* every remaining incoming edge of u is now priced relative to e */
            mpq_sub(top->delta, top->delta, top->key);
            hpop(&ctx, &heap[u]);
            Q[qi] = e; path[qi++] = u; seen[u] = s;
            u = ruf_find(&uf, src[e]);
            if (seen[u] == s) {                /* found a cycle: contract it */
                HNode* cyc = NULL;
                int end = qi, t = uf.sp, x;
                do { x = path[--qi]; cyc = hmerge(&ctx, cyc, heap[x]); }
                while (ruf_join(&uf, u, x));
                u = ruf_find(&uf, u);
                heap[u] = cyc; seen[u] = -1;
                cy_u[ncy] = u; cy_t[ncy] = t; cy_lo[ncy] = (int)clen;
                for (int i = qi; i < end; i++) cbuf[clen++] = Q[i];
                cy_hi[ncy++] = (int)clen;
            }
        }
        for (int i = 0; i < qi; i++) in_edge[ruf_find(&uf, dst[Q[i]])] = Q[i];
    }

    /* Expand the cycles, newest first. */
    for (int c = ncy - 1; c >= 0; c--) {
        ruf_rollback(&uf, cy_t[c]);
        int in_e = in_edge[cy_u[c]];
        for (int i = cy_lo[c]; i < cy_hi[c]; i++) in_edge[ruf_find(&uf, dst[cbuf[i]])] = cbuf[i];
        if (in_e >= 0) in_edge[ruf_find(&uf, dst[in_e])] = in_e;
    }
    in_edge[root] = -1;
    for (int v = 0; v < N; v++) if (skip[v]) in_edge[v] = -1;
    status = 1;

done:
    if (pool) for (int k = 0; k < built; k++) { mpq_clear(pool[k].key); mpq_clear(pool[k].delta); }
    free(pool); free(heap); free(seen); free(path); free(Q);
    free(uf.e); free(uf.st_a); free(uf.st_v);
    free(cy_u); free(cy_t); free(cy_lo); free(cy_hi); free(cbuf);
    return status;
}

/* Minimum branching of the directed graph (eu/ev/w) with the fewest roots, or
 * the minimum arborescence rooted at `root` over the vertices it reaches (mask)
 * when root >= 0. Appends the chosen edge ids. */
static int edmonds(int n, size_t ne, const int* eu, const int* ev, const mpq_t* w,
                   int root, const char* mask, int* pick, int* npick) {
    int rooted = root >= 0;
    int N = rooted ? n : n + 1;                   /* + a super root */
    int m = rooted ? (int)ne : (int)ne + n;
    size_t mm = (size_t)(m > 0 ? m : 1);
    int* src = malloc(mm * sizeof(int));
    int* dst = malloc(mm * sizeof(int));
    mpq_t* ww = malloc(mm * sizeof(mpq_t));
    char* skip = calloc((size_t)N, 1);
    int* in_edge = malloc((size_t)N * sizeof(int));
    int ok = 0, inited = 0;
    if (!src || !dst || !ww || !skip || !in_edge) goto out;
    for (int k = 0; k < m; k++) { mpq_init(ww[k]); inited++; }
    for (size_t k = 0; k < ne; k++) { src[k] = eu[k]; dst[k] = ev[k]; mpq_set(ww[k], w[k]); }
    if (rooted) {
        for (int v = 0; v < n; v++) skip[v] = !mask[v];
    } else {
        /* Super-root edges cost M = 2 * sum|w| + 1, more than any two
         * branchings' real weights can differ by, so the optimum first has the
         * fewest roots and then the least weight. */
        mpq_t M, a;
        mpq_init(M); mpq_init(a);
        for (size_t k = 0; k < ne; k++) { mpq_abs(a, w[k]); mpq_add(M, M, a); }
        mpq_add(M, M, M);
        mpq_set_ui(a, 1, 1); mpq_add(M, M, a);
        for (int v = 0; v < n; v++) {
            src[ne + (size_t)v] = n; dst[ne + (size_t)v] = v;
            mpq_set(ww[ne + (size_t)v], M);
        }
        mpq_clear(M); mpq_clear(a);
        root = n;
    }
    int st = dmst(N, root, m, src, dst, (const mpq_t*)ww, skip, in_edge);
    if (st == 1) {
        for (int v = 0; v < n; v++)
            if (in_edge[v] >= 0 && (size_t)in_edge[v] < ne) pick[(*npick)++] = in_edge[v];
        ok = 1;
    }
out:
    if (ww) for (int k = 0; k < inited; k++) mpq_clear(ww[k]);
    free(src); free(dst); free(ww); free(skip); free(in_edge);
    return ok;
}

/* ---- Result assembly ------------------------------------------------------- */

typedef struct { const int* a; const int* b; } PairCtx;

static int pair_cmp(const void* vctx, int x, int y) {
    const PairCtx* c = vctx;
    if (c->a[x] != c->a[y]) return c->a[x] < c->a[y] ? -1 : 1;
    if (c->b[x] != c->b[y]) return c->b[x] < c->b[y] ? -1 : 1;
    return 0;
}

/* Graph[verts in mask, tree edges (normalized, sorted), per-edge options]. */
static Expr* build_tree(const Expr* g, const int* eu, const int* ev, const unsigned char* edir,
                        const char* mask, int* pick, int npick) {
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t np = (size_t)(npick > 0 ? npick : 1);

    /* Endpoint positions as written in the result: an undirected edge lower
     * position first. Sort by them. */
    int* ta = malloc(np * sizeof(int));
    int* tb = malloc(np * sizeof(int));
    int* tmp = malloc(np * sizeof(int));
    int* ord = malloc(np * sizeof(int));
    int* newpos = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!ta || !tb || !tmp || !ord || !newpos) {
        free(ta); free(tb); free(tmp); free(ord); free(newpos);
        return NULL;
    }
    for (int i = 0; i < npick; i++) {
        int k = pick[i], a = eu[k], b = ev[k];
        if (!edir[k] && a > b) { int t = a; a = b; b = t; }
        ta[i] = a; tb[i] = b; ord[i] = i;
    }
    PairCtx pc = { ta, tb };
    merge_sort_ids(ord, tmp, (size_t)npick, pair_cmp, &pc);

    /* Vertex list and its index (keys borrow the copied vertices, which the
     * result keeps alive). */
    int nv = 0;
    for (int i = 0; i < n; i++) newpos[i] = (!mask || mask[i]) ? nv++ : -1;
    Expr** vs = malloc((size_t)(nv > 0 ? nv : 1) * sizeof(Expr*));
    Expr** es = malloc(np * sizeof(Expr*));
    int* neu = malloc(np * sizeof(int));
    int* nev = malloc(np * sizeof(int));
    unsigned char* ned = malloc(np);
    GraphVIdx* ix = graph_vidx_new((size_t)nv);
    const Expr* wl = graph_edge_weight_list(g);
    const Expr* cl = graph_edge_capacity_list(g);
    Expr** ws = wl ? malloc(np * sizeof(Expr*)) : NULL;
    Expr** cs = cl ? malloc(np * sizeof(Expr*)) : NULL;
    if (!vs || !es || !neu || !nev || !ned || !ix || (wl && !ws) || (cl && !cs)) {
        free(vs); free(es); free(neu); free(nev); free(ned); graph_vidx_free(ix);
        free(ws); free(cs);
        free(ta); free(tb); free(tmp); free(ord); free(newpos);
        return NULL;
    }
    for (int i = 0, j = 0; i < n; i++) {
        if (newpos[i] < 0) continue;
        vs[j] = expr_copy(verts->data.function.args[i]);
        graph_vidx_put(ix, vs[j], j);
        j++;
    }
    Expr* uhead = NULL;
    for (int j = 0; j < npick; j++) {
        int i = ord[j], k = pick[i];
        Expr* e = edges->data.function.args[k];
        if (!edir[k] && ta[i] != eu[k]) {
            /* written the other way round: re-emit lower position first */
            if (!uhead) uhead = expr_new_symbol(SYM_UndirectedEdge);
            Expr* ab[2] = { expr_copy(e->data.function.args[1]), expr_copy(e->data.function.args[0]) };
            es[j] = expr_new_function(expr_copy(uhead), ab, 2);
        } else {
            es[j] = expr_copy(e);
        }
        neu[j] = newpos[ta[i]]; nev[j] = newpos[tb[i]]; ned[j] = edir[k];
        if (ws) ws[j] = expr_copy(wl->data.function.args[k]);
        if (cs) cs[j] = expr_copy(cl->data.function.args[k]);
    }
    if (uhead) expr_free(uhead);

    Expr* gargs[4];
    size_t gargc = 0;
    gargs[gargc++] = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)nv);
    gargs[gargc++] = expr_new_function(expr_new_symbol(SYM_List), es, (size_t)npick);
    Expr* lists[2] = { ws ? expr_new_function(expr_new_symbol(SYM_List), ws, (size_t)npick) : NULL,
                       cs ? expr_new_function(expr_new_symbol(SYM_List), cs, (size_t)npick) : NULL };
    for (int r = 0; r < 2; r++) {
        if (!lists[r]) continue;
        Expr* rargs[2] = { expr_new_symbol(graph_edge_option_key(r)), lists[r] };
        gargs[gargc++] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_Graph), gargs, gargc);
    free(vs); free(es); free(ws); free(cs);
    free(ta); free(tb); free(tmp); free(ord); free(newpos);
    /* A subforest of a valid graph is valid; seeding the memo just saves the
     * next consumer a validation pass. It owns ix/neu/nev/ned either way. */
    (void)graph_memo_seed(out, ix, neu, nev, ned);
    return out;
}

/* ---- The builtin ----------------------------------------------------------- */

Expr* builtin_find_spanning_tree(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    for (size_t i = 1; i < argc; i++) {                /* options: accepted, unused */
        const Expr* o = res->data.function.args[i];
        if (!(head_is(o, SYM_Rule) || head_is(o, SYM_RuleDelayed))
            || o->data.function.arg_count != 2)
            return NULL;
    }

    const Expr* g = res->data.function.args[0];
    const Expr* rootv = NULL;
    if (graph_is_list(g)) {                            /* {g, v} */
        if (g->data.function.arg_count != 2) return NULL;
        rootv = g->data.function.args[1];
        g = g->data.function.args[0];
    }
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return NULL;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    long ndir = graph_directed_edge_count(g);
    if (ndir > 0 && (size_t)ndir < ne) return NULL;   /* mixed: not implemented */
    int directed = ndir > 0;

    int root = -1;
    if (rootv) {
        root = graph_vertex_position(g, rootv);
        if (root < 0) {
            char* s = expr_to_string((Expr*)rootv);
            char* call = expr_to_string(res);
            mth_message("FindSpanningTree", "inv",
                        "The argument %s in %s is not a valid vertex.",
                        s ? s : "?", call ? call : "FindSpanningTree");
            free(s); free(call);
            return NULL;
        }
    }

    /* Exact weights, before anything else can evict g from the memo. */
    const Expr* wl = graph_edge_weight_list(g);
    mpq_t* w = NULL;
    size_t winit = 0;
    int ok = 1;
    if (wl) {
        w = malloc((ne ? ne : 1) * sizeof(mpq_t));
        if (!w) return NULL;
        for (; winit < ne && ok; winit++) {
            mpq_init(w[winit]);
            ok = weight_to_mpq(wl->data.function.args[winit], w[winit]);
        }
    }
    /* weight_to_mpq may evaluate N[], which can touch other graphs; re-read
     * the (still valid, g being alive) memo arrays afterwards. */
    if (ok && wl && !graph_edge_indices(g, &eu, &ev, &edir)) ok = 0;

    size_t nn = (size_t)(n > 0 ? n : 1);
    int* pick = ok ? malloc(nn * sizeof(int)) : NULL;   /* a forest: < n edges */
    char* mask = NULL;
    int npick = 0;
    Inc inc;
    memset(&inc, 0, sizeof(inc));
    Expr* out = NULL;
    if (!pick || !inc_build(&inc, n, ne, eu, ev, edir)) goto done;
    if (root >= 0) {
        mask = calloc(nn, 1);
        if (!mask || !reach_from(&inc, n, root, mask)) goto done;
    }

    if (!wl) {
        int* starts;
        int nstarts;
        if (root >= 0) { starts = malloc(sizeof(int)); if (!starts) goto done; starts[0] = root; nstarts = 1; }
        else { starts = dfs_finish_desc(&inc, n); if (!starts) goto done; nstarts = n; }
        int bok = bfs_forest(&inc, n, starts, nstarts, pick, &npick);
        free(starts);
        if (!bok) goto done;
    } else if (!directed) {
        if (!kruskal(n, ne, eu, ev, (const mpq_t*)w, mask, pick, &npick)) goto done;
    } else {
        if (!edmonds(n, ne, eu, ev, (const mpq_t*)w, root, mask, pick, &npick)) goto done;
    }
    out = build_tree(g, eu, ev, edir, mask, pick, npick);

done:
    if (w) { for (size_t k = 0; k < winit; k++) mpq_clear(w[k]); free(w); }
    free(pick); free(mask); inc_free(&inc);
    return out;
}
