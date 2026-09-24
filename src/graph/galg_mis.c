/* galg_mis.c - independent sets and vertex covers.
 *
 *   FindIndependentVertexSet[g]              {s}: s a MAXIMUM independent set
 *   FindIndependentVertexSet[g, spec(, n)]   maximal independent sets by size
 *                                            spec (k, {k}, {kmin, kmax}; n or
 *                                            All of them), as for FindClique
 *   FindVertexCover[g]                       a MINIMUM vertex cover
 *   IndependentVertexSetQ[g, vs]             no two of vs adjacent
 *   VertexCoverQ[g, vs]                      every edge has an end in vs
 *   IndependentEdgeSetQ[g, es]               es edges of g, pairwise disjoint
 *   EdgeCoverQ[g, es]                        es edges of g touching every vertex
 *
 * Edge direction is ignored for independence and covering (Wolfram's
 * semantics). The *Q predicates give False for a non-graph, for list elements
 * that are not vertices / edges of g (an UndirectedEdge matches in either
 * orientation, a DirectedEdge only as given), and never stay unevaluated.
 * Repeated vertices are allowed in IndependentVertexSetQ / VertexCoverQ.
 *
 * Exactness. A minimum vertex cover is the complement of a maximum independent
 * set, so both heads share one exact solver (gm_mis). It returns an optimum or
 * gives up -- after MIS_MAX_NODES search nodes, or when TimeConstrained's
 * deadline passes (galg_poll) -- and then the head stays unevaluated. It never
 * returns a merely-maximal set.
 *
 * The solver:
 *   - Connected components are solved independently, and a component that is
 *     not sparse (average degree >= 8 or density >= 0.05, at most 3000
 *     vertices; the measured crossover) goes to the bitset
 *     maximum-clique branch and bound on its complement (galg_clique.c).
 *   - Sparse components use branch and reduce, seeded with a greedy
 *     minimum-degree solution as the incumbent. At every node the reductions
 *       degree 0 / degree 1   take the vertex (some optimum contains it),
 *       degree 2, triangle    take it,
 *       degree 2, folding     v with non-adjacent neighbours a, b: replace
 *                             {v, a, b} by one vertex w adjacent to
 *                             N(a) u N(b); alpha drops by exactly 1, and the
 *                             optimum lifts back ({a, b} if w was chosen,
 *                             else v),
 *       domination            if N[v] is inside N[u] for adjacent u, v, drop u,
 *     run to a fixed point (all undone through one trail); then the bound
 *     "solution so far + greedy clique cover of what is left" prunes (an
 *     independent set meets each clique at most once); a disconnected
 *     remainder splits into independent subproblems; otherwise the search
 *     branches on a maximum-degree vertex v: drop v together with its mirrors
 *     (Fomin-Grandoni-Kratsch), or take v.
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

#define MIS_MAX_NODES 20000000L
#define MIS_DENSE_MAX 3000
#define MIS_DOM_DEG 8           /* domination is tested from vertices this small */

/* A degree-2 fold: v's neighbours a, b (non-adjacent) and v collapse into the
 * new vertex w, adjacent to N(a) u N(b) - v; alpha(G) = alpha(G') + 1, and an
 * optimum of G' lifts back by taking {a, b} if w was chosen, else v. */
typedef struct { int v, a, b, w; } GmFold;

typedef struct {
    int cap;                     /* vertex-id capacity (n + folds)            */
    int nid;                     /* next free id (ids are stack-allocated)    */
    int** nb; int* nbn; int* nbcap;   /* growable adjacency lists             */
    char* alive;
    int* deg;                    /* live degree                               */
    int* tr; int tlen;           /* trail: v >= 0 removed v; -(w+1) fold w    */
    int* chosen; int clen;       /* vertices taken on the current path        */
    GmFold* folds; int flen;     /* folds on the current path                 */
    int* out; int olen; int ocap;/* solution stack (see gm_rec), growable     */
    int* wq; char* inq; int wqn; /* reduction worklist                        */
    int* mark; int stamp;        /* scratch stamps                            */
    int* mark2; int stamp2;
    int* cid; int* csz; int* ccnt; int* ctouch;   /* clique-cover scratch     */
    char* insol;
    long nodes;
    int aborted;
} GmMis;

/* Room for `extra` more entries on the solution stack. Nested search levels
 * each keep their incumbent there, so it can outgrow n; returns 0 (and marks
 * the search aborted) on allocation failure. M->out may move. */
static int gm_reserve(GmMis* M, int extra) {
    if (M->olen + extra <= M->ocap) return 1;
    long nc = M->ocap ? 2L * M->ocap : 1024;
    while (nc < (long)M->olen + extra) nc *= 2;
    int* t = realloc(M->out, (size_t)nc * sizeof(int));
    if (!t) { M->aborted = 1; return 0; }
    M->out = t; M->ocap = (int)nc;
    return 1;
}

static int gm_nb_push(GmMis* M, int x, int w) {
    if (M->nbn[x] == M->nbcap[x]) {
        int nc = M->nbcap[x] ? 2 * M->nbcap[x] : 4;
        int* t = realloc(M->nb[x], (size_t)nc * sizeof(int));
        if (!t) { M->aborted = 1; return 0; }
        M->nb[x] = t; M->nbcap[x] = nc;
    }
    M->nb[x][M->nbn[x]++] = w;
    return 1;
}

static void gm_push_wq(GmMis* M, int v) {
    if (!M->inq[v]) { M->inq[v] = 1; M->wq[M->wqn++] = v; }
}

static void gm_remove(GmMis* M, int v) {
    M->alive[v] = 0;
    M->tr[M->tlen++] = v;
    const int* a = M->nb[v];
    for (int j = 0, e = M->nbn[v]; j < e; j++) {
        int x = a[j];
        if (M->alive[x]) { M->deg[x]--; gm_push_wq(M, x); }
    }
}

static void gm_undo(GmMis* M, int mark) {
    while (M->tlen > mark) {
        int e = M->tr[--M->tlen];
        if (e >= 0) {
            const int* a = M->nb[e];
            for (int j = 0, k = M->nbn[e]; j < k; j++) if (M->alive[a[j]]) M->deg[a[j]]++;
            M->alive[e] = 1;
        } else {
            int w = -e - 1;
            const int* a = M->nb[w];
            for (int j = 0, k = M->nbn[w]; j < k; j++) { M->nbn[a[j]]--; M->deg[a[j]]--; }
            M->alive[w] = 0;
            M->nid--;
        }
    }
}

/* Take v: it joins the solution; v and its live neighbours leave. */
static void gm_take(GmMis* M, int v) {
    M->chosen[M->clen++] = v;
    const int* a = M->nb[v];
    for (int j = 0, e = M->nbn[v]; j < e; j++) if (M->alive[a[j]]) gm_remove(M, a[j]);
    gm_remove(M, v);
}

static int gm_adjacent(GmMis* M, int a, int b) {
    if (M->nbn[a] > M->nbn[b]) { int t = a; a = b; b = t; }
    const int* l = M->nb[a];
    for (int j = 0, e = M->nbn[a]; j < e; j++) if (l[j] == b) return M->alive[b];
    return 0;
}

static void gm_fold(GmMis* M, int v, int a, int b) {
    if (M->nid >= M->cap) { M->aborted = 1; return; }
    int w = M->nid++;
    M->nbn[w] = 0;
    M->stamp++;
    M->mark[v] = M->mark[a] = M->mark[b] = M->stamp;
    for (int s = 0; s < 2; s++) {
        int x = s ? b : a;
        const int* l = M->nb[x];
        for (int j = 0, e = M->nbn[x]; j < e; j++) {
            int y = l[j];
            if (!M->alive[y] || M->mark[y] == M->stamp) continue;
            M->mark[y] = M->stamp;
            if (!gm_nb_push(M, w, y)) return;
        }
    }
    gm_remove(M, v); gm_remove(M, a); gm_remove(M, b);
    M->alive[w] = 1;
    M->deg[w] = M->nbn[w];
    for (int j = 0, e = M->nbn[w]; j < e; j++) {
        int x = M->nb[w][j];
        if (!gm_nb_push(M, x, w)) return;
        M->deg[x]++;
        gm_push_wq(M, x);
    }
    gm_push_wq(M, w);
    M->tr[M->tlen++] = -(w + 1);
    GmFold f; f.v = v; f.a = a; f.b = b; f.w = w;
    M->folds[M->flen++] = f;
}

/* Is N[v] contained in N[x] (both alive, adjacent)? */
static int gm_dominates(GmMis* M, int x, int v) {
    M->stamp++;
    M->mark[x] = M->stamp;
    const int* l = M->nb[x];
    for (int j = 0, e = M->nbn[x]; j < e; j++) if (M->alive[l[j]]) M->mark[l[j]] = M->stamp;
    l = M->nb[v];
    for (int j = 0, e = M->nbn[v]; j < e; j++) {
        int y = l[j];
        if (M->alive[y] && M->mark[y] != M->stamp) return 0;
    }
    return 1;
}

/* Run the reductions to a fixed point over the worklist. */
static void gm_reduce(GmMis* M) {
    while (M->wqn > 0 && !M->aborted) {
        int v = M->wq[--M->wqn];
        M->inq[v] = 0;
        if (!M->alive[v]) continue;
        int d = M->deg[v];
        if (d <= 1) { gm_take(M, v); continue; }
        if (d == 2) {
            int a = -1, b = -1;
            const int* l = M->nb[v];
            for (int j = 0, e = M->nbn[v]; j < e; j++) {
                int w = l[j];
                if (M->alive[w]) { if (a < 0) a = w; else { b = w; break; } }
            }
            if (gm_adjacent(M, a, b)) gm_take(M, v);
            else gm_fold(M, v, a, b);
            continue;
        }
        if (d <= MIS_DOM_DEG) {
            const int* l = M->nb[v];
            for (int j = 0, e = M->nbn[v]; j < e; j++) {
                int w = l[j];
                if (M->alive[w] && M->deg[w] >= d && gm_dominates(M, w, v)) {
                    gm_remove(M, w);
                    gm_push_wq(M, v);
                    break;
                }
            }
        }
    }
}

/* Greedy clique cover size of the live vertices of vs: an upper bound on the
 * independence number. Each vertex joins the largest neighbouring clique it
 * is completely adjacent to. */
static int gm_clique_cover(GmMis* M, const int* vs, int k) {
    int nc = 0;
    for (int i = 0; i < k; i++) {
        int v = vs[i];
        if (!M->alive[v]) continue;
        int nt = 0, best = -1;
        const int* l = M->nb[v];
        for (int j = 0, e = M->nbn[v]; j < e; j++) {
            int w = l[j];
            if (!M->alive[w] || M->cid[w] < 0) continue;
            int c = M->cid[w];
            if (M->ccnt[c] == 0) M->ctouch[nt++] = c;
            M->ccnt[c]++;
        }
        for (int t = 0; t < nt; t++) {
            int c = M->ctouch[t];
            if (M->ccnt[c] == M->csz[c] && (best < 0 || M->csz[c] > M->csz[best])) best = c;
            M->ccnt[c] = 0;
        }
        if (best < 0) { best = nc++; M->csz[best] = 0; }
        M->cid[v] = best;
        M->csz[best]++;
    }
    for (int i = 0; i < k; i++) M->cid[vs[i]] = -1;
    return nc;
}

/* Mirrors of v: vertices u at distance 2 with N(v) - N(u) a clique. Some
 * optimum either contains v or avoids v and all its mirrors (Fomin, Grandoni,
 * Kratsch). Writes them to buf, returns the count. */
static int gm_mirrors(GmMis* M, int v, int* buf) {
    int cnt = 0;
    M->stamp++;
    int sv = M->stamp;
    M->mark[v] = sv;
    const int* lv = M->nb[v];
    int nv = M->nbn[v];
    for (int j = 0; j < nv; j++) if (M->alive[lv[j]]) M->mark[lv[j]] = sv;
    M->stamp2++;
    int sc = M->stamp2;
    for (int j = 0; j < nv; j++) {
        int x = lv[j];
        if (!M->alive[x]) continue;
        const int* lx = M->nb[x];
        for (int t = 0, e = M->nbn[x]; t < e; t++) {
            int u = lx[t];
            if (!M->alive[u] || M->mark[u] == sv || M->mark2[u] == sc) continue;
            M->mark2[u] = sc;
            buf[cnt++] = u;
        }
    }
    /* filter: N(v) - N(u) must be a clique */
    int keep = 0;
    int* dset = buf + cnt;       /* scratch after the candidates (<= deg v) */
    for (int i = 0; i < cnt; i++) {
        int u = buf[i];
        M->stamp2++;
        int su = M->stamp2;
        const int* lu = M->nb[u];
        for (int t = 0, e = M->nbn[u]; t < e; t++) if (M->alive[lu[t]]) M->mark2[lu[t]] = su;
        int nd = 0;
        for (int j = 0; j < nv; j++) {
            int y = lv[j];
            if (M->alive[y] && M->mark2[y] != su) dset[nd++] = y;
        }
        int clique = 1;
        for (int p = 0; p < nd && clique; p++)
            for (int q = p + 1; q < nd && clique; q++)
                if (!gm_adjacent(M, dset[p], dset[q])) clique = 0;
        if (clique) buf[keep++] = u;
    }
    return keep;
}

/* Assemble this level's solution: the entries out[obase..olen) returned by
 * subproblems, plus chosen[cfrom .. cfrom+ccount), lifted through the folds
 * folds[fbase..flen) in reverse. Leaves it at out[obase..], returns its size. */
static int gm_assemble(GmMis* M, int obase, int cfrom, int ccount, int fbase) {
    int nf = M->flen - fbase;
    if (!gm_reserve(M, ccount + 2 * nf)) return -1;
    memcpy(M->out + M->olen, M->chosen + cfrom, (size_t)ccount * sizeof(int));
    M->olen += ccount;
    for (int i = obase; i < M->olen; i++) M->insol[M->out[i]] = 1;
    for (int f = M->flen - 1; f >= fbase; f--) {
        GmFold* F = &M->folds[f];
        if (M->insol[F->w]) {
            M->insol[F->w] = 0; M->insol[F->a] = 1; M->insol[F->b] = 1;
            M->out[M->olen++] = F->a; M->out[M->olen++] = F->b;
        } else {
            M->insol[F->v] = 1;
            M->out[M->olen++] = F->v;
        }
    }
    int j = obase;
    for (int i = obase; i < M->olen; i++) {
        int x = M->out[i];
        if (M->insol[x]) { M->out[j++] = x; M->insol[x] = 0; }
    }
    M->olen = j;
    return j - obase;
}

/* Exact independence number of the live part of vs if it exceeds lb: pushes
 * the optimal set (in terms of the vertices alive on entry) onto M->out and
 * returns its size; otherwise returns a value <= lb and pushes nothing. The
 * live state is unchanged on return. */
static int gm_rec(GmMis* M, const int* vs, int k, int lb) {
    if (M->aborted) return lb;
    if ((++M->nodes & 1023) == 0) {
        galg_poll();
        if (M->nodes > MIS_MAX_NODES) { M->aborted = 1; return lb; }
    }
    int tmark = M->tlen, cbase = M->clen, fbase = M->flen, obase = M->olen;
    gm_reduce(M);
    int taken = (M->clen - cbase) + (M->flen - fbase);     /* folds count 1 */
    int maxlive = k + (M->flen - fbase);
    int* live = malloc((size_t)(maxlive > 0 ? maxlive : 1) * sizeof(int));
    int result = lb;
    if (!live) { M->aborted = 1; goto done; }
    int nl = 0;
    for (int i = 0; i < k; i++) if (M->alive[vs[i]]) live[nl++] = vs[i];
    for (int f = fbase; f < M->flen; f++) if (M->alive[M->folds[f].w]) live[nl++] = M->folds[f].w;
    if (M->aborted) goto done;
    if (nl == 0) {
        if (taken > lb) result = gm_assemble(M, obase, cbase, M->clen - cbase, fbase);
        goto done;
    }
    if (taken + gm_clique_cover(M, live, nl) <= lb) goto done;

    /* components of the live part */
    {
        M->stamp++;
        int s0 = M->stamp;
        int ncomp = 0, h = 0, t = 0;
        int* order = malloc((size_t)nl * sizeof(int));
        int* cstart = malloc((size_t)(nl + 1) * sizeof(int));
        if (!order || !cstart) { free(order); free(cstart); M->aborted = 1; goto done; }
        for (int i = 0; i < nl; i++) {
            int s = live[i];
            if (M->mark[s] == s0) continue;
            cstart[ncomp++] = t;
            M->mark[s] = s0; order[t++] = s;
            while (h < t) {
                int x = order[h++];
                const int* l = M->nb[x];
                for (int j = 0, e = M->nbn[x]; j < e; j++) {
                    int y = l[j];
                    if (M->alive[y] && M->mark[y] != s0) { M->mark[y] = s0; order[t++] = y; }
                }
            }
        }
        cstart[ncomp] = t;
        if (ncomp > 1) {
            /* independent subproblems, each solved exactly (lb -1) */
            int total = taken;
            for (int c = 0; c < ncomp && !M->aborted; c++)
                total += gm_rec(M, order + cstart[c], cstart[c + 1] - cstart[c], -1);
            free(order); free(cstart);
            if (!M->aborted && total > lb) result = gm_assemble(M, obase, cbase, M->clen - cbase, fbase);
            goto done;
        }
        free(order); free(cstart);
    }

    /* branch on a maximum-degree vertex: take it, or drop it with its mirrors */
    {
        int v = live[0];
        for (int i = 1; i < nl; i++) if (M->deg[live[i]] > M->deg[v]) v = live[i];
        int best = lb, bstart = obase, blen = 0;
        int nlc = M->clen - cbase;          /* chosen at this level so far */
        /* drop v and its mirrors */
        {
            int* mir = malloc((size_t)(M->nid + M->deg[v] + 1) * sizeof(int));
            if (!mir) { M->aborted = 1; goto done; }
            int nm = gm_mirrors(M, v, mir);
            int m2 = M->tlen;
            gm_remove(M, v);
            for (int i = 0; i < nm; i++) if (M->alive[mir[i]]) gm_remove(M, mir[i]);
            free(mir);
            int need = best - taken;
            int r = gm_rec(M, live, nl, need);
            if (!M->aborted && r > need) {
                /* child solution sits at out[bstart..]; lift it */
                int sz = gm_assemble(M, bstart, cbase, nlc, fbase);
                if (sz >= 0) { best = sz; blen = sz; }
            }
            gm_undo(M, m2);
        }
        /* take v */
        if (!M->aborted) {
            int m2 = M->tlen, c2 = M->clen;
            gm_take(M, v);
            int need = best - taken - 1;
            int base2 = M->olen;
            int r = gm_rec(M, live, nl, need);
            if (!M->aborted && r > need) {
                /* move the child's solution down over the old incumbent */
                memmove(M->out + bstart, M->out + base2, (size_t)r * sizeof(int));
                M->olen = bstart + r;
                int sz = gm_assemble(M, bstart, cbase, nlc + 1, fbase);
                if (sz >= 0) { best = sz; blen = sz; }
            }
            gm_undo(M, m2);
            M->clen = c2;
        }
        if (!M->aborted && best > lb) { result = best; M->olen = bstart + blen; }
    }
done:
    free(live);
    while (M->wqn > 0) M->inq[M->wq[--M->wqn]] = 0;
    gm_undo(M, tmark);
    M->clen = cbase;
    M->flen = fbase;
    if (result <= lb || M->aborted) { M->olen = obase; if (result > lb) result = lb; }
    return result;
}

/* Greedy minimum-degree independent set of the component vs (static graph):
 * a quick lower bound. Writes it to out, returns its size. */
static int gm_greedy(const GalgUG* u, const int* vs, int k, int* out, char* dead, int* degw) {
    int cnt = 0;
    /* lazy binary min-heap of (deg, v) */
    long hc = (long)k + 1, hn = 0;
    for (int i = 0; i < k; i++) { int v = vs[i]; hc += u->off[v + 1] - u->off[v]; }
    long* hk = malloc((size_t)hc * sizeof(long));
    if (!hk) return 0;
    for (int i = 0; i < k; i++) { int v = vs[i]; degw[v] = u->off[v + 1] - u->off[v]; dead[v] = 0; }
#define GM_HPUSH(key) do { long x_ = (key), i_ = hn++; \
        while (i_ > 0 && hk[(i_ - 1) / 2] > x_) { hk[i_] = hk[(i_ - 1) / 2]; i_ = (i_ - 1) / 2; } \
        hk[i_] = x_; } while (0)
    for (int i = 0; i < k; i++) GM_HPUSH(((long)degw[vs[i]] << 32) | vs[i]);
    while (hn > 0) {
        long top = hk[0], last = hk[--hn];
        long i = 0;
        for (;;) {
            long c = 2 * i + 1;
            if (c >= hn) break;
            if (c + 1 < hn && hk[c + 1] < hk[c]) c++;
            if (hk[c] >= last) break;
            hk[i] = hk[c]; i = c;
        }
        if (hn > 0) hk[i] = last;
        int v = (int)(top & 0xffffffffL), d = (int)(top >> 32);
        if (dead[v] || d != degw[v]) continue;
        out[cnt++] = v;
        dead[v] = 1;
        for (int j = u->off[v]; j < u->off[v + 1]; j++) {
            int x = u->adj[j];
            if (dead[x]) continue;
            dead[x] = 1;
            for (int t = u->off[x]; t < u->off[x + 1]; t++) {
                int y = u->adj[t];
                if (!dead[y]) { degw[y]--; GM_HPUSH(((long)degw[y] << 32) | y); }
            }
        }
    }
#undef GM_HPUSH
    free(hk);
    return cnt;
}

static void gm_free(GmMis* M) {
    if (M->nb) for (int i = 0; i < M->cap; i++) free(M->nb[i]);
    free(M->nb); free(M->nbn); free(M->nbcap); free(M->alive); free(M->deg); free(M->tr);
    free(M->chosen); free(M->folds); free(M->out); free(M->wq); free(M->inq); free(M->mark);
    free(M->mark2); free(M->cid); free(M->csz); free(M->ccnt); free(M->ctouch); free(M->insol);
}

/* Maximum independent set of the simple graph u: fills out[] (vertex ids,
 * any order) and returns its size, or -1 on budget / allocation failure. */
int galg_max_independent_set(const GalgUG* u, int* out) {
    int n = u->n;
    if (n == 0) return 0;
    GmMis M;
    memset(&M, 0, sizeof(M));
    int cap = n + n / 2 + 4;
    size_t cc = (size_t)cap;
    M.cap = cap; M.nid = n;
    M.nb = calloc(cc, sizeof(int*));
    M.nbn = calloc(cc, sizeof(int));
    M.nbcap = calloc(cc, sizeof(int));
    M.alive = calloc(cc, 1);
    M.deg = calloc(cc, sizeof(int));
    M.tr = malloc(2 * cc * sizeof(int));
    M.chosen = malloc(cc * sizeof(int));
    M.folds = malloc(cc * sizeof(GmFold));
    M.wq = malloc(cc * sizeof(int));
    M.inq = calloc(cc, 1);
    M.mark = calloc(cc, sizeof(int));
    M.mark2 = calloc(cc, sizeof(int));
    M.cid = malloc(cc * sizeof(int));
    M.csz = malloc(cc * sizeof(int));
    M.ccnt = calloc(cc, sizeof(int));
    M.ctouch = malloc(cc * sizeof(int));
    M.insol = calloc(cc, 1);
    int* comp = malloc((size_t)n * sizeof(int));
    int* cst = malloc(((size_t)n + 1) * sizeof(int));
    int* tmp = malloc((size_t)n * sizeof(int));
    int* degw = malloc((size_t)n * sizeof(int));
    char* dead = malloc((size_t)n);
    int total = -1;
    if (!M.nb || !M.nbn || !M.nbcap || !M.alive || !M.deg || !M.tr || !M.chosen || !M.folds
        || !M.wq || !M.inq || !M.mark || !M.mark2 || !M.cid || !M.csz || !M.ccnt || !M.ctouch
        || !M.insol || !comp || !cst || !tmp || !degw || !dead || !gm_reserve(&M, n)) goto out;
    for (int v = 0; v < cap; v++) M.cid[v] = -1;
    for (int v = 0; v < n; v++) {
        int d = u->off[v + 1] - u->off[v];
        M.nb[v] = malloc((size_t)(d > 0 ? d : 1) * sizeof(int));
        if (!M.nb[v]) goto out;
        memcpy(M.nb[v], u->adj + u->off[v], (size_t)d * sizeof(int));
        M.nbn[v] = M.nbcap[v] = d;
        M.alive[v] = 1; M.deg[v] = d;
    }
    /* top-level components (static graph) */
    int nc = 0;
    {
        char* seen = calloc((size_t)n, 1);
        if (!seen) goto out;
        int t = 0;
        for (int s = 0; s < n; s++) {
            if (seen[s]) continue;
            cst[nc++] = t;
            int h = t;
            seen[s] = 1; comp[t++] = s;
            while (h < t) {
                int x = comp[h++];
                for (int j = u->off[x]; j < u->off[x + 1]; j++)
                    if (!seen[u->adj[j]]) { seen[u->adj[j]] = 1; comp[t++] = u->adj[j]; }
            }
        }
        cst[nc] = t;
        free(seen);
    }
    total = 0;
    int olen = 0;
    for (int c = 0; c < nc; c++) {
        int* vs = comp + cst[c];
        int k = cst[c + 1] - cst[c];
        long edges2 = 0;
        for (int i = 0; i < k; i++) edges2 += u->off[vs[i] + 1] - u->off[vs[i]];
        double dens = k > 1 ? (double)edges2 / ((double)k * (k - 1)) : 0;
        int r;
        /* measured crossover: the complement-clique bound wins from about
         * average degree 8 up; branch and reduce below */
        double avgdeg = k > 0 ? (double)edges2 / k : 0;
        if (k > 1 && k <= MIS_DENSE_MAX && (avgdeg >= 8 || dens >= 0.05)) {
            r = galg_mis_dense(u, vs, k, tmp);
            if (r < 0) { total = -1; break; }
        } else {
            int g0 = gm_greedy(u, vs, k, tmp, dead, degw);
            for (int i = 0; i < k; i++) gm_push_wq(&M, vs[i]);
            M.olen = 0;
            int rr = gm_rec(&M, vs, k, g0);
            if (M.aborted) { total = -1; break; }
            if (rr > g0) { memcpy(tmp, M.out, (size_t)rr * sizeof(int)); r = rr; }
            else r = g0;
        }
        memcpy(out + olen, tmp, (size_t)r * sizeof(int));
        olen += r; total += r;
    }
out:
    gm_free(&M);
    free(comp); free(cst); free(tmp); free(degw); free(dead);
    return total;
}

/* ---- Builtins --------------------------------------------------------------- */

Expr* builtin_find_independent_vertex_set(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int kmin = 1, kmax = 0x3fffffff;
    long count = 1;
    if (argc >= 2 && !galg_parse_size_spec(res->data.function.args[1], &kmin, &kmax)) return NULL;
    if (argc == 3 && !galg_parse_count(res->data.function.args[2], &count)) return NULL;
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    Expr* out = NULL;
    if (u->n == 0) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else if (argc == 1) {
        int* s = malloc((size_t)u->n * sizeof(int));
        if (s) {
            int k = galg_max_independent_set(u, s);
            if (k >= 0) {
                Expr* one = galg_vertex_list(g, s, k, 1);
                out = expr_new_function(expr_new_symbol(SYM_List), &one, 1);
            }
        }
        free(s);
    } else {
        out = galg_clique_spec_query(g, u, 1, kmin, kmax, count);
    }
    galg_ug_free(u);
    return out;
}

Expr* builtin_find_vertex_cover(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    GalgUG* u = galg_ug_from_graph(g);
    if (!u) return NULL;
    int n = u->n;
    int* s = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    char* in = calloc((size_t)(n > 0 ? n : 1), 1);
    Expr* out = NULL;
    if (s && in) {
        int k = galg_max_independent_set(u, s);
        if (k >= 0) {
            for (int i = 0; i < k; i++) in[s[i]] = 1;
            int c = 0;
            for (int v = 0; v < n; v++) if (!in[v]) s[c++] = v;
            out = galg_vertex_list(g, s, c, 0);
        }
    }
    free(s); free(in); galg_ug_free(u);
    return out;
}

/* Vertex positions of the elements of list (a vertex list argument): returns
 * the count, or -1 if list is not a list or some element is not a vertex. */
static long gm_vertex_positions(const Expr* g, const Expr* arg, int** pos) {
    Expr* owned = NULL;
    const Expr* lst = galg_plain_list(arg, &owned);
    *pos = NULL;
    if (!lst) return -1;
    long k = (long)lst->data.function.arg_count;
    *pos = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
    if (!*pos) { expr_free(owned); return -1; }
    for (long i = 0; i < k; i++) {
        (*pos)[i] = galg_vertex_arg(g, lst->data.function.args[i]);
        if ((*pos)[i] < 0) { free(*pos); *pos = NULL; expr_free(owned); return -1; }
    }
    expr_free(owned);
    return k;
}

static Expr* gm_vertex_set_q(Expr* res, int cover) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return galg_truth(0);
    int* pos;
    long k = gm_vertex_positions(g, res->data.function.args[1], &pos);
    if (k < 0) return galg_truth(0);
    int n = galg_nv(g);
    long m = galg_ne(g);
    char* in = calloc((size_t)(n > 0 ? n : 1), 1);
    if (!in) { free(pos); return NULL; }
    for (long i = 0; i < k; i++) in[pos[i]] = 1;
    int ok = 1;
    for (long e = 0; e < m && ok; e++) {
        int a = in[eu[e]], b = in[ev[e]];
        if (cover ? !(a || b) : (a && b)) ok = 0;
    }
    free(in); free(pos);
    return galg_truth(ok);
}

Expr* builtin_independent_vertex_set_q(Expr* res) { return gm_vertex_set_q(res, 0); }
Expr* builtin_vertex_cover_q(Expr* res) { return gm_vertex_set_q(res, 1); }

/* Endpoint positions of an edge element if it is an edge of g: accepts
 * DirectedEdge / Rule (directed) and UndirectedEdge / TwoWayRule. */
static int gm_edge_of(const Expr* g, const Expr* e, int* a, int* b) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2
        || e->data.function.head->type != EXPR_SYMBOL) return 0;
    const char* h = e->data.function.head->data.symbol.name;
    int directed;
    if (h == SYM_DirectedEdge || h == SYM_Rule) directed = 1;
    else if (h == SYM_UndirectedEdge || galg_is_symbol(e->data.function.head, "TwoWayRule")) directed = 0;
    else return 0;
    if (graph_has_edge(g, e->data.function.args[0], e->data.function.args[1], directed) != 1) return 0;
    *a = galg_vertex_arg(g, e->data.function.args[0]);
    *b = galg_vertex_arg(g, e->data.function.args[1]);
    return *a >= 0 && *b >= 0;
}

static Expr* gm_edge_set_q(Expr* res, int cover) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return galg_truth(0);
    if (!graph_is_list(res->data.function.args[1])) return galg_truth(0);
    const Expr* lst = res->data.function.args[1];
    int n = galg_nv(g);
    char* hit = calloc((size_t)(n > 0 ? n : 1), 1);
    if (!hit) return NULL;
    int ok = 1;
    for (size_t i = 0; i < lst->data.function.arg_count && ok; i++) {
        int a, b;
        if (!gm_edge_of(g, lst->data.function.args[i], &a, &b)) { ok = 0; break; }
        if (!cover && (hit[a] || hit[b])) ok = 0;
        hit[a] = hit[b] = 1;
    }
    if (ok && cover) for (int v = 0; v < n; v++) if (!hit[v]) { ok = 0; break; }
    free(hit);
    return galg_truth(ok);
}

Expr* builtin_independent_edge_set_q(Expr* res) { return gm_edge_set_q(res, 0); }
Expr* builtin_edge_cover_q(Expr* res) { return gm_edge_set_q(res, 1); }
