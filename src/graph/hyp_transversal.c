/* hyp_transversal.c - hypergraph transversals (hitting sets).
 *
 *   TransversalHypergraph[h]      every MINIMAL transversal of h
 *   FindMinimumTransversal[h]     one transversal of minimum cardinality
 *
 * A transversal (hitting set, vertex cover) of h is a vertex set meeting every
 * hyperedge; it is minimal when no proper subset is one. The minimal ones form
 * the transversal hypergraph Tr(h) (Berge), the object the Function
 * Repository's TransversalHypergraph computes. Both heads accept a Hypergraph
 * or, for FR compatibility, a plain List of hyperedge Lists; the result's shape
 * follows the input (Hypergraph[VertexList[h], Tr] vs. the List Tr).
 *
 * Output order (deterministic, and Sort's order when the vertices are integers
 * in increasing VertexList order): each transversal lists its vertices in
 * VertexList order; transversals are sorted by size, then lexicographically by
 * VertexList positions.
 *
 * Conventions: no hyperedges -> Tr = {{}} (the empty set hits everything); an
 * empty hyperedge -> Tr = {} (nothing hits it), and FindMinimumTransversal is
 * then left unevaluated.
 *
 * ALGORITHMS.
 *   Tr(h): MMCS (Murakami & Uno, "Efficient algorithms for dualizing large-
 *   scale hypergraphs", 2014), the practical state of the art. Depth-first over
 *   partial solutions S that are always minimal (every s in S has a CRITICAL
 *   hyperedge -- one met by S only in s); at each node the uncovered hyperedge
 *   with the fewest remaining candidates is branched on. Crit bookkeeping is
 *   O(deg v) per add/remove: per hyperedge we keep |F ∩ S| and the SUM of S's
 *   members in F, so when |F ∩ S| == 1 the sum IS the unique member. Run
 *   iteratively (explicit frames), so a huge transversal cannot overflow the C
 *   stack.
 *   Minimum: branch and bound on the uncovered hyperedge with fewest open
 *   vertices (include v; then exclude it for later siblings), pruned by a
 *   greedy disjoint-hyperedge packing lower bound, seeded by the greedy
 *   max-coverage upper bound.
 *
 * BUDGETS (correct-or-unevaluated). Tr(h) can be exponentially large, and the
 * minimum hitting set is NP-hard, so both searches count nodes and give up --
 * returning unevaluated, never a partial or unproven answer -- past
 * HYP_TR_MAX_NODES search nodes (or HYP_TR_MAX_OUT output entries for Tr).
 * Node counts, not wall-clock, keep the verdict machine-independent; both poll
 * tc_check_deadline() every 4096 nodes so TimeConstrained interrupts them.
 */

#include "graph_hyper.h"
#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "sym_names.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define HYP_TR_MAX_NODES 50000000L
#define HYP_TR_MAX_OUT   20000000L
#define HYP_MIN_MAX_NODES 20000000L

/* ---- shared input handling ------------------------------------------------ */

static const Expr* tr_arg(const Expr* a, Expr** owned, int* was_list) {
    *owned = NULL;
    *was_list = 0;
    if (hypergraph_is_valid(a)) return a;
    if (!graph_is_list(a)) return NULL;
    /* FR form: a List of hyperedges; a non-List element (an "isolated vertex")
     * constrains nothing and is skipped. */
    size_t m = a->data.function.arg_count, c = 0;
    Expr** es = malloc((m ? m : 1) * sizeof(Expr*));
    if (!es) return NULL;
    for (size_t j = 0; j < m; j++)
        if (graph_is_list(a->data.function.args[j])) es[c++] = expr_copy(a->data.function.args[j]);
    Expr* el = expr_new_function(expr_new_symbol(SYM_List), es, c);
    free(es);
    Expr* call = expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), &el, 1);
    Expr* h = evaluate(call);
    expr_free(call);
    if (!hypergraph_is_valid(h)) { expr_free(h); return NULL; }
    *owned = h;
    *was_list = 1;
    return h;
}

/* ---- MMCS state ------------------------------------------------------------ */

typedef struct {
    const HypView* v;
    int* cnt;          /* |F ∩ S|                                           */
    int64_t* sum;      /* sum of S's members in F                            */
    int* crit;         /* crit[u]: # hyperedges where u is S's only member   */
    int zc;            /* # members of S with crit == 0                      */
    int* uncov; int* upos; int nu;   /* uncovered hyperedges, swap-removed  */
    unsigned char* incand;
    int* S; int ns;
} TrState;

static void tr_add(TrState* t, int x) {
    const HypView* v = t->v;
    t->S[t->ns++] = x;
    t->crit[x] = 0; t->zc++;
    for (int q = v->voff[x]; q < v->voff[x + 1]; q++) {
        int f = v->ve[q];
        if (t->cnt[f] == 0) {
            int p = t->upos[f], last = t->uncov[--t->nu];
            t->uncov[p] = last; t->upos[last] = p; t->upos[f] = -1;
            if (t->crit[x]++ == 0) t->zc--;
        } else if (t->cnt[f] == 1) {
            int u = (int)t->sum[f];
            if (--t->crit[u] == 0) t->zc++;
        }
        t->cnt[f]++; t->sum[f] += x;
    }
}

static void tr_remove(TrState* t, int x) {
    const HypView* v = t->v;
    for (int q = v->voff[x + 1] - 1; q >= v->voff[x]; q--) {
        int f = v->ve[q];
        t->cnt[f]--; t->sum[f] -= x;
        if (t->cnt[f] == 0) {
            t->upos[f] = t->nu; t->uncov[t->nu++] = f;
            if (--t->crit[x] == 0) t->zc++;
        } else if (t->cnt[f] == 1) {
            int u = (int)t->sum[f];
            if (t->crit[u]++ == 0) t->zc--;
        }
    }
    t->ns--;
    t->zc--;                                  /* crit[x] is 0 again */
}

/* Uncovered hyperedge with fewest candidate vertices; *nc gets that count. */
static int tr_choose(const TrState* t, int* nc) {
    const HypView* v = t->v;
    int best = -1, bc = INT32_MAX;
    for (int a = 0; a < t->nu && bc > 0; a++) {
        int f = t->uncov[a], c = 0;
        for (int q = v->soff[f]; q < v->soff[f + 1] && c < bc; q++) c += t->incand[v->sv[q]];
        if (c < bc) { bc = c; best = f; }
    }
    *nc = bc;
    return best;
}

typedef struct { int* a; size_t n, cap; } IVec;
static int iv_push(IVec* b, int x) {
    if (b->n == b->cap) {
        size_t nc = b->cap ? b->cap * 2 : 256;
        int* na = realloc(b->a, nc * sizeof(int));
        if (!na) return 0;
        b->a = na; b->cap = nc;
    }
    b->a[b->n++] = x;
    return 1;
}

static int cmp_int(const void* a, const void* b) {
    int x = *(const int*)a, y = *(const int*)b;
    return (x > y) - (x < y);
}

/* Sorting transversals: records (offset, length) into one int pool. */
static const int* g_pool;
typedef struct { size_t off; int len; } TrRec;
static int cmp_rec(const void* a, const void* b) {
    const TrRec* x = a; const TrRec* y = b;
    if (x->len != y->len) return x->len < y->len ? -1 : 1;
    for (int i = 0; i < x->len; i++) {
        int p = g_pool[x->off + (size_t)i], q = g_pool[y->off + (size_t)i];
        if (p != q) return p < q ? -1 : 1;
    }
    return 0;
}

/* Run MMCS. Returns 1 with pool/recs filled, 0 on budget/allocation failure. */
static int mmcs(const HypView* v, IVec* pool, TrRec** recs_out, size_t* nrec_out) {
    int n = v->n, m = v->m;
    size_t nn = (size_t)(n > 0 ? n : 1), mm = (size_t)(m > 0 ? m : 1);
    TrState t;
    memset(&t, 0, sizeof(t));
    t.v = v;
    t.cnt = calloc(mm, sizeof(int));
    t.sum = calloc(mm, sizeof(int64_t));
    t.crit = calloc(nn, sizeof(int));
    t.uncov = malloc(mm * sizeof(int));
    t.upos = malloc(mm * sizeof(int));
    t.incand = calloc(nn, 1);
    t.S = malloc(nn * sizeof(int));
    /* Frames: per depth, the start/length of its candidate list in cstack and
     * the next index to try. Depth <= |S| + 1 <= n + 1. */
    int* fstart = malloc((nn + 1) * sizeof(int));
    int* flen = malloc((nn + 1) * sizeof(int));
    int* fidx = malloc((nn + 1) * sizeof(int));
    IVec cstack = {0};
    size_t nrec = 0, caprec = 0;
    TrRec* recs = NULL;
    int ok = t.cnt && t.sum && t.crit && t.uncov && t.upos && t.incand && t.S
             && fstart && flen && fidx;
    if (ok) {
        for (int f = 0; f < m; f++) { t.uncov[f] = f; t.upos[f] = f; }
        t.nu = m;
        for (int x = 0; x < n; x++) t.incand[x] = (unsigned char)(v->voff[x + 1] > v->voff[x]);
    }
    long nodes = 0;
    int depth = 0;
    /* push_frame: pick the branching hyperedge for the current S. Returns 1
     * if a frame was pushed, 0 if S is a leaf (emitted or dead). */
#define EMIT_S() do {                                                        \
        if (nrec == caprec) {                                                \
            size_t nc = caprec ? caprec * 2 : 256;                           \
            TrRec* nr = realloc(recs, nc * sizeof(TrRec));                   \
            if (!nr) { ok = 0; break; }                                      \
            recs = nr; caprec = nc;                                          \
        }                                                                    \
        recs[nrec].off = pool->n; recs[nrec].len = t.ns;                     \
        for (int _i = 0; _i < t.ns && ok; _i++) ok = iv_push(pool, t.S[_i]); \
        if (ok) qsort(pool->a + recs[nrec].off, (size_t)t.ns, sizeof(int), cmp_int); \
        nrec++;                                                              \
        if ((long)pool->n > HYP_TR_MAX_OUT) ok = 0;                          \
    } while (0)

    if (ok) {
        if (t.nu == 0) {
            EMIT_S();
        } else {
            int nc, f = tr_choose(&t, &nc);
            if (nc > 0) {
                fstart[0] = (int)cstack.n; flen[0] = 0; fidx[0] = 0;
                for (int q = v->soff[f]; q < v->soff[f + 1] && ok; q++)
                    if (t.incand[v->sv[q]]) { ok = iv_push(&cstack, v->sv[q]); flen[0]++; }
                for (int q = 0; q < flen[0]; q++) t.incand[cstack.a[fstart[0] + q]] = 0;
                depth = 1;
            }
        }
    }
    while (ok && depth > 0) {
        int d = depth - 1;
        if (fidx[d] > 0) {                     /* undo the previous choice */
            int x = cstack.a[fstart[d] + fidx[d] - 1];
            tr_remove(&t, x);
            t.incand[x] = 1;
        }
        if (fidx[d] == flen[d]) {              /* frame exhausted */
            cstack.n = (size_t)fstart[d];
            depth--;
            continue;
        }
        int x = cstack.a[fstart[d] + fidx[d]];
        fidx[d]++;
        if (++nodes > HYP_TR_MAX_NODES) { ok = 0; break; }
        if ((nodes & 0xFFF) == 0) tc_check_deadline();
        tr_add(&t, x);
        if (t.zc != 0) continue;               /* not minimal: prune */
        if (t.nu == 0) { EMIT_S(); continue; }
        int nc, f = tr_choose(&t, &nc);
        if (nc == 0) continue;                 /* some hyperedge unhittable */
        fstart[depth] = (int)cstack.n; flen[depth] = 0; fidx[depth] = 0;
        for (int q = v->soff[f]; q < v->soff[f + 1] && ok; q++)
            if (t.incand[v->sv[q]]) { ok = iv_push(&cstack, v->sv[q]); flen[depth]++; }
        for (int q = 0; q < flen[depth]; q++) t.incand[cstack.a[fstart[depth] + q]] = 0;
        depth++;
    }
#undef EMIT_S
    free(t.cnt); free(t.sum); free(t.crit); free(t.uncov); free(t.upos);
    free(t.incand); free(t.S); free(fstart); free(flen); free(fidx); free(cstack.a);
    if (!ok) { free(recs); return 0; }
    g_pool = pool->a;
    if (nrec > 1) qsort(recs, nrec, sizeof(TrRec), cmp_rec);
    *recs_out = recs; *nrec_out = nrec;
    return 1;
}

Expr* builtin_transversal_hypergraph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned; int was_list;
    const Expr* h = tr_arg(res->data.function.args[0], &owned, &was_list);
    HypView v;
    if (!h || !hyp_view(h, &v, 1)) { expr_free(owned); return NULL; }
    IVec pool = {0};
    TrRec* recs = NULL; size_t nrec = 0;
    if (!mmcs(&v, &pool, &recs, &nrec)) { free(pool.a); expr_free(owned); return NULL; }
    Expr** out = malloc((nrec ? nrec : 1) * sizeof(Expr*));
    Expr** row = malloc((size_t)(v.n > 0 ? v.n : 1) * sizeof(Expr*));
    Expr* lh = expr_new_symbol(SYM_List);
    Expr* const* vx = v.verts->data.function.args;
    for (size_t r = 0; r < nrec; r++) {
        for (int i = 0; i < recs[r].len; i++) row[i] = expr_copy(vx[pool.a[recs[r].off + (size_t)i]]);
        out[r] = expr_new_function(expr_copy(lh), row, (size_t)recs[r].len);
    }
    Expr* tr = expr_new_function(lh, out, nrec);
    free(out); free(row); free(recs); free(pool.a);
    Expr* result = tr;
    if (!was_list) {
        Expr* args[2] = { expr_copy((Expr*)v.verts), tr };
        result = expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), args, 2);
    }
    expr_free(owned);
    return result;
}

/* ---- minimum transversal (branch and bound) -------------------------------- */

typedef struct {
    const HypView* v;
    int* cnt; int* uncov; int* upos; int nu;
    unsigned char* excl;
    int* S; int ns;
    int* best; int nbest;
    int* mark; int stamp;
    int* odeg; int* touched; int* dcount;   /* lower-bound scratch */
    int* ekey; int* eord;
    double* ey; double* slack;
    long nodes; int aborted;
} MinState;

static void mn_add(MinState* t, int x) {
    const HypView* v = t->v;
    t->S[t->ns++] = x;
    for (int q = v->voff[x]; q < v->voff[x + 1]; q++) {
        int f = v->ve[q];
        if (t->cnt[f]++ == 0) {
            int p = t->upos[f], last = t->uncov[--t->nu];
            t->uncov[p] = last; t->upos[last] = p; t->upos[f] = -1;
        }
    }
}
static void mn_remove(MinState* t, int x) {
    const HypView* v = t->v;
    for (int q = v->voff[x + 1] - 1; q >= v->voff[x]; q--) {
        int f = v->ve[q];
        if (--t->cnt[f] == 0) { t->upos[f] = t->nu; t->uncov[t->nu++] = f; }
    }
    t->ns--;
}

/* Lower bound on the vertices still needed to hit the uncovered hyperedges U,
 * the larger of two:
 *   degree bound  -- with open-vertex degrees d_1 >= d_2 >= ... counted within
 *                    U, the least k with d_1 + ... + d_k >= |U|;
 *   packing bound -- a set of pairwise-disjoint hyperedges of U (over open
 *                    vertices) needs one vertex each. Built greedily with the
 *                    hyperedges taken in increasing order of their total open
 *                    degree, since a hyperedge of rarely-used vertices blocks
 *                    the fewest others -- the ordering that makes the greedy
 *                    packing close to maximum in practice.
 * Scratch (odeg, dcount) is left zeroed. */
static int g_mn_key_base;
static const int* g_mn_key;
static int cmp_by_key(const void* a, const void* b) {
    int x = g_mn_key[*(const int*)a - g_mn_key_base], y = g_mn_key[*(const int*)b - g_mn_key_base];
    if (x != y) return x < y ? -1 : 1;
    return (*(const int*)a > *(const int*)b) - (*(const int*)a < *(const int*)b);
}

static int mn_lower_bound(MinState* t) {
    const HypView* v = t->v;
    int nt = 0, maxd = 0;
    for (int a = 0; a < t->nu; a++) {
        int f = t->uncov[a];
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++) {
            int x = v->sv[q];
            if (t->excl[x]) continue;
            if (t->odeg[x]++ == 0) t->touched[nt++] = x;
        }
    }
    for (int i = 0; i < nt; i++) {
        int d = t->odeg[t->touched[i]];
        if (d > maxd) maxd = d;
        t->dcount[d]++;
    }
    long need = t->nu;
    int kdeg = 0;
    for (int d = maxd; d >= 1 && need > 0; d--) {
        long c = t->dcount[d];
        if (c * d >= need) { kdeg += (int)((need + d - 1) / d); need = 0; }
        else { need -= c * d; kdeg += (int)c; }
    }
    for (int d = 1; d <= maxd; d++) t->dcount[d] = 0;
    if (need > 0) {                                    /* an unhittable edge */
        for (int i = 0; i < nt; i++) t->odeg[t->touched[i]] = 0;
        return INT32_MAX / 2;
    }
    /* Packing over U sorted by total open degree. */
    for (int a = 0; a < t->nu; a++) {
        int f = t->uncov[a], w = 0;
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++)
            if (!t->excl[v->sv[q]]) w += t->odeg[v->sv[q]];
        t->ekey[f] = w;
        t->eord[a] = f;
    }
    g_mn_key = t->ekey; g_mn_key_base = 0;
    qsort(t->eord, (size_t)t->nu, sizeof(int), cmp_by_key);
    for (int i = 0; i < nt; i++) t->odeg[t->touched[i]] = 0;
    if (++t->stamp == INT32_MAX) {
        for (int i = 0; i < v->n; i++) t->mark[i] = 0;
        t->stamp = 1;
    }
    int p = 0;
    for (int a = 0; a < t->nu; a++) {
        int f = t->eord[a], hit = 0;
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++) {
            int x = v->sv[q];
            if (!t->excl[x] && t->mark[x] == t->stamp) { hit = 1; break; }
        }
        if (hit) continue;
        p++;
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++) t->mark[v->sv[q]] = t->stamp;
    }
    int best = p > kdeg ? p : kdeg;

    /* Fractional packing (LP dual) bound: y_F = min over open v in F of
     * 1/d(v) is dual-feasible (each vertex's load sum_{F ni v} y_F <= 1), then
     * one dual-ascent pass raises each y_F by its vertices' least slack, in the
     * same low-degree-first order. Any feasible y bounds the optimum below by
     * sum y_F. */
    nt = 0;
    for (int a = 0; a < t->nu; a++) {
        int f = t->uncov[a];
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++) {
            int x = v->sv[q];
            if (!t->excl[x] && t->odeg[x]++ == 0) t->touched[nt++] = x;
        }
    }
    double tot = 0.0;
    for (int a = 0; a < t->nu; a++) {
        int f = t->eord[a];
        double y = 2.0;
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++) {
            int x = v->sv[q];
            if (!t->excl[x]) { double c = 1.0 / t->odeg[x]; if (c < y) y = c; }
        }
        t->ey[f] = y;
        tot += y;
    }
    for (int i = 0; i < nt; i++) t->slack[t->touched[i]] = 1.0;
    for (int a = 0; a < t->nu; a++) {
        int f = t->eord[a];
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++)
            if (!t->excl[v->sv[q]]) t->slack[v->sv[q]] -= t->ey[f];
    }
    for (int a = 0; a < t->nu; a++) {
        int f = t->eord[a];
        double sl = 2.0;
        for (int q = v->soff[f]; q < v->soff[f + 1]; q++)
            if (!t->excl[v->sv[q]] && t->slack[v->sv[q]] < sl) sl = t->slack[v->sv[q]];
        if (sl > 1e-12) {
            tot += sl;
            for (int q = v->soff[f]; q < v->soff[f + 1]; q++)
                if (!t->excl[v->sv[q]]) t->slack[v->sv[q]] -= sl;
        }
    }
    for (int i = 0; i < nt; i++) t->odeg[t->touched[i]] = 0;
    int kf = (int)(tot - 1e-9);
    if ((double)kf < tot - 1e-9) kf++;                 /* ceil, float-safe */
    return kf > best ? kf : best;
}

static void mn_search(MinState* t) {
    const HypView* v = t->v;
    if (t->aborted) return;
    if (++t->nodes > HYP_MIN_MAX_NODES) { t->aborted = 1; return; }
    if ((t->nodes & 0xFFF) == 0) tc_check_deadline();
    if (t->nu == 0) {
        memcpy(t->best, t->S, (size_t)t->ns * sizeof(int));
        t->nbest = t->ns;
        return;
    }
    if (t->ns + 1 >= t->nbest) return;                 /* need >= 1 more   */
    if (t->ns + mn_lower_bound(t) >= t->nbest) return;
    /* Branch on the uncovered hyperedge with the fewest open vertices. */
    int best = -1, bc = INT32_MAX;
    for (int a = 0; a < t->nu && bc > 1; a++) {
        int f = t->uncov[a], c = 0;
        for (int q = v->soff[f]; q < v->soff[f + 1] && c < bc; q++) c += !t->excl[v->sv[q]];
        if (c < bc) { bc = c; best = f; }
    }
    if (bc == 0) return;                               /* dead branch */
    int nx = 0;
    int* xs = malloc((size_t)bc * sizeof(int));
    if (!xs) { t->aborted = 1; return; }
    for (int q = v->soff[best]; q < v->soff[best + 1]; q++)
        if (!t->excl[v->sv[q]]) xs[nx++] = v->sv[q];
    for (int i = 0; i < nx && !t->aborted; i++) {
        mn_add(t, xs[i]);
        mn_search(t);
        mn_remove(t, xs[i]);
        t->excl[xs[i]] = 1;                            /* later siblings skip it */
    }
    for (int i = 0; i < nx; i++) t->excl[xs[i]] = 0;
    free(xs);
}

Expr* builtin_find_minimum_transversal(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* owned; int was_list;
    const Expr* h = tr_arg(res->data.function.args[0], &owned, &was_list);
    HypView v;
    if (!h || !hyp_view(h, &v, 1)) { expr_free(owned); return NULL; }
    int n = v.n, m = v.m;
    for (int f = 0; f < m; f++)
        if (v.soff[f + 1] == v.soff[f]) { expr_free(owned); return NULL; }  /* unhittable */
    size_t nn = (size_t)(n > 0 ? n : 1), mm = (size_t)(m > 0 ? m : 1);
    MinState t;
    memset(&t, 0, sizeof(t));
    t.v = &v;
    t.cnt = calloc(mm, sizeof(int));
    t.uncov = malloc(mm * sizeof(int));
    t.upos = malloc(mm * sizeof(int));
    t.excl = calloc(nn, 1);
    t.S = malloc(nn * sizeof(int));
    t.best = malloc(nn * sizeof(int));
    t.mark = calloc(nn, sizeof(int));
    t.odeg = calloc(nn, sizeof(int));
    t.touched = malloc(nn * sizeof(int));
    t.dcount = calloc(mm + 2, sizeof(int));
    t.ekey = malloc(mm * sizeof(int));
    t.eord = malloc(mm * sizeof(int));
    t.ey = malloc(mm * sizeof(double));
    t.slack = malloc(nn * sizeof(double));
    int* gain = malloc(nn * sizeof(int));
    Expr* out = NULL;
    if (!t.cnt || !t.uncov || !t.upos || !t.excl || !t.S || !t.best || !t.mark || !gain
        || !t.odeg || !t.touched || !t.dcount || !t.ekey || !t.eord || !t.ey || !t.slack) goto done;
    for (int f = 0; f < m; f++) { t.uncov[f] = f; t.upos[f] = f; }
    t.nu = m;

    /* Upper bound: greedy max-coverage, with a lazy max-heap on gains
     * (a stale entry is re-pushed with its current gain), so the whole greedy
     * pass is O(sum|e| log n). gain[x] = uncovered hyperedges containing x. */
    {
        int64_t* heap = malloc((nn + (size_t)v.soff[m] + 1) * sizeof(int64_t));
        unsigned char* covered = calloc(mm, 1);
        if (!heap || !covered) { free(heap); free(covered); goto done; }
        size_t hn = 0;
        /* key = gain * 2^32 + (2^32 - 1 - x): max gain, ties to smallest x. */
#define HKEY(g, x) (((int64_t)(g) << 32) | (int64_t)(0xFFFFFFFFu - (uint32_t)(x)))
        for (int x = 0; x < n; x++) {
            gain[x] = v.voff[x + 1] - v.voff[x];
            if (!gain[x]) continue;
            int64_t key = HKEY(gain[x], x);
            size_t i = hn++;
            while (i > 0 && heap[(i - 1) / 2] < key) { heap[i] = heap[(i - 1) / 2]; i = (i - 1) / 2; }
            heap[i] = key;
        }
        while (t.nu > 0 && hn > 0) {
            int64_t top = heap[0];
            int64_t last = heap[--hn];
            size_t i = 0;
            for (;;) {                                       /* sift down */
                size_t c = 2 * i + 1;
                if (c >= hn) break;
                if (c + 1 < hn && heap[c + 1] > heap[c]) c++;
                if (heap[c] <= last) break;
                heap[i] = heap[c]; i = c;
            }
            if (hn > 0) heap[i] = last;
            int x = (int)(0xFFFFFFFFu - (uint32_t)(top & 0xFFFFFFFF));
            int g = (int)(top >> 32);
            if (g != gain[x]) {                              /* stale */
                if (gain[x] > 0) {
                    int64_t key = HKEY(gain[x], x);
                    size_t j = hn++;
                    while (j > 0 && heap[(j - 1) / 2] < key) { heap[j] = heap[(j - 1) / 2]; j = (j - 1) / 2; }
                    heap[j] = key;
                }
                continue;
            }
            mn_add(&t, x);
            for (int q = v.voff[x]; q < v.voff[x + 1]; q++) {
                int f = v.ve[q];
                if (covered[f]) continue;
                covered[f] = 1;
                for (int r = v.soff[f]; r < v.soff[f + 1]; r++) gain[v.sv[r]]--;
            }
        }
#undef HKEY
        free(heap); free(covered);
    }
    memcpy(t.best, t.S, (size_t)t.ns * sizeof(int));
    t.nbest = t.ns;
    while (t.ns > 0) mn_remove(&t, t.S[t.ns - 1]);

    /* A packing matching the greedy bound proves it optimal without search;
     * otherwise search, whose recursion depth is < nbest. */
    if (mn_lower_bound(&t) < t.nbest) {
        if (t.nbest > 20000) goto done;                /* C-stack depth guard */
        mn_search(&t);
        if (t.aborted) goto done;
    }
    qsort(t.best, (size_t)t.nbest, sizeof(int), cmp_int);
    {
        Expr** it = malloc((size_t)(t.nbest > 0 ? t.nbest : 1) * sizeof(Expr*));
        if (!it) goto done;
        for (int i = 0; i < t.nbest; i++) it[i] = expr_copy(v.verts->data.function.args[t.best[i]]);
        out = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)t.nbest);
        free(it);
    }
done:
    free(t.cnt); free(t.uncov); free(t.upos); free(t.excl); free(t.S);
    free(t.best); free(t.mark); free(gain);
    free(t.odeg); free(t.touched); free(t.dcount); free(t.ekey); free(t.eord); free(t.ey); free(t.slack);
    expr_free(owned);
    return out;
}
