/* gmet_core.c - shared kernels for the graph metrics builtins.
 *
 * Everything the distance / centrality / clustering heads have in common:
 *
 *   - CSR adjacency built straight from the validated-graph memo's endpoint
 *     arrays (graph_edge_indices), so no vertex hashing happens per call;
 *   - edge-weight resolution to doubles, with Wolfram's rule that a symbolic,
 *     complex or negative weight leaves the call unevaluated;
 *   - a small thread team (pthreads, MATHILDA_THREADS builds only) that hands
 *     out chunks of sources dynamically;
 *   - bit-parallel multi-source BFS (MS-BFS, Then et al., VLDB 2015): 256
 *     sources advance together, one bit per source per vertex, so each level
 *     walks the adjacency once for the whole batch. On small-world graphs this
 *     is ~10-30x less adjacency traffic than 256 separate BFS runs, and it is
 *     never worse than them;
 *   - lazy-deletion binary-heap Dijkstra for non-negative weights;
 *   - the per-source distance summary (reach / sum / eccentricity) that
 *     Closeness, EccentricityCentrality, VertexEccentricity and the
 *     diameter/radius/center/periphery/mean-distance heads all reduce from,
 *     cached per graph node so a sequence of those heads on one graph pays for
 *     one all-pairs pass;
 *   - packed result builders, and a per-graph result cache.
 *
 * Thread safety: worker functions touch only plain C arrays. Workers are
 * started with every signal blocked, and the parallel region is bracketed by
 * tc_async_defer_push/pop, so a TimeConstrained SIGPROF is taken on the main
 * thread after the region instead of unwinding across live workers.
 */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE 1   /* _SC_NPROCESSORS_ONLN under _POSIX_C_SOURCE */
#endif

#include "graph_metrics.h"
#include "graph.h"
#include "expr.h"
#include "pack.h"
#include "core.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef MATHILDA_THREADS
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#endif

/* ---- Graph kind & weights ------------------------------------------------- */

int gmet_graph_kind(const Expr* g) {
    long nd = graph_directed_edge_count(g);
    if (nd < 0) return -1;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    if (nd == 0) return GMET_KIND_UNDIRECTED;
    if ((size_t)nd == ne) return GMET_KIND_DIRECTED;
    return GMET_KIND_MIXED;
}

int gmet_edge_weights(const Expr* g, double** w_out) {
    *w_out = NULL;
    if (g->data.function.arg_count != 3) return 0;
    const Expr* wl = g->data.function.args[2]->data.function.args[1];
    size_t ne = wl->data.function.arg_count;
    double* w = malloc((ne ? ne : 1) * sizeof(double));
    if (!w) return -1;
    for (size_t k = 0; k < ne; k++) {
        const Expr* x = wl->data.function.args[k];
        double d;
        if (x->type == EXPR_INTEGER) d = (double)x->data.integer;
        else if (x->type == EXPR_REAL) d = x->data.real;
        else {
            d = graph_weight_to_double(x);   /* BigInt, MPFR, Rational; NAN else */
        }
        if (!(d >= 0.0) || isinf(d)) { free(w); return -1; }
        w[k] = d;
    }
    *w_out = w;
    return 1;
}

/* ---- CSR ------------------------------------------------------------------ */

void gmet_csr_free(GmetCSR* c) {
    if (!c) return;
    free(c->off); free(c->adj); free(c->w); free(c->eid);
    memset(c, 0, sizeof(*c));
}

int gmet_csr_build(const Expr* g, int mode, const double* ew, int want_eid,
                   GmetCSR* c) {
    memset(c, 0, sizeof(*c));
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return 0;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;

    c->n = n;
    c->off = calloc((size_t)n + 1, sizeof(int64_t));
    if (!c->off) return 0;
    int64_t na = 0;
    for (size_t k = 0; k < ne; k++) {
        int both = !edir[k] || mode == GMET_UND;
        if (both) { c->off[eu[k] + 1]++; c->off[ev[k] + 1]++; na += 2; }
        else if (mode == GMET_OUT) { c->off[eu[k] + 1]++; na++; }
        else { c->off[ev[k] + 1]++; na++; }
    }
    for (int i = 0; i < n; i++) c->off[i + 1] += c->off[i];
    c->narcs = na;
    size_t nb = (size_t)(na > 0 ? na : 1);
    c->adj = malloc(nb * sizeof(int));
    if (ew) c->w = malloc(nb * sizeof(double));
    if (want_eid) c->eid = malloc(nb * sizeof(int));
    int64_t* fill = malloc(((size_t)n + 1) * sizeof(int64_t));
    if (!c->adj || (ew && !c->w) || (want_eid && !c->eid) || !fill) {
        free(fill); gmet_csr_free(c); return 0;
    }
    memcpy(fill, c->off, ((size_t)n + 1) * sizeof(int64_t));
    for (size_t k = 0; k < ne; k++) {
        int a = eu[k], b = ev[k];
        int both = !edir[k] || mode == GMET_UND;
        int64_t p;
        if (both || mode == GMET_OUT) {
            p = fill[a]++; c->adj[p] = b;
            if (c->w) c->w[p] = ew[k];
            if (c->eid) c->eid[p] = (int)k;
        }
        if (both || mode == GMET_IN) {
            p = fill[b]++; c->adj[p] = a;
            if (c->w) c->w[p] = ew[k];
            if (c->eid) c->eid[p] = (int)k;
        }
    }
    free(fill);
    return 1;
}

/* ---- Thread team ---------------------------------------------------------- */

#define GMET_MAX_THREADS 32

int gmet_thread_count(double work) {
#ifdef MATHILDA_THREADS
    long cap = -1;
    const char* env = getenv("MATHILDA_GRAPH_THREADS");
    if (env && *env) cap = strtol(env, NULL, 10);
    if (cap == 0 || cap == 1) return 1;
    if (work < 4e6) return 1;
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpu < 1) ncpu = 1;
    if (cap > 1 && ncpu > cap) ncpu = cap;
    if (ncpu > GMET_MAX_THREADS) ncpu = GMET_MAX_THREADS;
    double by_work = work / 2e6;
    if (by_work < (double)ncpu) ncpu = (long)(by_work < 1 ? 1 : by_work);
    return (int)ncpu;
#else
    (void)work;
    return 1;
#endif
}

int gmet_cap_threads(int nt, double per_thread_bytes) {
    if (per_thread_bytes > 0) {
        double cap = GMET_SCRATCH_BUDGET / per_thread_bytes;
        if (cap < (double)nt) nt = cap < 1.0 ? 1 : (int)cap;
    }
    return nt < 1 ? 1 : nt;
}

#ifdef MATHILDA_THREADS
typedef struct {
    gmet_work_fn fn;
    void* ctx;
    int64_t n, chunk, next;
    pthread_mutex_t mu;
} GmetTeam;

typedef struct { GmetTeam* team; int tid; } GmetWorker;

static void team_loop(GmetTeam* t, int tid) {
    for (;;) {
        pthread_mutex_lock(&t->mu);
        int64_t lo = t->next;
        t->next = lo + t->chunk;
        pthread_mutex_unlock(&t->mu);
        if (lo >= t->n) break;
        int64_t hi = lo + t->chunk;
        if (hi > t->n) hi = t->n;
        t->fn(t->ctx, tid, lo, hi);
    }
}

static void* team_thread(void* p) {
    GmetWorker* w = (GmetWorker*)p;
    team_loop(w->team, w->tid);
    return NULL;
}
#endif

void gmet_parallel_for(int64_t n, int64_t chunk, int nthreads,
                       gmet_work_fn fn, void* ctx) {
    if (n <= 0) return;
    if (chunk < 1) chunk = 1;
#ifdef MATHILDA_THREADS
    if (nthreads > GMET_MAX_THREADS) nthreads = GMET_MAX_THREADS;
    if (nthreads > 1 && n > chunk) {
        GmetTeam t;
        t.fn = fn; t.ctx = ctx; t.n = n; t.chunk = chunk; t.next = 0;
        pthread_mutex_init(&t.mu, NULL);
        pthread_t th[GMET_MAX_THREADS];
        GmetWorker wk[GMET_MAX_THREADS];
        int started[GMET_MAX_THREADS];
        /* Workers inherit the creator's signal mask: block everything while
         * spawning so no asynchronous signal (SIGPROF from TimeConstrained,
         * SIGINT) is ever delivered on a worker. */
        sigset_t all, old;
        sigfillset(&all);
        pthread_sigmask(SIG_SETMASK, &all, &old);
        tc_async_defer_push();
        for (int i = 1; i < nthreads; i++) {
            wk[i].team = &t; wk[i].tid = i;
            started[i] = pthread_create(&th[i], NULL, team_thread, &wk[i]) == 0;
        }
        pthread_sigmask(SIG_SETMASK, &old, NULL);
        team_loop(&t, 0);
        for (int i = 1; i < nthreads; i++) if (started[i]) pthread_join(th[i], NULL);
        tc_async_defer_pop();
        pthread_mutex_destroy(&t.mu);
        return;
    }
#else
    (void)nthreads;
#endif
    for (int64_t lo = 0; lo < n; lo += chunk) {
        int64_t hi = lo + chunk;
        if (hi > n) hi = n;
        fn(ctx, 0, lo, hi);
    }
}

/* ---- Multi-source BFS ----------------------------------------------------- */

#define W GMET_MSBFS_WORDS

struct GmetMSBFS {
    int n;
    uint64_t* seen;      /* n * W: lanes that have reached v                  */
    uint64_t* front;     /* n * W: lanes whose frontier contains v            */
    uint64_t* next;      /* n * W: lanes arriving at v this level             */
    int* active;         /* vertices with a non-zero frontier                 */
    int* active2;
    int* touched;        /* vertices whose next[] was written this level      */
    unsigned char* intouch;
    unsigned char* done; /* seen[v] == every lane in use                      */
};

void gmet_msbfs_free(GmetMSBFS* s) {
    if (!s) return;
    free(s->seen); free(s->front); free(s->next);
    free(s->active); free(s->active2); free(s->touched);
    free(s->intouch); free(s->done);
    free(s);
}

GmetMSBFS* gmet_msbfs_new(int n) {
    GmetMSBFS* s = calloc(1, sizeof(GmetMSBFS));
    if (!s) return NULL;
    size_t nn = (size_t)(n > 0 ? n : 1);
    s->n = n;
    s->seen = calloc(nn * W, sizeof(uint64_t));
    s->front = calloc(nn * W, sizeof(uint64_t));
    s->next = calloc(nn * W, sizeof(uint64_t));
    s->active = malloc(nn * sizeof(int));
    s->active2 = malloc(nn * sizeof(int));
    s->touched = malloc(nn * sizeof(int));
    s->intouch = calloc(nn, 1);
    s->done = calloc(nn, 1);
    if (!s->seen || !s->front || !s->next || !s->active || !s->active2
        || !s->touched || !s->intouch || !s->done) {
        gmet_msbfs_free(s);
        return NULL;
    }
    return s;
}

void gmet_msbfs_run(GmetMSBFS* s, const GmetCSR* g, const int* sources,
                    int nsrc, int maxlevel, gmet_msbfs_visit visit, void* ctx) {
    int n = g->n;
    uint64_t full[W];
    memset(full, 0, sizeof(full));
    for (int i = 0; i < nsrc; i++) full[i >> 6] |= (uint64_t)1 << (i & 63);
    memset(s->seen, 0, (size_t)n * W * sizeof(uint64_t));
    memset(s->done, 0, (size_t)n);

    int nact = 0;
    for (int i = 0; i < nsrc; i++) {
        int v = sources[i];
        uint64_t bit = (uint64_t)1 << (i & 63);
        if (!(s->seen[(size_t)v * W] | s->seen[(size_t)v * W + 1]
              | s->seen[(size_t)v * W + 2] | s->seen[(size_t)v * W + 3]))
            s->active[nact++] = v;
        s->seen[(size_t)v * W + (i >> 6)] |= bit;
        s->front[(size_t)v * W + (i >> 6)] |= bit;
    }
    const int64_t* off = g->off;
    const int* adj = g->adj;
    for (int level = 1; nact > 0 && (maxlevel < 0 || level <= maxlevel); level++) {
        int nt = 0;
        for (int a = 0; a < nact; a++) {
            int v = s->active[a];
            const uint64_t* fv = s->front + (size_t)v * W;
            uint64_t f0 = fv[0], f1 = fv[1], f2 = fv[2], f3 = fv[3];
            for (int64_t p = off[v]; p < off[v + 1]; p++) {
                int w = adj[p];
                if (s->done[w]) continue;
                uint64_t* nw = s->next + (size_t)w * W;
                if (!s->intouch[w]) { s->intouch[w] = 1; s->touched[nt++] = w; }
                nw[0] |= f0; nw[1] |= f1; nw[2] |= f2; nw[3] |= f3;
            }
        }
        for (int a = 0; a < nact; a++) {
            uint64_t* fv = s->front + (size_t)s->active[a] * W;
            fv[0] = fv[1] = fv[2] = fv[3] = 0;
        }
        int nn2 = 0;
        for (int t = 0; t < nt; t++) {
            int w = s->touched[t];
            s->intouch[w] = 0;
            uint64_t* nw = s->next + (size_t)w * W;
            uint64_t* sw = s->seen + (size_t)w * W;
            uint64_t x0 = nw[0] & ~sw[0], x1 = nw[1] & ~sw[1];
            uint64_t x2 = nw[2] & ~sw[2], x3 = nw[3] & ~sw[3];
            nw[0] = nw[1] = nw[2] = nw[3] = 0;
            if (!(x0 | x1 | x2 | x3)) continue;
            sw[0] |= x0; sw[1] |= x1; sw[2] |= x2; sw[3] |= x3;
            if (sw[0] == full[0] && sw[1] == full[1] && sw[2] == full[2]
                && sw[3] == full[3]) s->done[w] = 1;
            uint64_t* fw = s->front + (size_t)w * W;
            fw[0] = x0; fw[1] = x1; fw[2] = x2; fw[3] = x3;
            s->active2[nn2++] = w;
            uint64_t lanes[W] = { x0, x1, x2, x3 };
            visit(ctx, level, w, lanes);
        }
        int* tmp = s->active; s->active = s->active2; s->active2 = tmp;
        nact = nn2;
    }
    /* Leave front[] zeroed for the next run (a maxlevel cut can strand bits). */
    for (int a = 0; a < nact; a++) {
        uint64_t* fv = s->front + (size_t)s->active[a] * W;
        fv[0] = fv[1] = fv[2] = fv[3] = 0;
    }
}

/* ---- Dijkstra ------------------------------------------------------------- */

struct GmetHeap { double* key; int* val; int64_t size, cap; };

GmetHeap* gmet_heap_new(int64_t cap) {
    GmetHeap* h = calloc(1, sizeof(GmetHeap));
    if (!h) return NULL;
    if (cap < 16) cap = 16;
    h->key = malloc((size_t)cap * sizeof(double));
    h->val = malloc((size_t)cap * sizeof(int));
    h->cap = cap;
    if (!h->key || !h->val) { gmet_heap_free(h); return NULL; }
    return h;
}

void gmet_heap_free(GmetHeap* h) {
    if (!h) return;
    free(h->key); free(h->val); free(h);
}

static void heap_push(GmetHeap* h, double k, int v) {
    if (h->size == h->cap) {
        int64_t nc = h->cap * 2;
        double* nk = realloc(h->key, (size_t)nc * sizeof(double));
        if (nk) h->key = nk;
        int* nv = realloc(h->val, (size_t)nc * sizeof(int));
        if (nv) h->val = nv;
        if (!nk || !nv) return;          /* allocation failure: drop (cannot happen in practice) */
        h->cap = nc;
    }
    int64_t i = h->size++;
    while (i > 0) {
        int64_t p = (i - 1) >> 1;
        if (h->key[p] <= k) break;
        h->key[i] = h->key[p]; h->val[i] = h->val[p];
        i = p;
    }
    h->key[i] = k; h->val[i] = v;
}

static void heap_pop(GmetHeap* h, double* k, int* v) {
    *k = h->key[0]; *v = h->val[0];
    int64_t n = --h->size;
    if (n == 0) return;
    double lk = h->key[n]; int lv = h->val[n];
    int64_t i = 0;
    for (;;) {
        int64_t c = 2 * i + 1;
        if (c >= n) break;
        if (c + 1 < n && h->key[c + 1] < h->key[c]) c++;
        if (h->key[c] >= lk) break;
        h->key[i] = h->key[c]; h->val[i] = h->val[c];
        i = c;
    }
    h->key[i] = lk; h->val[i] = lv;
}

void gmet_dijkstra(const GmetCSR* g, int src, double* dist, int* order,
                   int* nsettled, GmetHeap* h) {
    int n = g->n;
    for (int i = 0; i < n; i++) dist[i] = INFINITY;
    h->size = 0;
    dist[src] = 0.0;
    heap_push(h, 0.0, src);
    int ns = 0;
    while (h->size > 0) {
        double d; int u;
        heap_pop(h, &d, &u);
        if (d > dist[u]) continue;           /* stale entry */
        if (order) order[ns] = u;
        ns++;
        for (int64_t p = g->off[u]; p < g->off[u + 1]; p++) {
            int v = g->adj[p];
            double nd = d + g->w[p];
            if (nd < dist[v]) { dist[v] = nd; heap_push(h, nd, v); }
        }
    }
    if (nsettled) *nsettled = ns;
}

/* ---- Distance summary ----------------------------------------------------- */

typedef struct {
    const GmetCSR* g;
    GmetDistSummary* out;
    GmetMSBFS* bfs[GMET_MAX_THREADS];
    GmetHeap* heap[GMET_MAX_THREADS];
    double* dist[GMET_MAX_THREADS];
} SummaryCtx;

typedef struct {
    GmetDistSummary* out;
    int base;
    int64_t reach[GMET_MSBFS_LANES];
    double sum[GMET_MSBFS_LANES];
    int ecc[GMET_MSBFS_LANES];
} SummaryBatch;

static void summary_visit(void* ctx, int level, int v, const uint64_t* lanes) {
    (void)v;
    SummaryBatch* b = (SummaryBatch*)ctx;
    for (int k = 0; k < W; k++) {
        uint64_t x = lanes[k];
        while (x) {
            int i = (k << 6) + GMET_CTZ64(x);
            x &= x - 1;
            b->reach[i]++;
            b->sum[i] += level;
            b->ecc[i] = level;
        }
    }
}

static void summary_bfs_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    SummaryCtx* c = (SummaryCtx*)vctx;
    int n = c->g->n;
    SummaryBatch* b = malloc(sizeof(SummaryBatch));
    int* src = malloc(GMET_MSBFS_LANES * sizeof(int));
    if (!b || !src) { free(b); free(src); return; }
    for (int64_t blk = lo; blk < hi; blk++) {
        int base = (int)(blk * GMET_MSBFS_LANES);
        int ns = n - base < GMET_MSBFS_LANES ? n - base : GMET_MSBFS_LANES;
        for (int i = 0; i < ns; i++) src[i] = base + i;
        memset(b, 0, sizeof(*b));
        gmet_msbfs_run(c->bfs[tid], c->g, src, ns, -1, summary_visit, b);
        for (int i = 0; i < ns; i++) {
            c->out->reach[base + i] = b->reach[i];
            c->out->sum[base + i] = b->sum[i];
            c->out->ecc[base + i] = b->ecc[i];
        }
    }
    free(b); free(src);
}

static void summary_dijkstra_work(void* vctx, int tid, int64_t lo, int64_t hi) {
    SummaryCtx* c = (SummaryCtx*)vctx;
    int n = c->g->n;
    double* dist = c->dist[tid];
    for (int64_t s = lo; s < hi; s++) {
        gmet_dijkstra(c->g, (int)s, dist, NULL, NULL, c->heap[tid]);
        int64_t r = 0; double sum = 0.0, ecc = 0.0;
        for (int v = 0; v < n; v++) {
            if (v == s || isinf(dist[v])) continue;
            r++; sum += dist[v];
            if (dist[v] > ecc) ecc = dist[v];
        }
        c->out->reach[s] = r; c->out->sum[s] = sum; c->out->ecc[s] = ecc;
    }
}

static void summary_free(GmetDistSummary* s) {
    if (!s) return;
    free(s->reach); free(s->sum); free(s->ecc); free(s);
}

static GmetDistSummary* summary_compute(const GmetCSR* g) {
    int n = g->n;
    size_t nn = (size_t)(n > 0 ? n : 1);
    GmetDistSummary* out = calloc(1, sizeof(GmetDistSummary));
    if (!out) return NULL;
    out->n = n;
    out->weighted = g->w != NULL;
    out->reach = calloc(nn, sizeof(int64_t));
    out->sum = calloc(nn, sizeof(double));
    out->ecc = calloc(nn, sizeof(double));
    if (!out->reach || !out->sum || !out->ecc) { summary_free(out); return NULL; }
    if (n == 0) return out;

    SummaryCtx c;
    memset(&c, 0, sizeof(c));
    c.g = g; c.out = out;
    int ok = 1;
    if (!g->w) {
        int64_t nblk = (n + GMET_MSBFS_LANES - 1) / GMET_MSBFS_LANES;
        double work = (double)nblk * (double)(g->narcs + n) * 8.0;
        int nt = gmet_thread_count(work);
        if (nt > nblk) nt = (int)nblk;
        nt = gmet_cap_threads(nt, (double)n * (3 * 8 * GMET_MSBFS_WORDS + 14));
        for (int t = 0; t < nt; t++) if (!(c.bfs[t] = gmet_msbfs_new(n))) ok = 0;
        if (ok) gmet_parallel_for(nblk, 1, nt, summary_bfs_work, &c);
        for (int t = 0; t < nt; t++) gmet_msbfs_free(c.bfs[t]);
    } else {
        double work = (double)n * (double)(g->narcs * 4 + n);
        int nt = gmet_cap_threads(gmet_thread_count(work), (double)(g->narcs + n) * 12 + 8.0 * n);
        for (int t = 0; t < nt; t++) {
            c.heap[t] = gmet_heap_new(g->narcs + n + 1);
            c.dist[t] = malloc(nn * sizeof(double));
            if (!c.heap[t] || !c.dist[t]) ok = 0;
        }
        if (ok) gmet_parallel_for(n, 16, nt, summary_dijkstra_work, &c);
        for (int t = 0; t < nt; t++) { gmet_heap_free(c.heap[t]); free(c.dist[t]); }
    }
    if (!ok) { summary_free(out); return NULL; }
    return out;
}

/* Two-slot cache of summaries keyed on the graph node (a held reference keeps
 * the node alive, hence immutable -- see the memo comment in graph_util.c). */
#define SUMMARY_SLOTS 2
static Expr* g_sum_graph[SUMMARY_SLOTS];
static GmetDistSummary* g_sum_val[SUMMARY_SLOTS];
static int g_sum_next = 0;

static void cache_sweep(void);

const GmetDistSummary* gmet_dist_summary(const Expr* g, int* status) {
    *status = -1;
    cache_sweep();
    for (int i = 0; i < SUMMARY_SLOTS; i++)
        if (g_sum_graph[i] == g && g_sum_val[i]) { *status = 0; return g_sum_val[i]; }
    if (!graph_is_valid(g)) return NULL;
    double* ew = NULL;
    int wk = gmet_edge_weights(g, &ew);
    if (wk < 0) return NULL;
    GmetCSR csr;
    if (!gmet_csr_build(g, GMET_OUT, ew, 0, &csr)) { free(ew); return NULL; }
    free(ew);
    GmetDistSummary* s = summary_compute(&csr);
    gmet_csr_free(&csr);
    if (!s) return NULL;
    int slot = -1;
    for (int i = 0; i < SUMMARY_SLOTS && slot < 0; i++) if (!g_sum_graph[i]) slot = i;
    if (slot < 0) { slot = g_sum_next; g_sum_next = (g_sum_next + 1) % SUMMARY_SLOTS; }
    if (g_sum_graph[slot]) expr_free(g_sum_graph[slot]);
    summary_free(g_sum_val[slot]);
    g_sum_graph[slot] = expr_copy((Expr*)g);
    g_sum_val[slot] = s;
    *status = 0;
    return s;
}

/* ---- Strongly connected components (iterative Tarjan) --------------------- */

int gmet_scc(const GmetCSR* g, int* comp) {
    int n = g->n;
    if (n == 0) return 0;
    size_t nn = (size_t)n;
    int* index = malloc(nn * sizeof(int));
    int* low = malloc(nn * sizeof(int));
    int* stack = malloc(nn * sizeof(int));
    int* cs_v = malloc(nn * sizeof(int));        /* call stack: vertex         */
    int64_t* cs_p = malloc(nn * sizeof(int64_t));/* call stack: next arc       */
    unsigned char* on = calloc(nn, 1);
    if (!index || !low || !stack || !cs_v || !cs_p || !on) {
        free(index); free(low); free(stack); free(cs_v); free(cs_p); free(on);
        return -1;
    }
    for (int i = 0; i < n; i++) index[i] = -1;
    int idx = 0, sp = 0, ncomp = 0;
    for (int r = 0; r < n; r++) {
        if (index[r] >= 0) continue;
        int top = 0;
        cs_v[0] = r; cs_p[0] = g->off[r];
        index[r] = low[r] = idx++;
        stack[sp++] = r; on[r] = 1;
        while (top >= 0) {
            int v = cs_v[top];
            if (cs_p[top] < g->off[v + 1]) {
                int w = g->adj[cs_p[top]++];
                if (index[w] < 0) {
                    index[w] = low[w] = idx++;
                    stack[sp++] = w; on[w] = 1;
                    top++;
                    cs_v[top] = w; cs_p[top] = g->off[w];
                } else if (on[w] && index[w] < low[v]) {
                    low[v] = index[w];
                }
            } else {
                if (low[v] == index[v]) {
                    int w;
                    do { w = stack[--sp]; on[w] = 0; comp[w] = ncomp; } while (w != v);
                    ncomp++;
                }
                top--;
                if (top >= 0) {
                    int u = cs_v[top];
                    if (low[v] < low[u]) low[u] = low[v];
                }
            }
        }
    }
    free(index); free(low); free(stack); free(cs_v); free(cs_p); free(on);
    return ncomp;
}

/* ---- Result builders ------------------------------------------------------ */

Expr* gmet_real_vector(const double* v, int64_t n) {
    double* buf = NULL;
    Expr* nd = ndbuild_open_f64(n, &buf);
    if (nd) {
        if (n > 0) memcpy(buf, v, (size_t)n * sizeof(double));
        return nd;
    }
    Expr** items = n > 0 ? malloc((size_t)n * sizeof(Expr*)) : NULL;
    for (int64_t i = 0; i < n; i++) items[i] = expr_new_real(v[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}

Expr* gmet_int_vector(const int64_t* v, int64_t n) {
    int64_t* buf = NULL;
    Expr* nd = ndbuild_open_i64(n, &buf);
    if (nd) {
        if (n > 0) memcpy(buf, v, (size_t)n * sizeof(int64_t));
        return nd;
    }
    Expr** items = n > 0 ? malloc((size_t)n * sizeof(Expr*)) : NULL;
    for (int64_t i = 0; i < n; i++) items[i] = expr_new_integer(v[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}

Expr* gmet_vertex_subset(const Expr* g, const unsigned char* flag) {
    const Expr* verts = g->data.function.args[0];
    size_t n = verts->data.function.arg_count, k = 0;
    for (size_t i = 0; i < n; i++) k += flag[i] != 0;
    Expr** items = k ? malloc(k * sizeof(Expr*)) : NULL;
    k = 0;
    for (size_t i = 0; i < n; i++)
        if (flag[i]) items[k++] = expr_copy(verts->data.function.args[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, k);
    free(items);
    return out;
}

Expr* gmet_distance_value(double d, int weighted) {
    if (isinf(d)) return expr_new_symbol(SYM_Infinity);
    if (weighted) return expr_new_real(d);
    return expr_new_integer((int64_t)d);
}

/* ---- Per-graph result cache ----------------------------------------------- */

#define RESULT_SLOTS 16
typedef struct { const char* head; Expr* call; Expr* value; } ResultSlot;
static ResultSlot g_results[RESULT_SLOTS];
static int g_results_next = 0;

static int call_matches(const Expr* a, const Expr* b) {
    if (a == b) return 1;
    size_t na = a->data.function.arg_count;
    if (na != b->data.function.arg_count || na == 0) return 0;
    if (a->data.function.args[0] != b->data.function.args[0]) return 0;  /* same graph node */
    for (size_t i = 1; i < na; i++)
        if (!expr_eq(a->data.function.args[i], b->data.function.args[i])) return 0;
    return 1;
}

/* Drop entries whose graph nobody else references any more: the cache must
 * not keep a large graph alive after the user has let go of it (it would
 * otherwise pin up to RESULT_SLOTS graphs). A graph still held by the
 * validated-graph memo stays until the memo evicts it, then goes here. */
static void cache_sweep(void) {
    for (int i = 0; i < RESULT_SLOTS; i++) {
        ResultSlot* s = &g_results[i];
        if (!s->call) continue;
        const Expr* g = s->call->data.function.args[0];
        if (s->call->refcount == 1 && g->refcount == 1) {
            expr_free(s->call); expr_free(s->value);
            s->call = NULL; s->value = NULL; s->head = NULL;
        }
    }
    for (int i = 0; i < SUMMARY_SLOTS; i++) {
        if (g_sum_graph[i] && g_sum_graph[i]->refcount == 1) {
            expr_free(g_sum_graph[i]); g_sum_graph[i] = NULL;
            summary_free(g_sum_val[i]); g_sum_val[i] = NULL;
        }
    }
}

Expr* gmet_cache_get(const char* head, const Expr* res) {
    cache_sweep();
    for (int i = 0; i < RESULT_SLOTS; i++) {
        ResultSlot* s = &g_results[i];
        if (s->call && (s->head == head || strcmp(s->head, head) == 0)
            && call_matches(s->call, res))
            return expr_copy(s->value);
    }
    return NULL;
}

static size_t result_size(const Expr* e) {
    if (e->type == EXPR_NDARRAY) return 1;   /* one buffer, not a node per leaf */
    if (e->type != EXPR_FUNCTION) return 1;
    size_t k = 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* a = e->data.function.args[i];
        k += (a->type == EXPR_FUNCTION) ? 1 + a->data.function.arg_count : 1;
    }
    return k;
}

void gmet_cache_put(const char* head, const Expr* res, Expr* value) {
    if (!value) return;
    size_t sz = result_size(value);
    if (value->type == EXPR_NDARRAY) {
        size_t elems = 1;
        for (int r = 0; r < value->data.ndarray.rank; r++)
            elems *= (size_t)value->data.ndarray.dims[r];
        sz = elems / 8 + 1;                  /* buffers are ~8x denser than nodes */
    }
    if (sz > GMET_CACHE_MAX_ELEMS / 4) return;
    cache_sweep();
    ResultSlot* s = NULL;
    for (int i = 0; i < RESULT_SLOTS && !s; i++) if (!g_results[i].call) s = &g_results[i];
    if (!s) s = &g_results[g_results_next];
    g_results_next = (g_results_next + 1) % RESULT_SLOTS;
    if (s->call) expr_free(s->call);
    if (s->value) expr_free(s->value);
    s->head = head;
    s->call = expr_copy((Expr*)res);
    s->value = expr_copy(value);
}
