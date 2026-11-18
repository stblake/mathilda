/*
 * fi_gf2.c — relations, GF(2) linear algebra and the square-root step shared
 * by the congruence-of-squares methods (CFRAC, Dixon, the quadratic sieve).
 *
 * A relation is X^2 == sign * prod p_i^e_i * E^2 (mod n) over a factor base,
 * where E collects squared large primes from combined partial relations.
 * Partial relations carry one or two large primes outside the factor base
 * and are kept as edges of a graph on the large primes (vertex 0 stands for
 * "no large prime", so a one-prime partial is an edge 0 - L).  Union-find
 * counts the independent cycles as edges arrive; every cycle multiplies out
 * to a full relation, since each of its vertices' primes appears exactly
 * twice (E = the product of the cycle's primes).  Two one-prime partials
 * with the same L are the 2-cycle 0 - L - 0: the classical pairing.
 *
 * Linear algebra finds subsets of relations whose exponent vectors sum to
 * zero mod 2.  The matrix is first pruned (singleton columns removed
 * repeatedly, since a relation holding the only odd occurrence of a prime
 * can never be in a dependency), then the survivors are solved:
 *   - dense Gaussian elimination on bit-packed columns below DENSE_MAX
 *     columns: simple, cache friendly, cost ~ c^3/64 word ops;
 *   - block Lanczos (Montgomery 1995, 64-bit blocks) above, O(w c^2 / 64)
 *     for w nonzeros, the only practical solver for the 10^4..10^5-column
 *     matrices of a 70..100 digit quadratic sieve.
 * Each dependency gives X = prod X_r, Y = sign-free sqrt of the product of
 * the right-hand sides (exponents halved) times prod E_r, and gcd(X - Y, n)
 * splits n with probability 1/2 per prime pair.
 *
 * Relation logs (FiRelLog) carry relations from worker threads to the main
 * thread, which alone touches the store: QS replays them per worker, CFRAC
 * and Dixon merge them by candidate index to keep the serial order.
 */
#include "fi_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Relation store                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t     X, E;
    uint32_t  off, cnt;      /* into the shared (col, exp) pool */
} Rel;

typedef struct {
    mpz_t     X;
    uint32_t  u, v;          /* vertices (large primes) */
    uint32_t  off, cnt;
    uint8_t   tree;          /* spanning-forest edge (else it closed a cycle) */
    uint8_t   done;          /* cycle already turned into a full relation */
} Edge;

struct FiRels {
    size_t    ncols;         /* factor-base columns, column 0 = sign */
    const uint32_t* fb;      /* fb[c] = prime for column c (fb[0] unused) */
    Rel*      rel;
    size_t    nrel, caprel;
    uint32_t* pcol;          /* pool of column indices */
    uint32_t* pexp;          /* matching exponents */
    size_t    npool, cappool;
    Edge*     edge;          /* partial relations */
    size_t    nedge, capedge;
    uint64_t* vkey;          /* open-addressed map large prime -> vertex */
    uint32_t* vval;
    size_t    vcap;
    uint64_t* vprime;        /* vertex -> large prime (vprime[0] = 1) */
    uint32_t* uf;            /* union-find parent */
    size_t    nvert, capvert;
    size_t    pending;       /* cycles not yet materialised */
    uint32_t* ex;            /* scratch: exponent per column */
    unsigned char* colseen;  /* column ever odd in a full relation */
    size_t    nactive;       /* number of such columns */
};

static uint32_t pool_push(FiRels* R, const uint32_t* col, const uint32_t* exps, size_t cnt) {
    if (R->npool + cnt > R->cappool) {
        while (R->npool + cnt > R->cappool) R->cappool = R->cappool ? 2 * R->cappool : 4096;
        R->pcol = realloc(R->pcol, R->cappool * sizeof(uint32_t));
        R->pexp = realloc(R->pexp, R->cappool * sizeof(uint32_t));
    }
    uint32_t off = (uint32_t)R->npool;
    memcpy(R->pcol + off, col, cnt * sizeof(uint32_t));
    memcpy(R->pexp + off, exps, cnt * sizeof(uint32_t));
    R->npool += cnt;
    return off;
}

FiRels* fi_rels_new(const uint32_t* fb, size_t ncols) {
    FiRels* R = calloc(1, sizeof *R);
    R->fb = fb;
    R->ncols = ncols;
    R->vcap = 1 << 12;
    R->vkey = calloc(R->vcap, sizeof(uint64_t));
    R->vval = calloc(R->vcap, sizeof(uint32_t));
    R->capvert = 1024;
    R->vprime = malloc(R->capvert * sizeof(uint64_t));
    R->uf = malloc(R->capvert * sizeof(uint32_t));
    R->vprime[0] = 1;
    R->uf[0] = 0;
    R->nvert = 1;
    R->ex = calloc(ncols + 1, sizeof(uint32_t));
    R->colseen = calloc(ncols + 1, 1);
    return R;
}

void fi_rels_free(FiRels* R) {
    if (!R) return;
    for (size_t i = 0; i < R->nrel; i++) mpz_clears(R->rel[i].X, R->rel[i].E, NULL);
    for (size_t i = 0; i < R->nedge; i++) mpz_clear(R->edge[i].X);
    free(R->rel); free(R->pcol); free(R->pexp); free(R->edge);
    free(R->vkey); free(R->vval); free(R->vprime); free(R->uf);
    free(R->ex); free(R->colseen);
    free(R);
}

size_t fi_rels_count(const FiRels* R) { return R->nrel + R->pending; }
size_t fi_rels_active(const FiRels* R) { return R->nactive; }
size_t fi_rels_partials(const FiRels* R) { return R->nedge; }

static void add_full(FiRels* R, const mpz_t X, const mpz_t E, const uint32_t* col,
                     const uint32_t* exps, size_t cnt) {
    if (R->nrel == R->caprel) {
        R->caprel = R->caprel ? 2 * R->caprel : 1024;
        R->rel = realloc(R->rel, R->caprel * sizeof(Rel));
    }
    Rel* r = &R->rel[R->nrel++];
    mpz_init_set(r->X, X);
    mpz_init_set(r->E, E);
    r->cnt = (uint32_t)cnt;
    r->off = pool_push(R, col, exps, cnt);
    for (size_t k = 0; k < cnt; k++)
        if ((exps[k] & 1) && !R->colseen[col[k]]) { R->colseen[col[k]] = 1; R->nactive++; }
}

/* Vertex id of large prime L (1 -> vertex 0), creating it if new. */
static uint32_t vertex(FiRels* R, uint64_t L) {
    if (L <= 1) return 0;
    if (2 * R->nvert >= R->vcap) {                      /* grow the map */
        size_t oc = R->vcap;
        uint64_t* ok = R->vkey;
        uint32_t* ov = R->vval;
        R->vcap *= 2;
        R->vkey = calloc(R->vcap, sizeof(uint64_t));
        R->vval = calloc(R->vcap, sizeof(uint32_t));
        for (size_t i = 0; i < oc; i++) {
            if (!ok[i]) continue;
            size_t h = (size_t)((ok[i] * UINT64_C(0x9E3779B97F4A7C15)) >> 20) & (R->vcap - 1);
            while (R->vkey[h]) h = (h + 1) & (R->vcap - 1);
            R->vkey[h] = ok[i];
            R->vval[h] = ov[i];
        }
        free(ok); free(ov);
    }
    size_t h = (size_t)((L * UINT64_C(0x9E3779B97F4A7C15)) >> 20) & (R->vcap - 1);
    while (R->vkey[h] && R->vkey[h] != L) h = (h + 1) & (R->vcap - 1);
    if (R->vkey[h]) return R->vval[h];
    if (R->nvert == R->capvert) {
        R->capvert *= 2;
        R->vprime = realloc(R->vprime, R->capvert * sizeof(uint64_t));
        R->uf = realloc(R->uf, R->capvert * sizeof(uint32_t));
    }
    uint32_t id = (uint32_t)R->nvert++;
    R->vprime[id] = L;
    R->uf[id] = id;
    R->vkey[h] = L;
    R->vval[h] = id;
    return id;
}

static uint32_t uf_find(uint32_t* uf, uint32_t x) {
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
}

/* Add a relation X^2 == prod fb[col]^exp * L1 * L2 (mod n).  col must be
 * sorted ascending (column 0 is the sign, -1).  L1, L2 are 1 or primes
 * outside the factor base.  Returns 1 if the count of usable relations
 * (fulls plus cycles) grew. */
int fi_rels_add2(FiRels* R, const mpz_t X, const uint32_t* col, const uint32_t* exps,
                 size_t cnt, uint64_t L1, uint64_t L2, const mpz_t n) {
    (void)n;
    if (L1 > L2) { uint64_t t = L1; L1 = L2; L2 = t; }
    if (L2 <= 1 || L1 == L2) {
        /* Smooth, or the cofactor was a square L^2: a full relation. */
        mpz_t E;
        mth_mpz_init_set_u64(E, (uint64_t)(L2 <= 1 ? 1 : L2));
        add_full(R, X, E, col, exps, cnt);
        mpz_clear(E);
        return 1;
    }
    uint32_t u = vertex(R, L1), v = vertex(R, L2);
    if (R->nedge == R->capedge) {
        R->capedge = R->capedge ? 2 * R->capedge : 4096;
        R->edge = realloc(R->edge, R->capedge * sizeof(Edge));
    }
    Edge* e = &R->edge[R->nedge++];
    mpz_init_set(e->X, X);
    e->u = u; e->v = v;
    e->cnt = (uint32_t)cnt;
    e->off = pool_push(R, col, exps, cnt);
    e->done = 0;
    uint32_t ru = uf_find(R->uf, u), rv = uf_find(R->uf, v);
    if (ru != rv) { R->uf[ru] = rv; e->tree = 1; return 0; }
    e->tree = 0;
    R->pending++;
    return 1;
}

int fi_rels_add(FiRels* R, const mpz_t X, const uint32_t* col, const uint32_t* exps,
                size_t cnt, uint64_t L, const mpz_t n) {
    return fi_rels_add2(R, X, col, exps, cnt, 1, L, n);
}

/* ------------------------------------------------------------------ */
/* Relation logs (worker threads -> main thread)                       */
/* ------------------------------------------------------------------ */

/* A log's records keep their mpz_t once initialised, so a worker that fills
 * the same log region after region allocates nothing in steady state. */
void fi_rellog_push(FiRelLog* lg, uint64_t tag, const mpz_t X, const uint32_t* col,
                    const uint32_t* exps, size_t cnt, uint64_t L1, uint64_t L2) {
    if (lg->n == lg->cap) {
        lg->cap = lg->cap ? 2 * lg->cap : 256;
        lg->r = realloc(lg->r, lg->cap * sizeof(FiLogRec));
    }
    if (lg->npool + cnt > lg->cappool) {
        while (lg->npool + cnt > lg->cappool) lg->cappool = lg->cappool ? 2 * lg->cappool : 4096;
        lg->col = realloc(lg->col, lg->cappool * sizeof(uint32_t));
        lg->exps = realloc(lg->exps, lg->cappool * sizeof(uint32_t));
    }
    FiLogRec* r = &lg->r[lg->n];
    if (lg->n == lg->ninit) { mpz_init_set(r->X, X); lg->ninit++; }
    else mpz_set(r->X, X);
    lg->n++;
    r->tag = tag;
    r->off = (uint32_t)lg->npool;
    r->cnt = (uint32_t)cnt;
    r->L1 = L1;
    r->L2 = L2;
    memcpy(lg->col + lg->npool, col, cnt * sizeof(uint32_t));
    memcpy(lg->exps + lg->npool, exps, cnt * sizeof(uint32_t));
    lg->npool += cnt;
}

void fi_rellog_replay(FiRelLog* lg, FiRels* R, const mpz_t n) {
    for (size_t i = 0; i < lg->n; i++) {
        FiLogRec* r = &lg->r[i];
        fi_rels_add2(R, r->X, lg->col + r->off, lg->exps + r->off, r->cnt, r->L1, r->L2, n);
    }
    lg->n = 0;
    lg->npool = 0;
}

int fi_rellog_merge(FiRelLog* logs, int nlogs,
                    int (*visit)(void* user, const FiLogRec* r, const uint32_t* col,
                                 const uint32_t* exps),
                    void* user) {
    size_t pos[FI_MAXTHREADS] = { 0 };
    int out = 0;
    if (nlogs > FI_MAXTHREADS) nlogs = FI_MAXTHREADS;
    /* At most FI_MAXTHREADS logs: a linear scan for the minimum is cheaper
     * than a heap, and the records are few (smooth candidates only). */
    for (;;) {
        int best = -1;
        for (int w = 0; w < nlogs; w++)
            if (pos[w] < logs[w].n && (best < 0 || logs[w].r[pos[w]].tag < logs[best].r[pos[best]].tag))
                best = w;
        if (best < 0) break;
        FiRelLog* lg = &logs[best];
        const FiLogRec* r = &lg->r[pos[best]++];
        out = visit(user, r, lg->col + r->off, lg->exps + r->off);
        if (out) break;
    }
    for (int w = 0; w < nlogs; w++) { logs[w].n = 0; logs[w].npool = 0; }
    return out;
}

void fi_rellog_free(FiRelLog* lg) {
    for (size_t i = 0; i < lg->ninit; i++) mpz_clear(lg->r[i].X);
    free(lg->r); free(lg->col); free(lg->exps);
    lg->r = NULL; lg->col = lg->exps = NULL;
    lg->n = lg->ninit = lg->cap = lg->npool = lg->cappool = 0;
}

/* Turn every pending cycle into a full relation.  The tree edges form a
 * spanning forest; a closing edge (u, v) plus the forest path v -> u is the
 * cycle. */
static void materialise_cycles(FiRels* R, const mpz_t n) {
    if (!R->pending) return;
    size_t V = R->nvert, E = R->nedge;
    uint32_t* head = malloc(V * sizeof(uint32_t));
    uint32_t* nxt = malloc(2 * E * sizeof(uint32_t) + 1);
    uint32_t* to = malloc(2 * E * sizeof(uint32_t) + 1);
    uint32_t* eid = malloc(2 * E * sizeof(uint32_t) + 1);
    for (size_t i = 0; i < V; i++) head[i] = UINT32_MAX;
    size_t na = 0;
    for (size_t k = 0; k < E; k++) {
        if (!R->edge[k].tree) continue;
        uint32_t a = R->edge[k].u, b = R->edge[k].v;
        to[na] = b; eid[na] = (uint32_t)k; nxt[na] = head[a]; head[a] = (uint32_t)na++;
        to[na] = a; eid[na] = (uint32_t)k; nxt[na] = head[b]; head[b] = (uint32_t)na++;
    }
    uint32_t* par = malloc(V * sizeof(uint32_t));      /* parent vertex */
    uint32_t* pe = malloc(V * sizeof(uint32_t));       /* edge to parent */
    uint32_t* dep = malloc(V * sizeof(uint32_t));
    uint32_t* queue = malloc(V * sizeof(uint32_t));
    for (size_t i = 0; i < V; i++) dep[i] = UINT32_MAX;
    for (size_t r = 0; r < V; r++) {
        if (dep[r] != UINT32_MAX) continue;
        size_t qh = 0, qt = 0;
        queue[qt++] = (uint32_t)r;
        dep[r] = 0; par[r] = (uint32_t)r; pe[r] = UINT32_MAX;
        while (qh < qt) {
            uint32_t x = queue[qh++];
            for (uint32_t a = head[x]; a != UINT32_MAX; a = nxt[a]) {
                uint32_t y = to[a];
                if (dep[y] != UINT32_MAX) continue;
                dep[y] = dep[x] + 1; par[y] = x; pe[y] = eid[a];
                queue[qt++] = y;
            }
        }
    }
    size_t ncols = R->ncols;
    uint32_t* touched = malloc((ncols + 1) * sizeof(uint32_t));
    uint32_t* oc = malloc((ncols + 1) * sizeof(uint32_t));
    uint32_t* oe = malloc((ncols + 1) * sizeof(uint32_t));
    mpz_t X, Ev;
    mpz_inits(X, Ev, NULL);
    for (size_t k = 0; k < E; k++) {
        Edge* c = &R->edge[k];
        if (c->tree || c->done) continue;
        c->done = 1;
        size_t nt = 0;
        mpz_set(X, c->X);
        mpz_set_ui(Ev, 1);
#define ACCUM(edge) do { \
            const Edge* e_ = (edge); \
            for (uint32_t q = 0; q < e_->cnt; q++) { \
                uint32_t cc = R->pcol[e_->off + q]; \
                if (!R->ex[cc]) touched[nt++] = cc; \
                R->ex[cc] += R->pexp[e_->off + q]; \
            } \
        } while (0)
        ACCUM(c);
        /* Walk u and v up to their common ancestor. */
        uint32_t a = c->u, b = c->v;
        while (a != b) {
            if (dep[a] >= dep[b]) {
                const Edge* t = &R->edge[pe[a]];
                ACCUM(t);
                mpz_mul(X, X, t->X); mpz_mod(X, X, n);
                mth_mpz_mul_u64(Ev, Ev, (uint64_t)R->vprime[a]); mpz_mod(Ev, Ev, n);
                a = par[a];
            } else {
                const Edge* t = &R->edge[pe[b]];
                ACCUM(t);
                mpz_mul(X, X, t->X); mpz_mod(X, X, n);
                mth_mpz_mul_u64(Ev, Ev, (uint64_t)R->vprime[b]); mpz_mod(Ev, Ev, n);
                b = par[b];
            }
        }
        mth_mpz_mul_u64(Ev, Ev, (uint64_t)R->vprime[a]); mpz_mod(Ev, Ev, n);   /* the apex */
#undef ACCUM
        /* Collect the summed exponents in column order. */
        size_t m = 0;
        for (size_t t = 0; t < nt; t++) { oc[m++] = touched[t]; }
        for (size_t t = 1; t < m; t++) {                   /* insertion sort */
            uint32_t x = oc[t];
            size_t y = t;
            while (y > 0 && oc[y - 1] > x) { oc[y] = oc[y - 1]; y--; }
            oc[y] = x;
        }
        for (size_t t = 0; t < m; t++) { oe[t] = R->ex[oc[t]]; R->ex[oc[t]] = 0; }
        add_full(R, X, Ev, oc, oe, m);
    }
    R->pending = 0;
    mpz_clears(X, Ev, NULL);
    free(touched); free(oc); free(oe);
    free(head); free(nxt); free(to); free(eid); free(par); free(pe); free(dep); free(queue);
}

/* ------------------------------------------------------------------ */
/* Dense Gaussian elimination                                          */
/* ------------------------------------------------------------------ */

/* Columns of the GF(2) system are the relations (nr of them); rows are the
 * factor-base primes (nc).  Returns dependencies as row-bit vectors over
 * relations. */
static size_t dense_deps(uint64_t* rows, size_t nc, size_t nr, uint64_t* deps, size_t max_deps) {
    size_t W = fi_gf2_words(nr);
    uint32_t* pivcol = malloc(nc * sizeof(uint32_t));
    unsigned char* is_piv = calloc(nr, 1);
    size_t rank = 0;
    for (size_t j = 0; j < nr && rank < nc; j++) {
        size_t w = j / 64;
        uint64_t bit = UINT64_C(1) << (j % 64);
        size_t r = rank;
        while (r < nc && !(rows[r * W + w] & bit)) r++;
        if (r == nc) continue;
        if (r != rank)
            for (size_t k = 0; k < W; k++) {
                uint64_t t = rows[r * W + k]; rows[r * W + k] = rows[rank * W + k]; rows[rank * W + k] = t;
            }
        uint64_t* pr = rows + rank * W;
        for (size_t i = 0; i < nc; i++) {
            if (i == rank) continue;
            uint64_t* ri = rows + i * W;
            if (ri[w] & bit) for (size_t k = w; k < W; k++) ri[k] ^= pr[k];
        }
        pivcol[rank++] = (uint32_t)j;
        is_piv[j] = 1;
    }
    size_t nd = 0;
    for (size_t f = 0; f < nr && nd < max_deps; f++) {
        if (is_piv[f]) continue;
        uint64_t* x = deps + nd * W;
        memset(x, 0, W * sizeof(uint64_t));
        x[f / 64] |= UINT64_C(1) << (f % 64);
        for (size_t i = 0; i < rank; i++)
            if ((rows[i * W + f / 64] >> (f % 64)) & 1)
                x[pivcol[i] / 64] |= UINT64_C(1) << (pivcol[i] % 64);
        nd++;
    }
    free(pivcol);
    free(is_piv);
    return nd;
}

/* ------------------------------------------------------------------ */
/* Sparse front end: prune singletons, dispatch to a solver            */
/* ------------------------------------------------------------------ */

size_t fi_gf2_lanczos(const FiSparseMat* m, uint64_t* deps, size_t max_deps);

/* Lanczos needs room for 64-wide blocks; tiny systems go dense. */
#define DENSE_MAX 1000

size_t fi_gf2_dependencies(const FiSparseMat* m, uint64_t* deps, size_t max_deps) {
    size_t nr = m->nrows, nc = m->ncols;
    /* Prune: drop rows containing a column of weight 1, repeatedly. */
    unsigned char* alive = malloc(nr);
    memset(alive, 1, nr);
    uint32_t* w = calloc(nc, sizeof(uint32_t));
    for (size_t r = 0; r < nr; r++)
        for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) w[m->col[k]]++;
    for (int changed = 1; changed; ) {
        changed = 0;
        for (size_t r = 0; r < nr; r++) {
            if (!alive[r]) continue;
            for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) {
                if (w[m->col[k]] == 1) {
                    alive[r] = 0;
                    for (uint32_t q = m->start[r]; q < m->start[r + 1]; q++) w[m->col[q]]--;
                    changed = 1;
                    break;
                }
            }
        }
    }
    /* Compact columns and rows. */
    uint32_t* cmap = malloc(nc * sizeof(uint32_t));
    size_t nc2 = 0;
    for (size_t c = 0; c < nc; c++) cmap[c] = w[c] ? (uint32_t)nc2++ : UINT32_MAX;
    uint32_t* rmap = malloc(nr * sizeof(uint32_t));
    size_t nr2 = 0;
    for (size_t r = 0; r < nr; r++) if (alive[r]) rmap[nr2++] = (uint32_t)r;
    size_t nd = 0;
    if (nr2 > nc2 || (nr2 > 0 && nc2 == 0)) {
        size_t W2 = fi_gf2_words(nr2);
        uint64_t* d2 = calloc(max_deps * W2 + 1, sizeof(uint64_t));
        size_t dense_max = DENSE_MAX;
        FI_TUNE_U("FI_DENSE_MAX", dense_max);
        if (nc2 <= dense_max) {
            uint64_t* rows = calloc(nc2 * W2 + 1, sizeof(uint64_t));
            for (size_t i = 0; i < nr2; i++) {
                size_t r = rmap[i];
                for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) {
                    uint32_t c = cmap[m->col[k]];
                    rows[c * W2 + i / 64] ^= UINT64_C(1) << (i % 64);
                }
            }
            nd = dense_deps(rows, nc2, nr2, d2, max_deps);
            free(rows);
        } else {
            FiSparseMat sm;
            sm.nrows = (uint32_t)nr2;
            sm.ncols = (uint32_t)nc2;
            sm.start = malloc((nr2 + 1) * sizeof(uint32_t));
            size_t tot = 0;
            for (size_t i = 0; i < nr2; i++) tot += m->start[rmap[i] + 1] - m->start[rmap[i]];
            sm.col = malloc((tot + 1) * sizeof(uint32_t));
            size_t q = 0;
            for (size_t i = 0; i < nr2; i++) {
                sm.start[i] = (uint32_t)q;
                size_t r = rmap[i];
                for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) sm.col[q++] = cmap[m->col[k]];
            }
            sm.start[nr2] = (uint32_t)q;
            nd = fi_gf2_lanczos(&sm, d2, max_deps);
            free(sm.start);
            free(sm.col);
            if (nd == 0 && nc2 <= 4 * DENSE_MAX) {
                /* Lanczos can fail on small or degenerate systems: fall back. */
                uint64_t* rows = calloc(nc2 * W2 + 1, sizeof(uint64_t));
                for (size_t i = 0; i < nr2; i++) {
                    size_t r = rmap[i];
                    for (uint32_t k = m->start[r]; k < m->start[r + 1]; k++) {
                        uint32_t c = cmap[m->col[k]];
                        rows[c * W2 + i / 64] ^= UINT64_C(1) << (i % 64);
                    }
                }
                nd = dense_deps(rows, nc2, nr2, d2, max_deps);
                free(rows);
            }
        }
        /* Map dependencies back to original row numbering. */
        size_t W = fi_gf2_words(nr);
        for (size_t k = 0; k < nd; k++) {
            uint64_t* out = deps + k * W;
            memset(out, 0, W * sizeof(uint64_t));
            for (size_t i = 0; i < nr2; i++)
                if ((d2[k * W2 + i / 64] >> (i % 64)) & 1)
                    out[rmap[i] / 64] |= UINT64_C(1) << (rmap[i] % 64);
        }
        free(d2);
    }
    free(alive); free(w); free(cmap); free(rmap);
    return nd;
}

/* ------------------------------------------------------------------ */
/* Square root step                                                    */
/* ------------------------------------------------------------------ */

/* Try every dependency until one splits n.  Returns 1 with d set. */
int fi_rels_solve(FiRels* R, const mpz_t n, mpz_t d) {
    materialise_cycles(R, n);
    size_t nr = R->nrel;
    if (nr == 0) return 0;
    FiSparseMat m;
    m.nrows = (uint32_t)nr;
    m.ncols = (uint32_t)R->ncols;
    m.start = malloc((nr + 1) * sizeof(uint32_t));
    m.col = malloc((R->npool + 1) * sizeof(uint32_t));
    size_t q = 0;
    for (size_t r = 0; r < nr; r++) {
        m.start[r] = (uint32_t)q;
        const Rel* x = &R->rel[r];
        for (uint32_t k = 0; k < x->cnt; k++)
            if (R->pexp[x->off + k] & 1) m.col[q++] = R->pcol[x->off + k];
    }
    m.start[nr] = (uint32_t)q;
    const size_t MAXD = 64;
    size_t W = fi_gf2_words(nr);
    uint64_t* deps = calloc(MAXD * W, sizeof(uint64_t));
    size_t nd = fi_gf2_dependencies(&m, deps, MAXD);
    free(m.start);
    free(m.col);

    int found = 0;
    uint64_t* ex = calloc(R->ncols, sizeof(uint64_t));
    mpz_t X, Y, t;
    mpz_inits(X, Y, t, NULL);
    for (size_t k = 0; k < nd && !found; k++) {
        memset(ex, 0, R->ncols * sizeof(uint64_t));
        mpz_set_ui(X, 1);
        mpz_set_ui(Y, 1);
        for (size_t r = 0; r < nr; r++) {
            if (!((deps[k * W + r / 64] >> (r % 64)) & 1)) continue;
            const Rel* x = &R->rel[r];
            mpz_mul(X, X, x->X); mpz_mod(X, X, n);
            mpz_mul(Y, Y, x->E); mpz_mod(Y, Y, n);
            for (uint32_t j = 0; j < x->cnt; j++) ex[R->pcol[x->off + j]] += R->pexp[x->off + j];
        }
        int bad = 0;
        for (size_t c = 1; c < R->ncols; c++) {
            if (ex[c] & 1) { bad = 1; break; }
            if (ex[c]) {
                mpz_set_ui(t, R->fb[c]);
                mpz_powm_ui(t, t, (unsigned long)(ex[c] / 2), n);
                mpz_mul(Y, Y, t); mpz_mod(Y, Y, n);
            }
        }
        if (bad || (ex[0] & 1)) continue;
        mpz_sub(t, X, Y);
        mpz_gcd(t, t, n);
        if (mpz_cmp_ui(t, 1) > 0 && mpz_cmp(t, n) < 0) { mpz_set(d, t); found = 1; }
    }
    mpz_clears(X, Y, t, NULL);
    free(ex);
    free(deps);
    return found;
}
