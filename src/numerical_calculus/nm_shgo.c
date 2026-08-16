/* nm_shgo.c — NMinimize SHGO (Simplicial Homology Global Optimization) engine.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ============================================================================
 *  SHGO — Simplicial Homology Global Optimization (Endres, Sandrock & Focke,
 *  J. Glob. Optim. 72 (2018) 181-217), mirroring scipy.optimize.shgo.
 *
 *  The method samples the bounded search box, builds a graph on the samples, and
 *  identifies the "minimizer pool": vertices that are strictly better (by Deb's
 *  feasibility rules) than every graph-neighbour. This is the discrete
 *  sub-level-set homology construction — each pool vertex sits in a distinct
 *  locally-convex sub-domain — so a single local minimization started from each
 *  pool vertex reaches every basin with far fewer local searches than blind
 *  multi-start. The number of distinct local minima (the homology group order)
 *  gives a principled stopping criterion under refinement.
 *
 *  Sampling (the "SamplingMethod" sub-option):
 *    "Simplicial" (default) — the exact 1-skeleton of the Freudenthal/Kuhn
 *        triangulation of the box, refined by longest-edge bisection. Convergence
 *        guaranteed; practical for low dimension (n <= SHGO_SIMPLICIAL_MAXDIM,
 *        else it falls back to Sobol, matching scipy's practical guidance).
 *    "Sobol" / "Halton" — a low-discrepancy point set whose graph connectivity is
 *        an adaptive k-nearest-neighbour directed graph. This kNN graph is a
 *        deliberate, documented approximation of scipy's Delaunay-based pool: the
 *        pool it yields seeds the same local searches and the returned optimum is
 *        validated against scipy (benchmarks/79-shgo), only the internal graph
 *        differs.
 *
 *  Everything downstream is shared with the other NMinimize engines: nm_eval
 *  (compiled-objective fast path + constraint penalty), nm_better (Deb selection)
 *  and nm_local_polish (BFGS / augmented-Lagrangian / integer descent). Only the
 *  bounded box (D->reg_lo/reg_hi) and those helpers are used, so constrained and
 *  mixed-integer problems flow through unchanged.
 * ==========================================================================*/

#define SHGO_SIMPLICIAL_MAXDIM 7      /* 2^7 = 128 corners, 7! = 5040 simplices */
#define SHGO_MAX_VERTICES      8000   /* hard ceiling on the complex size        */
#define SHGO_SOBOL_MAXDIM      21     /* dims served by the direction-number table */

/* Joe & Kuo (2008) Sobol' direction numbers for dimensions 2..SHGO_SOBOL_MAXDIM.
 * Row j (0-based, dim = j + 2): polynomial degree s, coefficient bits a, and the
 * s initial direction integers m_1..m_s (padded to the max degree 7). Dimension 1
 * is the van der Corput base-2 sequence and needs no table. The m-values are
 * validated at build time (each must be odd and < 2^k); a row that fails
 * validation degrades that coordinate to Halton, so a transcription slip can only
 * lower discrepancy, never corrupt the search. */
static const uint32_t SHGO_SOBOL_S[SHGO_SOBOL_MAXDIM - 1] = {
    1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7
};
static const uint32_t SHGO_SOBOL_A[SHGO_SOBOL_MAXDIM - 1] = {
    0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,19,22,25,1,4
};
static const uint32_t SHGO_SOBOL_M[SHGO_SOBOL_MAXDIM - 1][7] = {
    {1,0,0,0,0,0,0},      {1,3,0,0,0,0,0},      {1,3,1,0,0,0,0},
    {1,1,1,0,0,0,0},      {1,1,3,3,0,0,0},      {1,3,5,13,0,0,0},
    {1,1,5,5,17,0,0},     {1,1,5,5,5,0,0},      {1,1,7,11,19,0,0},
    {1,1,5,1,1,0,0},      {1,1,1,3,11,0,0},     {1,3,5,5,31,0,0},
    {1,3,3,9,7,49,0},     {1,1,1,15,21,21,0},   {1,3,1,13,27,49,0},
    {1,1,1,15,7,5,0},     {1,3,1,15,13,25,0},   {1,1,5,5,19,61,0},
    {1,3,7,11,23,15,103}, {1,3,7,9,9,25,29}
};

/* SplitMix64 mixer for the QMC digital/fractional shift (the "RandomSeed"). */
static uint64_t shgo_mix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/* Radical-inverse (van der Corput) in the given base — the Halton coordinate. */
static double shgo_radical_inverse(long i, int base) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= base; r += f * (double)(i % base); i /= base; }
    return r;
}

/* The (j+1)-th prime (0-based: j=0 -> 2), by trial division. j is small (one per
 * dimension), so the cost is negligible and no fixed table caps the dimension. */
static int shgo_nth_prime(size_t j) {
    size_t count = 0;
    int cand = 1;
    for (;;) {
        cand++;
        bool prime = true;
        for (int d = 2; d * d <= cand; d++) if (cand % d == 0) { prime = false; break; }
        if (prime) { if (count == j) return cand; count++; }
    }
}

/* Fill pts[n_samp * n] with Sobol' points (indices offset..offset+n_samp-1) in
 * [0,1)^n, digitally shifted by seed. Dimensions beyond the table (or that fail
 * validation) fall back to a seed-shifted Halton coordinate. */
static void shgo_sobol_fill(double* pts, size_t n_samp, size_t offset, size_t n,
                            uint64_t seed) {
    const uint32_t L = 32;
    size_t sd = n < (size_t)SHGO_SOBOL_MAXDIM ? n : (size_t)SHGO_SOBOL_MAXDIM;
    uint32_t*      Vdir  = (uint32_t*)malloc(sizeof(uint32_t) * sd * (L + 1));
    unsigned char* okdim = (unsigned char*)malloc(sd ? sd : 1);
    uint32_t*      shift = (uint32_t*)malloc(sizeof(uint32_t) * (sd ? sd : 1));
    uint32_t*      X     = (uint32_t*)calloc(sd ? sd : 1, sizeof(uint32_t));
    if (!Vdir || !okdim || !shift || !X) { free(Vdir); free(okdim); free(shift); free(X); return; }

    /* dimension 1: van der Corput (V_k = 2^(L-k)) */
    for (uint32_t k = 1; k <= L; k++) Vdir[0 * (L + 1) + k] = (uint32_t)1 << (L - k);
    okdim[0] = 1;
    for (size_t j = 1; j < sd; j++) {
        uint32_t s = SHGO_SOBOL_S[j - 1], a = SHGO_SOBOL_A[j - 1];
        const uint32_t* m = SHGO_SOBOL_M[j - 1];
        unsigned char good = 1;
        for (uint32_t k = 1; k <= s; k++)
            if ((m[k - 1] & 1u) == 0 || m[k - 1] >= ((uint32_t)1 << k)) { good = 0; break; }
        okdim[j] = good;
        if (!good) continue;
        uint32_t* Vj = &Vdir[j * (L + 1)];
        if (L <= s) {
            for (uint32_t k = 1; k <= L; k++) Vj[k] = m[k - 1] << (L - k);
        } else {
            for (uint32_t k = 1; k <= s; k++) Vj[k] = m[k - 1] << (L - k);
            for (uint32_t k = s + 1; k <= L; k++) {
                Vj[k] = Vj[k - s] ^ (Vj[k - s] >> s);
                for (uint32_t i = 1; i <= s - 1; i++)
                    if ((a >> (s - 1 - i)) & 1u) Vj[k] ^= Vj[k - i];
            }
        }
    }
    for (size_t j = 0; j < sd; j++)
        shift[j] = (uint32_t)shgo_mix(seed ^ (0xD1B54A32D192ED03ULL * (uint64_t)(j + 1)));

    for (size_t idx = 0; idx < offset + n_samp; idx++) {
        if (idx > 0) {
            uint32_t c = 1, val = (uint32_t)(idx - 1);
            while (val & 1u) { val >>= 1; c++; }
            if (c > L) c = L;
            for (size_t j = 0; j < sd; j++) if (okdim[j]) X[j] ^= Vdir[j * (L + 1) + c];
        }
        if (idx >= offset) {
            size_t row = idx - offset;
            for (size_t j = 0; j < n; j++) {
                if (j < sd && okdim[j]) {
                    pts[row * n + j] = (double)(X[j] ^ shift[j]) * (1.0 / 4294967296.0);
                } else {
                    double h = shgo_radical_inverse((long)(idx + 1), shgo_nth_prime(j));
                    double sf = (double)(shgo_mix(seed ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(j + 3))) >> 11)
                                * (1.0 / 9007199254740992.0);
                    double v = h + sf; v -= floor(v);
                    pts[row * n + j] = v;
                }
            }
        }
    }
    free(Vdir); free(okdim); free(shift); free(X);
}

/* Fill pts[n_samp * n] with Halton points (indices offset..) in [0,1)^n, each
 * coordinate fractionally shifted by seed (a Cranley-Patterson rotation). */
static void shgo_halton_fill(double* pts, size_t n_samp, size_t offset, size_t n,
                             uint64_t seed) {
    for (size_t j = 0; j < n; j++) {
        int base = shgo_nth_prime(j);
        double sf = (double)(shgo_mix(seed ^ (0x2545F4914F6CDD1DULL * (uint64_t)(j + 1))) >> 11)
                    * (1.0 / 9007199254740992.0);
        for (size_t i = 0; i < n_samp; i++) {
            double h = shgo_radical_inverse((long)(offset + i + 1), base);
            double v = h + sf; v -= floor(v);
            pts[i * n + j] = v;
        }
    }
}

static void shgo_cx_free(ShgoCx* c) { free(c->V); free(c->fv); free(c->pv); free(c->E);
                                      c->V = NULL; c->fv = NULL; c->pv = NULL; c->E = NULL; }

/* Append a vertex at coords (evaluated via nm_eval); returns its index, or
 * SIZE_MAX on allocation failure or when the vertex ceiling is reached. */
static size_t shgo_cx_add_vertex(ShgoCx* c, NmDriver* D, const double* coords) {
    if (c->nv >= SHGO_MAX_VERTICES) return SIZE_MAX;
    if (c->nv == c->vcap) {
        size_t nc = c->vcap ? c->vcap * 2 : 64;
        double* nV = (double*)realloc(c->V, sizeof(double) * nc * c->n);
        double* nf = (double*)realloc(c->fv, sizeof(double) * nc);
        double* np = (double*)realloc(c->pv, sizeof(double) * nc);
        if (!nV || !nf || !np) { c->V = nV ? nV : c->V; c->fv = nf ? nf : c->fv;
                                 c->pv = np ? np : c->pv; return SIZE_MAX; }
        c->V = nV; c->fv = nf; c->pv = np; c->vcap = nc;
    }
    size_t idx = c->nv;
    for (size_t j = 0; j < c->n; j++) c->V[idx * c->n + j] = coords[j];
    nm_eval(D, &c->V[idx * c->n], &c->fv[idx], &c->pv[idx]);
    c->nv++;
    return idx;
}

static void shgo_cx_add_edge(ShgoCx* c, size_t a, size_t b) {
    if (a == b) return;
    if (c->ne == c->ecap) {
        size_t nc = c->ecap ? c->ecap * 2 : 256;
        size_t* nE = (size_t*)realloc(c->E, sizeof(size_t) * nc * 2);
        if (!nE) return;
        c->E = nE; c->ecap = nc;
    }
    c->E[2 * c->ne] = a; c->E[2 * c->ne + 1] = b; c->ne++;
}

/* Recursively enumerate permutations of {0..n-1}; each permutation is a simplex
 * of the Kuhn triangulation whose vertices are the corners along the monotone
 * path 0 -> ... -> all-ones, and every pair among them is a 1-skeleton edge. The
 * adjacency bitset (nc*nc, nc = 2^n <= 128) dedups. */
static void shgo_kuhn_mark(unsigned char* adj, size_t nc, size_t n,
                           int* perm, unsigned used, size_t depth) {
    if (depth == n) {
        size_t path[SHGO_SIMPLICIAL_MAXDIM + 1];
        size_t cur = 0; path[0] = 0;
        for (size_t t = 0; t < n; t++) { cur |= ((size_t)1 << perm[t]); path[t + 1] = cur; }
        for (size_t a = 0; a <= n; a++)
            for (size_t b = a + 1; b <= n; b++) {
                size_t x = path[a], y = path[b];
                size_t lo = x < y ? x : y, hi = x < y ? y : x;
                adj[lo * nc + hi] = 1;
            }
        return;
    }
    for (int v = 0; v < (int)n; v++)
        if (!(used & (1u << v))) { perm[depth] = v; shgo_kuhn_mark(adj, nc, n, perm, used | (1u << v), depth + 1); }
}

/* Rebuild the kNN graph (used by the QMC samplers): each vertex is connected to
 * its k = min(nv-1, 2n+1) nearest neighbours (Euclidean), matching the ~2n
 * simplex-vertex degree so the minimizer-pool test behaves like the simplicial
 * one. Directed duplicates are harmless. O(nv^2 * k). */
static void shgo_knn_rebuild(ShgoCx* c) {
    c->ne = 0;
    size_t nv = c->nv, n = c->n;
    if (nv < 2) return;
    size_t k = (2 * n + 1 < nv - 1) ? (2 * n + 1) : (nv - 1);
    double* dist = (double*)malloc(sizeof(double) * nv);
    if (!dist) return;
    for (size_t i = 0; i < nv; i++) {
        for (size_t j = 0; j < nv; j++) {
            if (j == i) { dist[j] = INFINITY; continue; }
            double d2 = 0.0;
            for (size_t t = 0; t < n; t++) {
                double dd = c->V[i * n + t] - c->V[j * n + t];
                d2 += dd * dd;
            }
            dist[j] = d2;
        }
        for (size_t t = 0; t < k; t++) {
            size_t mn = nv; double mv = INFINITY;
            for (size_t j = 0; j < nv; j++) if (dist[j] < mv) { mv = dist[j]; mn = j; }
            if (mn == nv) break;
            shgo_cx_add_edge(c, i, mn);
            dist[mn] = INFINITY;
        }
    }
    free(dist);
}

/* Build the initial simplicial complex on the box: 2^n corners + the Kuhn
 * 1-skeleton, then longest-edge bisection until nv >= target. */
static void shgo_build_simplicial(ShgoCx* c, NmDriver* D, size_t target) {
    size_t n = c->n;
    size_t ncrn = (size_t)1 << n;
    double* coord = (double*)calloc(n ? n : 1, sizeof(double));
    if (!coord) return;
    for (size_t mask = 0; mask < ncrn; mask++) {
        for (size_t j = 0; j < n; j++)
            coord[j] = ((mask >> j) & 1) ? D->reg_hi[j] : D->reg_lo[j];
        shgo_cx_add_vertex(c, D, coord);      /* vertex index == corner mask */
    }
    free(coord);
    unsigned char* adj = (unsigned char*)calloc(ncrn * ncrn, 1);
    if (adj) {
        int perm[SHGO_SIMPLICIAL_MAXDIM];
        shgo_kuhn_mark(adj, ncrn, n, perm, 0u, 0);
        for (size_t a = 0; a < ncrn; a++)
            for (size_t b = a + 1; b < ncrn; b++)
                if (adj[a * ncrn + b]) shgo_cx_add_edge(c, a, b);
        free(adj);
    }
    /* longest-edge bisection: split the globally longest edge, connecting the new
     * midpoint to the edge endpoints and their common neighbours (the local
     * star), until the vertex target is met. */
    if (target > c->nv) {
        double* mid     = (double*)malloc(sizeof(double) * n);
        unsigned char* mark = (unsigned char*)calloc(SHGO_MAX_VERTICES, 1);
        size_t* commons = (size_t*)malloc(sizeof(size_t) * SHGO_MAX_VERTICES);
        while (mid && mark && commons && c->nv < target && c->nv < SHGO_MAX_VERTICES) {
            double best = -1.0; size_t be = SIZE_MAX;
            for (size_t e = 0; e < c->ne; e++) {
                size_t a = c->E[2 * e], b = c->E[2 * e + 1];
                double d2 = 0.0;
                for (size_t t = 0; t < n; t++) { double dd = c->V[a * n + t] - c->V[b * n + t]; d2 += dd * dd; }
                if (d2 > best) { best = d2; be = e; }
            }
            if (be == SIZE_MAX || best <= 0.0) break;
            size_t a = c->E[2 * be], b = c->E[2 * be + 1];
            /* collect common neighbours of a and b over the pre-insertion edges */
            memset(mark, 0, c->nv);
            size_t ne0 = c->ne, ncom = 0;
            for (size_t e = 0; e < ne0; e++) {
                if (c->E[2 * e] == a) mark[c->E[2 * e + 1]] = 1;
                else if (c->E[2 * e + 1] == a) mark[c->E[2 * e]] = 1;
            }
            for (size_t e = 0; e < ne0; e++) {
                size_t nb = SIZE_MAX;
                if (c->E[2 * e] == b) nb = c->E[2 * e + 1];
                else if (c->E[2 * e + 1] == b) nb = c->E[2 * e];
                if (nb != SIZE_MAX && nb != a && mark[nb] && ncom < SHGO_MAX_VERTICES) commons[ncom++] = nb;
            }
            for (size_t j = 0; j < n; j++) mid[j] = 0.5 * (c->V[a * n + j] + c->V[b * n + j]);
            size_t m = shgo_cx_add_vertex(c, D, mid);
            if (m == SIZE_MAX) break;
            /* remove the bisected edge (swap-with-last) and wire the new star */
            c->E[2 * be] = c->E[2 * (c->ne - 1)];
            c->E[2 * be + 1] = c->E[2 * (c->ne - 1) + 1];
            c->ne--;
            shgo_cx_add_edge(c, m, a);
            shgo_cx_add_edge(c, m, b);
            for (size_t t = 0; t < ncom; t++) shgo_cx_add_edge(c, m, commons[t]);
        }
        free(mid); free(mark); free(commons);
    }
}

/* Append a QMC batch (indices offset..offset+n_samp-1) and rebuild the kNN graph. */
static void shgo_build_qmc(ShgoCx* c, NmDriver* D, size_t n_samp, int method,
                           uint64_t seed, size_t offset) {
    size_t n = c->n;
    double* pts   = (double*)malloc(sizeof(double) * n_samp * n);
    double* coord = (double*)malloc(sizeof(double) * n);
    if (!pts || !coord) { free(pts); free(coord); return; }
    if (method == 1) shgo_sobol_fill(pts, n_samp, offset, n, seed);
    else             shgo_halton_fill(pts, n_samp, offset, n, seed);
    for (size_t i = 0; i < n_samp; i++) {
        for (size_t j = 0; j < n; j++)
            coord[j] = D->reg_lo[j] + pts[i * n + j] * (D->reg_hi[j] - D->reg_lo[j]);
        shgo_cx_add_vertex(c, D, coord);
    }
    free(pts); free(coord);
    shgo_knn_rebuild(c);
}

/* SHGO engine. Signature identical to the other NMinimize engines. */
void nm_shgo(NmDriver* D, const NmConfig* nc, NmRng* rng,
                    double* xbest, double* fbest, double* penbest) {
    (void)rng;                    /* SHGO is deterministic QMC + nc->seed         */
    size_t n = D->n;
    *fbest = 1e300; *penbest = 1e300;
    bool have = false;
    if (n == 0) { double f, p; double d0 = 0.0; nm_eval(D, &d0, &f, &p);
                  *fbest = f; *penbest = p; return; }

    int    sampling = nc->shgo_sampling;                        /* 0 simpl,1 sobol,2 halton */
    size_t n_samp   = nc->search_points > 0 ? (size_t)nc->search_points : 100;
    int    iters    = nc->shgo_iters > 0 ? nc->shgo_iters : 1;

    /* Simplicial is combinatorial in 2^n; above the practical ceiling fall back to
     * Sobol, warning once (scipy gives the same guidance). */
    if (sampling == 0 && n > (size_t)SHGO_SIMPLICIAL_MAXDIM) {
        fm_warn("NMinimize", "shgodim",
                "SHGO \"Simplicial\" sampling is impractical for %zu variables; "
                "using \"Sobol\" instead", n);
        sampling = 1;
    }

    ShgoCx cx; memset(&cx, 0, sizeof(cx)); cx.n = n;
    if (sampling == 0) shgo_build_simplicial(&cx, D, n_samp);
    else               shgo_build_qmc(&cx, D, n_samp, sampling, nc->seed, 0);
    if (cx.nv == 0) { shgo_cx_free(&cx); return; }             /* allocation failure */

    /* Basin-separation tolerance per coordinate (as in nm_de's multi-start). */
    double* btol = (double*)malloc(sizeof(double) * n);
    if (!btol) { shgo_cx_free(&cx); return; }
    for (size_t j = 0; j < n; j++) {
        double t = 1e-3 * (D->reg_hi[j] - D->reg_lo[j]);
        btol[j] = t > 1e-6 ? t : 1e-6;
    }

    /* per-iteration polish budget: the best K pool vertices (as in the plan) */
    size_t cap_polish = n_samp / 10; size_t two_n = 2 * n;
    if (cap_polish < two_n) cap_polish = two_n;
    if (cap_polish < 1) cap_polish = 1;

    /* The seed/minimum buffers hold at most cap_polish distinct entries per
     * iteration (each polish adds at most one); size to that bound, capped, and
     * guard the writes so a high-dimensional many-iteration run can never
     * overflow them (nseeds/nmins are the STORED counts, so the dedup scans stay
     * in bounds). */
    size_t slot_cap = cap_polish * (size_t)(iters > 0 ? iters : 1);
    if (slot_cap > SHGO_MAX_VERTICES) slot_cap = SHGO_MAX_VERTICES;
    if (slot_cap < 1) slot_cap = 1;

    /* Scratch reused across iterations. */
    unsigned char* is_min = (unsigned char*)malloc(SHGO_MAX_VERTICES);
    size_t*  order  = (size_t*)malloc(sizeof(size_t) * SHGO_MAX_VERTICES);
    unsigned char* used = (unsigned char*)malloc(SHGO_MAX_VERTICES);
    double*  seeds  = (double*)malloc(sizeof(double) * slot_cap * n);
    double*  minima = (double*)malloc(sizeof(double) * slot_cap * n);
    double*  xp     = (double*)malloc(sizeof(double) * n);
    if (!is_min || !order || !used || !seeds || !minima || !xp) {
        free(is_min); free(order); free(used); free(seeds); free(minima); free(xp);
        free(btol); shgo_cx_free(&cx); return;
    }
    size_t nseeds = 0, nmins = 0, prev_nmins = 0;

    for (int it = 0; it < iters; it++) {
        size_t nv = cx.nv;
        /* minimizer pool: is_min[i] iff i is strictly better than every neighbour */
        for (size_t i = 0; i < nv; i++) is_min[i] = 1;
        for (size_t e = 0; e < cx.ne; e++) {
            size_t a = cx.E[2 * e], b = cx.E[2 * e + 1];
            if (!nm_better(cx.fv[a], cx.pv[a], cx.fv[b], cx.pv[b])) is_min[a] = 0;
            if (!nm_better(cx.fv[b], cx.pv[b], cx.fv[a], cx.pv[a])) is_min[b] = 0;
        }
        /* global-best vertex is always a seed (covers ties / an empty pool) */
        size_t gb = 0;
        for (size_t i = 1; i < nv; i++)
            if (nm_better(cx.fv[i], cx.pv[i], cx.fv[gb], cx.pv[gb])) gb = i;
        is_min[gb] = 1;

        /* candidate list = pool vertices, ordered Deb-best first (selection) */
        size_t ncand = 0;
        for (size_t i = 0; i < nv; i++) if (is_min[i]) order[ncand++] = i;
        memset(used, 0, ncand);

        if (nc->post_process == 0) {
            /* PostProcess -> False: no local search; keep the Deb-best pool point. */
            for (size_t r = 0; r < ncand; r++) {
                size_t i = order[r];
                if (!have || nm_better(cx.fv[i], cx.pv[i], *fbest, *penbest)) {
                    for (size_t j = 0; j < n; j++) xbest[j] = cx.V[i * n + j];
                    *fbest = cx.fv[i]; *penbest = cx.pv[i]; have = true;
                }
            }
            break;                              /* no refinement without polishing */
        }

        size_t done = 0;
        while (done < cap_polish) {
            /* next Deb-best not-yet-used candidate */
            size_t bi = ncand;
            for (size_t r = 0; r < ncand; r++) {
                if (used[r]) continue;
                if (bi == ncand ||
                    nm_better(cx.fv[order[r]], cx.pv[order[r]], cx.fv[order[bi]], cx.pv[order[bi]]))
                    bi = r;
            }
            if (bi == ncand) break;
            used[bi] = 1;
            size_t i = order[bi];
            /* skip a candidate already inside a polished basin (seed dedup) */
            bool dup = false;
            for (size_t s = 0; s < nseeds && !dup; s++) {
                bool same = true;
                for (size_t j = 0; j < n; j++)
                    if (fabs(cx.V[i * n + j] - seeds[s * n + j]) > btol[j]) { same = false; break; }
                if (same) dup = true;
            }
            if (dup) continue;                  /* does not count against cap_polish */
            if (nseeds < slot_cap) {            /* record the basin seed (bounded) */
                for (size_t j = 0; j < n; j++) seeds[nseeds * n + j] = cx.V[i * n + j];
                nseeds++;
            }

            for (size_t j = 0; j < n; j++) xp[j] = cx.V[i * n + j];
            double fp = cx.fv[i], pp = cx.pv[i];
            nm_local_polish(D, xp, &fp, &pp);
            if (nm_better(cx.fv[i], cx.pv[i], fp, pp)) {   /* overshoot guard */
                for (size_t j = 0; j < n; j++) xp[j] = cx.V[i * n + j];
                fp = cx.fv[i]; pp = cx.pv[i];
            }
            /* record a distinct local minimum (homology group order) */
            bool known = false;
            for (size_t s = 0; s < nmins && !known; s++) {
                bool same = true;
                for (size_t j = 0; j < n; j++)
                    if (fabs(xp[j] - minima[s * n + j]) > btol[j]) { same = false; break; }
                if (same) known = true;
            }
            if (!known && nmins < slot_cap) {
                for (size_t j = 0; j < n; j++) minima[nmins * n + j] = xp[j];
                nmins++;
            }
            if (!have || nm_better(fp, pp, *fbest, *penbest)) {
                for (size_t j = 0; j < n; j++) xbest[j] = xp[j];
                *fbest = fp; *penbest = pp; have = true;
            }
            done++;
        }

        /* homology-growth stop: no new distinct minimum this round */
        if (it > 0 && nmins == prev_nmins) break;
        prev_nmins = nmins;
        if (it + 1 < iters && cx.nv < SHGO_MAX_VERTICES) {
            if (sampling == 0) shgo_build_simplicial(&cx, D, cx.nv + n_samp);
            else               shgo_build_qmc(&cx, D, n_samp, sampling, nc->seed, (size_t)(it + 1) * n_samp);
        }
    }

    free(is_min); free(order); free(used); free(seeds); free(minima); free(xp);
    free(btol); shgo_cx_free(&cx);
}
