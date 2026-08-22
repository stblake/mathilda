/* ============================================================================
 *  nm_direct.c — DIRECT (DIviding RECTangles) global optimizer for NMinimize.
 *
 *  Deterministic, derivative-free global search over a bounded hyperrectangle,
 *  after
 *    D. R. Jones, C. D. Perttunen & B. E. Stuckman, "Lipschitzian optimization
 *      without the Lipschitz constant", J. Optim. Theory Appl. 79 (1993) 157-181,
 *  with the locally-biased variant DIRECT-L of
 *    J. M. Gablonsky & C. T. Kelley, "A locally-biased form of the DIRECT
 *      algorithm", J. Global Optim. 21 (2001) 27-37.
 *
 *  This is the same algorithm the reference Fortran DIRECTv2.04 and
 *  scipy.optimize.direct implement, and the structure mirrors theirs:
 *
 *    - The search box is normalized to the unit hypercube [0,1]^n. A cell's
 *      side length along dimension j is 3^{-t[j]}, where t[j] counts how many
 *      times that dimension has been trisected; the cell center maps to a real
 *      point x[j] = lo[j] + c[j]*(hi[j]-lo[j]).
 *    - Cells are grouped by an integer "level" (size class): for DIRECT-L the
 *      level is the number of trisections of the longest side (min_j t[j]); for
 *      original DIRECT it is k*n+q with k = min_j t[j] and q = #{j : t[j] > k},
 *      which is a strictly size-decreasing key (see dc_size_of_level). Each
 *      level keeps a singly linked list of its cells ordered by value, so the
 *      list head is that level's minimum — the only cell in the class that can
 *      be "potentially optimal".
 *    - Each iteration selects the potentially-optimal cells (lower-right convex
 *      hull of the (size, value) level-heads, plus Jones's epsilon improvement
 *      test against the incumbent), then trisects each along its longest sides,
 *      sampling 2 new centers per longest dimension and dividing so the best
 *      children get the largest sub-cells.
 *
 *  Constrained / integer problems flow through the shared evaluation helpers:
 *  every center is scored with nm_eval (compiled fast path + integer rounding +
 *  constraint penalty) and the incumbent is chosen by Deb feasibility rules
 *  (nm_better); the value the algorithm minimizes is f + NM_PENALTY_MU*penalty
 *  (== f for a box-only problem, giving exact scipy parity). The engine returns
 *  the raw DIRECT incumbent; the NMinimize driver applies the exact local polish
 *  afterwards unless "PostProcess" -> False, so a raw-vs-raw race against
 *  scipy's direct() is available while the default result is polished.
 *
 *  Memory: contiguous structure-of-arrays stores grown by realloc doubling; no
 *  per-iteration reallocation of the whole store; all buffers are released on
 *  every exit path (including allocation failure). See nm_shgo.c for the sibling
 *  space-partitioning engine this follows.
 * ==========================================================================*/

#include "findmin_internal.h"

/* ------------------------------------------------------------------ *
 *  Engine context                                                    *
 * ------------------------------------------------------------------ */
typedef struct {
    NmDriver*     D;
    size_t        n;
    const double* lo;          /* box lower bounds (D->reg_lo)                 */
    const double* hi;          /* box upper bounds (D->reg_hi)                 */

    /* Cell store (structure of arrays), grown by realloc doubling. */
    size_t  nr;                /* number of cells                              */
    size_t  cap;               /* store capacity                               */
    double* c;                 /* cell centers in [0,1]^n, row-major nr*n      */
    int*    t;                 /* per-dim trisection counts, row-major nr*n    */
    double* fval;              /* objective value per cell                     */
    double* pval;              /* constraint penalty per cell                  */
    int*    next;              /* per-level linked-list successor (-1 = end)   */

    /* Per-level list heads: anchor[L] = index of the min-value cell at level L
     * (-1 = empty). Grown as deeper levels appear. */
    int*    anchor;
    size_t  anchor_cap;

    /* Algorithm parameters (resolved from NmConfig). */
    int     jones;             /* 0 original DIRECT, 1 DIRECT-L (locally biased)*/
    double  eps;               /* potentially-optimal slack                    */
    double  vol_tol;           /* stop when incumbent-cell volume < this       */
    double  len_tol;           /* stop when incumbent-cell size < this         */
    double  fmin;              /* known global minimum, -inf ⇒ inactive        */
    double  fmin_rtol;         /* rel. tol for the MinValue early stop         */

    /* Incumbent (best-so-far). best_rect is a stable store index (the store
     * never shrinks and a divided cell keeps its index as the central child). */
    bool    have;
    double  fbest, penbest;
    long    best_rect;
    long    nfev;

    bool    oom;               /* sticky allocation-failure flag               */

    /* Scratch reused across divisions (all sized n). */
    int*    dims;              /* indices of the longest dimensions (I in refs) */
    int*    order;            /* I sorted by best child value                  */
    double* wv;                /* min child value per longest dimension         */
    double* fp; double* penp;  /* + side child objective / penalty              */
    double* fm; double* penm;  /* - side child objective / penalty              */
    double* cbuf;              /* scratch child center                          */
    double* ccopy;             /* immutable copy of the divided cell's center   */
    int*    tbuf;              /* scratch child trisection counts               */
    int*    tcur;              /* mutable copy of the divided cell's counts     */
    double* xreal;             /* scratch mapped real point for nm_eval         */
} DirectCtx;

/* Scalar the algorithm minimizes: objective plus penalized infeasibility. */
static double dc_val(const DirectCtx* dc, size_t ri) {
    return dc->fval[ri] + NM_PENALTY_MU * dc->pval[ri];
}

/* Half the size measure of a level. DIRECT-L uses half the longest side
 * (0.5*3^{-k}); original DIRECT uses the half-diameter of the cell shape (k,q):
 *   0.5 * 3^{-k} * sqrt((n-q) + q/9)  = 0.5 * 3^{-k} * sqrt(n - (8/9)q).
 * Strictly decreasing in the level key, so distinct levels ⇒ distinct sizes. */
static double dc_size_of_level(const DirectCtx* dc, int L) {
    if (dc->jones) return 0.5 * pow(3.0, -(double)L);
    int n = (int)dc->n;
    int k = L / n, q = L % n;
    return 0.5 * pow(3.0, -(double)k) * sqrt((double)n - (8.0 / 9.0) * (double)q);
}

/* Integer size-class key of a cell (see file header). */
static int dc_level(const DirectCtx* dc, size_t ri) {
    size_t n = dc->n;
    const int* t = &dc->t[ri * n];
    int k = t[0];
    for (size_t j = 1; j < n; j++) if (t[j] < k) k = t[j];
    if (dc->jones) return k;
    int q = 0;
    for (size_t j = 0; j < n; j++) if (t[j] > k) q++;
    return k * (int)n + q;
}

/* ------------------------------------------------------------------ *
 *  Store / anchor growth                                             *
 * ------------------------------------------------------------------ */
static bool dc_reserve(DirectCtx* dc, size_t need) {
    if (need <= dc->cap) return true;
    size_t nc = dc->cap ? dc->cap * 2 : 64;
    while (nc < need) nc *= 2;
    size_t n = dc->n;
    double* c2 = (double*)realloc(dc->c,    sizeof(double) * nc * n);
    int*    t2 = (int*)   realloc(dc->t,    sizeof(int)    * nc * n);
    double* f2 = (double*)realloc(dc->fval, sizeof(double) * nc);
    double* p2 = (double*)realloc(dc->pval, sizeof(double) * nc);
    int*    x2 = (int*)   realloc(dc->next, sizeof(int)    * nc);
    if (c2) dc->c = c2;
    if (t2) dc->t = t2;
    if (f2) dc->fval = f2;
    if (p2) dc->pval = p2;
    if (x2) dc->next = x2;
    if (!c2 || !t2 || !f2 || !p2 || !x2) { dc->oom = true; return false; }
    dc->cap = nc;
    return true;
}

static bool dc_ensure_anchor(DirectCtx* dc, int L) {
    if ((size_t)L < dc->anchor_cap) return true;
    size_t nc = dc->anchor_cap ? dc->anchor_cap : 64;
    while ((size_t)L >= nc) nc *= 2;
    int* a = (int*)realloc(dc->anchor, sizeof(int) * nc);
    if (!a) { dc->oom = true; return false; }
    for (size_t i = dc->anchor_cap; i < nc; i++) a[i] = -1;
    dc->anchor = a;
    dc->anchor_cap = nc;
    return true;
}

/* Append a cell; update the incumbent by Deb rules. Returns the new index, or
 * -1 on allocation failure. Callers must not hold pointers into the store
 * across this call (it may realloc). */
static long dc_new_rect(DirectCtx* dc, const double* center, const int* t,
                        double f, double pen) {
    if (!dc_reserve(dc, dc->nr + 1)) return -1;
    size_t i = dc->nr++;
    size_t n = dc->n;
    memcpy(&dc->c[i * n], center, sizeof(double) * n);
    memcpy(&dc->t[i * n], t,      sizeof(int)    * n);
    dc->fval[i] = f;
    dc->pval[i] = pen;
    dc->next[i] = -1;
    if (!dc->have || nm_better(f, pen, dc->fbest, dc->penbest)) {
        dc->have = true;
        dc->fbest = f;
        dc->penbest = pen;
        dc->best_rect = (long)i;
    }
    return (long)i;
}

/* Insert a cell into its level's list, keeping the list ascending by value so
 * the head is the level minimum. O(list length); matches the reference. */
static bool dc_insert(DirectCtx* dc, size_t ri) {
    int L = dc_level(dc, ri);
    if (!dc_ensure_anchor(dc, L)) return false;
    double v = dc_val(dc, ri);
    int head = dc->anchor[L];
    if (head < 0 || v <= dc_val(dc, (size_t)head)) {
        dc->next[ri] = head;
        dc->anchor[L] = (int)ri;
        return true;
    }
    int prev = head;
    while (dc->next[prev] >= 0 && dc_val(dc, (size_t)dc->next[prev]) < v)
        prev = dc->next[prev];
    dc->next[ri] = dc->next[prev];
    dc->next[prev] = (int)ri;
    return true;
}

/* Map a unit-cube center to the real box and evaluate; bump the eval counter. */
static void dc_eval(DirectCtx* dc, const double* uc, double* f, double* pen) {
    size_t n = dc->n;
    for (size_t j = 0; j < n; j++)
        dc->xreal[j] = dc->lo[j] + uc[j] * (dc->hi[j] - dc->lo[j]);
    nm_eval(dc->D, dc->xreal, f, pen);
    dc->nfev++;
}

/* ------------------------------------------------------------------ *
 *  Potentially-optimal selection (DIRchoose)                        *
 * ------------------------------------------------------------------ *
 * Collect the min-value head of each non-empty level, keep those on the
 * lower-right convex hull of the (size, value) points that also pass Jones's
 * epsilon improvement test, DETACH the selected heads from their level lists,
 * and return their store indices. *out_n is the count; the caller frees the
 * result. Returns NULL when nothing is selected or on allocation failure. */
static long* dc_choose(DirectCtx* dc, size_t* out_n) {
    *out_n = 0;
    size_t K = 0;
    for (size_t L = 0; L < dc->anchor_cap; L++) if (dc->anchor[L] >= 0) K++;
    if (K == 0) return NULL;

    double* csz = (double*)malloc(sizeof(double) * K);
    double* cvl = (double*)malloc(sizeof(double) * K);
    int*    clv = (int*)   malloc(sizeof(int)    * K);
    long*   chd = (long*)  malloc(sizeof(long)   * K);
    long*   sel = (long*)  malloc(sizeof(long)   * K);
    int*    slv = (int*)   malloc(sizeof(int)    * K);
    if (!csz || !cvl || !clv || !chd || !sel || !slv) {
        dc->oom = true;
        free(csz); free(cvl); free(clv); free(chd); free(sel); free(slv);
        return NULL;
    }

    /* Candidates in size-descending order (level ascending). */
    size_t idx = 0;
    for (size_t L = 0; L < dc->anchor_cap; L++) {
        int h = dc->anchor[L];
        if (h < 0) continue;
        clv[idx] = (int)L;
        chd[idx] = h;
        csz[idx] = dc_size_of_level(dc, (int)L);
        cvl[idx] = dc_val(dc, (size_t)h);
        idx++;
    }

    double minf = dc->fbest + NM_PENALTY_MU * dc->penbest;
    double amp = fabs(minf);
    double rhs = minf - dc->eps * (amp > 1.0 ? amp : 1.0);

    size_t ns = 0;
    for (size_t j = 0; j < K; j++) {
        double helplower = HUGE_VAL;   /* min positive slope to a larger cell  */
        double helpgreater = 0.0;      /* max positive slope to a smaller cell */
        for (size_t i = 0; i < K; i++) {
            if (i == j) continue;
            double ds = csz[i] - csz[j];
            if (ds == 0.0) continue;
            double slope = (cvl[i] - cvl[j]) / ds;
            if (ds > 0.0) { if (slope > 0.0 && slope < helplower) helplower = slope; }
            else          { if (slope > helpgreater) helpgreater = slope; }
        }
        bool accept = false;
        if (helpgreater <= helplower) {
            if (helplower == HUGE_VAL) accept = true;   /* largest cell: always PO */
            else if (cvl[j] - helplower * csz[j] <= rhs) accept = true;
        }
        if (accept) { sel[ns] = chd[j]; slv[ns] = clv[j]; ns++; }
    }

    long* result = NULL;
    if (ns > 0) {
        result = (long*)malloc(sizeof(long) * ns);
        if (!result) {
            dc->oom = true;
        } else {
            for (size_t s = 0; s < ns; s++) {
                int L = slv[s];
                long h = sel[s];
                dc->anchor[L] = dc->next[h];   /* detach the head (O(1))        */
                dc->next[h] = -1;
                result[s] = h;
            }
            *out_n = ns;
        }
    }
    free(csz); free(cvl); free(clv); free(chd); free(sel); free(slv);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Division (DIRsamplepoints + DIRdivide)                           *
 * ------------------------------------------------------------------ *
 * Trisect a (detached) cell along its longest sides: sample c ± delta e_j for
 * each longest dimension j, then divide in ascending order of the better child
 * value so the best children keep the largest sub-cells. The pair split off at
 * rank r inherits the parent counts with dims order[0..r] incremented; the
 * central sub-cell (all longest dims incremented) reuses this cell's index. */
static void dc_divide(DirectCtx* dc, size_t r) {
    size_t n = dc->n;
    memcpy(dc->ccopy, &dc->c[r * n], sizeof(double) * n);
    memcpy(dc->tcur,  &dc->t[r * n], sizeof(int)    * n);

    int k = dc->tcur[0];
    for (size_t j = 1; j < n; j++) if (dc->tcur[j] < k) k = dc->tcur[j];
    double delta = pow(3.0, -(double)(k + 1));

    int m = 0;
    for (size_t j = 0; j < n; j++) if (dc->tcur[j] == k) dc->dims[m++] = (int)j;

    /* Sample both side centers per longest dimension. */
    for (int a = 0; a < m; a++) {
        int j = dc->dims[a];
        memcpy(dc->cbuf, dc->ccopy, sizeof(double) * n);
        dc->cbuf[j] = dc->ccopy[j] + delta;
        dc_eval(dc, dc->cbuf, &dc->fp[a], &dc->penp[a]);
        dc->cbuf[j] = dc->ccopy[j] - delta;
        dc_eval(dc, dc->cbuf, &dc->fm[a], &dc->penm[a]);
        double vp = dc->fp[a] + NM_PENALTY_MU * dc->penp[a];
        double vm = dc->fm[a] + NM_PENALTY_MU * dc->penm[a];
        dc->wv[a] = vp < vm ? vp : vm;
    }

    /* order = ranks 0..m-1 sorted by wv ascending (insertion sort; m ≤ n). */
    for (int a = 0; a < m; a++) dc->order[a] = a;
    for (int a = 1; a < m; a++) {
        int key = dc->order[a];
        double kw = dc->wv[key];
        int b = a - 1;
        while (b >= 0 && dc->wv[dc->order[b]] > kw) { dc->order[b + 1] = dc->order[b]; b--; }
        dc->order[b + 1] = key;
    }

    for (int rank = 0; rank < m; rank++) {
        int a = dc->order[rank];
        int j = dc->dims[a];
        memcpy(dc->tbuf, dc->tcur, sizeof(int) * n);
        dc->tbuf[j] = dc->tcur[j] + 1;

        memcpy(dc->cbuf, dc->ccopy, sizeof(double) * n);
        dc->cbuf[j] = dc->ccopy[j] + delta;
        long ip = dc_new_rect(dc, dc->cbuf, dc->tbuf, dc->fp[a], dc->penp[a]);
        if (ip < 0) return;
        if (!dc_insert(dc, (size_t)ip)) return;

        memcpy(dc->cbuf, dc->ccopy, sizeof(double) * n);
        dc->cbuf[j] = dc->ccopy[j] - delta;
        long im = dc_new_rect(dc, dc->cbuf, dc->tbuf, dc->fm[a], dc->penm[a]);
        if (im < 0) return;
        if (!dc_insert(dc, (size_t)im)) return;

        dc->tcur[j] += 1;   /* commit: the central column keeps dividing along j */
    }

    /* This cell is now the central sub-cell: write back the deepened counts and
     * re-file it (its center and value are unchanged, so the incumbent, if it is
     * this cell, still points at the same real point — now in a smaller box). */
    memcpy(&dc->t[r * n], dc->tcur, sizeof(int) * n);
    dc_insert(dc, r);
}

/* ------------------------------------------------------------------ *
 *  Driver                                                           *
 * ------------------------------------------------------------------ */
static void dc_free(DirectCtx* dc) {
    free(dc->c); free(dc->t); free(dc->fval); free(dc->pval); free(dc->next);
    free(dc->anchor);
    free(dc->dims); free(dc->order); free(dc->wv);
    free(dc->fp); free(dc->penp); free(dc->fm); free(dc->penm);
    free(dc->cbuf); free(dc->ccopy); free(dc->tbuf); free(dc->tcur); free(dc->xreal);
}

void nm_direct(NmDriver* D, const NmConfig* nc, NmRng* rng,
                      double* xbest, double* fbest, double* penbest) {
    (void)rng;                          /* DIRECT is deterministic              */
    size_t n = D->n;
    *fbest = 1e300;
    *penbest = 1e300;
    if (n == 0) {                       /* degenerate: nothing to divide        */
        double d0 = 0.0, f, p;
        nm_eval(D, &d0, &f, &p);
        *fbest = f;
        *penbest = p;
        return;
    }

    DirectCtx dc;
    memset(&dc, 0, sizeof dc);
    dc.D = D;
    dc.n = n;
    dc.lo = D->reg_lo;
    dc.hi = D->reg_hi;
    dc.best_rect = -1;
    dc.jones     = (nc->direct_locally_biased == 0) ? 0 : 1;   /* default DIRECT-L */
    dc.eps       = nc->direct_eps      >= 0.0 ? nc->direct_eps      : NM_DIRECT_EPS;
    dc.vol_tol   = nc->direct_vol_tol  >= 0.0 ? nc->direct_vol_tol  : NM_DIRECT_VOLTOL;
    dc.len_tol   = nc->direct_len_tol  >= 0.0 ? nc->direct_len_tol  : NM_DIRECT_LENTOL;
    dc.fmin      = nc->direct_fmin;
    dc.fmin_rtol = nc->direct_fmin_rtol >= 0.0 ? nc->direct_fmin_rtol : NM_DIRECT_FMINRTOL;

    long maxfun = nc->direct_max_fun > 0
                ? (long)nc->direct_max_fun
                : (long)NM_DIRECT_MAXFUN_PER_DIM * (long)n;
    if (maxfun > NM_DIRECT_MAXFUN_CAP) maxfun = NM_DIRECT_MAXFUN_CAP;
    if (maxfun < 2 * (long)n + 1) maxfun = 2 * (long)n + 1;   /* room for first divide */
    int maxiter = nc->direct_max_iter > 0 ? nc->direct_max_iter : NM_DIRECT_MAXITER;

    dc.dims  = (int*)   malloc(sizeof(int)    * n);
    dc.order = (int*)   malloc(sizeof(int)    * n);
    dc.wv    = (double*)malloc(sizeof(double) * n);
    dc.fp    = (double*)malloc(sizeof(double) * n);
    dc.penp  = (double*)malloc(sizeof(double) * n);
    dc.fm    = (double*)malloc(sizeof(double) * n);
    dc.penm  = (double*)malloc(sizeof(double) * n);
    dc.cbuf  = (double*)malloc(sizeof(double) * n);
    dc.ccopy = (double*)malloc(sizeof(double) * n);
    dc.tbuf  = (int*)   malloc(sizeof(int)    * n);
    dc.tcur  = (int*)   malloc(sizeof(int)    * n);
    dc.xreal = (double*)malloc(sizeof(double) * n);
    if (!dc.dims || !dc.order || !dc.wv || !dc.fp || !dc.penp || !dc.fm || !dc.penm ||
        !dc.cbuf || !dc.ccopy || !dc.tbuf || !dc.tcur || !dc.xreal) {
        dc_free(&dc);
        return;
    }

    /* Initial cell: center of the unit cube. */
    for (size_t j = 0; j < n; j++) { dc.cbuf[j] = 0.5; dc.tbuf[j] = 0; }
    {
        double f, p;
        dc_eval(&dc, dc.cbuf, &f, &p);
        long i0 = dc_new_rect(&dc, dc.cbuf, dc.tbuf, f, p);
        if (i0 < 0 || !dc_insert(&dc, (size_t)i0)) { dc_free(&dc); return; }
    }

    for (int it = 0; it < maxiter && dc.nfev < maxfun; it++) {
        size_t nsel = 0;
        long* sel = dc_choose(&dc, &nsel);
        if (dc.oom) { free(sel); break; }
        if (nsel == 0) { free(sel); break; }

        for (size_t s = 0; s < nsel; s++) {
            dc_divide(&dc, (size_t)sel[s]);
            if (dc.oom || dc.nfev >= maxfun) break;
        }
        free(sel);
        if (dc.oom) break;

        /* Resolution stops, measured on the cell that currently holds the best
         * point (it shrinks as that lineage is subdivided). */
        if (dc.best_rect >= 0) {
            int L = dc_level(&dc, (size_t)dc.best_rect);
            double half = dc_size_of_level(&dc, L);
            const int* bt = &dc.t[(size_t)dc.best_rect * n];
            long sumt = 0;
            for (size_t j = 0; j < n; j++) sumt += bt[j];
            double vol = pow(3.0, -(double)sumt);
            if (half < dc.len_tol) break;
            if (vol < dc.vol_tol) break;
        }
        if (dc.fmin > -HUGE_VAL && dc.have) {
            double den = fabs(dc.fmin);
            if (den < 1.0) den = 1.0;
            if ((dc.fbest - dc.fmin) / den < dc.fmin_rtol) break;
        }
    }

    if (dc.have && dc.best_rect >= 0) {
        const double* bc = &dc.c[(size_t)dc.best_rect * n];
        for (size_t j = 0; j < n; j++)
            xbest[j] = dc.lo[j] + bc[j] * (dc.hi[j] - dc.lo[j]);
        nm_project(D, xbest);            /* clamp to box, snap integer coords    */
        nm_eval(D, xbest, fbest, penbest);
    }

    dc_free(&dc);
}
