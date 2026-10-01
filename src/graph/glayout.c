/* glayout.c - deterministic graph layout engine (see glayout.h).
 *
 * METHODS
 *   Stress      SMACOF stress majorization on BFS graph distances (Gansner,
 *               Koren & North 2004), with weights d^-2, initialised by Pivot MDS
 *               (Brandes & Pich 2006). Components above STRESS_EXACT_MAX
 *               vertices skip the O(n^2)-per-sweep majorization and keep the
 *               Pivot MDS layout, which is O(k (n + m)) for k pivots.
 *   Spring      Hu's spring-electrical model (attraction d^2/K, repulsion
 *               C K^2 / d) with adaptive step length, initialised by Pivot MDS.
 *               Above SPRING_EXACT_MAX vertices repulsion is cut off at a
 *               radius and evaluated over a uniform grid (Fruchterman-Reingold).
 *   Layered     Trees: a tidy tree (leaves in DFS order, parents centred over
 *               their children) from the tree centre, or from the source of an
 *               arborescence. DAGs: longest-path layering, then barycentre
 *               sweeps to reduce crossings (best ordering kept), then x
 *               positions by isotonic regression (pool-adjacent-violators)
 *               towards neighbour barycentres with unit minimum separation.
 *               Other graphs: BFS layers from the graph centre.
 *   Circular, Bipartite, Grid   closed forms, barycentre-ordered for Bipartite.
 *
 * Every component is laid out on its own, normalised to unit mean edge length,
 * rotated to a canonical orientation (principal axis horizontal; snapped to the
 * axes when the edges are mostly axis-parallel, as in grids), and the
 * components are then shelf-packed, largest first.
 *
 * DETERMINISM. No random numbers anywhere: pivots, orders, starts and
 * tie-breaks are all fixed functions of the vertex indices, so a graph always
 * gets bit-identical coordinates.
 */

#include "glayout.h"
#include "sym_names.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define STRESS_EXACT_MAX 1000   /* full SMACOF up to this component size       */
#define STRESS_MAX_ITERS 400
#define STRESS_MULTISTART 60    /* components this small try several starts   */
#define SPRING_EXACT_MAX 1000   /* exact O(n^2) repulsion up to this size       */
#define SPRING_MAX_ITERS 500
#define PIVOTS 50
#define CENTER_EXACT_MAX 2000   /* exact graph centre (all BFS) up to this size */

/* ------------------------------------------------------------ CSR graph --- */

typedef struct {
    int n, m;                 /* vertices, edges                              */
    int *eu, *ev;             /* edge endpoints                               */
    unsigned char* dir;       /* directed flags (never NULL inside)           */
    int *off, *adj;           /* undirected adjacency, CSR                    */
} Sub;

static void sub_free(Sub* s) {
    free(s->eu); free(s->ev); free(s->dir); free(s->off); free(s->adj);
    memset(s, 0, sizeof(*s));
}

/* Builds the symmetric CSR adjacency of s's edge list. 0 on OOM. */
static int sub_csr(Sub* s) {
    s->off = calloc((size_t)s->n + 1, sizeof(int));
    s->adj = malloc(sizeof(int) * (size_t)(2 * s->m + 1));
    if (!s->off || !s->adj) return 0;
    for (int k = 0; k < s->m; k++) { s->off[s->eu[k] + 1]++; s->off[s->ev[k] + 1]++; }
    for (int i = 0; i < s->n; i++) s->off[i + 1] += s->off[i];
    int* fill = malloc(sizeof(int) * (size_t)(s->n + 1));
    if (!fill) return 0;
    memcpy(fill, s->off, sizeof(int) * (size_t)s->n);
    for (int k = 0; k < s->m; k++) {
        s->adj[fill[s->eu[k]]++] = s->ev[k];
        s->adj[fill[s->ev[k]]++] = s->eu[k];
    }
    free(fill);
    return 1;
}

/* BFS from src over s; d[] gets hop distances (-1 unreachable). queue is
 * scratch of size n. */
static void bfs(const Sub* s, int src, int* d, int* queue) {
    for (int i = 0; i < s->n; i++) d[i] = -1;
    int h = 0, t = 0;
    d[src] = 0; queue[t++] = src;
    while (h < t) {
        int u = queue[h++];
        for (int a = s->off[u]; a < s->off[u + 1]; a++) {
            int w = s->adj[a];
            if (d[w] < 0) { d[w] = d[u] + 1; queue[t++] = w; }
        }
    }
}

/* ---------------------------------------------------- small linear algebra */

/* Cyclic Jacobi eigen-decomposition of the symmetric k x k matrix A (row
 * major, destroyed). Eigenvalues to w[], eigenvectors as COLUMNS of V. */
static void jacobi_eigen(int k, double* A, double* w, double* V) {
    for (int i = 0; i < k * k; i++) V[i] = 0.0;
    for (int i = 0; i < k; i++) V[i * k + i] = 1.0;
    for (int sweep = 0; sweep < 60; sweep++) {
        double off = 0.0, diag = 0.0;
        for (int p = 0; p < k; p++)
            for (int q = 0; q < k; q++) {
                if (p == q) diag += A[p * k + q] * A[p * k + q];
                else off += A[p * k + q] * A[p * k + q];
            }
        if (off <= 1e-22 * (diag + 1e-300)) break;
        for (int p = 0; p < k - 1; p++)
            for (int q = p + 1; q < k; q++) {
                double apq = A[p * k + q];
                if (fabs(apq) < 1e-300) continue;
                double app = A[p * k + p], aqq = A[q * k + q];
                double theta = (aqq - app) / (2.0 * apq);
                double t = (theta >= 0 ? 1.0 : -1.0)
                         / (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0), sn = t * c;
                for (int r = 0; r < k; r++) {           /* columns p, q */
                    double arp = A[r * k + p], arq = A[r * k + q];
                    A[r * k + p] = c * arp - sn * arq;
                    A[r * k + q] = sn * arp + c * arq;
                }
                for (int r = 0; r < k; r++) {           /* rows p, q */
                    double apr = A[p * k + r], aqr = A[q * k + r];
                    A[p * k + r] = c * apr - sn * aqr;
                    A[q * k + r] = sn * apr + c * aqr;
                }
                for (int r = 0; r < k; r++) {
                    double vrp = V[r * k + p], vrq = V[r * k + q];
                    V[r * k + p] = c * vrp - sn * vrq;
                    V[r * k + q] = sn * vrp + c * vrq;
                }
            }
    }
    for (int i = 0; i < k; i++) w[i] = A[i * k + i];
}

/* ------------------------------------------------------------ utilities --- */

/* Scales xy so the mean edge length is 1 (no-op without edges or when all
 * edges are degenerate). */
static void normalise_edge_length(const Sub* s, double* xy) {
    if (s->m == 0) return;
    double tot = 0.0;
    for (int k = 0; k < s->m; k++) {
        double dx = xy[2 * s->eu[k]] - xy[2 * s->ev[k]];
        double dy = xy[2 * s->eu[k] + 1] - xy[2 * s->ev[k] + 1];
        tot += sqrt(dx * dx + dy * dy);
    }
    if (tot <= 1e-300) return;
    double f = (double)s->m / tot;
    for (int i = 0; i < 2 * s->n; i++) xy[i] *= f;
}

/* Tiny deterministic perturbation that breaks exact symmetries (coincident or
 * collinear starts) without any randomness. */
static void jitter(int n, double* xy, double amp) {
    for (int i = 0; i < n; i++) {
        xy[2 * i]     += amp * sin(1.0 + 2.399963 * (double)i);
        xy[2 * i + 1] += amp * cos(1.0 + 2.399963 * (double)i);
    }
}

/* Rotates to the canonical orientation: principal axis horizontal, then
 * snapped to the axes when the edge directions are strongly four-fold (grids,
 * cubes), then reflected so that vertex 0 sits left of and above the centre. */
static void canonical_orientation(const Sub* s, double* xy) {
    int n = s->n;
    if (n < 2) return;
    double cx = 0, cy = 0;
    for (int i = 0; i < n; i++) { cx += xy[2 * i]; cy += xy[2 * i + 1]; }
    cx /= n; cy /= n;
    double sxx = 0, syy = 0, sxy = 0;
    for (int i = 0; i < n; i++) {
        double dx = xy[2 * i] - cx, dy = xy[2 * i + 1] - cy;
        sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
    }
    double th = 0.5 * atan2(2.0 * sxy, sxx - syy);   /* principal axis angle */
    /* Only rotate when there is a genuine principal axis. */
    double aniso = sqrt((sxx - syy) * (sxx - syy) + 4 * sxy * sxy) / (sxx + syy + 1e-300);
    double rot = aniso > 0.05 ? -th : 0.0;
    /* Four-fold snap, measured after the principal rotation. */
    int snapped = 0;
    if (s->m > 0) {
        double S = 0, C = 0;
        for (int k = 0; k < s->m; k++) {
            double dx = xy[2 * s->ev[k]] - xy[2 * s->eu[k]];
            double dy = xy[2 * s->ev[k] + 1] - xy[2 * s->eu[k] + 1];
            double ph = atan2(dy, dx) + rot;
            S += sin(4 * ph); C += cos(4 * ph);
        }
        double R = sqrt(S * S + C * C) / s->m;
        if (R > 0.6) { rot -= atan2(S, C) / 4.0; snapped = 1; }
    }
    if (aniso <= 0.05 && !snapped) {
        /* No preferred axis (cycles, K_n, vertex-transitive graphs): stand
         * vertex 0 straight above the centre, as a textbook polygon does. */
        double dx = xy[0] - cx, dy = xy[1] - cy;
        if (dx * dx + dy * dy > 1e-18) rot = M_PI / 2 - atan2(dy, dx);
    }
    double c = cos(rot), sn = sin(rot);
    for (int i = 0; i < n; i++) {
        double dx = xy[2 * i] - cx, dy = xy[2 * i + 1] - cy;
        xy[2 * i] = c * dx - sn * dy;
        xy[2 * i + 1] = sn * dx + c * dy;
    }
    /* Clean up rounding noise so an axis-aligned layout is exactly aligned. */
    for (int i = 0; i < 2 * n; i++) if (fabs(xy[i]) < 1e-12) xy[i] = 0.0;
    if (xy[0] > 1e-9) for (int i = 0; i < n; i++) xy[2 * i] = -xy[2 * i];
    if (xy[1] < -1e-9) for (int i = 0; i < n; i++) xy[2 * i + 1] = -xy[2 * i + 1];
}

/* Union-find root with path halving. */
static int uf_find(int* p, int i) {
    while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
    return i;
}

/* Orthogonal straightening. When nearly every edge of an oriented drawing is
 * within 15 degrees of an axis (grids, ladders, tori drawn flat), stress leaves
 * the rows gently bowed, because hop distance is not Euclidean distance along
 * a diagonal. Vertices joined by near-horizontal edges then share one y (their
 * mean), and by near-vertical edges one x, so a grid comes out as a grid. */
static void straighten(const Sub* s, double* xy) {
    int n = s->n, m = s->m, aligned = 0;
    if (m < 4) return;
    const double tol = tan(15.0 * M_PI / 180.0);
    for (int k = 0; k < m; k++) {
        double dx = fabs(xy[2 * s->ev[k]] - xy[2 * s->eu[k]]);
        double dy = fabs(xy[2 * s->ev[k] + 1] - xy[2 * s->eu[k] + 1]);
        if (dy <= tol * dx || dx <= tol * dy) aligned++;
    }
    if (aligned < m) return;              /* every edge must be axis-like */
    int* p = malloc(sizeof(int) * (size_t)n);
    double* sum = malloc(sizeof(double) * (size_t)n);
    int* cnt = malloc(sizeof(int) * (size_t)n);
    if (p && sum && cnt) {
        for (int axis = 0; axis < 2; axis++) {  /* 0: rows share y; 1: columns share x */
            for (int i = 0; i < n; i++) { p[i] = i; sum[i] = 0; cnt[i] = 0; }
            for (int k = 0; k < m; k++) {
                double dx = fabs(xy[2 * s->ev[k]] - xy[2 * s->eu[k]]);
                double dy = fabs(xy[2 * s->ev[k] + 1] - xy[2 * s->eu[k] + 1]);
                int horiz = dy <= tol * dx;
                if (horiz == (axis == 0)) {
                    int a = uf_find(p, s->eu[k]), b = uf_find(p, s->ev[k]);
                    if (a != b) p[a < b ? b : a] = a < b ? a : b;
                }
            }
            int co = axis == 0 ? 1 : 0;
            for (int i = 0; i < n; i++) { int r = uf_find(p, i); sum[r] += xy[2 * i + co]; cnt[r]++; }
            for (int i = 0; i < n; i++) { int r = uf_find(p, i); xy[2 * i + co] = sum[r] / cnt[r]; }
        }
    }
    free(p); free(sum); free(cnt);
}

/* --------------------------------------------------------------- circle --- */

static void layout_circle(int n, const int* order, double* xy) {
    if (n == 1) { xy[0] = xy[1] = 0; return; }
    double R = (n == 2) ? 0.5 : 0.5 / sin(M_PI / n);    /* unit chord */
    for (int r = 0; r < n; r++) {
        int i = order ? order[r] : r;
        double t = M_PI / 2 + 2.0 * M_PI * r / n;
        xy[2 * i] = R * cos(t); xy[2 * i + 1] = R * sin(t);
    }
}

/* ------------------------------------------------------------ Pivot MDS --- */

/* Pivot MDS of s into xy. D (n x n, may be NULL) supplies distances when the
 * full matrix is already known. 0 on OOM. */
static int pivot_mds(const Sub* s, const int* D, double* xy) {
    int n = s->n;
    if (n <= 2) {
        xy[0] = 0; xy[1] = 0;
        if (n == 2) { xy[2] = 1; xy[3] = 0; }
        return 1;
    }
    int k = n < PIVOTS ? n : PIVOTS;
    double* C = malloc(sizeof(double) * (size_t)n * (size_t)k);
    int* piv = malloc(sizeof(int) * (size_t)k);
    int* mind = malloc(sizeof(int) * (size_t)n);
    int* d = malloc(sizeof(int) * (size_t)n);
    int* q = malloc(sizeof(int) * (size_t)n);
    double* M = malloc(sizeof(double) * (size_t)k * (size_t)k);
    double* V = malloc(sizeof(double) * (size_t)k * (size_t)k);
    double* w = malloc(sizeof(double) * (size_t)k);
    double* cm = calloc((size_t)k, sizeof(double));
    int ok = C && piv && mind && d && q && M && V && w && cm;
    if (ok) {
        /* First pivot: highest degree (lowest index on ties); then max-min. */
        int best = 0;
        for (int i = 1; i < n; i++)
            if (s->off[i + 1] - s->off[i] > s->off[best + 1] - s->off[best]) best = i;
        for (int i = 0; i < n; i++) mind[i] = 1 << 29;
        for (int j = 0; j < k; j++) {
            piv[j] = best;
            if (D) for (int i = 0; i < n; i++) d[i] = D[(size_t)best * n + i];
            else bfs(s, best, d, q);
            for (int i = 0; i < n; i++) {
                double dd = d[i] < 0 ? n : d[i];
                C[(size_t)i * k + j] = dd * dd;
                if (d[i] >= 0 && d[i] < mind[i]) mind[i] = d[i];
            }
            best = 0;
            for (int i = 1; i < n; i++) if (mind[i] > mind[best]) best = i;
        }
        /* Double centring. */
        double gm = 0;
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < n; i++) cm[j] += C[(size_t)i * k + j];
            cm[j] /= n; gm += cm[j];
        }
        gm /= k;
        for (int i = 0; i < n; i++) {
            double rm = 0;
            for (int j = 0; j < k; j++) rm += C[(size_t)i * k + j];
            rm /= k;
            for (int j = 0; j < k; j++)
                C[(size_t)i * k + j] = -0.5 * (C[(size_t)i * k + j] - rm - cm[j] + gm);
        }
        for (int a = 0; a < k; a++)
            for (int b = a; b < k; b++) {
                double acc = 0;
                for (int i = 0; i < n; i++) acc += C[(size_t)i * k + a] * C[(size_t)i * k + b];
                M[a * k + b] = M[b * k + a] = acc;
            }
        jacobi_eigen(k, M, w, V);
        int e1 = 0;
        for (int i = 1; i < k; i++) if (w[i] > w[e1]) e1 = i;
        int e2 = (e1 == 0) ? 1 : 0;
        for (int i = 0; i < k; i++) if (i != e1 && w[i] > w[e2]) e2 = i;
        for (int i = 0; i < n; i++) {
            double x = 0, y = 0;
            for (int j = 0; j < k; j++) {
                x += C[(size_t)i * k + j] * V[j * k + e1];
                y += C[(size_t)i * k + j] * V[j * k + e2];
            }
            xy[2 * i] = x; xy[2 * i + 1] = y;
        }
        normalise_edge_length(s, xy);
    }
    free(C); free(piv); free(mind); free(d); free(q); free(M); free(V); free(w); free(cm);
    return ok;
}

/* --------------------------------------------------------------- stress --- */

static double stress_value(int n, const int* D, const double* xy) {
    double st = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            double dij = D[(size_t)i * n + j];
            double dx = xy[2 * i] - xy[2 * j], dy = xy[2 * i + 1] - xy[2 * j + 1];
            double e = sqrt(dx * dx + dy * dy) - dij;
            st += e * e / (dij * dij);
        }
    return st;
}

/* Localised SMACOF (in-place majorization sweeps) from the start in xy.
 * Returns the final stress. */
static void smacof_sweep(int n, const int* D, double* xy) {
        for (int i = 0; i < n; i++) {
            double sx = 0, sy = 0, sw = 0;
            const int* Di = D + (size_t)i * n;
            double xi = xy[2 * i], yi = xy[2 * i + 1];
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                double dij = Di[j];
                double w = 1.0 / (dij * dij);
                double dx = xi - xy[2 * j], dy = yi - xy[2 * j + 1];
                double dist = sqrt(dx * dx + dy * dy);
                sw += w;
                if (dist > 1e-12) {
                    sx += w * (xy[2 * j] + dij * dx / dist);
                    sy += w * (xy[2 * j + 1] + dij * dy / dist);
                } else {
                    sx += w * xy[2 * j]; sy += w * xy[2 * j + 1];
                }
            }
            if (sw > 0) { xy[2 * i] = sx / sw; xy[2 * i + 1] = sy / sw; }
        }
}

static double smacof(int n, const int* D, double* xy) {
    double prev = stress_value(n, D, xy);
    for (int it = 0; it < STRESS_MAX_ITERS; it++) {
        smacof_sweep(n, D, xy);
        if (it % 4 == 3) {
            double cur = stress_value(n, D, xy);
            if (prev - cur <= 1e-6 * prev + 1e-12) { prev = cur; break; }
            prev = cur;
        }
    }
    return stress_value(n, D, xy);
}

/* All-pairs BFS distances (n x n ints), or NULL on OOM. */
static int* all_pairs(const Sub* s) {
    int n = s->n;
    int* D = malloc(sizeof(int) * (size_t)n * (size_t)n);
    int* q = malloc(sizeof(int) * (size_t)n);
    if (!D || !q) { free(D); free(q); return NULL; }
    for (int i = 0; i < n; i++) bfs(s, i, D + (size_t)i * n, q);
    free(q);
    return D;
}

/* Proper crossings between non-adjacent edges of a straight-line drawing. */
static long drawing_crossings(const Sub* s, const double* xy) {
    long c = 0;
    for (int a = 0; a < s->m; a++) {
        int p = s->eu[a], q = s->ev[a];
        double ax = xy[2 * p], ay = xy[2 * p + 1], bx = xy[2 * q], by = xy[2 * q + 1];
        for (int b = a + 1; b < s->m; b++) {
            int r = s->eu[b], t = s->ev[b];
            if (r == p || r == q || t == p || t == q) continue;
            double cx = xy[2 * r], cy = xy[2 * r + 1], dx = xy[2 * t], dy = xy[2 * t + 1];
            double d1 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
            double d2 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax);
            double d3 = (dx - cx) * (ay - cy) - (dy - cy) * (ax - cx);
            double d4 = (dx - cx) * (by - cy) - (dy - cy) * (bx - cx);
            if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0))
                && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) c++;
        }
    }
    return c;
}

/* A shortest cycle through v (BFS with branch labels): writes its vertices in
 * order to cyc and returns its length, or 0 if v lies on no cycle. */
static int shortest_cycle_through(const Sub* s, int v, int* cyc,
                                  int* d, int* par, int* br, int* q) {
    int n = s->n;
    for (int i = 0; i < n; i++) { d[i] = -1; par[i] = -1; br[i] = -1; }
    int h = 0, t = 0, best = 1 << 30, ba = -1, bb = -1;
    d[v] = 0; q[t++] = v;
    while (h < t) {
        int u = q[h++];
        if (2 * d[u] + 1 >= best) break;
        for (int a = s->off[u]; a < s->off[u + 1]; a++) {
            int w = s->adj[a];
            if (d[w] < 0) {
                d[w] = d[u] + 1; par[w] = u; br[w] = (u == v) ? w : br[u];
                q[t++] = w;
            } else if (w != par[u] && u != v && w != v && br[w] != br[u]) {
                int len = d[u] + d[w] + 1;
                if (len < best) { best = len; ba = u; bb = w; }
            } else if (w == v && u != v && par[u] != v) {
                /* closing edge straight back to v */
                int len = d[u] + 1;
                if (len < best) { best = len; ba = u; bb = v; }
            }
        }
    }
    if (ba < 0) return 0;
    int k = 0;
    for (int x = ba; x != v; x = par[x]) cyc[k++] = x;
    cyc[k++] = v;
    /* reverse so the cycle reads v .. ba, then append bb .. (towards v) */
    for (int i = 0; i < k / 2; i++) { int tmp = cyc[i]; cyc[i] = cyc[k - 1 - i]; cyc[k - 1 - i] = tmp; }
    if (bb != v) for (int x = bb; x != v; x = par[x]) cyc[k++] = x;
    return k;
}

/* Tutte-style start: the cycle cyc on a circle, every other vertex at the
 * barycentre of its neighbours (Gauss-Seidel). Crossing-free for a
 * 3-connected planar graph whose cycle is a face. */
static void tutte_start(const Sub* s, const int* cyc, int k, double* xy, unsigned char* pinned) {
    int n = s->n;
    memset(pinned, 0, (size_t)n);
    for (int i = 0; i < 2 * n; i++) xy[i] = 0;
    double R = (0.5 / sin(M_PI / k)) * (n > k ? sqrt((double)n / k) : 1.0);
    for (int i = 0; i < k; i++) {
        double t = M_PI / 2 + 2.0 * M_PI * i / k;
        xy[2 * cyc[i]] = R * cos(t); xy[2 * cyc[i] + 1] = R * sin(t); pinned[cyc[i]] = 1;
    }
    for (int it = 0; it < 300; it++)
        for (int i = 0; i < n; i++) {
            if (pinned[i] || s->off[i + 1] == s->off[i]) continue;
            double ax = 0, ay = 0; int na = s->off[i + 1] - s->off[i];
            for (int a = s->off[i]; a < s->off[i + 1]; a++) { ax += xy[2 * s->adj[a]]; ay += xy[2 * s->adj[a] + 1]; }
            xy[2 * i] = ax / na; xy[2 * i + 1] = ay / na;
        }
}

/* Rescales xy by the stress-optimal factor and returns the resulting stress,
 * so drawings of different scale compare fairly. */
static double scaled_stress(int n, const int* D, double* xy) {
    double num = 0, den = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            double dij = D[(size_t)i * n + j], w = 1.0 / (dij * dij);
            double dx = xy[2 * i] - xy[2 * j], dy = xy[2 * i + 1] - xy[2 * j + 1];
            double e = sqrt(dx * dx + dy * dy);
            num += w * dij * e; den += w * e * e;
        }
    if (den > 0) for (int i = 0; i < 2 * n; i++) xy[i] *= num / den;
    return stress_value(n, D, xy);
}

/* Planarity-preserving refinement of a Tutte drawing: Jacobi SMACOF sweeps
 * that move only the unpinned (interior) vertices -- the outer cycle stays a
 * regular polygon, and Jacobi updates keep every symmetry of the start --
 * stopped at the last crossing-free iterate. `nxt` is scratch of 2n doubles.
 * Returns the stress after optimal rescaling (so it compares with others). */
static double smacof_planar(const Sub* s, const int* D, double* xy,
                            const unsigned char* pinned, double* nxt) {
    int n = s->n;
    double prev = stress_value(n, D, xy);
    for (int it = 0; it < STRESS_MAX_ITERS; it++) {
        for (int i = 0; i < n; i++) {
            nxt[2 * i] = xy[2 * i]; nxt[2 * i + 1] = xy[2 * i + 1];
            if (pinned[i]) continue;
            double sx = 0, sy = 0, sw = 0;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                double dij = D[(size_t)i * n + j], w = 1.0 / (dij * dij);
                double dx = xy[2 * i] - xy[2 * j], dy = xy[2 * i + 1] - xy[2 * j + 1];
                double dist = sqrt(dx * dx + dy * dy);
                sw += w;
                sx += w * (xy[2 * j] + (dist > 1e-12 ? dij * dx / dist : 0));
                sy += w * (xy[2 * j + 1] + (dist > 1e-12 ? dij * dy / dist : 0));
            }
            if (sw > 0) { nxt[2 * i] = sx / sw; nxt[2 * i + 1] = sy / sw; }
        }
        if (drawing_crossings(s, nxt) > 0) break;
        memcpy(xy, nxt, sizeof(double) * 2 * (size_t)n);
        double cur = stress_value(n, D, xy);
        if (prev - cur <= 1e-6 * prev + 1e-12) break;
        prev = cur;
    }
    return scaled_stress(n, D, xy);
}

static int layout_stress(const Sub* s, double* xy) {
    int n = s->n;
    if (n <= 2) return pivot_mds(s, NULL, xy);
    if (n > STRESS_EXACT_MAX) return pivot_mds(s, NULL, xy);
    int* D = all_pairs(s);
    if (!D) return 0;
    int ok = pivot_mds(s, D, xy);
    if (ok) {
        jitter(n, xy, 1e-3);
        double best = smacof(n, D, xy);
        if (n <= STRESS_MULTISTART && s->m <= 400) {
            /* Extra deterministic starts -- a circle in VertexList order, and
             * circles in BFS order from several roots -- and keep, among the
             * layouts within 10% of the least stress, the one with the fewest
             * edge crossings (least stress on ties). */
            int ncirc = 1 + (n < 8 ? n : 8);
            int ntutte = 3;
            int nstart = ncirc + ntutte;
            int nall = nstart + ntutte;     /* + the Tutte drawings unrefined */
            double* alt = malloc(sizeof(double) * 2 * (size_t)n);
            double* cand = malloc(sizeof(double) * 2 * (size_t)n * (size_t)(nall + 1));
            double* st = malloc(sizeof(double) * (size_t)(nall + 1));
            long* cr = malloc(sizeof(long) * (size_t)(nall + 1));
            int* order = malloc(sizeof(int) * (size_t)n);
            int* sc = malloc(sizeof(int) * 5 * (size_t)n);
            unsigned char* pin = malloc((size_t)n);
            double* pin_xy = malloc(sizeof(double) * 2 * (size_t)n);
            if (alt && cand && st && cr && order && sc && pin && pin_xy) {
                memcpy(cand, xy, sizeof(double) * 2 * (size_t)n);
                st[0] = best; cr[0] = drawing_crossings(s, xy);
                for (int k = 0; k < nstart; k++) {
                    if (k >= ncirc) {
                        /* Tutte starts from shortest cycles through vertices
                         * spread over the index range. */
                        int v = (int)(((long)(k - ncirc) * n) / ntutte);
                        int len = shortest_cycle_through(s, v, order, sc, sc + n, sc + 2 * n, sc + 3 * n);
                        int raw = nstart + 1 + (k - ncirc);
                        if (len < 3) {
                            st[k + 1] = st[raw] = 1e300; cr[k + 1] = cr[raw] = 1L << 40;
                            continue;
                        }
                        tutte_start(s, order, len, alt, pin);
                        {   /* keep the unrefined drawing as a candidate too */
                            double* c = cand + 2 * (size_t)n * (size_t)raw;
                            memcpy(c, alt, sizeof(double) * 2 * (size_t)n);
                            cr[raw] = drawing_crossings(s, c);
                            st[raw] = cr[raw] == 0 ? smacof_planar(s, D, c, pin, pin_xy)
                                                   : scaled_stress(n, D, c);
                        }
                    } else if (k == 0) layout_circle(n, NULL, alt);
                    else {
                        int root = (int)(((long)(k - 1) * n) / (ncirc - 1));
                        const int* Dr = D + (size_t)root * n;
                        int t = 0;
                        for (int lev = 0; t < n; lev++)
                            for (int i = 0; i < n; i++) if (Dr[i] == lev) order[t++] = i;
                        layout_circle(n, order, alt);
                    }
                    jitter(n, alt, 1e-3);
                    st[k + 1] = smacof(n, D, alt);
                    cr[k + 1] = drawing_crossings(s, alt);
                    memcpy(cand + 2 * (size_t)n * (size_t)(k + 1), alt, sizeof(double) * 2 * (size_t)n);
                    if (st[k + 1] < best) best = st[k + 1];
                }
                /* Crossings of the least-stress drawing: a crossing-free
                 * (Tutte-derived) drawing may cost up to 2.5x the stress, but
                 * only when it removes many crossings (a dodecahedron's ten,
                 * not a cube's two -- the Necker cube is the textbook cube). */
                long crbest = 1L << 40;
                for (int k = 0; k <= nall; k++)
                    if (st[k] <= best * (1 + 1e-9) && cr[k] < crbest) crbest = cr[k];
                int planar_ok = crbest >= 4 && crbest * 6 >= s->m;
                int pick = -1;
                for (int k = 0; k <= nall; k++) {
                    double win = (cr[k] == 0 && planar_ok) ? 2.5 : 1.10;
                    if (st[k] > best * win + 1e-12) continue;
                    if (pick < 0 || cr[k] < cr[pick]
                        || (cr[k] == cr[pick] && st[k] < st[pick] * (1 - 1e-9))) pick = k;
                }
                if (pick > 0) memcpy(xy, cand + 2 * (size_t)n * (size_t)pick, sizeof(double) * 2 * (size_t)n);
            }
            free(alt); free(cand); free(st); free(cr); free(order); free(sc); free(pin); free(pin_xy);
        }
    }
    free(D);
    return ok;
}

/* ----------------------------------------------------- spring-electrical --- */

static int layout_spring(const Sub* s, double* xy) {
    int n = s->n;
    if (!pivot_mds(s, NULL, xy)) return 0;
    if (n <= 2) return 1;
    jitter(n, xy, 1e-3);
    const double K = 1.0, Cr = 0.2, cool = 0.9;
    double* f = malloc(sizeof(double) * 2 * (size_t)n);
    int* head = NULL; int* next = NULL;
    int use_grid = n > SPRING_EXACT_MAX;
    if (use_grid) next = malloc(sizeof(int) * (size_t)n);
    if (!f || (use_grid && !next)) { free(f); free(next); return 0; }
    double step = K, energy0 = 1e300;
    int progress = 0;
    const double R = 4.0 * K;               /* grid repulsion cut-off */
    for (int it = 0; it < SPRING_MAX_ITERS; it++) {
        memset(f, 0, sizeof(double) * 2 * (size_t)n);
        if (!use_grid) {
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++) {
                    double dx = xy[2 * i] - xy[2 * j], dy = xy[2 * i + 1] - xy[2 * j + 1];
                    double d2 = dx * dx + dy * dy + 1e-12;
                    double g = Cr * K * K / d2;
                    f[2 * i] += g * dx; f[2 * i + 1] += g * dy;
                    f[2 * j] -= g * dx; f[2 * j + 1] -= g * dy;
                }
        } else {
            double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
            for (int i = 0; i < n; i++) {
                if (xy[2 * i] < x0) x0 = xy[2 * i];
                if (xy[2 * i] > x1) x1 = xy[2 * i];
                if (xy[2 * i + 1] < y0) y0 = xy[2 * i + 1];
                if (xy[2 * i + 1] > y1) y1 = xy[2 * i + 1];
            }
            int gx = (int)((x1 - x0) / R) + 1, gy = (int)((y1 - y0) / R) + 1;
            if ((double)gx * gy > 4.0 * n) { gx = gy = (int)sqrt(4.0 * n) + 1; }
            double cwx = (x1 - x0) / gx + 1e-9, cwy = (y1 - y0) / gy + 1e-9;
            free(head);
            head = malloc(sizeof(int) * (size_t)gx * (size_t)gy);
            if (!head) break;
            for (int c = 0; c < gx * gy; c++) head[c] = -1;
            for (int i = 0; i < n; i++) {
                int cx = (int)((xy[2 * i] - x0) / cwx), cy = (int)((xy[2 * i + 1] - y0) / cwy);
                if (cx >= gx) cx = gx - 1;
                if (cy >= gy) cy = gy - 1;
                int c = cy * gx + cx;
                next[i] = head[c]; head[c] = i;
            }
            for (int i = 0; i < n; i++) {
                int cx = (int)((xy[2 * i] - x0) / cwx), cy = (int)((xy[2 * i + 1] - y0) / cwy);
                if (cx >= gx) cx = gx - 1;
                if (cy >= gy) cy = gy - 1;
                for (int ox = -1; ox <= 1; ox++)
                    for (int oy = -1; oy <= 1; oy++) {
                        int ax = cx + ox, ay = cy + oy;
                        if (ax < 0 || ay < 0 || ax >= gx || ay >= gy) continue;
                        for (int j = head[ay * gx + ax]; j >= 0; j = next[j]) {
                            if (j == i) continue;
                            double dx = xy[2 * i] - xy[2 * j], dy = xy[2 * i + 1] - xy[2 * j + 1];
                            double d2 = dx * dx + dy * dy + 1e-12;
                            if (d2 > R * R) continue;
                            double g = Cr * K * K / d2;
                            f[2 * i] += g * dx; f[2 * i + 1] += g * dy;
                        }
                    }
            }
        }
        for (int k = 0; k < s->m; k++) {
            int a = s->eu[k], b = s->ev[k];
            double dx = xy[2 * a] - xy[2 * b], dy = xy[2 * a + 1] - xy[2 * b + 1];
            double d = sqrt(dx * dx + dy * dy);
            double g = d / K;                          /* |F| = d^2 / K */
            f[2 * a] -= g * dx; f[2 * a + 1] -= g * dy;
            f[2 * b] += g * dx; f[2 * b + 1] += g * dy;
        }
        double energy = 0, moved = 0;
        for (int i = 0; i < n; i++) {
            double fx = f[2 * i], fy = f[2 * i + 1];
            double fl = sqrt(fx * fx + fy * fy);
            energy += fl * fl;
            if (fl > 1e-300) {
                xy[2 * i] += step * fx / fl; xy[2 * i + 1] += step * fy / fl;
                moved += step;
            }
        }
        if (energy < energy0) {
            if (++progress >= 5) { progress = 0; step /= cool; }
        } else { progress = 0; step *= cool; }
        energy0 = energy;
        if (moved < 1e-3 * K * n) break;
    }
    free(f); free(head); free(next);
    normalise_edge_length(s, xy);
    return 1;
}

/* -------------------------------------------------------------- layered --- */

/* Pool-adjacent-violators: the least-squares non-decreasing fit of y[0..n-1],
 * in place. blk/sum/cnt are scratch of size n. */
static void pava(int n, double* y, int* start, double* sum, int* cnt) {
    int nb = 0;
    for (int i = 0; i < n; i++) {
        start[nb] = i; sum[nb] = y[i]; cnt[nb] = 1; nb++;
        while (nb > 1 && sum[nb - 2] / cnt[nb - 2] > sum[nb - 1] / cnt[nb - 1]) {
            sum[nb - 2] += sum[nb - 1]; cnt[nb - 2] += cnt[nb - 1]; nb--;
        }
    }
    for (int b = 0; b < nb; b++) {
        double v = sum[b] / cnt[b];
        for (int i = start[b]; i < start[b] + cnt[b]; i++) y[i] = v;
    }
}

/* Longest-path layering over directed edges. Returns 1 and fills layer[] if
 * every edge is directed and they form a DAG; 0 otherwise (or on OOM). */
static int dag_layers(const Sub* s, int* layer) {
    int n = s->n;
    for (int k = 0; k < s->m; k++) if (!s->dir[k]) return 0;
    int* indeg = calloc((size_t)n, sizeof(int));
    int* off = calloc((size_t)n + 1, sizeof(int));
    int* out = malloc(sizeof(int) * (size_t)(s->m + 1));
    int* q = malloc(sizeof(int) * (size_t)n);
    int ok = indeg && off && out && q;
    if (ok) {
        for (int k = 0; k < s->m; k++) { indeg[s->ev[k]]++; off[s->eu[k] + 1]++; }
        for (int i = 0; i < n; i++) off[i + 1] += off[i];
        int* fill = malloc(sizeof(int) * (size_t)(n + 1));
        if (!fill) ok = 0;
        else {
            memcpy(fill, off, sizeof(int) * (size_t)n);
            for (int k = 0; k < s->m; k++) out[fill[s->eu[k]]++] = s->ev[k];
            free(fill);
            int h = 0, t = 0;
            for (int i = 0; i < n; i++) { layer[i] = 0; if (!indeg[i]) q[t++] = i; }
            while (h < t) {
                int u = q[h++];
                for (int a = off[u]; a < off[u + 1]; a++) {
                    int w = out[a];
                    if (layer[u] + 1 > layer[w]) layer[w] = layer[u] + 1;
                    if (--indeg[w] == 0) q[t++] = w;
                }
            }
            ok = (t == n);
        }
    }
    free(indeg); free(off); free(out); free(q);
    return ok;
}

/* A root for BFS layering: the graph centre (minimum eccentricity, lowest
 * index on ties) when affordable, else the vertex of maximum degree. */
static int layout_root(const Sub* s) {
    int n = s->n, best = 0;
    if (n <= CENTER_EXACT_MAX) {
        int* d = malloc(sizeof(int) * (size_t)n);
        int* q = malloc(sizeof(int) * (size_t)n);
        if (d && q) {
            int bestecc = 1 << 30;
            for (int i = 0; i < n; i++) {
                bfs(s, i, d, q);
                int ecc = 0;
                for (int j = 0; j < n; j++) if (d[j] > ecc) ecc = d[j];
                if (ecc < bestecc) { bestecc = ecc; best = i; }
            }
            free(d); free(q);
            return best;
        }
        free(d); free(q);
    }
    for (int i = 1; i < n; i++)
        if (s->off[i + 1] - s->off[i] > s->off[best + 1] - s->off[best]) best = i;
    return best;
}

/* Tidy tree from root: leaves at consecutive x in DFS order, each parent
 * centred over its first and last child; y = -depth. */
static int layout_tree(const Sub* s, int root, double* xy) {
    int n = s->n;
    int* parent = malloc(sizeof(int) * (size_t)n);
    int* depth = malloc(sizeof(int) * (size_t)n);
    int* pre = malloc(sizeof(int) * (size_t)n);
    int* stack = malloc(sizeof(int) * (size_t)(n + 1));
    int* first = malloc(sizeof(int) * (size_t)n);
    int* last = malloc(sizeof(int) * (size_t)n);
    int ok = parent && depth && pre && stack && first && last;
    if (ok) {
        for (int i = 0; i < n; i++) { parent[i] = -2; first[i] = last[i] = -1; }
        int sp = 0, np = 0;
        stack[sp++] = root; parent[root] = -1; depth[root] = 0;
        while (sp > 0) {
            int u = stack[--sp];
            pre[np++] = u;
            /* push children in reverse adjacency order so they pop in order */
            for (int a = s->off[u + 1] - 1; a >= s->off[u]; a--) {
                int w = s->adj[a];
                if (parent[w] != -2) continue;
                parent[w] = u; depth[w] = depth[u] + 1;
                stack[sp++] = w;
            }
        }
        /* first/last child in preorder */
        for (int t = 0; t < np; t++) {
            int u = pre[t], p = parent[u];
            if (p < 0) continue;
            if (first[p] < 0) first[p] = u;
            last[p] = u;
        }
        double leaf = 0;
        for (int t = 0; t < np; t++) {
            int u = pre[t];
            if (first[u] < 0) { xy[2 * u] = leaf; leaf += 1.0; }
            xy[2 * u + 1] = -(double)depth[u];
        }
        for (int t = np - 1; t >= 0; t--) {
            int u = pre[t];
            if (first[u] >= 0) xy[2 * u] = 0.5 * (xy[2 * first[u]] + xy[2 * last[u]]);
        }
    }
    free(parent); free(depth); free(pre); free(stack); free(first); free(last);
    return ok;
}

static int cmp_key_idx(const void* a, const void* b) {
    const double* x = (const double*)a; const double* y = (const double*)b;
    if (x[0] < y[0]) return -1;
    if (x[0] > y[0]) return 1;
    return (x[1] < y[1]) ? -1 : (x[1] > y[1]);
}

/* Crossings between edges joining consecutive layers (O(E^2) per layer pair,
 * called only for modest edge counts). */
static long count_crossings(const Sub* s, const int* layer, const int* pos) {
    long c = 0;
    for (int a = 0; a < s->m; a++) {
        int u1 = s->eu[a], v1 = s->ev[a];
        if (layer[u1] > layer[v1]) { int t = u1; u1 = v1; v1 = t; }
        if (layer[v1] - layer[u1] != 1) continue;
        for (int b = a + 1; b < s->m; b++) {
            int u2 = s->eu[b], v2 = s->ev[b];
            if (layer[u2] > layer[v2]) { int t = u2; u2 = v2; v2 = t; }
            if (layer[u2] != layer[u1] || layer[v2] != layer[v1]) continue;
            if ((long)(pos[u1] - pos[u2]) * (pos[v1] - pos[v2]) < 0) c++;
        }
    }
    return c;
}

/* Generic layered drawing given layer[]. */
static int layout_layers(const Sub* s, const int* layer, double* xy) {
    int n = s->n, nl = 0;
    for (int i = 0; i < n; i++) if (layer[i] + 1 > nl) nl = layer[i] + 1;
    int* lcount = calloc((size_t)nl + 1, sizeof(int));
    int* lstart = calloc((size_t)nl + 1, sizeof(int));
    int* order = malloc(sizeof(int) * (size_t)n);    /* vertices by layer, ranked */
    int* pos = malloc(sizeof(int) * (size_t)n);
    int* bestpos = malloc(sizeof(int) * (size_t)n);
    double* key = malloc(sizeof(double) * 2 * (size_t)n);
    double* y = malloc(sizeof(double) * (size_t)n);
    double* sum = malloc(sizeof(double) * (size_t)n);
    int* st = malloc(sizeof(int) * (size_t)n);
    int* cnt = malloc(sizeof(int) * (size_t)n);
    int ok = lcount && lstart && order && pos && bestpos && key && y && sum && st && cnt;
    if (ok) {
        for (int i = 0; i < n; i++) lcount[layer[i]]++;
        for (int l = 0; l < nl; l++) lstart[l + 1] = lstart[l] + lcount[l];
        { int* fill = calloc((size_t)nl, sizeof(int));
          if (!fill) ok = 0;
          else {
              for (int i = 0; i < n; i++) {
                  int l = layer[i];
                  order[lstart[l] + fill[l]] = i; pos[i] = fill[l]++;
              }
              free(fill);
          } }
    }
    if (ok) {
        int track = s->m <= 3000;
        long best = track ? count_crossings(s, layer, pos) : 0;
        memcpy(bestpos, pos, sizeof(int) * (size_t)n);
        for (int sweep = 0; sweep < 12 && (!track || best > 0); sweep++) {
            int down = (sweep % 2 == 0);
            for (int li = 1; li < nl; li++) {
                int l = down ? li : nl - 1 - li;
                int cntl = lcount[l];
                for (int r = 0; r < cntl; r++) {
                    int v = order[lstart[l] + r];
                    double acc = 0; int na = 0;
                    for (int a = s->off[v]; a < s->off[v + 1]; a++) {
                        int w = s->adj[a];
                        if (down ? layer[w] < l : layer[w] > l) { acc += pos[w]; na++; }
                    }
                    key[2 * r] = na ? acc / na : (double)pos[v];
                    key[2 * r + 1] = (double)pos[v];
                    (void)v;
                }
                /* sort (key, current pos); recover vertices via old ranks */
                for (int r = 0; r < cntl; r++) st[r] = order[lstart[l] + r];
                for (int r = 0; r < cntl; r++) key[2 * r + 1] = (double)r;
                qsort(key, (size_t)cntl, 2 * sizeof(double), cmp_key_idx);
                for (int r = 0; r < cntl; r++) {
                    int v = st[(int)key[2 * r + 1]];
                    order[lstart[l] + r] = v; pos[v] = r;
                }
            }
            if (track) {
                long c = count_crossings(s, layer, pos);
                if (c < best) { best = c; memcpy(bestpos, pos, sizeof(int) * (size_t)n); }
            }
        }
        if (track) {
            memcpy(pos, bestpos, sizeof(int) * (size_t)n);
            for (int i = 0; i < n; i++) order[lstart[layer[i]] + pos[i]] = i;
        }
        /* x coordinates: start at centred ranks, then pull each layer towards
         * its neighbours' barycentres subject to unit separation. */
        for (int i = 0; i < n; i++) xy[2 * i] = pos[i] - 0.5 * (lcount[layer[i]] - 1);
        for (int pass = 0; pass < 8; pass++) {
            int down = (pass % 2 == 0);
            for (int li = 0; li < nl; li++) {
                int l = down ? li : nl - 1 - li;
                int cntl = lcount[l];
                for (int r = 0; r < cntl; r++) {
                    int v = order[lstart[l] + r];
                    double acc = 0; int na = 0;
                    for (int a = s->off[v]; a < s->off[v + 1]; a++) {
                        int w = s->adj[a];
                        if (layer[w] != l) { acc += xy[2 * w]; na++; }
                    }
                    double want = na ? acc / na : xy[2 * v];
                    y[r] = want - r;
                }
                pava(cntl, y, st, sum, cnt);
                for (int r = 0; r < cntl; r++) xy[2 * order[lstart[l] + r]] = y[r] + r;
            }
        }
        /* Unit separation within a layer reads cramped against unit layer
         * spacing once arrowheads are drawn: widen the layers a little. */
        for (int i = 0; i < n; i++) { xy[2 * i] *= 1.3; xy[2 * i + 1] = -(double)layer[i]; }
    }
    free(lcount); free(lstart); free(order); free(pos); free(bestpos);
    free(key); free(y); free(sum); free(st); free(cnt);
    return ok;
}

/* Sugiyama-style: an edge spanning several layers is routed through one
 * dummy vertex per intermediate layer, so crossing reduction and x placement
 * see it, and the straight edge then clears the vertices in between. */
static int layout_layers_dummy(const Sub* s, const int* layer, double* xy) {
    long extra = 0;
    for (int k = 0; k < s->m; k++) {
        int sp = abs(layer[s->eu[k]] - layer[s->ev[k]]);
        if (sp > 1) extra += sp - 1;
    }
    if (extra == 0 || extra > 20000) return layout_layers(s, layer, xy);
    Sub a; memset(&a, 0, sizeof(a));
    a.n = s->n + (int)extra; a.m = s->m + (int)extra;
    a.eu = malloc(sizeof(int) * (size_t)a.m);
    a.ev = malloc(sizeof(int) * (size_t)a.m);
    a.dir = calloc((size_t)a.m, 1);
    int* L = malloc(sizeof(int) * (size_t)a.n);
    double* axy = malloc(sizeof(double) * 2 * (size_t)a.n);
    int ok = a.eu && a.ev && a.dir && L && axy;
    if (ok) {
        memcpy(L, layer, sizeof(int) * (size_t)s->n);
        int nv = s->n, ne = 0;
        for (int k = 0; k < s->m; k++) {
            int u = s->eu[k], v = s->ev[k];
            if (layer[u] > layer[v]) { int t = u; u = v; v = t; }
            int prev = u;
            for (int l = layer[u] + 1; l < layer[v]; l++) {
                L[nv] = l;
                a.eu[ne] = prev; a.ev[ne] = nv; ne++;
                prev = nv++;
            }
            a.eu[ne] = prev; a.ev[ne] = v; ne++;
        }
        ok = sub_csr(&a) && layout_layers(&a, L, axy);
        if (ok) memcpy(xy, axy, sizeof(double) * 2 * (size_t)s->n);
    }
    sub_free(&a); free(L); free(axy);
    return ok;
}

static int layout_layered(const Sub* s, double* xy) {
    int n = s->n;
    if (n == 1) { xy[0] = xy[1] = 0; return 1; }
    int* layer = malloc(sizeof(int) * (size_t)n);
    if (!layer) return 0;
    int ok;
    int dag = dag_layers(s, layer);
    if (s->m == n - 1) {
        /* A tree: tidy drawing. An arborescence hangs from its source. */
        int root = -1;
        if (dag) {
            int nsrc = 0;
            for (int i = 0; i < n; i++) if (layer[i] == 0) { nsrc++; root = i; }
            if (nsrc != 1) root = -1;
        }
        if (root < 0) root = layout_root(s);
        ok = layout_tree(s, root, xy);
    } else {
        if (!dag) {
            int* q = malloc(sizeof(int) * (size_t)n);
            if (!q) { free(layer); return 0; }
            bfs(s, layout_root(s), layer, q);
            free(q);
        }
        ok = layout_layers_dummy(s, layer, xy);
    }
    if (ok) {
        /* Wide, shallow drawings (big trees) get taller layer spacing so the
         * picture is not a flat strip: height at least 0.35 x width. */
        double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
        for (int i = 0; i < n; i++) {
            if (xy[2 * i] < x0) x0 = xy[2 * i];
            if (xy[2 * i] > x1) x1 = xy[2 * i];
            if (xy[2 * i + 1] < y0) y0 = xy[2 * i + 1];
            if (xy[2 * i + 1] > y1) y1 = xy[2 * i + 1];
        }
        if (y1 - y0 > 0 && 0.35 * (x1 - x0) > (y1 - y0)) {
            double f = 0.35 * (x1 - x0) / (y1 - y0);
            for (int i = 0; i < n; i++) xy[2 * i + 1] *= f;
        }
    }
    free(layer);
    return ok;
}

/* ------------------------------------------------------ whole-graph forms */

static void layout_grid(int n, double* xy) {
    int cols = (int)ceil(sqrt((double)n));
    if (cols < 1) cols = 1;
    for (int i = 0; i < n; i++) { xy[2 * i] = i % cols; xy[2 * i + 1] = -(double)(i / cols); }
}

/* Two columns, parts ordered by barycentre. 0 if s is not bipartite. */
static int layout_bipartite(const Sub* s, double* xy) {
    int n = s->n;
    int* col = malloc(sizeof(int) * (size_t)n);
    int* q = malloc(sizeof(int) * (size_t)n);
    int* layer = calloc((size_t)n + 1, sizeof(int));
    int ok = col && q && layer;
    if (ok) {
        for (int i = 0; i < n; i++) col[i] = -1;
        for (int r = 0; r < n && ok; r++) {
            if (col[r] >= 0) continue;
            int h = 0, t = 0;
            col[r] = 0; q[t++] = r;
            while (h < t && ok) {
                int u = q[h++];
                for (int a = s->off[u]; a < s->off[u + 1]; a++) {
                    int w = s->adj[a];
                    if (col[w] < 0) { col[w] = 1 - col[u]; q[t++] = w; }
                    else if (col[w] == col[u]) { ok = 0; break; }
                }
            }
        }
    }
    if (ok) {
        /* Reuse the layered machinery with two layers, then turn it sideways. */
        for (int i = 0; i < n; i++) layer[i] = col[i];
        ok = layout_layers(s, layer, xy);
        if (ok) {
            int nl = 0, nr = 0;
            for (int i = 0; i < n; i++) { if (col[i]) nr++; else nl++; }
            double sep = 0.45 * (nl > nr ? nl : nr);
            if (sep < 1.5) sep = 1.5;
            for (int i = 0; i < n; i++) {
                double x = xy[2 * i];
                xy[2 * i] = col[i] ? sep : 0.0;
                xy[2 * i + 1] = -x;
            }
        }
    }
    free(col); free(q); free(layer);
    return ok;
}

/* -------------------------------------------------------------- packing --- */

typedef struct { int comp; int size; int minv; } CompKey;

static int cmp_comp(const void* a, const void* b) {
    const CompKey* x = (const CompKey*)a; const CompKey* y = (const CompKey*)b;
    if (x->size != y->size) return y->size - x->size;
    return x->minv - y->minv;
}

/* ----------------------------------------------------------- dispatcher --- */

int glayout_parse_method(const Expr* v, GLMethod* out) {
    if (!v) return 0;
    if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic) {
        *out = GL_AUTOMATIC; return 1;
    }
    if (v->type != EXPR_STRING) return 0;
    const char* s = v->data.string;
    static const struct { const char* name; GLMethod m; } T[] = {
        {"CircularEmbedding", GL_CIRCULAR},
        {"SpringElectricalEmbedding", GL_SPRING},
        {"SpringEmbedding", GL_SPRING},
        {"StressEmbedding", GL_STRESS},
        {"LayeredEmbedding", GL_LAYERED},
        {"LayeredDigraphEmbedding", GL_LAYERED},
        {"TreeEmbedding", GL_LAYERED},
        {"BipartiteEmbedding", GL_BIPARTITE},
        {"GridEmbedding", GL_GRID},
    };
    for (size_t i = 0; i < sizeof(T) / sizeof(T[0]); i++)
        if (strcmp(s, T[i].name) == 0) { *out = T[i].m; return 1; }
    return 0;
}

/* Automatic choice: forests with a branch vertex -> layered (tidy trees);
 * all-directed DAGs -> layered; everything else (including paths) -> stress. */
static GLMethod choose_method(const Sub* g) {
    if (g->m == 0) return GL_STRESS;
    int anydir = 0, alldir = 1;
    for (int k = 0; k < g->m; k++) { if (g->dir[k]) anydir = 1; else alldir = 0; }
    if (!anydir) {
        /* forest <=> m == n - #components (no parallel edges in a Graph) */
        int* c = malloc(sizeof(int) * (size_t)g->n);
        int* q = malloc(sizeof(int) * (size_t)g->n);
        GLMethod res = GL_STRESS;
        if (c && q) {
            int nc = 0, maxdeg = 0;
            for (int i = 0; i < g->n; i++) c[i] = -1;
            for (int r = 0; r < g->n; r++) {
                if (c[r] >= 0) continue;
                nc++;
                int h = 0, t = 0;
                c[r] = 1; q[t++] = r;
                while (h < t) {
                    int u = q[h++];
                    for (int a = g->off[u]; a < g->off[u + 1]; a++)
                        if (c[g->adj[a]] < 0) { c[g->adj[a]] = 1; q[t++] = g->adj[a]; }
                }
            }
            for (int i = 0; i < g->n; i++)
                if (g->off[i + 1] - g->off[i] > maxdeg) maxdeg = g->off[i + 1] - g->off[i];
            if (g->m == g->n - nc && maxdeg >= 3) res = GL_LAYERED;
        }
        free(c); free(q);
        return res;
    }
    if (alldir) {
        int* layer = malloc(sizeof(int) * (size_t)g->n);
        int dag = layer ? dag_layers(g, layer) : 0;
        free(layer);
        if (dag) return GL_LAYERED;
    }
    return GL_STRESS;
}

int glayout_compute(int n, int m, const int* eu, const int* ev,
                    const unsigned char* dir, GLMethod method, double* xy) {
    if (n <= 0) return method == GL_AUTOMATIC ? GL_STRESS : (int)method;
    Sub g; memset(&g, 0, sizeof(g));
    g.n = n; g.m = m;
    g.eu = malloc(sizeof(int) * (size_t)(m + 1));
    g.ev = malloc(sizeof(int) * (size_t)(m + 1));
    g.dir = calloc((size_t)m + 1, 1);
    if (!g.eu || !g.ev || !g.dir) { sub_free(&g); return -1; }
    for (int k = 0; k < m; k++) {
        g.eu[k] = eu[k]; g.ev[k] = ev[k]; g.dir[k] = dir ? dir[k] : 0;
    }
    if (!sub_csr(&g)) { sub_free(&g); return -1; }

    if (method == GL_AUTOMATIC) method = choose_method(&g);

    /* Whole-graph embeddings. */
    if (method == GL_CIRCULAR) { layout_circle(n, NULL, xy); sub_free(&g); return method; }
    if (method == GL_GRID)     { layout_grid(n, xy);         sub_free(&g); return method; }
    if (method == GL_BIPARTITE) {
        if (layout_bipartite(&g, xy)) { sub_free(&g); return method; }
        method = GL_STRESS;
    }

    /* Per-component embeddings. */
    int* comp = malloc(sizeof(int) * (size_t)n);
    int* loc = malloc(sizeof(int) * (size_t)n);
    int* q = malloc(sizeof(int) * (size_t)n);
    int* members = malloc(sizeof(int) * (size_t)n);
    int* cstart = calloc((size_t)n + 1, sizeof(int));
    double* cxy = malloc(sizeof(double) * 2 * (size_t)n);
    double* bbox = malloc(sizeof(double) * 4 * (size_t)n);   /* per component */
    CompKey* keys = malloc(sizeof(CompKey) * (size_t)n);
    int ok = comp && loc && q && members && cstart && cxy && bbox && keys;
    int ncomp = 0;
    if (ok) {
        for (int i = 0; i < n; i++) comp[i] = -1;
        for (int r = 0; r < n; r++) {
            if (comp[r] >= 0) continue;
            int h = 0, t = 0;
            comp[r] = ncomp; q[t++] = r;
            while (h < t) {
                int u = q[h++];
                for (int a = g.off[u]; a < g.off[u + 1]; a++)
                    if (comp[g.adj[a]] < 0) { comp[g.adj[a]] = ncomp; q[t++] = g.adj[a]; }
            }
            ncomp++;
        }
        /* members grouped by component, increasing vertex index within */
        for (int i = 0; i < n; i++) cstart[comp[i] + 1]++;
        for (int c = 0; c < ncomp; c++) cstart[c + 1] += cstart[c];
        { int* fill = malloc(sizeof(int) * (size_t)(ncomp + 1));
          if (!fill) ok = 0;
          else {
              memcpy(fill, cstart, sizeof(int) * (size_t)ncomp);
              for (int i = 0; i < n; i++) { loc[i] = fill[comp[i]] - cstart[comp[i]]; members[fill[comp[i]]++] = i; }
              free(fill);
          } }
    }
    /* edges grouped by component */
    int* eorder = ok ? malloc(sizeof(int) * (size_t)(m + 1)) : NULL;
    int* estart = ok ? calloc((size_t)ncomp + 1, sizeof(int)) : NULL;
    if (ok && (!eorder || !estart)) ok = 0;
    if (ok) {
        for (int k = 0; k < m; k++) estart[comp[eu[k]] + 1]++;
        for (int c = 0; c < ncomp; c++) estart[c + 1] += estart[c];
        int* fill = malloc(sizeof(int) * (size_t)(ncomp + 1));
        if (!fill) ok = 0;
        else {
            memcpy(fill, estart, sizeof(int) * (size_t)ncomp);
            for (int k = 0; k < m; k++) eorder[fill[comp[eu[k]]]++] = k;
            free(fill);
        }
    }
    for (int c = 0; ok && c < ncomp; c++) {
        int nc = cstart[c + 1] - cstart[c];
        int mc = estart[c + 1] - estart[c];
        double* L = cxy + 2 * (size_t)cstart[c];
        if (nc == 1) { L[0] = L[1] = 0; }
        else {
            Sub s; memset(&s, 0, sizeof(s));
            s.n = nc; s.m = mc;
            s.eu = malloc(sizeof(int) * (size_t)(mc + 1));
            s.ev = malloc(sizeof(int) * (size_t)(mc + 1));
            s.dir = calloc((size_t)mc + 1, 1);
            if (!s.eu || !s.ev || !s.dir) { sub_free(&s); ok = 0; break; }
            for (int t = 0; t < mc; t++) {
                int k = eorder[estart[c] + t];
                s.eu[t] = loc[eu[k]]; s.ev[t] = loc[ev[k]]; s.dir[t] = g.dir[k];
            }
            if (!sub_csr(&s)) { sub_free(&s); ok = 0; break; }
            int r;
            if (method == GL_LAYERED) r = layout_layered(&s, L);
            else if (method == GL_SPRING) r = layout_spring(&s, L);
            else r = layout_stress(&s, L);
            if (r && method != GL_LAYERED) { canonical_orientation(&s, L); straighten(&s, L); }
            sub_free(&s);
            if (!r) { ok = 0; break; }
        }
        double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
        for (int i = 0; i < nc; i++) {
            if (L[2 * i] < x0) x0 = L[2 * i];
            if (L[2 * i] > x1) x1 = L[2 * i];
            if (L[2 * i + 1] < y0) y0 = L[2 * i + 1];
            if (L[2 * i + 1] > y1) y1 = L[2 * i + 1];
        }
        bbox[4 * c] = x0; bbox[4 * c + 1] = x1; bbox[4 * c + 2] = y0; bbox[4 * c + 3] = y1;
        keys[c].comp = c; keys[c].size = nc; keys[c].minv = members[cstart[c]];
    }
    if (ok) {
        /* Shelf packing, largest component first, rows of bounded width.
         * Isolated vertices are gathered into one square block at the end
         * rather than strung out one per slot. */
        const double gap = 1.0;
        qsort(keys, (size_t)ncomp, sizeof(CompKey), cmp_comp);
        int nsingle = 0;
        for (int t = 0; t < ncomp; t++) if (keys[t].size == 1) nsingle++;
        int nitems = ncomp - nsingle + (nsingle > 1 ? 1 : nsingle);
        int bcols = (int)ceil(sqrt((double)(nsingle > 0 ? nsingle : 1)));
        if (nsingle > 1) {
            /* block: singletons (in key order) on a grid; its bbox is the grid */
            int first = ncomp - nsingle;
            for (int t = first; t < ncomp; t++) {
                int c = keys[t].comp, r = t - first;
                double* L = cxy + 2 * (size_t)cstart[c];
                L[0] = r % bcols; L[1] = -(double)(r / bcols);
            }
        }
        /* item u: a component (keys[u].comp) or, last, the singleton block */
        double* iw = malloc(sizeof(double) * (size_t)(nitems + 1));
        double* ih = malloc(sizeof(double) * (size_t)(nitems + 1));
        double* ix0 = malloc(sizeof(double) * (size_t)(nitems + 1));
        double* iy1 = malloc(sizeof(double) * (size_t)(nitems + 1));
        double* ity = malloc(sizeof(double) * (size_t)(nitems + 1));
        double* itx = malloc(sizeof(double) * (size_t)(nitems + 1));
        if (!iw || !ih || !ix0 || !iy1 || !ity || !itx) ok = 0;
        int blockitem = (nsingle > 1) ? nitems - 1 : -1;
        for (int u = 0; ok && u < nitems; u++) {
            if (u == blockitem) {
                int rows = (nsingle + bcols - 1) / bcols;
                iw[u] = bcols - 1; ih[u] = rows - 1; ix0[u] = 0; iy1[u] = 0;
            } else {
                int c = keys[u].comp;
                iw[u] = bbox[4 * c + 1] - bbox[4 * c]; ih[u] = bbox[4 * c + 3] - bbox[4 * c + 2];
                ix0[u] = bbox[4 * c]; iy1[u] = bbox[4 * c + 3];
            }
        }
        if (ok) {
            double area = 0, maxw = 0;
            for (int u = 0; u < nitems; u++) {
                area += (iw[u] + gap) * (ih[u] + gap);
                if (iw[u] + gap > maxw) maxw = iw[u] + gap;
            }
            double rowmax = sqrt(area) * 1.4;
            if (rowmax < maxw) rowmax = maxw;
            double cx = 0, cy = 0, rowh = 0;
            int top_align = (method == GL_LAYERED);
            int rowbeg = 0;
            for (int u = 0; u <= nitems; u++) {
                if (u == nitems || (cx > 0 && cx + iw[u] > rowmax)) {
                    for (int t = rowbeg; t < u; t++)          /* centre in the row */
                        if (!top_align) ity[t] -= (rowh - ih[t]) / 2.0;
                    cy -= rowh + gap; cx = 0; rowh = 0; rowbeg = u;
                    if (u == nitems) break;
                }
                itx[u] = cx - ix0[u]; ity[u] = cy - iy1[u];
                cx += iw[u] + gap;
                if (ih[u] > rowh) rowh = ih[u];
            }
            for (int t = 0; t < ncomp; t++) {
                int c = keys[t].comp;
                int u = (blockitem >= 0 && t >= ncomp - nsingle) ? blockitem : t;
                double* L = cxy + 2 * (size_t)cstart[c];
                int nc = cstart[c + 1] - cstart[c];
                for (int q2 = 0; q2 < nc; q2++) { L[2 * q2] += itx[u]; L[2 * q2 + 1] += ity[u]; }
            }
        }
        free(iw); free(ih); free(ix0); free(iy1); free(ity); free(itx);
        for (int c = 0; c < ncomp; c++)
            for (int t = cstart[c]; t < cstart[c + 1]; t++) {
                int v = members[t];
                xy[2 * v] = cxy[2 * t]; xy[2 * v + 1] = cxy[2 * t + 1];
            }
    }
    free(comp); free(loc); free(q); free(members); free(cstart); free(cxy);
    free(bbox); free(keys); free(eorder); free(estart);
    sub_free(&g);
    return ok ? (int)method : -1;
}
