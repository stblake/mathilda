/* layout.c - vertex-coordinate layouts for GraphPlot / HighlightGraph.
 *
 * graph_compute_layout() fills caller-allocated x[]/y[] (length n) with 2D
 * coordinates for a validated graph under a named layout, roughly normalized
 * to the [-1, 1] box. Every kernel is fully deterministic (no RNG state, no
 * Date/random) so a notebook re-run reproduces the same picture.
 *
 * We implement a compact set of geometric + force-directed kernels and map the
 * full Wolfram-Language GraphLayout name list onto them (see layout_kernel()).
 * Members of the energy-minimization family (SpringElectrical, Gravity,
 * Spectral, Tutte, Planar, ...) are all served by a Fruchterman-Reingold
 * spring-electrical solver; tree/layer families share a BFS-layered kernel;
 * partite families share a two-column kernel. Names we approximate rather than
 * reproduce exactly are documented in docs/spec/builtins/graphs.md.
 *
 * Memory: read-only over g; allocates a transient GraphAdj for the kernels that
 * need adjacency and frees it before returning.
 */

#include "graph.h"
#include "expr.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum {
    LAYOUT_CIRCULAR,
    LAYOUT_SPRING,      /* Fruchterman-Reingold spring-electrical             */
    LAYOUT_GRAVITY,     /* FR with an added central gravity well              */
    LAYOUT_HIGHDIM,     /* two-pivot BFS-distance embedding projected to 2D   */
    LAYOUT_HYPERBOLIC,  /* FR then radial warp onto a Poincare-like disk      */
    LAYOUT_SPIRAL,
    LAYOUT_LINEAR,
    LAYOUT_GRID,
    LAYOUT_RANDOM,
    LAYOUT_STAR,
    LAYOUT_RADIAL,      /* concentric BFS shells from a root                  */
    LAYOUT_LAYERED,     /* stacked BFS layers                                 */
    LAYOUT_BIPARTITE    /* two columns from a BFS 2-coloring                  */
} LayoutKernel;

/* Case-sensitive exact match against the Wolfram GraphLayout name list, mapping
 * each to the kernel that best serves it. Unknown / None -> circular. */
static LayoutKernel layout_kernel(const char* name) {
    if (!name || !*name) return LAYOUT_CIRCULAR;
    struct { const char* n; LayoutKernel k; } table[] = {
        { "CircularEmbedding",            LAYOUT_CIRCULAR  },
        { "SpringElectricalEmbedding",    LAYOUT_SPRING    },
        { "SpringEmbedding",              LAYOUT_SPRING    },
        { "GravityEmbedding",             LAYOUT_GRAVITY   },
        { "HighDimensionalEmbedding",     LAYOUT_HIGHDIM   },
        { "SpectralEmbedding",            LAYOUT_HIGHDIM   },
        { "SphericalEmbedding",           LAYOUT_HYPERBOLIC},
        { "HyperbolicSpringEmbedding",    LAYOUT_HYPERBOLIC},
        { "TutteEmbedding",               LAYOUT_SPRING    },
        { "PlanarEmbedding",              LAYOUT_SPRING    },
        { "SpiralEmbedding",              LAYOUT_SPIRAL    },
        { "DiscreteSpiralEmbedding",      LAYOUT_SPIRAL    },
        { "LinearEmbedding",              LAYOUT_LINEAR    },
        { "GridEmbedding",                LAYOUT_GRID      },
        { "RandomEmbedding",              LAYOUT_RANDOM    },
        { "StarEmbedding",                LAYOUT_STAR      },
        { "RadialEmbedding",              LAYOUT_RADIAL    },
        { "BalloonEmbedding",             LAYOUT_RADIAL    },
        { "HyperbolicRadialEmbedding",    LAYOUT_RADIAL    },
        { "LayeredEmbedding",             LAYOUT_LAYERED   },
        { "LayeredDigraphEmbedding",      LAYOUT_LAYERED   },
        { "SymmetricLayeredEmbedding",    LAYOUT_LAYERED   },
        { "BipartiteEmbedding",           LAYOUT_BIPARTITE },
        { "MultipartiteEmbedding",        LAYOUT_BIPARTITE },
        { "CircularMultipartiteEmbedding",LAYOUT_BIPARTITE },
    };
    for (size_t i = 0; i < sizeof(table) / sizeof(table[0]); i++)
        if (strcmp(name, table[i].n) == 0) return table[i].k;
    return LAYOUT_CIRCULAR;
}

/* Rescale coordinates to fit the [-1, 1] box while preserving aspect ratio. */
static void normalize(double* x, double* y, int n) {
    if (n <= 0) return;
    double minx = x[0], maxx = x[0], miny = y[0], maxy = y[0];
    for (int i = 1; i < n; i++) {
        if (x[i] < minx) minx = x[i];
        if (x[i] > maxx) maxx = x[i];
        if (y[i] < miny) miny = y[i];
        if (y[i] > maxy) maxy = y[i];
    }
    double cx = 0.5 * (minx + maxx), cy = 0.5 * (miny + maxy);
    double span = fmax(maxx - minx, maxy - miny);
    double s = (span > 1e-9) ? 2.0 / span : 1.0;   /* map larger span to [-1,1] */
    for (int i = 0; i < n; i++) {
        x[i] = (x[i] - cx) * s;
        y[i] = (y[i] - cy) * s;
    }
}

static void layout_circular(double* x, double* y, int n) {
    for (int i = 0; i < n; i++) {
        if (n == 1) { x[i] = 0.0; y[i] = 0.0; continue; }
        double t = 2.0 * M_PI * (double)i / (double)n;
        x[i] = cos(t); y[i] = sin(t);
    }
}

static void layout_linear(double* x, double* y, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = (n == 1) ? 0.0 : (-1.0 + 2.0 * (double)i / (double)(n - 1));
        y[i] = 0.0;
    }
}

static void layout_spiral(double* x, double* y, int n) {
    /* Archimedean spiral: radius grows linearly with a constant angular step. */
    for (int i = 0; i < n; i++) {
        double t = (double)i * 0.6;          /* angular step (radians)         */
        double r = (n <= 1) ? 0.0 : (double)i / (double)(n - 1);
        x[i] = r * cos(t); y[i] = r * sin(t);
    }
}

static void layout_grid(double* x, double* y, int n) {
    int cols = (int)ceil(sqrt((double)n));
    if (cols < 1) cols = 1;
    for (int i = 0; i < n; i++) {
        int r = i / cols, c = i % cols;
        x[i] = (double)c; y[i] = -(double)r;   /* row 0 on top                 */
    }
}

static void layout_random(double* x, double* y, int n) {
    /* Deterministic pseudo-random: a fixed LCG folded over the vertex index,
     * so the picture is stable across runs (no global RNG dependence). */
    unsigned long s = 0x9E3779B9UL;
    for (int i = 0; i < n; i++) {
        s = s * 1103515245UL + 12345UL + (unsigned long)(i + 1) * 2654435761UL;
        double ux = (double)((s >> 16) & 0x7FFF) / 32767.0;
        s = s * 1103515245UL + 12345UL;
        double uy = (double)((s >> 16) & 0x7FFF) / 32767.0;
        x[i] = 2.0 * ux - 1.0; y[i] = 2.0 * uy - 1.0;
    }
}

static void layout_star(const GraphAdj* a, double* x, double* y, int n) {
    /* Highest-degree vertex at the center; the rest evenly on a circle. */
    int hub = 0, best = -1;
    for (int i = 0; i < n; i++) {
        int deg = a ? (a->outdeg[i] + a->indeg[i]) : 0;
        if (deg > best) { best = deg; hub = i; }
    }
    int placed = 0, rim = (n > 1) ? n - 1 : 1;
    for (int i = 0; i < n; i++) {
        if (i == hub) { x[i] = 0.0; y[i] = 0.0; continue; }
        double t = 2.0 * M_PI * (double)placed / (double)rim;
        x[i] = cos(t); y[i] = sin(t); placed++;
    }
}

/* BFS distance labels from the highest-degree vertex, over the underlying
 * undirected graph. dist[i] = shell index; unreached vertices get maxd+1. */
static int* bfs_levels(const GraphAdj* a, int n, int* maxd_out) {
    int* dist = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* queue = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!dist || !queue) { free(dist); free(queue); *maxd_out = 0; return NULL; }
    for (int i = 0; i < n; i++) dist[i] = -1;
    int root = 0, best = -1;
    for (int i = 0; i < n; i++) {
        int deg = a->outdeg[i] + a->indeg[i];
        if (deg > best) { best = deg; root = i; }
    }
    int head = 0, tail = 0, maxd = 0;
    dist[root] = 0; queue[tail++] = root;
    while (head < tail) {
        int u = queue[head++];
        if (dist[u] > maxd) maxd = dist[u];
        for (int e = 0; e < a->outdeg[u]; e++) {
            int v = a->out[u][e];
            if (dist[v] < 0) { dist[v] = dist[u] + 1; queue[tail++] = v; }
        }
        for (int e = 0; e < a->indeg[u]; e++) {
            int v = a->in[u][e];
            if (dist[v] < 0) { dist[v] = dist[u] + 1; queue[tail++] = v; }
        }
    }
    for (int i = 0; i < n; i++) if (dist[i] < 0) dist[i] = maxd + 1; /* disconnected */
    int any_disc = 0;
    for (int i = 0; i < n; i++) if (dist[i] == maxd + 1) { any_disc = 1; break; }
    if (any_disc) maxd += 1;
    free(queue);
    *maxd_out = maxd;
    return dist;
}

static void layout_radial(const GraphAdj* a, double* x, double* y, int n) {
    int maxd = 0;
    int* dist = bfs_levels(a, n, &maxd);
    if (!dist) { layout_circular(x, y, n); return; }
    int* count = calloc((size_t)(maxd + 2), sizeof(int));
    int* seen  = calloc((size_t)(maxd + 2), sizeof(int));
    for (int i = 0; i < n; i++) count[dist[i]]++;
    for (int i = 0; i < n; i++) {
        int d = dist[i];
        double r = (maxd == 0) ? 0.0 : (double)d / (double)maxd;
        int c = count[d];
        double t = (c > 0) ? 2.0 * M_PI * (double)seen[d] / (double)c : 0.0;
        seen[d]++;
        x[i] = r * cos(t); y[i] = r * sin(t);
    }
    free(count); free(seen); free(dist);
}

static void layout_layered(const GraphAdj* a, double* x, double* y, int n) {
    int maxd = 0;
    int* dist = bfs_levels(a, n, &maxd);
    if (!dist) { layout_circular(x, y, n); return; }
    int* count = calloc((size_t)(maxd + 2), sizeof(int));
    int* seen  = calloc((size_t)(maxd + 2), sizeof(int));
    for (int i = 0; i < n; i++) count[dist[i]]++;
    for (int i = 0; i < n; i++) {
        int d = dist[i], c = count[d];
        x[i] = (c <= 1) ? 0.0 : (-1.0 + 2.0 * (double)seen[d] / (double)(c - 1));
        y[i] = -(double)d;
        seen[d]++;
    }
    free(count); free(seen); free(dist);
}

static void layout_bipartite(const GraphAdj* a, double* x, double* y, int n) {
    /* 2-color the underlying undirected graph by BFS parity; place color 0 in
     * the left column and color 1 in the right. Non-bipartite graphs still get
     * a legible two-column split (parity is well defined per BFS tree). */
    int* color = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* queue = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!color || !queue) { free(color); free(queue); layout_circular(x, y, n); return; }
    for (int i = 0; i < n; i++) color[i] = -1;
    for (int s = 0; s < n; s++) {
        if (color[s] >= 0) continue;
        int head = 0, tail = 0;
        color[s] = 0; queue[tail++] = s;
        while (head < tail) {
            int u = queue[head++];
            for (int e = 0; e < a->outdeg[u]; e++) {
                int v = a->out[u][e];
                if (color[v] < 0) { color[v] = 1 - color[u]; queue[tail++] = v; }
            }
            for (int e = 0; e < a->indeg[u]; e++) {
                int v = a->in[u][e];
                if (color[v] < 0) { color[v] = 1 - color[u]; queue[tail++] = v; }
            }
        }
    }
    int c0 = 0, c1 = 0;
    for (int i = 0; i < n; i++) { if (color[i] == 0) c0++; else c1++; }
    int s0 = 0, s1 = 0;
    for (int i = 0; i < n; i++) {
        if (color[i] == 0) {
            x[i] = -1.0;
            y[i] = (c0 <= 1) ? 0.0 : (-1.0 + 2.0 * (double)s0++ / (double)(c0 - 1));
        } else {
            x[i] = 1.0;
            y[i] = (c1 <= 1) ? 0.0 : (-1.0 + 2.0 * (double)s1++ / (double)(c1 - 1));
        }
    }
    free(color); free(queue);
}

/* Fruchterman-Reingold spring-electrical core. Vertices repel each other
 * (~k^2/d); adjacent vertices attract (~d^2/k). When `gravity` > 0 an extra
 * central pull (~gravity*distance-from-origin, scaled by 1+degree) compacts the
 * graph and draws well-connected vertices inward — this is what distinguishes
 * GravityEmbedding from the plain spring-electrical layout. Initialized from the
 * circular layout (deterministic) and cooled linearly, so the result is stable.
 */
static void fr_core(const GraphAdj* a, double* x, double* y, int n, double gravity) {
    if (n <= 2) { layout_circular(x, y, n); return; }
    layout_circular(x, y, n);                     /* deterministic seed         */
    double area = 4.0;                            /* the [-1,1]^2 box           */
    double k = sqrt(area / (double)n);            /* natural edge length        */
    double* dx = calloc((size_t)n, sizeof(double));
    double* dy = calloc((size_t)n, sizeof(double));
    if (!dx || !dy) { free(dx); free(dy); return; }
    int iters = 300;
    double t = 0.20;                              /* initial max displacement   */
    double cool = t / (double)(iters + 1);
    for (int it = 0; it < iters; it++) {
        for (int i = 0; i < n; i++) { dx[i] = 0.0; dy[i] = 0.0; }
        /* Repulsion between every pair. */
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double ddx = x[i] - x[j], ddy = y[i] - y[j];
                double d2 = ddx * ddx + ddy * ddy;
                if (d2 < 1e-9) { ddx = 1e-3 * (i - j); ddy = 1e-3; d2 = ddx*ddx+ddy*ddy; }
                double d = sqrt(d2);
                double f = (k * k) / d;
                double ux = ddx / d, uy = ddy / d;
                dx[i] += ux * f; dy[i] += uy * f;
                dx[j] -= ux * f; dy[j] -= uy * f;
            }
        }
        /* Attraction along edges (undirected: use out-adjacency, each edge once
         * per stored direction — symmetric contributions balance out). */
        for (int u = 0; u < n; u++) {
            for (int e = 0; e < a->outdeg[u]; e++) {
                int v = a->out[u][e];
                double ddx = x[u] - x[v], ddy = y[u] - y[v];
                double d = sqrt(ddx * ddx + ddy * ddy);
                if (d < 1e-9) continue;
                double f = (d * d) / k;
                double ux = ddx / d, uy = ddy / d;
                dx[u] -= ux * f; dy[u] -= uy * f;
                dx[v] += ux * f; dy[v] += uy * f;
            }
        }
        /* Optional central gravity: pull toward the origin, stronger for
         * higher-degree vertices, so hubs settle in the middle. */
        if (gravity > 0.0) {
            for (int i = 0; i < n; i++) {
                double mass = 1.0 + 0.5 * (double)(a->outdeg[i] + a->indeg[i]);
                dx[i] -= gravity * mass * x[i];
                dy[i] -= gravity * mass * y[i];
            }
        }
        /* Displace, capped by the current temperature. */
        for (int i = 0; i < n; i++) {
            double d = sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
            if (d > 1e-9) {
                double cap = fmin(d, t);
                x[i] += (dx[i] / d) * cap;
                y[i] += (dy[i] / d) * cap;
            }
        }
        t -= cool;
    }
    free(dx); free(dy);
}

static void layout_spring(const GraphAdj* a, double* x, double* y, int n) {
    fr_core(a, x, y, n, 0.0);
}

static void layout_gravity(const GraphAdj* a, double* x, double* y, int n) {
    fr_core(a, x, y, n, 0.12);
}

/* HyperbolicSpring / Spherical: run the spring layout, then warp radii outward
 * so vertices crowd toward the boundary of a disk (a Poincare-disk feel). */
static void layout_hyperbolic(const GraphAdj* a, double* x, double* y, int n) {
    fr_core(a, x, y, n, 0.04);
    if (n <= 1) return;
    double cx = 0.0, cy = 0.0;
    for (int i = 0; i < n; i++) { cx += x[i]; cy += y[i]; }
    cx /= n; cy /= n;
    double maxr = 1e-9;
    for (int i = 0; i < n; i++) {
        double r = hypot(x[i] - cx, y[i] - cy);
        if (r > maxr) maxr = r;
    }
    for (int i = 0; i < n; i++) {
        double ddx = x[i] - cx, ddy = y[i] - cy;
        double r = hypot(ddx, ddy);
        if (r < 1e-9) continue;
        double rn = r / maxr;                 /* in [0,1]                       */
        double warped = pow(rn, 0.45);        /* push outward toward the rim    */
        double s = warped / rn;
        x[i] = cx + ddx * s;
        y[i] = cy + ddy * s;
    }
}

/* BFS distances from `src` over the underlying undirected graph; unreached
 * vertices get `n` (a value larger than any real distance). */
static void bfs_from(const GraphAdj* a, int n, int src, int* dist) {
    int* queue = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!queue) { for (int i = 0; i < n; i++) dist[i] = 0; return; }
    for (int i = 0; i < n; i++) dist[i] = -1;
    int head = 0, tail = 0;
    dist[src] = 0; queue[tail++] = src;
    while (head < tail) {
        int u = queue[head++];
        for (int e = 0; e < a->outdeg[u]; e++) {
            int v = a->out[u][e];
            if (dist[v] < 0) { dist[v] = dist[u] + 1; queue[tail++] = v; }
        }
        for (int e = 0; e < a->indeg[u]; e++) {
            int v = a->in[u][e];
            if (dist[v] < 0) { dist[v] = dist[u] + 1; queue[tail++] = v; }
        }
    }
    for (int i = 0; i < n; i++) if (dist[i] < 0) dist[i] = n;
    free(queue);
}

static int farthest(const int* dist, int n) {
    int best = 0, bd = -1;
    for (int i = 0; i < n; i++) {
        int d = (dist[i] >= n) ? -1 : dist[i];   /* ignore unreached           */
        if (d > bd) { bd = d; best = i; }
    }
    return best;
}

/* HighDimensional / Spectral (approx): embed each vertex by its BFS distance to
 * two far-apart pivots (found by a double sweep) and use those two distances as
 * the 2D coordinates — a cheap pivot-MDS that lays the graph out along its
 * "diameter". A tiny index-based nudge separates coincident vertices. */
static void layout_highdim(const GraphAdj* a, double* x, double* y, int n) {
    if (n <= 2) { layout_circular(x, y, n); return; }
    int* d = malloc((size_t)n * sizeof(int));
    int* d1 = malloc((size_t)n * sizeof(int));
    int* d2 = malloc((size_t)n * sizeof(int));
    if (!d || !d1 || !d2) { free(d); free(d1); free(d2); layout_circular(x, y, n); return; }
    bfs_from(a, n, 0, d);          int p1 = farthest(d, n);
    bfs_from(a, n, p1, d1);        int p2 = farthest(d1, n);
    bfs_from(a, n, p2, d2);
    for (int i = 0; i < n; i++) {
        double a1 = (d1[i] >= n) ? 0.0 : (double)d1[i];
        double a2 = (d2[i] >= n) ? 0.0 : (double)d2[i];
        x[i] = a1 + 0.001 * (double)(i % 7);
        y[i] = a2 + 0.001 * (double)(i % 5);
    }
    free(d); free(d1); free(d2);
}

void graph_compute_layout(const Expr* g, const char* layout, double* x, double* y) {
    const Expr* verts = g->data.function.args[0];
    int n = (int)verts->data.function.arg_count;
    if (n <= 0) return;

    LayoutKernel kern = layout_kernel(layout);

    /* Kernels that need adjacency build it once; geometric kernels don't. */
    GraphAdj* a = NULL;
    if (kern == LAYOUT_SPRING || kern == LAYOUT_GRAVITY || kern == LAYOUT_HIGHDIM
        || kern == LAYOUT_HYPERBOLIC || kern == LAYOUT_STAR || kern == LAYOUT_RADIAL
        || kern == LAYOUT_LAYERED || kern == LAYOUT_BIPARTITE) {
        a = graph_build_adj(g);
        if (!a) kern = LAYOUT_CIRCULAR;   /* defensive: fall back if build fails */
    }

    switch (kern) {
        case LAYOUT_CIRCULAR:  layout_circular(x, y, n);       break;
        case LAYOUT_LINEAR:    layout_linear(x, y, n);         break;
        case LAYOUT_SPIRAL:    layout_spiral(x, y, n);         break;
        case LAYOUT_GRID:      layout_grid(x, y, n);           break;
        case LAYOUT_RANDOM:    layout_random(x, y, n);         break;
        case LAYOUT_STAR:      layout_star(a, x, y, n);        break;
        case LAYOUT_RADIAL:    layout_radial(a, x, y, n);      break;
        case LAYOUT_LAYERED:   layout_layered(a, x, y, n);     break;
        case LAYOUT_BIPARTITE: layout_bipartite(a, x, y, n);   break;
        case LAYOUT_SPRING:    layout_spring(a, x, y, n);      break;
        case LAYOUT_GRAVITY:   layout_gravity(a, x, y, n);     break;
        case LAYOUT_HIGHDIM:   layout_highdim(a, x, y, n);     break;
        case LAYOUT_HYPERBOLIC:layout_hyperbolic(a, x, y, n);  break;
    }

    if (a) graph_adj_free(a);
    normalize(x, y, n);
}

/* ---- 3D layout (Graph3D) --------------------------------------------------- */

static void normalize3d(double* x, double* y, double* z, int n) {
    if (n <= 0) return;
    double mn[3] = { x[0], y[0], z[0] }, mx[3] = { x[0], y[0], z[0] };
    for (int i = 1; i < n; i++) {
        if (x[i] < mn[0]) mn[0] = x[i]; if (x[i] > mx[0]) mx[0] = x[i];
        if (y[i] < mn[1]) mn[1] = y[i]; if (y[i] > mx[1]) mx[1] = y[i];
        if (z[i] < mn[2]) mn[2] = z[i]; if (z[i] > mx[2]) mx[2] = z[i];
    }
    double span = fmax(mx[0] - mn[0], fmax(mx[1] - mn[1], mx[2] - mn[2]));
    double s = (span > 1e-9) ? 2.0 / span : 1.0;
    double cx = 0.5*(mn[0]+mx[0]), cy = 0.5*(mn[1]+mx[1]), cz = 0.5*(mn[2]+mx[2]);
    for (int i = 0; i < n; i++) {
        x[i] = (x[i]-cx)*s; y[i] = (y[i]-cy)*s; z[i] = (z[i]-cz)*s;
    }
}

/* Deterministic near-uniform points on the unit sphere (golden-angle spiral). */
static void fibonacci_sphere(double* x, double* y, double* z, int n) {
    double golden = M_PI * (3.0 - sqrt(5.0));   /* golden angle                 */
    for (int i = 0; i < n; i++) {
        double zt = (n == 1) ? 0.0 : 1.0 - 2.0 * (double)i / (double)(n - 1);
        double r = sqrt(fmax(0.0, 1.0 - zt * zt));
        double th = golden * (double)i;
        x[i] = r * cos(th); y[i] = r * sin(th); z[i] = zt;
    }
}

/* 3D Fruchterman-Reingold, seeded on a sphere (deterministic). */
static void fr3d(const GraphAdj* a, double* x, double* y, double* z, int n) {
    fibonacci_sphere(x, y, z, n);
    if (n <= 2) return;
    double k = pow(8.0 / (double)n, 1.0 / 3.0);    /* natural length in a cube   */
    double* dx = calloc((size_t)n, sizeof(double));
    double* dy = calloc((size_t)n, sizeof(double));
    double* dz = calloc((size_t)n, sizeof(double));
    if (!dx || !dy || !dz) { free(dx); free(dy); free(dz); return; }
    int iters = 250;
    double t = 0.20, cool = t / (double)(iters + 1);
    for (int it = 0; it < iters; it++) {
        for (int i = 0; i < n; i++) { dx[i]=dy[i]=dz[i]=0.0; }
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double ax = x[i]-x[j], ay = y[i]-y[j], az = z[i]-z[j];
                double d2 = ax*ax + ay*ay + az*az;
                if (d2 < 1e-9) { ax = 1e-3*(i-j); ay = 1e-3; az = 1e-3*(j-i); d2 = ax*ax+ay*ay+az*az; }
                double d = sqrt(d2), f = (k*k)/d;
                double ux=ax/d, uy=ay/d, uz=az/d;
                dx[i]+=ux*f; dy[i]+=uy*f; dz[i]+=uz*f;
                dx[j]-=ux*f; dy[j]-=uy*f; dz[j]-=uz*f;
            }
        }
        for (int u = 0; u < n; u++) {
            for (int e = 0; e < a->outdeg[u]; e++) {
                int v = a->out[u][e];
                double ax = x[u]-x[v], ay = y[u]-y[v], az = z[u]-z[v];
                double d = sqrt(ax*ax+ay*ay+az*az);
                if (d < 1e-9) continue;
                double f = (d*d)/k, ux=ax/d, uy=ay/d, uz=az/d;
                dx[u]-=ux*f; dy[u]-=uy*f; dz[u]-=uz*f;
                dx[v]+=ux*f; dy[v]+=uy*f; dz[v]+=uz*f;
            }
        }
        for (int i = 0; i < n; i++) {
            double d = sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);
            if (d > 1e-9) {
                double cap = fmin(d, t);
                x[i]+=(dx[i]/d)*cap; y[i]+=(dy[i]/d)*cap; z[i]+=(dz[i]/d)*cap;
            }
        }
        t -= cool;
    }
    free(dx); free(dy); free(dz);
}

void graph_compute_layout3d(const Expr* g, const char* layout,
                            double* x, double* y, double* z) {
    const Expr* verts = g->data.function.args[0];
    int n = (int)verts->data.function.arg_count;
    if (n <= 0) return;
    int spherical = (layout && strcmp(layout, "SphericalEmbedding") == 0);
    GraphAdj* a = graph_build_adj(g);
    if (!a || spherical) {
        fibonacci_sphere(x, y, z, n);
    } else {
        fr3d(a, x, y, z, n);
    }
    if (a) graph_adj_free(a);
    normalize3d(x, y, z, n);
}
