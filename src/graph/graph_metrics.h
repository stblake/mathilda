#ifndef GRAPH_METRICS_H
#define GRAPH_METRICS_H

/* graph_metrics.h - distance, centrality, clustering and generator builtins.
 *
 * Implementation lives in src/graph/gmet_*.c:
 *   gmet_core.c        shared kernels: CSR adjacency, weights, the thread
 *                      team, bit-parallel multi-source BFS, Dijkstra, packed
 *                      result builders, and a per-graph result cache
 *   gmet_distance.c    GraphDistanceMatrix, GraphDistance[g, s], eccentricity
 *                      and the diameter/radius/center/periphery family,
 *                      MeanGraphDistance, GraphDensity, KirchhoffMatrix
 *   gmet_centrality.c  Degree/Closeness/Eccentricity/Betweenness/
 *                      EdgeBetweenness centralities
 *   gmet_spectral.c    PageRank/Eigenvector/Katz/HITS centralities
 *   gmet_cluster.c     GraphTriangleCount and the clustering coefficients
 *   gmet_generators.c  WheelGraph, HypercubeGraph, GridGraph, KaryTree,
 *                      CompleteKaryTree, CirculantGraph, PetersenGraph,
 *                      TuranGraph, HararyGraph, CompleteGraph[{n1, n2, ...}]
 *   gmet_init.c        graph_metrics_init: registration, attributes, docs
 *
 * Every builtin follows the SPEC section 4 contract (NULL = unevaluated). The
 * internal kernels below operate on plain C arrays only, so they are safe to
 * run on worker threads; nothing in a worker touches an Expr.
 */

#include "expr.h"
#include <stdint.h>
#include <stddef.h>

/* ---- Builtins ------------------------------------------------------------ */
Expr* builtin_graph_distance_matrix(Expr* res);   /* GraphDistanceMatrix     */
Expr* builtin_gmet_graph_distance(Expr* res);     /* GraphDistance (all forms) */
Expr* builtin_vertex_eccentricity(Expr* res);     /* VertexEccentricity      */
Expr* builtin_graph_diameter(Expr* res);          /* GraphDiameter           */
Expr* builtin_graph_radius(Expr* res);            /* GraphRadius             */
Expr* builtin_graph_center(Expr* res);            /* GraphCenter             */
Expr* builtin_graph_periphery(Expr* res);         /* GraphPeriphery          */
Expr* builtin_mean_graph_distance(Expr* res);     /* MeanGraphDistance       */
Expr* builtin_graph_density(Expr* res);           /* GraphDensity            */
Expr* builtin_kirchhoff_matrix(Expr* res);        /* KirchhoffMatrix         */

Expr* builtin_degree_centrality(Expr* res);       /* DegreeCentrality        */
Expr* builtin_closeness_centrality(Expr* res);    /* ClosenessCentrality     */
Expr* builtin_eccentricity_centrality(Expr* res); /* EccentricityCentrality  */
Expr* builtin_betweenness_centrality(Expr* res);  /* BetweennessCentrality   */
Expr* builtin_edge_betweenness_centrality(Expr* res); /* EdgeBetweenness...  */
Expr* builtin_pagerank_centrality(Expr* res);     /* PageRankCentrality      */
Expr* builtin_eigenvector_centrality(Expr* res);  /* EigenvectorCentrality   */
Expr* builtin_katz_centrality(Expr* res);         /* KatzCentrality          */
Expr* builtin_hits_centrality(Expr* res);         /* HITSCentrality          */

Expr* builtin_graph_triangle_count(Expr* res);    /* GraphTriangleCount      */
Expr* builtin_local_clustering_coefficient(Expr* res);
Expr* builtin_global_clustering_coefficient(Expr* res);
Expr* builtin_mean_clustering_coefficient(Expr* res);

Expr* builtin_wheel_graph(Expr* res);             /* WheelGraph[n]           */
Expr* builtin_hypercube_graph(Expr* res);         /* HypercubeGraph[n]       */
Expr* builtin_grid_graph(Expr* res);              /* GridGraph[{n1, ...}]    */
Expr* builtin_kary_tree(Expr* res);               /* KaryTree[n] / [n, k]    */
Expr* builtin_complete_kary_tree(Expr* res);      /* CompleteKaryTree[n(,k)] */
Expr* builtin_circulant_graph(Expr* res);         /* CirculantGraph[n, js]   */
Expr* builtin_petersen_graph(Expr* res);          /* PetersenGraph[(n, k)]   */
Expr* builtin_turan_graph(Expr* res);             /* TuranGraph[n, k]        */
Expr* builtin_harary_graph(Expr* res);            /* HararyGraph[k, n]       */
Expr* builtin_gmet_complete_graph(Expr* res);     /* CompleteGraph (all forms) */

/* Registers every builtin above. Called once, last, from graph_init(). */
void graph_metrics_init(void);

/* ---- Internal kernels (gmet_core.c) --------------------------------------- */

/* Compressed sparse rows over vertex indices 0..n-1. Row v lists the arcs
 * leaving v: adj[off[v] .. off[v+1]-1]. `w` (arc weight) and `eid` (index of
 * the source edge in EdgeList order) are parallel to adj and may be NULL. */
typedef struct {
    int      n;
    int64_t  narcs;
    int64_t* off;       /* n + 1 entries */
    int*     adj;
    double*  w;         /* NULL unless requested */
    int*     eid;       /* NULL unless requested */
} GmetCSR;

/* Which arcs an edge contributes. An UndirectedEdge always contributes both
 * orientations; a DirectedEdge u->v contributes u->v (OUT), v->u (IN), or
 * both (UND, the underlying undirected graph). */
enum { GMET_OUT = 0, GMET_IN = 1, GMET_UND = 2 };

/* Kind of a valid graph, from its directed-edge count. */
enum { GMET_KIND_UNDIRECTED = 0, GMET_KIND_DIRECTED = 1, GMET_KIND_MIXED = 2 };

/* Graph kind, or -1 when g is not a valid graph. An edgeless graph is
 * undirected. */
int gmet_graph_kind(const Expr* g);

/* Edge-weight resolution. Returns 0 for an unweighted graph (*w_out = NULL),
 * 1 when every weight is a non-negative real number (*w_out = malloc'd array of
 * edge weights as doubles, EdgeList order), and -1 when the graph carries a
 * weight that is symbolic, complex or negative (the caller leaves the call
 * unevaluated, as Wolfram does). */
int gmet_edge_weights(const Expr* g, double** w_out);

/* Build the CSR of g in the given mode. `ew` (edge weights, may be NULL) fills
 * csr->w; want_eid fills csr->eid. Arcs of a row keep EdgeList order. Returns
 * 0 if g is not valid (or allocation failed); csr is then zeroed. */
int  gmet_csr_build(const Expr* g, int mode, const double* ew, int want_eid,
                    GmetCSR* csr);
void gmet_csr_free(GmetCSR* csr);

/* ---- Thread team ----------------------------------------------------------
 * gmet_parallel_for runs fn(ctx, tid, lo, hi) over [0, n) in chunks of `chunk`
 * items handed out dynamically to `nthreads` workers (tid in 0..nthreads-1).
 * nthreads <= 1, or a build without MATHILDA_THREADS, runs serially on the
 * calling thread with tid 0. */
typedef void (*gmet_work_fn)(void* ctx, int tid, int64_t lo, int64_t hi);
int  gmet_thread_count(double work);   /* threads worth using for `work` ops */
/* nt reduced so that nt * per_thread_bytes stays under a fixed scratch budget
 * (GMET_SCRATCH_BUDGET); never below 1. */
#define GMET_SCRATCH_BUDGET (512.0 * 1024 * 1024)
int  gmet_cap_threads(int nt, double per_thread_bytes);
void gmet_parallel_for(int64_t n, int64_t chunk, int nthreads,
                       gmet_work_fn fn, void* ctx);

/* ---- Multi-source BFS -----------------------------------------------------
 * Bit-parallel BFS (MS-BFS): up to GMET_MSBFS_LANES sources advance together,
 * each vertex carrying one bit per source, so the adjacency is walked once per
 * level for the whole batch instead of once per source. `visit` is called for
 * every (level, vertex, lane-mask) newly reached, level >= 1; lane i is
 * sources[i]. The sources themselves are level 0 and are not reported. */
#define GMET_MSBFS_WORDS 4
#if defined(__GNUC__) || defined(__clang__)
#define GMET_CTZ64(x) __builtin_ctzll(x)
#else
static inline int gmet_ctz64_fallback(uint64_t x) {
    int k = 0; while (!(x & 1)) { x >>= 1; k++; } return k;
}
#define GMET_CTZ64(x) gmet_ctz64_fallback(x)
#endif
#define GMET_MSBFS_LANES (64 * GMET_MSBFS_WORDS)
typedef struct GmetMSBFS GmetMSBFS;
typedef void (*gmet_msbfs_visit)(void* ctx, int level, int v, const uint64_t* lanes);
GmetMSBFS* gmet_msbfs_new(int n);
void       gmet_msbfs_free(GmetMSBFS* s);
void       gmet_msbfs_run(GmetMSBFS* s, const GmetCSR* g, const int* sources,
                          int nsrc, int maxlevel, gmet_msbfs_visit visit, void* ctx);

/* Single-source Dijkstra over non-negative arc weights, binary heap. dist[] is
 * filled (INFINITY = unreachable). order[] (optional) receives the vertices in
 * settling order and *nsettled their count. Scratch is caller-provided via
 * gmet_dijkstra_scratch so repeated calls do not allocate. */
typedef struct GmetHeap GmetHeap;
GmetHeap* gmet_heap_new(int64_t cap);
void      gmet_heap_free(GmetHeap* h);
void      gmet_dijkstra(const GmetCSR* g, int src, double* dist, int* order,
                        int* nsettled, GmetHeap* h);

/* Per-source distance summary over the OUT arcs: for every vertex s, the
 * number of other vertices reachable from s, the sum of their distances, and
 * the largest of them (0 when none). Weighted when g->w != NULL. */
typedef struct {
    int      n;
    int      weighted;
    int64_t* reach;     /* reachable vertices other than s                  */
    double*  sum;       /* sum of distances to them (exact for BFS < 2^53)   */
    double*  ecc;       /* max distance to a reachable vertex                */
} GmetDistSummary;

/* Summary for g in OUT mode, cached per graph node (see gmet_core.c). The
 * returned pointer is owned by the cache and stays valid until the next call
 * that could evict it; copy what must outlive that. NULL when g is invalid or
 * carries unusable weights (*status = -1) -- *status is 0 on success. */
const GmetDistSummary* gmet_dist_summary(const Expr* g, int* status);

/* Strongly connected components (Tarjan, iterative) of the OUT arcs of csr.
 * comp[v] receives a component id in 0..k-1; returns k. */
int gmet_scc(const GmetCSR* csr, int* comp);

/* ---- Result builders ------------------------------------------------------ */
/* Packed (when above the packing threshold) vector of doubles / int64s. The
 * data is copied; the caller keeps ownership of `v`. */
Expr* gmet_real_vector(const double* v, int64_t n);
Expr* gmet_int_vector(const int64_t* v, int64_t n);
/* List of the vertices whose flag is set, in VertexList order. */
Expr* gmet_vertex_subset(const Expr* g, const unsigned char* flag);
/* Exact number from a value known to be an integer or a real distance. */
Expr* gmet_distance_value(double d, int weighted);

/* ---- Per-graph result cache ----------------------------------------------
 * A small cache of finished results keyed on (head, graph node, remaining
 * arguments). Holding a reference to the graph node keeps it alive and hence
 * immutable, the same argument the validated-graph memo in graph_util.c rests
 * on, so a hit is always the answer the head would compute. Only results of
 * at most GMET_CACHE_MAX_ELEMS leaves are stored. */
#define GMET_CACHE_MAX_ELEMS 4000000
Expr* gmet_cache_get(const char* head, const Expr* res);
void  gmet_cache_put(const char* head, const Expr* res, Expr* value);

#endif /* GRAPH_METRICS_H */
