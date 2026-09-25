/* generators.c - standard graph constructors.
 *
 *   CompleteGraph[n]          - undirected K_n (all n(n-1)/2 edges)
 *   CycleGraph[n]             - undirected cycle on 1..n
 *   PathGraph[n]              - undirected path 1-2-...-n
 *   PathGraph[{v1,...,vk}]    - undirected path over the given vertices
 *   StarGraph[n]              - undirected star: hub 1 joined to leaves 2..n
 *   RandomGraph[{n, m}]       - undirected graph with n vertices, m random edges
 *   RandomGraph[{n, m}, k]    - a list of k such graphs
 *
 * Each assembles a Graph[List verts, List edges] expression and returns it; the
 * evaluator canonicalizes and validates it via builtin_graph. Vertices are the
 * integers 1..n (except the explicit PathGraph[{...}] form). RandomGraph draws
 * from the system RNG exactly as RandomSample over the candidate edges would
 * (random_sample_indices), so it honors SeedRandom.
 *
 * Memory (SPEC section 4): returns freshly-allocated trees; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "random.h"
#include <stdlib.h>
#include <stdint.h>   /* SIZE_MAX */

/* Small integer argument as a nonnegative long, or -1 if not a suitable int. */
static long as_count(const Expr* e) {
    if (!e || e->type != EXPR_INTEGER || e->data.integer < 0) return -1;
    return (long)e->data.integer;
}

static Expr* undirected_edge(long a, long b) {
    Expr* ea[2] = { expr_new_integer(a), expr_new_integer(b) };
    return expr_new_function(expr_new_symbol(SYM_UndirectedEdge), ea, 2);
}

/* Wrap a calloc'd Expr* array into a List[...]: the elements move into the new
 * node (expr_new_function copies the pointers, not the array), so the array
 * itself is ours to free. */
static Expr* make_list_owning(Expr** items, size_t n) {
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return list;
}

/* Wrap vertex/edge C-arrays into a Graph[...] (moves ownership, frees both
 * arrays). */
static Expr* make_graph(Expr** verts, size_t nv, Expr** edges, size_t ne) {
    Expr* gargs[2] = { make_list_owning(verts, nv), make_list_owning(edges, ne) };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}

static Expr** int_vertices(long n) {
    Expr** v = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (long i = 0; i < n; i++) v[i] = expr_new_integer(i + 1);
    return v;
}

Expr* builtin_complete_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    long n = as_count(res->data.function.args[0]);
    if (n < 0) return NULL;
    size_t ne = (size_t)n * (size_t)(n - 1) / 2;
    Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    size_t k = 0;
    for (long i = 1; i <= n; i++)
        for (long j = i + 1; j <= n; j++)
            edges[k++] = undirected_edge(i, j);
    return make_graph(int_vertices(n), (size_t)n, edges, ne);
}

Expr* builtin_cycle_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    long n = as_count(res->data.function.args[0]);
    if (n < 0) return NULL;
    /* Path edges 1-2-...-n, plus the wrap edge n-1 when n >= 3 (for n <= 2 the
     * wrap edge would duplicate an existing one). */
    size_t ne = (n >= 3) ? (size_t)n : (n >= 2 ? 1u : 0u);
    Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    size_t k = 0;
    for (long i = 1; i < n; i++) edges[k++] = undirected_edge(i, i + 1);
    if (n >= 3) edges[k++] = undirected_edge(n, 1);
    return make_graph(int_vertices(n), (size_t)n, edges, ne);
}

Expr* builtin_path_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* arg = res->data.function.args[0];

    if (graph_is_list(arg)) {
        /* PathGraph[{v1,...,vk}] over explicit vertices. */
        size_t nv = arg->data.function.arg_count;
        Expr** verts = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
        for (size_t i = 0; i < nv; i++) verts[i] = expr_copy(arg->data.function.args[i]);
        size_t ne = (nv > 0) ? nv - 1 : 0;
        Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
        for (size_t i = 0; i + 1 < nv; i++) {
            Expr* ea[2] = { expr_copy(arg->data.function.args[i]),
                            expr_copy(arg->data.function.args[i + 1]) };
            edges[i] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), ea, 2);
        }
        return make_graph(verts, nv, edges, ne);
    }

    long n = as_count(arg);
    if (n < 0) return NULL;
    size_t ne = (n > 0) ? (size_t)n - 1 : 0;
    Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    for (long i = 1; i < n; i++) edges[i - 1] = undirected_edge(i, i + 1);
    return make_graph(int_vertices(n), (size_t)n, edges, ne);
}

Expr* builtin_star_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    long n = as_count(res->data.function.args[0]);
    if (n < 0) return NULL;
    /* Hub is vertex 1, joined to each of 2..n: exactly n-1 edges, with no
     * duplicate possible at any size (unlike CycleGraph's wrap edge). */
    size_t ne = (n > 0) ? (size_t)n - 1 : 0;
    Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    for (long i = 2; i <= n; i++) edges[i - 2] = undirected_edge(1, i);
    return make_graph(int_vertices(n), (size_t)n, edges, ne);
}

/* One random undirected graph: n vertices, m of the n(n-1)/2 candidate edges.
 * Returns NULL if the sampler declines or an allocation fails. Caller has
 * already validated n >= 0, m >= 0, maxe representable, m <= maxe.
 */
static Expr* one_random_graph(long n, unsigned long long maxe, long m) {
    /* n <= 1 leaves no candidates: the empty graph. */
    if (maxe == 0) {
        Expr* gargs[2] = { make_list_owning(int_vertices(n), (size_t)n),
                           expr_new_function(expr_new_symbol(SYM_List), NULL, 0) };
        return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
    }

    /* Sample m of the maxe candidate edges without replacement. The candidate
     * list {1<->2, 1<->3, ..., (n-1)<->n} (row-major) is never built:
     * random_sample_indices makes exactly RandomSample[cand, m]'s draws, so a
     * seeded RandomGraph is unchanged, and each index decodes to its pair --
     * O(m) memory rather than O(n^2). */
    size_t* idx = random_sample_indices((size_t)maxe, (size_t)m);
    if (!idx) return NULL;
    Expr** es = malloc((size_t)(m > 0 ? m : 1) * sizeof(Expr*));
    if (!es) { free(idx); return NULL; }
    for (long e = 0; e < m; e++) {
        /* row i (1-based) holds pairs (i, i+1..n) and starts at
         * off(i) = (i-1)(2n-i)/2; find the last row with off(i) <= idx */
        unsigned long long x = idx[e];
        long lo = 1, hi = n - 1;
        while (lo < hi) {
            long mid = lo + (hi - lo + 1) / 2;
            unsigned long long off = (unsigned long long)(mid - 1)
                * (unsigned long long)(2 * n - mid) / 2;
            if (off <= x) lo = mid; else hi = mid - 1;
        }
        unsigned long long off = (unsigned long long)(lo - 1)
            * (unsigned long long)(2 * n - lo) / 2;
        es[e] = undirected_edge(lo, lo + 1 + (long)(x - off));
    }
    free(idx);
    Expr* sampled = expr_new_function(expr_new_symbol(SYM_List), es, (size_t)m);
    free(es);

    Expr* gargs[2] = { make_list_owning(int_vertices(n), (size_t)n), sampled };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}

Expr* builtin_random_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* spec = res->data.function.args[0];
    if (!graph_is_list(spec) || spec->data.function.arg_count != 2) return NULL;
    long n = as_count(spec->data.function.args[0]);
    long m = as_count(spec->data.function.args[1]);
    if (n < 0 || m < 0) return NULL;

    /* n(n-1)/2 overflows signed long for large n, and a negative maxe would let
     * the m gate below pass. Bound n BEFORE multiplying: checking the product
     * afterwards cannot catch a wrap, since n = 2^32 + 1 wraps n(n-1) to 2^32
     * and yields a small, plausible-looking maxe while the candidate loop still
     * runs to the true n(n-1)/2 and writes past the buffer. With n < 2^31 the
     * product cannot exceed 2^62, so the unsigned arithmetic below is exact and
     * the SIZE_MAX bound is meaningful. An n past that is unevaluated, per the
     * head's error channel. */
    if (n > 2147483647L) return NULL;
    unsigned long long maxe = (n < 2) ? 0ULL
        : (unsigned long long)n * (unsigned long long)(n - 1) / 2;
    if (maxe > (unsigned long long)(SIZE_MAX / sizeof(Expr*))) return NULL;
    if ((unsigned long long)m > maxe) return NULL;  /* more edges than a simple graph allows */

    if (argc == 1) return one_random_graph(n, maxe, m);

    /* RandomGraph[{n,m}, k]: k independent graphs. k = 0 gives {}; a negative,
     * non-integer, or symbolic k is silently unevaluated, matching the five
     * count-taking heads in src/random.c (SPEC section 4: NULL, no Message). */
    long kcount = as_count(res->data.function.args[1]);
    if (kcount < 0) return NULL;
    if ((unsigned long long)kcount > (unsigned long long)(SIZE_MAX / sizeof(Expr*)))
        return NULL;
    Expr** gs = (kcount > 0) ? calloc((size_t)kcount, sizeof(Expr*)) : NULL;
    if (kcount > 0 && !gs) return NULL;
    for (long i = 0; i < kcount; i++) {
        gs[i] = one_random_graph(n, maxe, m);
        if (!gs[i]) {                        /* never return a partial list */
            for (long j = 0; j < i; j++) expr_free(gs[j]);
            free(gs);
            return NULL;
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), gs, (size_t)kcount);
    free(gs);
    return out;
}
