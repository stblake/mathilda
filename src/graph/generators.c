/* generators.c - standard graph constructors.
 *
 *   CompleteGraph[n]          - undirected K_n (all n(n-1)/2 edges)
 *   CycleGraph[n]             - undirected cycle on 1..n
 *   PathGraph[n]              - undirected path 1-2-...-n
 *   PathGraph[{v1,...,vk}]    - undirected path over the given vertices
 *   RandomGraph[{n, m}]       - undirected graph with n vertices, m random edges
 *   RandomGraph[{n, m}, k]    - a list of k such graphs
 *
 * Each assembles a Graph[List verts, List edges] expression and returns it; the
 * evaluator canonicalizes and validates it via builtin_graph. Vertices are the
 * integers 1..n (except the explicit PathGraph[{...}] form). RandomGraph reuses
 * the system RNG by evaluating RandomSample over the candidate edges, so it
 * honors SeedRandom.
 *
 * Memory (SPEC section 4): returns freshly-allocated trees; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
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

/* Wrap vertex/edge C-arrays into a Graph[...] (moves ownership). */
static Expr* make_graph(Expr** verts, size_t nv, Expr** edges, size_t ne) {
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), edges, ne);
    Expr* gargs[2] = { vlist, elist };
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

/* Vertices 1..n wrapped as a List; frees the intermediate C array, which
 * expr_new_function memcpys rather than adopting (src/expr.c:257). */
static Expr* vertex_list(long n) {
    Expr** verts = int_vertices(n);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, (size_t)n);
    free(verts);
    return vlist;
}

/* One random undirected graph: n vertices, m of the n(n-1)/2 candidate edges.
 * Returns NULL if the sampler declines or an allocation fails. Caller has
 * already validated n >= 0, m >= 0, maxe representable, m <= maxe.
 *
 * Assembly is inline rather than via make_graph (:36-41) on purpose: that
 * helper frees none of the arrays it is handed (its "moves ownership" comment
 * predates the memcpy) and takes an Expr** edge array, not a built List. */
static Expr* one_random_graph(long n, unsigned long long maxe, long m) {
    /* n <= 1 leaves no candidates, and RandomSample[{}, m] is itself
     * unevaluated (is_nonempty_list, src/random.c:1707) — which is why
     * RandomGraph[{0,0}] failed before RG-1. Answer directly instead.
     *
     * Deliberately NOT extended to m == 0 with n >= 2: that case works via
     * RandomSample[cand, 0], and skipping the call would shift the RNG stream
     * for every later draw. */
    if (maxe == 0) {
        Expr* gargs[2] = { vertex_list(n),
                           expr_new_function(expr_new_symbol(SYM_List), NULL, 0) };
        return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
    }

    /* All candidate undirected edges. */
    size_t ncand = (size_t)maxe;
    Expr** cand = calloc(ncand, sizeof(Expr*));
    if (!cand) return NULL;                  /* absurd n: unevaluated, not a crash */
    size_t k = 0;
    for (long i = 1; i <= n; i++)
        for (long j = i + 1; j <= n; j++)
            cand[k++] = undirected_edge(i, j);
    Expr* cand_list = expr_new_function(expr_new_symbol(SYM_List), cand, ncand);
    free(cand);

    /* Sample m of them without replacement, via the seeded system RNG. */
    Expr* sample_args[2] = { cand_list, expr_new_integer(m) };
    Expr* sample_call = expr_new_function(expr_new_symbol("RandomSample"),
                                          sample_args, 2);
    Expr* sampled = evaluate(sample_call);   /* consumes sample_call */
    if (!graph_is_list(sampled)) { expr_free(sampled); return NULL; }

    Expr* gargs[2] = { vertex_list(n), sampled };
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
