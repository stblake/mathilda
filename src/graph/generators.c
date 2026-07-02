/* generators.c - standard graph constructors.
 *
 *   CompleteGraph[n]          - undirected K_n (all n(n-1)/2 edges)
 *   CycleGraph[n]             - undirected cycle on 1..n
 *   PathGraph[n]              - undirected path 1-2-...-n
 *   PathGraph[{v1,...,vk}]    - undirected path over the given vertices
 *   RandomGraph[{n, m}]       - undirected graph with n vertices, m random edges
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

Expr* builtin_random_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* spec = res->data.function.args[0];
    if (!graph_is_list(spec) || spec->data.function.arg_count != 2) return NULL;
    long n = as_count(spec->data.function.args[0]);
    long m = as_count(spec->data.function.args[1]);
    if (n < 0 || m < 0) return NULL;
    long maxe = n * (n - 1) / 2;
    if (m > maxe) return NULL;   /* more edges than a simple graph allows */

    /* All candidate undirected edges. */
    size_t ncand = (size_t)maxe;
    Expr** cand = (ncand > 0) ? calloc(ncand, sizeof(Expr*)) : NULL;
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

    Expr* gargs[2] = { expr_new_function(expr_new_symbol(SYM_List),
                                         int_vertices(n), (size_t)n),
                       sampled };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
