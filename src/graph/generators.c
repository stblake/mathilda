/* generators.c - standard graph constructors.
 *
 *   CompleteGraph[n]          - undirected K_n (all n(n-1)/2 edges)
 *   CycleGraph[n]             - undirected cycle on 1..n
 *   PathGraph[n]              - undirected path 1-2-...-n
 *   PathGraph[{v1,...,vk}]    - undirected path over the given vertices
 *   StarGraph[n]              - undirected star: hub 1 joined to leaves 2..n
 *   RandomGraph[{n, m}]       - undirected graph with n vertices, m random edges
 *   RandomGraph[{n, m}, k]    - a list of k such graphs
 *   RandomGraph[{n, m}, DirectedEdges -> True] - m random directed edges
 *   RandomGraph[BernoulliGraphDistribution[n, p]] - each edge with prob. p
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
#include <stdbool.h>
#include "checked_alloc.h"

/* Small integer argument as a nonnegative int64_t, or -1 if not a suitable int. */
static int64_t as_count(const Expr* e) {
    if (!e || e->type != EXPR_INTEGER || e->data.integer < 0) return -1;
    return (int64_t)e->data.integer;
}

static Expr* undirected_edge(int64_t a, int64_t b) {
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

/* Resource limit shared with the other family generators (gmet_generators.c):
 * a family past either bound is left unevaluated, as the spec states, rather
 * than allocating gigabytes of edge nodes. */
#define GEN_FAMILY_MAX_VERTICES INT64_C(100000000)
#define GEN_FAMILY_MAX_EDGES    INT64_C(50000000)

/* Whether a family on n vertices with ne edges is within the resource limit. */
static bool family_fits(int64_t n, int64_t ne) {
    return n >= 0 && ne >= 0 &&
           n <= GEN_FAMILY_MAX_VERTICES && ne <= GEN_FAMILY_MAX_EDGES;
}

/* Allocate the vertex array 1..n and an (unfilled) edge array of ne slots.
 * False (both freed, General::nomem issued) when memory cannot be had. */
static bool family_alloc(int64_t n, int64_t ne, Expr*** verts, Expr*** edges) {
    *verts = mth_malloc_n(n, sizeof(Expr*));
    *edges = *verts ? mth_malloc_n(ne, sizeof(Expr*)) : NULL;
    if (!*verts || !*edges) {
        free(*verts); free(*edges);
        *verts = *edges = NULL;
        return false;
    }
    for (int64_t i = 0; i < n; i++) (*verts)[i] = expr_new_integer(i + 1);
    return true;
}

/* True iff e is a function node whose head is the interned symbol h. */
static bool head_is_sym_name(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

/* True iff e is Rule[sym, value] (an option). */
static bool is_option(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Rule &&
           e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* See graph.h. */
Expr* graph_generator_with_options(Expr* res, size_t max_pos,
                                   Expr* (*core)(Expr*, int directed)) {
    size_t argc = res->data.function.arg_count;
    size_t npos = argc;
    while (npos > 0 && is_option(res->data.function.args[npos - 1])) npos--;
    if (npos == argc) return core(res, -1);
    if (npos == 0 || npos > max_pos) return NULL;

    /* Classify the options before doing any work. */
    int directed = -1;                       /* -1 unset, 0 False, 1 True */
    size_t nkeep = 0;
    for (size_t i = npos; i < argc; i++) {
        const Expr* o = res->data.function.args[i];
        const char* key = o->data.function.args[0]->data.symbol.name;
        const Expr* val = o->data.function.args[1];
        if (key == SYM_DirectedEdges) {
            if (directed >= 0 || val->type != EXPR_SYMBOL) return NULL;
            if (val->data.symbol.name == SYM_True) directed = 1;
            else if (val->data.symbol.name == SYM_False) directed = 0;
            else return NULL;
        } else if (graph_edge_option_rank(key) >= 0) {
            nkeep++;
        } else {
            return NULL;                     /* unsupported option */
        }
    }

    Expr** pa = mth_malloc_n((int64_t)npos, sizeof(Expr*));
    if (!pa) return mth_nomem_abort();
    for (size_t i = 0; i < npos; i++) pa[i] = expr_copy(res->data.function.args[i]);
    Expr* sub = expr_new_function(expr_copy(res->data.function.head), pa, npos);
    free(pa);
    Expr* g = core(sub, directed);
    expr_free(sub);
    if (!g) return NULL;
    if (g->type != EXPR_FUNCTION || g->data.function.arg_count != 2 ||
        !graph_is_list(g->data.function.args[1])) {
        expr_free(g);
        return NULL;
    }

    /* Rebuild the edge list with the requested orientation. */
    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;
    Expr** ea = mth_malloc_n((int64_t)(ne ? ne : 1), sizeof(Expr*));
    if (!ea) { expr_free(g); return mth_nomem_abort(); }
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        if (directed == 0 && kind == SYM_DirectedEdge) {
            for (size_t j = 0; j < k; j++) expr_free(ea[j]);
            free(ea);
            expr_free(g);
            return NULL;                     /* cannot drop a direction */
        }
        if (directed == 1 && kind == SYM_UndirectedEdge) {
            Expr* ab[2] = { expr_copy(e->data.function.args[0]),
                            expr_copy(e->data.function.args[1]) };
            ea[k] = expr_new_function(expr_new_symbol(SYM_DirectedEdge), ab, 2);
        } else {
            ea[k] = expr_copy((Expr*)e);
        }
    }
    Expr* elist = make_list_owning(ea, ne);

    size_t gargc = 2 + nkeep, gi = 0;
    Expr** ga = mth_malloc_n((int64_t)gargc, sizeof(Expr*));
    if (!ga) { expr_free(elist); expr_free(g); return mth_nomem_abort(); }
    ga[gi++] = expr_copy(g->data.function.args[0]);
    ga[gi++] = elist;
    for (size_t i = npos; i < argc; i++) {
        const Expr* o = res->data.function.args[i];
        if (o->data.function.args[0]->data.symbol.name != SYM_DirectedEdges)
            ga[gi++] = expr_copy((Expr*)o);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_Graph), ga, gargc);
    free(ga);
    expr_free(g);
    return out;
}

/* K_n: n(n-1)/2 edges i < j, in lexicographic order. */
Expr* builtin_complete_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int64_t n = as_count(res->data.function.args[0]);
    if (n < 0 || n > GEN_FAMILY_MAX_VERTICES) return NULL;
    int64_t ne = (n < 2) ? 0 : n * (n - 1) / 2;     /* n <= 1e8: no overflow */
    if (!family_fits(n, ne)) return NULL;
    Expr **verts, **edges;
    if (!family_alloc(n, ne, &verts, &edges)) return mth_nomem_abort();
    int64_t k = 0;
    for (int64_t i = 1; i <= n; i++)
        for (int64_t j = i + 1; j <= n; j++)
            edges[k++] = undirected_edge(i, j);
    return make_graph(verts, (size_t)n, edges, (size_t)ne);
}

static Expr* cycle_graph_core(Expr* res, int directed) {
    (void)directed;   /* applied by graph_generator_with_options */
    if (res->data.function.arg_count != 1) return NULL;
    int64_t n = as_count(res->data.function.args[0]);
    /* Path edges 1-2-...-n, plus the wrap edge n-1 when n >= 3 (for n <= 2 the
     * wrap edge would duplicate an existing one). */
    int64_t ne = (n >= 3) ? n : (n >= 2 ? 1 : 0);
    if (!family_fits(n, ne)) return NULL;
    Expr **verts, **edges;
    if (!family_alloc(n, ne, &verts, &edges)) return mth_nomem_abort();
    int64_t k = 0;
    for (int64_t i = 1; i < n; i++) edges[k++] = undirected_edge(i, i + 1);
    if (n >= 3) edges[k++] = undirected_edge(n, 1);
    return make_graph(verts, (size_t)n, edges, (size_t)ne);
}

static Expr* path_graph_core(Expr* res, int directed) {
    (void)directed;
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* arg = res->data.function.args[0];

    if (graph_is_list(arg)) {
        /* PathGraph[{v1,...,vk}] over explicit vertices (the list already
         * exists, so its size is already paid for). */
        size_t nv = arg->data.function.arg_count;
        Expr** verts = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
        size_t ne = (nv > 0) ? nv - 1 : 0;
        Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
        if ((nv > 0 && !verts) || (ne > 0 && !edges)) {
            free(verts); free(edges);
            return mth_nomem_abort();
        }
        for (size_t i = 0; i < nv; i++) verts[i] = expr_copy(arg->data.function.args[i]);
        for (size_t i = 0; i + 1 < nv; i++) {
            Expr* ea[2] = { expr_copy(arg->data.function.args[i]),
                            expr_copy(arg->data.function.args[i + 1]) };
            edges[i] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), ea, 2);
        }
        return make_graph(verts, nv, edges, ne);
    }

    int64_t n = as_count(arg);
    int64_t ne = (n > 0) ? n - 1 : 0;
    if (!family_fits(n, ne)) return NULL;
    Expr **verts, **edges;
    if (!family_alloc(n, ne, &verts, &edges)) return mth_nomem_abort();
    for (int64_t i = 1; i < n; i++) edges[i - 1] = undirected_edge(i, i + 1);
    return make_graph(verts, (size_t)n, edges, (size_t)ne);
}

static Expr* star_graph_core(Expr* res, int directed) {
    (void)directed;
    if (res->data.function.arg_count != 1) return NULL;
    int64_t n = as_count(res->data.function.args[0]);
    /* Hub is vertex 1, joined to each of 2..n: exactly n-1 edges, with no
     * duplicate possible at any size (unlike CycleGraph's wrap edge). */
    int64_t ne = (n > 0) ? n - 1 : 0;
    if (!family_fits(n, ne)) return NULL;
    Expr **verts, **edges;
    if (!family_alloc(n, ne, &verts, &edges)) return mth_nomem_abort();
    for (int64_t i = 2; i <= n; i++) edges[i - 2] = undirected_edge(1, i);
    return make_graph(verts, (size_t)n, edges, (size_t)ne);
}

/* Public heads: the cores above, plus trailing options (DirectedEdges,
 * EdgeWeight, EdgeCapacity) via graph_generator_with_options. */
Expr* builtin_cycle_graph(Expr* res) {
    return graph_generator_with_options(res, 1, cycle_graph_core);
}
Expr* builtin_path_graph(Expr* res) {
    return graph_generator_with_options(res, 1, path_graph_core);
}
Expr* builtin_star_graph(Expr* res) {
    return graph_generator_with_options(res, 1, star_graph_core);
}

/* List[1, ..., n], or NULL (General::nomem issued) when memory cannot be had. */
static Expr* int_vertex_list(int64_t n) {
    Expr** v = mth_malloc_n(n, sizeof(Expr*));
    if (!v) return NULL;
    for (int64_t i = 0; i < n; i++) v[i] = expr_new_integer(i + 1);
    return make_list_owning(v, (size_t)n);
}

/* One random undirected graph: n vertices, m of the n(n-1)/2 candidate edges.
 * Returns NULL if the sampler declines or an allocation fails. Caller has
 * already validated n >= 0, m >= 0, maxe representable, m <= maxe.
 */
static Expr* directed_edge(int64_t a, int64_t b) {
    Expr* ea[2] = { expr_new_integer(a), expr_new_integer(b) };
    return expr_new_function(expr_new_symbol(SYM_DirectedEdge), ea, 2);
}

static Expr* one_random_graph(int64_t n, unsigned long long maxe, int64_t m, bool directed) {
    /* n <= 1 leaves no candidates: the empty graph. */
    if (maxe == 0) {
        Expr* vl = int_vertex_list(n);
        if (!vl) return NULL;
        Expr* gargs[2] = { vl, expr_new_function(expr_new_symbol(SYM_List), NULL, 0) };
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
    for (int64_t e = 0; e < m; e++) {
        unsigned long long x = idx[e];
        if (directed) {
            /* the n(n-1) ordered pairs, row-major: row i holds (i, j != i) */
            int64_t i = (int64_t)(x / (unsigned long long)(n - 1)) + 1;
            int64_t k = (int64_t)(x % (unsigned long long)(n - 1));
            es[e] = directed_edge(i, k + 1 < i ? k + 1 : k + 2);
            continue;
        }
        /* row i (1-based) holds pairs (i, i+1..n) and starts at
         * off(i) = (i-1)(2n-i)/2; find the last row with off(i) <= idx */
        int64_t lo = 1, hi = n - 1;
        while (lo < hi) {
            int64_t mid = lo + (hi - lo + 1) / 2;
            unsigned long long off = (unsigned long long)(mid - 1)
                * (unsigned long long)(2 * n - mid) / 2;
            if (off <= x) lo = mid; else hi = mid - 1;
        }
        unsigned long long off = (unsigned long long)(lo - 1)
            * (unsigned long long)(2 * n - lo) / 2;
        es[e] = undirected_edge(lo, lo + 1 + (int64_t)(x - off));
    }
    free(idx);
    Expr* sampled = expr_new_function(expr_new_symbol(SYM_List), es, (size_t)m);
    free(es);

    Expr* vl = int_vertex_list(n);
    if (!vl) { expr_free(sampled); return NULL; }
    Expr* gargs[2] = { vl, sampled };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}

/* RandomGraph[BernoulliGraphDistribution[n, p]]: each of the n(n-1)/2
 * possible undirected edges present independently with probability p, one
 * uniform draw per candidate in row-major order (so it honours SeedRandom).
 * NULL for a malformed distribution or past the resource limit. */
static Expr* bernoulli_random_graph(const Expr* dist) {
    if (dist->data.function.arg_count != 2) return NULL;
    int64_t n = as_count(dist->data.function.args[0]);
    const Expr* pe = dist->data.function.args[1];
    double p;
    if (pe->type == EXPR_REAL) p = pe->data.real;
    else if (pe->type == EXPR_INTEGER) p = (double)pe->data.integer;
    else if (head_is_sym_name(pe, SYM_Rational) && pe->data.function.arg_count == 2 &&
             pe->data.function.args[0]->type == EXPR_INTEGER &&
             pe->data.function.args[1]->type == EXPR_INTEGER)
        p = (double)pe->data.function.args[0]->data.integer /
            (double)pe->data.function.args[1]->data.integer;
    else return NULL;
    if (n < 0 || !(p >= 0.0 && p <= 1.0)) return NULL;
    int64_t maxe = n < 2 ? 0 : n * (n - 1) / 2;
    if (n > GEN_FAMILY_MAX_VERTICES || !family_fits(n, maxe)) return NULL;
    Expr** es = mth_malloc_n(maxe > 0 ? maxe : 1, sizeof(Expr*));
    if (!es) return mth_nomem_abort();
    size_t m = 0;
    for (int64_t i = 1; i <= n; i++)
        for (int64_t j = i + 1; j <= n; j++)
            if (random_uniform_01() < p) es[m++] = undirected_edge(i, j);
    Expr* el = expr_new_function(expr_new_symbol(SYM_List), es, m);
    free(es);
    Expr* vl = int_vertex_list(n);
    if (!vl) { expr_free(el); return NULL; }
    Expr* gargs[2] = { vl, el };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}

Expr* builtin_random_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    /* a trailing DirectedEdges -> True | False */
    bool directed = false;
    if (argc >= 2 && is_option(res->data.function.args[argc - 1])) {
        const Expr* o = res->data.function.args[argc - 1];
        const Expr* v = o->data.function.args[1];
        if (o->data.function.args[0]->data.symbol.name != SYM_DirectedEdges ||
            v->type != EXPR_SYMBOL ||
            (v->data.symbol.name != SYM_True && v->data.symbol.name != SYM_False))
            return NULL;
        directed = v->data.symbol.name == SYM_True;
        argc--;
    }
    if (argc != 1 && argc != 2) return NULL;
    const Expr* spec = res->data.function.args[0];
    if (head_is_sym_name(spec, SYM_BernoulliGraphDistribution)) {
        if (directed) return NULL;            /* undirected distribution only */
        if (argc == 1) return bernoulli_random_graph(spec);
        int64_t kc = as_count(res->data.function.args[1]);
        if (kc < 0 || kc > GEN_FAMILY_MAX_VERTICES) return NULL;
        Expr** gs = kc > 0 ? calloc((size_t)kc, sizeof(Expr*)) : NULL;
        if (kc > 0 && !gs) return NULL;
        for (int64_t i = 0; i < kc; i++) {
            gs[i] = bernoulli_random_graph(spec);
            if (!gs[i]) {
                for (int64_t j = 0; j < i; j++) expr_free(gs[j]);
                free(gs);
                return NULL;
            }
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), gs, (size_t)kc);
        free(gs);
        return out;
    }
    if (!graph_is_list(spec) || spec->data.function.arg_count != 2) return NULL;
    int64_t n = as_count(spec->data.function.args[0]);
    int64_t m = as_count(spec->data.function.args[1]);
    if (n < 0 || m < 0) return NULL;
    if (!family_fits(n, m)) return NULL;          /* the shared resource limit */

    /* n(n-1)/2 overflows signed int64_t for large n, and a negative maxe would let
     * the m gate below pass. Bound n BEFORE multiplying: checking the product
     * afterwards cannot catch a wrap, since n = 2^32 + 1 wraps n(n-1) to 2^32
     * and yields a small, plausible-looking maxe while the candidate loop still
     * runs to the true n(n-1)/2 and writes past the buffer. With n < 2^31 the
     * product cannot exceed 2^62, so the unsigned arithmetic below is exact and
     * the SIZE_MAX bound is meaningful. An n past that is unevaluated, per the
     * head's error channel. */
    if (n > 2147483647L) return NULL;
    unsigned long long maxe = (n < 2) ? 0ULL
        : (unsigned long long)n * (unsigned long long)(n - 1) / (directed ? 1 : 2);
    if (maxe > (unsigned long long)(SIZE_MAX / sizeof(Expr*))) return NULL;
    if ((unsigned long long)m > maxe) return NULL;  /* more edges than a simple graph allows */

    if (argc == 1) return one_random_graph(n, maxe, m, directed);

    /* RandomGraph[{n,m}, k]: k independent graphs. k = 0 gives {}; a negative,
     * non-integer, or symbolic k is silently unevaluated, matching the five
     * count-taking heads in src/random.c (SPEC section 4: NULL, no Message). */
    int64_t kcount = as_count(res->data.function.args[1]);
    if (kcount < 0) return NULL;
    if ((unsigned long long)kcount > (unsigned long long)(SIZE_MAX / sizeof(Expr*)))
        return NULL;
    Expr** gs = (kcount > 0) ? calloc((size_t)kcount, sizeof(Expr*)) : NULL;
    if (kcount > 0 && !gs) return NULL;
    for (int64_t i = 0; i < kcount; i++) {
        gs[i] = one_random_graph(n, maxe, m, directed);
        if (!gs[i]) {                        /* never return a partial list */
            for (int64_t j = 0; j < i; j++) expr_free(gs[j]);
            free(gs);
            return NULL;
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), gs, (size_t)kcount);
    free(gs);
    return out;
}
