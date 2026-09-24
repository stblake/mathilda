#ifndef GRAPH_ALGOS_H
#define GRAPH_ALGOS_H

/* graph_algos.h - the "algos" stream of the graph subsystem: flows and cuts,
 * matchings and covers, cliques and independent sets, Hamiltonian cycles and
 * paths, isomorphism and canonical forms, and planarity.
 *
 * Every builtin here is registered by graph_algos_init() (galg_init.c), which
 * graph_init() calls last. The implementations live in src/graph/galg_*.c:
 *
 *   galg_common.c    CSR builders over the validated-graph memo, result
 *                    builders, the deadline/budget poll
 *   galg_flow.c      FindMaximumFlow, FindMinimumCut, FindEdgeCut,
 *                    FindVertexCut, EdgeConnectivity (+ galg_vertex_connectivity)
 *   galg_matching.c  FindIndependentEdgeSet, FindEdgeCover (Hopcroft-Karp,
 *                    Edmonds blossom)
 *   galg_mis.c       FindVertexCover, FindIndependentVertexSet (exact branch
 *                    and reduce) and the four *SetQ / *CoverQ predicates
 *   galg_clique.c    FindClique, FindKClique (bitset branch and bound)
 *   galg_hamilton.c  FindHamiltonianCycle, FindHamiltonianPath, HamiltonianGraphQ
 *   galg_iso.c       individualization-refinement engine (canonical labeling,
 *                    isomorphism search) -- no Expr dependency
 *   galg_isoheads.c  IsomorphicGraphQ, FindGraphIsomorphism, CanonicalGraph,
 *                    GraphAutomorphismGroup
 *   galg_planar.c    left-right planarity test -- no Expr dependency
 *
 * Ownership follows the SPEC section 4 builtin contract: a builtin returns a
 * fresh Expr (the evaluator frees res) or NULL to stay unevaluated.
 *
 * NP-hard heads (vertex cover, independent set, clique, Hamiltonian) are exact:
 * they answer with a proven optimum or stay unevaluated. They poll
 * galg_poll() (which calls tc_check_deadline) so TimeConstrained interrupts
 * them, and they carry a deterministic node budget as a backstop.
 */

#include "expr.h"

void graph_algos_init(void);

/* ---- Builtins ------------------------------------------------------------- */
Expr* builtin_find_maximum_flow(Expr* res);      /* FindMaximumFlow            */
Expr* builtin_find_minimum_cut(Expr* res);       /* FindMinimumCut             */
Expr* builtin_find_edge_cut(Expr* res);          /* FindEdgeCut                */
Expr* builtin_find_vertex_cut(Expr* res);        /* FindVertexCut              */
Expr* builtin_edge_connectivity(Expr* res);      /* EdgeConnectivity           */

Expr* builtin_find_independent_edge_set(Expr* res); /* FindIndependentEdgeSet  */
Expr* builtin_find_edge_cover(Expr* res);        /* FindEdgeCover              */
Expr* builtin_find_vertex_cover(Expr* res);      /* FindVertexCover            */
Expr* builtin_find_independent_vertex_set(Expr* res); /* FindIndependentVertexSet */
Expr* builtin_independent_vertex_set_q(Expr* res);    /* IndependentVertexSetQ */
Expr* builtin_vertex_cover_q(Expr* res);         /* VertexCoverQ               */
Expr* builtin_independent_edge_set_q(Expr* res); /* IndependentEdgeSetQ        */
Expr* builtin_edge_cover_q(Expr* res);           /* EdgeCoverQ                 */

Expr* builtin_find_clique(Expr* res);            /* FindClique                 */
Expr* builtin_find_k_clique(Expr* res);          /* FindKClique                */

Expr* builtin_find_hamiltonian_cycle(Expr* res); /* FindHamiltonianCycle       */
Expr* builtin_find_hamiltonian_path(Expr* res);  /* FindHamiltonianPath        */
Expr* builtin_hamiltonian_graph_q(Expr* res);    /* HamiltonianGraphQ          */

Expr* builtin_isomorphic_graph_q(Expr* res);     /* IsomorphicGraphQ           */
Expr* builtin_find_graph_isomorphism(Expr* res); /* FindGraphIsomorphism       */
Expr* builtin_canonical_graph(Expr* res);        /* CanonicalGraph             */
Expr* builtin_graph_automorphism_group(Expr* res); /* GraphAutomorphismGroup   */

Expr* builtin_planar_graph_q(Expr* res);         /* PlanarGraphQ               */

/* ---- Shared internals (galg_common.c) ------------------------------------- */

/* Cooperative interruption: calls tc_check_deadline() (which siglongjmps out of
 * a timed-out TimeConstrained body). Call it every few thousand units of work
 * in any search whose cost is not linear. */
void galg_poll(void);

/* Simple undirected graph in CSR form: the neighbours of v are
 * adj[off[v] .. off[v+1]-1], sorted ascending, without duplicates or v itself.
 * m is the number of undirected edges (adj holds 2m entries). */
typedef struct {
    int  n;
    long m;
    int* off;   /* n + 1 entries */
    int* adj;   /* 2m entries    */
} GalgUG;

/* The underlying simple undirected graph of a valid Graph expression: every
 * edge, directed or not, joins its endpoints, and an anti-parallel pair or a
 * directed edge alongside an undirected one collapses to one neighbour.
 * NULL if g is not a valid graph (or on allocation failure). */
GalgUG* galg_ug_from_graph(const Expr* g);

/* The same from raw endpoint arrays over vertices 0..n-1 (self-loops dropped,
 * duplicates merged). */
GalgUG* galg_ug_from_edges(int n, long m, const int* eu, const int* ev);
void    galg_ug_free(GalgUG* u);

/* Vertex count and edge count of a valid graph (no validation). */
int  galg_nv(const Expr* g);
long galg_ne(const Expr* g);

/* List[verts[idx[0]], ...] (sorted by index first when sort != 0). */
Expr* galg_vertex_list(const Expr* g, const int* idx, int k, int sort);

/* List[edges[k0], edges[k1], ...]: copies of g's own edge expressions. */
Expr* galg_edge_list(const Expr* g, const long* eidx, long k);

/* A fresh edge head[verts[u], verts[v]]; `head` (a DirectedEdge/UndirectedEdge
 * symbol node built once by the caller) is shared by reference, since interning
 * a symbol per edge is measurably slow on large results. */
Expr* galg_make_edge(const Expr* g, const Expr* head, int u, int v);

/* Symbol / string tests by name (no interning needed). */
int galg_is_symbol(const Expr* e, const char* name);
int galg_is_string(const Expr* e, const char* s);

/* A plain List view of e: e itself if it is a List, or a materialised copy
 * (returned through *owned, which the caller frees) if e is a packed list.
 * NULL if e is not a list at all. */
const Expr* galg_plain_list(const Expr* e, Expr** owned);

/* Rule / RuleDelayed with two arguments. */
int galg_is_rule(const Expr* e);

/* Index of the first Rule argument of res at or after `from` (arg_count if
 * none). */
size_t galg_first_option(const Expr* res, size_t from);

/* Scan res's arguments from `from` on: every one must be a rule whose key is
 * one of the NULL-terminated `names` (symbol or string); values[i] receives the
 * last value given for names[i] (NULL when absent). Returns 0 on any other
 * argument (the caller then stays unevaluated). */
int galg_split_options(const Expr* res, size_t from, const char* const* names,
                       const Expr** values);

/* Position of vertex v in valid graph g, or -1. */
int galg_vertex_arg(const Expr* g, const Expr* v);

/* True / False symbol. */
Expr* galg_truth(int b);

/* ---- Max-flow engine (galg_flow.c) -----------------------------------------
 * Internal entry point for a faster VertexConnectivity (not wired into
 * connectivity.c, which this stream does not own): the vertex connectivity of
 * the underlying undirected graph of g, via Even's algorithm on unit-capacity
 * split-vertex max flows. With s, t >= 0, the minimum number of vertices whose
 * removal separates s from t (0 when adjacent, matching Wolfram). -1 on error. */
long galg_vertex_connectivity(const Expr* g, int s, int t);

/* ---- Maximum matching (galg_matching.c) ------------------------------------
 * mate[v] = partner of v or -1. Returns the matching size, -1 on allocation
 * failure. Exact maximum cardinality: Hopcroft-Karp when the graph is
 * bipartite, Edmonds' blossom algorithm otherwise. */
long galg_max_matching(const GalgUG* u, int* mate);

/* ---- Cliques (galg_clique.c) ------------------------------------------------
 * Maximum clique of a simple graph: returns its size and fills out[] (sorted
 * ascending vertex ids), or -1 on budget exhaustion / allocation failure. */
int galg_max_clique(const GalgUG* u, int* out);

/* Size spec k / {k} / {kmin, kmax} (k an Integer >= 0 or Infinity) and count
 * (positive Integer, or All -> 0). Return 1 when well-formed. */
int galg_parse_size_spec(const Expr* spec, int* kmin, int* kmax);
int galg_parse_count(const Expr* e, long* count);

/* The Wolfram answer to a spec query: up to count (0: all) MAXIMAL cliques of
 * u -- of its complement when complement != 0, i.e. maximal independent sets --
 * with sizes in [kmin, kmax], as a List of vertex lists of g, largest first.
 * NULL (unevaluated) on budget exhaustion or when u is too large for the
 * bitset enumeration. */
Expr* galg_clique_spec_query(const Expr* g, const GalgUG* u, int complement,
                             int kmin, int kmax, long count);

/* Maximum independent set of the subgraph induced by vs (dense components;
 * complement-bitset branch and bound). Returns the size or -1. */
int galg_mis_dense(const GalgUG* u, const int* vs, int k, int* out);

/* ---- Independent sets (galg_mis.c) -------------------------------------------
 * Exact maximum independent set: fills out[] (vertex ids) and returns the
 * size, or -1 on budget exhaustion / allocation failure. */
int galg_max_independent_set(const GalgUG* u, int* out);

/* ---- Planarity (galg_planar.c) ---------------------------------------------
 * 1 if the simple undirected graph is planar, 0 if not, -1 on allocation
 * failure. Linear time (left-right planarity test), iterative. */
int galg_planar_test(const GalgUG* u);

/* ---- Isomorphism engine (galg_iso.c) ---------------------------------------
 * A graph for the engine: up to three adjacency relations in CSR form over the
 * vertices 0..n-1. Undirected neighbours (uoff/uadj), directed out-neighbours
 * (ooff/oadj) and directed in-neighbours (ioff/iadj); any relation may be
 * absent (NULL offsets). Neighbour lists are duplicate-free. vcol, when
 * non-NULL, is an initial vertex colouring that an isomorphism must preserve. */
typedef struct {
    int n;
    const int* uoff; const int* uadj;
    const int* ooff; const int* oadj;
    const int* ioff; const int* iadj;
    const int* vcol;
} GalgIsoGraph;

/* Search for an isomorphism g -> h. Returns 1 and fills map[0..n-1]
 * (map[v] = image of v) if one exists, 0 if the graphs are not isomorphic,
 * -1 if the search exceeded `budget` search nodes (budget <= 0: unlimited) or
 * an allocation failed. */
int galg_iso_find(const GalgIsoGraph* g, const GalgIsoGraph* h, int* map, long budget);

/* Enumerate isomorphisms g -> h, calling cb(map, ctx) for each (cb returns 0
 * to stop). Stops after max_count (max_count <= 0: all). Returns the number
 * reported, or -1 when the budget was exceeded / allocation failed before the
 * enumeration was complete. */
long galg_iso_enumerate(const GalgIsoGraph* g, const GalgIsoGraph* h, long max_count,
                        int (*cb)(const int* map, void* ctx), void* ctx, long budget);

/* Canonical labeling: fills lab[0..n-1] with a permutation (lab[i] = the
 * vertex placed at canonical position i) such that relabeling g and h by their
 * canonical labelings gives identical graphs iff g and h are isomorphic.
 * Returns 1, or -1 on budget exhaustion / allocation failure. */
int galg_iso_canon(const GalgIsoGraph* g, int* lab, long budget);

/* Automorphism group generators: fills *gens (malloc'd, ngens * n ints, each a
 * permutation image array) with a generating set of Aut(g). Returns the number
 * of generators (0 for the trivial group), or -1 on budget/allocation failure. */
long galg_iso_automorphisms(const GalgIsoGraph* g, int** gens, long budget);

#endif /* GRAPH_ALGOS_H */
