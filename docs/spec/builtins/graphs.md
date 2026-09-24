# Graphs

A graph subsystem modeled on the Wolfram Language's, implemented in
`src/graph/` (one builtin per translation unit, mirroring `src/linalg/`).
Graphs are ordinary `Expr` trees — there is **no new `EXPR_*` tag** — of the
canonical form

```
Graph[ List[v1, v2, ...], List[edge1, edge2, ...] ]
```

where each edge is `DirectedEdge[u, v]` or `UndirectedEdge[u, v]`. On
construction, `Rule`/`u -> v` is accepted as shorthand for `DirectedEdge`, and
`TwoWayRule`/`u <-> v` for `UndirectedEdge`. Vertices are arbitrary
expressions. Because graphs are plain expressions, generic tools (`Part`,
`Map`, `ReplaceAll`, …) work on them, and `AdjacencyMatrix[g]` returns a dense
`List`-of-`List`s consumable directly by `Det`, `Tr`, and `Eigenvalues`.

A graph may optionally carry per-edge weights via a third constructor
argument, `Graph[v, e, EdgeWeight -> {w1, ..., wm}]` — see `EdgeWeight` and
`WeightedAdjacencyMatrix` below.

**MVP scope (locked):** simple graphs only — no parallel edges, no self-loops,
no edge tags beyond `EdgeWeight`, no multigraphs, no hypergraphs, and no
vertex weights. Weighted shortest-path/distance and derived-vertex weighted
construction (`Graph[e, EdgeWeight -> {...}]`, no explicit vertex list) remain
out of scope.

### Performance model: the validated-graph memo

Because graphs are plain expressions, every accessor must know its argument is
a valid graph. Rather than re-validate (a vertex hash index plus a parallel-edge
check, `O(V + E)`) on every call, `src/graph/graph_util.c` keeps a small memo of
the last few validated graph **nodes**, holding for each one a reference, the
vertex index, the edge-key set, and every edge's endpoint indices. `Graph[...]`
seeds it during construction, from the index it builds anyway. So on a graph
you hold:

- `VertexQ`, `EdgeQ`, `EmptyGraphQ`, `UndirectedGraphQ`, `DirectedGraphQ`,
  `CompleteGraphQ` (single-kind graphs) are `O(1)`;
- adjacency-based algorithms skip validation and build a CSR adjacency with no
  hashing;
- `AcyclicGraphQ`, `TreeGraphQ`, `BipartiteGraphQ` and `TopologicalSort` cache
  their answer on the graph, so repeating a query is `O(1)`. Mathematica's atomic
  `Graph` object caches properties the same way.

Keying on the node pointer is sound because the memo holds a reference: a
referenced node can't be freed (so its address can't be reused), and a node with
more than one reference is immutable (mutators `expr_unshare` first). A
structurally equal graph at a different address just misses and is validated
again. The memo keeps at most 4 graphs alive past their last user reference.

## Graph
A graph value.
- `Graph[v, e]`: a graph with vertex list `v` and edge list `e`.
- `Graph[e]`: derives the vertex set from the edges, in first-appearance order
  (directed by default).
- `Graph[v, e, EdgeWeight -> {w1, ..., wm}]`: a weighted graph — `wi` is the
  weight of `e[[i]]`, matched by position. Requires the explicit-vertex form;
  `Graph[e, EdgeWeight -> {...}]` (derived vertices) is not accepted. A weight
  list whose length doesn't match `e` is malformed, same as any other
  rejection below.

On construction the edge list is normalized and validated, producing the
canonical `Graph[List[verts], List[edges]]` (or, when weighted,
`Graph[List[verts], List[edges], EdgeWeight -> List[weights]]`):
- `u -> v` (`Rule`) and `DirectedEdge[u, v]` become `DirectedEdge[u, v]`.
- `u <-> v` (`TwoWayRule`) and `UndirectedEdge[u, v]` become
  `UndirectedEdge[u, v]`.

Malformed input is left unevaluated: self-loops, parallel/duplicate edges,
3-argument edges, an edge endpoint absent from an explicit vertex list, or (for
a weighted graph) an `EdgeWeight` list whose length doesn't match the edge
list. (Anti-parallel directed edges `u -> v` and `v -> u` are distinct and
allowed.)

Printing: in standard output a graph shows a terse summary,
`Graph[<n vertices, m edges>]`. `InputForm` and `FullForm` print the literal
constructor, which round-trips through the parser.

```
Graph[{1,2,3,4}, {1->2, 2->3, 3->4, 4->1}]   (* Graph[<4 vertices, 4 edges>] *)
InputForm[Graph[{1,2}, {1<->2}]]              (* Graph[{1, 2}, {1 <-> 2}]      *)
```

## GraphQ
`GraphQ[g]` gives `True` if `g` is a valid graph, and `False` otherwise. A graph
is valid when it is the canonical `Graph[List, List]` with every edge a
2-argument `DirectedEdge`/`UndirectedEdge`, no self-loops, no parallel edges,
and every endpoint present in the vertex list.

```
GraphQ[Graph[{1,2}, {1->2}]]   (* True  *)
GraphQ[Graph[{1},   {1->1}]]   (* False -- self-loop *)
GraphQ[5]                      (* False *)
```

## Query / representation

All are thin readers over the canonical form and return unevaluated on a
non-graph argument.

- `VertexList[g]` — the vertices, in canonical order.
- `EdgeList[g]` — the edges (canonical `Directed`/`UndirectedEdge` form).
- `VertexCount[g]` / `EdgeCount[g]` — cardinalities.
- `AdjacencyList[g]` — `{neighbors(v1), …}` in vertex order;
  `AdjacencyList[g, v]` — neighbors of `v`. Directed edges contribute
  successors (`v -> u` makes `u` a neighbor of `v`); undirected edges go both
  ways.
- `VertexDegree[g]` / `VertexDegree[g, v]` — total degree (incident edges).
  `VertexInDegree` / `VertexOutDegree` give in-/out-degrees: a `DirectedEdge`
  adds to the source's out-degree and target's in-degree; an `UndirectedEdge`
  adds to both in- and out-degree of each endpoint.
- `DirectedGraphQ[g]` — `True` iff `g` is a valid graph with at least one
  edge, all of them directed. An edgeless graph counts as undirected (as in the
  Wolfram Language), so `DirectedGraphQ` and `UndirectedGraphQ` are never both
  `True`.
- `EdgeWeight[g]` — the weights of `g`'s edges, in `EdgeList` order. Defaults
  to all `1`s when `g` was built without an `EdgeWeight` option.

```
VertexList[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]   (* {1, 2, 3, 4}        *)
EdgeCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]    (* 4                   *)
VertexDegree[Graph[{1,2,3},{1<->2,2<->3}]]           (* {1, 2, 1}           *)
AdjacencyList[Graph[{1,2,3},{1<->2,2<->3}], 2]       (* {1, 3}              *)
EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]   (* {5, 7}        *)
EdgeWeight[Graph[{1,2,3},{1->2,2->3}]]                     (* {1, 1}        *)
```

## Matrix views (linear-algebra interop)

- `AdjacencyMatrix[g]` — the dense 0/1 adjacency matrix (`n x n`, canonical
  vertex order), symmetric for undirected graphs. It is an ordinary matrix, so
  `Det`, `Tr`, `Eigenvalues`, etc. apply directly.
- `IncidenceMatrix[g]` — the `|V| x |E|` incidence matrix; undirected edges mark
  both endpoints with `1`, directed edges are oriented (`-1` tail, `+1` head).
- `AdjacencyGraph[m]` — the inverse of `AdjacencyMatrix`: builds a graph on
  vertices `1..n` from a 0/1 matrix (undirected if `m` is symmetric, else
  directed). `AdjacencyGraph[AdjacencyMatrix[g]]` reproduces `g`'s edges.
- `WeightedAdjacencyMatrix[g]` — like `AdjacencyMatrix[g]`, but each nonzero
  entry is the corresponding edge's weight instead of `1` (`0` where there is
  no edge). Equal to `AdjacencyMatrix[g]` exactly when `g` has no
  `EdgeWeight` (every weight defaults to `1`).

```
AdjacencyMatrix[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
    (* {{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,0}} *)
Det[AdjacencyMatrix[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]]   (* -1 *)
WeightedAdjacencyMatrix[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]
    (* {{0,5,0},{0,0,7},{0,0,0}} *)
WeightedAdjacencyMatrix[CycleGraph[4]] == AdjacencyMatrix[CycleGraph[4]]  (* True *)
```

`FindShortestPath`/`GraphDistance` are weight-aware — see Search & computation below.

## Generators

Each builds a canonical graph (vertices `1..n`, undirected edges) via the
constructor path:

- `CompleteGraph[n]` — `K_n`, all `n(n-1)/2` edges.
- `CycleGraph[n]` — the cycle on `1..n`.
- `PathGraph[n]` — the path `1-2-...-n`; `PathGraph[{v1,...}]` uses the given
  vertices.
- `StarGraph[n]` — the star on `1..n`: the hub `1` joined to each of the `n-1`
  leaves `2..n`.
- `RandomGraph[{n, m}]` — a random undirected graph with `n` vertices and `m`
  distinct edges (uses the seeded system RNG, so `SeedRandom` makes it
  reproducible). Returns unevaluated if `m` exceeds `n(n-1)/2`. For `n <= 1`
  the only possible graph is the edgeless one, which is what is returned.
- `RandomGraph[{n, m}, k]` — a list of `k` independently sampled such graphs.
  `k = 0` gives `{}`; `k = 1` gives a one-element list, **not** a bare `Graph`.
  A negative, non-integer, or symbolic `k` leaves the expression unevaluated,
  silently — the convention the count-taking `Random*` heads share. Each
  element costs a full `O(n^2)` candidate materialisation and its own
  `RandomSample`, so time and peak memory scale as `O(k n^2)`. The candidate
  list is freed once sampled, so a call retains only the graphs it returns
  (until v0.183 it leaked roughly 364 KB per element at `n = 50`).

```
EdgeCount[CompleteGraph[5]]      (* 10                        *)
EdgeList[CycleGraph[4]]          (* {1<->2, 2<->3, 3<->4, 4<->1} *)
VertexDegree[PathGraph[5]]       (* {1, 2, 2, 2, 1}           *)
Length[RandomGraph[{6, 5}, 3]]   (* 3                         *)
EdgeList[StarGraph[4]]           (* {1<->2, 1<->3, 1<->4}     *)
VertexDegree[StarGraph[5]]       (* {4, 1, 1, 1, 1}           *)
```

## Search & computation

All build an integer-indexed adjacency on demand; all but `FindShortestPath`/`GraphDistance`
are unweighted.

- `FindShortestPath[g, s, t]` — a shortest path from `s` to `t` as a vertex
  list; `{}` if `t` is unreachable. **Weight-aware**: if `g` carries an
  `EdgeWeight` and every weight is non-negative and numeric, uses Dijkstra
  (minimum total weight); otherwise (unweighted, a symbolic weight, or a
  negative weight present) uses unweighted BFS (minimum hop count),
  following edge direction for directed graphs either way.
- `GraphDistance[g, s, t]` — the length/total weight of that path;
  `Infinity` if unreachable. Same weight-aware dispatch as `FindShortestPath`,
  and returns an exact value (`Integer`/`Rational`) whenever the weights are
  exact — never a `Real` artifact of the internal algorithm.
- `ConnectedComponents[g]` / `WeaklyConnectedComponents[g]` — components of the
  underlying undirected graph.
- `StronglyConnectedComponents[g]` — components following edge directions
  (Tarjan). For undirected graphs this coincides with the weak components.
- `FindSpanningTree[g]` — a spanning tree/forest as a graph (`VertexCount - 1`
  edges when connected); tree edges keep their original direction.
- `ConnectedGraphQ[g]` — `True` iff `g` is a single connected component.
- `VertexConnectivity[g]` — the minimum number of vertices whose removal
  disconnects `g` (`n-1` for `K_n`, `0` if already disconnected). Exact
  brute-force over vertex subsets, intended for small graphs.

```
FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4}], 1, 4]   (* {1, 2, 3, 4} *)
GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4}], 4, 1]      (* Infinity     *)
StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3}]]     (* {{1},{2},{3}} *)
VertexConnectivity[CycleGraph[5]]                           (* 2            *)

(* Weighted: the direct 1->4 edge (weight 10) loses to the longer,
   cheaper 1->2->3->4 route (weight 3). *)
FindShortestPath[
  Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}], 1, 4]
    (* {1, 2, 3, 4} *)
GraphDistance[
  Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}], 1, 4]
    (* 3 *)
```

Weighted `FindShortestPath`/`GraphDistance` use a plain O(V²) Dijkstra (no priority queue —
consistent with `VertexConnectivity`'s own small-graph exact-algorithm precedent above), and
fall back to unweighted BFS rather than erroring whenever a weight isn't usable for it (not
present, symbolic, or negative). No Bellman-Ford / negative-weight support.

## Structural predicates & ordering

Every `*Q` predicate here gives `False` — never unevaluated — for an argument
that is not a valid graph. All run in linear time on first use and `O(1)` when
repeated on the same graph (see *Performance model* above).

- `UndirectedGraphQ[g]` — `True` iff every edge is undirected. Edgeless graphs
  are undirected; a mixed graph is neither directed nor undirected.
- `EmptyGraphQ[g]` — `True` iff `g` has no edges (any number of vertices,
  including none).
- `CompleteGraphQ[g]` — `True` iff every ordered pair of distinct vertices
  `(u, v)` is joined by an edge usable from `u` to `v`: an undirected edge, or a
  directed `u -> v`. A complete directed graph therefore needs both directions.
  Graphs with 0 or 1 vertices are complete. `CompleteGraphQ[g, vlist]` tests the
  subgraph induced by `vlist` (`False` if some element is not a vertex of `g`;
  repeats are ignored; `{}` is complete).
- `BipartiteGraphQ[g]` — `True` iff the vertices split into two sets with every
  edge running between them. Edge direction is ignored; edgeless graphs are
  bipartite.
- `VertexQ[g, v]` — `True` iff `v` is a vertex of `g`, compared structurally (as
  by `SameQ`): the vertex `2` is not matched by `2.0`.
- `EdgeQ[g, e]` — `True` iff `e` is an edge of `g`. Accepts the constructor's
  sugar (`u -> v` is `DirectedEdge[u, v]`, `u <-> v` is `UndirectedEdge[u, v]`).
  An undirected edge matches in either orientation; direction must agree, so
  `u -> v` is not an edge of `Graph[{u, v}, {u <-> v}]`.
- `AcyclicGraphQ[g]` — `True` iff `g` has no cycle, where a cycle follows
  directed edges forwards and undirected edges either way, never reusing an
  edge. An undirected graph is acyclic iff it is a forest, a directed graph iff
  it is a DAG; `u -> v` with `v -> u` is a 2-cycle. Mixed graphs are decided
  exactly: undirected components are contracted (a repeated union-find join is
  an undirected cycle, and a directed edge inside one component closes a cycle
  through it), then the contracted digraph is checked with Kahn's algorithm.
- `TreeGraphQ[g]` — `True` iff `g` has at least one vertex, is connected, and has
  `VertexCount - 1` edges, i.e. is a tree when edge direction is ignored. An
  out-tree such as `1 -> 2, 1 -> 3` is a tree; a disconnected forest, the
  anti-parallel pair `1 -> 2, 2 -> 1`, and the null graph are not.
- `TopologicalSort[g]` — the vertices of a directed acyclic graph, ordered so
  that `u` precedes `v` for every edge `u -> v`. Kahn's algorithm; among the
  vertices ready at each step, the one earliest in `VertexList[g]` goes first,
  so the order is deterministic. `TopologicalSort[{v -> w, ...}]` uses the rules
  as the graph (built through `Graph`). Left unevaluated for a cyclic graph, a
  graph with any undirected edge, or a non-graph; an edgeless graph sorts to its
  `VertexList`.

```
UndirectedGraphQ[Graph[{1,2},{}]]                        (* True           *)
CompleteGraphQ[Graph[{1->2,2->1}]]                       (* True           *)
CompleteGraphQ[Graph[{1->2}]]                            (* False          *)
BipartiteGraphQ[CycleGraph[5]]                           (* False          *)
EdgeQ[CycleGraph[3], 2<->1]                              (* True           *)
EdgeQ[CycleGraph[3], 1->2]                               (* False          *)
AcyclicGraphQ[Graph[{1,2,3},{1<->2,2->3,3->1}]]          (* False          *)
TreeGraphQ[Graph[{1->2,1->3}]]                           (* True           *)
TopologicalSort[{1->3,1->4,2->1,2->4,3->4,5->2,5->3}]    (* {5,2,1,3,4}    *)
TopologicalSort[CycleGraph[3]]                           (* unevaluated    *)
```

## Visualization

- `GraphPlot[g]` — a `Graphics[...]` object drawing `g`: vertices are laid out
  on a circle (one `Disk` and one `Text` label each), edges are `Line`s. It
  renders through the standard graphics path (a window when `USE_GRAPHICS=1`,
  the text placeholder otherwise). Directed edges are drawn as plain lines in
  the MVP (no arrowheads yet); a force-directed layout is a future hook.

```
Head[GraphPlot[CycleGraph[8]]]                 (* Graphics *)
Count[GraphPlot[CompleteGraph[6]], _Line, Infinity]   (* 15 edges *)
```

## Graph algorithms: flows, matchings, cliques, Hamiltonian cycles, isomorphism, planarity

Implemented in `src/graph/galg_*.c` (declared in `src/graph/graph_algos.h`,
registered by `graph_algos_init()`). All read graphs through the
validated-graph memo, so on a graph built by `Graph[...]` they start from
pre-resolved integer endpoints. Every head is `Protected`, stays unevaluated on
a non-graph (the `*Q` predicates give `False`), and returns a fresh value.
Semantics and output forms were checked against Mathematica 15 by a randomized
differential test (`benchmarks/95-graph-algorithms/diff_mathematica.py`); where
the answer is unique the outputs are identical.

**Exactness policy.** Every NP-hard head (`FindVertexCover`,
`FindIndependentVertexSet`, `FindClique`, `FindKClique`, the Hamiltonian heads,
and the isomorphism search) returns a *proven* optimum / a complete answer, or
stays unevaluated when its deterministic node budget is exhausted. All of them
poll the `TimeConstrained` deadline, so `TimeConstrained[FindClique[g], 1]`
returns `$Aborted` rather than hanging. They never return a merely-good answer.

### Flows and cuts

- `FindMaximumFlow[g, s, t]` — the maximum flow value from `s` to `t` (`s`, `t`
  may be lists of sources/sinks; `s == t` gives `0`).
  `FindMaximumFlow[g, s, t, "prop"]` with `"FlowValue"`, `"FlowMatrix"` (a dense
  `n x n` matrix of edge flows, packed; Mathematica returns a `SparseArray`,
  which Mathilda does not have) or `"EdgeList"` (the edges carrying flow,
  oriented along it, in flow-matrix row order). Options `EdgeCapacity -> {c1,
  ...}` (EdgeList order) and `VertexCapacity -> {c1, ...}` (VertexList order;
  sources and sinks are uncapped). **Capacities ignore `EdgeWeight`**, exactly
  as Mathematica does; without `EdgeCapacity` every edge has capacity 1. An
  undirected edge carries flow either way.
- `FindMinimumCut[g]` — `{value, {part1, part2}}`, a global minimum cut
  weighted by `EdgeWeight` (else 1). For a graph with directed edges, the edges
  counted run from `part1` (the source side) to `part2`. Unevaluated for fewer
  than 2 vertices. Among equal cuts the shore without `VertexList[g][[1]]` is
  listed first for undirected graphs (Mathematica's tie choice is not
  reproducible; the value always agrees).
- `FindEdgeCut[g]` / `FindEdgeCut[g, s, t]` — the edges of a minimum (s-t) edge
  cut (weighted by `EdgeWeight`), in EdgeList order; the s-t cut is the one
  closest to `s`, as in Mathematica. Directed graphs use strong connectivity.
- `FindVertexCut[g]` / `FindVertexCut[g, s, t]` — a minimum vertex separator of
  the underlying undirected graph, in VertexList order (the s-t separator
  closest to `t`, `{}` for adjacent `s`, `t`). Mathematica's conventions: a
  complete undirected graph gives its first `n-1` vertices, a graph with a
  directed edge whose underlying graph is complete gives `{}`.
- `EdgeConnectivity[g]` / `EdgeConnectivity[g, s, t]` — the (s-t) edge
  connectivity, weighted by `EdgeWeight`; strong for directed graphs.

Numbers: integer capacities/weights give exact Integers; Rational or Real ones
give a Real (as Mathematica); `Infinity` is an allowed capacity; a negative or
symbolic one leaves the call unevaluated. Internally everything is int64: reals
are scaled by a common power of two, so the max flow of machine-real capacities
is computed exactly.

Algorithms: Dinic (BFS levels truncated at the sink, iterative blocking flow
with current-arc pointers) on a CSR residual network; Nagamochi-Ibaraki for
undirected global minimum cuts (maximum-adjacency orders that contract every
edge whose attachment reaches the current bound, not one pair per phase as in
Stoer-Wagner); `2(n-1)` bounded flows for directed global cuts; Even's
split-vertex network with Esfahanian-Hakimi pair selection for vertex
separators. `galg_vertex_connectivity(g, s, t)` exposes the same machinery as a
fast replacement for the brute-force `VertexConnectivity` in `connectivity.c`
(not wired in by this stream).

### Matchings, covers, independent sets

- `FindIndependentEdgeSet[g]` — a maximum matching (edge direction ignored),
  edges in EdgeList order. Karp-Sipser greedy start, then Hopcroft-Karp when
  the graph is bipartite and Edmonds' blossom algorithm otherwise (union-find
  blossom bases; failed searches retire their Hungarian trees).
- `FindEdgeCover[g]` — a minimum edge cover (maximum matching plus one edge per
  exposed vertex); `{}` when `g` has an isolated vertex, as in Mathematica.
- `FindVertexCover[g]` — a minimum vertex cover (complement of a maximum
  independent set), in VertexList order.
- `FindIndependentVertexSet[g]` — `{s}` with `s` a maximum independent set.
  `FindIndependentVertexSet[g, k]`, `[g, {k}]`, `[g, {kmin, kmax}]` and a third
  argument `n` / `All` enumerate MAXIMAL independent sets by size, exactly like
  the `FindClique` spec forms below (these forms are limited to 8192 vertices).
- `IndependentVertexSetQ[g, vs]`, `VertexCoverQ[g, vs]`,
  `IndependentEdgeSetQ[g, es]`, `EdgeCoverQ[g, es]` — membership predicates; an
  element that is not a vertex/edge of `g` gives `False` (an `UndirectedEdge`
  matches either orientation, a `DirectedEdge` only as given), repeats are
  allowed in the vertex forms.

Independent sets are exact: connected components are solved separately; a
component with average degree >= 8 (or density >= 0.05, up to 3000 vertices)
goes to the bitset maximum-clique search on its complement; sparser ones to
branch and reduce (degree-0/1 and triangle reductions, degree-2 folding,
domination, a greedy clique-cover upper bound, component splitting at every
node, and branching on a maximum-degree vertex with its mirrors).

Weighted graphs: these heads optimize cardinality, as the Wolfram
documentation states. (The differential test found Mathematica returning
non-maximum matchings and non-minimum edge covers on weighted graphs.)

### Cliques

- `FindClique[g]` — `{c}` with `c` a maximum clique (a directed graph needs
  edges both ways between clique members).
- `FindClique[g, k]` (largest maximal clique with at most `k` vertices; `k` may
  be `Infinity`), `FindClique[g, {k}]`, `FindClique[g, {kmin, kmax}]`, and
  `FindClique[g, spec, n]` / `FindClique[g, spec, All]` — maximal cliques by
  size, largest first; within one size Mathematica's order is reproduced
  (`FindClique[CycleGraph[5], Infinity, All]` gives
  `{{4,5},{3,4},{2,3},{1,5},{1,2}}`). `FindClique[g, 2]` is `{}` when every
  maximal clique is larger.
- `FindKClique[g, k]` — `{c}`, a largest set of vertices pairwise within
  distance `k` (maximum clique of the `k`-th power graph).

Maximum clique: bitset branch and bound with greedy-colouring bounds (Tomita's
MCQ/MCS family in San Segundo's BBMC form), vertices numbered by reverse
degeneracy order; graphs above 3000 vertices are decomposed by degeneracy
(each vertex's later neighbourhood is a small local problem, skipped when its
core number cannot beat the incumbent). Enumeration: pivoted Bron-Kerbosch on
bitsets with size pruning.

### Hamiltonian cycles and paths

- `FindHamiltonianCycle[g]` — `{c}` with `c` a Hamiltonian cycle as a list of
  edges (starting at `VertexList[g][[1]]`, each edge written in traversal
  order), or `{}`. `FindHamiltonianCycle[g, n]` / `[g, All]` give up to `n` /
  all of them, each once. The one-vertex graph gives `{{}}`; `K2` has none; a
  directed 2-cycle is one.
- `FindHamiltonianPath[g]` / `FindHamiltonianPath[g, s, t]` — a vertex list or
  `{}` (`{}` for the one-vertex graph, as in Mathematica).
- `HamiltonianGraphQ[g]` — `True` iff a Hamiltonian cycle exists (`True` for
  one vertex, `False` for no vertices).

Search over edge decisions with constraint propagation: each vertex needs
exactly two chosen edges (directed: one in, one out), chosen edges form path
fragments whose end-to-end links forbid short cycles, and the remaining graph
must stay biconnected (directed: strongly connected) at every node. A random
400-vertex cubic graph takes about 1 ms. Paths reduce to cycles through an
added vertex.

### Isomorphism and canonical forms

- `IsomorphicGraphQ[g1, g2, ...]` — `True` iff all the graphs are isomorphic
  (`False` if any argument is not a graph; one argument stays unevaluated, as
  in Mathematica).
- `FindGraphIsomorphism[g1, g2]` — `{assoc}` with `assoc` an Association
  `v -> image` over `VertexList[g1]` in order, or `{}`.
  `FindGraphIsomorphism[g1, g2, n]` / `[g1, g2, All]` give up to `n` / all of
  them (the list order is the engine's search order; Mathematica's differs).
  Two empty graphs give `{<||>}` (Mathematica gives `{}`, which contradicts its
  own `IsomorphicGraphQ` answer `True`).
- `CanonicalGraph[g]` — a graph on `1..n` with edges sorted, such that
  `CanonicalGraph[g] === CanonicalGraph[h]` iff `g` and `h` are isomorphic.
  The particular representative differs from Mathematica's (both are
  arbitrary); properties such as `EdgeWeight` are dropped, as in Mathematica.
- `GraphAutomorphismGroup[g]` — `PermutationGroup[{Cycles[...], ...}]` acting
  on vertex positions, given by a generating set (the group is Mathematica's;
  the generators may differ). Mathilda has no permutation-group functions, so
  the result is an inert expression.

Directed and mixed graphs are supported (Mathematica leaves mixed graphs
unevaluated). The reduction to the engine also folds self-loops into vertex
colours and turns an edge of multiplicity k > 1 into a coloured subdivision
vertex, so loops and multigraphs will work unchanged once `Graph` accepts
them (today's validator rejects both). Edge weights are ignored, as in
Mathematica.

Engine (`galg_iso.c`): individualization-refinement in the nauty / bliss /
Traces family. Hopcroft-style equitable refinement (U, out and in relations
counted separately, largest fragment skipped, every split decision a function
of cell positions, sizes and counts only) emits a 64-bit trace that is
compared event by event, so a branch that cannot match dies at its first
deviating split. The search tree is iterative with exact undo (no recursion
at depth n). Canonical labeling keeps the best leaf by (trace, relabeled
graph) with automorphism pruning (leaf and internal-node automorphisms,
jump-back, first-path orbits); the automorphisms found generate the whole
group. `g -> h` searches run in lockstep against `h`'s trace, and every map
returned is verified edge by edge. On an exactly 3-regular random graph with
10^4 vertices (where colour refinement alone learns nothing)
`IsomorphicGraphQ` takes about 4 ms against about 0.8 s in Mathematica 15; at
10^5 vertices about 0.2 s against about 60 s.
