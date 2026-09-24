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
again. The memo keeps at most 8 graphs alive past their last user reference. A fresh wrapper node around a memoized graph's own argument Lists (what the evaluator produces when re-evaluating a stored graph) re-keys the entry instead of re-validating.

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
  and, on a weighted graph, returns a machine real as the Wolfram Language does
  (weights `{5, 7}` give `12.`), identical to `GraphDistance[g, s]` and
  `GraphDistanceMatrix`. (Until v0.184 it returned an exact `Integer`/`Rational`,
  which disagreed with both Mathematica and the single-source form.)
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

## Editing, transforms, set operations, cycles and paths

Implemented in `src/graph/gops_*.c` (header `src/graph/graph_ops.h`). Every head
leaves a non-graph argument unevaluated (the `*Q` predicates give `False`).
Orders — of vertices, of edges, of cycle edges — follow Mathematica 15 unless a
deviation is listed.

**Performance.** Every edit is an integer pass over the validated-graph memo's
endpoint arrays: `O(V + E)` plus one hash per argument item, with vertex, edge
and weight nodes shared into the result. Results are registered with the memo
(`graph_memo_seed`) from the endpoint arrays already computed and stamped as
evaluated, so returning a graph costs one vertex hash per vertex — not a full
re-validation — and the first accessor on the result is a memo hit. The
cycle/path finders reuse a small per-graph cache of incidence lists (holding a
reference to the graph, like the memo), so repeating a query skips the CSR
build. At `10^5` vertices the edits run in 4–16 ms: 8–200x faster than
Mathematica 15 (warm and cold) and 20–60x faster than networkx; the transforms
`UndirectedGraph`, `DirectedGraph`, `FindEulerianCycle` and `FindPath` are
1.1–1.6x faster; `GraphComplement` of a 1000-cycle (498500 new edges) is at
parity, its time dominated by allocating and freeing edge expressions
(`benchmarks/93-graph-ops-editing`).

### Editing

- `VertexAdd[g, v]`, `VertexAdd[g, {v1, ...}]` — appends vertices not already
  present (repeats ignored). A list always means a list of vertices.
- `VertexDelete[g, v | {v1, ...} | patt]` — removes vertices and their incident
  edges; every listed vertex must exist (else unevaluated). A pattern removes the
  matching vertices. Orders and weights of the survivors are kept.
- `EdgeAdd[g, e | {e1, ...}]` — appends edges; endpoints not in `g` become new
  vertices (appended in order). `u -> v` takes the graph's kind: undirected in an
  undirected (or edgeless) graph, directed otherwise; `DirectedEdge` is always
  directed. A new edge has weight 1 in a weighted graph.
- `EdgeDelete[g, e | {e1, ...} | patt]` — removes edges; each must exist (an
  undirected edge matches either orientation; `1 -> 2` is not an edge of an
  undirected graph). Orders and remaining weights are kept.
- `Subgraph[g, {v1, ...} | patt]` — the induced subgraph. Vertices in the given
  order (non-vertices ignored, repeats dropped); edges emitted at their later
  endpoint in that order, following each vertex's incidence order (out-edges and
  undirected edges, then in-edges). Weights kept. An edge list (edge-induced
  subgraph) is not supported (unevaluated).
- `NeighborhoodGraph[g, v | {v1, ...}, k]` — subgraph induced by the vertices
  within distance `k` (default 1) of the centres, direction ignored. Vertices:
  the centres, then each centre's new vertices in `VertexList` order; edges as
  for `Subgraph`. Non-vertex centres are ignored. `k = Infinity` is accepted
  (Mathematica leaves it unevaluated).
- `VertexReplace[g, rules]` — renames vertices by `Replace` (patterns and
  `RuleDelayed` work); vertices mapped together merge.
- `EdgeRules[g]` — edges as `u -> v` rules.
- `VertexIndex[g, v]`, `EdgeIndex[g, e]` — 1-based positions (list arguments
  give lists). `EdgeIndex` matches `u -> v` against an undirected edge of an
  undirected graph, as Mathematica.
- `IndexGraph[g]`, `IndexGraph[g, r]` — vertices renamed `r, r+1, ...`
  (`r = 1`); weights kept.

Deviation: Mathilda graphs are simple, so an edit whose result would have a
self-loop or parallel edges — `EdgeAdd` of an existing edge, `VertexReplace`
merging two adjacent vertices — is left unevaluated (Mathematica returns a
multigraph).

### Transforms

- `GraphComplement[g]` — undirected `g`: every non-adjacent pair `i < j`;
  directed or mixed `g`: every ordered pair with no edge usable from `i` to `j`,
  as a directed edge. Row-major `VertexList` order; weights dropped.
- `ReverseGraph[g]` — directed edges reversed; undirected edges, order and
  weights kept.
- `UndirectedGraph[g]` — `u -> v` and `v -> u` merge into one edge whose weight
  is the sum (`Plus`) of theirs; edges oriented and ordered by `VertexList`
  position (upper triangle, row-major). An undirected `g` is returned unchanged.
- `DirectedGraph[g]` — each `u <-> v` becomes `u -> v, v -> u` in place (weight
  duplicated). `DirectedGraph[g, "Acyclic"]` orients each undirected edge from
  the earlier to the later vertex in `VertexList` (a DAG for undirected `g`) and
  sorts by (tail, head) position; a mixed `g` keeps its edge order. Other
  methods (`"Random"`, ...) are left unevaluated.
- `LineGraph[g]` — vertices `1..m` (EdgeList positions). Undirected: `j <-> i`
  for `i < j` sharing an endpoint, listed by `j`, then shared endpoint, then `i`.
  Directed: `i -> j` when edge `i` ends where edge `j` starts, in `(i, j)` order.
  Mixed graphs are unevaluated. Deviation: Mathematica numbers directed
  line-graph vertices in a traversal order of its own; Mathilda always uses
  EdgeList position (an isomorphic graph).

### Set operations

- `GraphUnion[g1, g2, ...]` — vertices: the union in canonical order. Edges: the
  distinct edges (an undirected edge equals its reversal) — all undirected:
  first-appearance order, each oriented by canonical vertex order; all directed:
  first-appearance order; mixed: canonical (`Sort`) order. `GraphUnion[g]` is `g`.
- `GraphIntersection[g1, g2, ...]` — vertices: the union (canonical order);
  edges: those of `g1` present in every graph, canonical order.
- `GraphDifference[g1, g2]` — vertices: the union (canonical order); edges:
  those of `g1` not in `g2`, canonical order.
- `GraphDisjointUnion[g1, g2, ...]` — vertices relabelled `1..n` (`g1`'s first),
  edges translated in order; `GraphDisjointUnion[g]` is `g`.

Weights are dropped by all four, as in Mathematica. Implementation: vertex lists
equal to the first graph's are mapped by an `O(V)` elementwise `SameQ` check (no
hashing); otherwise through one hash index over the union. The union is sorted
with `expr_compare` only when not already sorted (machine integers are sorted as
such), and edge keys over result positions are deduplicated in an integer hash
set and ordered by stable counting sorts.

### Predicates

- `SimpleGraphQ[g]`, `LoopFreeGraphQ[g]` — `True` for every valid graph.
- `MixedGraphQ[g]` — both directed and undirected edges.
- `WeightedGraphQ[g]`, `EdgeWeightedGraphQ[g]` — `g` carries `EdgeWeight`.
- `PathGraphQ[g]` — Mathematica's definition: at least one vertex, connected,
  and every degree `<= 2` (undirected) or every in/out-degree `<= 1` (directed);
  so cycles count (`PathGraphQ[CycleGraph[3]]` is `True`, as in Mathematica).
  Mixed graphs are never paths.
- `EulerianGraphQ[g]` — all degrees even (undirected) / in = out (directed) and
  all edges in one connected component; edgeless graphs with a vertex are
  Eulerian, the null graph is not. Mixed graphs are left unevaluated.

### Cycles and paths

- `FindEulerianCycle[g]`, `FindEulerianCycle[g, 1]` — `{cycle}` as a list of
  edges, `{}` if none, `{{}}` for an edgeless graph with a vertex. Hierholzer's
  algorithm, iterative, linear: starts at the first vertex with an edge, takes
  edges in EdgeList order, reports undirected cycles in pop order and directed
  ones forwards; undirected edges are written in the direction walked. This
  reproduces Mathematica's cycle in most cases (not all). `n > 1` / `All` and
  mixed graphs are left unevaluated.
- `FindCycle[g]` — `{cycle}` or `{}`, in linear time, by Mathematica's own
  search order (a stack DFS that scans a vertex's edges when visiting it), so the
  reported cycle matches Mathematica's. `FindCycle[g, k]` (length `<= k`, `k` may
  be `Infinity`), `FindCycle[g, {k}]`, `FindCycle[g, {kmin, kmax}]`,
  `FindCycle[g, kspec, n]` (`n` or `All`) and `FindCycle[{g, v}, ...]` (cycles
  through `v`; the plain form by BFS from `v`, linear) are supported. Each cycle
  is reported once. The length-bounded / enumerating forms backtrack from each
  vertex as the cycle's lowest vertex — exponential in the worst case (a cycle of
  exact length `n` is a Hamiltonian cycle) — poll `TimeConstrained`, and give up
  (unevaluated) after `5*10^7` steps. Their choice and order of cycles is
  Mathilda's own (Mathematica's differs). Cycle length counts edges; undirected
  cycles have length `>= 3`, directed `>= 2`. Mixed graphs, and weighted graphs
  with a length spec, are left unevaluated.
- `FindPath[g, s, t]` — `{path}` (a vertex list) or `{}`: the first path a DFS
  meets, neighbours in EdgeList order — linear time and the same path as
  Mathematica. `s == t` gives `{}`. `FindPath[g, s, t, kspec(, n)]` enumerates
  simple paths depth-first within the length bounds and reports the first `n`
  (or `All`) shortest first — Mathematica's order. Weighted graphs with a kspec
  are left unevaluated (Mathematica measures the kspec in total weight there).

```
VertexDelete[Graph[{1,2,3,4},{1<->2,2<->3,3<->4},EdgeWeight->{5,6,7}], 2]
                                     (* Graph[{1,3,4},{3<->4},EdgeWeight->{7}] *)
EdgeList[EdgeAdd[CycleGraph[3], {1->4, 4<->5}]]  (* {1<->2,2<->3,3<->1,1<->4,4<->5} *)
EdgeList[Subgraph[Graph[{1,2,3,4},{3<->4,1<->2,2<->3}], {3,2,4}]]  (* {2<->3,3<->4} *)
VertexList[NeighborhoodGraph[PathGraph[Range[6]], {1,6}]]     (* {1,6,2,5} *)
EdgeList[GraphUnion[Graph[{2<->1}], Graph[{3<->4}]]]           (* {1<->2,3<->4} *)
EdgeList[UndirectedGraph[Graph[{1->2,2->1,2->3}]]]             (* {1<->2,2<->3} *)
FindEulerianCycle[CycleGraph[4]]            (* {{1<->4,4<->3,3<->2,2<->1}} *)
FindCycle[Graph[{1->2,2->3,3->4,4->2,3->1}]]   (* {{1->2,2->3,3->1}} *)
FindPath[CompleteGraph[4], 1, 4, 2, All]       (* {{1,4},{1,3,4},{1,2,4}} *)
```

## Distances, centralities, clustering and graph families

Implemented in `src/graph/gmet_*.c` (header `src/graph/graph_metrics.h`,
registered by `graph_metrics_init()` at the end of `graph_init()`). Every
convention below was checked against Mathematica 15 with `wolframscript` on
small graphs, including the undocumented ones (noted *reverse-engineered*).
Numeric vector/matrix results are **packed** (`NDArrayQ` is `True`) whenever
they are uniform machine numbers; a result containing `Infinity` or exact
rationals is an ordinary list, as in Wolfram. Edge direction is followed; an
undirected edge is usable both ways.

**Weights.** Heads marked *(w)* use `EdgeWeight` as edge lengths and then
answer in machine reals (Wolfram converts even integer weights); the others
ignore weights exactly as Wolfram does. A symbolic, complex or negative weight
leaves a *(w)* head unevaluated (Wolfram also refuses symbolic weights;
negative weights, which Wolfram routes to Bellman–Ford, are not supported).

### Distances

- `GraphDistanceMatrix[g]` *(w)* — all-pairs distances, rows/columns in
  `VertexList` order, `Infinity` when unreachable. Integers unweighted (packed
  when every pair is reachable), reals weighted. `GraphDistanceMatrix[g, d]`
  keeps distances `<= d` (others `Infinity`).
- `GraphDistance[g, s]` *(w)* — distances from `s` to every vertex (the new
  single-source form; `GraphDistance[g, s, t]` is unchanged). Weighted: reals,
  except the source's own entry, which is an exact `0` as in Wolfram.
- `VertexEccentricity[g, v]` *(w)* — unweighted: the largest distance to a
  vertex `v` reaches (finite on disconnected graphs); weighted: `Infinity` if
  `v` does not reach every vertex (Wolfram's weighted rule).
- `GraphDiameter`, `GraphRadius`, `GraphCenter`, `GraphPeriphery`,
  `MeanGraphDistance` *(w)* — unweighted and not strongly connected (not
  connected, if undirected): `Infinity` / `{}`. Weighted: taken over the
  weighted eccentricities, so a weighted digraph that is not strongly connected
  can have a finite radius and a non-empty center (*reverse-engineered*).
  `MeanGraphDistance` averages over ordered pairs of distinct vertices; exact
  when unweighted; `Infinity` when some pair is unreachable; unevaluated for a
  single vertex. Graphs with no vertices give diameter/radius `0`.
- `GraphDensity[g]` — `(directed edges + 2 undirected edges)/(n (n - 1))`,
  exact; unevaluated for `n < 2`.
- `KirchhoffMatrix[g]` — `D - A`, `D` = number of incident edges, `A` the
  (directed) adjacency matrix; weights ignored. **Deviation:** a dense packed
  Integer matrix (Wolfram returns a `SparseArray`; this is its `Normal`).

### Centralities

- `DegreeCentrality[g]`, `[g, "In"]`, `[g, "Out"]` — exact degrees. On a mixed
  graph an undirected edge counts once in each direction (so it adds 2 to the
  default total), as in Wolfram.
- `ClosenessCentrality[g]` *(w)* — `r/s`: `r` vertices reachable from `v`,
  `s` the sum of their distances; 0 if none.
- `EccentricityCentrality[g]` *(w)* — `1/e(v)` with `e` measured over the
  reachable vertices (also when weighted); 0 when `e(v) = 0`.
- `BetweennessCentrality[g]` — unnormalized; unordered pairs on undirected
  graphs, ordered pairs on directed ones; weights ignored. **Deviation:** mixed
  graphs are left unevaluated (Wolfram's values there match no standard
  definition, e.g. `{1<->2, 2->3}` gives vertex 2 a betweenness of 0).
- `EdgeBetweennessCentrality[g]` *(w)* — per edge in `EdgeList` order, summed
  over **ordered** pairs even on undirected graphs (the edge of a `K2` scores
  2); weighted ties within a relative `1e-12` share paths.
- `PageRankCentrality[g]`, `[g, a]` — `x = a P^T x + (1 - a)/n`, dangling
  vertices jump uniformly, `Total[x] = 1`, default `a = 0.85`, `0 <= a <= 1`.
  Converged to `1e-14` (Wolfram stops near `1e-9`, so they agree to ~9 digits).
- `EigenvectorCentrality[g]`, `[g, "In"]` (default), `[g, "Out"]` — per
  strongly connected component: each component `C` with `|C| > 1` gets its
  Perron vector scaled to total `(|C| - 1)/Σ(|C'| - 1)`; single-vertex
  components get 0, so a DAG gives all zeros (*reverse-engineered*: reproduces
  Wolfram exactly on every disconnected / non-strongly-connected case tried,
  e.g. `K4 ⊔ K3 ⊔ K2` totals 3:2:1).
- `KatzCentrality[g, a]`, `[g, a, b]` — `x = a A^T x + b` (`b` = 1, a number,
  or a list). Answered even beyond the convergence radius, like Wolfram (`K4`
  with `a = 1/2` gives `-2`; dense LAPACK solve, `n <= 4000`); a singular
  system is unevaluated. No edges or an exact `a = 0` return `b` exactly.
- `HITSCentrality[g]` — `{h, A.h}`; `h` is built like `EigenvectorCentrality`
  but for `A^T A`, whose blocks are the classes of vertices sharing an
  in-neighbour (*reverse-engineered*; reproduces Wolfram exactly, including
  degenerate spectra such as `PathGraph[5]` and the all-zero answer for mixed
  graphs whose classes are singletons).

### Clustering

All exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
Undirected: `t(v)` triangles at `v`; local `t(v)/C(d(v), 2)`, global
`3T/Σ C(d, 2)`, mean of the locals. Directed (*reverse-engineered*): triangles
are directed 3-cycles; local `c(v)/(in(v) out(v) - r(v))` with `r(v)` the
reciprocally linked neighbours.

- `GraphTriangleCount[g]`, `LocalClusteringCoefficient[g]` / `[g, v]`,
  `GlobalClusteringCoefficient[g]`, `MeanClusteringCoefficient[g]`.

### Graph families

Undirected on `1..N`; the edge list is the sorted list of pairs `{i, j}`,
`i < j` — identical to Wolfram's `EdgeList` for each family.

- `WheelGraph[n]` (hub 1; `n = 1` or `n >= 4` — 2 and 3 would be multigraphs),
  `HypercubeGraph[n]`, `GridGraph[{n1, ..., nk}]` (first coordinate fastest),
  `KaryTree[n]`, `KaryTree[n, k]`, `CompleteKaryTree[n]`, `CompleteKaryTree[n, k]`
  (`n` levels), `CirculantGraph[n, j]`, `CirculantGraph[n, {j1, ...}]`,
  `PetersenGraph[]`, `PetersenGraph[n, k]` (inner star `1..n`, outer cycle
  `n+1..2n`), `TuranGraph[n, k]` (larger parts first),
  `HararyGraph[k, n]` (`k >= 2`, `n > k`),
  `CompleteGraph[{n1, n2, ...}]` (complete multipartite; `CompleteGraph[{n}]`
  is `K_n`). `CompleteGraph` and `GraphDistance` are re-registered by wrappers
  that delegate their pre-existing forms to the original builtins.
- Options (`DirectedEdges`, layout options) are not supported.

### Algorithms and performance

- **Bit-parallel multi-source BFS** (MS-BFS): 256 sources advance together,
  one bit per source per vertex, so each BFS level walks the adjacency once
  for the whole batch; batches run on a pthread team (`MATHILDA_THREADS`
  builds; `MATHILDA_GRAPH_THREADS=1` forces serial). Used for all-pairs
  distances and for the per-source summary (reach, distance sum, eccentricity)
  from which closeness, eccentricity centrality and the diameter family reduce;
  that summary is cached per graph node, so a sequence of those heads on one
  graph pays for one all-pairs pass. Weighted: binary-heap Dijkstra per source.
- **Brandes** betweenness (O(nm); weighted O(nm + n² log n) for edges),
  sources spread over threads with per-thread accumulators.
- **Triangle listing** with an acyclic orientation (degree order, or index
  order when `maxdeg² <= 4m`), O(m^1.5) worst case; directed 3-cycles via
  two-bit arc flags per edge.
- **Spectral**: restarted Arnoldi (Krylov dimension 40, BLAS re-orthogonal-
  ization, LAPACK on the Hessenberg matrix) per block, residual-tested to
  `1e-13`; PageRank/Katz by (Jacobi) iteration with a correct stopping rule.
- Finished results are cached per `(head, graph node, other arguments)` in a
  16-slot cache holding references, the same soundness argument as the
  validated-graph memo.

Benchmarks: `benchmarks/94-graph-metrics` (warm and cold cases, vs Mathematica
and networkx).

```
GraphDistanceMatrix[Graph[{1->2, 2->3, 3->1, 3->4}]]
      (* {{0,1,2,3},{2,0,1,2},{1,2,0,1},{Infinity,Infinity,Infinity,0}} *)
BetweennessCentrality[StarGraph[5]]          (* {6., 0., 0., 0., 0.} *)
EdgeBetweennessCentrality[PathGraph[{1,2,3,4}]] (* {6., 8., 6.} *)
LocalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]] (* {1, 1, 1/3, 0} *)
EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]  (* {0.333333, 0.333333, 0.333333, 0.} *)
EdgeList[WheelGraph[5]]
      (* {1<->2,1<->3,1<->4,1<->5,2<->3,2<->5,3<->4,4<->5} *)
```
