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
- `DirectedGraphQ[g]` — `True` iff `g` is a valid graph whose edges are all
  directed.
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
  reproducible). Returns unevaluated if `m` exceeds `n(n-1)/2`.

```
EdgeCount[CompleteGraph[5]]      (* 10                        *)
EdgeList[CycleGraph[4]]          (* {1<->2, 2<->3, 3<->4, 4<->1} *)
VertexDegree[PathGraph[5]]       (* {1, 2, 2, 2, 1}           *)
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
