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
again. The memo keeps at most 8 graphs alive past their last user reference
(plus up to 4 more held by the incidence cache behind `FindCycle`/`FindPath`/
`FindEulerianCycle`, one per traversal mode); each is released as soon as its
slot is reused. A fresh wrapper node around a memoized graph's own argument Lists (what the evaluator produces when re-evaluating a stored graph) re-keys the entry instead of re-validating.

## Graph

- `Graph[v, e]`: a graph with vertex list `v` and edge list `e`.
- `Graph[e]`: derives the vertex set from the edges, in first-appearance order
  (directed by default).
- `Graph[v, e, EdgeWeight -> {w1, ..., wm}]`: a weighted graph — `wi` is the
  weight of `e[[i]]`, matched by position.

**Features**:
- `Protected`. A graph is a value: the constructor normalizes and validates its
  input and returns the canonical `Graph[List[verts], List[edges]]` (or, when
  weighted, `Graph[List[verts], List[edges], EdgeWeight -> List[weights]]`).
- Edge normalization: `u -> v` (`Rule`) and `DirectedEdge[u, v]` become
  `DirectedEdge[u, v]`; `u <-> v` (`TwoWayRule`) and `UndirectedEdge[u, v]`
  become `UndirectedEdge[u, v]`. Directed and undirected edges may be mixed.
- Malformed input is left unevaluated: self-loops, parallel/duplicate edges,
  3-argument edges, an edge endpoint absent from an explicit vertex list, or
  (for a weighted graph) an `EdgeWeight` list whose length doesn't match the
  edge list. Anti-parallel directed edges `u -> v` and `v -> u` are distinct and
  allowed.
- The weighted form requires the explicit-vertex form;
  `Graph[e, EdgeWeight -> {...}]` (derived vertices) is not accepted and stays
  unevaluated. A weight list whose length doesn't match `e` is malformed, like
  any other rejection above. Read the weights back with `EdgeWeight`.
- Printing: in standard output a graph shows a terse summary,
  `Graph[<n vertices, m edges>]`. `InputForm` and `FullForm` print the literal
  constructor, which round-trips through the parser.
- Validation is memoized per graph node (see the *Performance model* section of
  this file's preamble), so repeated queries on the same graph do not re-check
  it. `GraphQ` tests validity.

### Basic Examples

```mathematica
In[1]:= Graph[{1,2,3,4}, {1->2, 2->3, 3->4, 4->1}]
Out[1]= Graph[<4 vertices, 4 edges>]

In[2]:= InputForm[Graph[{1,2}, {1<->2}]]
Out[2]= Graph[{1, 2}, {1 <-> 2}]

In[3]:= InputForm[Graph[{1->2, 2->3, 3->1}]]
Out[3]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 3 -> 1}]
```

### Scope

```mathematica
In[4]:= InputForm[Graph[{a,b,c}, {DirectedEdge[a,b], UndirectedEdge[b,c]}]]
Out[4]= Graph[{a, b, c}, {a -> b, b <-> c}]

In[5]:= InputForm[Graph[{1,2,3}, {1->2, 2->3}, EdgeWeight -> {5, 7}]]
Out[5]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5, 7}]

In[6]:= FullForm[Graph[{1,2},{1<->2}]]
Out[6]= Graph[List[1, 2], List[UndirectedEdge[1, 2]]]
```

**Malformed input** is returned unevaluated (a self-loop, a duplicate edge, an
endpoint missing from the vertex list, a weight list of the wrong length):

```mathematica
In[7]:= Graph[{1,2}, {1->1}]
Out[7]= Graph[{1, 2}, {1 -> 1}]

In[8]:= Graph[{1,2}, {1->2, 1->2}]
Out[8]= Graph[{1, 2}, {1 -> 2, 1 -> 2}]

In[9]:= Graph[{1,2}, {1->3}]
Out[9]= Graph[{1, 2}, {1 -> 3}]

In[10]:= Graph[{1,2,3}, {1->2, 2->3}, EdgeWeight -> {5}]
Out[10]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5}]
```

## GraphQ

- `GraphQ[g]`: gives `True` if `g` is a valid graph, and `False` otherwise.

**Features**:
- `Protected`. Never left unevaluated: any non-graph gives `False`.
- A graph is valid when it is the canonical `Graph[List, List]` with every edge
  a 2-argument `DirectedEdge`/`UndirectedEdge`, no self-loops, no parallel
  edges, and every endpoint present in the vertex list (the same conditions the
  `Graph` constructor enforces).
- The null graph `Graph[{}, {}]` is valid.

```mathematica
In[1]:= GraphQ[Graph[{1,2}, {1->2}]]
Out[1]= True

In[2]:= GraphQ[Graph[{1}, {1->1}]]
Out[2]= False

In[3]:= GraphQ[5]
Out[3]= False

In[4]:= GraphQ[Graph[{}, {}]]
Out[4]= True

In[5]:= GraphQ[CycleGraph[4]]
Out[5]= True
```

## VertexList

- `VertexList[g]`: the vertices of `g`, in canonical order.

**Features**:
- `Protected`. The query/representation heads (`VertexList`, `EdgeList`,
  `VertexCount`, `EdgeCount`, `AdjacencyList`, `VertexDegree`,
  `VertexInDegree`, `VertexOutDegree`, `EdgeWeight`) are all thin readers over
  the canonical form and return unevaluated on a non-graph argument.
- Canonical order is the explicit vertex list, or first-appearance order when
  `Graph[e]` derived the vertices from the edges.

```mathematica
In[1]:= VertexList[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= {1, 2, 3, 4}

In[2]:= VertexList[Graph[{c->a, a->b}]]
Out[2]= {c, a, b}

In[3]:= VertexList[Graph[{}, {}]]
Out[3]= {}

In[4]:= VertexList[5]
Out[4]= VertexList[5]
```

## EdgeList

- `EdgeList[g]`: the edges of `g`, in canonical `DirectedEdge`/`UndirectedEdge`
  form.

**Features**:
- `Protected`. A thin reader over the canonical form; unevaluated on a non-graph
  (see `VertexList`).
- Edges print in operator form (`1 -> 2`, `2 <-> 3`), but are
  `DirectedEdge`/`UndirectedEdge` internally.

```mathematica
In[1]:= EdgeList[Graph[{1,2,3},{1->2, 2<->3}]]
Out[1]= {1 -> 2, 2 <-> 3}

In[2]:= EdgeList[CycleGraph[4]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}

In[3]:= EdgeList[Graph[{1,2},{}]]
Out[3]= {}

In[4]:= EdgeList[5]
Out[4]= EdgeList[5]
```

## VertexCount / EdgeCount

- `VertexCount[g]`: the number of vertices of `g`.
- `EdgeCount[g]`: the number of edges of `g`.

**Features**:
- `Protected`. Cardinalities read from the canonical form; unevaluated on a
  non-graph (see `VertexList`).

```mathematica
In[1]:= VertexCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= 4

In[2]:= EdgeCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[2]= 4

In[3]:= VertexCount[Graph[{}, {}]]
Out[3]= 0

In[4]:= EdgeCount[Graph[{1,2,3},{}]]
Out[4]= 0

In[5]:= EdgeCount[x]
Out[5]= EdgeCount[x]
```

## AdjacencyList

- `AdjacencyList[g]`: `{neighbors(v1), ...}`, one list per vertex in vertex
  order.
- `AdjacencyList[g, v]`: the neighbors of `v`.

**Features**:
- `Protected`. Directed edges contribute successors (`v -> u` makes `u` a
  neighbor of `v`, not the reverse); undirected edges go both ways.
- Unevaluated on a non-graph (see `VertexList`) or when `v` is not a vertex.

```mathematica
In[1]:= AdjacencyList[Graph[{1,2,3},{1<->2,2<->3}]]
Out[1]= {{2}, {1, 3}, {2}}

In[2]:= AdjacencyList[Graph[{1,2,3},{1<->2,2<->3}], 2]
Out[2]= {1, 3}

In[3]:= AdjacencyList[Graph[{1,2,3},{1->2,3->1}]]
Out[3]= {{2}, {}, {1}}

In[4]:= AdjacencyList[Graph[{1,2,3},{1->2,3->1}], 1]
Out[4]= {2}

In[5]:= AdjacencyList[5]
Out[5]= AdjacencyList[5]
```

## VertexDegree

- `VertexDegree[g]`: the total degree (number of incident edges) of each vertex,
  in vertex order.
- `VertexDegree[g, v]`: the degree of vertex `v`.

**Features**:
- `Protected`. Counts every incident edge regardless of direction; for the
  directed split see `VertexInDegree` / `VertexOutDegree`.
- Unevaluated on a non-graph (see `VertexList`) or when `v` is not a vertex.

```mathematica
In[1]:= VertexDegree[Graph[{1,2,3},{1<->2,2<->3}]]
Out[1]= {1, 2, 1}

In[2]:= VertexDegree[Graph[{1,2,3},{1<->2,2<->3}], 2]
Out[2]= 2

In[3]:= VertexDegree[Graph[{1,2,3},{1->2,1->3}]]
Out[3]= {2, 1, 1}

In[4]:= VertexDegree[Graph[{1,2},{}]]
Out[4]= {0, 0}
```

## VertexInDegree / VertexOutDegree

- `VertexInDegree[g]` / `VertexInDegree[g, v]`: in-degrees of all vertices, or
  of `v`.
- `VertexOutDegree[g]` / `VertexOutDegree[g, v]`: out-degrees of all vertices,
  or of `v`.

**Features**:
- `Protected`. A `DirectedEdge` adds to the source's out-degree and the target's
  in-degree; an `UndirectedEdge` adds to both the in- and out-degree of each
  endpoint.
- Unevaluated on a non-graph (see `VertexList`).

```mathematica
In[1]:= VertexInDegree[Graph[{1,2,3},{1->2,1->3}]]
Out[1]= {0, 1, 1}

In[2]:= VertexOutDegree[Graph[{1,2,3},{1->2,1->3}]]
Out[2]= {2, 0, 0}

In[3]:= VertexInDegree[Graph[{1,2,3},{1->2,1->3}], 2]
Out[3]= 1

In[4]:= VertexOutDegree[Graph[{1,2,3},{1->2,1->3}], 1]
Out[4]= 2

In[5]:= VertexInDegree[Graph[{1,2,3},{1<->2,2->3}]]
Out[5]= {1, 1, 1}

In[6]:= VertexOutDegree[Graph[{1,2,3},{1<->2,2->3}]]
Out[6]= {1, 2, 0}
```

## DirectedGraphQ

- `DirectedGraphQ[g]`: `True` iff `g` is a valid graph with at least one edge,
  all of them directed.

**Features**:
- `Protected`. `False` (never unevaluated) for a non-graph.
- An edgeless graph counts as undirected (as in the Wolfram Language), so
  `DirectedGraphQ` and `UndirectedGraphQ` are never both `True`. A mixed graph
  is neither.

```mathematica
In[1]:= DirectedGraphQ[Graph[{1->2,2->3}]]
Out[1]= True

In[2]:= DirectedGraphQ[CycleGraph[3]]
Out[2]= False

In[3]:= DirectedGraphQ[Graph[{1,2},{}]]
Out[3]= False

In[4]:= DirectedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]
Out[4]= False

In[5]:= DirectedGraphQ[5]
Out[5]= False
```

## UndirectedGraphQ

- `UndirectedGraphQ[g]`: `True` iff every edge of `g` is undirected.

**Features**:
- `Protected`. Edgeless graphs are undirected; a mixed graph is neither directed
  nor undirected (see `DirectedGraphQ`).
- Shared by all the structural predicates (`UndirectedGraphQ`, `EmptyGraphQ`,
  `CompleteGraphQ`, `BipartiteGraphQ`, `VertexQ`, `EdgeQ`, `AcyclicGraphQ`,
  `TreeGraphQ`): each gives `False` — never unevaluated — for an argument that
  is not a valid graph, and each runs in linear time on first use and `O(1)`
  when repeated on the same graph (see the *Performance model* section of this
  file's preamble).

```mathematica
In[1]:= UndirectedGraphQ[Graph[{1,2},{}]]
Out[1]= True

In[2]:= UndirectedGraphQ[CycleGraph[3]]
Out[2]= True

In[3]:= UndirectedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]
Out[3]= False

In[4]:= UndirectedGraphQ[x]
Out[4]= False
```

## EdgeWeight

- `EdgeWeight[g]`: the weights of `g`'s edges, in `EdgeList` order.
- `EdgeWeight -> {w1, ...}`: the `Graph` option that attaches the weights (see
  `Graph`).

**Features**:
- `Protected`. Defaults to all `1`s when `g` was built without an `EdgeWeight`
  option. Weights may be symbolic or exact; they are returned as given.
- Unevaluated on a non-graph (see `VertexList`).
- Weight-aware consumers: `WeightedAdjacencyMatrix`, `FindShortestPath` and
  `GraphDistance` (see `FindShortestPath`). The other search/computation heads
  in this section ignore weights.

```mathematica
In[1]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]
Out[1]= {5, 7}

In[2]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3}]]
Out[2]= {1, 1}

In[3]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{a,1/2}]]
Out[3]= {a, 1/2}

In[4]:= EdgeWeight[Graph[{1,2},{}]]
Out[4]= {}

In[5]:= EdgeWeight[5]
Out[5]= EdgeWeight[5]
```

## AdjacencyMatrix

- `AdjacencyMatrix[g]`: the dense 0/1 adjacency matrix of `g` (`n x n`,
  canonical vertex order).

**Features**:
- `Protected`. Symmetric for undirected graphs; entry `(i, j)` is `1` for a
  directed edge `vi -> vj`.
- Linear-algebra interop: the result is an ordinary matrix, so `Det`, `Tr`,
  `Eigenvalues`, `MatrixPower`, etc. apply directly. Inverse: `AdjacencyGraph`.
  Weighted variant: `WeightedAdjacencyMatrix`.
- Mathematica returns a `SparseArray`; Mathilda returns a dense list of lists.
- Unevaluated on a non-graph.

```mathematica
In[1]:= AdjacencyMatrix[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}

In[2]:= Det[AdjacencyMatrix[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]]
Out[2]= -1

In[3]:= AdjacencyMatrix[PathGraph[3]]
Out[3]= {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}

In[4]:= Eigenvalues[AdjacencyMatrix[CycleGraph[4]]]
Out[4]= {-2, 2, 0, 0}

In[5]:= AdjacencyMatrix[Graph[{},{}]]
Out[5]= {}
```

## IncidenceMatrix

- `IncidenceMatrix[g]`: the `|V| x |E|` incidence matrix of `g`.

**Features**:
- `Protected`. Rows follow `VertexList`, columns follow `EdgeList`. Undirected
  edges mark both endpoints with `1`; directed edges are oriented (`-1` at the
  tail, `+1` at the head).
- Unevaluated on a non-graph.

```mathematica
In[1]:= IncidenceMatrix[PathGraph[3]]
Out[1]= {{1, 0}, {1, 1}, {0, 1}}

In[2]:= IncidenceMatrix[Graph[{1,2,3},{1->2,2->3}]]
Out[2]= {{-1, 0}, {1, -1}, {0, 1}}

In[3]:= IncidenceMatrix[Graph[{1,2,3},{1->2,2<->3}]]
Out[3]= {{-1, 0}, {1, 1}, {0, 1}}

In[4]:= IncidenceMatrix[Graph[{1,2},{}]]
Out[4]= {{}, {}}
```

## AdjacencyGraph

- `AdjacencyGraph[m]`: the graph on vertices `1..n` whose adjacency matrix is
  the 0/1 matrix `m`.

**Features**:
- `Protected`. The inverse of `AdjacencyMatrix`: undirected if `m` is symmetric,
  else directed. `AdjacencyGraph[AdjacencyMatrix[g]]` reproduces `g`'s edges
  (as a set — edges are regenerated in row-major order, so `EdgeList` order and
  undirected-edge orientation may differ from `g`'s; the adjacency matrices
  agree).
- A non-square (ragged) matrix, or entries other than 0/1, leave the call
  unevaluated. A `1` on the diagonal is currently dropped silently rather than
  rejected (the graph model has no self-loops).

```mathematica
In[1]:= AdjacencyGraph[{{0,1},{1,0}}]
Out[1]= Graph[<2 vertices, 1 edge>]

In[2]:= EdgeList[AdjacencyGraph[{{0,1,0},{0,0,1},{1,0,0}}]]
Out[2]= {1 -> 2, 2 -> 3, 3 -> 1}

In[3]:= EdgeList[AdjacencyGraph[{{0,1,1},{1,0,0},{1,0,0}}]]
Out[3]= {1 <-> 2, 1 <-> 3}

In[4]:= AdjacencyMatrix[AdjacencyGraph[AdjacencyMatrix[CycleGraph[4]]]] == AdjacencyMatrix[CycleGraph[4]]
Out[4]= True

In[5]:= AdjacencyGraph[{{1,2},{3}}]
Out[5]= AdjacencyGraph[{{1, 2}, {3}}]
```

## WeightedAdjacencyMatrix

- `WeightedAdjacencyMatrix[g]`: like `AdjacencyMatrix[g]`, but each nonzero
  entry is the corresponding edge's weight instead of `1` (`0` where there is
  no edge).

**Features**:
- `Protected`. Equal to `AdjacencyMatrix[g]` exactly when `g` has no
  `EdgeWeight` (every weight defaults to `1`). Undirected edges put their weight
  in both symmetric positions.
- Mathematica returns a `SparseArray`; Mathilda returns a dense matrix.
- Unevaluated on a non-graph.

```mathematica
In[1]:= WeightedAdjacencyMatrix[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]
Out[1]= {{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}

In[2]:= WeightedAdjacencyMatrix[Graph[{1,2,3},{1<->2,2<->3},EdgeWeight->{5,7}]]
Out[2]= {{0, 5, 0}, {5, 0, 7}, {0, 7, 0}}

In[3]:= WeightedAdjacencyMatrix[CycleGraph[4]] == AdjacencyMatrix[CycleGraph[4]]
Out[3]= True

In[4]:= WeightedAdjacencyMatrix[5]
Out[4]= WeightedAdjacencyMatrix[5]
```

## CycleGraph

- `CycleGraph[n]`: the cycle on vertices `1..n`.

**Features**:
- `Protected`. Like all the graph generators (`CompleteGraph`, `CycleGraph`,
  `PathGraph`, `StarGraph`, `RandomGraph`), it builds a canonical graph with
  vertices `1..n` and undirected edges via the `Graph` constructor path.
- Related generator `CompleteGraph[n]` — `K_n`, with all `n(n-1)/2` edges.
- Small cases: `CycleGraph[2]` is the single edge `1 <-> 2` (no parallel edges),
  `CycleGraph[1]` is one isolated vertex and `CycleGraph[0]` the null graph. A
  symbolic `n` is left unevaluated.

```mathematica
In[1]:= EdgeList[CycleGraph[4]]
Out[1]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}

In[2]:= VertexCount[CycleGraph[10]]
Out[2]= 10

In[3]:= EdgeList[CycleGraph[2]]
Out[3]= {1 <-> 2}

In[4]:= EdgeCount[CompleteGraph[5]]
Out[4]= 10

In[5]:= CycleGraph[x]
Out[5]= CycleGraph[x]
```

## PathGraph

- `PathGraph[n]`: the path `1-2-...-n`.
- `PathGraph[{v1, ...}]`: the path through the given vertices, in order.

**Features**:
- `Protected`. Undirected edges, built through the `Graph` constructor (see
  `CycleGraph`). `PathGraph[1]` is a single vertex; a symbolic argument is left
  unevaluated.

```mathematica
In[1]:= VertexDegree[PathGraph[5]]
Out[1]= {1, 2, 2, 2, 1}

In[2]:= EdgeList[PathGraph[4]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 4}

In[3]:= EdgeList[PathGraph[{a,b,c}]]
Out[3]= {a <-> b, b <-> c}

In[4]:= PathGraph[1]
Out[4]= Graph[<1 vertex, 0 edges>]

In[5]:= PathGraph[x]
Out[5]= PathGraph[x]
```

## StarGraph

- `StarGraph[n]`: the star on `1..n` — the hub `1` joined to each of the `n-1`
  leaves `2..n`.

**Features**:
- `Protected`. Undirected edges, built through the `Graph` constructor (see
  `CycleGraph`). `StarGraph[1]` is a single vertex; a symbolic argument is left
  unevaluated.

```mathematica
In[1]:= EdgeList[StarGraph[4]]
Out[1]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[2]:= VertexDegree[StarGraph[5]]
Out[2]= {4, 1, 1, 1, 1}

In[3]:= StarGraph[1]
Out[3]= Graph[<1 vertex, 0 edges>]

In[4]:= StarGraph[x]
Out[4]= StarGraph[x]
```

## RandomGraph

- `RandomGraph[{n, m}]`: a random undirected graph with `n` vertices and `m`
  distinct edges.
- `RandomGraph[{n, m}, k]`: a list of `k` independently sampled such graphs.

**Features**:
- `Protected`. Uses the seeded system RNG, so `SeedRandom` makes it
  reproducible. Vertices are `1..n`, edges undirected (see `CycleGraph`).
- Returns unevaluated if `m` exceeds `n(n-1)/2`. For `n <= 1` the only possible
  graph is the edgeless one, which is what is returned.
- `k = 0` gives `{}`; `k = 1` gives a one-element list, **not** a bare `Graph`.
  A negative, non-integer, or symbolic `k` leaves the expression unevaluated,
  silently — the convention the count-taking `Random*` heads share.
- Algorithm: the `n(n-1)/2` candidate edges are never materialised. Each graph
  draws exactly what `RandomSample` over the row-major candidate list would (so
  seeded output is that of `RandomSample`), in `O(m)` time and memory per graph
  — `RandomGraph[{200000, 300000}]` takes about 80 ms. (Until v0.185 each
  element built the full `O(n^2)` candidate list, and leaked it.)

```mathematica
In[1]:= SeedRandom[42]; EdgeList[RandomGraph[{5, 4}]]
Out[1]= {3 <-> 5, 1 <-> 5, 2 <-> 5, 3 <-> 4}

In[2]:= Length[RandomGraph[{6, 5}, 3]]
Out[2]= 3

In[3]:= RandomGraph[{3, 4}]
Out[3]= RandomGraph[{3, 4}]

In[4]:= RandomGraph[{6, 5}, 0]
Out[4]= {}

In[5]:= Head[RandomGraph[{6, 5}, 1]]
Out[5]= List

In[6]:= RandomGraph[{6, 5}, -1]
Out[6]= RandomGraph[{6, 5}, -1]
```

## FindShortestPath

- `FindShortestPath[g, s, t]`: a shortest path from `s` to `t` as a vertex
  list; `{}` if `t` is unreachable.

**Features**:
- `Protected`. Shared by the search & computation heads (`FindShortestPath`,
  `GraphDistance`, `ConnectedComponents`, `WeaklyConnectedComponents`,
  `StronglyConnectedComponents`, `FindSpanningTree`, `ConnectedGraphQ`,
  `VertexConnectivity`): all build an integer-indexed adjacency on demand, and
  all but `FindShortestPath`/`GraphDistance` are unweighted.
- **Weight-aware**: if `g` carries an `EdgeWeight` and every weight is
  non-negative and numeric, uses Dijkstra (minimum total weight); otherwise
  (unweighted, a symbolic weight, or a negative weight present) uses unweighted
  BFS (minimum hop count), following edge direction for directed graphs either
  way.
- Weighted search is a plain `O(V^2)` Dijkstra (no priority queue — consistent
  with `VertexConnectivity`'s original small-graph exact-algorithm precedent),
  and it falls back to unweighted BFS rather than erroring whenever a weight
  isn't usable for it (not present, symbolic, or negative). There is no
  Bellman-Ford / negative-weight support.
- Unevaluated when `s` or `t` is not a vertex.
- Companion `GraphDistance[g, s, t]` gives the length/total weight of that path;
  `Infinity` if unreachable. Same weight-aware dispatch as `FindShortestPath`,
  and, on a weighted graph, returns a machine real as the Wolfram Language does
  (weights `{5, 7}` give `12.`), identical to `GraphDistance[g, s]` and
  `GraphDistanceMatrix`. (Until v0.186 it returned an exact
  `Integer`/`Rational`, which disagreed with both Mathematica and the
  single-source form.)

```mathematica
In[1]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4}], 1, 4]
Out[1]= {1, 2, 3, 4}

In[2]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4}], 4, 1]
Out[2]= {}

In[3]:= FindShortestPath[CycleGraph[6], 1, 4]
Out[3]= {1, 2, 3, 4}
```

**Weighted**: the direct `1 -> 4` edge (weight 10) loses to the longer, cheaper
`1 -> 2 -> 3 -> 4` route (weight 3). A negative weight falls back to hop count.

```mathematica
In[4]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}], 1, 4]
Out[4]= {1, 2, 3, 4}

In[5]:= FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,-10}], 1, 4]
Out[5]= {1, 4}

In[6]:= GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}], 1, 4]
Out[6]= 3.0

In[7]:= GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4}], 4, 1]
Out[7]= Infinity
```

## ConnectedComponents / WeaklyConnectedComponents

- `ConnectedComponents[g]`: the connected components of the underlying
  undirected graph, as lists of vertices.
- `WeaklyConnectedComponents[g]`: the same components.

**Features**:
- `Protected`. Edge direction is ignored, so both heads agree on every graph.
  Unweighted; see `FindShortestPath` for the shared search machinery.
- Mathematica's `ConnectedComponents` on a directed graph gives the *strongly*
  connected components; Mathilda's gives the weak ones (use
  `StronglyConnectedComponents` for the directed notion).
- The null graph has no components; unevaluated on a non-graph.

```mathematica
In[1]:= ConnectedComponents[Graph[{1,2,3,4,5},{1<->2,3<->4}]]
Out[1]= {{1, 2}, {3, 4}, {5}}

In[2]:= ConnectedComponents[Graph[{1,2,3},{1->2,3->2}]]
Out[2]= {{1, 2, 3}}

In[3]:= WeaklyConnectedComponents[Graph[{1,2,3,4},{1->2,3->2}]]
Out[3]= {{1, 2, 3}, {4}}

In[4]:= ConnectedComponents[Graph[{},{}]]
Out[4]= {}

In[5]:= ConnectedComponents[5]
Out[5]= ConnectedComponents[5]
```

## StronglyConnectedComponents

- `StronglyConnectedComponents[g]`: the components of `g` following edge
  directions.

**Features**:
- `Protected`. Tarjan's algorithm. For undirected graphs this coincides with the
  weak components (see `ConnectedComponents`).
- Unevaluated on a non-graph.

```mathematica
In[1]:= StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3}]]
Out[1]= {{1}, {2}, {3}}

In[2]:= StronglyConnectedComponents[Graph[{1,2,3,4},{1->2,2->1,2->3,3->4,4->3}]]
Out[2]= {{1, 2}, {3, 4}}

In[3]:= StronglyConnectedComponents[PathGraph[3]]
Out[3]= {{1, 2, 3}}
```

## FindSpanningTree

- `FindSpanningTree[g]`: a spanning tree (or forest) of `g`, as a graph.

**Features**:
- `Protected`. Has `VertexCount - 1` edges when `g` is connected; a
  disconnected graph gives a spanning forest. Tree edges keep their original
  direction. Unweighted (not a minimum-weight tree); see `FindShortestPath` for
  the shared machinery.
- Unevaluated on a non-graph.

```mathematica
In[1]:= EdgeList[FindSpanningTree[CycleGraph[4]]]
Out[1]= {1 <-> 2, 4 <-> 1, 2 <-> 3}

In[2]:= EdgeCount[FindSpanningTree[CompleteGraph[6]]]
Out[2]= 5

In[3]:= EdgeList[FindSpanningTree[Graph[{1,2,3,4},{1->2,2->3,3->1}]]]
Out[3]= {1 -> 2, 3 -> 1}

In[4]:= EdgeList[FindSpanningTree[Graph[{1,2,3,4},{1<->2,3<->4}]]]
Out[4]= {1 <-> 2, 3 <-> 4}
```

## ConnectedGraphQ

- `ConnectedGraphQ[g]`: `True` iff `g` is a single connected component.

**Features**:
- `Protected`. Connectivity of the underlying undirected graph (weak
  connectivity for directed graphs), matching `ConnectedComponents`.
- The null graph is not connected.
- Unlike the `*Q` structural predicates (see `UndirectedGraphQ`), a non-graph
  argument leaves `ConnectedGraphQ` unevaluated; Mathematica gives `False`.

```mathematica
In[1]:= ConnectedGraphQ[CycleGraph[5]]
Out[1]= True

In[2]:= ConnectedGraphQ[Graph[{1,2,3},{1<->2}]]
Out[2]= False

In[3]:= ConnectedGraphQ[Graph[{1,2,3},{1->2,3->2}]]
Out[3]= True

In[4]:= ConnectedGraphQ[Graph[{},{}]]
Out[4]= False

In[5]:= ConnectedGraphQ[5]
Out[5]= ConnectedGraphQ[5]
```

## VertexConnectivity

- `VertexConnectivity[g]`: the minimum number of vertices whose removal
  disconnects `g`.

**Features**:
- `Protected`. Gives `n-1` for `K_n`, and `0` if `g` is already disconnected.
- Algorithm: originally an exact brute-force search over vertex subsets,
  intended for small graphs. It now uses the max-flow machinery shared with
  `FindVertexCut` (via `galg_vertex_connectivity`: Even's split-vertex network
  with Esfahanian-Hakimi pair selection), replacing the former exponential
  subset search.
- Edge direction is currently ignored (the underlying undirected graph is
  used), so a directed cycle scores like an undirected one. Mathematica uses
  strong connectivity for directed graphs (a directed cycle has vertex
  connectivity 1).

```mathematica
In[1]:= VertexConnectivity[CycleGraph[5]]
Out[1]= 2

In[2]:= VertexConnectivity[CompleteGraph[5]]
Out[2]= 4

In[3]:= VertexConnectivity[StarGraph[5]]
Out[3]= 1

In[4]:= VertexConnectivity[Graph[{1,2,3,4},{1<->2,3<->4}]]
Out[4]= 0

In[5]:= VertexConnectivity[Graph[{1->2,2->3,3->4,4->1}]]
Out[5]= 2
```

## EmptyGraphQ

- `EmptyGraphQ[g]`: `True` iff `g` has no edges.

**Features**:
- `Protected`. Any number of vertices, including none. `False` for a non-graph
  (see `UndirectedGraphQ`).

```mathematica
In[1]:= EmptyGraphQ[Graph[{1,2,3},{}]]
Out[1]= True

In[2]:= EmptyGraphQ[Graph[{},{}]]
Out[2]= True

In[3]:= EmptyGraphQ[PathGraph[2]]
Out[3]= False

In[4]:= EmptyGraphQ[5]
Out[4]= False
```

## CompleteGraphQ

- `CompleteGraphQ[g]`: `True` iff every ordered pair of distinct vertices is
  joined by an edge usable from the first to the second.
- `CompleteGraphQ[g, vlist]`: tests the subgraph induced by `vlist`.

**Features**:
- `Protected`. A pair `(u, v)` is covered by an undirected edge or a directed
  `u -> v`, so a complete directed graph needs both directions.
- Graphs with 0 or 1 vertices are complete.
- `CompleteGraphQ[g, vlist]` gives `False` if some element is not a vertex of
  `g`; repeats are ignored; `{}` is complete.
- `False` for a non-graph (see `UndirectedGraphQ`).

```mathematica
In[1]:= CompleteGraphQ[Graph[{1->2,2->1}]]
Out[1]= True

In[2]:= CompleteGraphQ[Graph[{1->2}]]
Out[2]= False

In[3]:= CompleteGraphQ[CompleteGraph[5]]
Out[3]= True

In[4]:= CompleteGraphQ[CycleGraph[4], {1,2}]
Out[4]= True

In[5]:= CompleteGraphQ[CycleGraph[4], {1,3}]
Out[5]= False

In[6]:= CompleteGraphQ[CycleGraph[4], {1,9}]
Out[6]= False
```

## BipartiteGraphQ

- `BipartiteGraphQ[g]`: `True` iff the vertices split into two sets with every
  edge running between them.

**Features**:
- `Protected`. Edge direction is ignored; edgeless graphs are bipartite.
  `False` for a non-graph (see `UndirectedGraphQ`).

```mathematica
In[1]:= BipartiteGraphQ[CycleGraph[5]]
Out[1]= False

In[2]:= BipartiteGraphQ[CycleGraph[6]]
Out[2]= True

In[3]:= BipartiteGraphQ[Graph[{1,2,3},{}]]
Out[3]= True

In[4]:= BipartiteGraphQ[x]
Out[4]= False
```

## VertexQ

- `VertexQ[g, v]`: `True` iff `v` is a vertex of `g`.

**Features**:
- `Protected`. Vertices are compared structurally (as by `SameQ`): the vertex
  `2` is not matched by `2.0`. `False` for a non-graph (see
  `UndirectedGraphQ`).

```mathematica
In[1]:= VertexQ[CycleGraph[3], 2]
Out[1]= True

In[2]:= VertexQ[CycleGraph[3], 2.0]
Out[2]= False

In[3]:= VertexQ[CycleGraph[3], 7]
Out[3]= False

In[4]:= VertexQ[5, 1]
Out[4]= False
```

## EdgeQ

- `EdgeQ[g, e]`: `True` iff `e` is an edge of `g`.

**Features**:
- `Protected`. Accepts the constructor's sugar (`u -> v` is
  `DirectedEdge[u, v]`, `u <-> v` is `UndirectedEdge[u, v]`).
- An undirected edge matches in either orientation; direction must agree, so
  `u -> v` is not an edge of `Graph[{u, v}, {u <-> v}]`.
- `False` for a non-graph (see `UndirectedGraphQ`).

```mathematica
In[1]:= EdgeQ[CycleGraph[3], 2<->1]
Out[1]= True

In[2]:= EdgeQ[CycleGraph[3], 1->2]
Out[2]= False

In[3]:= EdgeQ[Graph[{1->2}], DirectedEdge[1,2]]
Out[3]= True

In[4]:= EdgeQ[Graph[{1->2}], 2->1]
Out[4]= False

In[5]:= EdgeQ[5, 1->2]
Out[5]= False
```

## AcyclicGraphQ

- `AcyclicGraphQ[g]`: `True` iff `g` has no cycle.

**Features**:
- `Protected`. A cycle follows directed edges forwards and undirected edges
  either way, never reusing an edge. An undirected graph is acyclic iff it is a
  forest, a directed graph iff it is a DAG; `u -> v` with `v -> u` is a 2-cycle.
- Mixed graphs are decided exactly: undirected components are contracted (a
  repeated union-find join is an undirected cycle, and a directed edge inside
  one component closes a cycle through it), then the contracted digraph is
  checked with Kahn's algorithm.
- `False` for a non-graph (see `UndirectedGraphQ`).

```mathematica
In[1]:= AcyclicGraphQ[Graph[{1,2,3},{1<->2,2->3,3->1}]]
Out[1]= False

In[2]:= AcyclicGraphQ[PathGraph[4]]
Out[2]= True

In[3]:= AcyclicGraphQ[CycleGraph[4]]
Out[3]= False

In[4]:= AcyclicGraphQ[Graph[{1->2,2->1}]]
Out[4]= False

In[5]:= AcyclicGraphQ[Graph[{1->2,1->3,2->4,3->4}]]
Out[5]= True
```

## TreeGraphQ

- `TreeGraphQ[g]`: `True` iff `g` is a tree when edge direction is ignored.

**Features**:
- `Protected`. Requires at least one vertex, connectivity, and exactly
  `VertexCount - 1` edges.
- An out-tree such as `1 -> 2, 1 -> 3` is a tree; a disconnected forest, the
  anti-parallel pair `1 -> 2, 2 -> 1`, and the null graph are not. A single
  vertex is a tree.
- `False` for a non-graph (see `UndirectedGraphQ`).

```mathematica
In[1]:= TreeGraphQ[Graph[{1->2,1->3}]]
Out[1]= True

In[2]:= TreeGraphQ[Graph[{1,2,3,4},{1<->2,3<->4}]]
Out[2]= False

In[3]:= TreeGraphQ[Graph[{1->2,2->1}]]
Out[3]= False

In[4]:= TreeGraphQ[Graph[{},{}]]
Out[4]= False

In[5]:= TreeGraphQ[Graph[{1},{}]]
Out[5]= True
```

## TopologicalSort

- `TopologicalSort[g]`: the vertices of a directed acyclic graph, ordered so
  that `u` precedes `v` for every edge `u -> v`.
- `TopologicalSort[{v -> w, ...}]`: uses the rules as the graph (built through
  `Graph`).

**Features**:
- `Protected`. Kahn's algorithm; among the vertices ready at each step, the one
  earliest in `VertexList[g]` goes first, so the order is deterministic.
- Left unevaluated for a cyclic graph, a graph with any undirected edge, or a
  non-graph. An edgeless graph sorts to its `VertexList`.

```mathematica
In[1]:= TopologicalSort[{1->3,1->4,2->1,2->4,3->4,5->2,5->3}]
Out[1]= {5, 2, 1, 3, 4}

In[2]:= TopologicalSort[Graph[{a,b,c},{c->a}]]
Out[2]= {b, c, a}

In[3]:= TopologicalSort[Graph[{3,1,2},{}]]
Out[3]= {3, 1, 2}

In[4]:= TopologicalSort[Graph[{1->2,2->3,3->1}]]
Out[4]= TopologicalSort[Graph[<3 vertices, 3 edges>]]

In[5]:= TopologicalSort[CycleGraph[3]]
Out[5]= TopologicalSort[Graph[<3 vertices, 3 edges>]]
```

## GraphPlot

- `GraphPlot[g]`: a `Graphics[...]` object drawing `g`.

**Features**:
- `Protected`. Vertices are laid out on a circle, each drawn as a `Disk`; edges
  are `Line`s. The specification calls for one `Text` label per vertex as well,
  but the current binary emits no `Text` primitives (see the example below).
- Renders through the standard graphics path (a window when `USE_GRAPHICS=1`,
  the text placeholder otherwise).
- MVP limitations: directed edges are drawn as plain lines (no arrowheads yet);
  a force-directed layout is a future hook. Mathematica's `GraphPlot` uses a
  spring-electrical layout.
- Unevaluated on a non-graph.

```mathematica
In[1]:= Head[GraphPlot[CycleGraph[8]]]
Out[1]= Graphics

In[2]:= Count[GraphPlot[CompleteGraph[6]], _Line, Infinity]
Out[2]= 15

In[3]:= Count[GraphPlot[CycleGraph[5]], _Disk, Infinity]
Out[3]= 5

In[4]:= Count[GraphPlot[CycleGraph[5]], _Text, Infinity]
Out[4]= 0

In[5]:= GraphPlot[5]
Out[5]= GraphPlot[5]
```

## FindVertexColoring

- `FindVertexColoring[g]`: a minimal vertex colouring of `g`, as a list of
  integer colours in `VertexList` order.

**Features**:
- `Protected`. The number of distinct colours equals the chromatic number, and
  no edge joins two vertices of equal colour. Edge direction is ignored.
- Minimality is proven by exact search (Wolfram's `"BacktrackingDS"` method)
  seeded by DSATUR upper and clique lower bounds.
- Exact colouring is NP-hard, so there are two guards, and exceeding either
  returns the expression unevaluated — never a valid-but-larger colouring:
  graphs of more than 128 vertices are refused outright, and the search gives
  up after 8 million nodes (on the order of 100 seconds). The limit is a node
  count rather than a clock, so the answer does not depend on host speed.
- To bound how long a call may take, wrap it in `TimeConstrained` — that is the
  intended lever, and the search polls for the deadline. The node ceiling is a
  last-resort backstop for unattended runs.
- Cost is driven by density, not vertex count: a sparse 128-vertex graph
  answers instantly, while a dense one may exhaust the budget and refuse.
- `FindVertexColoring[g, {c1, ...}]` and `FindVertexColoring[g, l]`
  (Mathematica's colour-list forms) are not implemented and stay unevaluated.

```mathematica
In[1]:= FindVertexColoring[CycleGraph[5]]
Out[1]= {1, 2, 1, 2, 3}

In[2]:= FindVertexColoring[CycleGraph[4]]
Out[2]= {1, 2, 1, 2}

In[3]:= FindVertexColoring[CompleteGraph[4]]
Out[3]= {1, 2, 3, 4}

In[4]:= FindVertexColoring[Graph[{1,2,3},{}]]
Out[4]= {1, 1, 1}

In[5]:= FindVertexColoring[Graph[{},{}]]
Out[5]= {}

In[6]:= FindVertexColoring[5]
Out[6]= FindVertexColoring[5]
```

## VertexAdd

- `VertexAdd[g, v]`: appends the vertex `v` to the graph `g` if it is not already present.
- `VertexAdd[g, {v1, ...}]`: appends each listed vertex not already present (repeats ignored).

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- A list always means a list of vertices; to add a vertex that is itself a list,
  wrap it in another list (`VertexAdd[g, {{1, 2}}]`).
- Vertex and edge orders follow Mathematica 15. Existing weights are kept.
- Shared conventions of the graph-editing family (`VertexAdd`, `VertexDelete`,
  `EdgeAdd`, `EdgeDelete`, `Subgraph`, `NeighborhoodGraph`, `VertexReplace`,
  `EdgeRules`, `VertexIndex`, `EdgeIndex`, `IndexGraph`), the transforms, set
  operations, predicates and cycle/path finders: implemented in
  `src/graph/gops_*.c` (header `src/graph/graph_ops.h`). Every head leaves a
  non-graph argument unevaluated (the `*Q` predicates give `False`). Orders —
  of vertices, of edges, of cycle edges — follow Mathematica 15 unless a
  deviation is listed.
- Performance: every edit is an integer pass over the validated-graph memo's
  endpoint arrays: `O(V + E)` plus one hash per argument item, with vertex, edge
  and weight nodes shared into the result. Results are registered with the memo
  (`graph_memo_seed`) from the endpoint arrays already computed and stamped as
  evaluated, so returning a graph costs one vertex hash per vertex — not a full
  re-validation — and the first accessor on the result is a memo hit. At `10^5`
  vertices the edits run in 4–16 ms: 8–200x faster than Mathematica 15 (warm
  and cold) and 20–60x faster than networkx (`benchmarks/93-graph-ops-editing`).
- Deviation (whole editing family): Mathilda graphs are simple, so an edit whose
  result would have a self-loop or parallel edges — `EdgeAdd` of an existing
  edge, `VertexReplace` merging two adjacent vertices — is left unevaluated
  (Mathematica returns a multigraph).

```mathematica
In[1]:= VertexList[VertexAdd[CycleGraph[3], 4]]
Out[1]= {1, 2, 3, 4}

In[2]:= VertexList[VertexAdd[CycleGraph[3], {4, 5, 4, 1}]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= VertexList[VertexAdd[Graph[{a<->b}], {{1,2}}]]
Out[3]= {a, b, {1, 2}}

In[4]:= VertexAdd[x, 4]
Out[4]= VertexAdd[x, 4]
```

## VertexDelete

- `VertexDelete[g, v]`: removes the vertex `v` and its incident edges.
- `VertexDelete[g, {v1, ...}]`: removes every listed vertex and their incident edges.
- `VertexDelete[g, patt]`: removes the vertices matching the pattern `patt`.

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- Every listed vertex must exist; otherwise the call is left unevaluated.
- Orders and weights of the surviving vertices and edges are kept.
- `O(V + E)` integer pass over the memoized endpoint arrays plus one hash per
  argument item; the result is seeded into the graph memo (see `VertexAdd`).

```mathematica
In[1]:= InputForm[VertexDelete[Graph[{1,2,3,4},{1<->2,2<->3,3<->4},EdgeWeight->{5,6,7}], 2]]
Out[1]= Graph[{1, 3, 4}, {3 <-> 4}, EdgeWeight -> {7}]

In[2]:= EdgeList[VertexDelete[CycleGraph[5], {1, 3}]]
Out[2]= {4 <-> 5}

In[3]:= VertexList[VertexDelete[PathGraph[Range[6]], _?EvenQ]]
Out[3]= {1, 3, 5}

In[4]:= VertexDelete[CycleGraph[3], 7]
Out[4]= VertexDelete[Graph[<3 vertices, 3 edges>], 7]
```

## EdgeAdd

- `EdgeAdd[g, e]`: appends the edge `e` to `g`.
- `EdgeAdd[g, {e1, ...}]`: appends each listed edge in order.

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- Endpoints not in `g` become new vertices, appended in order.
- `u -> v` takes the graph's kind: undirected in an undirected (or edgeless)
  graph, directed otherwise; `DirectedEdge` is always directed.
- A new edge has weight 1 in a weighted graph.
- Deviation: Mathilda graphs are simple, so adding an edge that already exists
  (which would create parallel edges) or a self-loop is left unevaluated;
  Mathematica returns a multigraph.
- `O(V + E)` plus one hash per argument item; result seeded into the graph memo.

```mathematica
In[1]:= EdgeList[EdgeAdd[CycleGraph[3], {1->4, 4<->5}]]
Out[1]= {1 <-> 2, 2 <-> 3, 3 <-> 1, 1 <-> 4, 4 <-> 5}

In[2]:= EdgeList[EdgeAdd[Graph[{1->2}], 2->3]]
Out[2]= {1 -> 2, 2 -> 3}

In[3]:= EdgeList[EdgeAdd[Graph[{1->2}], DirectedEdge[3,1]]]
Out[3]= {1 -> 2, 3 -> 1}

In[4]:= EdgeAdd[CycleGraph[3], 1<->2]
Out[4]= EdgeAdd[Graph[<3 vertices, 3 edges>], TwoWayRule[1, 2]]
```

## EdgeDelete

- `EdgeDelete[g, e]`: removes the edge `e` from `g`.
- `EdgeDelete[g, {e1, ...}]`: removes each listed edge.
- `EdgeDelete[g, patt]`: removes the edges matching the pattern `patt`.

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- Each listed edge must exist, otherwise the call is left unevaluated. An
  undirected edge matches either orientation; `1 -> 2` is not an edge of an
  undirected graph.
- Orders and the remaining weights are kept.
- `O(V + E)` plus one hash per argument item; result seeded into the graph memo.

```mathematica
In[1]:= EdgeList[EdgeDelete[CycleGraph[4], 1<->2]]
Out[1]= {2 <-> 3, 3 <-> 4, 4 <-> 1}

In[2]:= EdgeList[EdgeDelete[CycleGraph[4], {2<->1, 3<->4}]]
Out[2]= {2 <-> 3, 4 <-> 1}

In[3]:= EdgeList[EdgeDelete[CycleGraph[4], _[1, _]]]
Out[3]= {2 <-> 3, 3 <-> 4, 4 <-> 1}

In[4]:= EdgeDelete[CycleGraph[4], 1->2]
Out[4]= EdgeDelete[Graph[<4 vertices, 4 edges>], 1 -> 2]

In[5]:= EdgeDelete[CycleGraph[4], 1<->3]
Out[5]= EdgeDelete[Graph[<4 vertices, 4 edges>], TwoWayRule[1, 3]]
```

## Subgraph

- `Subgraph[g, {v1, ...}]`: the subgraph of `g` induced by the listed vertices.
- `Subgraph[g, patt]`: the subgraph induced by the vertices matching `patt`.

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- Vertices appear in the given order; non-vertices are ignored and repeats dropped.
- Edges are emitted at their later endpoint in that order, following each
  vertex's incidence order (out-edges and undirected edges, then in-edges).
- Weights are kept.
- An edge list (the edge-induced subgraph form) is not supported; the call is
  left unevaluated.

```mathematica
In[1]:= EdgeList[Subgraph[Graph[{1,2,3,4},{3<->4,1<->2,2<->3}], {3,2,4}]]
Out[1]= {2 <-> 3, 3 <-> 4}

In[2]:= VertexList[Subgraph[CycleGraph[5], {3,9,1,3,2}]]
Out[2]= {3, 1, 2}

In[3]:= EdgeList[Subgraph[CycleGraph[6], {1,2,3}]]
Out[3]= {1 <-> 2, 2 <-> 3}

In[4]:= EdgeList[Subgraph[CycleGraph[6], _?OddQ]]
Out[4]= {}

In[5]:= Subgraph[CycleGraph[4], {1<->2}]
Out[5]= Subgraph[Graph[<4 vertices, 4 edges>], {TwoWayRule[1, 2]}]
```

## NeighborhoodGraph

- `NeighborhoodGraph[g, v]`: the subgraph induced by `v` and the vertices within distance 1 of it.
- `NeighborhoodGraph[g, {v1, ...}]`: the subgraph induced by the vertices within distance 1 of any centre.
- `NeighborhoodGraph[g, v | {v1, ...}, k]`: uses distance `k` instead of 1.

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- Distance ignores edge direction.
- Vertices: the centres, then each centre's new vertices in `VertexList` order;
  edges are ordered as for `Subgraph`.
- Non-vertex centres are ignored.
- `k = Infinity` is accepted (Mathematica leaves it unevaluated).

```mathematica
In[1]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], {1,6}]]
Out[1]= {1, 6, 2, 5}

In[2]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 3, 2]]
Out[2]= {3, 1, 2, 4, 5}

In[3]:= EdgeList[NeighborhoodGraph[Graph[{1->2,3->1,2->4}], 1]]
Out[3]= {1 -> 2, 3 -> 1}

In[4]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 2, Infinity]]
Out[4]= {2, 1, 3, 4, 5, 6}

In[5]:= VertexList[NeighborhoodGraph[PathGraph[Range[6]], 9]]
Out[5]= {}
```

## VertexReplace

- `VertexReplace[g, rules]`: renames the vertices of `g` by applying `rules` with `Replace`.

**Features**:
- `Protected`. A non-graph first argument is left unevaluated.
- Patterns and `RuleDelayed` work, as in `Replace`.
- Vertices mapped to the same name merge.
- Deviation: Mathilda graphs are simple, so merging two adjacent vertices (which
  would create a self-loop) is left unevaluated; Mathematica returns a
  multigraph.

```mathematica
In[1]:= EdgeList[VertexReplace[CycleGraph[3], {1->a, 2->b}]]
Out[1]= {a <-> b, b <-> 3, 3 <-> a}

In[2]:= VertexList[VertexReplace[PathGraph[Range[4]], x_?EvenQ :> x^2]]
Out[2]= {1, 4, 3, 16}

In[3]:= EdgeList[VertexReplace[Graph[{1<->2, 3<->4}], {3->1}]]
Out[3]= {1 <-> 2, 1 <-> 4}

In[4]:= VertexReplace[PathGraph[Range[3]], 2->1]
Out[4]= VertexReplace[Graph[<3 vertices, 2 edges>], 2 -> 1]
```

## EdgeRules

- `EdgeRules[g]`: the edges of `g` as a list of `u -> v` rules.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Undirected and directed edges both become rules, in `EdgeList` order.

```mathematica
In[1]:= EdgeRules[CycleGraph[3]]
Out[1]= {1 -> 2, 2 -> 3, 3 -> 1}

In[2]:= EdgeRules[Graph[{1->2, 2->3}]]
Out[2]= {1 -> 2, 2 -> 3}

In[3]:= EdgeRules[5]
Out[3]= EdgeRules[5]
```

## VertexIndex / EdgeIndex

- `VertexIndex[g, v]`: the 1-based position of vertex `v` in `VertexList[g]`.
- `EdgeIndex[g, e]`: the 1-based position of edge `e` in `EdgeList[g]`.
- A list as the second argument gives a list of positions.

**Features**:
- `Protected`. A non-graph first argument, or a vertex/edge not in the graph, is
  left unevaluated.
- `EdgeIndex` matches `u -> v` against an undirected edge of an undirected
  graph, as Mathematica does; an undirected edge matches either orientation.

```mathematica
In[1]:= VertexIndex[Graph[{a<->b, b<->c}], b]
Out[1]= 2

In[2]:= VertexIndex[Graph[{a<->b, b<->c}], {c, a}]
Out[2]= {3, 1}

In[3]:= EdgeIndex[CycleGraph[4], 3<->2]
Out[3]= 2

In[4]:= EdgeIndex[CycleGraph[4], 3->4]
Out[4]= 3

In[5]:= EdgeIndex[CycleGraph[4], {1<->2, 4<->1}]
Out[5]= {1, 4}

In[6]:= EdgeIndex[Graph[{1->2}], 2->1]
Out[6]= EdgeIndex[Graph[<2 vertices, 1 edge>], 2 -> 1]
```

## IndexGraph

- `IndexGraph[g]`: renames the vertices of `g` to `1, 2, ...` in `VertexList` order.
- `IndexGraph[g, r]`: renames the vertices to `r, r+1, ...`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- The default start is `r = 1`. Edge order and weights are kept.

```mathematica
In[1]:= EdgeList[IndexGraph[Graph[{a<->b, b<->c}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= EdgeList[IndexGraph[Graph[{a<->b, b<->c}], 10]]
Out[2]= {10 <-> 11, 11 <-> 12}

In[3]:= IndexGraph[{1,2}]
Out[3]= IndexGraph[{1, 2}]
```

## GraphComplement

- `GraphComplement[g]`: the graph on the vertices of `g` whose edges are the pairs not adjacent in `g`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Undirected `g`: every non-adjacent pair `i < j`, as an undirected edge.
- Directed or mixed `g`: every ordered pair with no edge usable from `i` to `j`,
  as a directed edge.
- Edges in row-major `VertexList` order; weights are dropped.
- Performance: the complement of a 1000-cycle (498500 new edges) is at parity
  with Mathematica 15, its time dominated by allocating and freeing edge
  expressions (`benchmarks/93-graph-ops-editing`).

```mathematica
In[1]:= EdgeList[GraphComplement[CycleGraph[5]]]
Out[1]= {1 <-> 3, 1 <-> 4, 2 <-> 4, 2 <-> 5, 3 <-> 5}

In[2]:= EdgeList[GraphComplement[Graph[{1->2,2->3}]]]
Out[2]= {1 -> 3, 2 -> 1, 3 -> 1, 3 -> 2}

In[3]:= InputForm[GraphComplement[Graph[{1,2,3},{1<->2},EdgeWeight->{4}]]]
Out[3]= Graph[{1, 2, 3}, {1 <-> 3, 2 <-> 3}]

In[4]:= EdgeList[GraphComplement[CompleteGraph[4]]]
Out[4]= {}

In[5]:= GraphComplement[{1}]
Out[5]= GraphComplement[{1}]
```

## ReverseGraph

- `ReverseGraph[g]`: reverses every directed edge of `g`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Undirected edges are kept as they are; edge order and weights are kept.

```mathematica
In[1]:= EdgeList[ReverseGraph[Graph[{1->2,2->3,3<->4}]]]
Out[1]= {2 -> 1, 3 -> 2, 3 <-> 4}

In[2]:= InputForm[ReverseGraph[Graph[{1,2,3},{1->2,3->2},EdgeWeight->{4,5}]]]
Out[2]= Graph[{1, 2, 3}, {2 -> 1, 2 -> 3}, EdgeWeight -> {4, 5}]

In[3]:= EdgeList[ReverseGraph[CycleGraph[3]]]
Out[3]= {1 <-> 2, 2 <-> 3, 3 <-> 1}
```

## UndirectedGraph

- `UndirectedGraph[g]`: the undirected graph obtained by replacing each directed edge of `g` with an undirected one.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- `u -> v` and `v -> u` merge into one edge whose weight is the sum (`Plus`) of theirs.
- Edges are oriented and ordered by `VertexList` position (upper triangle, row-major).
- An undirected `g` is returned unchanged.
- Performance: 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

```mathematica
In[1]:= EdgeList[UndirectedGraph[Graph[{1->2,2->1,2->3}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= InputForm[UndirectedGraph[Graph[{1,2,3},{1->2,2->1,3->1},EdgeWeight->{2,3,4}]]]
Out[2]= Graph[{1, 2, 3}, {1 <-> 2, 1 <-> 3}, EdgeWeight -> {5, 4}]

In[3]:= EdgeList[UndirectedGraph[Graph[{3->1, 2->3}]]]
Out[3]= {3 <-> 1, 3 <-> 2}
```

## DirectedGraph

- `DirectedGraph[g]`: replaces each undirected edge `u <-> v` of `g` by the pair `u -> v, v -> u`.
- `DirectedGraph[g, "Acyclic"]`: orients each undirected edge from the earlier to the later vertex in `VertexList`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Plain form: the two directed edges replace the undirected one in place; the
  weight is duplicated onto both.
- `"Acyclic"`: gives a DAG for undirected `g`, with edges sorted by (tail, head)
  `VertexList` position; a mixed `g` keeps its edge order.
- Other methods (`"Random"`, ...) are left unevaluated.
- Performance: 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

```mathematica
In[1]:= EdgeList[DirectedGraph[PathGraph[Range[3]]]]
Out[1]= {1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2}

In[2]:= InputForm[DirectedGraph[Graph[{1,2},{1<->2},EdgeWeight->{7}]]]
Out[2]= Graph[{1, 2}, {1 -> 2, 2 -> 1}, EdgeWeight -> {7, 7}]

In[3]:= EdgeList[DirectedGraph[Graph[{3<->1, 2<->3, 1<->2}], "Acyclic"]]
Out[3]= {3 -> 1, 3 -> 2, 1 -> 2}

In[4]:= EdgeList[DirectedGraph[Graph[{1->2, 3<->1}], "Acyclic"]]
Out[4]= {1 -> 2, 1 -> 3}

In[5]:= DirectedGraph[CycleGraph[3], "Random"]
Out[5]= DirectedGraph[Graph[<3 vertices, 3 edges>], Random]
```

## LineGraph

- `LineGraph[g]`: the line graph of `g`, whose vertices are the edges of `g`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Vertices are `1..m`, the `EdgeList` positions.
- Undirected `g`: `j <-> i` for `i < j` sharing an endpoint, listed by `j`, then
  by shared endpoint, then by `i`.
- Directed `g`: `i -> j` when edge `i` ends where edge `j` starts, in `(i, j)` order.
- Mixed graphs are left unevaluated.
- Deviation: Mathematica numbers directed line-graph vertices in a traversal
  order of its own; Mathilda always uses `EdgeList` position (an isomorphic graph).

```mathematica
In[1]:= EdgeList[LineGraph[PathGraph[Range[4]]]]
Out[1]= {2 <-> 1, 3 <-> 2}

In[2]:= EdgeList[LineGraph[StarGraph[4]]]
Out[2]= {2 <-> 1, 3 <-> 1, 3 <-> 2}

In[3]:= EdgeList[LineGraph[Graph[{1->2,2->3,3->1}]]]
Out[3]= {1 -> 2, 2 -> 3, 3 -> 1}

In[4]:= LineGraph[Graph[{1->2,2<->3}]]
Out[4]= LineGraph[Graph[<3 vertices, 2 edges>]]
```

## GraphUnion / GraphIntersection / GraphDifference

- `GraphUnion[g1, g2, ...]`: the graph with the union of the vertices and of the edges.
- `GraphIntersection[g1, g2, ...]`: the union of the vertices, with the edges of `g1` present in every graph.
- `GraphDifference[g1, g2]`: the union of the vertices, with the edges of `g1` not in `g2`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Vertices are always the union of the inputs' vertices, in canonical order.
- `GraphUnion` edges: the distinct edges (an undirected edge equals its
  reversal). All undirected: first-appearance order, each oriented by canonical
  vertex order. All directed: first-appearance order. Mixed: canonical (`Sort`)
  order. `GraphUnion[g]` is `g`.
- `GraphIntersection` / `GraphDifference` edges are in canonical order.
  `GraphDifference` takes exactly two graphs.
- Weights are dropped, as in Mathematica (the one-argument `GraphUnion[g]`
  returns `g` itself, so its weights survive).
- Implementation: vertex lists equal to the first graph's are mapped by an
  `O(V)` elementwise `SameQ` check (no hashing); otherwise through one hash index
  over the union. The union is sorted with `expr_compare` only when not already
  sorted (machine integers are sorted as such), and edge keys over result
  positions are deduplicated in an integer hash set and ordered by stable
  counting sorts. The same machinery serves `GraphDisjointUnion`.

```mathematica
In[1]:= EdgeList[GraphUnion[Graph[{2<->1}], Graph[{3<->4}]]]
Out[1]= {1 <-> 2, 3 <-> 4}

In[2]:= EdgeList[GraphUnion[CycleGraph[3], Graph[{2<->1, 3<->4}]]]
Out[2]= {1 <-> 2, 2 <-> 3, 1 <-> 3, 3 <-> 4}

In[3]:= EdgeList[GraphUnion[Graph[{2->1}], Graph[{1<->3}]]]
Out[3]= {2 -> 1, 1 <-> 3}

In[4]:= EdgeList[GraphIntersection[CompleteGraph[4], CycleGraph[4], PathGraph[Range[4]]]]
Out[4]= {1 <-> 2, 2 <-> 3, 3 <-> 4}

In[5]:= EdgeList[GraphDifference[CompleteGraph[4], CycleGraph[4]]]
Out[5]= {1 <-> 3, 2 <-> 4}

In[6]:= GraphUnion[CycleGraph[3], 5]
Out[6]= GraphUnion[Graph[<3 vertices, 3 edges>], 5]
```

## GraphDisjointUnion

- `GraphDisjointUnion[g1, g2, ...]`: the disjoint union of the graphs, with vertices relabelled `1..n`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Vertices are relabelled `1..n`, `g1`'s first; edges are translated in order.
- `GraphDisjointUnion[g]` is `g` (not relabelled).
- Weights are dropped, as in Mathematica. Shares the set-operation machinery
  described under `GraphUnion`.

```mathematica
In[1]:= EdgeList[GraphDisjointUnion[Graph[{a<->b}], Graph[{a<->b, b<->c}]]]
Out[1]= {1 <-> 2, 3 <-> 4, 4 <-> 5}

In[2]:= VertexList[GraphDisjointUnion[CycleGraph[3], PathGraph[{x,y}]]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= InputForm[GraphDisjointUnion[Graph[{a<->b}]]]
Out[3]= Graph[{a, b}, {a <-> b}]
```

## SimpleGraphQ / LoopFreeGraphQ

- `SimpleGraphQ[g]`: `True` if `g` is a simple graph (no self-loops, no parallel edges).
- `LoopFreeGraphQ[g]`: `True` if `g` has no self-loops.

**Features**:
- `Protected`. A non-graph argument gives `False`.
- Mathilda graphs are always simple, so both give `True` for every valid graph.

```mathematica
In[1]:= {SimpleGraphQ[CycleGraph[3]], LoopFreeGraphQ[CycleGraph[3]], SimpleGraphQ[5], LoopFreeGraphQ[{1}]}
Out[1]= {True, True, False, False}
```

## MixedGraphQ

- `MixedGraphQ[g]`: `True` if `g` has both directed and undirected edges.

**Features**:
- `Protected`. A non-graph argument gives `False`.

```mathematica
In[1]:= {MixedGraphQ[Graph[{1->2, 2<->3}]], MixedGraphQ[Graph[{1->2}]], MixedGraphQ[CycleGraph[3]], MixedGraphQ[x]}
Out[1]= {True, False, False, False}
```

## WeightedGraphQ / EdgeWeightedGraphQ

- `WeightedGraphQ[g]`: `True` if `g` carries `EdgeWeight`.
- `EdgeWeightedGraphQ[g]`: `True` if `g` carries `EdgeWeight`.

**Features**:
- `Protected`. A non-graph argument gives `False`.

```mathematica
In[1]:= {WeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], WeightedGraphQ[CycleGraph[3]], EdgeWeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], EdgeWeightedGraphQ[CycleGraph[3]], WeightedGraphQ[x]}
Out[1]= {True, False, True, False, False}
```

## PathGraphQ

- `PathGraphQ[g]`: `True` if `g` is a path graph in Mathematica's sense.

**Features**:
- `Protected`. A non-graph argument gives `False`.
- Mathematica's definition: at least one vertex, connected, and every degree
  `<= 2` (undirected) or every in/out-degree `<= 1` (directed). So cycles count:
  `PathGraphQ[CycleGraph[3]]` is `True`, as in Mathematica.
- Mixed graphs are never paths.

```mathematica
In[1]:= {PathGraphQ[PathGraph[Range[4]]], PathGraphQ[CycleGraph[3]], PathGraphQ[StarGraph[4]], PathGraphQ[Graph[{1->2,2->3}]], PathGraphQ[Graph[{1->2,3->2}]]}
Out[1]= {True, True, False, True, False}

In[2]:= {PathGraphQ[Graph[{1->2,2<->3}]], PathGraphQ[Graph[{1<->2,3<->4}]], PathGraphQ[Graph[{1},{}]], PathGraphQ[Graph[{},{}]], PathGraphQ[x]}
Out[2]= {False, False, True, False, False}
```

## EulerianGraphQ

- `EulerianGraphQ[g]`: `True` if `g` has an Eulerian cycle.

**Features**:
- `Protected`. A non-graph argument gives `False`.
- Tests that all degrees are even (undirected) or in-degree equals out-degree
  (directed), and that all edges lie in one connected component (isolated
  vertices are allowed).
- An edgeless graph with at least one vertex is Eulerian; the null graph is not.
- Mixed graphs are left unevaluated.

```mathematica
In[1]:= {EulerianGraphQ[CycleGraph[4]], EulerianGraphQ[PathGraph[Range[3]]], EulerianGraphQ[Graph[{1->2,2->3,3->1}]], EulerianGraphQ[Graph[{1},{}]], EulerianGraphQ[Graph[{},{}]]}
Out[1]= {True, False, True, True, False}

In[2]:= {EulerianGraphQ[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]], EulerianGraphQ[Graph[{1,2,3,4},{1<->2,2<->3,3<->1}]], EulerianGraphQ[x]}
Out[2]= {False, True, False}

In[3]:= EulerianGraphQ[Graph[{1->2,2<->3}]]
Out[3]= EulerianGraphQ[Graph[<3 vertices, 2 edges>]]
```

## FindEulerianCycle

- `FindEulerianCycle[g]`: a list `{cycle}` containing one Eulerian cycle of `g` as a list of edges, or `{}` if none exists.
- `FindEulerianCycle[g, 1]`: the same.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Gives `{{}}` for an edgeless graph with a vertex.
- Hierholzer's algorithm, iterative, linear time: starts at the first vertex with
  an edge, takes edges in `EdgeList` order, reports undirected cycles in pop
  order and directed ones forwards; undirected edges are written in the
  direction walked. This reproduces Mathematica's cycle in most cases (not all).
- `n > 1` / `All` and mixed graphs are left unevaluated.
- The cycle/path finders reuse a small per-graph cache of incidence lists
  (holding a reference to the graph, like the memo), so repeating a query skips
  the CSR build. 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

```mathematica
In[1]:= FindEulerianCycle[CycleGraph[4]]
Out[1]= {{1 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}

In[2]:= FindEulerianCycle[Graph[{1->2,2->3,3->1}]]
Out[2]= {{1 -> 2, 2 -> 3, 3 -> 1}}

In[3]:= FindEulerianCycle[Graph[{1<->2,2<->3,3<->1,3<->4,4<->5,5<->3}]]
Out[3]= {{1 <-> 3, 3 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}

In[4]:= FindEulerianCycle[PathGraph[Range[3]]]
Out[4]= {}

In[5]:= FindEulerianCycle[Graph[{1},{}]]
Out[5]= {{}}

In[6]:= FindEulerianCycle[CycleGraph[3], All]
Out[6]= FindEulerianCycle[Graph[<3 vertices, 3 edges>], All]
```

## FindCycle

- `FindCycle[g]`: a list `{cycle}` containing one cycle of `g` as a list of edges, or `{}` if `g` is acyclic.
- `FindCycle[g, k]`: a cycle of length at most `k` (`k` may be `Infinity`).
- `FindCycle[g, {k}]`: a cycle of length exactly `k`.
- `FindCycle[g, {kmin, kmax}]`: a cycle with length between `kmin` and `kmax`.
- `FindCycle[g, kspec, n]`: at most `n` cycles (`n` may be `All`).
- `FindCycle[{g, v}, ...]`: cycles through the vertex `v`.

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Plain form: linear time, by Mathematica's own search order (a stack DFS that
  scans a vertex's edges when visiting it), so the reported cycle matches
  Mathematica's.
- `FindCycle[{g, v}]` (plain form) is found by BFS from `v`, in linear time.
- Each cycle is reported once. Cycle length counts edges; undirected cycles have
  length `>= 3`, directed `>= 2`.
- The length-bounded / enumerating forms backtrack from each vertex as the
  cycle's lowest vertex — exponential in the worst case (a cycle of exact length
  `n` is a Hamiltonian cycle) — poll `TimeConstrained`, and give up (unevaluated)
  after `5*10^7` steps. Their choice and order of cycles is Mathilda's own
  (Mathematica's differs).
- Mixed graphs, and weighted graphs with a length spec, are left unevaluated.
- Reuses the per-graph incidence-list cache shared with `FindPath` and
  `FindEulerianCycle`.

```mathematica
In[1]:= FindCycle[Graph[{1->2,2->3,3->4,4->2,3->1}]]
Out[1]= {{1 -> 2, 2 -> 3, 3 -> 1}}

In[2]:= FindCycle[PathGraph[Range[4]]]
Out[2]= {}

In[3]:= FindCycle[CompleteGraph[4], 3, All]
Out[3]= {{1 <-> 2, 2 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 4, 4 <-> 1}, {1 <-> 3, 3 <-> 4, 4 <-> 1}, {2 <-> 3, 3 <-> 4, 4 <-> 2}}

In[4]:= FindCycle[CompleteGraph[4], {3,4}, 2]
Out[4]= {{1 <-> 2, 2 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}

In[5]:= FindCycle[{Graph[{1<->2,2<->3,3<->1,3<->4,4<->5,5<->3}], 4}]
Out[5]= {{4 <-> 3, 3 <-> 5, 5 <-> 4}}

In[6]:= FindCycle[Graph[{1->2,2<->3}]]
Out[6]= FindCycle[Graph[<3 vertices, 2 edges>]]
```

## FindPath

- `FindPath[g, s, t]`: a list `{path}` containing one path from `s` to `t` as a list of vertices, or `{}` if none exists.
- `FindPath[g, s, t, kspec]`: a path whose length satisfies `kspec` (`k`, `{k}`, `{kmin, kmax}`, `Infinity`).
- `FindPath[g, s, t, kspec, n]`: at most `n` paths (`n` may be `All`).

**Features**:
- `Protected`. A non-graph argument is left unevaluated.
- Plain form: the first path a DFS meets, taking neighbours in `EdgeList` order
  — linear time and the same path as Mathematica.
- `s == t` gives `{}`.
- With a `kspec`, simple paths are enumerated depth-first within the length
  bounds and the first `n` (or `All`) are reported shortest first —
  Mathematica's order.
- Weighted graphs with a `kspec` are left unevaluated (Mathematica measures the
  `kspec` in total weight there).
- Reuses the per-graph incidence-list cache; 1.1–1.6x faster than Mathematica 15
  at `10^5` vertices (`benchmarks/93-graph-ops-editing`).

```mathematica
In[1]:= FindPath[CycleGraph[5], 1, 3]
Out[1]= {{1, 2, 3}}

In[2]:= FindPath[CompleteGraph[4], 1, 4, 2, All]
Out[2]= {{1, 4}, {1, 3, 4}, {1, 2, 4}}

In[3]:= FindPath[CycleGraph[5], 1, 3, Infinity, All]
Out[3]= {{1, 2, 3}, {1, 5, 4, 3}}

In[4]:= FindPath[Graph[{1->2,3->2}], 1, 3]
Out[4]= {}

In[5]:= FindPath[CycleGraph[5], 2, 2]
Out[5]= {}

In[6]:= FindPath[Graph[{1,2,3},{1<->2,2<->3,3<->1},EdgeWeight->{1,2,3}], 1, 3, 2]
Out[6]= FindPath[Graph[<3 vertices, 3 edges>], 1, 3, 2]
```

## GraphDistanceMatrix

- `GraphDistanceMatrix[g]`: the matrix of shortest-path distances between all pairs of vertices of `g`, rows and columns in `VertexList` order.
- `GraphDistanceMatrix[g, d]`: keeps only distances `<= d`; every other entry is `Infinity`.

**Features**:
- Part of the distances / centralities / clustering / graph-families module,
  implemented in `src/graph/gmet_*.c` (header `src/graph/graph_metrics.h`,
  registered by `graph_metrics_init()` at the end of `graph_init()`). Every
  convention of this module was checked against Mathematica 15 with
  `wolframscript` on small graphs, including the undocumented ones (noted
  *reverse-engineered* in the sections concerned).
- Numeric vector/matrix results of the module are **packed** (`NDArrayQ` is
  `True`) whenever they are uniform machine numbers; a result containing
  `Infinity` or exact rationals is an ordinary list, as in Wolfram.
- Edge direction is followed; an undirected edge is usable both ways.
  `Infinity` marks an unreachable pair.
- **Weights** *(w)*: uses `EdgeWeight` as edge lengths and then answers in
  machine reals (Wolfram converts even integer weights). Unweighted results are
  Integers (packed when every pair is reachable). A symbolic, complex or
  negative weight leaves the call unevaluated (Wolfram also refuses symbolic
  weights; negative weights, which Wolfram routes to Bellman–Ford, are not
  supported). The same weight rule applies to every head marked *(w)* in this
  module (`GraphDistance`, `VertexEccentricity`, the diameter family,
  `MeanGraphDistance`, `ClosenessCentrality`, `EccentricityCentrality`,
  `EdgeBetweennessCentrality`); the other heads ignore weights exactly as
  Wolfram does.
- **Algorithm — bit-parallel multi-source BFS** (MS-BFS): 256 sources advance
  together, one bit per source per vertex, so each BFS level walks the
  adjacency once for the whole batch; batches run on a pthread team
  (`MATHILDA_THREADS` builds; `MATHILDA_GRAPH_THREADS=1` forces serial). Used
  for all-pairs distances and for the per-source summary (reach, distance sum,
  eccentricity) from which closeness, eccentricity centrality and the diameter
  family reduce; that summary is cached per graph node, so a sequence of those
  heads on one graph pays for one all-pairs pass. Weighted: binary-heap
  Dijkstra per source.
- **Result cache**: finished results of every head in this module are cached
  per `(head, graph node, other arguments)` in a 16-slot cache holding
  references, the same soundness argument as the validated-graph memo.
- Benchmarks: `benchmarks/94-graph-metrics` (warm and cold cases, vs
  Mathematica and networkx).

```mathematica
In[1]:= GraphDistanceMatrix[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {{0, 1, 2, 3}, {2, 0, 1, 2}, {1, 2, 0, 1}, {Infinity, Infinity, Infinity, 0}}

In[2]:= GraphDistanceMatrix[PathGraph[4]]
Out[2]= {{0, 1, 2, 3}, {1, 0, 1, 2}, {2, 1, 0, 1}, {3, 2, 1, 0}}

In[3]:= NDArrayQ[GraphDistanceMatrix[PathGraph[4]]]
Out[3]= True

In[4]:= GraphDistanceMatrix[PathGraph[5], 2]
Out[4]= {{0, 1, 2, Infinity, Infinity}, {1, 0, 1, 2, Infinity}, {2, 1, 0, 1, 2}, {Infinity, 2, 1, 0, 1}, {Infinity, Infinity, 2, 1, 0}}

In[5]:= GraphDistanceMatrix[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[5]= {{0.0, 1.0, 2.0}, {1.0, 0.0, 1.0}, {2.0, 1.0, 0.0}}

In[6]:= GraphDistanceMatrix[x]
Out[6]= GraphDistanceMatrix[x]
```

## GraphDistance

- `GraphDistance[g, s, t]`: the length (or total weight) of a shortest path from `s` to `t` in `g`; `Infinity` if `t` is unreachable.
- `GraphDistance[g, s]`: the distances from `s` to every vertex, in `VertexList` order.

**Features**:
- *(w)* weight-aware. `GraphDistance[g, s]` is the single-source form added
  by the graph-metrics module (`src/graph/gmet_*.c`); `GraphDistance[g, s, t]`
  is unchanged. `GraphDistance` is re-registered by a wrapper that delegates
  its pre-existing form (`[g, s, t]`) to the original builtin.
- `[g, s, t]` uses the same weight-aware dispatch as `FindShortestPath`
  (Dijkstra when every weight is non-negative and numeric, BFS otherwise) and
  on a weighted graph returns a machine real, as the Wolfram Language does
  (weights `{5, 7}` give `12.`), identical to `GraphDistance[g, s]` and
  `GraphDistanceMatrix`.
- Single-source form, weighted: machine reals, except the source's own entry,
  which is an exact `0` as in Wolfram. A symbolic, complex or negative weight
  leaves the single-source form unevaluated.
- Unweighted single-source results are Integers, packed when every vertex is
  reachable. Edge direction is followed; an undirected edge is usable both
  ways.

```mathematica
In[1]:= GraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}], 1]
Out[1]= {0, 1, 2, 3}

In[2]:= GraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}], 4, 1]
Out[2]= Infinity

In[3]:= GraphDistance[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}], 1]
Out[3]= {0, 1.0, 2.0}

In[4]:= GraphDistance[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{5,7}], 1, 3]
Out[4]= 12.0

In[5]:= GraphDistance[Graph[{1,2},{1<->2}, EdgeWeight->{a}], 1]
Out[5]= GraphDistance[Graph[<2 vertices, 1 edge>], 1]
```

## VertexEccentricity

- `VertexEccentricity[g, v]`: the largest distance from `v` to a vertex of `g`.

**Features**:
- *(w)* weight-aware (machine reals when weighted; symbolic, complex or
  negative weights leave it unevaluated).
- Unweighted: the largest distance to a vertex that `v` reaches, so it is
  finite on disconnected graphs.
- Weighted: `Infinity` if `v` does not reach every vertex (Wolfram's weighted
  rule).
- Computed from the cached per-source MS-BFS summary (see
  `GraphDistanceMatrix`).

```mathematica
In[1]:= VertexEccentricity[PathGraph[5], 3]
Out[1]= 2

In[2]:= VertexEccentricity[Graph[{1->2, 2->3, 3->1, 3->4}], 1]
Out[2]= 3

In[3]:= VertexEccentricity[Graph[{1<->2,3<->4}], 1]
Out[3]= 1

In[4]:= VertexEccentricity[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}], 1]
Out[4]= 3.0

In[5]:= VertexEccentricity[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}], 2]
Out[5]= Infinity
```

## GraphDiameter / GraphRadius

- `GraphDiameter[g]`: the largest vertex eccentricity of `g`.
- `GraphRadius[g]`: the smallest vertex eccentricity of `g`.

**Features**:
- *(w)* weight-aware (machine reals when weighted).
- Unweighted and not strongly connected (not connected, if undirected):
  `Infinity`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a finite radius
  (*reverse-engineered*).
- Graphs with no vertices give diameter/radius `0`.
- Reduced from the cached MS-BFS per-source summary, so calling several of
  `GraphDiameter`, `GraphRadius`, `GraphCenter`, `GraphPeriphery`,
  `ClosenessCentrality`, `EccentricityCentrality` on one graph pays for one
  all-pairs pass.

```mathematica
In[1]:= GraphDiameter[PathGraph[5]]
Out[1]= 4

In[2]:= GraphRadius[PathGraph[5]]
Out[2]= 2

In[3]:= GraphDiameter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= Infinity

In[4]:= GraphRadius[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[4]= 3.0

In[5]:= GraphDiameter[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[5]= 2.0

In[6]:= GraphDiameter[Graph[{},{}]]
Out[6]= 0
```

## GraphCenter / GraphPeriphery

- `GraphCenter[g]`: the vertices of `g` whose eccentricity equals the radius.
- `GraphPeriphery[g]`: the vertices of `g` whose eccentricity equals the diameter.

**Features**:
- *(w)* weight-aware.
- Unweighted and not strongly connected (not connected, if undirected): `{}`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a non-empty center
  (*reverse-engineered*).
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

```mathematica
In[1]:= GraphCenter[PathGraph[5]]
Out[1]= {3}

In[2]:= GraphPeriphery[PathGraph[5]]
Out[2]= {1, 5}

In[3]:= GraphCenter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {}

In[4]:= GraphCenter[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[4]= {1}

In[5]:= GraphPeriphery[Graph[{1<->2,3<->4}]]
Out[5]= {}
```

## MeanGraphDistance

- `MeanGraphDistance[g]`: the mean of the distances between all ordered pairs of distinct vertices of `g`.

**Features**:
- *(w)* weight-aware (machine real when weighted).
- Averages over ordered pairs of distinct vertices; exact when unweighted.
- `Infinity` when some pair is unreachable (in particular, unweighted and not
  strongly connected).
- Unevaluated for a single vertex.

```mathematica
In[1]:= MeanGraphDistance[PathGraph[4]]
Out[1]= 5/3

In[2]:= MeanGraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= Infinity

In[3]:= MeanGraphDistance[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[3]= 1.33333

In[4]:= MeanGraphDistance[Graph[{1},{}]]
Out[4]= MeanGraphDistance[Graph[<1 vertex, 0 edges>]]
```

## GraphDensity

- `GraphDensity[g]`: the ratio of the number of edges of `g` to the maximum possible number of edges.

**Features**:
- Computed as `(directed edges + 2 undirected edges)/(n (n - 1))`, exact.
- Weights ignored.
- Unevaluated for `n < 2`.

```mathematica
In[1]:= GraphDensity[PathGraph[4]]
Out[1]= 1/2

In[2]:= GraphDensity[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= 1/3

In[3]:= GraphDensity[CompleteGraph[5]]
Out[3]= 1

In[4]:= GraphDensity[Graph[{1<->2, 2->3}]]
Out[4]= 1/2

In[5]:= GraphDensity[Graph[{1},{}]]
Out[5]= GraphDensity[Graph[<1 vertex, 0 edges>]]
```

## KirchhoffMatrix

- `KirchhoffMatrix[g]`: the Kirchhoff (Laplacian) matrix `D - A` of `g`.

**Features**:
- `D` is the diagonal matrix of the number of incident edges of each vertex,
  `A` the (directed) adjacency matrix; weights ignored.
- **Deviation:** returns a dense packed Integer matrix (Wolfram returns a
  `SparseArray`; this is its `Normal`).

```mathematica
In[1]:= KirchhoffMatrix[PathGraph[3]]
Out[1]= {{1, -1, 0}, {-1, 2, -1}, {0, -1, 1}}

In[2]:= KirchhoffMatrix[Graph[{1->2,2->3}]]
Out[2]= {{1, -1, 0}, {0, 2, -1}, {0, 0, 1}}

In[3]:= NDArrayQ[KirchhoffMatrix[PathGraph[3]]]
Out[3]= True

In[4]:= KirchhoffMatrix[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[4]= {{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}}

In[5]:= KirchhoffMatrix[5]
Out[5]= KirchhoffMatrix[5]
```

## DegreeCentrality

- `DegreeCentrality[g]`: the list of vertex degrees of `g`.
- `DegreeCentrality[g, "In"]`: the in-degrees.
- `DegreeCentrality[g, "Out"]`: the out-degrees.

**Features**:
- Exact degrees, in `VertexList` order; weights ignored.
- On a mixed graph an undirected edge counts once in each direction (so it
  adds 2 to the default total), as in Wolfram.
- An unknown direction string leaves the call unevaluated.

```mathematica
In[1]:= DegreeCentrality[StarGraph[5]]
Out[1]= {4, 1, 1, 1, 1}

In[2]:= DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "In"]
Out[2]= {1, 1, 1, 1}

In[3]:= DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Out"]
Out[3]= {1, 1, 2, 0}

In[4]:= DegreeCentrality[Graph[{1<->2, 2->3}]]
Out[4]= {2, 3, 1}

In[5]:= DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Foo"]
Out[5]= DegreeCentrality[Graph[<4 vertices, 4 edges>], Foo]
```

## ClosenessCentrality

- `ClosenessCentrality[g]`: the closeness centrality of each vertex of `g`.

**Features**:
- *(w)* weight-aware; machine reals, packed.
- For a vertex `v`: `r/s`, with `r` the number of vertices reachable from `v`
  and `s` the sum of their distances; 0 if none are reachable.
- Reduced from the cached MS-BFS per-source summary (weighted: Dijkstra per
  source); see `GraphDistanceMatrix`.

```mathematica
In[1]:= ClosenessCentrality[PathGraph[4]]
Out[1]= {0.5, 0.75, 0.75, 0.5}

In[2]:= ClosenessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.5, 0.6, 0.75, 0.0}

In[3]:= ClosenessCentrality[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[3]= {0.5, 0.5, 0.0}

In[4]:= ClosenessCentrality[Graph[{1,2,3},{}]]
Out[4]= {0.0, 0.0, 0.0}
```

## EccentricityCentrality

- `EccentricityCentrality[g]`: the eccentricity centrality `1/e(v)` of each vertex of `g`.

**Features**:
- *(w)* weight-aware; machine reals.
- `e(v)` is measured over the vertices reachable from `v` (also when
  weighted); the centrality is 0 when `e(v) = 0`.
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

```mathematica
In[1]:= EccentricityCentrality[PathGraph[5]]
Out[1]= {0.25, 0.333333, 0.5, 0.333333, 0.25}

In[2]:= EccentricityCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.333333, 0.5, 0.5, 0.0}

In[3]:= EccentricityCentrality[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[3]= {0.333333, 0.5, 0.0}

In[4]:= EccentricityCentrality[Graph[{1,2},{}]]
Out[4]= {0.0, 0.0}
```

## BetweennessCentrality

- `BetweennessCentrality[g]`: the betweenness centrality of each vertex of `g`.

**Features**:
- Unnormalized; summed over unordered pairs on undirected graphs, ordered
  pairs on directed ones; weights ignored. Machine reals.
- **Deviation:** mixed graphs are left unevaluated (Wolfram's values there
  match no standard definition, e.g. `{1<->2, 2->3}` gives vertex 2 a
  betweenness of 0).
- **Algorithm:** Brandes, O(nm), with the sources spread over threads with
  per-thread accumulators.

```mathematica
In[1]:= BetweennessCentrality[StarGraph[5]]
Out[1]= {6.0, 0.0, 0.0, 0.0, 0.0}

In[2]:= BetweennessCentrality[PathGraph[4]]
Out[2]= {0.0, 2.0, 2.0, 0.0}

In[3]:= BetweennessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {1.0, 2.0, 3.0, 0.0}

In[4]:= BetweennessCentrality[Graph[{1<->2, 2->3}]]
Out[4]= BetweennessCentrality[Graph[<3 vertices, 2 edges>]]
```

## EdgeBetweennessCentrality

- `EdgeBetweennessCentrality[g]`: the betweenness centrality of each edge of `g`, in `EdgeList` order.

**Features**:
- *(w)* weight-aware; machine reals.
- Summed over **ordered** pairs even on undirected graphs (the edge of a `K2`
  scores 2).
- Weighted: ties within a relative `1e-12` share paths.
- **Algorithm:** Brandes, O(nm) unweighted, O(nm + n² log n) weighted, sources
  spread over threads with per-thread accumulators.

```mathematica
In[1]:= EdgeBetweennessCentrality[PathGraph[{1,2,3,4}]]
Out[1]= {6.0, 8.0, 6.0}

In[2]:= EdgeBetweennessCentrality[PathGraph[2]]
Out[2]= {2.0}

In[3]:= EdgeBetweennessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {4.0, 5.0, 3.0, 3.0}

In[4]:= EdgeBetweennessCentrality[Graph[{1,2,3,4},{1<->2,2<->3,3<->4,1<->4}, EdgeWeight->{1,1,1,5}]]
Out[4]= {6.0, 8.0, 6.0, 0.0}
```

## PageRankCentrality

- `PageRankCentrality[g]`: the PageRank centrality of each vertex of `g`, with damping `0.85`.
- `PageRankCentrality[g, a]`: uses damping factor `a`.

**Features**:
- Solves `x = a P^T x + (1 - a)/n`; dangling vertices jump uniformly;
  normalized to `Total[x] = 1`. Default `a = 0.85`; requires `0 <= a <= 1`
  (otherwise unevaluated). Weights ignored.
- Computed by (Jacobi) iteration with a correct stopping rule, converged to
  `1e-14` (Wolfram stops near `1e-9`, so they agree to ~9 digits).

```mathematica
In[1]:= PageRankCentrality[CycleGraph[4]]
Out[1]= {0.25, 0.25, 0.25, 0.25}

In[2]:= PageRankCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.213762, 0.264622, 0.307853, 0.213762}

In[3]:= PageRankCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], 0.5]
Out[3]= {0.22449, 0.265306, 0.285714, 0.22449}

In[4]:= PageRankCentrality[StarGraph[5]]
Out[4]= {0.475676, 0.131081, 0.131081, 0.131081, 0.131081}

In[5]:= PageRankCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], 2]
Out[5]= PageRankCentrality[Graph[<4 vertices, 4 edges>], 2]
```

## EigenvectorCentrality

- `EigenvectorCentrality[g]`: the eigenvector centrality of each vertex of `g`, following in-edges.
- `EigenvectorCentrality[g, "In"]`: the same (the default).
- `EigenvectorCentrality[g, "Out"]`: follows out-edges.

**Features**:
- Computed per strongly connected component: each component `C` with
  `|C| > 1` gets its Perron vector scaled to total
  `(|C| - 1)/Σ(|C'| - 1)`; single-vertex components get 0, so a DAG gives all
  zeros (*reverse-engineered*: reproduces Wolfram exactly on every
  disconnected / non-strongly-connected case tried, e.g. `K4 ⊔ K3 ⊔ K2` totals
  3:2:1).
- Weights ignored.
- **Algorithm:** restarted Arnoldi (Krylov dimension 40, BLAS
  re-orthogonalization, LAPACK on the Hessenberg matrix) per block,
  residual-tested to `1e-13`.

```mathematica
In[1]:= EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {0.333333, 0.333333, 0.333333, 0.0}

In[2]:= EigenvectorCentrality[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[2]= {0.269594, 0.269594, 0.315449, 0.145362}

In[3]:= EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Out"]
Out[3]= {0.333333, 0.333333, 0.333333, 0.0}

In[4]:= EigenvectorCentrality[Graph[{1->2,2->3}]]
Out[4]= {0.0, 0.0, 0.0}

In[5]:= EigenvectorCentrality[Graph[Join[EdgeList[CompleteGraph[4]], {5<->6,6<->7,5<->7}, {8<->9}]]]
Out[5]= {0.125, 0.125, 0.125, 0.125, 0.111111, 0.111111, 0.111111, 0.0833333, 0.0833333}
```

## KatzCentrality

- `KatzCentrality[g, a]`: the Katz centrality of each vertex of `g` with attenuation factor `a`.
- `KatzCentrality[g, a, b]`: uses the constant term `b` (a number or a list of per-vertex values).

**Features**:
- Solves `x = a A^T x + b`, with `b` = 1 by default, a number, or a list.
  Weights ignored.
- Answered even beyond the convergence radius, like Wolfram (`K4` with
  `a = 1/2` gives `-2`); uses a dense LAPACK solve, `n <= 4000`. A system
  singular up to rounding (LU pivot ratio below `1e-12`, e.g. `K4` with
  `a = 1/3`) is left unevaluated. Mathematica 15 returns rounding noise there
  (about `1.8*10^16` per vertex for `K4`, `a = 1/3`).
- No edges, or an exact `a = 0`, return `b` exactly.

```mathematica
In[1]:= KatzCentrality[PathGraph[3], 0.1]
Out[1]= {1.12245, 1.22449, 1.12245}

In[2]:= KatzCentrality[CompleteGraph[4], 1/2]
Out[2]= {-2.0, -2.0, -2.0, -2.0}

In[3]:= KatzCentrality[PathGraph[3], 0.1, 2]
Out[3]= {2.2449, 2.44898, 2.2449}

In[4]:= KatzCentrality[PathGraph[3], 0.1, {1, 2, 3}]
Out[4]= {1.2449, 2.44898, 3.2449}

In[5]:= KatzCentrality[PathGraph[3], 0]
Out[5]= {1, 1, 1}

In[6]:= KatzCentrality[Graph[{1,2,3},{}], 0.3]
Out[6]= {1, 1, 1}

In[7]:= KatzCentrality[CompleteGraph[4], 1/3]
Out[7]= KatzCentrality[Graph[<4 vertices, 6 edges>], 1/3]
```

## HITSCentrality

- `HITSCentrality[g]`: the list `{authorities, hubs}` of HITS centralities of the vertices of `g`.

**Features**:
- Returns `{h, A.h}`; `h` is built like `EigenvectorCentrality` but for
  `A^T A`, whose blocks are the classes of vertices sharing an in-neighbour
  (*reverse-engineered*; reproduces Wolfram exactly, including degenerate
  spectra such as `PathGraph[5]` and the all-zero answer for mixed graphs whose
  classes are singletons).
- Weights ignored. Spectral blocks use the same restarted Arnoldi solver as
  `EigenvectorCentrality`.

```mathematica
In[1]:= HITSCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {{0.5, 0.0, 0.0, 0.5}, {0.0, 0.0, 1.0, 0.0}}

In[2]:= HITSCentrality[Graph[{1->2,1->3,2->3}]]
Out[2]= {{0.0, 0.381966, 0.618034}, {1.0, 0.618034, 0.0}}

In[3]:= HITSCentrality[PathGraph[5]]
Out[3]= {{0.166667, 0.166667, 0.333333, 0.166667, 0.166667}, {0.166667, 0.5, 0.333333, 0.5, 0.166667}}

In[4]:= HITSCentrality[Graph[{1<->2, 3->4}]]
Out[4]= {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
```

## GraphTriangleCount

- `GraphTriangleCount[g]`: the number of triangles in `g`.

**Features**:
- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- Undirected: the number of 3-cliques. Directed (*reverse-engineered*):
  triangles are directed 3-cycles, so a transitive triple `1->2, 2->3, 1->3`
  does not count.
- **Algorithm:** triangle listing with an acyclic orientation (degree order,
  or index order when `maxdeg² <= 4m`), O(m^1.5) worst case; directed 3-cycles
  via two-bit arc flags per edge.

```mathematica
In[1]:= GraphTriangleCount[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 1

In[2]:= GraphTriangleCount[CompleteGraph[5]]
Out[2]= 10

In[3]:= GraphTriangleCount[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 1

In[4]:= GraphTriangleCount[Graph[{1->2,2->3,1->3}]]
Out[4]= 0

In[5]:= GraphTriangleCount[Graph[{1<->2, 2->3}]]
Out[5]= GraphTriangleCount[Graph[<3 vertices, 2 edges>]]
```

## LocalClusteringCoefficient

- `LocalClusteringCoefficient[g]`: the local clustering coefficient of each vertex of `g`.
- `LocalClusteringCoefficient[g, v]`: the local clustering coefficient of vertex `v`.

**Features**:
- All clustering heads are exact, ignore weights, and leave mixed graphs
  unevaluated (as in Wolfram).
- Undirected: `t(v)/C(d(v), 2)`, with `t(v)` the triangles at `v` and `d(v)`
  its degree.
- Directed (*reverse-engineered*): triangles are directed 3-cycles; local
  coefficient `c(v)/(in(v) out(v) - r(v))`, with `c(v)` the directed 3-cycles
  through `v` and `r(v)` the reciprocally linked neighbours.
- Triangles are found by the oriented triangle-listing algorithm of
  `GraphTriangleCount`.

```mathematica
In[1]:= LocalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= {1, 1, 1/3, 0}

In[2]:= LocalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}], 3]
Out[2]= 1/3

In[3]:= LocalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {1, 1, 1/2, 0}

In[4]:= LocalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}, EdgeWeight->{1,2,3,4}]]
Out[4]= {1, 1, 1/3, 0}

In[5]:= LocalClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[5]= LocalClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

## GlobalClusteringCoefficient

- `GlobalClusteringCoefficient[g]`: the global clustering coefficient (transitivity) of `g`.

**Features**:
- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- Undirected: `3T/Σ C(d, 2)`, with `T` the number of triangles and the sum
  over the vertex degrees `d`.
- Directed (*reverse-engineered*): triangles are directed 3-cycles, with the
  same in/out/reciprocal counting as `LocalClusteringCoefficient`.

```mathematica
In[1]:= GlobalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 3/5

In[2]:= GlobalClusteringCoefficient[PathGraph[4]]
Out[2]= 0

In[3]:= GlobalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 3/4

In[4]:= GlobalClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= GlobalClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

## MeanClusteringCoefficient

- `MeanClusteringCoefficient[g]`: the mean of the local clustering coefficients of the vertices of `g`.

**Features**:
- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- The arithmetic mean of `LocalClusteringCoefficient[g]` (undirected or
  directed definition accordingly).

```mathematica
In[1]:= MeanClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 7/12

In[2]:= MeanClusteringCoefficient[CompleteGraph[4]]
Out[2]= 1

In[3]:= MeanClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 5/8

In[4]:= MeanClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= MeanClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

## CompleteGraph

- `CompleteGraph[n]`: the complete graph `K_n` on `1..n`, with all `n(n-1)/2` edges.
- `CompleteGraph[{n1, n2, ...}]`: the complete multipartite graph with parts of sizes `n1, n2, ...`.

**Features**:
- All graph families of this module are undirected on `1..N`; the edge list is
  the sorted list of pairs `{i, j}`, `i < j` — identical to Wolfram's
  `EdgeList` for each family.
- `CompleteGraph[{n}]` is `K_n`.
- `CompleteGraph` is re-registered by a wrapper that delegates its
  pre-existing form (`CompleteGraph[n]`) to the original builtin.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: a family with more than 10^8 vertices or 5×10^7 edges is
  left unevaluated rather than allocating gigabytes.

```mathematica
In[1]:= CompleteGraph[5]
Out[1]= Graph[<5 vertices, 10 edges>]

In[2]:= CompleteGraph[{2, 3}]
Out[2]= Graph[<5 vertices, 6 edges>]

In[3]:= EdgeList[CompleteGraph[{2, 3}]]
Out[3]= {1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 4, 2 <-> 5}

In[4]:= EdgeList[CompleteGraph[{1, 1, 2}]]
Out[4]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4}

In[5]:= CompleteGraph[{4}]
Out[5]= Graph[<4 vertices, 6 edges>]

In[6]:= CompleteGraph[x]
Out[6]= CompleteGraph[x]
```

## WheelGraph

- `WheelGraph[n]`: the wheel graph on `n` vertices: hub `1` joined to every vertex of the cycle `2..n`.

**Features**:
- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Accepts `n = 1` or `n >= 4`; `n = 2` and `3` would be multigraphs and are
  left unevaluated.
- Options (`DirectedEdges`, layout options) are not supported: a call with
  options is left unevaluated.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= WheelGraph[5]
Out[1]= Graph[<5 vertices, 8 edges>]

In[2]:= EdgeList[WheelGraph[5]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}

In[3]:= WheelGraph[1]
Out[3]= Graph[<1 vertex, 0 edges>]

In[4]:= WheelGraph[3]
Out[4]= WheelGraph[3]

In[5]:= WheelGraph[5, DirectedEdges -> True]
Out[5]= WheelGraph[5, DirectedEdges -> True]
```

## HypercubeGraph

- `HypercubeGraph[n]`: the `n`-dimensional hypercube graph, with `2^n` vertices.

**Features**:
- Undirected on `1..2^n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- `HypercubeGraph[0]` is a single vertex.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= HypercubeGraph[3]
Out[1]= Graph[<8 vertices, 12 edges>]

In[2]:= EdgeList[HypercubeGraph[2]]
Out[2]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4}

In[3]:= VertexDegree[HypercubeGraph[4]]
Out[3]= {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}

In[4]:= HypercubeGraph[0]
Out[4]= Graph[<1 vertex, 0 edges>]
```

## GridGraph

- `GridGraph[{n1, ..., nk}]`: the `k`-dimensional grid graph with `n1 × ... × nk` vertices.

**Features**:
- Undirected on `1..N`, vertices numbered with the first coordinate varying
  fastest; the edge list is the sorted list of pairs `{i, j}`, `i < j`,
  identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limits (shared by every graph family of this module): a family
  with more than 10^8 vertices or 5×10^7 edges is left unevaluated rather than
  allocating gigabytes (e.g. `GridGraph[Table[2, {26}]]`, ~8.7×10^8 edges).

```mathematica
In[1]:= GridGraph[{3, 2}]
Out[1]= Graph[<6 vertices, 7 edges>]

In[2]:= EdgeList[GridGraph[{3, 2}]]
Out[2]= {1 <-> 2, 1 <-> 4, 2 <-> 3, 2 <-> 5, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[3]:= GridGraph[{2, 2, 2}]
Out[3]= Graph[<8 vertices, 12 edges>]

In[4]:= GridGraph[{5}]
Out[4]= Graph[<5 vertices, 4 edges>]

In[5]:= GridGraph[Table[2, {26}]]
Out[5]= GridGraph[{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}]
```

## KaryTree / CompleteKaryTree

- `KaryTree[n]`: the binary tree with `n` vertices.
- `KaryTree[n, k]`: the `k`-ary tree with `n` vertices.
- `CompleteKaryTree[n]`: the complete binary tree with `n` levels.
- `CompleteKaryTree[n, k]`: the complete `k`-ary tree with `n` levels.

**Features**:
- Undirected on `1..N`, root `1`, vertices numbered level by level; the edge
  list is the sorted list of pairs `{i, j}`, `i < j`, identical to Wolfram's
  `EdgeList`.
- `KaryTree` counts vertices; `CompleteKaryTree` counts levels.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= EdgeList[KaryTree[7]]
Out[1]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 2 <-> 5, 3 <-> 6, 3 <-> 7}

In[2]:= EdgeList[KaryTree[5, 3]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 5}

In[3]:= CompleteKaryTree[3]
Out[3]= Graph[<7 vertices, 6 edges>]

In[4]:= CompleteKaryTree[3, 3]
Out[4]= Graph[<13 vertices, 12 edges>]

In[5]:= EdgeList[CompleteKaryTree[2, 3]]
Out[5]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[6]:= KaryTree[0]
Out[6]= KaryTree[0]
```

## CirculantGraph

- `CirculantGraph[n, j]`: the circulant graph on `n` vertices in which each vertex `i` is joined to `i ± j (mod n)`.
- `CirculantGraph[n, {j1, j2, ...}]`: joins each vertex to `i ± j1, i ± j2, ...`.

**Features**:
- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= CirculantGraph[5, 1]
Out[1]= Graph[<5 vertices, 5 edges>]

In[2]:= EdgeList[CirculantGraph[6, 2]]
Out[2]= {1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 5, 4 <-> 6}

In[3]:= EdgeList[CirculantGraph[6, {1, 3}]]
Out[3]= {1 <-> 2, 1 <-> 4, 1 <-> 6, 2 <-> 3, 2 <-> 5, 3 <-> 4, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[4]:= VertexDegree[CirculantGraph[8, {1, 2}]]
Out[4]= {4, 4, 4, 4, 4, 4, 4, 4}
```

## PetersenGraph

- `PetersenGraph[]`: the Petersen graph (10 vertices, 15 edges).
- `PetersenGraph[n, k]`: the generalized Petersen graph `GP(n, k)`.

**Features**:
- Undirected on `1..2n`: the inner star is `1..n`, the outer cycle `n+1..2n`;
  the edge list is the sorted list of pairs `{i, j}`, `i < j`, identical to
  Wolfram's `EdgeList`.
- `PetersenGraph[]` is `PetersenGraph[5, 2]`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= PetersenGraph[]
Out[1]= Graph[<10 vertices, 15 edges>]

In[2]:= VertexDegree[PetersenGraph[]]
Out[2]= {3, 3, 3, 3, 3, 3, 3, 3, 3, 3}

In[3]:= GraphDiameter[PetersenGraph[]]
Out[3]= 2

In[4]:= EdgeList[PetersenGraph[4, 1]]
Out[4]= {1 <-> 2, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 8, 6 <-> 7, 7 <-> 8}
```

## TuranGraph

- `TuranGraph[n, k]`: the Turán graph `T(n, k)`: the complete `k`-partite graph on `n` vertices with part sizes as equal as possible.

**Features**:
- Undirected on `1..n`, larger parts first; the edge list is the sorted list
  of pairs `{i, j}`, `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= TuranGraph[5, 2]
Out[1]= Graph[<5 vertices, 6 edges>]

In[2]:= EdgeList[TuranGraph[5, 2]]
Out[2]= {1 <-> 4, 1 <-> 5, 2 <-> 4, 2 <-> 5, 3 <-> 4, 3 <-> 5}

In[3]:= TuranGraph[6, 3]
Out[3]= Graph[<6 vertices, 12 edges>]

In[4]:= EdgeList[TuranGraph[4, 3]]
Out[4]= {1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4, 3 <-> 4}
```

## HararyGraph

- `HararyGraph[k, n]`: the Harary graph `H(k, n)`: the minimal `k`-connected graph on `n` vertices.

**Features**:
- Requires `k >= 2` and `n > k`; otherwise unevaluated.
- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

```mathematica
In[1]:= HararyGraph[2, 5]
Out[1]= Graph[<5 vertices, 5 edges>]

In[2]:= EdgeList[HararyGraph[3, 6]]
Out[2]= {1 <-> 2, 1 <-> 4, 1 <-> 6, 2 <-> 3, 2 <-> 5, 3 <-> 4, 3 <-> 6, 4 <-> 5, 5 <-> 6}

In[3]:= EdgeList[HararyGraph[3, 5]]
Out[3]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}

In[4]:= HararyGraph[1, 5]
Out[4]= HararyGraph[1, 5]

In[5]:= HararyGraph[4, 4]
Out[5]= HararyGraph[4, 4]
```

## FindMaximumFlow

- `FindMaximumFlow[g, s, t]`: the maximum flow value from source `s` to sink `t` in the graph `g`.
- `FindMaximumFlow[g, s, t, "prop"]`: the flow property `"prop"` — `"FlowValue"`, `"FlowMatrix"` or `"EdgeList"`.

**Features**:
- Part of the graph-algorithm family (flows, cuts, matchings, covers, cliques,
  independent sets, Hamiltonian cycles, isomorphism, planarity) implemented in
  `src/graph/galg_*.c`, declared in `src/graph/graph_algos.h` and registered by
  `graph_algos_init()`. Every head in the family reads graphs through the
  validated-graph memo, so on a graph built by `Graph[...]` it starts from
  pre-resolved integer endpoints. Every head is `Protected`, stays unevaluated on
  a non-graph (the `*Q` predicates give `False`), and returns a fresh value.
  Semantics and output forms were checked against Mathematica 15 by a randomized
  differential test (`benchmarks/95-graph-algorithms/diff_mathematica.py`, about
  12k cases, 0 failures); where the answer is unique the outputs are identical.
- `s` and `t` may be lists of sources / sinks; `s == t` gives `0`.
- `"FlowValue"` is the default. `"FlowMatrix"` is a dense `n x n` matrix of edge
  flows (packed); Mathematica returns a `SparseArray`, which Mathilda does not
  have. `"EdgeList"` gives the edges carrying flow, oriented along the flow, in
  flow-matrix row order.
- Options: `EdgeCapacity -> {c1, ...}` (in `EdgeList` order) and
  `VertexCapacity -> {c1, ...}` (in `VertexList` order). Sources and sinks are
  capped by their own vertex capacity too, as in Mathematica. The two lists are
  parsed together onto one exact scale, so an exact `EdgeCapacity` mixes with a
  Rational/Real `VertexCapacity` (Mathematica leaves a non-integer
  `VertexCapacity` unevaluated). Matches Mathematica on 160 random
  vertex-capacitated instances (single and multiple terminals, directed and
  undirected). Where every capacity on the source side is `Infinity`,
  Mathematica 15 answers `0`; Mathilda gives the true value (`Infinity`, or the
  finite bottleneck).
- **Capacities ignore `EdgeWeight`**, exactly as Mathematica does; without
  `EdgeCapacity` every edge has capacity 1. An undirected edge carries flow
  either way.
- Numbers (shared by the whole flow/cut family): integer capacities/weights give
  exact Integers; Rational or Real ones give a Real (as Mathematica); `Infinity`
  is an allowed capacity; a negative or symbolic one leaves the call
  unevaluated. Internally everything is int64: reals are scaled by a common power
  of two, so the max flow of machine-real capacities is computed exactly.
- Algorithm: Dinic (BFS levels truncated at the sink, iterative blocking flow
  with current-arc pointers) on a CSR residual network.

```mathematica
In[1]:= FindMaximumFlow[CompleteGraph[4], 1, 4]
Out[1]= 3

In[2]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {2, 3, 1}]
Out[2]= 3

In[3]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, "FlowMatrix", EdgeCapacity -> {2, 3, 1}]
Out[3]= {{0, 2, 1}, {0, 0, 2}, {0, 0, 0}}

In[4]:= FindMaximumFlow[CycleGraph[6], 1, 4, "EdgeList"]
Out[4]= {1 <-> 2, 1 <-> 6, 2 <-> 3, 3 <-> 4, 5 <-> 4, 6 <-> 5}

In[5]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {1/2, 3, 1}]
Out[5]= 1.5

In[6]:= FindMaximumFlow[CycleGraph[6], {1, 2}, {4, 5}]
Out[6]= 2
```

`EdgeWeight` is not a capacity, an infinite capacity is allowed, `s == t` gives
`0`, and a negative capacity or a non-graph leaves the call unevaluated:

```mathematica
In[7]:= FindMaximumFlow[Graph[{1,2,3},{1->2,2->3,1->3}, EdgeWeight->{5,5,5}], 1, 3]
Out[7]= 2

In[8]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {Infinity, 3, 1}]
Out[8]= 4

In[9]:= FindMaximumFlow[CycleGraph[4], 2, 2]
Out[9]= 0

In[10]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {-1, 3, 1}]
Out[10]= FindMaximumFlow[Graph[<3 vertices, 3 edges>], 1, 3, EdgeCapacity -> {-1, 3, 1}]

In[11]:= FindMaximumFlow[PathGraph[{1, 2, 3}], 1, 3, EdgeCapacity -> {10, 10}, VertexCapacity -> {2, 10, 10}]
Out[11]= 2

In[12]:= FindMaximumFlow[CompleteGraph[4], {1, 2}, {3, 4}, VertexCapacity -> {3, 1, 1, 5}]
Out[12]= 3

In[13]:= FindMaximumFlow[PathGraph[{1, 2, 3}], 1, 3, VertexCapacity -> {1, 1/2, 1}]
Out[13]= 0.5
```

## FindMinimumCut

- `FindMinimumCut[g]`: `{value, {part1, part2}}`, a global minimum cut of `g` weighted by `EdgeWeight` (else 1).

**Features**:
- `Protected`; unevaluated on a non-graph and for fewer than 2 vertices.
- For a graph with directed edges, the edges counted run from `part1` (the
  source side) to `part2`.
- Among equal cuts the shore without `VertexList[g][[1]]` is listed first for
  undirected graphs (Mathematica's tie choice is not reproducible; the value
  always agrees).
- Integer weights give an Integer; Rational or Real weights give a Real (as
  Mathematica); a negative or symbolic weight leaves the call unevaluated.
  Computation is exact in int64 (reals scaled by a common power of two).
- Algorithms: Nagamochi-Ibaraki for undirected global minimum cuts
  (maximum-adjacency orders that contract every edge whose attachment reaches
  the current bound, not one pair per phase as in Stoer-Wagner); `2(n-1)`
  bounded flows (Dinic) for directed global cuts.

```mathematica
In[1]:= FindMinimumCut[CycleGraph[5]]
Out[1]= {2, {{2, 3, 4, 5}, {1}}}

In[2]:= FindMinimumCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3]}, EdgeWeight->{5,2}]]
Out[2]= {2, {{3}, {1, 2}}}

In[3]:= FindMinimumCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3]}, EdgeWeight->{1/2,2}]]
Out[3]= {0.5, {{2, 3}, {1}}}

In[4]:= FindMinimumCut[Graph[{1->2,2->3,3->1,1->3}]]
Out[4]= {1, {{1, 3}, {2}}}

In[5]:= FindMinimumCut[Graph[{1},{}]]
Out[5]= FindMinimumCut[Graph[<1 vertex, 0 edges>]]
```

## FindEdgeCut

- `FindEdgeCut[g]`: the edges of a minimum edge cut of `g`.
- `FindEdgeCut[g, s, t]`: the edges of a minimum `s`-`t` edge cut.

**Features**:
- `Protected`; unevaluated on a non-graph.
- Cuts are weighted by `EdgeWeight`; the edges are returned in `EdgeList` order.
- The `s`-`t` cut is the one closest to `s`, as in Mathematica.
- Directed graphs use strong connectivity.
- Uses the same exact int64 flow machinery as `FindMaximumFlow` /
  `FindMinimumCut` (Dinic; Nagamochi-Ibaraki for undirected global cuts).

```mathematica
In[1]:= FindEdgeCut[CycleGraph[5]]
Out[1]= {1 <-> 2, 5 <-> 1}

In[2]:= FindEdgeCut[CycleGraph[6], 1, 4]
Out[2]= {1 <-> 2, 6 <-> 1}

In[3]:= FindEdgeCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[3,1]}, EdgeWeight->{2,3,4}]]
Out[3]= {1 <-> 2, 2 <-> 3}

In[4]:= FindEdgeCut[Graph[{1->2,2->3,3->1}]]
Out[4]= {1 -> 2}

In[5]:= FindEdgeCut[x]
Out[5]= FindEdgeCut[x]
```

## FindVertexCut

- `FindVertexCut[g]`: a minimum vertex separator of `g`.
- `FindVertexCut[g, s, t]`: a minimum vertex separator between `s` and `t`.

**Features**:
- `Protected`; unevaluated on a non-graph.
- Works on the underlying undirected graph; vertices are returned in
  `VertexList` order.
- The `s`-`t` separator is the one closest to `t`; `{}` for adjacent `s`, `t`.
- Mathematica's conventions: a complete undirected graph gives its first `n-1`
  vertices; a graph with a directed edge whose underlying graph is complete
  gives `{}`.
- Algorithm: Even's split-vertex network with Esfahanian-Hakimi pair selection,
  solved by Dinic max flow. `VertexConnectivity` uses the same machinery (via
  `galg_vertex_connectivity`), replacing its former exponential subset search
  (identical answers on 300 random graphs).

```mathematica
In[1]:= FindVertexCut[CycleGraph[6]]
Out[1]= {2, 4}

In[2]:= FindVertexCut[CycleGraph[6], 1, 4]
Out[2]= {3, 5}

In[3]:= FindVertexCut[CycleGraph[6], 1, 2]
Out[3]= {}

In[4]:= FindVertexCut[CompleteGraph[4]]
Out[4]= {1, 2, 3}

In[5]:= FindVertexCut[Graph[{1->2,UndirectedEdge[1,3],UndirectedEdge[2,3]}]]
Out[5]= {}

In[6]:= VertexConnectivity[CycleGraph[6]]
Out[6]= 2
```

## EdgeConnectivity

- `EdgeConnectivity[g]`: the edge connectivity of `g`.
- `EdgeConnectivity[g, s, t]`: the `s`-`t` edge connectivity of `g`.

**Features**:
- `Protected`; unevaluated on a non-graph.
- Weighted by `EdgeWeight`; strong connectivity for directed graphs.
- Integer weights give an Integer; Rational or Real ones give a Real; computed
  exactly in int64 by the flow/cut engine (Dinic; Nagamochi-Ibaraki for the
  undirected global value; `2(n-1)` bounded flows for directed graphs).

```mathematica
In[1]:= EdgeConnectivity[CompleteGraph[5]]
Out[1]= 4

In[2]:= EdgeConnectivity[CycleGraph[6], 1, 4]
Out[2]= 2

In[3]:= EdgeConnectivity[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[3,1]}, EdgeWeight->{2,3,4}]]
Out[3]= 5

In[4]:= EdgeConnectivity[Graph[{1->2,2->3,3->1}]]
Out[4]= 1
```

## FindIndependentEdgeSet

- `FindIndependentEdgeSet[g]`: a maximum matching of `g` (a largest set of pairwise disjoint edges).

**Features**:
- `Protected`; unevaluated on a non-graph.
- Edge direction is ignored; edges are returned in `EdgeList` order.
- Algorithm: Karp-Sipser greedy start, then Hopcroft-Karp when the graph is
  bipartite and Edmonds' blossom algorithm otherwise (union-find blossom bases;
  failed searches retire their Hungarian trees). Matching on a triangulated grid:
  about 2.8 ms against 190 ms in Mathematica.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.
  (The differential test found Mathematica returning non-maximum matchings on
  weighted graphs.)

```mathematica
In[1]:= FindIndependentEdgeSet[PetersenGraph[]]
Out[1]= {1 <-> 3, 2 <-> 4, 5 <-> 10, 6 <-> 7, 8 <-> 9}

In[2]:= FindIndependentEdgeSet[StarGraph[4]]
Out[2]= {1 <-> 2}

In[3]:= FindIndependentEdgeSet[Graph[{1->2,3->2,3->4}]]
Out[3]= {1 -> 2, 3 -> 4}

In[4]:= FindIndependentEdgeSet[x]
Out[4]= FindIndependentEdgeSet[x]
```

## FindEdgeCover

- `FindEdgeCover[g]`: a minimum edge cover of `g` (a smallest set of edges touching every vertex).

**Features**:
- `Protected`; unevaluated on a non-graph.
- Computed as a maximum matching (`FindIndependentEdgeSet`) plus one edge per
  exposed vertex.
- `{}` when `g` has an isolated vertex, as in Mathematica.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.
  (The differential test found Mathematica returning non-minimum edge covers on
  weighted graphs.)

```mathematica
In[1]:= FindEdgeCover[PathGraph[{1,2,3,4,5}]]
Out[1]= {1 <-> 2, 2 <-> 3, 4 <-> 5}

In[2]:= FindEdgeCover[StarGraph[4]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[3]:= FindEdgeCover[Graph[{1,2,3},{UndirectedEdge[1,2]}]]
Out[3]= {}
```

## FindVertexCover

- `FindVertexCover[g]`: a minimum vertex cover of `g`, in `VertexList` order.

**Features**:
- `Protected`; unevaluated on a non-graph.
- Computed as the complement of a maximum independent set, so it shares the
  exact independent-set engine and its budgets (see `FindIndependentVertexSet`).
- **Exactness policy** (all NP-hard heads: `FindVertexCover`,
  `FindIndependentVertexSet`, `FindClique`, `FindKClique`, the Hamiltonian heads,
  and the isomorphism search): the result is a *proven* optimum / a complete
  answer, or the head stays unevaluated when its deterministic node budget is
  exhausted. All of them poll the `TimeConstrained` deadline, so e.g.
  `TimeConstrained[FindClique[g], 1]` returns `$Aborted` rather than hanging.
  They never return a merely-good answer.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.

```mathematica
In[1]:= FindVertexCover[PetersenGraph[]]
Out[1]= {1, 4, 5, 7, 8, 10}

In[2]:= FindVertexCover[StarGraph[5]]
Out[2]= {1}

In[3]:= FindVertexCover[CycleGraph[5]]
Out[3]= {1, 2, 4}

In[4]:= FindVertexCover[x]
Out[4]= FindVertexCover[x]
```

## FindIndependentVertexSet

- `FindIndependentVertexSet[g]`: `{s}` with `s` a maximum independent vertex set of `g`.
- `FindIndependentVertexSet[g, k]`: the largest maximal independent set with at most `k` vertices.
- `FindIndependentVertexSet[g, {k}]`: maximal independent sets with exactly `k` vertices.
- `FindIndependentVertexSet[g, {kmin, kmax}]`: maximal independent sets with between `kmin` and `kmax` vertices.
- `FindIndependentVertexSet[g, spec, n]` / `FindIndependentVertexSet[g, spec, All]`: up to `n` / all such sets.

**Features**:
- `Protected`; unevaluated on a non-graph.
- The size-spec forms enumerate MAXIMAL independent sets by size, largest
  first, exactly like the `FindClique` spec forms; they are limited to 8192
  vertices.
- Exact (see the exactness policy under `FindVertexCover`): connected
  components are solved separately.
  - A bipartite component of at least 256 vertices (grids, meshes, trees, even
    cycles) is solved in O(m sqrt n) by König's theorem (Hopcroft-Karp matching,
    then alternating reachability); `GridGraph[{1000, 1000}]` takes about 1.2 s.
  - Otherwise a component with average degree >= 8 (or density >= 0.05, up to
    3000 vertices) goes to the bitset maximum-clique search on its complement.
  - Sparser ones go to branch and reduce: degree-0/1 and triangle reductions,
    degree-2 folding, domination, a greedy clique-cover upper bound, component
    splitting at every node, and branching on a maximum-degree vertex with its
    mirrors.
- Budgets: the search gives up -- the head stays unevaluated, never a merely
  maximal set -- after 20 million nodes, 6 x 10^8 units of aggregate work
  (live-set size summed over nodes), 64M ints of per-level scratch or depth
  20000, so a huge non-bipartite sparse input (e.g.
  `RandomGraph[{200000, 300000}]`) returns unevaluated in bounded memory rather
  than exhausting it. The `TimeConstrained` deadline is polled.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.

```mathematica
In[1]:= FindIndependentVertexSet[PetersenGraph[]]
Out[1]= {{2, 3, 6, 9}}

In[2]:= FindIndependentVertexSet[CycleGraph[6], {2}, All]
Out[2]= {{3, 6}, {2, 5}, {1, 4}}

In[3]:= FindIndependentVertexSet[CycleGraph[5], Infinity, All]
Out[3]= {{3, 5}, {2, 5}, {2, 4}, {1, 4}, {1, 3}}

In[4]:= FindIndependentVertexSet[PathGraph[{1,2,3,4}], {1,2}, 2]
Out[4]= {{1, 4}, {1, 3}}

In[5]:= FindIndependentVertexSet[CompleteGraph[4], 3]
Out[5]= {{1}}

In[6]:= FindIndependentVertexSet[x]
Out[6]= FindIndependentVertexSet[x]
```

## IndependentVertexSetQ / VertexCoverQ

- `IndependentVertexSetQ[g, vs]`: `True` if no two vertices of `vs` are adjacent in `g`.
- `VertexCoverQ[g, vs]`: `True` if every edge of `g` has an endpoint in `vs`.

**Features**:
- `Protected` membership predicates; give `False` (never stay unevaluated) when
  `g` is not a graph.
- An element that is not a vertex of `g` gives `False`.
- Repeated vertices are allowed in `vs`.

```mathematica
In[1]:= IndependentVertexSetQ[CycleGraph[6], {1, 3, 5}]
Out[1]= True

In[2]:= IndependentVertexSetQ[CycleGraph[6], {1, 2}]
Out[2]= False

In[3]:= IndependentVertexSetQ[CycleGraph[6], {1, 1, 3}]
Out[3]= True

In[4]:= IndependentVertexSetQ[CycleGraph[6], {1, 7}]
Out[4]= False

In[5]:= VertexCoverQ[CycleGraph[4], {1, 3}]
Out[5]= True

In[6]:= VertexCoverQ[x, {1}]
Out[6]= False
```

## IndependentEdgeSetQ / EdgeCoverQ

- `IndependentEdgeSetQ[g, es]`: `True` if the edges `es` of `g` are pairwise disjoint (a matching).
- `EdgeCoverQ[g, es]`: `True` if every vertex of `g` is an endpoint of some edge in `es`.

**Features**:
- `Protected` membership predicates; give `False` when `g` is not a graph.
- An element that is not an edge of `g` gives `False`. An `UndirectedEdge`
  matches either orientation; a `DirectedEdge` matches only as given.

```mathematica
In[1]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[4,3]}]
Out[1]= True

In[2]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[2]= False

In[3]:= IndependentEdgeSetQ[Graph[{1->2,3->4}], {DirectedEdge[2,1]}]
Out[3]= False

In[4]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[3,4]}]
Out[4]= True

In[5]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[5]= False

In[6]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,3]}]
Out[6]= False
```

## FindClique

- `FindClique[g]`: `{c}` with `c` a maximum clique of `g`.
- `FindClique[g, k]`: the largest maximal clique with at most `k` vertices (`k` may be `Infinity`).
- `FindClique[g, {k}]`: maximal cliques with exactly `k` vertices.
- `FindClique[g, {kmin, kmax}]`: maximal cliques with between `kmin` and `kmax` vertices.
- `FindClique[g, spec, n]` / `FindClique[g, spec, All]`: up to `n` / all such maximal cliques.

**Features**:
- `Protected`; unevaluated on a non-graph.
- A directed graph needs edges both ways between clique members.
- Maximal cliques are listed by size, largest first; within one size
  Mathematica's order is reproduced
  (`FindClique[CycleGraph[5], Infinity, All]` gives
  `{{4,5},{3,4},{2,3},{1,5},{1,2}}`).
- `FindClique[g, 2]` is `{}` when every maximal clique is larger.
- Exact (see the exactness policy under `FindVertexCover`): a proven maximum or
  unevaluated when the deterministic node budget runs out; the
  `TimeConstrained` deadline is polled, so a timed-out search gives `$Aborted`.
- Maximum clique: bitset branch and bound with greedy-colouring bounds (Tomita's
  MCQ/MCS family in San Segundo's BBMC form), vertices numbered by reverse
  degeneracy order; graphs above 3000 vertices are decomposed by degeneracy
  (each vertex's later neighbourhood is a small local problem, skipped when its
  core number cannot beat the incumbent).
- Enumeration: pivoted Bron-Kerbosch on bitsets with size pruning.

```mathematica
In[1]:= FindClique[CompleteGraph[5]]
Out[1]= {{1, 2, 3, 4, 5}}

In[2]:= FindClique[CycleGraph[5], Infinity, All]
Out[2]= {{4, 5}, {3, 4}, {2, 3}, {1, 5}, {1, 2}}

In[3]:= FindClique[WheelGraph[6], {3}, 2]
Out[3]= {{1, 2, 6}, {1, 2, 3}}

In[4]:= FindClique[Graph[{1,2,3,4},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[1,3],UndirectedEdge[3,4]}], {2,3}, All]
Out[4]= {{1, 2, 3}, {3, 4}}

In[5]:= FindClique[CompleteGraph[4], 2]
Out[5]= {}

In[6]:= FindClique[Graph[{1->2,2->1,2->3,3->2,1->3}]]
Out[6]= {{1, 2}}
```

A search that exceeds its `TimeConstrained` limit aborts instead of hanging:

```mathematica
In[7]:= TimeConstrained[FindClique[RandomGraph[{400, 40000}]], 0.001]
Out[7]= $Aborted
```

## FindKClique

- `FindKClique[g, k]`: `{c}`, a largest set of vertices of `g` pairwise within distance `k`.

**Features**:
- `Protected`; unevaluated on a non-graph.
- Computed as a maximum clique of the `k`-th power graph, using the same exact
  BBMC maximum-clique search as `FindClique` (proven optimum or unevaluated on
  budget exhaustion; polls `TimeConstrained`).

```mathematica
In[1]:= FindKClique[CycleGraph[8], 2]
Out[1]= {{1, 2, 3}}

In[2]:= FindKClique[PathGraph[{1,2,3,4,5}], 1]
Out[2]= {{1, 2}}

In[3]:= FindKClique[x, 2]
Out[3]= FindKClique[x, 2]
```

## FindHamiltonianCycle

- `FindHamiltonianCycle[g]`: `{c}` with `c` a Hamiltonian cycle of `g` as a list of edges, or `{}` if there is none.
- `FindHamiltonianCycle[g, n]`: up to `n` Hamiltonian cycles.
- `FindHamiltonianCycle[g, All]`: all Hamiltonian cycles, each once.

**Features**:
- `Protected`; unevaluated on a non-graph.
- The cycle starts at `VertexList[g][[1]]`, each edge written in traversal
  order.
- The one-vertex graph gives `{{}}`; `K2` has none; a directed 2-cycle is one.
- Exact (see the exactness policy under `FindVertexCover`): complete answer or
  unevaluated on budget exhaustion; polls `TimeConstrained`.
- Algorithm: search over edge decisions with constraint propagation: each vertex
  needs exactly two chosen edges (directed: one in, one out), chosen edges form
  path fragments whose end-to-end links forbid short cycles, and the remaining
  graph must stay biconnected (directed: strongly connected) at every node.
- Performance: a random 400-vertex cubic graph takes about 1 ms. Large meshes
  are not a fast case: the linear per-node biconnectivity check makes
  `GridGraph[{n, n}]` roughly quadratic in the vertex count (100 x 100 about
  1.5 s, 150 x 150 about 7 s), so wrap far larger instances in
  `TimeConstrained`.

```mathematica
In[1]:= FindHamiltonianCycle[CycleGraph[4]]
Out[1]= {{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}

In[2]:= FindHamiltonianCycle[CompleteGraph[4], All]
Out[2]= {{1 <-> 2, 2 <-> 4, 4 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}, {1 <-> 3, 3 <-> 2, 2 <-> 4, 4 <-> 1}}

In[3]:= FindHamiltonianCycle[CompleteGraph[5], 2]
Out[3]= {{1 <-> 2, 2 <-> 4, 4 <-> 5, 5 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 1}}

In[4]:= FindHamiltonianCycle[PetersenGraph[]]
Out[4]= {}

In[5]:= FindHamiltonianCycle[Graph[{1->2,2->1}]]
Out[5]= {{1 -> 2, 2 -> 1}}

In[6]:= FindHamiltonianCycle[Graph[{1},{}]]
Out[6]= {{}}
```

```mathematica
In[7]:= FindHamiltonianCycle[CompleteGraph[2]]
Out[7]= {}

In[8]:= TimeConstrained[FindHamiltonianCycle[GridGraph[{150, 150}]], 0.05]
Out[8]= $Aborted
```

## FindHamiltonianPath

- `FindHamiltonianPath[g]`: a Hamiltonian path of `g` as a vertex list, or `{}` if there is none.
- `FindHamiltonianPath[g, s, t]`: a Hamiltonian path from `s` to `t`, or `{}`.

**Features**:
- `Protected`; unevaluated on a non-graph.
- The one-vertex graph gives `{}`, as in Mathematica.
- Paths reduce to cycles through an added vertex and use the
  `FindHamiltonianCycle` engine (exact; polls `TimeConstrained`).

```mathematica
In[1]:= FindHamiltonianPath[PetersenGraph[]]
Out[1]= {6, 7, 2, 4, 1, 3, 5, 10, 9, 8}

In[2]:= FindHamiltonianPath[CycleGraph[5], 1, 2]
Out[2]= {1, 5, 4, 3, 2}

In[3]:= FindHamiltonianPath[CycleGraph[5], 1, 3]
Out[3]= {}

In[4]:= FindHamiltonianPath[Graph[{1->2,2->3}]]
Out[4]= {1, 2, 3}

In[5]:= FindHamiltonianPath[StarGraph[4]]
Out[5]= {}

In[6]:= FindHamiltonianPath[Graph[{1},{}]]
Out[6]= {}
```

## HamiltonianGraphQ

- `HamiltonianGraphQ[g]`: `True` if `g` has a Hamiltonian cycle, `False` otherwise.

**Features**:
- `Protected`; gives `False` for a non-graph.
- `True` for the one-vertex graph, `False` for the graph with no vertices.
- Uses the `FindHamiltonianCycle` constraint-propagation search (exact;
  polls `TimeConstrained`).

```mathematica
In[1]:= HamiltonianGraphQ[CompleteGraph[5]]
Out[1]= True

In[2]:= HamiltonianGraphQ[PetersenGraph[]]
Out[2]= False

In[3]:= HamiltonianGraphQ[Graph[{1},{}]]
Out[3]= True

In[4]:= HamiltonianGraphQ[Graph[{},{}]]
Out[4]= False

In[5]:= HamiltonianGraphQ[x]
Out[5]= False
```

## IsomorphicGraphQ

- `IsomorphicGraphQ[g1, g2, ...]`: `True` if all the graphs `gi` are isomorphic.

**Features**:
- `Protected`. `False` if any argument is not a graph; one argument stays
  unevaluated, as in Mathematica.
- Directed and mixed graphs are supported (Mathematica leaves mixed graphs
  unevaluated). Edge weights are ignored, as in Mathematica.
- The reduction to the engine also folds self-loops into vertex colours and
  turns an edge of multiplicity k > 1 into a coloured subdivision vertex, so
  loops and multigraphs will work unchanged once `Graph` accepts them (today's
  validator rejects both). This applies to the whole isomorphism family
  (`FindGraphIsomorphism`, `CanonicalGraph`, `GraphAutomorphismGroup`).
- Engine (`galg_iso.c`): individualization-refinement in the nauty / bliss /
  Traces family. Hopcroft-style equitable refinement (U, out and in relations
  counted separately, largest fragment skipped, every split decision a function
  of cell positions, sizes and counts only) emits a 64-bit trace that is
  compared event by event, so a branch that cannot match dies at its first
  deviating split. The search tree is iterative with exact undo (no recursion at
  depth n). `g -> h` searches run in lockstep against `h`'s trace.
- Exact (see the exactness policy under `FindVertexCover`): complete answer or
  unevaluated on budget exhaustion; polls `TimeConstrained`.
- Performance: on an exactly 3-regular random graph with 10^4 vertices (where
  colour refinement alone learns nothing) `IsomorphicGraphQ` takes about 4 ms
  against about 0.8 s in Mathematica 15; at 10^5 vertices about 0.2 s against
  about 60 s.

```mathematica
In[1]:= IsomorphicGraphQ[CycleGraph[5], Graph[{UndirectedEdge[1,3],UndirectedEdge[3,5],UndirectedEdge[5,2],UndirectedEdge[2,4],UndirectedEdge[4,1]}]]
Out[1]= True

In[2]:= IsomorphicGraphQ[CycleGraph[6], PathGraph[Range[6]]]
Out[2]= False

In[3]:= IsomorphicGraphQ[CycleGraph[4], CycleGraph[4], CycleGraph[4]]
Out[3]= True

In[4]:= IsomorphicGraphQ[Graph[{1->2,2->3}], Graph[{1->2,1->3}]]
Out[4]= False

In[5]:= IsomorphicGraphQ[CycleGraph[4], x]
Out[5]= False

In[6]:= IsomorphicGraphQ[CycleGraph[4]]
Out[6]= IsomorphicGraphQ[Graph[<4 vertices, 4 edges>]]
```

## FindGraphIsomorphism

- `FindGraphIsomorphism[g1, g2]`: `{assoc}` with `assoc` an Association `v -> image` giving an isomorphism from `g1` to `g2`, or `{}` if none exists.
- `FindGraphIsomorphism[g1, g2, n]`: up to `n` isomorphisms.
- `FindGraphIsomorphism[g1, g2, All]`: all isomorphisms.

**Features**:
- `Protected`; unevaluated on non-graph arguments.
- The Association runs over `VertexList[g1]` in order.
- The order of the list of isomorphisms is the engine's search order;
  Mathematica's differs.
- Two empty graphs give `{<||>}` (Mathematica gives `{}`, which contradicts its
  own `IsomorphicGraphQ` answer `True`).
- Uses the `IsomorphicGraphQ` individualization-refinement engine; every map
  returned is verified edge by edge. Directed and mixed graphs are supported;
  edge weights are ignored.

```mathematica
In[1]:= FindGraphIsomorphism[CycleGraph[4], Graph[{a,b,c,d},{UndirectedEdge[a,c],UndirectedEdge[c,b],UndirectedEdge[b,d],UndirectedEdge[d,a]}]]
Out[1]= {<|1 -> a, 2 -> c, 3 -> b, 4 -> d|>}

In[2]:= FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], 2]
Out[2]= {<|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>, <|1 -> 1, 2 -> 4, 3 -> 3, 4 -> 2|>}

In[3]:= Length[FindGraphIsomorphism[CycleGraph[4], CycleGraph[4], All]]
Out[3]= 8

In[4]:= FindGraphIsomorphism[CycleGraph[4], StarGraph[4]]
Out[4]= {}

In[5]:= FindGraphIsomorphism[Graph[{},{}], Graph[{},{}]]
Out[5]= {<||>}
```

## CanonicalGraph

- `CanonicalGraph[g]`: a canonical representative of the isomorphism class of `g`, a graph on vertices `1..n` with sorted edges.

**Features**:
- `Protected`; unevaluated on a non-graph.
- `CanonicalGraph[g] === CanonicalGraph[h]` iff `g` and `h` are isomorphic.
- The particular representative differs from Mathematica's (both are
  arbitrary).
- Properties such as `EdgeWeight` are dropped, as in Mathematica.
- Canonical labeling keeps the best leaf by (trace, relabeled graph) with
  automorphism pruning (leaf and internal-node automorphisms, jump-back,
  first-path orbits) on the `IsomorphicGraphQ` engine. Directed and mixed graphs
  are supported.

```mathematica
In[1]:= EdgeList[CanonicalGraph[Graph[{UndirectedEdge[3,1],UndirectedEdge[1,2]}]]]
Out[1]= {1 <-> 3, 2 <-> 3}

In[2]:= CanonicalGraph[PathGraph[{1,2,3}]] === CanonicalGraph[Graph[{UndirectedEdge[2,1],UndirectedEdge[3,1]}]]
Out[2]= True

In[3]:= CanonicalGraph[CycleGraph[4]] === CanonicalGraph[StarGraph[4]]
Out[3]= False

In[4]:= CanonicalGraph[Graph[{1,2,3},{UndirectedEdge[1,2]}, EdgeWeight->{7}]] === CanonicalGraph[Graph[{1,2,3},{UndirectedEdge[1,3]}]]
Out[4]= True

In[5]:= CanonicalGraph[x]
Out[5]= CanonicalGraph[x]
```

## GraphAutomorphismGroup

- `GraphAutomorphismGroup[g]`: the automorphism group of `g` as `PermutationGroup[{Cycles[...], ...}]` acting on vertex positions.

**Features**:
- `Protected`; unevaluated on a non-graph.
- The group is given by a generating set: the group is Mathematica's, but the
  generators may differ.
- Mathilda has no permutation-group functions, so the result is an inert
  expression.
- The automorphisms found during canonical labeling (leaf and internal-node
  automorphisms, with jump-back and first-path orbit pruning) generate the whole
  group. Directed and mixed graphs are supported; edge weights are ignored.

```mathematica
In[1]:= GraphAutomorphismGroup[CycleGraph[4]]
Out[1]= PermutationGroup[{Cycles[{{2, 4}}], Cycles[{{1, 3}}], Cycles[{{1, 2}, {3, 4}}]}]

In[2]:= GraphAutomorphismGroup[PathGraph[{1,2,3}]]
Out[2]= PermutationGroup[{Cycles[{{1, 3}}]}]

In[3]:= GraphAutomorphismGroup[Graph[{1->2,2->3}]]
Out[3]= PermutationGroup[{}]

In[4]:= GraphAutomorphismGroup[x]
Out[4]= GraphAutomorphismGroup[x]
```

## PlanarGraphQ

- `PlanarGraphQ[g]`: `True` if `g` can be drawn in the plane without edge crossings, `False` otherwise.

**Features**:
- `Protected`; gives `False` for a non-graph. Edge directions are ignored.
- Algorithm (`galg_planar.c`): the linear-time left-right planarity test of de
  Fraysseix and Rosenstiehl in Brandes' formulation (the same algorithm as
  networkx's `check_planarity`), without the embedding phase. Two iterative DFS
  passes per connected component (orientation with lowpoints and nesting
  depths, then conflict-pair testing), so path-like graphs with 10^6 vertices
  need no call stack. The Euler bound rejects `m > 3n - 6` up front.
- O(n + m) time and memory.

```mathematica
In[1]:= PlanarGraphQ[CompleteGraph[4]]
Out[1]= True

In[2]:= PlanarGraphQ[CompleteGraph[5]]
Out[2]= False

In[3]:= PlanarGraphQ[CompleteGraph[{3,3}]]
Out[3]= False

In[4]:= PlanarGraphQ[PetersenGraph[]]
Out[4]= False

In[5]:= PlanarGraphQ[Graph[{1->2,2->3,3->1}]]
Out[5]= True

In[6]:= PlanarGraphQ[x]
Out[6]= False
```

