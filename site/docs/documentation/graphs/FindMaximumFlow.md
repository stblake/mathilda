# FindMaximumFlow

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMaximumFlow[g, s, t] gives the value of a maximum flow from s to t (s and t may be lists of sources and sinks). FindMaximumFlow[g, s, t, "prop"] gives "FlowValue", "FlowMatrix" (dense n x n matrix of edge flows) or "EdgeList" (edges carrying flow, oriented along it). Capacities come from the EdgeCapacity -> {c1, ...} option (EdgeList order; default 1, EdgeWeight is ignored); VertexCapacity -> {c1, ...} caps the flow through each vertex. Undirected edges carry flow either way. Dinic's algorithm; exact for integer capacities.`**

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindMaximumFlow[CompleteGraph[4], 1, 4]
Out[1]= 3

In[2]:= FindMaximumFlow[CycleGraph[6], 1, 4, "EdgeList"]
Out[2]= {1 <-> 2, 1 <-> 6, 2 <-> 3, 3 <-> 4, 5 <-> 4, 6 <-> 5}

In[3]:= FindMaximumFlow[CycleGraph[6], {1, 2}, {4, 5}]
Out[3]= 2
```

### Scope (1)

```mathematica
In[4]:= FindMaximumFlow[CycleGraph[4], 2, 2]
Out[4]= 0
```

### Options (9)

```mathematica
In[5]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {2, 3, 1}]
Out[5]= 3

In[6]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, "FlowMatrix", EdgeCapacity -> {2, 3, 1}]
Out[6]= {{0, 2, 1}, {0, 0, 2}, {0, 0, 0}}

In[7]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {1/2, 3, 1}]
Out[7]= 1.5

In[8]:= FindMaximumFlow[Graph[{1,2,3},{1->2,2->3,1->3}, EdgeWeight->{5,5,5}], 1, 3]
Out[8]= 2

In[9]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {Infinity, 3, 1}]
Out[9]= 4

In[10]:= FindMaximumFlow[Graph[{1->2,2->3,1->3}], 1, 3, EdgeCapacity -> {-1, 3, 1}]
Out[10]= FindMaximumFlow[Graph[<3 vertices, 3 edges>], 1, 3, EdgeCapacity -> {-1, 3, 1}]

In[11]:= FindMaximumFlow[PathGraph[{1, 2, 3}], 1, 3, EdgeCapacity -> {10, 10}, VertexCapacity -> {2, 10, 10}]
Out[11]= 2

In[12]:= FindMaximumFlow[CompleteGraph[4], {1, 2}, {3, 4}, VertexCapacity -> {3, 1, 1, 5}]
Out[12]= 3

In[13]:= FindMaximumFlow[PathGraph[{1, 2, 3}], 1, 3, VertexCapacity -> {1, 1/2, 1}]
Out[13]= 0.5
```

## Options & behaviour

`EdgeWeight` is not a capacity, an infinite capacity is allowed, `s == t` gives
`0`, and a negative capacity or a non-graph leaves the call unevaluated:

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/), [VertexList](../../graphs/VertexList/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
