# GraphDistanceMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDistanceMatrix[g] gives the matrix of shortest-path distances between all pairs of vertices of g, rows and columns in VertexList order; entry {i, j} is the distance from vertex i to vertex j (Infinity if unreachable). GraphDistanceMatrix[g, d] keeps only distances at most d. Integers for unweighted graphs, machine reals when g has EdgeWeight (edge weights as lengths).`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphDistanceMatrix[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {{0, 1, 2, 3}, {2, 0, 1, 2}, {1, 2, 0, 1}, {Infinity, Infinity, Infinity, 0}}

In[2]:= GraphDistanceMatrix[PathGraph[4]]
Out[2]= {{0, 1, 2, 3}, {1, 0, 1, 2}, {2, 1, 0, 1}, {3, 2, 1, 0}}

In[3]:= NDArrayQ[GraphDistanceMatrix[PathGraph[4]]]
Out[3]= True

In[4]:= GraphDistanceMatrix[PathGraph[5], 2]
Out[4]= {{0, 1, 2, Infinity, Infinity}, {1, 0, 1, 2, Infinity}, {2, 1, 0, 1, 2}, {Infinity, 2, 1, 0, 1}, {Infinity, Infinity, 2, 1, 0}}

In[5]:= GraphDistanceMatrix[x]
Out[5]= GraphDistanceMatrix[x]
```

### Options (1)

```mathematica
In[6]:= GraphDistanceMatrix[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[6]= {{0.0, 1.0, 2.0}, {1.0, 0.0, 1.0}, {2.0, 1.0, 0.0}}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [NDArrayQ](../../other-advanced/NDArrayQ/), [EdgeWeight](../../graphs/EdgeWeight/), [GraphDistance](../../graphs/GraphDistance/), [VertexEccentricity](../../graphs/VertexEccentricity/), [MeanGraphDistance](../../graphs/MeanGraphDistance/), [ClosenessCentrality](../../graphs/ClosenessCentrality/), [EccentricityCentrality](../../graphs/EccentricityCentrality/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
