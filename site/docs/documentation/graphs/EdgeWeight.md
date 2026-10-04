# EdgeWeight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeWeight[g] gives the weights of g's edges, in EdgeList order. Defaults to all 1s if g was built without an EdgeWeight option.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3}]]
Out[1]= {1, 1}

In[2]:= EdgeWeight[Graph[{1,2},{}]]
Out[2]= {}

In[3]:= EdgeWeight[5]
Out[3]= EdgeWeight[5]
```

### Options (3)

```mathematica
In[4]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]
Out[4]= {5, 7}

In[5]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{a,1/2}]]
Out[5]= {a, 1/2}

In[6]:= EdgeWeight[Graph[{1->2, 2->3}, EdgeWeight -> {5, 7}]]
Out[6]= {5, 7}
```

### Applications (2)

An unweighted graph defaults every weight to 1

```mathematica
In[7]:= EdgeWeight[CycleGraph[3]]
Out[7]= {1, 1, 1}
```

The stored weights, in EdgeList order

```mathematica
In[8]:= EdgeWeight[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 7}]]
Out[8]= {5, 7}
```

## Algorithm

edgeweight.c - EdgeWeight[g]: the graph's per-edge weights, in EdgeList order. Defaults to List[1, 1, ..., 1] (one per edge) when g carries no EdgeWeight -- matching Wolfram Language's own behavior for an unweighted graph, and giving WeightedAdjacencyMatrix[g] a well-defined answer for every valid graph, not just ones explicitly built with weights.

Memory (SPEC section 4): returns a fresh list; the evaluator frees res.

## Implementation notes

**Algorithm.** `builtin_edge_weight` returns `g`'s per-edge weights in `EdgeList`
order. It delegates to the shared `graph_resolve_edge_weights`
(`src/graph/graph_util.c`): when `g` carries an `EdgeWeight -> {w1, ...}` option
(the 3-argument canonical form) it returns a copy of that list, and otherwise
`{1, 1, ..., 1}` — one `1` per edge. Treating an unweighted edge as weight `1`
matches the Wolfram Language and means every valid graph has a well-defined
weight list, not just those built with explicit weights.

**Data structures.** The weight list lives inside the canonical `Graph` as a
borrowed `EdgeWeight` option list (read via `graph_edge_weight_list`, never by
argument position — which option is present decides where it sits, per
`src/graph/graph.h`). The returned list is a fresh copy so the caller owns it.

**Complexity / limits.** `O(E)`. Sharing the resolver with
`WeightedAdjacencyMatrix` guarantees the two builtins never disagree on what
"unweighted" defaults to. Returns `NULL` (unevaluated) when `g` is not a valid
graph.

- `Protected`. Defaults to all `1`s when `g` was built without an `EdgeWeight`
  option. Weights may be symbolic or exact; they are returned as given.
- Unevaluated on a non-graph (see `VertexList`).
- Weight-aware consumers: `WeightedAdjacencyMatrix`, `FindShortestPath` and
  `GraphDistance` (see `FindShortestPath`), `FindSpanningTree` (a minimum
  spanning tree), and the cut family. The option is accepted by both
  `Graph[e, EdgeWeight -> w]` and `Graph[v, e, EdgeWeight -> w]`.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/), [Graph](../../graphs/Graph/), [VertexList](../../graphs/VertexList/), [WeightedAdjacencyMatrix](../../graphs/WeightedAdjacencyMatrix/), [FindShortestPath](../../graphs/FindShortestPath/), [GraphDistance](../../graphs/GraphDistance/), [FindSpanningTree](../../graphs/FindSpanningTree/)

- Source: [`src/graph/edgeweight.c`](https://github.com/stblake/mathilda/blob/main/src/graph/edgeweight.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeWeight[g]` always returns a weight for every edge, in `EdgeList` order. A
graph built without weights reports `1` for each edge, matching the Wolfram
convention that an unweighted edge has weight `1` — so `WeightedAdjacencyMatrix`
and the weighted shortest-path routines have a defined weight to use on any
graph.

The weight list is matched to the edges by position: `EdgeWeight[g][[k]]` is the
weight of `EdgeList[g][[k]]`.
