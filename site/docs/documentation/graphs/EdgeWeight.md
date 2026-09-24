# EdgeWeight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeWeight[g] gives the weights of g's edges, in EdgeList order. Defaults to all 1s if g was built without an EdgeWeight option.`**

## Examples (5)

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

### Options (2)

```mathematica
In[4]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]
Out[4]= {5, 7}

In[5]:= EdgeWeight[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{a,1/2}]]
Out[5]= {a, 1/2}
```

## Algorithm

edgeweight.c - EdgeWeight[g]: the graph's per-edge weights, in EdgeList order. Defaults to List[1, 1, ..., 1] (one per edge) when g carries no EdgeWeight -- matching Wolfram Language's own behavior for an unweighted graph, and giving WeightedAdjacencyMatrix[g] a well-defined answer for every valid graph, not just ones explicitly built with weights.

Memory (SPEC section 4): returns a fresh list; the evaluator frees res.

## Implementation notes

- `Protected`. Defaults to all `1`s when `g` was built without an `EdgeWeight`
  option. Weights may be symbolic or exact; they are returned as given.
- Unevaluated on a non-graph (see `VertexList`).
- Weight-aware consumers: `WeightedAdjacencyMatrix`, `FindShortestPath` and
  `GraphDistance` (see `FindShortestPath`). The other search/computation heads
  in this section ignore weights.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/), [Graph](../../graphs/Graph/), [VertexList](../../graphs/VertexList/), [WeightedAdjacencyMatrix](../../graphs/WeightedAdjacencyMatrix/), [FindShortestPath](../../graphs/FindShortestPath/), [GraphDistance](../../graphs/GraphDistance/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
