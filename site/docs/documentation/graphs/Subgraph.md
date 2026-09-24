# Subgraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Subgraph[g, {v1, v2, ...}] gives the subgraph of g induced by the listed vertices (elements that are not vertices of g are ignored); Subgraph[g, patt] uses the vertices matching patt. Vertices come in the order given; edges in lower-triangular order of that vertex order. Edge weights are kept.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

## Implementation notes

- `Protected`. A non-graph first argument is left unevaluated.
- Vertices appear in the given order; non-vertices are ignored and repeats dropped.
- Edges are emitted at their later endpoint in that order, following each
  vertex's incidence order (out-edges and undirected edges, then in-edges).
- Weights are kept.
- An edge list (the edge-induced subgraph form) is not supported; the call is
  left unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
