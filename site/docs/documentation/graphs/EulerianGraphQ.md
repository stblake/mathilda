# EulerianGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EulerianGraphQ[g] gives True if g has a cycle using every edge exactly once: all degrees even (undirected) or in-degree = out-degree (directed), with all edges in one connected component. Left unevaluated for mixed graphs.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= {EulerianGraphQ[CycleGraph[4]], EulerianGraphQ[PathGraph[Range[3]]], EulerianGraphQ[Graph[{1->2,2->3,3->1}]], EulerianGraphQ[Graph[{1},{}]], EulerianGraphQ[Graph[{},{}]]}
Out[1]= {True, False, True, True, False}

In[2]:= {EulerianGraphQ[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]], EulerianGraphQ[Graph[{1,2,3,4},{1<->2,2<->3,3<->1}]], EulerianGraphQ[x]}
Out[2]= {False, True, False}

In[3]:= EulerianGraphQ[Graph[{1->2,2<->3}]]
Out[3]= EulerianGraphQ[Graph[<3 vertices, 2 edges>]]
```

## Implementation notes

- `Protected`. A non-graph argument gives `False`.
- Tests that all degrees are even (undirected) or in-degree equals out-degree
  (directed), and that all edges lie in one connected component (isolated
  vertices are allowed).
- An edgeless graph with at least one vertex is Eulerian; the null graph is not.
- Mixed graphs are left unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
