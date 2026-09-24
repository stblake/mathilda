# EdgeAdd

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeAdd[g, e] adds the edge e to g; EdgeAdd[g, {e1, ...}] adds several. Edges may be written u->v / DirectedEdge[u,v] or u<->v / UndirectedEdge[u,v]; endpoints not yet in g become new vertices (appended). A new edge gets weight 1 in a weighted graph. Adding a self-loop or an edge g already has (a multigraph) is left unevaluated.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

## Implementation notes

- `Protected`. A non-graph first argument is left unevaluated.
- Endpoints not in `g` become new vertices, appended in order.
- `u -> v` takes the graph's kind: undirected in an undirected (or edgeless)
  graph, directed otherwise; `DirectedEdge` is always directed.
- A new edge has weight 1 in a weighted graph.
- Deviation: Mathilda graphs are simple, so adding an edge that already exists
  (which would create parallel edges) or a self-loop is left unevaluated;
  Mathematica returns a multigraph.
- `O(V + E)` plus one hash per argument item; result seeded into the graph memo.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
