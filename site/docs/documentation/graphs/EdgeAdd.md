# EdgeAdd

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeAdd[g, e] adds the edge e to g; EdgeAdd[g, {e1, ...}] adds several. Edges may be written u->v / DirectedEdge[u,v] or u<->v / UndirectedEdge[u,v]; endpoints not yet in g become new vertices (appended). A new edge gets weight 1 in a weighted graph. Adding a self-loop or an edge g already has (a multigraph) is left unevaluated.`**

## Examples (10)

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

### Applications (6)

Closing a path into a triangle

```mathematica
In[5]:= EdgeList[EdgeAdd[PathGraph[{1, 2, 3}], 1 <-> 3]]
Out[5]= {1 <-> 2, 2 <-> 3, 1 <-> 3}
```

A list of edges extends the path

```mathematica
In[6]:= EdgeList[EdgeAdd[PathGraph[{1, 2, 3}], {3 <-> 4, 4 <-> 5}]]
Out[6]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}
```

An unseen endpoint becomes a new vertex

```mathematica
In[7]:= VertexList[EdgeAdd[PathGraph[{1, 2, 3}], 3 <-> 9]]
Out[7]= {1, 2, 3, 9}
```

In a directed graph the arrow keeps its direction

```mathematica
In[8]:= EdgeList[EdgeAdd[Graph[{1 -> 2, 2 -> 3}], 3 -> 1]]
Out[8]= {1 -> 2, 2 -> 3, 3 -> 1}
```

The arrow sugar is read as undirected here

```mathematica
In[9]:= EdgeList[EdgeAdd[Graph[{1 <-> 2}], 2 -> 3]]
Out[9]= {1 <-> 2, 2 <-> 3}
```

Edges can be added to an edgeless graph

```mathematica
In[10]:= EdgeList[EdgeAdd[Graph[{1, 2, 3}, {}], 1 <-> 2]]
Out[10]= {1 <-> 2}
```

## Implementation notes

**Algorithm.** `builtin_edge_add` accepts `EdgeAdd[g, e]` or `EdgeAdd[g, {e1, e2, ...}]`, where each item must parse as an edge (`gops_parse_edge`: `DirectedEdge`, `UndirectedEdge`, or the `->`/`<->` sugar); anything else leaves the call unevaluated. It copies the vertex, edge and endpoint arrays with room for the new edges, then walks the new edges in order. An endpoint already in the graph is found through the memoized vertex index; an unseen endpoint is appended to the vertex list once, in first-appearance order, via a small scratch hash. `u -> v` sugar takes the graph's kind: undirected when the graph has no directed edges, directed otherwise, while an explicit `DirectedEdge` stays directed. A weighted graph keeps its `EdgeWeight` list aligned, and each new edge gets weight 1.

**Data structures.** The graph is an ordinary `Graph[List[vertices], List[edges]]` expression tree with no dedicated `EXPR_*` tag. The edit works on the per-graph memo (`gops_view`: endpoint arrays `eu[k]`/`ev[k]`/`directed[k]`) and builds the result through `gops_graph_new`, sharing existing vertex and edge nodes and rebuilding only sugared edges.

**Complexity / limits.** `O(V + E)` for the copy plus one hash probe per new endpoint. Mathilda graphs are simple, so an edit that would create a self-loop or a parallel edge, such as adding an edge that already exists, is rejected when the result is validated and the call stays unevaluated instead of returning a multigraph.

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

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

Graphs here are simple, so adding an edge that already exists (or a self-loop) is not accepted and `EdgeAdd` returns unevaluated. Endpoints that are not yet vertices are appended to the vertex list in order of first appearance.

Each edge may be written `u <-> v`, `u -> v`, `UndirectedEdge[u, v]` or `DirectedEdge[u, v]`. In a graph with no directed edges, `u -> v` is read as undirected.
