# AdjacencyList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyList[g] gives the adjacency list of g; AdjacencyList[g,v] gives the vertices adjacent to v (successors for directed edges).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

## Algorithm

adjlist.c - AdjacencyList[g] and AdjacencyList[g, v].

Neighbors of v: successors for directed edges (v -> u yields u), and both endpoints for undirected edges. This matches the row convention of AdjacencyMatrix (a 1 in row v, column u means v is adjacent to u). Neighbors are returned in first-appearance order, de-duplicated.

```text
AdjacencyList[g]     -> {neighbors(v1), neighbors(v2), ...} in vertex order.
AdjacencyList[g, v]  -> neighbors(v).
```

Memory (SPEC section 4): returns freshly-allocated lists; evaluator frees res.

## Implementation notes

- `Protected`. Directed edges contribute successors (`v -> u` makes `u` a
  neighbor of `v`, not the reverse); undirected edges go both ways.
- Unevaluated on a non-graph (see `VertexList`) or when `v` is not a vertex.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
