# DegreeCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DegreeCentrality[g] gives the list of vertex degrees of g; DegreeCentrality[g, "In"] and [g, "Out"] count incoming and outgoing edges. In a mixed graph an undirected edge counts as one edge in each direction.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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
Out[5]= DegreeCentrality[Graph[<4 vertices, 4 edges>], "Foo"]
```

### Applications (5)

The hub has degree 4, each leaf degree 1

```mathematica
In[6]:= DegreeCentrality[StarGraph[5]]
Out[6]= {4, 1, 1, 1, 1}
```

Count incoming arcs only

```mathematica
In[7]:= DegreeCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[1, 3], DirectedEdge[3, 2]}], "In"]
Out[7]= {0, 2, 1}
```

Count outgoing arcs only

```mathematica
In[8]:= DegreeCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[1, 3], DirectedEdge[3, 2]}], "Out"]
Out[8]= {2, 0, 1}
```

The default is in-degree plus out-degree

```mathematica
In[9]:= DegreeCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[1, 3], DirectedEdge[3, 2]}]]
Out[9]= {2, 2, 2}
```

Each vertex touches the other four

```mathematica
In[10]:= DegreeCentrality[CompleteGraph[5]]
Out[10]= {4, 4, 4, 4, 4}
```

## Implementation notes

**Algorithm.** `builtin_degree_centrality` accepts `[g]`, `[g, "In"]` or `[g, "Out"]`. One pass over the edge list increments `out[u]` and `in[v]` for each edge, and for an undirected edge also `out[v]` and `in[u]`, so an undirected edge counts as one arc each way. `"In"` and `"Out"` return the respective counts. The default returns their sum, except on a purely undirected graph, where `out` already equals the degree. In a mixed graph an undirected edge therefore adds 2 to the total at each endpoint, as in Mathematica. Weights are ignored and the result is exact integers.

**Data structures.** The edge-index view from `graph_edge_indices` supplies the parallel `eu`, `ev` and `directed` arrays. Two `int64_t` count arrays of length `n` are the only working storage, and the answer is a packed integer vector (`gmet_int_vector`).

**Complexity / limits.** `O(n + m)` time and `O(n)` space. Any mode string other than `"In"` or `"Out"`, or a non-string second argument, leaves the call unevaluated.

- Exact degrees, in `VertexList` order; weights ignored.
- On a mixed graph an undirected edge counts once in each direction (so it
  adds 2 to the default total), as in Wolfram.
- An unknown direction string leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gmet_centrality.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_centrality.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The result is a list of exact integers in `VertexList` order. On an undirected graph every mode equals the ordinary vertex degree. On a mixed graph an undirected edge counts as one arc in each direction, so it adds 2 to the default total at each endpoint. Edge weights are ignored.
