# StronglyConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StronglyConnectedComponents[g] gives the strongly connected components of g (following edge directions).`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3}]]
Out[1]= {{1}, {2}, {3}}

In[2]:= StronglyConnectedComponents[Graph[{1,2,3,4},{1->2,2->1,2->3,3->4,4->3}]]
Out[2]= {{1, 2}, {3, 4}}

In[3]:= StronglyConnectedComponents[PathGraph[3]]
Out[3]= {{1, 2, 3}}
```

## Algorithm

components.c - connected-component builtins.

```text
  ConnectedComponents[g]          weak components (underlying undirected)
  WeaklyConnectedComponents[g]    same as ConnectedComponents
  StronglyConnectedComponents[g]  strong components (Tarjan) over directed
                                  adjacency; for undirected graphs this
                                  coincides with the weak components.
```

Each returns a List of Lists of vertices, components in first-appearance order, vertices within a component in canonical index order.

Memory (SPEC section 4): returns freshly-allocated lists; frees res.

## Implementation notes

- `Protected`. Tarjan's algorithm. For undirected graphs this coincides with the
  weak components (see `ConnectedComponents`).
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [ConnectedComponents](../../graphs/ConnectedComponents/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
