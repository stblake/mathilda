# ConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedComponents[g] gives the connected components of g (weak, on the underlying undirected graph).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ConnectedComponents[Graph[{1,2,3,4,5},{1<->2,3<->4}]]
Out[1]= {{1, 2}, {3, 4}, {5}}

In[2]:= ConnectedComponents[Graph[{1,2,3},{1->2,3->2}]]
Out[2]= {{1, 2, 3}}

In[3]:= WeaklyConnectedComponents[Graph[{1,2,3,4},{1->2,3->2}]]
Out[3]= {{1, 2, 3}, {4}}

In[4]:= ConnectedComponents[Graph[{},{}]]
Out[4]= {}

In[5]:= ConnectedComponents[5]
Out[5]= ConnectedComponents[5]
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

- `Protected`. Edge direction is ignored, so both heads agree on every graph.
  Unweighted; see `FindShortestPath` for the shared search machinery.
- Mathematica's `ConnectedComponents` on a directed graph gives the *strongly*
  connected components; Mathilda's gives the weak ones (use
  `StronglyConnectedComponents` for the directed notion).
- The null graph has no components; unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [WeaklyConnectedComponents](../../graphs/WeaklyConnectedComponents/), [FindShortestPath](../../graphs/FindShortestPath/), [StronglyConnectedComponents](../../graphs/StronglyConnectedComponents/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
