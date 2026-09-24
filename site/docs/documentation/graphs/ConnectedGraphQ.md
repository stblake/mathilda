# ConnectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedGraphQ[g] gives True if g is connected.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ConnectedGraphQ[CycleGraph[5]]
Out[1]= True

In[2]:= ConnectedGraphQ[Graph[{1,2,3},{1<->2}]]
Out[2]= False

In[3]:= ConnectedGraphQ[Graph[{1,2,3},{1->2,3->2}]]
Out[3]= True

In[4]:= ConnectedGraphQ[Graph[{},{}]]
Out[4]= False

In[5]:= ConnectedGraphQ[5]
Out[5]= ConnectedGraphQ[5]
```

## Algorithm

connectivity.c - ConnectedGraphQ[g] and VertexConnectivity[g].

Both operate on the underlying undirected graph.

```text
  ConnectedGraphQ[g]    True iff g has >= 1 vertex and forms a single
                        connected component.
  VertexConnectivity[g] the least number of vertices whose removal
                        disconnects g (n-1 for a complete graph, 0 if already
                        disconnected or trivial). Computed by brute-force
                        search over vertex subsets -- exact, intended for the
                        small graphs of a pico-CAS.
```

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

- `Protected`. Connectivity of the underlying undirected graph (weak
  connectivity for directed graphs), matching `ConnectedComponents`.
- The null graph is not connected.
- Unlike the `*Q` structural predicates (see `UndirectedGraphQ`), a non-graph
  argument leaves `ConnectedGraphQ` unevaluated; Mathematica gives `False`.

**Attributes:** `Protected`.

## References

**See also:** [ConnectedComponents](../../graphs/ConnectedComponents/), [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
