# ConnectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedGraphQ[g] gives True if g is connected (strongly connected when g has directed edges).`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= ConnectedGraphQ[CycleGraph[5]]
Out[1]= True

In[2]:= ConnectedGraphQ[Graph[{1,2,3},{1<->2}]]
Out[2]= False

In[3]:= ConnectedGraphQ[Graph[{1,2,3},{1->2,3->2}]]
Out[3]= False

In[4]:= ConnectedGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]
Out[4]= True

In[5]:= ConnectedGraphQ[Graph[{},{}]]
Out[5]= False

In[6]:= ConnectedGraphQ[5]
Out[6]= ConnectedGraphQ[5]
```

## Algorithm

connectivity.c - ConnectedGraphQ[g] and VertexConnectivity[g].

```text
  ConnectedGraphQ[g]    True iff g has >= 1 vertex and forms a single
                        connected component: strongly connected when g has
                        a directed edge (Mathematica's rule), connected
                        otherwise.
  VertexConnectivity[g] the least number of vertices whose removal
                        disconnects g (n-1 for a complete graph, 0 if already
                        disconnected or trivial). Computed by brute-force
                        search over vertex subsets -- exact, intended for the
                        small graphs of a pico-CAS.
```

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

- `Protected`. As in Mathematica, a graph with a directed edge must be
  **strongly** connected (every vertex reaches every other along the edges'
  directions; an undirected edge goes both ways); an undirected graph must be
  connected. So it agrees with `Length[ConnectedComponents[g]] == 1`.
- The null graph is not connected.
- Unlike the `*Q` structural predicates (see `UndirectedGraphQ`), a non-graph
  argument leaves `ConnectedGraphQ` unevaluated; Mathematica gives `False`.

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
