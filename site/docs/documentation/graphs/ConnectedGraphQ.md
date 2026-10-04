# ConnectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedGraphQ[g] gives True if g is connected (strongly connected when g has directed edges).`**

## Examples (9)

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

### Applications (3)

A path is connected

```mathematica
In[7]:= ConnectedGraphQ[PathGraph[{1, 2, 3}]]
Out[7]= True
```

Vertex 3 is isolated

```mathematica
In[8]:= ConnectedGraphQ[Graph[{1, 2, 3}, {1 <-> 2}]]
Out[8]= False
```

A cycle is connected

```mathematica
In[9]:= ConnectedGraphQ[CycleGraph[4]]
Out[9]= True
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

**Algorithm.** `builtin_connected_graph_q` follows Mathematica's rule: a graph
carrying any directed edge is "connected" only when **strongly** connected, and a
purely undirected graph when it has a single component. It builds the integer
adjacency and branches: with at least one directed edge it labels strongly
connected components (`graph_strong_label`, iterative Tarjan over the
out-adjacency, where an undirected edge links both ways) and requires exactly one
component; otherwise it counts undirected components (`graph_count_components`,
union/BFS). The empty graph is not connected (`n >= 1` is required).

**Data structures.** The `GraphAdj` CSR successor/predecessor lists
(`src/graph/graph.h`) plus a `comp[]` labelling array for the strong-component
pass. Everything runs on integer vertex indices from the validated-graph memo;
no expression is hashed during the scan.

**Complexity / limits.** `O(V + E)`. `NULL` (unevaluated) only on allocation
failure; a non-graph argument also returns `NULL`, since the builtin requires a
buildable adjacency.

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

- Source: [`src/graph/connectivity.c`](https://github.com/stblake/mathilda/blob/main/src/graph/connectivity.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

For an undirected graph this is ordinary connectivity: one component containing
every vertex. For a graph with directed edges it is *strong* connectivity —
every vertex must reach every other along the edge directions — matching
Mathematica, where a directed graph is "connected" only when strongly connected.

The empty graph (no vertices) is not connected.
