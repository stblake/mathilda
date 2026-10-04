# ConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConnectedComponents[g] gives the connected components of g: the strongly connected components when g has directed edges (listed with no edge from a component to a later one), else the components, largest first. ConnectedComponents[g, {v1, ...}] keeps only those containing some vi.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= ConnectedComponents[Graph[{1,2,3,4,5},{1<->2,3<->4,4<->5}]]
Out[1]= {{3, 4, 5}, {1, 2}}

In[2]:= ConnectedComponents[Graph[{1->2,2->3,3->1,3->4}]]
Out[2]= {{4}, {1, 2, 3}}

In[3]:= ConnectedComponents[Graph[{3->1, 1->5, 2->4, 2->6, 3->5, 4->6}]]
Out[3]= {{5}, {1}, {3}, {6}, {4}, {2}}

In[4]:= ConnectedComponents[Graph[{1->2,2->3,3->1,3->4}], {4}]
Out[4]= {{4}}

In[5]:= WeaklyConnectedComponents[Graph[{1,2,3,4},{1->2,3->2}]]
Out[5]= {{1, 2, 3}, {4}}

In[6]:= ConnectedComponents[Graph[{},{}]]
Out[6]= {}

In[7]:= ConnectedComponents[5]
Out[7]= ConnectedComponents[5]
```

### Applications (2)

Undirected: largest component first

```mathematica
In[8]:= ConnectedComponents[Graph[{1 <-> 2, 3 <-> 4, 4 <-> 5}]]
Out[8]= {{3, 4, 5}, {1, 2}}
```

Directed: strongly connected components

```mathematica
In[9]:= ConnectedComponents[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]
Out[9]= {{4}, {1, 2, 3}}
```

## Algorithm

components.c - connected-component builtins.

```text
  ConnectedComponents[g]          Mathematica's semantics, which depend on g:
                                    - g has a directed edge: the STRONGLY
                                      connected components, in an order with no
                                      edge from c_i to any later c_j (sinks
                                      first) -- Tarjan's completion order;
                                    - g is undirected: the components, largest
                                      first (ties by first appearance).
  ConnectedComponents[g, {v...}]  only the components containing some v.
  WeaklyConnectedComponents[g]    components of the underlying undirected
                                  graph, largest first; [g, {v...}] as above.
  StronglyConnectedComponents[g]  strong components (Tarjan) in the same
                                  sinks-first order as the directed
                                  ConnectedComponents (a Mathilda extension;
                                  Mathematica spells it ConnectedComponents).
```

Every form lists the vertices of a component in VertexList order.

Ordering, checked against Mathematica 15: Tarjan run from the vertices in VertexList order, following out-edges in EdgeList order, emits the strong components exactly in the order Mathematica returns them -- e.g. ConnectedComponents[{3->1, 1->5, 2->4, 2->6, 3->5, 4->6}] is {{5}, {1}, {3}, {6}, {4}, {2}} in both. Mathematica's order WITHIN an undirected component, and between equal-sized undirected components, follows no documented rule (it is not stable across such ties), so those use VertexList order and first appearance respectively.

Memory (SPEC section 4): results are freshly allocated; res is never touched, so the evaluator frees it on success and keeps it on NULL.

## Implementation notes

**Algorithm.** `builtin_connected_components` dispatches on direction. If the graph has any
directed edge it computes the **strongly connected components** with an iterative Tarjan
scan (`graph_strong_label`); otherwise it labels the **undirected components** by an iterative
DFS flood-fill over the combined out+in adjacency, then stably counting-sorts the components so
the largest comes first (ties broken by first appearance). Within each component the vertices
are kept in `VertexList` order.

**Data structures.** A CSR `GraphAdj` built from the validated graph. Tarjan uses `index[]`,
`low[]`, an `onstack[]` byte mask, an explicit vertex stack, a per-node child cursor and an
explicit recursion stack (no C recursion); the undirected path uses a `comp[]` label array with
an explicit DFS stack. A `keep[]` mask supports the selection form `ConnectedComponents[g, {v,
...}]`.

**Complexity / limits.** `O(V + E)` for either labelling, plus an `O(n+k)` counting sort. No
cap. The directed (Tarjan) components come out in completion order — a reverse topological order
of the condensation, so there is no edge from a component to a later one, matching Mathematica.
A non-graph argument returns unevaluated.

- `Protected`. Directed and mixed graphs: Tarjan's algorithm, the components
  listed so that no edge runs from a component to a later one (sinks first) —
  exactly Mathematica's order, e.g. `ConnectedComponents[Graph[{3->1, 1->5,
  2->4, 2->6, 3->5, 4->6}]]` is `{{5}, {1}, {3}, {6}, {4}, {2}}` in both. An
  undirected edge links its endpoints both ways.
- Undirected graphs, and `WeaklyConnectedComponents` always: the largest
  component first, as Mathematica documents; equal-sized components keep
  first-appearance order.
- Vertices within a component are in `VertexList` order. (Mathematica's order
  inside an undirected component, and between equal-sized ones, follows no
  documented rule and is not reproduced.)
- Linear time. The null graph has no components; unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [WeaklyConnectedComponents](../../graphs/WeaklyConnectedComponents/), [VertexList](../../graphs/VertexList/)

- R. E. Tarjan, *Depth-first search and linear graph algorithms*, SIAM J. Comput. **1** (1972) 146-160.
- Source: [`src/graph/components.c`](https://github.com/stblake/mathilda/blob/main/src/graph/components.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

For an undirected graph these are the connected components, listed largest first (ties by first
appearance). For a graph with directed edges they are the strongly connected components, listed
in reverse-topological order of the condensation — no edge runs from a component to a later one,
so sinks come first. Vertices inside a component keep `VertexList` order.

`ConnectedComponents[g, {v1, ...}]` keeps only the components containing one of the listed
vertices. `WeaklyConnectedComponents` ignores edge directions, and `StronglyConnectedComponents`
forces the directed interpretation even on an undirected graph.
