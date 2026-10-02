# WeaklyConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeaklyConnectedComponents[g] gives the weakly connected components of g (edge directions ignored), largest first. WeaklyConnectedComponents[g, {v1, ...}] keeps only those containing some vi.`**

## Examples (7)

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

**See also:** [ConnectedComponents](../../graphs/ConnectedComponents/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
