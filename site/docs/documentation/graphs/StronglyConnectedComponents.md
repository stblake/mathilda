# StronglyConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StronglyConnectedComponents[g] gives the strongly connected components of g (following edge directions), in the same order as ConnectedComponents on a directed graph.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3}]]
Out[1]= {{3}, {2}, {1}}

In[2]:= StronglyConnectedComponents[Graph[{1,2,3,4},{1->2,2->1,2->3,3->4,4->3}]]
Out[2]= {{3, 4}, {1, 2}}

In[3]:= StronglyConnectedComponents[PathGraph[3]]
Out[3]= {{1, 2, 3}}
```

### Applications (1)

The 3-cycle, then its sink

```mathematica
In[4]:= StronglyConnectedComponents[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]
Out[4]= {{4}, {1, 2, 3}}
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

**Algorithm.** `builtin_strongly_connected_components` always computes strongly connected
components with an iterative Tarjan scan (`graph_strong_label`), never the undirected flood-fill
that `ConnectedComponents` falls back to. Two vertices share a component when each is reachable
from the other following edge directions; an undirected edge is traversed both ways. The
components are returned in Tarjan completion order — the same order the directed form of
`ConnectedComponents` uses.

**Data structures.** A CSR `GraphAdj` and Tarjan's `index[]`, `low[]`, `onstack[]` byte mask,
an explicit vertex stack, a per-node child cursor and an explicit recursion stack (no C
recursion). The component order array is the identity, since Tarjan's order is already the
answer order.

**Complexity / limits.** `O(V + E)`. Exactly one argument; a non-graph argument returns
unevaluated. (Mathematica spells this operation `ConnectedComponents`; Mathilda keeps the
explicit name as well.)

- `Protected`. Tarjan's algorithm, listed in the same sinks-first order as
  `ConnectedComponents` on a directed graph (Mathematica has no separate head:
  there `ConnectedComponents` is the strong notion). On an undirected graph the
  strong components are the ordinary ones.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [ConnectedComponents](../../graphs/ConnectedComponents/)

- R. E. Tarjan, *Depth-first search and linear graph algorithms*, SIAM J. Comput. **1** (1972) 146-160.
- Source: [`src/graph/components.c`](https://github.com/stblake/mathilda/blob/main/src/graph/components.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

Two vertices lie in the same strongly connected component when each can reach the other
following edge directions. The components come out in reverse-topological order of the
condensation (sinks first), the same order the directed `ConnectedComponents` uses.

An undirected edge is traversed in both directions, so on an undirected graph the strongly and
weakly connected components coincide.
