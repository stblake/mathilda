# EdgeRules

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeRules[g] gives the edges of g as a list of rules u -> v (undirected edges too), in EdgeList order.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeRules[CycleGraph[3]]
Out[1]= {1 -> 2, 2 -> 3, 3 -> 1}

In[2]:= EdgeRules[Graph[{1->2, 2->3}]]
Out[2]= {1 -> 2, 2 -> 3}

In[3]:= EdgeRules[5]
Out[3]= EdgeRules[5]
```

### Applications (1)

Edges as u -> v rules

```mathematica
In[4]:= EdgeRules[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]
Out[4]= {1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}
```

## Implementation notes

**Algorithm.** `builtin_edge_rules` maps each edge of the graph to a `Rule` `u -> v`, in
`EdgeList` order. It reads the edge list directly and, for each edge, copies its two endpoints
and wraps them in `Rule` — regardless of whether the edge was a `DirectedEdge` or an
`UndirectedEdge`. Direction is therefore *not* encoded in the output: every edge becomes a plain
`->` rule, which is the form the `Graph` constructor accepts as input sugar.

**Data structures.** None beyond the result `List` of `Rule` nodes; no adjacency or vertex index
is built, and the edge order is the graph's stored order.

**Complexity / limits.** `O(E)`. Exactly one argument; a non-graph argument returns unevaluated.
The result is a plain list of rules, not a `Graph`.

- `Protected`. A non-graph argument is left unevaluated.
- Undirected and directed edges both become rules, in `EdgeList` order.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeRules[g]` gives the edges of `g` as a list of rules `u -> v`, in `EdgeList` order. Both
directed and undirected edges become plain `->` rules — the direction is not preserved in the
output, since the point of this form is to feed the edges back into `Graph` or `ReplaceAll`.

The result is an ordinary list of `Rule`s, not a graph; use `EdgeList` if you need the edges with
their `DirectedEdge` / `UndirectedEdge` heads intact.
