# DirectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DirectedGraphQ[g] gives True if g has at least one edge and all edges of g are directed. An edgeless graph counts as undirected.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= DirectedGraphQ[Graph[{1->2,2->3}]]
Out[1]= True

In[2]:= DirectedGraphQ[CycleGraph[3]]
Out[2]= False

In[3]:= DirectedGraphQ[Graph[{1,2},{}]]
Out[3]= False

In[4]:= DirectedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]
Out[4]= False

In[5]:= DirectedGraphQ[5]
Out[5]= False
```

### Applications (3)

Every edge is a DirectedEdge

```mathematica
In[6]:= DirectedGraphQ[Graph[{1, 2}, {1 -> 2}]]
Out[6]= True
```

An undirected graph is not directed

```mathematica
In[7]:= DirectedGraphQ[CycleGraph[3]]
Out[7]= False
```

An undirected edge, so False

```mathematica
In[8]:= DirectedGraphQ[Graph[{1, 2}, {1 <-> 2}]]
Out[8]= False
```

## Algorithm

directedq.c - DirectedGraphQ[g]: True iff g is a valid graph with at least one edge, all of them DirectedEdge. False otherwise -- including for an edgeless graph, which (as in the Wolfram Language) counts as undirected, so DirectedGraphQ and UndirectedGraphQ are never both True. Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res.

## Implementation notes

**Algorithm.** `builtin_directed_graph_q` is `True` exactly when the graph is
valid, has at least one edge, and *every* edge is a `DirectedEdge`. It reads the
directed-edge count from the validated-graph memo (`graph_directed_edge_count`,
`O(1)` after the first query on the node) and compares it with the total edge
count: equal and nonzero means fully directed. An edgeless graph counts as
undirected (as in the Wolfram Language), so `DirectedGraphQ` and
`UndirectedGraphQ` are never both `True`; a mixed graph fails both.

**Data structures.** None of its own — the whole test is two integers read off
the memo entry that `graph_is_valid` populates (`src/graph/graph_util.c`). The
result is a fresh `True`/`False` symbol; the evaluator frees the argument.

**Complexity / limits.** `O(1)` given a memoized graph (`O(V + E)` to populate
the memo on first contact). A non-graph argument gives `False`, never
unevaluated.

- `Protected`. `False` (never unevaluated) for a non-graph.
- An edgeless graph counts as undirected (as in the Wolfram Language), so
  `DirectedGraphQ` and `UndirectedGraphQ` are never both `True`. A mixed graph
  is neither.

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/directedq.c`](https://github.com/stblake/mathilda/blob/main/src/graph/directedq.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`DirectedGraphQ` is `True` only when the graph has at least one edge and all of
its edges are directed. An edgeless graph is treated as undirected, so it gives
`False` — and `DirectedGraphQ` and `UndirectedGraphQ` are never both `True` on
the same graph. A mixed graph (some directed, some undirected edges) gives
`False`; test for that with `MixedGraphQ`.
