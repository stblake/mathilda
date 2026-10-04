# SimpleGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SimpleGraphQ[g] gives True if g is a graph with no self-loops or parallel edges -- every valid Mathilda graph -- and False otherwise.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {SimpleGraphQ[CycleGraph[3]], LoopFreeGraphQ[CycleGraph[3]], SimpleGraphQ[5], LoopFreeGraphQ[{1}]}
Out[1]= {True, True, False, False}
```

### Applications (4)

Every valid Mathilda graph is simple

```mathematica
In[2]:= SimpleGraphQ[CompleteGraph[4]]
Out[2]= True
```

A directed triangle is still simple

```mathematica
In[3]:= SimpleGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]
Out[3]= True
```

An undirected cycle

```mathematica
In[4]:= SimpleGraphQ[CycleGraph[6]]
Out[4]= True
```

A non-graph argument is False

```mathematica
In[5]:= SimpleGraphQ[5]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_simple_graph_q` returns `gops_truth(graph_is_valid(g))` — it is
`True` for every valid Mathilda graph and `False` otherwise. The reason is structural:
`Graph[...]` construction rejects self-loops and parallel/duplicate edges, so a canonical
graph is *by construction* simple. There is no separate scan for multi-edges or loops to
run; validity is the whole question.

**Data structures.** `graph_is_valid` is the memoized validator from `graph_util.c`: on
`g`'s first use it checks the canonical shape (`Graph[List verts, List edges]`, every edge
a 2-argument `DirectedEdge`/`UndirectedEdge`, no self-loops, no parallel edges, endpoints
all in `verts`) and caches the result on the node, so repeat calls are `O(1)`. `gops_truth`
maps the `int` to a fresh `True`/`False` symbol.

**Complexity / limits.** `O(1)` on a memo hit, one `O(V + E)` validation pass otherwise.
A non-graph argument is simply not valid, so the head gives `False`; it never stays
unevaluated.

- `Protected`. A non-graph argument gives `False`.
- Mathilda graphs are always simple, so both give `True` for every valid graph.

**Attributes:** `Protected`.

## References

**See also:** [LoopFreeGraphQ](../../graphs/LoopFreeGraphQ/)

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A simple graph has no self-loops and no parallel (duplicate) edges. Mathilda's `Graph[...]`
constructor rejects both at build time, so every canonical graph is simple — `SimpleGraphQ`
is therefore `True` for any valid graph and reduces to a validity check. It coincides with
`LoopFreeGraphQ` for the same reason.

The predicate is memoized through the validated-graph cache, so it is `O(1)` after the
first query. A non-graph argument gives `False` rather than staying unevaluated.
