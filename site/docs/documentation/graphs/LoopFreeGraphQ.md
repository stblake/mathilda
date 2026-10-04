# LoopFreeGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LoopFreeGraphQ[g] gives True if g is a graph with no self-loops -- every valid Mathilda graph -- and False otherwise.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {SimpleGraphQ[CycleGraph[3]], LoopFreeGraphQ[CycleGraph[3]], SimpleGraphQ[5], LoopFreeGraphQ[{1}]}
Out[1]= {True, True, False, False}
```

### Applications (4)

No edge joins a vertex to itself

```mathematica
In[2]:= LoopFreeGraphQ[CycleGraph[3]]
Out[2]= True
```

Complete graphs have no self-loops

```mathematica
In[3]:= LoopFreeGraphQ[CompleteGraph[4]]
Out[3]= True
```

A directed path is loop-free

```mathematica
In[4]:= LoopFreeGraphQ[Graph[{1 -> 2, 2 -> 3}]]
Out[4]= True
```

A non-graph argument is False

```mathematica
In[5]:= LoopFreeGraphQ["x"]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_loop_free_graph_q` returns `gops_truth(graph_is_valid(g))`: it is
`True` for every valid Mathilda graph and `False` otherwise. A self-loop is an edge whose
two endpoints coincide, and `Graph[...]` construction rejects self-loops outright, so a
canonical graph has none — loop-freeness is guaranteed by validity, with no edge scan of
its own. (This makes `LoopFreeGraphQ` and `SimpleGraphQ` the same test in Mathilda.)

**Data structures.** The work is entirely in `graph_is_valid`, the memoized validator in
`graph_util.c`, which caches the canonical-shape / no-self-loop / no-parallel-edge result
on `g`'s node so repeat calls are `O(1)`. The graph is the ordinary
`Graph[List verts, List edges]` `Expr` tree; `gops_truth` builds the `True`/`False` symbol.

**Complexity / limits.** `O(1)` on a memo hit, one `O(V + E)` validation pass otherwise.
A non-graph argument is not valid, so the head gives `False` and never stays unevaluated.

- `Protected`. A non-graph argument gives `False`.
- Mathilda graphs are always simple, so both give `True` for every valid graph.

**Attributes:** `Protected`.

## References

**See also:** [SimpleGraphQ](../../graphs/SimpleGraphQ/)

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A loop-free graph has no self-loops — no edge from a vertex to itself. Mathilda's `Graph[...]`
constructor refuses self-loops, so every valid graph is loop-free and the predicate reduces
to a validity check. For the same structural reason it agrees with `SimpleGraphQ` on every
input.

The answer is memoized through the validated-graph cache (`O(1)` after the first query),
and a non-graph argument yields `False` rather than remaining unevaluated.
