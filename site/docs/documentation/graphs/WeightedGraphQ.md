# WeightedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeightedGraphQ[g] gives True if g carries edge weights (EdgeWeight).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Options (1)

```mathematica
In[1]:= {WeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], WeightedGraphQ[CycleGraph[3]], EdgeWeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], EdgeWeightedGraphQ[CycleGraph[3]], WeightedGraphQ[x]}
Out[1]= {True, False, True, False, False}
```

### Applications (4)

Carries an EdgeWeight list

```mathematica
In[2]:= WeightedGraphQ[Graph[{1 -> 2, 2 -> 3}, EdgeWeight -> {5, 7}]]
Out[2]= True
```

An undirected weighted edge

```mathematica
In[3]:= WeightedGraphQ[Graph[{1 <-> 2}, EdgeWeight -> {2}]]
Out[3]= True
```

A generator graph carries no weights

```mathematica
In[4]:= WeightedGraphQ[CompleteGraph[3]]
Out[4]= False
```

Unweighted, so False

```mathematica
In[5]:= WeightedGraphQ[CycleGraph[4]]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_weighted_graph_q` returns `True` iff `g` is a valid graph that
carries an `EdgeWeight` option list — the single test
`graph_is_valid(g) && graph_edge_weight_list(g) != NULL`. A graph built with
`EdgeWeight -> {w1, ...}` stores that list in the canonical 3-argument form; an unweighted
graph has no such list and gives `False`. `EdgeWeightedGraphQ` is defined to call this same
function, because Mathilda has no vertex weights and so the edge-weighted and weighted
questions coincide.

**Data structures.** `graph_edge_weight_list` returns a borrowed pointer to the
`EdgeWeight` `List` stored inside the canonical `Graph[List verts, List edges, opts...]`
tree (per-edge options are read by key, not position — see `graph.h`). `graph_is_valid` is
the memoized validator, so the validity half is `O(1)` after the first query. `gops_truth`
builds the result symbol.

**Complexity / limits.** `O(1)` apart from the first validation pass. The test is purely
presence-of-option: it does not inspect whether the stored weights are numeric or
non-negative (that stricter check is `graph_weights_usable`, used by the shortest-path
code). A non-graph argument gives `False` and never stays unevaluated.

- `Protected`. A non-graph argument gives `False`.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeightedGraphQ](../../graphs/EdgeWeightedGraphQ/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`WeightedGraphQ[g]` is `True` exactly when `g` is a valid graph built with an `EdgeWeight`
option. It coincides with `EdgeWeightedGraphQ` — Mathilda has no vertex weights, so a graph
is weighted iff its edges are. The test checks only for the *presence* of the weight list,
not that the weights are numeric or usable by a shortest-path algorithm.

An unweighted graph (every accessor then treats each edge's weight as `1`) gives `False`,
as does any non-graph argument.
