# EdgeWeightedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeWeightedGraphQ[g] gives True if g carries edge weights (EdgeWeight).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Options (1)

```mathematica
In[1]:= {WeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], WeightedGraphQ[CycleGraph[3]], EdgeWeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], EdgeWeightedGraphQ[CycleGraph[3]], WeightedGraphQ[x]}
Out[1]= {True, False, True, False, False}
```

### Applications (4)

One weighted undirected edge

```mathematica
In[2]:= EdgeWeightedGraphQ[Graph[{1 <-> 2}, EdgeWeight -> {2}]]
Out[2]= True
```

A weighted directed path

```mathematica
In[3]:= EdgeWeightedGraphQ[Graph[{1 -> 2, 2 -> 3}, EdgeWeight -> {3, 4}]]
Out[3]= True
```

No EdgeWeight, so False

```mathematica
In[4]:= EdgeWeightedGraphQ[CycleGraph[4]]
Out[4]= False
```

A non-graph argument is False

```mathematica
In[5]:= EdgeWeightedGraphQ[0]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_edge_weighted_graph_q` is defined as a direct call to
`builtin_weighted_graph_q`, so it answers `True` iff `g` is a valid graph carrying an
`EdgeWeight` option list (`graph_is_valid(g) && graph_edge_weight_list(g) != NULL`). The
two predicates are deliberately identical: Mathilda has no vertex weights, so "edge
weighted" and "weighted" are the same property.

**Data structures.** Delegation means there is no separate state — `graph_edge_weight_list`
returns a borrowed pointer to the `EdgeWeight` `List` inside the canonical
`Graph[List verts, List edges, opts...]` tree, and `graph_is_valid` is the memoized
validator, so the validity half is `O(1)` after the first query on a node.

**Complexity / limits.** `O(1)` apart from the first validation pass. Like `WeightedGraphQ`
it tests only that an `EdgeWeight` list is present, not that the stored weights are numeric
or non-negative. A non-graph argument gives `False` and never stays unevaluated.

- `Protected`. A non-graph argument gives `False`.

**Attributes:** `Protected`.

## References

**See also:** [WeightedGraphQ](../../graphs/WeightedGraphQ/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeWeightedGraphQ[g]` is `True` exactly when `g` carries an `EdgeWeight` option — the
same test as `WeightedGraphQ`, since Mathilda associates weights only with edges. It checks
for the presence of the weight list, not that the weights are numeric.

A graph built without `EdgeWeight` is unweighted (each edge's weight defaults to `1` in the
accessors) and gives `False`, as does any non-graph argument.
