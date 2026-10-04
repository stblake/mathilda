# HypergraphEdgeDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphEdgeDelete[h, e] or [h, {e1, ...}] deletes every hyperedge identical (SameQ) to one given.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {7, 8}]]
Out[1]= Hypergraph[{1, 2, 3, 4, 5, 6, 7, 8}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}, {7, 8}}]

In[2]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1,2}}], {1,2}]]
Out[2]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}}]

In[3]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {{3,4},{7}}]]
Out[3]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {4, 5, 6}}]

In[4]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1,2},{2,3},{1,2}}], {1,2}]]
Out[4]= Hypergraph[{1, 2, 3}, {{2, 3}}]

In[5]:= HypergraphEdgeDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {4, 3}]
Out[5]= HypergraphEdgeDelete[Hypergraph[<7 vertices, 4 hyperedges>], {4, 3}]
```

### Applications (4)

```mathematica
In[6]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[6]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Delete the hyperedge {3, 4}

```mathematica
In[7]:= InputForm[HypergraphEdgeDelete[h, {3, 4}]]
Out[7]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {4, 5, 6}, {7}}]
```

Delete several hyperedges

```mathematica
In[8]:= InputForm[HypergraphEdgeDelete[h, {{3, 4}, {7}}]]
Out[8]= Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {4, 5, 6}}]
```

Removes all copies of a repeated hyperedge

```mathematica
In[9]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1, 2}, {1, 2}}], {1, 2}]]
Out[9]= Hypergraph[{1, 2}, {}]
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_edge_delete` removes every hyperedge that is
`SameQ` to a named one — compared **as written**, so `{2, 1}` does not delete
`{1, 2}`, and all copies of a repeated hyperedge are removed. `edit_items`
(`edges = 1`) reads the argument as one hyperedge or a list of hyperedges. The
named hyperedges are put into a `GraphVIdx del` keyed on the whole `List` node;
each of `h`'s hyperedges is looked up there, its keep-bit set to miss, and a `hit`
bit recorded. If some named hyperedge never matched (its first slot's `hit` is
unset) the call is left unevaluated; otherwise `rebuild` emits the survivors.

**Data structures.** A `GraphVIdx del` over the named hyperedge Lists, a `hit`
mask, and vertex/hyperedge keep-masks (vertices are all kept); `rebuild` copies the
surviving hyperedges into a fresh `Hypergraph`.

**Complexity / limits.** `O(m + k)` hyperedge comparisons (each an `expr_hash` /
`SameQ` on the List). The vertex set is unchanged — deleting a hyperedge does not
remove vertices that became isolated.

- In the Edge heads, a List whose every element is a List is a list of
  hyperedges; otherwise it is one hyperedge.
- `HypergraphEdgeAdd` allows repeats (a multi-hypergraph).
- `HypergraphEdgeDelete` compares hyperedges as written (`{2, 1}` does not
  delete `{1, 2}`) and removes all copies of a repeated hyperedge. It is
  unevaluated if a named hyperedge does not occur.

**Attributes:** `Protected`.

## References

**See also:** [HypergraphEdgeAdd](../../hypergraphs/HypergraphEdgeAdd/), [SameQ](../../comparisons/SameQ/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphEdgeDelete[h, e]` removes every hyperedge identical (`SameQ`) to `e`;
the list form removes every hyperedge matching one of several. The vertices are
left in place, so deleting a hyperedge can leave its vertices isolated.

This is one of the order-sensitive heads: hyperedges are compared **as written**,
so `{2, 1}` does not delete `{1, 2}`. All copies of a repeated hyperedge are
removed together. A named hyperedge that does not occur in `h` leaves the call
unevaluated. The disambiguation rule matches `HypergraphEdgeAdd`: a List of Lists
is a list of hyperedges, anything else a single hyperedge.
