# HypergraphEdgeAdd

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphEdgeAdd[h, e] appends the hyperedge e (a List); HypergraphEdgeAdd[h, {e1, ...}] appends several. New vertices are added.`**

## Examples (8)

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

### Applications (3)

One hyperedge; vertex 3 is new

```mathematica
In[6]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1, 2}}], {2, 3}]]
Out[6]= Hypergraph[{1, 2, 3}, {{1, 2}, {2, 3}}]
```

Several hyperedges

```mathematica
In[7]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1, 2}}], {{2, 3}, {3, 4}}]]
Out[7]= Hypergraph[{1, 2, 3, 4}, {{1, 2}, {2, 3}, {3, 4}}]
```

Repeats are allowed: a multi-hypergraph

```mathematica
In[8]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1, 2}}], {1, 2}]]
Out[8]= Hypergraph[{1, 2}, {{1, 2}, {1, 2}}]
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_edge_add` appends hyperedges to `h`. `edit_items`
with `edges = 1` reads the second argument as a *list of hyperedges* when it is a
List whose every element is itself a List, otherwise as a single hyperedge. Each
appended hyperedge's vertices that are new (not in the memo's vertex index, and
de-duplicated among the additions via a temporary `GraphVIdx extra`) are added to
the vertex set in first-appearance order; the hyperedges are then appended after
the existing ones. Repeats are allowed — the result is a multi-hypergraph.

**Data structures.** The memo's vertex `GraphVIdx` plus a temporary `extra` index;
the vertex and hyperedge Lists are rebuilt with `mk_hyp`/`mk_list` (existing
hyperedges shared by `expr_copy`, appended ones copied from the argument).

**Complexity / limits.** `O(n + Σ|existing| + Σ|added|)`. A plain non-List element
inside a would-be hyperedge leaves the call unevaluated.

- In the Edge heads, a List whose every element is a List is a list of
  hyperedges; otherwise it is one hyperedge.
- `HypergraphEdgeAdd` allows repeats (a multi-hypergraph).
- `HypergraphEdgeDelete` compares hyperedges as written (`{2, 1}` does not
  delete `{1, 2}`) and removes all copies of a repeated hyperedge. It is
  unevaluated if a named hyperedge does not occur.

**Attributes:** `Protected`.

## References

**See also:** [HypergraphEdgeDelete](../../hypergraphs/HypergraphEdgeDelete/), [SameQ](../../comparisons/SameQ/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphEdgeAdd[h, e]` appends the hyperedge `e`, adding any new vertices in
first-appearance order; `HypergraphEdgeAdd[h, {e1, ...}]` appends several. Repeats
are permitted, so the result can be a multi-hypergraph.

The disambiguation rule for the Edge heads: a `List` second argument is a *list of
hyperedges* when every one of its elements is itself a List, otherwise it is a
single hyperedge. So `{2, 3}` adds one hyperedge, while `{{2, 3}, {3, 4}}` adds
two. The edit builds a fresh object. To remove hyperedges, use
`HypergraphEdgeDelete`.
