# HypergraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphQ[h] gives True if h is a valid Hypergraph, and False otherwise.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HypergraphQ[Hypergraph[{{1,2,3},{3,4}}]]
Out[1]= True

In[2]:= HypergraphQ[Graph[{1<->2}]]
Out[2]= False

In[3]:= GraphQ[Hypergraph[{{1,2,3},{3,4}}]]
Out[3]= False

In[4]:= HypergraphQ[Hypergraph[{1,2},{{1,5}}]]
Out[4]= False
```

### Applications (3)

A repeated vertex inside a hyperedge is allowed

```mathematica
In[5]:= HypergraphQ[Hypergraph[{{1, 1, 2}, {2, 3}}]]
Out[5]= True
```

A hypergraph is not a Graph

```mathematica
In[6]:= GraphQ[Hypergraph[{{1, 2}}]]
Out[6]= False
```

```mathematica
In[7]:= HypergraphQ["not a hypergraph"]
Out[7]= False
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_q` is a one-argument predicate that returns
`True`/`False` from `hypergraph_is_valid(h)`. That first checks the head is the
interned `Hypergraph` symbol, then calls `hyp_memo(h)`, which returns a memo slot
(validating and memoizing on a miss) or `NULL`. Validation requires the
`Hypergraph[List, List]` shape, pairwise-distinct vertices (a repeated vertex in
the vertex List is not canonical), and every hyperedge a `List` of declared
vertices. A `Graph`, a bare List of hyperedges, and a malformed — hence
unevaluated — `Hypergraph[...]` all give `False`; correspondingly `GraphQ` of a
hypergraph is `False`.

**Data structures.** No result tree beyond the `True`/`False` symbol; the work is
the shared validated-hypergraph memo (vertex `GraphVIdx` index plus the
vertex-index CSR of the hyperedges) described on the `Hypergraph` page.

**Complexity / limits.** `O(1)` on a memoized object (pointer-keyed slot hit);
`O(Σ|e|)` on the first validation of a new object, which also fills the memo. The
predicate has no side effects and never mutates its argument.

- `False` for a `Graph`, for a bare List of hyperedges, and for a malformed
  (unevaluated) `Hypergraph[...]`. Conversely `GraphQ` of a hypergraph is `False`.
- `O(1)` on a memoized hypergraph object.

**Attributes:** `Protected`.

## References

**See also:** [Hypergraph](../../hypergraphs/Hypergraph/), [Graph](../../graphs/Graph/), [GraphQ](../../graphs/GraphQ/)

- Source: [`src/graph/hyp_util.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_util.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphQ` is `True` only for a canonical, valid `Hypergraph` object: a
`Hypergraph[List, List]` with pairwise-distinct vertices and every hyperedge a
`List` of declared vertices. It is `False` for a `Graph`, for a bare List of
hyperedges, and for a malformed `Hypergraph[...]` that the constructor left
unevaluated (for instance one naming a vertex that is not in its vertex List).

The check is `O(1)` once the object has been seen, because validation is the same
pointer-keyed memo lookup every hypergraph head shares. Like the other `*Q`
predicates it always returns a Boolean and never leaves itself unevaluated.
