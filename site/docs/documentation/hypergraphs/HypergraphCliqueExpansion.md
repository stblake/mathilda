# HypergraphCliqueExpansion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphCliqueExpansion[h] gives the 2-section of h: the undirected Graph on VertexList[h] with u<->v whenever u and v share a hyperedge.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= InputForm[HypergraphCliqueExpansion[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Graph[{1, 2, 3, 4, 5, 6, 7}, {1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 4, 4 <-> 5, 4 <-> 6, 5 <-> 6}]

In[2]:= InputForm[HypergraphCliqueExpansion[Hypergraph[{{1,1,2},{2,3}}]]]
Out[2]= Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]

In[3]:= HypergraphCliqueExpansion[{{1, 2}}]
Out[3]= Graph[<2 vertices, 1 edge>]
```

### Applications (3)

```mathematica
In[4]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[4]= Hypergraph[<7 vertices, 4 hyperedges>]
```

An edge for every pair of vertices sharing a hyperedge

```mathematica
In[5]:= EdgeList[HypergraphCliqueExpansion[h]]
Out[5]= {1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 4, 4 <-> 5, 4 <-> 6, 5 <-> 6}
```

The 2-section keeps every vertex

```mathematica
In[6]:= VertexCount[HypergraphCliqueExpansion[h]]
Out[6]= 7
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_clique_expansion` gives the 2-section: the
undirected `Graph` on `VertexList[h]` with `u <-> v` whenever `u` and `v` share a
hyperedge. Each hyperedge is read as its distinct-vertex set `sv[soff[j]..]`; for
every pair `(a, b)` with `a < b` in that set, the oriented key `pair_key(min, max)`
is inserted into an open-addressing `U64Set`, and a first insertion emits an
`UndirectedEdge` in first-co-occurrence order. A vertex repeated inside a
hyperedge yields no self-loop (the distinct set collapses it). The argument may be
a `Hypergraph` or, for Function-Repository compatibility, a bare List of
hyperedges (`hyp_arg` wraps and evaluates it).

**Data structures.** The distinct-vertex CSR `soff/sv` from the memo; a growable
`U64Set` (64-bit packed-pair keys, multiplicative hash) for edge de-duplication;
a growable `EVec` of `UndirectedEdge` expressions. The result is a `Graph`, itself
validated and memoized by the Graph subsystem.

**Complexity / limits.** `O(Σ_v deg(v)²)` in the worst case — the number of
co-occurring pairs, which bounds the output edge count anyway. This is also the
graph whose shortest-path distance `HypergraphDistance` measures.

- Returns a simple `Graph` (validated and memoized as usual) on `VertexList[h]`;
  edges in order of first co-occurrence. A vertex repeated inside a hyperedge
  gives no self-loop.
- Also equals the graph whose distance `HypergraphDistance` measures.
- Benchmark (experiment 96): 10⁵ hyperedges 32.3 ms warm / 33.8 ms cold
  (Mathematica 278 ms, xgi 729 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [HypergraphDistance](../../hypergraphs/HypergraphDistance/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

The clique expansion (2-section) replaces each hyperedge by a clique: `u <-> v` in
the resulting `Graph` exactly when `u` and `v` lie together in some hyperedge. A
vertex repeated inside a hyperedge produces no self-loop, since the hyperedge is
read as its distinct-vertex set, and every vertex of `h` is kept even if isolated.

This is precisely the graph whose shortest-path distance `HypergraphDistance`
measures. Edges are listed in order of first co-occurrence. The head also accepts a
bare List of hyperedges, as the Function Repository function does.
