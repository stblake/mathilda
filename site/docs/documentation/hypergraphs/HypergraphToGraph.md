# HypergraphToGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphToGraph[h] converts h, read as an ordered hypergraph, to the directed Graph with v_a->v_b for every a<b in each hyperedge {v_1, ..., v_k} (the Function Repository's HypergraphToGraph). Self-loops are dropped and parallel edges merged. h may be a plain List of hyperedges.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= InputForm[HypergraphToGraph[{{1,2,3},{3,4}}]]
Out[1]= Graph[{1, 2, 3, 4}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4}]

In[2]:= InputForm[HypergraphToGraph[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[2]= Graph[{1, 2, 3, 4, 5, 6, 7}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4, 4 -> 5, 4 -> 6, 5 -> 6}]

In[3]:= InputForm[HypergraphToGraph[{{1,1,2},{1,2},{5}}]]
Out[3]= Graph[{1, 2, 5}, {1 -> 2}]

In[4]:= HypergraphToGraph[5]
Out[4]= HypergraphToGraph[5]
```

### Applications (3)

```mathematica
In[5]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[5]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Directed edges v_a -> v_b for every a < b in a hyperedge

```mathematica
In[6]:= EdgeCount[HypergraphToGraph[h]]
Out[6]= 7
```

A repeated vertex would be a self-loop, which is dropped

```mathematica
In[7]:= InputForm[HypergraphToGraph[{{1, 1, 2}}]]
Out[7]= Graph[{1, 2}, {1 -> 2}]
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_to_graph` reads `h` as an **ordered**
hypergraph (the Wolfram Function Repository convention): hyperedge
`{v1, ..., vk}` contributes the directed edge `v_a -> v_b` for every position
`a < b`. It therefore reads the **raw** hyperedge CSR `ev` (order and repeats
preserved), not the distinct set. Because Mathilda graphs are simple, a self-loop
(from a repeated vertex, `x == y`) is dropped and parallel copies are merged via a
`U64Set` on `pair_key(x, y)`; every vertex of `h` is kept, including one that
occurs only in a unary hyperedge (the FR function drops those). Accepts a
`Hypergraph` or a bare List of hyperedges (`hyp_arg`).

**Data structures.** The raw hyperedge CSR `eoff/ev`; a `U64Set` for
directed-edge de-duplication; an `EVec` of `DirectedEdge` expressions. The result
is a `Graph` on `VertexList[h]`.

**Complexity / limits.** `O(Σ_j |e_j|²)` directed pairs, de-duplicated in
amortised `O(1)` each. The simple-graph normalisation (dropped self-loops, merged
parallels, retained vertices) is the documented deviation from the FR multigraph.

- Wolfram Function Repository name and semantics: hyperedge `{v1, ..., vk}`
  contributes `v_a -> v_b` for every `a < b`.
- Accepts a bare List of hyperedges, as the FR function does.
- *Deviation:* Mathilda graphs are simple, so self-loops (from a repeated
  vertex) are dropped and parallel copies merged, and every vertex is kept (the
  FR function returns a multigraph and drops vertices that occur only in unary
  hyperedges).
- Benchmark (experiment 96): 10⁵ hyperedges 59.3 ms (FR function 233 ms, xgi
  376 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphToGraph` reads `h` as an **ordered** hypergraph (the Function
Repository convention): a hyperedge `{v1, ..., vk}` contributes the directed edge
`v_a -> v_b` for every position `a < b`. This is the one head besides the arity
family that cares about hyperedge order and repeats.

Because Mathilda graphs are simple, the result deviates from the FR multigraph: a
self-loop from a repeated vertex is dropped, parallel copies are merged, and every
vertex is kept — including one occurring only in a unary hyperedge, which the FR
function drops. The head also accepts a bare List of hyperedges.
