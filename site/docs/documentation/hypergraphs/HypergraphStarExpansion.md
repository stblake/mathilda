# HypergraphStarExpansion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphStarExpansion[h] gives the incidence (star) expansion of h: the bipartite Graph on VertexList[h] and nodes Hyperedge[1], ..., Hyperedge[m], with v<->Hyperedge[j] whenever v lies in hyperedge j.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= InputForm[HypergraphStarExpansion[Hypergraph[{{1,2},{2,3}}]]]
Out[1]= Graph[{1, 2, 3, Hyperedge[1], Hyperedge[2]}, {1 <-> Hyperedge[1], 2 <-> Hyperedge[1], 2 <-> Hyperedge[2], 3 <-> Hyperedge[2]}]

In[2]:= EdgeCount[HypergraphStarExpansion[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[2]= 9

In[3]:= InputForm[HypergraphStarExpansion[Hypergraph[{Hyperedge[1], 2},{{Hyperedge[1],2}}]]]
Out[3]= HypergraphStarExpansion[Hypergraph[{Hyperedge[1], 2}, {{Hyperedge[1], 2}}]]
```

### Applications (3)

```mathematica
In[4]:= g = HypergraphStarExpansion[Hypergraph[{{1, 2}, {2, 3}}]]
Out[4]= Graph[<5 vertices, 4 edges>]
```

The original vertices, then one Hyperedge[j] node per hyperedge

```mathematica
In[5]:= VertexList[g]
Out[5]= {1, 2, 3, Hyperedge[1], Hyperedge[2]}
```

Each vertex joins the hyperedge nodes it lies in

```mathematica
In[6]:= EdgeList[g]
Out[6]= {1 <-> Hyperedge[1], 2 <-> Hyperedge[1], 2 <-> Hyperedge[2], 3 <-> Hyperedge[2]}
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_star_expansion` gives the incidence (star)
bipartite `Graph`. Its vertices are `VertexList[h]` followed by one node
`Hyperedge[j]` per hyperedge; for each vertex `v` in the distinct-vertex set of
hyperedge `j` it emits `v <-> Hyperedge[j]`. Before building, each `Hyperedge[j]`
node is looked up in the hypergraph's own vertex index; if some vertex of `h` *is*
such a node the construction would be ambiguous, so the call is left unevaluated
(`clash`). The argument may be a `Hypergraph` or a bare List of hyperedges
(`hyp_arg`).

**Data structures.** The distinct-vertex CSR `soff/sv` and the memo's vertex
`GraphVIdx` (for the clash check). Vertices and edges are gathered into `Expr**`
arrays; the result is a `Graph`. Edge count equals the total incidence
`Σ|distinct e_j|`.

**Complexity / limits.** `O(n + Σ|e|)` — one `Hyperedge[j]` node per hyperedge and
one edge per incidence. The clash guard is the only reason the head declines a
well-formed hypergraph.

- Vertices are `VertexList[h]` followed by nodes `Hyperedge[1], ...,
  Hyperedge[m]`, with an edge `v <-> Hyperedge[j]` for each `v ∈ e_j`.
- Returns a simple `Graph`, validated and memoized as usual.
- Unevaluated if some vertex of `h` is itself such a `Hyperedge[j]` node.
- Benchmark (experiment 96): 10⁵ hyperedges 38.5 ms warm / 38.3 ms cold
  (Mathematica 313 ms, xgi 466 ms).

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

The star (incidence) expansion is the bipartite `Graph` of the
vertex–hyperedge incidence: its vertices are `VertexList[h]` followed by the nodes
`Hyperedge[1], ..., Hyperedge[m]`, with an edge `v <-> Hyperedge[j]` for each
vertex `v` of hyperedge `j`. Its edge count is the total incidence `Σ|e|`.

The star expansion is also the layout skeleton `HypergraphPlot` draws from. The one
case it declines: if a vertex of `h` happens to *be* a `Hyperedge[j]` node the
construction would be ambiguous, so the call is left unevaluated. It accepts a bare
List of hyperedges as well as a `Hypergraph`.
