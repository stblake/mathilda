# FindEdgeCover

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindEdgeCover[g] gives a minimum edge cover of g: a smallest set of edges touching every vertex. Gives {} when g has an isolated vertex (no edge cover exists).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindEdgeCover[PathGraph[{1,2,3,4,5}]]
Out[1]= {1 <-> 2, 2 <-> 3, 4 <-> 5}

In[2]:= FindEdgeCover[StarGraph[4]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[3]:= FindEdgeCover[Graph[{1,2,3},{UndirectedEdge[1,2]}]]
Out[3]= {}
```

### Applications (2)

Two edges cover all four vertices

```mathematica
In[4]:= FindEdgeCover[PathGraph[4]]
Out[4]= {1 <-> 2, 3 <-> 4}
```

An isolated vertex has no edge cover

```mathematica
In[5]:= FindEdgeCover[Graph[{1, 2, 3}, {2 <-> 3}]]
Out[5]= {}
```

## Implementation notes

**Algorithm.** `builtin_find_edge_cover` computes a minimum edge cover — a smallest set of edges
touching every vertex — by Gallai's reduction: take a maximum matching, then add one arbitrary
incident edge for each still-uncovered vertex, giving a cover of size `n - nu(g)`. The maximum
matching starts from a Karp–Sipser greedy phase (repeatedly match a degree-1 vertex to its
neighbour) and is completed exactly by **Hopcroft–Karp** when the graph is bipartite and by
**Edmonds' blossom algorithm** otherwise (union-find blossom bases, Hungarian-tree retirement).
Edge direction is ignored, matching Mathematica's semantics.

**Data structures.** A `GalgUG` CSR, a `mate[]` matching array, Karp–Sipser degree/queue arrays,
and for the blossom path a union-find `dsu[]` with tree labels, a `dead[]` retirement mask and
LCA marking. The cover is assembled in `EdgeList` order from the original edge expressions.

**Complexity / limits.** Polynomial — `O(E sqrt(V))` on the bipartite path, and each blossom
search `O(E·alpha(V))` over the part it reaches. No hard node cap, only a `TimeConstrained`
poll. If the graph has an isolated vertex no edge cover exists and the result is `{}`. A
non-graph argument returns unevaluated. The cover returned *is* a minimum.

- `Protected`; unevaluated on a non-graph.
- Computed as a maximum matching (`FindIndependentEdgeSet`) plus one edge per
  exposed vertex.
- `{}` when `g` has an isolated vertex, as in Mathematica.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.
  (The differential test found Mathematica returning non-minimum edge covers on
  weighted graphs.)

**Attributes:** `Protected`.

## References

**See also:** [FindIndependentEdgeSet](../../graphs/FindIndependentEdgeSet/)

- T. Gallai, *Über extreme Punkt- und Kantenmengen*, Ann. Univ. Sci. Budapest. Eötvös Sect. Math. **2** (1959) 133-138.
- J. Edmonds, *Paths, trees, and flowers*, Canad. J. Math. **17** (1965) 449-467.
- J. E. Hopcroft and R. M. Karp, *An n^{5/2} algorithm for maximum matchings in bipartite graphs*, SIAM J. Comput. **2** (1973) 225-231.
- Source: [`src/graph/galg_matching.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_matching.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

A minimum edge cover is a smallest set of edges that touches every vertex; its size is `n` minus
the size of a maximum matching. The edges are returned in `EdgeList` order and the set is
guaranteed minimum.

A graph with any isolated vertex has no edge cover at all (nothing can touch that vertex), so the
result is `{}`.
