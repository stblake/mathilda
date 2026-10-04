# FindIndependentEdgeSet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindIndependentEdgeSet[g] gives a maximum independent edge set (maximum matching) of g, ignoring edge direction. Hopcroft-Karp for bipartite graphs, Edmonds' blossom algorithm otherwise.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FindIndependentEdgeSet[PetersenGraph[]]
Out[1]= {1 <-> 3, 2 <-> 4, 5 <-> 10, 6 <-> 7, 8 <-> 9}

In[2]:= FindIndependentEdgeSet[StarGraph[4]]
Out[2]= {1 <-> 2}

In[3]:= FindIndependentEdgeSet[Graph[{1->2,3->2,3->4}]]
Out[3]= {1 -> 2, 3 -> 4}

In[4]:= FindIndependentEdgeSet[x]
Out[4]= FindIndependentEdgeSet[x]
```

### Applications (3)

The complete graph on 4 vertices has a perfect matching of two edges

```mathematica
In[5]:= FindIndependentEdgeSet[CompleteGraph[4]]
Out[5]= {1 <-> 2, 3 <-> 4}
```

Alternate edges of the path

```mathematica
In[6]:= FindIndependentEdgeSet[PathGraph[{1, 2, 3, 4}]]
Out[6]= {1 <-> 2, 3 <-> 4}
```

An odd cycle matches all but one vertex

```mathematica
In[7]:= FindIndependentEdgeSet[CycleGraph[5]]
Out[7]= {1 <-> 2, 4 <-> 5}
```

## Implementation notes

**Algorithm.** `FindIndependentEdgeSet[g]` returns a maximum **matching** (a
largest set of pairwise non-adjacent edges), ignoring edge direction as Wolfram
does. `galg_max_matching` starts from a Karp-Sipser greedy matching
(`gm_karp_sipser`): it repeatedly matches a degree-1 vertex to its only neighbour
— always safe, some maximum matching contains that edge — and otherwise an
arbitrary remaining edge, which on sparse graphs is already optimal or within a
few edges. It then completes to an exact maximum: a bipartiteness BFS
(`gm_bipartite`) routes bipartite graphs to Hopcroft-Karp (`gm_hopcroft_karp`,
phases of shortest vertex-disjoint augmenting paths) and everything else to
Edmonds' blossom algorithm (`gm_edmonds`), one alternating-tree search per
exposed vertex with odd cycles contracted into blossoms.

**Data structures.** The graph is a CSR `GalgUG`; the matching is a single
`mate[]` array (`mate[v]` = partner or `-1`). Blossom contraction keeps bases in
a union-find forest (`gm_find` with path halving), so one search costs
`O(E·α(V))` over the part of the graph it reaches; only the vertices a search
touched are reset afterwards, and a search that fails leaves a Hungarian tree
whose vertices are retired for the rest of the run (Edmonds' lemma). The matched
edges are emitted as `g`'s own edge expressions in `EdgeList` order by
`gm_matched_edges`.

**Complexity / limits.** Hopcroft-Karp is `O(E√V)`; the blossom search is
`O(V·E·α)` overall. Both are exact maximum-cardinality algorithms, so the result
size is the matching number `ν(g)`; the particular maximum matching is not
specified (the greedy start and search order fix it). The companion
`FindEdgeCover` reuses the matching through the Gallai identity (a maximum
matching plus one incident edge per exposed vertex, size `n − ν(g)`).

- `Protected`; unevaluated on a non-graph.
- Edge direction is ignored; edges are returned in `EdgeList` order.
- Algorithm: Karp-Sipser greedy start, then Hopcroft-Karp when the graph is
  bipartite and Edmonds' blossom algorithm otherwise (union-find blossom bases;
  failed searches retire their Hungarian trees). Matching on a triangulated grid:
  about 2.8 ms against 190 ms in Mathematica.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.
  (The differential test found Mathematica returning non-maximum matchings on
  weighted graphs.)

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- J. E. Hopcroft and R. M. Karp, *An n^{5/2} algorithm for maximum matchings in bipartite graphs*, SIAM J. Comput. **2** (1973) 225-231.
- J. Edmonds, *Paths, trees, and flowers*, Canad. J. Math. **17** (1965) 449-467.
- R. M. Karp and M. Sipser, *Maximum matchings in sparse random graphs*, FOCS 1981, 364-375.
- Source: [`src/graph/galg_matching.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_matching.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`FindIndependentEdgeSet[g]` returns a maximum **matching** — a largest set of
pairwise non-adjacent edges — as `g`'s own edge expressions in `EdgeList` order.
Edge direction is ignored. The size of the result is the matching number
`ν(g)`; the particular maximum matching returned is not specified.

The solver is exact: a Karp-Sipser greedy start, then Hopcroft-Karp for
bipartite graphs and Edmonds' blossom algorithm for general ones. Because K4 and
the 4-path both admit a perfect (resp. near-perfect) matching, the odd 5-cycle
can cover only four of its five vertices, leaving one exposed.
