# FindIndependentVertexSet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindIndependentVertexSet[g] gives {s} with s a maximum independent vertex set of g. FindIndependentVertexSet[g, k] finds a maximal independent set of at most k vertices, [g, {k}] of exactly k, [g, {kmin, kmax}] within the range; a third argument n (or All) gives up to n such sets. Exact; unevaluated if the search budget is exhausted.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindIndependentVertexSet[PetersenGraph[]]
Out[1]= {{2, 3, 6, 9}}

In[2]:= FindIndependentVertexSet[CycleGraph[6], {2}, All]
Out[2]= {{3, 6}, {2, 5}, {1, 4}}

In[3]:= FindIndependentVertexSet[CycleGraph[5], Infinity, All]
Out[3]= {{3, 5}, {2, 5}, {2, 4}, {1, 4}, {1, 3}}

In[4]:= FindIndependentVertexSet[PathGraph[{1,2,3,4}], {1,2}, 2]
Out[4]= {{1, 4}, {1, 3}}

In[5]:= FindIndependentVertexSet[CompleteGraph[4], 3]
Out[5]= {{1}}

In[6]:= FindIndependentVertexSet[x]
Out[6]= FindIndependentVertexSet[x]
```

### Applications (5)

One maximum set, wrapped in a list

```mathematica
In[7]:= FindIndependentVertexSet[PathGraph[{a, b, c, d, e}]]
Out[7]= {{a, c, e}}
```

The leaves

```mathematica
In[8]:= FindIndependentVertexSet[StarGraph[5]]
Out[8]= {{2, 3, 4, 5}}
```

Alpha of a five-cycle is two

```mathematica
In[9]:= FindIndependentVertexSet[CycleGraph[5]]
Out[9]= {{3, 5}}
```

Alpha of the petersen graph is four

```mathematica
In[10]:= FindIndependentVertexSet[PetersenGraph[]]
Out[10]= {{2, 3, 6, 9}}
```

Maximal sets of size two, all of them

```mathematica
In[11]:= FindIndependentVertexSet[CycleGraph[6], {2}, All]
Out[11]= {{3, 6}, {2, 5}, {1, 4}}
```

## Implementation notes

**Algorithm.** `builtin_find_independent_vertex_set` ignores edge direction. With one argument it calls the exact solver `galg_max_independent_set` and returns `{s}`, a single *maximum* independent set. With a size spec (`k`, `{k}`, `{kmin, kmax}`) and an optional count (`n` or `All`) it instead enumerates *maximal* independent sets by size through `galg_clique_spec_query`, the same engine `FindClique` uses on the complement. The solver works per connected component:

1. a bipartite component of at least 256 vertices (grids, trees, even cycles) is solved in `O(m sqrt n)` by König's theorem (`gm_konig`);
2. a dense component (average degree >= 8 or density >= 0.05, at most 3000 vertices) goes to the bitset maximum-clique branch and bound on its complement;
3. a sparse component uses branch and reduce, seeded with a greedy minimum-degree incumbent. Reductions run to a fixed point: degree 0 and 1, degree-2 triangle, degree-2 *folding* (replace `{v, a, b}` by one vertex adjacent to `N(a) u N(b)`, alpha drops by exactly 1), and domination (`N[v]` inside `N[u]` drops `u`). The bound "taken so far + greedy clique cover of the rest" prunes, a disconnected remainder splits, and otherwise the search branches on a maximum-degree vertex: drop it with its mirrors (Fomin-Grandoni-Kratsch) or take it.

**Data structures.** A `GalgUG` CSR view built from the `Graph[List, List]` expression tree (edges memoized by `graph_edge_indices`), then a `GmMis` search state: growable adjacency lists, a removal/fold trail so every reduction is undone exactly, a solution stack, and stamped scratch arrays. Dense components use bitsets. Results are rebuilt as vertex expressions from the integer positions.

**Complexity / limits.** Maximum independent set is NP-hard; the worst case is exponential. The solver gives up, and the head stays unevaluated, after 2e7 search nodes, 6e8 aggregate work, 64M ints of per-level scratch, 20000 levels of recursion, or when `TimeConstrained`'s deadline passes. It never returns a merely maximal set for the one-argument form. Which optimum is returned among ties is fixed by the vertex order and the reductions, so repeated calls agree.

- `Protected`; unevaluated on a non-graph.
- The size-spec forms enumerate MAXIMAL independent sets by size, largest
  first, exactly like the `FindClique` spec forms; they are limited to 8192
  vertices.
- Exact (see the exactness policy under `FindVertexCover`): connected
  components are solved separately.
  - A bipartite component of at least 256 vertices (grids, meshes, trees, even
    cycles) is solved in O(m sqrt n) by König's theorem (Hopcroft-Karp matching,
    then alternating reachability); `GridGraph[{1000, 1000}]` takes about 1.2 s.
  - Otherwise a component with average degree >= 8 (or density >= 0.05, up to
    3000 vertices) goes to the bitset maximum-clique search on its complement.
  - Sparser ones go to branch and reduce: degree-0/1 and triangle reductions,
    degree-2 folding, domination, a greedy clique-cover upper bound, component
    splitting at every node, and branching on a maximum-degree vertex with its
    mirrors.
- Budgets: the search gives up -- the head stays unevaluated, never a merely
  maximal set -- after 20 million nodes, 6 x 10^8 units of aggregate work
  (live-set size summed over nodes), 64M ints of per-level scratch or depth
  20000, so a huge non-bipartite sparse input (e.g.
  `RandomGraph[{200000, 300000}]`) returns unevaluated in bounded memory rather
  than exhausting it. The `TimeConstrained` deadline is polled.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.

**Attributes:** `Protected`.

## References

**See also:** [FindClique](../../graphs/FindClique/), [FindVertexCover](../../graphs/FindVertexCover/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- F. V. Fomin, F. Grandoni and D. Kratsch, *A measure and conquer approach for the analysis of exact algorithms*, J. ACM **56**(5) (2009) Art. 25.
- Source: [`src/graph/galg_mis.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_mis.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

With only a graph the head returns `{s}`, one maximum independent set. With a size specification (`k`, `{k}` or `{kmin, kmax}`) and an optional count (`n` or `All`) it lists *maximal* independent sets of those sizes, as `FindClique` does.

Which maximum set is returned among ties is deterministic for a given graph. The search is exact but NP-hard in general; past its node/work budget, or a `TimeConstrained` limit, the call stays unevaluated.
