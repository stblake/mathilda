# FindIndependentVertexSet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindIndependentVertexSet[g] gives {s} with s a maximum independent vertex set of g. FindIndependentVertexSet[g, k] finds a maximal independent set of at most k vertices, [g, {k}] of exactly k, [g, {kmin, kmax}] within the range; a third argument n (or All) gives up to n such sets. Exact; unevaluated if the search budget is exhausted.`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
