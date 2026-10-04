# FindClique

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindClique[g] gives {c} with c a maximum clique of g. FindClique[g, k] finds a maximal clique of at most k vertices, [g, {k}] of exactly k, [g, {kmin, kmax}] within the range; a third argument n (or All) gives up to n such cliques. For directed graphs a clique needs edges both ways. Exact branch and bound with colouring bounds.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindClique[CompleteGraph[5]]
Out[1]= {{1, 2, 3, 4, 5}}

In[2]:= FindClique[CycleGraph[5], Infinity, All]
Out[2]= {{4, 5}, {3, 4}, {2, 3}, {1, 5}, {1, 2}}

In[3]:= FindClique[WheelGraph[6], {3}, 2]
Out[3]= {{1, 2, 6}, {1, 2, 3}}

In[4]:= FindClique[Graph[{1,2,3,4},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[1,3],UndirectedEdge[3,4]}], {2,3}, All]
Out[4]= {{1, 2, 3}, {3, 4}}

In[5]:= FindClique[CompleteGraph[4], 2]
Out[5]= {}

In[6]:= FindClique[Graph[{1->2,2->1,2->3,3->2,1->3}]]
Out[6]= {{1, 2}}
```

### Scope (1)

```mathematica
In[7]:= TimeConstrained[FindClique[RandomGraph[{400, 40000}]], 0.001]
Out[7]= $Aborted
```

### Applications (5)

The whole graph is one clique

```mathematica
In[8]:= FindClique[CompleteGraph[4]]
Out[8]= {{1, 2, 3, 4}}
```

The triangle, not the pendant edge

```mathematica
In[9]:= FindClique[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}]]
Out[9]= {{1, 2, 3}}
```

Triangle-free: the maximum clique is a single edge

```mathematica
In[10]:= FindClique[CycleGraph[5]]
Out[10]= {{1, 2}}
```

Maximal cliques of <= 2 vertices -- none, every maximal clique is K4

```mathematica
In[11]:= FindClique[CompleteGraph[4], 2]
Out[11]= {}
```

Every maximal clique of exactly 2 vertices

```mathematica
In[12]:= FindClique[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}], {2}, All]
Out[12]= {{3, 4}}
```

## Options & behaviour

A search that exceeds its `TimeConstrained` limit aborts instead of hanging:

## Implementation notes

**Algorithm.** `builtin_find_clique` first reduces `g` to the simple
mutual-adjacency graph `gc_mutual_graph` — undirected edges, and directed pairs
joined both ways, so a clique in a digraph needs the arcs in both directions.
With no size spec it calls `galg_max_clique`, a bitset branch-and-bound
(Tomita's MCS with San Segundo's BBMC bit-parallel colouring bound): vertices are
numbered by reverse degeneracy order (`gc_degeneracy`, an `O(n+m)` bucket peel),
greedy colouring supplies the pruning bound `depth + colour > incumbent`, and the
incumbent clique is returned sorted. A size spec (`k`, `{k}`, `{kmin,kmax}`, with
an optional count or `All`) switches to pivoted Bron-Kerbosch enumeration
(`gc_bk`), which lists **maximal** cliques in the size window largest-first; so
`FindClique[g, 2]` is `{}` when every maximal clique is larger. For `count == 1`
the enumeration tries sizes downward so the `kmin` prune bites hardest.

**Data structures.** Adjacency is a row-per-vertex bitset (`GcBits`, 64-bit
words), built induced on the working vertex set by `gc_bits_induced` (with an
optional complement for the independent-set sibling). Branch-and-bound keeps a
per-depth `P`-stack of candidate bitsets plus lazily allocated colouring scratch
(`O(ω·n)`, not `O(n²)`); Bron-Kerbosch concatenates each reported clique into a
growable `res`/`rstart` buffer and sorts the runs by `(size, lexicographic)` to
reproduce Wolfram's listing order. Above `GC_GLOBAL_MAX = 3000` vertices the
maximum-clique search is decomposed by degeneracy — every clique lies in
`{v} ∪ N⁺(v)` for its earliest vertex, a small local bitset problem skipped
outright when its core number cannot beat the incumbent.

**Complexity / limits.** Exponential in the worst case but pruned hard in
practice; capped at `GC_MAX_NODES = 5·10⁷` search nodes and `GC_BITSET_MAX = 8192`
vertices for the spec/complement paths, polled against the `TimeConstrained`
deadline every 4096 nodes. An exhausted budget leaves the call **unevaluated** —
never a non-maximum answer. `FindClique[g]` returns `{c}` (a one-element list
holding one maximum clique); the branch-and-bound proves it maximum but, like
Mathematica, does not specify which maximum clique.

- `Protected`; unevaluated on a non-graph.
- A directed graph needs edges both ways between clique members.
- Maximal cliques are listed by size, largest first; within one size
  Mathematica's order is reproduced
  (`FindClique[CycleGraph[5], Infinity, All]` gives
  `{{4,5},{3,4},{2,3},{1,5},{1,2}}`).
- `FindClique[g, 2]` is `{}` when every maximal clique is larger.
- Exact (see the exactness policy under `FindVertexCover`): a proven maximum or
  unevaluated when the deterministic node budget runs out; the
  `TimeConstrained` deadline is polled, so a timed-out search gives `$Aborted`.
- Maximum clique: bitset branch and bound with greedy-colouring bounds (Tomita's
  MCQ/MCS family in San Segundo's BBMC form), vertices numbered by reverse
  degeneracy order; graphs above 3000 vertices are decomposed by degeneracy
  (each vertex's later neighbourhood is a small local problem, skipped when its
  core number cannot beat the incumbent).
- Enumeration: pivoted Bron-Kerbosch on bitsets with size pruning.

**Attributes:** `Protected`.

## References

**See also:** [FindVertexCover](../../graphs/FindVertexCover/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- E. Tomita, Y. Sutani, T. Higashi, S. Takahashi and M. Wakatsuki, *A simple and faster branch-and-bound algorithm for finding a maximum clique*, WALCOM 2010, LNCS 5942, 191-203.
- P. San Segundo, D. Rodríguez-Losada and A. Jiménez, *An exact bit-parallel algorithm for the maximum clique problem*, Comput. Oper. Res. **38** (2011) 571-581.
- C. Bron and J. Kerbosch, *Algorithm 457: finding all cliques of an undirected graph*, Comm. ACM **16** (1973) 575-577.
- Source: [`src/graph/galg_clique.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_clique.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`FindClique[g]` returns a one-element list `{c}` holding one maximum clique, each
clique in `VertexList` order. The maximum is proved by an exact bitset
branch-and-bound (Tomita MCS with a bit-parallel colouring bound); which maximum
clique is returned is not specified, as in Mathematica.

The size-spec forms report **maximal** cliques — not extendable — so
`FindClique[g, 2]` is empty whenever every maximal clique is larger than two
vertices. A spec may be `k` (at most `k`), `{k}` (exactly `k`) or `{kmin, kmax}`,
and a trailing count or `All` asks for several, largest first. In a directed
graph a clique needs the arcs in both directions (an `UndirectedEdge` counts as
both). The search is budgeted and leaves the call unevaluated on exhaustion
rather than returning a non-maximum clique.
