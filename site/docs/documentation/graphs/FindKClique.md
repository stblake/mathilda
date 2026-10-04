# FindKClique

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindKClique[g, k] gives {c} with c a largest k-clique of g: a vertex set whose members are pairwise within distance k.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindKClique[CycleGraph[8], 2]
Out[1]= {{1, 2, 3}}

In[2]:= FindKClique[PathGraph[{1,2,3,4,5}], 1]
Out[2]= {{1, 2}}

In[3]:= FindKClique[x, 2]
Out[3]= FindKClique[x, 2]
```

### Applications (2)

K = 1 is an ordinary maximum clique

```mathematica
In[4]:= FindKClique[Graph[{1 <-> 2, 2 <-> 3, 1 <-> 3, 3 <-> 4}], 1]
Out[4]= {{1, 2, 3}}
```

K = 2: vertices pairwise within distance 2

```mathematica
In[5]:= FindKClique[PathGraph[5], 2]
Out[5]= {{3, 4, 5}}
```

## Implementation notes

**Algorithm.** `builtin_find_k_clique` finds a largest set of vertices that are pairwise within
graph distance `k` — a "`k`-clique" in the distance sense, which for `k = 1` is an ordinary
maximum clique. It first builds the `k`-th power graph (an edge between every pair at distance
`<= k`) by a truncated BFS to depth `k` from each vertex, then runs a **maximum clique** search
on that power graph. The clique engine is a bitset branch-and-bound in the Tomita/San-Segundo
(BBMC) style: vertices are ordered by reverse degeneracy, and a greedy graph-colouring of the
remaining candidates gives the bound that prunes a branch (`depth + colour[i] <= best`).

**Data structures.** A `GalgUG` CSR for the input; the power graph as a `GcBits` bitset matrix
(`n` rows of `w` 64-bit words); BFS `dist[]`/queue arrays to generate the power-graph edges; and
the branch-and-bound candidate sets as bitsets.

**Complexity / limits.** NP-hard in the worst case. The bitset engine is capped at `n <=
GC_BITSET_MAX = 8192` vertices and a node budget `GC_MAX_NODES = 5·10^7`, polled against
`TimeConstrained`. If the budget is exhausted the head returns unevaluated — never a
non-maximum clique. It returns `{c}` (one largest vertex set) or `{}` for the empty graph; a
non-graph argument, or a `k` that is not a positive integer, returns unevaluated.

- `Protected`; unevaluated on a non-graph.
- Computed as a maximum clique of the `k`-th power graph, using the same exact
  BBMC maximum-clique search as `FindClique` (proven optimum or unevaluated on
  budget exhaustion; polls `TimeConstrained`).

**Attributes:** `Protected`.

## References

**See also:** [FindClique](../../graphs/FindClique/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- E. Tomita, A. Tanaka and H. Takahashi, *The worst-case time complexity for generating all maximal cliques and computational experiments*, Theoret. Comput. Sci. **363** (2006) 28-42.
- P. San Segundo, D. Rodríguez-Losada and A. Jiménez, *An exact bit-parallel algorithm for the maximum clique problem*, Comput. Oper. Res. **38** (2011) 571-581.
- Source: [`src/graph/galg_clique.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_clique.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`FindKClique[g, k]` returns `{c}`, where `c` is a largest set of vertices that are pairwise
within graph distance `k`. With `k = 1` this is a maximum clique of `g`; larger `k` relaxes
adjacency to "close enough", so on a path the three consecutive vertices within distance 2 of
each other form the answer.

The empty graph gives `{}`. The search is exact: if the node budget is exhausted the call stays
unevaluated rather than return a set that might not be largest.
