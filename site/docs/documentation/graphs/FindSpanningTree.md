# FindSpanningTree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindSpanningTree[g] gives a spanning tree (forest) of g as a graph: a minimum spanning tree when g has EdgeWeight (an arborescence for a directed graph), else a breadth-first one. FindSpanningTree[{g, v}] gives the tree grown from vertex v.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[FindSpanningTree[CycleGraph[4]]]
Out[1]= {1 <-> 2, 1 <-> 4, 2 <-> 3}

In[2]:= EdgeCount[FindSpanningTree[CompleteGraph[6]]]
Out[2]= 5

In[3]:= EdgeList[FindSpanningTree[Graph[{1,2,3,4},{1->2,2->3,3->1}]]]
Out[3]= {1 -> 2, 2 -> 3}

In[4]:= InputForm[FindSpanningTree[{Graph[{1,2,3,4},{1<->2,3<->4}], 4}]]
Out[4]= Graph[{3, 4}, {3 <-> 4}]
```

### Options (3)

```mathematica
In[5]:= InputForm[FindSpanningTree[Graph[{a<->b, a<->c, b<->c, b<->d, c<->d, d<->e}, EdgeWeight -> {4,1,2,5,8,3}]]]
Out[5]= Graph[{a, b, c, d, e}, {a <-> c, b <-> c, b <-> d, d <-> e}, EdgeWeight -> {1, 2, 5, 3}]

In[6]:= InputForm[FindSpanningTree[Graph[{1->2, 2->3, 1->3, 3->4, 4->1}, EdgeWeight -> {5,1,2,3,4}]]]
Out[6]= Graph[{1, 2, 3, 4}, {2 -> 3, 3 -> 4, 4 -> 1}, EdgeWeight -> {1, 3, 4}]

In[7]:= InputForm[FindSpanningTree[{Graph[{1->2, 2->3, 1->3, 3->4, 4->1}, EdgeWeight -> {5,1,2,3,4}], 1}]]
Out[7]= Graph[{1, 2, 3, 4}, {1 -> 2, 2 -> 3, 3 -> 4}, EdgeWeight -> {5, 1, 3}]
```

### Applications (2)

Drops one edge to break the cycle

```mathematica
In[8]:= EdgeList[FindSpanningTree[CycleGraph[4]]]
Out[8]= {1 <-> 2, 1 <-> 4, 2 <-> 3}
```

A spanning tree has n-1 edges

```mathematica
In[9]:= EdgeCount[FindSpanningTree[CompleteGraph[5]]]
Out[9]= 4
```

## Algorithm

spanningtree.c - FindSpanningTree[g], FindSpanningTree[{g, v}].

Mathematica's semantics, checked case by case against Mathematica 15:

```text
  undirected, unweighted   a BFS spanning forest;
  undirected, EdgeWeight   a MINIMUM spanning forest (Kruskal): ties between
                           equal weights are broken by the edge's vertex
                           positions, lowest first, which reproduces
                           Mathematica's choice on tied inputs;
  directed,   unweighted   a BFS branching (edges followed forwards), its
                           roots taken in decreasing DFS finishing time, so
                           the forest needs as few roots as possible;
  directed,   EdgeWeight   a MINIMUM-WEIGHT spanning branching (Chu-Liu /
                           Edmonds) among those with the fewest roots -- a
                           minimum spanning arborescence whenever one vertex
                           reaches all;
  mixed                    left unevaluated, as Mathematica does ("not
                           implemented").
```

{g, v} restricts the answer to the tree grown from v: v's component (undirected) or the vertices v reaches (directed), rooted at v; the result's VertexList is just those vertices. A v that is not a vertex of g emits FindSpanningTree::inv and stays unevaluated. Trailing options (Method -> ...) are accepted and ignored: every method yields the same optimum.

Weights are compared EXACTLY. Integer, Rational, Real and arbitrary-precision weights become GMP rationals (a Real's exact binary value), so 1/3 and 0.3333333333333333 are told apart and nothing is rounded; other numeric expressions (Sqrt[2], Pi) are compared through N[w, 40]. A weight that is not a real number (a symbol, a Complex) leaves the call unevaluated, as in Mathematica.

The result is Graph[verts, treeEdges(, EdgeWeight -> w)(, EdgeCapacity -> c)] with each tree edge's own weight and capacity carried over. As in Mathematica, an undirected tree edge is written lower VertexList position first and the edges are sorted by the positions of their endpoints.

Everything runs on the validated-graph memo's integer endpoint arrays: no expression is hashed or compared after the weights are read.

Memory (SPEC section 4): returns a freshly allocated Graph built from copies of g's parts; res is never modified (the evaluator frees it on success and keeps it on NULL).

## Implementation notes

**Algorithm.** `builtin_find_spanning_tree` returns a spanning tree/forest as a
`Graph`, choosing the method from the graph's kind (checked case by case against
Mathematica 15): an **undirected unweighted** graph gives a BFS spanning forest;
an **undirected weighted** graph a minimum spanning forest by **Kruskal** (ties
broken by endpoint positions, reproducing Mathematica's choice); a **directed
unweighted** graph a BFS branching with roots taken in decreasing DFS finishing
time (fewest roots); a **directed weighted** graph a minimum-weight spanning
branching by **Chu-Liu/Edmonds** (the `O(E log V)` contraction formulation with a
rollback union-find and leftist mergeable heaps carrying a lazy additive delta).
A mixed graph is left unevaluated. `FindSpanningTree[{g, v}]` restricts the tree
to `v`'s component/reach, rooted at `v`.

**Data structures.** An out-incidence CSR with edge ids (`Inc`), and for the
weighted paths exact **GMP rational** weights (`mpq_t`): Integer/Rational/Real/
MPFR become exact binary rationals so `1/3` and `0.3333333333333333` are told
apart, other numerics go through `N[w, 40]`. Everything else runs on the memo's
integer endpoint arrays; `build_tree` re-emits edges in sorted endpoint order and
seeds the result into the graph memo.

**Complexity / limits.** `O(V + E)` unweighted, `O(E log V)` weighted. Weights
are compared exactly; a non-real weight (symbol, Complex) leaves the call
unevaluated. Trailing `Method -> ...` options are accepted and ignored — every
method yields the same optimum.

- `Protected`. Matches Mathematica case by case:
  - undirected, unweighted: a BFS spanning forest;
  - undirected, weighted: a minimum spanning forest (Kruskal with union-find),
    ties between equal weights broken by the edge's vertex positions, lowest
    first — which reproduces Mathematica's choice on tied inputs;
  - directed, unweighted: a BFS *branching* (edges followed forwards), roots
    taken in decreasing DFS finishing time so as few roots as possible are used;
  - directed, weighted: a minimum-weight spanning branching (Chu–Liu/Edmonds,
    `O(E log V)`) among those with the fewest roots — the minimum spanning
    arborescence whenever some vertex reaches all;
  - mixed graphs are left unevaluated, as in Mathematica.
- Weights are compared **exactly**: integers, rationals and reals become GMP
  rationals (a real's exact binary value), so `1/3` and `0.3333333333333333`
  are told apart; other numeric weights (`Sqrt[2]`, `Pi`) are compared via
  `N[w, 40]`. A weight that is not a real number (a symbol, a complex number)
  leaves the call unevaluated. Negative weights are fine.
- The result carries each tree edge's own `EdgeWeight` (and `EdgeCapacity`).
  As in Mathematica, an undirected tree edge is written lower vertex position
  first, and the edges are sorted by their endpoints' positions. The result's
  `VertexList` is `g`'s, or for `{g, v}` just the vertices of `v`'s tree.
- A `v` that is not a vertex of `g` emits `FindSpanningTree::inv` and leaves the
  call unevaluated; so does a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/), [Pi](../../mathematical-constants/Pi/), [VertexList](../../graphs/VertexList/)

- J. B. Kruskal, *On the shortest spanning subtree of a graph and the traveling salesman problem*, Proc. Amer. Math. Soc. **7** (1956) 48-50.
- Y. J. Chu and T. H. Liu, *On the shortest arborescence of a directed graph*, Sci. Sinica **14** (1965) 1396-1400; J. Edmonds, *Optimum branchings*, J. Res. NBS **71B** (1967) 233-240.
- Source: [`src/graph/spanningtree.c`](https://github.com/stblake/mathilda/blob/main/src/graph/spanningtree.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

The result is a `Graph` on the same vertices whose edges form a spanning tree
(or forest, if the graph is disconnected). A tree on `n` vertices has `n - 1`
edges.

When the graph carries `EdgeWeight`, the tree is a *minimum* spanning tree —
Kruskal for undirected graphs, Chu-Liu/Edmonds for directed ones — and weights
are compared exactly, so an integer and the nearest machine real are never
confused. A mixed graph is left unevaluated.
