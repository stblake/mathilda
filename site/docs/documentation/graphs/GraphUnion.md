# GraphUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphUnion[g1, g2, ...] gives the graph whose vertices and edges are the unions of those of the gi (vertices in canonical order; an undirected edge equals its reversal). Weights are dropped.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= EdgeList[GraphUnion[Graph[{2<->1}], Graph[{3<->4}]]]
Out[1]= {1 <-> 2, 3 <-> 4}

In[2]:= EdgeList[GraphUnion[CycleGraph[3], Graph[{2<->1, 3<->4}]]]
Out[2]= {1 <-> 2, 2 <-> 3, 1 <-> 3, 3 <-> 4}

In[3]:= EdgeList[GraphUnion[Graph[{2->1}], Graph[{1<->3}]]]
Out[3]= {2 -> 1, 1 <-> 3}

In[4]:= EdgeList[GraphIntersection[CompleteGraph[4], CycleGraph[4], PathGraph[Range[4]]]]
Out[4]= {1 <-> 2, 2 <-> 3, 3 <-> 4}

In[5]:= EdgeList[GraphDifference[CompleteGraph[4], CycleGraph[4]]]
Out[5]= {1 <-> 3, 2 <-> 4}

In[6]:= GraphUnion[CycleGraph[3], 5]
Out[6]= GraphUnion[Graph[<3 vertices, 3 edges>], 5]
```

### Applications (4)

The vertices are the set-union 1..5

```mathematica
In[7]:= VertexCount[GraphUnion[CycleGraph[3], PathGraph[{3, 4, 5}]]]
Out[7]= 5
```

Every distinct edge of either graph

```mathematica
In[8]:= EdgeList[GraphUnion[CycleGraph[3], PathGraph[{3, 4, 5}]]]
Out[8]= {1 <-> 2, 2 <-> 3, 1 <-> 3, 3 <-> 4, 4 <-> 5}
```

The shared edge 2 <-> 3 appears once

```mathematica
In[9]:= EdgeList[GraphUnion[PathGraph[{1, 2, 3}], PathGraph[{2, 3, 4}]]]
Out[9]= {1 <-> 2, 2 <-> 3, 3 <-> 4}
```

A graph unioned with itself is itself

```mathematica
In[10]:= EdgeCount[GraphUnion[CompleteGraph[3], CompleteGraph[3]]]
Out[10]= 3
```

## Implementation notes

**Algorithm.** `builtin_graph_union` requires every argument to be a valid graph;
`GraphUnion[g]` returns a copy of `g`. The vertex set of the result is the union of
all the graphs' vertices in canonical (`Sort`) order. `union_build` maps each
graph's vertices to union ids — by an `O(V)` elementwise `SameQ` check when a
vertex list equals the first graph's (the common case of graphs sharing one vertex
set: no hashing at all), else through a single `GraphVIdx` hash over the union —
then orders the union with `expr_compare`, skipping the sort when it is already
sorted and running on machine integers when every vertex is one. The edges are the
**distinct** edges of all the graphs (an undirected edge equals its reversal), with
duplicates found in a `GopsKeySet` over integer edge keys. Edge order follows
Mathematica 15: all-undirected edges are oriented by canonical vertex order in
first-appearance order, all-directed edges keep first-appearance order, and a mixed
set is put in canonical order by a counting sort on `(kind, first, second)`.

**Data structures.** A `Union` struct holds the per-graph vertex-id maps and the
ordered union vertex array. Each kept edge is an `OutEdge` of result-position
endpoints plus a direction flag and the source edge node it may share. Deduplication
and canonical ordering are integer passes — a `GopsKeySet` for distinctness and a
stable counting sort (`gops_csort`) for order — so the whole assembly is `O(V + E)`
apart from the vertex sort. `union_graph` builds the result with `gops_graph_new`,
sharing a source edge node by `expr_copy` whenever its orientation is unchanged.

**Complexity / limits.** `O(V + E)` plus the `O(V log V)` vertex sort, which drops
to `O(V)` when the union is already sorted or every vertex is an integer. The result
is a simple graph; edge **weights are dropped**, as in Mathematica. `res` is
borrowed and never modified.

- `Protected`. A non-graph argument is left unevaluated.
- Vertices are always the union of the inputs' vertices, in canonical order.
- `GraphUnion` edges: the distinct edges (an undirected edge equals its
  reversal). All undirected: first-appearance order, each oriented by canonical
  vertex order. All directed: first-appearance order. Mixed: canonical (`Sort`)
  order. `GraphUnion[g]` is `g`.
- `GraphIntersection` / `GraphDifference` edges are in canonical order.
  `GraphDifference` takes exactly two graphs.
- Weights are dropped, as in Mathematica (the one-argument `GraphUnion[g]`
  returns `g` itself, so its weights survive).
- Implementation: vertex lists equal to the first graph's are mapped by an
  `O(V)` elementwise `SameQ` check (no hashing); otherwise through one hash index
  over the union. The union is sorted with `expr_compare` only when not already
  sorted (machine integers are sorted as such), and edge keys over result
  positions are deduplicated in an integer hash set and ordered by stable
  counting sorts. The same machinery serves `GraphDisjointUnion`.

**Attributes:** `Protected`.

## References

**See also:** [GraphIntersection](../../graphs/GraphIntersection/), [GraphDifference](../../graphs/GraphDifference/), [Sort](../../data-structures/Sort/), [SameQ](../../comparisons/SameQ/), [GraphDisjointUnion](../../graphs/GraphDisjointUnion/)

- Source: [`src/graph/gops_setops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_setops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`GraphUnion[g1, g2, ...]` is the graph whose vertices are the union of all the
graphs' vertices (in canonical `Sort` order) and whose edges are every **distinct**
edge of any of them — an undirected edge and its reversal count as one, so a shared
edge appears a single time. `GraphUnion[g]` is just `g`. The result is an opaque
`Graph`, read through `VertexCount`, `EdgeCount`, `EdgeList`, and the other
accessors.

Edge order matches Mathematica 15 (first-appearance for a uniform direction,
canonical order when directed and undirected edges are mixed), and edge **weights
are dropped** — the union is an unweighted simple graph.
