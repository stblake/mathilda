# GraphIntersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphIntersection[g1, g2, ...] gives the graph on the union of the vertex sets whose edges are those common to all the gi, in canonical order. Weights are dropped.`**

## Examples (9)

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

### Applications (3)

```mathematica
In[7]:= g1 = Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}];

In[8]:= g2 = Graph[{2 <-> 3, 3 <-> 4, 4 <-> 1}];
```

The edges shared by both

```mathematica
In[9]:= EdgeList[GraphIntersection[g1, g2]]
Out[9]= {2 <-> 3, 3 <-> 4}
```

## Implementation notes

**Algorithm.** `builtin_graph_intersection` builds the graph on the union of the argument graphs'
vertex sets whose edges are those common to *every* graph. Vertices are merged in canonical
(`expr_compare`) order. The edges of the first graph are tested for membership in each other
graph — an edge is kept when it is present in all of them (`hits == ng-1`) — and the kept edges
are counting-sorted into canonical order. The first graph's orientation is preserved, and edge
weights are dropped. `GraphIntersection[g]` of a single graph returns `g`.

**Data structures.** A `GopsView` per graph with integer endpoints; a `GopsKeySet` hash of the
first graph's edges, probed by each other graph's edges mapped onto the shared vertex indices;
the vertex union sort uses radix/counting for all-integer vertices and a comparison sort
otherwise.

**Complexity / limits.** `O(V + E)` apart from the `O(V log V)` vertex-union sort (`O(V)` when
vertices are already ordered integers). Variadic, at least one argument; a non-graph argument
returns unevaluated. An undirected edge equals its reversal under the edge key, so orientation
does not block a match.

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

**See also:** [GraphUnion](../../graphs/GraphUnion/), [GraphDifference](../../graphs/GraphDifference/), [Sort](../../data-structures/Sort/), [SameQ](../../comparisons/SameQ/), [GraphDisjointUnion](../../graphs/GraphDisjointUnion/)

- Source: [`src/graph/gops_setops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_setops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`GraphIntersection[g1, g2, ...]` gives the graph on the union of the vertex sets whose edges are
those common to all of the `gi`, in canonical order. Weights are dropped.

An undirected edge matches its reversal, so `1 <-> 2` and `2 <-> 1` are the same edge for the
set operation. The result is a canonical `Graph`; query it with `EdgeList` / `EdgeCount`.
