# GraphIntersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphIntersection[g1, g2, ...] gives the graph on the union of the vertex sets whose edges are those common to all the gi, in canonical order. Weights are dropped.`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
