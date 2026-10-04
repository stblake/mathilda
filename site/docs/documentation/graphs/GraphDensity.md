# GraphDensity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDensity[g] gives the number of edges of g divided by the number of possible edges: (directed edges + 2 undirected edges)/(n(n-1)).`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphDensity[PathGraph[4]]
Out[1]= 1/2

In[2]:= GraphDensity[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= 1/3

In[3]:= GraphDensity[CompleteGraph[5]]
Out[3]= 1

In[4]:= GraphDensity[Graph[{1<->2, 2->3}]]
Out[4]= 1/2

In[5]:= GraphDensity[Graph[{1},{}]]
Out[5]= GraphDensity[Graph[<1 vertex, 0 edges>]]
```

### Applications (4)

A complete graph has density 1

```mathematica
In[6]:= GraphDensity[CompleteGraph[4]]
Out[6]= 1
```

Five edges out of ten possible

```mathematica
In[7]:= GraphDensity[CycleGraph[5]]
Out[7]= 1/2
```

Directed arcs each count once

```mathematica
In[8]:= GraphDensity[Graph[{1 -> 2, 2 -> 3}]]
Out[8]= 1/3
```

A complete bipartite graph

```mathematica
In[9]:= GraphDensity[CompleteGraph[{2, 3}]]
Out[9]= 3/5
```

## Implementation notes

**Algorithm.** `builtin_graph_density` reads three numbers straight off the `Graph`
expression and does no traversal. It counts the directed arcs `nd` with
`graph_directed_edge_count`, reads the total edge count `ne` and the vertex count
`n` from the lengths of the graph's vertex-list and edge-list arguments, and
returns the exact rational `(nd + 2(ne − nd)) / (n(n−1))`. Each directed arc
contributes one to the numerator and each undirected edge two, so the numerator is
the number of ordered endpoint pairs joined by an edge; the denominator
`n(n−1)` is the number of ordered pairs of distinct vertices, i.e. the maximum
possible directed arcs on `n` labelled vertices. A complete undirected graph
therefore has density `1`.

**Data structures.** None beyond the input `Expr` tree: the arg counts are read
directly, and `make_rational` builds the result, reducing it to lowest terms.

**Complexity / limits.** `O(E)` to classify the edges as directed or undirected,
then `O(1)` arithmetic; the answer is an exact `Rational` with no floating point.
Graphs with fewer than two vertices leave the call unevaluated (the density is
undefined — `n(n−1) = 0`), matching Wolfram. Self-loops, if present, are counted
among the edges by `ne`.

- Computed as `(directed edges + 2 undirected edges)/(n (n - 1))`, exact.
- Weights ignored.
- Unevaluated for `n < 2`.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The density is the number of edges present divided by the number that could be
present: `(directed arcs + 2 × undirected edges) / (n(n−1))`. An undirected edge
counts twice because it fills both ordered endpoint pairs, so a complete undirected
graph scores exactly `1` and a sparse graph scores near `0`.

The result is an exact rational, computed from the vertex and edge counts alone
with no traversal. A graph on fewer than two vertices is left unevaluated, since
`n(n−1)` would be zero and the density is undefined.
