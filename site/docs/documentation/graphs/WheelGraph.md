# WheelGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WheelGraph[n] gives the wheel graph with n vertices: vertex 1 joined to every vertex of the cycle 2, ..., n.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= WheelGraph[5]
Out[1]= Graph[<5 vertices, 8 edges>]

In[2]:= EdgeList[WheelGraph[5]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}

In[3]:= WheelGraph[1]
Out[3]= Graph[<1 vertex, 0 edges>]

In[4]:= WheelGraph[3]
Out[4]= WheelGraph[3]
```

### Options (1)

```mathematica
In[5]:= WheelGraph[5, DirectedEdges -> True]
Out[5]= WheelGraph[5, DirectedEdges -> True]
```

### Applications (6)

The smallest proper wheel is the complete graph on 4 vertices

```mathematica
In[6]:= EdgeList[WheelGraph[4]]
Out[6]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4, 3 <-> 4}
```

A hub joined to a 4-cycle rim

```mathematica
In[7]:= EdgeList[WheelGraph[5]]
Out[7]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}
```

Hub degree n - 1, every rim vertex degree 3

```mathematica
In[8]:= VertexDegree[WheelGraph[6]]
Out[8]= {5, 3, 3, 3, 3, 3}
```

2 (n - 1) edges: spokes plus rim

```mathematica
In[9]:= EdgeCount[WheelGraph[10]]
Out[9]= 18
```

A single vertex and no edges

```mathematica
In[10]:= EdgeList[WheelGraph[1]]
Out[10]= {}
```

The hub plus a rim of six vertices

```mathematica
In[11]:= VertexCount[WheelGraph[7]]
Out[11]= 7
```

## Implementation notes

**Algorithm.** `builtin_wheel_graph[n]` joins hub vertex 1 to every other vertex (`n - 1` spokes) and, for `n >= 4`, adds the rim cycle on vertices `2..n`: edges `2-3`, `3-4`, ..., `(n-1)-n`, and the closing edge `2-n`. `WheelGraph[1]` is a single vertex. `n = 2` and `n = 3` would be multigraphs and are left unevaluated, as in Mathematica.

**Data structures.** Spokes and rim edges are appended to the packed-64-bit `Pairs` buffer in already-sorted order, then `pairs_graph` builds `Graph[Range[n], {UndirectedEdge[i, j], ...}]` with a shared `UndirectedEdge` head node and shared vertex integers.

**Complexity / limits.** `O(n)` time, `2(n - 1)` edges for `n >= 4`. Non-integer, non-positive, or over 10^8 vertex arguments leave the call unevaluated; options are not supported.

- Undirected on `1..n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- Accepts `n = 1` or `n >= 4`; `n = 2` and `3` would be multigraphs and are
  left unevaluated.
- Options (`DirectedEdges`, layout options) are not supported: a call with
  options is left unevaluated.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`WheelGraph[n]` has hub vertex 1 joined to every vertex of the cycle on `2..n`. For `n = 1` it is a single vertex; `n = 2` and `n = 3` would need repeated edges and are left unevaluated, as in Mathematica.

Vertex 1 always has degree `n - 1` and every rim vertex has degree 3.
