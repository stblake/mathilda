# HypercubeGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypercubeGraph[n] gives the n-dimensional hypercube graph on 2^n vertices.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HypercubeGraph[3]
Out[1]= Graph[<8 vertices, 12 edges>]

In[2]:= EdgeList[HypercubeGraph[2]]
Out[2]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4}

In[3]:= VertexDegree[HypercubeGraph[4]]
Out[3]= {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}

In[4]:= HypercubeGraph[0]
Out[4]= Graph[<1 vertex, 0 edges>]
```

### Applications (6)

The zero-dimensional cube is a single vertex

```mathematica
In[5]:= VertexCount[HypercubeGraph[0]]
Out[5]= 1
```

The cube of dimension one is one edge

```mathematica
In[6]:= EdgeList[HypercubeGraph[1]]
Out[6]= {1 <-> 2}
```

The square, a 4-cycle

```mathematica
In[7]:= EdgeList[HypercubeGraph[2]]
Out[7]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4}
```

Every vertex of the d-cube has degree d

```mathematica
In[8]:= VertexDegree[HypercubeGraph[3]]
Out[8]= {3, 3, 3, 3, 3, 3, 3, 3}
```

D times 2^(d-1) edges

```mathematica
In[9]:= EdgeCount[HypercubeGraph[4]]
Out[9]= 32
```

Vertices adjacent when their labels differ in one bit

```mathematica
In[10]:= AdjacencyMatrix[HypercubeGraph[2]]
Out[10]= {{0, 1, 1, 0}, {1, 0, 0, 1}, {1, 0, 0, 1}, {0, 1, 1, 0}}
```

## Implementation notes

**Algorithm.** `builtin_hypercube_graph[d]` builds the `d`-dimensional cube graph on `2^d` vertices. Vertex numbers are `1 + v` for bit patterns `v`; for each `v` and each bit `b` that is clear in `v`, the edge `v ~ v | (1 << b)` is emitted, so two vertices are adjacent exactly when their 0-based labels differ in one bit.

**Data structures.** Edges are written to the packed-64-bit `Pairs` buffer (one `uint64_t` key per edge) and `pairs_graph` sorts them lexicographically and constructs `Graph[Range[2^d], {UndirectedEdge[i, j], ...}]` as an ordinary `Expr` tree.

**Complexity / limits.** `O(d * 2^d)` time and `d * 2^(d-1)` edges. The dimension must be a machine integer in `0..24` (`HypercubeGraph[0]` is a single vertex); larger or non-integer arguments, and anything exceeding the 5 x 10^7 edge cap, are left unevaluated.

- Undirected on `1..2^n`; the edge list is the sorted list of pairs `{i, j}`,
  `i < j`, identical to Wolfram's `EdgeList`.
- `HypercubeGraph[0]` is a single vertex.
- Options (`DirectedEdges`, layout options) are not supported.
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

`HypercubeGraph[d]` has `2^d` vertices; vertices `i` and `j` are adjacent exactly when `i - 1` and `j - 1` differ in a single binary digit. It is `d`-regular and bipartite.

The dimension must be an integer from 0 to 24; anything else leaves the call unevaluated.
