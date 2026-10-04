# StarGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StarGraph[n] gives the star on n vertices: the hub 1 joined to each of the n-1 leaves 2..n.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[StarGraph[4]]
Out[1]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[2]:= VertexDegree[StarGraph[5]]
Out[2]= {4, 1, 1, 1, 1}

In[3]:= StarGraph[1]
Out[3]= Graph[<1 vertex, 0 edges>]

In[4]:= StarGraph[x]
Out[4]= StarGraph[x]
```

### Applications (6)

Hub 1 joined to four leaves

```mathematica
In[5]:= EdgeList[StarGraph[5]]
Out[5]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5}
```

The hub has degree n - 1, each leaf degree 1

```mathematica
In[6]:= VertexDegree[StarGraph[6]]
Out[6]= {5, 1, 1, 1, 1, 1}
```

A single vertex has no edges

```mathematica
In[7]:= EdgeCount[StarGraph[1]]
Out[7]= 0
```

Only the first row and column carry ones

```mathematica
In[8]:= AdjacencyMatrix[StarGraph[4]]
Out[8]= {{0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 0}}
```

The result is an ordinary Graph expression

```mathematica
In[9]:= GraphQ[StarGraph[8]]
Out[9]= True
```

Cheap even for large n, with n - 1 edges

```mathematica
In[10]:= VertexCount[StarGraph[100]]
Out[10]= 100
```

## Implementation notes

**Algorithm.** `builtin_star_graph[n]` reads a non-negative integer count via `as_count` and joins hub vertex 1 to each of vertices `2..n`, giving exactly `n - 1` edges `UndirectedEdge[1, i]`. `StarGraph[1]` is one isolated vertex and `StarGraph[0]` is the empty graph.

**Data structures.** A `calloc`'d array of `n - 1` edge nodes built by `undirected_edge` and an integer vertex list from `int_vertices`, passed to `make_graph`, which assembles the canonical `Graph[List[1..n], List[edges]]` `Expr` tree. No duplicate edge is possible, so no dedup pass is needed (unlike `CycleGraph`'s wrap edge).

**Complexity / limits.** `O(n)` time and memory. A non-integer or negative argument, or a wrong argument count, leaves the call unevaluated.

- `Protected`. Undirected edges, built through the `Graph` constructor (see
  `CycleGraph`). `StarGraph[1]` is a single vertex; a symbolic argument is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [CycleGraph](../../graphs/CycleGraph/)

- Source: [`src/graph/generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`StarGraph[n]` has `n` vertices with vertex 1 as the hub joined to each of `2..n`; it is the complete bipartite graph `K(1, n-1)`. `StarGraph[1]` is one isolated vertex.

The argument must be a non-negative integer; anything else leaves the call unevaluated.
