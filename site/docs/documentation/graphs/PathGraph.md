# PathGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PathGraph[n] gives the path on n vertices; PathGraph[{v1,...}] the path over the given vertices.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexDegree[PathGraph[5]]
Out[1]= {1, 2, 2, 2, 1}

In[2]:= EdgeList[PathGraph[4]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 4}

In[3]:= EdgeList[PathGraph[{a,b,c}]]
Out[3]= {a <-> b, b <-> c}

In[4]:= PathGraph[1]
Out[4]= Graph[<1 vertex, 0 edges>]

In[5]:= PathGraph[x]
Out[5]= PathGraph[x]
```

### Applications (3)

A chain 1 - 2 - 3 - 4 - 5

```mathematica
In[6]:= EdgeList[PathGraph[5]]
Out[6]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}
```

A path over the given vertices

```mathematica
In[7]:= VertexList[PathGraph[{a, b, c}]]
Out[7]= {a, b, c}
```

```mathematica
In[8]:= EdgeList[PathGraph[{a, b, c}]]
Out[8]= {a <-> b, b <-> c}
```

## Implementation notes

**Algorithm.** `builtin_path_graph` builds an undirected path. `PathGraph[n]` with an integer
argument gives the vertices `1..n` with edges `i <-> i+1` for `i = 1..n-1`. `PathGraph[{v1, ...,
vk}]` with a list argument uses the given expressions as vertices (copied) and joins consecutive
ones, `v_i <-> v_{i+1}`. The function assembles a raw `Graph[List[verts], List[edges]]`
expression, which the evaluator then re-runs through `builtin_graph` to canonicalise and
validate.

**Data structures.** The vertex and edge `List`s are built directly; all edges are
`UndirectedEdge` nodes. No weights.

**Complexity / limits.** `O(n)`. The integer form needs a non-negative integer (negative or
non-integer returns unevaluated); `PathGraph[0]` is the empty path.

- `Protected`. Undirected edges, built through the `Graph` constructor (see
  `CycleGraph`). `PathGraph[1]` is a single vertex; a symbolic argument is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [CycleGraph](../../graphs/CycleGraph/)

- Source: [`src/graph/generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`PathGraph[n]` is the undirected path on vertices `1..n`, with `n-1` edges joining consecutive
vertices. `PathGraph[{v1, ..., vk}]` builds the same chain over the given vertex expressions,
which may be any symbols or values.

The result is a canonical `Graph` object; query it with `VertexList`, `EdgeList`, `VertexCount`
or `EdgeCount` rather than relying on its printed form.
