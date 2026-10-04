# CompleteGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompleteGraph[n] gives the complete graph on n vertices. CompleteGraph[{n1, n2, ...}] gives the complete multipartite graph with parts of sizes n1, n2, ....`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= CompleteGraph[5]
Out[1]= Graph[<5 vertices, 10 edges>]

In[2]:= CompleteGraph[{2, 3}]
Out[2]= Graph[<5 vertices, 6 edges>]

In[3]:= EdgeList[CompleteGraph[{2, 3}]]
Out[3]= {1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 4, 2 <-> 5}

In[4]:= EdgeList[CompleteGraph[{1, 1, 2}]]
Out[4]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4}

In[5]:= CompleteGraph[{4}]
Out[5]= Graph[<4 vertices, 6 edges>]

In[6]:= CompleteGraph[x]
Out[6]= CompleteGraph[x]
```

### Applications (4)

K5 has n(n-1)/2 = 10 edges

```mathematica
In[7]:= EdgeCount[CompleteGraph[5]]
Out[7]= 10
```

Vertices are 1..n

```mathematica
In[8]:= VertexCount[CompleteGraph[4]]
Out[8]= 4
```

Every pair joined, in row-major order

```mathematica
In[9]:= EdgeList[CompleteGraph[3]]
Out[9]= {1 <-> 2, 1 <-> 3, 2 <-> 3}
```

Complete bipartite K(2,3): 2*3 = 6 edges

```mathematica
In[10]:= EdgeCount[CompleteGraph[{2, 3}]]
Out[10]= 6
```

## Implementation notes

**Algorithm.** `builtin_complete_graph` builds the undirected complete graph
`K_n` on the vertices `1..n`: it emits every pair `i < j` as an
`UndirectedEdge[i, j]`, giving exactly `n(n-1)/2` edges. The multipartite form
`CompleteGraph[{n1, n2, ...}]` is handled by a second registration,
`builtin_gmet_complete_graph` in `src/graph/gmet_generators.c`, which wraps the
base builtin: it lays the parts out as consecutive vertex blocks and joins every
pair of vertices in *different* parts, so `CompleteGraph[{2, 3}]` is the complete
bipartite graph on `2 + 3` vertices.

**Data structures.** Both forms assemble plain C arrays of `Expr*` for the
vertices and edges and hand them to `expr_new_function` as a
`Graph[List, List]`; the evaluator's `builtin_graph` canonicalizes and validates
the result and seeds the graph memo. Edges are generated in row-major pair order
(`undirected_edge(i, j)` shares a single interned `UndirectedEdge` head path).

**Complexity / limits.** `O(n^2)` — the output has `Theta(n^2)` edges, so the
cost is output-bound. A negative or non-integer `n` leaves the call unevaluated.

- All graph families of this module are undirected on `1..N`; the edge list is
  the sorted list of pairs `{i, j}`, `i < j` — identical to Wolfram's
  `EdgeList` for each family.
- `CompleteGraph[{n}]` is `K_n`.
- `CompleteGraph` is re-registered by a wrapper that delegates its
  pre-existing form (`CompleteGraph[n]`) to the original builtin.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: a family with more than 10^8 vertices or 5×10^7 edges is
  left unevaluated rather than allocating gigabytes.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`CompleteGraph[n]` joins every pair of the `n` vertices, so it has `n(n-1)/2`
edges and is the densest simple graph on `n` vertices.

`CompleteGraph[{n1, n2, ...}]` is the complete multipartite graph: the vertices
split into blocks of the given sizes and every pair in *different* blocks is
joined, none within a block. With two parts this is the complete bipartite
graph.
