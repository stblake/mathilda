# CycleGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CycleGraph[n] gives the cycle graph on n vertices.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EdgeList[CycleGraph[4]]
Out[1]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}

In[2]:= VertexCount[CycleGraph[10]]
Out[2]= 10

In[3]:= EdgeList[CycleGraph[2]]
Out[3]= {1 <-> 2}

In[4]:= EdgeCount[CompleteGraph[5]]
Out[4]= 10

In[5]:= CycleGraph[x]
Out[5]= CycleGraph[x]
```

### Applications (5)

N edges when n >= 3

```mathematica
In[6]:= EdgeCount[CycleGraph[5]]
Out[6]= 5
```

The rim, closing with 4 <-> 1

```mathematica
In[7]:= EdgeList[CycleGraph[4]]
Out[7]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}
```

Every vertex of a cycle has degree two

```mathematica
In[8]:= VertexDegree[CycleGraph[5]]
Out[8]= {2, 2, 2, 2, 2}
```

The integer vertices 1..n

```mathematica
In[9]:= VertexList[CycleGraph[6]]
Out[9]= {1, 2, 3, 4, 5, 6}
```

N <= 2 drops the duplicate wrap edge

```mathematica
In[10]:= EdgeCount[CycleGraph[2]]
Out[10]= 1
```

## Implementation notes

**Algorithm.** `builtin_cycle_graph` builds the undirected cycle `C_n` on the
integer vertices `1..n`. It emits the path edges `1 <-> 2, 2 <-> 3, ...,
(n-1) <-> n`, then adds the wrap edge `n <-> 1` only when `n >= 3`: for `n <= 2`
that edge would duplicate one already present (`n = 2` is a single edge, `n <= 1`
is edgeless), and Mathilda graphs are simple. The assembled
`Graph[List verts, List edges]` is returned and the evaluator canonicalises and
validates it through `builtin_graph`. A non-integer or negative argument, or any
arity other than one, returns `NULL` and the call is left unevaluated.

**Data structures.** Vertices are `EXPR_INTEGER` nodes `1..n` from `int_vertices`;
the edges are `UndirectedEdge` nodes in a `calloc`'d `Expr*` array. `make_graph`
wraps both arrays into the two `List`s of a `Graph[...]`, moving ownership of the
array contents into the new nodes (the arrays themselves are then freed). No
adjacency or incidence structure is built here — that is the validated-graph memo's
job, populated lazily once `builtin_graph` accepts the tree.

**Complexity / limits.** `O(n)` vertices and `O(n)` edges in a single allocation
pass. The output is a simple undirected graph — the `n <= 2` wrap-edge guard is
exactly what keeps it so — in which every vertex has degree two for `n >= 3`. `res`
is freed by the evaluator on success.

- `Protected`. Like all the graph generators (`CompleteGraph`, `CycleGraph`,
  `PathGraph`, `StarGraph`, `RandomGraph`), it builds a canonical graph with
  vertices `1..n` and undirected edges via the `Graph` constructor path.
- Related generator `CompleteGraph[n]` — `K_n`, with all `n(n-1)/2` edges.
- Small cases: `CycleGraph[2]` is the single edge `1 <-> 2` (no parallel edges),
  `CycleGraph[1]` is one isolated vertex and `CycleGraph[0]` the null graph. A
  symbolic `n` is left unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [CompleteGraph](../../graphs/CompleteGraph/), [PathGraph](../../graphs/PathGraph/), [StarGraph](../../graphs/StarGraph/), [RandomGraph](../../graphs/RandomGraph/), [Graph](../../graphs/Graph/)

- Source: [`src/graph/generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`CycleGraph[n]` is the undirected cycle `C_n` on the integer vertices `1..n`,
closing the path `1-2-...-n` with the edge `n <-> 1`. For `n >= 3` it has exactly
`n` edges and every vertex has degree two; `n = 2` is a single edge (the wrap edge
would duplicate it) and `n <= 1` is edgeless. The result is an opaque `Graph`, read
through accessors such as `EdgeList`, `VertexDegree`, and `EdgeCount`.
