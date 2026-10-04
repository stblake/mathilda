# LineGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LineGraph[g] gives the line graph of g: one vertex per edge (numbered by EdgeList position); undirected edges sharing an endpoint are joined, and for directed g, edge i -> edge j when i ends where j starts. Mixed graphs are left unevaluated.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[LineGraph[PathGraph[Range[4]]]]
Out[1]= {2 <-> 1, 3 <-> 2}

In[2]:= EdgeList[LineGraph[StarGraph[4]]]
Out[2]= {2 <-> 1, 3 <-> 1, 3 <-> 2}

In[3]:= EdgeList[LineGraph[Graph[{1->2,2->3,3->1}]]]
Out[3]= {1 -> 2, 2 -> 3, 3 -> 1}

In[4]:= LineGraph[Graph[{1->2,2<->3}]]
Out[4]= LineGraph[Graph[<3 vertices, 2 edges>]]
```

### Applications (2)

```mathematica
In[5]:= lg = LineGraph[PathGraph[4]];
```

Three edges become three vertices along a path

```mathematica
In[6]:= {VertexList[lg], EdgeList[lg]}
Out[6]= {{1, 2, 3}, {2 <-> 1, 3 <-> 2}}
```

## Implementation notes

**Algorithm.** `builtin_line_graph` builds the line graph: one vertex per edge of the input,
labelled by its `EdgeList` position `1..m` (rather than by the original edge expression). For an
undirected graph two line-graph vertices are joined when the corresponding edges share an
endpoint: scanning the incidence of each edge's endpoints, an undirected edge `j <-> i` is
emitted once, at the later edge. For a fully directed graph, edge `i -> edge j` whenever edge
`i` ends at the vertex where edge `j` starts. A mixed directed/undirected graph is left
unevaluated.

**Data structures.** A CSR incidence structure (`gops_inc_build`, modes `GOPS_INC_ALL` for
undirected and `GOPS_INC_OUT` for directed); the result is assembled as a canonical `Graph` on
integer vertices.

**Complexity / limits.** `O(E + sum of deg^2)` — the per-endpoint neighbour scans dominate, so a
high-degree vertex of degree `d` contributes about `d^2` line-graph edges. The result is a
canonical `Graph`; a non-graph or mixed argument returns unevaluated.

- `Protected`. A non-graph argument is left unevaluated.
- Vertices are `1..m`, the `EdgeList` positions.
- Undirected `g`: `j <-> i` for `i < j` sharing an endpoint, listed by `j`, then
  by shared endpoint, then by `i`.
- Directed `g`: `i -> j` when edge `i` ends where edge `j` starts, in `(i, j)` order.
- Mixed graphs are left unevaluated.
- Deviation: Mathematica numbers directed line-graph vertices in a traversal
  order of its own; Mathilda always uses `EdgeList` position (an isomorphic graph).

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gops_transform.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_transform.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

The line graph has one vertex per edge of the original graph, labelled by its `EdgeList`
position `1, 2, ..., m`. Two of these vertices are adjacent when the corresponding edges share an
endpoint (for a directed graph, when the head of one edge is the tail of the next).

So the line graph of a path `P_n` is a path `P_{n-1}`. A mixed graph (both directed and
undirected edges) is left unevaluated. The result is a canonical `Graph`; inspect it with
`VertexList` / `EdgeList`.
