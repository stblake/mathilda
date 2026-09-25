# LineGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LineGraph[g] gives the line graph of g: one vertex per edge (numbered by EdgeList position); undirected edges sharing an endpoint are joined, and for directed g, edge i -> edge j when i ends where j starts. Mixed graphs are left unevaluated.`**

## Examples (4)

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

## Implementation notes

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

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
