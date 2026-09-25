# DirectedGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DirectedGraph[g] replaces each undirected edge u<->v of g by the pair u->v, v->u (weights duplicated). DirectedGraph[g, "Acyclic"] instead orients each undirected edge from the vertex earlier in VertexList to the later one, giving a DAG for undirected g.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[DirectedGraph[PathGraph[Range[3]]]]
Out[1]= {1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2}

In[2]:= EdgeList[DirectedGraph[Graph[{3<->1, 2<->3, 1<->2}], "Acyclic"]]
Out[2]= {3 -> 1, 3 -> 2, 1 -> 2}

In[3]:= EdgeList[DirectedGraph[Graph[{1->2, 3<->1}], "Acyclic"]]
Out[3]= {1 -> 2, 1 -> 3}

In[4]:= DirectedGraph[CycleGraph[3], "Random"]
Out[4]= DirectedGraph[Graph[<3 vertices, 3 edges>], "Random"]
```

### Options (1)

```mathematica
In[5]:= InputForm[DirectedGraph[Graph[{1,2},{1<->2},EdgeWeight->{7}]]]
Out[5]= Graph[{1, 2}, {1 -> 2, 2 -> 1}, EdgeWeight -> {7, 7}]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Plain form: the two directed edges replace the undirected one in place; the
  weight is duplicated onto both.
- `"Acyclic"`: gives a DAG for undirected `g`, with edges sorted by (tail, head)
  `VertexList` position; a mixed `g` keeps its edge order.
- Other methods (`"Random"`, ...) are left unevaluated.
- Performance: 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
