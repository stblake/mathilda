# GraphComplement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphComplement[g] gives the complement of g: the same vertices, with an edge wherever g has none. For a directed or mixed g the complement is directed (u->v present iff no edge of g leads from u to v). Edges are in row-major VertexList order; weights are dropped.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[GraphComplement[CycleGraph[5]]]
Out[1]= {1 <-> 3, 1 <-> 4, 2 <-> 4, 2 <-> 5, 3 <-> 5}

In[2]:= EdgeList[GraphComplement[Graph[{1->2,2->3}]]]
Out[2]= {1 -> 3, 2 -> 1, 3 -> 1, 3 -> 2}

In[3]:= EdgeList[GraphComplement[CompleteGraph[4]]]
Out[3]= {}

In[4]:= GraphComplement[{1}]
Out[4]= GraphComplement[{1}]
```

### Options (1)

```mathematica
In[5]:= InputForm[GraphComplement[Graph[{1,2,3},{1<->2},EdgeWeight->{4}]]]
Out[5]= Graph[{1, 2, 3}, {1 <-> 3, 2 <-> 3}]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Undirected `g`: every non-adjacent pair `i < j`, as an undirected edge.
- Directed or mixed `g`: every ordered pair with no edge usable from `i` to `j`,
  as a directed edge.
- Edges in row-major `VertexList` order; weights are dropped.
- Performance: the complement of a 1000-cycle (498500 new edges) is at parity
  with Mathematica 15, its time dominated by allocating and freeing edge
  expressions (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
