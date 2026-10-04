# GraphComplement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphComplement[g] gives the complement of g: the same vertices, with an edge wherever g has none. For a directed or mixed g the complement is directed (u->v present iff no edge of g leads from u to v). Edges are in row-major VertexList order; weights are dropped.`**

## Examples (8)

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

### Applications (3)

The two missing diagonals

```mathematica
In[6]:= EdgeList[GraphComplement[CycleGraph[4]]]
Out[6]= {1 <-> 3, 2 <-> 4}
```

The complement of K_n is edgeless

```mathematica
In[7]:= EdgeCount[GraphComplement[CompleteGraph[5]]]
Out[7]= 0
```

Only the 1-3 pair is missing

```mathematica
In[8]:= EdgeCount[GraphComplement[PathGraph[{1, 2, 3}]]]
Out[8]= 1
```

## Implementation notes

**Algorithm.** `builtin_graph_complement` gives the graph on the same vertices
with an edge wherever `g` has none. For an undirected `g` every non-adjacent pair
`i < j` becomes `i <-> j`; for a directed or mixed `g` every ordered pair `(i,
j)`, `i != j`, with no edge usable from `i` to `j` becomes the directed edge `i
-> j`. For each source vertex `i` it stamps `i`'s existing out-neighbours into a
scratch array and then emits an edge to every unstamped `j`, in row-major
`VertexList` order. Weights are dropped.

**Data structures.** A `GopsView` of the graph, an out-incidence CSR
(`gops_inc_build`), a per-row `stamp[]` array marking existing neighbours, and a
growable `EdgeBuf` sized up front to the exact complement edge count. Result
edges share a single interned edge head (`GopsHeads`) and the finished graph is
seeded into the memo (`gops_graph_new`).

**Complexity / limits.** Output-bound: `O(V^2)` since the complement of a sparse
graph is dense. The `tc_check_deadline` poll every 64 rows keeps it interruptible
under `TimeConstrained`.

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

- Source: [`src/graph/gops_transform.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_transform.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`GraphComplement[g]` keeps the vertices and flips adjacency: an edge appears
exactly where `g` has none. The complement of `CompleteGraph[n]` has no edges,
and the complement of the empty graph on `n` vertices is `CompleteGraph[n]`.

For a directed or mixed graph the complement is directed: `i -> j` is present
precisely when no edge of `g` leads from `i` to `j`. Edge weights are dropped.
