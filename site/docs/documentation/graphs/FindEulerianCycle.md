# FindEulerianCycle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindEulerianCycle[g] gives {c}, an Eulerian cycle c of g as a list of edges (Hierholzer's algorithm, linear time), or {} if there is none. Only FindEulerianCycle[g] and FindEulerianCycle[g, 1] are supported; mixed graphs are left unevaluated.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindEulerianCycle[CycleGraph[4]]
Out[1]= {{1 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}

In[2]:= FindEulerianCycle[Graph[{1->2,2->3,3->1}]]
Out[2]= {{1 -> 2, 2 -> 3, 3 -> 1}}

In[3]:= FindEulerianCycle[Graph[{1<->2,2<->3,3<->1,3<->4,4<->5,5<->3}]]
Out[3]= {{1 <-> 3, 3 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}

In[4]:= FindEulerianCycle[PathGraph[Range[3]]]
Out[4]= {}

In[5]:= FindEulerianCycle[Graph[{1},{}]]
Out[5]= {{}}

In[6]:= FindEulerianCycle[CycleGraph[3], All]
Out[6]= FindEulerianCycle[Graph[<3 vertices, 3 edges>], All]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Gives `{{}}` for an edgeless graph with a vertex.
- Hierholzer's algorithm, iterative, linear time: starts at the first vertex with
  an edge, takes edges in `EdgeList` order, reports undirected cycles in pop
  order and directed ones forwards; undirected edges are written in the
  direction walked. This reproduces Mathematica's cycle in most cases (not all).
- `n > 1` / `All` and mixed graphs are left unevaluated.
- The cycle/path finders reuse a small per-graph cache of incidence lists
  (holding a reference to the graph, like the memo), so repeating a query skips
  the CSR build. 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
