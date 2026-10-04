# EulerianGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EulerianGraphQ[g] gives True if g has a cycle using every edge exactly once: all degrees even (undirected) or in-degree = out-degree (directed), with all edges in one connected component. Left unevaluated for mixed graphs.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= {EulerianGraphQ[CycleGraph[4]], EulerianGraphQ[PathGraph[Range[3]]], EulerianGraphQ[Graph[{1->2,2->3,3->1}]], EulerianGraphQ[Graph[{1},{}]], EulerianGraphQ[Graph[{},{}]]}
Out[1]= {True, False, True, True, False}

In[2]:= {EulerianGraphQ[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]], EulerianGraphQ[Graph[{1,2,3,4},{1<->2,2<->3,3<->1}]], EulerianGraphQ[x]}
Out[2]= {False, True, False}

In[3]:= EulerianGraphQ[Graph[{1->2,2<->3}]]
Out[3]= EulerianGraphQ[Graph[<3 vertices, 2 edges>]]
```

### Applications (4)

Every degree is 2, so a cycle is Eulerian

```mathematica
In[4]:= EulerianGraphQ[CycleGraph[5]]
Out[4]= True
```

The two ends have odd degree

```mathematica
In[5]:= EulerianGraphQ[PathGraph[{1, 2, 3}]]
Out[5]= False
```

Every vertex has even degree 4

```mathematica
In[6]:= EulerianGraphQ[CompleteGraph[5]]
Out[6]= True
```

Degree 3 is odd, so not Eulerian

```mathematica
In[7]:= EulerianGraphQ[CompleteGraph[4]]
Out[7]= False
```

## Implementation notes

**Algorithm.** `builtin_eulerian_graph_q` (via `gops_eulerian`) is `True` when
`g` has a closed walk using every edge exactly once. The degree test is the
classical one: every vertex of even degree for an undirected graph, or
`in-degree = out-degree` for a directed graph, **and** all edges lying in one
connected component. One pass over the edges accumulates a per-vertex balance
(parity XOR for undirected, `+1`/`-1` for directed) and marks the touched
vertices; a union-find then checks that the non-isolated vertices form a single
component. For a directed graph, balanced plus weakly connected implies strongly
connected, so weak connectivity suffices. An edgeless graph with at least one
vertex is Eulerian; the null graph is not; a mixed graph is left unevaluated.

**Data structures.** A `GopsView` over the graph plus three `O(V)` scratch
arrays: the balance counters `bal[]`, a `touched[]` bitmap, and the union-find
`parent[]`. All work is on the view's integer endpoint arrays.

**Complexity / limits.** `O(V + E)`. A non-graph argument gives `False`; a mixed
graph returns `NULL` (unevaluated), as Mathematica leaves it. `FindEulerianCycle`
constructs an actual tour (Hierholzer) when this predicate holds.

- `Protected`. A non-graph argument gives `False`.
- Tests that all degrees are even (undirected) or in-degree equals out-degree
  (directed), and that all edges lie in one connected component (isolated
  vertices are allowed).
- An edgeless graph with at least one vertex is Eulerian; the null graph is not.
- Mixed graphs are left unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A graph is Eulerian when it has a closed walk traversing every edge once. The
condition is: all edges in one component, and every vertex of even degree
(undirected) or with equal in- and out-degree (directed). `K_n` is Eulerian
exactly when `n` is odd.

Mixed graphs (both directed and undirected edges) are left unevaluated. Use
`FindEulerianCycle` to obtain an actual tour.
