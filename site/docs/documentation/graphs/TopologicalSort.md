# TopologicalSort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TopologicalSort[g] gives the vertices of the directed acyclic graph g ordered so that u precedes v for every edge u->v; ties go to the vertex earlier in VertexList[g]. TopologicalSort[{v->w, ...}] uses the rules as the graph. Left unevaluated if g has a cycle or an undirected edge.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= TopologicalSort[{1->3,1->4,2->1,2->4,3->4,5->2,5->3}]
Out[1]= {5, 2, 1, 3, 4}

In[2]:= TopologicalSort[Graph[{a,b,c},{c->a}]]
Out[2]= {b, c, a}

In[3]:= TopologicalSort[Graph[{3,1,2},{}]]
Out[3]= {3, 1, 2}

In[4]:= TopologicalSort[Graph[{1->2,2->3,3->1}]]
Out[4]= TopologicalSort[Graph[<3 vertices, 3 edges>]]

In[5]:= TopologicalSort[CycleGraph[3]]
Out[5]= TopologicalSort[Graph[<3 vertices, 3 edges>]]
```

### Applications (5)

A DAG: u before v for every edge u -> v

```mathematica
In[6]:= TopologicalSort[Graph[{1 -> 2, 2 -> 3, 1 -> 3}]]
Out[6]= {1, 2, 3}
```

The rule form builds Graph[...] first, then sorts

```mathematica
In[7]:= TopologicalSort[{4 -> 2, 2 -> 1, 4 -> 3}]
Out[7]= {4, 2, 1, 3}
```

An edgeless graph sorts to its VertexList

```mathematica
In[8]:= TopologicalSort[Graph[{1, 2, 3}, {}]]
Out[8]= {1, 2, 3}
```

A cycle has no topological order -- left unevaluated

```mathematica
In[9]:= TopologicalSort[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]
Out[9]= TopologicalSort[Graph[<3 vertices, 3 edges>]]
```

An undirected edge imposes no order, so this is not a DAG

```mathematica
In[10]:= TopologicalSort[CycleGraph[4]]
Out[10]= TopologicalSort[Graph[<4 vertices, 4 edges>]]
```

## Algorithm

acyclic.c - cycle structure and ordering.

```text
  AcyclicGraphQ[g]           True iff g has no cycle, respecting direction
  TreeGraphQ[g]              True iff g is a tree: connected, >= 1 vertex,
                             and exactly VertexCount - 1 edges
  TopologicalSort[g]         vertices of a directed acyclic graph with u
                             before v for every edge u -> v
  TopologicalSort[{rules}]   the same, for Graph[{rules}]
```

AcyclicGraphQ -- a cycle is a closed walk that uses no edge twice, following directed edges forwards only and undirected edges either way. So an undirected graph is acyclic iff it is a forest, a directed graph iff it is a DAG, and u -> v with v -> u is a (2-)cycle. Mixed graphs are decided exactly by contraction:

```text
  1. Union-find over the undirected edges. Joining two vertices already in
     one set closes an undirected cycle.
  2. Otherwise every undirected component is a tree, so any vertex in it can
     reach any other without reusing an edge. Collapse each to one node. A
     directed edge inside a component closes a cycle with the tree path back
     to its tail; the rest form a directed multigraph on the components.
  3. g is acyclic iff that contracted digraph is (Kahn's algorithm). A simple
     cycle there visits each component once, so it lifts to a cycle of g
     using one tree path per component; conversely any cycle of g that uses
     a directed edge projects to a closed walk in the contraction.
```

Linear time, O(V + E alpha(V)).

All three cache their answer on the graph's validated-graph memo entry (graph_prop_*, graph_cached_*), so repeating a query on the same graph node is O(1), as with Mathematica's atomic Graph object.

TreeGraphQ follows Wolfram's "a simple connected graph with no cycles" on the underlying undirected graph, where "connected with n - 1 edges" is the equivalent test (connectivity by union-find over the edges). Direction is ignored, so an out-tree like 1->2, 1->3 is a tree; the anti-parallel pair 1->2, 2->1 is two edges on two vertices, not a tree. The null graph (no vertices) is not a tree, matching ConnectedGraphQ.

TopologicalSort is Kahn's algorithm with ties broken by VertexList position (the smallest ready index goes first, via a binary min-heap), so the order is deterministic and an already-sorted vertex list is returned unchanged. It is left unevaluated for a graph with a cycle, or with any undirected edge (an undirected edge imposes no order, so the graph is not a DAG); an edgeless graph sorts to its VertexList.

AcyclicGraphQ/TreeGraphQ give False for a non-graph; TopologicalSort is left unevaluated.

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

**Algorithm.** `builtin_topological_sort` orders the vertices of a directed
acyclic graph so that `u` precedes `v` for every edge `u -> v`, by **Kahn's
algorithm**: repeatedly remove a vertex of in-degree 0 and decrement its
successors' in-degrees. Ties — several vertices ready at once — are broken by
**VertexList position** through a binary min-heap keyed on vertex index, so the
smallest ready index always goes first; the order is therefore deterministic and
an already-sorted vertex list comes back unchanged. The call is left
**unevaluated** (returns `NULL`) for a graph that has a cycle (Kahn's queue
drains before all vertices are emitted) or that carries any **undirected** edge
(an undirected edge imposes no order, so the graph is not a DAG); an edgeless
graph sorts to its `VertexList`. `TopologicalSort[{rules}]` builds `Graph[{rules}]`
through the ordinary constructor first — so the edge sugar and validation are
exactly `Graph`'s — then sorts that; a non-graph, non-list argument declines.

**Data structures.** The adjacency comes from `graph_build_adj` as a `GraphAdj`
(forward out-neighbour lists plus per-vertex in-degrees), read straight off the
validated-graph memo. Working state is an `indeg[]` counter array and the
integer min-heap; the result is assembled as a `List` of `expr_copy`'d vertex
expressions. That list is cached on the graph node's memo slot
(`GRAPH_CACHED_TOPOSORT`): a `List` is an immutable value, so handing out
references to the one cached copy is safe, and a repeated query on the same graph
node is free.

**Complexity / limits.** Linear in the graph, `O(V + E)` plus the `O(V log V)`
heap work the deterministic tie-break costs; `O(1)` on a memo hit. The output is a
plain vertex `List`, not a graph. There is no partial-order output for a cyclic or
mixed graph — the head declines rather than return a best effort, matching
Mathematica's leaving `TopologicalSort` unevaluated there.

- `Protected`. Kahn's algorithm; among the vertices ready at each step, the one
  earliest in `VertexList[g]` goes first, so the order is deterministic.
- Left unevaluated for a cyclic graph, a graph with any undirected edge, or a
  non-graph. An edgeless graph sorts to its `VertexList`.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/), [VertexList](../../graphs/VertexList/)

- A. B. Kahn, *Topological sorting of large networks*, Communications of the ACM **5** (1962) 558-562.
- Source: [`src/graph/acyclic.c`](https://github.com/stblake/mathilda/blob/main/src/graph/acyclic.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`TopologicalSort` is Kahn's algorithm: it repeatedly emits a vertex with no
remaining incoming edges. When several vertices are ready at once the tie is broken
by `VertexList` position (smallest index first), so the order is deterministic and
an already-sorted vertex list is returned unchanged.

The head applies only to a **directed acyclic graph**. A graph with a cycle, or one
carrying any undirected edge (an undirected edge imposes no `u` before `v`
ordering), is left unevaluated rather than given a best-effort order — so the
cyclic and undirected examples above print back as the input. A plain list of rules
is accepted as shorthand and built into a `Graph` with exactly the constructor's
edge sugar before sorting.
