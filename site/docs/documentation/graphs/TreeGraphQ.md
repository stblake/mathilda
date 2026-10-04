# TreeGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TreeGraphQ[g] gives True if g is a tree: at least one vertex, connected, and with no cycles, ignoring edge direction. A disconnected forest is not a tree.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= TreeGraphQ[Graph[{1->2,1->3}]]
Out[1]= True

In[2]:= TreeGraphQ[Graph[{1,2,3,4},{1<->2,3<->4}]]
Out[2]= False

In[3]:= TreeGraphQ[Graph[{1->2,2->1}]]
Out[3]= False

In[4]:= TreeGraphQ[Graph[{},{}]]
Out[4]= False

In[5]:= TreeGraphQ[Graph[{1},{}]]
Out[5]= True
```

### Applications (6)

A star is a tree

```mathematica
In[6]:= TreeGraphQ[StarGraph[5]]
Out[6]= True
```

A cycle has as many edges as vertices

```mathematica
In[7]:= TreeGraphQ[CycleGraph[4]]
Out[7]= False
```

A path is a tree

```mathematica
In[8]:= TreeGraphQ[PathGraph[{1, 2, 3, 4}]]
Out[8]= True
```

Direction is ignored, so an out-tree counts

```mathematica
In[9]:= TreeGraphQ[Graph[{1 -> 2, 1 -> 3}]]
Out[9]= True
```

Two components, so too few edges to be a tree

```mathematica
In[10]:= TreeGraphQ[Graph[{1, 2, 3, 4}, {1 <-> 2, 3 <-> 4}]]
Out[10]= False
```

Anti-parallel pair is two edges on two vertices

```mathematica
In[11]:= TreeGraphQ[Graph[{1 -> 2, 2 -> 1}]]
Out[11]= False
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

**Algorithm.** `builtin_tree_graph_q` returns `False` for a non-graph, then reads the cached `GRAPH_PROP_TREE` flag. On a miss it tests that the graph has at least one vertex and exactly `V - 1` edges, then runs union-find over all edges. It answers `True` when every edge joins two different components. With `V - 1` joins the graph is one component, so it is connected and acyclic. Direction is ignored, so the out-tree `1 -> 2, 1 -> 3` is a tree. The anti-parallel pair `1 -> 2, 2 -> 1` is two edges on two vertices and is not. The null graph is not a tree.

**Data structures.** It reads the memoized `graph_edge_indices` views (`eu[k]`, `ev[k]`, `directed[k]`) over the `Graph[List, List]` expression. Union-find uses two `int` arrays (`parent` and `size`) with path halving and union by size. The answer is stored with `graph_prop_set`.

**Complexity / limits.** `O(V + E alpha(V))` on the first query, and `O(1)` on a repeat for the same graph node. The edge count test rejects most non-trees before any union-find work. The result is always `True` or `False`: an argument that is not a valid graph gives `False`, not an unevaluated form.

- `Protected`. Requires at least one vertex, connectivity, and exactly
  `VertexCount - 1` edges.
- An out-tree such as `1 -> 2, 1 -> 3` is a tree; a disconnected forest, the
  anti-parallel pair `1 -> 2, 2 -> 1`, and the null graph are not. A single
  vertex is a tree.
- `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/acyclic.c`](https://github.com/stblake/mathilda/blob/main/src/graph/acyclic.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`TreeGraphQ` ignores edge direction, so an out-tree counts, while the anti-parallel pair `1 -> 2, 2 -> 1` does not. The null graph is not a tree. The verdict is cached on the graph, so repeating the query on the same graph is `O(1)`.
