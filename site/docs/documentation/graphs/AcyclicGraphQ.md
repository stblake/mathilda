# AcyclicGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AcyclicGraphQ[g] gives True if g has no cycle, following directed edges forwards and undirected edges either way: a forest when undirected, a DAG when directed. u->v together with v->u is a cycle.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= AcyclicGraphQ[Graph[{1,2,3},{1<->2,2->3,3->1}]]
Out[1]= False

In[2]:= AcyclicGraphQ[PathGraph[4]]
Out[2]= True

In[3]:= AcyclicGraphQ[CycleGraph[4]]
Out[3]= False

In[4]:= AcyclicGraphQ[Graph[{1->2,2->1}]]
Out[4]= False

In[5]:= AcyclicGraphQ[Graph[{1->2,1->3,2->4,3->4}]]
Out[5]= True
```

### Applications (4)

A directed chain is a DAG

```mathematica
In[6]:= AcyclicGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 4}]]
Out[6]= True
```

The 1-2-3 cycle makes it False

```mathematica
In[7]:= AcyclicGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]
Out[7]= False
```

An undirected tree is a forest

```mathematica
In[8]:= AcyclicGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}]]
Out[8]= True
```

Anti-parallel edges are a 2-cycle

```mathematica
In[9]:= AcyclicGraphQ[Graph[{1 -> 2, 2 -> 1}]]
Out[9]= False
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

**Algorithm.** `builtin_acyclic_graph_q` decides whether the graph has a cycle, treating
directed edges as forward-only and undirected edges as traversable either way, by an exact
three-step contraction. (1) A **union-find** (path halving, union by size) joins the endpoints
of every undirected edge; a join of two already-connected vertices closes an undirected cycle,
so the graph is cyclic. (2) Each directed edge is mapped onto the contracted components; a
directed edge inside one component also closes a cycle. (3) The contracted digraph is run
through **Kahn's topological sort** — the graph is acyclic iff every vertex is eventually
removed. So an undirected graph is acyclic exactly when it is a forest, and a directed graph
exactly when it is a DAG; an anti-parallel pair `u -> v, v -> u` is a 2-cycle. The 0/1 answer is
cached on the graph node (`GRAPH_PROP_ACYCLIC`).

**Data structures.** The union-find `parent[]`/`size[]`; pre-resolved endpoint index arrays
(`eu`/`ev`/`edir`) from the validated-graph memo; staged directed tails/heads and in/out degree
arrays; a CSR `start[]`/`succ[]` for the contracted digraph and a Kahn work queue.

**Complexity / limits.** `O(V + E·alpha(V))`, no cap. Since the `Graph` constructor rejects
self-loops, such input is never a valid graph and `AcyclicGraphQ` returns `False`; any non-graph
argument also returns `False` (it is a predicate), never unevaluated.

- `Protected`. A cycle follows directed edges forwards and undirected edges
  either way, never reusing an edge. An undirected graph is acyclic iff it is a
  forest, a directed graph iff it is a DAG; `u -> v` with `v -> u` is a 2-cycle.
- Mixed graphs are decided exactly: undirected components are contracted (a
  repeated union-find join is an undirected cycle, and a directed edge inside
  one component closes a cycle through it), then the contracted digraph is
  checked with Kahn's algorithm.
- `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- A. B. Kahn, *Topological sorting of large networks*, Comm. ACM **5** (1962) 558-562.
- Source: [`src/graph/acyclic.c`](https://github.com/stblake/mathilda/blob/main/src/graph/acyclic.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

A cycle is a closed walk using no edge twice, following directed edges forwards and undirected
edges either way. So an undirected graph is acyclic exactly when it is a forest, and a directed
graph exactly when it is a DAG. A pair `u -> v` together with `v -> u` counts as a cycle.

Being a predicate, it returns `False` (never stays unevaluated) for anything that is not an
acyclic graph, including a non-graph argument.
