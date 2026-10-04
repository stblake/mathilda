# FindEulerianCycle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindEulerianCycle[g] gives {c}, an Eulerian cycle c of g as a list of edges (Hierholzer's algorithm, linear time), or {} if there is none. Only FindEulerianCycle[g] and FindEulerianCycle[g, 1] are supported; mixed graphs are left unevaluated.`**

## Examples (11)

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

### Applications (5)

Every vertex has even degree, so a closed Euler tour exists

```mathematica
In[7]:= FindEulerianCycle[CycleGraph[4]]
Out[7]= {{1 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}
```

The triangle is Eulerian

```mathematica
In[8]:= FindEulerianCycle[CompleteGraph[3]]
Out[8]= {{1 <-> 3, 3 <-> 2, 2 <-> 1}}
```

Directed: in-degree equals out-degree at each vertex

```mathematica
In[9]:= FindEulerianCycle[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]
Out[9]= {{1 -> 2, 2 -> 3, 3 -> 1}}
```

Odd-degree endpoints, so no Euler cycle -- an empty result

```mathematica
In[10]:= FindEulerianCycle[PathGraph[{1, 2, 3}]]
Out[10]= {}
```

Every vertex has odd degree 3, so again none exists

```mathematica
In[11]:= FindEulerianCycle[CompleteGraph[4]]
Out[11]= {}
```

## Implementation notes

**Algorithm.** `builtin_find_eulerian_cycle` returns `{cycle}` — a closed walk
using every edge exactly once, as a list of edges — or `{}` when the graph has no
Eulerian cycle. It first calls `gops_eulerian`, which decides the precondition
from degrees and connectivity alone: for a directed graph every vertex must have
in-degree equal to out-degree, for an undirected graph every vertex must have even
degree, and the edge-touched vertices must form one connected block (balance plus
weak connectivity gives strong connectivity in the directed case). A non-Eulerian
graph gives `{}`; an edgeless graph with at least one vertex gives `{{}}`; the null
graph `{}`. When the graph is Eulerian the walk is built by **Hierholzer's
algorithm**, run iteratively: starting at the first vertex in `VertexList` order
that has an incident edge, it follows unused edges (taken in `EdgeList` order)
until it returns stuck, popping vertices to record the tour. Directed cycles are
reported forwards; undirected cycles in Hierholzer's pop order, each undirected
edge written in the direction it was actually walked, as Mathematica does. Only
`FindEulerianCycle[g]` and `FindEulerianCycle[g, 1]` are supported — mixed graphs,
and a second argument other than `1` (`n > 1` or `All`), leave the call unevaluated.

**Data structures.** The graph is read through a `GopsView`; traversal uses a
cached incidence CSR (`GOPS_INC_OUT` for a directed graph, `GOPS_INC_ALL` for an
undirected one). Scratch is a `used[]` edge-visited bitmap, a per-vertex incidence
cursor `ptr[]`, an explicit DFS stack (`sv`/`se` for vertex and entering edge), and
pop-order arrays (`cv`/`ce`) — no recursion, so there is no call-stack bound. The
result is `{List of edges}`, each edge `expr_copy`'d from the original (or
re-oriented for a reversed undirected traversal via `walked_edge`).

**Complexity / limits.** `O(V + E)` in both the Eulerian test and the walk, with
one incidence cursor advance per edge. A `TimeConstrained` abort can `longjmp` out
and leak the scratch arrays, the standard tradeoff across `src/graph/`. The output
is a concrete edge list, not a graph.

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

- C. Hierholzer, *Über die Möglichkeit, einen Linienzug ohne Wiederholung und ohne Unterbrechung zu umfahren*, Mathematische Annalen **6** (1873) 30-32.
- Source: [`src/graph/gops_cycles.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_cycles.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`FindEulerianCycle[g]` returns `{cycle}`, where the cycle is the list of edges of a
closed walk that uses every edge of `g` exactly once, or `{}` when no such walk
exists. A graph is Eulerian when it is connected and every vertex has even degree
(undirected) or equal in- and out-degree (directed); the head checks that first and
returns `{}` immediately otherwise.

The tour is built by Hierholzer's algorithm, starting from the first vertex (in
`VertexList` order) that has an incident edge and taking edges in `EdgeList` order,
so the walk is deterministic. Only `FindEulerianCycle[g]` and `FindEulerianCycle[g,
1]` are supported; a mixed directed/undirected graph, or a request for more than
one tour, is left unevaluated.
