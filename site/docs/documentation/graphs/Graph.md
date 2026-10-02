# Graph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Graph[v, e] represents a graph with vertices v and edges e. Graph[e] derives the vertices from the edge list. Edges are DirectedEdge[u,v] or UndirectedEdge[u,v]; u->v and u<->v are accepted as shorthand. Graph[e, opts] and Graph[v, e, opts] attach per-edge lists, matched to e by position: EdgeWeight -> {w1, ...} and EdgeCapacity -> {c1, ...}. Simple graphs only: no self-loops or parallel edges.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (3)

```mathematica
In[1]:= Graph[{1,2,3,4}, {1->2, 2->3, 3->4, 4->1}]
Out[1]= Graph[<4 vertices, 4 edges>]

In[2]:= InputForm[Graph[{1,2}, {1<->2}]]
Out[2]= Graph[{1, 2}, {1 <-> 2}]

In[3]:= InputForm[Graph[{1->2, 2->3, 3->1}]]
Out[3]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 3 -> 1}]
```

### Scope (9)

```mathematica
In[4]:= InputForm[Graph[{a,b,c}, {DirectedEdge[a,b], UndirectedEdge[b,c]}]]
Out[4]= Graph[{a, b, c}, {a -> b, b <-> c}]

In[5]:= InputForm[Graph[{1,2,3}, {1->2, 2->3}, EdgeWeight -> {5, 7}]]
Out[5]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5, 7}]

In[6]:= FullForm[Graph[{1,2},{1<->2}]]
Out[6]= Graph[List[1, 2], List[UndirectedEdge[1, 2]]]

In[7]:= InputForm[Graph[{a<->b, b<->c, a<->c}, EdgeWeight -> {1.5, 2, 1}]]
Out[7]= Graph[{a, b, c}, {a <-> b, b <-> c, a <-> c}, EdgeWeight -> {1.5, 2, 1}]

In[8]:= InputForm[Graph[{1->2, 2->3}, EdgeCapacity -> {4, 5}, EdgeWeight -> {1, 2}]]
Out[8]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {1, 2}, EdgeCapacity -> {4, 5}]

In[9]:= Graph[{1,2}, {1->1}]
Out[9]= Graph[{1, 2}, {1 -> 1}]

In[10]:= Graph[{1,2}, {1->2, 1->2}]
Out[10]= Graph[{1, 2}, {1 -> 2, 1 -> 2}]

In[11]:= Graph[{1,2}, {1->3}]
Out[11]= Graph[{1, 2}, {1 -> 3}]

In[12]:= Graph[{1,2,3}, {1->2, 2->3}, EdgeWeight -> {5}]
Out[12]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5}]
```

## Options & behaviour

### Scope

**Malformed input** is returned unevaluated (a self-loop, a duplicate edge, an
endpoint missing from the vertex list, a weight list of the wrong length):

## Algorithm

construct.c - builtin_graph: normalize, derive, validate, canonicalize.

Accepts:

```text
  Graph[edges]                        -- vertices derived from the edges (directed default)
  Graph[verts, edges]                 -- explicit vertex list
  Graph[edges, opt...]
  Graph[verts, edges, opt...]         -- either form followed by per-edge options:
                                           EdgeWeight   -> {w1, ..., wm}
                                           EdgeCapacity -> {c1, ..., cm}
                                         each matched to `edges` by position. A list of
                                         the wrong length, a repeated option, or any
                                         other option leaves Graph[...] unevaluated,
                                         like every other rejection below.
```

The two forms are told apart by the second argument: a List is the edge list of the explicit-vertex form, a Rule the first option of the edges-only form. The canonical result stores the options after the edge list in their canonical order (EdgeWeight, then EdgeCapacity; see graph.h), whatever order they were given in, so every spelling of one graph canonicalizes to the same tree.

Edge sugar is normalized on construction:

```text
  Rule[u,v]        / u -> v    ->  DirectedEdge[u, v]
  TwoWayRule[u,v]  / u <-> v   ->  UndirectedEdge[u, v]
  DirectedEdge[u,v] / UndirectedEdge[u,v]   pass through unchanged
```

The result is the canonical Graph[List[verts], List[edges]] with vertices in first-appearance order (when derived). Malformed input -- 3-arg edges, self-loops, parallel edges, or an edge endpoint absent from an explicit vertex list -- leaves Graph[...] unevaluated (returns NULL).

Memory (SPEC section 4): the canonical tree is built entirely from expr_copy of the argument's parts, so `res` is never cannibalized. On success the evaluator frees `res`; on NULL it retains it. The "already canonical" case returns NULL so evaluation reaches a fixed point.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Graph construction, 40000 edges | 5.81 s | 3.98 s | 14.7 s |
| VertexDegree, 20000 vertices | 0 s | 0.045 s | 0.981 s |
| EdgeCount | 0 s | 0 s | 1.07 s |
| ConnectedComponents | 0 s | 0.001 s | 6.39 s |
| GraphDistance from vertex 1 | 0 s | 1.05e+04 s | 3.91 s |
| FindShortestPath 1 to 10000 | 0 s | 0.692 s | 0.057 s |

## Implementation notes

- `Protected`. A graph is a value: the constructor normalizes and validates its
  input and returns the canonical `Graph[List[verts], List[edges]]`, followed
  by `EdgeWeight -> List[weights]` and then `EdgeCapacity -> List[caps]` when
  given (always in that order, so every spelling of one graph is the same
  expression).
- Edge normalization: `u -> v` (`Rule`) and `DirectedEdge[u, v]` become
  `DirectedEdge[u, v]`; `u <-> v` (`TwoWayRule`) and `UndirectedEdge[u, v]`
  become `UndirectedEdge[u, v]`. Directed and undirected edges may be mixed.
- Malformed input is left unevaluated: self-loops, parallel/duplicate edges,
  3-argument edges, an edge endpoint absent from an explicit vertex list, an
  `EdgeWeight`/`EdgeCapacity` list whose length doesn't match the edge list, a
  repeated option, or any other option. Anti-parallel directed edges `u -> v` and `v -> u` are distinct and
  allowed.
- The two forms are told apart by the second argument: a `List` is the edge list
  of `Graph[v, e, ...]`, a rule the first option of `Graph[e, ...]`. Read the
  weights back with `EdgeWeight`; every weight-aware head (`GraphDistance`,
  `FindShortestPath`, `FindSpanningTree`, the cut family, ...) sees them the
  same way whichever form built the graph.
- Graph-editing heads (`EdgeDelete`, `Subgraph`, ...) carry `EdgeWeight` over;
  they do not yet carry `EdgeCapacity`.
- Printing: in standard output a graph shows a terse summary,
  `Graph[<n vertices, m edges>]`. `InputForm` and `FullForm` print the literal
  constructor, which round-trips through the parser.
- Validation is memoized per graph node (see the *Performance model* section of
  this file's preamble), so repeated queries on the same graph do not re-check
  it. `GraphQ` tests validity.

**Attributes:** `Protected`.

## References

**See also:** [FindMaximumFlow](../../graphs/FindMaximumFlow/), [Rule](../../assignment-and-rules/Rule/), [EdgeWeight](../../graphs/EdgeWeight/), [List](../../other-advanced/List/), [GraphDistance](../../graphs/GraphDistance/), [FindShortestPath](../../graphs/FindShortestPath/), [FindSpanningTree](../../graphs/FindSpanningTree/), [EdgeDelete](../../graphs/EdgeDelete/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
