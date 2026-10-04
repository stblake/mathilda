# GraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphQ[g] gives True if g is a valid graph, and False otherwise.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphQ[Graph[{1,2}, {1->2}]]
Out[1]= True

In[2]:= GraphQ[Graph[{1}, {1->1}]]
Out[2]= False

In[3]:= GraphQ[5]
Out[3]= False

In[4]:= GraphQ[Graph[{}, {}]]
Out[4]= True

In[5]:= GraphQ[CycleGraph[4]]
Out[5]= True
```

### Applications (6)

A generated graph is valid

```mathematica
In[6]:= GraphQ[CycleGraph[4]]
Out[6]= True
```

An atom is not a graph

```mathematica
In[7]:= GraphQ[5]
Out[7]= False
```

Explicit vertices and edges

```mathematica
In[8]:= GraphQ[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]]
Out[8]= True
```

An edge endpoint missing from the vertex list

```mathematica
In[9]:= GraphQ[Graph[{1, 2}, {1 <-> 3}]]
Out[9]= False
```

A bare list of edges is not yet a graph

```mathematica
In[10]:= GraphQ[{1 <-> 2, 2 <-> 3}]
Out[10]= False
```

Directed graphs are valid too

```mathematica
In[11]:= GraphQ[Graph[{1 -> 2, 2 -> 3}]]
Out[11]= True
```

## Algorithm

graphq.c - GraphQ[g]: is g a valid graph?

A thin wrapper over graph_is_valid (graph_util.c): returns the symbol True when the (already-evaluated) argument is a canonical, valid graph, and False otherwise. Non-unary calls are left unevaluated (NULL).

Memory (SPEC section 4): returns a freshly-allocated symbol; the evaluator frees `res`.

## Implementation notes

**Algorithm.** `builtin_graph_q` is a thin wrapper over `graph_is_valid`. It returns `True` for a canonical, valid `Graph[List[vertices], List[edges]]` and `False` for any other expression, including atoms. A call with other than one argument is left unevaluated. Validity means that the edge endpoints are vertices of the graph, that each edge is a `DirectedEdge` or `UndirectedEdge`, and that there are no duplicate or self-loop edges.

**Data structures.** Validation lives in `graph_util.c`, which memoizes the verdict on the graph node together with a vertex hash index, the edge-key set, and the `eu`/`ev`/`directed` edge-index views. Every other graph predicate and algorithm reuses that entry, so `GraphQ` is the call that fills it.

**Complexity / limits.** The first call on a graph is `O(V + E)`. Repeats on the same graph node are an `O(1)` memo hit. It is purely structural and never evaluates or rewrites its argument. Hypergraphs and other graph-like objects are not `Graph` expressions, so `GraphQ` gives `False` for them.

- `Protected`. Never left unevaluated: any non-graph gives `False`.
- A graph is valid when it is the canonical `Graph[List, List]` with every edge
  a 2-argument `DirectedEdge`/`UndirectedEdge`, no self-loops, no parallel
  edges, and every endpoint present in the vertex list (the same conditions the
  `Graph` constructor enforces).
- The null graph `Graph[{}, {}]` is valid.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- Source: [`src/graph/graphq.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graphq.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`GraphQ` is a structural test: it checks the `Graph[vertices, edges]` shape, that every edge endpoint is a listed vertex, and that edges are well formed. It never raises a message, and any non-graph gives `False`. The first call on a graph fills the validation memo that the other graph heads reuse.
