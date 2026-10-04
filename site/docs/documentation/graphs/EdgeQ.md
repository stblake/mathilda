# EdgeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeQ[g, e] gives True if e is an edge of the graph g. u->v means DirectedEdge[u,v] and u<->v UndirectedEdge[u,v]; an undirected edge matches in either orientation, and direction must agree.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EdgeQ[CycleGraph[3], 2<->1]
Out[1]= True

In[2]:= EdgeQ[CycleGraph[3], 1->2]
Out[2]= False

In[3]:= EdgeQ[Graph[{1->2}], DirectedEdge[1,2]]
Out[3]= True

In[4]:= EdgeQ[Graph[{1->2}], 2->1]
Out[4]= False

In[5]:= EdgeQ[5, 1->2]
Out[5]= False
```

### Applications (7)

An undirected edge present

```mathematica
In[6]:= EdgeQ[CycleGraph[4], 1 <-> 2]
Out[6]= True
```

Undirected edges match in either orientation

```mathematica
In[7]:= EdgeQ[CycleGraph[4], 2 <-> 1]
Out[7]= True
```

No such edge in the cycle

```mathematica
In[8]:= EdgeQ[CycleGraph[4], 1 <-> 3]
Out[8]= False
```

A directed edge in its own orientation

```mathematica
In[9]:= EdgeQ[Graph[{1 -> 2, 2 -> 3}], 1 -> 2]
Out[9]= True
```

Direction is never blurred

```mathematica
In[10]:= EdgeQ[Graph[{1 -> 2, 2 -> 3}], 2 -> 1]
Out[10]= False
```

An undirected query does not match a directed edge

```mathematica
In[11]:= EdgeQ[Graph[{1 -> 2}], 1 <-> 2]
Out[11]= False
```

A non-graph gives False

```mathematica
In[12]:= EdgeQ[5, 1 <-> 2]
Out[12]= False
```

## Algorithm

membership.c - VertexQ[g, v] and EdgeQ[g, e].

```text
  VertexQ[g, v]   True iff v is a vertex of g
  EdgeQ[g, e]     True iff e is an edge of g
```

Membership is structural (SameQ, via expr_eq), as in the Wolfram Language: a vertex 1 is not matched by 1.0.

EdgeQ accepts the same edge sugar the constructor does -- u -> v means DirectedEdge[u, v], u <-> v means UndirectedEdge[u, v]. An undirected query matches an UndirectedEdge in either orientation; a directed query matches only a DirectedEdge with the same ordered endpoints. Direction is never blurred: u -> v is not an edge of Graph[{u, v}, {u <-> v}].

Both give False for anything that is not a valid graph. Both are O(1) hash probes into the validated-graph memo (graph_util.c), which already holds the vertex index and edge-key set that validating g built.

Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res.

## Implementation notes

**Algorithm.** `builtin_edge_q` takes `EdgeQ[g, e]`. It returns `False` when `g` is not a valid graph, or when `e` is not a two-argument edge. `query_kind` accepts `DirectedEdge` and `Rule` (`u -> v`) as directed, and `UndirectedEdge` and `TwoWayRule` (`u <-> v`) as undirected. It then calls `graph_has_edge(g, u, v, directed)`. An undirected query matches an `UndirectedEdge` in either orientation. A directed query matches only a `DirectedEdge` with the same ordered endpoints, so `1 -> 2` is not an edge of a graph whose only edge is `1 <-> 2`.

**Data structures.** The graph is `Graph[List, List]`. Membership is structural (`expr_eq`), so vertex `1` does not match `1.0`. The probe is an `O(1)` hash lookup in the validated-graph memo, which already holds the vertex index and edge-key set that validating `g` built.

**Complexity / limits.** `O(1)` per query once `g` has been validated. The first query on a new graph pays the validation cost, which is `O(V + E)`. The result is always `True` or `False`, never unevaluated, except for a wrong argument count. An edge whose endpoints are not vertices of `g` is `False`.

- `Protected`. Accepts the constructor's sugar (`u -> v` is
  `DirectedEdge[u, v]`, `u <-> v` is `UndirectedEdge[u, v]`).
- An undirected edge matches in either orientation; direction must agree, so
  `u -> v` is not an edge of `Graph[{u, v}, {u <-> v}]`.
- `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/membership.c`](https://github.com/stblake/mathilda/blob/main/src/graph/membership.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeQ` accepts the edge sugar `u -> v` (directed) and `u <-> v` (undirected). An undirected query matches either orientation, but direction is never blurred between the two kinds. Membership is structural, so vertex `1` does not match `1.0`. Anything that is not a valid graph gives `False`.
