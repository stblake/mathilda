# DirectedGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DirectedGraph[g] replaces each undirected edge u<->v of g by the pair u->v, v->u (weights duplicated). DirectedGraph[g, "Acyclic"] instead orients each undirected edge from the vertex earlier in VertexList to the later one, giving a DAG for undirected g.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[DirectedGraph[PathGraph[Range[3]]]]
Out[1]= {1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2}

In[2]:= EdgeList[DirectedGraph[Graph[{3<->1, 2<->3, 1<->2}], "Acyclic"]]
Out[2]= {3 -> 1, 3 -> 2, 1 -> 2}

In[3]:= EdgeList[DirectedGraph[Graph[{1->2, 3<->1}], "Acyclic"]]
Out[3]= {1 -> 2, 1 -> 3}

In[4]:= DirectedGraph[CycleGraph[3], "Random"]
Out[4]= DirectedGraph[Graph[<3 vertices, 3 edges>], "Random"]
```

### Options (1)

```mathematica
In[5]:= InputForm[DirectedGraph[Graph[{1,2},{1<->2},EdgeWeight->{7}]]]
Out[5]= Graph[{1, 2}, {1 -> 2, 2 -> 1}, EdgeWeight -> {7, 7}]
```

### Applications (2)

Each undirected edge becomes a pair of arcs

```mathematica
In[6]:= EdgeList[DirectedGraph[Graph[{1 <-> 2, 2 <-> 3}]]]
Out[6]= {1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2}
```

Oriented low-to-high: a DAG

```mathematica
In[7]:= EdgeList[DirectedGraph[Graph[{1 <-> 2, 2 <-> 3}], "Acyclic"]]
Out[7]= {1 -> 2, 2 -> 3}
```

## Implementation notes

**Algorithm.** `builtin_directed_graph` converts a graph to a directed one. By default each
undirected edge `u <-> v` becomes the two directed edges `u -> v` and `v -> u` (weight
duplicated), while existing directed edges are kept in place; a fully directed graph is returned
unchanged. `DirectedGraph[g, "Acyclic"]` instead *orients* each undirected edge from the vertex
earlier in `VertexList` to the one later, yielding a DAG when the input is undirected (the edges
are then counting-sorted by tail then head); a mixed graph keeps its existing directed edges and
orders the rest.

**Data structures.** A `GopsView` with integer endpoints and a directed-edge count `ndir`; the
default form allocates `ne + (ne - ndir)` result edges, the "Acyclic" form reuses the edge count.
The result is a canonical `Graph` built by `gops_graph_new`.

**Complexity / limits.** `O(V + E)` (plus a counting sort in the undirected "Acyclic" case). A
second argument other than the string `"Acyclic"`, or a non-graph first argument, returns
unevaluated.

- `Protected`. A non-graph argument is left unevaluated.
- Plain form: the two directed edges replace the undirected one in place; the
  weight is duplicated onto both.
- `"Acyclic"`: gives a DAG for undirected `g`, with edges sorted by (tail, head)
  `VertexList` position; a mixed `g` keeps its edge order.
- Other methods (`"Random"`, ...) are left unevaluated.
- Performance: 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_transform.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_transform.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

By default `DirectedGraph[g]` replaces each undirected edge `u <-> v` by the two arcs `u -> v`
and `v -> u` (weights duplicated), so the result has the same reachability as `g`.

`DirectedGraph[g, "Acyclic"]` instead orients each undirected edge from the earlier to the later
vertex in `VertexList`, producing a DAG. A graph that is already fully directed is returned
unchanged. The result is a canonical `Graph` — inspect it with `EdgeList`.
