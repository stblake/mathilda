### Worked examples

```mathematica
In[1]:= VertexList[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* a directed triangle; vertices in first-appearance order *)
```

```mathematica
In[1]:= EdgeList[Graph[{1 <-> 2, 2 <-> 3}]]  (* the <-> sugar normalises to UndirectedEdge *)
```

```mathematica
In[1]:= DirectedGraphQ[Graph[{1 -> 2, 2 -> 3}]]  (* arrow edges build a directed graph *)
```

```mathematica
In[1]:= EdgeList[Graph[{1, 2, 3}, {1 <-> 2}]]  (* an explicit vertex list keeps the isolated vertex 3 *)
```

```mathematica
In[1]:= EdgeWeight[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 7}]]  (* per-edge weights, matched to edges by position *)
```

### Notes

`Graph` is the canonicalising constructor. Its result is an opaque object that
prints as `Graph[<n vertices, m edges>]` and is read through the accessors —
`VertexList`, `EdgeList`, `EdgeWeight`, `VertexCount`, `DirectedGraphQ`, and the
rest. `Graph[edges]` derives the vertices from the edges in first-appearance order;
`Graph[{verts}, {edges}]` gives them explicitly, so isolated vertices survive.
`Rule`/`->` and `TwoWayRule`/`<->` edge sugar normalise to `DirectedEdge` /
`UndirectedEdge`, and the options are stored in a fixed canonical order, so every
spelling of one graph is the same tree.

Mathilda graphs are **simple**: a self-loop (`Graph[{1 -> 1}]`), a parallel edge, a
wrong-length `EdgeWeight` list, or an edge endpoint missing from an explicit vertex
list leaves `Graph[...]` unevaluated rather than building a multigraph.
