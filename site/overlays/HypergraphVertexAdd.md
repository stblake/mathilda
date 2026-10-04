### Worked examples

```mathematica
In[1]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1, 2}}], 3]]  (* a single new vertex *)
```

```mathematica
In[1]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1, 2}}], {3, 4}]]  (* several at once *)
```

```mathematica
In[1]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1, 2}}], 2]]  (* an existing vertex is left alone *)
```

### Notes

`HypergraphVertexAdd[h, v]` appends the vertex `v`, and `HypergraphVertexAdd[h, {v1, ...}]`
appends several; the hyperedges are untouched, so this introduces isolated
vertices. A vertex already present is silently skipped, and the new vertices keep
their given order after the existing ones.

In the Vertex heads any `List` second argument is read as a *list of vertices*, so
to add a single List-valued vertex wrap it as `{{...}}`. Because a memoized
hypergraph is held immutably, the edit builds a fresh object and never changes the
original. To add hyperedges instead, use `HypergraphEdgeAdd`.
