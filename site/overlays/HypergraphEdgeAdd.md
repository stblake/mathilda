### Worked examples

```mathematica
In[1]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1, 2}}], {2, 3}]]  (* one hyperedge; vertex 3 is new *)
```

```mathematica
In[1]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1, 2}}], {{2, 3}, {3, 4}}]]  (* several hyperedges *)
```

```mathematica
In[1]:= InputForm[HypergraphEdgeAdd[Hypergraph[{{1, 2}}], {1, 2}]]  (* repeats are allowed: a multi-hypergraph *)
```

### Notes

`HypergraphEdgeAdd[h, e]` appends the hyperedge `e`, adding any new vertices in
first-appearance order; `HypergraphEdgeAdd[h, {e1, ...}]` appends several. Repeats
are permitted, so the result can be a multi-hypergraph.

The disambiguation rule for the Edge heads: a `List` second argument is a *list of
hyperedges* when every one of its elements is itself a List, otherwise it is a
single hyperedge. So `{2, 3}` adds one hyperedge, while `{{2, 3}, {3, 4}}` adds
two. The edit builds a fresh object. To remove hyperedges, use
`HypergraphEdgeDelete`.
