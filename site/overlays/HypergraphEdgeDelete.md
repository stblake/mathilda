### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= InputForm[HypergraphEdgeDelete[h, {3, 4}]]  (* delete the hyperedge {3, 4} *)
```

```mathematica
In[1]:= InputForm[HypergraphEdgeDelete[h, {{3, 4}, {7}}]]  (* delete several hyperedges *)
```

```mathematica
In[1]:= InputForm[HypergraphEdgeDelete[Hypergraph[{{1, 2}, {1, 2}}], {1, 2}]]  (* removes all copies of a repeated hyperedge *)
```

### Notes

`HypergraphEdgeDelete[h, e]` removes every hyperedge identical (`SameQ`) to `e`;
the list form removes every hyperedge matching one of several. The vertices are
left in place, so deleting a hyperedge can leave its vertices isolated.

This is one of the order-sensitive heads: hyperedges are compared **as written**,
so `{2, 1}` does not delete `{1, 2}`. All copies of a repeated hyperedge are
removed together. A named hyperedge that does not occur in `h` leaves the call
unevaluated. The disambiguation rule matches `HypergraphEdgeAdd`: a List of Lists
is a list of hyperedges, anything else a single hyperedge.
