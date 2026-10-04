### Worked examples

```mathematica
In[1]:= ConnectedHypergraphQ[Hypergraph[{{1, 2}, {2, 3}}]]  (* a single component *)
```

```mathematica
In[1]:= ConnectedHypergraphQ[Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]]  (* vertex 7 is isolated *)
```

```mathematica
In[1]:= ConnectedHypergraphQ[{{1, 2}, {2, 3}, {3, 4}}]  (* accepts a bare list of hyperedges *)
```

### Notes

`ConnectedHypergraphQ[h]` is `True` when `h` has at least one vertex and exactly
one connected component — the Boolean companion of
`HypergraphConnectedComponents`. It is the Wolfram Function Repository name and,
like that function, accepts a bare List of hyperedges as well as a `Hypergraph`
object.

It gives `False` for any non-hypergraph and for the empty List `{}` (which the FR
function leaves unevaluated). Internally it reuses the same vertex union–find and
simply checks for a single root, so it is near-linear in the total incidence.
