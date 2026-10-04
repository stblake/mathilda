### Worked examples

```mathematica
In[1]:= TransversalHypergraph[{{1, 2}, {2, 3}}]  (* {2} hits both; so does the minimal {1, 3} *)
```

```mathematica
In[1]:= InputForm[TransversalHypergraph[Hypergraph[{{a, b}, {b, c}}]]]  (* returned as a Hypergraph on h's vertices *)
```

```mathematica
In[1]:= TransversalHypergraph[{{1, 2}, {3, 4}}]  (* disjoint hyperedges: pick one vertex from each *)
```

### Notes

A transversal (hitting set) is a vertex set meeting every hyperedge; it is minimal
when no proper subset still does. `TransversalHypergraph[h]` returns all minimal
transversals — the transversal hypergraph `Tr(h)` of Berge — as
`Hypergraph[VertexList[h], Tr]`, or, for a bare List of hyperedges, as the List
`Tr`. Each transversal lists its vertices in `VertexList` order, and the
transversals are sorted by size then lexicographically.

The algorithm is MMCS (Murakami & Uno, 2014), the practical state of the art, run
iteratively so a large output cannot overflow the C stack. Because `Tr(h)` can be
exponentially large, the search is correct-or-unevaluated: it counts nodes and
returns unevaluated — never a partial answer — past its budget, and
`TimeConstrained` can interrupt it. Degenerate inputs the FR function declines are
answered here: no hyperedges gives `{{}}`, an empty hyperedge gives `{}`.
