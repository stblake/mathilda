### Worked examples

```mathematica
In[1]:= EdgeList[CanonicalGraph[CompleteGraph[3]]]  (* relabelled onto 1..n, edges sorted *)
```

```mathematica
In[1]:= EdgeList[CanonicalGraph[Graph[{2 <-> 3, 3 <-> 4, 4 <-> 2}]]]  (* a triangle on {2,3,4} canonicalises to one on {1,2,3} *)
```

```mathematica
In[1]:= CanonicalGraph[Graph[{1 <-> 2, 2 <-> 3}]] === CanonicalGraph[Graph[{5 <-> 9, 9 <-> 7}]]  (* isomorphic graphs share a canonical form *)
```

### Notes

`CanonicalGraph[g]` returns a relabelling of `g` onto the vertices `1..n` that is
an isomorphism invariant: two graphs are isomorphic exactly when their canonical
graphs are identical. The last example is the defining property — two differently
labelled paths of length two collapse to the same canonical graph, so `===`
(verbatim equality) is `True`.

The canonical labelling comes from an individualization-refinement engine (the
nauty / Traces family) run on a coloured reduction of `g`, so directed graphs,
self-loops and multigraphs are all supported; edge weights and other properties
are ignored. Since the result is an opaque `Graph` object, the examples read it
back through `EdgeList`. The search is budgeted and leaves the call unevaluated
on exhaustion.
