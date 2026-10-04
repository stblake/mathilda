### Worked examples

```mathematica
In[1]:= EdgeList[GraphComplement[CycleGraph[4]]]  (* the two missing diagonals *)
```

```mathematica
In[1]:= EdgeCount[GraphComplement[CompleteGraph[5]]]  (* the complement of K_n is edgeless *)
```

```mathematica
In[1]:= EdgeCount[GraphComplement[PathGraph[{1, 2, 3}]]]  (* only the 1-3 pair is missing *)
```

### Notes

`GraphComplement[g]` keeps the vertices and flips adjacency: an edge appears
exactly where `g` has none. The complement of `CompleteGraph[n]` has no edges,
and the complement of the empty graph on `n` vertices is `CompleteGraph[n]`.

For a directed or mixed graph the complement is directed: `i -> j` is present
precisely when no edge of `g` leads from `i` to `j`. Edge weights are dropped.
