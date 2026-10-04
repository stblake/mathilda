### Worked examples

```mathematica
In[1]:= FindClusters[{1, 2, 3, 10, 11, 12, 20}]  (* the gaps choose three clusters *)
```

```mathematica
In[1]:= FindClusters[Join[Range[5], Range[20, 25]]]  (* two well-separated runs *)
```

```mathematica
In[1]:= FindClusters[{1.0, 1.1, 5.0, 5.1, 5.2}, 2]  (* ask for exactly two clusters *)
```

### Notes

`FindClusters[list]` partitions `list` into clusters of nearby elements, choosing
the number automatically from the gaps; `FindClusters[list, n]` forces exactly
`n`, and `FindClusters[list, UpTo[n]]` gives at most `n`. Clusters appear in order
of the first occurrence of any member, and elements keep their input order within
a cluster.

All elements must be of one kind: real numbers (distance on the line),
equal-length numeric vectors (squared Euclidean), colours whose arguments are
coordinates, or strings (`EditDistance`). One dimension has no size cap; above it
the partition is built from a minimum spanning tree and is capped (20000 machine
points, 2000 exact points or strings). The result is **not** intended to match
Mathematica's — it auto-selects an unpublished metric — so this implements the
textbook algorithm for each named `Method`.
