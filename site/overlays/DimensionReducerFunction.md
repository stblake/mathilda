### Worked examples

```mathematica
In[1]:= Head[DimensionReduction[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}, 1]]  (* returns a reusable reducer *)
```

```mathematica
In[1]:= r = DimensionReduction[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}, 1]; r["ReducedDimension"]  (* query a property *)
```

### Notes

`DimensionReducerFunction[method, {means, loadings...}, featureCount, reducedDimension]`
is the projector that `DimensionReduction[data, k]` returns — a reusable object, not a
user-written head. Apply it to a feature vector, or to a matrix of vectors, to project
into the reduced `k`-dimensional space; it also answers `"Method"`, `"FeatureCount"` and
`"ReducedDimension"`.

New rows are centred on the stored training-column means before projection, which is
what makes projections comparable across batches. That is the difference from
`DimensionReduce[data, k]`, which returns the reduced *training* data: the
`DimensionReducerFunction` generalises to data it was not trained on. It is `Protected`.
