### Worked examples

```mathematica
In[1]:= tr = {{100., 200.}, {102., 205.}, {104., 210.}, {106., 215.}, {108., 220.}}  (* deliberately off-centre training data *)
```

```mathematica
In[1]:= r = DimensionReduction[tr, 1]  (* a reusable reducer, not reduced data *)
```

```mathematica
In[1]:= r[{110., 225.}]  (* project a point that was NOT in the training set *)
```

```mathematica
In[1]:= r[{104., 210.}]  (* exactly at the training mean, which projects to 0 *)
```

```mathematica
In[1]:= r["ReducedDimension"]  (* query a stored property *)
```

### Notes

`DimensionReduction` is the reusable counterpart to `DimensionReduce`: the latter gives
you reduced *training data*, this gives you a *reducer* you can apply to data it was
never trained on.

New rows are centred on the **training** column means, not the incoming batch's. That is
the entire point of a reusable reducer — it is what makes projections comparable across
batches — and the alternative bug is hard to see, since centring a single new point
against itself gives all zeros, which looks correct on data near the origin. The
off-centre training set above is chosen to expose it.

The reducer applies to a single feature vector or a matrix of them, reproduces
`DimensionReduce` on its own training data, and answers `"Method"`, `"FeatureCount"` and
`"ReducedDimension"`. Only `"PrincipalComponentsAnalysis"` is available as a reducer so
far.
