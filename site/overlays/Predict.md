### Worked examples

```mathematica
In[1]:= p = Predict[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]  (* fit a linear regression *)
```

```mathematica
In[1]:= {p[10.], p["Coefficients"], p["Method"]}  (* apply the model and read its properties *)
```

```mathematica
In[1]:= Predict[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]  (* a matrix whose last column is the response *)
```

```mathematica
In[1]:= knn = Predict[{1. -> 1., 2. -> 4., 3. -> 9., 4. -> 16.}, Method -> "NearestNeighbors"]  (* the training set is the model *)
```

```mathematica
In[1]:= knn[2.5]  (* the mean response of the three nearest rows *)
```

### Notes

Data is either a list of rules `{features -> value, …}` or a matrix whose last column is
the response — both give the same model, and which was used is not remembered. A
one-feature model accepts a bare scalar, so `p[10.]` works without wrapping it in a list.

`"LinearRegression"` is the default and the only parametric method; its coefficients are
`{intercept, c₁, …}`. A collinear feature set, or fewer observations than parameters,
has no unique fit and returns unevaluated rather than inventing one.

`"NearestNeighbors"` has no fitted parameters — the training set *is* the model, which is
why it costs nothing to fit and everything to apply. The prediction is the mean response
of the `k` nearest rows (`k = 3` by default, clamped to the number of rows), so for a
query of `2.5` the three nearest responses `4`, `9`, `1` average to `14/3`. It answers
`"NeighborCount"` and `"TrainingData"` as well as the shared properties.
