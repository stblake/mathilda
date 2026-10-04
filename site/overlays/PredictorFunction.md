### Worked examples

```mathematica
In[1]:= p = Predict[{1 -> 2, 2 -> 4, 3 -> 6}]  (* fit a predictor; it prints elided *)
```

```mathematica
In[1]:= p[10]  (* apply it to a new input *)
```

```mathematica
In[1]:= p["Coefficients"]  (* read its coefficients: {intercept, slope} *)
```

```mathematica
In[1]:= p["Method"]  (* and the method it was fitted with *)
```

### Notes

`PredictorFunction[method, coefficients, featureCount]` is what `Predict` and
`LinearModelFit` return. Apply it to a feature vector (or, for a one-feature model, a
bare scalar) to get a prediction, or to `"Method"`, `"Coefficients"` or
`"FeatureCount"` to read it back. It prints elided because its parameters are
derived. `Predict` implements only `"LinearRegression"`, whose coefficient list is
`{intercept, w1, ...}`; a collinear feature set has no unique fit and returns
unevaluated.
