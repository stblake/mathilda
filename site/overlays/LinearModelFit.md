### Worked examples

```mathematica
In[1]:= m = LinearModelFit[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]  (* an intercept plus a slope *)
```

```mathematica
In[1]:= m[10.]  (* apply the fitted line *)
```

```mathematica
In[1]:= m["Coefficients"]  (* the same coefficients Predict finds *)
```

```mathematica
In[1]:= LinearModelFit[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]  (* multiple regressors, last column the response *)
```

### Notes

`LinearModelFit[data]` fits a linear model with an intercept and returns a
`PredictorFunction` carrying its coefficients — the very ones `Predict` finds, through
the same least-squares kernel. It is applied and queried exactly like any other
`PredictorFunction`.

Wolfram's `LinearModelFit` returns a `FittedModel` whose properties are regression
diagnostics — `RSquared`, standard errors, ANOVA. Those are a separate piece of work and
are deliberately **not** approximated here; only the coefficients are offered, under a
name a reader of regression code will look for. A singular system returns unevaluated.
