# PredictorFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PredictorFunction[method, coefficients, featureCount] is the fitted object Predict returns. Apply it to a feature vector to get a prediction, or to "Method", "Coefficients" or "FeatureCount" to read it. A one-feature model also accepts a bare scalar.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

Fit a predictor; it prints elided

```mathematica
In[1]:= p = Predict[{1 -> 2, 2 -> 4, 3 -> 6}]
Out[1]= PredictorFunction["LinearRegression", <>]
```

Apply it to a new input

```mathematica
In[2]:= p[10]
Out[2]= 20.0
```

Read its coefficients: {intercept, slope}

```mathematica
In[3]:= p["Coefficients"]
Out[3]= {0.0, 2.0}
```

And the method it was fitted with

```mathematica
In[4]:= p["Method"]
Out[4]= "LinearRegression"
```

## Implementation notes

**Definition.** `PredictorFunction[method, coefficients, featureCount]` is the
fitted predictor object that `Predict` and `LinearModelFit` return. Apply it to a
feature vector to get a prediction, or to the strings `"Method"`, `"Coefficients"`
or `"FeatureCount"` to read it back; a one-feature model also accepts a bare scalar.
It is `Protected`, prints **elided** (its parameters are derived), has no builtin of
its own, and its docstring is in `predict.c`.

**Representation.** A positional inert head built by `ml_make_model`
(`src/ml/predict.c`): a method-tag string (`"LinearRegression"`,
`"NearestNeighbors"`), a `parameters` slot whose shape the method chooses (a
coefficient `List` for linear regression, the training matrix for nearest
neighbours), the feature count, and an extra count. Applying the object routes
through `ml_model_apply`, which confirms the head is `PredictorFunction`, reads the
method tag, and evaluates the prediction; for linear regression the coefficient list
is `{intercept, w1, ...}` and the result is its dot product with the feature vector.

**Usage & limits.** Produced by fitting, not written by hand. `Predict` implements
only `"LinearRegression"` (the default); `LinearModelFit` returns the same object
without regression diagnostics (RSquared, standard errors, ANOVA are not
implemented). A collinear feature set has no unique fit, so `Predict` returns
unevaluated rather than inventing one.

**Attributes:** `Protected`.

## References

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)

## Notes & additional examples

### Notes

`PredictorFunction[method, coefficients, featureCount]` is what `Predict` and
`LinearModelFit` return. Apply it to a feature vector (or, for a one-feature model, a
bare scalar) to get a prediction, or to `"Method"`, `"Coefficients"` or
`"FeatureCount"` to read it back. It prints elided because its parameters are
derived. `Predict` implements only `"LinearRegression"`, whose coefficient list is
`{intercept, w1, ...}`; a collinear feature set has no unique fit and returns
unevaluated.
