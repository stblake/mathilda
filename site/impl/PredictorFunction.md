---
source: src/ml/predict.c
---
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
