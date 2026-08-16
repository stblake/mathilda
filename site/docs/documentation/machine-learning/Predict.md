# Predict

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Predict[data] fits a predictor to data and returns a PredictorFunction, which can be stored and applied to new inputs. Data is either a list of rules {features -> value, ...} or a matrix whose last column is the response. Method -> "LinearRegression" is the only method implemented and is the default; any other declines rather than silently linear-regressing. The returned object also answers "Method", "Coefficients" and "FeatureCount". A collinear feature set has no unique fit and returns unevaluated.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= p = Predict[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]
Out[1]= PredictorFunction["LinearRegression", <>]

In[2]:= {p[10.], p["Coefficients"], p["Method"]}
Out[2]= {21.0, {1.0, 2.0}, "LinearRegression"}

In[3]:= Predict[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]
Out[3]= {1.0, 2.0, 3.0}
```

## Options & behaviour

### On verification

the 1-neighbour predictor is cross-checked against the existing
independent `Nearest` builtin, on both sides of a midpoint, so neighbour *selection* is
validated against separate code. `Nearest` supports only the one-neighbour scalar form
here (its `k` form and point form decline), so that check covers selection at `k = 1`
and says nothing about the averaging, which is asserted directly instead.

## Implementation notes

- Data is either a list of rules `{features -> value, …}` or a **matrix whose last
  column is the response** — the shape data usually arrives in. Both give the same
  model; which was used is not remembered.
- A one-feature model accepts a bare scalar as well as a one-element list, so
  `p[3.]` works for a single-variable regression.
- **A singular system returns unevaluated.** Perfectly collinear features, or fewer
  observations than parameters, have no unique fit — a pseudo-inverse would return one
  of infinitely many answers and *look* like a successful fit.
- `"LinearRegression"` is the only method implemented and the default. Any other
  declines rather than silently linear-regressing.
- A wrongly-shaped input leaves the application unevaluated rather than guessing.
- **`"NearestNeighbors"` has no fitted parameters — the training set *is* the model.**
  That is why it costs nothing to fit and everything to apply. The prediction is the
  mean response of the `k` nearest training rows; `k` defaults to 3, enough to average
  away one noisy response and few enough to stay local, and is clamped to the number of
  training rows. It answers `"NeighborCount"` and `"TrainingData"` as well as the shared
  properties.
- `"NeighborsNumber"` on a `"LinearRegression"` is **refused, not ignored** — silently
  accepting a meaningless option would hide a real mistake. So is an unrecognised
  sub-option, and a non-positive `k`.

**Attributes:** `Protected`.

## References

**See also:** [PredictorFunction](../../other-advanced/PredictorFunction/), [Nearest](../../lists-and-iteration/Nearest/)

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)
