# Predict

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Predict[data] fits a predictor to data and returns a PredictorFunction, which can be stored and applied to new inputs. Data is either a list of rules {features -> value, ...} or a matrix whose last column is the response. Method -> "LinearRegression" is the only method implemented and is the default; any other declines rather than silently linear-regressing. The returned object also answers "Method", "Coefficients" and "FeatureCount". A collinear feature set has no unique fit and returns unevaluated.`**

## Examples (8)

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

### Applications (5)

Fit a linear regression

```mathematica
In[4]:= p = Predict[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]
Out[4]= PredictorFunction["LinearRegression", <>]
```

Apply the model and read its properties

```mathematica
In[5]:= {p[10.], p["Coefficients"], p["Method"]}
Out[5]= {21.0, {1.0, 2.0}, "LinearRegression"}
```

A matrix whose last column is the response

```mathematica
In[6]:= Predict[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]
Out[6]= {1.0, 2.0, 3.0}
```

The training set is the model

```mathematica
In[7]:= knn = Predict[{1. -> 1., 2. -> 4., 3. -> 9., 4. -> 16.}, Method -> "NearestNeighbors"]
Out[7]= PredictorFunction["NearestNeighbors", <>]
```

The mean response of the three nearest rows

```mathematica
In[8]:= knn[2.5]
Out[8]= 4.66667
```

## Options & behaviour

### On verification

the 1-neighbour predictor is cross-checked against the existing
independent `Nearest` builtin, on both sides of a midpoint, so neighbour *selection* is
validated against separate code. `Nearest` supports only the one-neighbour scalar form
here (its `k` form and point form decline), so that check covers selection at `k = 1`
and says nothing about the averaging, which is asserted directly instead.

## Implementation notes

**Algorithm.** `builtin_predict` reads its training data with `ml_read_training`, which
accepts either a list of rules `{features -> value, …}` or a plain matrix whose **last
column is the response** — both give the same model, and which was used is not
remembered.

For the default `"LinearRegression"`, `ml_ols` solves the normal equations
`(AᵀA) c = Aᵀy` with `A = [1 | x]` (the intercept column formed on the fly) by Gaussian
elimination with partial pivoting. Normal equations rather than LAPACK's `dgels`
deliberately: `AᵀA` is only `(dim+1) × (dim+1)`, and a zero pivot (tested scaled against
the largest diagonal, `> 1e-12`) says "these features are collinear" *directly*, so the
singular case — which is meant to decline rather than return one of infinitely many
fits — is legible. Fewer observations than parameters (`n < dim + 1`) also declines.

`"NearestNeighbors"` fits nothing: the training set *is* the model, stored as one matrix
with the response appended to each row so features and answers cannot drift apart. `k`
defaults to 3 (enough to average away one noisy response, few enough to stay local) and
is clamped to `n`. Application (`ml_model_apply`, via the evaluator's composite-head
dispatch) returns the mean response of the `k` nearest rows, found by a partial
insertion selection into a `k`-sized list rather than a full sort, ties keeping the
earlier row. A `"NeighborsNumber"` sub-option on a `"LinearRegression"` is refused, not
ignored.

**Data structures.** A `PredictorFunction[method, params, dim, extra]` whose `params`
is a coefficient `List` `[intercept, c₁, …, c_d]` for linear regression, or the joint
training matrix for nearest neighbours — a positional, method-tagged `EXPR_FUNCTION`,
not a new node type (see `src/ml/predict.h`). Named properties (`"Coefficients"`,
`"Method"`, `"FeatureCount"`, and for k-NN `"NeighborCount"` / `"TrainingData"`) are
read through the application.

**Complexity / limits.** OLS is `O(n · p² + p³)` with `p = dim + 1`; a nearest-neighbour
prediction is `O(n · dim)`. A one-feature model accepts a bare scalar as well as a
one-element list. A wrongly-shaped input leaves the application unevaluated.

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

- T. Hastie, R. Tibshirani and J. Friedman, *The Elements of Statistical Learning*, 2nd ed. (Springer, 2009), §3.2 (least squares).
- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)

## Notes & additional examples

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
