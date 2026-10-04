---
references:
  - "T. Hastie, R. Tibshirani and J. Friedman, *The Elements of Statistical Learning*, 2nd ed. (Springer, 2009), §3.2 (least squares)."
source: src/ml/predict.c
---
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
