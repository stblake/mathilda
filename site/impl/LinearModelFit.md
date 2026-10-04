---
references:
  - "T. Hastie, R. Tibshirani and J. Friedman, *The Elements of Statistical Learning*, 2nd ed. (Springer, 2009), §3.2 (least squares)."
source: src/ml/predict.c
---
**Algorithm.** `builtin_linear_model_fit` fits an ordinary least-squares line with an
intercept through the *same* `ml_ols` kernel `Predict` uses — the normal equations
`(AᵀA) c = Aᵀy` with `A = [1 | x]`, solved by Gaussian elimination with partial
pivoting — and returns a `PredictorFunction["LinearRegression", {intercept, c₁, …}, dim,
0]`. Data is read by `ml_read_training`, so a list of rules `{features -> value, …}` or
a matrix whose last column is the response both work. A singular system (collinear
features, or fewer observations than parameters) declines, since there is no unique fit.

The object is applied and queried exactly like any other `PredictorFunction`:
`m[x]` evaluates `intercept + c · x`, and `m["Coefficients"]` reads the fitted vector —
the same coefficients `Predict` finds.

**Data structures.** The shared positional, method-tagged model representation (see
`src/ml/predict.h`); the coefficient `List` is the only payload.

**Complexity / limits.** `O(n · p² + p³)` with `p = dim + 1`. Wolfram's `LinearModelFit`
returns a `FittedModel` whose properties are regression diagnostics (`RSquared`,
standard errors, ANOVA); those are a separate piece of work and are **not** approximated
here — only the coefficients, under a name a reader of regression code will look for.
