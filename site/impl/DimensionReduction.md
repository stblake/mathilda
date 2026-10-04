---
references:
  - "I. T. Jolliffe, *Principal Component Analysis*, 2nd ed. (Springer, 2002)."
source: src/ml/predict.c
---
**Algorithm.** `builtin_dimension_reduction` is the reusable counterpart to
`DimensionReduce`: where the latter returns reduced *data*, this returns a reusable
*reducer*. It computes the training column means (`ml_column_mean`) and the full
eigenvector matrix (`ml_pca`, with `correlation = false`), then stores a
`DimensionReducerFunction["PrincipalComponentsAnalysis", {means, loading₁, …,
loading_k}, dim, k]` — the means as the first payload row and the leading `k` loadings
as the rest. Only PCA is offered as a reducer (MDS has no out-of-sample extension
without an explicit Nyström step; LSA would need the term-document vocabulary carried
along).

Applying the reducer (`ml_reducer_apply`, reached through the evaluator's composite-head
dispatch) projects a point with `ml_project_row`: it subtracts the **training** means —
not the incoming batch's — and contracts the shifted vector with each stored loading.
Using the training means is the entire point of a reusable reducer, and the subtle bug
it avoids is that centring a lone new point against itself gives all zeros, which looks
correct on data near the origin. A single feature vector or a whole matrix of them is
accepted, so a batch needs no `Map`. The object also answers the string queries
`"Method"`, `"FeatureCount"` and `"ReducedDimension"`.

**Data structures.** The payload is an ordinary `List`; a reducer is the *same*
positional, method-tagged representation as a `PredictorFunction`, which is why it needed
no new node type and no change to `eval.c` beyond `ml_model_apply_probe` recognising one
more head.

**Complexity / limits.** Fitting is a PCA: `O(n · dim² + dim³)`. Applying is
`O(target · dim)` per point. A target dimension larger than `dim`, or a flat-list
input, declines. On its own training data the reducer reproduces `DimensionReduce`
exactly.
