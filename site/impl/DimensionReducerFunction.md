---
source: src/ml/predict.c
---
**Definition.** `DimensionReducerFunction[method, {means, loadings...}, featureCount,
reducedDimension]` is the reusable projector that `DimensionReduction[data, k]` returns.
It is an **inert head** (no builtin, just a docstring and the `Protected` attribute),
produced by `DimensionReduction` and consumed by `predict.c` when applied.

**Representation.** A four-argument `EXPR_FUNCTION`: the `method`, a parameter block
holding the training-column `means` and the projection `loadings` (the retained
principal axes), the input `featureCount`, and `reducedDimension` — the target
dimension `k`. New rows are centred on the stored training means before projection,
which is what makes projections comparable across batches.

**Usage & limits.** `Protected`. Apply it to one feature vector, or to a matrix of
vectors, to project into the reduced space; it also answers `"Method"`,
`"FeatureCount"` and `"ReducedDimension"`. The difference from `DimensionReduce[data,
k]` — which returns the reduced *training* data — is that this object generalises to
data it was not trained on. It holds only the means and loadings, so it cannot
reconstruct the original rows.
