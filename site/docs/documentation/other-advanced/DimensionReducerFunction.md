# DimensionReducerFunction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DimensionReducerFunction[method, {means, loadings...}, featureCount, reducedDimension] is the reusable reducer DimensionReduction returns. Apply it to a feature vector, or to a matrix of them, to project into the reduced space.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Returns a reusable reducer

```mathematica
In[1]:= Head[DimensionReduction[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}, 1]]
Out[1]= DimensionReducerFunction
```

Query a property

```mathematica
In[2]:= r = DimensionReduction[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}, 1]; r["ReducedDimension"]
Out[2]= 1
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)

## Notes & additional examples

### Notes

`DimensionReducerFunction[method, {means, loadings...}, featureCount, reducedDimension]`
is the projector that `DimensionReduction[data, k]` returns — a reusable object, not a
user-written head. Apply it to a feature vector, or to a matrix of vectors, to project
into the reduced `k`-dimensional space; it also answers `"Method"`, `"FeatureCount"` and
`"ReducedDimension"`.

New rows are centred on the stored training-column means before projection, which is
what makes projections comparable across batches. That is the difference from
`DimensionReduce[data, k]`, which returns the reduced *training* data: the
`DimensionReducerFunction` generalises to data it was not trained on. It is `Protected`.
