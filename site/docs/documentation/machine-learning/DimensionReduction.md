# DimensionReduction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DimensionReduction[data, k] returns a DimensionReducerFunction projecting into k dimensions, applicable to data it was NOT trained on -- the difference from DimensionReduce[data, k], which returns the reduced training data. New rows are centred on the TRAINING column means, which is what makes projections comparable across batches. Accepts one point or a matrix of points, and answers "Method", "FeatureCount" and "ReducedDimension".`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

A point NOT in the training set

```mathematica
In[1]:= r[{110., 225.}]
Out[1]= r[{110.0, 225.0}]
```

Exactly at the training mean

```mathematica
In[2]:= r[{104., 210.}]
Out[2]= r[{104.0, 210.0}]
```

## Options & behaviour

A reducer is the *same* representation as a `PredictorFunction` — positional and
method-tagged — which is why adding it needed no new evaluation machinery.

## Implementation notes

- **New rows are centred on the training column means**, not on the incoming batch's.
  That is the entire point of a reusable reducer — it is what makes projections
  comparable across batches — and it is worth stating because the alternative bug is
  hard to see: centring a single new point against itself gives all zeros, which looks
  correct on data that happens to sit near the origin.
- Applies to a single feature vector or to a matrix of them, so a batch needs no `Map`.
- On its own training data it reproduces `DimensionReduce` exactly.
- Answers `"Method"`, `"FeatureCount"` and `"ReducedDimension"`.
- Only `"PrincipalComponentsAnalysis"` is available as a reducer so far; MDS has no
  out-of-sample extension without an explicit one (Nyström), and LSA's would need the
  term-document vocabulary carried along.

**Attributes:** `Protected`.

## References

**See also:** [DimensionReducerFunction](../../other-advanced/DimensionReducerFunction/), [DimensionReduce](../../machine-learning/DimensionReduce/), [Map](../../data-structures/Map/), [PredictorFunction](../../other-advanced/PredictorFunction/)

- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)
