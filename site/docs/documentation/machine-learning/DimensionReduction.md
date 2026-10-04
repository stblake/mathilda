# DimensionReduction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DimensionReduction[data, k] returns a DimensionReducerFunction projecting into k dimensions, applicable to data it was NOT trained on -- the difference from DimensionReduce[data, k], which returns the reduced training data. New rows are centred on the TRAINING column means, which is what makes projections comparable across batches. Accepts one point or a matrix of points, and answers "Method", "FeatureCount" and "ReducedDimension".`**

## Examples (7)

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

### Applications (5)

Deliberately off-centre training data

```mathematica
In[3]:= tr = {{100., 200.}, {102., 205.}, {104., 210.}, {106., 215.}, {108., 220.}}
Out[3]= {{100.0, 200.0}, {102.0, 205.0}, {104.0, 210.0}, {106.0, 215.0}, {108.0, 220.0}}
```

A reusable reducer, not reduced data

```mathematica
In[4]:= r = DimensionReduction[tr, 1]
Out[4]= DimensionReducerFunction["PrincipalComponentsAnalysis", <>]
```

Project a point that was NOT in the training set

```mathematica
In[5]:= r[{110., 225.}]
Out[5]= {16.1555}
```

Exactly at the training mean, which projects to 0

```mathematica
In[6]:= r[{104., 210.}]
Out[6]= {0.0}
```

Query a stored property

```mathematica
In[7]:= r["ReducedDimension"]
Out[7]= 1
```

## Options & behaviour

A reducer is the *same* representation as a `PredictorFunction` — positional and
method-tagged — which is why adding it needed no new evaluation machinery.

## Implementation notes

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

- I. T. Jolliffe, *Principal Component Analysis*, 2nd ed. (Springer, 2002).
- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)

## Notes & additional examples

### Notes

`DimensionReduction` is the reusable counterpart to `DimensionReduce`: the latter gives
you reduced *training data*, this gives you a *reducer* you can apply to data it was
never trained on.

New rows are centred on the **training** column means, not the incoming batch's. That is
the entire point of a reusable reducer — it is what makes projections comparable across
batches — and the alternative bug is hard to see, since centring a single new point
against itself gives all zeros, which looks correct on data near the origin. The
off-centre training set above is chosen to expose it.

The reducer applies to a single feature vector or a matrix of them, reproduces
`DimensionReduce` on its own training data, and answers `"Method"`, `"FeatureCount"` and
`"ReducedDimension"`. Only `"PrincipalComponentsAnalysis"` is available as a reducer so
far.
