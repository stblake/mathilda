# Machine Learning

Kernels for this category live in [`src/ml/`](../../../src/ml/) rather than inside the
builtin that first needed them, and that is deliberate. Each is written against a
row-major `n × dim` buffer of machine doubles rather than against an `Expr`, so a
second consumer can reach it: `src/ml/gmm.c`'s EM fit serves both
`FindClusters[…, Method -> "GaussianMixture"]` and (in prospect) `LearnDistribution`,
and `src/ml/pca.c`'s column statistics serve `Standardize`, `PrincipalComponents`,
and any later feature scaling. `src/list/find_clusters.c`, by contrast, holds 55
static functions and exports exactly one symbol, which is why `ClusteringComponents`,
`ClusteringTree`, `Dendrogram` and `NearestNeighborGraph` cannot yet be built on it.

**Results are plain `List`s.** The machine bridge's `na_build_matrix` returns a
*visible* `NDArray`, whose head is `NDArray`, so a result built that way compares
`False` against the literal list a user would write — while `Inverse`, `Dot` and
`LinearSolve` compare `True`. These builtins construct `List`s and let the evaluator's
own packing gate decide whether the result is held as a buffer, keeping the surface
consistent with the rest of the system.

## Standardize

Shifts each column to zero mean and rescales it to unit standard deviation.
Attributes: `Protected`.

- `Standardize[data]`

**Features**:
- **Columns are variables, rows are observations.** A flat list is treated as `n`
  observations of *one* variable, not one observation of `n`.
- The divisor is `n - 1` (the sample standard deviation), matching
  `StandardDeviation` — so `Standardize[x]` agrees with
  `(x - Mean[x])/StandardDeviation[x]` written out by hand. A mismatch here would be
  invisible on the mean but not on the scale.
- **A constant column becomes exactly `0`, not `Indeterminate`.** Zero variance
  carries no information, so "no deviation from the mean" is the honest value;
  dividing by the zero standard deviation would propagate `Indeterminate` through
  every reduction over the row.

```mathematica
In[1]:= Standardize[{1., 2., 3., 4.}]
Out[1]= {-1.1619, -0.387298, 0.387298, 1.1619}

In[2]:= Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]
Out[2]= {{-1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}}

In[3]:= Standardize[{{1., 5.}, {2., 5.}, {3., 5.}}]
Out[3]= {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}
```

## PrincipalComponents

Gives the rows of a matrix in principal-component coordinates, components ordered by
decreasing variance. Attributes: `Protected`.

- `PrincipalComponents[matrix]`
- `PrincipalComponents[matrix, Method -> "Covariance" | "Correlation"]`

**Features**:
- Rows are observations, columns are variables. A flat list declines: one variable has
  no components to rotate.
- **The transform is an orthogonal rotation, so total variance is preserved** — it is
  redistributed into the leading components, not created or destroyed. Rank-deficient
  input therefore puts exactly zero variance in the trailing components: five points
  on a line give a second coordinate of `0`.
- `Method -> "Correlation"` standardises each variable to unit variance first, which
  is what you want when the columns have incommensurable units — otherwise the
  variable with the largest raw scale dominates for no statistical reason. The default
  `"Covariance"` does not. Note that "variance explained" means a different thing
  under each, so eigenvalues are not comparable across the two.
- An unrecognised `Method` returns unevaluated rather than silently choosing, so a
  typo is visible instead of quietly changing the statistics.
- **Eigenvector signs are canonical.** An eigenvector is defined only up to sign, and
  LAPACK's `dsyev` and the in-house Jacobi fallback do not agree on which they return;
  each component is flipped so its largest-magnitude loading is positive, so the
  output does not depend on whether the binary was linked against LAPACK.

```mathematica
In[1]:= PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]
Out[1]= {{-2.82843, 0.0}, {-1.41421, 0.0}, {0.0, 0.0}, {1.41421, 0.0}, {2.82843, 0.0}}

In[2]:= p = PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}];
        {Variance[Map[First, p]], Variance[Map[Last, p]]}
Out[2]= {13.4817, 0.518295}
```

## DimensionReduce

Reduces each row of a matrix to a lower-dimensional representation.
Attributes: `Protected`.

- `DimensionReduce[data, k]`
- `DimensionReduce[data, k, Method -> m]`

**Features**:
- The three methods are **one algorithm with three ways of forming the symmetric
  matrix to decompose**, which is why they share the eigendecomposition rather than
  each carrying its own linear algebra:

| Method | What it decomposes |
|---|---|
| `"PrincipalComponentsAnalysis"` (default) | the covariance of the **centred** columns |
| `"LatentSemanticAnalysis"` | the Gram matrix `X'X`, **without** centring — a truncated SVD |
| `"MultidimensionalScaling"` | the double-centred squared-distance matrix (classical Torgerson scaling) |

- **Skipping the centring is the entire difference between PCA and LSA.** A
  term-document matrix is sparse and non-negative, and centring destroys both
  properties along with the meaning of a zero entry — which is why LSA does not.
- **Reducing to `k` gives exactly the first `k` principal components**, not a
  separately-fitted `k`-component model: `DimensionReduce[data, 2]` equals
  `Map[Take[#, 2] &, PrincipalComponents[data]]`.
- **Classical MDS on Euclidean distances is the same embedding as PCA**, reached by a
  different route (an `n × n` double-centred distance matrix rather than a
  `dim × dim` covariance). That agreement is used as a cross-check on both in the test
  suite; it also means MDS earns its keep only when the distances come from somewhere
  other than the coordinates.
- `"MultidimensionalScaling"` is **capped at 2000 rows**, its matrix being `n × n` —
  the same order of ceiling, for the same reason, as `FindClusters`' `"Spectral"`.
- **Asking for more dimensions than the data supports returns unevaluated** rather
  than padding with zeros, since padding would look like a successful reduction to a
  caller checking only the shape. An unknown `Method`, a non-positive `k`, an omitted
  `k`, and a flat list all decline too.

**Not implemented**: Wolfram's `DimensionReduce` can also choose `k` itself and can
return a `DimensionReducerFunction` applicable to *new* data. The second is the
substantive gap — a reusable reducer is a trained model, and that representation is
being designed with the `Predict` family rather than invented twice.

```mathematica
In[1]:= d = {{1., 2., 3.}, {3., 5., 4.}, {4., 4., 8.}, {6., 9., 2.}, {7., 8., 9.}};
        DimensionReduce[d, 2]
Out[1]= {{-5.28583, 0.302136}, {-1.6734, -0.58098}, {0.141634, 3.22424},
         {1.80921, -4.66364}, {5.00839, 1.71825}}

In[2]:= DimensionReduce[d, 2, Method -> "LatentSemanticAnalysis"]
Out[2]= {{3.53625, 1.0399}, {7.03311, -0.290957}, {9.24375, 3.24698},
         {9.903, -4.78819}, {13.8764, 1.13662}}
```

## Clustering

`FindClusters` and its ten methods are documented under
[`lists-and-iteration.md`](lists-and-iteration.md), where the rest of the list
operations live. All ten cluster n-dimensional points; only string input is declined
above one dimension, having no coordinates.
