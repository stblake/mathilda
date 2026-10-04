# DimensionReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DimensionReduce[data, k] reduces each row of data to k dimensions. Method -> "PrincipalComponentsAnalysis" (default) centres the columns and projects onto the leading eigenvectors of the covariance; "LatentSemanticAnalysis" skips the centring, giving a truncated SVD, which is what a sparse non-negative term-document matrix wants; "MultidimensionalScaling" double-centres the squared distance matrix (classical Torgerson scaling) and is capped at 2000 rows, its matrix being n x n. Asking for more dimensions than the data supports returns unevaluated rather than padding with zeros.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= d = {{1., 2., 3.}, {3., 5., 4.}, {4., 4., 8.}, {6., 9., 2.}, {7., 8., 9.}}; DimensionReduce[d, 2]
Out[1]= {{-5.28583, 0.302136}, {-1.6734, -0.58098}, {0.141634, 3.22424}, {1.80921, -4.66364}, {5.00839, 1.71825}}
```

### Options (1)

```mathematica
In[2]:= DimensionReduce[d, 2, Method -> "LatentSemanticAnalysis"]
Out[2]= {{3.53625, 1.0399}, {7.03311, -0.290957}, {9.24375, 3.24698}, {9.903, -4.78819}, {13.8764, 1.13662}}
```

### Applications (5)

A 5x3 data matrix

```mathematica
In[3]:= d = {{1., 2., 3.}, {3., 5., 4.}, {4., 4., 8.}, {6., 9., 2.}, {7., 8., 9.}}
Out[3]= {{1.0, 2.0, 3.0}, {3.0, 5.0, 4.0}, {4.0, 4.0, 8.0}, {6.0, 9.0, 2.0}, {7.0, 8.0, 9.0}}
```

Default PCA: reduce each row to two dimensions

```mathematica
In[4]:= DimensionReduce[d, 2]
Out[4]= {{-5.28583, 0.302136}, {-1.6734, -0.58098}, {0.141634, 3.22424}, {1.80921, -4.66364}, {5.00839, 1.71825}}
```

Reducing to k is exactly the first k principal components

```mathematica
In[5]:= DimensionReduce[d, 2] == Map[Take[#, 2] &, PrincipalComponents[d]]
Out[5]= True
```

LSA skips the centring: a truncated SVD

```mathematica
In[6]:= DimensionReduce[d, 2, Method -> "LatentSemanticAnalysis"]
Out[6]= {{3.53625, 1.0399}, {7.03311, -0.290957}, {9.24375, 3.24698}, {9.903, -4.78819}, {13.8764, 1.13662}}
```

Classical MDS, the PCA embedding up to a reflection

```mathematica
In[7]:= DimensionReduce[d, 2, Method -> "MultidimensionalScaling"]
Out[7]= {{5.28583, -0.302136}, {1.6734, 0.58098}, {-0.141634, -3.22424}, {-1.80921, 4.66364}, {-5.00839, -1.71825}}
```

## Options & behaviour

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

### Not implemented

Wolfram's `DimensionReduce` can also choose `k` itself and can
return a `DimensionReducerFunction` applicable to *new* data. The second is the
substantive gap — a reusable reducer is a trained model, and that representation is
being designed with the `Predict` family rather than invented twice.

## Implementation notes

**Algorithm.** `builtin_dimension_reduce` validates `k` (a positive integer; omitting
it, or asking for more dimensions than the data supports, declines rather than padding
with zeros) and dispatches to `ml_reduce`, where the three methods are **one algorithm
with three ways of forming the symmetric matrix to decompose**, sharing
`ml_sym_eigen_desc`:

- `"PrincipalComponentsAnalysis"` (default) centres the columns (`ml_standardize`, no
  rescale), forms the `dim × dim` scatter matrix, eigendecomposes, and projects each
  centred row onto the leading `k` eigenvectors.
- `"LatentSemanticAnalysis"` does the same **without** centring — the Gram matrix
  `XᵀX`, i.e. a truncated SVD. Skipping the one `if (method == ML_REDUCE_PCA)
  ml_standardize(...)` is the entire mathematical difference, and it is what keeps a
  sparse, non-negative term-document matrix sparse and non-negative.
- `"MultidimensionalScaling"` is classical (Torgerson) scaling: build the `n × n`
  squared-distance matrix, double-centre it as `B = −½ J D² J` with `J = I − 11ᵀ/n`
  (turning squared distances back into inner products), eigendecompose, and take
  coordinate `t` of row `i` as `eigenvectorₜ[i] · √λₜ` (a non-positive eigenvalue
  contributes `0`, not the root of a negative number). On Euclidean distances this is
  the same embedding as PCA, which the test suite uses as a mutual cross-check.

**Data structures.** `double` buffers; PCA/LSA work with a `dim × dim` scatter matrix,
MDS with an `n × n` matrix `B`. Rows are assembled into a plain `List`.

**Complexity / limits.** PCA/LSA are `O(n · dim² + dim³)`; MDS is `O(n² · dim + n³)`,
which is why it is capped at `ML_MDS_MAX_N = 2000` rows — the same order of ceiling, for
the same quadratic-memory reason, as `FindClusters`' `"Spectral"` method. Reducing to
`k` gives exactly the first `k` principal components; there is no separately-fitted
`k`-component model.

- The three methods are **one algorithm with three ways of forming the symmetric
  matrix to decompose**, which is why they share the eigendecomposition rather than
  each carrying its own linear algebra:

**Attributes:** `Protected`.

## References

**See also:** [FindClusters](../../lists-and-iteration/FindClusters/), [DimensionReducerFunction](../../other-advanced/DimensionReducerFunction/), [Predict](../../machine-learning/Predict/)

- W. S. Torgerson, *Multidimensional scaling: I. Theory and method*, Psychometrika **17** (1952) 401-419 (classical / Torgerson scaling).
- T. F. Cox and M. A. A. Cox, *Multidimensional Scaling*, 2nd ed. (Chapman & Hall/CRC, 2001).
- S. Deerwester et al., *Indexing by latent semantic analysis*, J. Amer. Soc. Inf. Sci. **41** (1990) 391-407.
- Source: [`src/ml/pca.c`](https://github.com/stblake/mathilda/blob/main/src/ml/pca.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)

## Notes & additional examples

### Notes

The three methods are one algorithm with three ways of forming the symmetric matrix to
decompose: PCA uses the covariance of the centred columns, LSA the Gram matrix `XᵀX`
without centring, and MDS the double-centred squared-distance matrix. Skipping the
centring is the *entire* difference between PCA and LSA — a term-document matrix is
sparse and non-negative, and centring would destroy both.

Classical MDS on Euclidean distances is the same embedding as PCA, reached by a
different route (an `n × n` distance matrix rather than a `dim × dim` covariance); it is
identical up to the sign of each axis, and earns its keep only when the distances come
from somewhere other than the coordinates. MDS is capped at 2000 rows because its matrix
is `n × n`.

Asking for more dimensions than the data supports returns unevaluated rather than
padding with zeros, which would look like a successful reduction to a caller checking
only the shape.
