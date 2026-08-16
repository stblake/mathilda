# DimensionReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DimensionReduce[data, k] reduces each row of data to k dimensions. Method -> "PrincipalComponentsAnalysis" (default) centres the columns and projects onto the leading eigenvectors of the covariance; "LatentSemanticAnalysis" skips the centring, giving a truncated SVD, which is what a sparse non-negative term-document matrix wants; "MultidimensionalScaling" double-centres the squared distance matrix (classical Torgerson scaling) and is capped at 2000 rows, its matrix being n x n. Asking for more dimensions than the data supports returns unevaluated rather than padding with zeros.`**

## Examples

_No verified examples yet for this function._

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

- The three methods are **one algorithm with three ways of forming the symmetric
  matrix to decompose**, which is why they share the eigendecomposition rather than
  each carrying its own linear algebra:

**Attributes:** `Protected`.

## References

**See also:** [FindClusters](../../lists-and-iteration/FindClusters/), [DimensionReducerFunction](../../other-advanced/DimensionReducerFunction/), [Predict](../../machine-learning/Predict/)

- Source: [`src/ml/pca.c`](https://github.com/stblake/mathilda/blob/main/src/ml/pca.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)
