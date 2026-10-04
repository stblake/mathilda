---
references:
  - "W. S. Torgerson, *Multidimensional scaling: I. Theory and method*, Psychometrika **17** (1952) 401-419 (classical / Torgerson scaling)."
  - "T. F. Cox and M. A. A. Cox, *Multidimensional Scaling*, 2nd ed. (Chapman & Hall/CRC, 2001)."
  - "S. Deerwester et al., *Indexing by latent semantic analysis*, J. Amer. Soc. Inf. Sci. **41** (1990) 391-407."
source: src/ml/pca.c
---
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
