---
references:
  - "I. T. Jolliffe, *Principal Component Analysis*, 2nd ed. (Springer, 2002)."
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013), §8.5 (the cyclic Jacobi eigenvalue method)."
source: src/ml/pca.c
---
**Algorithm.** `builtin_principal_components` reads an `n × dim` matrix (a flat list
declines — one variable has nothing to rotate) and calls `ml_pca`. That routine centres
the columns with `ml_standardize`, and for `Method -> "Correlation"` *also* rescales
each column to unit variance — the single flag that is the whole difference between a
covariance and a correlation analysis. It forms the symmetric `dim × dim` covariance
(`n − 1` divisor) and diagonalises it with `ml_sym_eigen_desc`, which calls LAPACK's
`dsyev` when linked and otherwise falls back to an in-house cyclic Jacobi sweep
(`ml_jacobi`, written here because the linalg symmetric solvers are static to their
units). Eigenpairs are then selection-sorted into **descending** eigenvalue order (so
the leading component carries the most variance regardless of which backend ran), and
each eigenvector is **sign-canonicalised** — flipped so its largest-magnitude loading is
positive — because an eigenvector is defined only up to sign and `dsyev` and Jacobi
disagree on which they return. Finally each centred row is projected onto the ordered
eigenvectors. An unrecognised `Method` returns unevaluated rather than silently choosing.

**Data structures.** Plain `double` buffers throughout: a centred copy of the data, the
`dim × dim` covariance, an eigenvalue vector and a row-major eigenvector matrix. The
transformed rows are returned as a plain `List` (`ml_list_matrix`).

**Complexity / limits.** `O(n · dim²)` to form the covariance and `O(dim³)` to
diagonalise (the Jacobi fallback runs up to 100 sweeps). Because the transform is an
orthogonal rotation, total variance is preserved — rank-deficient input puts exactly
zero variance in the trailing components. Eigenvalues are *not* comparable between the
Covariance and Correlation methods, since "variance explained" means a different thing
under each.
