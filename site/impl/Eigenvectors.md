---
references:
  - "Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013)."
  - "J. H. Wilkinson, *The Algebraic Eigenvalue Problem* (Oxford University Press, 1965)."
source: src/linalg/eigen.c
---
**Algorithm.** `builtin_eigenvectors` shares its eigenvalue computation and `Method` dispatch with `Eigenvalues` (see that page: Faddeev–Leverrier characteristic polynomial for the exact/symbolic path; Householder tridiagonalisation + Wilkinson-shift symmetric QR, Hessenberg + Francis QR, and Banded/Arnoldi/FEAST kernels for the inexact path, with `LAPACK-HOOK` sites for `dtrevc`-style vector recovery). Once the eigenvalues are known and sorted by descending `Abs`, eigenvectors are obtained per eigenvalue by **null-space computation**: equal eigenvalues are collapsed into runs, the residual matrix `m − λI` (or `m − λa` for the generalised pencil; `a` itself for `Infinity` eigenvalues) is substituted and row-reduced, and `eigen_null_space` returns up to `multiplicity` basis vectors. A defective eigenvalue whose geometric multiplicity is short is padded in place with zero vectors so the `i`-th eigenvector lines up positionally with the `i`-th eigenvalue.

For inexact input the matrix is first **rationalised** (`common_rationalize_input` at the minimum precision present) so the rank defect needed to expose the eigenvector is not destroyed by floating-point round-off; the null-space and normalisation are done in exact rational arithmetic, then the result is numericalised back (`common_numericalize_result`) and each vector normalised to unit `Norm`. A `{k}`/`-k`/`UpTo[k]` spec trims the result.

**Data structures.** `Expr` matrices/vectors throughout the symbolic path; the residual substitution uses `ReplaceAll` of the internal `λ` symbol, and `eigen_null_space` drives the exact Gauss-Jordan row reducer. Numerical kernels use dense row-major `double`/MPFR buffers (shared with `Eigenvalues`).

**Complexity / limits.** Dominated by the eigenvalue solve plus one null-space (row-reduction) per distinct eigenvalue. Defective matrices yield fewer independent vectors than the eigenvalue multiplicity, made explicit by zero-vector padding; the generalised case is restricted to small `n`.
