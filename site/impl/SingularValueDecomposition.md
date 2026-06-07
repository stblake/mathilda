---
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013)."
  - "J. Demmel and K. Veselić, \"Jacobi's Method is More Accurate than QR\", SIAM J. Matrix Anal. Appl. 13 (1992)."
source: src/linalg/svdecomp.c
---
**Algorithm.** `builtin_singularvaluedecomposition` returns `{u, sigma, v}` with `m == u.sigma.ConjugateTranspose[v]`, supporting `SingularValueDecomposition[m, k]`, `UpTo[k]`, the generalized `{m, a}` form, and `Tolerance`/`TargetStructure` options. `svd_dispatch` selects a kernel per numeric domain; all three feed `svd_apply_postprocess`, which centralises truncation, tolerance, and `TargetStructure` handling.

- **Exact / symbolic** (`svd_symbolic_dispatch` → `svd_symbolic_core`): forms the smaller Gram matrix `B = mᴴm` (p × p) when `n >= p`, else `B = m mᴴ` (n × n), and eigendecomposes it through the evaluator's `Eigenvalues`/`Eigenvectors`, with a 2×2 closed-form fast path. The singular values are `Sqrt[lambda]`; the eigenvectors give the primary factor (`v` if `mᴴm`, else `u`); the other factor is recovered as `m.v.Sigma⁻¹` (resp. `mᴴ.u/sigma_i`) for the non-zero singular values, with the remaining columns filled by `qr_symbolic_core`'s orthogonal completion to span the null space. If the eigendecomposition has no closed form, this path returns NULL.
- **Inexact, `min_bits <= 53`** (`svd_machine_dispatch`): LAPACK `dgesdd`/`zgesdd` (divide-and-conquer) for the standard form, `dggsvd3`/`zggsvd3` for the generalized `{m, a}` form.
- **Inexact, `min_bits > 53`** (`svd_mpfr_dispatch`): one-sided Jacobi SVD over MPFR arrays (Demmel-Veselić), preceded by a QR/Paige-Van Loan reduction.

**Data structures.** `Expr*` `List`-of-`List` matrices for the symbolic path (Gram matrix, eigenpairs, `Sqrt`-valued `sigma`); dense `double` / interleaved-complex LAPACK buffers for the machine path; arbitrary-precision MPFR arrays for the high-precision path. The exact path picks the smaller of `mᴴm` and `m mᴴ` to keep the eigenproblem small. Inexact input uses the rationalise / numericalise round-trip at the minimum input precision.

**Complexity / limits.** Exact path cost is dominated by the symbolic eigendecomposition of a min(n,p)-square matrix and only succeeds when that closes; the generalized `{m, a}` form has no symbolic kernel and requires the machine path. Machine path is LAPACK-bound (~O(n p · min(n,p))). The `TargetStructure -> "Structured"` head is not fully realised (results are returned dense).
