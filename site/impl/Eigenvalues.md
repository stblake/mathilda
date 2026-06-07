---
references:
  - "Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013)."
  - "J. H. Wilkinson, *The Algebraic Eigenvalue Problem* (Oxford University Press, 1965)."
source: src/linalg/eigen.c
---
**Algorithm.** `builtin_eigenvalues` (dispatcher in `eigen.c`; kernels in `eigen_common.c`, `eigen_direct.c`, `eigen_banded.c`, `eigen_arnoldi.c`, `eigen_feast.c`) chooses between an exact/symbolic path and numerical kernels based on whether the matrix is inexact.

*Exact / symbolic path* (`eigen_compute_eigenvalues_full`): the characteristic polynomial is formed by the **Faddeev–Leverrier–Souriau** recurrence (`eigen_char_poly_faddeev`), which builds the coefficients in `O(n^4)` matrix multiplications — far cheaper than Laplace expansion of `det(m − λI)` over a polynomial-entry matrix (`O(n!)`) once `n` grows. Its roots are found by `eigen_solve_poly` (radical formulas for cubics/quartics controlled by `Cubics`/`Quartics`, otherwise `Root`/`Solve`), extracted with multiplicity, and sorted by descending `Abs`. The **generalised** problem `Eigenvalues[{m, a}]` instead forms `det(m − λa)` by Laplace expansion (only used for small pencils) and pads the short result with `Infinity` for the degree-drop branch.

*Numerical path* (inexact input): dispatched by `Method`. `Direct` (`direct_dispatch`) uses **Householder reduction** to tridiagonal (symmetric, Golub & Van Loan Alg. 8.3.1) followed by **implicit-shift symmetric tridiagonal QR with Wilkinson shift** (`direct_symtridiag_qr`) for symmetric matrices, and Householder reduction to upper **Hessenberg + Francis implicit double-shift QR** for non-symmetric matrices. `Automatic` prefers `Banded` for narrow-band Hermitian input and `Arnoldi` when only a small `k` is requested, with FEAST available for Hermitian interval problems; each dispatcher returns `NULL` for shapes it doesn't support and falls through to `Direct`, then ultimately to the symbolic pipeline. Numerical-noise imaginary parts are chopped (`eigen_chop`). A `{k}`/`-k`/`UpTo[k]` spec trims the sorted result (`eigen_apply_k_spec`).

**Data structures.** Symbolic side: `Expr` trees for the polynomial and roots. Numerical kernels operate on dense row-major `double` (and MPFR) matrices; the source carries `LAPACK-HOOK` annotations marking where `dsytrd`/`dsteqr`/`dgehrd`/`dhseqr` would drop in under a `USE_LAPACK` build (the hooks are present but no LAPACK backend is currently wired).

**Complexity / limits.** Faddeev–Leverrier is `O(n^4)`; symmetric QR is `O(n^3)` with cubic per-eigenvalue convergence under the Wilkinson shift. The generalised path is restricted to small `n` (Laplace expansion). Closed-form eigenvalues are limited by the degree-≤4 radical solver; higher degrees come back as `Root` objects.
