---
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013)."
  - "P. A. Businger and G. H. Golub, \"Linear Least Squares Solutions by Householder Transformations\", Numer. Math. 7 (1965)."
source: src/linalg/qrdecomp.c
---
**Algorithm.** `builtin_qrdecomposition` returns the thin QR `{q, r}` (or `{q, r, p}` with `Pivoting -> True`) satisfying `m == ConjugateTranspose[q].r`. `qr_dispatch` routes inexact `min_bits <= 53` input to `qr_machine_dispatch` (LAPACK `dgeqrf`/`dgeqp3`, or `zgeqrf`/`zgeqp3` for complex, with the standard rank-revealing cutoff), higher precision to `qr_mpfr_dispatch` (textbook Householder QR over MPFR), and exact integer/rational/complex/symbolic input to `qr_symbolic_dispatch`. Any fast-path failure (LAPACK off, non-coercible leaf, non-zero `info`, rank-deficient with no pivoting) falls through to the symbolic kernel.

The symbolic core `qr_symbolic_core` is **Modified Gram-Schmidt** on the columns of `A`, evaluated symbolically: for each column it subtracts the projections `coeff = <Q[:,j], v>` onto each existing orthonormal column (storing `coeff` into `R[j,k]`), takes `norm = Sqrt[<v,v>]`, and appends `v/norm` as a new orthonormal column; a zero residual norm marks a dependent column that is skipped (so `r = rank`). `q` is built as the conjugate transpose of `Q` so the result is `Q^T` for real input and the proper Hermitian transpose for complex. With `Pivoting -> True` the next column chosen at each step is the remaining one whose residual projection has the largest squared norm — Businger-Golub column pivoting expressed in MGS form, making the diagonal of `R` decrease in magnitude — and the order is inflated into a `p × p` permutation matrix with `p[perm[j], j] = 1`.

**Data structures.** Flat row-major `Expr**` working buffers for `Q` (n × r) and `R` (r × p); inner products, norms, and scalings go through the evaluator via `eval_and_free`, leaving exact `Sqrt[...]` / rational / symbolic entries. The `is_definitely_zero` test (`Together` then `is_zero_poly`) detects dependent columns. Inexact input uses the `common_rationalize_input` / `common_numericalize_result` round-trip at the minimum input precision.

**Complexity / limits.** O(n·p·r) evaluator-level vector operations for the symbolic path; symbolic norms accumulate `Sqrt` nestings. The machine path is LAPACK-bound.
