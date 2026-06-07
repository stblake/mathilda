---
references:
  - "C. C. Paige, M. A. Saunders, \"LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares\", ACM TOMS 8 (1982)."
  - "Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013)."
source: src/linalg/lstsq.c
---
**Algorithm.** `builtin_leastsquares` returns an `x` minimising `Norm[m.x - b]`. The default `"Direct"` method computes the Moore–Penrose solution `x = PseudoInverse[m] . b` (reusing the `PseudoInverse` builtin, which itself does an exact full-rank `B·C` decomposition via `RowReduce`), so for full-column-rank `m` it gives the unique minimiser and for rank-deficient `m` the minimum-norm minimiser. The `LeastSquares == PseudoInverse . b` identity is the fundamental specification, and `"Direct"` works on every input family (integer, rational, symbolic, Real/MPFR, complex). Optional iterative methods are also provided and dispatched by input type:

- `"IterativeRefinement"` — residual-correction loop on top of Direct (one pass for exact input; drives round-off to `Tolerance` for inexact).
- `"LSQR"` — Paige–Saunders LSQR via Lanczos bidiagonalisation with Givens rotation updates, using their `|φ̄·α_{k+1}|` estimate of `‖Aᵀr‖`; symbolic input falls back to Direct, exact/complex input to CGLS, pure-real inexact to the canonical double-precision algorithm.
- `"Krylov"` — Conjugate-Gradient on the normal equations (CGLS / Hestenes–Stiefel), restricted to numeric inputs; symbolic falls back to Direct.

`Method` and `Tolerance` options may appear in either order.

**Data structures.** Solutions are built from the generic `PseudoInverse`/`Dot`/`Plus`/`Times` evaluator pipeline; the iterative kernels operate on the numeric tower (Real/MPFR/Complex) with exact-arithmetic variants for Integer/Rational inputs to avoid square-root growth.

**Complexity / limits.** Direct cost is dominated by `PseudoInverse` (exact Gauss-Jordan, `O(mn·rank)` style). Krylov/LSQR iterate to a `2·cols + O(1)` cap with `Tolerance`-based stopping; for symbolic inputs only Direct is well-defined since the iterative stopping tests are undecidable.
