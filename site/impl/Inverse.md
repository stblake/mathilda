---
references:
  - "Erwin H. Bareiss, \"Sylvester's Identity and Multistep Integer-Preserving Gaussian Elimination\", Mathematics of Computation 22 (1968)."
source: src/linalg/inv.c
---
**Algorithm.** `builtin_inverse` validates that the argument is a non-empty square matrix, parses an optional `Method`, and dispatches (`MatsolMethod`) to one of three exact workers:

- `inverse_divfree` (default, `Method -> Automatic` / `"DivisionFreeRowReduction"`): Bareiss-like **fraction-free Gauss-Jordan** elimination on the augmented matrix `[A | I]`. A running pivot product `P` divides each updated entry so the elimination stays division-free (exact, no rational/GCD blow-up); the right half becomes `A^{-1}` once the left half is reduced to (a scalar multiple of) the identity. Singular matrices emit `Inverse::sing` and return `NULL`.
- `inverse_onestep` (`"OneStepRowReduction"`): classical Gauss-Jordan with one division per pivot, each entry canonicalised via `Together` so symbolic cancellations are still detected.
- `inverse_cofactor` (`"CofactorExpansion"`): the adjugate/determinant formula `A^{-1}[i,j] = (-1)^{i+j} det(M_{j,i})/det(A)`, with each minor computed by the same Laplace expansion used by `Det` — `O(n!·n^2)`, for tiny `n` only.

Inexact (`Real`/MPFR) matrices are routed through the standard `common_scan_inexact` → `common_rationalize_input` pipeline: the matrix is rationalised at the minimum precision present, inverted exactly, then numericalised back to that precision, so the rank/pivot decisions are exact.

**Data structures.** Dense flat `Expr**` augmented matrix of `n × 2n` element pointers, row-major; the pivot product is a single shared `Expr*`. The same `inv.c` module also implements `PseudoInverse` via a full-rank `B·C` decomposition from `RowReduce`.

**Complexity / limits.** Fraction-free Gauss-Jordan is `O(n^3)` arithmetic ops with controlled intermediate growth; cofactor expansion is `O(n!)`. There is no machine-precision LU `dgetrf`-style kernel — machine matrices are handled by exact rationalisation rather than floating-point LU.
