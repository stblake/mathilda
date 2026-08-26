---
references:
  - "Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013)."
  - "D. K. Faddeev and V. N. Faddeeva, *Computational Methods of Linear Algebra* (W. H. Freeman, 1963) — the Faddeev–LeVerrier–Souriau recurrence."
source: src/linalg/charpoly.c
---
**Algorithm.** `builtin_characteristicpolynomial` (`src/linalg/charpoly.c`) is a thin builtin over the eigen module's characteristic-polynomial machinery — the same polynomial `Eigenvalues` solves for. `eigen_extract_matrix_pair` classifies the argument as an ordinary square matrix `m` or a generalised pencil `{m, a}` and validates squareness.

*Ordinary case* `CharacteristicPolynomial[m, x]`: the polynomial is formed by the **Faddeev–LeVerrier–Souriau** recurrence (`eigen_char_poly_faddeev`), which builds the coefficients in `O(n^4)` matrix multiplications — far cheaper than the naïve `Det[m - x I]`, which would face an `O(n!)` Laplace expansion of a symbolic-in-`x` matrix (so a 100×100 machine matrix is sub-second). That routine returns `det(λI − m)`, so for odd `n` the result is negated to match `Det[m − x I] = (−1)^n det(x I − m)` — hence odd-degree characteristic polynomials are monic-negative.

*Generalised case* `CharacteristicPolynomial[{m, a}, x]`: `eigen_build_lambda_matrix` forms the entry-wise matrix `m − λa` and `eigen_compute_det` takes its determinant by Laplace expansion, which already carries the correct sign. A null space shared by `m` and `a` drops the leading term(s), so an infinite generalised eigenvalue shows as a degree deficit.

The polynomial is built in a private internal lambda symbol; the user's second argument (a symbol, a number, or any expression) is then substituted for it via `ReplaceAll` and the result is expanded (`Det`/Laplace do not multiply out). Entries may be integer, rational, machine- or arbitrary-precision real, complex, or symbolic.

**Complexity / limits.** Faddeev–LeVerrier is `O(n^4)` for the ordinary case; the generalised case is Laplace expansion (`O(n!)`, in practice used only for small pencils, matching the generalised `Eigenvalues` path). The result is a symbolic `Plus`, so the head is exempt from the packed-array/`Compile` surfaces, but it accepts a packed or visible `NDArray` matrix as input (the buffer is materialised before the symbolic recurrence runs).
