---
references:
  - "Erwin H. Bareiss, \"Sylvester's Identity and Multistep Integer-Preserving Gaussian Elimination\", Mathematics of Computation 22 (1968)."
source: src/linalg/linsolve.c
---
**Algorithm.** `builtin_linearsolve` solves `m . x = b` and shares its method dispatcher (`MatsolMethod`, `matsol_parse_method_option`) with `RowReduce` and `Inverse`. An optional `Method` argument selects:

- `"DivisionFreeRowReduction"` (default / `Automatic`): `linearsolve_divfree` performs Bareiss-like **fraction-free Gauss-Jordan** elimination on the augmented matrix `[m | b]` with a running pivot product to keep arithmetic exact (no GCD blow-up). After reduction it reads off a particular solution; the implementation handles rectangular/under-determined systems (pivot column bookkeeping per row) and inconsistent systems (emitting the appropriate diagnostic), and supports vector, matrix, and higher-rank RHS by tracking the trailing dimensions of `b`.
- `"OneStepRowReduction"`: classical Gauss-Jordan, one division per pivot.
- `"CofactorExpansion"`: Cramer's rule via Laplace cofactor expansion (small square non-singular systems).

The RHS shape is normalised: a rank-1 `b` against a rank-2 `m` returns a flat solution vector; a matrix RHS returns shape `{c, k}`; higher-rank `b` returns `{c, trail_dims…}`. Non-rectangular `m` emits a `LinearSolve` shape diagnostic.

**Data structures.** Dense flat `Expr**` augmented matrix `r × (c + k)`, row-major; `int* pivot_col_for_row` records pivot columns; the pivot product is a single `Expr*`. All arithmetic is symbolic/exact via the evaluator, so integer, rational, and symbolic systems solve exactly.

**Complexity / limits.** Fraction-free elimination is `O(r·c·(c+k))` arithmetic ops with bounded intermediate-coefficient growth (Bareiss). Inexact input is handled by the same exact pipeline used elsewhere in `linalg` (rationalise → solve → numericalise) rather than a floating-point LU.
