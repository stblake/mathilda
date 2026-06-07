---
references:
  - "W. S. Brown and J. F. Traub, \"On Euclid's Algorithm and the Theory of Subresultants\", J. ACM 18(4), 1971."
  - "M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005)."
source: src/poly/subresultantpoly.c
---
**Algorithm.** `builtin_subresultantpolynomials` returns the list of subresultant *polynomials* `{S_0, ..., S_m}` of `poly1`, `poly2` w.r.t. `var`, where `m = Exponent[poly2, var]` (requires `deg p1 ≥ deg p2`, exact coefficients). `S_0` is `Resultant[p1, p2, var]`, and the coefficient of `var^j` in `S_j` is the j-th principal subresultant coefficient (matching `Subresultants[...][[j+1]]`).

By the fundamental theorem of subresultants each `S_j` is either zero or a scalar multiple of a single member of the subresultant polynomial remainder sequence. The implementation reuses the same Bronstein gamma/beta/delta PRS as `Subresultants` and `Resultant`, then classifies each output index: a *regular* index (`j == deg(R_p)` at a strict-drop step) gives `S_j = (psc_j / lc(R_p)) · R_p` — the chain member rescaled so its leading coefficient is `psc_j = clc_p^{δ_p}` (leaving `S_j = R_p` when `δ_p == 1`); a *defective* index (across a degree gap, `δ_p > 1`) is computed directly from the determinant-polynomial definition (small Sylvester minor since these sit high in the chain); all other indices are zero. For algebraic-number coefficients the whole list is built from the determinant-polynomial definition, mirroring `Resultant`/`Subresultants`.

**Data structures.** Descending coefficient arrays `Expr**` (`desc_coeffs`); the determinant-polynomial path builds truncated shifted polynomials via `trunc_shift_poly`. Output is a `List` of polynomial `Expr*`.

**Complexity / limits.** Same coefficient-growth control as the subresultant PRS; defective and algebraic-coefficient indices fall back to small-minor determinant evaluation.
