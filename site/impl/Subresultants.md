---
references:
  - "W. S. Brown and J. F. Traub, \"On Euclid's Algorithm and the Theory of Subresultants\", J. ACM 18(4), 1971."
  - "M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005)."
source: src/poly/subresultants.c
---
**Algorithm.** `builtin_subresultants` returns the list of principal subresultant coefficients (PSCs) of `poly1`, `poly2` w.r.t. `var`. The list has length `min(deg p1, deg p2) + 1`; element 0 equals `Resultant[p1, p2, var]`, and the first k PSCs vanish exactly when the polynomials share k roots (with multiplicity). Both arguments are checked polynomial in `var` (`internal_polynomialq`), then expanded.

The efficient path is `subresultants_prs`: the Bronstein subresultant polynomial remainder sequence with the gamma/beta/delta recurrence (`gamma_i = (-lc_{i-1})^{δ_{i-1}} γ_{i-1}^{1-δ_{i-1}}`, `beta_i = -lc_{i-1} γ_i^{δ_i}`, `r_{i+1} = prem(r_{i-1}, r_i)/β_i`) — the same recurrence as `Resultant` in `src/poly/poly.c`. It orients so `deg A ≥ deg B`, runs the chain via `pseudo_rem_standard` and `poly_divide_by_scalar`, retains each member's degree and leading coefficient, and extracts PSCs by the fundamental theorem of subresultants using the cumulative recurrence `s_p = lc(R_p)^{δ_p} / s_{p-1}^{δ_p-1}`, placing each at its degree index. Equal input degrees set `psc_L = 1` (empty top minor); a swap-correction sign `(-1)^{(n-j)(m-j)}` is applied when the inputs were oriented. Intermediate gamma/beta and PSC arithmetic goes through `internal_cancel` + `expr_expand`.

For algebraic-number coefficients (Sqrt, cube roots — `subres_has_algebraic`), where the pseudo-remainder chain bloats, it falls back to `subresultants_determinant`: `PSC_j = Det(M_j)`, where `M_j` is the Sylvester matrix restricted to the first `m-j` `poly1`-shift rows, `n-j` `poly2`-shift rows, and `n+m-2j` columns, evaluated via the `Det` builtin and `expr_expand`. The same determinant path is the fallback when the PRS path returns NULL.

**Data structures.** Descending coefficient arrays `Expr**` (`desc_coeffs`, leading coeff at index 0), built via the bulk `get_all_coeffs_expanded` extractor when possible. The PRS keeps parallel growable arrays `R[]` (chain members), `deg[]`, `lc[]`. Matrices are `List`-of-`List` `Expr*`.

**Complexity / limits.** The subresultant PRS keeps coefficient growth polynomial (the whole point of the Bronstein recurrence vs. naive PRS). The determinant fallback is `O(min(n,m))` minors, each a Sylvester-minor determinant. Identically-zero input is left unevaluated.
