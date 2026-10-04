---
source: src/stats/corrcov.c
---
**Algorithm.** `builtin_correlation` dispatches through `corrcov_dispatch(res, 1,
"Correlation")` — the same engine as `Covariance` with the `correlate` flag set —
accepting one to three arguments (`Correlation::argb` otherwise) and routing an
`NDArray` / packed argument to `nd_correlation` (`src/linalg/ndcorrcov.c`). A
correlation is a covariance normalized by the standard deviations:

- **`Correlation[v, w]`** — `g_corr_vectors` computes `Covariance[v, w]` and
  divides by `StandardDeviation[v] · StandardDeviation[w]` (the `n−1` factors
  cancel). For real data `−1 ≤ ρ ≤ 1`; exact / complex / symbolic data are carried
  through exactly as in `Covariance`.
- **`Correlation[a, b]`** — the `p×q` cross-correlation of the columns of two
  matrices; `g_matrix` divides each covariance cell by the two column standard
  deviations (precomputed per column in `sds_a` / `sds_b`).
- **`Correlation[a]`** — the `p×p` auto-correlation, symmetric with a unit
  diagonal. The diagonal is set outright rather than computed: it is the `Real`
  `1.` when the data has a machine real (`common_has_machine_real`), keeping the
  matrix one uniform type so `SymmetricMatrixQ` and `==` stay decidable, and the
  exact `Integer` `1` for exact / symbolic data.

**Data structures.** Shares `Covariance`'s `g_cov_vectors` / `g_matrix`; adds
`Expr**` arrays of per-column standard deviations for the normalization. Results
are `List`-of-`List` for the matrix forms.

**Complexity / limits.** `O(n)` for a vector pair; `O(p·q·n)` for the matrix forms,
plus `O(p + q)` standard-deviation reductions. On the `AWARE` list in `src/pack.c`
but not `INT64_OK` (shares `Covariance`'s `Rational`-degradation on an integer
buffer). Lowered inside `Compile[]`. Stays unevaluated for a single vector,
mismatched shapes, or fewer than two observations.
