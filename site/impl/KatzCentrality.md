---
references:
  - "L. Katz, *A new status index derived from sociometric analysis*, Psychometrika **18** (1953) 39-43."
source: src/graph/gmet_spectral.c
---
**Algorithm.** `builtin_katz_centrality` solves the Katz fixed point
`x = b + α Aᵀx` (Katz 1953), where `α` is the attenuation factor (second argument)
and `b` the baseline (the third argument — a scalar or a length-`n` list —
defaulting to the all-ones vector). The transpose `Aᵀ` is realised directly: the
graph is built as an in-arc CSR (`GMET_IN`), so each vertex accumulates over its
in-neighbours. The primary path is a fixed-point iteration
`xₙₑₓₜ = b + α Aᵀx`, up to 20000 sweeps, declared converged when the L1 change
falls to `≤ 1e-15` of the L1 norm; a growth guard aborts after 50 consecutive
increases (a divergent `α`). On non-convergence with `n ≤ 4000` it falls back to a
dense LAPACK solve (`mat_lapack_dgesv`) of `(I − α Aᵀ) x = b` in column-major form.

**Data structures.** A `GMET_IN` `GmetCSR`, two `double` work vectors for the
iteration, and — only on the dense fallback — an `n × n` column-major matrix `M`
holding `I − α Aᵀ` (then its LU factors) plus a pivot array. Results are memoised
through `gmet_cache`, keyed on the whole call.

**Complexity / limits.** `O(iter · E)` for the iteration, `O(n³)` for the dense
fallback (hence the `n ≤ 4000` cap). Edge weights are ignored. Two cases short-
circuit to `b` exactly (with exact entries preserved): a zero `α`, and a graph with
no edges. The fallback carries a singularity guard — if the smallest LU pivot is
below `1e-12` of the largest, the system is treated as singular (an `α` at a
reciprocal eigenvalue of `A`) and the call is left unevaluated rather than
returning a `~1e16` garbage answer.
