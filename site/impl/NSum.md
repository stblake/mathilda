---
source: src/numerical_calculus/nsum.c
references:
  - "H. Cohen, F. Rodriguez Villegas and D. Zagier, *Convergence acceleration of alternating series*, Experiment. Math. **9** (2000) 3–12."
  - "D. Levin, *Development of non-linear transformations for improving convergence of sequences*, Internat. J. Comput. Math. **3** (1973) 371–388."
  - "P. Wynn, *On a device for computing the e_m(S_n) transformation*, MTAC **10** (1956) 91–96 — the epsilon algorithm."
  - "J. B. Keiper, *Numerical computation of infinite products*, Wolfram Research tech. report (1992) — the Exp[NSum[Log]] reduction used by NProduct."
---
**Algorithm.** `builtin_nsum` (HoldAll) parses the iterator `{i, imin, imax,
di}`, reindexes terms to `k = 0, 1, 2, …` with `x_k = imin + k·di` built and
evaluated under a Block-style index binding, and sums the first `NSumTerms`
(default 15) explicitly before accelerating the tail. `ns_choose_method` picks
from a sampled profile:

- **Euler–Maclaurin** (`ns_em_*`) for a monotone summand defined off the
  integers: `Σ ≈ (1/di)∫_N^∞ f + f(N)/2 − Σ_j B_{2j}/(2j)! · di^{2j-1}
  f^{(2j-1)}(N)`. The tail integral is a dedicated double-exponential (exp-sinh)
  quadrature (`dequad_halfline_*`), **not** `NIntegrate`; the derivative
  corrections are a hybrid — symbolic `D` while the derivative tree stays small,
  switching to Taylor coefficients recovered from a circle DFT (Cauchy's
  formula, as `NSeries` does) once it balloons — and the correction series is
  truncated asymptotically at its smallest term.
- **Cohen–Villegas–Zagier** (`ns_cvz_*`) for a strictly alternating real series:
  Chebyshev weights `d_n = ((3+√8)^n + (3+√8)^{-n})/2` in a single linear pass
  giving ≈ 2.54 n bits.
- **Wynn's epsilon** (`seqaccel.c`, the iterated Shanks transform) otherwise,
  with **Levin's u/t/v transform** as a last resort for logarithmically/
  algebraically convergent tails.

An integer-only summand such as `1/Prime[n]` is detected *behaviourally* — the
summand is probed at two non-integer points, and if neither is finite,
Euler–Maclaurin is forbidden (no continuous tail) and the series is extrapolated
instead (with a larger head-term count). Multidimensional sums nest an inner
`NSum` as the summand (dependent inner bounds see the outer index via HoldAll).
A large finite sum of decaying terms is computed as the difference of two
infinite tails. With `VerifyConvergence -> True` (default), a summand whose
magnitude is still rising far into the sampled tail emits `NSum::div` and
returns `ComplexInfinity`.

`NProduct` (`nprod.c`) is evaluated as **`Exp[NSum[Log[f], …]]`** (Keiper 1992):
it maps `NProductFactors -> NSumTerms` (so the default factor count is also 15),
`NProductExtraFactors -> NSumExtraTerms`, runs the inner `NSum` at ten guard
digits above the request (because `Exp` turns the exponent's absolute error into
the product's relative error), and rounds the final `Exp` back. A divergent
log-sum propagates as `ComplexInfinity`.

**Data structures.** Two parallel implementations gated on `USE_MPFR`. The
machine path is `double _Complex` throughout — partial-sum sequence `P[]`, the
Wynn ε-table `(terms+1)²`, Levin's `a`/`omega` arrays. The MPFR path keeps split
`mpfr_t` real/imag buffers at `2·target + 32` internal bits (absorbing
cancellation in near-1 Euler–Maclaurin samples), rounded back to the target.
`WorkingPrecision` (`MachinePrecision` or a digit count) selects between them;
the extrapolation sequence length scales with the bit count. A machine request
draws every term through one of three cached compiled programs (real, at
precision, and complex-input for the contour), so a single compile serves all
methods; the MPFR path stays on the interpreter per term.

**Complexity / limits.** Linear in the head terms plus `O(seq²)` for the
ε-table (`NS_MAX_SEQ = 64`). Options: `Method` (`Automatic` | `EulerMaclaurin` |
`AlternatingSigns` | `WynnEpsilon` | `"Levin"`/`"LevinU"`/`"LevinT"`/`"LevinV"`),
`WorkingPrecision`, `NSumTerms`, `NSumExtraTerms`, `WynnDegree`,
`VerifyConvergence`, `AccuracyGoal` (default `MachinePrecision`), `PrecisionGoal`
(default Automatic); a shortfall against the combined tolerance warns
`NSum::ncvg`. NSum declines (stays unevaluated) on a non-numeric finite bound, a
summand that never numericalises, or a non-decaying finite sum beyond
2 000 000 terms. Per-term internal probes are muted via
`arith_warnings_mute`, and all diagnostics route through `mth_message`.
