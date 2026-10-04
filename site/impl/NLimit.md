---
source: src/numerical_calculus/nlimit.c
references:
  - "P. Wynn, *On a device for computing the e_m(S_n) transformation*, MTAC **10** (1956) 91–96 — the epsilon algorithm (SequenceLimit)."
  - "D. Levin, *Development of non-linear transformations for improving convergence of sequences*, Internat. J. Comput. Math. **3** (1973) 371–388."
---
**Algorithm.** `builtin_nlimit` parses `NLimit[expr, z -> z0]` and constructs a
**geometric sequence of sample points** approaching `z0`: for a finite target,
`z_k = z0 − d·Scale·2^{−k}` (the default `Direction -> Automatic == −1`
approaches from larger values); for an infinite target (`Infinity`,
`ComplexInfinity`, `DirectedInfinity[d]`, `I Infinity`, …), `z_k = u·Scale·2^k`
marches outward along the ray. Each sample binds `z` Block-style and bumps the
evaluation clock. The limit is recovered by **sequence acceleration**; `Method
-> Automatic` runs three accelerators and keeps the one with the smallest
internal convergence residual:

- **Richardson/Romberg** (`EulerSum`, the `2^j − 1` tableau shared with `ND`) —
  trusted even at an exact-zero residual (integer-power tails converge exactly);
- **Wynn's epsilon** (`SequenceLimit`, the iterated Shanks transform) at degrees
  `1 … (terms−1)/2` — admitted only with a strictly positive residual;
- **Levin's u-transform** — admitted only when the sample increments are
  contracting (otherwise it can collapse to a spurious finite value with a
  deceptively small residual).

`Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`) forces the chosen Levin
variant; the kernels live in `seqaccel.c` (shared with `NSum`). Two robustness
gates both return the form unevaluated rather than a meaningless extrapolant: a
**non-decaying oscillation** screen (`NLimit::osc`) that counts increment
reversals, then — only if the screen fires — resamples on a 20-octave ladder at
two offsets (one scaled by the golden ratio to defeat power-of-two aliasing) and
refuses when the far/near envelope ratio stays above 0.6; and a **noise/
divergence** gate (`NLimit::noise`) that rejects a non-finite result, a result
far larger than the sample scale, or successive extrapolates that do not settle.
Spurious tiny residuals are *not* chopped. Each method is also exposed directly
as `NLimit`Automatic`/`EulerSum`/`SequenceLimit`/`Levin`/`LevinU`/`LevinT`/`LevinV`.

**Data structures.** The sample sequence is a `double _Complex*` of the sample
count on the machine path, or paired `mpfr_t*` real/imag buffers under MPFR; the
oscillation ladder is a small `double` magnitude array. `WorkingPrecision`
(`MachinePrecision` or a digit count) selects the path. The limit point is
classified and its direction normalised to a unit vector before sampling.

**Complexity / limits.** `[start, 50]` symbolic samples with an `O(terms²)`
Richardson tableau, plus Wynn over `~terms/2` degrees and one Levin pass; the
oscillation stage 2 (only when triggered) adds 40 samples. Options: `Method`
(`Automatic` | `EulerSum` | `SequenceLimit` | `"Levin"` variants), `Direction`
(`Automatic == −1` or a complex vector), `Scale` (default 1), `Terms` (default
**13**, grown adaptively — the tableau depth, not the arithmetic precision, sets
the accuracy of a branch-point/fractional-power approach), `WynnDegree`,
`WorkingPrecision`, `AccuracyGoal` (default `MachinePrecision`), `PrecisionGoal`.
Diagnostics (`badmeth`, `notnum`, `osc`, `noise`, `nc_warn_goal`) route through
`mth_message`.
