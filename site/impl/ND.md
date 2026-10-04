---
source: src/numerical_calculus/nderiv.c
references:
  - "W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §5.7 — Richardson extrapolation of finite differences."
---
**Algorithm.** `builtin_nd` handles `ND[expr, x, x0]` and `ND[expr, {x, n}, x0]`
(threading manually over a `List` first argument). Each sample binds `x` to a
numeric value as a temporary OwnValue (Block-style, saving and restoring the
symbol's OwnValues and attributes), bumps the evaluation clock so memoization
cannot return a stale value, and evaluates then numericalises the held
expression. Two methods:

- **`Method -> EulerSum`** (default): Richardson (Romberg/Neville) extrapolation
  of the n-th *forward* finite difference taken along the direction `Scale`,
  `D(h) = (1/(s h)^n) Σ_{k=0}^n (−1)^{n−k} C(n,k) f(x0 + k s h)` with step
  sequence `h_i = s·2^{−i}` and tableau `T(i, j) = T(i, j−1) + (T(i, j−1) −
  T(i−1, j−1))/(2^j − 1)`. The `2^j − 1` (all-powers) denominator — not
  Neville's `4^j − 1` — is correct because a one-sided difference's error
  expansion runs in every power of `h`. The one-sided stencil is what lets a
  complex or directional `Scale` give directional derivatives, and it works for
  non-analytic `f`; it requires integer order `n ≥ 1`. Binomials are exact
  (GMP `mpz_bin_uiui`).
- **`Method -> NIntegrate`**: Cauchy's differentiation formula via the existing
  `NResidue`, `f^{(n)}(x0) = Gamma(n+1)·NResidue[expr/(x−x0)^{n+1}, {x, x0},
  Radius -> Scale]`. `Radius` is always passed explicitly (NResidue's tiny
  default would cause heavy `1/r^n` cancellation), and using `Gamma(n+1)` rather
  than `n!` admits **fractional or complex order**; it needs `expr` analytic
  near `x0`.

**Data structures.** The extrapolation tableau is a `double _Complex*` of
`terms²` on the machine path, or two separate `mpfr_t*` real/imag arrays under
MPFR (complex multiply/divide open-coded, binomials still via GMP). `Scale` may
be complex. `WorkingPrecision` (`MachinePrecision` or a digit count) selects the
path; the NIntegrate method's precision is NResidue's.

**Complexity / limits.** EulerSum does `(n+1)·T` symbolic samples plus an
`O(T²)` tableau. `Terms` (default 7) is the *starting* depth, grown adaptively up
to 50 to meet the goal, keeping the smallest-residual estimate — refinements
that move the estimate by more than `8·best_err` are rejected as
cancellation-driven noise, and halving past the round-off floor amplifies
cancellation (so e.g. a high-order derivative can drift). Options: `Method`
(`EulerSum` | `NIntegrate`, `Automatic` → `EulerSum`), `Scale` (step / contour
radius / complex direction, default 1), `Terms`, `WorkingPrecision`,
`AccuracyGoal` (default `MachinePrecision`), `PrecisionGoal`, `MaxRecursion`
(NIntegrate only). `ND` is not `Listable` (it threads manually to protect the
`{x, n}` spec) and does not chop spurious tiny values; a shortfall warns via
`nc_warn_goal`, and `ND::nnum`/`ivar`/`ord`/`badscl` route through `mth_message`.
