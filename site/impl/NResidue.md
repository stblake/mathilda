---
source: src/numerical_calculus/nresidue.c
references:
  - "L. N. Trefethen and J. A. C. Weideman, *The exponentially convergent trapezoidal rule*, SIAM Rev. **56** (2014) 385–458 — periodic-trapezoidal contour integration."
---
**Algorithm.** `builtin_nresidue` evaluates `NResidue[expr, {z, z0}]`, the
residue at `z = z0` (the coefficient of `(z − z0)^{−1}`), by the **periodic-
trapezoidal Cauchy integral** `(1/2πi) ∮ f dz` around a small circle — the
trapezoidal rule on a circle, which converges exponentially for an analytic
integrand. Unlike the symbolic `Residue` (which needs a Laurent series), this
works for essential singularities (`Exp[1/x]`, `Sin[1/x]`). It threads manually
over a `List` first argument, and each sample binds `z` Block-style and bumps
the evaluation clock. `NResidue` supplies only the sampler callbacks and option
plumbing; the actual trapezoidal evaluation, N-doubling refinement, and
branch-cut detection live in the shared quadrature layer
(`quadrature.c`/`.h`, `qd_contour_residue_machine`/`_mpfr`), which also provides
the `Radius -> Automatic` adaptive-radius walk.

**Data structures.** The machine path uses a `double _Complex` sampler callback,
the MPFR path a pair of `mpfr_t` real/imag callbacks; the trapezoidal sample
buffers (and the Aitken/Shanks Δ² extrapolation and radius search) belong to the
quadrature layer, not this file. `NumericSpec` selects precision;
`WorkingPrecision` above machine routes to the MPFR quadrature twin.

**Complexity / limits.** Trapezoidal sampling with geometric N-doubling bounded
by `MaxRecursion` (default 10 refinements), each doubling the full-symbolic
sample count. `Radius` defaults to `1/100` (or `Automatic`). The quadrature
status is reported through `nr_report`: a non-convergent result is returned
without warning when it already meets the absolute `AccuracyGoal` floor, else
warns `NResidue::ncvi`; a detected branch cut on the contour warns
`NResidue::bcut` (the result is unreliable — the contour must enclose only the
one singularity and cross no cut); a non-numeric sample warns `NResidue::nnum`
and returns `NULL`. NResidue cannot distinguish a tiny spurious residual from a
true zero — `Chop` when needed. Options: `Radius`, `WorkingPrecision`,
`AccuracyGoal` (default `MachinePrecision`), `PrecisionGoal` (unset → two guard
digits below working precision), `MaxRecursion`, `Method` (`"Trapezoidal"`, the
only one). All diagnostics route through `mth_message`.
