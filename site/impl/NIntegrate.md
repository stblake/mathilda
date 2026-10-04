---
source: src/numerical_calculus/nint.c
references:
  - "R. Piessens, E. de Doncker-Kapenga, C. W. Überhuber and D. K. Kahaner, *QUADPACK: A Subroutine Package for Automatic Integration* (Springer, 1983) — the QAG/QAGS adaptive Gauss–Kronrod strategy."
  - "A. C. Genz and A. A. Malik, *An adaptive algorithm for numerical integration over an N-dimensional rectangular region*, J. Comput. Appl. Math. **6** (1980) 295–302."
  - "H. Takahasi and M. Mori, *Double exponential formulas for numerical integration*, Publ. RIMS Kyoto **9** (1974) 721–741."
  - "T. Ooura and M. Mori, *A robust double exponential formula for Fourier-type integrals*, J. Comput. Appl. Math. **112** (1999) 229–241."
  - "D. Levin, *Fast integration of rapidly oscillatory functions*, J. Comput. Appl. Math. **67** (1996) 95–101."
---
**Algorithm.** `builtin_nintegrate` peels options into an `NiOpts`, validates the
`{x, a, b}` specs, and dispatches. A finite real interval under `Method ->
Automatic` runs a best-of cascade in `ni_core_finite` — adaptive Gauss–Kronrod,
then tanh–sinh, then Levin collocation (if an `f·{cos|sin|e^{i}}` kernel is
detected), then an exponential endpoint-singularity map, then
integration-between-the-zeros — each tried only while no earlier engine has
converged (`ni_consider` keeps the converged result, else the smallest error
estimate). The engines:

- **Gauss–Kronrod** (`gkadapt.c`): the 7-point Gauss rule embedded in the
  15-point Kronrod rule (G7–K15, QUADPACK abscissae), with subintervals in a
  binary max-heap keyed by local error so the worst panel bisects first (QAG);
  running totals feed Wynn's ε (shared `seqaccel`) for endpoint singularities
  (QAGS). Error = `|resk − resg|` refined by the QUADPACK `resasc` heuristic.
- **Double-exponential** (`denint.c`, `dequad.c`): tanh–sinh on a finite
  interval, exp–sinh on a half line (`x = a + exp((π/2)sinh t)`), sinh–sinh on
  the whole line; trapezoidal on the transformed variable, step-halved level by
  level, converging double-exponentially.
- **Oscillatory**: Levin collocation (`levincoll.c`) solves the Levin ODE
  `p' + i g' p = f` in a Chebyshev basis at Chebyshev–Gauss–Lobatto nodes via a
  LAPACK complex LU (accuracy *improves* with oscillation rate; the factored
  matrix is reused across right-hand sides); `oscint.c` locates half-periods by
  sign changes and sums/Wynn-extrapolates one lobe per panel; `oscde.c` is the
  MPFR Ooura–Mori Fourier double-exponential for an aligned `amp·{sin|cos}(ωx)`
  half-line integrand.
- **Multidimensional**: adaptive Genz–Malik degree-7 (embedded degree-5 error)
  box cubature (`cubature.c`, region max-heap), or iterated 1-D quadrature
  (dependent inner bounds via `HoldAll` + `Block` localisation), or a 2-D Levin
  reduction, or quasi-Monte-Carlo (Halton radical-inverse, `mcint.c`) chosen
  automatically for a `Boole`/region integrand or ≥ 6 specs.

Infinite and complex endpoints give exp-sinh rays and piecewise-linear contours
(`ni_run_contour`, each segment `z0 + t(z1−z0)`); `PrincipalValue` uses a
symmetric mirror sample so a simple pole cancels; `Exclusions` are solved with
`Solve` and split the interval so the tanh–sinh rule straddles each singularity.

**Data structures.** Machine path is `double _Complex` throughout, with a
Block-style snapshot/restore of the variable's OwnValues and
`arith_warnings_mute` around each sample. Arbitrary precision is **MPFR**
(`mpfr_t`, run at `target + 64` guard bits and rounded back), with complex
values carried as hand-rolled `mpfr_t` pairs (`qcx`/`ncpx`) — there is no
Arb/MPC dependency; every engine has an MPFR twin. The integrand is
**auto-compiled to bytecode** lazily on the first sample (`autocompile`), with
the interpreter as the exact fallback at any point the compiled real program
reports non-real/singular. Panels (GK) and regions (cubature) live in binary
max-heaps; Levin caches its LU factors.

**Complexity / limits.** GK is adaptive QAG/QAGS; DE converges
double-exponentially (hundreds of digits in a few thousand smooth samples);
Levin is `O(n³)` per solve (`O(n²)` prepared) and cheaper the faster the
oscillation; MC is `O(1/√N)`. Options: `WorkingPrecision`, `AccuracyGoal`
(default `MachinePrecision`), `PrecisionGoal` (default `WorkingPrecision/2`),
`MaxRecursion`, `MinRecursion`, `MaxPoints`, `Exclusions`. An explicitly
requested but unimplemented method (`ClenshawCurtisRule`, `LobattoKronrodRule`,
…) warns `NIntegrate::method` and stays unevaluated rather than silently
approximating; all diagnostics route through `mth_message` (gated on the
arithmetic-warning mute so an inner iterated integral stays quiet).
