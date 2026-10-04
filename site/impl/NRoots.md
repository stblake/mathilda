---
source: src/numerical_roots/nroots.c
references:
  - "D. A. Bini, *Numerical computation of polynomial zeros by means of Aberth's method*, Numer. Algorithms **13** (1996) 179–200."
  - "O. Aberth, *Iteration methods for finding all zeros of a polynomial simultaneously*, Math. Comp. **27** (1973) 339–344."
  - "M. A. Jenkins and J. F. Traub, *Algorithm 419: Zeros of a complex polynomial*, Comm. ACM **15** (1972) 97–99 (CPOLY, ACM TOMS 419)."
---
**Algorithm.** `builtin_nroots` parses `Equal[lhs, rhs]` and the variable, forms
`poly = lhs − rhs`, expands it, extracts coefficients (each numericalised to a
complex MPFR value at working precision), strips a trailing `x^m` factor (exact
zero roots), and dispatches the reduced polynomial. The output is a disjunction
`x == r1 || … || x == rd` (a bare `Equal` for one root), with multiple roots
emitted as repeated identical equations. Three engines:

- **Aberth–Ehrlich** (`nroots_aberth.c`), the default: all `d` roots iterated
  together from **Bini's convex-hull (Newton-polygon) initial placement** — the
  upper hull of `(k, log|a_k|)` seeds starting points on circles — using the
  Aberth correction `w_i = N_i / (1 − N_i S_i)` with `N_i = p(z_i)/p'(z_i)` and
  `S_i = Σ_{j≠i} 1/(z_i − z_j)`, applied Gauss–Seidel in place. Cubically
  convergent at simple roots; a multiplicity-`m` root emerges as a tight
  cluster. Every root is Newton-polished against the original polynomial.
- **CompanionMatrix** (`"CompanionMatrix"`, also the Automatic machine fast
  path): the Frobenius companion's eigenvalues — real coefficients via a real
  QR (LAPACK `dgeev` on the machine path, the MPFR real-QR kernel otherwise,
  "`numpy.roots` exactly"), complex coefficients via the `ℂ^{n×n} → ℝ^{2n×2n}`
  embedding `[[Re, −Im], [Im, Re]]` whose spectrum doubles the roots, filtered
  by residual and read for multiplicity by synthetic division.
- **JenkinsTraub** (`nroots_jt.c`): the three-stage shifted-deflation CPOLY —
  Stage 1 no-shift, Stage 2 fixed-shift near a Cauchy-circle zero (sweeping up
  to 12 angles), Stage 3 variable-shift (cubically convergent) — the found zero
  polished against the original polynomial, then deflated; a multiplicity-`m`
  zero is found `m` times.

An exact-integer polynomial of degree ≥ 2 is first **Yun squarefree-
decomposed** so a high multiplicity like `(x²−2)^30` stays well conditioned.

**Data structures.** `NrPoly { int deg; ncpx* c; int is_real; mpfr_prec_t prec }`,
where `ncpx` is the MPFR complex type (a pair of `mpfr_t`). An Automatic machine
request works in `double`/`double _Complex` (real/interleaved companion arrays
for LAPACK), rounding the final roots back to `double`/`Complex`; a
`PrecisionGoal` digit count runs the *same* Aberth loop entirely in MPFR complex
arithmetic at `target_bits + max(48, target_bits/2)` bits (the convergence
tolerance `2^{−(wp−8)}` scales with it). Post-processing chops noise, single-
linkage clusters near-equal roots to a shared value (so multiplicities become
numerically identical), symmetrises conjugate pairs for real coefficients, and
canonically sorts. The whole numeric core is `#ifdef USE_MPFR`; without MPFR,
`NRoots` emits `nompfr` and stays unevaluated.

**Complexity / limits.** Aberth is `O(d²)` per sweep over `MaxIterations`
(default `100 + 20d`); companion is the `O(d³)` QR; Jenkins–Traub is one root at
a time with fixed inner caps (`MaxIterations` is ignored by companion and JT).
`PrecisionGoal` (`Automatic`/`Infinity` → machine; a digit count → MPFR) selects
precision; `AccuracyGoal` (default `MachinePrecision`) drives only the post-solve
residual check — a root whose Newton correction exceeds the combined tolerance
triggers `NRoots::accgl` (clustered/stationary roots are skipped in that test).
Declines (stay unevaluated) on a non-equation, non-polynomial, non-numeric
coefficient, non-variable, unknown method, or engine non-convergence (`conv`);
all diagnostics route through `mth_message`. `NRoots` is `Protected`, not
`Listable`.
