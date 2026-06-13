# ND — Numerical Derivative builtin

`ND[expr, x, x0]` / `ND[expr, {x, n}, x0]` — numerical approximation to the
(n-th) derivative of `expr` w.r.t. `x` at `x = x0`. Mirrors the legacy
Mathematica `NumericalCalculus`ND`.

## Design summary

Two methods (reuse `src/numerical_calculus/quadrature.c`; new module
`src/numerical_calculus/nderiv.{c,h}`):

- **Method -> EulerSum (DEFAULT).** Richardson extrapolation of *forward*
  (one-sided, direction-`s`) finite differences. The n-th forward difference
  `Delta^n_h f / (s h)^n` with `Delta^n_h f(x0) = sum_{k=0}^n (-1)^{n-k} C(n,k)
  f(x0 + k s h)`, sampled on the geometric ladder `h_i = h0 / 2^i` (h0 = Scale),
  extrapolated h->0 with the Romberg/Neville tableau
  `T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1)`.
  Forward (not central) stencil is what reproduces the documented directional
  results, e.g. `ND[Abs[x],{x,1},0,Scale->1+I] = 0.707107 - 0.707107 I`.
  Requires integer `n >= 1`. Works for non-analytic `f` (only samples along `s`).

- **Method -> NIntegrate.** Cauchy integral formula:
  `f^(n)(x0) = n! * Res_{z=x0}[ f(z) / (z-x0)^(n+1) ]`. **Reuse the existing
  NResidue builtin wholesale** — evaluate
  `Gamma(n+1) * NResidue[expr/(x-x0)^(n+1), {x, x0}, Radius->Scale, <mapped opts>]`.
  Crucial: pass `Radius -> Scale` (default 1) to override NResidue's tiny 1/100
  default radius — that small radius (not the method) is what makes the doc's
  `4!NResidue[...]` example lose ~9 digits; at radius 1, `1/r^(n+1)=1` so the
  pointwise `(z-x0)^(n+1)` division is well-conditioned and recovers the
  documented ~1e-16 accuracy. Forming `f(z)/(z-x0)^(n+1)` pointwise is exactly
  the direct Cauchy `e^{-i n theta}` weight on the contour, so no accuracy is
  lost. `n!` -> `Gamma(n+1)` so fractional/complex order works (e.g.
  `ND[x,{x,-1/2},1,Method->NIntegrate] = 4/(3 Sqrt[Pi]) ~ 0.752253`); accept a
  QD_BRANCHCUT warning for fractional order (best-effort). Requires `f` analytic
  near `x0` (gives wrong answers on `Re`/`Im`/`Abs`/...). Reuses all of
  NResidue: option parsing, machine/MPFR dispatch, Block binding, reporting.

Options (trailing `Rule`s, NResidue-style peeling): `Method` (EulerSum |
NIntegrate, default EulerSum), `Scale` (default 1; step for EulerSum / radius
for NIntegrate; may be complex), `Terms` (default 7; EulerSum tableau depth),
`WorkingPrecision` (MachinePrecision | digits -> MPFR), `PrecisionGoal`,
`MaxRecursion` (NIntegrate doublings).

Precision: machine `double _Complex` and MPFR complex paths, exactly like
NResidue. Variable bound Block-style via temporary OwnValues; `numericalize`
the integrand at each sample. Manual list-threading on arg0 (like Series),
`ATTR_PROTECTED` only (NO `ATTR_LISTABLE` — would split the `{x,n}` spec).

## Tasks

- [ ] Add interned symbols `ND`, `Scale`, `Terms`, `EulerSum`, `NIntegrate` in
      `src/sym_names.{c,h}`.
- [ ] `src/numerical_calculus/nderiv.h` — declare `builtin_nd`, `nd_init`.
- [ ] `src/numerical_calculus/nderiv.c`:
      - arg/spec parse (`x` or `{x, n}`), x0 numeric-eval, option peeling.
      - Block-style binding + sampler helpers (mirror nresidue.c).
      - EulerSum: machine + MPFR forward-difference Richardson tableau (file-local
        complex toolkit, numeric_complex.c style; `(s h)^n` via complex pow).
      - NIntegrate: build `Gamma(n+1) * NResidue[expr/(x-x0)^(n+1), {x,x0},
        Radius->Scale, <mapped opts>]` and evaluate (reuse NResidue; no new
        sampler / quadrature wiring needed in nderiv.c).
      - status reporting + Expr result (NULL-on-fail; never free `res`).
      - manual threading over a List in arg0.
- [ ] Register in `src/core.c` (forward decl + `nd_init()` call).
- [ ] Docstring in `src/info.c` (terse, no examples).
- [ ] `tests/CMakeLists.txt`: add `nderiv.c` to COMMON_SRC + `nderiv_tests`
      executable.
- [ ] `tests/test_nderiv.c`: cover the documented examples as oracle
      (Exp/Cos^3/Sin@(Pi I)/list-threading/Scale/Terms/WorkingPrecision/
      NIntegrate/fractional/non-analytic-uses-EulerSum), plus error paths.
- [ ] Docs: `docs/spec/builtins/` (numerical-calculus / calculus) + changelog
      `docs/spec/changelog/2026-06-08.md`; overview row if needed.
- [ ] Build clean (`-std=c99 -Wall -Wextra`), run `nderiv_tests` + `nresidue_tests`,
      valgrind a representative call (diff against Sin[1.0] baseline noise).

## Review

Done. ND implemented in `src/numerical_calculus/nderiv.{c,h}` (~640 lines),
wired into `core.c`, `info.c`, `sym_names.c`, and the test CMake.

- **EulerSum (default)**: forward-difference Richardson tableau, machine
  `double _Complex` + MPFR complex paths. Reproduced every documented example
  exactly (Exp/Cos^3/Sin@PiI/Cos[Ix]@1+I, list threading, non-analytic
  Re[Cos[I y]], Abs directional incl. Scale->1+I, Sin[100x] Scale/Terms, the
  WorkingPrecision->40 cases). The documented instability cases (default
  `ND[Exp[x],{x,10},0]`, `Terms->20`) reproduce the documented "nonsense"
  behaviour, with higher precision recovering the right value.
- **NIntegrate**: reuses NResidue as planned — `Gamma(n+1) *
  NResidue[expr/(x-x0)^(n+1), {x,x0}, Radius->Scale]`. Matched the documented
  analytic, fractional (`x^{-1/2}` = 4/(3√π)), and complex-order (`x^4` at `I`)
  results to displayed precision; `Terms->20,WorkingPrecision->40` matched the
  doc's 23-digit value.
- 18 unit tests in `tests/test_nderiv.c`, all passing; `nresidue_tests` still
  green. Compiles clean under `-std=c99 -Wall -Wextra`. Valgrind: no Mathilda
  frames in any leak stack; `definitely lost` equals the macOS baseline noise
  (12,800 B / 400 blocks), i.e. no new leaks.

Notes worth remembering: the EulerSum stencil is *forward* (not central) — that
is what reproduces directional/one-sided derivatives; and the test harness
compares inside the language (`N[Abs[Re d] + Abs[Im d]] < tol`) because the
printer rounds machine reals to ~6 sig figs and `Abs` of a machine complex
leaves a spurious tiny imaginary residual.
