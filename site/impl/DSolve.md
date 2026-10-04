---
references:
  - "E. S. Cheb-Terrab and A. D. Roche, *Integrating factors for second-order ODEs* (1999) — the ReducibleIntegratingFactor / ReducibleFirstIntegral classes, cited in the dispatcher."
source: src/calculus/dsolve.c
---
**Algorithm.** `builtin_dsolve` is a cascade polyalgorithm that mirrors
`Integrate`: a `Method`-option enum (`ds_method_from_string`) selects either the
automatic cascade or a single pinned method (strict, no fallback). After
`dsolve_parse` builds a shared `DSolveProblem`, the dispatcher branches on the
problem shape — PDE, a system (`nfun > 1`), or a scalar ODE — and in each branch
tries methods in a fixed order with `if (!result) result = dsolve_run(&P,
<method>_try)` until one succeeds. The scalar `DS_AUTOMATIC` chain runs the
specialists roughly front-to-back by specificity: Factorable and NthAlgebraic
split products/powers of the top derivative first; then first-order named classes
(Quadrature, LinearFirstOrder, Bernoulli, Homogeneous, Separable, Exact,
Clairaut, Lagrange); then constant-coefficient and Euler–Cauchy linear ODEs
(UndeterminedCoefficients, LinearConstantCoefficients); then the heavier
second-order machinery (Kovacic, special-function / change-of-variable,
variation-of-parameters); then substitution and reduction methods (Riccati,
Chini, Abel, Lie point symmetry); and finally the always-available Frobenius /
power-series fallbacks. Each method file (`src/calculus/dsolve_<method>.c`)
returns an array of solution branches or `NULL` to fall through.

**Data structures.** The problem substrate (parse / verify / fit / assemble)
lives in `dsolve_common.c`; the dispatcher only sequences method `*_try`
functions and wraps them in `dsolve_run`, `dsolve_run_implicit`,
`dsolve_run_parametric`, `dsolve_run_first_integral`, `dsolve_run_system` and
`dsolve_run_pde`. A per-command fail-memo keyed on `eval_toplevel_id()`
(`ds_fail_tab`, 32 slots) records each `(equation, variable, method)` the
deterministic cascade already declined, so the fixed-point loop does not re-run
the whole cascade on re-entry; `g_dsolve_depth` separates the outermost user
call from internal recursions. For the whole cascade the cosmetic
`Power::infy` / `Infinity::indet` warnings are muted
(`arith_warnings_mute_push`), because speculative probes legitimately form `1/0`
while classifying an equation — the back-substitution verifier is the real gate.

**Complexity / limits.** Cost is the sum of the declining probes up to the first
method that claims the equation, so ordering is tuned to reach each family
before an earlier general method spins on it (several methods carry per-attempt
`TimeConstrained` deadlines and decline memos for exactly this reason). Options
default to `GeneratedParameters -> C`, `Assumptions -> True`,
`Method -> Automatic`, `IncludeSingularSolutions -> False`. Every returned branch
is verified by back-substitution before it is kept. A few methods are
pinned-only (not in the automatic chain) — `FirstOrderPowerSeries`, `SolvableForX`,
the eigenvalue problem — so by default a first-order ODE with no closed form
stays unevaluated rather than returning a truncated series.
