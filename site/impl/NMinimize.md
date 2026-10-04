---
source: src/numerical_calculus/nm_driver.c
references:
  - "R. Storn and K. Price, *Differential Evolution — a simple and efficient heuristic for global optimization over continuous spaces*, J. Global Optim. **11** (1997) 341–359."
  - "K. Deb, *An efficient constraint handling method for genetic algorithms*, Comput. Methods Appl. Mech. Engrg. **186** (2000) 311–338 — the feasibility-rule selection."
  - "M. J. D. Powell, *A fast algorithm for nonlinearly constrained optimization calculations*, Lecture Notes in Math. **630** (1978) — the PHR augmented-Lagrangian polish."
---
**Algorithm.** `nm_minimize_driver` (shared with `NMaximize`, registered in
`findmin.c`) parses the objective and variables: bare symbols, `{x, y, …}`
lists, `{x, lo, hi}` search boxes, indexed `x[i]` (rewritten to fresh scalar
symbols so the symbol-keyed solver runs unchanged), and `Element[x, Integers]`
domain declarations. A `{f, c1, c2, …}` first argument is the objective with its
constraints implicitly `And`-ed; `fm_collect_constraints` recurses through
`And`, expands chained inequalities, and routes each leaf to a simple-variable
**box**, a general constraint `g(x) ≤ 0` / `h(x) == 0`, or an `Or` disjunction.
Constraints are enforced two ways at once: a quadratic exterior **penalty**
`Σ max(0, gᵢ)² + Σ hⱼ²` (a disjunction contributes its minimum branch, zero iff
one branch holds), and **Deb's feasibility rules** for *selection* (feasible
beats infeasible; then smaller objective; then smaller violation), with a loose
threshold for search ranking and a tight one for what the caller receives.

The default `Method -> Automatic` is **differential evolution** (`nm_de.c`),
DE/rand/1/bin with Deb-rule selection. Automatic mode strengthens it to
current-to-best/1/bin with a Latin-hypercube initial population (`15n`, clamped),
per-generation dithered `F ∈ [0.5, 1.0]`, `CR = 0.9`, and bounce-back bound
handling; an explicit `"DifferentialEvolution"` keeps the classic `F = 0.6`,
`NP = 10n` path. It stops when the feasible sub-population's objective spread
collapses, then polishes up to `Min[2n, 50]` distinct basins. The other engines
(`nm_neldermead.c`, `nm_randomsearch.c`, `nm_sa.c`, `nm_shgo.c`,
`nm_dual_annealing.c`, `nm_direct.c`, `nm_basin_hopping.c`) implement
restart Nelder–Mead, multi-start local search, geometric-cooling simulated
annealing, SHGO, Tsallis dual annealing, deterministic DIRECT/DIRECT-L, and
basin hopping respectively. The global best is then refined by the **reused
FindMinimum local solver**: plain BFGS (box bounds only) or BFGS inside a PHR
augmented-Lagrangian outer loop when there are general constraints — not SLSQP
or interior-point, which are selectable only for `FindMinimum`. Mixed-integer
problems add integer coordinate descent plus a continuous-relaxation/round/pin
recovery. `NMaximize` builds a synthetic `NMinimize` of `−f` (negating only the
objective element, leaving constraints intact) and negates the reported optimum.

**Data structures.** Population and working buffers are flat row-major `double*`
(`NP×n` in DE, the `(n+1)×n` simplex in Nelder–Mead, chain buffers in SA); the
context `NmDriver` carries the objective program, variables, region bounds,
boxes, general constraints, disjunctions, and integer flags. At
`MachinePrecision` the objective and every general constraint are
**auto-compiled to bytecode** (`compile_expr_ex`, `COMPILE_FOLD_GLOBALS`) and
evaluated by the register machine per trial point, falling back to the
interpreter where a body cannot be lowered. Variables are bound and restored
`Block`-style (`FmVarBind` saves OwnValues + attributes, clears them so the
symbol is a free argument during the solve, and restores on every exit). A
`WorkingPrecision > machine` request refines a continuous, general-constraint-free
result with an MPFR BFGS at the requested bits.

**Complexity / limits.** The RNG is SplitMix64, so the search is **deterministic
for a fixed `RandomSeed`** (default seed fixed); per-attempt seeds derive from
it. Options: `Method` (with per-method sub-options), `MaxIterations` (default
100; setting it disables the Automatic-only best-of-K runs and region growth),
`WorkingPrecision`, `AccuracyGoal`/`PrecisionGoal` (default `WorkingPrecision/2`,
feeding the DE convergence tolerance), `EvaluationMonitor`/`StepMonitor` (fire
inside the local polish, not the global loop). An empty or unmet feasible set
returns `{Infinity, {x -> Indeterminate, …}}`; a success returns
`{fmin, {x -> xmin, …}}` with integer coordinates as exact `Integer`s. Per-trial
evaluation runs under `arith_warnings_mute` and the local solvers under
`g_fm_quiet`, so numeric-domain chatter from probing bad points is suppressed
while driver diagnostics still honour `Quiet[]`/`Check[]`. `NMinimize` is
`Protected` but **not** `HoldAll`.
