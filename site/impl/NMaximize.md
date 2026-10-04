---
source: src/numerical_calculus/nm_driver.c
references:
  - "R. Storn and K. Price, *Differential Evolution — a simple and efficient heuristic for global optimization over continuous spaces*, J. Global Optim. **11** (1997) 341–359."
  - "K. Deb, *An efficient constraint handling method for genetic algorithms*, Comput. Methods Appl. Mech. Engrg. **186** (2000) 311–338."
---
**Algorithm.** `builtin_nmaximize` is a thin wrapper over the `NMinimize`
driver. It constructs a synthetic `NMinimize` call that **minimises `−f`**: for a
bare objective it negates `f`; for a `{f, c1, c2, …}` list it negates only the
objective element (`Times[-1, f]`) and copies every constraint unchanged — a
fix for an earlier bug where negating the whole list threaded `Times[-1, …]` over
the constraints and corrupted them. It then calls `nm_minimize_driver(...,
"NMaximize")` and **negates the reported optimum** back (`mpfr_neg` or `-fmin`,
preserving the numeric type). Everything else — the global-method catalogue
(differential evolution by default, with Deb's feasibility rules), the variable
and constraint grammar, integer domains, the BFGS/augmented-Lagrangian local
polish, auto-compilation at `MachinePrecision`, and fixed-seed determinism — is
`NMinimize`'s, unchanged. See the `NMinimize` implementation notes for the full
account.

**Data structures.** Identical to `NMinimize`: the synthetic objective is
compiled to bytecode and evaluated per trial point; variables are bound and
restored `Block`-style; populations and simplices are flat `double*` buffers; a
`WorkingPrecision > machine` request uses the MPFR BFGS refinement.

**Complexity / limits.** As for `NMinimize`. Returns `{fmax, {x -> xmax, …}}`;
an empty or unmet feasible set returns `{-Infinity, {x -> Indeterminate, …}}`
(the negation of `NMinimize`'s infeasible `Infinity`). `NMaximize` is `Protected`
but **not** `HoldAll`, and shares `NMinimize`'s options
(`Method`, `MaxIterations`, `WorkingPrecision`, `AccuracyGoal`, `PrecisionGoal`,
`EvaluationMonitor`, `StepMonitor`) and message routing.
