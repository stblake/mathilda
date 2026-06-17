---
references:
  - "W. H. Press et al., *Numerical Recipes*, 3rd ed. (Cambridge, 2007) — Newton, secant, Brent's method."
source: src/numerical_roots/findroot.c
---
**Algorithm.** `FindRoot` (`HoldAll | Protected`) does iterative numerical
root-finding (src/numerical_roots/findroot.c). It Block-binds the search variables via temporary
OwnValues so the global symbol table is unperturbed, restoring them on every
exit path. The variable spec selects the method:

- `{x, x0}` → **Newton** from one start. The derivative is obtained
  symbolically by evaluating `D[f, x]` (`fr_compute_derivative`) and re-evaluated
  numerically each step; the update is `x -= damping·f/f'`. Real and complex
  variants exist (`fr_run_newton_real`, `fr_run_newton_complex`).
- `{x, x0, x1}` → **secant** (`fr_run_secant_real`), used as the derivative-free
  fallback.
- `{x, xstart, xmin, xmax}` → **Brent** bracketing (`fr_run_brent_real`).
- `{{x, x0}, {y, y0}, ...}` with a list of equations → **multivariate Newton**
  (`fr_run_newton_system_real`): the Jacobian `J[i][j] = D[f_i, var_j]` is
  precomputed symbolically (or taken from a user `Jacobian` option),
  re-evaluated each step, and the Newton step solves `J·dx = f` by **Gaussian
  elimination with partial pivoting** (inline LU-style forward elimination +
  back substitution), then `x -= damping·dx`.

Equation forms `lhs == rhs` are normalised to `lhs - rhs`. Options include
`Method` (Automatic/"Newton"/"Secant"/"Brent"), `WorkingPrecision`
(MachinePrecision or an MPFR bit count), `MaxIterations` (default 100),
`AccuracyGoal`/`PrecisionGoal`, `DampingFactor`, `Jacobian`, and held
`StepMonitor`/`EvaluationMonitor`. When `WorkingPrecision` requests extended
precision and MPFR is built in, dedicated MPFR Newton paths
(`fr_run_newton_mpfr_real`, …) run the iteration at the requested bit width.

**Data structures.** `double` / C99 `double complex` for machine precision;
`mpfr_t` scalars for extended precision; a flat row-major `double` augmented
matrix for the system solve. Function/derivative evaluation re-enters the
Mathilda evaluator with the current numeric binding installed.

**Complexity / limits.** Per-iteration cost is dominated by symbolic-expression
evaluation (and `O(n^3)` Gaussian elimination for an n-variable system).
Convergence is local — quadratic for Newton near a simple root. Emits
diagnostics (`dsing`, `noconv`, `nlnum`) and returns NULL/unevaluated when the
derivative vanishes, the Jacobian is singular, evaluation is non-numeric, or it
fails to converge.
