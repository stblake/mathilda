---
references:
  - "J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006)."
source: src/findmin.c
---
**Algorithm.** `FindMaximum` (`HoldAll | Protected`) is a thin wrapper over
`FindMinimum` (src/findmin.c, `builtin_findmaximum`): it negates the objective,
runs the same local optimizer, and negates the first component of the resulting
`{f_min, {x -> x_min, ...}}` pair to report `{f_max, {x -> x_max, ...}}`. All
machinery — Brent in 1-D, BFGS quasi-Newton / conjugate-gradient / Newton in
n-D, symbolic gradients/Hessian with a central-difference fallback, Armijo line
search, box-projection and quadratic-penalty constraint handling, MPFR extended
precision — is inherited unchanged from `FindMinimum`. See `FindMinimum` for the
full description.

**Complexity / limits.** Same as `FindMinimum`: local search only. The negation
is precision-aware (`mpfr_neg` for `EXPR_MPFR` results, plain real otherwise).
