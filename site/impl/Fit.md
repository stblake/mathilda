---
source: src/fit.c
---
**Algorithm.** `builtin_fit` fits a linear combination `a_1 f_1 + … + a_n f_n` of basis functions to data and returns the symbolic fit expression. `Fit[data, {f1,…,fn}, vars]` first forms the design matrix `m` (the same construction as `DesignMatrix`) and response vector `v`, solves the linear least-squares problem for the coefficient vector `a`, and reassembles `Σ a_i f_i`. `Fit[{m, v}]` solves directly given a design matrix and response. The least-squares solve is reuse-first:

- Plain L2 and **ridge/Tikhonov** (`FitRegularization -> {"L2"|"Tikhonov"|"RidgeRegression", λ}`) route through the `LeastSquares` builtin (`PseudoInverse . v`); ridge is reduced to ordinary least squares on the augmented system `[m; √λ I]`, `[v; 0]`.
- **LASSO** (`{"L1"|"LASSO", λ}`) uses cyclic coordinate descent with soft-thresholding (machine precision).
- `NormFunction -> Norm[#,1]` (least absolute deviations) uses iteratively reweighted least squares (IRLS).
- Any other norm, or a norm combined with regularisation, falls back to `FindMinimum`, warm-started from the L2 solution.

`WorkingPrecision` selects exact rational (`Infinity`), machine (`Automatic`), or `n`-digit MPFR arithmetic.

**Data structures.** Design matrix and response are `List`s of `List`s / `List`; coefficients come back as a vector that is recombined with the basis-function expressions. Data shapes are normalised exactly as in `DesignMatrix` (implicit `{1,2,…}` abscissae, univariate, or multivariate coordinate rows).
