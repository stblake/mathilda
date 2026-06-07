---
source: src/fit.c
---
**Algorithm.** `builtin_designmatrix` (in `fit.c`, the same module as `Fit`) builds the design matrix `m_{ij} = f_i(coords_j)` for `DesignMatrix[data, {f1,…,fn}, vars]`. It normalises the data shape (`{v1,…}` → `{{1,v1},…}`; `{{x,v},…}` univariate; `{{x1,…,xk,v},…}` multivariate, dropping the trailing response value to obtain the coordinate vector for each row), then evaluates each basis function `f_i` at each data point by substituting the `vars` symbols with that row's coordinates (pattern/`ReplaceAll`-style substitution followed by `evaluate`). The result is a `List` of rows, one per data point, each the vector `{f_1(coords), …, f_n(coords)}` — exactly the matrix that `Fit`/`LeastSquares` solve against.

**Data structures.** A `List` of `List`s of evaluated basis-function values; entry representation (exact vs. machine vs. MPFR) follows the same `WorkingPrecision` handling as `Fit`. This is the Vandermonde-like basis matrix for a polynomial basis but works for any list of basis functions.
