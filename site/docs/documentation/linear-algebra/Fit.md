# Fit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Fit[data, {f1, ..., fn}, vars] finds a least-squares fit a1 f1 + ... + an fn to data for functions of the variables vars (a symbol or list of symbols).`**

**`Fit[{m, v}] gives the coefficient vector minimizing ||m.a - v|| for design matrix m and response vector v.`**

<details>
<summary>Notes</summary>

Data may be a list of values {v1, ...} (coordinates 1, 2, ...), univariate pairs {{x, v}, ...}, or multivariate rows {{x, y, ..., v}, ...}. Options: WorkingPrecision (Automatic | n | Infinity), FitRegularization ({"Tikhonov"|"L2"|"RidgeRegression"|"LASSO"|"L1", lambda}), NormFunction.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= 0.186441 + 0.694915 x

In[2]:= Fit[{N[HilbertMatrix[4]], Range[4]}]
Out[2]= {-64.0, 900.0, -2520.0, 1820.0}

In[3]:= Fit[{{0,0,0},{1,0,1},{0,1,2},{1,1,0},{1/2,1/2,1}}, {1, x, y}, {x, y}]
Out[3]= 0.8 - 0.5 x + 0.5 y
```

### Options (3)

```mathematica
In[4]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x, x^2}, x, WorkingPrecision -> Infinity]
Out[4]= 135/199 - 53/199 x + 38/199 x^2

In[5]:= Fit[{{0.,0.},{0.001,1},{0.01,1}}, {1, x, x^2}, x, FitRegularization -> {"Tikhonov", 1}]
Out[5]= 0.499985 + 0.00549961 x + 5.0496e-05 x^2

In[6]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x, NormFunction -> Function[Norm[#, 1]]]
Out[6]= -1.0 + 1.0 x
```

### Applications (4)

```mathematica
In[1]:= Fit[{1, 2, 1.3, 3.75, 2.25}, {1, x}, x]
Out[1]= 0.785 + 0.425 x
```

```mathematica
In[1]:= Fit[{1, 4, 9, 16}, {1, x, x^2}, x]
Out[1]= 0.0 + 0.0 x + 1.0 x^2
```

```mathematica
In[1]:= Fit[{{0, 1}, {1, 2.7}, {2, 7.4}, {3, 20.1}}, {1, x, x^2}, x]
Out[1]= 1.25 - 2.05 x + 2.75 x^2
```

```mathematica
In[1]:= Fit[{{0, 0, 1}, {1, 0, 2}, {0, 1, 3}, {1, 1, 5}}, {1, x, y}, {x, y}]
Out[1]= 0.75 + 1.5 x + 2.5 y
```

## Algorithm

fit.c — linear least-squares regression (Fit) and DesignMatrix.

This module implements Mathematica's `Fit` builtin: it fits a linear

```text
combination  a1 f1 + ... + an fn  of basis functions to data, plus the
```

companion `DesignMatrix` (the matrix of basis functions evaluated at the data coordinates).

Call forms ----------

```text
  Fit[data, {f1,...,fn}, vars]
      Fits a1 f1 + ... + an fn to `data`.  `vars` is a single symbol `x`
      or a list {x, y, ...}.  Returns the symbolic fit expression.
  Fit[{m, v}]
      Given a design matrix `m` and response vector `v`, returns the
      coefficient vector a minimising ||m.a - v||.
  DesignMatrix[data, {f1,...,fn}, vars]
      Returns the design matrix m_ij = f_i(coords_j).
```

Data shapes (3-argument form) -----------------------------

```text
  {v1,...,vn}            equivalent to {{1,v1},...,{n,vn}}.
  {{x1,v1},...}          univariate: coordinate x_i, response v_i.
  {{x1,...,xk,v1},...}   multivariate: leading k coordinates, last value.
```

Options -------

```text
  WorkingPrecision -> Automatic  (default; exact input -> machine reals)
                    -> n         (n-digit MPFR arithmetic)
                    -> Infinity  (exact rational arithmetic)
  FitRegularization -> {"Tikhonov"|"L2"|"RidgeRegression", lambda}
                         minimise ||m.a-v||^2 + lambda ||a||^2 (ridge).
                    -> {"LASSO"|"L1", lambda}
                         minimise ||m.a-v||^2 + lambda ||a||_1.
  NormFunction -> Function[Norm[#,p]]
                    minimise normf[m.a - v] instead of the 2-norm.
```

Solvers (reuse-first) ---------------------

```text
  * Plain L2 and ridge route through the existing LeastSquares builtin,
    which already supports exact (rational), machine (Real) and MPFR
    arithmetic.  Ridge is reduced to ordinary least squares on the
    augmented system [m; sqrt(lambda) I] / [v; 0].
  * LASSO uses cyclic coordinate descent with soft-thresholding
    (machine precision).
  * NormFunction -> Norm[#,1] (least absolute deviations) uses iteratively
    reweighted least squares (IRLS, machine precision).
  * Any other norm (or a norm combined with regularisation) falls back to
    the FindMinimum builtin, warm-started from the L2 solution.
```

Memory ownership follows the standard builtin contract: the evaluator owns `res` and frees it; this file never frees `res` or its argument

```text
subtrees.  Every intermediate built here is freed on every return path.
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Fit degree-8 polynomial, 5000 points | 0.999 s | 1.29 s | 0.327 s |
| LeastSquares 2000x20 | 0.496 s | 0.612 s | 0.293 s |
| Fit trig basis, 5000 points | 0.224 s | 0.409 s | 0.092 s |
| Fit quadratic, 5000 points | 0.223 s | 0.413 s | 0.141 s |
| FindMinimum Rosenbrock 2-D | 0.19 s | 0.158 s | 3.57 s |
| Fit linear, 5000 points | 0.163 s | 0.274 s | 0.099 s |

## Implementation notes

**Algorithm.** `builtin_fit` fits a linear combination `a_1 f_1 + … + a_n f_n` of basis functions to data and returns the symbolic fit expression. `Fit[data, {f1,…,fn}, vars]` first forms the design matrix `m` (the same construction as `DesignMatrix`) and response vector `v`, solves the linear least-squares problem for the coefficient vector `a`, and reassembles `Σ a_i f_i`. `Fit[{m, v}]` solves directly given a design matrix and response. The least-squares solve is reuse-first:

- Plain L2 and **ridge/Tikhonov** (`FitRegularization -> {"L2"|"Tikhonov"|"RidgeRegression", λ}`) route through the `LeastSquares` builtin (`PseudoInverse . v`); ridge is reduced to ordinary least squares on the augmented system `[m; √λ I]`, `[v; 0]`.
- **LASSO** (`{"L1"|"LASSO", λ}`) uses cyclic coordinate descent with soft-thresholding (machine precision).
- `NormFunction -> Norm[#,1]` (least absolute deviations) uses iteratively reweighted least squares (IRLS).
- Any other norm, or a norm combined with regularisation, falls back to `FindMinimum`, warm-started from the L2 solution.

`WorkingPrecision` selects exact rational (`Infinity`), machine (`Automatic`), or `n`-digit MPFR arithmetic.

**Data structures.** Design matrix and response are `List`s of `List`s / `List`; coefficients come back as a vector that is recombined with the basis-function expressions. Data shapes are normalised exactly as in `DesignMatrix` (implicit `{1,2,…}` abscissae, univariate, or multivariate coordinate rows).

**Attributes:** `Protected`.

## See also

[FindMinimum](../../calculus/FindMinimum/), [LeastSquares](../../linear-algebra/LeastSquares/), [NDArray](../../linear-algebra/NDArray/)

## References

- Source: [`src/fit.c`](https://github.com/stblake/mathilda/blob/main/src/fit.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_fit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fit.c)

## Notes & additional examples

### Notes

`Fit[data, {f1, ..., fn}, x]` returns the least-squares linear combination
`a1 f1 + ... + an fn` of the basis functions. Plain `{v1, v2, ...}` data is
taken at abscissae `1, 2, ...`, while `{{x, v}, ...}` pairs supply explicit
abscissae; the perfect-square data recovers `x^2` exactly. The third example
fits a quadratic trend to noisy data, and the last shows a multivariate fit
in two predictors `{x, y}` from `{x, y, v}` rows — the design matrix is
solved by normal equations in either case.
