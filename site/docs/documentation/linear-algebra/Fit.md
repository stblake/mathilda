# Fit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Fit[data, {f1, ..., fn}, vars] finds a least-squares fit a1 f1 + ... + an fn to data for functions of the variables vars (a symbol or list of symbols).
Fit[{m, v}] gives the coefficient vector minimizing ||m.a - v|| for design matrix m and response vector v.
Data may be a list of values {v1, ...} (coordinates 1, 2, ...), univariate pairs {{x, v}, ...}, or multivariate rows {{x, y, ..., v}, ...}.
Options: WorkingPrecision (Automatic | n | Infinity), FitRegularization ({"Tikhonov"|"L2"|"RidgeRegression"|"LASSO"|"L1", lambda}), NormFunction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= 0.186441 + 0.694915 x

In[2]:= Fit[{{0,1},{1,0},{3,2},{5,4}}, {1, x, x^2}, x, WorkingPrecision -> Infinity]
Out[2]= 135/199 - 53/199 x + 38/199 x^2

In[3]:= Fit[{N[HilbertMatrix[4]], Range[4]}]
Out[3]= {-64.0, 900.0, -2520.0, 1820.0}

In[4]:= Fit[{{0,0,0},{1,0,1},{0,1,2},{1,1,0},{1/2,1/2,1}}, {1, x, y}, {x, y}]
Out[4]= 0.8 - 0.5 x + 0.5 y
```

## Implementation notes

**Algorithm.** `builtin_fit` fits a linear combination `a_1 f_1 + … + a_n f_n` of basis functions to data and returns the symbolic fit expression. `Fit[data, {f1,…,fn}, vars]` first forms the design matrix `m` (the same construction as `DesignMatrix`) and response vector `v`, solves the linear least-squares problem for the coefficient vector `a`, and reassembles `Σ a_i f_i`. `Fit[{m, v}]` solves directly given a design matrix and response. The least-squares solve is reuse-first:

- Plain L2 and **ridge/Tikhonov** (`FitRegularization -> {"L2"|"Tikhonov"|"RidgeRegression", λ}`) route through the `LeastSquares` builtin (`PseudoInverse . v`); ridge is reduced to ordinary least squares on the augmented system `[m; √λ I]`, `[v; 0]`.
- **LASSO** (`{"L1"|"LASSO", λ}`) uses cyclic coordinate descent with soft-thresholding (machine precision).
- `NormFunction -> Norm[#,1]` (least absolute deviations) uses iteratively reweighted least squares (IRLS).
- Any other norm, or a norm combined with regularisation, falls back to `FindMinimum`, warm-started from the L2 solution.

`WorkingPrecision` selects exact rational (`Infinity`), machine (`Automatic`), or `n`-digit MPFR arithmetic.

**Data structures.** Design matrix and response are `List`s of `List`s / `List`; coefficients come back as a vector that is recombined with the basis-function expressions. Data shapes are normalised exactly as in `DesignMatrix` (implicit `{1,2,…}` abscissae, univariate, or multivariate coordinate rows).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/fit.c`](https://github.com/stblake/mathilda/blob/main/src/fit.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

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

### Notes

`Fit[data, {f1, ..., fn}, x]` returns the least-squares linear combination
`a1 f1 + ... + an fn` of the basis functions. Plain `{v1, v2, ...}` data is
taken at abscissae `1, 2, ...`, while `{{x, v}, ...}` pairs supply explicit
abscissae; the perfect-square data recovers `x^2` exactly. The third example
fits a quadratic trend to noisy data, and the last shows a multivariate fit
in two predictors `{x, y}` from `{x, y, v}` rows — the design matrix is
solved by normal equations in either case.
