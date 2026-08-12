# DesignMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DesignMatrix[data, {f1, ..., fn}, vars] gives the design matrix with entries f_i evaluated at the data coordinates.`**

<details>
<summary>Notes</summary>

Data shapes match Fit. The WorkingPrecision option converts entries to machine or n-digit reals; otherwise they are exact.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= DesignMatrix[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= {{1, 0}, {1, 1}, {1, 3}, {1, 5}}

In[2]:= DesignMatrix[{{0,0,0},{1,0,1},{0,1,2}}, {1, x, y}, {x, y}]
Out[2]= {{1, 0, 0}, {1, 1, 0}, {1, 0, 1}}
```

### Applications (4)

```mathematica
In[1]:= DesignMatrix[{{0, 1}, {1, 0}, {3, 2}, {5, 4}}, {1, x}, x]
Out[1]= {{1, 0}, {1, 1}, {1, 3}, {1, 5}}
```

```mathematica
In[1]:= DesignMatrix[{{1, 1}, {2, 8}, {3, 27}}, {1, x, x^2, x^3}, x]
Out[1]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}}
```

```mathematica
In[1]:= DesignMatrix[{{1, 1, 5}, {2, 4, 6}, {3, 9, 2}}, {1, x, y, x*y}, {x, y}]
Out[1]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}}
```

```mathematica
In[1]:= DesignMatrix[{{1, 2}, {2, 5}}, {1, Sin[x]}, x, WorkingPrecision -> 40]
Out[1]= {{1.0, 0.84147098480789650665250232163029899962254}, {1.0, 0.90929742682568169539601986591174484270222}}
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

## Implementation notes

**Algorithm.** `builtin_designmatrix` (in `fit.c`, the same module as `Fit`) builds the design matrix `m_{ij} = f_i(coords_j)` for `DesignMatrix[data, {f1,…,fn}, vars]`. It normalises the data shape (`{v1,…}` → `{{1,v1},…}`; `{{x,v},…}` univariate; `{{x1,…,xk,v},…}` multivariate, dropping the trailing response value to obtain the coordinate vector for each row), then evaluates each basis function `f_i` at each data point by substituting the `vars` symbols with that row's coordinates (pattern/`ReplaceAll`-style substitution followed by `evaluate`). The result is a `List` of rows, one per data point, each the vector `{f_1(coords), …, f_n(coords)}` — exactly the matrix that `Fit`/`LeastSquares` solve against.

**Data structures.** A `List` of `List`s of evaluated basis-function values; entry representation (exact vs. machine vs. MPFR) follows the same `WorkingPrecision` handling as `Fit`. This is the Vandermonde-like basis matrix for a polynomial basis but works for any list of basis functions.

**Attributes:** `Protected`.

## See also

[Fit](../../linear-algebra/Fit/), [List](../../other-advanced/List/)

## References

- Source: [`src/fit.c`](https://github.com/stblake/mathilda/blob/main/src/fit.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_fit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fit.c)

## Notes & additional examples

### Notes

`DesignMatrix[data, {f1, ..., fn}, vars]` builds the matrix whose rows are the
basis functions `f_i` evaluated at each data coordinate — exactly the matrix
`Fit` assembles internally before solving the normal equations. The data shapes
match `Fit`: with a single variable, each row is `{x_k, y_k}` (or just `{x_k}`)
and only the leading coordinate(s) are substituted, so the response column is
ignored when forming the design entries. A polynomial basis `{1, x, x^2, x^3}`
therefore produces a Vandermonde matrix. For several predictors the basis may mix
the variables freely (`{1, x, y, x*y}`). Entries are kept exact unless
`WorkingPrecision -> MachinePrecision` or a digit count is supplied, in which
case each entry is converted to an approximate number — useful when the basis
contains transcendental functions such as `Sin`.
