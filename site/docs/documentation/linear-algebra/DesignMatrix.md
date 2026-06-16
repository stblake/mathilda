# DesignMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DesignMatrix[data, {f1, ..., fn}, vars] gives the design matrix with entries f_i evaluated at the data coordinates.
Data shapes match Fit. The WorkingPrecision option converts entries to machine or n-digit reals; otherwise they are exact.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DesignMatrix[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= {{1, 0}, {1, 1}, {1, 3}, {1, 5}}

In[2]:= DesignMatrix[{{0,0,0},{1,0,1},{0,1,2}}, {1, x, y}, {x, y}]
Out[2]= {{1, 0, 0}, {1, 1, 0}, {1, 0, 1}}
```

## Implementation notes

**Algorithm.** `builtin_designmatrix` (in `fit.c`, the same module as `Fit`) builds the design matrix `m_{ij} = f_i(coords_j)` for `DesignMatrix[data, {f1,…,fn}, vars]`. It normalises the data shape (`{v1,…}` → `{{1,v1},…}`; `{{x,v},…}` univariate; `{{x1,…,xk,v},…}` multivariate, dropping the trailing response value to obtain the coordinate vector for each row), then evaluates each basis function `f_i` at each data point by substituting the `vars` symbols with that row's coordinates (pattern/`ReplaceAll`-style substitution followed by `evaluate`). The result is a `List` of rows, one per data point, each the vector `{f_1(coords), …, f_n(coords)}` — exactly the matrix that `Fit`/`LeastSquares` solve against.

**Data structures.** A `List` of `List`s of evaluated basis-function values; entry representation (exact vs. machine vs. MPFR) follows the same `WorkingPrecision` handling as `Fit`. This is the Vandermonde-like basis matrix for a polynomial basis but works for any list of basis functions.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/fit.c`](https://github.com/stblake/mathilda/blob/main/src/fit.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

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
