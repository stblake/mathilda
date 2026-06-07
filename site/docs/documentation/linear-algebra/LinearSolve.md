# LinearSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LinearSolve[m, b]
    finds an x that solves the matrix equation m . x == b.
LinearSolve[m, b, Method -> "<name>"]
    runs a specific elimination algorithm.

LinearSolve works on both numerical and symbolic matrices.
The matrix m may be square or rectangular.
The argument b may be a vector or a matrix; when b is a matrix
(one column per RHS) LinearSolve returns a matrix of solutions.
Higher-rank tensor inputs are also supported: when m has
dimensions {d1, ..., d(N-1), n}, b may have dimensions
{d1, ..., d(N-1), e1, ..., ep} and the result has dimensions
{n, e1, ..., ep}.

For under-determined systems LinearSolve returns a particular
solution in which the free (non-pivot) variables are taken to be
0; Solve returns the general solution.  When the equation has no
solution LinearSolve emits LinearSolve::nosol and returns
unevaluated.

Accepted method names:
  "Automatic"                 — alias for "DivisionFreeRowReduction" (default)
  "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan on [m | b]
  "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot
  "CofactorExpansion"         — Cramer's rule via Laplace cofactor expansion
                                 (square non-singular m only; LinearSolve::cofnsq
                                 / ::cofsng on shape / singularity errors)

Default implementation: fraction-free Gauss-Jordan elimination on
the augmented matrix [m | b] (the Bareiss-like algorithm shared
with RowReduce and Inverse), so exact integer / rational /
symbolic inputs flow through without any spurious denominator
blow-up.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LinearSolve[{{r, s}, {t, u}}, {y, z}]
Out[1]= {(u y)/(-s t + r u) - (s z)/(-s t + r u), -(t y)/(-s t + r u) + (r z)/(-s t + r u)}

In[2]:= LinearSolve[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]
Out[2]= {{-3, -4}, {4, 5}}

In[3]:= LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, {9, 10, 11, 12}]
Out[3]= {-1, 2}

In[4]:= LinearSolve[{{1, 2, 3}, {4, 5, 6}}, {6, 15}]
Out[4]= {0, 3, 0}

In[5]:= LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1, 2}]
Out[5]= LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1, 2}]

In[6]:= LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Method -> "CofactorExpansion"]
Out[6]= {-4, 9/2}

In[7]:= LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, {9, 10, 11, 12}, Method -> "OneStepRowReduction"]
Out[7]= {-1, 2}
```

## Implementation notes

- `Protected`.
- The matrix `m` may be square or rectangular.
- The argument `b` may be a vector (in which case the result is a

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — solving linear systems by Gaussian elimination.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — LU factorisation and triangular solves.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
Out[1]= {-4, 9/2}
```

```mathematica
In[1]:= LinearSolve[{{2, 0}, {0, 4}}, {6, 8}]
Out[1]= {3, 2}
```

```mathematica
In[1]:= LinearSolve[{{1, 1}, {0, 1}}, {3, 1}]
Out[1]= {2, 1}
```

```mathematica
In[1]:= LinearSolve[{{1, 2}, {2, 4}}, {1, 3}]
Out[1]= LinearSolve[{{1, 2}, {2, 4}}, {1, 3}]
```

### Notes

`LinearSolve[m, b]` returns an `x` satisfying `m . x == b`; over exact integer inputs the solution is exact rational, as in the first example. The default algorithm is fraction-free (Bareiss-like) Gauss-Jordan elimination on the augmented matrix `[m | b]`, so no spurious denominators appear during elimination. When the system is inconsistent (last example) `LinearSolve::nosol` is emitted on stderr and the call is returned unevaluated; for under-determined systems a particular solution with free variables set to 0 is returned. The right-hand side may also be a matrix (one column per system), giving a matrix of solutions.
