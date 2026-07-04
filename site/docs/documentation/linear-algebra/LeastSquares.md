# LeastSquares

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LeastSquares[m, b]
    finds an x that solves the linear least-squares problem
    for the matrix equation m . x == b, i.e. an x minimising
    Norm[m . x - b].
LeastSquares[m, b, Method -> "<name>"]
    selects an explicit solver.
LeastSquares[m, b, Tolerance -> t]
    specifies the singular-value truncation tolerance forwarded
    to the underlying PseudoInverse call (Tolerance -> Automatic
    by default).  Method and Tolerance options may appear in any
    order.

LeastSquares works on every input family supported by
PseudoInverse: exact (Integer / Rational), symbolic, inexact
(Real / MPFR), and complex.  The matrix m may be square or
rectangular and of any rank.  When m is rank-deficient the
result is the minimum-norm minimiser -- the Moore-Penrose
pseudoinverse solution PseudoInverse[m] . b.

The right-hand side b may be a vector or a matrix.  When b is
a matrix (one column per RHS), LeastSquares returns a matrix
of solutions, the j-th column of which is the least-squares
solution for the j-th column of b -- minimising
Norm[m . x - b, "Frobenius"] over the multi-RHS system.

Accepted Method names:
  "Automatic"           — alias for "Direct" (default)
  "Direct"              — PseudoInverse[m] . b; works for all
                          input families (dense or sparse,
                          exact or numeric, real or complex)
  "IterativeRefinement" — residual-correction loop on Direct,
                          x <- x + PseudoInverse[m] . (b - m.x),
                          terminated when ||dx||^2 <= Tolerance^2
                          or at a 50-iteration cap.  Exact inputs
                          converge in one pass; inexact inputs
                          drive round-off down to Tolerance.
  "Krylov"              — Conjugate-Gradient-on-Least-Squares
                          (Hestenes-Stiefel CG on the normal
                          equations) with x_0 = 0.  Converges
                          to the minimum-norm LS solution for
                          rank-deficient m, capped at 2 cols(m)
                          + 10 iterations.  Free symbols fall
                          back to Direct.
  "LSQR"                — Paige-Saunders LSQR: Lanczos
                          bidiagonalisation with Givens rotations.
                          Free symbols fall back to Direct; exact
                          rationals and complex inputs fall back
                          to Krylov / CGLS (equivalent without
                          square-root growth); pure-real inputs
                          with at least one Real entry run the
                          double-precision algorithm.

An unknown Method or Tolerance value leaves the call
unevaluated.  When m . x == b has an exact solution,
LeastSquares[m, b] coincides with LinearSolve[m, b].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}]
Out[1]= {19/3, 1/2}

In[2]:= LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {2, -4, 2}]
Out[2]= {0, 0, 0}

In[3]:= LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}, {1, 2, 4, 8}]
Out[3]= {157/180, 23/90, -13/36}

In[4]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {{7, 1}, {7, 2}, {8, 3}}]
Out[4]= {{19/3, 0}, {1/2, 1}}

In[5]:= LeastSquares[IdentityMatrix[4], {1, 2, 3, 4}]
Out[5]= {1, 2, 3, 4}

In[6]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Method -> "Krylov"]
Out[6]= {19/3, 1/2}
```

## Implementation notes

**Algorithm.** `builtin_leastsquares` returns an `x` minimising `Norm[m.x - b]`. The default `"Direct"` method computes the Moore–Penrose solution `x = PseudoInverse[m] . b` (reusing the `PseudoInverse` builtin, which itself does an exact full-rank `B·C` decomposition via `RowReduce`), so for full-column-rank `m` it gives the unique minimiser and for rank-deficient `m` the minimum-norm minimiser. The `LeastSquares == PseudoInverse . b` identity is the fundamental specification, and `"Direct"` works on every input family (integer, rational, symbolic, Real/MPFR, complex). Optional iterative methods are also provided and dispatched by input type:

- `"IterativeRefinement"` — residual-correction loop on top of Direct (one pass for exact input; drives round-off to `Tolerance` for inexact).
- `"LSQR"` — Paige–Saunders LSQR via Lanczos bidiagonalisation with Givens rotation updates, using their `|φ̄·α_{k+1}|` estimate of `‖Aᵀr‖`; symbolic input falls back to Direct, exact/complex input to CGLS, pure-real inexact to the canonical double-precision algorithm.
- `"Krylov"` — Conjugate-Gradient on the normal equations (CGLS / Hestenes–Stiefel), restricted to numeric inputs; symbolic falls back to Direct.

`Method` and `Tolerance` options may appear in either order.

**Data structures.** Solutions are built from the generic `PseudoInverse`/`Dot`/`Plus`/`Times` evaluator pipeline; the iterative kernels operate on the numeric tower (Real/MPFR/Complex) with exact-arithmetic variants for Integer/Rational inputs to avoid square-root growth.

**Complexity / limits.** Direct cost is dominated by `PseudoInverse` (exact Gauss-Jordan, `O(mn·rank)` style). Krylov/LSQR iterate to a `2·cols + O(1)` cap with `Tolerance`-based stopping; for symbolic inputs only Direct is well-defined since the iterative stopping tests are undecidable.

- `Protected`.
- The matrix `m` may be square or rectangular and of any rank. When `m`

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- C. C. Paige, M. A. Saunders, "LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares", ACM TOMS 8 (1982).
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- Source: [`src/linalg/lstsq.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/lstsq.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {1, 2, 2}]
Out[1]= {2/3, 1/2}

In[2]:= LeastSquares[{{1, 0}, {0, 1}, {1, 1}}, {1, 1, 3}]
Out[2]= {4/3, 4/3}
```

A least-squares quadratic fit `{c0, c1, c2}` of `y = {6, 5, 7, 10}` at `x = 1, 2, 3, 4` (Vandermonde design matrix), returned as exact rationals:

```mathematica
In[1]:= LeastSquares[{{1, 1, 1}, {1, 2, 4}, {1, 3, 9}, {1, 4, 16}}, {6, 5, 7, 10}]
Out[1]= {17/2, -18/5, 1}
```

### Notes

`LeastSquares[m, b]` returns the `x` minimising `Norm[m . x - b]` for the
overdetermined system `m . x == b`. With exact (rational) input it gives an
exact rational answer via the pseudoinverse; pass `Method ->` or
`Tolerance ->` to control the solver and singular-value truncation.
