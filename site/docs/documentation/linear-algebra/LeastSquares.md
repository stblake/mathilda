# LeastSquares

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LeastSquares[m, b]`**

finds an x that solves the linear least-squares problem for the matrix equation m . x == b, i.e. an x minimising Norm\[m . x - b\].

**`LeastSquares[m, b, Method -> "<name>"]`**

selects an explicit solver.

**`LeastSquares[m, b, Tolerance -> t]`**

specifies the singular-value truncation tolerance forwarded to the underlying PseudoInverse call (Tolerance -\> Automatic by default).  Method and Tolerance options may appear in any order.

**`Norm[m . x - b, "Frobenius"] over the multi-RHS system.`**

**`LeastSquares[m, b] coincides with LinearSolve[m, b].`**

<details>
<summary>Notes</summary>

LeastSquares works on every input family supported by PseudoInverse: exact (Integer / Rational), symbolic, inexact (Real / MPFR), and complex.  The matrix m may be square or rectangular and of any rank.  When m is rank-deficient the result is the minimum-norm minimiser -- the Moore-Penrose pseudoinverse solution PseudoInverse\[m\] . b. The right-hand side b may be a vector or a matrix.  When b is a matrix (one column per RHS), LeastSquares returns a matrix of solutions, the j-th column of which is the least-squares solution for the j-th column of b -- minimising Accepted Method names: "Automatic"           — alias for "Direct" (default) "Direct"              — PseudoInverse\[m\] . b; works for all input families (dense or sparse, exact or numeric, real or complex) "IterativeRefinement" — residual-correction loop on Direct, x \<- x + PseudoInverse\[m\] . (b - m.x), terminated when ||dx||^2 \<= Tolerance^2 or at a 50-iteration cap.  Exact inputs converge in one pass; inexact inputs drive round-off down to Tolerance. "Krylov"              — Conjugate-Gradient-on-Least-Squares (Hestenes-Stiefel CG on the normal equations) with x\_0 = 0.  Converges to the minimum-norm LS solution for rank-deficient m, capped at 2 cols(m) + 10 iterations.  Free symbols fall back to Direct. "LSQR"                — Paige-Saunders LSQR: Lanczos bidiagonalisation with Givens rotations. Free symbols fall back to Direct; exact rationals and complex inputs fall back to Krylov / CGLS (equivalent without square-root growth); pure-real inputs with at least one Real entry run the double-precision algorithm. An unknown Method or Tolerance value leaves the call unevaluated.  When m . x == b has an exact solution,

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

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

In[6]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 7}] == LinearSolve[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 7}]
Out[6]= True
```

### Options (4)

```mathematica
In[7]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Method -> "IterativeRefinement", Tolerance -> 1/100]
Out[7]= {19/3, 1/2}

In[8]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Method -> "Krylov"]
Out[8]= {19/3, 1/2}

In[9]:= LeastSquares[{{3.2, 2.2, 1.2}, {2.1, 7.1, 8.5}, {9.5, 6.7, 3.7}}, {7., 8., 9.}, Method -> "LSQR"]
Out[9]= {73.9499, -174.379, 128.329}

In[10]:= LeastSquares[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}, {2., -4., 2.}, Method -> "LSQR"]
Out[10]= {0.0, 0.0, 0.0}
```

### Applications (3)

```mathematica
In[11]:= LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {1, 2, 2}]
Out[11]= {2/3, 1/2}

In[12]:= LeastSquares[{{1, 0}, {0, 1}, {1, 1}}, {1, 1, 3}]
Out[12]= {4/3, 4/3}

In[13]:= LeastSquares[{{1, 1, 1}, {1, 2, 4}, {1, 3, 9}, {1, 4, 16}}, {6, 5, 7, 10}]
Out[13]= {17/2, -18/5, 1}
```

## Algorithm

matlstsq.c

LeastSquares[m, b] -- linear least-squares solver.

```text
  LeastSquares[m, b]
      Returns an x that minimises Norm[m . x - b].  When m has full
      column rank the minimiser is unique; when m is rank-deficient
      the result is the minimum-norm minimiser
      (PseudoInverse[m] . b).

  LeastSquares[m, b, Method -> "<name>"]
  LeastSquares[m, b, Tolerance -> t]
      Optional Method and Tolerance arguments may appear in either
      order and both may be given.
```

Method names (parsed case-sensitively; "Automatic" symbol also accepted, matching the Mathematica grammar):

```text
  "Automatic"           -- alias for "Direct" (default).
  "Direct"              -- Moore-Penrose solve PseudoInverse[m] . b.
                           Works on every input family: exact
                           (Integer / Rational), symbolic, inexact
                           (Real / MPFR), and complex.  This is the
                           workhorse method; the LeastSquares ==
                           PseudoInverse . b identity is the
                           fundamental specification.

  "IterativeRefinement" -- residual-correction loop on top of Direct.
                           Starting from x = PseudoInverse[m] . b we
                           repeatedly compute r = b - m . x,
                           dx = PseudoInverse[m] . r, x <- x + dx
                           until ||dx||^2 <= tol^2 or a 50-iteration
                           cap is hit.  For exact inputs the first
                           correction is exactly zero (the pseudoinverse
                           identity x = pinv b implies dx = pinv (I - m
                           pinv) b = 0), so the loop converges in a
                           single pass; for inexact inputs the loop
                           drives round-off down to Tolerance.

  "LSQR"                -- Paige-Saunders LSQR.  Lanczos
                           bidiagonalisation of A combined with a
                           Givens-rotation update of the upper
                           triangular factor R, exactly the algorithm
                           from Paige & Saunders, ACM TOMS 1982.
                           Convergence test uses their |phi_bar *
                           alpha_{k+1}| estimate of ||A^T r||, scaled
                           against the initial gradient ||A^T b||.
                           Dispatches by input type: free symbols go
                           to Direct (the iteration's stopping test is
                           undecidable); exact (Integer / Rational)
                           and Complex inputs go to CGLS (Krylov is
                           mathematically equivalent and avoids the
                           square-root growth in exact arithmetic);
                           pure-real inexact inputs run the canonical
                           double-precision algorithm.

  "Krylov"              -- Conjugate-Gradient-on-Least-Squares
                           (Hestenes-Stiefel CG applied to the normal
                           equations).  Iterates
                               q = A p,  alpha = |s|^2 / |q|^2,
                               x <- x + alpha p,  r <- r - alpha q,
                               s = A^T r,  beta = |s_new|^2 / |s|^2,
                               p <- s_new + beta p
                           with x_0 = 0, r_0 = b, s_0 = A^T b, p_0 = s_0.
                           Stops on |s|^2 <= tol^2 (or exact zero), a
                           null search direction |q|^2 = 0, or a
                           2 cols(A) + 10 iteration cap.  Restricted to
                           numeric inputs (Integer / Rational / Real /
                           Complex with numeric components); symbolic
                           inputs fall back to Direct to avoid running
                           the loop with an undecidable convergence
                           test.  Matrix RHS are solved column by
                           column and recombined via Transpose.
```

Tolerance:

```text
  Tolerance -> Automatic (default), or a non-negative number.
  Forwarded verbatim as the Tolerance option of the underlying
  PseudoInverse call so a future singular-value-truncation pass in
  PseudoInverse picks it up automatically.  The iterative methods
  will also use Tolerance as a convergence threshold once they have
  dedicated implementations.
```

Memory ownership follows the standard builtin contract: this file owns `res` on success and frees it; on failure (returning NULL) the caller (the evaluator) retains ownership and the expression remains

```text
unevaluated.  Every intermediate PseudoInverse / Dot / Plus / Times
```

call goes through eval_and_free so its argument tree is consumed and its return value is owned by this function.

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
  has full column rank the minimiser is unique; when `m` is
  rank-deficient the result is the minimum-norm minimiser
  (`PseudoInverse[m] . b`).
- The right-hand side `b` may be a vector or a matrix. When `b` is a
  matrix, the result has one column per RHS — column `j` is the
  least-squares solution for `b[[All, j]]`, i.e. the `x` minimising
  `Norm[m . x - b, "Frobenius"]` over the multi-RHS system.
- Works on every input family supported by `PseudoInverse`:
  integer / rational, symbolic, machine-precision real / MPFR,
  exact and inexact complex.
- Method and Tolerance options may appear in any order; each may
  appear at most once. Duplicates or unknown option names leave the
  call unevaluated.
- Accepted Method names:
  - `Method -> Automatic` or `Method -> "Automatic"` — alias for `"Direct"` (default).
  - `Method -> "Direct"` — Moore-Penrose solve `PseudoInverse[m] . b`. Works on every input family (dense or sparse, exact or numeric, real or complex). The workhorse method. A **machine-real** coefficient matrix takes the same thin SVD `PseudoInverse` does, applying `V`, `S^+` and `U^T` to `b` in turn rather than forming the pseudo-inverse and multiplying: `LeastSquares[A500, b500]` did not finish in 180 s on the exact pipeline and costs 77 ms here. The iterative methods below are deliberately left on the exact path — they are chosen for their iteration, not their answer.
  - `Method -> "IterativeRefinement"` — residual-correction loop on top of Direct: `x <- PseudoInverse[m] . b`, then repeatedly `r = b - m . x`, `dx = PseudoInverse[m] . r`, `x <- x + dx`, capped at 50 iterations and terminated when `Total[Flatten[dx]^2] <= Tolerance^2`. For exact inputs the first correction is exactly zero (Moore-Penrose identity) so the loop converges in one pass; for inexact inputs the loop drives round-off down to Tolerance.
  - `Method -> "Krylov"` — Conjugate-Gradient-on-Least-Squares (Hestenes-Stiefel CG applied to the normal equations). Iterates `q = A p`, `alpha = |s|^2 / |q|^2`, `x <- x + alpha p`, `r <- r - alpha q`, `s = A^T r`, `beta = |s_new|^2 / |s|^2`, `p <- s_new + beta p` from `x_0 = 0` (so the iterate stays in `range(A^T)` and converges to the minimum-norm LS solution for rank-deficient `m`). Capped at `2 cols(m) + 10` iterations. Matrix RHS are solved column-by-column and recombined via `Transpose`. Symbolic inputs (the convergence test is undecidable for them) fall back to Direct.
  - `Method -> "LSQR"` — Paige-Saunders LSQR (ACM TOMS 1982): Lanczos bidiagonalisation of `m` with a Givens-rotation update of the resulting upper triangular factor. Dispatches by input grammar: free-symbol inputs go to Direct; exact (Integer / Rational) and Complex inputs go to Krylov / CGLS (mathematically equivalent without the square-root growth in exact arithmetic); pure-real inputs with at least one Real entry run the canonical double-precision algorithm. Uses the Paige-Saunders estimate `|phi_bar * alpha_{k+1}|` of `||A^T r||` for the convergence test, scaled against the initial gradient `||A^T b||`. Cap is `2 cols(m) + 10` iterations. Detects rank deficiency by monitoring `alpha_new / max(|A_ij|)` and `beta_new / max(|A_ij|)`, avoiding the catastrophic blowup of dividing by a near-zero `alpha`.
- `Tolerance -> Automatic` (default), or a non-negative integer / real /
  `Rational`. Forwarded verbatim as the Tolerance option of the
  underlying `PseudoInverse` call so any future singular-value
  truncation pass in `PseudoInverse` is picked up automatically.
- When `m . x == b` is consistent, `LeastSquares[m, b]` coincides with
  `LinearSolve[m, b]`.
- Satisfies the identity `LeastSquares[m, b] == PseudoInverse[m] . b`.
- Issues `LeastSquares::matrix` / `::lvec` / `::lvec1` and returns
  unevaluated for shape errors.
- Lives in `src/linalg/lstsq.c`.

**Attributes:** `Protected`.

## References

**See also:** [PseudoInverse](../../linear-algebra/PseudoInverse/), [Transpose](../../structural-manipulation/Transpose/), [Rational](../../arithmetic/Rational/)

- C. C. Paige, M. A. Saunders, "LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares", ACM TOMS 8 (1982).
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- Source: [`src/linalg/lstsq.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/lstsq.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_matlstsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matlstsq.c)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`LeastSquares[m, b]` returns the `x` minimising `Norm[m . x - b]` for the
overdetermined system `m . x == b`. With exact (rational) input it gives an
exact rational answer via the pseudoinverse; pass `Method ->` or
`Tolerance ->` to control the solver and singular-value truncation.
