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

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- The matrix `m` may be square or rectangular and of any rank. When `m`

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
