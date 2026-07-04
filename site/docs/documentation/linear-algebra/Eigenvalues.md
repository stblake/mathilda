# Eigenvalues

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Eigenvalues[m]
    gives a list of the eigenvalues of the square matrix m.
Eigenvalues[{m, a}]
    gives the generalized eigenvalues of m with respect to a.
Eigenvalues[m, k]
    gives the first k eigenvalues (largest by absolute value).
Eigenvalues[m, -k]
    gives the k eigenvalues smallest in absolute value.
Eigenvalues[m, UpTo[k]]
    gives k eigenvalues, or as many as are available.

Eigenvalues are computed from the roots of the characteristic
polynomial Det[m - lambda I] (or Det[m - lambda a] for the
generalised case). Approximate (Real / MPFR) matrices flow through
the Solve rationalise -> solve -> numericalize pipeline and yield
numerical eigenvalues sorted in order of decreasing absolute value.
Repeated eigenvalues appear with their algebraic multiplicity.

Options:
    Cubics    -> True       (use radicals to solve cubics)
    Quartics  -> True       (use radicals to solve quartics)
    Method    -> Automatic  (numeric-matrix method dispatch)

Method values for approximate-numeric matrices:
  Automatic    selects Direct unless k is small (-> Arnoldi)
               or the matrix is Hermitian-banded (-> Banded).
  "Direct"     Hessenberg + implicit shifted QR (general); for
               Hermitian inputs tridiagonalisation + Wilkinson-
               shift symmetric QR.  Returns all eigenvalues.
  "Arnoldi"    Krylov-subspace iteration for the k extreme
               eigenvalues; accepts Method -> {"Arnoldi",
               MaxIterations -> n, Tolerance -> t, BasisSize -> m}.
  "Banded"     Hermitian only; auto-detects band structure and
               reduces to tridiagonal before symmetric QR.
  "FEAST"      Hermitian only; eigenvalues in a user-specified
               Interval -> {a, b}; accepts Method ->
               {"FEAST", "Interval" -> {a, b},
                "ContourPoints" -> Ne, "SubspaceSize" -> m0,
                "MaxIterations" -> k, "Tolerance" -> t}.
Non-numeric matrices ignore Method and use the symbolic
characteristic-polynomial pipeline.

Implementation status: "Direct" runs the hand-rolled Householder
tridiagonalisation + Wilkinson-shift symmetric QR kernel at
machine precision for real symmetric matrices, the Hessenberg
+ implicit double-shift Francis QR kernel for real non-symmetric
matrices, a complex Householder tridiagonalisation + diagonal-
phase correction + symmetric QR kernel for complex Hermitian
matrices (returns real eigenvalues sorted by |lambda|
descending), and a real-block-embedding kernel for complex
non-Hermitian matrices (M = [[Re A, -Im A], [Im A, Re A]]
routed through real Hessenberg + Francis QR with grouped
complex Gram-Schmidt disambiguation of M's spec to recover
spec(A)).  Automatic routes here too.  Arbitrary-precision
(MPFR) inputs go through a parallel "Direct" kernel at the
input's combined precision: all four shapes -- real symmetric
(step 2d-A), real non-symmetric (step 2d-B), complex Hermitian
(step 2d-C), and complex non-Hermitian (step 2d-D) -- return
eigenvalues / eigenvectors carrying full input precision.
"Arnoldi" is implemented in Phase 3 at both machine and MPFR
precision: m-step classical Gram-Schmidt with one re-orthog-
onalisation pass builds the Krylov basis V_m and the m x m
upper Hessenberg H_m; H_m is diagonalised by reusing the
"Direct" Francis QR pipeline (real machine, real MPFR, or via
a 2mu x 2mu real-block embedding for complex H_m), and Ritz
vectors V_m y_i lift back to A-eigenvectors.  Complex inputs
use paired re/im storage for V_m and H_m.  Automatic routes to
Arnoldi when n > 32 and k_spec is given with k <= max(20, n/10);
small matrices always go through Direct (faster + exact).
Default BasisSize is max(2k, 20) capped at n; on lucky breakdown
(||w|| below tolerance at some step j) Arnoldi terminates early
with j+1 exact eigenpairs.  MPFR Arnoldi carries through the
input's combined precision via the same scratch-pool discipline
as the Direct MPFR kernels.
"Banded" (Phase 4, machine + MPFR) handles real symmetric and
complex Hermitian matrices.  It auto-detects the half-bandwidth
and reduces to symmetric tridiagonal form via Schwarz-style
two-sided Givens rotations with bulge chasing (one off-band
entry zeroed per Givens; the introduced bulge is chased b
columns at a time until it falls past the matrix edge); the
resulting tridiagonal eigenproblem reuses the Phase 2 symmetric
QR.  Complex Hermitian banded uses paired re/im Givens with a
real-c / complex-s parameterisation and the same phase-
correction step as the Direct Hermitian kernel.  Banded refuses
(returns NULL, falls back to Direct) when the matrix isn't
Hermitian or when it's fully dense (b == n - 1).  Automatic
routes here when the matrix is Hermitian, n > 8, and the half-
bandwidth is at most max(8, n/4); narrower bands save more flops
than wider ones.
"FEAST" (Phase 5, machine + MPFR) handles Hermitian (real
symmetric or complex Hermitian) input and returns only the
eigenvalues in a user-supplied real Interval -> {a, b} -- a
spectral-slice query rather than a full decomposition.  Uses
Ne-point Gauss-Legendre quadrature (default Ne = 8; supported:
2, 4, 8, 16) on the upper half of the elliptic contour through
(a, 0) and (b, 0) to approximate the spectral projector
P_[a,b](A); Schwarz symmetry halves the number of complex
linear solves.  A Rayleigh-Ritz reduction with Cholesky
B_q = L L^* extracts the in-interval eigenpairs by reusing the
Direct symmetric / Hermitian kernel on L^-1 A_q L^-*.  Output
is sorted by |lambda| descending so k_spec composes naturally
with the in-interval filter.  Automatic never routes to FEAST
(it requires the user to commit to an interval).  FEAST falls
back to Direct (NULL return + one-shot stderr warning tagged
with the reason) on: non-Hermitian input, missing Interval,
degenerate or invalid {a, b} (interval_high <= interval_low),
generalised eigenproblem, Cholesky failure on B_q (subspace
too small for the spectral count), LU singular at any quad-
rature node, or non-convergence within MaxIterations.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {1/2 (15 + 3 Sqrt[33]), 1/2 (15 - 3 Sqrt[33]), 0}

In[2]:= Eigenvalues[{{a, b}, {c, d}}]
Out[2]= {1/2 (a + d + Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1/2 (a + d - Sqrt[(-a - d)^2 - 4 (-b c + a d)])}

In[3]:= Eigenvalues[IdentityMatrix[12]]
Out[3]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
```

## Implementation notes

**Algorithm.** `builtin_eigenvalues` (dispatcher in `eigen.c`; kernels in `eigen_common.c`, `eigen_direct.c`, `eigen_banded.c`, `eigen_arnoldi.c`, `eigen_feast.c`) chooses between an exact/symbolic path and numerical kernels based on whether the matrix is inexact.

*Exact / symbolic path* (`eigen_compute_eigenvalues_full`): the characteristic polynomial is formed by the **Faddeev–Leverrier–Souriau** recurrence (`eigen_char_poly_faddeev`), which builds the coefficients in `O(n^4)` matrix multiplications — far cheaper than Laplace expansion of `det(m − λI)` over a polynomial-entry matrix (`O(n!)`) once `n` grows. Its roots are found by `eigen_solve_poly` (radical formulas for cubics/quartics controlled by `Cubics`/`Quartics`, otherwise `Root`/`Solve`), extracted with multiplicity, and sorted by descending `Abs`. The **generalised** problem `Eigenvalues[{m, a}]` instead forms `det(m − λa)` by Laplace expansion (only used for small pencils) and pads the short result with `Infinity` for the degree-drop branch.

*Numerical path* (inexact input): dispatched by `Method`. `Direct` (`direct_dispatch`) uses **Householder reduction** to tridiagonal (symmetric, Golub & Van Loan Alg. 8.3.1) followed by **implicit-shift symmetric tridiagonal QR with Wilkinson shift** (`direct_symtridiag_qr`) for symmetric matrices, and Householder reduction to upper **Hessenberg + Francis implicit double-shift QR** for non-symmetric matrices. `Automatic` prefers `Banded` for narrow-band Hermitian input and `Arnoldi` when only a small `k` is requested, with FEAST available for Hermitian interval problems; each dispatcher returns `NULL` for shapes it doesn't support and falls through to `Direct`, then ultimately to the symbolic pipeline. Numerical-noise imaginary parts are chopped (`eigen_chop`). A `{k}`/`-k`/`UpTo[k]` spec trims the sorted result (`eigen_apply_k_spec`).

**Data structures.** Symbolic side: `Expr` trees for the polynomial and roots. Numerical kernels operate on dense row-major `double` (and MPFR) matrices; the source carries `LAPACK-HOOK` annotations marking where `dsytrd`/`dsteqr`/`dgehrd`/`dhseqr` would drop in under a `USE_LAPACK` build (the hooks are present but no LAPACK backend is currently wired).

**Complexity / limits.** Faddeev–Leverrier is `O(n^4)`; symmetric QR is `O(n^3)` with cubic per-eigenvalue convergence under the Wilkinson shift. The generalised path is restricted to small `n` (Laplace expansion). Closed-form eigenvalues are limited by the degree-≤4 radical solver; higher degrees come back as `Root` objects.

- `Protected`.
- Implemented via the characteristic polynomial `Det[m - lambda I]`

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the symmetric and unsymmetric eigenvalue problems.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — eigenvalue algorithms and the QR iteration.
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- J. H. Wilkinson, *The Algebraic Eigenvalue Problem* (Oxford University Press, 1965).
- Source: [`src/linalg/eigen.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/eigen.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Eigenvalues[{{2, 1}, {0, 3}}]
Out[1]= {3, 2}
```

```mathematica
In[1]:= Eigenvalues[{{2, 0}, {0, 5}}]
Out[1]= {5, 2}
```

```mathematica
In[1]:= Eigenvalues[{{0, 1}, {-1, 0}}]
Out[1]= {-I, I}
```

```mathematica
In[1]:= Eigenvalues[{{a, b}, {c, d}}]
Out[1]= {1/2 (a + d + Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1/2 (a + d - Sqrt[(-a - d)^2 - 4 (-b c + a d)])}
```

```mathematica
In[1]:= Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}]
Out[1]= {1, -(-1)^(1/3), (-1)^(2/3)}
```

```mathematica
In[1]:= Eigenvalues[{{5, 4, 2}, {4, 5, 2}, {2, 2, 2}}]
Out[1]= {10, 1, 1}
```

### Notes

Eigenvalues of an exact matrix are computed as the roots of the characteristic polynomial `Det[m - lambda I]`, so triangular matrices return their diagonal entries directly and a rotation matrix returns the complex conjugate pair `{-I, I}`. Results are ordered by decreasing absolute value, and repeated eigenvalues appear with their full algebraic multiplicity. Symbolic matrices return closed-form roots (the `2x2` case uses the quadratic formula). Approximate Real / MPFR matrices flow through dedicated numerical kernels (Householder + QR for the symmetric path, Hessenberg + Francis QR for the general path) selectable through the `Method` option.
