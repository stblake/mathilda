# Eigenvalues

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Eigenvalues[m]`**

gives a list of the eigenvalues of the square matrix m.

**`Eigenvalues[{m, a}]`**

gives the generalized eigenvalues of m with respect to a.

**`Eigenvalues[m, k]`**

gives the first k eigenvalues (largest by absolute value).

**`Eigenvalues[m, -k]`**

gives the k eigenvalues smallest in absolute value.

**`Eigenvalues[m, UpTo[k]]`**

gives k eigenvalues, or as many as are available.

<details>
<summary>Notes</summary>

Eigenvalues are computed from the roots of the characteristic polynomial Det\[m - lambda I\] (or Det\[m - lambda a\] for the generalised case). Approximate (Real / MPFR) matrices flow through the Solve rationalise -\> solve -\> numericalize pipeline and yield numerical eigenvalues sorted in order of decreasing absolute value. Repeated eigenvalues appear with their algebraic multiplicity. Options: Cubics    -\> False      (radicals for a general cubic; True forces them) Quartics  -\> False      (radicals for a general quartic; True forces them) Method    -\> Automatic  (numeric-matrix method dispatch) With Cubics/Quartics -\> False (the default) a general irreducible cubic or quartic characteristic polynomial is returned as held Root\[\] objects. Special always-solvable families -- binomials, quadratic-in-x^m, and biquadratic-after-depression quartics -- are always returned in radical form regardless of these options. Method values for approximate-numeric matrices: Automatic    selects Direct unless k is small (-\> Arnoldi) or the matrix is Hermitian-banded (-\> Banded). "Direct"     Hessenberg + implicit shifted QR (general); for Hermitian inputs tridiagonalisation + Wilkinson- shift symmetric QR.  Returns all eigenvalues. "Arnoldi"    Krylov-subspace iteration for the k extreme eigenvalues; accepts Method -\> {"Arnoldi", MaxIterations -\> n, Tolerance -\> t, BasisSize -\> m}. "Banded"     Hermitian only; auto-detects band structure and reduces to tridiagonal before symmetric QR. "FEAST"      Hermitian only; eigenvalues in a user-specified Interval -\> {a, b}; accepts Method -\> {"FEAST", "Interval" -\> {a, b}, "ContourPoints" -\> Ne, "SubspaceSize" -\> m0, "MaxIterations" -\> k, "Tolerance" -\> t}. Non-numeric matrices ignore Method and use the symbolic characteristic-polynomial pipeline. Implementation status: "Direct" runs the hand-rolled Householder tridiagonalisation + Wilkinson-shift symmetric QR kernel at machine precision for real symmetric matrices, the Hessenberg + implicit double-shift Francis QR kernel for real non-symmetric matrices, a complex Householder tridiagonalisation + diagonal- phase correction + symmetric QR kernel for complex Hermitian matrices (returns real eigenvalues sorted by |lambda| descending), and a real-block-embedding kernel for complex non-Hermitian matrices (M = \[\[Re A, -Im A\], \[Im A, Re A\]\] routed through real Hessenberg + Francis QR with grouped complex Gram-Schmidt disambiguation of M's spec to recover spec(A)).  Automatic routes here too.  Arbitrary-precision (MPFR) inputs go through a parallel "Direct" kernel at the input's combined precision: all four shapes -- real symmetric (step 2d-A), real non-symmetric (step 2d-B), complex Hermitian (step 2d-C), and complex non-Hermitian (step 2d-D) -- return eigenvalues / eigenvectors carrying full input precision. "Arnoldi" is implemented in Phase 3 at both machine and MPFR precision: m-step classical Gram-Schmidt with one re-orthog- onalisation pass builds the Krylov basis V\_m and the m x m upper Hessenberg H\_m; H\_m is diagonalised by reusing the "Direct" Francis QR pipeline (real machine, real MPFR, or via a 2mu x 2mu real-block embedding for complex H\_m), and Ritz vectors V\_m y\_i lift back to A-eigenvectors.  Complex inputs use paired re/im storage for V\_m and H\_m.  Automatic routes to Arnoldi when n \> 32 and k\_spec is given with k \<= max(20, n/10); small matrices always go through Direct (faster + exact). Default BasisSize is max(2k, 20) capped at n; on lucky breakdown (||w|| below tolerance at some step j) Arnoldi terminates early with j+1 exact eigenpairs.  MPFR Arnoldi carries through the input's combined precision via the same scratch-pool discipline as the Direct MPFR kernels. "Banded" (Phase 4, machine + MPFR) handles real symmetric and complex Hermitian matrices.  It auto-detects the half-bandwidth and reduces to symmetric tridiagonal form via Schwarz-style two-sided Givens rotations with bulge chasing (one off-band entry zeroed per Givens; the introduced bulge is chased b columns at a time until it falls past the matrix edge); the resulting tridiagonal eigenproblem reuses the Phase 2 symmetric QR.  Complex Hermitian banded uses paired re/im Givens with a real-c / complex-s parameterisation and the same phase- correction step as the Direct Hermitian kernel.  Banded refuses (returns NULL, falls back to Direct) when the matrix isn't Hermitian or when it's fully dense (b == n - 1).  Automatic routes here when the matrix is Hermitian, n \> 8, and the half- bandwidth is at most max(8, n/4); narrower bands save more flops than wider ones. "FEAST" (Phase 5, machine + MPFR) handles Hermitian (real symmetric or complex Hermitian) input and returns only the eigenvalues in a user-supplied real Interval -\> {a, b} -- a spectral-slice query rather than a full decomposition.  Uses Ne-point Gauss-Legendre quadrature (default Ne = 8; supported: 2, 4, 8, 16) on the upper half of the elliptic contour through (a, 0) and (b, 0) to approximate the spectral projector P\_\[a,b\](A); Schwarz symmetry halves the number of complex linear solves.  A Rayleigh-Ritz reduction with Cholesky B\_q = L L^\* extracts the in-interval eigenpairs by reusing the Direct symmetric / Hermitian kernel on L^-1 A\_q L^-\*.  Output is sorted by |lambda| descending so k\_spec composes naturally with the in-interval filter.  Automatic never routes to FEAST (it requires the user to commit to an interval).  FEAST falls back to Direct (NULL return + one-shot stderr warning tagged with the reason) on: non-Hermitian input, missing Interval, degenerate or invalid {a, b} (interval\_high \<= interval\_low), generalised eigenproblem, Cholesky failure on B\_q (subspace too small for the spectral count), LU singular at any quad- rature node, or non-convergence within MaxIterations.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Scope (4)

```mathematica
In[1]:= Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {1/2 (15 + 3 Sqrt[33]), 1/2 (15 - 3 Sqrt[33]), 0}

In[2]:= Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}]
Out[2]= {4, 4, 3, 2}

In[3]:= Eigenvalues[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}}, {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}]
Out[3]= {Infinity, 1/2 (1 + Sqrt[5]), 1/2 (1 - Sqrt[5])}

In[4]:= Eigenvalues[IdentityMatrix[12]]
Out[4]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
```

### Options (2)

```mathematica
In[5]:= Eigenvalues[N[{{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0}, {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1}, {0, 0, 0, -1, 2}}], Method -> {"FEAST", "Interval" -> {2.5, 4}}]
Out[5]= {3.73205, 3.0}

In[6]:= Eigenvalues[N[{{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0}, {0, -1, 2, -1, 0}, {0, 0, -1, 2, -1}, {0, 0, 0, -1, 2}}], 1, Method -> {"FEAST", "Interval" -> {0, 4}}]
Out[6]= {3.73205}
```

### Worked examples (1)

```mathematica
In[7]:= Eigenvalues[{{0,-1},{1,0}}]
Out[7]= {I, -I}
```

### Applications (6)

```mathematica
In[8]:= Eigenvalues[{{2, 1}, {0, 3}}]
Out[8]= {3, 2}

In[9]:= Eigenvalues[{{2, 0}, {0, 5}}]
Out[9]= {5, 2}

In[10]:= Eigenvalues[{{0, 1}, {-1, 0}}]
Out[10]= {-I, I}

In[11]:= Eigenvalues[{{a, b}, {c, d}}]
Out[11]= {1/2 (a + d + Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1/2 (a + d - Sqrt[(-a - d)^2 - 4 (-b c + a d)])}

In[12]:= Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}]
Out[12]= {1, -(-1)^(1/3), (-1)^(2/3)}

In[13]:= Eigenvalues[{{5, 4, 2}, {4, 5, 2}, {2, 2, 2}}]
Out[13]= {10, 1, 1}
```

## Options & behaviour

### `Method -> "FEAST"` (interval slice)

For Hermitian (real symmetric or complex Hermitian) inexact numeric
input, `Method -> "FEAST"` returns only the eigenpairs whose
eigenvalues lie in a user-supplied real interval `[a, b]` — a
spectral-slice query rather than a full decomposition.  Implements
the contour-integral spectral-projector algorithm of Polizzi (2009):
`Ne`-point Gauss-Legendre quadrature on the upper half of the
elliptic contour through `(a, 0)` and `(b, 0)` approximates
`P_[a, b](A) Y` (Schwarz symmetry halves the number of complex
linear solves), then a Rayleigh–Ritz reduction with Cholesky
`B_q = L L^*` extracts the in-interval eigenpairs.

Sub-options:

| Key | Default | Meaning |
|-----|---------|---------|
| `"Interval" -> {a, b}` | required | Real interval to slice.  Auto-swapped if `a > b`. |
| `"ContourPoints" -> Ne` | `8` | Gauss-Legendre order.  Supported: `2`, `4`, `8`, `16`. |
| `"SubspaceSize" -> m0` | `Max[20, n/4]` (capped at `n`) | Working subspace dimension; must be `>=` spectral count in `[a, b]`. |
| `"MaxIterations" -> k` | `20` | Outer iteration cap. |
| `"Tolerance" -> t` | precision-aware | Residual stopping criterion. |

Automatic never routes to FEAST; explicit `Method -> "FEAST"` is the
only way in.  The output is sorted by `|lambda|` descending so the
optional `k_spec` (`Eigenvalues[m, k]`, `Eigenvalues[m, -k]`,
`Eigenvalues[m, UpTo[k]]`) composes naturally with the in-interval
filter.  MPFR inputs run a parallel kernel at the input's combined
precision.

### Fail-soft cascade

FEAST returns `NULL` (and the call falls
through to Direct) on any of: non-Hermitian input, missing
`"Interval"`, `interval_high <= interval_low` (catches degenerate
`{c, c}` and NaN coercion failures), generalised eigenproblem
(`Eigenvalues[{m, a}]`), Cholesky failure on `B_q` (subspace too
small for the spectral count in the interval), LU singular at any
quadrature node, or non-convergence within `MaxIterations`.  The
first such fall-back per process emits a single
`Eigenvalues::feast: ... falling back to the Direct method.` stderr
line tagged with the reason, so an explicit FEAST call always
returns *some* sensible answer — at worst the full Direct spectrum.

## Implementation notes

**Algorithm.** `builtin_eigenvalues` (dispatcher in `eigen.c`; kernels in `eigen_common.c`, `eigen_direct.c`, `eigen_banded.c`, `eigen_arnoldi.c`, `eigen_feast.c`) chooses between an exact/symbolic path and numerical kernels based on whether the matrix is inexact.

*Exact / symbolic path* (`eigen_compute_eigenvalues_full`): the characteristic polynomial is formed by the **Faddeev–Leverrier–Souriau** recurrence (`eigen_char_poly_faddeev`), which builds the coefficients in `O(n^4)` matrix multiplications — far cheaper than Laplace expansion of `det(m − λI)` over a polynomial-entry matrix (`O(n!)`) once `n` grows. Its roots are found by `eigen_solve_poly` (radical formulas for cubics/quartics controlled by `Cubics`/`Quartics`, otherwise `Root`/`Solve`), extracted with multiplicity, and sorted by descending `Abs`. The **generalised** problem `Eigenvalues[{m, a}]` instead forms `det(m − λa)` by Laplace expansion (only used for small pencils) and pads the short result with `Infinity` for the degree-drop branch.

*Numerical path* (inexact input): dispatched by `Method`. `Direct` (`direct_dispatch`) uses **Householder reduction** to tridiagonal (symmetric, Golub & Van Loan Alg. 8.3.1) followed by **implicit-shift symmetric tridiagonal QR with Wilkinson shift** (`direct_symtridiag_qr`) for symmetric matrices, and Householder reduction to upper **Hessenberg + Francis implicit double-shift QR** for non-symmetric matrices. `Automatic` prefers `Banded` for narrow-band Hermitian input and `Arnoldi` when only a small `k` is requested, with FEAST available for Hermitian interval problems; each dispatcher returns `NULL` for shapes it doesn't support and falls through to `Direct`, then ultimately to the symbolic pipeline. Numerical-noise imaginary parts are chopped (`eigen_chop`). A `{k}`/`-k`/`UpTo[k]` spec trims the sorted result (`eigen_apply_k_spec`).

**Data structures.** Symbolic side: `Expr` trees for the polynomial and roots. Numerical kernels operate on dense row-major `double` (and MPFR) matrices; the source carries `LAPACK-HOOK` annotations marking where `dsytrd`/`dsteqr`/`dgehrd`/`dhseqr` would drop in under a `USE_LAPACK` build (the hooks are present but no LAPACK backend is currently wired).

**Complexity / limits.** Faddeev–Leverrier is `O(n^4)`; symmetric QR is `O(n^3)` with cubic per-eigenvalue convergence under the Wilkinson shift. The generalised path is restricted to small `n` (Laplace expansion). Closed-form eigenvalues are limited by the degree-≤4 radical solver; higher degrees come back as `Root` objects.

- `Protected`.
- Implemented via the characteristic polynomial `Det[m - lambda I]`
  (or `Det[m - lambda a]`) followed by `Solve`. The ordinary case uses a
  Faddeev–Leverrier–Souriau fast path in `O(n^4)` matrix multiplications,
  so eigenvalues of large rational / diagonal matrices return instantly.
- Approximate (`Real` / MPFR) input flows through Solve's rationalise →
  solve → numericalize pipeline; tiny imaginary noise introduced by the
  Cardano formula on real cubics is chopped automatically.
- Repeated eigenvalues appear with their algebraic multiplicity.
- Eigenvalues are sorted in order of decreasing absolute value. On an exact
  modulus tie — the case that matters is a complex-conjugate pair `a ± b I` —
  the `+imag` member is listed first, matching Mathematica
  (`Eigenvalues[{{0,-1},{1,0}}]` → `{I, -I}`). This holds for both the numeric
  Direct kernel and the symbolic (`Root[]` / radical) path, so
  `N[Eigenvalues[m]]` and `Eigenvalues[N[m]]` agree on the sign order
  (e.g. `Eigenvalues[{{0,1,0},{0,0,1},{1,1,0}}]` returns the characteristic
  roots as `{Root[…,1], Root[…,3], Root[…,2]}`, whose `N` is
  `{1.32472, -0.662+0.562 I, -0.662-0.562 I}`). Symbolic eigenvalues that carry
  free variables (no concrete modulus) retain Solve's natural order.
- When `m, a` have a shared null space, `Eigenvalues[{m, a}]` returns
  `Infinity` for each degree drop in the characteristic polynomial.
- Options: `Cubics -> False`, `Quartics -> False` (defaults), `Method`.  With
  the defaults a general irreducible cubic/quartic characteristic polynomial is
  returned as held `Root[]` objects (matching Mathematica); `Cubics -> True` /
  `Quartics -> True` force explicit radicals for the general case.  The special
  always-solvable families — binomials (`a x^n + b`), quadratic-in-`x^m`, and
  biquadratic-after-depression quartics — are always returned in radical form
  regardless of these options, so a matrix whose eigenvalues are compact nested
  radicals keeps its closed form.  `Root[]` objects numericalize (including
  after numeric substitution), so `N` and the numeric optimizers work either
  way.
- For approximate-numeric (`Real` / MPFR) matrices `Method` selects the
  numerical kernel.  Accepted values: `Automatic` (default; routes
  among Direct / Arnoldi / Banded based on shape and `k`), `"Direct"`,
  `"Arnoldi"`, `"Banded"`, and `"FEAST"`.  Each accepts a sub-option
  list form (`Method -> {"Name", "Key" -> value, ...}`); see the
  per-method sections in the `?Eigenvalues` REPL docstring and the
  weekly changelog entries under `docs/spec/changelog/`.

**Attributes:** `Protected`.

## References

**See also:** [Solve](../../solutions-of-equations/Solve/), [N](../../arithmetic/N/)

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the symmetric and unsymmetric eigenvalue problems.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — eigenvalue algorithms and the QR iteration.
- Gene H. Golub, Charles F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins University Press, 2013).
- J. H. Wilkinson, *The Algebraic Eigenvalue Problem* (Oxford University Press, 1965).
- Source: [`src/linalg/eigen.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/eigen.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_lapack_builtin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lapack_builtin.c)
- Tests: [`tests/test_mateigen_arnoldi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_arnoldi.c)
- Tests: [`tests/test_mateigen_banded.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_banded.c)

## Notes & additional examples

### Notes

Eigenvalues of an exact matrix are computed as the roots of the characteristic polynomial `Det[m - lambda I]`, so triangular matrices return their diagonal entries directly and a rotation matrix returns the complex conjugate pair `{-I, I}`. Results are ordered by decreasing absolute value, and repeated eigenvalues appear with their full algebraic multiplicity. Symbolic matrices return closed-form roots (the `2x2` case uses the quadratic formula). Approximate Real / MPFR matrices flow through dedicated numerical kernels (Householder + QR for the symmetric path, Hessenberg + Francis QR for the general path) selectable through the `Method` option.
