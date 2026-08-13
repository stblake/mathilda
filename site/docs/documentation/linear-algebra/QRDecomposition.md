# QRDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`QRDecomposition[m]`**

gives the QR decomposition of m as a list {q, r}, where q is row-orthonormal (row-unitary in the complex case) and r is upper triangular.  The original matrix satisfies m == ConjugateTranspose\[q\] . r.

**`QRDecomposition[m, Pivoting -> True]`**

gives a list {q, r, p} where p is a p x p permutation matrix such that m . p == ConjugateTranspose\[q\] . r.  With pivoting the diagonal of r appears in order of decreasing magnitude.

<details>
<summary>Notes</summary>

QRDecomposition computes the "thin" QR factorisation: when m has rank r, both q and r have r rows.  For an n x p input, q has dimensions r x n and r has dimensions r x p, so q's rows live in the column space of m and r encodes the original columns in that basis. QRDecomposition works on every input family supported by the rest of the linear-algebra builtins: - exact integer / rational matrices (output stays exact, with Sqrt\[...\] in the column norms) - complex matrices (q's rows are unitary in the Hermitian inner product) - machine-precision Real matrices (output is Real at machine precision, matching the inexact-in / inexact-out contract) - arbitrary-precision MPFR matrices (output at the input precision) - free-symbolic matrices (output in closed symbolic form) The algorithm is Modified Gram-Schmidt on the columns of m, applied through the evaluator so symbolic, exact, and inexact inputs share one code path.  Rank-deficient inputs (columns in the span of earlier columns) produce a shorter q / r without any error. A non-rank-2 or empty matrix emits QRDecomposition::matrix and the call is left unevaluated.  Unknown option keys or values emit QRDecomposition::opts and the call is left unevaluated. TargetStructure -\> "Structured" is reserved for a future release and currently leaves the call unevaluated.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= {q, r} = QRDecomposition[{{1, 2}, {3, 4}, {5, 6}}]; Transpose[q] . r
Out[1]= {{1, 26/35 + (2 Sqrt[7/5] + 6 Sqrt[5/7])/Sqrt[35]}, {3, 8/35 + 3 (2 Sqrt[7/5] + 6 Sqrt[5/7])/Sqrt[35]}, {5, -2/7 + 5 (2 Sqrt[7/5] + 6 Sqrt[5/7])/Sqrt[35]}}

In[2]:= {q, r} = QRDecomposition[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]; {Length[q], Length[r]}
Out[2]= {2, 2}

In[3]:= {q, r} = QRDecomposition[{{1.2, 2.3, 3.4}, {2.3, 4.5, 5.6}, {3.2, 7.6, 6.5}}]; Chop[Transpose[q] . r - {{1.2, 2.3, 3.4}, {2.3, 4.5, 5.6}, {3.2, 7.6, 6.5}}]
Out[3]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
```

### Options (1)

```mathematica
In[4]:= QRDecomposition[{{1, 2}, {3, 4}}, Pivoting -> True]
Out[4]= {{{1/Sqrt[5], 2/Sqrt[5]}, {-2/Sqrt[5], 1/Sqrt[5]}}, {{2 Sqrt[5], 7/Sqrt[5]}, {0, 1/Sqrt[5]}}, {{0, 1}, {1, 0}}}
```

### Applications (6)

```mathematica
In[5]:= QRDecomposition[{{3, 0}, {0, 4}}]
Out[5]= {{{1, 0}, {0, 1}}, {{3, 0}, {0, 4}}}

In[6]:= QRDecomposition[{{1, 1}, {0, 1}}]
Out[6]= {{{1, 0}, {0, 1}}, {{1, 1}, {0, 1}}}

In[7]:= q = QRDecomposition[{{1, 1}, {0, 1}}][[1]]; r = QRDecomposition[{{1, 1}, {0, 1}}][[2]]; ConjugateTranspose[q] . r
Out[7]= {{1, 1}, {0, 1}}

In[8]:= QRDecomposition[{{12, -51}, {6, 167}}]
Out[8]= {{{2/Sqrt[5], 1/Sqrt[5]}, {-1/Sqrt[5], 2/Sqrt[5]}}, {{6 Sqrt[5], 13 Sqrt[5]}, {0, 77 Sqrt[5]}}}

In[9]:= q = QRDecomposition[{{2, -1}, {1, 2}}][[1]]; q . ConjugateTranspose[q]
Out[9]= {{1, 0}, {0, 1}}

In[10]:= QRDecomposition[{{1, 2}, {3, 4}}, Pivoting -> True]
Out[10]= {{{1/Sqrt[5], 2/Sqrt[5]}, {-2/Sqrt[5], 1/Sqrt[5]}}, {{2 Sqrt[5], 7/Sqrt[5]}, {0, 1/Sqrt[5]}}, {{0, 1}, {1, 0}}}
```

## Options & behaviour

> Implementation is split across:
> - `src/linalg/qrdecomp.c` -- builtin entry, option parsing, and
>   `qr_dispatch` which routes to the precision-matched kernel.
>   Hosts the symbolic / fallback MGS pipeline.
> - `src/linalg/qrdecomp_machine.c` -- LAPACK fast path
>   (`qr_machine_dispatch`).  Loads the matrix into a column-major
>   double buffer (interleaved re/im pairs for complex), calls the
>   wrappers in `lapack.c` (`mat_lapack_dgeqp3`, etc.), then
>   reconstructs `q` / `r` / `p` as Mathilda lists.  Numerical rank
>   uses LAPACK's standard cutoff `max(m, n) * eps * |R[0,0]|`.
> - `src/linalg/qrdecomp_mpfr.c` -- MPFR Householder fast path
>   (`qr_mpfr_dispatch`).  Loads the matrix into column-major MPFR
>   arrays at `min_bits` precision (paired re/im planes for complex),
>   runs Householder reflections in place with Businger-Golub
>   pivoting, then reconstructs `q` / `r` / `p` as MPFR-precision
>   Mathilda lists.  Updates already-stored R rows in-step with
>   column swaps so R's column ordering stays consistent with the
>   pivoted A.  Reconstruction residuals scale as `2^(-bits)`.
> - `src/linalg/lapack.h` / `lapack.c` -- platform-papering Fortran
>   ABI wrappers shared across machine-precision linalg kernels.
>
> The MGS loop allocates a column-major Q buffer and a row-major R
> buffer of size `min(n,p) x max(n,p)`, frees the unused tail when
> the rank turns out to be smaller, and steals the in-use cells
> into the final `List[List[...]]` result.  For complex inputs the
> `q` entries are conjugated at construction time; real inputs (no
> Complex head, no `I` leaf) skip the conjugation to keep the
> printed form free of `Conjugate[Sqrt[...]]` residues that
> Mathilda's simplifier does not reduce.

## Algorithm

qrdecomp.c

```text
QRDecomposition[m]                   -- {q, r} thin-QR factorisation.
```

QRDecomposition[m, Pivoting -> True] -- {q, r, p}, m . p == q^H . r.

```text
Strategy.  One algorithmic core - Modified Gram-Schmidt on the
```

columns of m, driven through the Mathilda evaluator - serves every input family:

```text
  - exact integer / rational / complex / free-symbolic matrices
    run the pipeline as-is.  The output is exact (Sqrt[...] in
    the norms, Rational / symbolic entries elsewhere).

  - inexact matrices (Real or MPFR leaves) follow the
    rationalise -> exact pipeline -> numericalise round-trip used
    by PseudoInverse / Eigenvalues / Solve.  The output precision
    matches the minimum precision present in the input, mirroring
    the inexact-in / inexact-out contract advertised across the
    rest of the system.

Conventions.  We work internally with a standard "thin" QR

    A = Q . R          Q n x r orthonormal-columns,  R r x p upper
                       trapezoidal,  r = MatrixRank[A]

and at the end return  q = ConjugateTranspose[Q],  r = R  so the
Mathematica identity  m == ConjugateTranspose[q] . r  holds.
```

Because ConjugateTranspose is involutive this matches the spec convention exactly: Length[q] == Length[r] == r (the rank), the rows of q are orthonormal in the complex inner product, and r has zeros below the leading-diagonal echelon.

Modified Gram-Schmidt loop, column k = 0 .. p-1:

```text
    v = A[:, k]
    for each existing orthonormal column Q[:, j]:
        coeff = <Q[:, j], v> = Sum_i Conjugate[Q[i, j]] * v[i]
        R[j, k] = coeff
        v -= coeff * Q[:, j]
    norm_sq = <v, v>
    if norm_sq == 0: column is dependent, skip (no new orthonormal row)
    norm = Sqrt[norm_sq]
    R[rk, k] = norm
    Q[:, rk] = v / norm
    rk += 1

After the loop q is built as  q[j, i] = Conjugate[Q[i, j]] - this
```

collapses to Q^T for real matrices and gives the proper conjugate-transpose for complex matrices.

```text
Pivoting (when Pivoting -> True).  At the start of each step we
```

pick, among the remaining columns of A, the one whose residual orthogonal projection (after subtracting components along the

```text
already-built Q columns) has the largest squared norm.  This is
```

exactly Householder column pivoting expressed in MGS form and makes the diagonal of R appear in order of decreasing magnitude,

```text
matching the Mathematica example.  The permutation array is then
```

inflated into a p x p permutation matrix p with

```text
    p[perm[j], j] = 1
```

so that A . p picks the columns in the chosen order, satisfying A . p == ConjugateTranspose[q] . r.

```text
Memory contract.  Standard builtin contract.  This file does NOT
```

call expr_free(res) - the evaluator owns `res` and frees it on a

```text
non-NULL return (MEMORY.md / SPEC.md §4.1).  Every intermediate
```

allocation is tracked: the Q and R working buffers are freed after the q / r List wrappers have stolen / copied their entries.

## Implementation notes

**Algorithm.** `builtin_qrdecomposition` returns the thin QR `{q, r}` (or `{q, r, p}` with `Pivoting -> True`) satisfying `m == ConjugateTranspose[q].r`. `qr_dispatch` routes inexact `min_bits <= 53` input to `qr_machine_dispatch` (LAPACK `dgeqrf`/`dgeqp3`, or `zgeqrf`/`zgeqp3` for complex, with the standard rank-revealing cutoff), higher precision to `qr_mpfr_dispatch` (textbook Householder QR over MPFR), and exact integer/rational/complex/symbolic input to `qr_symbolic_dispatch`. Any fast-path failure (LAPACK off, non-coercible leaf, non-zero `info`, rank-deficient with no pivoting) falls through to the symbolic kernel.

The symbolic core `qr_symbolic_core` is **Modified Gram-Schmidt** on the columns of `A`, evaluated symbolically: for each column it subtracts the projections `coeff = <Q[:,j], v>` onto each existing orthonormal column (storing `coeff` into `R[j,k]`), takes `norm = Sqrt[<v,v>]`, and appends `v/norm` as a new orthonormal column; a zero residual norm marks a dependent column that is skipped (so `r = rank`). `q` is built as the conjugate transpose of `Q` so the result is `Q^T` for real input and the proper Hermitian transpose for complex. With `Pivoting -> True` the next column chosen at each step is the remaining one whose residual projection has the largest squared norm — Businger-Golub column pivoting expressed in MGS form, making the diagonal of `R` decrease in magnitude — and the order is inflated into a `p × p` permutation matrix with `p[perm[j], j] = 1`.

**Data structures.** Flat row-major `Expr**` working buffers for `Q` (n × r) and `R` (r × p); inner products, norms, and scalings go through the evaluator via `eval_and_free`, leaving exact `Sqrt[...]` / rational / symbolic entries. The `is_definitely_zero` test (`Together` then `is_zero_poly`) detects dependent columns. Inexact input uses the `common_rationalize_input` / `common_numericalize_result` round-trip at the minimum input precision.

**Complexity / limits.** O(n·p·r) evaluator-level vector operations for the symbolic path; symbolic norms accumulate `Sqrt` nestings. The machine path is LAPACK-bound.

- `Protected`.
- Computes the "thin" QR factorisation: when `m` has rank `r`, both
  `q` and `r` have `r` rows. For an `n x p` input, `q` has dimensions
  `r x n` and `r` has dimensions `r x p`.
- Works on every input family:
  - exact integer / rational matrices (output stays exact with
    `Sqrt[...]` in the column norms);
  - complex matrices (rows of `q` are unitary in the Hermitian
    inner product);
  - machine-precision Real matrices (output Real at machine
    precision);
  - arbitrary-precision MPFR matrices (output at the input
    precision via the shared rationalise → exact → numericalise
    pipeline);
  - free-symbolic matrices (closed-form symbolic output).
- Algorithm: dispatched on leaf precision.
  - **MachinePrecision inputs (`min_bits <= 53`)** use a LAPACK
    Householder kernel (`dgeqrf` / `dgeqp3` for real, `zgeqrf` /
    `zgeqp3` for complex, plus `dorgqr` / `zungqr` to form `q`).
    Wired through the four-tier autodetection ladder described in
    `src/linalg/lapack.h` (Apple Accelerate → pkg-config lapacke →
    system lapacke → graceful fall-back).
  - **MPFR inputs (`min_bits > 53`)** use a hand-rolled Householder
    kernel over column-major MPFR arrays (paired re/im planes for
    complex; no MPC dependency, same convention as the eigen
    kernels).  Column pivoting follows Businger-Golub; numerical
    rank uses the cutoff `|R[i,i]| < 2^(-bits/2) * |R[0,0]|`.
    Reconstruction error scales as `2^(-bits)`, matching the
    requested working precision.
  - **Exact / symbolic inputs** stay on the Modified Gram-Schmidt
    kernel, driven through the evaluator so symbolic-real, exact
    rational, and free-variable inputs share one code path.
  - Rank-deficient inputs produce a shorter `q` / `r` without error.
    With `Pivoting -> False` on a rank-deficient input the MPFR
    kernel bails to symbolic (which handles mid-stream rank
    deficiency cleanly); with `Pivoting -> True` the MPFR kernel
    truncates the output at the numerical rank.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`),
    machine-precision inputs transparently route to the symbolic
    kernel; similarly `USE_MPFR=0` routes MPFR inputs to symbolic.
    Correctness is preserved across every combination; only
    performance changes.
- Issues `QRDecomposition::matrix` and returns unevaluated if the
  argument is not a non-empty rank-2 tensor.
- Issues `QRDecomposition::opts` and returns unevaluated for an
  unknown option key or value. `TargetStructure -> "Structured"` is
  reserved for a future release and currently leaves the call
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [I](../../mathematical-constants/I/)

- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gram-Schmidt orthogonalisation and the QR factorisation.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — QR factorisation algorithms.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013).
- P. A. Businger and G. H. Golub, "Linear Least Squares Solutions by Householder Transformations", Numer. Math. 7 (1965).
- Source: [`src/linalg/qrdecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/qrdecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)
- Tests: [`tests/test_qrdecomposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_qrdecomposition.c)
- Tests: [`tests/test_qrdecomposition_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_qrdecomposition_machine.c)
- Tests: [`tests/test_qrdecomposition_mpfr.c`](https://github.com/stblake/mathilda/blob/main/tests/test_qrdecomposition_mpfr.c)

## Notes & additional examples

### Notes

`QRDecomposition` returns the **thin** factorisation `{q, r}` where `q` is row-orthonormal and `r` is upper triangular, with the original matrix recovered as `ConjugateTranspose[q] . r` (note the conjugate-transpose convention: `q`'s *rows* are the orthonormal basis). When the columns are already orthogonal — as for a diagonal matrix or the unit-column examples above — `q` is the identity and `r` reproduces the input. The algorithm is Modified Gram-Schmidt applied through the evaluator, so exact integer inputs keep `Sqrt[...]` column norms in closed form while machine-precision Real inputs return Real factors. `Pivoting -> True` adds a permutation matrix `p` so that `m . p == ConjugateTranspose[q] . r`.
