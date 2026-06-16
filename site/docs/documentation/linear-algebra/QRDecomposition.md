# QRDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
QRDecomposition[m]
    gives the QR decomposition of m as a list {q, r}, where q is
    row-orthonormal (row-unitary in the complex case) and r is
    upper triangular.  The original matrix satisfies
    m == ConjugateTranspose[q] . r.

QRDecomposition computes the "thin" QR factorisation: when m
has rank r, both q and r have r rows.  For an n x p input,
q has dimensions r x n and r has dimensions r x p, so q's
rows live in the column space of m and r encodes the original
columns in that basis.

QRDecomposition[m, Pivoting -> True]
    gives a list {q, r, p} where p is a p x p permutation matrix
    such that m . p == ConjugateTranspose[q] . r.  With pivoting
    the diagonal of r appears in order of decreasing magnitude.

QRDecomposition works on every input family supported by the
rest of the linear-algebra builtins:
    - exact integer / rational matrices (output stays exact,
      with Sqrt[...] in the column norms)
    - complex matrices (q's rows are unitary in the Hermitian
      inner product)
    - machine-precision Real matrices (output is Real at machine
      precision, matching the inexact-in / inexact-out contract)
    - arbitrary-precision MPFR matrices (output at the input
      precision)
    - free-symbolic matrices (output in closed symbolic form)

The algorithm is Modified Gram-Schmidt on the columns of m,
applied through the evaluator so symbolic, exact, and inexact
inputs share one code path.  Rank-deficient inputs (columns in
the span of earlier columns) produce a shorter q / r without
any error.

A non-rank-2 or empty matrix emits QRDecomposition::matrix and
the call is left unevaluated.  Unknown option keys or values
emit QRDecomposition::opts and the call is left unevaluated.
TargetStructure -> "Structured" is reserved for a future
release and currently leaves the call unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_qrdecomposition` returns the thin QR `{q, r}` (or `{q, r, p}` with `Pivoting -> True`) satisfying `m == ConjugateTranspose[q].r`. `qr_dispatch` routes inexact `min_bits <= 53` input to `qr_machine_dispatch` (LAPACK `dgeqrf`/`dgeqp3`, or `zgeqrf`/`zgeqp3` for complex, with the standard rank-revealing cutoff), higher precision to `qr_mpfr_dispatch` (textbook Householder QR over MPFR), and exact integer/rational/complex/symbolic input to `qr_symbolic_dispatch`. Any fast-path failure (LAPACK off, non-coercible leaf, non-zero `info`, rank-deficient with no pivoting) falls through to the symbolic kernel.

The symbolic core `qr_symbolic_core` is **Modified Gram-Schmidt** on the columns of `A`, evaluated symbolically: for each column it subtracts the projections `coeff = <Q[:,j], v>` onto each existing orthonormal column (storing `coeff` into `R[j,k]`), takes `norm = Sqrt[<v,v>]`, and appends `v/norm` as a new orthonormal column; a zero residual norm marks a dependent column that is skipped (so `r = rank`). `q` is built as the conjugate transpose of `Q` so the result is `Q^T` for real input and the proper Hermitian transpose for complex. With `Pivoting -> True` the next column chosen at each step is the remaining one whose residual projection has the largest squared norm — Businger-Golub column pivoting expressed in MGS form, making the diagonal of `R` decrease in magnitude — and the order is inflated into a `p × p` permutation matrix with `p[perm[j], j] = 1`.

**Data structures.** Flat row-major `Expr**` working buffers for `Q` (n × r) and `R` (r × p); inner products, norms, and scalings go through the evaluator via `eval_and_free`, leaving exact `Sqrt[...]` / rational / symbolic entries. The `is_definitely_zero` test (`Together` then `is_zero_poly`) detects dependent columns. Inexact input uses the `common_rationalize_input` / `common_numericalize_result` round-trip at the minimum input precision.

**Complexity / limits.** O(n·p·r) evaluator-level vector operations for the symbolic path; symbolic norms accumulate `Sqrt` nestings. The machine path is LAPACK-bound.

- `Protected`.
- Computes the "thin" QR factorisation: when `m` has rank `r`, both

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gram-Schmidt orthogonalisation and the QR factorisation.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — QR factorisation algorithms.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013).
- P. A. Businger and G. H. Golub, "Linear Least Squares Solutions by Householder Transformations", Numer. Math. 7 (1965).
- Source: [`src/linalg/qrdecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/qrdecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= QRDecomposition[{{3, 0}, {0, 4}}]
Out[1]= {{{1, 0}, {0, 1}}, {{3, 0}, {0, 4}}}
```

```mathematica
In[1]:= QRDecomposition[{{1, 1}, {0, 1}}]
Out[1]= {{{1, 0}, {0, 1}}, {{1, 1}, {0, 1}}}
```

Reconstructing the original matrix from `{q, r}` uses `ConjugateTranspose[q] . r`:

```mathematica
In[1]:= q = QRDecomposition[{{1, 1}, {0, 1}}][[1]]; r = QRDecomposition[{{1, 1}, {0, 1}}][[2]]; ConjugateTranspose[q] . r
Out[1]= {{1, 1}, {0, 1}}
```

For a non-orthogonal integer matrix the factors stay exact, with `Sqrt[...]` column norms carried symbolically (the classic Householder textbook matrix):

```mathematica
In[1]:= QRDecomposition[{{12, -51}, {6, 167}}]
Out[1]= {{{2/Sqrt[5], 1/Sqrt[5]}, {-1/Sqrt[5], 2/Sqrt[5]}}, {{6 Sqrt[5], 13 Sqrt[5]}, {0, 77 Sqrt[5]}}}
```

The rows of `q` are genuinely orthonormal, so `q . ConjugateTranspose[q]` is the identity even though every entry is an irrational radical:

```mathematica
In[1]:= q = QRDecomposition[{{2, -1}, {1, 2}}][[1]]; q . ConjugateTranspose[q]
Out[1]= {{1, 0}, {0, 1}}
```

With `Pivoting -> True` the columns are reordered by decreasing norm and a permutation matrix `p` is returned alongside `{q, r}`:

```mathematica
In[1]:= QRDecomposition[{{1, 2}, {3, 4}}, Pivoting -> True]
Out[1]= {{{1/Sqrt[5], 2/Sqrt[5]}, {-2/Sqrt[5], 1/Sqrt[5]}}, {{2 Sqrt[5], 7/Sqrt[5]}, {0, 1/Sqrt[5]}}, {{0, 1}, {1, 0}}}
```

### Notes

`QRDecomposition` returns the **thin** factorisation `{q, r}` where `q` is row-orthonormal and `r` is upper triangular, with the original matrix recovered as `ConjugateTranspose[q] . r` (note the conjugate-transpose convention: `q`'s *rows* are the orthonormal basis). When the columns are already orthogonal — as for a diagonal matrix or the unit-column examples above — `q` is the identity and `r` reproduces the input. The algorithm is Modified Gram-Schmidt applied through the evaluator, so exact integer inputs keep `Sqrt[...]` column norms in closed form while machine-precision Real inputs return Real factors. `Pivoting -> True` adds a permutation matrix `p` so that `m . p == ConjugateTranspose[q] . r`.
