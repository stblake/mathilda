# LUDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LUDecomposition[m]
    gives the LU decomposition of a square matrix m as a list
    {lu, p, c}.  The first element lu is the combined Doolittle
    factor matrix: its strictly-lower triangle is L (with an
    implicit unit diagonal) and its upper triangle is U.  The
    second element p is a 1-indexed row-permutation vector such
    that m[[p]] == l . u where l = LowerTriangularize[lu, -1] +
    IdentityMatrix[n] and u = UpperTriangularize[lu].  The third
    element c is an L-infinity condition-number estimate for
    approximate numerical matrices, or 0 for exact / symbolic m.

LUDecomposition works on every input family supported by the
rest of the linear-algebra builtins:
    - exact integer / rational matrices (output stays exact)
    - complex matrices
    - machine-precision Real matrices (LAPACK dgetrf / zgetrf
      with dgecon / zgecon for the condition estimate)
    - arbitrary-precision MPFR matrices (Householder-free
      Doolittle at the input precision; condition estimate via
      the explicit inverse)
    - free-symbolic matrices (output in closed symbolic form)

The algorithm is Doolittle's elimination with partial row
pivoting.  Numerical inputs use largest-|pivot| selection;
symbolic / exact inputs advance to the next non-zero pivot
only when the natural choice is provably zero.

A singular m emits LUDecomposition::sing and the factorisation
completes with a zero on U's diagonal at the singular step.

A non-square or empty matrix emits LUDecomposition::matsq and
the call is left unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LUDecomposition[{{1, 1, 1}, {2, 4, 8}, {3, 9, 27}}]
Out[1]= {{{1, 1, 1}, {2, 2, 6}, {3, 3, 6}}, {1, 2, 3}, 0}

In[2]:= LUDecomposition[{{a, b}, {c, d}}]
Out[2]= {{{a, b}, {c/a, (-b c + a d)/a}}, {1, 2}, 0}
```

## Implementation notes

**Algorithm.** `builtin_ludecomposition` returns `{lu, p, c}` — the combined Doolittle factorisation, the row-permutation list, and a condition-number estimate. A top-level router, `lu_dispatch`, inspects the input with `common_scan_inexact`: machine-precision inexact matrices (`min_bits <= 53`) go to `lu_machine_dispatch` (a LAPACK fast path), higher precision goes to `lu_mpfr_dispatch`, and everything else — including exact integer, rational, complex, and free-symbolic matrices — goes to `lu_symbolic_dispatch`. Any fast-path failure falls through to the symbolic kernel, which absorbs inexact input via the standard `common_rationalize_input` → exact-pipeline → `common_numericalize_result` round-trip.

The symbolic core `lu_symbolic_core` is Doolittle's algorithm (Gaussian elimination producing unit-lower `L` and upper `U` stored in one buffer) with partial pivoting, driven through the Mathilda evaluator so every divide/multiply/subtract works on symbolic, rational, `Sqrt`, or complex entries. Pivoting has two regimes: when a column (from the current row down) is entirely exact-numeric (Integer/BigInt/Rational/Complex of those), `lu_column_all_numeric` is true and the pivot is the entry of *smallest* `|entry|^2` (computed exactly via `numeric_abs_sq_as_mpq`), matching Mathematica's exact-numeric behaviour; otherwise it falls back to "first structurally non-zero" (`is_definitely_zero`, which runs `Together` then `is_zero_poly`). A fully-zero column flags the matrix singular (one-shot `LUDecomposition::sing` warning) but elimination still completes.

**Data structures.** A flat row-major `Expr**` buffer of size `rows*cols` (strict lower triangle = `L`, upper triangle = `U`, unit diagonal implicit), and an `int*` 1-indexed `perm` of length `rows` so `m[[perm]] == l.u`. Exact magnitudes use GMP `mpq_t`. The final `lu` matrix is passed through element-wise `Together` (`tidy_matrix`) to collapse cancellations. The condition slot `c` is exact Integer `0` for the symbolic/exact path (the LAPACK/MPFR kernels supply a real L-infinity estimate).

**Complexity / limits.** O(min(m,n)·m·n) evaluator-level arithmetic operations; symbolic entries can grow in size as elimination proceeds. Rectangular `m × n` input is accepted (elimination stops at step `min(m,n)−1`); only non-list, empty, or higher-rank input is rejected with `LUDecomposition::matsq`.

- `Protected`.
- Any non-empty rectangular `rows x cols` matrix is accepted.  Empty

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — LU factorisation with partial pivoting.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gaussian elimination and LU factorisation.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013).
- Source: [`src/linalg/ludecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/ludecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LUDecomposition[{{4, 3}, {6, 3}}]
Out[1]= {{{4, 3}, {3/2, -3/2}}, {1, 2}, 0}
```

```mathematica
In[1]:= LUDecomposition[{{2, 1}, {4, 1}}]
Out[1]= {{{2, 1}, {2, -1}}, {1, 2}, 0}
```

The combined factor reconstructs `m` via `L . U` (here `L = {{1, 0}, {3/2, 1}}`, `U = {{4, 3}, {0, -3/2}}`):

```mathematica
In[1]:= {{1, 0}, {3/2, 1}} . {{4, 3}, {0, -3/2}}
Out[1]= {{4, 3}, {6, 3}}
```

A non-trivial `3x3` integer matrix factors exactly with no row swap (`p = {1, 2, 3}`), the multipliers `3, 2, 1` packed into the lower triangle:

```mathematica
In[1]:= LUDecomposition[{{1, 2, 4}, {3, 8, 14}, {2, 6, 13}}]
Out[1]= {{{1, 2, 4}, {3, 2, 2}, {2, 1, 3}}, {1, 2, 3}, 0}
```

### Notes

`LUDecomposition` returns a list `{lu, p, c}`: `lu` is the combined Doolittle factor whose strictly-lower triangle is `L` (with an implicit unit diagonal) and whose upper triangle is `U`; `p` is the 1-indexed row-permutation vector (here `{1, 2}`, i.e. no swap was needed); and `c` is an `L`-infinity condition estimate that is `0` for exact or symbolic input. The relation is `m[[p]] == L . U`, as the manual reconstruction above confirms. The algorithm is Doolittle elimination with partial pivoting; exact integer inputs keep exact rational factors. A singular matrix emits `LUDecomposition::sing` and completes with a zero pivot on `U`'s diagonal; a non-square or empty matrix emits `LUDecomposition::matsq`.
