# LUDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LUDecomposition[m]`**

gives the LU decomposition of a square matrix m as a list {lu, p, c}.  The first element lu is the combined Doolittle factor matrix: its strictly-lower triangle is L (with an implicit unit diagonal) and its upper triangle is U.  The second element p is a 1-indexed row-permutation vector such that m\[\[p\]\] == l . u where l = LowerTriangularize\[lu, -1\] + IdentityMatrix\[n\] and u = UpperTriangularize\[lu\].  The third element c is an L-infinity condition-number estimate for approximate numerical matrices, or 0 for exact / symbolic m.

<details>
<summary>Notes</summary>

LUDecomposition works on every input family supported by the rest of the linear-algebra builtins: - exact integer / rational matrices (output stays exact) - complex matrices - machine-precision Real matrices (LAPACK dgetrf / zgetrf with dgecon / zgecon for the condition estimate) - arbitrary-precision MPFR matrices (Householder-free Doolittle at the input precision; condition estimate via the explicit inverse) - free-symbolic matrices (output in closed symbolic form) The algorithm is Doolittle's elimination with partial row pivoting.  Numerical inputs use largest-|pivot| selection; symbolic / exact inputs advance to the next non-zero pivot only when the natural choice is provably zero. A singular m emits LUDecomposition::sing and the factorisation completes with a zero on U's diagonal at the singular step. A non-square or empty matrix emits LUDecomposition::matsq and the call is left unevaluated.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= LUDecomposition[{{1, 1, 1}, {2, 4, 8}, {3, 9, 27}}]
Out[1]= {{{1, 1, 1}, {2, 2, 6}, {3, 3, 6}}, {1, 2, 3}, 0}

In[2]:= LUDecomposition[{{a, b}, {c, d}}]
Out[2]= {{{a, b}, {c/a, (-b c + a d)/a}}, {1, 2}, 0}

In[3]:= {lu, p, c} = LUDecomposition[{{1.6, 2.7, 3.6}, {1.2, 3.2, 5.2}, {3.3, 3.4, 6.5}}]; c
Out[3]= 20.8391
```

### Applications (4)

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

## Options & behaviour

> Implementation is split across:
> - `src/linalg/ludecomp.c` -- builtin entry, top-level dispatcher,
>   and the symbolic Doolittle core driven through the evaluator.
> - `src/linalg/ludecomp_machine.c` -- LAPACK fast path; loads to a
>   column-major double buffer, runs `dgetrf` / `zgetrf` then
>   `dgecon` / `zgecon`, builds `{lu, p, c}` as Mathilda lists.
> - `src/linalg/ludecomp_mpfr.c` -- MPFR Doolittle kernel; row-major
>   MPFR arrays (paired re/im for complex), largest-magnitude pivot
>   selection, ‖A‖∞ * ‖A⁻¹‖∞ condition estimate.

## Algorithm

ludecomp.c

LUDecomposition[m] -- {lu, p, c} Doolittle factorisation with row pivoting.

```text
Strategy.  One algorithmic core - Doolittle's algorithm with partial
```

pivoting driven through the Mathilda evaluator - serves every input family:

```text
  - exact integer / rational / complex / free-symbolic matrices
    run the pipeline as-is.  Pivoting is zero-only (advance to the
    next non-zero pivot if the natural choice is provably zero).

  - inexact matrices (Real or MPFR leaves) follow the
    rationalise -> exact pipeline -> numericalise round-trip used by
    PseudoInverse / Eigenvalues / QRDecomposition.  Output precision
    matches the minimum precision present in the input.

Conventions.  We work internally on a flat n x n row-major buffer:

    LU[i * n + k] = (i, k) entry of the combined factorisation
```

with the strict lower triangle storing L (unit diagonal implicit)

```text
and the upper triangle storing U.  perm[k] is the 1-indexed original
```

row that was placed at row k by the pivoting steps, so the public identity m[[perm]] == l . u holds:

```text
    l = LowerTriangularize[lu, -1] + IdentityMatrix[n]
    u = UpperTriangularize[lu]
```

For exact / symbolic inputs the condition-number slot is the exact

```text
Integer 0 (matching the Mathematica example).  The machine and MPFR
```

kernels emit a Real / MPFR L-infinity estimate; in this top-level file we only need the exact-zero default.

```text
Memory contract.  Standard builtin contract.  This file does NOT
```

call expr_free(res) - the evaluator owns `res` and frees it on a

```text
non-NULL return (MEMORY.md / SPEC.md §4.1).  Every intermediate
```

allocation is matched by a free along every exit path.

## Implementation notes

**Algorithm.** `builtin_ludecomposition` returns `{lu, p, c}` — the combined Doolittle factorisation, the row-permutation list, and a condition-number estimate. A top-level router, `lu_dispatch`, inspects the input with `common_scan_inexact`: machine-precision inexact matrices (`min_bits <= 53`) go to `lu_machine_dispatch` (a LAPACK fast path), higher precision goes to `lu_mpfr_dispatch`, and everything else — including exact integer, rational, complex, and free-symbolic matrices — goes to `lu_symbolic_dispatch`. Any fast-path failure falls through to the symbolic kernel, which absorbs inexact input via the standard `common_rationalize_input` → exact-pipeline → `common_numericalize_result` round-trip.

The symbolic core `lu_symbolic_core` is Doolittle's algorithm (Gaussian elimination producing unit-lower `L` and upper `U` stored in one buffer) with partial pivoting, driven through the Mathilda evaluator so every divide/multiply/subtract works on symbolic, rational, `Sqrt`, or complex entries. Pivoting has two regimes: when a column (from the current row down) is entirely exact-numeric (Integer/BigInt/Rational/Complex of those), `lu_column_all_numeric` is true and the pivot is the entry of *smallest* `|entry|^2` (computed exactly via `numeric_abs_sq_as_mpq`), following the exact-numeric convention; otherwise it falls back to "first structurally non-zero" (`is_definitely_zero`, which runs `Together` then `is_zero_poly`). A fully-zero column flags the matrix singular (one-shot `LUDecomposition::sing` warning) but elimination still completes.

**Data structures.** A flat row-major `Expr**` buffer of size `rows*cols` (strict lower triangle = `L`, upper triangle = `U`, unit diagonal implicit), and an `int*` 1-indexed `perm` of length `rows` so `m[[perm]] == l.u`. Exact magnitudes use GMP `mpq_t`. The final `lu` matrix is passed through element-wise `Together` (`tidy_matrix`) to collapse cancellations. The condition slot `c` is exact Integer `0` for the symbolic/exact path (the LAPACK/MPFR kernels supply a real L-infinity estimate).

**Complexity / limits.** O(min(m,n)·m·n) evaluator-level arithmetic operations; symbolic entries can grow in size as elimination proceeds. Rectangular `m × n` input is accepted (elimination stops at step `min(m,n)−1`); only non-list, empty, or higher-rank input is rejected with `LUDecomposition::matsq`.

- `Protected`.
- Any non-empty rectangular `rows x cols` matrix is accepted.  Empty
  matrices and non-matrix arguments emit `LUDecomposition::matsq` and
  the call is left unevaluated.
- Algorithm: Doolittle's elimination with partial row pivoting.
  - **MachinePrecision inputs (`min_bits <= 53`)** use the LAPACK
    fast path: `dgetrf` / `zgetrf` for the factorisation, plus
    `dgecon` / `zgecon` (with `dlange` / `zlange` for ‖A‖∞) for the
    condition estimate.
  - **MPFR inputs (`min_bits > 53`)** use a hand-rolled Doolittle
    kernel over row-major MPFR arrays (paired re/im for complex).
    For real matrices the condition number is estimated by the
    Hager-Higham one-norm iteration on `A^{-T}` (LAPACK's `*lacn2`
    strategy; typically 2-5 triangular-solve pairs, each `O(n^2)`).
    For complex matrices the kernel falls back to the explicit
    inverse via `n` back-substitution pairs (`O(n^3)`).
  - **Exact / symbolic inputs** run Doolittle through the Mathilda
    evaluator, so integer / rational / Complex / Sqrt-bearing /
    free-symbolic entries all share one code path.  Pivot selection:
    when every entry of the active column is an exact numeric
    (`Integer` / `BigInt` / `Rational` / `Complex` of those) the pivot
    with the **smallest absolute value** is chosen — matching
    Mathematica (e.g. `LUDecomposition[{{1/2, 1/3}, {1/5, 1/7}}]`
    picks the `1/5` pivot, keeping intermediate `L` entries small).
    For any column containing a free symbol, `Sqrt`, or other
    non-exact-numeric leaf, the rule falls back to "first non-zero"
    — matching the Mathematica example
    `LUDecomposition[{{a, b}, {c, d}}]` returning `p = {1, 2}`.
  - When BLAS/LAPACK is unavailable at build time (`USE_LAPACK=0`)
    or MPFR is unavailable (`USE_MPFR=0`) the corresponding fast
    path transparently routes to the symbolic kernel, which itself
    understands inexact input via the standard rationalise → exact
    → numericalise round-trip.
- Singular inputs emit `LUDecomposition::sing` and the factorisation
  completes with a zero on U's diagonal at the singular step
  (matching Mathematica's behaviour).
- Ill-conditioned numerical inputs emit `LUDecomposition::luc`: the
  factorisation completes but the reported `c` exceeds the
  precision-loss threshold (`1 / $MachineEpsilon` for machine input,
  `2^bits` for MPFR input), matching Mathematica's
  `LUDecomposition::luc` behaviour.  Both warnings are one-shot per
  process; subsequent calls are silent.

**Attributes:** `Protected`.

## See also

[Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [Sqrt](../../arithmetic/Sqrt/)

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — LU factorisation with partial pivoting.
- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — Gaussian elimination and LU factorisation.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Johns Hopkins, 2013).
- Source: [`src/linalg/ludecomp.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/ludecomp.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_ludecomposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ludecomposition.c)
- Tests: [`tests/test_ludecomposition_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ludecomposition_machine.c)
- Tests: [`tests/test_ludecomposition_mpfr.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ludecomposition_mpfr.c)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)

## Notes & additional examples

### Notes

`LUDecomposition` returns a list `{lu, p, c}`: `lu` is the combined Doolittle factor whose strictly-lower triangle is `L` (with an implicit unit diagonal) and whose upper triangle is `U`; `p` is the 1-indexed row-permutation vector (here `{1, 2}`, i.e. no swap was needed); and `c` is an `L`-infinity condition estimate that is `0` for exact or symbolic input. The relation is `m[[p]] == L . U`, as the manual reconstruction above confirms. The algorithm is Doolittle elimination with partial pivoting; exact integer inputs keep exact rational factors. A singular matrix emits `LUDecomposition::sing` and completes with a zero pivot on `U`'s diagonal; a non-square or empty matrix emits `LUDecomposition::matsq`.
