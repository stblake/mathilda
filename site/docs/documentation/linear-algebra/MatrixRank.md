# MatrixRank

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MatrixRank[m]`**

gives the rank of the matrix m -- the number of linearly independent rows (equivalently, of linearly independent columns).

**`MatrixRank[m, Method -> "<name>"]`**

runs a specific elimination algorithm for the exact path. Accepted method names match NullSpace / RowReduce / LinearSolve / Inverse: "Automatic"                 -- alias for "DivisionFreeRowReduction" (default) "DivisionFreeRowReduction"  -- Bareiss-like fraction-free Gauss-Jordan "OneStepRowReduction"       -- classical Gauss-Jordan with division per pivot "CofactorExpansion"         -- identity-if-invertible (falls back to DivisionFreeRowReduction on singular / rectangular m)

**`MatrixRank[m, Tolerance -> t]`**

treats |entry| \<= t as zero during pivot selection.  With Tolerance -\> 0 even arbitrarily small entries count; the default, Tolerance -\> Automatic, applies max(rows, cols) \* MachineEpsilon \* Max\[|entries|\] for machine-precision (Real / MPFR) matrices and 0 otherwise.

<details>
<summary>Notes</summary>

MatrixRank works on both numerical and symbolic matrices and on square or rectangular matrices.  The default exact path routes through RowReduce and counts the non-zero rows; the numerical path (triggered by inexact leaves or an explicit Tolerance) runs partial-pivot Gaussian elimination over double-precision complex. An unknown Method value or Tolerance form emits MatrixRank::opt and leaves the call unevaluated.  A non-rank-2 or empty matrix emits MatrixRank::matrix and the call is left unevaluated.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 2

In[2]:= MatrixRank[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[2]= 3

In[3]:= MatrixRank[{{0, 5, 2, 4, 4}, {2, 5, 0, 4, 0}, {5, 1, 5, 4, 5}}]
Out[3]= 3

In[4]:= MatrixRank[{{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {1.1, 2.5, -1.5}}]
Out[4]= 3

In[5]:= MatrixRank[{{a, b}, {2 a, 2 b}}]
Out[5]= 1

In[6]:= m = {{1, 1, 1}, {0, 10^-10, 0}, {0, 0, 10^-20}}; MatrixRank[m]
Out[6]= 3

In[7]:= MatrixRank[N[m]]
Out[7]= 2
```

### Options (1)

```mathematica
In[8]:= MatrixRank[N[m], Tolerance -> 0]
Out[8]= 3
```

### Applications (5)

```mathematica
In[1]:= MatrixRank[{{1, 2}, {2, 4}}]
Out[1]= 1
```

```mathematica
In[1]:= MatrixRank[{{1, 2}, {3, 4}}]
Out[1]= 2
```

```mathematica
In[1]:= MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 2
```

```mathematica
In[1]:= MatrixRank[{{a, b}, {2 a, 2 b}}]
Out[1]= 1
```

```mathematica
In[1]:= MatrixRank[Table[1/(i + j - 1), {i, 4}, {j, 4}]]
Out[1]= 4
```

## Options & behaviour

> Implementation lives in `src/linalg/matrank.c` (registered by
> `matrank_init`). The exact path delegates to RowReduce, so any
> improvement to RowReduce / its `Method` kernels propagates to
> MatrixRank.  The numerical path is self-contained: a portable
> `cplx_t = {re, im}` struct stands in for `double _Complex` to keep
> the build strictly C99.

## Algorithm

matrank.c

```text
MatrixRank[m]                       -- rank of m (exact path through
                                      RowReduce when m has no Real /
                                      MPFR / Complex-of-Real entries).
MatrixRank[m, Method  -> "<name>"]  -- explicit RREF method dispatch
                                      (exact path only).
MatrixRank[m, Tolerance -> t]       -- treat |entry| <= t as zero
                                      during pivot selection.
```

MatrixRank[m, Method->..., Tolerance->...]

```text
                                    -- both options simultaneously.
```

Algorithm dispatch:

```text
  - Numerical path (matrix has any Real / MPFR leaf, or every leaf
    converts to a complex double, OR the user supplies Tolerance):
    run a partial-pivot Gaussian forward-elimination over
    `double _Complex` (here modelled as a {re, im} struct so we stay
    portable to compilers without _Complex).  A column is skipped
    when its largest sub-pivot |entry| is <= tolerance; the rank is
    the number of accepted pivots.

  - Exact path (every leaf is exact AND no Tolerance supplied):
    route through `RowReduce[m, Method->...]` and count the number
    of RREF rows whose entries are not all structurally zero.  This
    gives the rank with no precision concerns.
```

Default tolerance (when Tolerance -> Automatic, the default):

```text
  - If the matrix contains any Real / MPFR leaf:
      tol = max(rows, cols) * DBL_EPSILON * max(|entries|)
    (the standard "rank-by-SVD" surrogate; we substitute max(|entries|)
    for the largest singular value, which is exact in the row /
    column-scaled cases and at most a small constant factor off in
    general).
  - Otherwise tol = 0, so the numerical path agrees with the exact
    path on integer / rational matrices.

Memory ownership: standard builtin contract.  This file does NOT
```

free `res` on success or failure -- the evaluator owns it (see MEMORY.md / SPEC.md §4.1).

## Implementation notes

**Algorithm.** `builtin_matrixrank` returns the rank as the number of pivots, choosing between two paths. The **numerical path** is taken when the user supplies a finite `Tolerance` or the matrix has any inexact (`Real`/`MPFR`) leaf: every entry is coerced to a `cplx_t` `{re, im}` struct via `entry_to_cplx` (handling Integer, BigInt, Real, MPFR, `Rational`, `Complex`, `I`, and a `N[expr]` fallback for `Pi`, `Sqrt[2]`, etc.), then `gauss_rank_cplx` runs partial-pivot Gaussian forward elimination over `double`-complex and counts accepted pivots — a column is skipped when its largest sub-pivot magnitude is `<= tol`. The default tolerance for inexact input is `max(rows,cols) · DBL_EPSILON · max(|entries|)` (the standard rank-by-SVD surrogate, with `max(|entries|)` substituted for the top singular value); for exact input the default is `0`.

The **exact path** (every leaf exact, no `Tolerance`, or the numeric coercion failed because of symbolic entries) calls `RowReduce[m, Method -> ...]` through the evaluator (`call_rowreduce`) and counts RREF rows that are not all structurally zero (`count_nonzero_rows`, using `is_zero_poly`). An optional `Method` option (`DivisionFreeRowReduction`, `OneStepRowReduction`, `CofactorExpansion`, `Automatic`) is forwarded to `RowReduce`.

**Data structures / limits.** Numerical path: a flat `cplx_t` array of size `rows*cols` with hand-rolled complex add/sub/mul/div. Exact path: defers entirely to the `RowReduce` dispatcher. Bad options emit `MatrixRank::opt`; non-rectangular input emits `MatrixRank::matrix`.

- `Protected`.
- Returns a non-negative `Integer` equal to the number of linearly
  independent rows of `m` (equivalently, of linearly independent
  columns).
- Works on numerical (Integer / Rational / Real / MPFR / Complex),
  big-integer, and symbolic matrices, square or rectangular.
- **Two execution paths**:
  - *Exact path* (every leaf is exact, no `Tolerance`): an all-integer/rational
    matrix gets its rank directly from FLINT (`fmpq_mat_rref`) in polynomial
    time when built with FLINT (rank is basis-independent, so the value matches
    the classical count); otherwise routes
    through `RowReduce[m, Method -> "<name>"]` and counts the
    non-zero rows of the RREF, using `is_zero_poly` for structural
    zero. Honors the same `Method` grammar as NullSpace / RowReduce /
    LinearSolve / Inverse (`Automatic` / `"DivisionFreeRowReduction"`
    / `"OneStepRowReduction"` / `"CofactorExpansion"`). The FLINT kernel is
    also exposed directly as `` FLINT`MatrixRank `` (see *Structural
    Manipulation*).
  - *Numerical path* (any inexact leaf, or any explicit `Tolerance`):
    runs partial-pivot Gaussian forward-elimination over a portable
    `double`-complex kernel with tolerance-aware pivot selection: a
    column is skipped when its largest sub-pivot `|entry|` is `<= t`.
    `Method` does not affect this path.
- **Tolerance** accepted forms:
  - `Tolerance -> Automatic` (default) — `max(rows, cols) *
    MachineEpsilon * Max[|entries|]` for inexact matrices; `0`
    otherwise (so integer / rational input agrees with the exact
    path).
  - `Tolerance -> 0` — no tolerance; even arbitrarily small entries
    count.
  - `Tolerance -> <non-negative number>` — Integer, Real, Rational,
    or a `Power[10, k]`-style expression resolved via `N[...]`.
- Issues `MatrixRank::matrix` and returns unevaluated if the
  argument is not a non-empty rank-2 tensor.
- Issues `MatrixRank::opt` and returns unevaluated for an unknown
  `Method` value, an unknown option key, or a negative / symbolic
  `Tolerance`.

**Attributes:** `Protected`.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — rank and the row echelon form.
- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — rank and linear independence.
- Source: [`src/linalg/matrank.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/matrank.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_diagonal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
- Tests: [`tests/test_matrank.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matrank.c)

## Notes & additional examples

### Notes

`MatrixRank` returns the number of linearly independent rows, which equals the number of independent columns. The first example is rank 1 because its rows are proportional, and the `3x3` numeric example is rank 2 because its three rows satisfy a linear relation. The default exact path routes through `RowReduce` and counts the non-zero rows; rectangular matrices are accepted. Working over the exact field means it also handles purely symbolic matrices — `{{a, b}, {2 a, 2 b}}` is recognised as rank 1 for *every* value of `a` and `b` because the second row is a symbolic multiple of the first. The `4x4` Hilbert matrix is famously ill-conditioned, yet exact rational arithmetic confirms it has full rank 4 with no round-off ambiguity. For inexact (Real / MPFR) input, or with an explicit `Tolerance`, a numerical partial-pivot elimination is used and entries with `|entry| <= t` are treated as zero.
