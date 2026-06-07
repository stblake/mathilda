# MatrixRank

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MatrixRank[m]
    gives the rank of the matrix m -- the number of linearly
    independent rows (equivalently, of linearly independent
    columns).
MatrixRank[m, Method -> "<name>"]
    runs a specific elimination algorithm for the exact path.
    Accepted method names match NullSpace / RowReduce /
    LinearSolve / Inverse:
      "Automatic"                 -- alias for "DivisionFreeRowReduction" (default)
      "DivisionFreeRowReduction"  -- Bareiss-like fraction-free Gauss-Jordan
      "OneStepRowReduction"       -- classical Gauss-Jordan with division per pivot
      "CofactorExpansion"         -- identity-if-invertible (falls back to
                                     DivisionFreeRowReduction on singular /
                                     rectangular m)
MatrixRank[m, Tolerance -> t]
    treats |entry| <= t as zero during pivot selection.  With
    Tolerance -> 0 even arbitrarily small entries count; the
    default, Tolerance -> Automatic, applies
    max(rows, cols) * MachineEpsilon * Max[|entries|] for
    machine-precision (Real / MPFR) matrices and 0 otherwise.

MatrixRank works on both numerical and symbolic matrices and
on square or rectangular matrices.  The default exact path
routes through RowReduce and counts the non-zero rows; the
numerical path (triggered by inexact leaves or an explicit
Tolerance) runs partial-pivot Gaussian elimination over
double-precision complex.

An unknown Method value or Tolerance form emits
MatrixRank::opt and leaves the call unevaluated.  A non-rank-2
or empty matrix emits MatrixRank::matrix and the call is left
unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

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

In[6]:= MatrixRank[N[m]]
Out[6]= 2

In[7]:= MatrixRank[N[m], Tolerance -> 0]
Out[7]= 3
```

## Implementation notes

- `Protected`.
- Returns a non-negative `Integer` equal to the number of linearly

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — rank and the row echelon form.
- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — rank and linear independence.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

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

### Notes

`MatrixRank` returns the number of linearly independent rows, which equals the number of independent columns. The first example is rank 1 because its rows are proportional, and the `3x3` example is rank 2 because its three rows satisfy a linear relation. The default exact path routes through `RowReduce` and counts the non-zero rows; rectangular matrices are accepted. For inexact (Real / MPFR) input, or with an explicit `Tolerance`, a numerical partial-pivot elimination is used and entries with `|entry| <= t` are treated as zero.
