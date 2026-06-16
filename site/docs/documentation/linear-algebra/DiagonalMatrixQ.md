# DiagonalMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DiagonalMatrixQ[m]
    gives True if m is diagonal, and False otherwise.
DiagonalMatrixQ[m, k]
    gives True if m has nonzero elements only on the k-th diagonal,
    and False otherwise.  Positive k refers to superdiagonals above
    the main diagonal; negative k refers to subdiagonals below it.
    Works for rectangular as well as square matrices.

Option:
    Tolerance -> Automatic   numeric tolerance for approximate matrices.

With Tolerance -> t, off-diagonal entries e are taken to be zero
when Abs[e] <= t evaluates to True.  Without a tolerance the test
is structural: only literal numeric zeros (Integer 0, Real 0.0,
BigInt 0) count as zero.  Returns False on non-matrix, ragged,
empty (i.e. {}), or higher-rank tensor inputs; an n-by-0 matrix
(e.g. {{}, {}}) is vacuously diagonal.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_diagonal_matrix_q` tests whether a matrix has nonzero entries only on the `k`-th diagonal (default `k = 0`). It accepts an optional integer `k` at position 2 and a `Tolerance` option; an empty or malformed argument list yields a `DiagonalMatrixQ::argt`/`::nonopt` diagnostic, and missing args return `False`. After validating that the matrix is a rectangular `List` of `List`s, it returns `True`/`False` according to whether every off-`k`-diagonal entry is structurally (or within `Tolerance`) zero.

- `Protected`.
- Works for rectangular matrices, not only square -- only the entry-zero

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= DiagonalMatrixQ[{{1, 0}, {0, 2}}]
Out[1]= True

In[2]:= DiagonalMatrixQ[{{1, 2}, {0, 3}}]
Out[2]= False
```

```mathematica
In[1]:= DiagonalMatrixQ[DiagonalMatrix[{a, b, c}]]
Out[1]= True
```

```mathematica
In[1]:= DiagonalMatrixQ[{{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}, 1]
Out[1]= True
```

```mathematica
In[1]:= DiagonalMatrixQ[{{0.0, 1.0*10^-15}, {0, 0.0}}, Tolerance -> 10^-10]
Out[1]= True
```

### Notes

Without a `Tolerance` option the test is structural: only literal numeric zeros off the main diagonal count as zero. Returns `False` on non-matrix, ragged, or higher-rank inputs.
