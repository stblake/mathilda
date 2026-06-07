# UpperTriangularMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
UpperTriangularMatrixQ[m]
    gives True if m is upper triangular, and False otherwise.
UpperTriangularMatrixQ[m, k]
    gives True if m is upper triangular starting from the k-th
    diagonal (every entry m[i,j] with j - i < k is zero), and
    False otherwise.  Positive k refers to superdiagonals above
    the main diagonal; negative k refers to subdiagonals below it.
    Works for rectangular as well as square matrices.

Option:
    Tolerance -> Automatic   numeric tolerance for approximate matrices.

With Tolerance -> t, sub-diagonal entries e are taken to be zero
when Abs[e] <= t evaluates to True.  Without a tolerance the test
is structural: only literal numeric zeros (Integer 0, Real 0.0,
BigInt 0) count as zero.  Returns False on non-matrix, ragged,
empty (i.e. {}), or higher-rank tensor inputs; an n-by-0 matrix
(e.g. {{}, {}}) is vacuously upper triangular.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Works for rectangular matrices, not only square -- only the entry-zero

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
