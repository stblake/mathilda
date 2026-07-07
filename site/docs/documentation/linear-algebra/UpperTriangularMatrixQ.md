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

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= UpperTriangularMatrixQ[{{a, b, c}, {0, e, f}, {0, 0, g}}]
Out[1]= True

In[2]:= UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[2]= False

In[3]:= UpperTriangularMatrixQ[{{0, 1, 2}, {0, 0, 3}, {0, 0, 0}}, 1]
Out[3]= True

In[4]:= UpperTriangularMatrixQ[{{1, 2, 3}, {4, 5, 6}, {0, 7, 9}}, -1]
Out[4]= True

In[5]:= UpperTriangularMatrixQ[{{1, 2, 3}, {0, 4, 5}}]
Out[5]= True

In[6]:= UpperTriangularMatrixQ[{{1, 2}, {0, 4}, {0, 0}}]
Out[6]= True

In[7]:= UpperTriangularMatrixQ[IdentityMatrix[5]]
Out[7]= True
```

## Implementation notes

`builtin_upper_triangular_matrix_q` validates that the argument is a non-empty rectangular `List` of equal-length `List`s (no deeper nesting), then returns `True` iff every entry strictly below the `k`-th diagonal (column−row `< k`, default `k = 0`) is zero. The optional second Integer argument selects the diagonal `k`; a `Tolerance -> t` option relaxes the zero test (otherwise a structural zero check is used). Bad arguments/options emit `UpperTriangularMatrixQ::nonopt` / `::argt`; shape rejections return `False`.

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
In[1]:= UpperTriangularMatrixQ[{{1, 2}, {0, 3}}]
Out[1]= True

In[2]:= UpperTriangularMatrixQ[{{1, 0}, {2, 3}}]
Out[2]= False
```

The test is purely structural, so it works on a symbolically generated triangular matrix just as well as on numeric data:

```mathematica
In[1]:= UpperTriangularMatrixQ[Table[If[j >= i, a[i, j], 0], {i, 4}, {j, 4}]]
Out[1]= True
```

The two-argument form shifts the reference diagonal. A strictly upper-triangular matrix (zero main diagonal) passes the `k = 1` test:

```mathematica
In[1]:= UpperTriangularMatrixQ[{{0, 1, 2}, {0, 0, 3}, {0, 0, 0}}, 1]
Out[1]= True
```

For approximate matrices, a `Tolerance` lets tiny round-off entries below the diagonal count as zero:

```mathematica
In[1]:= UpperTriangularMatrixQ[{{1.0, 2}, {1*10^-15, 3}}, Tolerance -> 10^-10]
Out[1]= True
```

### Notes

A matrix is upper triangular when every entry below the main diagonal is zero. The two-argument form `UpperTriangularMatrixQ[m, k]` tests against the k-th diagonal (positive `k` for superdiagonals, negative for subdiagonals), and a `Tolerance` option relaxes the zero test for approximate matrices. Rectangular matrices are supported.
