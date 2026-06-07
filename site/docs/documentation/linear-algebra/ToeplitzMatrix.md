# ToeplitzMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToeplitzMatrix[n] gives the n x n Toeplitz matrix with first row and column the integers 1..n.
ToeplitzMatrix[{c1, ..., cn}] gives the n x n symmetric Toeplitz matrix with first column the given list.
ToeplitzMatrix[{c1, ..., cm}, {r1, ..., rn}] gives the m x n Toeplitz matrix with first column the first list and first row the second.
A Toeplitz matrix is constant along its diagonals; entries are copied verbatim.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ToeplitzMatrix[4]
Out[1]= {{1, 2, 3, 4}, {2, 1, 2, 3}, {3, 2, 1, 2}, {4, 3, 2, 1}}

In[2]:= ToeplitzMatrix[{a, b, c, d}]
Out[2]= {{a, b, c, d}, {b, a, b, c}, {c, b, a, b}, {d, c, b, a}}

In[3]:= ToeplitzMatrix[{1, 2, 3, 4, 5}, {1, 6, 7}]
Out[3]= {{1, 6, 7}, {2, 1, 6}, {3, 2, 1}, {4, 3, 2}, {5, 4, 3}}

In[4]:= ToeplitzMatrix[{1, 2, 3}, {1, 4, 5, 6, 7}]
Out[4]= {{1, 4, 5, 6, 7}, {2, 1, 4, 5, 6}, {3, 2, 1, 4, 5}}

In[5]:= N[ToeplitzMatrix[3]]
Out[5]= {{1.0, 2.0, 3.0}, {2.0, 1.0, 2.0}, {3.0, 2.0, 1.0}}
```

## Implementation notes

**Algorithm.** `builtin_toeplitzmatrix` constructs a matrix constant along diagonals. Three forms dispatch on argument shape: `ToeplitzMatrix[n]` builds the symmetric `n × n` integer matrix with entry `(i,j) = |i−j| + 1`; `ToeplitzMatrix[{c}]` builds the symmetric matrix whose first column and row are the list `c`; `ToeplitzMatrix[{c}, {r}]` builds the `m × n` matrix with `c` down the first column and `r` across the first row. The builder `tz_build` sets entry `(i,j)` to `cvals[i−j]` when `i >= j` and `rvals[j−i]` otherwise, deep-copying source entries verbatim so symbolic/complex/exact/inexact entries flow through unchanged.

**Limits.** The shared corner reads `c_1` on the diagonal; if `c_1 != r_1` it warns `ToeplitzMatrix::crs` and uses the column element. Zero arguments emit `ToeplitzMatrix::argb`; empty lists or any other shape leave the call unevaluated.

- `Protected`.
- Entries are copied verbatim, so symbolic, exact, complex, machine and

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/toeplitzmat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/toeplitzmat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
