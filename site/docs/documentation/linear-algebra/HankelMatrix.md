# HankelMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HankelMatrix[n] gives the n x n Hankel matrix with first row and column the integers 1..n.
HankelMatrix[{c1, ..., cm}] gives the m x m Hankel matrix with first column the given list.
HankelMatrix[{c1, ..., cm}, {r1, ..., rn}] gives the m x n Hankel matrix with first column the first list and last row the second.
A Hankel matrix is constant along its antidiagonals; entries are copied verbatim.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HankelMatrix[4]
Out[1]= {{1, 2, 3, 4}, {2, 3, 4, 0}, {3, 4, 0, 0}, {4, 0, 0, 0}}

In[2]:= HankelMatrix[{a, b, c, d}]
Out[2]= {{a, b, c, d}, {b, c, d, 0}, {c, d, 0, 0}, {d, 0, 0, 0}}

In[3]:= HankelMatrix[{x, y, z}, {z, a, b, c, d}]
Out[3]= {{x, y, z, a, b}, {y, z, a, b, c}, {z, a, b, c, d}}

In[4]:= HankelMatrix[{1, 1 + 2 I, 3 + 4 I}]
Out[4]= {{1, 1 + 2*I, 3 + 4*I}, {1 + 2*I, 3 + 4*I, 0}, {3 + 4*I, 0, 0}}

In[5]:= N[HankelMatrix[3]]
Out[5]= {{1.0, 2.0, 3.0}, {2.0, 3.0, 0.0}, {3.0, 0.0, 0.0}}
```

## Implementation notes

**Algorithm.** `builtin_hankelmatrix` builds a matrix constant along antidiagonals, where entry `(i,j)` depends only on `s = i+j-1`. Three forms are handled: `HankelMatrix[n]` (square, antidiagonal index `s` for `s ≤ n` else `0`, the integer form); `HankelMatrix[{c1,…,cm}]` (square, first column `c`, zeros below the antidiagonal); and `HankelMatrix[{c…}, {r…}]` (`m×n`, first column `c` and last row `r`, with `(i,j) = c_s` when `s ≤ m` and `r_{s-m}` otherwise). The shared corner `c_m` must equal `r_1`; if not, it warns via `HankelMatrix::crs` and uses the column element. Zero arguments emit `HankelMatrix::argb`; any other shape (non-list, empty list) is left unevaluated.

**Data structures.** A `List` of `List`s built by `hk_build`; source entries are deep-copied (`expr_copy`) so symbolic/complex/exact/inexact entries pass through verbatim — arbitrary precision comes from the entries themselves. Complexity `O(mn)`.

- `Protected`.
- Entries are copied verbatim, so symbolic, exact, complex, machine and

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/hankelmat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/hankelmat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= HankelMatrix[3]
Out[1]= {{1, 2, 3}, {2, 3, 0}, {3, 0, 0}}
```

Give an explicit first column; the matrix is constant along every
antidiagonal and zero-padded below the secondary diagonal:

```mathematica
In[1]:= HankelMatrix[{a, b, c}]
Out[1]= {{a, b, c}, {b, c, 0}, {c, 0, 0}}
```

A first column plus a last row builds a rectangular catalecticant — here the
shared corner is taken from the column:

```mathematica
In[1]:= HankelMatrix[{1, 2, 3, 4}, {4, 5, 6}]
Out[1]= {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {4, 5, 6}}
```

Hankel determinants encode sequence properties. The catalecticant of the
Fibonacci numbers is a perfect power of 2:

```mathematica
In[1]:= Det[HankelMatrix[{1, 1, 2, 3, 5, 8}]]
Out[1]= -262144
```

### Notes

A Hankel matrix is constant along its antidiagonals. With a single integer `n`, the first row and column are the integers `1..n` and the lower-right triangle is zero-filled. Supplying a first column (and optionally a last row) lets you build the catalecticant of any sequence; `Det[HankelMatrix[...]]` then gives that sequence's Hankel determinant.
