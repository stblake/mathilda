# HilbertMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HilbertMatrix[n] gives the n x n Hilbert matrix with entries 1/(i + j - 1).
HilbertMatrix[{m, n}] gives the m x n Hilbert matrix.
Entries are exact Rationals unless the WorkingPrecision option requests MachinePrecision or a digit count.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HilbertMatrix[3]
Out[1]= {{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}

In[2]:= HilbertMatrix[{3, 5}]
Out[2]= {{1, 1/2, 1/3, 1/4, 1/5}, {1/2, 1/3, 1/4, 1/5, 1/6}, {1/3, 1/4, 1/5, 1/6, 1/7}}

In[3]:= HilbertMatrix[3, WorkingPrecision -> MachinePrecision]
Out[3]= {{1.0, 0.5, 0.333333}, {0.5, 0.333333, 0.25}, {0.333333, 0.25, 0.2}}

In[4]:= Det[HilbertMatrix[3]]
Out[4]= 1/2160

In[5]:= Inverse[HilbertMatrix[3]]
Out[5]= {{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}
```

## Implementation notes

**Algorithm.** `builtin_hilbertmatrix` constructs the `m×n` Hilbert matrix with entry `(i,j) = 1/(i+j-1)`. The dimension spec (`hm_parse_dims`) is a positive integer `n` (square) or a pair `{m, n}` of positive integers; bad specs emit `HilbertMatrix::dims`, zero arguments emit `HilbertMatrix::argx`. The only recognised option is `WorkingPrecision` (`hm_parse_working_precision`, last-valid-setting-wins): `Infinity` (default) yields exact `Rational` entries via `make_rational`; `MachinePrecision` (or a digit count at/below machine precision) yields machine-precision `Real`s; a larger digit count yields MPFR reals (`mpfr_div_ui`) when built with `USE_MPFR`, degrading to machine reals otherwise (`HilbertMatrix::wprec`). Any trailing non-option argument triggers `HilbertMatrix::nonopt`.

**Data structures.** A `List` of `List`s built row by row; each entry is created by `hm_entry` according to the selected `hm_prec_mode` (`EXACT`/`MACHINE`/`MPFR`). Complexity is `O(mn)` entry constructions.

- `Protected`.
- Entries are exact `Rational`s by default. The matrix is symmetric and

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/hilbertmat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/hilbertmat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
