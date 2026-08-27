# HilbertMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HilbertMatrix[n] gives the n x n Hilbert matrix with entries 1/(i + j - 1).`**

**`HilbertMatrix[{m, n}] gives the m x n Hilbert matrix.`**

<details>
<summary>Notes</summary>

Entries are exact Rationals unless the WorkingPrecision option requests MachinePrecision or a digit count.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HilbertMatrix[3]
Out[1]= {{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}

In[2]:= HilbertMatrix[{3, 5}]
Out[2]= {{1, 1/2, 1/3, 1/4, 1/5}, {1/2, 1/3, 1/4, 1/5, 1/6}, {1/3, 1/4, 1/5, 1/6, 1/7}}

In[3]:= Det[HilbertMatrix[3]]
Out[3]= 1/2160

In[4]:= Inverse[HilbertMatrix[3]]
Out[4]= {{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}
```

### Options (1)

```mathematica
In[5]:= HilbertMatrix[3, WorkingPrecision -> MachinePrecision]
Out[5]= {{1.0, 0.5, 0.333333}, {0.5, 0.333333, 0.25}, {0.333333, 0.25, 0.2}}
```

### Applications (5)

```mathematica
In[6]:= HilbertMatrix[3]
Out[6]= {{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}

In[7]:= HilbertMatrix[{2, 4}]
Out[7]= {{1, 1/2, 1/3, 1/4}, {1/2, 1/3, 1/4, 1/5}}

In[8]:= Inverse[HilbertMatrix[3]]
Out[8]= {{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}

In[9]:= Det[HilbertMatrix[5]]
Out[9]= 1/266716800000

In[10]:= Eigenvalues[HilbertMatrix[2]]
Out[10]= {1/24 (16 + 4 Sqrt[13]), 1/24 (16 - 4 Sqrt[13])}
```

## Options & behaviour

| Option | Default | Meaning |
|--------|---------|---------|
| `WorkingPrecision` | `Infinity` | precision at which to create entries |

**Diagnostics** (the call is returned unevaluated):

## Algorithm

HilbertMatrix — the m x n Hilbert matrix with entries 1/(i + j - 1).

```text
  HilbertMatrix[n]        n x n Hilbert matrix.
  HilbertMatrix[{m, n}]   m x n Hilbert matrix.
```

Entries are exact Rationals by default (WorkingPrecision -> Infinity). The WorkingPrecision option selects the entry representation:

```text
  WorkingPrecision -> Infinity          exact Rationals (default)
  WorkingPrecision -> MachinePrecision  machine-precision Reals
  WorkingPrecision -> d                  d-digit MPFR Reals (d above
                                         machine precision; otherwise
                                         machine Reals, matching the
                                         rest of Mathilda's numeric
                                         tower).
```

Diagnostics mirror Wolfram's surface text:

```text
  - zero arguments               -> HilbertMatrix::argx
  - bad dimension specification  -> HilbertMatrix::dims
  - non-option trailing argument -> HilbertMatrix::nonopt
```

## Implementation notes

**Algorithm.** `builtin_hilbertmatrix` constructs the `m×n` Hilbert matrix with entry `(i,j) = 1/(i+j-1)`. The dimension spec (`hm_parse_dims`) is a positive integer `n` (square) or a pair `{m, n}` of positive integers; bad specs emit `HilbertMatrix::dims`, zero arguments emit `HilbertMatrix::argx`. The only recognised option is `WorkingPrecision` (`hm_parse_working_precision`, last-valid-setting-wins): `Infinity` (default) yields exact `Rational` entries via `make_rational`; `MachinePrecision` (or a digit count at/below machine precision) yields machine-precision `Real`s; a larger digit count yields MPFR reals (`mpfr_div_ui`) when built with `USE_MPFR`, degrading to machine reals otherwise (`HilbertMatrix::wprec`). Any trailing non-option argument triggers `HilbertMatrix::nonopt`.

**Data structures.** A `List` of `List`s built row by row; each entry is created by `hm_entry` according to the selected `hm_prec_mode` (`EXACT`/`MACHINE`/`MPFR`). Complexity is `O(mn)` entry constructions.

- `Protected`.
- Entries are exact `Rational`s by default. The matrix is symmetric and
  notoriously ill-conditioned, making it a standard test case for numeric
  linear-algebra routines.
- The `WorkingPrecision` option chooses the entry representation:
  - `WorkingPrecision -> Infinity` (default): exact `Rational`s.
  - `WorkingPrecision -> MachinePrecision`: machine-precision `Real`s.
  - `WorkingPrecision -> d`: `d`-digit arbitrary-precision (MPFR) `Real`s.
    A digit count at or below machine precision (or a build without MPFR)
    degrades to machine `Real`s.

**Attributes:** `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/)

- Source: [`src/linalg/hilbertmat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/hilbertmat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_fit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fit.c)
- Tests: [`tests/test_hilbertmatrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hilbertmatrix.c)
- Tests: [`tests/test_ludecomposition_mpfr.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ludecomposition_mpfr.c)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)

## Notes & additional examples

### Notes

`HilbertMatrix[n]` has entries `1/(i + j - 1)` and is the canonical ill-conditioned test matrix. Entries are exact `Rational`s, so `Det`, `Inverse`, and `Eigenvalues` return exact answers; request `WorkingPrecision -> MachinePrecision` (or a digit count) for an inexact version.
