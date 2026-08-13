# ToeplitzMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToeplitzMatrix[n] gives the n x n Toeplitz matrix with first row and column the integers 1..n.`**

**`ToeplitzMatrix[{c1, ..., cn}] gives the n x n symmetric Toeplitz matrix with first column the given list.`**

**`ToeplitzMatrix[{c1, ..., cm}, {r1, ..., rn}] gives the m x n Toeplitz matrix with first column the first list and first row the second.`**

<details>
<summary>Notes</summary>

A Toeplitz matrix is constant along its diagonals; entries are copied verbatim.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (4)

```mathematica
In[6]:= ToeplitzMatrix[3]
Out[6]= {{1, 2, 3}, {2, 1, 2}, {3, 2, 1}}

In[7]:= ToeplitzMatrix[{a, b, c}, {a, x, y}]
Out[7]= {{a, x, y}, {b, a, x}, {c, b, a}}

In[8]:= Det[ToeplitzMatrix[{2, -1, 0, 0, 0}, {2, -1, 0, 0, 0}]]
Out[8]= 6

In[9]:= Eigenvalues[ToeplitzMatrix[{2, -1, 0}, {2, -1, 0}]]
Out[9]= {1/2 (4 + 2 Sqrt[2]), 2, 1/2 (4 - 2 Sqrt[2])}
```

## Options & behaviour

### Diagnostics

A first argument that is neither a positive integer nor a list (and any
over-arity call) is returned unevaluated.

## Algorithm

ToeplitzMatrix — a matrix that is constant along its diagonals.

```text
  ToeplitzMatrix[n]            n x n Toeplitz matrix whose first row and
                               first column are the successive integers
                               1..n (symmetric: entry (i, j) is |i - j| + 1).
  ToeplitzMatrix[{c1,...,cn}]  n x n symmetric Toeplitz matrix whose first
                               column (and first row) is the given list.
  ToeplitzMatrix[{c1,...,cm},  m x n Toeplitz matrix with the first list
                 {r1,...,rn}]  down the first column and the second list
                               across the first row.
```

The entry (i, j) is c_{i-j+1} when i >= j, and r_{j-i+1} otherwise (1-based

```text
indices).  The shared corner c_1 and r_1 should be equal; if they differ,
```

the column element is used and a ToeplitzMatrix::crs warning is emitted (the

```text
formula always reads c_1 on the diagonal, never r_1).  Entries are copied
```

verbatim, so symbolic, complex, exact and inexact entries all flow through unchanged; arbitrary precision comes from the entries themselves (e.g.

```text
`1`20`) or from wrapping in N.
```

Diagnostics mirror Wolfram's surface text:

```text
  - zero arguments  ->  ToeplitzMatrix::argb
  - mismatched corner element  ->  ToeplitzMatrix::crs (warning; still builds)
```

## Implementation notes

**Algorithm.** `builtin_toeplitzmatrix` constructs a matrix constant along diagonals. Three forms dispatch on argument shape: `ToeplitzMatrix[n]` builds the symmetric `n × n` integer matrix with entry `(i,j) = |i−j| + 1`; `ToeplitzMatrix[{c}]` builds the symmetric matrix whose first column and row are the list `c`; `ToeplitzMatrix[{c}, {r}]` builds the `m × n` matrix with `c` down the first column and `r` across the first row. The builder `tz_build` sets entry `(i,j)` to `cvals[i−j]` when `i >= j` and `rvals[j−i]` otherwise, deep-copying source entries verbatim so symbolic/complex/exact/inexact entries flow through unchanged.

**Limits.** The shared corner reads `c_1` on the diagonal; if `c_1 != r_1` it warns `ToeplitzMatrix::crs` and uses the column element. Zero arguments emit `ToeplitzMatrix::argb`; empty lists or any other shape leave the call unevaluated.

- `Protected`.
- Entries are copied verbatim, so symbolic, exact, complex, machine and
  arbitrary-precision entries all flow through unchanged; precision comes from
  the entries themselves (e.g. ``1`20``) or from wrapping the result in `N`.
  The single-list form is plain symmetric (no conjugation).
- For `m = n` the matrix is symmetric, and has real eigenvalues when the
  entries are real.
- The shared corner element `c_1` must equal `r_1`. If they differ, the
  column element `c_1` is used (the formula never reads `r_1`, which sits on
  the diagonal as `c_1`) and a `ToeplitzMatrix::crs` warning is emitted; the
  matrix is still produced.
- An all-integer or all-machine-real argument writes the rank-2 result buffer
  directly (`ndbuild_open`); both the single-vector form
  (`Compile[{{v, _Real, 1}}, ToeplitzMatrix[v]]`, rank 1 → rank 2) and the
  two-vector form (`ToeplitzMatrix[c, r]`, rank 1 × rank 1 → rank 2)
  **compile**. See [`packed-arrays.md`](../packed-arrays/index.md).

**Attributes:** `Protected`.

## References

**See also:** [N](../../arithmetic/N/)

- Source: [`src/linalg/toeplitzmat.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/toeplitzmat.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_toeplitzmatrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_toeplitzmatrix.c)
