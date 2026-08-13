# DiagonalMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DiagonalMatrix[list] gives a matrix with the elements of list on the leading diagonal, and zero elsewhere.`**

**`DiagonalMatrix[list, k] gives a matrix with the elements of list on the k-th diagonal.`**

**`DiagonalMatrix[list, k, n] pads with zeros to create an n x n matrix.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= DiagonalMatrix[{a, b, c}]
Out[1]= {{a, 0, 0}, {0, b, 0}, {0, 0, c}}

In[2]:= DiagonalMatrix[{a, b}, 1]
Out[2]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}

In[3]:= DiagonalMatrix[{1, 2, 3}, 0, {3, 5}]
Out[3]= {{1, 0, 0, 0, 0}, {0, 2, 0, 0, 0}, {0, 0, 3, 0, 0}}

In[4]:= DiagonalMatrix[{1., 2., 3.}]
Out[4]= {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}}
```

### Applications (4)

```mathematica
In[5]:= DiagonalMatrix[{1, 2, 3}]
Out[5]= {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}

In[6]:= DiagonalMatrix[{a, b}, 1]
Out[6]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}

In[7]:= DiagonalMatrix[{x, y, z}, -1, 4]
Out[7]= {{0, 0, 0, 0}, {x, 0, 0, 0}, {0, y, 0, 0}, {0, 0, z, 0}}

In[8]:= DiagonalMatrix[{1, 1, 1, 1}, 2]
Out[8]= {{0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}
```

## Implementation notes

`builtin_diagonalmatrix` builds a matrix placing the entries of the given `List` on the `k`-th diagonal. With one argument the entries go on the main diagonal (`k = 0`); a second integer argument `k` selects a super-/sub-diagonal (`j - i == k`), sizing the matrix to `(s + |k|) × (s + |k|)` where `s` is the list length; an optional third argument fixes the output dimensions as `n` or `{m, n}`. Off-diagonal cells are `Integer` `0`; diagonal cells are deep-copied verbatim from the input list, so symbolic, exact, and inexact entries flow through unchanged. Malformed `k`/dimension specs return the call unevaluated. The result is a `List` of `List`s.

- `Protected`.
- For `k > 0`, places elements `k` positions above the leading diagonal.
- For `k < 0`, places elements `k` positions below the leading diagonal.
- By default, size is optimally bounded to fit the full array cleanly. Extraneous elements are dropped if manual constraints fall short of required lengths.
- **The zeros take the diagonal's exactness.** A machine `Real` anywhere on the
  diagonal makes the whole matrix machine-real — the invented zeros and the
  exact entries alike — so `DiagonalMatrix[{1, 2, 3.}]` is
  `{{1., 0., 0.}, {0., 2., 0.}, {0., 0., 3.}}` and `DiagonalMatrix[{1/2, 1.}]`
  is `{{0.5, 0.}, {0., 1.}}`. An exact diagonal keeps exact zeros. Only
  *machine* `Real` is contagious: an MPFR entry keeps them exact
  (``DiagonalMatrix[{1.`30, 2}]``), and a symbolic entry stays symbolic while
  the zeros around it still turn `Real` (`DiagonalMatrix[{a, 1.}]` is
  `{{a, 0.}, {0., 1.}}`).
- A uniformly exact-integer or uniformly machine-real result is a
  [packed list](../packed-arrays/index.md) written directly, without building the
  elements — which is what the exactness rule above buys: a matrix of two heads
  fits no buffer, and the exact zeros used to cost a `Real` diagonal 320×
  against NumPy's `np.diag`. The single-vector form **compiles**
  (`Compile[{{v, _Real, 1}}, DiagonalMatrix[v]]`, rank 1 → rank 2). See
  [`packed-arrays.md`](../packed-arrays/index.md).

**Attributes:** `Protected`.

## References

- Source: [`src/linalg/construct.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/construct.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_diagonal_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal_matrix_q.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_lapack_builtin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lapack_builtin.c)

## Notes & additional examples

### Notes

`DiagonalMatrix[list]` places `list` on the main diagonal of an otherwise zero
square matrix. The two-argument form `DiagonalMatrix[list, k]` shifts the band to
the `k`-th diagonal — positive `k` lies above the main diagonal (a superdiagonal),
negative `k` below it (a subdiagonal) — and the matrix grows just large enough to
hold that band, so `DiagonalMatrix[{a, b}, 1]` is `3 x 3`. The three-argument form
`DiagonalMatrix[list, k, n]` pads with zeros to force an explicit `n x n` size, as
in the `2 x` Jordan-style superdiagonal block of ones above.
