# DiagonalMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DiagonalMatrix[list] gives a matrix with the elements of list on the leading diagonal, and zero elsewhere.
DiagonalMatrix[list, k] gives a matrix with the elements of list on the k-th diagonal.
DiagonalMatrix[list, k, n] pads with zeros to create an n x n matrix.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DiagonalMatrix[{a, b, c}]
Out[1]= {{a, 0, 0}, {0, b, 0}, {0, 0, c}}

In[2]:= DiagonalMatrix[{a, b}, 1]
Out[2]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}

In[3]:= DiagonalMatrix[{1, 2, 3}, 0, {3, 5}]
Out[3]= {{1, 0, 0, 0, 0}, {0, 2, 0, 0, 0}, {0, 0, 3, 0, 0}}
```

## Implementation notes

`builtin_diagonalmatrix` builds a matrix placing the entries of the given `List` on the `k`-th diagonal. With one argument the entries go on the main diagonal (`k = 0`); a second integer argument `k` selects a super-/sub-diagonal (`j - i == k`), sizing the matrix to `(s + |k|) × (s + |k|)` where `s` is the list length; an optional third argument fixes the output dimensions as `n` or `{m, n}`. Off-diagonal cells are `Integer` `0`; diagonal cells are deep-copied verbatim from the input list, so symbolic, exact, and inexact entries flow through unchanged. Malformed `k`/dimension specs return the call unevaluated. The result is a `List` of `List`s.

- `Protected`.
- For `k > 0`, places elements `k` positions above the leading diagonal.
- For `k < 0`, places elements `k` positions below the leading diagonal.
- By default, size is optimally bounded to fit the full array cleanly. Extraneous elements are dropped if manual constraints fall short of required lengths.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/construct.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/construct.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
