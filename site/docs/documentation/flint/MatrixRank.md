# FLINT`MatrixRank

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`MatrixRank[m] gives the rank of the matrix m when every entry is an integer or rational, computed exactly via FLINT (fmpq_mat_rref). Returns unevaluated for a matrix with any non-rational entry.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`MatrixRank[{{1, 2}, {2, 4}}]
Out[1]= 1
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/flint_mat_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/flint_mat_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
