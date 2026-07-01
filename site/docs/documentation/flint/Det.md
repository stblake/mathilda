# FLINT`Det

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`Det[m] gives the determinant of the square matrix m when every entry is an integer or rational, computed exactly and directly via FLINT (fmpq_mat_det). Returns unevaluated for a matrix with any non-rational entry.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`Det[{{1, 2}, {3, 4}}]
Out[1]= -2
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/flint_mat_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/flint_mat_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
