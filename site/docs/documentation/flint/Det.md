# FLINT`Det

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

FLINT\`Det\[m\] gives the determinant of the square matrix m when every entry is an integer or rational, computed exactly and directly via FLINT (fmpq\_mat\_det). Returns unevaluated for a matrix with any non-rational entry.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FLINT`Det[{{1, 2}, {3, 4}}]
Out[1]= -2
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/linalg/flint_mat_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/flint_mat_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
