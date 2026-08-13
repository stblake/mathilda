# FLINT`RowReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

FLINT\`RowReduce\[m\] gives the reduced row echelon form of the matrix m when every entry is an integer or rational, computed exactly via FLINT (fmpq\_mat\_rref). Returns unevaluated for a matrix with any non-rational entry.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FLINT`RowReduce[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{1, 0, -1}, {0, 1, 2}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/linalg/flint_mat_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/flint_mat_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
