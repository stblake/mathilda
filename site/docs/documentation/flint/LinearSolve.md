# FLINT`LinearSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`LinearSolve[m, b] solves the square system m.x == b exactly via FLINT (fmpq_mat_solve) when m is a nonsingular rational matrix and b a rational vector or matrix. Returns unevaluated for a non-square, singular, or non-rational system.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
Out[1]= {-4, 9/2}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/flint_mat_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/flint_mat_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
