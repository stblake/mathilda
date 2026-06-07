# MatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MatrixQ[expr]
    gives True if expr is a list of lists that can represent a matrix, and gives False otherwise.
MatrixQ[expr, test]
    gives True only if test yields True when applied to each of the matrix elements in expr.

MatrixQ[expr] gives True only if expr is a list and each of its elements is a list of the same length,
containing no elements that are themselves lists.
MatrixQ[expr, NumberQ] tests whether expr is a numerical matrix.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
