# VectorQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
VectorQ[expr]
    gives True if expr is a list, none of whose elements are themselves lists, and gives False otherwise.
VectorQ[expr, test]
    gives True only if test yields True when applied to each of the elements in expr.

VectorQ[expr, NumberQ] tests whether expr is a vector of numbers.
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
