# Most

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Most[expr] gives all but the last element of expr.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_most` returns a copy of the input with its last element dropped: it
copies args `0 .. n−2` into a new function node with the same head. Returns `NULL` (unevaluated)
for atoms or empty expressions.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
