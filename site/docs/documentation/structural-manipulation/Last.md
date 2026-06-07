# Last

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Last[expr] gives the last element of expr.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_last` (in `src/part.c`) takes a single argument and returns a deep copy of its final element (`args[arg_count - 1]`). It returns `NULL` (unevaluated) when the argument is atomic or empty.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
