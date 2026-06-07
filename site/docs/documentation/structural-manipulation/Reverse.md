# Reverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Reverse[expr] reverses the order of elements in expr.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_reverse` recurses through the expression with `reverse_rec`, reversing
the argument order at the levels selected by an optional level spec (`should_reverse_at_level`
matches an integer level, or any level listed in a `List` spec; default is level 1). At each
visited node it builds a new function with the same head, drawing children either forward or
reversed depending on whether the current level is selected, and recursing into each child.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
