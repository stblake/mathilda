# Position

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Position[expr, pattern] gives a list of the positions at which objects matching pattern appear in expr.
Position[expr, pattern, levelspec] finds only objects that appear on levels specified by levelspec.
Position[expr, pattern, levelspec, n] gives the positions of the first n objects found.
Position[pattern] represents an operator form of Position that can be applied to an expression.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Defaults to levels `{0, Infinity}` with `Heads -> True`.
- Yields lists of indices in lexicographic order.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
