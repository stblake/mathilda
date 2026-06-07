# DeleteCases

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DeleteCases[expr, pattern] removes all elements of expr that match pattern.
DeleteCases[expr, pattern, levelspec] removes all parts of expr on levels specified by levelspec that match pattern.
DeleteCases[expr, pattern, levelspec, n] removes the first n parts of expr that match pattern.
DeleteCases[pattern] represents an operator form of DeleteCases that can be applied to an expression.
The default levelspec is {1}. With Heads -> True, the heads of expressions are also tested; deleting a head is equivalent to applying FlattenAt at that location.
DeleteCases traverses expr in depth-first post-order (leaves before roots).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Uses standard level specifications, defaulting to level `{1}`.
- Option `Heads -> True` tests heads as well; deleting a head is equivalent to applying `FlattenAt` at that point, splicing the remaining arguments into the parent.
- Traverses `expr` in depth-first post-order (leaves before roots) so that the `n` budget is spent on deeper matches before shallower ones.
- The match test is applied to the original subexpression (not the version with children already deleted), matching Mathematica semantics.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
