# Extract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Extract[expr, pos]
    extracts the part of expr at the position specified by pos.
Extract[expr, {pos1, pos2, ...}]
    extracts a list of parts of expr.
Extract[expr, pos, h]
    extracts parts of expr, wrapping each of them with head h before evaluation.
Extract[pos]
    represents an operator form of Extract that can be applied to an expression.

The position pos has the same form as a Position result: a list of indices {i1, i2, ...} that descends i1 into expr, then i2, etc. A scalar index n is treated as the path {n}, so Extract[expr, n] is equivalent to Extract[expr, {n}]; in particular Extract[expr, 0] gives Head[expr].
Indices are 1-based and may be negative; index 0 selects the head. Extract is treated as atomic on Integer, Real, String, Symbol, Rational[n, d], and Complex[re, im].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Position specifications have the same form as those returned by `Position`.
- `Extract[expr, {i, j, ...}]` is equivalent to `Part[expr, i, j, ...]`.
- `pos` can be of the more general form `{part1, part2, ...}` where `parti` are `Part` specifications such as an integer `i`, `All` or `Span`.
- You can use `Extract[expr, ..., Hold]` to extract parts without evaluation.

**Attributes:** `NHoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
