# Take

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Take[list, n]
    gives the first n elements of list.
Take[list, -n]
    gives the last n elements.
Take[list, {m, n}]
    gives elements m through n.
Take[list, {m, n, s}]
    gives elements m through n in steps of s.
Take[list, {m}]
    gives the single element at position m (wrapped in the head of list).
Take[list, spec1, spec2, ...]
    takes elements at successive levels, e.g. a sub-block of a matrix.

Negative indices count from the end; UpTo[n], All, and None are also accepted
as specifications. Indices are 1-based; out-of-range requests leave the
expression unevaluated. Take operates on any expression, not just List.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `NHoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
