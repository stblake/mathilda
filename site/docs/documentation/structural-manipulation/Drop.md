# Drop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Drop[list, n]
    gives list with its first n elements dropped.
Drop[list, -n]
    drops the last n elements.
Drop[list, {m, n}]
    drops elements m through n.
Drop[list, {m, n, s}]
    drops elements m through n in steps of s.
Drop[list, {m}]
    drops the single element at position m.
Drop[list, spec1, spec2, ...]
    drops elements at successive levels.

Negative indices count from the end; UpTo[n], All, and None are also accepted.
Indices are 1-based; out-of-range requests leave the expression unevaluated.
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
