# FixedPoint

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FixedPoint[f, expr]
    starts with expr and applies f repeatedly until the result no longer
    changes, returning the final value.
FixedPoint[f, expr, n]
    stops after at most n applications of f, returning the last value
    obtained even if a fixed point has not been reached.
FixedPoint[f, expr, SameTest -> s]
FixedPoint[f, expr, n, SameTest -> s]
    uses the binary predicate s instead of SameQ to test successive pairs.

FixedPoint[f, expr] gives the last element of FixedPointList[f, expr].
Throw can be used inside f to exit early.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
