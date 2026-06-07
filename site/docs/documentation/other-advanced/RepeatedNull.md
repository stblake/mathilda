# RepeatedNull

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
p... or RepeatedNull[p] is a pattern object that represents a sequence of zero or more expressions, each matching p.
RepeatedNull[p, max] represents from 0 to max expressions matching p.
RepeatedNull[p, {min, max}] represents between min and max expressions matching p.
RepeatedNull[p, {n}] represents exactly n expressions matching p.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
