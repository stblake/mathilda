# Hold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Hold[expr]
    maintains expr in an unevaluated form.
Hold has attribute HoldAll: its arguments are not evaluated.
Evaluate[expr] inside Hold overrides the hold and evaluates expr once.
Sequence expressions inside Hold are flattened; use HoldComplete to prevent this.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
