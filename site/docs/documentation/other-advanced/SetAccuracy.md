# SetAccuracy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SetAccuracy[x, n]
    Returns an expression equivalent to x with numeric values
    re-rounded or promoted to n decimal digits of accuracy.
    Requires a USE_MPFR build for high-accuracy outputs.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
