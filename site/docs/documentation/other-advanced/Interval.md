# Interval

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Interval[{min, max}] represents the range of real values between min and max, inclusive. Interval[{a1,b1}, {a2,b2}, ...] is the union of the ranges. Arithmetic and elementary functions thread through intervals, producing rigorous enclosures; exact endpoints are kept exact.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/interval.c`](https://github.com/stblake/mathilda/blob/main/src/interval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)
