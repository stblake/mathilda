# SessionTime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SessionTime[] gives the total number of seconds of wall-clock time elapsed since the beginning of the current Mathilda session.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Measured from a monotonic clock captured at kernel start-up; includes time
  spent in `Pause`.

**Attributes:** `Protected`.

## References

**See also:** [Pause](../../time-and-date/Pause/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)
