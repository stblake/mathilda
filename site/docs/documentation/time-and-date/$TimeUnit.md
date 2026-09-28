# $TimeUnit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

$TimeUnit gives the minimum time interval in seconds recorded on the computer system.

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected` (read-only system constant).
- A real equal to the resolution of the `clock()`-based timers
  (`1 / CLOCKS_PER_SEC`), the granularity `Pause` documents itself against.

**Attributes:** `Protected`.

## References

**See also:** [Pause](../../time-and-date/Pause/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
