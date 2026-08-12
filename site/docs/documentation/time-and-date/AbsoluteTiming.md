# AbsoluteTiming

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
AbsoluteTiming[expr] evaluates expr, and returns a list of the absolute number of seconds of elapsed wall-clock time, together with the result obtained.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, `Protected`, `SequenceHold`.
- Returns `{seconds, result}`.
- Elapsed real time from a monotonic clock, so a clock adjustment during a long

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
